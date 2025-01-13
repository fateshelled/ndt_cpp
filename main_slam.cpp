#include <numeric>
#include <chrono>
#include <string>
#include <math.h>
#include <cstddef>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include <unordered_set>

#include "type.hpp"
#include "ndt-cpu-single.hpp"

struct Polar2
{
    float angle;
    float range;

    static ndtcpp::point2 to_cart(const Polar2& polar) {
        const float x = polar.range * std::cos(polar.angle);
        const float y = polar.range * std::sin(polar.angle);
        return {x, y};
    }

    static std::vector<ndtcpp::point2> to_carts(const std::vector<Polar2>& polars) {
        std::vector<ndtcpp::point2> points;
        std::transform(polars.begin(), polars.end(),  std::back_inserter(points), [](const Polar2& p) { return to_cart(p); });
        return points;
    }
};


std::vector<std::vector<Polar2>> load_dataset(const std::string& file_path, float min_dist, float max_dist) {

    std::ifstream file(file_path);
    if (!file.is_open()) {
        std::cerr << "File could not be opened." << std::endl;
        return {};
    }

    std::vector<std::vector<Polar2>> dataset;
    std::string line_str;
    while(std::getline(file, line_str)){
        std::istringstream iss(line_str);

        std::string type_str;
        iss >> type_str;
        if (type_str == "LASERSCAN") {
            std::string tmp;
            iss >> tmp >> tmp >> tmp;
            int point_num;
            iss >> point_num;

            float angle, range;
            std::vector<Polar2> points;
            points.reserve(point_num);
            for (int i = 0; i < point_num; ++i) {
                iss >> angle >> range;
                angle *= (M_PI / 180.0f);
                if (range <= min_dist || max_dist <= range) continue;
                points.push_back({angle, range});
            }
            std::sort(points.begin(), points.end(), [](const Polar2& a, const Polar2& b) { return a.angle < b.angle; });
            dataset.push_back(points);
        }
    }

    return dataset;
}

inline std::tuple<ndtcpp::point2, float> to_se2(const ndtcpp::mat3x3& trans) {
    return {{trans.c, trans.f}, std::asin(trans.c)};
}


struct InterporationParam
{
    bool enable = true;
    float point_interval = 0.03f;
    float point_far_threshold = 0.3f;
};

inline std::vector<ndtcpp::ndtpoint2> preprocess(
    std::vector<Polar2>& points, float voxel_size, size_t voxel_min_count, size_t neighbor_n,
    const InterporationParam& interporation_param = InterporationParam()
    // , float error_ellipse_area_thredhold=1.0f
) {

    // ndtcpp::compute_ndt_points(dataset[0], target_points);
    // ndtcpp::compute_ndt_points_downsampling(dataset[0], target_points, voxel_size, voxel_min_count);
    auto points_carts = Polar2::to_carts(points);

    if (interporation_param.enable) {
        std::vector<ndtcpp::point2> eq_interval_points;
        {
            const float interval_squared = interporation_param.point_interval * interporation_param.point_interval;
            ndtcpp::point2& pt0 = points_carts.front();
            eq_interval_points.push_back(pt0);
            for (size_t i = 1; i < points_carts.size(); ++i) {
                const auto dx = pt0.x - points_carts[i].x;
                const auto dy = pt0.y - points_carts[i].y;
                const auto dist = dx * dx + dy * dy;
                if (dist < interval_squared) {
                    continue;
                }
                if (dist > interporation_param.point_far_threshold) {
                    eq_interval_points.push_back(points_carts[i]);
                    pt0 = points_carts[i];
                    continue;
                }
                const float ratio = interporation_param.point_interval / std::sqrt(dist);
                const ndtcpp::point2 new_point = {
                    ratio * dx + points_carts[i].x,
                    ratio * dy + points_carts[i].y,
                };
                eq_interval_points.push_back(new_point);

                pt0 = new_point;
            }
        }

        points_carts = eq_interval_points;
    }

    std::vector<ndtcpp::point2> downsampled;
    ndtcpp::compute_voxel_downsampling(points_carts, downsampled, voxel_size, voxel_min_count);

    // ndtcpp::writePointsToSVG(points_carts, downsampled, "slam_output/interpolate_ds.svg");

    std::vector<ndtcpp::ndtpoint2> result;
    const auto point_size = downsampled.size();
    result.reserve(point_size);

    kdtree::construct(points_carts.begin(), points_carts.end());

    std::vector<ndtcpp::point2> result_points(neighbor_n);
    std::vector<float> result_distances(neighbor_n);

    for(std::size_t i = 0; i < point_size; i++) {
        kdtree::search_knn(
            points_carts.begin(), points_carts.end(),
            result_points.begin(), result_distances.begin(), neighbor_n,
            downsampled[i]);
        // const auto mean = ndtcpp::compute_mean(result_points);
        const auto cov = ndtcpp::compute_covariance(result_points, downsampled[i]);

        // const float a = 0.5f * (cov.a + cov.d);
        // const float b = (cov.a - cov.d);
        // const float c = 0.5f * std::sqrt(b * b + 4.0f * cov.b * cov.b);
        // const float u = a + c;
        // const float v = a - c;
        // // u * v * PI = 誤差楕円の大きさ
        // if (u * v * M_PI > error_ellipse_area_thredhold) continue;

        result.push_back({downsampled[i], cov});
    }
    return result;
}

inline std::vector<ndtcpp::ndtpoint2> remove_large_covariance_points(const std::vector<ndtcpp::ndtpoint2>& points, float error_ellipse_area_thredhold=0.01f) {
    std::vector<ndtcpp::ndtpoint2> ret;
    ret.reserve(points.size());

    for (const auto& pt: points) {
        const auto& cov = pt.cov;
        const float a = 0.5f * (cov.a + cov.d);
        const float b = (cov.a - cov.d);
        const float c = 0.5f * std::sqrt(b * b + 4.0f * cov.b * cov.b);
        const float u = a + c;
        const float v = a - c;
        if (u * v * M_PI > error_ellipse_area_thredhold) continue;
        ret.push_back(pt);
    }
    return ret;
}


class VoxelMap {
private:
    static const int STATE_UNKOWN = 0;
    static const int STATE_OCCUPIED = 1;
    static const int STATE_EMPTY = 2;
    struct Voxel {
        size_t count = 0;
        ndtcpp::point2 mean = {0.f, 0.f};
        float occupancy = 0.5f;
        float probability = std::exp(0.5f) / (1.0f + std::exp(0.5f));
        int state = VoxelMap::STATE_UNKOWN;
    };
    float voxel_size_;
    float voxel_size_inv_;
    float log_odds_;
    float occupied_threshold_ = 0.8f;
    float empty_threshold_ = 0.2f;
    std::unordered_map<std::tuple<int, int>, Voxel, ndtcpp::tuple_int_hash> occupancy_;

    std::unordered_set<std::tuple<int, int>, ndtcpp::tuple_int_hash> line(int sx, int sy, int ex, int ey) {
        std::unordered_set<std::tuple<int, int>, ndtcpp::tuple_int_hash> ret;
        const int dx = std::abs(ex - sx);
        const int dy = std::abs(ey - sy);
        const int x = sx < ex ? 1 : -1;
        const int y = sy < ey ? 1 : -1;
        int error = dx - dy;
        int x0 = sx;
        int y0 = sy;
        const int x1 = ex;
        const int y1 = ey;

        while (true)
        {
            if (x0 == x1 && y0 == y1) break;

            ret.emplace(std::tuple<int, int>(x0, y0));

            const int error2 = 2 * error;
            if (error2 > -dy) {
                error -= dy;
                x0 += x;
            }
            if (error2 < dx) {
                error += dx;
                y0 += y;
            }

        }
        return ret;
    }

    static float to_probability(const Voxel& voxel) {
        const float exp = std::exp(voxel.occupancy);
        return exp / (1.0f + exp);
    }

public:
    VoxelMap(float voxel_size, float odds=0.4f, float occupied_threshold=0.8f, float empty_threshold=0.2f)
    : voxel_size_(voxel_size), voxel_size_inv_(1.f / voxel_size), log_odds_(std::log((1.0f - odds)/odds)),
      occupied_threshold_(occupied_threshold), empty_threshold_(empty_threshold) {
    }
    ~VoxelMap(){}

    void set_occupied_threshold(float value) {this->occupied_threshold_ = value;}
    void set_empty_threshold(float value) {this->empty_threshold_ = value;}

    void addPoints(const std::vector<ndtcpp::point2> points_no_trans, const ndtcpp::mat3x3& odom) {
        const auto pos = std::get<0>(to_se2(odom));
        const int pos_voxel_x = std::floor(pos.x * voxel_size_inv_);
        const int pos_voxel_y = std::floor(pos.y * voxel_size_inv_);

        std::unordered_map<std::tuple<int, int>, Voxel, ndtcpp::tuple_int_hash> new_occupied;
        std::unordered_set<std::tuple<int, int>, ndtcpp::tuple_int_hash> new_empty;
        for (const auto& pt: points_no_trans) {
            const auto transformed = ndtcpp::transformPointCopy(odom, pt);
            // occupied
            const std::tuple<int, int> index = {
                std::floor(transformed.x * voxel_size_inv_),
                std::floor(transformed.y * voxel_size_inv_)
            };
            if (new_occupied.count(index) == 0) {
                new_occupied[index] = Voxel();
            }
            new_occupied[index].mean = transformed;
            ++new_occupied[index].count;

            // empty
            const auto empty_cells = line(pos_voxel_x, pos_voxel_y, std::get<0>(index), std::get<1>(index));
            for (const auto& c: empty_cells) {
                if (new_occupied.count(c) == 0) {
                    new_empty.emplace(c);
                }
            }
        }

        for (const auto& [index, voxel]: new_occupied) {
            ++this->occupancy_[index].count;
            this->occupancy_[index].mean.x += (voxel.mean.x - this->occupancy_[index].mean.x) / this->occupancy_[index].count;
            this->occupancy_[index].mean.y += (voxel.mean.y - this->occupancy_[index].mean.y) / this->occupancy_[index].count;
            this->occupancy_[index].occupancy += this->log_odds_;
        }

        for (const auto& index: new_empty) {
            this->occupancy_[index].occupancy -= this->log_odds_;
        }
    }

    bool updateStatus() {
        bool updated = false;
        for (auto& [index, voxel]: this->occupancy_) {
            // const auto old_prob = voxel.probability;
            int old_state = voxel.state;

            voxel.probability = this->to_probability(voxel);
            if (voxel.probability >= occupied_threshold_) {
                voxel.state = STATE_OCCUPIED;
            } else if (voxel.probability <= empty_threshold_) {
                voxel.state = STATE_EMPTY;
            }

            if (old_state != voxel.state) {
                updated = true;
            }
        }
        return updated;
    }

    std::vector<ndtcpp::point2> to_point_cloud() {
        std::vector<ndtcpp::point2> ret;
        for (const auto&[index, voxel]: this->occupancy_) {
            if (voxel.state == STATE_OCCUPIED) {
                ret.push_back(voxel.mean);
            }
        }
        return ret;
    }

    size_t saveAsSVG(const std::string& file_name) {
        ndtcpp::writeSVGSetting setting;
        setting.size = 250;
        setting.point1_pt_color = "black";
        setting.point2_pt_color = "white";
        setting.flip_y = true;

        std::ofstream file(file_name);
        if (!file.is_open()) {
            std::cerr << "Cannot open file for writing." << std::endl;
            return 0;
        }
        const int size = setting.size;
        const float scale = setting.scale;
        // const float ellipse_scale = setting.ellipse_scale;
        const float offset = size / 2.0f;
        // const std::string point1_ellipse_color = setting.point1_ellipse_color;
        // const std::string point2_ellipse_color = setting.point2_ellipse_color;
        const std::string point1_pt_color = setting.point1_pt_color;
        const std::string point2_pt_color = setting.point2_pt_color;
        // const float voxel_size = setting.voxel_size;
        const std::string bg_color = "gray";
        const int sign = setting.flip_y ? -1: 1;

        file << "<svg xmlns='http://www.w3.org/2000/svg' width='" << size << "' height='" << size << "'>\n";
        file << "<rect width='" << size << "' height='" << size << "' x='0' y='0' fill='" << bg_color << "' stroke='#000' />\n";

        size_t count = 0;
        for (const auto& [index, voxel] : this->occupancy_) {
            const auto prob = this->to_probability(voxel);
            // const auto prob = voxel.probability;
            if (prob >= this->occupied_threshold_) {
                file << "<rect width='1' height='1' x='" << std::get<0>(index) + offset << "' y='" << sign * std::get<1>(index) + offset  << "' fill='" << point1_pt_color << "'/>\n";
                ++count;
            } else if (prob <= this->empty_threshold_) {
                file << "<rect width='1' height='1' x='" << std::get<0>(index) + offset << "' y='" << sign * std::get<1>(index) + offset  << "' fill='" << point2_pt_color << "'/>\n";
                ++count;
            }
        }
        file << "</svg>\n";
        file.close();
        return count;

    }
};

int main(void) {
    // std::string dataset_path = "dataset/corridor.lsc";
    std::string dataset_path = "dataset/hall.lsc";

    const float min_dist = 0.01f;
    const float max_dist = 20.0f;
    auto dataset = load_dataset(dataset_path, min_dist, max_dist);

    std::vector<double> durations_preprocess;
    std::vector<double> durations_scan_matching;
    std::vector<double> durations_map_matching;
    std::vector<double> durations_mapping;

    const size_t start_index = 0;
    const size_t warmup_num = 40;
    // const size_t N = dataset.size();
    const size_t N = std::min(static_cast<size_t>(400 + 1), dataset.size());
    const float voxel_size = 0.2f;
    // const float voxel_size = 0.3f;
    const size_t voxel_min_count = 1;
    const size_t neighbor_n = 10;

    const bool is_gicp = true;
    // const float scan_error_threshold = 0.1f;
    // const float map_register_error_threshold = 0.4f;
    const float scan_error_threshold = 0.5f;
    const float map_register_error_threshold = 1.0f;
    const bool verbose = true;

    // debug
    ndtcpp::writeSVGSetting setting;
    setting.voxel_size = voxel_size;
    setting.size = 500;

    const auto init_pose = ndtcpp::mat3x3::eye();
    ndtcpp::mat3x3 odometry = init_pose;
    std::vector<ndtcpp::mat3x3> odom_trajectory;

    const auto target_points_raw = Polar2::to_carts(dataset[start_index]);
    auto target_points = preprocess(dataset[start_index], voxel_size, voxel_min_count, neighbor_n);

    VoxelMap map(voxel_size);
    // map.set_occupied_threshold(0.75f);
    // map.set_empty_threshold(0.25f);
    map.set_occupied_threshold(0.7f);
    map.set_empty_threshold(0.3f);
    map.addPoints(target_points_raw, odometry);
    map.updateStatus();
    ndtcpp::GICP_PARAMS map_gicp_param;
    map_gicp_param.max_iter_num = 30;
    // map_gicp_param.max_correspondence_distance *= 0.5f;
    // map_gicp_param.min_correspondence = 20;

    std::vector<ndtcpp::point2> simple_map_points; // debug
    std::vector<ndtcpp::ndtpoint2> map_points;
    map_points = target_points;


    for (size_t i = start_index + 1; i < N; ++i) {
        std::cout << "[" << i << "]" << std::endl;

        std::vector<ndtcpp::point2> source_points_raw;
        std::vector<ndtcpp::ndtpoint2> source_points;
        {
            auto start_time = std::chrono::high_resolution_clock::now();

            source_points_raw = Polar2::to_carts(dataset[i]);
            source_points = preprocess(dataset[i], voxel_size, voxel_min_count, neighbor_n);

            auto end_time = std::chrono::high_resolution_clock::now();
            auto microsec = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() / 1e6;
            durations_preprocess.push_back(microsec);
        }

        ndtcpp::scan_matching_result scan2scan_result;
        auto scan2scan_trans_mat = ndtcpp::mat3x3::eye();
        /* scan-to-scan matching */
        {
            auto start_time = std::chrono::high_resolution_clock::now();

            /* scan matching */
            if (is_gicp) {
                scan2scan_result = ndtcpp::gicp_scan_matching(scan2scan_trans_mat, source_points, target_points, verbose);
            } else {
                scan2scan_result = ndtcpp::ndt_scan_matching(scan2scan_trans_mat, source_points_raw, target_points, verbose);
            }

            auto end_time = std::chrono::high_resolution_clock::now();

            auto microsec = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() / 1e6;
            durations_scan_matching.push_back(microsec);

            // // debug
            // {
            //     std::string output_path = "slam_output/";
            //     if (is_gicp) {
            //         output_path += "gicp_scan2scan[" + std::to_string(i) + "]_";
            //         if (scan2scan_result.converged) output_path += "conv_";
            //         output_path += std::to_string(scan2scan_result.error);
            //         {
            //             auto source = source_points;
            //             for (auto& pt: source) {
            //                 pt.mean = ndtcpp::transformPointCopy(scan2scan_trans_mat, pt.mean);
            //             }
            //             ndtcpp::writePointsToSVG(source, target_points, output_path + ".svg", setting);
            //         }
            //         {
            //             auto source = source_points;
            //             for (auto& pt: source) {
            //                 pt.mean = ndtcpp::transformPointCopy(odometry * scan2scan_trans_mat, pt.mean);
            //             }
            //             ndtcpp::writePointsToSVG(source, map_points, output_path + "_map.svg", setting);
            //         }
            //     } else {
            //         output_path += "ndt_scan2scan[" + std::to_string(i) + "]_";
            //         if (scan2scan_result.converged) output_path += "conv_";
            //         output_path += std::to_string(scan2scan_result.error) + ".svg";
            //         {
            //             auto source = source_points_raw;
            //             ndtcpp::transformPointsZeroCopy(scan2scan_trans_mat, source);
            //             ndtcpp::writePointsToSVG(source, target_points, output_path + ".svg", setting);
            //         }
            //         {
            //             auto source = source_points_raw;
            //             ndtcpp::transformPointsZeroCopy(odometry * scan2scan_trans_mat, source);
            //             ndtcpp::writePointsToSVG(source, map_points, output_path + "_map.svg", setting);
            //         }
            //     }
            // }
        }

        target_points = source_points;
        const ndtcpp::mat3x3 scan2scan_odom = odometry * scan2scan_trans_mat;

        auto map_update_function = [&map, &map_points](const std::vector<ndtcpp::point2>& points, const ndtcpp::mat3x3& odom, bool update_always){
            map.addPoints(points, odom);
            bool updated = map.updateStatus();

            // update map_points
            if (update_always) updated = true;
            if (updated) {
                auto cloud = map.to_point_cloud();
                ndtcpp::compute_ndt_points(cloud, map_points);
                // map_points = remove_large_covariance_points(map_points);
                // std::cout << "REGISTER MAP: " << cloud.size() << std::endl;
            }
            return updated;
        };

        bool map_update = false;
        if (i < start_index + warmup_num) {
        // if (i < N) {
            odometry = scan2scan_odom;

            auto start_time = std::chrono::high_resolution_clock::now();

            // std::vector<ndtcpp::point2> points;
            // std::transform(
            //     source_points.begin(), source_points.end(),
            //     std::back_insert_iterator<std::vector<ndtcpp::point2>>(points),
            //     [](const ndtcpp::ndtpoint2& pt){return pt.mean;});
            // map_update = map_update_function(points, scan2scan_odom, true);
            map_update = map_update_function(source_points_raw, scan2scan_odom, true);

            auto end_time = std::chrono::high_resolution_clock::now();
            auto microsec = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() / 1e6;
            durations_mapping.push_back(microsec);

        } else if (scan2scan_result.error < scan_error_threshold) {

            /* scan-to-map matching */
            ndtcpp::scan_matching_result scan2map_result;
            ndtcpp::mat3x3 scan2map_odom = scan2scan_odom;
            // ndtcpp::mat3x3 scan2map_odom = odometry * scan2scan_trans_mat;
            {
                auto start_time = std::chrono::high_resolution_clock::now();
                if (is_gicp) {
                    scan2map_result = ndtcpp::gicp_scan_matching(scan2map_odom, source_points, map_points, verbose, map_gicp_param);
                } else {
                    scan2map_result = ndtcpp::ndt_scan_matching(scan2map_odom, source_points_raw, map_points, verbose);
                }
                auto end_time = std::chrono::high_resolution_clock::now();

                auto microsec = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() / 1e6;
                durations_map_matching.push_back(microsec);

                // debug
                {
                    std::string output_path = "slam_output/";
                    if (is_gicp) {
                        output_path += "gicp_scan2map[" + std::to_string(i) + "]_";
                        if (scan2map_result.converged) output_path += "conv_";
                        output_path += std::to_string(scan2map_result.error) + ".svg";
                        auto source = source_points;
                        for (auto& pt: source) {
                            pt.mean = ndtcpp::transformPointCopy(scan2map_odom, pt.mean);
                        }
                        ndtcpp::writePointsToSVG(source, map_points, output_path, setting);
                    } else {
                        output_path += "ndt_scan2map[" + std::to_string(i) + "]_";
                        if (scan2map_result.converged) output_path += "conv_";
                        output_path += std::to_string(scan2map_result.error) + ".svg";
                        auto source = source_points_raw;
                        ndtcpp::transformPointsZeroCopy(scan2map_odom, source);
                        ndtcpp::writePointsToSVG(source, map_points, output_path, setting);
                    }
                }
            }

            // mapping
            // if (scan2map_result.error < map_register_error_threshold || (scan2map_result.converged && scan2map_result.error < map_register_error_threshold)) {
            if (scan2map_result.error < map_register_error_threshold) {
                odometry = scan2map_odom; // ?

                auto start_time = std::chrono::high_resolution_clock::now();

                // std::vector<ndtcpp::point2> points;
                // std::transform(
                //     source_points.begin(), source_points.end(),
                //     std::back_insert_iterator<std::vector<ndtcpp::point2>>(points),
                //     [](const ndtcpp::ndtpoint2& pt){return pt.mean;});
                // map_update = map_update_function(points, odometry, false);
                map_update = map_update_function(source_points_raw, odometry, false);

                auto end_time = std::chrono::high_resolution_clock::now();
                auto microsec = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() / 1e6;
                durations_mapping.push_back(microsec);
            }
            else {
                odometry = scan2scan_odom;
            }
        }


        //debug
        {
            // std::vector<ndtcpp::point2> cur_points;
            // for (const auto& pt: source_points_raw) {
            //     const auto transformed = ndtcpp::transformPointCopy(odometry, pt);
            //     cur_points.push_back(transformed);
            //     simple_map_points.push_back(transformed);
            // }
            // ndtcpp::writePointsToSVG(simple_map_points, cur_points, "slam_output/simple_map_[" + std::to_string(i) + "].svg", setting);
            map.saveAsSVG("slam_output/map_only[" + std::to_string(i) + "].svg");

            // if (map_update) {
            //     const auto map_count = map.saveAsSVG("slam_output/map_" + std::to_string(i) + ".svg");
            //     std::cout << " map[" << i << "]: " << map_count << std::endl;
            // }

            // std::cout << std::asin(-scan2scan_odom.b) << ", " << scan2scan_odom.c << ", " << scan2scan_odom.f << std::endl;
            // // std::cout << "|" << scan2scan_odom.a << ", " << scan2scan_odom.b << ", " << scan2scan_odom.c << "|" << std::endl;
            // // std::cout << "|" << scan2scan_odom.d << ", " << scan2scan_odom.e << ", " << scan2scan_odom.f << "|" << std::endl;
            // // std::cout << "|" << scan2scan_odom.g << ", " << scan2scan_odom.h << ", " << scan2scan_odom.i << "|" << std::endl;

            // std::cout << std::asin(-scan2map_odom.b) << ", " << scan2map_odom.c << ", " << scan2map_odom.f << std::endl;
            // // std::cout << "|" << scan2map_odom.a << ", " << scan2map_odom.b << ", " << scan2map_odom.c << "|" << std::endl;
            // // std::cout << "|" << scan2map_odom.d << ", " << scan2map_odom.e << ", " << scan2map_odom.f << "|" << std::endl;
            // // std::cout << "|" << scan2map_odom.g << ", " << scan2map_odom.h << ", " << scan2map_odom.i << "|" << std::endl;

            // std::vector<ndtcpp::point2> mean;
            // const ndtcpp::mat2x2 trans2x2 = {trans_mat.a, trans_mat.b, trans_mat.d, trans_mat.e};
            // const ndtcpp::mat2x2 trans2x2_T = {trans2x2.a, trans2x2.c, trans2x2.b, trans2x2.d};
            // std::vector<ndtcpp::ndtpoint2> source_transformed;
            // for (const auto& pt: source_points) {
            //     const auto p = ndtcpp::transformPointCopy(trans_mat, pt.mean);
            //     source_transformed.push_back({
            //         p, trans2x2 * pt.cov * trans2x2_T
            //     });
            //     mean.push_back(p);
            // }
        }
    }

    {
        const double mean = std::accumulate(durations_preprocess.begin(), durations_preprocess.end(), 0.0) / durations_preprocess.size();
        std::cout << "PREPROCESS MEAN: " << mean << " mill sec" << std::endl;
    }
    {
        const double mean = std::accumulate(durations_scan_matching.begin(), durations_scan_matching.end(), 0.0) / durations_scan_matching.size();
        std::cout << "SCAN-TO-SCAN MATCHING MEAN: " << mean << " mill sec" << std::endl;
    }
    {
        const double mean = std::accumulate(durations_map_matching.begin(), durations_map_matching.end(), 0.0) / durations_map_matching.size();
        std::cout << "SCAN-TO-MAP MATCHING MEAN: " << mean << " mill sec" << std::endl;
    }
    {
        const double mean = std::accumulate(durations_mapping.begin(), durations_mapping.end(), 0.0) / durations_mapping.size();
        std::cout << "MAPPING MEAN: " << mean << " mill sec" << std::endl;
    }

}
