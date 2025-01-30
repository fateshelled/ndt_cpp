#include <numeric>
#include <chrono>
#include <string>
#include <math.h>
#include <cstddef>
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>

#include "type.hpp"
#include "ndt-cpu-single.hpp"
#include "voxel_map.hpp"

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

struct InterporationParam
{
    bool enable = true;
    float point_interval = 0.03f;
    float point_far_threshold = 0.3f;
};

inline std::vector<ndtcpp::ndtpoint2> preprocess(
    std::vector<Polar2>& points, float voxel_size, size_t voxel_min_count, size_t neighbor_n,
    const InterporationParam& interporation_param = InterporationParam()
) {

    auto points_carts = Polar2::to_carts(points);

    if (interporation_param.enable) {
        std::vector<ndtcpp::point2> eq_interval_points;
        {
            const float interval_squared = interporation_param.point_interval * interporation_param.point_interval;
            ndtcpp::point2 pt0 = points_carts.front();
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
        // const auto cov = ndtcpp::compute_covariance(result_points, downsampled[i]);
        const auto cov = ndtcpp::compute_covariance_line(result_points, downsampled[i]);
        result.push_back({downsampled[i], cov});
    }
    return result;
}


int main(void) {
    // std::string dataset_path = "dataset/corridor.lsc";
    std::string dataset_path = "dataset/hall.lsc";

    const float min_dist = 0.3f;
    const float max_dist = 20.0f;
    auto dataset = load_dataset(dataset_path, min_dist, max_dist);

    std::vector<double> durations_preprocess;
    std::vector<double> durations_scan_matching;
    std::vector<double> durations_map_matching;
    std::vector<double> durations_mapping;

    const size_t start_index = 0;
    const size_t N = dataset.size();
    const float voxel_size = 0.2f;
    const size_t voxel_min_count = 1;
    const size_t neighbor_n = 10;

    const bool is_gicp = true;
    const float scan_error_threshold = 0.5f;
    const float map_register_error_threshold = 1.0f;
    const bool verbose = true;

    const auto init_pose = ndtcpp::mat3x3::eye();
    ndtcpp::mat3x3 odometry = init_pose;
    std::vector<ndtcpp::mat3x3> odom_trajectory;

    const auto target_points_raw = Polar2::to_carts(dataset[start_index]);
    auto target_points = preprocess(dataset[start_index], voxel_size, voxel_min_count, neighbor_n);

    VoxelMap map(voxel_size);
    map.set_occupied_threshold(0.7f);
    map.set_empty_threshold(0.3f);
    map.addPoints(target_points_raw, odometry);
    map.updateStatus();

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
        }

        target_points = source_points;
        const ndtcpp::mat3x3 scan2scan_odom = odometry * scan2scan_trans_mat;

        ndtcpp::scan_matching_result scan2map_result;
        const bool do_scan2map = scan2scan_result.error < scan_error_threshold;
        bool success = false;
        if (do_scan2map) {
            /* scan-to-map matching */
            ndtcpp::mat3x3 scan2map_odom = scan2scan_odom;
            {
                auto start_time = std::chrono::high_resolution_clock::now();
                if (is_gicp) {
                    scan2map_result = ndtcpp::gicp_scan_matching(scan2map_odom, source_points, map_points, verbose);
                } else {
                    scan2map_result = ndtcpp::ndt_scan_matching(scan2map_odom, source_points_raw, map_points, verbose);
                }
                auto end_time = std::chrono::high_resolution_clock::now();

                auto microsec = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() / 1e6;
                durations_map_matching.push_back(microsec);
            }

            // mapping
            success = scan2map_result.error < map_register_error_threshold;
            if (success) {
                odometry = scan2map_odom;

                auto start_time = std::chrono::high_resolution_clock::now();

                /* add raw points */
                // {
                //     map.addPoints(source_points_raw, odometry);
                // }

                /* add downsampled points */
                {
                    std::vector<ndtcpp::point2> points;
                    std::transform(
                        source_points.begin(), source_points.end(),
                        std::back_insert_iterator<std::vector<ndtcpp::point2>>(points),
                        [](const ndtcpp::ndtpoint2& pt){return pt.mean;});
                    map.addPoints(points, odometry);
                }

                const bool updated = map.updateStatus();

                // update map_points
                if (updated) {
                    auto cloud = map.to_point_cloud();
                    ndtcpp::compute_ndt_points(cloud, map_points);
                }

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
            static ndtcpp::writeSVGSetting setting;
            setting.voxel_size = voxel_size;
            setting.size = 500;

            // save map
            map.saveAsSVG("slam_output/map[" + std::to_string(i) + "].svg");

            std::string output_path = "slam_output/";
            if (is_gicp) {
                output_path += "gicp_scan2map[" + std::to_string(i) + "]_";
                if (scan2map_result.converged) output_path += "conv_";
                output_path += std::to_string(scan2map_result.error) + ".svg";
                auto source = source_points;

                const ndtcpp::mat2x2 trans2x2 = {odometry.a, odometry.b, odometry.d, odometry.e};
                const ndtcpp::mat2x2 trans2x2_T = trans2x2.transpose();
                for (auto& pt: source) {
                    pt.mean = ndtcpp::transformPointCopy(odometry, pt.mean);
                    pt.cov = trans2x2 * pt.cov * trans2x2_T;
                }
                // ndtcpp::writePointsToSVG(source, map_points, output_path, setting);
                ndtcpp::writePointsToSVG(source, map_points, odometry, scan2map_result.H, output_path, setting);
            } else {
                output_path += "ndt_scan2map[" + std::to_string(i) + "]_";
                if (scan2map_result.converged) output_path += "conv_";
                output_path += std::to_string(scan2map_result.error) + ".svg";
                auto source = source_points_raw;
                ndtcpp::transformPointsZeroCopy(odometry, source);
                ndtcpp::writePointsToSVG(source, map_points, output_path, setting);
            }
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
