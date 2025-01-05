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
                if (range <= min_dist || max_dist < range) continue;
                points.push_back({angle, range});
            }
            dataset.push_back(points);
        }
    }

    return dataset;
}

inline std::tuple<ndtcpp::point2, float> to_se2(const ndtcpp::mat3x3& trans) {
    return {{trans.c, trans.f}, std::asin(trans.c)};
}

inline std::vector<ndtcpp::ndtpoint2> preprocess(
    std::vector<ndtcpp::point2>& points, float voxel_size, size_t voxel_min_count, size_t neighbor_n, float error_ellipse_area_thredhold=0.01f) {

    // ndtcpp::compute_ndt_points(dataset[0], target_points);
    // ndtcpp::compute_ndt_points_downsampling(dataset[0], target_points, voxel_size, voxel_min_count);

    auto downsampled = std::vector<ndtcpp::point2>();
    ndtcpp::compute_voxel_downsampling(points, downsampled, voxel_size, voxel_min_count);

    auto result = std::vector<ndtcpp::ndtpoint2>();
    const auto point_size = downsampled.size();
    result.reserve(point_size);

    kdtree::construct(points.begin(), points.end());

    std::vector<ndtcpp::point2> result_points(neighbor_n);
    std::vector<float> result_distances(neighbor_n);

    const float area_thresh = error_ellipse_area_thredhold;
    for(std::size_t i = 0; i < point_size; i++) {
        kdtree::search_knn(
            points.begin(), points.end(),
            result_points.begin(), result_distances.begin(), neighbor_n,
            downsampled[i]);
        // const auto mean = ndtcpp::compute_mean(result_points);
        const auto cov = ndtcpp::compute_covariance(result_points, downsampled[i]);

        const float a = 0.5f * (cov.a + cov.d);
        const float b = 0.5f * std::sqrt((cov.a - cov.d) * (cov.a - cov.d) + 4.0f * cov.b * cov.b);
        const float u = a + b;
        const float v = a - b;
        // u * v * PI = 誤差楕円の大きさ
        if (u * v * M_PI > area_thresh) continue;

        result.push_back({downsampled[i], cov});
    }
    return result;
}

int main(void) {
    // std::string dataset_path = "dataset/corridor.lsc";
    std::string dataset_path = "dataset/hall.lsc";

    const float min_dist = 0.01f;
    const float max_dist = 20.0f;
    auto dataset = load_dataset(dataset_path, min_dist, max_dist);

    std::vector<double> durations_scan_matching;
    std::vector<double> durations_map_matching;
    // std::vector<double> durations_mapping;
    // std::vector<double> durations_svg;

    const size_t start_index = 0;
    const size_t N = dataset.size();
    // const size_t N = std::min(static_cast<size_t>(130 + 2), dataset.size());
    const float voxel_size = 0.2f;
    // const float voxel_size = 0.3f;
    const size_t voxel_min_count = 1;
    const size_t neighbor_n = 10;

    const float map_voxel_size = 0.4f;
    const size_t map_voxel_min_count = 1;
    const size_t map_neighbor_n = 6;
    const bool verbose = true;
    const bool is_gicp = true;
    const bool use_scan2scan = true;

    // debug
    ndtcpp::writeSVGSetting setting;
    setting.voxel_size = voxel_size;
    setting.size = 200;

    const auto init_pose = ndtcpp::mat3x3::eye();
    ndtcpp::mat3x3 odometry = init_pose;
    std::vector<ndtcpp::mat3x3> odom_trajectory;

    auto target_points_raw = Polar2::to_carts(dataset[start_index]);
    auto target_points = preprocess(target_points_raw, voxel_size, voxel_min_count, neighbor_n);

    auto map_points = preprocess(target_points_raw, map_voxel_size, voxel_min_count, neighbor_n);
    const float keyframe_register_threshold_dist = 0.1f;
    const float keyframe_register_threshold_angle = (10.0f) * (M_PI / 180.0f);
    std::vector<ndtcpp::point2> keyframes;
    keyframes.push_back(std::get<0>(to_se2(odometry)));

    for (size_t i = start_index + 1; i < N; ++i) {

        auto source_points_raw = Polar2::to_carts(dataset[i]);
        auto source_points = preprocess(source_points_raw, voxel_size, voxel_min_count, neighbor_n);

        ndtcpp::scan_matching_result scan2scan_result;
        auto trans_mat = ndtcpp::mat3x3::eye();
        /* scan-to-scan matching */
        if (use_scan2scan) {
            auto start_time = std::chrono::high_resolution_clock::now();

            {
                /* scan matching */
                if (is_gicp) {
                    scan2scan_result = ndtcpp::gicp_scan_matching(trans_mat, source_points, target_points, verbose);
                } else {
                    scan2scan_result = ndtcpp::ndt_scan_matching(trans_mat, source_points_raw, target_points, verbose);
                }
            }

            auto end_time = std::chrono::high_resolution_clock::now();

            auto microsec = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() / 1e6;
            durations_scan_matching.push_back(microsec);
        }

        /* scan-to-map matching */
        ndtcpp::scan_matching_result scan2map_result;
        ndtcpp::mat3x3 new_odom;
        {
            auto start_time = std::chrono::high_resolution_clock::now();

            {
                // if (use_scan2scan && scan2scan_result.converged) {
                if (use_scan2scan) {
                    new_odom = odometry * trans_mat;
                } else {
                    new_odom = odometry;
                }

                if (is_gicp) {
                    scan2map_result = ndtcpp::gicp_scan_matching(new_odom, source_points, map_points, verbose);
                } else {
                    scan2map_result = ndtcpp::ndt_scan_matching(new_odom, source_points_raw, map_points, verbose);
                }
            }

            auto end_time = std::chrono::high_resolution_clock::now();

            auto microsec = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() / 1e6;
            durations_map_matching.push_back(microsec);
        }

        //debug
        {
            // const auto scan2scan_odom = odometry * trans_mat;
            // std::cout << std::asin(-scan2scan_odom.b) << ", " << scan2scan_odom.c << ", " << scan2scan_odom.f << std::endl;
            // // std::cout << "|" << scan2scan_odom.a << ", " << scan2scan_odom.b << ", " << scan2scan_odom.c << "|" << std::endl;
            // // std::cout << "|" << scan2scan_odom.d << ", " << scan2scan_odom.e << ", " << scan2scan_odom.f << "|" << std::endl;
            // // std::cout << "|" << scan2scan_odom.g << ", " << scan2scan_odom.h << ", " << scan2scan_odom.i << "|" << std::endl;

            // const auto scan2map_odom = new_odom;
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

            std::string output_path = "slam_output/";
            if (is_gicp) {
                output_path += "gicp_" + std::to_string(i) + ".svg";
            } else {
                output_path += "ndt_" + std::to_string(i) + ".svg";
            }
            // ndtcpp::writePointsToSVG(source_transformed, map_points, output_path, setting);
            // ndtcpp::writePointsToSVG(source_transformed, target_points, output_path, setting);
            // ndtcpp::writePointsToSVG(source_points_raw, target_points, output_path, setting);
            {
                auto source = source_points_raw;
                ndtcpp::transformPointsZeroCopy(odometry, source);
                ndtcpp::writePointsToSVG(source, map_points, output_path, setting);
            }
            // {
            //     std::vector<ndtcpp::point2> map_pts;
            //     std::vector<ndtcpp::point2> empty;
            //     for (const auto& pt: map_points) {
            //         map_pts.push_back(pt.mean);
            //     }
            //     ndtcpp::writePointsToSVG(map_pts, empty, output_path);
            // }
            std::cout << output_path << std::endl;
        }

        if (scan2map_result.converged && scan2map_result.error < 100.0f) {
        // if (result.converged) {

            odometry = new_odom;
            target_points = source_points;
            const auto& [pos, rot] = to_se2(odometry);

            const float dx = pos.x - keyframes[keyframes.size() - 1].x;
            const float dy = pos.y - keyframes[keyframes.size() - 1].y;
            const float dist = std::sqrt(dx * dx + dy * dy);
            if (dist >= keyframe_register_threshold_dist || std::fabs(rot) >= keyframe_register_threshold_angle) {
                const ndtcpp::mat2x2 trans2x2 = {odometry.a, odometry.b, odometry.d, odometry.e};
                const ndtcpp::mat2x2 trans2x2_T = {trans2x2.a, trans2x2.c, trans2x2.b, trans2x2.d};

                // // 単純追加
                // map_points.reserve(map_points.size() + source_points.size());
                // for (auto& pt: source_points) {
                //     auto mean = ndtcpp::transformPointCopy(odometry, pt.mean);
                //     auto cov = trans2x2 * pt.cov * trans2x2_T;
                //     map_points.push_back({mean, cov});
                // }

                // 単純追加しつつ、ダウンサンプリング
                std::vector<ndtcpp::point2> points(map_points.size() + source_points.size());
                for (const auto& pt: map_points) {
                    points.push_back(pt.mean);
                }
                for (const auto& pt: source_points) {
                    points.push_back(ndtcpp::transformPointCopy(odometry, pt.mean));
                }
                // const auto tmp_map = preprocess(points, map_voxel_size, 1, map_neighbor_n);
                const auto tmp_map = preprocess(points, map_voxel_size, map_voxel_min_count, map_neighbor_n, 0.02f);
                // const auto tmp_map = preprocess(points, map_voxel_size, map_voxel_min_count, map_neighbor_n, 0.05f);
                // const auto tmp_map = preprocess(points, map_voxel_size, map_voxel_min_count, map_neighbor_n, 0.1f);

                // // 点群数が増えていたらマップを更新
                // if (tmp_map.size() >= map_points.size()) {
                //     std::cout << "REGISTER KEYFRAME" << std::endl;
                //     map_points = tmp_map;
                //     keyframes.push_back(pos);
                // }

                // マップを更新
                {
                    map_points = tmp_map;
                    keyframes.push_back(pos);
                    std::cout << "REGISTER KEYFRAME" << std::endl;
                }
            }
        }
    }

    {
        const double mean = std::accumulate(durations_scan_matching.begin(), durations_scan_matching.end(), 0.0) / durations_scan_matching.size();
        std::cout << "SCAN-TO-SCAN MATCHING MEAN: " << mean << " mill sec" << std::endl;
    }
    {
        const double mean = std::accumulate(durations_map_matching.begin(), durations_map_matching.end(), 0.0) / durations_map_matching.size();
        std::cout << "SCAN-TO-MAP MATCHING MEAN: " << mean << " mill sec" << std::endl;
    }

}
