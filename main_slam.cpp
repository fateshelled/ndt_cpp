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


std::vector<std::vector<ndtcpp::point2>> load_dataset(const std::string& file_path, float max_dist) {

    std::ifstream file(file_path);
    if (!file.is_open()) {
        std::cerr << "File could not be opened." << std::endl;
        return {};
    }

    std::vector<std::vector<ndtcpp::point2>> dataset;
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
            std::vector<ndtcpp::point2> points;
            points.reserve(point_num);
            for (int i = 0; i < point_num; ++i) {
                iss >> angle >> range;
                angle *= (M_PI / 180.0f);
                const float x = range * std::cos(angle);
                const float y = range * std::sin(angle);
                if (range <= 0.0f || max_dist < range) continue;
                points.push_back({x, y});
            }
            dataset.push_back(points);
        }
    }

    return dataset;
}

inline std::tuple<ndtcpp::point2, float> to_se2(const ndtcpp::mat3x3& trans) {
    return {{trans.c, trans.f}, std::acos(trans.a)};
}
inline std::vector<ndtcpp::ndtpoint2> preprocess(
    std::vector<ndtcpp::point2>& points, float voxel_size, size_t voxel_min_count, size_t neighbor_n) {

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

    for(std::size_t i = 0; i < point_size; i++) {
        kdtree::search_knn(
            points.begin(), points.end(),
            result_points.begin(), result_distances.begin(), neighbor_n,
            downsampled[i]);
        // const auto mean = ndtcpp::compute_mean(result_points);
        const auto cov = ndtcpp::compute_covariance(result_points, downsampled[i]);

        const float u = 0.5f * ((cov.a + cov.d) + std::sqrt((cov.a - cov.d) * (cov.a - cov.d) + 4.0f * cov.b * cov.b));
        const float v = 0.5f * ((cov.a + cov.d) - std::sqrt((cov.a - cov.d) * (cov.a - cov.d) + 4.0f * cov.b * cov.b));
        const float cov_size = std::max(u, v);
        if (cov_size > 0.5f) continue;

        result.push_back({downsampled[i], cov});
    }
    return result;
}

int main(void) {
    // std::string dataset_path = "dataset/corridor.lsc";
    std::string dataset_path = "dataset/hall.lsc";

    const float max_dist = 5.0f;
    auto dataset = load_dataset(dataset_path, max_dist);

    std::vector<double> durations_scan_matching;
    std::vector<double> durations_map_matching;
    // std::vector<double> durations_mapping;
    // std::vector<double> durations_svg;

    const size_t start_index = 0;
    // const size_t N = dataset.size();
    const size_t N = std::min(static_cast<size_t>(62 + 2), dataset.size());
    const float voxel_size = 0.2f;
    // const float voxel_size = 0.3f;
    const size_t voxel_min_count = 1;
    const float map_voxel_size = 0.4f;
    const bool verbose = true;
    const size_t neighbor_n = 10;
    const bool is_gicp = true;

    // debug
    ndtcpp::writeSVGSetting setting;
    setting.voxel_size = voxel_size;

    const auto init_pose = ndtcpp::mat3x3::eye();
    ndtcpp::mat3x3 odometry = init_pose;
    std::vector<ndtcpp::mat3x3> odom_trajectory;

    auto target_points = preprocess(dataset[start_index], voxel_size, voxel_min_count, neighbor_n);

    auto map_points = preprocess(dataset[start_index], map_voxel_size, voxel_min_count, neighbor_n);
    const float keyframe_register_threshold = 0.2f;
    std::vector<ndtcpp::point2> keyframes;
    keyframes.push_back(std::get<0>(to_se2(odometry)));

    for (size_t i = start_index + 1; i < N; ++i) {

        auto source_points = preprocess(dataset[i], voxel_size, voxel_min_count, neighbor_n);

        auto trans_mat = ndtcpp::mat3x3::eye();
        // /* scan-to-scan matching */
        // {
        //     auto start_time = std::chrono::high_resolution_clock::now();

        //     {
        //         /* scan matching */
        //         if (is_gicp) {
        //             ndtcpp::gicp_scan_matching(trans_mat, source_points, target_points, verbose);
        //         } else {
        //             ndtcpp::ndt_scan_matching(trans_mat, dataset[i], target_points, verbose);
        //         }
        //     }

        //     auto end_time = std::chrono::high_resolution_clock::now();

        //     auto microsec = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() / 1e6;
        //     durations_scan_matching.push_back(microsec);
        // }

        /* scan-to-map matching */
        ndtcpp::scan_matching_result result;
        ndtcpp::mat3x3 new_odom;
        {
            auto start_time = std::chrono::high_resolution_clock::now();

            {
                // new_odom = trans_mat * odometry;
                new_odom = odometry * trans_mat;

                if (is_gicp) {
                    result = ndtcpp::gicp_scan_matching(new_odom, source_points, map_points, verbose);
                } else {
                    result = ndtcpp::ndt_scan_matching(new_odom, dataset[i], map_points, verbose);
                }
            }

            auto end_time = std::chrono::high_resolution_clock::now();

            auto microsec = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() / 1e6;
            durations_map_matching.push_back(microsec);
        }

        //debug
        {
            // std::cout << std::acos(trans_mat.a) << ", " << std::asin(-trans_mat.b) << ", " << trans_mat.c << ", " << trans_mat.f << std::endl;
            // std::cout << "|" << trans_mat.a << ", " << trans_mat.b << ", " << trans_mat.c << "|" << std::endl;
            // std::cout << "|" << trans_mat.d << ", " << trans_mat.e << ", " << trans_mat.f << "|" << std::endl;
            // std::cout << "|" << trans_mat.g << ", " << trans_mat.h << ", " << trans_mat.i << "|" << std::endl;

            // std::cout << std::acos(odometry.a) << ", " << std::asin(-odometry.b) << ", " << odometry.c << ", " << odometry.f << std::endl;
            // std::cout << "|" << odometry.a << ", " << odometry.b << ", " << odometry.c << "|" << std::endl;
            // std::cout << "|" << odometry.d << ", " << odometry.e << ", " << odometry.f << "|" << std::endl;
            // std::cout << "|" << odometry.g << ", " << odometry.h << ", " << odometry.i << "|" << std::endl;

            std::vector<ndtcpp::point2> mean;
            const ndtcpp::mat2x2 trans2x2 = {trans_mat.a, trans_mat.b, trans_mat.d, trans_mat.e};
            const ndtcpp::mat2x2 trans2x2_T = {trans2x2.a, trans2x2.c, trans2x2.b, trans2x2.d};
            for (auto& pt: source_points) {
                pt.mean = ndtcpp::transformPointCopy(trans_mat, pt.mean);
                pt.cov = trans2x2 * pt.cov * trans2x2_T;
                mean.push_back(pt.mean);
            }

            std::string output_path = "slam_output/";
            if (is_gicp) {
                output_path += "gicp_" + std::to_string(i) + ".svg";
            } else {
                output_path += "ndt_" + std::to_string(i) + ".svg";
            }
            // ndtcpp::writePointsToSVG(source_points, map_points, output_path, setting);
            // ndtcpp::writePointsToSVG(source_points, target_points, output_path, setting);
            // ndtcpp::writePointsToSVG(dataset[i], target_points, output_path, setting);
            {

                auto source = dataset[i];
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

        if (result.converged && result.error < 100.0f) {
        // if (result.converged) {

            odometry = new_odom;
            target_points = source_points;
            // const auto& [pose, rot] = to_se2(odometry);

            const ndtcpp::point2 pos = std::get<0>(to_se2(odometry));
            const float dx = pos.x - keyframes[keyframes.size() - 1].x;
            const float dy = pos.y - keyframes[keyframes.size() - 1].y;
            const float dist = std::sqrt(dx * dx + dy * dy);
            if (dist >= keyframe_register_threshold) {
                const ndtcpp::mat2x2 trans2x2 = {odometry.a, odometry.b, odometry.d, odometry.e};
                const ndtcpp::mat2x2 trans2x2_T = {trans2x2.a, trans2x2.c, trans2x2.b, trans2x2.d};
                // map_points.reserve(map_points.size() + source_points.size());
                // for (auto& pt: source_points) {
                //     auto mean = ndtcpp::transformPointCopy(odometry, pt.mean);
                //     auto cov = trans2x2 * pt.cov * trans2x2_T;
                //     map_points.push_back({mean, cov});
                // }

                std::vector<ndtcpp::point2> points(map_points.size() + source_points.size());
                for (const auto& pt: map_points) {
                    points.push_back(pt.mean);
                }
                for (const auto& pt: source_points) {
                    points.push_back(ndtcpp::transformPointCopy(odometry, pt.mean));
                }
                // map_points = preprocess(points, voxel_size, 1, neighbor_n);
                const auto tmp_map = preprocess(points, map_voxel_size, keyframes.size() / 2, neighbor_n);
                if (tmp_map.size() >= map_points.size() * 0.7) {
                    map_points = tmp_map;
                    keyframes.push_back(pos);
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
