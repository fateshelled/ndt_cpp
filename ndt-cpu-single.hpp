// Copyright 2024 Taiga Takano
//
// Licensed under the Apache License, Version 2.0 (the "License");
// you may not use this file except in compliance with the License.
// You may obtain a copy of the License at
//
//     http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing, software
// distributed under the License is distributed on an "AS IS" BASIS,
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
// See the License for the specific language governing permissions and
// limitations under the License.

#ifndef NDTCPP_NDT_CPU_SINGLE_HPP_
#define NDTCPP_NDT_CPU_SINGLE_HPP_

#include <cstddef>
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <unordered_map>
#include <cmath>
#include <limits>
#include <numeric>

#include "flatkdtree.h"
#include "type.hpp"
#include "ndtcpputil.hpp"
#include "matrixutil.hpp"
#include "util.hpp"


namespace ndtcpp {
struct ndtpoint2 {
    ndtcpp::point2 mean;
    ndtcpp::mat2x2 cov;
};

struct scan_matching_result {
    bool converged = false;
    float error = std::numeric_limits<float>::max();
    size_t iter = 0;
    size_t correspondence_num = 0;
    ndtcpp::mat3x3 H = ndtcpp::mat3x3::zeros();
    ndtcpp::point3 b = ndtcpp::point3::zeros();
};

} // namespace ndtcpp

template <std::size_t I>
struct kdtree::trait::access<ndtcpp::point2, I> {
    static auto get(const ndtcpp::point2 &p) -> float
    {
        return I == 0 ? p.x : p.y;
    }
};

template <>
struct kdtree::trait::dimension<ndtcpp::point2> {
    static constexpr std::size_t value = 2;
};

template <std::size_t I>
struct kdtree::trait::access<ndtcpp::ndtpoint2, I> {
    static auto get(const ndtcpp::ndtpoint2 &p) -> float
    {
        return I == 0 ? p.mean.x : p.mean.y;
    }
};

template <>
struct kdtree::trait::dimension<ndtcpp::ndtpoint2> {
    static constexpr std::size_t value = 2;
};

namespace ndtcpp {

inline std::vector<point2> extract_points(const std::vector<ndtpoint2> &points) {
    std::vector<point2> result;
    std::transform(points.begin(), points.end(), std::back_inserter(result), [](const ndtpoint2 &pt) { return pt.mean; });
    return result;
}

inline ndtcpp::mat3x3 makeTransformationMatrix(const float& tx, const float& ty, const float& theta) {
    ndtcpp::mat3x3 mat = {
        cosf(theta), sinf(theta) * -1.0f, tx,
        sinf(theta), cosf(theta)        , ty,
        0.0f, 0.0f, 1.0f
    };
    return mat;
}

inline void transformPointsZeroCopy(const ndtcpp::mat3x3& mat, std::vector<ndtcpp::point2>& points) {
    ndtcpp::point2 transformedPoint;

    for (auto& point : points) {
        transformedPoint.x = mat.a * point.x + mat.b * point.y + mat.c;
        transformedPoint.y = mat.d * point.x + mat.e * point.y + mat.f;
        point.x = transformedPoint.x;
        point.y = transformedPoint.y;
    }
}

inline float multiplyPowPoint3(const ndtcpp::point3& vec){
    return vec.x * vec.x + vec.y * vec.y + vec.z * vec.z;
}

inline std::vector<float> compute_eigen_values(const ndtcpp::mat2x2& m) {
    const auto a_d = m.a + m.d;
    const auto root = a_d * a_d - 4.0f * (m.a * m.d - m.b * m.c);
    if (root < 0) {
        return {};
    }
    const auto val0 = 0.5f * (a_d - std::sqrt(root));
    const auto val1 = 0.5f * (a_d + std::sqrt(root));
    if (std::abs(val0 - val1) < 1e-5f) {
        return {val0};
    }
    if (val0 < val1) {
        return {val0, val1};
    }
    return {val1, val0};
}

inline ndtcpp::point2 compute_eigen_vector(const ndtcpp::mat2x2& m, const float eigen_value) {
    const auto M = m - eigen_value * ndtcpp::mat2x2::eye();

    // 逆べき乗法
    constexpr size_t k = 5;
    ndtcpp::point2 eigen_vector = {2.f, 1.f};
    for (size_t i = 0; i < k; ++i) {
        eigen_vector = M * eigen_vector;
    }
    eigen_vector *= (1.0f / eigen_vector.norm());
    return eigen_vector;
}

inline std::vector<std::tuple<float, ndtcpp::point2>> compute_eigen(const ndtcpp::mat2x2& m) {
    const auto eigen_values = compute_eigen_values(m);
    std::vector<std::tuple<float, ndtcpp::point2>> result;
    for (const auto& val: eigen_values) {
        const auto eigen_vector = compute_eigen_vector(m, val);
        result.push_back({val, eigen_vector});
    }
    return result;
}

inline ndtcpp::point3 solve3x3(const ndtcpp::mat3x3& m, const ndtcpp::point3& p) {
    float A[3][4] = {
        {m.a, m.b, m.c, p.x},
        {m.d, m.e, m.f, p.y},
        {m.g, m.h, m.i, p.z}
    };

    constexpr int n = 3;

    for (int i = 0; i < n; i++) {
        // Pivot選択
        float maxEl = std::abs(A[i][i]);
        int maxRow = i;
        for (int k = i+1; k < n; k++) {
            const auto el = std::abs(A[k][i]);
            if (el > maxEl) {
                maxEl = el;
                maxRow = k;
            }
        }

        // Pivotのある行を交換
        for (int k = i; k < n+1;k++) {
            float tmp = A[maxRow][k];
            A[maxRow][k] = A[i][k];
            A[i][k] = tmp;
        }

        // すべての行について消去を行う
        for (int k = i+1; k < n; k++) {
            const float c = -A[k][i] / A[i][i];
            for (int j = i; j < n+1; j++) {
                if (i == j) {
                    A[k][j] = 0.0f;
                } else {
                    A[k][j] += c * A[i][j];
                }
            }
        }
    }

    // 解の計算 (後退代入)
    ndtcpp::point3 solution;
    solution.z = A[2][3] / A[2][2];
    solution.y = (A[1][3] - A[1][2] * solution.z) / A[1][1];
    solution.x = (A[0][3] - A[0][2] * solution.z - A[0][1] * solution.y) / A[0][0];

    return solution;
}

inline ndtcpp::point2 solve2x2(const ndtcpp::mat2x2& m, const ndtcpp::point2& p) {
    float A[2][3] = {
        {m.a, m.b, p.x},
        {m.c, m.d, p.y},
    };

    constexpr int n = 2;

    for (int i = 0; i < n; i++) {
        // Pivot選択
        float maxEl = std::abs(A[i][i]);
        int maxRow = i;
        for (int k = i+1; k < n; k++) {
            const auto el = std::abs(A[k][i]);
            if (el > maxEl) {
                maxEl = el;
                maxRow = k;
            }
        }

        // Pivotのある行を交換
        for (int k = i; k < n+1;k++) {
            float tmp = A[maxRow][k];
            A[maxRow][k] = A[i][k];
            A[i][k] = tmp;
        }

        // すべての行について消去を行う
        for (int k = i+1; k < n; k++) {
            const float c = -A[k][i] / A[i][i];
            for (int j = i; j < n+1; j++) {
                if (i == j) {
                    A[k][j] = 0.0f;
                } else {
                    A[k][j] += c * A[i][j];
                }
            }
        }
    }

    // 解の計算 (後退代入)
    ndtcpp::point2 solution;
    solution.y = A[1][2] / A[1][1];
    solution.x = (A[0][2] - A[0][1] * solution.y) / A[0][0];

    return solution;
}

inline ndtcpp::mat3x3 expmap(const ndtcpp::point3& point){
    auto t = point.z;
    auto c = cosf(t);
    auto s = sinf(t);

    ndtcpp::mat2x2 R {
        c, s * -1.0f,
        s, c
    };

    ndtcpp::mat3x3 T {
        R.a, R.b, point.x,
        R.c, R.d, point.y,
        0.0f, 0.0f, 1.0f
    };

    return T;
}

inline ndtcpp::point2 transformPointCopy(const ndtcpp::mat3x3& mat, const ndtcpp::point2& point) {
    ndtcpp::point2 transformedPoint;

    transformedPoint.x = mat.a * point.x + mat.b * point.y + mat.c;
    transformedPoint.y = mat.d * point.x + mat.e * point.y + mat.f;

    return transformedPoint;
}

inline ndtcpp::point2 skewd(const ndtcpp::point2& input_point){
    const ndtcpp::point2 skewd_point {
        input_point.y,
        input_point.x * -1.0f
    };
    return skewd_point;
}

inline ndtcpp::point2 compute_mean(const std::vector<ndtcpp::point2>& points){
    ndtcpp::point2 mean;
    mean.x = 0.0f;
    mean.y = 0.0f;
    for(const auto& point : points){
        mean.x += point.x;
        mean.y += point.y;
    }
    mean.x = mean.x / (float)points.size();
    mean.y = mean.y / (float)points.size();
    return mean;
}

inline ndtcpp::mat2x2 compute_covariance(const std::vector<ndtcpp::point2>& points, const ndtcpp::point2& mean){
    auto point_size = points.size();
    auto vxx = 0.0f;
    auto vxy = 0.0f;
    auto vyy = 0.0f;

    for(const auto& point : points){
        const auto dx = point.x - mean.x;
        const auto dy = point.y - mean.y;
        vxx += dx * dx;
        vxy += dx * dy;
        vyy += dy * dy;
    }

    ndtcpp::mat2x2 cov;
    cov.a = vxx / point_size;
    cov.b = vxy / point_size;
    cov.c = cov.b;
    cov.d = vyy / point_size;
    return cov;
}

inline ndtcpp::mat2x2 update_covariance(const ndtcpp::mat2x2& covariance, const ndtcpp::point2& val){

    auto ret = compute_eigen(covariance);
    auto eig_vec0 = std::get<1>(ret[0]);
    auto eig_vec1 = std::get<1>(ret[1]);
    ndtcpp::mat2x2 mat;
    mat.a = eig_vec0.x;
    mat.b = eig_vec1.x;
    mat.c = eig_vec0.y;
    mat.d = eig_vec1.y;

    auto vals = ndtcpp::mat2x2::diagonal(val.x, val.y);
    return mat * vals * mat.transpose();
}

inline void update_covariance_line(ndtcpp::ndtpoint2& point){
    point.cov = update_covariance(point.cov, {1.0f, 0.1f});
}

inline void update_covariances_line(std::vector<ndtcpp::ndtpoint2>& points){
    for(auto& point : points){
        update_covariance_line(point);
    }
}

inline void compute_ndt_points(std::vector<ndtcpp::point2>& points, std::vector<ndtpoint2> &results){
    auto N = 10;

    const auto point_size = points.size();

    kdtree::construct(points.begin(), points.end());
    std::vector<ndtcpp::point2> result_points(N);
    std::vector<float> result_distances(N);

    std::vector<ndtcpp::mat2x2> covs(point_size);
    results.resize(point_size);

    for(std::size_t i = 0; i < point_size; i++) {
        kdtree::search_knn(points.begin(), points.end(), result_points.begin(), result_distances.begin(), N, points[i]);
        const auto mean = compute_mean(result_points);
        const auto cov = compute_covariance(result_points, mean);
        results[i] = {mean, cov};
    }
}

inline void compute_voxel_downsampling(
    const std::vector<ndtcpp::point2>& points, std::vector<ndtcpp::point2> &results,
    float voxel_size = 1.0f, std::size_t voxel_min_count = 4) {

    const auto point_size = points.size();

    std::unordered_map<std::tuple<int, int>, std::vector<size_t>, tuple_int_hash> voxel_indices;
    const float voxel_size_inv = 1.0f / voxel_size;

    for(size_t i = 0; i < point_size; i++) {
        const auto& pt0 = points[i];
        const std::tuple<int, int> voxel = {
            std::floor(pt0.x * voxel_size_inv),
            std::floor(pt0.y * voxel_size_inv)
        };

        voxel_indices[voxel].push_back(i);
    }

    results.clear();
    results.reserve(voxel_indices.size());
    for (const auto& [voxel, indices]: voxel_indices) {
        if (indices.size() < voxel_min_count) continue;

        std::vector<ndtcpp::point2> result_points;
        result_points.reserve(indices.size());
        for (const auto& i: indices) {
            result_points.push_back(points[i]);
        }
        const auto mean = compute_mean(result_points);
        results.push_back(mean);
    }
}

inline void compute_ndt_points_downsampling(
    const std::vector<ndtcpp::point2>& points, std::vector<ndtpoint2> &results,
    float voxel_size = 1.0f, std::size_t voxel_min_count = 4) {

    const auto point_size = points.size();

    std::unordered_map<std::tuple<int, int>, std::vector<size_t>, tuple_int_hash> voxel_indices;
    const float voxel_size_inv = 1.0f / voxel_size;

    for(size_t i = 0; i < point_size; i++) {
        const auto& pt0 = points[i];
        const std::tuple<int, int> voxel = {
            std::floor(pt0.x * voxel_size_inv),
            std::floor(pt0.y * voxel_size_inv)
        };

        voxel_indices[voxel].push_back(i);
    }

    results.clear();
    results.reserve(voxel_indices.size());
    for (const auto& [voxel, indices]: voxel_indices) {
        if (indices.size() < voxel_min_count) continue;

        std::vector<ndtcpp::point2> result_points;
        result_points.reserve(indices.size());
        for (const auto& i: indices) {
            result_points.push_back(points[i]);
        }
        const auto mean = compute_mean(result_points);
        const auto cov = compute_covariance(result_points, mean);
        results.push_back({mean, cov});
    }
}

inline scan_matching_result ndt_scan_matching(
    ndtcpp::mat3x3& trans_mat,
    const std::vector<ndtcpp::point2>& source_points,
    std::vector<ndtpoint2>& target_points, bool verbose = false
) {
    const size_t max_iter_num = 20;
    const float max_correspondence_distance = 3.0f;
    const float max_distance2 = max_correspondence_distance * max_correspondence_distance;
    const size_t point_step = 10;

    const size_t target_points_size = target_points.size();
    const size_t source_points_size = source_points.size();

    bool is_converged = false;
    ndtcpp::point3 prev_delta;
    float min_error = std::numeric_limits<float>::max();
    ndtcpp::mat3x3 min_trans_mat;

    kdtree::construct(target_points.begin(), target_points.end());
    size_t iter = 0;
    for(iter = 0; iter < max_iter_num; iter++){
        auto H_Mat = ndtcpp::mat3x3::zeros();
        auto b_Point = ndtcpp::point3::zeros();

        for(auto point_iter = 0; point_iter < source_points_size; point_iter += point_step){
            ndtpoint2 query_point = {transformPointCopy(trans_mat, source_points[point_iter]), {}};
            ndtpoint2 target_point;
            float target_distance;
            kdtree::search_knn(target_points.begin(), target_points.end(), &target_point, &target_distance, 1, query_point);

            if(target_distance > max_distance2){continue;}

            const auto identity_plus_cov = ndtcpp::mat3x3{
                target_point.cov.a + 1.0f, target_point.cov.b + 1.0f, 0.0f,
                target_point.cov.c + 1.0f, target_point.cov.d + 1.0f, 0.0f,
                0.0f, 0.0f, 1.0f
            };

            const ndtcpp::mat3x3 target_cov_inv = identity_plus_cov.inv(); //IM


            const auto error = ndtcpp::point3{
                target_point.mean.x - query_point.mean.x,
                target_point.mean.y - query_point.mean.y,
                0.0f
            };

            const ndtcpp::point2 v_point = transformPointCopy(trans_mat, skewd(source_points[point_iter]));

            const auto mat_J = ndtcpp::mat3x3{
                trans_mat.a * -1.0f, trans_mat.b * -1.0f, v_point.x,
                trans_mat.d * -1.0f, trans_mat.e * -1.0f, v_point.y,
                trans_mat.g * -1.0f, trans_mat.h * -1.0f, trans_mat.i * -1.0f
            };

            const ndtcpp::mat3x3 mat_J_T = mat_J.transpose();

            H_Mat += (mat_J_T * (target_cov_inv * mat_J));

            b_Point += (mat_J_T * (target_cov_inv * error));

        }
        b_Point *= -1.0f;

        const ndtcpp::point3 delta = solve3x3(H_Mat + 1e-6f * ndtcpp::mat3x3::eye(), b_Point);
        trans_mat = trans_mat * expmap(delta);

        const float error = multiplyPowPoint3(delta);
        if(error < 1e-4){
            is_converged = true;
        }

        if (iter > 0) {
            const float dx = prev_delta.x - delta.x;
            const float dy = prev_delta.y - delta.y;
            const float dz = prev_delta.z - delta.z;
            const auto d = std::max(std::max(std::fabs(dx), std::fabs(dy)), std::fabs(dz));
            if (d < 1e-4) {
                is_converged = true;
            }
        }

        if (is_converged) {
            min_error = error;
            if (verbose) {
                std::cout << "END NDT. ITER: " << iter;
                std::cout << ", ERROR VALUE: " << error << std::endl;
            }
            break;
        }

        prev_delta = delta;

        if (min_error > error) {
            min_error = error;
            min_trans_mat = trans_mat;
        }

        if (iter == max_iter_num - 1) {
            if (verbose) {
                std::cout << "END NDT NOT CONVERGED. ERROR VALUE: " << min_error << std::endl;
            }
            trans_mat = min_trans_mat;
        }
    }
    return {is_converged, min_error};
}

namespace {
inline float calc_gicp_error(
    const ndtcpp::point2& trans_source, const ndtcpp::point2& target, const ndtcpp::mat3x3& IM)
{
    const auto residual = ndtcpp::point3{
        target.x - trans_source.x,
        target.y - trans_source.y,
        0.0f
    };
    const ndtcpp::point3 IM_residual = IM * residual;
    return 0.5f * (residual.x * IM_residual.x + residual.y * IM_residual.y);
}
}

struct GICP_PARAMS {
    size_t max_iter_num = 20;
    float max_correspondence_distance = 0.5f;
    // float max_correspondence_distance = 3.0f;
    size_t point_step = 1;
    size_t min_correspondence = 10;
    float converged_error_th = 1e-4f;
    float converged_delta_xy_th = 1e-4f;
    float converged_delta_rot_th = 1e-6f;

    // for Levenberg-Marquardt
    // if max_inner_iter_num == 1 and lambda_factor == 1.0f -> Gauss-Newton
    size_t max_inner_iter_num = 10;
    float init_lambda = 1e-6f;
    float lambda_factor = 10.0f;
};

inline scan_matching_result gicp_scan_matching(
    ndtcpp::mat3x3& trans_mat,
    const std::vector<ndtpoint2>& source_points,
    std::vector<ndtpoint2>& target_points, bool verbose = false,
    const GICP_PARAMS& param = GICP_PARAMS()
) {
    auto max_distance = param.max_correspondence_distance;

    float lambda = param.init_lambda;

    ndtcpp::point3 prev_delta;
    scan_matching_result result = {};
    ndtcpp::mat3x3 min_trans_mat = trans_mat;

    const size_t target_points_size = target_points.size();
    const size_t source_points_size = source_points.size();

    kdtree::construct(target_points.begin(), target_points.end());
    size_t iter = 0;
    for(iter = 0; iter < param.max_iter_num; ++iter){
        const float max_distance2 = max_distance * max_distance;
        auto H_Mat = ndtcpp::mat3x3::zeros();
        auto b_Point = ndtcpp::point3::zeros();
        float error = 0.0f;

        const ndtcpp::mat2x2 trans_mat2x2 = {trans_mat.a, trans_mat.b, trans_mat.d, trans_mat.e};
        const auto trans_mat2x2_T = trans_mat2x2.transpose();

        std::vector<std::tuple<ndtcpp::mat3x3, ndtcpp::point2, int>> IMs;

        for(size_t point_iter = 0; point_iter < source_points_size; point_iter += param.point_step){
            ndtpoint2 query_point = {
                transformPointCopy(trans_mat, source_points[point_iter].mean),
                {}
            };
            ndtpoint2 target_point;
            float target_sq_distance;
            kdtree::search_knn(target_points.begin(), target_points.end(), &target_point, &target_sq_distance, 1, query_point);

            if(target_sq_distance > max_distance2){continue;}

            const auto identity_plus_target_cov = ndtcpp::mat3x3{
                target_point.cov.a + 1.0f, target_point.cov.b       , 0.0f,
                target_point.cov.c       , target_point.cov.d + 1.0f, 0.0f,
                0.0f, 0.0f, 1.0f
            };

            query_point.cov = trans_mat2x2 * source_points[point_iter].cov * trans_mat2x2_T;
            const auto identity_plus_query_cov = ndtcpp::mat3x3{
                query_point.cov.a + 1.0f, query_point.cov.b       , 0.0f,
                query_point.cov.c       , query_point.cov.d + 1.0f, 0.0f,
                0.0f, 0.0f, 1.0f
            };

            // Information Matrix
            const ndtcpp::mat3x3 IM = (identity_plus_target_cov + identity_plus_query_cov).inv();

            const auto residual = ndtcpp::point3{
                target_point.mean.x - query_point.mean.x,
                target_point.mean.y - query_point.mean.y,
                0.0f
            };

            const ndtcpp::point2 v_point = transformPointCopy(trans_mat, skewd(source_points[point_iter].mean));

            const auto mat_J = ndtcpp::mat3x3{
                trans_mat.a * -1.0f, trans_mat.b * -1.0f, v_point.x,
                trans_mat.d * -1.0f, trans_mat.e * -1.0f, v_point.y,
                trans_mat.g * -1.0f, trans_mat.h * -1.0f, trans_mat.i * -1.0f
            };

            const ndtcpp::mat3x3 mat_J_TxIM = mat_J.transpose() * IM;

            H_Mat += (mat_J_TxIM * mat_J);      // J.T * IM * J
            b_Point += (mat_J_TxIM * residual); // J.T * IM * residual

            error += calc_gicp_error(query_point.mean, target_point.mean, IM);
            IMs.emplace_back(IM, target_point.mean, point_iter);
        }
        b_Point *= -1.0f;
        if (IMs.size() > 0) {
            error /= IMs.size();
        }

        if (IMs.size() < param.min_correspondence) {
            result.error = error;
            result.correspondence_num = IMs.size();
            result.H = H_Mat;
            result.b = b_Point;
            break;
        }

        ndtcpp::point3 delta;
        ndtcpp::point3 prev_delta_inner;
        for (size_t inner_iter = 0; inner_iter < param.max_inner_iter_num; ++inner_iter) {
            // damping
            const auto H = H_Mat + lambda * ndtcpp::mat3x3::eye();

            delta = solve3x3(H, b_Point);
            trans_mat = trans_mat * expmap(delta);

            float new_error = 0.0f;
            for (const auto& [IM, target, point_iter]: IMs) {
                const auto trans_source = transformPointCopy(trans_mat, source_points[point_iter].mean);
                new_error += calc_gicp_error(trans_source, target, IM);
            }
            new_error /= IMs.size();
            if (new_error <= error) {
                error = new_error;
                if(error < param.converged_error_th){
                    result.converged = true;
                }
                lambda /= param.lambda_factor;
                break;
            }
            else {
                lambda *= param.lambda_factor;
            }

            if (inner_iter > 0) {
                const float dx = std::fabs(prev_delta_inner.x - delta.x);
                const float dy = std::fabs(prev_delta_inner.y - delta.y);
                const float dz = std::fabs(prev_delta_inner.z - delta.z);
                if (std::max(dx, dy) < param.converged_delta_xy_th && dz < param.converged_delta_rot_th) {
                    result.converged = true;
                    break;
                }
            }
            prev_delta_inner = delta;

            if (result.converged) {
                error = new_error;
            }
        }

        if (iter > 0) {
            const float dx = std::fabs(prev_delta_inner.x - delta.x);
            const float dy = std::fabs(prev_delta_inner.y - delta.y);
            const float dz = std::fabs(prev_delta_inner.z - delta.z);
            if (std::max(dx, dy) < param.converged_delta_xy_th && dz < param.converged_delta_rot_th) {
                result.converged = true;
            }
        }

        if (result.converged) {
            result.error = error;
            result.correspondence_num = IMs.size();
            result.H = H_Mat;
            result.b = b_Point;
            break;
        }

        prev_delta = delta;

        if (result.error > error) {
            result.error = error;
            result.correspondence_num = IMs.size();
            result.H = H_Mat;
            result.b = b_Point;
            min_trans_mat = trans_mat;
        }

        if (iter == param.max_iter_num - 1) {
            trans_mat = min_trans_mat;
        }

    }
    if (verbose) {
        if (result.converged) {
            std::cout << "END GICP. ";
        } else {
            std::cout << "END GICP NOT CONVERGED. ";
        }
        std::cout << "ITER: " << iter;
        std::cout << ", ERROR VALUE: " << result.error;
        std::cout << ", CORRESPONDENCE: " << result.correspondence_num << "/" << source_points_size << std::endl;
    }
    result.iter = iter;
    return result;
}

//debug
struct writeSVGSetting {
    float voxel_size = 1.0f;
    int size = 300;
    float scale = 10.0f;
    float ellipse_scale = 3.0f;
    std::string point1_ellipse_color = "pink";
    std::string point2_ellipse_color = "green";
    std::string point1_pt_color = "red";
    std::string point2_pt_color = "darkgreen";
    bool draw_point_covariance = false;
    bool draw_odom_covariance = false;
    bool flip_y = true;
};

inline void writePointsToSVG(const std::vector<ndtcpp::point2>& point_1, const std::vector<ndtcpp::point2>& point_2, const std::string& file_name, writeSVGSetting setting={}) {
    std::ofstream file(file_name);
    if (!file.is_open()) {
        std::cerr << "Cannot open file for writing." << std::endl;
        return;
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
    const float sign = setting.flip_y ? -1.0f: 1.0f;

    file << "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"500\" height=\"500\">\n";
    file << "<rect width=\"100%\" height=\"100%\" fill=\"white\"/>\n";
    for (const auto& point : point_1) {
        file << "<circle cx='" << point.x * scale + offset << "' cy='" << sign * point.y * scale + offset << "' r='1' fill='" << point1_pt_color << "' />\n";
    }

    for (const auto& point : point_2) {
        file << "<circle cx='" << point.x * scale + offset << "' cy='" << sign * point.y * scale + offset << "' r='1' fill='" << point2_pt_color << "' />\n";
    }

    file << "</svg>\n";
    file.close();
}


inline void writePointsToSVG(const std::vector<ndtcpp::point2>& point_1, const std::vector<ndtpoint2>& point_2, const std::string& file_name, writeSVGSetting setting={}) {
    std::ofstream file(file_name);
    if (!file.is_open()) {
        std::cerr << "Cannot open file for writing." << std::endl;
        return;
    }
    const int size = setting.size;
    const float scale = setting.scale;
    const float ellipse_scale = setting.ellipse_scale;
    const float offset = size / 2.0f;
    // const std::string point1_ellipse_color = setting.point1_ellipse_color;
    const std::string point2_ellipse_color = setting.point2_ellipse_color;
    const std::string point1_pt_color = setting.point1_pt_color;
    const std::string point2_pt_color = setting.point2_pt_color;
    const float voxel_size = setting.voxel_size;
    const float sign = setting.flip_y ? -1.0f: 1.0f;

    file << "<svg xmlns='http://www.w3.org/2000/svg' width='" << size << "' height='" << size << "'>\n";
    file << "<g fill='#fff' stroke='#ddd' stroke-width='1'>\n";
    const int voxel_interval = static_cast<int>(std::floor(voxel_size * scale));
    for (size_t i = 0; i < size + voxel_interval; i+=voxel_interval) {
        file << "<path d='M" << i << ",0 L" << i << "," << size << "' />\n";
        file << "<path d='M0," << i << " L" << size << "," << i << "' />\n";
    }
    file << "</g>\n";
    file << "<g fill='#fff' stroke='#000' stroke-width='1'>\n";
    file << "<path d='M0,0 L0," << size << "' />\n";
    file << "<path d='M0,0 L" << size << ",0' />\n";
    file << "<path d='M0," << size << " L" << size << "," << size << "' />\n";
    file << "<path d='M" << size << ",0 L" << size << "," << size << "' />\n";
    file << "</g>\n";

    for (const auto& point : point_1) {
        file << "<circle cx='" << point.x * scale + offset << "' cy='" << sign * point.y * scale + offset << "' r='1' fill='" << point1_pt_color << "' />\n";
    }

    for (const auto& point : point_2) {
        const auto cx = point.mean.x * scale + offset;
        const auto cy = sign * point.mean.y * scale + offset;
        if (setting.draw_point_covariance) {
            const auto& cov = point.cov;
            const float u = 0.5f * ((cov.a + cov.d) + std::sqrt((cov.a - cov.d) * (cov.a - cov.d) + 4.0f * cov.b * cov.b));
            const float v = 0.5f * ((cov.a + cov.d) - std::sqrt((cov.a - cov.d) * (cov.a - cov.d) + 4.0f * cov.b * cov.b));
            const float e1 = (u - cov.a) / cov.b;
            // const float e2 = (v - cov.a) / cov.b;
            // 95%
            const float rx = 2.0f * 2.448f * std::sqrt(u) * ellipse_scale;
            const float ry = 2.0f * 2.448f * std::sqrt(v) * ellipse_scale;
            const auto rot = std::atan(e1) * (180.0f / M_PI);

            file << "<ellipse cx='" << cx << "' cy='" << cy << "' rx='" << rx << "' ry='" << ry << "' fill='" << point2_ellipse_color << "' fill-opacity='0.5' transform='rotate(" << rot << ", " << cx << ", " << cy << ")'/>\n";
        }
        file << "<circle cx='" << cx << "' cy='" << cy << "' r='1' fill='" << point2_pt_color << "' />\n";
    }

    file << "</svg>\n";
    file.close();
}

inline void writePointsToSVG(const std::vector<ndtpoint2>& point_1, const std::vector<ndtpoint2>& point_2, const std::string& file_name, writeSVGSetting setting={}) {
    std::ofstream file(file_name);
    if (!file.is_open()) {
        std::cerr << "Cannot open file for writing." << std::endl;
        return;
    }
    const int size = setting.size;
    const float scale = setting.scale;
    const float ellipse_scale = setting.ellipse_scale;
    const float offset = size / 2.0f;
    const std::string point1_ellipse_color = setting.point1_ellipse_color;
    const std::string point2_ellipse_color = setting.point2_ellipse_color;
    const std::string point1_pt_color = setting.point1_pt_color;
    const std::string point2_pt_color = setting.point2_pt_color;
    const float voxel_size = setting.voxel_size;
    const float sign = setting.flip_y ? -1.0f: 1.0f;

    file << "<svg xmlns='http://www.w3.org/2000/svg' width='" << size << "' height='" << size << "'>\n";
    file << "<g fill='#fff' stroke='#ddd' stroke-width='1'>\n";
    const int voxel_interval = static_cast<int>(std::floor(voxel_size * scale));
    for (size_t i = 0; i < size + voxel_interval; i+=voxel_interval) {
        file << "<path d='M" << i << ",0 L" << i << "," << size << "' />\n";
        file << "<path d='M0," << i << " L" << size << "," << i << "' />\n";
    }
    file << "</g>\n";
    file << "<g fill='#fff' stroke='#000' stroke-width='1'>\n";
    file << "<path d='M0,0 L0," << size << "' />\n";
    file << "<path d='M0,0 L" << size << ",0' />\n";
    file << "<path d='M0," << size << " L" << size << "," << size << "' />\n";
    file << "<path d='M" << size << ",0 L" << size << "," << size << "' />\n";
    file << "</g>\n";

    for (const auto& point : point_1) {
        const auto cx = point.mean.x * scale + offset;
        const auto cy = sign * point.mean.y * scale + offset;

        if (setting.draw_point_covariance) {
            const auto& cov = point.cov;
            const float u = 0.5f * ((cov.a + cov.d) + std::sqrt((cov.a - cov.d) * (cov.a - cov.d) + 4.0f * cov.b * cov.b));
            const float v = 0.5f * ((cov.a + cov.d) - std::sqrt((cov.a - cov.d) * (cov.a - cov.d) + 4.0f * cov.b * cov.b));
            const float e1 = (u - cov.a) / cov.b;
            // const float e2 = (v - cov.a) / cov.b;
            // 95%
            const float rx = 2.0f * 2.448f * std::sqrt(u) * ellipse_scale;
            const float ry = 2.0f * 2.448f * std::sqrt(v) * ellipse_scale;
            const auto rot = std::atan(e1) * (180.0f / M_PI);

            file << "<ellipse cx='" << cx << "' cy='" << cy << "' rx='" << rx << "' ry='" << ry << "' fill='" << point1_ellipse_color << "' fill-opacity='0.5' transform='rotate(" << rot << ", " << cx << ", " << cy << ")'/>\n";
        }
        file << "<circle cx='" << cx << "' cy='" << cy << "' r='1' fill='" << point1_pt_color << "' />\n";
    }

    for (const auto& point : point_2) {
        const auto cx = point.mean.x * scale + offset;
        const auto cy = sign * point.mean.y * scale + offset;

        if (setting.draw_point_covariance) {
            const auto& cov = point.cov;
            const float u = 0.5f * ((cov.a + cov.d) + std::sqrt((cov.a - cov.d) * (cov.a - cov.d) + 4.0f * cov.b * cov.b));
            const float v = 0.5f * ((cov.a + cov.d) - std::sqrt((cov.a - cov.d) * (cov.a - cov.d) + 4.0f * cov.b * cov.b));
            const float e1 = (u - cov.a) / cov.b;
            // const float e2 = (v - cov.a) / cov.b;
            // 95%
            const float rx = 2.0f * 2.448f * std::sqrt(u) * ellipse_scale;
            const float ry = 2.0f * 2.448f * std::sqrt(v) * ellipse_scale;
            const auto rot = std::atan(e1) * (180.0f / M_PI);

            file << "<ellipse cx='" << cx << "' cy='" << cy << "' rx='" << rx << "' ry='" << ry << "' fill='" << point2_ellipse_color << "' fill-opacity='0.5' transform='rotate(" << rot << ", " << cx << ", " << cy << ")'/>\n";
        }
        file << "<circle cx='" << cx << "' cy='" << cy << "' r='1' fill='" << point2_pt_color << "' />\n";
    }

    file << "</svg>\n";
    file.close();
}

inline void writePointsToSVG(const std::vector<ndtpoint2>& point_1, const std::vector<ndtpoint2>& point_2, const ndtcpp::mat3x3& odom, const ndtcpp::mat3x3& H, const std::string& file_name, writeSVGSetting setting={}) {
    std::ofstream file(file_name);
    if (!file.is_open()) {
        std::cerr << "Cannot open file for writing." << std::endl;
        return;
    }
    const int size = setting.size;
    const float scale = setting.scale;
    const float ellipse_scale = setting.ellipse_scale;
    const float offset = size / 2.0f;
    const auto point1_ellipse_color = setting.point1_ellipse_color;
    const auto point2_ellipse_color = setting.point2_ellipse_color;
    const auto point1_pt_color = setting.point1_pt_color;
    const auto point2_pt_color = setting.point2_pt_color;
    const float odom_scale = 5.0f;
    const auto odom_color = "blue";
    const float voxel_size = setting.voxel_size;
    const float sign = setting.flip_y ? -1.0f: 1.0f;

    file << "<svg xmlns='http://www.w3.org/2000/svg' width='" << size << "' height='" << size << "'>\n";
    file << "<g fill='#fff' stroke='#ddd' stroke-width='1'>\n";
    const int voxel_interval = static_cast<int>(std::floor(voxel_size * scale));
    for (size_t i = 0; i < size + voxel_interval; i+=voxel_interval) {
        file << "<path d='M" << i << ",0 L" << i << "," << size << "' />\n";
        file << "<path d='M0," << i << " L" << size << "," << i << "' />\n";
    }
    file << "</g>\n";
    file << "<g fill='#fff' stroke='#000' stroke-width='1'>\n";
    file << "<path d='M0,0 L0," << size << "' />\n";
    file << "<path d='M0,0 L" << size << ",0' />\n";
    file << "<path d='M0," << size << " L" << size << "," << size << "' />\n";
    file << "<path d='M" << size << ",0 L" << size << "," << size << "' />\n";
    file << "</g>\n";

    for (const auto& point : point_1) {
        const auto cx = point.mean.x * scale + offset;
        const auto cy = sign * point.mean.y * scale + offset;
        if (setting.draw_point_covariance) {
            const auto& cov = point.cov;
            const float u = 0.5f * ((cov.a + cov.d) + std::sqrt((cov.a - cov.d) * (cov.a - cov.d) + 4.0f * cov.b * cov.b));
            const float v = 0.5f * ((cov.a + cov.d) - std::sqrt((cov.a - cov.d) * (cov.a - cov.d) + 4.0f * cov.b * cov.b));
            const float e1 = (u - cov.a) / cov.b;
            // const float e2 = (v - cov.a) / cov.b;
            // 95%
            const float rx = 2.0f * 2.448f * std::sqrt(u) * ellipse_scale;
            const float ry = 2.0f * 2.448f * std::sqrt(v) * ellipse_scale;
            const auto rot = std::atan(e1) * (180.0f / M_PI);

            file << "<ellipse cx='" << cx << "' cy='" << cy << "' rx='" << rx << "' ry='" << ry << "' fill='" << point1_ellipse_color << "' fill-opacity='0.5' transform='rotate(" << rot << ", " << cx << ", " << cy << ")'/>\n";
        }
        file << "<circle cx='" << cx << "' cy='" << cy << "' r='1' fill='" << point1_pt_color << "' />\n";
    }

    for (const auto& point : point_2) {
        const auto cx = point.mean.x * scale + offset;
        const auto cy = sign * point.mean.y * scale + offset;
        if (setting.draw_point_covariance) {
            const auto& cov = point.cov;
            const float u = 0.5f * ((cov.a + cov.d) + std::sqrt((cov.a - cov.d) * (cov.a - cov.d) + 4.0f * cov.b * cov.b));
            const float v = 0.5f * ((cov.a + cov.d) - std::sqrt((cov.a - cov.d) * (cov.a - cov.d) + 4.0f * cov.b * cov.b));
            const float e1 = (u - cov.a) / cov.b;
            // const float e2 = (v - cov.a) / cov.b;
            // 95%
            const float rx = 2.0f * 2.448f * std::sqrt(u) * ellipse_scale;
            const float ry = 2.0f * 2.448f * std::sqrt(v) * ellipse_scale;
            const auto rot = std::atan(e1) * (180.0f / M_PI);

            file << "<ellipse cx='" << cx << "' cy='" << cy << "' rx='" << rx << "' ry='" << ry << "' fill='" << point2_ellipse_color << "' fill-opacity='0.5' transform='rotate(" << rot << ", " << cx << ", " << cy << ")'/>\n";
        }
        file << "<circle cx='" << cx << "' cy='" << cy << "' r='1' fill='" << point2_pt_color << "' />\n";
    }
    // odom
    {
        const auto x = odom.c;
        const auto y = odom.f;
        const auto cx = x * scale + offset;
        const auto cy = sign * y * scale + offset;
        if (setting.draw_odom_covariance) {
            const auto cov = H.inv();
            const float u = 0.5f * ((cov.a + cov.d) + std::sqrt((cov.a - cov.d) * (cov.a - cov.d) + 4.0f * cov.b * cov.b));
            const float v = 0.5f * ((cov.a + cov.d) - std::sqrt((cov.a - cov.d) * (cov.a - cov.d) + 4.0f * cov.b * cov.b));
            const float e1 = (u - cov.a) / cov.b;
            const float rx = 2.0f * 2.448f * std::sqrt(u) * ellipse_scale * odom_scale;
            const float ry = 2.0f * 2.448f * std::sqrt(v) * ellipse_scale * odom_scale;
            const auto rot = std::atan(e1) * (180.0f / M_PI);
            file << "<ellipse cx='" << cx << "' cy='" << cy << "' rx='" << rx << "' ry='" << ry << "' fill='" << odom_color << "' fill-opacity='0.5' transform='rotate(" << rot << ", " << cx << ", " << cy << ")'/>\n";
        }
        file << "<circle cx='" << cx << "' cy='" << cy << "' r='1' fill='" << odom_color << "' />\n";

    }

    file << "</svg>\n";
    file.close();
}

} // namespace ndt_cpp

#endif // NDTCPP_NDT_CPU_SINGLE_HPP_
