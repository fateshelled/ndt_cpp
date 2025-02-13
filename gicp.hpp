#pragma once

#include "ndt-cpu-single.hpp"
#include <numeric>


namespace ndtcpp {

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
            const float dx = std::fabs(prev_delta.x - delta.x);
            const float dy = std::fabs(prev_delta.y - delta.y);
            const float dz = std::fabs(prev_delta.z - delta.z);
            result.converged = std::max(dx, dy) < param.converged_delta_xy_th && dz < param.converged_delta_rot_th;
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

} // namespace ndtcpp
