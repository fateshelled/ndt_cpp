#pragma once

#include "ndt-cpu-single.hpp"
#include "gicp.hpp"
#include <numeric>


namespace ndtcpp {

void extract_feature_points(const std::vector<ndtcpp::ndtpoint2>& points, const float corner_threshold, std::vector<ndtcpp::ndtpoint2>& lines, std::vector<ndtcpp::ndtpoint2>& corners) {
    lines.clear();
    corners.clear();
    lines.reserve(points.size());
    corners.reserve(points.size());
    for (const auto& pt: points) {
        const auto eigen = ndtcpp::compute_eigen(pt.cov);
        const auto& [eigen_val0, eigen_vec0] = eigen[0];
        const auto& [eigen_val1, eigen_vec1] = eigen[1];
        const float ratio = eigen_val0 / eigen_val1;
        ndtcpp::mat2x2 mat = {eigen_vec0.x, eigen_vec1.x, eigen_vec0.y, eigen_vec1.y};
        if (ratio < corner_threshold) {
            mat = mat * ndtcpp::mat2x2::diagonal(5.0f, 5.0f) * mat.transpose();
            corners.emplace_back(pt.mean, mat);
        } else {
            mat = mat * ndtcpp::mat2x2::diagonal(1.0f, .1f) * mat.transpose();
            lines.emplace_back(pt.mean, mat);
        }
    }
}

inline ndtcpp::scan_matching_result loam_scan_matching(
    ndtcpp::mat3x3& trans_mat,
    const std::vector<ndtcpp::ndtpoint2>& source_lines,
    const std::vector<ndtcpp::ndtpoint2>& source_corners,
    std::vector<ndtcpp::ndtpoint2>& target_lines,
    std::vector<ndtcpp::ndtpoint2>& target_corners,
    bool verbose = false,
    const ndtcpp::GICP_PARAMS& param = ndtcpp::GICP_PARAMS()
) {
    auto max_distance = param.max_correspondence_distance;

    float lambda = param.init_lambda;

    ndtcpp::point3 prev_delta;
    ndtcpp::scan_matching_result result = {};
    ndtcpp::mat3x3 min_trans_mat = trans_mat;

    const size_t source_lines_size = source_lines.size();
    const size_t target_lines_size = target_lines.size();

    const size_t source_corners_size = source_corners.size();
    const size_t target_corners_size = target_corners.size();

    kdtree::construct(target_lines.begin(), target_lines.end());
    kdtree::construct(target_corners.begin(), target_corners.end());

    size_t iter = 0;
    for(iter = 0; iter < param.max_iter_num; ++iter){
        const float max_distance2 = max_distance * max_distance;
        auto H_Mat = ndtcpp::mat3x3::zeros();
        auto b_Point = ndtcpp::point3::zeros();
        float error = 0.0f;

        const ndtcpp::mat2x2 trans_mat2x2 = {trans_mat.a, trans_mat.b, trans_mat.d, trans_mat.e};
        const auto trans_mat2x2_T = trans_mat2x2.transpose();

        std::vector<std::tuple<ndtcpp::mat3x3, ndtcpp::point2, int>> IMs_line;
        std::vector<std::tuple<ndtcpp::mat3x3, ndtcpp::point2, int>> IMs_corner;

        for(size_t point_iter = 0; point_iter < source_lines_size; point_iter += param.point_step){
            ndtcpp::ndtpoint2 query_point = {
                transformPointCopy(trans_mat, source_lines[point_iter].mean),
                trans_mat2x2 * source_lines[point_iter].cov * trans_mat2x2_T
            };
            ndtcpp::ndtpoint2 target_point;
            float target_sq_distance;
            kdtree::search_knn(target_lines.begin(), target_lines.end(), &target_point, &target_sq_distance, 1, query_point);

            if(target_sq_distance > max_distance2){ continue; }

            const auto identity_plus_target_cov = ndtcpp::mat3x3{
                target_point.cov.a + 1.0f, target_point.cov.b       , 0.0f,
                target_point.cov.c       , target_point.cov.d + 1.0f, 0.0f,
                0.0f, 0.0f, 1.0f
            };

            query_point.cov = trans_mat2x2 * source_lines[point_iter].cov * trans_mat2x2_T;
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

            const ndtcpp::point2 v_point = transformPointCopy(trans_mat, skewd(source_lines[point_iter].mean));

            const auto mat_J = ndtcpp::mat3x3{
                trans_mat.a * -1.0f, trans_mat.b * -1.0f, v_point.x,
                trans_mat.d * -1.0f, trans_mat.e * -1.0f, v_point.y,
                trans_mat.g * -1.0f, trans_mat.h * -1.0f, trans_mat.i * -1.0f
            };

            const ndtcpp::mat3x3 mat_J_TxIM = mat_J.transpose() * IM;

            H_Mat += (mat_J_TxIM * mat_J);      // J.T * IM * J
            b_Point += (mat_J_TxIM * residual); // J.T * IM * residual

            error += ndtcpp::calc_gicp_error(query_point.mean, target_point.mean, IM);
            IMs_line.emplace_back(IM, target_point.mean, point_iter);
        }

        for(size_t point_iter = 0; point_iter < source_corners_size; point_iter += param.point_step){
            ndtcpp::ndtpoint2 query_point = {
                transformPointCopy(trans_mat, source_corners[point_iter].mean),
                trans_mat2x2 * source_corners[point_iter].cov * trans_mat2x2_T
            };
            ndtcpp::ndtpoint2 target_point;
            float target_sq_distance;
            kdtree::search_knn(target_corners.begin(), target_corners.end(), &target_point, &target_sq_distance, 1, query_point);

            if(target_sq_distance > max_distance2){ continue; }

            const auto identity_plus_target_cov = ndtcpp::mat3x3{
                target_point.cov.a + 1.0f, target_point.cov.b       , 0.0f,
                target_point.cov.c       , target_point.cov.d + 1.0f, 0.0f,
                0.0f, 0.0f, 1.0f
            };

            query_point.cov = trans_mat2x2 * source_corners[point_iter].cov * trans_mat2x2_T;
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

            const ndtcpp::point2 v_point = transformPointCopy(trans_mat, skewd(source_corners[point_iter].mean));

            const auto mat_J = ndtcpp::mat3x3{
                trans_mat.a * -1.0f, trans_mat.b * -1.0f, v_point.x,
                trans_mat.d * -1.0f, trans_mat.e * -1.0f, v_point.y,
                trans_mat.g * -1.0f, trans_mat.h * -1.0f, trans_mat.i * -1.0f
            };

            const ndtcpp::mat3x3 mat_J_TxIM = mat_J.transpose() * IM;

            H_Mat += (mat_J_TxIM * mat_J);      // J.T * IM * J
            b_Point += (mat_J_TxIM * residual); // J.T * IM * residual

            error += ndtcpp::calc_gicp_error(query_point.mean, target_point.mean, IM);
            IMs_corner.emplace_back(IM, target_point.mean, point_iter);
        }
        b_Point *= -1.0f;
        const size_t total_IMs = IMs_line.size() + IMs_corner.size();
        if (total_IMs > 0) {
            error /= IMs_line.size() + IMs_corner.size();
        }

        if (total_IMs < param.min_correspondence) {
            result.error = error;
            result.correspondence_num = total_IMs;
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
            for (const auto& [IM, target, point_iter]: IMs_line) {
                const auto trans_source = transformPointCopy(trans_mat, source_lines[point_iter].mean);
                new_error += ndtcpp::calc_gicp_error(trans_source, target, IM);
            }
            for (const auto& [IM, target, point_iter]: IMs_corner) {
                const auto trans_source = transformPointCopy(trans_mat, source_corners[point_iter].mean);
                new_error += ndtcpp::calc_gicp_error(trans_source, target, IM);
            }
            new_error /= total_IMs;
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
            result.correspondence_num = total_IMs;
            result.H = H_Mat;
            result.b = b_Point;
            break;
        }

        prev_delta = delta;

        if (result.error > error) {
            result.error = error;
            result.correspondence_num = total_IMs;
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
            std::cout << "END LOAM. ";
        } else {
            std::cout << "END LOAM NOT CONVERGED. ";
        }
        std::cout << "ITER: " << iter;
        std::cout << ", ERROR VALUE: " << result.error;
        std::cout << ", CORRESPONDENCE: " << result.correspondence_num << "/" << source_lines_size + source_corners_size << std::endl;
    }
    result.iter = iter;
    return result;
}

} // namespace ndtcpp
