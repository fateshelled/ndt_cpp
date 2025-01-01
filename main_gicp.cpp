#include <numeric>
#include <chrono>
#include "ndt-cpu-single.hpp"


int main(void){

    auto scan_points1 = ndtcpp::read_scan_points("./data/scan_1.txt");
    auto target_points = ndtcpp::read_scan_points("./data/scan_2.txt");

    std::vector<double> durations;
    const size_t N = 10;

    for (size_t i = 0; i < N; ++i) {

        auto source = scan_points1;
        auto target = target_points;
        auto trans_mat1 = ndtcpp::makeTransformationMatrix(1.0f, 0.0f, 0.5f);
        transformPointsZeroCopy(trans_mat1, source);

        auto source_ndt = std::vector<ndtcpp::ndtpoint2>();
        auto target_ndt = std::vector<ndtcpp::ndtpoint2>();
        auto start_time = std::chrono::high_resolution_clock::now();

        const bool verbose = true;
        ndtcpp::compute_ndt_points_downsampling(source, source_ndt);
        ndtcpp::compute_ndt_points_downsampling(target, target_ndt);
        ndtcpp::gicp_scan_matching(trans_mat1, source_ndt, target_ndt, verbose);

        auto end_time = std::chrono::high_resolution_clock::now();

        ndtcpp::transformPointsZeroCopy(trans_mat1, source);

        const ndtcpp::mat2x2 trans2x2 = {trans_mat1.a, trans_mat1.b, trans_mat1.d, trans_mat1.e};
        const ndtcpp::mat2x2 trans2x2_T = {trans2x2.a, trans2x2.c, trans2x2.b, trans2x2.d};
        for (auto& pt: source_ndt) {
            pt.mean = ndtcpp::transformPointCopy(trans_mat1, pt.mean);
            pt.cov = trans2x2 * pt.cov * trans2x2_T;
        }

        //debug
        auto microsec = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count() / 1e6;
        durations.push_back(microsec);
        if (i == N - 1) {
            ndtcpp::writePointsToSVG(source, target, "scan_points_gicp.svg");
            ndtcpp::writePointsToSVG(source_ndt, target_ndt, "scan_points_gicp_cov.svg");
        }
    }
    const double mean = std::accumulate(durations.begin(), durations.end(), 0.0) / durations.size();
    std::cout << "MEAN: " << mean << " mill sec" << std::endl;

}
