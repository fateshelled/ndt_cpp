#pragma once

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
#include <unordered_map>

#include "type.hpp"
#include "ndt-cpu-single.hpp"


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
        const ndtcpp::point2 pos = {odom.c, odom.f};
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
