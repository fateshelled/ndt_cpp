#pragma once

#include <cstddef>
#include <tuple>
#include <functional>

namespace ndtcpp {
struct tuple_int_hash {
    size_t operator()(const std::tuple<int, int>& v) const {
        const auto hash0 = std::hash<int>{}(std::get<0>(v));
        const auto hash1 = std::hash<int>{}(std::get<1>(v));
        size_t seed = 0;
        seed ^= hash0 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        seed ^= hash1 + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        return seed;
    }
};
} // namespace ndtcpp
