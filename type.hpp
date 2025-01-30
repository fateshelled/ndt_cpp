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

#ifndef NDTCPP_TYPE_H_
#define NDTCPP_TYPE_H_

namespace ndtcpp
{
    struct point2;
    struct point2_T;
    struct point3;
    struct point3_T;

    struct point2{
        float x, y;
        static point2 zeros() {
            return {0.f, 0.f};
        }
        point2_T transpose() const;

        float norm() const {
            return std::sqrt(x * x + y * y);
        }
    };

    struct point2_T {
        float x, y;
        static point2_T zeros() {
            return point2_T {0.f, 0.f};
        }
        point2 transpose() const;
    };

    point2_T point2::transpose() const {
        return point2_T {x, y};
    }
    point2 point2_T::transpose() const {
        return point2 {x, y};
    }

    struct point3{
        float x, y, z;
        static point3 zeros() {
            return {0.f, 0.f, 0.f};
        }
        point3_T transpose() const;
    };

    struct point3_T{
        float x, y, z;
        static point3 zeros() {
            return {0.f, 0.f, 0.f};
        }
        point3 transpose() const;
    };

    point3 point3_T::transpose() const {
        return {x, y, z};
    }
    point3_T point3::transpose() const {
        return {x, y, z};
    }

    struct mat2x2{
        float a, b;
        float c, d;
        static mat2x2 eye() {
            return {1.f, 0.f, 0.f, 1.f};
        }
        static mat2x2 zeros() {
            return {0.f, 0.f, 0.f, 0.f};
        }
        static mat2x2 diagonal(float a, float d) {
            return {a, 0.f, 0.f, d};
        }
        mat2x2 transpose() const {
            return {a, c, b, d};
        }
        mat2x2 diagonal() {
            return {a, 0.f, 0.f, d};
        }
        point2 diagonal_vector() const {
            return {a, d};
        }
        mat2x2 inv() const {
            const auto val = 1.0f / (
                a * d - b * c
            );
            return {
                val * d , -val * b,
                -val * c, val * a
            };
        }
        float det() const {
            return a * d - b * c;
        }
    };

    struct mat3x3{
        float a, b, c;
        float d, e, f;
        float g, h, i;
        static mat3x3 eye() {
            return {1.f, 0.f, 0.f, 0.f, 1.f, 0.f, 0.f, 0.f, 1.f};
        }
        static mat3x3 zeros() {
            return {0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f};
        }
        mat3x3 transpose() const {
            return {
                a, d, g,
                b, e, h,
                c, f, i,
            };
        }
        mat3x3 diagonal() const {
            return {
                a,    0.0f, 0.0f,
                0.0f,    e, 0.0f,
                0.0f, 0.0f, i,
            };
        }
        mat3x3 inv() const {

            const auto val = 1.0f / (
                a * e * i +
                b * f * g +
                c * d * h -
                c * e * g -
                b * d * i -
                a * f * h
            );

            ndtcpp::mat3x3 inv_mat;
            inv_mat.a = e * i - f * h;
            inv_mat.b = b * i - c * h;
            inv_mat.c = b * f - c * e;

            inv_mat.d = d * i - f * g;
            inv_mat.e = a * i - c * g;
            inv_mat.f = a * f - c * d;

            inv_mat.g = d * h - e * g;
            inv_mat.h = a * h - b * g;
            inv_mat.i = a * e - b * d;

            inv_mat.a = inv_mat.a * val;
            inv_mat.b = inv_mat.b * val * -1.0f;
            inv_mat.c = inv_mat.c * val;

            inv_mat.d = inv_mat.d * val * -1.0f;
            inv_mat.e = inv_mat.e * val;
            inv_mat.f = inv_mat.f * val * -1.0f;

            inv_mat.g = inv_mat.g * val;
            inv_mat.h = inv_mat.h * val * -1.0f;
            inv_mat.i = inv_mat.i * val;

            return inv_mat;
        }
        float det() const {
            return a * e * i + b * f * g + c * d * h - c * e * g - b * d * i - a * f * h;
        }

    };
} // namespace ndtcpp

#endif // NDTCPP_TYPE_H_
