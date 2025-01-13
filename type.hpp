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
    struct point2{
        float x, y;
        static point2 zeros() {
            return {0.f, 0.f};
        }
    };
    struct point3{
        float x, y, z;
        static point3 zeros() {
            return {0.f, 0.f, 0.f};
        }
    };

    struct mat2x2{
        float a, b;
        float c, d;
        static mat2x2 eye() {
            return {1.f, 0.f, 0.f, 1.f};
        }
        static mat2x2 zeros() {
            return {0.f, 0.f, 0.f, 0.f};
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
    };
} // namespace ndtcpp

#endif // NDTCPP_TYPE_H_
