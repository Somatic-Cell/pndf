#pragma once
#include <string>
#include <vector>
#include "pndf/math.hpp"

namespace pndf {

struct NormalMap {
    int width = 0;
    int height = 0;
    std::vector<Vec2> nxy;

    Vec2 at_wrapped(int x, int y) const;
    Vec2 sample_periodic(Vec2 uv) const;
};

NormalMap read_nxy_binary(const std::string& path);
void write_nxy_binary(const NormalMap& map, const std::string& path);

} // namespace pndf
