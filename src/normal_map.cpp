#include "pndf/normal_map.hpp"
#include <fstream>
#include <stdexcept>
#include <cstdint>

namespace pndf {

Vec2 NormalMap::at_wrapped(int x, int y) const {
    if (width <= 0 || height <= 0) return {};
    x %= width; if (x < 0) x += width;
    y %= height; if (y < 0) y += height;
    return nxy[static_cast<size_t>(y) * width + x];
}

Vec2 NormalMap::sample_periodic(Vec2 uv) const {
    const double fx = wrap01(uv.x) * width - 0.5;
    const double fy = wrap01(uv.y) * height - 0.5;
    const int x0 = static_cast<int>(std::floor(fx));
    const int y0 = static_cast<int>(std::floor(fy));
    const double tx = fx - x0;
    const double ty = fy - y0;
    const Vec2 v00 = at_wrapped(x0, y0);
    const Vec2 v10 = at_wrapped(x0+1, y0);
    const Vec2 v01 = at_wrapped(x0, y0+1);
    const Vec2 v11 = at_wrapped(x0+1, y0+1);
    Vec2 out = v00 * ((1-tx)*(1-ty)) + v10 * (tx*(1-ty)) + v01 * ((1-tx)*ty) + v11 * (tx*ty);
    return clamp_projected_normal(out);
}

NormalMap read_nxy_binary(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("Cannot open nxy binary: " + path);
    int32_t n = 0;
    f.read(reinterpret_cast<char*>(&n), sizeof(int32_t));
    if (n <= 0) throw std::runtime_error("Invalid nxy binary size");
    NormalMap map;
    map.width = n;
    map.height = n;
    const size_t total = static_cast<size_t>(n) * n;
    std::vector<float> nx(total), ny(total);
    f.read(reinterpret_cast<char*>(nx.data()), static_cast<std::streamsize>(total * sizeof(float)));
    f.read(reinterpret_cast<char*>(ny.data()), static_cast<std::streamsize>(total * sizeof(float)));
    if (!f) throw std::runtime_error("Unexpected EOF while reading nxy binary");
    map.nxy.resize(total);
    for (size_t i=0;i<total;++i) map.nxy[i] = clamp_projected_normal({nx[i], ny[i]});
    return map;
}

void write_nxy_binary(const NormalMap& map, const std::string& path) {
    if (map.width != map.height) throw std::runtime_error("Only square maps are supported by this binary format");
    std::ofstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("Cannot write nxy binary: " + path);
    int32_t n = map.width;
    f.write(reinterpret_cast<const char*>(&n), sizeof(int32_t));
    const size_t total = map.nxy.size();
    std::vector<float> nx(total), ny(total);
    for (size_t i=0;i<total;++i) { nx[i] = static_cast<float>(map.nxy[i].x); ny[i] = static_cast<float>(map.nxy[i].y); }
    f.write(reinterpret_cast<const char*>(nx.data()), static_cast<std::streamsize>(total * sizeof(float)));
    f.write(reinterpret_cast<const char*>(ny.data()), static_cast<std::streamsize>(total * sizeof(float)));
}

} // namespace pndf
