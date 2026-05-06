#include "pndf/eval.hpp"
#include <algorithm>
#include <fstream>
#include <numeric>
#include <stdexcept>

namespace pndf {

NormalErrorStats evaluate_vertex_sample_error(const NormalMap& gt, const TorusMesh& mesh) {
    std::vector<double> errors;
    errors.reserve(mesh.vertices.size());
    for (const Vertex& v : mesh.vertices) {
        if (!v.alive) continue;
        Vec2 ref = gt.sample_periodic(v.uv);
        errors.push_back(angular_error_deg(ref, v.nxy));
    }
    NormalErrorStats s;
    if (errors.empty()) return s;
    std::sort(errors.begin(), errors.end());
    const double sum = std::accumulate(errors.begin(), errors.end(), 0.0);
    s.mean_deg = sum / errors.size();
    auto pct = [&](double p) {
        size_t idx = static_cast<size_t>(p * (errors.size()-1));
        return errors[idx];
    };
    s.p95_deg = pct(0.95);
    s.p99_deg = pct(0.99);
    s.max_deg = errors.back();
    return s;
}

void write_stats_csv(const NormalErrorStats& s, const std::string& path) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot write stats: " + path);
    f << "mean_deg,p95_deg,p99_deg,max_deg\n";
    f << s.mean_deg << "," << s.p95_deg << "," << s.p99_deg << "," << s.max_deg << "\n";
}

} // namespace pndf
