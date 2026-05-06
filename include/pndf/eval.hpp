#pragma once
#include <string>
#include "pndf/normal_map.hpp"
#include "pndf/torus_mesh.hpp"

namespace pndf {

struct NormalErrorStats {
    double mean_deg = 0.0;
    double p95_deg = 0.0;
    double p99_deg = 0.0;
    double max_deg = 0.0;
};

NormalErrorStats evaluate_vertex_sample_error(const NormalMap& gt, const TorusMesh& mesh);
void write_stats_csv(const NormalErrorStats& s, const std::string& path);

} // namespace pndf
