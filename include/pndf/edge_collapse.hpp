#pragma once
#include <queue>
#include <string>
#include "pndf/torus_mesh.hpp"
#include "pndf/qem.hpp"

namespace pndf {

struct CollapseOptions {
    int target_vertices = 0;
    QemMode mode = QemMode::NormalOnly;
    bool lock_boundary = false;       // false = true torus; true = conservative rectangular debug mode
    int rebuild_interval = 50000;
    bool verbose = true;
};

struct CollapseStats {
    int accepted = 0;
    int rejected = 0;
    int final_vertices = 0;
    int final_faces = 0;
};

CollapseStats simplify_qem(TorusMesh& mesh, const CollapseOptions& opt);

} // namespace pndf
