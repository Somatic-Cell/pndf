#pragma once
#include <array>
#include "pndf/math.hpp"
#include "pndf/torus_mesh.hpp"

namespace pndf {

enum class QemMode {
    NormalOnly,
    Qem4D
};

struct Placement {
    Vec2 uv;
    Vec2 nxy;
    double cost = 0.0;
};

Mat5 triangle_quadric_normal_only(const TorusMesh& mesh, const Face& face);
Mat5 triangle_quadric_4d(const TorusMesh& mesh, const Face& face);
Mat5 point_normal_stabilizer(Vec2 nxy, double weight = 1.0);
Placement best_placement(const Mat5& q, Vec2 u0, Vec2 n0, Vec2 u1, Vec2 n1, QemMode mode);

} // namespace pndf
