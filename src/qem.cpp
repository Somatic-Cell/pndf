#include "pndf/qem.hpp"
#include <algorithm>
#include <cmath>

namespace pndf {

static bool affine_coefficients(Vec2 p0, Vec2 p1, Vec2 p2, double f0, double f1, double f2,
                                double& a, double& b, double& c) {
    const double det = p0.x*(p1.y-p2.y) + p1.x*(p2.y-p0.y) + p2.x*(p0.y-p1.y);
    if (std::abs(det) < 1e-20) return false;
    a = (f0*(p1.y-p2.y) + f1*(p2.y-p0.y) + f2*(p0.y-p1.y)) / det;
    b = (p0.x*(f1-f2) + p1.x*(f2-f0) + p2.x*(f0-f1)) / det;
    c = (p0.x*(p1.y*f2-p2.y*f1) + p1.x*(p2.y*f0-p0.y*f2) + p2.x*(p0.y*f1-p1.y*f0)) / det;
    return true;
}

Mat5 triangle_quadric_normal_only(const TorusMesh& mesh, const Face& face) {
    Mat5 q;
    const auto uv = unwrapped_face_uvs(mesh, face);
    const Vec2 n0 = mesh.vertices[face.v[0]].nxy;
    const Vec2 n1 = mesh.vertices[face.v[1]].nxy;
    const Vec2 n2 = mesh.vertices[face.v[2]].nxy;
    double ax, bx, cx, ay, by, cy;
    if (!affine_coefficients(uv[0], uv[1], uv[2], n0.x, n1.x, n2.x, ax, bx, cx)) return q;
    if (!affine_coefficients(uv[0], uv[1], uv[2], n0.y, n1.y, n2.y, ay, by, cy)) return q;
    q.add_outer({ax, bx, -1.0, 0.0, cx});
    q.add_outer({ay, by, 0.0, -1.0, cy});
    return q;
}

Mat5 triangle_quadric_4d(const TorusMesh& mesh, const Face& face) {
    Mat5 q;
    const auto uv = unwrapped_face_uvs(mesh, face);
    const Vec2 n0 = mesh.vertices[face.v[0]].nxy;
    const Vec2 n1 = mesh.vertices[face.v[1]].nxy;
    const Vec2 n2 = mesh.vertices[face.v[2]].nxy;
    const Vec4 p0{uv[0].x, uv[0].y, n0.x, n0.y};
    const Vec4 p1{uv[1].x, uv[1].y, n1.x, n1.y};
    const Vec4 p2{uv[2].x, uv[2].y, n2.x, n2.y};
    const Vec4 e1 = p1 - p0;
    const Vec4 e2 = p2 - p0;
    const double g00 = dot4(e1,e1);
    const double g01 = dot4(e1,e2);
    const double g11 = dot4(e2,e2);
    const double det = g00*g11 - g01*g01;
    if (std::abs(det) < 1e-20) return q;
    const double inv00 = g11 / det;
    const double inv01 = -g01 / det;
    const double inv11 = g00 / det;
    // P = I - E (E^T E)^-1 E^T, projection onto the 2D normal complement in R4.
    double P[4][4]{};
    for (int i=0;i<4;++i) P[i][i] = 1.0;
    for (int i=0;i<4;++i) {
        for (int j=0;j<4;++j) {
            const double term = e1[i]*(inv00*e1[j] + inv01*e2[j]) + e2[i]*(inv01*e1[j] + inv11*e2[j]);
            P[i][j] -= term;
        }
    }
    // Homogeneous Q for (p - p0)^T P (p - p0), variables [u,v,nx,ny,1].
    double Pp0[4]{};
    for (int i=0;i<4;++i) for (int j=0;j<4;++j) Pp0[i] += P[i][j] * p0[j];
    for (int i=0;i<4;++i) for (int j=0;j<4;++j) q(i,j) += P[i][j];
    for (int i=0;i<4;++i) {
        q(i,4) += -Pp0[i];
        q(4,i) += -Pp0[i];
    }
    double c = 0.0;
    for (int i=0;i<4;++i) c += p0[i] * Pp0[i];
    q(4,4) += c;
    return q;
}

Mat5 point_normal_stabilizer(Vec2 nxy, double weight) {
    Mat5 q;
    q.add_outer({0.0, 0.0, 1.0, 0.0, -nxy.x}, weight);
    q.add_outer({0.0, 0.0, 0.0, 1.0, -nxy.y}, weight);
    return q;
}

static Vec2 solve_normal_for_fixed_uv(const Mat5& q, Vec2 uv, Vec2 fallback) {
    const double A00 = q(2,2);
    const double A01 = q(2,3);
    const double A11 = q(3,3);
    const double b0 = q(2,0)*uv.x + q(2,1)*uv.y + q(2,4);
    const double b1 = q(3,0)*uv.x + q(3,1)*uv.y + q(3,4);
    double nx, ny;
    if (!solve2x2(A00, A01, A11, -b0, -b1, nx, ny)) return fallback;
    return clamp_projected_normal({nx, ny});
}

Placement best_placement_unwrapped(const Mat5& q, Vec2 u0, Vec2 n0, Vec2 u1_unwrapped, Vec2 n1, QemMode mode) {
    // u0 and u1_unwrapped must be in the same local torus chart.  The
    // returned uv is deliberately not wrapped back to [0,1): doing so would
    // make seam-crossing collapses attract to the arbitrary cut line.
    const Vec2 umid = (u0 + u1_unwrapped) * 0.5;
    const Vec2 nmid = clamp_projected_normal((n0+n1)*0.5);
    std::array<Vec2,3> uv_candidates{umid, u0, u1_unwrapped};
    std::array<Vec2,3> n_candidates{nmid, n0, n1};
    Placement best;
    best.cost = std::numeric_limits<double>::infinity();
    for (int i=0;i<3;++i) {
        Vec2 uv = uv_candidates[i];
        Vec2 n = (mode == QemMode::NormalOnly) ? solve_normal_for_fixed_uv(q, uv, n_candidates[i]) : n_candidates[i];
        const double cost = eval_quadric(q, {uv.x, uv.y, n.x, n.y, 1.0});
        if (cost < best.cost) best = {uv, n, cost};
    }
    return best;
}

Placement best_placement(const Mat5& q, Vec2 u0, Vec2 n0, Vec2 u1, Vec2 n1, QemMode mode) {
    return best_placement_unwrapped(q, u0, n0, unwrap_near(u1, u0), n1, mode);
}

} // namespace pndf
