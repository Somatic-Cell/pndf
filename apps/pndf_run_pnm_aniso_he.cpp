#include "pndf/edge_collapse.hpp"
#include "pndf/mesh_io.hpp"
#include "pndf/normal_map.hpp"
#include "pndf/torus_mesh.hpp"

#include <pmp/surface_mesh.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <limits>
#include <queue>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <vector>

using namespace pndf;

namespace {

struct Options {
    std::string input;
    std::string output;
    int target = 0;
    QemMode mode = QemMode::NormalOnly;
    int iterations = 5;
    int flip_passes = 1;
    int relocate_passes = 1;
    double relocate_step_pixels = 1.0;
    double min_area_ratio = 0.05;
    double flip_improvement_eps = 1e-9;
    int flip_max_accepted_per_pass = 0; // 0 = unlimited
    int rebuild_interval = 50000;
    int progress_interval = 50000;
    bool quiet = false;

    // Add:
    std::string checkpoint_dir;
    std::string checkpoint_prefix = "pnm_he";
};

void usage() {
    std::cerr
        << "Usage:\n"
        << "  pndf_run_pnm_aniso_he --input nxy.bin --output meshbin --target N [options]\n\n"
        << "Options:\n"
        << "  --mode normal|4d                 QEM mode used for staged collapse. Default: normal.\n"
        << "  --iterations N                   Staged collapse/alignment iterations. Default: 5.\n"
        << "  --flip-passes N                  Half-edge dirty-queue flip passes after each stage. Default: 1.\n"
        << "  --relocate-passes N              Vertex relocation passes after each stage. Default: 1.\n"
        << "  --relocate-step-pixels X         Candidate relocation step in texels. Default: 1.0.\n"
        << "  --min-area-ratio X               Reject/limit moves that shrink incident area below X. Default: 0.05.\n"
        << "  --flip-improvement-eps X         Required local energy decrease for flips. Default: 1e-9.\n"
        << "  --flip-max-accepted-per-pass N   Debug cap. 0 means unlimited.\n"
        << "  --rebuild-interval N             QEM heap rebuild interval. Default: 50000.\n"
        << "  --progress-interval N            QEM progress interval. Default: 50000.\n"
        << "  --quiet                          Reduce progress output.\n\n"
        << "This executable keeps the proven QEM collapse engine, but runs the\n"
        << "anisotropic alignment stage on a pmp::SurfaceMesh half-edge structure.\n"
        << "The flip stage follows a local dirty-queue update pattern, and the\n"
        << "relocation stage uses orientation-root step limiting.\n"
        << "  --checkpoint-dir DIR            Write stage checkpoint meshbins into DIR.\n"
        << "  --checkpoint-prefix NAME        Prefix for checkpoint files. Default: pnm_he.\n";
}

struct HeMesh {
    pmp::SurfaceMesh mesh;
    pmp::VertexProperty<Vec2> uv;
    pmp::VertexProperty<Vec2> nxy;
    int width = 0;
    int height = 0;
};

pmp::Point point_from_uv(Vec2 uv) {
    return pmp::Point(static_cast<pmp::Scalar>(uv.x),
                      static_cast<pmp::Scalar>(uv.y),
                      static_cast<pmp::Scalar>(0));
}

bool valid_vertex(const HeMesh& he, pmp::Vertex v) {
    return he.mesh.is_valid(v) && !he.mesh.is_deleted(v);
}

bool valid_edge(const HeMesh& he, pmp::Edge e) {
    return he.mesh.is_valid(e) && !he.mesh.is_deleted(e);
}

bool valid_face(const HeMesh& he, pmp::Face f) {
    return he.mesh.is_valid(f) && !he.mesh.is_deleted(f);
}

HeMesh to_halfedge_mesh(const TorusMesh& src) {
    HeMesh he;
    he.width = src.width;
    he.height = src.height;
    he.uv = he.mesh.add_vertex_property<Vec2>("v:uv");
    he.nxy = he.mesh.add_vertex_property<Vec2>("v:nxy");

    std::vector<pmp::Vertex> remap(src.vertices.size());
    for (std::size_t i = 0; i < src.vertices.size(); ++i) {
        if (!src.vertices[i].alive) continue;
        pmp::Vertex v = he.mesh.add_vertex(point_from_uv(src.vertices[i].uv));
        he.uv[v] = src.vertices[i].uv;
        he.nxy[v] = clamp_projected_normal(src.vertices[i].nxy);
        remap[i] = v;
    }

    for (const Face& f : src.faces) {
        if (!f.alive) continue;
        if (f.v[0] < 0 || f.v[1] < 0 || f.v[2] < 0) continue;
        if (!src.vertices[f.v[0]].alive || !src.vertices[f.v[1]].alive || !src.vertices[f.v[2]].alive) continue;
        if (f.v[0] == f.v[1] || f.v[1] == f.v[2] || f.v[2] == f.v[0]) continue;
        he.mesh.add_triangle(remap[f.v[0]], remap[f.v[1]], remap[f.v[2]]);
    }

    return he;
}

TorusMesh from_halfedge_mesh(HeMesh& he) {
    he.mesh.garbage_collection();

    TorusMesh out;
    out.width = he.width;
    out.height = he.height;

    std::vector<int> remap(he.mesh.vertices_size(), -1);
    for (pmp::Vertex v : he.mesh.vertices()) {
        if (he.mesh.is_deleted(v)) continue;
        Vertex pv;
        pv.uv = fract(he.uv[v]);
        pv.q_origin = pv.uv;
        pv.nxy = clamp_projected_normal(he.nxy[v]);
        pv.alive = true;
        pv.version = 0;
        remap[v.idx()] = static_cast<int>(out.vertices.size());
        out.vertices.push_back(pv);
    }

    for (pmp::Face f : he.mesh.faces()) {
        if (he.mesh.is_deleted(f)) continue;
        std::array<int, 3> tri{};
        int k = 0;
        for (pmp::Vertex v : he.mesh.vertices(f)) {
            if (k < 3) tri[k++] = remap[v.idx()];
        }
        if (k == 3 && tri[0] >= 0 && tri[1] >= 0 && tri[2] >= 0 &&
            tri[0] != tri[1] && tri[1] != tri[2] && tri[2] != tri[0]) {
            Face pf;
            pf.v = tri;
            pf.alive = true;
            out.faces.push_back(pf);
        }
    }

    out.rebuild_vertex_face_adjacency();
    return out;
}

TorusMesh snapshot_halfedge_mesh(const HeMesh& he) {
    TorusMesh out;
    out.width = he.width;
    out.height = he.height;

    std::vector<int> remap(he.mesh.vertices_size(), -1);

    for (pmp::Vertex v : he.mesh.vertices()) {
        if (!valid_vertex(he, v)) continue;

        Vertex pv;
        pv.uv = fract(he.uv[v]);
        pv.q_origin = pv.uv;
        pv.nxy = clamp_projected_normal(he.nxy[v]);
        pv.alive = true;
        pv.version = 0;

        remap[v.idx()] = static_cast<int>(out.vertices.size());
        out.vertices.push_back(pv);
    }

    for (pmp::Face f : he.mesh.faces()) {
        if (!valid_face(he, f)) continue;

        std::array<int, 3> tri{};
        int k = 0;
        for (pmp::Vertex v : he.mesh.vertices(f)) {
            if (k < 3) tri[k++] = remap[v.idx()];
        }

        if (
            k == 3 &&
            tri[0] >= 0 && tri[1] >= 0 && tri[2] >= 0 &&
            tri[0] != tri[1] &&
            tri[1] != tri[2] &&
            tri[2] != tri[0]
        ) {
            Face pf;
            pf.v = tri;
            pf.alive = true;
            out.faces.push_back(pf);
        }
    }

    out.rebuild_vertex_face_adjacency();
    return out;
}

std::array<pmp::Vertex, 3> face_vertices3(const HeMesh& he, pmp::Face f) {
    std::array<pmp::Vertex, 3> vs{};
    int k = 0;
    for (pmp::Vertex v : he.mesh.vertices(f)) {
        if (k < 3) vs[k++] = v;
    }
    if (k != 3) return {};
    return vs;
}

bool face_contains(const HeMesh& he, pmp::Face f, pmp::Vertex v) {
    if (!valid_face(he, f) || !valid_vertex(he, v)) return false;
    for (pmp::Vertex fv : he.mesh.vertices(f)) {
        if (fv == v) return true;
    }
    return false;
}

pmp::Vertex opposite_vertex(const HeMesh& he, pmp::Face f, pmp::Vertex a, pmp::Vertex b) {
    if (!valid_face(he, f)) return pmp::Vertex();
    for (pmp::Vertex v : he.mesh.vertices(f)) {
        if (v != a && v != b) return v;
    }
    return pmp::Vertex();
}

std::array<Vec2, 3> unwrapped_uvs(const HeMesh& he, pmp::Vertex a, pmp::Vertex b, pmp::Vertex c) {
    std::array<Vec2, 3> uv{};
    uv[0] = he.uv[a];
    uv[1] = unwrap_near(he.uv[b], uv[0]);
    uv[2] = unwrap_near(he.uv[c], uv[0]);
    return uv;
}

std::array<Vec2, 3> unwrapped_uvs_for_face(const HeMesh& he, pmp::Face f) {
    const auto vs = face_vertices3(he, f);
    return unwrapped_uvs(he, vs[0], vs[1], vs[2]);
}

double signed_area(const std::array<Vec2, 3>& uv) {
    const Vec2 e0 = uv[1] - uv[0];
    const Vec2 e1 = uv[2] - uv[0];
    return e0.x * e1.y - e0.y * e1.x;
}

double orient2d(Vec2 a, Vec2 b, Vec2 c) {
    const Vec2 e0 = b - a;
    const Vec2 e1 = c - a;
    return e0.x * e1.y - e0.y * e1.x;
}

struct LocalFlipPatchUv {
    Vec2 a;
    Vec2 b;
    Vec2 c;
    Vec2 d;
};

LocalFlipPatchUv unwrap_flip_patch_uv(
    const HeMesh& he,
    pmp::Vertex a,
    pmp::Vertex b,
    pmp::Vertex c,
    pmp::Vertex d
) {
    LocalFlipPatchUv p;

    p.a = he.uv[a];

    // Put the current diagonal endpoint b in the same local periodic chart as a.
    p.b = unwrap_near(he.uv[b], p.a);

    // Opposite vertices should be near the current edge, not necessarily near a only.
    const Vec2 edge_mid = (p.a + p.b) * 0.5;
    p.c = unwrap_near(he.uv[c], edge_mid);
    p.d = unwrap_near(he.uv[d], edge_mid);

    return p;
}

bool has_strict_opposite_signs(double x, double y, double eps) {
    return (x > eps && y < -eps) || (x < -eps && y > eps);
}

// Lightweight UV triangle degeneracy guard.
// These thresholds are intentionally weak: they reject zero-area / vertex-on-edge
// cases while still allowing the highly anisotropic, thin brush triangles.
constexpr double kPnmUvArea2AbsEps = 1e-10;
constexpr double kPnmUvMinSinAngle = 1e-7;
constexpr double kPnmUvLen2Eps = 1e-24;

double pnm_quality_norm2(Vec2 v) {
    return v.x * v.x + v.y * v.y;
}

double pnm_quality_cross(Vec2 a, Vec2 b) {
    return a.x * b.y - a.y * b.x;
}

bool pnm_triangle_quality_ok_unwrapped(
    const std::array<Vec2, 3>& uv,
    double expected_sign
) {
    const double area2 = signed_area(uv);
    if (!std::isfinite(area2)) return false;

    const double signed_area2 = expected_sign * area2;
    if (signed_area2 <= kPnmUvArea2AbsEps) {
        return false;
    }

    // Reject near-collinear triangles.  This is effectively a minimum sine-angle
    // test, but with a very small threshold so that anisotropic brush triangles
    // are still allowed.
    for (int i = 0; i < 3; ++i) {
        const Vec2 a = uv[(i + 1) % 3] - uv[i];
        const Vec2 b = uv[(i + 2) % 3] - uv[i];

        const double la2 = pnm_quality_norm2(a);
        const double lb2 = pnm_quality_norm2(b);
        if (la2 <= kPnmUvLen2Eps || lb2 <= kPnmUvLen2Eps) {
            return false;
        }

        const double denom = std::sqrt(la2 * lb2);
        if (denom <= 0.0 || !std::isfinite(denom)) {
            return false;
        }

        const double sin_angle = std::abs(pnm_quality_cross(a, b)) / denom;
        if (sin_angle <= kPnmUvMinSinAngle) {
            return false;
        }
    }

    return true;
}

bool pnm_triangle_quality_preserve_positive_orientation(
    const std::array<Vec2, 3>& old_uv,
    const std::array<Vec2, 3>& new_uv
) {
    const double old_area2 = signed_area(old_uv);
    if (!std::isfinite(old_area2)) return false;

    // Existing local chart must already be positively oriented.
    if (old_area2 <= kPnmUvArea2AbsEps) return false;

    // New candidate must also be positively oriented and non-degenerate.
    return pnm_triangle_quality_ok_unwrapped(new_uv, +1.0);
}

bool flip_candidate_triangle_quality_ok(
    const HeMesh& he,
    pmp::Vertex a,
    pmp::Vertex b,
    pmp::Vertex c,
    pmp::Vertex d
) {
    const LocalFlipPatchUv p = unwrap_flip_patch_uv(he, a, b, c, d);

    // evaluate_flip() constructs the post-flip faces as {c,d,a} and {d,c,b}.
    // Require positive UV orientation in the local unwrapped chart.
    const std::array<Vec2, 3> t0{p.c, p.d, p.a};
    const std::array<Vec2, 3> t1{p.d, p.c, p.b};

    if (!pnm_triangle_quality_ok_unwrapped(t0, +1.0)) return false;
    if (!pnm_triangle_quality_ok_unwrapped(t1, +1.0)) return false;

    return true;
}

bool is_periodic_uv_convex_flip_patch(
    const HeMesh& he,
    pmp::Vertex a,
    pmp::Vertex b,
    pmp::Vertex c,
    pmp::Vertex d
) {
    const LocalFlipPatchUv p = unwrap_flip_patch_uv(he, a, b, c, d);

    // Scale-independent small tolerance for UV-domain orientation tests.
    // The ordinary grid triangle has area2 around 1 / (W * H), so 1e-18 is only a degeneracy guard.
    constexpr double eps = kPnmUvArea2AbsEps;

    // Current diagonal is a-b. The two opposite vertices must lie on opposite sides.
    const double side_c_ab = orient2d(p.a, p.b, p.c);
    const double side_d_ab = orient2d(p.a, p.b, p.d);
    if (!has_strict_opposite_signs(side_c_ab, side_d_ab, eps)) {
        return false;
    }

    // New diagonal is c-d. The old diagonal endpoints must lie on opposite sides.
    const double side_a_cd = orient2d(p.c, p.d, p.a);
    const double side_b_cd = orient2d(p.c, p.d, p.b);
    if (!has_strict_opposite_signs(side_a_cd, side_b_cd, eps)) {
        return false;
    }

    // The two flipped triangles used by evaluate_flip() are {c,d,a} and {d,c,b}.
    // They should be non-degenerate and have the same orientation in the local chart.
    const double area0 = orient2d(p.c, p.d, p.a);
    const double area1 = orient2d(p.d, p.c, p.b);
    if (std::abs(area0) <= eps || std::abs(area1) <= eps) {
        return false;
    }
    if (area0 * area1 <= 0.0) {
        return false;
    }

    return true;
}

Vec2 lerp_nxy(const HeMesh& he, pmp::Vertex a, pmp::Vertex b, pmp::Vertex c,
              double l0, double l1, double l2) {
    return clamp_projected_normal(he.nxy[a] * l0 + he.nxy[b] * l1 + he.nxy[c] * l2);
}

double projected_error2(Vec2 a, Vec2 b) {
    a = clamp_projected_normal(a);
    b = clamp_projected_normal(b);
    const double dx = a.x - b.x;
    const double dy = a.y - b.y;
    return dx * dx + dy * dy;
}

double triangle_energy(const HeMesh& he, const NormalMap& normal_map,
                       pmp::Vertex a, pmp::Vertex b, pmp::Vertex c) {
    if (!valid_vertex(he, a) || !valid_vertex(he, b) || !valid_vertex(he, c)) return std::numeric_limits<double>::infinity();
    if (a == b || b == c || c == a) return std::numeric_limits<double>::infinity();
    const auto uv = unwrapped_uvs(he, a, b, c);
    const double area2 = std::abs(signed_area(uv));
    if (area2 <= kPnmUvArea2AbsEps) return std::numeric_limits<double>::infinity();

    static constexpr std::array<std::array<double, 3>, 7> samples{{
        {{1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0}},
        {{0.60, 0.20, 0.20}}, {{0.20, 0.60, 0.20}}, {{0.20, 0.20, 0.60}},
        {{0.80, 0.10, 0.10}}, {{0.10, 0.80, 0.10}}, {{0.10, 0.10, 0.80}},
    }};

    double e = 0.0;
    for (const auto& s : samples) {
        const Vec2 u = uv[0] * s[0] + uv[1] * s[1] + uv[2] * s[2];
        const Vec2 nhat = lerp_nxy(he, a, b, c, s[0], s[1], s[2]);
        const Vec2 gt = normal_map.sample_periodic(u);
        e += projected_error2(nhat, gt);
    }
    return area2 * e / static_cast<double>(samples.size());
}

double face_energy(const HeMesh& he, const NormalMap& normal_map, pmp::Face f) {
    if (!valid_face(he, f)) return 0.0;
    const auto vs = face_vertices3(he, f);
    return triangle_energy(he, normal_map, vs[0], vs[1], vs[2]);
}

double pair_energy(const HeMesh& he, const NormalMap& normal_map,
                   const std::array<pmp::Vertex, 3>& t0,
                   const std::array<pmp::Vertex, 3>& t1) {
    const double e0 = triangle_energy(he, normal_map, t0[0], t0[1], t0[2]);
    const double e1 = triangle_energy(he, normal_map, t1[0], t1[1], t1[2]);
    if (!std::isfinite(e0) || !std::isfinite(e1)) return std::numeric_limits<double>::infinity();
    return e0 + e1;
}

std::vector<pmp::Edge> local_edges_around_vertices(const HeMesh& he, const std::vector<pmp::Vertex>& vertices) {
    std::vector<pmp::Edge> edges;
    for (pmp::Vertex v : vertices) {
        if (!valid_vertex(he, v)) continue;
        for (pmp::Edge e : he.mesh.edges(v)) {
            if (valid_edge(he, e)) edges.push_back(e);
        }
    }
    std::sort(edges.begin(), edges.end(), [](pmp::Edge a, pmp::Edge b) { return a.idx() < b.idx(); });
    edges.erase(std::unique(edges.begin(), edges.end()), edges.end());
    return edges;
}

struct FlipCandidate {
    bool ok = false;
    double delta = 0.0;
    std::vector<pmp::Vertex> touched;
};

FlipCandidate evaluate_flip(const HeMesh& he, const NormalMap& normal_map, pmp::Edge e, double improvement_eps) {
    FlipCandidate out;
    if (!valid_edge(he, e)) return out;
    if (he.mesh.is_boundary(e)) return out;
    if (!he.mesh.is_flip_ok(e)) return out;

    const pmp::Vertex a = he.mesh.vertex(e, 0);
    const pmp::Vertex b = he.mesh.vertex(e, 1);
    if (!valid_vertex(he, a) || !valid_vertex(he, b)) return out;

    const pmp::Face f0 = he.mesh.face(e, 0);
    const pmp::Face f1 = he.mesh.face(e, 1);
    if (!valid_face(he, f0) || !valid_face(he, f1)) return out;
    if (!face_contains(he, f0, a) || !face_contains(he, f0, b)) return out;
    if (!face_contains(he, f1, a) || !face_contains(he, f1, b)) return out;

    const pmp::Vertex c = opposite_vertex(he, f0, a, b);
    const pmp::Vertex d = opposite_vertex(he, f1, a, b);
    if (!valid_vertex(he, c) || !valid_vertex(he, d)) return out;
    if (a == b || a == c || a == d || b == c || b == d || c == d) return out;

    const pmp::Edge cd = he.mesh.find_edge(c, d);
    if (he.mesh.is_valid(cd) && !he.mesh.is_deleted(cd)) return out;

    if (!is_periodic_uv_convex_flip_patch(he, a, b, c, d)) return out;

    // New local degeneracy guard: reject flips that create zero-area,
    // near-collinear, or wrong-orientation triangles in periodic UV.
    if (!flip_candidate_triangle_quality_ok(he, a, b, c, d)) return out;

    const double before = face_energy(he, normal_map, f0) + face_energy(he, normal_map, f1);
    const std::array<pmp::Vertex, 3> t0{c, d, a};
    const std::array<pmp::Vertex, 3> t1{d, c, b};
    const double after = pair_energy(he, normal_map, t0, t1);
    if (!std::isfinite(before) || !std::isfinite(after)) return out;
    if (after + improvement_eps >= before) return out;

    out.ok = true;
    out.delta = after - before;
    out.touched = {a, b, c, d};
    return out;
}

struct FlipItem {
    double delta = 0.0;
    pmp::Edge edge;
    int generation = 0;
    bool operator<(const FlipItem& rhs) const {
        return delta > rhs.delta; // min-heap behavior through std::priority_queue
    }
};

bool current_face_triangle_quality_ok(const HeMesh& he, pmp::Face f) {
    if (!valid_face(he, f)) return false;

    const auto uv = unwrapped_uvs_for_face(he, f);

    // Audit と同じく、最終 mesh は positive orientation を要求する。
    return pnm_triangle_quality_ok_unwrapped(uv, +1.0);
}

bool edge_adjacent_faces_quality_ok(const HeMesh& he, pmp::Edge e) {
    if (!valid_edge(he, e)) return false;
    if (he.mesh.is_boundary(e)) return false;

    const pmp::Face f0 = he.mesh.face(e, 0);
    const pmp::Face f1 = he.mesh.face(e, 1);

    if (!valid_face(he, f0) || !valid_face(he, f1)) return false;

    if (!current_face_triangle_quality_ok(he, f0)) return false;
    if (!current_face_triangle_quality_ok(he, f1)) return false;

    return true;
}

std::vector<pmp::Face> incident_faces(const HeMesh& he, pmp::Vertex v) {
    std::vector<pmp::Face> out;
    if (!valid_vertex(he, v)) return out;
    for (pmp::Face f : he.mesh.faces(v)) {
        if (valid_face(he, f)) out.push_back(f);
    }
    return out;
}

bool touched_vertex_star_triangle_quality_ok(
    const HeMesh& he,
    const std::vector<pmp::Vertex>& touched
) {
    std::vector<pmp::Face> faces;

    for (pmp::Vertex v : touched) {
        if (!valid_vertex(he, v)) continue;

        for (pmp::Face f : incident_faces(he, v)) {
            if (valid_face(he, f)) {
                faces.push_back(f);
            }
        }
    }

    std::sort(faces.begin(), faces.end(), [](pmp::Face a, pmp::Face b) {
        return a.idx() < b.idx();
    });
    faces.erase(std::unique(faces.begin(), faces.end()), faces.end());

    for (pmp::Face f : faces) {
        if (!current_face_triangle_quality_ok(he, f)) {
            return false;
        }
    }

    return true;
}

int flip_pass_dirty_queue(HeMesh& he, const NormalMap& normal_map,
                          double improvement_eps,
                          int max_accepted) {
    std::vector<int> generation(std::max<std::size_t>(1, he.mesh.edges_size()), 0);
    std::priority_queue<FlipItem> queue;

    auto ensure_generation = [&](pmp::Edge e) {
        if (e.idx() >= generation.size()) generation.resize(e.idx() + 1, 0);
    };

    auto push_if_good = [&](pmp::Edge e) {
        if (!valid_edge(he, e)) return;
        ensure_generation(e);
        const FlipCandidate cand = evaluate_flip(he, normal_map, e, improvement_eps);
        if (cand.ok) queue.push({cand.delta, e, generation[e.idx()]});
    };

    for (pmp::Edge e : he.mesh.edges()) push_if_good(e);

    int accepted = 0;
    while (!queue.empty()) {
        const std::size_t queue_limit = std::max<std::size_t>(1000000ull, 128ull * std::max<std::size_t>(1, he.mesh.edges_size()));
        if (queue.size() > queue_limit) {
            std::cerr << "WARNING: heFlip queue explosion: queue=" << queue.size()
                      << " edges=" << he.mesh.edges_size()
                      << ". Stopping this flip pass.\n";
            break;
        }
        const FlipItem item = queue.top();
        queue.pop();
        pmp::Edge e = item.edge;
        if (!valid_edge(he, e)) continue;
        ensure_generation(e);
        if (item.generation != generation[e.idx()]) continue;

        const FlipCandidate cand = evaluate_flip(he, normal_map, e, improvement_eps);
        if (!cand.ok) continue;

        const pmp::Vertex c = cand.touched.size() >= 3 ? cand.touched[2] : pmp::Vertex();
        const pmp::Vertex d = cand.touched.size() >= 4 ? cand.touched[3] : pmp::Vertex();

        he.mesh.flip(e);

        // Validate the actual post-flip local star, not only the predicted two triangles.
        if (!touched_vertex_star_triangle_quality_ok(he, cand.touched)) {
            pmp::Edge undo_e = pmp::Edge();
        
            if (valid_vertex(he, c) && valid_vertex(he, d)) {
                undo_e = he.mesh.find_edge(c, d);
            }
        
            if (!valid_edge(he, undo_e)) {
                undo_e = e;
            }
        
            if (valid_edge(he, undo_e) && he.mesh.is_flip_ok(undo_e)) {
                he.mesh.flip(undo_e);
            } else {
                std::cerr << "WARNING: could not undo invalid flip around edge="
                          << e.idx() << "\n";
            }
        
            continue;
        }

        ++accepted;

        if (accepted > 0 && accepted % 1000 == 0) {
        std::cerr << "heFlip accepted=" << accepted
              << " queue=" << queue.size()
              << " edges=" << he.mesh.edges_size()
              << "\n";
        }

        const auto dirty = local_edges_around_vertices(he, cand.touched);
        for (pmp::Edge de : dirty) {
            if (!valid_edge(he, de)) continue;
            ensure_generation(de);
            ++generation[de.idx()];
            push_if_good(de);
        }

        if (max_accepted > 0 && accepted >= max_accepted) break;
    }

    return accepted;
}

std::vector<pmp::Vertex> live_neighbors(const HeMesh& he, pmp::Vertex v) {
    std::vector<pmp::Vertex> out;
    if (!valid_vertex(he, v)) return out;
    for (pmp::Vertex w : he.mesh.vertices(v)) {
        if (valid_vertex(he, w)) out.push_back(w);
    }
    std::sort(out.begin(), out.end(), [](pmp::Vertex a, pmp::Vertex b) { return a.idx() < b.idx(); });
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

std::array<Vec2, 3> relocation_trial_face_uvs(
    const HeMesh& he,
    pmp::Face f,
    pmp::Vertex moving,
    Vec2 trial_uv
) {
    const auto vs = face_vertices3(he, f);

    std::array<Vec2, 3> uv{};

    // Important:
    // Trial check must mirror the final stored mesh and audit convention.
    // The actual accepted position is stored as fract(best_uv), and
    // unwrapped_uvs_for_face() unwraps relative to face vertex 0.
    const Vec2 trial_wrapped = fract(trial_uv);

    const Vec2 p0 = (vs[0] == moving) ? trial_wrapped : he.uv[vs[0]];
    uv[0] = p0;

    for (int i = 1; i < 3; ++i) {
        const Vec2 pi = (vs[i] == moving) ? trial_wrapped : he.uv[vs[i]];
        uv[i] = unwrap_near(pi, uv[0]);
    }

    return uv;
}

bool vertex_star_triangle_quality_ok(const HeMesh& he, pmp::Vertex v) {
    if (!valid_vertex(he, v)) return false;

    for (pmp::Face f : incident_faces(he, v)) {
        if (!valid_face(he, f)) continue;

        if (!current_face_triangle_quality_ok(he, f)) {
            return false;
        }
    }

    return true;
}

bool relocation_trial_triangle_quality_ok(
    HeMesh& he,
    pmp::Vertex moving,
    Vec2 trial_uv
) {
    if (!valid_vertex(he, moving)) return false;

    const Vec2 old_uv = he.uv[moving];
    const pmp::Point old_p = he.mesh.position(moving);

    // Mirror the actual accepted update exactly.
    he.uv[moving] = fract(trial_uv);
    he.mesh.position(moving) = point_from_uv(he.uv[moving]);

    const bool ok = vertex_star_triangle_quality_ok(he, moving);

    he.uv[moving] = old_uv;
    he.mesh.position(moving) = old_p;

    return ok;
}

double vertex_local_energy(const HeMesh& he, const NormalMap& normal_map, pmp::Vertex v) {
    double e = 0.0;
    for (pmp::Face f : incident_faces(he, v)) e += face_energy(he, normal_map, f);
    return e;
}

double vertex_local_energy_trial(HeMesh& he, const NormalMap& normal_map,
                                 pmp::Vertex v, Vec2 trial_uv, Vec2 trial_nxy) {
    const Vec2 old_uv = he.uv[v];
    const Vec2 old_n = he.nxy[v];
    const pmp::Point old_p = he.mesh.position(v);

    he.uv[v] = fract(trial_uv);
    he.nxy[v] = clamp_projected_normal(trial_nxy);
    he.mesh.position(v) = point_from_uv(he.uv[v]);
    const double e = vertex_local_energy(he, normal_map, v);

    he.uv[v] = old_uv;
    he.nxy[v] = old_n;
    he.mesh.position(v) = old_p;
    return e;
}

double area_with_vertex_delta(const HeMesh& he, pmp::Face f, pmp::Vertex moving, Vec2 delta, double alpha) {
    const auto vs = face_vertices3(he, f);
    Vec2 base = he.uv[moving];
    std::array<Vec2, 3> uv{};
    for (int i = 0; i < 3; ++i) {
        if (vs[i] == moving) uv[i] = base + delta * alpha;
        else uv[i] = unwrap_near(he.uv[vs[i]], base);
    }
    return signed_area(uv);
}

double orientation_root_limited_alpha(const HeMesh& he, pmp::Vertex v, Vec2 delta, double min_area_ratio) {
    double alpha = 1.0;
    for (pmp::Face f : incident_faces(he, v)) {
        const double c = area_with_vertex_delta(he, f, v, delta, 0.0);
        if (std::abs(c) < 1e-20) return 0.0;
        const double a1 = area_with_vertex_delta(he, f, v, delta, 1.0);
        const double m = a1 - c;
        const double sign = (c >= 0.0) ? 1.0 : -1.0;
        const double signed_slope = sign * m;
        const double abs_c = std::abs(c);
        const double min_abs = std::max(0.0, min_area_ratio) * abs_c;

        if (signed_slope < 0.0) {
            // sign*(c + alpha*m) >= min_abs
            const double max_alpha = (abs_c - min_abs) / (-signed_slope);
            if (max_alpha <= 0.0 || !std::isfinite(max_alpha)) return 0.0;
            alpha = std::min(alpha, 0.75 * max_alpha);
        }
    }
    return std::max(0.0, std::min(1.0, alpha));
}

int relocate_pass_root_limited(HeMesh& he, const NormalMap& normal_map,
                               double step_pixels,
                               double min_area_ratio) {
    const double sx = step_pixels / std::max(1, he.width);
    const double sy = step_pixels / std::max(1, he.height);
    int moved = 0;

    std::vector<pmp::Vertex> verts;
    for (pmp::Vertex v : he.mesh.vertices()) if (valid_vertex(he, v)) verts.push_back(v);

    for (pmp::Vertex v : verts) {
        if (!valid_vertex(he, v)) continue;
        const Vec2 base = he.uv[v];
        const double base_e = vertex_local_energy(he, normal_map, v);
        double best_e = base_e;
        Vec2 best_uv = base;

        std::vector<Vec2> candidates;
        candidates.push_back({base.x + sx, base.y});
        candidates.push_back({base.x - sx, base.y});
        candidates.push_back({base.x, base.y + sy});
        candidates.push_back({base.x, base.y - sy});
        candidates.push_back({base.x + sx, base.y + sy});
        candidates.push_back({base.x + sx, base.y - sy});
        candidates.push_back({base.x - sx, base.y + sy});
        candidates.push_back({base.x - sx, base.y - sy});

        const auto nbrs = live_neighbors(he, v);
        if (!nbrs.empty()) {
            Vec2 centroid{0.0, 0.0};
            for (pmp::Vertex w : nbrs) centroid += unwrap_near(he.uv[w], base);
            centroid = centroid / static_cast<double>(nbrs.size());
            candidates.push_back(base + (centroid - base) * 0.5);
        }

        for (Vec2 cand : candidates) {
            const Vec2 delta = cand - base;
            if (norm2(delta) < 1e-24) continue;
            const double alpha = orientation_root_limited_alpha(he, v, delta, min_area_ratio);
            if (alpha <= 1e-8) continue;
            const Vec2 trial = base + delta * alpha;

            // New local degeneracy guard.
            // Orientation-root limiting prevents crossing the zero-orientation root,
            // but it can still allow moves very close to zero area.  This rejects
            // zero-area / near-collinear incident triangles.
            if (!relocation_trial_triangle_quality_ok(he, v, trial)) continue;

            const Vec2 trial_n = normal_map.sample_periodic(trial);
            const double e = vertex_local_energy_trial(he, normal_map, v, trial, trial_n);

            if (std::isfinite(e) && e + 1e-15 < best_e) {
                best_e = e;
                best_uv = trial;
            }
        }

        if (norm2(wrap_delta(best_uv - base)) > 1e-20) {
            he.uv[v] = fract(best_uv);
            he.nxy[v] = normal_map.sample_periodic(best_uv);
            he.mesh.position(v) = point_from_uv(he.uv[v]);
            ++moved;
        }
    }

    return moved;
}

std::vector<int> make_stage_targets(int initial_vertices, int final_target, int iterations) {
    std::vector<int> stages;
    if (final_target <= 0 || final_target >= initial_vertices) return stages;
    iterations = std::max(1, iterations);
    const int start = std::min(initial_vertices, std::max(final_target, final_target * 10));
    if (start < initial_vertices) stages.push_back(start);
    for (int i = 1; i <= iterations; ++i) {
        const double t = static_cast<double>(i) / static_cast<double>(iterations);
        const double v = std::exp((1.0 - t) * std::log(static_cast<double>(start)) +
                                  t * std::log(static_cast<double>(final_target)));
        const int stage = std::max(final_target, static_cast<int>(std::llround(v)));
        if (stages.empty() || stage < stages.back()) stages.push_back(stage);
    }
    if (stages.empty() || stages.back() != final_target) stages.push_back(final_target);
    return stages;
}

bool checkpoint_enabled(const Options& opt) {
    return !opt.checkpoint_dir.empty();
}

std::filesystem::path checkpoint_path(
    const Options& opt,
    int stage,
    const std::string& tag
) {
    return std::filesystem::path(opt.checkpoint_dir) /
        (opt.checkpoint_prefix +
         "_stage" + std::to_string(stage) +
         "_" + tag +
         ".meshbin");
}

void write_torus_checkpoint(
    const Options& opt,
    const TorusMesh& mesh,
    int stage,
    const std::string& tag
) {
    if (!checkpoint_enabled(opt)) return;

    std::filesystem::create_directories(opt.checkpoint_dir);

    const auto path = checkpoint_path(opt, stage, tag);
    write_mesh_binary(mesh, path.string());

    if (!opt.quiet) {
        std::cerr << "checkpoint stage=" << stage
                  << " tag=" << tag
                  << " path=" << path.string()
                  << "\n";
    }
}

void write_he_checkpoint(
    const Options& opt,
    const HeMesh& he,
    int stage,
    const std::string& tag
) {
    if (!checkpoint_enabled(opt)) return;

    TorusMesh snap = snapshot_halfedge_mesh(he);
    write_torus_checkpoint(opt, snap, stage, tag);
}

Options parse(int argc, char** argv) {
    Options opt;
    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        auto next = [&]() -> std::string {
            if (i + 1 >= argc) {
                usage();
                std::exit(2);
            }
            return argv[++i];
        };
        if (a == "--input") opt.input = next();
        else if (a == "--output") opt.output = next();
        else if (a == "--target") opt.target = std::stoi(next());
        else if (a == "--mode") {
            const std::string m = next();
            if (m == "normal") opt.mode = QemMode::NormalOnly;
            else if (m == "4d") opt.mode = QemMode::Qem4D;
            else { usage(); std::exit(2); }
        } else if (a == "--iterations") opt.iterations = std::stoi(next());
        else if (a == "--flip-passes") opt.flip_passes = std::stoi(next());
        else if (a == "--relocate-passes") opt.relocate_passes = std::stoi(next());
        else if (a == "--relocate-step-pixels") opt.relocate_step_pixels = std::stod(next());
        else if (a == "--min-area-ratio") opt.min_area_ratio = std::stod(next());
        else if (a == "--flip-improvement-eps") opt.flip_improvement_eps = std::stod(next());
        else if (a == "--flip-max-accepted-per-pass") opt.flip_max_accepted_per_pass = std::stoi(next());
        else if (a == "--rebuild-interval") opt.rebuild_interval = std::stoi(next());
        else if (a == "--progress-interval") opt.progress_interval = std::stoi(next());
        else if (a == "--quiet") opt.quiet = true;
        else if (a == "--checkpoint-dir") opt.checkpoint_dir = next();
        else if (a == "--checkpoint-prefix") opt.checkpoint_prefix = next();
        else { usage(); std::exit(2); }
    }
    if (opt.input.empty() || opt.output.empty() || opt.target <= 0) {
        usage();
        std::exit(2);
    }
    return opt;
}

} // namespace

int main(int argc, char** argv) {
    try {
        Options opt = parse(argc, argv);
        NormalMap map = read_nxy_binary(opt.input);

        TorusMesh mesh;
        mesh.build_full_torus(map);
        const int initial_vertices = mesh.alive_vertex_count();
        const std::vector<int> stages = make_stage_targets(initial_vertices, opt.target, opt.iterations);

        if (!opt.quiet) {
            std::cerr << "pnm_aniso_he stages:";
            for (int s : stages) std::cerr << ' ' << s;
            std::cerr << "\n";
        }

        for (int stage : stages) {
            if (mesh.alive_vertex_count() > stage) {
                CollapseOptions copt;
                copt.target_vertices = stage;
                copt.mode = opt.mode;
                copt.verbose = !opt.quiet;
                copt.rebuild_interval = opt.rebuild_interval;
                copt.progress_interval = opt.progress_interval;
            
                auto stats = simplify_qem(mesh, copt);
            
                if (!opt.quiet) {
                    std::cerr << "stage collapse target=" << stage
                              << " finalV=" << stats.final_vertices
                              << " sec=" << stats.total_seconds << "\n";
                }
            }
        
            // Checkpoint after QEM collapse, before conversion to half-edge.
            write_torus_checkpoint(opt, mesh, stage, "after_collapse");
        
            if (opt.flip_passes > 0 || opt.relocate_passes > 0) {
                HeMesh he = to_halfedge_mesh(mesh);
            
                // Optional: verify whether conversion itself already changes anything.
                write_he_checkpoint(opt, he, stage, "after_to_halfedge");
            
                for (int p = 0; p < opt.flip_passes; ++p) {
                    const int n = flip_pass_dirty_queue(
                        he,
                        map,
                        opt.flip_improvement_eps,
                        opt.flip_max_accepted_per_pass
                    );
                
                    if (!opt.quiet) {
                        std::cerr << "stage=" << stage
                                  << " heFlipPass=" << p
                                  << " flipped=" << n << "\n";
                    }
                
                    write_he_checkpoint(
                        opt,
                        he,
                        stage,
                        "after_flip_p" + std::to_string(p)
                    );
                }
            
                for (int p = 0; p < opt.relocate_passes; ++p) {
                    const int n = relocate_pass_root_limited(
                        he,
                        map,
                        opt.relocate_step_pixels,
                        opt.min_area_ratio
                    );
                
                    if (!opt.quiet) {
                        std::cerr << "stage=" << stage
                                  << " heRelocatePass=" << p
                                  << " moved=" << n << "\n";
                    }
                
                    write_he_checkpoint(
                        opt,
                        he,
                        stage,
                        "after_reloc_p" + std::to_string(p)
                    );
                }
            
                mesh = from_halfedge_mesh(he);
            
                // Checkpoint after HE -> TorusMesh conversion.
                write_torus_checkpoint(opt, mesh, stage, "after_from_halfedge");
            }
        }

        const std::filesystem::path out_path(opt.output);
        if (out_path.parent_path() != std::filesystem::path{}) {
            std::filesystem::create_directories(out_path.parent_path());
        }
        write_mesh_binary(mesh, opt.output);
        std::cout << "final_vertices," << mesh.alive_vertex_count() << "\n";
        std::cout << "final_faces," << mesh.alive_face_count() << "\n";
        std::cout << "output," << opt.output << "\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
