#include "pndf/edge_collapse.hpp"
#include <algorithm>
#include <chrono>
#include <cstdint>
#include <iostream>
#include <set>
#include <unordered_set>

namespace pndf {

namespace {

struct QueueItem {
    double cost = 0.0;
    int u = -1;
    int v = -1;
    int version_u = 0;
    int version_v = 0;
    bool operator<(const QueueItem& b) const { return cost > b.cost; }
};

uint64_t edge_key(int a, int b) {
    if (a > b) std::swap(a,b);
    return (uint64_t(uint32_t(a)) << 32) | uint32_t(b);
}

bool face_has(const Face& f, int v) { return f.v[0] == v || f.v[1] == v || f.v[2] == v; }

bool shared_live_face(const TorusMesh& mesh, int u, int v) {
    const auto& fu = mesh.vertices[u].faces;
    const auto& fv = mesh.vertices[v].faces;
    const auto& small = fu.size() < fv.size() ? fu : fv;
    for (int fi : small) {
        const Face& f = mesh.faces[fi];
        if (!f.alive) continue;
        if (face_has(f, u) && face_has(f, v)) return true;
    }
    return false;
}

std::vector<int> live_neighbors(const TorusMesh& mesh, int u) {
    std::vector<int> out;
    for (int fi : mesh.vertices[u].faces) {
        const Face& f = mesh.faces[fi];
        if (!f.alive) continue;
        for (int k=0;k<3;++k) {
            int w = f.v[k];
            if (w != u && mesh.vertices[w].alive) out.push_back(w);
        }
    }
    std::sort(out.begin(), out.end());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    return out;
}

bool link_condition_closed_manifold(const TorusMesh& mesh, int u, int v) {
    auto a = live_neighbors(mesh, u);
    auto b = live_neighbors(mesh, v);
    std::vector<int> inter;
    std::set_intersection(a.begin(), a.end(), b.begin(), b.end(), std::back_inserter(inter));
    // In a closed triangle manifold, an edge contraction is safe only if exactly the two opposite vertices are shared.
    return inter.size() == 2;
}

bool candidate_non_flipping(const TorusMesh& mesh, int keep, int remove, Vec2 new_uv) {
    std::vector<int> affected = mesh.vertices[keep].faces;
    affected.insert(affected.end(), mesh.vertices[remove].faces.begin(), mesh.vertices[remove].faces.end());
    std::sort(affected.begin(), affected.end());
    affected.erase(std::unique(affected.begin(), affected.end()), affected.end());
    for (int fi : affected) {
        const Face& f = mesh.faces[fi];
        if (!f.alive) continue;
        const bool hk = face_has(f, keep);
        const bool hr = face_has(f, remove);
        if (hk && hr) continue; // face disappears
        if (!hk && !hr) continue;
        std::array<Vec2,3> olduv = unwrapped_face_uvs(mesh, f);
        const double old_area = oriented_area_unwrapped(olduv);
        if (std::abs(old_area) < 1e-18) return false;
        std::array<Vec2,3> newuv = olduv;
        // unwrap new_uv near first face vertex for a conservative local test
        for (int k=0;k<3;++k) {
            if (f.v[k] == keep || f.v[k] == remove) {
                Vec2 base = olduv[0];
                Vec2 d = wrap_delta(new_uv - base);
                newuv[k] = base + d;
            }
        }
        const double new_area = oriented_area_unwrapped(newuv);
        if (std::abs(new_area) < 1e-18) return false;
        if (old_area * new_area < -1e-18) return false;
    }
    return true;
}

QueueItem make_item(const TorusMesh& mesh, int u, int v, QemMode mode, bool lock_boundary) {
    if (u > v) std::swap(u,v);
    constexpr double kBig = 1e100;
    if (!mesh.vertices[u].alive || !mesh.vertices[v].alive) return {kBig,u,v,0,0};
    if (lock_boundary && (mesh.vertices[u].boundary_locked || mesh.vertices[v].boundary_locked)) return {kBig,u,v,0,0};
    if (!shared_live_face(mesh, u, v)) return {kBig,u,v,0,0};
    if (!link_condition_closed_manifold(mesh, u, v)) return {kBig,u,v,0,0};
    const Mat5 q = mesh.vertices[u].q + mesh.vertices[v].q;
    Placement p = best_placement(q, mesh.vertices[u].uv, mesh.vertices[u].nxy, mesh.vertices[v].uv, mesh.vertices[v].nxy, mode);
    if (!candidate_non_flipping(mesh, u, v, p.uv)) return {kBig,u,v,0,0};
    return {p.cost, u, v, mesh.vertices[u].version, mesh.vertices[v].version};
}

void initialize_quadrics(TorusMesh& mesh, QemMode mode) {
    for (auto& v : mesh.vertices) v.q = point_normal_stabilizer(v.nxy, 1.0);
    for (const Face& f : mesh.faces) {
        if (!f.alive) continue;
        Mat5 q = (mode == QemMode::NormalOnly) ? triangle_quadric_normal_only(mesh, f) : triangle_quadric_4d(mesh, f);
        for (int k=0;k<3;++k) mesh.vertices[f.v[k]].q += q;
    }
}

void rebuild_heap(const TorusMesh& mesh, std::priority_queue<QueueItem>& heap, QemMode mode, bool lock_boundary) {
    std::vector<uint64_t> keys;
    keys.reserve(mesh.faces.size()*3/2);
    for (const Face& f : mesh.faces) {
        if (!f.alive) continue;
        for (int k=0;k<3;++k) {
            int a = f.v[k], b = f.v[(k+1)%3];
            if (mesh.vertices[a].alive && mesh.vertices[b].alive) keys.push_back(edge_key(a,b));
        }
    }
    std::sort(keys.begin(), keys.end());
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
    heap = std::priority_queue<QueueItem>();
    for (uint64_t key : keys) {
        int a = int(key >> 32);
        int b = int(key & 0xffffffffu);
        QueueItem item = make_item(mesh, a, b, mode, lock_boundary);
        if (item.cost < 1e90) heap.push(item);
    }
}

void compact_vertex_faces(TorusMesh& mesh, int u) {
    auto& vf = mesh.vertices[u].faces;
    std::sort(vf.begin(), vf.end());
    vf.erase(std::unique(vf.begin(), vf.end()), vf.end());
    std::vector<int> out;
    for (int fi : vf) {
        if (fi < 0 || fi >= static_cast<int>(mesh.faces.size())) continue;
        const Face& f = mesh.faces[fi];
        if (f.alive && face_has(f, u)) out.push_back(fi);
    }
    vf.swap(out);
}

void collapse_edge(TorusMesh& mesh, int keep, int remove, QemMode mode) {
    Mat5 q = mesh.vertices[keep].q + mesh.vertices[remove].q;
    Placement p = best_placement(q, mesh.vertices[keep].uv, mesh.vertices[keep].nxy, mesh.vertices[remove].uv, mesh.vertices[remove].nxy, mode);
    mesh.vertices[keep].uv = p.uv;
    mesh.vertices[keep].nxy = p.nxy;
    mesh.vertices[keep].q = q;
    mesh.vertices[keep].version++;
    mesh.vertices[remove].alive = false;
    mesh.vertices[remove].version++;

    std::vector<int> affected = mesh.vertices[keep].faces;
    affected.insert(affected.end(), mesh.vertices[remove].faces.begin(), mesh.vertices[remove].faces.end());
    std::sort(affected.begin(), affected.end());
    affected.erase(std::unique(affected.begin(), affected.end()), affected.end());
    for (int fi : affected) {
        Face& f = mesh.faces[fi];
        if (!f.alive) continue;
        const bool hk = face_has(f, keep);
        const bool hr = face_has(f, remove);
        if (hk && hr) {
            f.alive = false;
            continue;
        }
        for (int k=0;k<3;++k) if (f.v[k] == remove) f.v[k] = keep;
        if (f.v[0] == f.v[1] || f.v[1] == f.v[2] || f.v[2] == f.v[0]) {
            f.alive = false;
        } else {
            mesh.vertices[keep].faces.push_back(fi);
        }
    }
    mesh.vertices[remove].faces.clear();
    compact_vertex_faces(mesh, keep);
}

} // namespace

CollapseStats simplify_qem(TorusMesh& mesh, const CollapseOptions& opt) {
    CollapseStats stats;
    initialize_quadrics(mesh, opt.mode);
    std::priority_queue<QueueItem> heap;
    rebuild_heap(mesh, heap, opt.mode, opt.lock_boundary);
    const auto t0 = std::chrono::steady_clock::now();
    int since_rebuild = 0;
    while (mesh.alive_vertex_count() > opt.target_vertices && !heap.empty()) {
        QueueItem item = heap.top(); heap.pop();
        if (item.cost > 1e90) { stats.rejected++; continue; }
        if (!mesh.vertices[item.u].alive || !mesh.vertices[item.v].alive) { stats.rejected++; continue; }
        if (item.version_u != mesh.vertices[item.u].version || item.version_v != mesh.vertices[item.v].version) { stats.rejected++; continue; }
        QueueItem fresh = make_item(mesh, item.u, item.v, opt.mode, opt.lock_boundary);
        if (fresh.cost > 1e90) { stats.rejected++; continue; }
        if (std::abs(fresh.cost - item.cost) > 1e-8 * (1.0 + std::abs(item.cost))) {
            heap.push(fresh);
            continue;
        }
        collapse_edge(mesh, item.u, item.v, opt.mode);
        stats.accepted++;
        since_rebuild++;
        if (since_rebuild >= opt.rebuild_interval) {
            since_rebuild = 0;
            rebuild_heap(mesh, heap, opt.mode, opt.lock_boundary);
            if (opt.verbose) {
                auto t1 = std::chrono::steady_clock::now();
                double sec = std::chrono::duration<double>(t1-t0).count();
                std::cerr << "accepted=" << stats.accepted << " aliveV=" << mesh.alive_vertex_count()
                          << " aliveF=" << mesh.alive_face_count() << " sec=" << sec << "\n";
            }
        }
    }
    stats.final_vertices = mesh.alive_vertex_count();
    stats.final_faces = mesh.alive_face_count();
    return stats;
}

} // namespace pndf
