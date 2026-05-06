#include "pndf/edge_collapse.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <set>
#include <unordered_set>

namespace pndf {

namespace {

using Clock = std::chrono::steady_clock;

static double seconds_since(Clock::time_point t0) {
    return std::chrono::duration<double>(Clock::now() - t0).count();
}

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

QueueItem make_item(const TorusMesh& mesh, int u, int v, QemMode mode, bool lock_boundary,
                    CollapseStats* stats = nullptr) {
    const auto t_make = Clock::now();
    auto finish = [&](QueueItem item) -> QueueItem {
        if (stats) {
            stats->make_item_calls++;
            stats->make_item_seconds += seconds_since(t_make);
        }
        return item;
    };

    if (u > v) std::swap(u,v);
    constexpr double kBig = 1e100;
    if (!mesh.vertices[u].alive || !mesh.vertices[v].alive) return finish({kBig,u,v,0,0});
    if (lock_boundary && (mesh.vertices[u].boundary_locked || mesh.vertices[v].boundary_locked)) return finish({kBig,u,v,0,0});
    if (!shared_live_face(mesh, u, v)) return finish({kBig,u,v,0,0});

    const auto t_link = Clock::now();
    const bool link_ok = link_condition_closed_manifold(mesh, u, v);
    if (stats) stats->link_condition_seconds += seconds_since(t_link);
    if (!link_ok) return finish({kBig,u,v,0,0});

    const Mat5 q = mesh.vertices[u].q + mesh.vertices[v].q;
    const auto t_place = Clock::now();
    Placement p = best_placement(q, mesh.vertices[u].uv, mesh.vertices[u].nxy,
                                 mesh.vertices[v].uv, mesh.vertices[v].nxy, mode);
    if (stats) stats->placement_seconds += seconds_since(t_place);

    const auto t_flip = Clock::now();
    const bool no_flip = candidate_non_flipping(mesh, u, v, p.uv);
    if (stats) stats->non_flip_seconds += seconds_since(t_flip);
    if (!no_flip) return finish({kBig,u,v,0,0});

    return finish({p.cost, u, v, mesh.vertices[u].version, mesh.vertices[v].version});
}

void initialize_quadrics(TorusMesh& mesh, QemMode mode) {
    for (auto& v : mesh.vertices) v.q = point_normal_stabilizer(v.nxy, 1.0);
    for (const Face& f : mesh.faces) {
        if (!f.alive) continue;
        Mat5 q = (mode == QemMode::NormalOnly) ? triangle_quadric_normal_only(mesh, f) : triangle_quadric_4d(mesh, f);
        for (int k=0;k<3;++k) mesh.vertices[f.v[k]].q += q;
    }
}

void rebuild_heap(const TorusMesh& mesh, std::priority_queue<QueueItem>& heap, QemMode mode,
                  bool lock_boundary, CollapseStats* stats = nullptr) {
    const auto t_rebuild = Clock::now();
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
        QueueItem item = make_item(mesh, a, b, mode, lock_boundary, stats);
        if (item.cost < 1e90) heap.push(item);
    }
    if (stats) {
        stats->rebuild_count++;
        stats->rebuild_seconds += seconds_since(t_rebuild);
        stats->max_heap_size = std::max<std::uint64_t>(stats->max_heap_size, heap.size());
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

void collapse_edge(TorusMesh& mesh, int keep, int remove, QemMode mode, CollapseStats* stats = nullptr) {
    const auto t_collapse = Clock::now();
    Mat5 q = mesh.vertices[keep].q + mesh.vertices[remove].q;
    Placement p = best_placement(q, mesh.vertices[keep].uv, mesh.vertices[keep].nxy,
                                 mesh.vertices[remove].uv, mesh.vertices[remove].nxy, mode);
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

    if (stats) {
        stats->collapse_updates++;
        stats->collapse_update_seconds += seconds_since(t_collapse);
    }
}

void push_local_edges(const TorusMesh& mesh, int keep, std::priority_queue<QueueItem>& heap,
                      QemMode mode, bool lock_boundary, CollapseStats* stats = nullptr) {
    std::vector<uint64_t> keys;
    keys.reserve(mesh.vertices[keep].faces.size() * 2);
    for (int fi : mesh.vertices[keep].faces) {
        if (fi < 0 || fi >= static_cast<int>(mesh.faces.size())) continue;
        const Face& f = mesh.faces[fi];
        if (!f.alive) continue;
        if (!face_has(f, keep)) continue;
        for (int k = 0; k < 3; ++k) {
            const int w = f.v[k];
            if (w == keep) continue;
            if (!mesh.vertices[w].alive) continue;
            keys.push_back(edge_key(keep, w));
        }
    }
    std::sort(keys.begin(), keys.end());
    keys.erase(std::unique(keys.begin(), keys.end()), keys.end());
    for (uint64_t key : keys) {
        int a = int(key >> 32);
        int b = int(key & 0xffffffffu);
        QueueItem item = make_item(mesh, a, b, mode, lock_boundary, stats);
        if (item.cost < 1e90) heap.push(item);
    }
}

void write_progress_header(std::ofstream& csv) {
    csv << "reason,accepted,rejected,alive_vertices,alive_faces,heap_size,elapsed_seconds,"
        << "heap_pops,stale_pops,dead_pops,invalid_cost_pops,fresh_recomputes,fresh_requeues,fresh_rejects,"
        << "rebuild_count,max_heap_size,initialize_seconds,rebuild_seconds,make_item_seconds,"
        << "link_condition_seconds,placement_seconds,non_flip_seconds,collapse_update_seconds\n";
}

void emit_progress(const char* reason, const TorusMesh& mesh, const std::priority_queue<QueueItem>& heap,
                   const CollapseStats& stats, Clock::time_point t0, bool verbose, std::ofstream* csv) {
    const double elapsed = seconds_since(t0);
    if (verbose) {
        std::cerr << "reason=" << reason
                  << " accepted=" << stats.accepted
                  << " rejected=" << stats.rejected
                  << " aliveV=" << mesh.alive_vertex_count()
                  << " aliveF=" << mesh.alive_face_count()
                  << " heap=" << heap.size()
                  << " sec=" << elapsed
                  << " rebuildSec=" << stats.rebuild_seconds
                  << " makeItemSec=" << stats.make_item_seconds
                  << " linkSec=" << stats.link_condition_seconds
                  << " placeSec=" << stats.placement_seconds
                  << " nonFlipSec=" << stats.non_flip_seconds
                  << " collapseSec=" << stats.collapse_update_seconds
                  << " heapPops=" << stats.heap_pops
                  << " stale=" << stats.stale_pops
                  << " dead=" << stats.dead_pops
                  << " fresh=" << stats.fresh_recomputes
                  << " requeues=" << stats.fresh_requeues
                  << "\n";
    }
    if (csv && csv->is_open()) {
        (*csv) << reason << ','
               << stats.accepted << ',' << stats.rejected << ','
               << mesh.alive_vertex_count() << ',' << mesh.alive_face_count() << ',' << heap.size() << ','
               << elapsed << ','
               << stats.heap_pops << ',' << stats.stale_pops << ',' << stats.dead_pops << ','
               << stats.invalid_cost_pops << ',' << stats.fresh_recomputes << ',' << stats.fresh_requeues << ','
               << stats.fresh_rejects << ',' << stats.rebuild_count << ',' << stats.max_heap_size << ','
               << stats.initialize_seconds << ',' << stats.rebuild_seconds << ',' << stats.make_item_seconds << ','
               << stats.link_condition_seconds << ',' << stats.placement_seconds << ',' << stats.non_flip_seconds << ','
               << stats.collapse_update_seconds << '\n';
        csv->flush();
    }
}

} // namespace

CollapseStats simplify_qem(TorusMesh& mesh, const CollapseOptions& opt) {
    CollapseStats stats;
    const auto t0 = Clock::now();

    const auto t_init = Clock::now();
    initialize_quadrics(mesh, opt.mode);
    stats.initialize_seconds = seconds_since(t_init);

    std::priority_queue<QueueItem> heap;
    rebuild_heap(mesh, heap, opt.mode, opt.lock_boundary, &stats);

    std::ofstream csv;
    if (!opt.progress_csv.empty()) {
        csv.open(opt.progress_csv);
        if (csv) write_progress_header(csv);
        else std::cerr << "warning: could not open progress csv: " << opt.progress_csv << "\n";
    }

    int since_rebuild = 0;
    emit_progress("initial", mesh, heap, stats, t0, opt.verbose, csv.is_open() ? &csv : nullptr);

    while (mesh.alive_vertex_count() > opt.target_vertices && !heap.empty()) {
        stats.max_heap_size = std::max<std::uint64_t>(stats.max_heap_size, heap.size());
        QueueItem item = heap.top();
        heap.pop();
        stats.heap_pops++;

        if (item.cost > 1e90) {
            stats.rejected++;
            stats.invalid_cost_pops++;
            continue;
        }
        if (!mesh.vertices[item.u].alive || !mesh.vertices[item.v].alive) {
            stats.rejected++;
            stats.dead_pops++;
            continue;
        }
        if (item.version_u != mesh.vertices[item.u].version || item.version_v != mesh.vertices[item.v].version) {
            stats.rejected++;
            stats.stale_pops++;
            continue;
        }

        stats.fresh_recomputes++;
        QueueItem fresh = make_item(mesh, item.u, item.v, opt.mode, opt.lock_boundary, &stats);
        if (fresh.cost > 1e90) {
            stats.rejected++;
            stats.fresh_rejects++;
            continue;
        }
        if (std::abs(fresh.cost - item.cost) > 1e-8 * (1.0 + std::abs(item.cost))) {
            heap.push(fresh);
            stats.fresh_requeues++;
            continue;
        }

        collapse_edge(mesh, item.u, item.v, opt.mode, &stats);
        // Critical for QEM correctness: after a collapse, all edges incident to
        // the surviving vertex have changed quadrics/topology and need fresh
        // queue items.  Periodic full rebuilds should be a cleanup tool, not the
        // only mechanism that updates edge priorities.
        push_local_edges(mesh, item.u, heap, opt.mode, opt.lock_boundary, &stats);
        stats.max_heap_size = std::max<std::uint64_t>(stats.max_heap_size, heap.size());
        stats.accepted++;
        since_rebuild++;

        if (opt.rebuild_interval > 0 && since_rebuild >= opt.rebuild_interval) {
            since_rebuild = 0;
            rebuild_heap(mesh, heap, opt.mode, opt.lock_boundary, &stats);
            emit_progress("rebuild", mesh, heap, stats, t0, opt.verbose, csv.is_open() ? &csv : nullptr);
        } else if (opt.progress_interval > 0 && stats.accepted % opt.progress_interval == 0) {
            emit_progress("progress", mesh, heap, stats, t0, opt.verbose, csv.is_open() ? &csv : nullptr);
        }
    }

    stats.final_vertices = mesh.alive_vertex_count();
    stats.final_faces = mesh.alive_face_count();
    stats.total_seconds = seconds_since(t0);
    emit_progress(heap.empty() ? "final_heap_empty" : "final", mesh, heap, stats, t0, opt.verbose, csv.is_open() ? &csv : nullptr);
    return stats;
}

} // namespace pndf
