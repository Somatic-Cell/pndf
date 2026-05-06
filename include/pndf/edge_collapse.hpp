#pragma once
#include <cstdint>
#include <queue>
#include <string>
#include "pndf/torus_mesh.hpp"
#include "pndf/qem.hpp"

namespace pndf {

struct CollapseOptions {
    int target_vertices = 0;
    QemMode mode = QemMode::NormalOnly;
    bool lock_boundary = false;       // false = true torus; true = conservative rectangular debug mode

    // 0 disables periodic full heap rebuilds.  The lazy priority queue still
    // remains valid because stale items are rejected by vertex versions.
    int rebuild_interval = 50000;

    // 0 disables progress messages.  This is independent of rebuild_interval,
    // so profiling can still be emitted when rebuilds are disabled.
    int progress_interval = 50000;

    // Optional CSV file for per-progress profiling samples.  Empty means disabled.
    std::string progress_csv;

    bool verbose = true;
};

struct CollapseStats {
    int accepted = 0;
    int rejected = 0;
    int final_vertices = 0;
    int final_faces = 0;

    std::uint64_t heap_pops = 0;
    std::uint64_t stale_pops = 0;
    std::uint64_t dead_pops = 0;
    std::uint64_t invalid_cost_pops = 0;
    std::uint64_t fresh_recomputes = 0;
    std::uint64_t fresh_requeues = 0;
    std::uint64_t fresh_rejects = 0;

    std::uint64_t rebuild_count = 0;
    std::uint64_t make_item_calls = 0;
    std::uint64_t collapse_updates = 0;
    std::uint64_t max_heap_size = 0;

    double initialize_seconds = 0.0;
    double rebuild_seconds = 0.0;
    double make_item_seconds = 0.0;
    double link_condition_seconds = 0.0;
    double placement_seconds = 0.0;
    double non_flip_seconds = 0.0;
    double collapse_update_seconds = 0.0;
    double total_seconds = 0.0;
};

CollapseStats simplify_qem(TorusMesh& mesh, const CollapseOptions& opt);

} // namespace pndf
