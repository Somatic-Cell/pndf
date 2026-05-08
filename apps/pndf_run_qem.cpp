#include "pndf/edge_collapse.hpp"
#include "pndf/mesh_io.hpp"
#include "pndf/normal_map.hpp"
#include "pndf/torus_mesh.hpp"
#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <functional>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using namespace pndf;

namespace {
static void usage() {
    std::cerr
        << "Usage:\n"
        << "  Single target:\n"
        << "    pndf_run_qem --input nxy.bin --output mesh.bin --target N --mode normal|4d [options]\n"
        << "  Multi target checkpoint sweep:\n"
        << "    pndf_run_qem --input nxy.bin --output-prefix prefix --targets N1,N2,... --mode normal|4d [options]\n"
        << "Options:\n"
        << "  --quiet                         Disable progress messages on stderr.\n"
        << "  --rebuild-interval N            Full heap rebuild interval in accepted collapses.\n"
        << "                                  Use 0 to disable periodic full rebuilds. Default: 50000.\n"
        << "  --progress-interval N           Progress print/CSV interval in accepted collapses.\n"
        << "                                  Use 0 to disable interval progress. Default: 50000.\n"
        << "  --profile-csv path.csv          Write progress/profiling samples to CSV.\n";
}
static std::vector<int> parse_targets(const std::string& text) {
    std::vector<int> out;
    std::stringstream ss(text);
    std::string item;
    while (std::getline(ss, item, ',')) {
        if (item.empty()) continue;
        out.push_back(std::stoi(item));
    }
    out.erase(std::remove_if(out.begin(), out.end(), [](int v) { return v <= 0; }), out.end());
    std::sort(out.begin(), out.end(), std::greater<int>());
    out.erase(std::unique(out.begin(), out.end()), out.end());
    if (out.empty()) throw std::runtime_error("--targets must contain at least one positive integer");
    return out;
}
static std::filesystem::path checkpoint_path(const std::string& prefix, int target) {
    return std::filesystem::path(prefix + "_v" + std::to_string(target) + ".meshbin");
}
static void print_stats(const CollapseStats& stats) {
    std::cout << "accepted," << stats.accepted << "\n";
    std::cout << "rejected," << stats.rejected << "\n";
    std::cout << "final_vertices," << stats.final_vertices << "\n";
    std::cout << "final_faces," << stats.final_faces << "\n";
    std::cout << "seconds," << stats.total_seconds << "\n";
    std::cout << "heap_pops," << stats.heap_pops << "\n";
    std::cout << "stale_pops," << stats.stale_pops << "\n";
    std::cout << "dead_pops," << stats.dead_pops << "\n";
    std::cout << "invalid_cost_pops," << stats.invalid_cost_pops << "\n";
    std::cout << "fresh_recomputes," << stats.fresh_recomputes << "\n";
    std::cout << "fresh_requeues," << stats.fresh_requeues << "\n";
    std::cout << "fresh_rejects," << stats.fresh_rejects << "\n";
    std::cout << "rebuild_count," << stats.rebuild_count << "\n";
    std::cout << "max_heap_size," << stats.max_heap_size << "\n";
    std::cout << "initialize_seconds," << stats.initialize_seconds << "\n";
    std::cout << "rebuild_seconds," << stats.rebuild_seconds << "\n";
    std::cout << "make_item_seconds," << stats.make_item_seconds << "\n";
    std::cout << "link_condition_seconds," << stats.link_condition_seconds << "\n";
    std::cout << "placement_seconds," << stats.placement_seconds << "\n";
    std::cout << "non_flip_seconds," << stats.non_flip_seconds << "\n";
    std::cout << "collapse_update_seconds," << stats.collapse_update_seconds << "\n";
}
} // namespace

int main(int argc, char** argv) {
    std::string input, output, output_prefix, mode_str = "normal";
    std::string targets_text;
    int target = 0;
    bool verbose = true;
    int rebuild_interval = 50000;
    int progress_interval = 50000;
    std::string profile_csv;

    for (int i=1;i<argc;++i) {
        std::string a = argv[i];
        auto next = [&]() -> std::string { if (i+1>=argc) { usage(); std::exit(2); } return argv[++i]; };
        if (a == "--input") input = next();
        else if (a == "--output") output = next();
        else if (a == "--output-prefix") output_prefix = next();
        else if (a == "--target") target = std::stoi(next());
        else if (a == "--targets") targets_text = next();
        else if (a == "--mode") mode_str = next();
        else if (a == "--quiet") verbose = false;
        else if (a == "--rebuild-interval") rebuild_interval = std::stoi(next());
        else if (a == "--progress-interval") progress_interval = std::stoi(next());
        else if (a == "--profile-csv" || a == "--progress-csv") profile_csv = next();
        else { usage(); return 2; }
    }

    const bool multi_target = !targets_text.empty();
    if (input.empty()) { usage(); return 2; }
    if (mode_str != "normal" && mode_str != "4d") { usage(); return 2; }
    if (multi_target) {
        if (!output.empty() || target > 0 || output_prefix.empty()) { usage(); return 2; }
    } else {
        if (output.empty() || target <= 0 || !output_prefix.empty()) { usage(); return 2; }
    }

    try {
        std::vector<int> targets;
        if (multi_target) {
            targets = parse_targets(targets_text);
            target = targets.back();
        }

        auto map = read_nxy_binary(input);
        TorusMesh mesh;
        mesh.build_full_torus(map);

        CollapseOptions opt;
        opt.target_vertices = target;
        opt.verbose = verbose;
        opt.mode = (mode_str == "4d") ? QemMode::Qem4D : QemMode::NormalOnly;
        opt.rebuild_interval = rebuild_interval;
        opt.progress_interval = progress_interval;
        opt.progress_csv = profile_csv;

        std::vector<int> saved_checkpoints;
        if (multi_target) {
            opt.checkpoint_vertices = targets;
            opt.checkpoint_callback = [&](int checkpoint_target, const TorusMesh& checkpoint_mesh, const CollapseStats&) {
                const auto path = checkpoint_path(output_prefix, checkpoint_target);
                if (path.parent_path() != std::filesystem::path{}) std::filesystem::create_directories(path.parent_path());
                write_mesh_binary(checkpoint_mesh, path.string());
                saved_checkpoints.push_back(checkpoint_target);
                if (verbose) std::cerr << "wrote checkpoint target=" << checkpoint_target << " path=" << path.string() << "\n";
            };
        }

        auto stats = simplify_qem(mesh, opt);

        if (multi_target) {
            if (saved_checkpoints.empty() || saved_checkpoints.back() != stats.final_vertices) {
                const auto path = std::filesystem::path(output_prefix + "_final_v" + std::to_string(stats.final_vertices) + ".meshbin");
                if (path.parent_path() != std::filesystem::path{}) std::filesystem::create_directories(path.parent_path());
                write_mesh_binary(mesh, path.string());
                if (verbose) std::cerr << "wrote final mesh path=" << path.string() << "\n";
            }
            std::cout << "checkpoints_written," << saved_checkpoints.size() << "\n";
            for (int t : saved_checkpoints) std::cout << "checkpoint," << t << "," << checkpoint_path(output_prefix, t).string() << "\n";
        } else {
            write_mesh_binary(mesh, output);
        }

        print_stats(stats);
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
