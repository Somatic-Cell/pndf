#include "pndf/edge_collapse.hpp"
#include "pndf/mesh_io.hpp"
#include "pndf/normal_map.hpp"
#include "pndf/torus_mesh.hpp"
#include <chrono>
#include <iostream>
#include <string>

using namespace pndf;

static void usage() {
    std::cerr
        << "Usage: pndf_run_qem --input nxy.bin --output mesh.bin --target N --mode normal|4d [options]\n"
        << "Options:\n"
        << "  --quiet                         Disable progress messages on stderr.\n"
        << "  --rebuild-interval N            Full heap rebuild interval in accepted collapses.\n"
        << "                                  Use 0 to disable periodic full rebuilds. Default: 50000.\n"
        << "  --progress-interval N           Progress print/CSV interval in accepted collapses.\n"
        << "                                  Use 0 to disable interval progress. Default: 50000.\n"
        << "  --profile-csv path.csv          Write progress/profiling samples to CSV.\n";
}

int main(int argc, char** argv) {
    std::string input, output, mode_str = "normal";
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
        else if (a == "--target") target = std::stoi(next());
        else if (a == "--mode") mode_str = next();
        else if (a == "--quiet") verbose = false;
        else if (a == "--rebuild-interval") rebuild_interval = std::stoi(next());
        else if (a == "--progress-interval") progress_interval = std::stoi(next());
        else if (a == "--profile-csv" || a == "--progress-csv") profile_csv = next();
        else { usage(); return 2; }
    }

    if (input.empty() || output.empty() || target <= 0) { usage(); return 2; }
    if (mode_str != "normal" && mode_str != "4d") { usage(); return 2; }

    try {
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

        auto stats = simplify_qem(mesh, opt);
        write_mesh_binary(mesh, output);

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
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
