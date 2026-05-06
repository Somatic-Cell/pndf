#include "pndf/eval.hpp"
#include "pndf/mesh_io.hpp"
#include "pndf/normal_map.hpp"
#include <iostream>
#include <string>

using namespace pndf;

static void usage() {
    std::cerr << "Usage: pndf_eval_mesh --input nxy.bin --mesh mesh.bin --csv stats.csv\n";
}

int main(int argc, char** argv) {
    std::string input, mesh_path, csv;
    for (int i=1;i<argc;++i) {
        std::string a = argv[i];
        auto next = [&]() -> std::string { if (i+1>=argc) { usage(); std::exit(2); } return argv[++i]; };
        if (a == "--input") input = next();
        else if (a == "--mesh") mesh_path = next();
        else if (a == "--csv") csv = next();
        else { usage(); return 2; }
    }
    if (input.empty() || mesh_path.empty() || csv.empty()) { usage(); return 2; }
    try {
        NormalMap gt = read_nxy_binary(input);
        TorusMesh mesh = read_mesh_binary(mesh_path);
        auto stats = evaluate_vertex_sample_error(gt, mesh);
        write_stats_csv(stats, csv);
        std::cout << "mean_deg," << stats.mean_deg << "\n";
        std::cout << "p95_deg," << stats.p95_deg << "\n";
        std::cout << "p99_deg," << stats.p99_deg << "\n";
        std::cout << "max_deg," << stats.max_deg << "\n";
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
