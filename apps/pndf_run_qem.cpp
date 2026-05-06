#include "pndf/edge_collapse.hpp"
#include "pndf/mesh_io.hpp"
#include "pndf/normal_map.hpp"
#include "pndf/torus_mesh.hpp"
#include <iostream>
#include <string>
#include <chrono>

using namespace pndf;

static void usage() {
    std::cerr << "Usage: pndf_run_qem --input nxy.bin --output mesh.bin --target N --mode normal|4d [--quiet]\n";
}

int main(int argc, char** argv) {
    std::string input, output, mode_str = "normal";
    int target = 0;
    bool verbose = true;
    for (int i=1;i<argc;++i) {
        std::string a = argv[i];
        auto next = [&]() -> std::string { if (i+1>=argc) { usage(); std::exit(2); } return argv[++i]; };
        if (a == "--input") input = next();
        else if (a == "--output") output = next();
        else if (a == "--target") target = std::stoi(next());
        else if (a == "--mode") mode_str = next();
        else if (a == "--quiet") verbose = false;
        else { usage(); return 2; }
    }
    if (input.empty() || output.empty() || target <= 0) { usage(); return 2; }
    try {
        auto map = read_nxy_binary(input);
        TorusMesh mesh;
        mesh.build_full_torus(map);
        CollapseOptions opt;
        opt.target_vertices = target;
        opt.verbose = verbose;
        opt.mode = (mode_str == "4d") ? QemMode::Qem4D : QemMode::NormalOnly;
        auto t0 = std::chrono::steady_clock::now();
        auto stats = simplify_qem(mesh, opt);
        auto t1 = std::chrono::steady_clock::now();
        write_mesh_binary(mesh, output);
        double sec = std::chrono::duration<double>(t1-t0).count();
        std::cout << "accepted," << stats.accepted << "\n";
        std::cout << "rejected," << stats.rejected << "\n";
        std::cout << "final_vertices," << stats.final_vertices << "\n";
        std::cout << "final_faces," << stats.final_faces << "\n";
        std::cout << "seconds," << sec << "\n";
    } catch (const std::exception& e) {
        std::cerr << "error: " << e.what() << "\n";
        return 1;
    }
}
