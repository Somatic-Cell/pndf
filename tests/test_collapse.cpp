#include "pndf/edge_collapse.hpp"
#include "pndf/normal_map.hpp"
#include "pndf/torus_mesh.hpp"
#include "test_common.hpp"


using namespace pndf;

TEST_CASE(test_small_torus_collapse_normal_qem) {
    NormalMap map;
    map.width = 8; map.height = 8; map.nxy.resize(64, {0.0, 0.0});
    TorusMesh mesh; mesh.build_full_torus(map);
    CollapseOptions opt;
    opt.target_vertices = 32;
    opt.mode = QemMode::NormalOnly;
    opt.verbose = false;
    opt.rebuild_interval = 16;
    auto stats = simplify_qem(mesh, opt);
    REQUIRE(stats.final_vertices <= 64);
    REQUIRE(stats.final_faces > 0);
}