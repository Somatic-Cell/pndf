#include "pndf/math.hpp"
#include "pndf/normal_map.hpp"
#include "pndf/torus_mesh.hpp"
#include <cmath>
#include "test_common.hpp"


using namespace pndf;

TEST_CASE(test_wrap_delta_midpoint) {
    REQUIRE_NEAR(wrap_delta(0.9), -0.1, 1e-12);
    REQUIRE_NEAR(wrap_delta(-0.9), 0.1, 1e-12);
    Vec2 m = periodic_midpoint({0.99, 0.4}, {0.01, 0.4});
    REQUIRE_NEAR(m.x, 0.0, 1e-12);
    REQUIRE_NEAR(m.y, 0.4, 1e-12);
}

TEST_CASE(test_torus_mesh_counts) {
    NormalMap map;
    map.width = 8; map.height = 8; map.nxy.resize(64, {0.0, 0.0});
    TorusMesh mesh; mesh.build_full_torus(map);
    REQUIRE(mesh.vertices.size() == 64);
    REQUIRE(mesh.faces.size() == 128);
    REQUIRE(mesh.vertex_id(8, 0) == mesh.vertex_id(0,0));
    REQUIRE(mesh.vertex_id(-1, 0) == mesh.vertex_id(7,0));
}