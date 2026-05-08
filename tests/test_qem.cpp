#include "pndf/qem.hpp"
#include "pndf/normal_map.hpp"
#include "pndf/torus_mesh.hpp"
#include <cmath>
#include "test_common.hpp"


using namespace pndf;

TEST_CASE(test_normal_qem_affine_field_low_cost) {
    NormalMap map;
    map.width = 4; map.height = 4; map.nxy.resize(16);
    for (int y=0;y<4;y++) for (int x=0;x<4;x++) {
        double u = (x+0.5)/4.0;
        double v = (y+0.5)/4.0;
        map.nxy[y*4+x] = {0.1*u + 0.05*v, -0.03*u + 0.07*v};
    }
    TorusMesh mesh; mesh.build_full_torus(map);
    const Face& f = mesh.faces[0];
    Mat5 q = triangle_quadric_normal_only(mesh, f);
    const Vertex& a = mesh.vertices[f.v[0]];
    double c = eval_quadric(q, {a.uv.x, a.uv.y, a.nxy.x, a.nxy.y, 1.0});
    REQUIRE(c < 1e-8);
}

TEST_CASE(test_uv_quadric_rebase_translation) {
    Mat5 q;
    // Plane: u_old - nx = 0.
    q.add_outer({1.0, 0.0, -1.0, 0.0, 0.0});
    Mat5 qr = translate_uv_quadric(q, {1.0, 0.0});
    // In the new chart, u_old = u_new + 1, so nx must equal u_new + 1.
    double c = eval_quadric(qr, {0.25, 0.0, 1.25, 0.0, 1.0});
    REQUIRE(c < 1e-10);
}

TEST_CASE(test_best_placement_does_not_wrap_seam_midpoint) {
    Mat5 q;
    // A zero quadric makes all candidates equal; this tests the candidate
    // construction itself.  For a seam-crossing edge 0.99 -> 0.01, the
    // unwrapped midpoint should be 1.0, not 0.0.
    Placement p = best_placement(q, {0.99, 0.5}, {0.0, 0.0}, {0.01, 0.5}, {0.0, 0.0}, QemMode::Qem4D);
    REQUIRE_NEAR(p.uv.x, 1.0, 1e-12);
    REQUIRE_NEAR(p.uv.y, 0.5, 1e-12);
}
