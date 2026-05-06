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