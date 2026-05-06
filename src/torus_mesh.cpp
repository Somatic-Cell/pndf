#include "pndf/torus_mesh.hpp"
#include <algorithm>

namespace pndf {

int TorusMesh::vertex_id(int x, int y) const {
    x %= width; if (x < 0) x += width;
    y %= height; if (y < 0) y += height;
    return y * width + x;
}

void TorusMesh::build_full_torus(const NormalMap& normal_map) {
    width = normal_map.width;
    height = normal_map.height;
    vertices.clear();
    faces.clear();
    vertices.resize(static_cast<size_t>(width) * height);
    for (int y=0;y<height;++y) {
        for (int x=0;x<width;++x) {
            const int id = vertex_id(x,y);
            vertices[id].uv = { (x + 0.5) / double(width), (y + 0.5) / double(height) };
            vertices[id].nxy = normal_map.nxy[static_cast<size_t>(y)*width + x];
            vertices[id].alive = true;
            vertices[id].version = 0;
            vertices[id].q = {};
        }
    }
    faces.reserve(static_cast<size_t>(2) * width * height);
    for (int y=0;y<height;++y) {
        for (int x=0;x<width;++x) {
            const int v00 = vertex_id(x,y);
            const int v10 = vertex_id(x+1,y);
            const int v01 = vertex_id(x,y+1);
            const int v11 = vertex_id(x+1,y+1);
            faces.push_back({{v00, v10, v01}, true});
            faces.push_back({{v01, v10, v11}, true});
        }
    }
    rebuild_vertex_face_adjacency();
}

void TorusMesh::rebuild_vertex_face_adjacency() {
    for (auto& v : vertices) v.faces.clear();
    for (int fi=0; fi<static_cast<int>(faces.size()); ++fi) {
        if (!faces[fi].alive) continue;
        for (int k=0;k<3;++k) {
            const int vi = faces[fi].v[k];
            if (vertices[vi].alive) vertices[vi].faces.push_back(fi);
        }
    }
}

int TorusMesh::alive_vertex_count() const {
    int n=0;
    for (const auto& v : vertices) if (v.alive) ++n;
    return n;
}

int TorusMesh::alive_face_count() const {
    int n=0;
    for (const auto& f : faces) if (f.alive) ++n;
    return n;
}

std::array<Vec2,3> unwrapped_face_uvs(const TorusMesh& mesh, const Face& face) {
    std::array<Vec2,3> uv{};
    uv[0] = mesh.vertices[face.v[0]].uv;
    for (int i=1;i<3;++i) {
        Vec2 raw = mesh.vertices[face.v[i]].uv;
        Vec2 d = wrap_delta(raw - uv[0]);
        uv[i] = uv[0] + d;
    }
    return uv;
}

double oriented_area_unwrapped(const std::array<Vec2,3>& uv) {
    const Vec2 a = uv[1] - uv[0];
    const Vec2 b = uv[2] - uv[0];
    return a.x*b.y - a.y*b.x;
}

} // namespace pndf
