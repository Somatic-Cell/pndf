#pragma once
#include <array>
#include <vector>
#include "pndf/math.hpp"
#include "pndf/normal_map.hpp"

namespace pndf {

struct Vertex {
    Vec2 uv;     // torus coordinate in [0,1)
    Vec2 nxy;    // projected normal
    bool alive = true;
    bool boundary_locked = false; // debug mode only; true torus does not need this
    int version = 0;
    Mat5 q;
    std::vector<int> faces;
};

struct Face {
    std::array<int,3> v{};
    bool alive = true;
};

struct TorusMesh {
    int width = 0;
    int height = 0;
    std::vector<Vertex> vertices;
    std::vector<Face> faces;

    int vertex_id(int x, int y) const;
    void build_full_torus(const NormalMap& normal_map);
    void rebuild_vertex_face_adjacency();
    int alive_vertex_count() const;
    int alive_face_count() const;
};

std::array<Vec2,3> unwrapped_face_uvs(const TorusMesh& mesh, const Face& face);

double oriented_area_unwrapped(const std::array<Vec2,3>& uv);

} // namespace pndf
