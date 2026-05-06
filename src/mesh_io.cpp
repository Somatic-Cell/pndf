#include "pndf/mesh_io.hpp"
#include <fstream>
#include <stdexcept>
#include <unordered_map>
#include <cstdint>

namespace pndf {

void write_mesh_binary(const TorusMesh& mesh, const std::string& path) {
    std::ofstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("Cannot write mesh: " + path);
    std::vector<int32_t> remap(mesh.vertices.size(), -1);
    int32_t vcount = 0;
    for (size_t i=0;i<mesh.vertices.size();++i) if (mesh.vertices[i].alive) remap[i] = vcount++;
    std::vector<std::array<int32_t,3>> faces;
    for (const Face& face : mesh.faces) {
        if (!face.alive) continue;
        int32_t a=remap[face.v[0]], b=remap[face.v[1]], c=remap[face.v[2]];
        if (a>=0 && b>=0 && c>=0 && a!=b && b!=c && c!=a) faces.push_back({a,b,c});
    }
    int32_t width = mesh.width, height = mesh.height, fcount = static_cast<int32_t>(faces.size());
    f.write(reinterpret_cast<const char*>(&width), sizeof(int32_t));
    f.write(reinterpret_cast<const char*>(&height), sizeof(int32_t));
    f.write(reinterpret_cast<const char*>(&vcount), sizeof(int32_t));
    f.write(reinterpret_cast<const char*>(&fcount), sizeof(int32_t));
    for (const Vertex& v : mesh.vertices) if (v.alive) {
        float rec[4] = {static_cast<float>(v.uv.x), static_cast<float>(v.uv.y), static_cast<float>(v.nxy.x), static_cast<float>(v.nxy.y)};
        f.write(reinterpret_cast<const char*>(rec), sizeof(rec));
    }
    for (auto tri : faces) f.write(reinterpret_cast<const char*>(tri.data()), 3*sizeof(int32_t));
}

TorusMesh read_mesh_binary(const std::string& path) {
    std::ifstream f(path, std::ios::binary);
    if (!f) throw std::runtime_error("Cannot read mesh: " + path);
    int32_t width=0,height=0,vcount=0,fcount=0;
    f.read(reinterpret_cast<char*>(&width), sizeof(int32_t));
    f.read(reinterpret_cast<char*>(&height), sizeof(int32_t));
    f.read(reinterpret_cast<char*>(&vcount), sizeof(int32_t));
    f.read(reinterpret_cast<char*>(&fcount), sizeof(int32_t));
    TorusMesh mesh;
    mesh.width = width; mesh.height = height;
    mesh.vertices.resize(vcount);
    for (int32_t i=0;i<vcount;++i) {
        float rec[4]; f.read(reinterpret_cast<char*>(rec), sizeof(rec));
        mesh.vertices[i].uv = {rec[0], rec[1]};
        mesh.vertices[i].nxy = clamp_projected_normal({rec[2], rec[3]});
        mesh.vertices[i].alive = true;
    }
    mesh.faces.resize(fcount);
    for (int32_t i=0;i<fcount;++i) {
        int32_t tri[3]; f.read(reinterpret_cast<char*>(tri), sizeof(tri));
        mesh.faces[i].v = {tri[0], tri[1], tri[2]};
        mesh.faces[i].alive = true;
    }
    mesh.rebuild_vertex_face_adjacency();
    return mesh;
}

void write_obj_debug(const TorusMesh& mesh, const std::string& path) {
    std::ofstream f(path);
    if (!f) throw std::runtime_error("Cannot write obj: " + path);
    std::vector<int> remap(mesh.vertices.size(), -1);
    int id = 1;
    for (size_t i=0;i<mesh.vertices.size();++i) if (mesh.vertices[i].alive) {
        remap[i] = id++;
        const auto& v = mesh.vertices[i];
        f << "v " << v.uv.x << " " << v.uv.y << " " << 0.0 << "\n";
        f << "vt " << v.uv.x << " " << v.uv.y << "\n";
    }
    for (const Face& face : mesh.faces) if (face.alive) {
        int a=remap[face.v[0]], b=remap[face.v[1]], c=remap[face.v[2]];
        if (a>0 && b>0 && c>0) f << "f " << a << "/" << a << " " << b << "/" << b << " " << c << "/" << c << "\n";
    }
}

} // namespace pndf
