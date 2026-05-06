#pragma once
#include <string>
#include "pndf/torus_mesh.hpp"

namespace pndf {

void write_mesh_binary(const TorusMesh& mesh, const std::string& path);
TorusMesh read_mesh_binary(const std::string& path);
void write_obj_debug(const TorusMesh& mesh, const std::string& path);

} // namespace pndf
