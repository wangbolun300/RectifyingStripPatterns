
#include <igl/read_triangle_mesh.h>
#include <igl/triangulated_grid.h>
#include <igl/heat_geodesics.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/avg_edge_length.h>
#include <igl/isolines_map.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/create_shader_program.h>
#include <igl/opengl/destroy_shader_program.h>
#include <iostream>
#include <lsc/tools.h>
int main(int argc, char *argv[])
{
  runRotateArrows(0.15);
}