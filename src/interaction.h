#include<igl/unproject_onto_mesh.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include<lsc/basic.h>
// to use sleep()
#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#endif

void draw_stroke_on_mesh(const CGMesh &mesh, igl::opengl::glfw::Viewer &viewer,
                         const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const std::vector<double> &xs,
                         const std::vector<double> &ys, Eigen::MatrixXd& Vers);