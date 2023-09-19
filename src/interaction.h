#pragma once
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
                         const std::vector<double> &ys, Eigen::MatrixXd& Vers, std::vector<int> &fout, 
                         std::vector<Eigen::Vector3f> &bcout);
class PartOfAutoRun{
    public:
    PartOfAutoRun(){};
    int iterations; // nbr of iterations
    double weight_gravity;
    double weight_lap ;
    double weight_bnd;
    double weight_pg;
    double weight_strip_width;
};


// the default arguments for running the pseudo-geodesics with drawing strokes
class AutoRunArgs{
    public:
    AutoRunArgs(){};
    std::vector<PartOfAutoRun> parts;
    bool compute_pg;
    double stop_step_length; // when the step length lower than this, stop
    double stop_energy_sqrt; // when the energy sqrt lower than this, stop

    
};
void AssignAutoRunDefaultArgs(AutoRunArgs &args, const bool compute_pg);