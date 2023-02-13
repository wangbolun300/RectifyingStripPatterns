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

// the default arguments for running the 
class AutoRunArgs{
    public:
    AutoRunArgs(){};
    bool compute_pg;
    double stop_step_length; // when the step length lower than this, stop

    std::vector<int> iterations; // nbr of iterations
    std::vector<double> weight_gravity;
    std::vector<double> weight_lap ;
    std::vector<double> weight_bnd;
    std::vector<double> weight_pg;
    std::vector<double> weight_strip_width;
};
void AssignAutoRunDefaultArgs(AutoRunArgs &args, const bool compute_pg);