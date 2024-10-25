#pragma once
#include <igl/readOFF.h>
#include <igl/slice.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/arap.h>
#include <igl/harmonic.h>
#include <igl/readOBJ.h>
#include <iostream>
#include <igl/unproject_ray.h>
#include <lsc/MeshProcessing.h>
#include <lsc/basic.h>
#include <lsc/tools.h>
#include <lsc/interaction.h>
// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <igl/parula.h>
#include <igl/isolines_map.h>
#include <igl/Timer.h>
#include <functional>

typedef OpenMesh::PolyMesh_ArrayKernelT<> CGMesh;

// the interface operators
class lscif
{
private:
#ifdef _WIN32
	std::string spliter = "\\";
#else
	std::string spliter = "/";
#endif
public:
	lscif();
	igl::Timer timer_global;
	double time_total = 0;
	int iteration_total = 0;
	std::vector<double> errorlist;
	std::vector<double> timelist; 

	Eigen::RowVector3d sea_green;
	Eigen::RowVector3d hot_red;
	Eigen::RowVector3d pure_blk;
	MeshProcessing MP;

	// This V and F are only for visulization of the default mesh
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;

	bool keyPress_1 = false;
	bool keyPress_2 = false;
	bool keyPress_d = false;

	// Optimization Parameters
	double weigth_closeness = 0.0;
	double refer_AveEL = 1;
	double weight_mass = 100.;	  // weight of the mass function to
	double weight_boundary = 100; // weight of the boundary condition
	double weight_laplacian = 0.001;
	double weight_pseudo_geodesic = 0.01;
	double weight_strip_width = 0.0001;
	double weight_geodesic = 3;
	// double weight_shading = 0.01; // For shading: the light is othogonal to the tangent direction

	double maximal_step_length = 0.5;
	bool enable_extreme_cases = false;
	bool Given_Const_Direction = false;
	int which_levelset = 0;
	bool fix_angle_of_two_levelsets = false;
	double angle_between_two_levelsets = 60;
	double weight_fix_two_ls_angle = 0;

	// Tracing Parameters
	int which_seg_id = 0; // the face id of which we will assign value to
	int nbr_of_intervals = 3;
	int start_bnd_he = 0; // start boundary halfedge, for init the values by assign values on some edges
	int nbr_edges = 20;
	int id_debug_global = 5;
	double target_angle = 60;
	double target_angle_2 = 60;
	double start_angle = 60;
	double threadshold_angel_degree = 150; // the threadshold for detecting boundary corners
	int dbg_int = 0;
	int dbg_int2 = 5;
	double dbg_dbl = 30;
	double dbg_dbl2 = 60;
	double vector_scaling = 1;
	double vector_scaling2 = 0;
	int extracted_nbr = 5;
	bool debug_flag = false;
	std::vector<int> fixedVertices;
	double ptInx = 0;
	double ptIny = 0;
	double ptInz = 0;
	int TracingType = 0;
	// int QuadOptType = 0;

	int OpIter = 10;
	CGMesh mesh;
	CGMesh RefMesh;
	bool enable_pg_energy_checkbox = false;
	bool enable_strip_width_checkbox = false;
	bool fixBoundary_checkbox = false;

	// mesh processing for level set
	int nbr_lines_first_ls = 30;
	int nbr_lines_second_ls = 30;
	Eigen::VectorXd readed_LS1;
	Eigen::VectorXd readed_LS2;
	Eigen::VectorXd readed_LS3;
	std::vector<CGMesh> readed_mesh1;
	CGMesh readed_mesh2;
	int Nbr_Iterations_Mesh_Opt = 10;
	double weight_Mesh_mass = 100;
	double weight_Mesh_smoothness = 0.001;
	double weight_Mesh_edgelength = 0.01;
	double weight_Mesh_pesudo_geodesic = 30;
	double Mesh_opt_max_step_length = 0.1;
	double weight_Mesh_approximation = 0.01;

	int update_ver_id = 0;
	int ring_nbr;
	std::vector<int> update_verlist;

	bool selectFEV[30] = {true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true};
	std::vector<std::string> meshFileName;
	std::vector<CGMesh> Meshes;
	lsTools tools;
	PolyOpt poly_tool;
	QuadOpt quad_tool;
	double weight_angle = 0;
	double InputPx = 0;	  // theta
	double InputPy = 180; // phi

	double InputPx1 = 0;
	double InputPy1 = 180;

	double InputThetaTol = 23.5;  // give 10 degrees of tolerance
	double InputThetaTol1 = 23.5; // give 10 degrees of tolerance
	double InputPhiTol = 45;	  // give 10 degrees of tolerance
	double InputPhiTol1 = 45;	  // give 10 degrees of tolerance
	int vinrow = -1;

	double Shading_Latitude = 34;
	bool shading_init = false;
	bool let_ray_through = false;
	bool let_ray_reflect = false;
	bool recompute_auxiliaries = false;

	double weight_binormal = 1;
	double weight_smt_binormal = 0;
	double weight_endpoint_ratio = 1;
	bool pick_single_ply = false;
	int pick_line_id = 0;
	std::vector<std::vector<Eigen::Vector3d>> Poly_readed;
	std::vector<std::vector<Eigen::Vector3d>> Bino_readed;
	std::vector<std::vector<Eigen::Vector3d>> Poly_opt;
	std::vector<std::vector<Eigen::Vector3d>> Bino_opt;

	std::vector<int> VertexType;

	// interaction stuff
	bool draw_strokes = false;
	bool left_button_down = false;
	std::vector<std::vector<double>> project_2dx;
	std::vector<std::vector<double>> project_2dy;
	std::vector<double> project_x_tmp;
	std::vector<double> project_y_tmp;
	std::vector<std::vector<int>> flist;
	std::vector<std::vector<Eigen::Vector3f>> bclist;
	std::vector<int> flist_tmp;
	std::vector<Eigen::Vector3f> bclist_tmp;
	std::vector<std::vector<Eigen::Vector3d>> project_pts;
	AutoRunArgs autorunner[2];

	// functional angles stuff.
	bool Disable_Changed_Angles = false;
	bool Use_Fitting_Angles = false;
	bool Use_Opt_Only_BNMS = false;

	bool AxisFixedForSlopes = false;
	bool AnglesFixedForSlopes = false;
	double AxisFixInX = 0;
	double AxisFixInY = 0;
	double AxisFixInZ = 1;
	double AngleFixIn = 45;

	bool CentroidFixedForNodes = false;
	bool RadiusFixedForNodes = false;
	bool FindCandRForNodes = false; // take C and R as variables
	double CentroidFixInX = 0;
	double CentroidFixInY = 0;
	double CentroidFixInZ = 1;
	double RadiusFixIn = 45;


	// add a new mesh into the mesh lists, and show the new mesh along with previous showed meshes
	void updateMeshViewer(igl::opengl::glfw::Viewer &viewer, CGMesh &mesh);

	void plot2dFittingResults(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd& xs, const Eigen::VectorXd& ys,
		const std::vector<double> &ply);

	bool pre_draw(igl::opengl::glfw::Viewer &viewer);
	bool key_up(igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods);
	bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods);

	// when holding "1" or "2", with mouse left button down, points can be choosen
	// TODO fix the bug when select one mesh but click on the other one, and when swich to
	// the other mesh, all the histories are recorded and reproduce points on the second mesh
	bool mouse_down(igl::opengl::glfw::Viewer &viewer, int button, int modifier);
	bool mouse_up(igl::opengl::glfw::Viewer &viewer, int button, int modifier);
	bool mouse_move(igl::opengl::glfw::Viewer &viewer, int mouse_x, int mouse_y);

	void show_axis(Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, double ratio);

	void viewer_launch_functions(igl::opengl::glfw::Viewer& viewer, igl::opengl::glfw::imgui::ImGuiPlugin &plugin, 
	igl::opengl::glfw::imgui::ImGuiMenu &menu, igl::opengl::glfw::imgui::ImGuiMenu &menu_mp);
	void draw_menu1(igl::opengl::glfw::Viewer& viewer, igl::opengl::glfw::imgui::ImGuiPlugin &plugin, 
	igl::opengl::glfw::imgui::ImGuiMenu &menu);
	void draw_menu2(igl::opengl::glfw::Viewer& viewer, igl::opengl::glfw::imgui::ImGuiPlugin &plugin, 
	igl::opengl::glfw::imgui::ImGuiMenu &menu);
};

void launch_viewer();