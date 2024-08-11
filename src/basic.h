﻿#pragma once
#include <lsc/MeshProcessing.h>
#include <lsc/igl_tool.h>
#include <lsc/basic.h>
#include <igl/AABB.h>

// Efunc represent a elementary value, which is the linear combination of
// some function values on their corresponding vertices, of this vertex.
// i.e. Efunc[0]*f_0+Efunc[1]*f_1+...+Efunc[n]*f_n.
typedef Eigen::SparseVector<double> Efunc;
typedef Eigen::SparseVector<int> SpVeci;
// typedef Eigen::SparseMatrix<double> spMat;
typedef Eigen::Triplet<double> Trip;
#define SCALAR_ZERO 1e-8
#define MERGE_VERTEX_RATIO 0.01
#define LSC_PI 3.14159265358979
#define ANGLE_TOLERANCE 2.
#define QUADRANT_TOLERANCE 0.04 // over 2 degree tolerance

static bool produce_no_smooth_oscu = false;// this is a flag to create wrong results when set as true.
static bool produce_small_search_range = false; // this is to enable small searching range for tracing.

enum BCtype{ // boundary condition type
    Traced, //0 // the traced values as the guide directions
    Opp_fixed// // fixing the vales of two opposite boundary edges
};

class EnergyPrepare{
    public:
    EnergyPrepare(){};
    double weight_gravity;  
    double weight_lap; 
    double weight_bnd;
    double weight_pg; 
    double weight_strip_width;
    bool solve_pseudo_geodesic;
    bool solve_strip_width_on_traced;
    bool enable_extreme_cases;
    double target_angle;
    double max_step_length;
    bool Given_Const_Direction;
    // Eigen::Vector3d Reference_ray;
};
class MeshEnergyPrepare{
    public:
    MeshEnergyPrepare(){};
    double weight_Mesh_smoothness;
    double weight_Mesh_pesudo_geodesic;
    double weight_Mesh_edgelength;
    double weight_mass;
    double Mesh_opt_max_step_length;
    bool enable_extreme_cases;
    double target_angle;

    bool Given_Const_Direction;
    // Eigen::Vector3d Reference_ray;
};
class TracingPrepare{
    public:
    TracingPrepare(){};
    int every_n_edges;
    int which_boundary_segment;
    double target_angle;
    double start_angle;
    double threadshold_angel_degree;
    int start_bnd_he;// start boundary halfedge, for init the values by assign values on some edges
    int nbr_edges;

};
class LSAnalizer {
public:
    LSAnalizer() {};
    Eigen::VectorXi LocalActInner;
    Eigen::VectorXi LocalActBpair;
    std::vector<bool> Correspondance;
    Eigen::VectorXi Special; // record if we apply the second condition in shading
    // Eigen::VectorXi HighEnergy; // detect high energy vertices
    std::vector<CGMesh::HalfedgeHandle> heh0; 
    std::vector<CGMesh::HalfedgeHandle> heh1;
    std::vector<double> t1s;
    std::vector<double> t2s;

};
class NeighbourInfo
{
public:
    NeighbourInfo()
    {
        is_vertex = 0;
        round = 0;
    };
    bool is_vertex = false;
    Eigen::Vector3d pnorm;
    CGMesh::VertexHandle center_handle;
    std::vector<CGMesh::HalfedgeHandle> edges;
    int round = 0;
};
// for polyline optimization
class PolyOpt
{
public:
    PolyOpt(){};
    // sample_nbr is the nbr of resampling points on the longest polyline. -1 means do not sample.
    void init(const std::vector<std::vector<Eigen::Vector3d>> &ply, const std::vector<std::vector<Eigen::Vector3d>> &bi, const int sample_nbr);
    void init(const Eigen::MatrixXd& ply_in);
    int VerNbr;
    Eigen::VectorXd PlyVars; // vars for the polylines
    Eigen::VectorXd OriVars; // original vars
    Eigen::VectorXd VinPly; // mark which polyline this vertex belong to
    // Eigen::MatrixXd OriV; // original vertex list of the (sampled) polylines 
    Eigen::VectorXi Front;   // the vertex id of the front point
    Eigen::VectorXi Back;    // the vertex id of the back point
    Eigen::MatrixXd norm_v;
    std::vector<std::vector<Eigen::Vector3d>> ply_extracted;
    std::vector<std::vector<Eigen::Vector3d>> bin_extracted;
    std::vector<std::vector<Eigen::Vector3d>> ply_rotated;
    std::vector<std::vector<Eigen::Vector3d>> bin_rotated;
    CGMesh RecMesh;// rectifying mesh
    // spMat endpts_signs;
    double weight_smooth;
    // double weight_plysmooth;
    double weight_mass;
    
    double weight_binormal; // the weight for correct binormals 
    double binormal_ratio; // the ratio added to binormal smoothness
    double max_step = 1;
    double strip_scale = 1;
    double weight_angle; // the weight for binormals having a expected angle with the surface normals
    double ratio_endpts;
    double target_angle;
    bool first_compute = true;
    bool pick_single_line = false;
    int pick_line_id;
    // if opt for polyline, then don't opt for crease.
    bool opt_for_polyline = false;
    bool opt_for_crease = false;
    Eigen::MatrixXd Vref;
    Eigen::MatrixXi Fref;
    Eigen::MatrixXd Nref;
    // int MaxNbrPinC = 0;
    bool OrientEndPts = true;
    bool singleLineProcessing = false;// if we only process single line, no need to glid on a surface.

    void opt();
    // make the values not far from original data
    void assemble_gravity(spMat& H, Eigen::VectorXd& B, Eigen::VectorXd &energy);
    void assemble_gravity_single_ply(spMat& H, Eigen::VectorXd& B, Eigen::VectorXd &energy);
    void assemble_polyline_smooth(const bool crease,  spMat& H, Eigen::VectorXd& B, Eigen::VectorXd &energy);
    void assemble_binormal_condition(spMat& H, Eigen::VectorXd& B, Eigen::VectorXd &energy);
    void assemble_angle_condition(spMat& H, Eigen::VectorXd& B, Eigen::VectorXd &energy);
    void extract_rectifying_plane_mesh();
    void extract_polylines_and_binormals();
    void rotate_back_the_model_to_horizontal_coordinates(const double latitude);
    void orient_binormals_of_plyline(const std::vector<std::vector<Eigen::Vector3d>> &bi, std::vector<std::vector<Eigen::Vector3d>> &bout);
    void save_polyline_and_binormals_as_files(const bool rotated = false);
    void save_polyline_and_binormals_as_files(const std::vector<std::vector<Eigen::Vector3d>> &ply, const std::vector<std::vector<Eigen::Vector3d>> &bi);
    void save_polyline_and_binormals_as_files(const std::string &fname,
                                              const std::vector<std::vector<Eigen::Vector3d>> &ply, const std::vector<std::vector<Eigen::Vector3d>> &bi);
    void polyline_to_matrix(const std::vector<std::vector<Eigen::Vector3d>> &ply, Eigen::MatrixXd &V);
    void vertex_matrix_to_polyline(std::vector<std::vector<Eigen::Vector3d>> &ply, Eigen::MatrixXd& V);
    void get_normal_vector_from_reference(const Eigen::MatrixXd& Vr, const Eigen::MatrixXi& Fr, const Eigen::MatrixXd& nr);
    void force_smoothing_binormals();
    void force_invert_binormals();
    void related_error();// the related error between the strip midlines and the surface


    // the crease opt framework
    Eigen::MatrixXd vertices_cp;
    Eigen::MatrixXd tangents_cp;
    Eigen::MatrixXd binormal_cp;
    Eigen::VectorXd AngCollector;
    Eigen::VectorXd Inflecs;
    std::vector<bool> FlatPts;
    Eigen::VectorXd FlatDeterm;
    // Eigen::MatrixXd normalve_cp;
void init_crease_opt(const std::vector<std::vector<Eigen::Vector3d>> &vertices,
                     const std::vector<std::vector<Eigen::Vector3d>> &tangents,
                     const std::vector<std::vector<Eigen::Vector3d>> &binormals);
void assemble_gravity_crease(spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy);
void assemble_crease_planarity(spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy);
void extract_rectifying_plane_mesh_from_crease();
void extract_polylines_and_binormals_from_creases();

void opt_planarity();
void draw_inflections(Eigen::MatrixXd& pts);
};


class QuadOpt{
public:
    QuadOpt(){};
    void init(CGMesh& mesh_in, const std::vector<std::vector<double>> rowinfo, const std::vector<std::vector<double>>colinfo);
    void init(CGMesh& mesh_in, const std::string& prefix);
    void init(CGMesh& mesh_in);
	void initISO(const std::vector<Eigen::Vector3d> &Vlist, const int rnbr);
    void initAAG(const std::vector<Eigen::Vector3d> &Vlist, const int rnbr);
    void initAGG(const std::vector<Eigen::Vector3d> &Vlist, const int rnbr);
	void initGGG(const std::vector<Eigen::Vector3d> &Vlist, const int rnbr);
    int vNbrInRow = -1;
    std::vector<std::vector<double>> rowinfo; // to record the adjancency info of rows  
    std::vector<std::vector<double>> colinfo; // to record the adjancency info of cols
    Eigen::MatrixXi F;
    Eigen::MatrixXd V;
    CGMesh mesh_original;
    CGMesh mesh_update;
    std::vector<Eigen::Vector3d> verOriginal;
    std::vector<Eigen::Vector3d> verUpdate;
    std::vector<Eigen::Vector3d> curveRef;
    bool AGGAAGApproInitStrip = false; // approximate the first strip instead of the first curve
	bool GggEnableEditing = false;
	int GggEditingVid = -1;
	Eigen::Vector3d GggTargetPosition;
    double real_step_length;
	// evaluation of the error and the smoothness. For now it is only for the propagation method.
	std::vector<double> error_eval; 
	std::vector<double> smt_eval;
    
    // input from outside
    double angle_degree0;
    double angle_degree1;
    int OptType = 0;// 0: AAG. 1: GGA. 2: PP. 3: PPG. 4: AAGG
    int WhichDiagonal = 0; // 0 or 1
    double weight_fairness;
    double weight_pg;
    double weight_gravity;
    double weight_mass = 1;
    double weight_curve = 0;
    double pg_ratio = 1;// the ratio of diagonal (geodesic) energy to weight_pg
    double max_step = 1;
    Eigen::MatrixXd B0;
    Eigen::MatrixXd B1;
    Eigen::MatrixXd N;
    // int pg1_type = 0;//by default disabled. 1 means asymptotic, 2 means geodesic, 3 pseudo-geodesic
    // int pg2_type = 0;
    
    void opt();
    void optAAG();
    void optAGG();
	void optGGG();
	void optISO();
    Eigen::MatrixXd propagateBoundary();
    void reset();
    void load_triangle_mesh_tree(const igl::AABB<Eigen::MatrixXd, 3> &tree, const Eigen::MatrixXd &Vt,
                                 const Eigen::MatrixXi &Ft, const Eigen::MatrixXd &Nt);
    void show_curve_families(std::array<Eigen::MatrixXd, 3>& edges); 
    void show_diagonals(Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, Eigen::MatrixXd &E2, Eigen::MatrixXd &E3); 
    void extract_binormals(const int family, const int bnm_start, const int vid, Eigen::Vector3d& bi);
    void extract_binormals_propagation(const int family, const int vid, Eigen::Vector3d& bi, int varsOff);
    void extract_diagonals(const int family, std::vector<std::vector<int>> &digs);
    void write_polyline_info();
    void write_polyline_info_propagation();
    // void evaluateGGGcosineConstraints();
    Eigen::MatrixXd Debugtool;
private:
    MeshProcessing MP;
    int varsize;
    int d0_type = 0; // for the diagonal 0, 0 is disabled, 1 is asymptotic, 2 is geodesic, 3 is pg.
    int d1_type = 0; // for the diagonal 1, 0 is disabled, 1 is asymptotic, 2 is geodesic, 3 is pg.
    Eigen::VectorXd GlobVars; // variables for optimization.
    Eigen::VectorXd OrigVars; // copied original vars.
    Eigen::VectorXi row_front;
    Eigen::VectorXi row_back;
    Eigen::VectorXi col_front;
    Eigen::VectorXi col_back;
    Eigen::VectorXi d0_front; // diagonal 0 front
    Eigen::VectorXi d1_front; // diagonal 1 front
    Eigen::VectorXi d0_back;
    Eigen::VectorXi d1_back;
    bool ComputeAuxiliaries = true;
    bool Estimate_PG_Angles = true;
    spMat gravity_matrix;
    igl::AABB<Eigen::MatrixXd, 3> aabbtree;
    Eigen::MatrixXd Vtri;
    Eigen::MatrixXi Ftri;
    Eigen::MatrixXd Ntri;
	// some variables for isometric deformations
	Eigen::VectorXd d0Lengths;
	Eigen::VectorXd d1Lengths;
	Eigen::VectorXd d0d1dots;
    
    std::vector<Eigen::Vector3d> OriginalCurve;
    void get_Bnd(Eigen::VectorXi& Bnd); // the boundary vertices: the net corners
    void assemble_fairness(spMat& H, Eigen::VectorXd& B, Eigen::VectorXd &energy, const int order = 0);
    void assemble_gravity(spMat& H, Eigen::VectorXd& B, Eigen::VectorXd &energy);
    void assemble_gravity_AAG_AGG(spMat& H, Eigen::VectorXd& B, Eigen::VectorXd &energy, int order = 0);
    void assemble_gravity_ApproOriginal(spMat& H, Eigen::VectorXd& B, Eigen::VectorXd &energy);
    // type: 0 disabled. 1 means asymptotic, 2 means geodesic, 3 pseudo-geodesic
    // family: 0 rows, 1 cols, 2 diagonal NO0, 3 diagonal NO1.
    void assemble_pg_extreme_cases(spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy,
                                   const int type, const int family, const int aux_start_location, const int order = 0,
                                   const int whichBnm = 0);
    void assemble_pg_cases(const double angle_radian, spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy,
                           const int family, const int aux_start_location);
    void assemble_binormal_conditions(spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy, int family,
                                      int bnm_start, const int order = 0, const int whichBnm = 0);
    void assemble_normal_conditions(spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy, const int order = 0);
    void assemble_approximate_curve_conditions(spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy, const int order);
	void assemble_isometric_condition(spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy);

public:
    // debug tools
    Eigen::MatrixXd Enm0;
    Eigen::MatrixXd Enm1;
};

// The basic tool of LSC. Please initialize it with a mesh
class lsTools
{
private:
    MeshProcessing MP;
    Eigen::MatrixXd norm_f;               // normal per face
    Eigen::MatrixXd angF;                 // angel per vertex of each F. nx3, n is the number of faces
    Eigen::VectorXd areaF;                // the area of each face.
    Eigen::VectorXd areaPF;               // the area of each face in parametric domain;
    std::vector<Eigen::Matrix3d> Rotate;  // the 90 degree rotation matrices for each face.
    std::vector<Eigen::Matrix3d> RotateV; // the 90 degree rotation matrices on vertices
    // the rotated 3 half edges for each face, in the plane of this face
    // the order of the 3 directions are: v0-v2, v1-v0, v2-v1
    std::vector<std::array<Eigen::Vector3d, 3>> Erotate;
    // 2d rotated half edges in parametric domain
    std::vector<std::array<Eigen::Vector2d, 3>> Erotate2d;
    std::array<spMat, 3> gradVF;                  // gradient of function in each face
    std::array<spMat, 3> gradV;                   // gradient of function in each vertex
    std::array<std::array<spMat, 3>, 3> HessianV; // Hessian (2 order deriavate) of function in each vertex
    std::array<Eigen::MatrixXd, 2> Deriv1;        // the 1 order derivates for each vertex; 
    std::array<Eigen::MatrixXd, 4> Deriv2;        // the 2 order derivates for each vertex; 
    std::vector<double> II_L;                     // second fundamental form: L 
    std::vector<double> II_M;                     // second fundamental form: M 
    std::vector<double> II_N;                     // second fundamental form: N 
    // spMat Dlpsqr;                       // the derivates of ||laplacian F||^2 for all the vertices.
    spMat QcH; // curved harmonic energy matrix
    spMat Dlps;       // the derivates of laplacian F for all the vertices.
    spMat mass;       // mass(i,i) is the area of the voronoi cell of vertex vi
    Eigen::MatrixXd norm_e; // normal directions on each edge


    /*
        Initialize mesh properties
    */
    void get_mesh_angles();
    void get_mesh_normals_per_face();
    void get_mesh_normals_per_ver();
    // This is to calculate the 90 degree rotation matrices in
    // each triangle face. It can be used to calculate gradients.
    void get_face_rotation_matices();
    void get_rotated_edges_for_each_face();
    void gradient_v2f(std::array<spMat, 3> &output); // calculate gradient in each face from vertex values
    void gradient_f2v(spMat &output);                // calculate gradient in each vertex by averging face values
    void get_function_gradient_vertex();             // get the gradient of f on vertices
    void get_function_hessian_vertex();              // get the hessian matrix of f on vertices
    // get the gradients and hessians on each vertices using the assigned level-set function values
    // CAUTION: assign function values before call this function.
    void get_gradient_hessian_values(const Eigen::VectorXd& func, Eigen::MatrixXd& Vgrad, Eigen::MatrixXd& Fgrad);

    void get_vertex_rotation_matices();

    void get_I_and_II_locally();
    void get_bnd_vers_and_handles();
    void get_all_the_edge_normals();// the edge normals and the active edges
    
    /*
        tracing params
    */
    std::vector<CGMesh::HalfedgeHandle> tracing_start_edges;
    std::vector<std::vector<Eigen::Vector3d>> trace_vers;        // traced vertices
    std::vector<std::vector<CGMesh::HalfedgeHandle>> trace_hehs; // traced half edge handles
    std::vector<double> assigned_trace_ls; // function values for each traced curve.
    std::vector<std::vector<Eigen::Vector2d>> traced_paras;//parameterizations for the traced curves. This is only useful for illustrations.
    double trace_start_angle_degree = 90; //the angle between the first segment and the given boundary

    // void assemble_solver_pseudo_geodesic_part(spMat &H, Efunc& B);
    bool find_geodesic_intersection_p1_is_NOT_ver(const Eigen::Vector3d &p0, const Eigen::Vector3d &p1,
                                                  const CGMesh::HalfedgeHandle &edge_middle, const Eigen::Vector3d &pnorm, CGMesh::HalfedgeHandle &edge_out, Eigen::Vector3d &p_end);
    bool find_osculating_plane_intersection_not_geodesic_p1_is_not_ver(const Eigen::Vector3d &p0, const Eigen::Vector3d &p1,
                                                                       const CGMesh::HalfedgeHandle &edge_middle, const Eigen::Vector3d &pnorm, const double angle, std::vector<CGMesh::HalfedgeHandle> &edge_out, std::vector<Eigen::Vector3d> &p_end);
    bool get_pseudo_vertex_and_trace_forward(
        QuadricCalculator &cc,
        const std::vector<Eigen::Vector3d> &curve, std::array<Eigen::Vector3d, 3> &pcurve_local, const double angle_degree,
        const CGMesh::HalfedgeHandle &edge_middle,
        const Eigen::Vector3d &point_in, const Eigen::Vector3d &point_middle,
        const bool calculate_pseudo_vertex, CGMesh::HalfedgeHandle &edge_out,
        Eigen::Vector3d &point_out, bool &generate_pseudo_vertex, Eigen::Vector3d &pseudo_vertex_out);
    bool trace_pseudo_geodesic_forward(
        const NeighbourInfo &ninfo,
        const std::vector<Eigen::Vector3d> &curve, const double angle_degree,
        const CGMesh::HalfedgeHandle &edge_middle,
        const Eigen::Vector3d &point_in, const Eigen::Vector3d &point_middle,
        CGMesh::HalfedgeHandle &edge_out,
        Eigen::Vector3d &point_out);
    bool init_pseudo_geodesic_first_segment(
        const CGMesh::HalfedgeHandle &start_boundary_edge_pre, const double &start_point_para,
        const double start_boundary_angle_degree,
        CGMesh::HalfedgeHandle &intersected_handle,
        Eigen::Vector3d &intersected_point);
    
    bool trace_single_pseudo_geodesic_curve(const double target_angle_degree,
                                            const CGMesh::HalfedgeHandle &start_boundary_edge, const double &start_point_para,
                                            const double start_boundary_angle_degree,
                                            std::vector<Eigen::Vector3d> &curve,
                                            std::vector<CGMesh::HalfedgeHandle> &handles);
    bool get_checking_edges(const std::vector<int> &start_point_ids, const CGMesh::HalfedgeHandle &edge_middle, const Eigen::Vector3d &point_middle,
                            Efunc &edges_checked, Efunc &points_checked, NeighbourInfo &ninfo, std::vector<int> &point_to_check);
    void estimate_strip_width_according_to_tracing();




    /*
        Interactive design
    */
    std::vector<std::vector<int>> interactive_flist;
    std::vector<std::vector<Eigen::Vector3f>> interactive_bclist;

    /*
        level set opt functions and variables
    */
    std::array<LSAnalizer, 3> analizers;
    Eigen::VectorXd Glob_lsvars; // global variables for levelset opt
    // normal level set analyzer
    void analysis_pseudo_geodesic_on_vertices(const Eigen::VectorXd& func_values, Eigen::VectorXi& LocalActInner,
        std::vector<CGMesh::HalfedgeHandle>& heh0, std::vector<CGMesh::HalfedgeHandle>& heh1,
        std::vector<double>& t1s, std::vector<double>& t2s);
    void analysis_pseudo_geodesic_on_vertices(const Eigen::VectorXd& func_values, LSAnalizer& analizer);
    
    
    // level set angle controller
    void calculate_pseudo_geodesic_opt_expanded_function_values(Eigen::VectorXd &vars, const std::vector<double> &angle_degree,
                                                                const LSAnalizer &analizer, const int vars_start_loc, const int aux_start_loc, std::vector<Trip> &tripletes, Eigen::VectorXd &Energy);
    void calculate_extreme_pseudo_geodesic_values(Eigen::VectorXd &vars, const bool asymptotic,
                                                  const LSAnalizer &analizer, const int vars_start_loc, std::vector<Trip> &tripletes, Eigen::VectorXd &Energy);
    
    
    void assemble_solver_laplacian_part(spMat &H, Efunc &B);
    void assemble_solver_biharmonic_smoothing(const Eigen::VectorXd &func, spMat &H, Eigen::VectorXd &B); // Biharmonic smoothing (natural curved Hessian boundary)
    void assemble_solver_boundary_condition_part(const Eigen::VectorXd &func, spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &bcfvalue);
    void assemble_solver_interactive_boundary_condition_part(const Eigen::VectorXd &func, spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &bcfvalue);
    void assemble_solver_strip_width_part(const Eigen::MatrixXd& GradValue,  spMat& H, Eigen::VectorXd& B);
    void assemble_solver_extreme_cases_part_vertex_based(Eigen::VectorXd &vars, const bool asymptotic, const bool use_given_direction, 
                                                         const LSAnalizer &analizer, const int vars_start_loc, spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy);

    void assemble_solver_pesudo_geodesic_energy_part_vertex_based(Eigen::VectorXd &vars, const std::vector<double> &angle_degree,
                                                                  const LSAnalizer &analizer, const int vars_start_loc, const int aux_start_loc, spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy);

    
    void assemble_solver_othogonal_to_given_face_directions(const Eigen::VectorXd &func, const Eigen::MatrixXd &directions,
                                                            const Eigen::VectorXi &fids, spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy);
    // this function takes the reference level set as constant
    void assemble_solver_fix_angle_to_given_face_directions(const Eigen::VectorXd &func, const Eigen::MatrixXd &directions,
                                                            const Eigen::MatrixXd &grads, const double angle_fix,
                                                            const Eigen::VectorXi &fids, spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy);
    
    
    /*
        Shading design
    */
    void calculate_shading_condition_inequivalent(Eigen::VectorXd &vars,
                                                  const LSAnalizer &analizer, const int vars_start_loc, const int aux_start_loc, std::vector<Trip> &tripletes, Eigen::VectorXd &Energy);
    void calculate_shading_init(Eigen::VectorXd &vars,
                                const LSAnalizer &analizer, const int vars_start_loc, const int aux_start_loc, std::vector<Trip> &tripletes, Eigen::VectorXd &Energy);
    void calculate_binormal_regulizer(Eigen::VectorXd& vars,
	const LSAnalizer &analizer, const int vars_start_loc, const int aux_start_loc,std::vector<Trip> &tripletes,  Eigen::VectorXd& Energy);
    void assemble_solver_binormal_regulizer(Eigen::VectorXd &vars,
												 const LSAnalizer &analizer, const int vars_start_loc,
												 const int aux_start_loc, spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy);
   

   /*
        Topology: stitching boundary
   */
  // level set analyzer with detecting seams.
      std::vector<std::array<int, 2>> BndPairs; // pairs of overlapping boundary vertices.
    void analysis_pseudo_geodesic_on_vertices(const Eigen::VectorXd &func_values,
                                              const std::vector<std::array<int, 2>> boundary_pairs,
                                              Eigen::VectorXi &LocalActInner,
                                              Eigen::VectorXi &LocalActBpair,
                                              std::vector<bool> &correspondance,
                                              std::vector<CGMesh::HalfedgeHandle> &heh0, std::vector<CGMesh::HalfedgeHandle> &heh1,
                                              std::vector<double> &t1s, std::vector<double> &t2s);
    void analysis_pseudo_geodesic_on_vertices(const Eigen::VectorXd& func_values,const std::vector<std::array<int, 2>> boundary_pairs, LSAnalizer& analizer);
    
   void assemble_solver_pesudo_geodesic_energy_part_stitch_boundary(Eigen::VectorXd &vars, const std::vector<double> &angle_degree,
                                                                  const LSAnalizer &analizer, const int vars_start_loc, const int aux_start_loc, spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy);                                                                  

    

    /*
        level set structures
    */

    // this function takes the reference level set also as variables
    // fix the angle between the first and the second levelset
    void assemble_solver_fix_two_ls_angle(const Eigen::VectorXd &vars, const double angle_fix,
                                          spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy);


    /*
        extract polylines
    */
    void extract_one_curve(const double value, Eigen::MatrixXd& E0, Eigen::MatrixXd &E1);
    


    /*
        Triangle Mesh Optimization Part
    */
    Eigen::VectorXd Glob_Vars;
    spMat MVLap; // mean value laplacian
    spMat Elmat; // the edge length constriant
    double weight_Mesh_smoothness;
    double weight_Mesh_pesudo_geodesic;
    double weight_Mesh_edgelength;
public:
    igl::AABB<Eigen::MatrixXd, 3> aabbtree; // tree of the original mesh
    Eigen::VectorXd ElStored;               // edge length of the original mesh
    Eigen::MatrixXd Vstored;                // the stored vertices of the original mesh.
    Eigen::MatrixXd Nstored;                // stored normal vectors of the original mesh.
    double weight_Mesh_approximation;
    double weight_Mesh_mass;
private:
    double Mesh_opt_max_step_length;
    bool Last_Opt_Mesh=false; // the last step was mesh optimization. need to update all the mesh properties
    void calculate_mesh_opt_expanded_function_values( Eigen::VectorXd& vars,
        const LSAnalizer &analizer,
        const std::vector<double>& angle_degree,
        const int aux_start_loc, std::vector<Trip>& tripletes, Eigen::VectorXd& MTenergy);
    void calculate_mesh_opt_extreme_values(Eigen::VectorXd &vars, const int aux_start_loc, const Eigen::VectorXd &func, const bool asymptotic, const bool use_given_direction, const Eigen::Vector3d &ray,
                                           const LSAnalizer &analizer, std::vector<Trip> &tripletes, Eigen::VectorXd &MTenergy);
    void assemble_solver_mesh_opt_part( Eigen::VectorXd& vars,
        const LSAnalizer &analizer,
        const std::vector<double>& angle_degrees, const int aux_start_loc, spMat& JTJ, Eigen::VectorXd& B, Eigen::VectorXd& MTEnergy);
    
    void assemble_solver_mesh_smoothing(const Eigen::VectorXd &vars, spMat &H, Eigen::VectorXd &B);
    void calculate_mesh_opt_shading_condition_values(const Eigen::VectorXd &func, const Eigen::Vector3d &ray,
                                                const LSAnalizer &analizer, std::vector<Trip> &tripletes, Eigen::VectorXd &MTenergy);
    void assemble_solver_mean_value_laplacian(const Eigen::VectorXd& vars, spMat& H, Eigen::VectorXd& B);
    void assemble_solver_mesh_edge_length_part(const Eigen::VectorXd vars, spMat& H, Eigen::VectorXd& B, 
    Eigen::VectorXd& energy);
    void assemble_solver_mesh_extreme(Eigen::VectorXd &vars, const int aux_start_loc, const Eigen::VectorXd &func, const bool asymptotic, const bool use_given_direction,
                                      const Eigen::Vector3d &ray, const LSAnalizer &analizer,
                                      spMat &JTJ, Eigen::VectorXd &B, Eigen::VectorXd &MTEnergy);
    void assemble_solver_shading_mesh_opt(const Eigen::VectorXd &func, const Eigen::Vector3d &ray,
                                           const LSAnalizer &analizer,
                                           spMat &JTJ, Eigen::VectorXd &B, Eigen::VectorXd &MTEnergy);
    void assemble_solver_curve_smooth_mesh_opt(const LSAnalizer &analizer,
                                               spMat &JTJ, Eigen::VectorXd &B, Eigen::VectorXd &energy);
    void assemble_solver_approximate_original(spMat &JTJ, Eigen::VectorXd &B, Eigen::VectorXd &energy);
    // void assemble_solver_mesh_strip_width(spMat &JTJ, Eigen::VectorXd &B, Eigen::VectorXd &energy);
    void update_mesh_properties();// the mesh properties (normals, laplacian, etc.) get updated before optimization



    /*
    ///////////////functional angle values part//////////////////////
    */
    

    bool Analyze_Optimized_LS_Angles = false; // analyze the pg angles from the variables.
    Eigen::VectorXd AnalizedAngelVector;

    std::vector<double> changed_angles;
    Eigen::VectorXd RawTargetAngles;
    std::vector<int> PosFids;// the positive curved fids
    std::vector<int> NegFids;// the nagative curved fids
    std::vector<Eigen::Vector3d> OrthMinK;// the direction orthogonal to minimal curvature direction
    std::vector<double> CosSqr; // the squared cosin value of the target direction to CurvatureDirections[i][0]
    void assemble_solver_follow_min_abs_curvature(spMat &H, Eigen::VectorXd &B, Eigen::VectorXd& Eng);
    // disable changed angles and use the current version.
public:
    bool Disable_Changed_Angles = false;
    // By default false, meaning we still change the angles to make it close to asymptotics. 
    // if true, fit the angles as a polynormial curve, and approximate the curve values.
    bool Use_Fitting_Angles = false; 
    bool Use_Opt_Only_BNMS = false;
    bool Init_ChagingAngles = false;
    std::vector<std::vector<Eigen::Vector3d>> CurvatureDirections;
    void write_fitting_data();
    void show_minimal_curvature_directions(Eigen::MatrixXd& E0, Eigen::MatrixXd& E1, const double scaling);


    /*
        constant slope
    */
private:
    std::vector<int> Fslope; // the faces of the marked slopes.
    std::vector<Eigen::Vector3d> OrthoSlope; // the orthogonal vectors of the desired slopes
    Eigen::MatrixXd bxis0;
    Eigen::MatrixXd bxis1;
    // The following is for constant slope optimization
    void assemble_solver_constant_slope(const Eigen::VectorXd &vars,
                                        const bool fix_axis, const Eigen::Vector3d &axis_in, const bool fix_angle, const double angle_degree,
                                        const LSAnalizer &analizer,
                                        spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy);
    void assemble_solver_constant_slope_vertex_based(const Eigen::VectorXd &vars,
                                        const bool fix_axis, const Eigen::Vector3d &axis_in, const bool fix_angle, const double angle_degree,
                                        const LSAnalizer &analizer,
                                        spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy);
public:
    bool AxisFixedForSlopes = false;
    bool AnglesFixedForSlopes = false;
    Eigen::Vector3d AxisFixIn;
    double AngleFixIn;
    void show_slopes(const double scaleing, Eigen::MatrixXd &E0, Eigen::MatrixXd &E1);
    
    /*
        node design
    */

   // the rulings passing through the same point.
    void assemble_solver_co_point_ruling_node(const Eigen::VectorXd &vars,
                                                   const bool fix_centroid, const Eigen::Vector3d &centroid_in, 
												   const bool fix_radius, const std::vector<double>& radius_input,
                                                   const LSAnalizer &analizer,
                                                   spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy);
    void Run_ruling_opt();
    void statistic_tangent_vector_to_pt_distance();
    void test_quadratic_compare_0(const double a, const double b, const double c, const double xmin, const double xmax);
    

    bool RulingFixCentroid = false;
    bool RulingFixRadius = false;
    Eigen::Vector3d RulingCentroid;
    bool RulingOnlyFindingSphere = false; // take only centroid and radius as variables, keep the level sets unchanged
    Eigen::MatrixXd TangentE0;
    Eigen::MatrixXd TangentE1;
    Eigen::MatrixXd NodeProjectPts;
    bool write_ls_radius_file = false;
    std::vector<double> RadiusCollected; 
    std::vector<double> NodeFitPly;
    std::vector<double> NodeFitTargetRadius;
    Eigen::VectorXd XsNode;
    Eigen::VectorXd YsNode;

    /*
    //////////////////////  Interfaces  ////////////////////
    */
public:
    // parametrization and find the boundary loop
    lsTools(CGMesh &mesh);
    lsTools(){};
    void init(CGMesh &mesh);
    CGMesh lsmesh;                        // the input mesh
    Eigen::MatrixXd V;
    Eigen::MatrixXd norm_v;               // normal perf vertex
    Eigen::MatrixXi F;
    Eigen::MatrixXi E;
    Eigen::VectorXi bnd;     // boundary loop
    Eigen::MatrixXd paras;   // parameters of the mesh vertices, nx2
    Eigen::VectorXd fvalues; // function values
    Eigen::VectorXi InnerV; // the vector show if it is a inner ver. size is vnbr
    std::vector<int> IVids; // the ids of the inner vers, size is ninner
    spMat mass_uniform; // the uniformed mass matrix to make sure the average area is 1.

    double weight_mass;              // weight of the mass function energy
    double weight_laplacian;              // weight of the laplacian energy
    double weight_boundary;              // weight of the boundary energy
    double weight_pseudo_geodesic_energy;// weight of the pseudo-geodesic energy
    double weight_strip_width;
    double weight_geodesic; // For AAG and shading
    double weight_smt_binormal = 0; // smooth the binormals. Please only use it for shading system
    // double weight_shading; // For shading: the light is othogonal to the tangent direction 
    double strip_width = 1;       // strip width, defined as h/w, h: level set function value difference. w: distance between two points on the surface
    double max_step_length;
    double step_length;
    double pgerror;
    
    bool enable_pseudo_geodesic_energy=false; // decide if we include pseudo-geodesic energy
    bool enable_strip_width_energy=false;
    
    bool enable_extreme_cases = false; // osculating plane othogonal to the normal direction, or contains the given direction
    bool Given_Const_Direction = false;// use the given direction for shading 
    bool enable_shading_init = false;
    bool enable_let_ray_through = false; // optimize that whole surface let ray through for shading
    bool enable_reflection = false;
    bool recompute_auxiliaries = false;// recompute auxiliaries to reduce the energy.
    bool fix_angle_of_two_levelsets = false;
    double angle_between_two_levelsets;
    double weight_fix_two_ls_angle;

    double Reference_theta;  // the ray feed to the optimization as a constant direction
    double Reference_phi;
    double Reference_theta1;  // the ray feed to the optimization as a constant direction
    double Reference_phi1;
    double Theta_tol; // the tolerances for theta and phi
    double Phi_tol;
    double ShadingLatitude;
    double Theta_tol1;
    double Phi_tol1;

    double weight_binormal;
    Eigen::VectorXi shading_condition_info;
    std::vector<std::vector<int>> PG_Vertex_Types;// record the vertices types used in variational angle
    Eigen::VectorXd PGGeoDistance;// TODO maybe not needed since this method is bad.
    Eigen::VectorXd PGVariationalAngles;// TODO maybe not needed since this method is bad.
    
    double pseudo_geodesic_target_angle_degree; // the target pseudo-geodesic angle
    double pseudo_geodesic_target_angle_degree_2; // the target pseudo-geodesic angle
    // bool enable_pseudo_vertex=false;

    // Eigen::MatrixXd vBinormal;

    std::vector<CGMesh::HalfedgeHandle> Boundary_Edges; // Caution: Either it or it's opposite handle is the boundary halfedge handle

    bool Compute_Auxiliaries=true;
    bool Compute_Auxiliaries_Mesh = true;
    Eigen::MatrixXd Binormals;
    Eigen::MatrixXd Lights;
    Eigen::VectorXd PGE;// pseudo geodesic energy

    std::vector<int> refids;                      // output the ids of current dealing points. just for debug purpose
    Eigen::MatrixXd Ppro0;// the projected corresponding points of the mesh to the reference mesh
    Eigen::MatrixXd Npro0;// the normal vector of the projected corresponding points on the reference mesh
    
    // Eigen::Vector3d ray_catcher(int vid);
    // boundary conditions
    // int boundary_type = 
    // 
    // this function should be calculated first once the class get constructed
    //  1. get face normals;
    //  2. get all the angles;
    //  3. get vertex normals using angles and face normals;
    //  4. get the face rotation matrices.
    //  5. get the rotated edges in each face.
    //  6. get the gradients on vertices
    //  7. get the Hessian of function f on vertices
    //  8. get rotated edges in parametric domain
    //  9. get the I and II fundamental forms
    void initialize_mesh_properties()
    {
        // 1
        get_mesh_normals_per_face();
        // 2
        get_mesh_angles();
        // // 3
        get_mesh_normals_per_ver();
        // 4
        get_face_rotation_matices();

        // 5
        get_rotated_edges_for_each_face();
        // 6
        get_function_gradient_vertex();

        // 7
        get_function_hessian_vertex();

        // 8
        // get_rotated_parameter_edges();
        // 9
        // get_I_and_II_locally();
        get_vertex_rotation_matices();
        get_bnd_vers_and_handles();
        get_all_the_edge_normals();
        
    }
    void prepare_level_set_solving(const EnergyPrepare &Energy_initializer);
    void prepare_mesh_optimization_solving(const MeshEnergyPrepare& initializer);
    
    void show_level_set(Eigen::VectorXd &val);

    void show_gradients(Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, double ratio);
    void show_gradient_scalar(Eigen::VectorXd &values, int which);
    void show_face_gradients(Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, double ratio);
    void show_current_reference_points(Eigen::MatrixXd &pts);
    void show_1_order_derivate(Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, Eigen::MatrixXd &E2, Eigen::MatrixXd &E3, double ratio);
    void show_vertex_normal(Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, double ratio);
    void show_binormals(const Eigen::VectorXd &func, Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, Eigen::MatrixXd &binormals, double ratio);
    void show_max_pg_energy(Eigen::VectorXd& e);
    void show_max_pg_energy_all(Eigen::MatrixXd &energy);
    void show_posi_negt_curvature(Eigen::MatrixXd& color);
    void extract_levelset_curves(const int nbr, std::vector<Eigen::MatrixXd> &E0, std::vector<Eigen::MatrixXd> &E1);
    // works only for const slope method
    void show_bxis(Eigen::MatrixXd& E0, Eigen::MatrixXd& E1);
   
   
    // show the traced pg curves
    void show_pseudo_geodesic_curve(std::vector<Eigen::MatrixXd> &E0, std::vector<Eigen::MatrixXd> &E1, Eigen::MatrixXd &vers);
    void show_traced_binormals(Eigen::MatrixXd &bE0, Eigen::MatrixXd &bE1, Eigen::MatrixXd &nE0, Eigen::MatrixXd &nE1, const double scale);
    void show_traced_curve_params(Eigen::MatrixXd &curve);
    void save_traced_curves(const std::string filename);
    void clear_traced_curves();
    void load_one_traced_curve();
    void Trace_One_Guide_Pseudo_Geodesic();
    
    
    // after tracing, use this function to get smooth level set
    void Run_Level_Set_Opt();
    void Run_Level_Set_Opt_interactive(const bool compute_pg);
    void Run_Level_Set_Opt_Angle_Variable();
    void Run_Mesh_Opt();
    void Run_Mesh_Smoothness();
    void Run_AAG(Eigen::VectorXd& func0, Eigen::VectorXd& func1, Eigen::VectorXd& func2);
    void Run_AAG_Mesh_Opt(Eigen::VectorXd& func0, Eigen::VectorXd& func1, Eigen::VectorXd& func2);
    void Run_AGG(Eigen::VectorXd& func0, Eigen::VectorXd& func1, Eigen::VectorXd& func2);
    void Run_AGG_Mesh_Opt(Eigen::VectorXd& func0, Eigen::VectorXd& func1, Eigen::VectorXd& func2);
    void Run_PPG(Eigen::VectorXd& func0, Eigen::VectorXd& func1, Eigen::VectorXd& func2);
    void Run_PPG_Mesh_Opt(Eigen::VectorXd& func0, Eigen::VectorXd& func1, Eigen::VectorXd& func2);
    void Run_Othogonal_Levelset(const Eigen::VectorXd &func_ref);
    void Run_AsOrthAsPossible_LS();
    void RunInitAOAP();
    void Run_ConstSlopeOpt();

    void initialize_level_set_accroding_to_parametrization();
    // the following two functions actually trace curves or segments and assign values on them
    void initialize_level_set_by_tracing(const TracingPrepare& Tracing_initializer);
    void initialize_level_set_by_select_boundary_segment(const TracingPrepare& Tracing_initializer, bool trace);
    bool receive_interactive_strokes_and_init_ls(const std::vector<std::vector<int>> &flist,
                                          const std::vector<std::vector<Eigen::Vector3f>> &bclist);

    CGMesh write_parameterization_mesh();
    void debug_tool();
    void print_info(const int vid);
};

// // partial derivate tools, regard level set function as variates.
// class PDtools{
//     PDtools(){};
// };
