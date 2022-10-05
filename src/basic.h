#pragma once
#include <lsc/MeshProcessing.h>
#include <lsc/igl_tool.h>
#include <lsc/basic.h>

// Efunc represent a elementary value, which is the linear combination of
// some function values on their corresponding vertices, of this vertex.
// i.e. Efunc[0]*f_0+Efunc[1]*f_1+...+Efunc[n]*f_n.
typedef Eigen::SparseVector<double> Efunc;
typedef Eigen::SparseVector<int> SpVeci;
// typedef Eigen::SparseMatrix<double> spMat;
typedef Eigen::Triplet<double> Trip;
#define SCALAR_ZERO 1e-8
#define MERGE_VERTEX_RATIO 0.01
#define LSC_PI 3.14159265
#define ANGLE_TOLERANCE 2.
#define QUADRANT_TOLERANCE 0.04 // over 2 degree tolerance
class TracCurve
{
public:
    TracCurve(){};

    std::vector<int> face_ids;
    std::vector<Eigen::Vector3d> dirc;        // the directions of the tracing curve
    std::vector<int> edge_ids;                // ids of the corresponding edges;
    std::vector<Eigen::Vector3d> edge_points; // intersection points on edges.
    std::vector<double> func_values;
    bool size_correct();
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
// The basic tool of LSC. Please initialize it with a mesh
class lsTools
{
private:
    MeshProcessing MP;
    CGMesh lsmesh;                        // the input mesh
    Eigen::MatrixXd norm_f;               // normal per face
    Eigen::MatrixXd norm_v;               // normal perf vertex
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
    Eigen::MatrixXd gvvalue;                      // the calculated gradient values of f for each vertex.
    Eigen::MatrixXd gfvalue;                      // the calculated gradient values of f for each face.
    std::vector<Eigen::Matrix3d> hfvalue;         // the calculated Hessian values of f for each vertex
    std::vector<int> refids;                      // output the ids of current dealing points. just for debug purpose
    // spMat Dlpsqr;                       // the derivates of ||laplacian F||^2 for all the vertices.
    spMat Dlps;       // the derivates of laplacian F for all the vertices.
    spMat Dgrad_norm; // the derivates of the norm of gradients for all the vertices.
    spMat mass;       // mass(i,i) is the area of the voronoi cell of vertex vi
    Efunc bcfvalue; // boundary condition function values
    
    // spMat F2V;                                       // the matrix averaging face values to vertices values, acorrding to the angels
    spMat ORB; // the matrix where the (i,i) element show if it is a boundary vertex or one-ring vertex of boundary
    std::vector<CGMesh::HalfedgeHandle> Boundary_Edges;
    Eigen::MatrixXd norm_e; // normal directions on each edge
    std::vector<TracCurve> trace;
    std::vector<std::vector<Eigen::Vector3d>> trace_vers;        // traced vertices
    std::vector<std::vector<CGMesh::HalfedgeHandle>> trace_hehs; // traced half edge handles
    std::vector<double> assigned_trace_ls;
    std::vector<double> func_values;// function values for each traced curve.
    
    Efunc ActE;// active edges, which means they are not co-planar edges, and not boundary edges.
    bool derivates_calculated = false;
    void
    get_mesh_angles();
    void get_mesh_normals_per_face();
    // void get_mesh_normals_per_ver();
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
    void get_gradient_hessian_values();

    void get_vertex_rotation_matices();

    void get_I_and_II_locally();

    // CAUTION: please call this after you initialize the function values.
    void get_all_the_derivate_matrices();
    // get the boundary vertices and one ring vertices from them
    void get_bnd_and_bnd_one_ring();
    void get_all_the_edge_normals();// the edge normals and the active edges
    void calculate_gradient_partial_parts();// calculate gradient partial derivatives. Edge based
    void assemble_solver_boundary_condition_part(spMat& H, Efunc& B);
    void assemble_solver_laplacian_part(spMat &H, Efunc &B);
    void assemble_solver_strip_width_part(spMat &H, Efunc &B);

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
    // start_point_para is the parameter t, start_point=from+t*(to-from).
    // start_angle_degree is the angle between the initial direction and the direction from->to
    bool trace_single_pseudo_geodesic_curve_pseudo_vertex_method(const double target_angle,
                                                                 const CGMesh::HalfedgeHandle &start_boundary_edge, const double &start_point_para,
                                                                 const double start_angle_degree,
                                                                 std::vector<Eigen::Vector3d> &curve,
                                                                 std::vector<CGMesh::HalfedgeHandle> &handles);
    bool trace_single_pseudo_geodesic_curve(const double target_angle_degree,
                                            const CGMesh::HalfedgeHandle &start_boundary_edge, const double &start_point_para,
                                            const double start_boundary_angle_degree,
                                            std::vector<Eigen::Vector3d> &curve,
                                            std::vector<CGMesh::HalfedgeHandle> &handles);

    bool get_checking_edges(const std::vector<int> &start_point_ids, const CGMesh::HalfedgeHandle &edge_middle, const Eigen::Vector3d &point_middle,
                            Efunc &edges_checked, Efunc &points_checked, NeighbourInfo &ninfo, std::vector<int> &point_to_check);

public:
    // parametrization and find the boundary loop
    lsTools(CGMesh &mesh);
    lsTools(){};
    void init(CGMesh &mesh);
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXi E;
    Eigen::VectorXi bnd;     // boundary loop
    Eigen::MatrixXd paras;   // parameters of the mesh vertices, nx2
    Eigen::VectorXd fvalues; // function values
    double level_set_step_length;
    double weight_assign_face_value; // weight of the assigned value shown in the solver
    double weight_mass;              // weight of the mass function energy
    double weight_laplacian;              // weight of the laplacian energy
    double weight_boundary;              // weight of the boundary energy
    int assign_face_id;              // the face id of which we will assign value to
    double assign_value[3];
    double strip_width = 0;       // strip width, defined as h/w.
    double pseudo_geodesic_ratio; // k_g/k_n of the pseudo geodesic
    // bool enable_pseudo_vertex=false;
    Eigen::MatrixXd ver_dbg;
    Eigen::MatrixXd ver_dbg1;
    bool flag_dbg = false;
    int id_dbg;
    Eigen::MatrixXd E0_dbg;
    Eigen::MatrixXd direction_dbg;
    Eigen::MatrixXd pnorm_dbg;
    std::vector<Eigen::Vector3d> pseudo_vers_dbg;
    std::vector<Eigen::Vector3d> pnorm_list_dbg;
    Eigen::MatrixXd visual_pts_dbg;

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
        // get_mesh_normals_per_ver();
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
        get_I_and_II_locally();
        get_vertex_rotation_matices();
        get_bnd_and_bnd_one_ring();
        get_all_the_edge_normals();
    }
    void show_level_set(Eigen::VectorXd &val);
    // convert the parameters into a mesh, for visulization purpose
    void convert_paras_as_meshes(CGMesh &output);

    void show_gradients(Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, double ratio);
    void show_hessian(Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, double ratio, int which);
    void show_gradient_scalar(Eigen::VectorXd &values, int which);
    void show_face_gradients(Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, double ratio);
    void show_current_reference_points(Eigen::MatrixXd &pts);
    void show_1_order_derivate(Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, Eigen::MatrixXd &E2, Eigen::MatrixXd &E3, double ratio);
    void show_face_1_order_derivate(Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, Eigen::MatrixXd &E2, Eigen::MatrixXd &E3, double ratio);
    void show_vertex_normal(Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, double ratio);
    void show_pseudo_geodesic_curve(std::vector<Eigen::MatrixXd> &E0, std::vector<Eigen::MatrixXd> &E1, Eigen::MatrixXd &vers);

    void initialize_and_smooth_level_set_by_laplacian();
    void initialize_and_optimize_strip_width();
    void initialize_and_optimize_pseudo_geodesic();
    // after tracing, use this function to get smooth level set
    void optimize_laplacian_with_traced_boundary_condition();
    void initialize_level_set_accroding_to_parametrization();
    void debug_tool(int id = 0, int id2 = 0, double value = 0);
    void debug_tool_v2(const std::vector<int> &ids, const std::vector<double> values);
    void debug_tool_v3(int id = 0, int id2 = 0, double value = 0);
    void debug_tool_v4(const std::vector<int> &ids, const std::vector<double> values);
    void initialize_level_set_by_tracing(const std::vector<int> &ids, const std::vector<double> values);
    void make_sphere_ls_example(int rowid);
};

// // partial derivate tools, regard level set function as variates.
// class PDtools{
//     PDtools(){};
// };