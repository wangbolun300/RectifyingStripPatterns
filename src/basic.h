#pragma once
#include <lsc/MeshProcessing.h>
// theta <= pi/2, phi<=pi
void sphere_example(double radius, double theta, double phi, int nt, int np);
// Efunc represent a elementary value, which is the linear combination of
// some function values on their corresponding vertices, of this vertex.
// i.e. Efunc[0]*f_0+Efunc[1]*f_1+...+Efunc[n]*f_n.
typedef Eigen::SparseVector<double> Efunc;
// typedef Eigen::SparseMatrix<double> spMat;
typedef Eigen::Triplet<double> Trip;
#define SCALAR_ZERO 1e-8
// class Vectorlf
// {
// private:
//     std::vector<Efunc> mat; // the vector of LCF, for each vertex there is a LCF

// public:
//     Vectorlf(){};
//     // int size() const
//     // {

//     //     return mat.size();
//     // }
//     // int size()
//     // {

//     //     return mat.size();
//     // }

//     // void resize(int sz)
//     // {
//     //     mat.resize(sz);
//     // }
//     // void resize(int sz1, int sz2)
//     // {
//     //     mat.resize(sz1);
//     //     for (int i = 0; i < sz1; i++)
//     //     {
//     //         mat[i].resize(sz2);
//     //     }
//     // }
//     // void clear()
//     // {
//     //     mat.clear();
//     // }

//     // Efunc operator()(int id) const
//     // {
//     //     assert(id >= 0 && id < mat.size());
//     //     return mat[id];
//     // }
//     // Efunc &operator()(int id)
//     // {
//     //     assert(id >= 0 && id < mat.size());
//     //     return mat[id];
//     // }
// };

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
    std::vector<Eigen::Matrix3d> RotateV; // on vertices
    // the rotated 3 half edges for each face, in the plane of this face
    // the order of the 3 directions are: v0-v2, v1-v0, v2-v1
    std::vector<std::array<Eigen::Vector3d, 3>> Erotate;
    // 2d rotated half edges in parametric domain
    std::vector<std::array<Eigen::Vector2d, 3>> Erotate2d;
    std::array<spMat, 3> gradVF;                  // gradient of function in each face
    std::array<spMat, 3> gradV;                   // gradient of function in each vertex
    std::array<std::array<spMat, 3>, 3> HessianV; // Hessian (2 order deriavate) of function in each vertex
    std::array<Eigen::MatrixXd, 2> Deriv1;           // the 1 order derivates for each vertex;
    std::array<Eigen::MatrixXd, 4> Deriv2;           // the 2 order derivates for each vertex;
    std::vector<double> II_L;                        // second fundamental form: L
    std::vector<double> II_M;                        // second fundamental form: M
    std::vector<double> II_N;                        // second fundamental form: N
    Eigen::MatrixXd gvvalue;                         // the calculated gradient values of f for each vertex.
    Eigen::MatrixXd gfvalue;                         // the calculated gradient values of f for each face.
    std::vector<Eigen::Matrix3d> hfvalue;            // the calculated Hessian values of f for each vertex
    std::vector<int> refids;                         // output the ids of current dealing points. just for debug purpose
    // spMat Dlpsqr;                       // the derivates of ||laplacian F||^2 for all the vertices.
    spMat Dlps;                         // the derivates of laplacian F for all the vertices.
    // std::array<spMat,3> Dgrad;                        // the derivates of gradients for all the vertices.
    // std::vector<spMat> Dhess;                        // the derivates of hessian for all the vertices.
    // spMat Dgrad_norm;                   // the derivates of the norm of gradients for all the vertices.
    spMat mass;                                      // mass(i,i) is the area of the voronoi cell of vertex vi
    // spMat F2V;                                       // the matrix averaging face values to vertices values, acorrding to the angels

    bool derivates_calculated = false;
    void
    get_mesh_angles();
    void get_mesh_normals_per_face();
    // void get_mesh_normals_per_ver();
    // This is to calculate the 90 degree rotation matrices in
    // each triangle face. It can be used to calculate gradients.
    void get_face_rotation_matices();
    void get_rotated_edges_for_each_face();
    void gradient_v2f(std::array<spMat, 3> &output);                // calculate gradient in each face from vertex values
    void gradient_f2v(spMat& output); // calculate gradient in each vertex by averging face values
    void get_function_gradient_vertex(); // get the gradient of f on vertices
    void get_function_hessian_vertex();  // get the hessian matrix of f on vertices
    // get the gradients and hessians on each vertices using the assigned level-set function values
    // CAUTION: assign function values before call this function.
    void get_gradient_hessian_values();
    void make_sphere_ls_example(int rowid);
    void get_vertex_rotation_matices();

    void get_I_and_II_locally();
    // get the partial derivate of gradient gradF on each vertex i
    void get_gradient_partial_cofficient_matrix(int i);
    // get the partial derivate cofficient matrix of the hessian on each vertex i
    void get_hessian_partial_cofficient_matrix(int i);
    // the partial derivate of ||gradF||: partial{||gradF||}{fj} on each vertex vi. takes partial matrix of gradient GP of vertex i as input
    void get_gradient_norm_partial_cofficient_matrix(int i, const spMat &GP);

    // the partial derivate of laplacian(F on vertex vi. takes partial matrix of hessian HP as input
    void get_laplacian_partial_cofficient_matrix(int i, const spMat &HP);
    // the partial derivate of ||laplacian(F)||^2 on vertex vi. takes partial matrix of hessian HP as input
    // void get_laplacian_square_partial_cofficient_matrix(int i, const spMat &HP);
    // CAUTION: please call this after you initialize the function values.
    void get_all_the_derivate_matrices();
    void assemble_solver_laplacian_part(spMat &H, Efunc &B);

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
    double weight_mass;              // weight of the mass function to
    int assign_face_id;              // the face id of which we will assign value to
    double assign_value[3];
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
    }
    void show_level_set(Eigen::VectorXd &val);
    // convert the parameters into a mesh, for visulization purpose
    void convert_paras_as_meshes(CGMesh &output);
    void show_gradients(Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, double ratio);
    void show_face_gradients(Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, double ratio);
    void show_current_reference_points(Eigen::MatrixXd &pts);
    void show_1_order_derivate(Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, Eigen::MatrixXd &E2, Eigen::MatrixXd &E3, double ratio);
    void show_face_1_order_derivate(Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, Eigen::MatrixXd &E2, Eigen::MatrixXd &E3, double ratio);
    void show_vertex_normal(Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, double ratio);
    void initialize_and_smooth_level_set_by_laplacian();
    void debug_tool(int id = 0, double value = 0)
    {
        // sphere_example(5, 80.0/180*3.1415926, 1.8, 17, 20);
        make_sphere_ls_example(id);
    }
};

std::vector<Trip> to_triplets(spMat &M);
// convert a vector to a matrix with size 1*n
spMat sparse_vec_to_sparse_maxrix(Efunc &vec);
Efunc sparse_mat_col_to_sparse_vec(const spMat &mat, const int col);
Efunc dense_vec_to_sparse_vec(const Eigen::VectorXd& vec);
// // partial derivate tools, regard level set function as variates.
// class PDtools{
//     PDtools(){};
// };