#pragma once
#include <lsc/MeshProcessing.h>
// theta <= pi/2, phi<=pi
 void sphere_example(double radius, double theta, double phi, int nt, int np);
// Efunc represent a elementary value, which is the linear combination of
// some function values on their corresponding vertices, of this vertex.
// i.e. Efunc[0]*f_0+Efunc[1]*f_1+...+Efunc[n]*f_n.
typedef Eigen::SparseVector<double> Efunc;
#define SCALAR_ZERO 1e-8
class Vectorlf
{
private:
    std::vector<Efunc> mat; // the vector of LCF, for each vertex there is a LCF

public:
    Vectorlf(){};
    int size() const 
    {
        
        return mat.size();
    }
    int size() 
    {
        
        return mat.size();
    }

    void resize(int sz)
    {
        mat.resize(sz);
    }
    void resize(int sz1, int sz2)
    {
        mat.resize(sz1);
        for (int i = 0; i < sz1; i++)
        {
            mat[i].resize(sz2);
        }
    }
    void clear()
    {
        mat.clear();
    }

    Efunc operator()(int id) const
    {
        assert(id >= 0 && id < mat.size());
        return mat[id];
    }
    Efunc &operator()(int id)
    {
        assert(id >= 0 && id < mat.size());
        return mat[id];
    }
};

// The basic tool of LSC. Please initialize it with a mesh
class lsTools
{
private:
    MeshProcessing MP;
    CGMesh lsmesh;                       // the input mesh
    Eigen::MatrixXd norm_f;              // normal per face
    Eigen::MatrixXd norm_v;              // normal perf vertex
    Eigen::MatrixXd angF;                // angel per vertex of each F. nx3, n is the number of faces
    Eigen::VectorXd areaF;               // the area of each face.
    Eigen::VectorXd areaPF;              // the area of each face in parametric domain;
    std::vector<Eigen::Matrix3d> Rotate; // the 90 degree rotation matrices for each face.
    std::vector<Eigen::Matrix3d> RotateV;// on vertices
    // the rotated 3 half edges for each face, in the plane of this face
    // the order of the 3 directions are: v0-v2, v1-v0, v2-v1
    std::vector<std::array<Eigen::Vector3d, 3>> Erotate;
    // 2d rotated half edges in parametric domain
    std::vector<std::array<Eigen::Vector2d, 3>> Erotate2d;
    std::array<Vectorlf, 3> gradVF;                  // gradient of function in each face
    std::array<Vectorlf, 3> gradV;                   // gradient of function in each vertex
    std::array<std::array<Vectorlf, 3>, 3> HessianV; // Hessian (2 order deriavate) of function in each vertex
    std::array<Eigen::MatrixXd, 2> Deriv1;           // the 1 order derivates for each vertex;
    std::array<Eigen::MatrixXd, 4> Deriv2;           // the 2 order derivates for each vertex;
    std::vector<double> II_L;                        // second fundamental form: L
    std::vector<double> II_M;                        // second fundamental form: M
    std::vector<double> II_N;                        // second fundamental form: N
    Eigen::MatrixXd gvvalue;                         // the calculated gradient values of f for each vertex.
    Eigen::MatrixXd gfvalue;                         // the calculated gradient values of f for each face.
    std::vector<Eigen::Matrix3d> hfvalue;             // the calculated Hessian values of f for each vertex
    std::vector<int> refids;// output the ids of current dealing points. just for debug purpose
    // Eigen::MatrixXd Fcenters;// face centers, used for debug;
        void
        get_mesh_angles();
    void get_mesh_normals_per_face();
    void get_mesh_normals_per_ver();
    // This is to calculate the 90 degree rotation matrices in
    // each triangle face. It can be used to calculate gradients.
    void get_face_rotation_matices();
    void get_rotated_edges_for_each_face();
    void gradient_v2f(Vectorlf &values, std::array<Vectorlf, 3> &output);                // calculate gradient in each face from vertex values
    void gradient_f2v(std::array<Vectorlf, 3> &values, std::array<Vectorlf, 3> &output); // calculate gradient in each vertex by averging face values
    void initialize_function_coefficient(Vectorlf &input);                               // initialize the coefficients as an identity matrix format.
    void gradient_easy_interface(Vectorlf &func, std::array<Vectorlf, 3> &gonf, std::array<Vectorlf, 3> &gonv);
    void get_function_gradient_vertex();                                                                         // get the gradient of f on vertices
    void get_function_hessian_vertex();                                                                          // get the hessian matrix of f on vertices
    void get_rotated_parameter_edges();                                                                          // get the rotated edges on the 2d parametric domain
    void surface_derivate_v2f(const Eigen::MatrixXd &vvalues, std::array<Eigen::MatrixXd, 2> &DonF);             // calculate gradient in each face from vertex values
    void surface_derivate_f2v(const std::array<Eigen::MatrixXd, 2> &DonF, std::array<Eigen::MatrixXd, 2> &DonV); // calculate gradient in each vertex by averging face values
    void surface_derivate_easy_interface(const Eigen::MatrixXd &vvalues, std::array<Eigen::MatrixXd, 2> &DonV);
    void get_surface_derivate();        // get the derivates
    void get_surface_II_each_ver();     // get the second fundamental forms L, M and N
    void get_gradient_hessian_values(); // get the gradients and hessians on each vertices using the assigned level-set function values
    void get_lf_value(const Vectorlf& coff, Eigen::VectorXd& res); // get the linear combinations of function values.
    void make_sphere_ls_example(int rowid);
    void get_vertex_rotation_matices();
    
        
public:
    // parametrization and find the boundary loop
    lsTools(CGMesh &mesh);

    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXi E;
    Eigen::VectorXi bnd;     // boundary loop
    Eigen::MatrixXd paras;   // parameters of the mesh vertices, nx2
    Eigen::VectorXd fvalues; // function values
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
        // 3
        get_mesh_normals_per_ver();
        // 4
        get_face_rotation_matices();
        get_vertex_rotation_matices();
        // 5
        get_rotated_edges_for_each_face();
        // 6
        get_function_gradient_vertex();
        // 7
        get_function_hessian_vertex();
        // 8
        get_rotated_parameter_edges();
        // 9
        get_surface_derivate();
        get_surface_II_each_ver();
    }
    void show_level_set(Eigen::VectorXd &val); 
    // convert the parameters into a mesh, for visulization purpose
    void convert_paras_as_meshes(CGMesh &output);
    void show_gradients(Eigen::MatrixXd& E0, Eigen::MatrixXd &E1, double ratio);
    void show_face_gradients(Eigen::MatrixXd& E0, Eigen::MatrixXd &E1, double ratio);
    void show_current_reference_points(Eigen::MatrixXd& pts);
    void show_face_grad_max_angle(Eigen::MatrixXd& points);
    void debug_tool(int id = 0, double value = 0)
    {
        make_sphere_ls_example(id);
        

        // show gradients as edges
        
        //sphere_example(value, 1.3, 1.8, id, id);
        // Efunc vec, vec1;
        // vec.resize(3);
        // vec1.resize(3);
        // std::cout << "print the test sparse vector\n"
        //           << vec << std::endl;
        // std::cout << "size " << vec.size() << std::endl;
        // std::cout << "values " << vec.coeffRef(0) << std::endl;
        // vec.coeffRef(0) = 1.0;
        // vec.coeffRef(2) = 2.0;
        // vec1.coeffRef(0) = 1.0;
        // vec = vec * 2;
        // vec += vec1;

        // std::cout << "print the test sparse vector\n"
        //           << vec << std::endl;
        // std::cout << "size " << vec.size() << std::endl;
        // std::cout << "values " << vec.coeffRef(0) << std::endl;

        // std::cout << " testing Vectorlf" << std::endl;
        // Vectorlf newvec;
        // newvec.resize(1, 3);
        // newvec(0) = vec;
        // std::cout << "Vectorlf value is " << newvec(0) << std::endl;
    }
};