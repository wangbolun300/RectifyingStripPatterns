#pragma once
#include <lsc/basic.h>
std::vector<Trip> to_triplets(spMat &M);
// convert a vector to a matrix with size 1*n
spMat sparse_vec_to_sparse_maxrix(Efunc &vec);
Efunc sparse_mat_col_to_sparse_vec(const spMat &mat, const int col);
Efunc dense_vec_to_sparse_vec(const Eigen::VectorXd &vec);
void mat_col_to_triplets(const spMat &mat, const int col, const int ref, const bool inverse, std::vector<Trip> &triplets);
spMat dense_mat_list_as_sparse_diagnal(const std::vector<Eigen::MatrixXd> &mlist);
Eigen::MatrixXd vec_list_to_matrix(const std::vector<Eigen::Vector3d> &vec);

bool triangles_coplanar(const Eigen::Vector3d &t0, const Eigen::Vector3d &t1, const Eigen::Vector3d &t2,
                        const Eigen::Vector3d &p0, const Eigen::Vector3d &p1, const Eigen::Vector3d &p2);

bool triangles_coplanar(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const int fid1, const int fid2);
bool segments_colinear(const Eigen::Vector3d& s1, const Eigen::Vector3d& s2);
// -1: boundary edge
// 0 not coplanar
// 1 co-planar
int edge_is_coplanar(const CGMesh &lsmesh, const CGMesh::EdgeHandle &edge_middle, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F);
std::vector<double> polynomial_simplify(const std::vector<double> &poly);
std::vector<double> polynomial_add(const std::vector<double> &poly1, const std::vector<double> &poly2);
std::vector<double> polynomial_times(const std::vector<double> &poly1, const std::vector<double> &poly2);
std::vector<double> polynomial_times(const std::vector<double> &poly1, const double &nbr);
double polynomial_value(const std::vector<double> &poly, const double para);
std::vector<double> polynomial_integration(const std::vector<double> &poly);
double polynomial_integration(const std::vector<double> &poly, const double lower, const double upper);

void cylinder_example(double radius, double height, int nr, int nh);
void cylinder_open_example(double radius, double height, int nr, int nh);
// theta <= pi/2, phi<=pi
void sphere_example(double radius, double theta, double phi, int nt, int np);

void split_mesh_boundary_by_corner_detection(CGMesh& lsmesh, const Eigen::MatrixXd& V, const double threadshold_angel_degree,
const std::vector<CGMesh::HalfedgeHandle> &Boundary_Edges, std::vector<std::vector<CGMesh::HalfedgeHandle>>& boundaries);

bool save_levelset(const Eigen::VectorXd &ls, const Eigen::MatrixXd& binormals);
bool save_levelset(const Eigen::VectorXd &ls);
bool read_levelset(Eigen::VectorXd &ls);
bool read_bi_normals(Eigen::MatrixXd &bn);
std::array<int,4> get_vers_around_edge(CGMesh& lsmesh, int edgeid, int& fid1, int &fid2);
double get_t_of_segment(const Eigen::Vector3d &ver, const Eigen::Vector3d &start, const Eigen::Vector3d &end);
double get_t_of_value(const double &ver, const double &start, const double &end);
Eigen::Vector2d get_2d_ver_from_t(const double t, const Eigen::Vector2d& start, const Eigen::Vector2d& end);
Eigen::Vector3d get_3d_ver_from_t(const double t, const Eigen::Vector3d& start, const Eigen::Vector3d& end);
bool angles_match(const double angle_degree1, const double angle_degree2);
void assign_angles_based_on_funtion_values(const Eigen::VectorXd &fvalues,
                                            const double angle_large_function, const double angle_small_function,
                                           std::vector<double>& angle_degree);
bool find_one_ls_segment_on_triangle(const double value, const Eigen::MatrixXi &F, const Eigen::MatrixXd &V,
                                      const Eigen::VectorXd &fvalues, const int fid, Eigen::Vector3d &E0, Eigen::Vector3d &E1);
void extract_levelset_web(const CGMesh &lsmesh, const Eigen::MatrixXd &V,
                          const Eigen::MatrixXi &F, const Eigen::VectorXd &ls0, const Eigen::VectorXd &ls1,
                          const int nbr_ls0, const int nbr_ls1, const int threadshold_nbr,
                          Eigen::MatrixXd &vers, Eigen::MatrixXi &Faces, bool even_pace);
void extract_Origami_web(const CGMesh &lsmesh, const Eigen::MatrixXd &V,const std::vector<CGMesh::HalfedgeHandle>& loop,
                         const Eigen::MatrixXi &F, const Eigen::VectorXd &ls,
                         const int expect_nbr_ls, const int expect_nbr_dis,
                         Eigen::MatrixXd &vers, Eigen::MatrixXi &Faces);
double get_mat_max_diag(spMat& M);
void solve_mean_value_laplacian_mat(CGMesh& lsmesh, const std::vector<int>& IVids, spMat& mat);
spMat sum_uneven_spMats(const spMat& mat_small, const spMat& mat_large);
Eigen::VectorXd sum_uneven_vectors(const Eigen::VectorXd& vsmall, const Eigen::VectorXd& vlarge);
spMat three_spmat_in_diag(const spMat& mat0, const spMat& mat1, const spMat& mat2, const int ntarget);
Eigen::VectorXd three_vec_in_row(const Eigen::VectorXd& ve0, const Eigen::VectorXd& ve1, const Eigen::VectorXd& ve2, const int ntarget);
// Ft, Vt belong to the triangle mesh,  
bool write_quad_mesh_with_binormal(const std::string & fname, const Eigen::MatrixXd& Vt, const Eigen::MatrixXi& Ft, const Eigen::MatrixXd& bi,
const Eigen::MatrixXd& Vq, const Eigen::MatrixXi& Fq);
void levelset_unit_scale(Eigen::VectorXd& func, Eigen::MatrixXd &GradValueF, const double length);
void extract_Quad_Mesh_Zigzag(const CGMesh &lsmesh,const std::vector<CGMesh::HalfedgeHandle>& loop, const Eigen::MatrixXd &V,
                         const Eigen::MatrixXi &F, const Eigen::VectorXd &ls,
                         const int expect_nbr_ls, const int expect_nbr_dis, const int threadshold_nbr,
                         Eigen::MatrixXd &vers, Eigen::MatrixXi &Faces);
Eigen::MatrixXd get_each_face_direction(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXd &func);