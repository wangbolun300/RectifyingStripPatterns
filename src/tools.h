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

bool save_levelset(const Eigen::VectorXd &ls);
bool read_levelset(Eigen::VectorXd &ls);
std::array<int,4> get_vers_around_edge(CGMesh& lsmesh, int edgeid, int& fid1, int &fid2);
double get_t_of_segment(const Eigen::Vector3d &ver, const Eigen::Vector3d &start, const Eigen::Vector3d &end);
double get_t_of_value(const double &ver, const double &start, const double &end);
Eigen::Vector2d get_2d_ver_from_t(const double t, const Eigen::Vector2d& start, const Eigen::Vector2d& end);
Eigen::Vector3d get_3d_ver_from_t(const double t, const Eigen::Vector3d& start, const Eigen::Vector3d& end);
bool angles_match(const double angle_degree1, const double angle_degree2);