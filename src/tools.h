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

bool triangles_coplanar(const Eigen::Vector3d& t0, const Eigen::Vector3d& t1, const Eigen::Vector3d& t2,
 const Eigen::Vector3d&p0, const Eigen::Vector3d&p1, const Eigen::Vector3d&p2);
 // -1: boundary edge
 // 0 not coplanar
 // 1 co-planar
int edge_is_coplanar(const CGMesh& lsmesh, const CGMesh::EdgeHandle& edge_middle, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);
std::vector<double> polynomial_simplify(const std::vector<double> &poly);
std::vector<double> polynomial_add(const std::vector<double> &poly1, const std::vector<double> &poly2);
std::vector<double> polynomial_times(const std::vector<double> &poly1, const std::vector<double> &poly2);
std::vector<double> polynomial_times(const std::vector<double> &poly1, const double &nbr);
double polynomial_value(const std::vector<double> &poly, const double para);
std::vector<double> polynomial_integration(const std::vector<double> &poly);
double polynomial_integration(const std::vector<double> &poly, const double lower, const double upper);

void cylinder_example(double radius, double height, int nr, int nh);
// theta <= pi/2, phi<=pi
void sphere_example(double radius, double theta, double phi, int nt, int np);