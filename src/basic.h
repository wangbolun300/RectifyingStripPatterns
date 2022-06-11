#pragma once
#include<lsc/MeshProcessing.h>
#include <igl/lscm.h>
#include <igl/harmonic.h>
#include <igl/boundary_loop.h>
class lsTools{
public:
// Eigen::SparseVector<double> a;
// void dddd(){
//     a.resize();

// }
void mesh_parameterization(
		const std::string &meshfile, Eigen::MatrixXd &V, Eigen::MatrixXd &param, Eigen::MatrixXi &F)
	{
		
		// /////////////////////////////////////////
		// // parameterization part
		// Eigen::VectorXi bnd;
		// igl::boundary_loop(F, bnd); // boundary vertices detection
		// Eigen::MatrixXd bnd_uv;
		// /*igl::map_vertices_to_circle(V, bnd, bnd_uv);*/
		// map_vertices_to_square(V, bnd, bnd_uv);

		// igl::harmonic(V, F, bnd, bnd_uv, 1, param);
	}



};