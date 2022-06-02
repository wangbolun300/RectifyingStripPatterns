#ifndef OPTIMIZER_CHECKBOARD_H
#define OPTIMIZER_CHECKBOARD_H

//#include <Eigen/CholmodSupport>
//#include <Eigen/UmfPackSupport>
//#include <Eigen/SPQRSupport>
#include <Eigen/Sparse>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/IterativeLinearSolvers>
#include <Eigen/SparseQR>
// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#include <ANN/ANN.h>

// ----------------------------------------------------------------------------
typedef OpenMesh::PolyMesh_ArrayKernelT<>  CGMesh;
typedef Eigen::SparseMatrix<double> spMat;
typedef Eigen::Triplet<double> T;

struct spTriple {
	int row;
	int col;
	double value;
};

class OptimizerCheckboard
{
public:
	OptimizerCheckboard() {};
	~OptimizerCheckboard();



private:
	void finalize();
	void updateNormalization(double* x, int startId, int num);

	double LinearEquation_closeness(double* x, double weight, int StId, bool FirstBuild);
	double LinearEquation_closenessReferPt(double* x, double weight, int StId, bool FirstBuild);
	double LinearEquation_polyhedron(double* x, double weight, int StId, bool FirstBuild);
	double LinearEquation_polyhedron_vertex(double* x, double weight, int StIdVn, int StIdPt, bool FirstBuild);
	double LinearEquation_Normalize_vnormal(double* x, double weight, int StId, bool FirstBuild);
	double LinearEquation_EqualDiagornals(double* x, double weight, int StId, bool FirstBuild);
	double LinearEquation_EqualDiagornals2D3D(double* x, double weight, int StId, bool FirstBuild);
	double LinearEquation_ConstDiagornals(double* x, double weight, int StId, bool FirstBuild);
	double LinearEquation_OrthDiagornals(double* x, double weight, int StId, bool FirstBuild);
	double LinearEquation_UniformDiagornals(double* x, double weight, int StId, bool FirstBuild);
	double LinearEquation_Fairness(double* x, double weight, int StId, bool FirstBuild);
	double LinearEquation_FairnessDiag(double* x, double weight, int StId, bool FirstBuild);


	


	double getAvgEdgeLength(double* x);



	void GetInitialValue(std::vector<double>& values);
	void GetInitialTriplet(double* x, bool firstBuild, bool PrintFlag, int Iter);



	double get_closest_point(const CGMesh::Point& querypoint, CGMesh::Point& closestpoint, CGMesh::Point& closestpointnormal);
	double get_closest_point_ID(const CGMesh::Point& querypoint, int& id);
	
	double closest_point_triangle(const CGMesh::Point& querypoint, std::vector<CGMesh::Point>& Tripts, CGMesh::Point& closestpoint);
	double closest_point_lineSegement(const CGMesh::Point& querypoint, CGMesh::Point& p1, CGMesh::Point& p2, CGMesh::Point& fn, CGMesh::Point& closestpoint);




	
	//
	void getDiagInfo(CGMesh& mesh, std::vector<OpenMesh::Vec3d>& Diag3);

	




public:
	std::vector<double> fillMatandSolve_CheckBoard();
	double closest_point_projection(const CGMesh::Point& querypoint, CGMesh::Point& closestpoint);


public:
	CGMesh coarseMesh;
	CGMesh referenceMesh;
	CGMesh CloseRefMesh;
	std::vector<OpenMesh::Vec3i> Diag3pts; // 3 pts along diag mesh
	int OpIter;
	// weight for different energy terms. 
	double weight_closeness;
	double weight_polyhedron;
	double weight_EqualDiagornals;
	double weight_OrthDiagornals;
	double weight_UniformDiagornals;
	double weight_Fairness;
	double weight_FairnessDiag;
	double weight_polyhedron_vertex;
	double weight_Nornalization_vertexNormal;

	double refAveEL; // reference average edge length;
	int fairnessType;
	bool relativeFairness;
	int fairnessChangeType;
	double weight_FairnessWeightScale;
	double weight_ForceBalanceDiagMesh;

	double weight_FaceNormal_QuadCurvature;
	double weight_SelfAdjoint_QuadCurvature;

	double weight_ElasticalEnergy;


	double aveEL; // average Edge length;
	double Lsqure;
	double shellthickness;
	double PoissonRatio;
	bool BendingFlag;
	bool StrechingFlag;
	double DampingPara;
	bool TripletFlag;


	// 
	std::vector<int> FixedVertices;






private:
	std::vector<T> tripletList;
	std::vector<spTriple> spTripletList;
	std::vector<int> tripleRow;
	std::vector<int> tripleCol;
	std::vector<double> tripleValue;
	std::vector<int> FixedVerticesFlag;
	int EqID;   // equation id;
	int TID;   // triplet id;
	std::vector<double> Rhs;

	std::vector<OpenMesh::Vec5i> V4DiagE;
	std::vector<OpenMesh::Vec5i> V4DiagV;
	std::vector<double> vLoad;

	std::vector<OpenMesh::Vec3d> DiagInfos;

	//! kd-tree
	ANNkd_tree* _kdTree;
	//! points in the kd-tree
	ANNpointArray _dataPts;
	//! normals of points in the kd-tree
	CGMesh::Point* _pointnormals;


};
#endif