#pragma once

#include <Eigen/Core>
#include <Eigen/Sparse>
// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>

#include <ANN/ANN.h>

// ----------------------------------------------------------------------------
typedef OpenMesh::PolyMesh_ArrayKernelT<>  CGMesh;
typedef Eigen::SparseMatrix<double> spMat;
typedef Eigen::Triplet<double> T;


class MeshProcessing
{
public:
	MeshProcessing() {};
	~MeshProcessing();



private:


public:
	//std::vector<double> fillMatandSolve_CheckBoard();
	void mesh2Matrix(CGMesh& mesh, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
	void matrix2Mesh(CGMesh& mesh, Eigen::MatrixXd& V, Eigen::MatrixXi& F);
	void meshEdges(CGMesh& mesh, Eigen::MatrixXi& E);
	void MeshUpdate(CGMesh& originMesh, CGMesh& updatedMesh, std::vector<double>& pts);
	void MeshUpdate(CGMesh& originMesh, CGMesh& updatedMesh, double* pts);
	std::vector<int> meshCorners(CGMesh& mesh);
	void MeshEdgeDual(CGMesh& coarseMesh, CGMesh& fineMesh, std::vector<int>& Vtype, double s, int type);

	spMat SubdivCoarseMesh_CaltmullClack(CGMesh& coarseMesh, CGMesh& fineMesh);
	std::vector<int> MeshVertexClassificationBlackWhite(CGMesh& mesh);
	std::vector<int> MeshFaceClassificationBlackWhite(CGMesh& mesh);
	std::vector<int> MeshFaceClassificationCube4types(CGMesh& mesh, int stFaceId, std::vector<int>& vType);
	void getEdgeLengthStatistics(CGMesh& mesh, double& totalL, double& minL, double& maxL, double& avgL);
	void MeshUnitScale(CGMesh& originMesh, CGMesh& updatedMesh);
	void MeshScale(CGMesh& originMesh, CGMesh& updatedMesh, double s);
	void MeshDiag(CGMesh& originMesh, CGMesh& updatedMesh, int type);


	//***************************************** mesh refinement *********************
	void subdivision_EdgeSplit(CGMesh& coarseMesh, CGMesh& fineMesh);
	
	//***************************************** mesh Planarity *********************
	void meshPlanarity(CGMesh& inputMesh, std::vector<double>& planarity);
	void meshFaceCenter(CGMesh& inputMesh, std::vector<OpenMesh::Vec3d>& fcen);



public:
	//CGMesh coarseMesh;




private:

};
