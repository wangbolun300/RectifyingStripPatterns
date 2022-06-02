#include "OptimizerCheckboard.h"

OptimizerCheckboard::~OptimizerCheckboard()
{
	finalize();
}

void OptimizerCheckboard::finalize()
{
	//if (_pointnormals)
	//	delete[] _pointnormals;
	//if (_kdTree)
	//	delete _kdTree;
	//if (_dataPts)
	//	annDeallocPts(_dataPts);

}

void OptimizerCheckboard::updateNormalization(double* x, int startId, int num)
{
	for (int i = 0; i < num; i++)
	{
		OpenMesh::Vec3d p;
		p = OpenMesh::Vec3d(x[startId + 3 * i], x[startId + 3 * i + 1], x[startId + 3 * i + 2]);
		p.normalize();
		x[startId + 3 * i] = p[0];
		x[startId + 3 * i + 1] = p[1];
		x[startId + 3 * i + 2] = p[2];
	}
}

double OptimizerCheckboard::LinearEquation_closeness(double* x, double weight, int StId, bool FirstBuild)
{
	int nv = coarseMesh.n_vertices();
	double f = 0;

	for (CGMesh::VertexIter v_it = coarseMesh.vertices_begin(); v_it != (coarseMesh.vertices_end()); ++v_it)
	{
		CGMesh::Point p, cpt, cfn;
		int vid = v_it.handle().idx();
		int m = 3 * vid + StId;
		p = OpenMesh::Vec3d(x[m], x[m + 1], x[m + 2]);
		get_closest_point(p, cpt, cfn);

		double residual = (p - cpt) | cfn;
		f += weight * weight * residual * residual;
		for (int j = 0; j < 3; j++)
		{
			if (FirstBuild)
				spTripletList.push_back(spTriple({ EqID, m + j, cfn[j] * weight }));

			else
			{
				spTripletList[TID].value = cfn[j] * weight; TID++;
			}
		}

		Rhs.push_back(-residual * weight);
		EqID++;
	}
	return f;
}

double OptimizerCheckboard::LinearEquation_closenessReferPt(double* x, double weight, int StId, bool FirstBuild)
{

	int nv = coarseMesh.n_vertices();
	double f = 0;

	for (CGMesh::VertexIter v_it = coarseMesh.vertices_begin(); v_it != (coarseMesh.vertices_end()); ++v_it)
	{
		OpenMesh::Vec3d p, cpt, cfn;
		int vid = v_it.handle().idx();
		int m = 3 * vid + StId;
		p = OpenMesh::Vec3d(x[m], x[m + 1], x[m + 2]);
		cpt = p;
		cfn = coarseMesh.normal(v_it);

		OpenMesh::Vec3d term = p - cpt;

		for (int j = 0; j < 3; j++)
		{
			double residual = term[j];
			if (FirstBuild)
				spTripletList.push_back(spTriple({ EqID, m + j, 1 * weight }));

			else
			{
				spTripletList[TID].value = 1 * weight; TID++;
			}
			f += weight * weight * residual * residual;
			Rhs.push_back(-residual * weight);
			EqID++;
		}
	}
	return f;

}

double OptimizerCheckboard::LinearEquation_polyhedron(double* x, double weight, int StId, bool FirstBuild)
{
	int nv = coarseMesh.n_vertices();
	double f = 0;

	for (CGMesh::FaceIter f_it = coarseMesh.faces_begin(); f_it != (coarseMesh.faces_end()); ++f_it)
	{
		std::vector <OpenMesh::VertexHandle> vecfvh;
		for (CGMesh::FaceVertexIter fv_it = coarseMesh.fv_begin(f_it.handle()); fv_it != coarseMesh.fv_end(f_it.handle()); ++fv_it)
			vecfvh.push_back(fv_it.handle());

		OpenMesh::Vec3d fn;
		int fid = f_it.handle().idx();
		int m = 3 * nv + 3 * fid + StId;
		fn = OpenMesh::Vec3d(x[m], x[m + 1], x[m + 2]);
		int numvalence = vecfvh.size();
		for (int i = 0; i < numvalence; i++)
		{
			int i1 = i;
			int i2 = (i + 1) % numvalence;

			OpenMesh::VertexHandle vh1, vh2;
			OpenMesh::Vec3d p1, p2;
			vh1 = vecfvh[i1]; vh2 = vecfvh[i2];
			int vid1, vid2;
			vid1 = vh1.idx();
			vid2 = vh2.idx();
			int m1, m2;
			m1 = 3 * vid1 + StId;
			m2 = 3 * vid2 + StId;

			p1 = OpenMesh::Vec3d(x[m1], x[m1 + 1], x[m1 + 2]);
			p2 = OpenMesh::Vec3d(x[m2], x[m2 + 1], x[m2 + 2]);

			double residual;
			residual = (p1 - p2) | fn;
			f += residual * residual * weight * weight;

			for (int j = 0; j < 3; j++)
			{
				if (FirstBuild)
				{
					spTripletList.push_back(spTriple({ EqID, m1 + j, fn[j] * weight }));
					spTripletList.push_back(spTriple({ EqID, m2 + j, -fn[j] * weight }));
					spTripletList.push_back(spTriple({ EqID, m + j, (p1 - p2)[j] * weight }));
				}
				else
				{
					spTripletList[TID].value = fn[j] * weight; TID++;
					spTripletList[TID].value = -fn[j] * weight; TID++;
					spTripletList[TID].value = (p1 - p2)[j] * weight; TID++;
				}
			}

			Rhs.push_back(-residual * weight);
			EqID++;
		}
	}
	return f;
}

double OptimizerCheckboard::LinearEquation_polyhedron_vertex(double* x, double weight, int StIdVn, int StIdPt, bool FirstBuild)
{
	int nv = coarseMesh.n_vertices();
	double f = 0;

	for (CGMesh::VertexIter v_it = coarseMesh.vertices_begin(); v_it != (coarseMesh.vertices_end()); ++v_it)
	{
		int nvv = coarseMesh.valence(v_it);
		if (coarseMesh.is_boundary(v_it) || nvv != 4)
			continue;

		std::vector <OpenMesh::VertexHandle> vecvvh;
		for (CGMesh::VertexVertexIter vv_it = coarseMesh.vv_begin(v_it.handle()); vv_it != coarseMesh.vv_end(v_it.handle()); ++vv_it)
			vecvvh.push_back(vv_it.handle());

		OpenMesh::Vec3d vn;
		int vid = v_it.handle().idx();
		int m = 3 * vid + StIdVn;
		vn = OpenMesh::Vec3d(x[m], x[m + 1], x[m + 2]);
		int numvalence = vecvvh.size();

		for (int i = 0; i < numvalence; i++)
		{
			int i1 = i;
			int i2 = (i + 1) % numvalence;

			OpenMesh::VertexHandle vh1, vh2;
			OpenMesh::Vec3d p1, p2;
			vh1 = vecvvh[i1]; vh2 = vecvvh[i2];
			int vid1, vid2;
			vid1 = vh1.idx();
			vid2 = vh2.idx();
			int m1, m2;
			m1 = 3 * vid1 + StIdPt;
			m2 = 3 * vid2 + StIdPt;

			p1 = OpenMesh::Vec3d(x[m1], x[m1 + 1], x[m1 + 2]);
			p2 = OpenMesh::Vec3d(x[m2], x[m2 + 1], x[m2 + 2]);

			double residual;
			residual = (p1 - p2) | vn;
			f += residual * residual * weight * weight;

			for (int j = 0; j < 3; j++)
			{
				if (FirstBuild)
				{
					spTripletList.push_back(spTriple({ EqID, m1 + j, vn[j] * weight }));
					spTripletList.push_back(spTriple({ EqID, m2 + j, -vn[j] * weight }));
					spTripletList.push_back(spTriple({ EqID, m + j, (p1 - p2)[j] * weight }));
				}
				else
				{
					spTripletList[TID].value = vn[j] * weight; TID++;
					spTripletList[TID].value = -vn[j] * weight; TID++;
					spTripletList[TID].value = (p1 - p2)[j] * weight; TID++;
				}
			}

			Rhs.push_back(-residual * weight);
			EqID++;
		}
	}
	return f;
}

double OptimizerCheckboard::LinearEquation_Normalize_vnormal(double* x, double weight, int StId, bool FirstBuild)
{
	int nv = coarseMesh.n_vertices();
	double f = 0;

	for (CGMesh::VertexIter v_it = coarseMesh.vertices_begin(); v_it != (coarseMesh.vertices_end()); ++v_it)
	{
		OpenMesh::Vec3d vn;
		int vid = v_it.handle().idx();
		int m = 3 * vid + StId;
		vn = OpenMesh::Vec3d(x[m], x[m + 1], x[m + 2]);
		double residual = (vn | vn) - 1;
		f += residual * residual * weight * weight;
		for (int j = 0; j < 3; j++)
		{
			if (FirstBuild)
			{
				spTripletList.push_back(spTriple({ EqID, m + j, 2 * vn[j] * weight }));
			}
			else
			{
				spTripletList[TID].value = 2 * vn[j] * weight; TID++;
			}
		}
		Rhs.push_back(-residual * weight);
		EqID++;
	}
	return f;
}

double OptimizerCheckboard::LinearEquation_EqualDiagornals(double* x, double weight, int StId, bool FirstBuild)
{
	double f = 0;

	for (CGMesh::FaceIter f_it = coarseMesh.faces_begin(); f_it != (coarseMesh.faces_end()); ++f_it)
	{
		std::vector <OpenMesh::VertexHandle> vecfvh;
		for (CGMesh::FaceVertexIter fv_it = coarseMesh.fv_begin(f_it.handle()); fv_it != coarseMesh.fv_end(f_it.handle()); ++fv_it)
			vecfvh.push_back(fv_it.handle());

		int m0, m1, m2, m3;
		OpenMesh::Vec3d p0, p1, p2, p3, p02, p13;
		if (vecfvh.size() == 4)
		{
			m0 = vecfvh[0].idx() * 3 + StId; m1 = vecfvh[1].idx() * 3 + StId; m2 = vecfvh[2].idx() * 3 + StId; m3 = vecfvh[3].idx() * 3 + StId;
			p0 = OpenMesh::Vec3d(x[m0], x[m0 + 1], x[m0 + 2]);
			p1 = OpenMesh::Vec3d(x[m1], x[m1 + 1], x[m1 + 2]);
			p2 = OpenMesh::Vec3d(x[m2], x[m2 + 1], x[m2 + 2]);
			p3 = OpenMesh::Vec3d(x[m3], x[m3 + 1], x[m3 + 2]);
			p02 = p0 - p2;
			p13 = p1 - p3;
			double residual = (p02 | p02) - (p13 | p13);
			f += residual * residual * weight * weight;

			for (int j = 0; j < 3; j++)
			{
				if (FirstBuild)
				{
					spTripletList.push_back(spTriple({ EqID, m0 + j, 2 * p02[j] * weight }));
					spTripletList.push_back(spTriple({ EqID, m1 + j, -2 * p13[j] * weight }));
					spTripletList.push_back(spTriple({ EqID, m2 + j, -2 * p02[j] * weight }));
					spTripletList.push_back(spTriple({ EqID, m3 + j, 2 * p13[j] * weight }));
				}
				else
				{
					spTripletList[TID].value = 2 * p02[j] * weight; TID++;
					spTripletList[TID].value = -2 * p13[j] * weight; TID++;
					spTripletList[TID].value = -2 * p02[j] * weight; TID++;
					spTripletList[TID].value = 2 * p13[j] * weight; TID++;
				}
			}
			Rhs.push_back(-residual * weight);
			EqID++;
		}
	}
	return f;
}

double OptimizerCheckboard::LinearEquation_EqualDiagornals2D3D(double* x, double weight, int StId, bool FirstBuild)
{

	double f = 0;
	if (DiagInfos.size() != coarseMesh.n_faces())
		return f;

	for (CGMesh::FaceIter f_it = coarseMesh.faces_begin(); f_it != (coarseMesh.faces_end()); ++f_it)
	{
		int fid = f_it.handle().idx();
		std::vector <OpenMesh::VertexHandle> vecfvh;
		for (CGMesh::FaceVertexIter fv_it = coarseMesh.fv_begin(f_it.handle()); fv_it != coarseMesh.fv_end(f_it.handle()); ++fv_it)
			vecfvh.push_back(fv_it.handle());

		int m0, m1, m2, m3;
		OpenMesh::Vec3d p0, p1, p2, p3, p02, p13;
		if (vecfvh.size() == 4)
		{
			OpenMesh::Vec3d Diag = DiagInfos[fid];
			m0 = vecfvh[0].idx() * 3 + StId; m1 = vecfvh[1].idx() * 3 + StId; m2 = vecfvh[2].idx() * 3 + StId; m3 = vecfvh[3].idx() * 3 + StId;
			p0 = OpenMesh::Vec3d(x[m0], x[m0 + 1], x[m0 + 2]);
			p1 = OpenMesh::Vec3d(x[m1], x[m1 + 1], x[m1 + 2]);
			p2 = OpenMesh::Vec3d(x[m2], x[m2 + 1], x[m2 + 2]);
			p3 = OpenMesh::Vec3d(x[m3], x[m3 + 1], x[m3 + 2]);
			p02 = p0 - p2;
			p13 = p1 - p3;

			//***** d0 ******************
			double residual = (p02 | p02) - Diag[0];
			f += residual * residual * weight * weight;
			for (int j = 0; j < 3; j++)
			{
				if (FirstBuild)
				{
					spTripletList.push_back(spTriple({ EqID, m0 + j, 2 * p02[j] * weight }));
					spTripletList.push_back(spTriple({ EqID, m2 + j, -2 * p02[j] * weight }));
				}
				else
				{
					spTripletList[TID].value = 2 * p02[j] * weight; TID++;
					spTripletList[TID].value = -2 * p02[j] * weight; TID++;
				}
			}
			Rhs.push_back(-residual * weight);
			EqID++;

			//***** d1 ****************
			residual = (p13 | p13) - Diag[1];
			f += residual * residual * weight * weight;
			for (int j = 0; j < 3; j++)
			{
				if (FirstBuild)
				{
					spTripletList.push_back(spTriple({ EqID, m1 + j, 2 * p13[j] * weight }));
					spTripletList.push_back(spTriple({ EqID, m3 + j, -2 * p13[j] * weight }));
				}
				else
				{
					spTripletList[TID].value = 2 * p13[j] * weight; TID++;
					spTripletList[TID].value = -2 * p13[j] * weight; TID++;
				}
			}
			Rhs.push_back(-residual * weight);
			EqID++;

			//***** d2 ******************


			residual = (p02 | p13) - Diag[2];
			f += residual * residual * weight * weight;

			for (int j = 0; j < 3; j++)
			{
				if (FirstBuild)
				{
					spTripletList.push_back(spTriple({ EqID, m0 + j,  p13[j] * weight }));
					spTripletList.push_back(spTriple({ EqID, m1 + j,  p02[j] * weight }));
					spTripletList.push_back(spTriple({ EqID, m2 + j, -1 * p13[j] * weight }));
					spTripletList.push_back(spTriple({ EqID, m3 + j,  -1 * p02[j] * weight }));
				}
				else
				{
					spTripletList[TID].value = p13[j] * weight; TID++;
					spTripletList[TID].value = p02[j] * weight; TID++;
					spTripletList[TID].value = -1 * p13[j] * weight; TID++;
					spTripletList[TID].value = -1 * p02[j] * weight; TID++;
				}
			}
			Rhs.push_back(-residual * weight);
			EqID++;
		}
	}
	return f;
}

double OptimizerCheckboard::LinearEquation_ConstDiagornals(double* x, double weight, int StId, bool FirstBuild)
{
	double f = 0;


	for (CGMesh::FaceIter f_it = coarseMesh.faces_begin(); f_it != (coarseMesh.faces_end()); ++f_it)
	{
		std::vector <OpenMesh::VertexHandle> vecfvh;
		for (CGMesh::FaceVertexIter fv_it = coarseMesh.fv_begin(f_it.handle()); fv_it != coarseMesh.fv_end(f_it.handle()); ++fv_it)
			vecfvh.push_back(fv_it.handle());

		int m0, m1, m2, m3;
		OpenMesh::Vec3d p0, p1, p2, p3, p02, p13;
		if (vecfvh.size() == 4)
		{
			for (int num = 0; num < 2; num++)
			{
				m0 = vecfvh[num].idx() * 3 + StId;  m2 = vecfvh[num + 2].idx() * 3 + StId;
				p0 = OpenMesh::Vec3d(x[m0], x[m0 + 1], x[m0 + 2]);
				p2 = OpenMesh::Vec3d(x[m2], x[m2 + 1], x[m2 + 2]);
				p02 = p0 - p2;
				double residual = (p02 | p02) - Lsqure;
				f += residual * residual * weight * weight;

				for (int j = 0; j < 3; j++)
				{
					if (FirstBuild)
					{
						spTripletList.push_back(spTriple({ EqID, m0 + j, 2 * p02[j] * weight }));
						spTripletList.push_back(spTriple({ EqID, m2 + j, -2 * p02[j] * weight }));
					}
					else
					{
						spTripletList[TID].value = 2 * p02[j] * weight; TID++;
						spTripletList[TID].value = -2 * p02[j] * weight; TID++;
					}
				}
				Rhs.push_back(-residual * weight);
				EqID++;
			}
		}
	}
	return f;
}

double OptimizerCheckboard::LinearEquation_OrthDiagornals(double* x, double weight, int StId, bool FirstBuild)
{
	double f = 0;

	for (CGMesh::FaceIter f_it = coarseMesh.faces_begin(); f_it != (coarseMesh.faces_end()); ++f_it)
	{
		std::vector <OpenMesh::VertexHandle> vecfvh;
		for (CGMesh::FaceVertexIter fv_it = coarseMesh.fv_begin(f_it.handle()); fv_it != coarseMesh.fv_end(f_it.handle()); ++fv_it)
			vecfvh.push_back(fv_it.handle());

		int m0, m1, m2, m3;
		OpenMesh::Vec3d p0, p1, p2, p3;
		if (vecfvh.size() == 4)
		{
			m0 = vecfvh[0].idx() * 3 + StId; m1 = vecfvh[1].idx() * 3 + StId; m2 = vecfvh[2].idx() * 3 + StId; m3 = vecfvh[3].idx() * 3 + StId;
			p0 = OpenMesh::Vec3d(x[m0], x[m0 + 1], x[m0 + 2]);
			p1 = OpenMesh::Vec3d(x[m1], x[m1 + 1], x[m1 + 2]);
			p2 = OpenMesh::Vec3d(x[m2], x[m2 + 1], x[m2 + 2]);
			p3 = OpenMesh::Vec3d(x[m3], x[m3 + 1], x[m3 + 2]);
			double residual = (p0 - p2) | (p1 - p3);
			f += residual * residual * weight * weight;

			for (int j = 0; j < 3; j++)
			{
				if (FirstBuild)
				{
					spTripletList.push_back(spTriple({ EqID, m0 + j,  (p1[j] - p3[j]) * weight }));
					spTripletList.push_back(spTriple({ EqID, m1 + j,  (p0[j] - p2[j]) * weight }));
					spTripletList.push_back(spTriple({ EqID, m2 + j, -(p1[j] - p3[j]) * weight }));
					spTripletList.push_back(spTriple({ EqID, m3 + j, -(p0[j] - p2[j]) * weight }));
				}
				else
				{
					spTripletList[TID].value = (p1[j] - p3[j]) * weight; TID++;
					spTripletList[TID].value = (p0[j] - p2[j]) * weight; TID++;
					spTripletList[TID].value = -(p1[j] - p3[j]) * weight; TID++;
					spTripletList[TID].value = -(p0[j] - p2[j]) * weight; TID++;
				}
			}
			Rhs.push_back(-residual * weight);
			EqID++;
		}
	}
	return f;
}

double OptimizerCheckboard::LinearEquation_UniformDiagornals(double* x, double weight, int StId, bool FirstBuild)
{
	CGMesh::FaceIter   f_it;
	CGMesh::FaceVertexIter fv_it;

	double f = 0;
	for (CGMesh::EdgeIter e_it = coarseMesh.edges_begin(); e_it != coarseMesh.edges_end(); ++e_it)
	{
		if (coarseMesh.is_boundary(e_it))
			continue;
		CGMesh::HalfedgeHandle hh1, hh2;
		/*	coarseMesh->mesh()->getEdgeFaceHandles(e_it, fh1, fh2);
			hh1 = coarseMesh->mesh()->getHalfedge(e_it, fh1);
			hh2 = coarseMesh->mesh()->getHalfedge(e_it, fh2);*/
		hh1 = coarseMesh.halfedge_handle(e_it, 0);
		hh2 = coarseMesh.halfedge_handle(e_it, 1);
		std::vector<OpenMesh::Vec3i> ids;
		int id0, id1, id2, id3, id4, id5;
		id0 = coarseMesh.from_vertex_handle(hh1).idx();
		id1 = coarseMesh.to_vertex_handle(hh1).idx();
		id2 = coarseMesh.to_vertex_handle(coarseMesh.next_halfedge_handle(hh1)).idx();
		id3 = coarseMesh.from_vertex_handle(coarseMesh.prev_halfedge_handle(hh1)).idx();
		id4 = coarseMesh.to_vertex_handle(coarseMesh.next_halfedge_handle(hh2)).idx();
		id5 = coarseMesh.from_vertex_handle(coarseMesh.prev_halfedge_handle(hh2)).idx();
		ids.push_back(OpenMesh::Vec3i(id0, id2, id5));
		ids.push_back(OpenMesh::Vec3i(id1, id3, id4));
		int m0, m1, m2;
		for (int i = 0; i < 2; i++)
		{
			m0 = ids[i][0] * 3 + StId;; m1 = ids[i][1] * 3 + StId; m2 = ids[i][2] * 3 + StId;;

			OpenMesh::Vec3d p0, p1, p2, p01, p02;


			p0 = OpenMesh::Vec3d(x[m0], x[m0 + 1], x[m0 + 2]);
			p1 = OpenMesh::Vec3d(x[m1], x[m1 + 1], x[m1 + 2]);
			p2 = OpenMesh::Vec3d(x[m2], x[m2 + 1], x[m2 + 2]);

			p01 = p0 - p1;
			p02 = p0 - p2;
			double residual = (p01 | p01) - (p02 | p02);
			f += residual * residual * weight * weight;

			for (int j = 0; j < 3; j++)
			{
				if (FirstBuild)
				{
					spTripletList.push_back(spTriple({ EqID, m0 + j, 2 * (p2[j] - p1[j]) * weight }));
					spTripletList.push_back(spTriple({ EqID, m1 + j, -2 * p01[j] * weight }));
					spTripletList.push_back(spTriple({ EqID, m2 + j, 2 * p02[j] * weight }));
				}
				else
				{
					spTripletList[TID].value = 2 * (p2[j] - p1[j]) * weight; TID++;
					spTripletList[TID].value = -2 * p01[j] * weight; TID++;
					spTripletList[TID].value = 2 * p02[j] * weight; TID++;
				}
			}


			Rhs.push_back(-residual * weight);
			EqID++;
		}
	}
	return f;
}

double OptimizerCheckboard::LinearEquation_Fairness(double* x, double weight, int StId, bool FirstBuild)
{

	CGMesh::VertexIter   v_it;
	CGMesh::VertexVertexIter vv_it;

	double f = 0;
	for (v_it = coarseMesh.vertices_begin(); v_it != (coarseMesh.vertices_end()); ++v_it)
	{

		if (!coarseMesh.is_boundary(v_it.handle()))
		//if (coarseMesh.valence(v_it.handle())==4)
		{
			std::vector<OpenMesh::VertexHandle> vvh;
			std::vector<int> vvindex;
			std::vector<double> diff;
			OpenMesh::Vec3d p0, p1, p2, p3, p4, vn;
			OpenMesh::VertexHandle vh;

			for (vv_it = coarseMesh.vv_begin(v_it.handle()); vv_it != coarseMesh.vv_end(v_it.handle()); ++vv_it)
			{
				vvh.push_back(vv_it);
				vvindex.push_back(vv_it.handle().idx());
			}
			int vnnum = v_it.handle().idx();
			int m0 = 3 * vnnum + StId;
			p0 = OpenMesh::Vec3d(x[m0], x[m0 + 1], x[m0 + 2]);

			int num = vvh.size();

			if (num % 2 == 0)
			{
				int pair = num / 2;
				for (int i = 0; i < pair; i++)
				{
					int m1 = vvindex[i] * 3 + StId;
					int m2 = vvindex[i + pair] * 3 + StId;


					p1 = OpenMesh::Vec3d(x[m1], x[m1 + 1], x[m1 + 2]);
					p2 = OpenMesh::Vec3d(x[m2], x[m2 + 1], x[m2 + 2]);
					OpenMesh::Vec3d term;
					term = (p1 + p2 - p0 * 2);

					if (relativeFairness)
						term = term - OpenMesh::Vec3d(referenceMesh.point(vvh[i]) + referenceMesh.point(vvh[i + pair]) - 2 * referenceMesh.point(v_it));



					double residual;

					for (int j = 0; j < 3; j++)
					{
						residual = term[j];
						/*if (relativeFairness)
							residual = 0;*/

						f += residual * residual * weight * weight;

						if (FirstBuild)
						{
							spTripletList.push_back(spTriple({ EqID, m1 + j, weight }));
							spTripletList.push_back(spTriple({ EqID, m2 + j, weight }));
							spTripletList.push_back(spTriple({ EqID, m0 + j, -2 * weight }));
						}
						else
						{
							spTripletList[TID].value = weight; TID++;
							spTripletList[TID].value = weight; TID++;
							spTripletList[TID].value = -2 * weight; TID++;
						}
						Rhs.push_back(-residual * weight);
						EqID++;
					}
				}
			}
			else
			{

				OpenMesh::Vec3d sumpt = OpenMesh::Vec3d(0.0, 0.0, 0.0);
				for (int i = 0; i < num; i++)
				{
					int mi = vvindex[i] * 3 + StId;
					OpenMesh::Vec3d pt = OpenMesh::Vec3d(x[mi], x[mi + 1], x[mi + 2]);
					sumpt = sumpt + pt;
				}
				OpenMesh::Vec3d term = sumpt - p0 * num;

				double residual;
				//double weight=mf_Fairness_Weight*scale;
				for (int j = 0; j < 3; j++)
				{
					residual = term[j];
					f += residual * residual * weight * weight;

					for (int i = 0; i < num; i++)
					{
						int mi = vvindex[i] * 3 + StId;
						if (FirstBuild)
							spTripletList.push_back(spTriple({ EqID, mi + j, weight }));
						else
							spTripletList[TID].value = weight; TID++;
					}
					if (FirstBuild)
						spTripletList.push_back(spTriple({ EqID, m0 + j, -num * weight }));
					else
						spTripletList[TID].value = -num * weight; TID++;
					
					Rhs.push_back(-residual * weight);
					EqID++;
				}
			}
		}
		else
		{
			if (coarseMesh.valence(v_it) > 2)
			{
				std::vector<int> vids;
				std::vector<OpenMesh::VertexHandle> vvh;
				for (vv_it = coarseMesh.vv_begin(v_it.handle()); vv_it != coarseMesh.vv_end(v_it.handle()); ++vv_it)
				{
					if (coarseMesh.is_boundary(vv_it))
					{
						vids.push_back(vv_it.handle().idx());
						vvh.push_back(vv_it.handle());
					}
				}
				if (vids.size() == 2)
				{
					int vnnum = v_it.handle().idx();
					int m0 = 3 * vnnum + StId;
					OpenMesh::Vec3d p0 = OpenMesh::Vec3d(x[m0], x[m0 + 1], x[m0 + 2]);
					int m1 = vids[0] * 3 + StId;
					int m2 = vids[1] * 3 + StId;
					OpenMesh::Vec3d p1 = OpenMesh::Vec3d(x[m1], x[m1 + 1], x[m1 + 2]);
					OpenMesh::Vec3d p2 = OpenMesh::Vec3d(x[m2], x[m2 + 1], x[m2 + 2]);
					OpenMesh::Vec3d term = (p1 + p2 - p0 * 2);
					if (relativeFairness)
						term = term - OpenMesh::Vec3d(referenceMesh.point(vvh[0]) + referenceMesh.point(vvh[1]) - 2 * referenceMesh.point(v_it));

					double residual;
					double s = 2;
					for (int j = 0; j < 3; j++)
					{
						residual = term[j];
						/*				if (relativeFairness)
											residual = 0;*/

						f += residual * residual * weight * weight * s * s;

						if (FirstBuild)
						{
							spTripletList.push_back(spTriple({ EqID, m1 + j, weight * s }));
							spTripletList.push_back(spTriple({ EqID, m2 + j, weight * s }));
							spTripletList.push_back(spTriple({ EqID, m0 + j, -2 * weight * s }));
						}
						else
						{
							spTripletList[TID].value = weight * s; TID++;
							spTripletList[TID].value = weight * s; TID++;
							spTripletList[TID].value = -2 * weight * s; TID++;
						}
						Rhs.push_back(-residual * weight * s);
						EqID++;
					}
				}
			}
		}
	}
	return f;
}


double OptimizerCheckboard::LinearEquation_FairnessDiag(double* x, double weight, int StId, bool FirstBuild)
{
	double f = 0;
	for (int i = 0; i < Diag3pts.size(); i++)
	{
		OpenMesh::Vec3i pts = Diag3pts[i];
		int m0 = 3 * pts[0] + StId;
		int m1 = 3 * pts[1] + StId;
		int m2 = 3 * pts[2] + StId;
		OpenMesh::Vec3d p0 = OpenMesh::Vec3d(x[m0], x[m0 + 1], x[m0 + 2]);
		OpenMesh::Vec3d p1 = OpenMesh::Vec3d(x[m1], x[m1 + 1], x[m1 + 2]);
		OpenMesh::Vec3d p2 = OpenMesh::Vec3d(x[m2], x[m2 + 1], x[m2 + 2]);
		OpenMesh::Vec3d term = (p0 + p2 - p1 * 2);
		double residual;

		for (int j = 0; j < 3; j++)
		{
			if (FirstBuild)
			{
				spTripletList.push_back(spTriple({ EqID, m0 + j, weight }));
				spTripletList.push_back(spTriple({ EqID, m2 + j, weight }));
				spTripletList.push_back(spTriple({ EqID, m1 + j, -2 * weight }));
			}
			else
			{
				spTripletList[TID].value = weight; TID++;
				spTripletList[TID].value = weight; TID++;
				spTripletList[TID].value = -2 * weight; TID++;
			}
			residual = term[j];
			f += residual * residual * weight * weight;
			Rhs.push_back(-residual * weight);
			EqID++;
		}
	}
	return f;
}

double OptimizerCheckboard::getAvgEdgeLength(double* x)
{
	int ne = coarseMesh.n_edges();
	double totalL = 0;
	double minL = std::numeric_limits<double>::max();
	double maxL = 0;
	for (CGMesh::EdgeIter e_it = coarseMesh.edges_begin(); e_it != coarseMesh.edges_end(); ++e_it)
	{
		OpenMesh::HalfedgeHandle hh = coarseMesh.halfedge_handle(e_it, 0);
		if (!hh.is_valid()) hh = coarseMesh.halfedge_handle(e_it, 1);
		int m1 = coarseMesh.from_vertex_handle(hh).idx();
		int m2 = coarseMesh.to_vertex_handle(hh).idx();
		OpenMesh::Vec3d p1 = OpenMesh::Vec3d(x[3 * m1], x[3 * m1 + 1], x[3 * m1 + 2]);
		OpenMesh::Vec3d p2 = OpenMesh::Vec3d(x[3 * m2], x[3 * m2 + 1], x[3 * m2 + 2]);

		double eL = (p1 - p2).norm();
		if (eL < minL)
			minL = eL;
		if (eL > maxL)
			maxL = eL;

		totalL += eL;
	}
	double avgL = totalL / ne;
	return avgL;
}


void OptimizerCheckboard::GetInitialValue(std::vector<double>& values)
{
	int nv = coarseMesh.n_vertices();
	int nf = coarseMesh.n_faces();
	values.resize(3 * nv, 0);	// adding 

	if (weight_Nornalization_vertexNormal > 0)
	{
		values.resize(6 * nv, 0);
	}

	std::cout << "Number of total para is :" << values.size() << std::endl;



	OpenMesh::Vec3d p;
	for (CGMesh::VertexIter v_it = coarseMesh.vertices_begin(); v_it != coarseMesh.vertices_end(); ++v_it)
	{
		p = coarseMesh.point(v_it);
		int id = v_it.handle().idx();
		for (int j = 0; j < 3; j++)
			values[3 * id + j] = p[j];
	}

	if (weight_Nornalization_vertexNormal > 0)
	{
		coarseMesh.request_vertex_normals();
		coarseMesh.request_face_normals();
		coarseMesh.update_face_normals();
		coarseMesh.update_vertex_normals();

		for (int i = 0; i < nv; i++)
		{
			CGMesh::Normal vn = coarseMesh.normal(coarseMesh.vertex_handle(i));
			values[3 * nv + 3 * i] = vn[0]; values[3 * nv + 3 * i + 1] = vn[1]; values[3 * nv + 3 * i + 2] = vn[2];
		}
	}


	Diag3pts.clear();
	for (CGMesh::VertexIter viter = coarseMesh.vertices_begin(); viter != coarseMesh.vertices_end(); ++viter) {
		if (coarseMesh.valence(viter) != 4)
			continue;
		int id0 = viter.handle().idx();
		std::vector<int> vecpts;
		for (CGMesh::VertexFaceIter vfit = coarseMesh.vf_begin(viter); vfit != coarseMesh.vf_end(viter); ++vfit)
		{
			OpenMesh::HalfedgeHandle hh;
			for (CGMesh::FaceHalfedgeIter fhhit = coarseMesh.fh_begin(vfit); fhhit != coarseMesh.fh_end(vfit); ++fhhit)
			{
				if (coarseMesh.from_vertex_handle(fhhit) == viter.handle())
				{
					hh = fhhit.handle();
					break;
				}
			}
			OpenMesh::VertexHandle vh = coarseMesh.to_vertex_handle(coarseMesh.next_halfedge_handle(hh));
			vecpts.push_back(vh.idx());
		}
		Diag3pts.push_back(OpenMesh::Vec3i(vecpts[0], id0, vecpts[2]));
		Diag3pts.push_back(OpenMesh::Vec3i(vecpts[1], id0, vecpts[3]));

	}

	// initialization of FixedVerticesFlag;
	FixedVerticesFlag.resize(nv, 0);
	if (FixedVertices.size() > 0)
	{
		for (int i = 0; i < FixedVertices.size(); i++)
			FixedVerticesFlag[FixedVertices[i]] = 1;
	}

	if (weight_closeness > 0)   // construct reference surfaces.
	{
		referenceMesh.request_vertex_normals();
		referenceMesh.request_face_normals();
		referenceMesh.update_face_normals();
		referenceMesh.update_vertex_normals();
		int _nPts = referenceMesh.n_vertices();

		_dataPts = annAllocPts(_nPts, CGMesh::Point::size());

		_pointnormals = new  CGMesh::Point[_nPts];

		for (unsigned int i = 0; i < _nPts; i++)
		{
			CGMesh::Point p = referenceMesh.point(referenceMesh.vertex_handle(i));
			for (unsigned int j = 0; j < CGMesh::Point::size(); j++)
				_dataPts[i][j] = p[j];

			p = referenceMesh.normal(referenceMesh.vertex_handle(i));
			_pointnormals[i] = p;
		}

		_kdTree = new ANNkd_tree(_dataPts, _nPts, CGMesh::Point::size());

	}
}



void OptimizerCheckboard::GetInitialTriplet(double* x, bool firstBuild, bool PrintFlag, int Iter)
{
	EqID = 0;
	TID = 0;
	Rhs.clear();
	int nv = coarseMesh.n_vertices();

	if (firstBuild)
		spTripletList.clear();

	double E_total = 0;

	double E_Closeness = 0.0;
	if (weight_closeness > 0)
	{
		double s;
		if (fairnessChangeType == 0)
			s = weight_closeness;
		else
			s = weight_closeness * (OpIter - 1 - Iter) / OpIter;


		if (referenceMesh.n_vertices() < 1)
			referenceMesh = coarseMesh;
		E_Closeness = LinearEquation_closeness(x, s, 0, firstBuild);
		E_total = E_Closeness + E_total;
	}

	double E_polyhedron = 0.0;
	if (weight_polyhedron > 0)
	{
		E_polyhedron = LinearEquation_polyhedron(x, weight_polyhedron, 0, firstBuild);
		E_total = E_polyhedron + E_total;
	}

	double E_polyhedron_vertex = 0.0;
	if (weight_polyhedron_vertex > 0)
	{
		E_polyhedron_vertex = LinearEquation_polyhedron_vertex(x, weight_polyhedron_vertex, 3 * nv, 0, firstBuild);
		E_total = E_polyhedron_vertex + E_total;
	}

	double E_VertexNormal = 0.0;
	if (weight_Nornalization_vertexNormal > 0)
	{
		E_VertexNormal = LinearEquation_Normalize_vnormal(x, weight_Nornalization_vertexNormal, 3 * nv, firstBuild);
		E_total = E_VertexNormal + E_total;
	}


	double E_OrthDiagornals = 0.0;
	if (weight_OrthDiagornals > 0)
	{
		E_OrthDiagornals = LinearEquation_OrthDiagornals(x, weight_OrthDiagornals, 0, firstBuild);
		E_total = E_OrthDiagornals + E_total;
	}

	double E_EqualDiagornals = 0.0;
	if (weight_EqualDiagornals > 0)
	{
		E_EqualDiagornals = LinearEquation_EqualDiagornals(x, weight_EqualDiagornals, 0, firstBuild);
		//E_EqualDiagornals = LinearEquation_ConstDiagornals(x, weight_EqualDiagornals, 0, firstBuild);
		//E_EqualDiagornals = LinearEquation_EqualDiagornals2D3D(x, weight_EqualDiagornals, 0, firstBuild);;
		E_total = E_EqualDiagornals + E_total;
	}

	double E_UniformDiagornals = 0.0;
	if (weight_UniformDiagornals > 0)
	{
		E_UniformDiagornals = LinearEquation_UniformDiagornals(x, weight_UniformDiagornals, 0, firstBuild);
		E_total = E_UniformDiagornals + E_total;
	}

	double E_Fairness = 0.0;
	if (weight_Fairness > 0)
	{
		double s;

		if (fairnessChangeType == 0)
			s = weight_Fairness;
		else
			s = weight_Fairness * (OpIter - 1 - Iter) / OpIter;


		E_Fairness = LinearEquation_Fairness(x, s, 0, firstBuild);



		E_total = E_Fairness + E_total;
	}

	double E_FairnessDiag = 0.0;
	if (weight_FairnessDiag > 0)
	{
		double s = weight_FairnessDiag * (OpIter - 1 - Iter) / OpIter;
		E_FairnessDiag = LinearEquation_FairnessDiag(x, s, 0, firstBuild);
		E_total = E_FairnessDiag + E_total;
	}


	if (PrintFlag)
	{
		std::cout << "************************** Iter: " << Iter << " ******************************" << std::endl;

		std::cout << "Total Enery is :: " << E_total << std::endl;

		if (weight_closeness > 0)
			std::cout << "E_Closeness is :: " << E_Closeness << std::endl;


		if (weight_polyhedron > 0)
			std::cout << "E_polyhedron is :: " << E_polyhedron << std::endl;

		if (weight_OrthDiagornals > 0)
			std::cout << "E_OrthDiagornals is :: " << E_OrthDiagornals << std::endl;

		if (weight_EqualDiagornals > 0)
			std::cout << "E_EqualDiagornals is :: " << E_EqualDiagornals << std::endl;


		if (weight_UniformDiagornals > 0)
			std::cout << "E_UniformDiagornals is :: " << E_UniformDiagornals << std::endl;

		if (weight_Fairness > 0)
			std::cout << "E_Fairness is :: " << E_Fairness << std::endl;

		if (weight_FairnessDiag > 0)
			std::cout << "E_FairnessDiag is :: " << E_FairnessDiag << std::endl;

		if (weight_polyhedron_vertex > 0)
			std::cout << "E_polyhedron_vertex is :: " << E_polyhedron_vertex << std::endl;

		if (weight_Nornalization_vertexNormal > 0)
			std::cout << "E_VertexNormal is :: " << E_VertexNormal << std::endl;

	}
}

std::vector<double> OptimizerCheckboard::fillMatandSolve_CheckBoard()
{
	clock_t tStart = clock();
	double time1 = 0;
	double time0 = 0;
	double firstBuildTime = 0;

	std::vector<double> x;

	GetInitialValue(x);

	if (referenceMesh.n_faces() > 0)
	{
		getDiagInfo(referenceMesh, DiagInfos);
	}

	int nv = coarseMesh.n_vertices();
	int	nf = coarseMesh.n_faces();

	clock_t tStart_firstbuild = clock();
	GetInitialTriplet(&x[0], true, true, 0);

	aveEL = getAvgEdgeLength(&x[0]);


	int nEqs = EqID;
	int nVars = x.size();
	spMat M(nEqs, nVars);


	std::vector<T> tripletList;
	tripletList.reserve(TID);
	for (int i = 0; i < spTripletList.size(); i++)
	{
		int vid = spTripletList[i].col / 3;
		if (vid < nv && FixedVerticesFlag[vid] == 1)
			continue;

		tripletList.push_back(T(spTripletList[i].row, spTripletList[i].col, 1.0));
	}
	M.setFromTriplets(tripletList.begin(), tripletList.end());



	firstBuildTime = (double)(clock() - tStart_firstbuild) / CLOCKS_PER_SEC;

	spMat MtM(nVars, nVars);
	MtM = M.transpose() * M;


	double MaxinDiag = 0;
	for (int i = 0; i < nVars; i++)
	{
		if (MtM.coeffRef(i, i) > MaxinDiag)
		{
			MaxinDiag = MtM.coeffRef(i, i);
		}
	}

	for (int i = 0; i < nVars; i++)
	{
		MtM.coeffRef(i, i) += 1e-8 * MaxinDiag;
	}
	MtM.makeCompressed();


	Eigen::SimplicialLLT<Eigen::SparseMatrix<double > > solver;
	
	double Engergy = 10e+10;
	for (int num = 0; num < OpIter; num++)
	{
		clock_t tStart0 = clock();

		GetInitialTriplet(&x[0], false, num == OpIter - 1, num);



		for (int i = 0; i < spTripletList.size(); i++)
		{
			int row = spTripletList[i].row;
			int col = spTripletList[i].col;
			M.coeffRef(row, col) = 0;

		}

		for (int i = 0; i < spTripletList.size(); i++)
		{
			int vid = spTripletList[i].col / 3;
			if (vid < nv && FixedVerticesFlag[vid] == 1)
				continue;

			int row = spTripletList[i].row;
			int col = spTripletList[i].col;
			M.coeffRef(row, col) += spTripletList[i].value;
		}
		M.makeCompressed();

		spMat MtM(nVars, nVars);
		MtM = M.transpose() * M;

	

		double MaxinDiag = 0;
		for (int i = 0; i < nVars; i++)
		{
			if (MtM.coeffRef(i, i) > MaxinDiag)
			{
				MaxinDiag = MtM.coeffRef(i, i);
			}
		}

		for (int i = 0; i < nVars; i++)
		{
			MtM.coeffRef(i, i) += 1e-6 * MaxinDiag;
		}
		MtM.makeCompressed();

		Eigen::VectorXd b(nEqs);
		for (uint i = 0; i < nEqs; ++i) b[i] = Rhs[i];
		Eigen::VectorXd Mtb = M.transpose() * b;

		time0 = time0 + (double)(clock() - tStart0) / CLOCKS_PER_SEC;

		clock_t tStart1 = clock();
		solver.compute(MtM);

		if (solver.info() != Eigen::Success) {
			// decomposition failed
			std::cout << "solver decomposition failed at iteration::::" << num << std::endl;
			break;
		}
		Eigen::VectorXd sol = solver.solve(Mtb);
		//Eigen::VectorXd sol = solver.solve(b);


		if (solver.info() != Eigen::Success) {
			// solving failed
			std::cout << "solver fail at iteration::::" << num << std::endl;
			break;
		}
		for (int i = 0; i < nVars; i++)
		{

			x[i] = x[i] + sol[i];
		}
		if (weight_Nornalization_vertexNormal > 0)
		{
			updateNormalization(&x[0], 3 * nv, nv);
		}

		time1 = time1 + (double)(clock() - tStart1) / CLOCKS_PER_SEC;

	}

	std::cout << "Total Time taken:" << (double)(clock() - tStart) / CLOCKS_PER_SEC << "\n";
	std::cout << "Solve Time taken:" << time1 << "\n";
	std::cout << "Later Fill Matrix Time:" << time0 << "\n";
	std::cout << "First build Time taken:" << firstBuildTime << "\n";
	std::cout << "Fill Matrix time:" << (double)(clock() - tStart) / CLOCKS_PER_SEC - time1 << "\n";

	return x;
}

void OptimizerCheckboard::getDiagInfo(CGMesh& mesh, std::vector<OpenMesh::Vec3d>& Diag3)
{
	for (CGMesh::FaceIter fiter = mesh.faces_begin(); fiter != mesh.faces_end(); ++fiter)
	{
		if (mesh.valence(fiter) != 4)
			continue;
		std::vector<OpenMesh::Vec3d> pts;
		for (CGMesh::FaceVertexIter fv_it = mesh.fv_begin(fiter.handle()); fv_it != mesh.fv_end(fiter.handle()); ++fv_it)
		{
			OpenMesh::VertexHandle vh = fv_it.handle();
			OpenMesh::Vec3d pt = OpenMesh::Vec3d(mesh.point(vh));
			pts.push_back(pt);
		}
		double d1 = (pts[2] - pts[0]) | (pts[2] - pts[0]);
		double d2 = (pts[3] - pts[1]) | (pts[3] - pts[1]);
		double d3 = (pts[3] - pts[1]) | (pts[2] - pts[0]);
		Diag3.push_back(OpenMesh::Vec3d(d1, d2, d3));
	}
}


double OptimizerCheckboard::get_closest_point(const CGMesh::Point& querypoint, CGMesh::Point& closestpoint, CGMesh::Point& closestpointnormal)
{
	ANNpoint queryPt;	// query point
	ANNidx nnIdx;		// near neigh index
	ANNdist dists;		// near neighbor distance

	queryPt = annAllocPt(CGMesh::Point::size());
	for (unsigned int i = 0; i < CGMesh::Point::size(); i++)
		queryPt[i] = querypoint[i];

	_kdTree->annkSearch(	// search
		queryPt, 			// query point
		1, 					// number of near neighbors
		&nnIdx, 			// nearest neighbors (returned)
		&dists); 			// distance (returned)

	annDeallocPt(queryPt);

	for (unsigned int i = 0; i < CGMesh::Point::size(); i++)
		closestpoint[i] = _dataPts[nnIdx][i];

	closestpointnormal = _pointnormals[nnIdx];

	return dists;
}
double OptimizerCheckboard::get_closest_point_ID(const CGMesh::Point& querypoint, int& id)
{
	ANNpoint queryPt;	// query point
	ANNidx nnIdx;		// near neigh index
	ANNdist dists;		// near neighbor distance

	queryPt = annAllocPt(CGMesh::Point::size());
	for (unsigned int i = 0; i < CGMesh::Point::size(); i++)
		queryPt[i] = querypoint[i];

	_kdTree->annkSearch(	// search
		queryPt, 			// query point
		1, 					// number of near neighbors
		&nnIdx, 			// nearest neighbors (returned)
		&dists); 			// distance (returned)

	annDeallocPt(queryPt);

	id = nnIdx;

	return dists;
}
double OptimizerCheckboard::closest_point_lineSegement(const CGMesh::Point& querypoint, CGMesh::Point& p1, CGMesh::Point& p2, CGMesh::Point& fn, CGMesh::Point& closestpoint)
{
	CGMesh::Point d1 = p2 - p1;
	d1.normalize();
	CGMesh::Point d2 = querypoint - p1;
	CGMesh::Point d3 = fn % d1;
	d3.normalize();
	CGMesh::Point e = (d2 - d3 * (d3 | d2));
	double s = (e | d1) / (p2 - p1).norm();
	if (s < 0)
	{
		closestpoint = p1;
	}
	else if (s > 1)
	{
		closestpoint = p2;
	}
	else
	{
		closestpoint = p1 + s * (p2 - p1);
	}
	double dis = (closestpoint - querypoint).norm();
	return dis;
}
double OptimizerCheckboard::closest_point_triangle(const CGMesh::Point& querypoint, std::vector<CGMesh::Point>& Tripts, CGMesh::Point& closestpoint)
{
	double dis = std::numeric_limits<double>::max();
	if (Tripts.size() == 3)
	{
		int closeID;
		double ppdis = std::numeric_limits<double>::max();
		for (int i = 0; i < 3; i++)
		{
			double d = (querypoint - Tripts[i]).norm();
			if (d < ppdis)
			{
				closeID = i;
				ppdis = d;
				closestpoint = Tripts[i];
			}
		}
		dis = ppdis;
		int i1 = (closeID + 1) % 3; int i2 = (closeID + 2) % 3;
		CGMesh::Point fn = (Tripts[i2] - Tripts[i1]) % (Tripts[closeID] - Tripts[i1]);

		fn.normalize();
		CGMesh::Point dir1 = querypoint - Tripts[closeID];
		CGMesh::Point closetPt1 = Tripts[closeID] + dir1 - (dir1 | fn) * fn;
		// if the projection point is located in the middle of the triangle
		if (((closetPt1 - Tripts[0]) % (closetPt1 - Tripts[1]) | fn) > 0 && ((closetPt1 - Tripts[1]) % (closetPt1 - Tripts[2]) | fn) > 0 && ((closetPt1 - Tripts[2]) % (closetPt1 - Tripts[0]) | fn) > 0)
		{
			closestpoint = closetPt1;
			dis = (closestpoint - querypoint).norm();
		}
		else
		{
			for (int i = 0; i < 3; i++)
			{
				int i1 = (i + 1) % 3;
				CGMesh::Point closetPt2;
				double d = closest_point_lineSegement(querypoint, Tripts[i], Tripts[i1], fn, closetPt2);
				if (d < ppdis)
				{
					ppdis = d;
					closestpoint = closetPt2;
				}
			}
			dis = (closestpoint - querypoint).norm();
		}

	}
	return dis;
}
double OptimizerCheckboard::closest_point_projection(const CGMesh::Point& querypoint, CGMesh::Point& closestpoint)
{
	double dis = std::numeric_limits<double>::max();
	if (referenceMesh.n_vertices() > 0)   // construct reference surfaces.
	{
		referenceMesh.request_vertex_normals();
		referenceMesh.request_face_normals();
		referenceMesh.update_face_normals();
		referenceMesh.update_vertex_normals();
		int _nPts = referenceMesh.n_vertices();

		_dataPts = annAllocPts(_nPts, CGMesh::Point::size());

		_pointnormals = new  CGMesh::Point[_nPts];

		for (unsigned int i = 0; i < _nPts; i++)
		{
			CGMesh::Point p = referenceMesh.point(referenceMesh.vertex_handle(i));
			for (unsigned int j = 0; j < CGMesh::Point::size(); j++)
				_dataPts[i][j] = p[j];

			p = referenceMesh.normal(referenceMesh.vertex_handle(i));
			_pointnormals[i] = p;
		}

		_kdTree = new ANNkd_tree(_dataPts, _nPts, CGMesh::Point::size());

		int vid;
		get_closest_point_ID(querypoint, vid);
		OpenMesh::VertexHandle vh = referenceMesh.vertex_handle(vid);
		CGMesh::Point Pt = referenceMesh.point(vh);

		CGMesh::Point closePt;
		for (CGMesh::VertexFaceIter vfit = referenceMesh.vf_begin(vh); vfit != referenceMesh.vf_end(vh); ++vfit)
		{
			std::vector<CGMesh::Point> fvpts;
			CGMesh::Point closePtfi;
			for (CGMesh::FaceHalfedgeIter fhhit = referenceMesh.fh_begin(vfit); fhhit != referenceMesh.fh_end(vfit); ++fhhit)
			{
				fvpts.push_back(referenceMesh.point(referenceMesh.from_vertex_handle(fhhit)));
			}
			double d = closest_point_triangle(querypoint, fvpts, closePtfi);
			if (d < dis)
			{
				dis = d;
				closePt = closePtfi;
			}
		}
		closestpoint = closePt;
		//closestpoint = Pt;
	}
	return dis;
}


