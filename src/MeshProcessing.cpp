#include "MeshProcessing.h"
#include <Eigen/Eigenvalues> 


MeshProcessing::~MeshProcessing()
{
	//finalize();
}

std::vector<int> MeshProcessing::meshCorners(CGMesh& mesh)
{
	std::vector<int> ids;
	for (CGMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		if (mesh.valence(v_it) == 2)
			ids.push_back(v_it.handle().idx());
	}
	return ids;
}

std::vector<int> MeshProcessing::MeshVertexClassificationBlackWhite(CGMesh& mesh)
{
	int nVerts = mesh.n_vertices();
	if (nVerts < 1)
	{
		std::cout << "Please Update the Reference Mesh" << std::endl;
		std::vector<int> none;
		return none;
	}
	std::vector<int> vType(nVerts, -1);

	std::vector<int> stVts;
	stVts.push_back(0);
	vType[0] = 0;

	while (stVts.size() > 0)
	{
		std::vector<int> Vts;
		for (int i = 0; i < stVts.size(); i++)
		{
			int Typei = vType[stVts[i]];
			CGMesh::VertexHandle vh = mesh.vertex_handle(stVts[i]);
			for (CGMesh::VertexVertexIter vv_it = mesh.vv_begin(vh); vv_it != mesh.vv_end(vh); ++vv_it)
			{
				int id = vv_it.handle().idx();
				if (vType[id] < 0)
				{
					Vts.push_back(id);
					vType[id] = 1 - Typei;
				}
			}
		}
		stVts.clear();
		for (int i = 0; i < Vts.size(); i++)
		{
			stVts.push_back(Vts[i]);
		}

	}
	return vType;
}

std::vector<int> MeshProcessing::MeshFaceClassificationBlackWhite(CGMesh& mesh)
{
	int nf = mesh.n_faces();
	if (nf < 1)
	{
		std::cout << "Please Update the Reference Mesh" << std::endl;
		std::vector<int> none;
		return none;
	}
	std::vector<int> fType(nf, -1);

	std::vector<int> stfs;
	stfs.push_back(0);
	fType[0] = 0;

	while (stfs.size() > 0)
	{
		std::vector<int> fs;
		for (int i = 0; i < stfs.size(); i++)
		{
			int Typei = fType[stfs[i]];
			CGMesh::FaceHandle fh = mesh.face_handle(stfs[i]);
			for (CGMesh::FaceFaceIter ff_it = mesh.ff_begin(fh); ff_it != mesh.ff_end(fh); ++ff_it)
			{
				int id = ff_it.handle().idx();
				if (fType[id] < 0)
				{
					fs.push_back(id);
					fType[id] = 1 - Typei;
				}
			}
		}
		stfs.clear();
		for (int i = 0; i < fs.size(); i++)
		{
			stfs.push_back(fs[i]);
		}

	}
	return fType;
}

std::vector<int> MeshProcessing::MeshFaceClassificationCube4types(CGMesh& mesh, int stFaceId, std::vector<int>& vType)
{
	int nf = mesh.n_faces();
	int nv = mesh.n_vertices();
	if (nf < 1)
	{
		std::cout << "Please Update the Reference Mesh" << std::endl;
		std::vector<int> none;
		return none;
	}
	std::vector<int> fType(nf, -1);
	vType.resize(nv, -1);
	//std::vector<int> vType(nv, -1);

	std::vector<int> stfs;
	stfs.push_back(stFaceId);
	fType[stFaceId] = 0;

	while (stfs.size() > 0)
	{
		//std::vector<int> fs;
		//for (int i = 0; i < stfs.size(); i++)
		//{
		//	int Typei = fType[stfs[i]];
		//	CGMesh::FaceHandle fh = mesh.face_handle(stfs[i]);
		//	for (CGMesh::FaceHalfedgeIter fh_it = mesh.fh_begin(fh); fh_it != mesh.fh_end(fh); ++fh_it)
		//	{
		//		if (!mesh.is_boundary(fh_it))
		//		{
		//			CGMesh::HalfedgeHandle oph = mesh.opposite_halfedge_handle(fh_it);
		//			CGMesh::HalfedgeHandle nexth = mesh.next_halfedge_handle(mesh.next_halfedge_handle(oph));
		//			if (!mesh.is_boundary(nexth))
		//			{
		//				oph = mesh.opposite_halfedge_handle(nexth);
		//				CGMesh::FaceHandle fh= mesh.face_handle(oph);
		//				int id = fh.idx();
		//				if (fType[id] < 0)
		//				{
		//					fs.push_back(id);
		//					fType[id] = 1 - Typei;
		//				}
		//			}
		//		}
		//	}
		//}
		//stfs.clear();
		//for (int i = 0; i < fs.size(); i++)
		//{
		//	stfs.push_back(fs[i]);
		//}

		std::vector<int> fs;
		for (int i = 0; i < stfs.size(); i++)
		{
			int Typei = fType[stfs[i]];
			CGMesh::FaceHandle fh = mesh.face_handle(stfs[i]);
			for (CGMesh::FaceHalfedgeIter fh_it = mesh.fh_begin(fh); fh_it != mesh.fh_end(fh); ++fh_it)
			{
				CGMesh::HalfedgeHandle oph = mesh.opposite_halfedge_handle(fh_it);
				if (!mesh.is_boundary(oph))
				{
					CGMesh::HalfedgeHandle nexth = mesh.next_halfedge_handle(mesh.next_halfedge_handle(oph));
					oph = mesh.opposite_halfedge_handle(nexth);
					if (!mesh.is_boundary(oph))
					{
						CGMesh::FaceHandle fh = mesh.face_handle(oph);
						int id = fh.idx();
						if (fType[id] < 0)
						{
							fs.push_back(id);
							fType[id] = 1 - Typei;
						}
					}
				}
			}
		}
		stfs.clear();
		for (int i = 0; i < fs.size(); i++)
		{
			stfs.push_back(fs[i]);
		}

	}

	for (CGMesh::FaceIter f_it = mesh.faces_begin(); f_it != (mesh.faces_end()); ++f_it)
	{
		int fid = f_it.handle().idx();
		int t = fType[fid];
		if (t > -1)
		{
			for (CGMesh::FaceVertexIter fv_it = mesh.fv_begin(f_it); fv_it != (mesh.fv_end(f_it)); ++fv_it)
			{
				int fvid = fv_it.handle().idx();
				vType[fvid] = t;
			}
		}
	}

	for (CGMesh::FaceIter f_it = mesh.faces_begin(); f_it != (mesh.faces_end()); ++f_it)
	{
		int fid = f_it.handle().idx();
		int t = fType[fid];
		if (t == 0 || t == 1)
		{
			for (CGMesh::FaceHalfedgeIter fh_it = mesh.fh_begin(f_it); fh_it != (mesh.fh_end(f_it)); ++fh_it)
			{
				CGMesh::HalfedgeHandle oph = mesh.opposite_halfedge_handle(fh_it);
				if (!mesh.is_boundary(oph))
				{
					int ffid = mesh.face_handle(oph).idx();
					fType[ffid] = 2;
					CGMesh::HalfedgeHandle nexth = mesh.next_halfedge_handle(oph);
					nexth = mesh.next_halfedge_handle(nexth);
					int vid1 = mesh.from_vertex_handle(nexth).idx();
					int vid2 = mesh.to_vertex_handle(nexth).idx();
					vType[vid1] = 1 - t; vType[vid2] = 1 - t;
				}
			}
		}
	}

	for (CGMesh::FaceIter f_it = mesh.faces_begin(); f_it != (mesh.faces_end()); ++f_it)
	{
		int fid = f_it.handle().idx();
		int t = fType[fid];
		if (t == -1)
			fType[fid] = 3;
	}

	/*for (CGMesh::FaceIter f_it = mesh.faces_begin(); f_it != (mesh.faces_end()); ++f_it)
	{
		int fid = f_it.handle().idx();
		int t = fType[fid];
		if (t == -1)
		{
			std::vector<int> fvtypes;
			int n = mesh.valence(f_it);
			for (CGMesh::FaceVertexIter fv_it = mesh.fv_begin(f_it); fv_it != (mesh.fv_end(f_it)); ++fv_it)
			{
				int fvid = fv_it.handle().idx();
				int vt=vType[fvid];
				if(vt>-1)
					fvtypes.push_back(vt);
			}
			if (fvtypes.size() == n)
				if (fvtypes[0] == fvtypes[2])
					fType[fid] = 3;
				else
					fType[fid] = 2;
			else
				std::cout << "face type is not determined" << std::endl;

		}
	}*/


	return fType;
}


void MeshProcessing::MeshUnitScale(CGMesh& originMesh, CGMesh& updatedMesh)
{
	double totalL, minL, maxL, avgL;
	getEdgeLengthStatistics(originMesh, totalL, minL, maxL, avgL);
	double s = 1 / avgL;
	updatedMesh = originMesh;
	for (CGMesh::VertexIter viter = updatedMesh.vertices_begin(); viter != updatedMesh.vertices_end(); ++viter) {
		CGMesh::Point pt = originMesh.point(viter);
		pt = pt * s;
		updatedMesh.set_point(viter, pt);
	}
}

void MeshProcessing::MeshScale(CGMesh& originMesh, CGMesh& updatedMesh, double s)
{
	updatedMesh = originMesh;
	for (CGMesh::VertexIter viter = updatedMesh.vertices_begin(); viter != updatedMesh.vertices_end(); ++viter) {
		CGMesh::Point pt = originMesh.point(viter);
		pt = pt * s;
		updatedMesh.set_point(viter, pt);
	}
}

void MeshProcessing::MeshDiag(CGMesh& originMesh, CGMesh& updatedMesh, int type)
{
	std::vector<int> Vtype = MeshVertexClassificationBlackWhite(originMesh);
	std::vector<CGMesh::VertexHandle> vhandles;
	for (CGMesh::VertexIter viter = originMesh.vertices_begin(); viter != originMesh.vertices_end(); ++viter) {
		CGMesh::Point pt = originMesh.point(viter);
		vhandles.push_back(updatedMesh.add_vertex(pt));
	}
	for (CGMesh::VertexIter viter = originMesh.vertices_begin(); viter != originMesh.vertices_end(); ++viter) {
		int vid = viter.handle().idx();
		if (Vtype[vid] == type)
		{
			std::vector<CGMesh::VertexHandle>vvhs;
			for (CGMesh::VertexVertexIter vvit = originMesh.vv_begin(viter); vvit != originMesh.vv_end(viter); ++vvit)
			{
				vvhs.push_back(vvit.handle());
			}
			if (vvhs.size() == 4)
				updatedMesh.add_face(vvhs);
		}
	}
}

spMat MeshProcessing::SubdivCoarseMesh_CaltmullClack(CGMesh& coarseMesh, CGMesh& fineMesh)
{
	std::vector<T> tripletList; // for sparse matrix connect the coarse and refine meshes.
	std::vector<OpenMesh::VertexHandle> vecvh;
	OpenMesh::VertexHandle vh;

	for (CGMesh::VertexIter v_it = coarseMesh.vertices_begin(); v_it != coarseMesh.vertices_end(); ++v_it) // new point for vertices
	{
		std::vector<OpenMesh::VertexHandle> vhs;
		std::vector<double> w;

		if (coarseMesh.is_boundary(v_it))
		{
			if (coarseMesh.valence(v_it) == 2)
			{
				vhs.push_back(v_it.handle()); w.push_back(1.0);
			}
			else
			{
				vhs.push_back(v_it.handle()); w.push_back(0.75);
				for (CGMesh::VertexVertexIter vv_it = coarseMesh.vv_begin(v_it); vv_it != coarseMesh.vv_end(v_it); ++vv_it)
				{
					if (coarseMesh.is_boundary(vv_it))
					{
						vhs.push_back(vv_it.handle()); w.push_back(0.125);
					}
				}
			}
		}
		else
		{
			// version 1
			vhs.push_back(v_it.handle()); w.push_back(9.0 / 16);
			for (CGMesh::VertexOHalfedgeIter voh_it = coarseMesh.voh_begin(v_it); voh_it != coarseMesh.voh_end(v_it); ++voh_it)
			{
				vhs.push_back(coarseMesh.to_vertex_handle(voh_it)); w.push_back(3.0 / 32);
				vhs.push_back(coarseMesh.to_vertex_handle(coarseMesh.next_halfedge_handle(voh_it))); w.push_back(1.0 / 64);
			}

			//// version 2
			//vhs.push_back(v_it.handle()); w.push_back(8.0/16);
			//for (CGMesh::VertexOHalfedgeIter voh_it = coarseMesh.voh_begin(v_it); voh_it != coarseMesh.voh_end(v_it); ++voh_it)
			//{
			//	vhs.push_back(coarseMesh.to_vertex_handle(voh_it)); w.push_back(3.0 / 32);
			//	vhs.push_back(coarseMesh.to_vertex_handle(coarseMesh.next_halfedge_handle(voh_it))); w.push_back(1.0 / 32);
			//}
		}
		OpenMesh::Vec3d p(0, 0, 0);
		for (int i = 0; i < vhs.size(); i++)
			p = coarseMesh.point(vhs[i]) * w[i] + p;

		vh = fineMesh.add_vertex(CGMesh::Point(p[0], p[1], p[2]));
		vecvh.push_back(vh);

		for (int i = 0; i < vhs.size(); i++)
			tripletList.push_back(T(vh.idx(), vhs[i].idx(), w[i]));
	}

	for (CGMesh::EdgeIter e_it = coarseMesh.edges_begin(); e_it != coarseMesh.edges_end(); ++e_it)   // new point for edge
	{
		std::vector<OpenMesh::VertexHandle> vhs;
		std::vector<double> w;
		OpenMesh::HalfedgeHandle hh = coarseMesh.halfedge_handle(e_it, 0);
		OpenMesh::VertexHandle vh1, vh2;
		vh1 = coarseMesh.from_vertex_handle(hh);
		vh2 = coarseMesh.to_vertex_handle(hh);
		if (coarseMesh.is_boundary(e_it))
		{
			vhs.push_back(vh1); vhs.push_back(vh2);
			w.push_back(0.5); w.push_back(0.5);
		}
		else
		{
			vhs.push_back(vh1); vhs.push_back(vh2);
			w.push_back(3.0 / 8); w.push_back(3.0 / 8);
			vhs.push_back(coarseMesh.to_vertex_handle(coarseMesh.next_halfedge_handle(hh))); w.push_back(1.0 / 16);
			vhs.push_back(coarseMesh.from_vertex_handle(coarseMesh.prev_halfedge_handle(hh))); w.push_back(1.0 / 16);

			hh = coarseMesh.opposite_halfedge_handle(hh);
			vhs.push_back(coarseMesh.to_vertex_handle(coarseMesh.next_halfedge_handle(hh))); w.push_back(1.0 / 16);
			vhs.push_back(coarseMesh.from_vertex_handle(coarseMesh.prev_halfedge_handle(hh))); w.push_back(1.0 / 16);
		}

		OpenMesh::Vec3d p(0, 0, 0);
		for (int i = 0; i < vhs.size(); i++)
			p = coarseMesh.point(vhs[i]) * w[i] + p;

		vh = fineMesh.add_vertex(CGMesh::Point(p[0], p[1], p[2]));
		vecvh.push_back(vh);

		for (int i = 0; i < vhs.size(); i++)
			tripletList.push_back(T(vh.idx(), vhs[i].idx(), w[i]));

	}

	for (CGMesh::FaceIter f_it = coarseMesh.faces_begin(); f_it != coarseMesh.faces_end(); ++f_it) // new point for face
	{
		std::vector<OpenMesh::VertexHandle> vhs;
		std::vector<double> w;

		int fvn = coarseMesh.valence(f_it);
		for (CGMesh::FaceVertexIter fv_it = coarseMesh.fv_begin(f_it); fv_it != coarseMesh.fv_end(f_it); ++fv_it)
		{
			vhs.push_back(fv_it.handle());
			w.push_back(1.0 / fvn);
		}
		OpenMesh::Vec3d p(0, 0, 0);
		for (int i = 0; i < vhs.size(); i++)
			p = coarseMesh.point(vhs[i]) * w[i] + p;

		vh = fineMesh.add_vertex(CGMesh::Point(p[0], p[1], p[2]));
		vecvh.push_back(vh);

		for (int i = 0; i < vhs.size(); i++)
			tripletList.push_back(T(vh.idx(), vhs[i].idx(), w[i]));
	}

	// addfaces
	int nv = coarseMesh.n_vertices();
	int ne = coarseMesh.n_edges();
	for (CGMesh::FaceIter f_it = coarseMesh.faces_begin(); f_it != coarseMesh.faces_end(); ++f_it) // new point for face
	{

		int fid = f_it.handle().idx();
		for (CGMesh::FaceHalfedgeIter fh_it = coarseMesh.fh_begin(f_it); fh_it != coarseMesh.fh_end(f_it); ++fh_it)
		{
			std::vector<OpenMesh::VertexHandle> vhs;
			int eid1 = coarseMesh.edge_handle(fh_it.current_halfedge_handle()).idx();
			int eid2 = coarseMesh.edge_handle(coarseMesh.next_halfedge_handle(fh_it)).idx();
			int vid = coarseMesh.to_vertex_handle(fh_it).idx();
			vhs.push_back(vecvh[fid + nv + ne]); vhs.push_back(vecvh[eid1 + nv]);  vhs.push_back(vecvh[vid]); vhs.push_back(vecvh[eid2 + nv]);
			fineMesh.add_face(vhs);
		}
	}

	int nvfine = fineMesh.n_vertices();
	spMat M(nvfine, nv);
	M.setFromTriplets(tripletList.begin(), tripletList.end());
	return M;
}

void MeshProcessing::MeshEdgeDual(CGMesh& coarseMesh, CGMesh& fineMesh, std::vector<int>& Vtype, double s, int type)
{
	CGMesh::VertexIter v_it;
	CGMesh::FaceIter f_it;
	CGMesh::EdgeIter e_it;
	CGMesh::FaceEdgeIter fe_it;
	CGMesh::VertexEdgeIter ve_it;

	std::vector<OpenMesh::VertexHandle> vecvh;


	OpenMesh::EdgeHandle eh;
	OpenMesh::VertexHandle vh1, vh2, vh;
	OpenMesh::Vec3d p1, p2, pev1, pev2;

	for (e_it = coarseMesh.edges_begin(); e_it != coarseMesh.edges_end(); ++e_it)
	{
		eh = e_it.handle();
		OpenMesh::HalfedgeHandle hh = coarseMesh.halfedge_handle(eh, 0);
		if (!hh.is_valid()) hh = coarseMesh.halfedge_handle(eh, 1);
		vh1 = coarseMesh.from_vertex_handle(hh);
		vh2 = coarseMesh.to_vertex_handle(hh);

		OpenMesh::Vec3d p1, p2, p12;
		p1 = coarseMesh.point(vh1);
		p2 = coarseMesh.point(vh2);
		int id1 = vh1.idx(); int id2 = vh2.idx();
		double w;
		if (Vtype[id1] == 0)
			w = s;
		else
			w = 1 - s;

		p12 = p1 * w + p2 * (1 - w);
		vh = fineMesh.add_vertex(CGMesh::Point(p12[0], p12[1], p12[2]));
		vecvh.push_back(vh);
	}

	bool flag1, flag2;
	flag1 = true;
	flag2 = true;
	if (type == 0)
	{
		flag1 = true; flag2 = false;
	}
	if (type == 1)
	{
		flag1 = false; flag2 = true;
	}

	if (flag1)
	{
		for (f_it = coarseMesh.faces_begin(); f_it != (coarseMesh.faces_end()); ++f_it)
		{
			std::vector<int> veceid;
			for (fe_it = coarseMesh.fe_begin(f_it); fe_it != (coarseMesh.fe_end(f_it)); ++fe_it)
			{
				int eid = fe_it.handle().idx();
				veceid.push_back(eid);
			}
			std::vector<OpenMesh::VertexHandle> vhs;

			for (int i = veceid.size() - 1; i > -1; i--)
			{
				vhs.push_back(vecvh[veceid[i]]);
			}
			if (vhs.size() == 4)
				fineMesh.add_face(vhs);
		}

	}

	if (flag2)
	{
		for (v_it = coarseMesh.vertices_begin(); v_it != (coarseMesh.vertices_end()); ++v_it)
		{
			//if (coarseMesh.is_boundary(v_it))
			//	continue;

			std::vector<int> veceid;
			for (ve_it = coarseMesh.ve_begin(v_it); ve_it != (coarseMesh.ve_end(v_it)); ++ve_it)
			{
				int eid = ve_it.handle().idx();
				veceid.push_back(eid);
			}
			if (veceid.size() == 2)
				continue;

			std::vector<OpenMesh::VertexHandle> vhs;
			if (veceid.size() == 2)
			{
				OpenMesh::VertexHandle vh = fineMesh.add_vertex(coarseMesh.point(v_it));
				vhs.push_back(vh);
			}

			for (int i = 0; i < veceid.size(); i++)
			{
				vhs.push_back(vecvh[veceid[i]]);
			}
			if (vhs.size() == 4)
				fineMesh.add_face(vhs);
		}
	}
}

void MeshProcessing::matrix2Mesh(CGMesh& mesh, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
	mesh.clear();
	int nv = V.rows();
	std::vector<CGMesh::VertexHandle> vhs;
	for (int i = 0; i < nv; i++)
		vhs.push_back(mesh.add_vertex(CGMesh::Point(V(i, 0), V(i, 1), V(i, 2))));

	int nf = F.rows();
	int fv = F.cols();
	for (int i = 0; i < nf; i++)
	{
		std::vector<CGMesh::VertexHandle>  face_vhandles;
		face_vhandles.push_back(vhs[F(i, 0)]);
		face_vhandles.push_back(vhs[F(i, 1)]);
		face_vhandles.push_back(vhs[F(i, 2)]);
		if (fv == 4)
		{
			face_vhandles.push_back(vhs[F(i, 3)]);
		}
		mesh.add_face(face_vhandles);
	}
}


void MeshProcessing::mesh2Matrix(CGMesh& mesh, Eigen::MatrixXd& V, Eigen::MatrixXi& F)
{
	int nv = mesh.n_vertices();
	V = Eigen::MatrixXd(nv, 3);
	for (CGMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
		CGMesh::Point pt = mesh.point(*v_it);
		int i = v_it.handle().idx();
		V.row(i) << pt[0], pt[1], pt[2];
	}

	int n = 0;
	std::vector<int> Fvids;
	for (CGMesh::FaceIter f_it = mesh.faces_begin(); f_it != mesh.faces_end(); ++f_it)
	{
		std::vector<int> vhs;
		for (CGMesh::FaceVertexIter fv_it = mesh.fv_begin(*f_it); fv_it != mesh.fv_end(*f_it); ++fv_it)
			vhs.push_back(fv_it.handle().idx());

		for (int i = 0; i < vhs.size() - 2; i++)
		{
			int i1 = i + 1;
			int i2 = i + 2;
			Fvids.push_back(vhs[0]); Fvids.push_back(vhs[i1]); Fvids.push_back(vhs[i2]);
			//F.row(n) << vhs[0], vhs[i1], vhs[i2]; 
			n++;
		}
	}
	F = Eigen::MatrixXi(n, 3);
	for (int i = 0; i < n; i++)
		F.row(i) << Fvids[3 * i], Fvids[3 * i + 1], Fvids[3 * i + 2];
}


void MeshProcessing::meshEdges(CGMesh& mesh, Eigen::MatrixXi& E)
{
	int ne = mesh.n_edges();
	E = Eigen::MatrixXi(ne, 2);
	int i = 0;
	for (CGMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it)
	{
		OpenMesh::HalfedgeHandle hh = mesh.halfedge_handle(e_it, 0);
		if (!hh.is_valid()) hh = mesh.halfedge_handle(e_it, 1);
		int i1 = mesh.from_vertex_handle(hh).idx();
		int i2 = mesh.to_vertex_handle(hh).idx();
		E.row(i) << i1, i2;
		i++;
	}
}


void MeshProcessing::MeshUpdate(CGMesh& originMesh, CGMesh& updatedMesh, std::vector<double>& pts)
{
	int nVerts = originMesh.n_vertices();
	for (int i = 0; i < nVerts; i++)
	{
		CGMesh::Point pt = CGMesh::Point(pts[3 * i], pts[3 * i + 1], pts[3 * i + 2]);
		updatedMesh.add_vertex(pt);
	}

	for (CGMesh::FaceIter fiter = originMesh.faces_begin(); fiter != originMesh.faces_end(); ++fiter) {
		std::vector<CGMesh::VertexHandle> faceVertexHandles;
		for (CGMesh::FaceVertexIter fv_it = originMesh.fv_begin(fiter.handle()); fv_it != originMesh.fv_end(fiter.handle()); ++fv_it)
			faceVertexHandles.push_back(fv_it.handle());

		updatedMesh.add_face(faceVertexHandles);
	}
}

void MeshProcessing::MeshUpdate(CGMesh& originMesh, CGMesh& updatedMesh, double* pts)
{
	int nVerts = originMesh.n_vertices();
	for (int i = 0; i < nVerts; i++)
	{
		CGMesh::Point pt = CGMesh::Point(pts[3 * i], pts[3 * i + 1], pts[3 * i + 2]);
		updatedMesh.add_vertex(pt);
	}

	for (CGMesh::FaceIter fiter = originMesh.faces_begin(); fiter != originMesh.faces_end(); ++fiter) {
		std::vector<CGMesh::VertexHandle> faceVertexHandles;
		for (CGMesh::FaceVertexIter fv_it = originMesh.fv_begin(fiter.handle()); fv_it != originMesh.fv_end(fiter.handle()); ++fv_it)
			faceVertexHandles.push_back(fv_it.handle());

		updatedMesh.add_face(faceVertexHandles);
	}
}


void MeshProcessing::getEdgeLengthStatistics(CGMesh& mesh, double& totalL, double& minL, double& maxL, double& avgL)
{
	int ne = mesh.n_edges();
	totalL = 0;
	minL = std::numeric_limits<double>::max();
	maxL = 0;
	for (CGMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it)
	{
		//OpenMesh::HalfedgeHandle hh = mesh.halfedge_handle(e_it, 0);
		//if (!hh.is_valid()) hh = mesh.halfedge_handle(e_it, 1);
		//OpenMesh::Vec3d p1 = mesh.point(mesh.from_vertex_handle(hh));
		//OpenMesh::Vec3d p2 = mesh.point(mesh.to_vertex_handle(hh));
		double eL = mesh.calc_edge_length(e_it);
		if (eL < minL)
			minL = eL;
		if (eL > maxL)
			maxL = eL;

		totalL += eL;
	}
	avgL = totalL / ne;
}


void MeshProcessing::subdivision_EdgeSplit(CGMesh& coarseMesh, CGMesh& fineMesh)
{
	int nVerts = coarseMesh.n_vertices();
	std::vector<CGMesh::VertexHandle> vhandles;
	for (CGMesh::VertexIter viter = coarseMesh.vertices_begin(); viter != coarseMesh.vertices_end(); ++viter) {
		CGMesh::Point pt = coarseMesh.point(viter);
		vhandles.push_back(fineMesh.add_vertex(pt));
	}

	int nFaces = coarseMesh.n_faces();
	for (CGMesh::FaceIter fiter = coarseMesh.faces_begin(); fiter != coarseMesh.faces_end(); ++fiter) {
		CGMesh::Point fpt(0, 0, 0);
		for (CGMesh::FaceVertexIter fv_iter = coarseMesh.fv_iter(fiter); fv_iter; ++fv_iter) {
			CGMesh::Point pt = coarseMesh.point(fv_iter);
			fpt += pt;
		}
		fpt /= coarseMesh.valence(fiter);
		vhandles.push_back(fineMesh.add_vertex(fpt));
	}

	int nEdges = coarseMesh.n_edges();
	for (CGMesh::EdgeIter eiter = coarseMesh.edges_begin(); eiter != coarseMesh.edges_end(); ++eiter) {
		CGMesh::HalfedgeHandle heHandle = coarseMesh.halfedge_handle(eiter, 0);
		CGMesh::Point pt1 = coarseMesh.point(coarseMesh.from_vertex_handle(heHandle));
		CGMesh::Point pt2 = coarseMesh.point(coarseMesh.to_vertex_handle(heHandle));
		CGMesh::Point ept = (pt1 + pt2) / 2;
		vhandles.push_back(fineMesh.add_vertex(ept));
	}

	for (CGMesh::FaceIter fiter = coarseMesh.faces_begin(); fiter != coarseMesh.faces_end(); ++fiter) {
		int fIdx = fiter.handle().idx();
		for (CGMesh::FaceHalfedgeIter fe_iter = coarseMesh.fh_iter(fiter); fe_iter; ++fe_iter) {

			std::vector<CGMesh::VertexHandle> faceVertexHandles;

			int vIdx = coarseMesh.from_vertex_handle(fe_iter).idx();
			int prevEdgeIdx = coarseMesh.edge_handle(coarseMesh.prev_halfedge_handle(fe_iter)).idx();
			int currEdgeIdx = coarseMesh.edge_handle(fe_iter.current_halfedge_handle()).idx();

			faceVertexHandles.push_back(vhandles[nVerts + fIdx]);
			faceVertexHandles.push_back(vhandles[nVerts + nFaces + prevEdgeIdx]);
			faceVertexHandles.push_back(vhandles[vIdx]);
			faceVertexHandles.push_back(vhandles[nVerts + nFaces + currEdgeIdx]);
			fineMesh.add_face(faceVertexHandles);
		}
	}
}

// mesh planarity
void MeshProcessing::meshPlanarity(CGMesh& inputMesh, std::vector<double>& planarity)
{
	planarity.clear();
	double p;
	for (CGMesh::FaceIter f_it = inputMesh.faces_begin(); f_it != (inputMesh.faces_end()); ++f_it)
	{
		if (inputMesh.valence(f_it) == 3)
			p = 0;
		else
		{

			std::vector<OpenMesh::Vec3d> pts;
			for (CGMesh::FaceVertexIter fv_it = inputMesh.fv_begin(f_it.handle()); fv_it != inputMesh.fv_end(f_it.handle()); ++fv_it)
				pts.push_back(OpenMesh::Vec3d(inputMesh.point(fv_it)));

			p = 0;
			int n = 0;
			for (int i = 0; i <= pts.size() - 4; i++)
			{
				OpenMesh::Vec3d pt1, pt2, pt3, pt4, c;
				double f, s;
				pt1 = pts[i]; pt2 = pts[i + 1]; pt3 = pts[i + 2]; pt4 = pts[i + 3];
				c = (pt1 - pt3) % (pt2 - pt4);
				c.normalize();
				f = fabs((pt2 - pt1) | c);
				s = ((pt1 - pt3).norm() + (pt2 - pt4).norm()) / 2;
				p = p + f / s;
				n++;
			}
			p = p / n;
		}
		planarity.push_back(p);
	}
}


void MeshProcessing::meshFaceCenter(CGMesh& inputMesh, std::vector<OpenMesh::Vec3d>& fcen)
{
	fcen.clear();
	OpenMesh::Vec3d fceni(0, 0, 0);
	for (CGMesh::FaceIter f_it = inputMesh.faces_begin(); f_it != (inputMesh.faces_end()); ++f_it)
	{
		int num = 0;
		for (CGMesh::FaceVertexIter fv_it = inputMesh.fv_begin(f_it.handle()); fv_it != inputMesh.fv_end(f_it.handle()); ++fv_it)
		{
			fceni = fceni + inputMesh.point(fv_it);
			num++;
		}
		fceni = fceni / num;
		fcen.push_back(fceni);
	}
}


