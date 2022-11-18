#include <lsc/basic.h>
#include <lsc/tools.h>

// directions[0] is the inward pointing to the point, and direction[1] is outward shooting from the point. 
// handles are the two halfedges opposite to the point. 
// l2s shows if the value is from large to small along the halfedge direction
bool find_active_faces_and_directions_around_ver(CGMesh& lsmesh, const Eigen::VectorXd& fvalues,
	const Eigen::MatrixXd& V,
	const Eigen::Vector3d& normal, const int vid,
	std::array<Eigen::Vector3d, 2>& directions,
	std::array<bool, 2>& l2s,// if the halfedge handle is from large to small
	std::array<CGMesh::HalfedgeHandle, 2>& handles)
{
	CGMesh::VertexHandle vh = lsmesh.vertex_handle(vid);
	std::vector<Eigen::Vector3d> dircs;
	std::vector<int> order;// 0 means the first segment, 1 means the second segment
	std::vector<CGMesh::HalfedgeHandle> hes;
	std::vector<bool> large_to_small;
	double value = fvalues[vid];
	for (CGMesh::VertexOHalfedgeIter voh = lsmesh.voh_begin(vh); voh != lsmesh.voh_end(vh); ++voh) {
		CGMesh::HalfedgeHandle next_he = lsmesh.next_halfedge_handle(voh);
		int id_from = lsmesh.from_vertex_handle(next_he).idx();
		int id_to = lsmesh.to_vertex_handle(next_he).idx();
		double value1 = fvalues[id_from];
		double value2 = fvalues[id_to];
		if (value1 == value2) {
			continue;
		}
		double value_min;
		double value_max;
		int idmin, idmax;
		bool tmp_l2s;
		if (value1 < value2) {
			value_min = value1;
			value_max = value2;
			idmin = id_from;
			idmax = id_to;
			tmp_l2s = false;
		}
		else {
			value_min = value2;
			value_max = value1;
			idmin = id_to;
			idmax = id_from;
			tmp_l2s = true;
		}
		if (value < value_min || value >= value_max) {
			continue;
		}
		double t;
		t = get_t_of_value(value, value1, value2);
		Eigen::Vector3d vedge = get_3d_ver_from_t(t, V.row(id_from), V.row(id_to));
		// Eigen::Vector3d big2small=V.row(idmin)-V.row(idmax);
		Eigen::Vector3d edge2ver = Eigen::Vector3d(V.row(vid)) - vedge;
		// double flag=big2small.cross(edge2ver).dot(normal);
		if (tmp_l2s) {// the halfedge is from large to small, it is the shoot-in direction
			dircs.push_back(edge2ver);
			order.push_back(0);
		}
		else {// it is the shoot-out direction
			dircs.push_back(-edge2ver);
			order.push_back(1);
		}
		int fid = lsmesh.face_handle(next_he).idx();
		hes.push_back(next_he);
		large_to_small.push_back(tmp_l2s);

	}
	if (dircs.size() != 2) {// it is a singularity
		return false;
	}
	assert(order[0] + order[1] == 1);
	for (int i = 0; i < 2; i++) {
		handles[order[i]] = hes[i];
		directions[order[i]] = dircs[i];
		l2s[order[i]] = large_to_small[i];
	}
	assert(lsmesh.face_handle(handles[0]).idx() != lsmesh.face_handle(handles[1]).idx() && "the two halfedges correspond to different faces");
	return true;

}
// TODO this can be pre-computed before iterations: for each triangle we compute 3 vectors
Eigen::Vector3d get_coff_vec_for_gradient(std::array<spMat, 3>& gradVF, const int fid, const int vid)
{
	Eigen::Vector3d result;
	result[0] = gradVF[0].coeffRef(fid, vid);
	result[1] = gradVF[1].coeffRef(fid, vid);
	result[2] = gradVF[2].coeffRef(fid, vid);
	return result;
}
// v1, v2, vM is one triangle. v3, v4, vM is one triangle. 
// This code works only for v1, v2, v3, v4 are not duplicated
std::vector<Trip> get_triplets_for_gradient_partial_parts_vertex_based(std::array<spMat, 3>& gradVF, const int v1, const int v2, const int v3, const int v4, const int vM,
	const int fid1, const int fid2, const Eigen::Vector3d& norm, const Eigen::Vector3d& norm1,
	const Eigen::Vector3d& norm2) {
	std::vector<Trip> triplets;
	triplets.reserve(17);
	Eigen::Vector3d d1, d2, dM, b3, b4, bM; // the coffecient vectors of each point for assembling the gradients.
	d1 = get_coff_vec_for_gradient(gradVF, fid1, v1);
	d2 = get_coff_vec_for_gradient(gradVF, fid1, v2);
	dM = get_coff_vec_for_gradient(gradVF, fid1, vM);
	b3 = get_coff_vec_for_gradient(gradVF, fid2, v3);
	b4 = get_coff_vec_for_gradient(gradVF, fid2, v4);
	bM = get_coff_vec_for_gradient(gradVF, fid2, vM);

	d1 = norm1.cross(d1);
	d2 = norm1.cross(d2);
	dM = norm1.cross(dM);
	b3 = norm2.cross(b3);
	b4 = norm2.cross(b4);
	bM = norm2.cross(bM);

	double cfmfm = dM.cross(bM).dot(norm);
	double cf1fm = d1.cross(bM).dot(norm);
	double cf2fm = d2.cross(bM).dot(norm);
	double cfmf3 = dM.cross(b3).dot(norm);
	double cf1f3 = d1.cross(b3).dot(norm);
	double cf2f3 = d2.cross(b3).dot(norm);
	double cfmf4 = dM.cross(b4).dot(norm);
	double cf1f4 = d1.cross(b4).dot(norm);
	double cf2f4 = d2.cross(b4).dot(norm);

	// derivate v1
	triplets.push_back(Trip(v1, v3, cf1f3));
	triplets.push_back(Trip(v1, v4, cf1f4));
	triplets.push_back(Trip(v1, vM, cf1fm));
	// derivate v2
	triplets.push_back(Trip(v2, v3, cf2f3));
	triplets.push_back(Trip(v2, v4, cf2f4));
	triplets.push_back(Trip(v2, vM, cf2fm));
	// derivate v3
	triplets.push_back(Trip(v3, v1, cf1f3));
	triplets.push_back(Trip(v3, v2, cf2f3));
	triplets.push_back(Trip(v3, vM, cfmf3));
	// derivate v4
	triplets.push_back(Trip(v4, v1, cf1f4));
	triplets.push_back(Trip(v4, v2, cf2f4));
	triplets.push_back(Trip(v4, vM, cfmf4));
	// derivate vM
	triplets.push_back(Trip(vM, v1, cf1fm));
	triplets.push_back(Trip(vM, v2, cf2fm));
	triplets.push_back(Trip(vM, v3, cfmf3));
	triplets.push_back(Trip(vM, v4, cfmf4));
	triplets.push_back(Trip(vM, vM, 2 * cfmfm));
	return triplets;
}

std::vector<Trip> get_triplets_for_gradient_partial_parts_edge_based(std::array<spMat, 3>& gradVF,
	const int v1, const int v2, const int vA, const int vB,
	const int fid1, const int fid2, const Eigen::Vector3d& norm, const Eigen::Vector3d& norm1,
	const Eigen::Vector3d& norm2)
{
	std::vector<Trip> triplets;
	triplets.reserve(16);

	Eigen::Vector3d d1, d2, dA, b1, b2, bB; // the coffecient vectors of each point for assembling the gradients.
	d1 = get_coff_vec_for_gradient(gradVF, fid1, v1);
	d2 = get_coff_vec_for_gradient(gradVF, fid1, v2);
	dA = get_coff_vec_for_gradient(gradVF, fid1, vA);
	b1 = get_coff_vec_for_gradient(gradVF, fid2, v1);
	b2 = get_coff_vec_for_gradient(gradVF, fid2, v2);
	bB = get_coff_vec_for_gradient(gradVF, fid2, vB);

	// rotate the gradient to present the iso-line direction
	d1 = norm1.cross(d1);
	d2 = norm1.cross(d2);
	dA = norm1.cross(dA);
	b1 = norm2.cross(b1);
	b2 = norm2.cross(b2);
	bB = norm2.cross(bB);

	double cf1f1 = d1.cross(b1).dot(norm);
	double cf1f2 = (d2.cross(b1) + d1.cross(b2)).dot(norm);
	double cf2f2 = d2.cross(b2).dot(norm);
	double cfaf1 = dA.cross(b1).dot(norm);
	double cfaf2 = dA.cross(b2).dot(norm);
	double cfafb = dA.cross(bB).dot(norm);
	double cf1fb = d1.cross(bB).dot(norm);
	double cf2fb = d2.cross(bB).dot(norm);
	// derivate v1
	triplets.push_back(Trip(v1, v1, 2 * cf1f1));
	triplets.push_back(Trip(v1, v2, cf1f2));
	triplets.push_back(Trip(v1, vA, cfaf1));
	triplets.push_back(Trip(v1, vB, cf1fb));

	// derivate v2
	triplets.push_back(Trip(v2, v1, cf1f2));
	triplets.push_back(Trip(v2, v2, 2 * cf2f2));
	triplets.push_back(Trip(v2, vA, cfaf2));
	triplets.push_back(Trip(v2, vB, cf2fb));

	// derivate vA
	triplets.push_back(Trip(vA, v1, cfaf1));
	triplets.push_back(Trip(vA, v2, cfaf2));
	triplets.push_back(Trip(vA, vB, cfafb));

	// derivate vB
	triplets.push_back(Trip(vB, v1, cf1fb));
	triplets.push_back(Trip(vB, v2, cf2fb));
	triplets.push_back(Trip(vB, vA, cfafb));
	return triplets;
}
// calculate vector of matrix A, A[i]*f is the jacobian of (g1xg2).dot(norm).
// here g1, g2 are the rotated gradients
void lsTools::calculate_gradient_partial_parts() {

	int act_size = Actid.size();
	PEJC.resize(act_size);
	for (int i = 0; i < act_size; i++) {
		PEJC[i].resize(V.rows(), V.rows());
	}
	// std::array<Eigen::MatrixXd,8> 
	for (int i = 0; i < act_size; i++) {
		int fid1, fid2;
		std::array<int, 4> vids;
		// get the v1, v2, vA, vB.
		vids = get_vers_around_edge(lsmesh, Actid[i], fid1, fid2);
		int v1 = vids[0];
		int v2 = vids[1];
		int vA = vids[2];
		int vB = vids[3];

		Eigen::Vector3d norm = norm_e.row(Actid[i]);
		Eigen::Vector3d norm1 = norm_f.row(fid1);
		Eigen::Vector3d norm2 = norm_f.row(fid2);
		std::vector<Trip> triplets = get_triplets_for_gradient_partial_parts_edge_based(gradVF, v1, v2, vA, vB,
			fid1, fid2, norm, norm1, norm2);

		PEJC[i].setFromTriplets(triplets.begin(), triplets.end());
	}
}
// calculate vector of matrix A, A[i]*f is the jacobian of (g1xg2).dot(norm).
// here g1, g2 are the rotated gradients
void lsTools::calculate_gradient_partial_parts_ver_based() {
	int ninner = IVids.size();
	int vnbr = V.rows();
	bool first_time_compute = false; // if it is first time, we need to compute all; otherwise we just need to update
	if (VPEJC.empty())
	{
		VPEJC.resize(ninner);
		Vheh0.resize(ninner);
		Vheh1.resize(ninner);
		Vdire0.resize(ninner);
		Vdire1.resize(ninner);
		ActInner.resize(ninner);
		first_time_compute = true;
	}
	if (Last_Opt_Mesh) {
		first_time_compute = true;
	}


	for (int i = 0; i < ninner; i++)
	{
		int vid = IVids[i];
		std::array<Eigen::Vector3d, 2> directions;
		std::array<bool, 2> l2s; // if the halfedge handle is from large to small
		std::array<CGMesh::HalfedgeHandle, 2> handles;
		bool active = find_active_faces_and_directions_around_ver(lsmesh, fvalues, V, norm_v.row(vid), vid, directions,
			l2s, handles);
		if (active)
		{
			ActInner[i] = true;
		}
		else
		{
			ActInner[i] = false;
			continue;
		}
		CGMesh::HalfedgeHandle hd1 = handles[0]; // the inward edge
		CGMesh::HalfedgeHandle hd2 = handles[1]; // the outward edge
		int v1 = lsmesh.from_vertex_handle(hd1).idx();
		int v2 = lsmesh.to_vertex_handle(hd1).idx();
		int v3 = lsmesh.from_vertex_handle(hd2).idx();
		int v4 = lsmesh.to_vertex_handle(hd2).idx();
		int fid1 = lsmesh.face_handle(hd1).idx();
		int fid2 = lsmesh.face_handle(hd2).idx();
		assert(fid1 != fid2);
		bool need_update = false;
		if (first_time_compute)
		{
			need_update = true;
		}
		else
		{
			if (hd1 == Vheh0[i] && hd2 == Vheh1[i])
			{ // if the active information is the same as the last iteration, no need to update the matrices
				continue;
			}
			need_update = true;
		}
		if (need_update)
		{
			Eigen::Vector3d norm = norm_v.row(vid);
			Eigen::Vector3d norm1 = norm_f.row(fid1);
			Eigen::Vector3d norm2 = norm_f.row(fid2);
			std::vector<Trip> triplets;
			if (v2 == v3)
			{
				triplets = get_triplets_for_gradient_partial_parts_edge_based(gradVF, vid, v2, v1, v4, fid1, fid2, norm, norm1, norm2);
			}
			else if (v1 == v4)
			{
				triplets = get_triplets_for_gradient_partial_parts_edge_based(gradVF, v1, vid, v2, v3, fid1, fid2, norm, norm1, norm2);
			}
			else
			{
				triplets = get_triplets_for_gradient_partial_parts_vertex_based(gradVF, v1, v2, v3, v4, vid, fid1, fid2, norm, norm1, norm2);
			}
			VPEJC[i].resize(vnbr, vnbr);
			VPEJC[i].setFromTriplets(triplets.begin(), triplets.end());
			Vheh0[i] = hd1;
			Vheh1[i] = hd2;
			Vdire0[i] = directions[0];
			Vdire1[i] = directions[1];
		}
	}
}
void lsTools::analysis_pseudo_geodesic_on_vertices(const Eigen::VectorXd& func_values, Eigen::VectorXd& LocalActInner,
	std::vector<CGMesh::HalfedgeHandle>& heh0, std::vector<CGMesh::HalfedgeHandle>& heh1,
	std::vector<double>& t1s, std::vector<double>& t2s) {
	int ninner = IVids.size();
	int vnbr = V.rows();
	if (LocalActInner.size() == 0 || heh0.empty() || t1s.empty()) {
		LocalActInner.resize(ninner);
		heh0.resize(ninner);
		heh1.resize(ninner);
		t1s.resize(ninner);
		t2s.resize(ninner);

	}
	for (int i = 0; i < ninner; i++)
	{
		int vid = IVids[i];
		std::array<Eigen::Vector3d, 2> directions;
		std::array<bool, 2> l2s; // if the halfedge handle is from large to small
		std::array<CGMesh::HalfedgeHandle, 2> handles;
		bool active = find_active_faces_and_directions_around_ver(lsmesh, func_values, V, norm_v.row(vid), vid, directions,
			l2s, handles);
		if (active)
		{
			LocalActInner[i] = true;
		}
		else
		{
			LocalActInner[i] = false;
			continue;
		}
		CGMesh::HalfedgeHandle hd1 = handles[0]; // the inward edge
		CGMesh::HalfedgeHandle hd2 = handles[1]; // the outward edge
		int v1 = lsmesh.from_vertex_handle(hd1).idx();
		int v2 = lsmesh.to_vertex_handle(hd1).idx();
		int v3 = lsmesh.from_vertex_handle(hd2).idx();
		int v4 = lsmesh.to_vertex_handle(hd2).idx();
		int fid1 = lsmesh.face_handle(hd1).idx();
		int fid2 = lsmesh.face_handle(hd2).idx();
		double t1 = get_t_of_value(func_values[vid], func_values[v1], func_values[v2]);
		double t2 = get_t_of_value(func_values[vid], func_values[v3], func_values[v4]);

		assert(func_values[v1] > func_values[v2]);
		assert(fid1 != fid2);
		heh0[i] = hd1;
		heh1[i] = hd2;
		t1s[i] = t1;
		t2s[i] = t2;
	}
}
void lsTools::analysis_pseudo_geodesic_on_vertices(const Eigen::VectorXd& func_values, LSAnalizer& ana) {
	analysis_pseudo_geodesic_on_vertices(func_values, ana.LocalActInner, ana.heh0, ana.heh1, ana.t1s, ana.t2s);
}
void lsTools::calculate_pseudo_geodesic_opt_expanded_function_values(const Eigen::VectorXd& vars, const std::vector<double>& angle_degree,
	const Eigen::VectorXd& LocalActInner, const std::vector<CGMesh::HalfedgeHandle>& heh0, const std::vector<CGMesh::HalfedgeHandle>& heh1,
	const std::vector<double>& t1s, const std::vector<double>& t2s, const bool first_compute, const int aux_start_loc, std::vector<Trip>& tripletes, Eigen::VectorXd& Energy) {
	double cos_angle;
	if (angle_degree.size() == 1) {
		double angle_radian = angle_degree[0] * LSC_PI / 180.; // the angle in radian
		cos_angle = cos(angle_radian);
	}
	int vnbr = V.rows();
	int ninner = LocalActInner.size();
	tripletes.clear();
	tripletes.reserve(ninner * 20); // the number of rows is ninner*4, the number of cols is aux_start_loc + ninner * 3 (all the function values and auxiliary vars)
	Energy = Eigen::VectorXd::Zero(ninner * 4); // mesh total energy values

	for (int i = 0; i < ninner; i++)
	{
		if (LocalActInner[i] == false) {
			continue;
		}
		int vm = IVids[i];
		if (angle_degree.size() == vnbr) {
			double angle_radian = angle_degree[vm] * LSC_PI / 180.; // the angle in radian
			cos_angle = cos(angle_radian);
		}
		CGMesh::HalfedgeHandle inhd = heh0[i], outhd = heh1[i];
		int v1 = lsmesh.from_vertex_handle(inhd).idx();
		int v2 = lsmesh.to_vertex_handle(inhd).idx();
		int v3 = lsmesh.from_vertex_handle(outhd).idx();
		int v4 = lsmesh.to_vertex_handle(outhd).idx();
		double t1 = t1s[i];
		double t2 = t2s[i];
		Eigen::Vector3d ver0 = V.row(v1) + (V.row(v2) - V.row(v1)) * t1;
		Eigen::Vector3d ver1 = V.row(vm);
		Eigen::Vector3d ver2 = V.row(v3) + (V.row(v4) - V.row(v3)) * t2;
		// the locations
		int lrx = i + aux_start_loc;
		int lry = i + aux_start_loc + ninner;
		int lrz = i + aux_start_loc + ninner * 2;
		if (first_compute) {
			Eigen::Vector3d real_r = (ver1 - ver0).cross(ver2 - ver1);
			real_r = real_r.normalized();
			vars[lrx] = real_r[0];
			vars[lry] = real_r[1];
			vars[lrz] = real_r[2];
		}
		Eigen::Vector3d r = Eigen::Vector3d(vars[lrx], vars[lry], vars[lrz]);
		// r dot (vm+(t1-1)*vf-t1*vt)
		// vf = v1, vt = v2 
		tripletes.push_back(Trip(i, lrx, (V(v1, 0) - V(v2, 0)) * vars[vm] + (V(v2, 0) - V(vm, 0)) * vars[v1] + (V(vm, 0) - V(v1, 0)) * vars[v2]));
		tripletes.push_back(Trip(i, lry, (V(v1, 1) - V(v2, 1)) * vars[vm] + (V(v2, 1) - V(vm, 1)) * vars[v1] + (V(vm, 1) - V(v1, 1)) * vars[v2]));
		tripletes.push_back(Trip(i, lrz, (V(v1, 2) - V(v2, 2)) * vars[vm] + (V(v2, 2) - V(vm, 2)) * vars[v1] + (V(vm, 2) - V(v1, 2)) * vars[v2]));

		double r12 = (V.row(v1) - V.row(v2)).dot(r);
		double rm1 = (V.row(vm) - V.row(v1)).dot(r);
		double r2m = (V.row(v2) - V.row(vm)).dot(r);
		tripletes.push_back(Trip(i, vm, r12));
		tripletes.push_back(Trip(i, v1, r2m));
		tripletes.push_back(Trip(i, v2, rm1));
		Energy[i] = r12 * vars[vm] + rm1 * vars[v2] + r2m * vars[v1];

		// vf = v3, vt = v4 
		tripletes.push_back(Trip(i + ninner, lrx, (V(v3, 0) - V(v4, 0)) * vars[vm] + (V(v4, 0) - V(vm, 0)) * vars[v1] + (V(vm, 0) - V(v3, 0)) * vars[v4]));
		tripletes.push_back(Trip(i + ninner, lry, (V(v3, 1) - V(v4, 1)) * vars[vm] + (V(v4, 1) - V(vm, 1)) * vars[v1] + (V(vm, 1) - V(v3, 1)) * vars[v4]));
		tripletes.push_back(Trip(i + ninner, lrz, (V(v3, 2) - V(v4, 2)) * vars[vm] + (V(v4, 2) - V(vm, 2)) * vars[v1] + (V(vm, 2) - V(v3, 2)) * vars[v4]));

		r12 = (V.row(v3) - V.row(v4)).dot(r);
		rm1 = (V.row(vm) - V.row(v3)).dot(r);
		r2m = (V.row(v4) - V.row(vm)).dot(r);
		tripletes.push_back(Trip(i + ninner, vm, r12));
		tripletes.push_back(Trip(i + ninner, v3, r2m));
		tripletes.push_back(Trip(i + ninner, v4, rm1));

		Energy[i + ninner] = r12 * vars[vm] + rm1 * vars[v4] + r2m * vars[v3];

		// r*r=1
		tripletes.push_back(Trip(i + ninner * 2, lrx, 2 * vars(lrx)));
		tripletes.push_back(Trip(i + ninner * 2, lry, 2 * vars(lry)));
		tripletes.push_back(Trip(i + ninner * 2, lrz, 2 * vars(lrz)));

		Energy[i + ninner * 2] = r.dot(r) - 1;

		// r*norm - cos = 0
		Eigen::Vector3d norm = norm_v.row(vm);
		tripletes.push_back(Trip(i + ninner * 3, lrx, norm(0)));
		tripletes.push_back(Trip(i + ninner * 3, lry, norm(1)));
		tripletes.push_back(Trip(i + ninner * 3, lrz, norm(2)));

		Energy[i + ninner * 3] = norm.dot(r) - cos_angle;
	}

}
// this function only need be called after initializing the level set
void lsTools::get_traced_boundary_triangle_direction_derivatives() {
	int size = tracing_start_edges.size();

	spMat jacobian;
	jacobian.resize(size, V.rows());
	std::vector<Trip> triplets;
	triplets.reserve(size * 3);
	refids.clear();
	refids.reserve(size * 3);
	for (int i = 0; i < size; i++) {
		CGMesh::HalfedgeHandle hd = tracing_start_edges[i];
		CGMesh::HalfedgeHandle op = lsmesh.opposite_halfedge_handle(hd);
		// get the face
		int fid = lsmesh.face_handle(hd).idx();
		if (fid < 0) {
			fid = lsmesh.face_handle(op).idx();
		}
		int idfrom = lsmesh.from_vertex_handle(hd).idx();
		int idto = lsmesh.to_vertex_handle(hd).idx();
		int idlarge = idto;
		int idsmall = idfrom;
		if (fvalues[idfrom] > fvalues[idto]) {
			idlarge = idfrom;
			idsmall = idto;
		}
		Eigen::Vector3d bdirec = (V.row(idsmall) - V.row(idlarge)).normalized();// the direction large to small
		Eigen::Vector3d norm = norm_f.row(fid);
		Eigen::Vector3d d1 = get_coff_vec_for_gradient(gradVF, fid, F(fid, 0));
		Eigen::Vector3d d2 = get_coff_vec_for_gradient(gradVF, fid, F(fid, 1));
		Eigen::Vector3d d3 = get_coff_vec_for_gradient(gradVF, fid, F(fid, 2));
		double c1 = norm.cross(-d1).dot(bdirec);
		double c2 = norm.cross(-d2).dot(bdirec);
		double c3 = norm.cross(-d3).dot(bdirec);
		triplets.push_back(Trip(i, F(fid, 0), c1));
		triplets.push_back(Trip(i, F(fid, 1), c2));
		triplets.push_back(Trip(i, F(fid, 2), c3));
		refids.push_back(F(fid, 0));
		refids.push_back(F(fid, 1));
		refids.push_back(F(fid, 2));
	}
	jacobian.setFromTriplets(triplets.begin(), triplets.end());
	DBdirections = jacobian;
}

void lsTools::calculate_boundary_direction_energy_function_values(const Eigen::MatrixXd& GradFValue,
	const std::vector<CGMesh::HalfedgeHandle>& edges, const Eigen::VectorXd& func, Eigen::VectorXd& lens, Eigen::VectorXd& energy) {
	int size = edges.size();
	energy.resize(size);
	lens.resize(size);
	for (int i = 0; i < size; i++) {
		CGMesh::HalfedgeHandle hd = edges[i];
		CGMesh::HalfedgeHandle op = lsmesh.opposite_halfedge_handle(hd);
		CGMesh::HalfedgeHandle next = lsmesh.next_halfedge_handle(hd);
		// get the face
		int fid = lsmesh.face_handle(hd).idx();
		if (fid < 0) {
			fid = lsmesh.face_handle(op).idx();
			next = lsmesh.next_halfedge_handle(op);
		}
		int idfrom = lsmesh.from_vertex_handle(hd).idx();
		int idto = lsmesh.to_vertex_handle(hd).idx();
		int idlarge = idto;
		int idsmall = idfrom;
		if (func[idfrom] > func[idto]) {
			idlarge = idfrom;
			idsmall = idto;
		}
		Eigen::Vector3d bdirec = (V.row(idsmall) - V.row(idlarge)).normalized();// the direction large to small
		Eigen::Vector3d norm = norm_f.row(fid);
		Eigen::Vector3d gradient = -GradFValue.row(fid);
		Eigen::Vector3d cross = norm.cross(gradient);
		lens[i] = cross.norm();

		double angle_radian = pseudo_geodesic_start_angle_degree * LSC_PI / 180.; // the angle in radian
		double value = cross.dot(bdirec) / lens[i] - cos(angle_radian);
		energy[i] = value;

		// check correctness
		CGMesh::VertexHandle vh1 = lsmesh.vertex_handle(idfrom);
		CGMesh::VertexHandle vh2 = lsmesh.vertex_handle(idto);

		if (!(lsmesh.is_boundary(vh1) && lsmesh.is_boundary(vh2)))
		{
			std::cout << "not boundary!" << std::endl;
		}

		if (gradient.dot(bdirec) < 0) {
			std::cout << "not correct gradient" << std::endl;
		}
		int vop = lsmesh.to_vertex_handle(next).idx();
		if (vop == idfrom || vop == idto) {
			std::cout << "wrong wrong" << std::endl;
		}
		Eigen::Vector3d rot = norm.cross(bdirec);
		Eigen::Vector3d other = V.row(vop) - V.row(idsmall);
		if (rot.dot(other) < 0) {
			std::cout << "the rot doesnot point to correct " << i << std::endl;
		}
		if (rot.dot(cross) < 0) {
			std::cout << "wrong cross!, " << i << " ang " << rot.dot(cross) << std::endl;
		}
	}
}
void lsTools::assemble_solver_fixed_boundary_direction_part(const Eigen::MatrixXd& GradFValue, const std::vector<CGMesh::HalfedgeHandle>& edges,
	const Eigen::VectorXd& func,
	spMat& H, Eigen::VectorXd& B, Eigen::VectorXd& energy) {
	Eigen::VectorXd lens;
	calculate_boundary_direction_energy_function_values(GradFValue, edges, func, lens, energy);
	spMat J = spMat(lens.asDiagonal().inverse()) * DBdirections;
	spMat JTJ = J.transpose() * J;
	B = -J.transpose() * BDEvalue;
	H = JTJ;
}
void lsTools::calculate_pseudo_energy_function_values_vertex_based(const double angle_degree, Eigen::VectorXd& lens)
{
	int ninner = IVids.size();
	VPEvalue.resize(ninner);
	LsOrient.resize(ninner);// first make g1 g2 same direction, then check if g1xg2 or g2xg1.
	PeWeight.resize(ninner);
	lens.resize(ninner);
	vBinormal.resize(ninner, 3);
	double angle_radian = angle_degree * LSC_PI / 180.; // the angle in radian
	double cos_angle = cos(angle_radian);

	for (int i = 0; i < ninner; i++)
	{
		if (ActInner[i] == false) {// this is a singularity
			continue;
		}
		int vid = IVids[i];
		CGMesh::HalfedgeHandle heh1 = Vheh0[i];
		CGMesh::HalfedgeHandle heh2 = Vheh1[i];
		int fid1, fid2;
		fid1 = lsmesh.face_handle(heh1).idx();
		fid2 = lsmesh.face_handle(heh2).idx();
		assert(fid1 != fid2);
		Eigen::Vector3d norm = norm_v.row(vid);
		Eigen::Vector3d norm1 = norm_f.row(fid1);
		Eigen::Vector3d norm2 = norm_f.row(fid2);
		Eigen::Vector3d g1 = gfvalue.row(fid1);
		Eigen::Vector3d g2 = gfvalue.row(fid2);
		// rotate gradients to get iso-curve directions
		g1 = norm1.cross(g1);
		g2 = norm2.cross(g2);
		Eigen::Vector3d g1xg2 = g1.cross(g2);
		// deal with the orientation of the level set: the left of the curve should be bigger value
		double flag1 = g1.dot(Vdire0[i]);
		double flag2 = g2.dot(Vdire1[i]);
		double flag = flag1 * flag2;
		if (flag >= 0) {
			LsOrient[i] = 1;
		}
		else {
			LsOrient[i] = -1;
		}

		// deal with PeWeight
		double lg1 = g1.norm();
		double lg2 = g2.norm();
		if (lg1 < 1e-16 || lg2 < 1e-16) {
			PeWeight[i] = 0;// it means the g1xg2 will not be accurate
		}
		else {
			Eigen::Vector3d g1n = g1.normalized();
			Eigen::Vector3d g2n = LsOrient[i] * g2.normalized();// get the oriented g2, to make g1 g2 goes to the same direction
			vBinormal.row(i) = g1n.cross(g2n).normalized();
			double cos_real = g1n.cross(g2n).normalized().dot(norm);// angle between the binormal and the surface normal
			if (cos_real < -1 || cos_real>1) {
				std::cout << "angle is wrong" << std::endl;
			}
			if (cos_angle < -1 || cos_angle>1) {
				std::cout << "angle is wrong1" << std::endl;
			}
			double cos_diff = fabs(cos_real - cos_angle) / 2;
			double angle_real = acos(g1n.dot(g2n)) * 180 / LSC_PI;// the angle between the two iso-lines
			// cos = 1, weight is 1; cos = -1, means it is very sharp turn, weight = 0
			// it is still resonable for vertex-based method. since there can be multiple intersections between the
			// osculating plane and the one-ring of vertex
			//std::cout<<"angle is "<<angle_real<<std::endl;
			// if the triangles are coplanar, or the segment is too straight, there is no room for optimization
			if (triangles_coplanar(V, F, fid1, fid2) || segments_colinear(g1n, g2n))
			{

				PeWeight[i] = 0;
			}
			else
			{
				PeWeight[i] = cos_diff;
			}
		}

		// if(fvalues[v1]<fvalues[v2]){ // orient g1xg2 or g2xg1
		//     LsOrient[i]*=-1;
		// }// we don't need it for vertex based method, since it is implicitly given by Vheh0 and Vheh1
		double length = g1xg2.norm();
		double value;
		if (length < 1e-16) {
			value = 1e6;// if it is really small, we set the value very big
		}
		else// divided by the length to avoid straight lines
		{
			value = LsOrient[i] * g1xg2.dot(norm) / length - cos_angle;
		}
		VPEvalue.coeffRef(i) = value;
		lens[i] = length;
	}
}
void lsTools::assemble_solver_boundary_condition_part(const Eigen::VectorXd& func, spMat& H, Eigen::VectorXd& B, Eigen::VectorXd &bcfvalue) {
	spMat bcJacobian;
	assert(trace_vers.size() == trace_hehs.size());
	int size = 0;// the number of constraints
	std::vector<double> vec_elements;
	std::vector<Trip> triplets;
	vec_elements.reserve(lsmesh.n_edges());
	triplets.reserve(lsmesh.n_edges());
	assert(assigned_trace_ls.size() == trace_vers.size());
	// assigned_trace_ls.resize(trace_vers.size());
	// for (int i = 0; i < trace_vers.size(); i++)// TODO assign assigned_trace_ls when tracing
	// {
	//     assigned_trace_ls[i]=i;
	// }
	// std::cout<<"before going into for loop"<<std::endl;
	for (int i = 0; i < trace_vers.size(); i++)
	{

		assert(trace_vers[i].size() == trace_hehs[i].size());
		for (int j = 0; j < trace_vers[i].size(); j++) {

			CGMesh::HalfedgeHandle edge = trace_hehs[i][j];
			Eigen::Vector3d point_middle = trace_vers[i][j];
			int id_from = lsmesh.from_vertex_handle(edge).idx();
			int id_to = lsmesh.to_vertex_handle(edge).idx();
			Eigen::Vector3d point_from = V.row(id_from);
			Eigen::Vector3d point_to = V.row(id_to);
			double d1 = (point_middle - point_from).norm();
			double d2 = (point_to - point_from).norm();
			triplets.push_back(Trip(size, id_from, d2 - d1));
			triplets.push_back(Trip(size, id_to, d1));
			double right_value = (d2 - d1) * func[id_from] + d1 * func[id_to] - d2 * assigned_trace_ls[i];
			vec_elements.push_back(right_value);
			size++;
		}
	}
	// std::cout<<"after loop"<<std::endl;
	// get boundary condition: the jacobian matrix
	bcJacobian.resize(size, V.rows());
	bcJacobian.setFromTriplets(triplets.begin(), triplets.end());// TODO calculate this before optimization
	bcfvalue.resize(size);
	// std::cout<<"calculated a and b"<<std::endl;
	for (int i = 0; i < size; i++) {// get the f(x)
		double value = vec_elements[i];
		bcfvalue(i) = value;
	}
	
	// get JTJ
	H = bcJacobian.transpose() * bcJacobian;
	// get the -(J^T)*f(x)
	B = -bcJacobian.transpose() * bcfvalue;
}
double get_mat_max_diag(spMat& M) {
	int nbr = M.rows();
	double value = 0;
	for (int i = 0; i < nbr; i++)
	{
		double dv = M.coeffRef(i, i);
		if (dv > value)
		{
			value = dv;
		}
	}
	return value;
}
void lsTools::assemble_solver_pesudo_geodesic_energy_part_vertex_based(const Eigen::VectorXd& vars, const std::vector<double>& angle_degree, Eigen::VectorXd& LocalActInner,
	std::vector<CGMesh::HalfedgeHandle>& heh0, std::vector<CGMesh::HalfedgeHandle>& heh1,
	std::vector<double>& t1s, std::vector<double>& t2s, const bool first_compute, const int aux_start_loc, spMat& H, Eigen::VectorXd& B, Eigen::VectorXd& energy)
{
	std::vector<Trip> tripletes;
	int vsize = V.rows();
	int ninner = LocalActInner.size();
	calculate_pseudo_geodesic_opt_expanded_function_values(vars, angle_degree,
		LocalActInner, heh0, heh1, t1s, t2s, first_compute, aux_start_loc, tripletes, energy);
	int nvars = vars.size();
	int ncondi = ninner * 4;
	spMat J;
	J.resize(ncondi, nvars);
	J.setFromTriplets(tripletes.begin(), tripletes.end());
	H = J.transpose() * J;
	B = -J.transpose() * energy;
}
// check if active faces exist around one vertex
// the information can be used to construct pseudo-geodesic energy


// before calling this function, please get the traced curves and assign values 
// for these curves
void lsTools::Run_Level_Set_Opt() {
	
	LSAnalizer ana;
	Eigen::MatrixXd GradValueF, GradValueV;
	std::array<spMat, 3> GradMatF;
	Eigen::VectorXd PGenergy;
	Eigen::VectorXd func = fvalues;
	std::vector<double> angle_degree;
	
	int vnbr = V.rows();
	int fnbr = F.rows();
	bool first_compute = false; // if we need initialize auxiliary vars
	analysis_pseudo_geodesic_on_vertices(func, ana);
	int ninner = ana.LocalActInner.size();
	int final_size = ninner * 3 + vnbr;// Change this when using more auxilary vars

	assert(assigned_trace_ls.size() && "Please trace the curves first before solve the energy");
	assert(assigned_trace_ls.size() == trace_vers.size());
	// initialize the level set with some number
	if (func.size() != vnbr)
	{
		initialize_level_set_accroding_to_parametrization();
		std::cout << "level set get initialized for smoothing" << std::endl;
		return;
	}
	if (Last_Opt_Mesh || func.size() != vnbr)
	{
		get_traced_boundary_triangle_direction_derivatives();
	}
	//  update quantities associated with level set values
	
	
	get_gradient_hessian_values(func, GradValueV, GradValueF);// TODO may not need this
	if (Glob_lsvars.size() == 0) {
		first_compute = true;
		std::cout << "Initializing Global Variable For LevelSet Opt ... " << std::endl;
		Glob_lsvars = Eigen::VectorXd::Zero(final_size);// We change the size if opt more than 1 level set
		Glob_lsvars.segment(0, vnbr) = func;
	}
	

	spMat H;
	H.resize(vnbr, vnbr);
	Eigen::VectorXd B=Eigen::VectorXd::Zero(vnbr);
	
	spMat LTL;  // left of laplacian
	Eigen::VectorXd mLTF; // right of laplacian
	assemble_solver_biharmonic_smoothing(func, LTL, mLTF);
	H = weight_laplacian * LTL;
	B = weight_laplacian * mLTF;
	assert(mass.rows() == vnbr);

	// fix inner vers and smooth boundary
	if (enable_inner_vers_fixed)
	{
		H += weight_mass * InnerV.asDiagonal();
	}
	Eigen::VectorXd bcfvalue;
	// boundary condition (traced as boundary condition)
	if (!enable_boundary_angles) {
		spMat bc_JTJ;
		Eigen::VectorXd bc_mJTF;
		assemble_solver_boundary_condition_part(func, bc_JTJ, bc_mJTF, bcfvalue);
		H += weight_boundary * bc_JTJ;
		B += weight_boundary * bc_mJTF;
	}

	// strip width condition
	if (enable_strip_width_energy)
	{
		spMat sw_JTJ;
		Eigen::VectorXd sw_mJTF;
		assemble_solver_strip_width_part(GradValueF, GradMatF, sw_JTJ, sw_mJTF);
		H += weight_strip_width * sw_JTJ;
		B += weight_strip_width * sw_mJTF;
	}

	// boundary edge angles 
	Eigen::VectorXd fbdenergy;
	if (enable_boundary_angles) {
		spMat JTJ;
		Eigen::VectorXd mJTF;
		assemble_solver_fixed_boundary_direction_part(GradValueF, tracing_start_edges, func, JTJ, mJTF, fbdenergy);
		H += weight_boundary * JTJ;
		B += weight_boundary * mJTF;
	}
	assert(H.rows() == vnbr);
	assert(H.cols() == vnbr);

	// pseudo geodesic
	spMat Hlarge;
	Eigen::VectorXd Blarge;
	Hlarge.resize(final_size, final_size);
	Blarge = Eigen::VectorXd::Zero(final_size);
	Eigen::VectorXd PGEnergy;
	if (enable_pseudo_geodesic_energy && weight_pseudo_geodesic_energy > 0)
	{
		spMat pg_JTJ;
		Efunc pg_mJTF;
		int aux_start_loc = vnbr;
		assemble_solver_pesudo_geodesic_energy_part_vertex_based(Glob_lsvars, angle_degree, ana.LocalActInner,
			ana.heh0, ana.heh1,
			ana.t1s, ana.t2s, first_compute, aux_start_loc, pg_JTJ, pg_mJTF, PGEnergy);
		Hlarge = sum_uneven_spMats(H, weight_pseudo_geodesic_energy * pg_JTJ);
		Blarge = sum_uneven_vectors(B, weight_pseudo_geodesic_energy * pg_mJTF);
		
	}
	
	Hlarge += 1e-6 * weight_mass * Eigen::VectorXd::Ones(vnbr).asDiagonal();
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver(Hlarge);
	// assert(solver.info() == Eigen::Success);
	if (solver.info() != Eigen::Success)
	{
		// solving failed
		std::cout << "solver fail" << std::endl;
		return;
	}
	// std::cout<<"solved successfully"<<std::endl;
	Eigen::VectorXd dx = solver.solve(Blarge).eval();
	dx *= 0.75;
	// std::cout << "step length " << dx.norm() << std::endl;
	double level_set_step_length = dx.norm();
	// double inf_norm=dx.cwiseAbs().maxCoeff();
	if (enable_boundary_angles || enable_pseudo_geodesic_energy)// only pg energy and boundary angles requires step length
	{
		if (level_set_step_length > max_step_length)
		{
			dx *= max_step_length / level_set_step_length;
		}
	}

	func += dx.topRows(vnbr);
	double energy_biharmonic = func.transpose() * QcH * func;
	double energy_boundary = (bcfvalue).norm();
	std::cout << "energy: harm " << energy_biharmonic << ", bnd " << energy_boundary << ", ";
	if (enable_pseudo_geodesic_energy)
	{
		double energy_pg = PGEnergy.norm();
		double max_energy_ls = PGEnergy.lpNorm<Eigen::Infinity>();
		std::cout << "pg, " << energy_pg << ", " << "AngleDiffMax," << max_energy_ls << ",";
		double max_ls_angle_energy = PGEnergy.bottomRows(ninner).norm();
		std::cout << "total angle energy, " << max_ls_angle_energy << ", ";
		// std::cout<<"PeWeight\n"<<PeWeight<<std::endl;
	}
	if (enable_strip_width_energy)
	{

		Eigen::VectorXd ener = GradValueF.rowwise().norm();
		ener = ener.asDiagonal() * ener;
		Eigen::VectorXd wds = Eigen::VectorXd::Ones(fnbr) * strip_width * strip_width;

		ener -= wds;
		double stp_energy = Eigen::VectorXd(ener).dot(ener);
		std::cout << "strip, " << stp_energy << ", ";
	}
	if (enable_boundary_angles) {
		double energy = fbdenergy.norm();
		std::cout << "bound_dirc, " << energy << ", ";

	}
	step_length = dx.norm();
	std::cout << "step " << step_length << std::endl;

	Last_Opt_Mesh = false;
}