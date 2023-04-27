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
void lsTools::analysis_pseudo_geodesic_on_vertices(const Eigen::VectorXd& func_values, Eigen::VectorXi& LocalActInner,
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
		assert(func_values[v4] > func_values[v3]);
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

// auxiliaries:
// r: the bi-normal. position: aux_start_loc ~ aux_start_loc + ninner * 3
// u: the side vector, position: aux_start_loc + ninner * 3 ~ aux_start_loc + ninner * 6
// s: Q-P （reducted without denominators, position: aux_start_loc + ninner * 6 ~ aux_start_loc + ninner * 9
// h: the r.dot(u). position: aux_start_loc + ninner * 9 ~ aux_start_loc + ninner * 10
// the nbr of vars: vnbr + ninner * 10
 
void lsTools::calculate_pseudo_geodesic_opt_expanded_function_values(Eigen::VectorXd& vars, const std::vector<double>& angle_degree,
	const LSAnalizer &analizer, const int vars_start_loc, const int aux_start_loc, std::vector<Trip>& tripletes, Eigen::VectorXd& Energy) {
	double cos_angle, sin_angle;
	int vnbr = V.rows();
	if(Binormals.rows()!=vnbr){
		Binormals = Eigen::MatrixXd::Zero(vnbr, 3);
	}
	int ninner = analizer.LocalActInner.size();
	if (angle_degree.size() == 1) {
		double angle_radian = angle_degree[0] * LSC_PI / 180.; // the angle in radian
		cos_angle = cos(angle_radian);
		sin_angle = sin(angle_radian);
	}
	if (Analyze_Optimized_LS_Angles)
	{
		// std::cout<<"Check 1"<<std::endl;
		AnalizedAngelVector.resize(ninner);
	}

	tripletes.clear();
	tripletes.reserve(ninner * 60); // 
	Energy = Eigen::VectorXd::Zero(ninner * 12); // mesh total energy values
#ifdef LSC_VARIABLE_CURVE_ANGLES
	std::vector<int> type_record = PG_Vertex_Types[0];
	for (int itr = 0; itr < type_record.size(); itr++)
	{
		int i = type_record[itr]; // the id in inner vers

#else
	assert(angle_degree.size() == 1 || angle_degree.size() == vnbr);
	for (int i = 0; i < ninner; i++)
	{

#endif
		int vm = IVids[i];
		if (analizer.LocalActInner[i] == false) {
			std::cout << "singularity, " << vm << std::endl;
			continue;
		}
#ifdef LSC_VARIABLE_CURVE_ANGLES
		double angle_radian = PGVariationalAngles[vm] * LSC_PI / 180.;
		cos_angle = cos(angle_radian);
		sin_angle = sin(angle_radian);
#else
		if (angle_degree.size() == vnbr)
		{
			double angle_radian = angle_degree[vm] * LSC_PI / 180.; // the angle in radian
			cos_angle = cos(angle_radian);
			sin_angle = sin(angle_radian);
		}
#endif

		CGMesh::HalfedgeHandle inhd = analizer.heh0[i], outhd = analizer.heh1[i];
		int v1 = lsmesh.from_vertex_handle(inhd).idx();
		int v2 = lsmesh.to_vertex_handle(inhd).idx();
		int v3 = lsmesh.from_vertex_handle(outhd).idx();
		int v4 = lsmesh.to_vertex_handle(outhd).idx();
		
		double t1 = analizer.t1s[i];
		double t2 = analizer.t2s[i];
		
		Eigen::Vector3d ver0 = V.row(v1) + (V.row(v2) - V.row(v1)) * t1;
		Eigen::Vector3d ver1 = V.row(vm);
		Eigen::Vector3d ver2 = V.row(v3) + (V.row(v4) - V.row(v3)) * t2;
		
		// the locations
		int lvm = vars_start_loc + vm;
		int lv1 = vars_start_loc + v1;
		int lv2 = vars_start_loc + v2;
		int lv3 = vars_start_loc + v3;
		int lv4 = vars_start_loc + v4;

		int lrx = i + aux_start_loc;
		int lry = i + aux_start_loc + ninner;
		int lrz = i + aux_start_loc + ninner * 2;
		int lux = i + aux_start_loc + ninner * 3;
		int luy = i + aux_start_loc + ninner * 4;
		int luz = i + aux_start_loc + ninner * 5;
		int lsx = i + aux_start_loc + ninner * 6;
		int lsy = i + aux_start_loc + ninner * 7;
		int lsz = i + aux_start_loc + ninner * 8;
		int lh = i + aux_start_loc + ninner * 9;
		Eigen::Vector3d norm = norm_v.row(vm);
		Eigen::Vector3d v31 = V.row(v3) - V.row(v1);
		Eigen::Vector3d v43 = V.row(v4) - V.row(v3);
		Eigen::Vector3d v21 = V.row(v2) - V.row(v1);
		double f1 = vars(lv1);
		double f2 = vars(lv2);
		double f3 = vars(lv3);
		double f4 = vars(lv4);
		double fm = vars(lvm);
		
		Eigen::Vector3d real_u = norm.cross(ver2 - ver0);
		real_u = real_u.normalized();
		if (Compute_Auxiliaries) {
			// std::cout<<"init auxiliaries"<<std::endl;
			Eigen::Vector3d real_s = v31 * (f4 - f3) * (f2 - f1) + v43 * (fm - f3) * (f2 - f1) - v21 * (fm - f1) * (f4 - f3);
			Eigen::Vector3d real_r = (ver1 - ver0).normalized().cross((ver2 - ver1).normalized());
			real_r = real_r.normalized();
			
			double real_h= real_r.dot(real_u);
			vars[lrx] = real_r[0];
			vars[lry] = real_r[1];
			vars[lrz] = real_r[2];
			vars[lux] = real_u[0];
			vars[luy] = real_u[1];
			vars[luz] = real_u[2];
			vars[lsx] = real_s[0];
			vars[lsy] = real_s[1];
			vars[lsz] = real_s[2];
			vars[lh] = real_h;
		}
		Eigen::Vector3d s = Eigen::Vector3d(vars(lsx), vars(lsy), vars(lsz));
		Eigen::Vector3d u = Eigen::Vector3d(vars(lux), vars(luy), vars(luz));
		double ns = s.norm();
		// r can flip to the opposite position, but u cannot, since it is the side vector.
		if (u.dot(real_u) < 0)
		{
			vars(lux) = real_u[0];
			vars(luy) = real_u[1];
			vars(luz) = real_u[2];
			u = real_u;
			std::cout<<"flip u"<<std::endl;
		}

		Eigen::Vector3d r = Eigen::Vector3d(vars[lrx], vars[lry], vars[lrz]);
		Binormals.row(vm) = r.dot(norm) < 0 ? -r : r;
		// the weights
		double dis0 = ((V.row(v1) - V.row(v2)) * vars[lvm] + (V.row(v2) - V.row(vm)) * vars[lv1] + (V.row(vm) - V.row(v1)) * vars[lv2]).norm();
		double dis1 = ((V.row(v3) - V.row(v4)) * vars[lvm] + (V.row(v4) - V.row(vm)) * vars[lv3] + (V.row(vm) - V.row(v3)) * vars[lv4]).norm();
		if (dis0 == 0 || dis1 == 0)
		{
			std::cout << "error, vars, " << vars[lvm] << ", " << vars[lv1] << ", " << vars[lv2] << ", " << vars[lv3] << ", " << vars[lv4] << std::endl;
		}
		assert(dis0 != 0 || dis1 != 0);
		
		double scale = mass_uniform.coeff(vm, vm);
		double EhanceBinormal = 1;
		if (Analyze_Optimized_LS_Angles && !Use_Fitting_Angles)
		{
			EhanceBinormal = 100;
		}
		assert(scale != 0);
		// r dot (vm+(t1-1)*vf-t1*vt)
		// vf = v1, vt = v2
		tripletes.push_back(Trip(i, lrx, ((V(v1, 0) - V(v2, 0)) * vars[lvm] + (V(v2, 0) - V(vm, 0)) * vars[lv1] + (V(vm, 0) - V(v1, 0)) * vars[lv2]) / dis0 * scale * EhanceBinormal));
		tripletes.push_back(Trip(i, lry, ((V(v1, 1) - V(v2, 1)) * vars[lvm] + (V(v2, 1) - V(vm, 1)) * vars[lv1] + (V(vm, 1) - V(v1, 1)) * vars[lv2]) / dis0 * scale * EhanceBinormal));
		tripletes.push_back(Trip(i, lrz, ((V(v1, 2) - V(v2, 2)) * vars[lvm] + (V(v2, 2) - V(vm, 2)) * vars[lv1] + (V(vm, 2) - V(v1, 2)) * vars[lv2]) / dis0 * scale * EhanceBinormal));

		double r12 = (V.row(v1) - V.row(v2)).dot(r);
		double rm1 = (V.row(vm) - V.row(v1)).dot(r);
		double r2m = (V.row(v2) - V.row(vm)).dot(r);
		tripletes.push_back(Trip(i, lvm, r12 / dis0 * scale * EhanceBinormal));
		tripletes.push_back(Trip(i, lv1, r2m / dis0 * scale * EhanceBinormal));
		tripletes.push_back(Trip(i, lv2, rm1 / dis0 * scale * EhanceBinormal));
		Energy[i] = (r12 * vars[lvm] + rm1 * vars[lv2] + r2m * vars[lv1]) / dis0 * scale * EhanceBinormal;

		// vf = v3, vt = v4
		tripletes.push_back(Trip(i + ninner, lrx, ((V(v3, 0) - V(v4, 0)) * vars[lvm] + (V(v4, 0) - V(vm, 0)) * vars[lv3] + (V(vm, 0) - V(v3, 0)) * vars[lv4]) / dis1 * scale * EhanceBinormal));
		tripletes.push_back(Trip(i + ninner, lry, ((V(v3, 1) - V(v4, 1)) * vars[lvm] + (V(v4, 1) - V(vm, 1)) * vars[lv3] + (V(vm, 1) - V(v3, 1)) * vars[lv4]) / dis1 * scale * EhanceBinormal));
		tripletes.push_back(Trip(i + ninner, lrz, ((V(v3, 2) - V(v4, 2)) * vars[lvm] + (V(v4, 2) - V(vm, 2)) * vars[lv3] + (V(vm, 2) - V(v3, 2)) * vars[lv4]) / dis1 * scale * EhanceBinormal));

		r12 = (V.row(v3) - V.row(v4)).dot(r);
		rm1 = (V.row(vm) - V.row(v3)).dot(r);
		r2m = (V.row(v4) - V.row(vm)).dot(r);
		tripletes.push_back(Trip(i + ninner, lvm, r12 / dis1 * scale * EhanceBinormal));
		tripletes.push_back(Trip(i + ninner, lv3, r2m / dis1 * scale * EhanceBinormal));
		tripletes.push_back(Trip(i + ninner, lv4, rm1 / dis1 * scale * EhanceBinormal));

		Energy[i + ninner] = (r12 * vars[lvm] + rm1 * vars[lv4] + r2m * vars[lv3]) / dis1 * scale * EhanceBinormal;

		// r*r=1
		tripletes.push_back(Trip(i + ninner * 2, lrx, 2 * vars(lrx) * scale * EhanceBinormal));
		tripletes.push_back(Trip(i + ninner * 2, lry, 2 * vars(lry) * scale * EhanceBinormal));
		tripletes.push_back(Trip(i + ninner * 2, lrz, 2 * vars(lrz) * scale * EhanceBinormal));

		Energy[i + ninner * 2] = (r.dot(r) - 1) * scale * EhanceBinormal;
		if (Analyze_Optimized_LS_Angles)
		{
			Eigen::Vector3d nu = u.normalized();
			Eigen::Vector3d nn = norm;
			Eigen::Vector3d nb = r.normalized();
			double cosina = nb.dot(nn);
			double sina = nb.dot(nu);
			if (sina < 0)
            {
                cosina *= -1;
            }
            double angle = acos(cosina);
			AnalizedAngelVector[i] = angle;
			if(Use_Opt_Only_BNMS){
				continue;
			}
		}
		// (r*norm)^2 - cos^2 = 0

		tripletes.push_back(Trip(i + ninner * 3, lrx,
								 (2 * norm(0) * norm(0) * vars(lrx) + 2 * norm(0) * norm(1) * vars(lry) + 2 * norm(0) * norm(2) * vars(lrz)) * scale));
		tripletes.push_back(Trip(i + ninner * 3, lry,
								 (2 * norm(1) * norm(1) * vars(lry) + 2 * norm(0) * norm(1) * vars(lrx) + 2 * norm(2) * norm(1) * vars(lrz)) * scale));
		tripletes.push_back(Trip(i + ninner * 3, lrz,
								 (2 * norm(2) * norm(2) * vars(lrz) + 2 * norm(0) * norm(2) * vars(lrx) + 2 * norm(1) * norm(2) * vars(lry)) * scale));

		double ndr = norm.dot(r);
		Energy[i + ninner * 3] = (ndr * ndr - cos_angle * cos_angle) * scale;

		// u*u=1
		
		tripletes.push_back(Trip(i + ninner * 4, lux, 2 * vars(lux) * scale));
		tripletes.push_back(Trip(i + ninner * 4, luy, 2 * vars(luy) * scale));
		tripletes.push_back(Trip(i + ninner * 4, luz, 2 * vars(luz) * scale));

		Energy[i + ninner * 4] = (u.dot(u) - 1) * scale;

		// u * norm = 0
		tripletes.push_back(Trip(i + ninner * 5, lux, norm(0) * scale));
		tripletes.push_back(Trip(i + ninner * 5, luy, norm(1) * scale));
		tripletes.push_back(Trip(i + ninner * 5, luz, norm(2) * scale));

		Energy[i + ninner * 5] = norm.dot(u) * scale;

		// u*s = 0

		tripletes.push_back(Trip(i + ninner * 6, lux, s(0) / ns * scale));
		tripletes.push_back(Trip(i + ninner * 6, luy, s(1) / ns * scale));
		tripletes.push_back(Trip(i + ninner * 6, luz, s(2) / ns * scale));

		tripletes.push_back(Trip(i + ninner * 6, lsx, u(0) / ns * scale));
		tripletes.push_back(Trip(i + ninner * 6, lsy, u(1) / ns * scale));
		tripletes.push_back(Trip(i + ninner * 6, lsz, u(2) / ns * scale));

		Energy[i + ninner * 6] = s.dot(u) / ns * scale;
		assert(ns != 0);
		// xxx-s = 0
		
		// x axis
		tripletes.push_back(Trip(i + ninner * 7, lsx, -1 * scale));
		tripletes.push_back(Trip(i + ninner * 7, lv1, (-v31(0) * (f4 - f3) - v43(0) * (fm - f3) + v21(0) * (f4 - f3)) * scale));
		tripletes.push_back(Trip(i + ninner * 7, lv2, (v31(0) * (f4 - f3) + v43(0) * (fm - f3)) * scale));
		tripletes.push_back(Trip(i + ninner * 7, lv3, (-v31(0) * (f2 - f1) - v43(0) * (f2 - f1) + v21(0) * (fm - f1)) * scale));
		tripletes.push_back(Trip(i + ninner * 7, lv4, (v31(0) * (f2 - f1) - v21(0) * (fm - f1)) * scale));
		tripletes.push_back(Trip(i + ninner * 7, lvm, (v43(0) * (f2 - f1) - v21(0) * (f4 - f3)) * scale));
		Energy[i + ninner * 7] = (v31(0) * (f4 - f3) * (f2 - f1) + v43(0) * (fm - f3) * (f2 - f1) - v21(0) * (fm - f1) * (f4 - f3) - s(0)) * scale;

		// y axis
		tripletes.push_back(Trip(i + ninner * 8, lsy, -1 * scale));
		tripletes.push_back(Trip(i + ninner * 8, lv1, (-v31(1) * (f4 - f3) - v43(1) * (fm - f3) + v21(1) * (f4 - f3)) * scale));
		tripletes.push_back(Trip(i + ninner * 8, lv2, (v31(1) * (f4 - f3) + v43(1) * (fm - f3)) * scale));
		tripletes.push_back(Trip(i + ninner * 8, lv3, (-v31(1) * (f2 - f1) - v43(1) * (f2 - f1) + v21(1) * (fm - f1)) * scale));
		tripletes.push_back(Trip(i + ninner * 8, lv4, (v31(1) * (f2 - f1) - v21(1) * (fm - f1)) * scale));
		tripletes.push_back(Trip(i + ninner * 8, lvm, (v43(1) * (f2 - f1) - v21(1) * (f4 - f3)) * scale));
		Energy[i + ninner * 8] = (v31(1) * (f4 - f3) * (f2 - f1) + v43(1) * (fm - f3) * (f2 - f1) - v21(1) * (fm - f1) * (f4 - f3) - s(1)) * scale;

		// z axis
		tripletes.push_back(Trip(i + ninner * 9, lsz, -1 * scale));
		tripletes.push_back(Trip(i + ninner * 9, lv1, (-v31(2) * (f4 - f3) - v43(2) * (fm - f3) + v21(2) * (f4 - f3)) * scale));
		tripletes.push_back(Trip(i + ninner * 9, lv2, (v31(2) * (f4 - f3) + v43(2) * (fm - f3)) * scale));
		tripletes.push_back(Trip(i + ninner * 9, lv3, (-v31(2) * (f2 - f1) - v43(2) * (f2 - f1) + v21(2) * (fm - f1)) * scale));
		tripletes.push_back(Trip(i + ninner * 9, lv4, (v31(2) * (f2 - f1) - v21(2) * (fm - f1)) * scale));
		tripletes.push_back(Trip(i + ninner * 9, lvm, (v43(2) * (f2 - f1) - v21(2) * (f4 - f3)) * scale));
		Energy[i + ninner * 9] = (v31(2) * (f4 - f3) * (f2 - f1) + v43(2) * (fm - f3) * (f2 - f1) - v21(2) * (fm - f1) * (f4 - f3) - s(2)) * scale;

		// r*u - h =0
		tripletes.push_back(Trip(i + ninner * 10, lrx, u(0) * scale));
		tripletes.push_back(Trip(i + ninner * 10, lry, u(1) * scale));
		tripletes.push_back(Trip(i + ninner * 10, lrz, u(2) * scale));

		tripletes.push_back(Trip(i + ninner * 10, lux, r(0) * scale));
		tripletes.push_back(Trip(i + ninner * 10, luy, r(1) * scale));
		tripletes.push_back(Trip(i + ninner * 10, luz, r(2) * scale));

		tripletes.push_back(Trip(i + ninner * 10, lh, -1 * scale));

		Energy[i + ninner * 10] = (r.dot(u) - vars[lh]) * scale;

		// h times (r.dot(n)) - sin*cos = 0
		double h = vars[lh];
		tripletes.push_back(Trip(i + ninner * 11, lh, r.dot(norm) * scale));
		tripletes.push_back(Trip(i + ninner * 11, lrx, norm(0) * h * scale));
		tripletes.push_back(Trip(i + ninner * 11, lry, norm(1) * h * scale));
		tripletes.push_back(Trip(i + ninner * 11, lrz, norm(2) * h * scale));

		Energy[i + ninner * 11] = (r.dot(norm) * h - sin_angle * cos_angle) * scale;

		//
		
	}
}
// convert angles (theta, phi) into a ray
Eigen::Vector3d angle_ray_converter(const double theta, const double phi)
{
	double theta_radian = theta * LSC_PI / 180.;
	double phi_radian = phi * LSC_PI / 180.;
	Eigen::Vector3d ray;
	ray[0] = cos(theta_radian) * sin(phi_radian);
	ray[1] = cos(theta_radian) * cos(phi_radian);
	ray[2] = sin(theta_radian);
	return ray;
}
// 
void angle_range_calculator(double theta, double phi, double theta_tol, double phi_tol, double &zmin, double &zmax, double &tanmin,
							double &tanmax)
{
	double target_theta = theta * LSC_PI / 180.;
	double target_phi = phi * LSC_PI / 180.;
	double target_theta_tol = theta_tol * LSC_PI / 180.;
	double target_phi_tol = phi_tol * LSC_PI / 180.;
	
	double theta_upper = target_theta + target_theta_tol;
	double theta_lower = target_theta - target_theta_tol;
	double phi_upper = target_phi + target_phi_tol;
	double phi_lower = target_phi - target_phi_tol;
	zmin = sin(theta_lower);
	zmax = sin(theta_upper);
	if(zmin>zmax){
		std::cout<<"Error in angle_range_calculator, zmin zmax"<<std::endl;
	}
	double tan1 = tan(phi_lower);
	double tan2 = tan(phi_upper);

	tanmin = std::min(tan1, tan2);
	tanmax = std::max(tan1, tan2);
}


// reformulate the problem as inequivalent problems without sin / cos
// the 14 auxiliaries: binormal r, variables for inequivalent: zl, zr, xl, xr, yr, the ray, the principle normal np.
void lsTools::calculate_shading_condition_inequivalent(Eigen::VectorXd &vars,
													   const LSAnalizer &analizer, const int vars_start_loc, const int aux_start_loc, std::vector<Trip> &tripletes, Eigen::VectorXd &Energy)
{
	int vnbr = V.rows();
	if (Binormals.rows() != vnbr)
	{
		Binormals = Eigen::MatrixXd::Zero(vnbr, 3);
	}
	int ninner = analizer.LocalActInner.size();

	tripletes.clear();
	tripletes.reserve(ninner * 60);				 //
	Energy = Eigen::VectorXd::Zero(ninner * 16); // mesh total energy values
	if (Lights.rows() != vnbr)
	{
		Lights = Eigen::MatrixXd::Zero(vnbr, 3);
	}
	double latitude_radian = ShadingLatitude * LSC_PI / 180.;
	double rho = (LSC_PI / 2 - latitude_radian);
	Eigen::Matrix3d rotation;
	rotation << 1, 0, 0,
		0, cos(rho), -sin(rho),
		0, sin(rho), cos(rho);
	Eigen::Vector3d direction_ground = Eigen::Vector3d(0, 0, -1); // the ground direction
	direction_ground = rotation * direction_ground;
	// std::cout<<"ground direction, "<<direction_ground.transpose()<<std::endl;
	for (int i = 0; i < ninner; i++)
	{
		if (analizer.LocalActInner[i] == false)
		{
			std::cout << "singularity" << std::endl;
			continue;
		}
		int vm = IVids[i];
		CGMesh::HalfedgeHandle inhd = analizer.heh0[i], outhd = analizer.heh1[i];
		int v1 = lsmesh.from_vertex_handle(inhd).idx();
		int v2 = lsmesh.to_vertex_handle(inhd).idx();
		int v3 = lsmesh.from_vertex_handle(outhd).idx();
		int v4 = lsmesh.to_vertex_handle(outhd).idx();

		double t1 = analizer.t1s[i];
		double t2 = analizer.t2s[i];

		Eigen::Vector3d ver0 = V.row(v1) + (V.row(v2) - V.row(v1)) * t1;
		Eigen::Vector3d ver1 = V.row(vm);
		Eigen::Vector3d ver2 = V.row(v3) + (V.row(v4) - V.row(v3)) * t2;

		// the locations
		int lvm = vars_start_loc + vm;
		int lv1 = vars_start_loc + v1;
		int lv2 = vars_start_loc + v2;
		int lv3 = vars_start_loc + v3;
		int lv4 = vars_start_loc + v4;

		int lrx = i + aux_start_loc;
		int lry = i + aux_start_loc + ninner;
		int lrz = i + aux_start_loc + ninner * 2;

		int lzl = i + aux_start_loc + ninner * 3;
		int lzr = i + aux_start_loc + ninner * 4;
		int lxl = i + aux_start_loc + ninner * 5;
		int lxr = i + aux_start_loc + ninner * 6;
		int lyr = i + aux_start_loc + ninner * 7;

		int lrayx = i + aux_start_loc + ninner * 8;
		int lrayy = i + aux_start_loc + ninner * 9;
		int lrayz = i + aux_start_loc + ninner * 10;
		// principle normals
		int lnpx = i + aux_start_loc + ninner * 11;
		int lnpy = i + aux_start_loc + ninner * 12;
		int lnpz = i + aux_start_loc + ninner * 13;
		assert(lnpz < vars.size());

		Eigen::Vector3d norm = norm_v.row(vm);
		Eigen::Vector3d v31 = V.row(v3) - V.row(v1);
		Eigen::Vector3d v43 = V.row(v4) - V.row(v3);
		Eigen::Vector3d v21 = V.row(v2) - V.row(v1);
		double f1 = vars(lv1);
		double f2 = vars(lv2);
		double f3 = vars(lv3);
		double f4 = vars(lv4);
		double fm = vars(lvm);
		double zmin, zmax, tanmin, tanmax;

		
		int condition_type = 1; // 0: through, 1: block, 2: transition, 3: reflection
		if(enable_reflection){
			condition_type = 3;
		}
		if(enable_let_ray_through){// if the surface totally let the ray through,
			condition_type = 0;
		}
		
		if (analizer.Special.size() > 0)
		{ // if we set multiple conditions on one surface
			// std::cout<<"special vers";
			if(!enable_let_ray_through){// partly shading partly transition.
				if (analizer.Special[i] == 1)
				{ // the unmarked vertices
					condition_type = 1;
				}
				else{
					condition_type = 2;
				}
			}
			else{ // shading, through and transition
				if (analizer.Special[i] == 0)
				{ // the unmarked vertices
					condition_type = 1;
				}
				if (analizer.Special[i] == 1)
				{ // the marked vertices
					condition_type = 0;
				}
				if (analizer.Special[i] == 2)
				{ // the transition part in between of marked/unmarked points
					condition_type = 2;
				}
			}
			
		}
		angle_range_calculator(Reference_theta, Reference_phi, Theta_tol, Phi_tol, zmin, zmax, tanmin, tanmax);
		if (analizer.Special.size() > 0 && condition_type == 0) // mix shading and through. set special condition for through
		{
			angle_range_calculator(Reference_theta1, Reference_phi1, Theta_tol1, Phi_tol1, zmin, zmax, tanmin, tanmax);
		}
		if (Compute_Auxiliaries || recompute_auxiliaries)
		{
			// std::cout<<"init auxiliaries"<<std::endl;
			// Eigen::Vector3d real_s = v31 * (f4 - f3) * (f2 - f1) + v43 * (fm - f3) * (f2 - f1) - v21 * (fm - f1) * (f4 - f3);
			Eigen::Vector3d real_r = (ver1 - ver0).normalized().cross((ver2 - ver1).normalized());

			real_r = real_r.normalized();

			// the binormal
			vars[lrx] = real_r[0];
			vars[lry] = real_r[1];
			vars[lrz] = real_r[2];

			// the ray and theta, phi
			Eigen::Vector3d real_ray = angle_ray_converter(Reference_theta, Reference_phi);
			if (analizer.Special.size() > 0 && condition_type == 0) // mix shading and through. set special condition for through
			{
				real_ray = angle_ray_converter(Reference_theta1, Reference_phi1);
			}
			double x = real_ray[0];
			double y = real_ray[1];
			double z = real_ray[2];
			double zl = sqrt(std::max(z - zmin, 0.0));
			double zr = sqrt(std::max(zmax - z, 0.0));
			double xl = sqrt(std::max(x / y - tanmin, 0.0));
			double xr = sqrt(std::max(tanmax - x / y, 0.0));
			double yr = sqrt(std::max(-y, 0.0));
			if (y > 0 && !recompute_auxiliaries)
			{
				std::cout << "Error, Please Check Here: calculate_shading_condition_inequivalent" << std::endl;
			}
			vars[lzl] = zl;
			vars[lzr] = zr;
			vars[lxl] = xl;
			vars[lxr] = xr;
			vars[lyr] = yr;

			vars[lrayx] = real_ray[0];
			vars[lrayy] = real_ray[1];
			vars[lrayz] = real_ray[2];

			Eigen::Vector3d real_np = real_r.cross(ver2 - ver0);
			real_np = real_np.normalized();
			vars[lnpx] = real_np[0];
			vars[lnpy] = real_np[1];
			vars[lnpz] = real_np[2];
		}
		

		// the weight of ray * tangent for the places who does not have solutions
		double weight_loose = 1;
		if (analizer.HighEnergy.size() == ninner)
		{
			if (analizer.HighEnergy[i] == false)
			{
				weight_loose = weight_geodesic;
			}
		}
		Eigen::Vector3d r = Eigen::Vector3d(vars[lrx], vars[lry], vars[lrz]);
		Eigen::Vector3d ray = Eigen::Vector3d(vars[lrayx], vars[lrayy], vars[lrayz]);
		Eigen::Vector3d np = Eigen::Vector3d(vars[lnpx], vars[lnpy], vars[lnpz]);
		if(ray[1]>0){ // flip y to make sure that y < 0
			ray[1] *= -1;
			vars[lrayy] *= -1;
		}

		double zl = vars[lzl];
		double zr = vars[lzr];
		double xl = vars[lxl];
		double xr = vars[lxr];
		double yr = vars[lyr];

		// std::cout<<"up_lower, "<<theta_upper<<", "<<theta_lower<<" "<<phi_upper<<" "<<phi_lower<<", ";
		Binormals.row(vm) = r.dot(norm) < 0 ? -r : r; // orient the binormal
		Lights.row(vm) = ray;
		// the weights
		double dis0 = ((V.row(v1) - V.row(v2)) * vars[lvm] + (V.row(v2) - V.row(vm)) * vars[lv1] + (V.row(vm) - V.row(v1)) * vars[lv2]).norm();
		double dis1 = ((V.row(v3) - V.row(v4)) * vars[lvm] + (V.row(v4) - V.row(vm)) * vars[lv3] + (V.row(vm) - V.row(v3)) * vars[lv4]).norm();
		double scale = mass_uniform.coeff(vm, vm);

		
		double weight_test = 1;
		if (condition_type == 1) // this part blocks light
		{
			// std::cout<<"blk,";
			// r * ray = 0
			tripletes.push_back(Trip(i + ninner * 12, lrx, weight_test * ray[0] * scale));
			tripletes.push_back(Trip(i + ninner * 12, lry, weight_test * ray[1] * scale));
			tripletes.push_back(Trip(i + ninner * 12, lrz, weight_test * ray[2] * scale));

			tripletes.push_back(Trip(i + ninner * 12, lrayx, weight_test * r[0] * scale));
			tripletes.push_back(Trip(i + ninner * 12, lrayy, weight_test * r[1] * scale));
			tripletes.push_back(Trip(i + ninner * 12, lrayz, weight_test * r[2] * scale));
			Energy[i + ninner * 12] = weight_test * r.dot(ray) * scale;

			// tangent * ray = 0
			Eigen::Vector3d v12 = V.row(v1) - V.row(v2);

			double d31 = ray.dot(v31);
			double d43 = ray.dot(v43);
			double d12 = ray.dot(v12);
			Eigen::Vector3d tangent = v31 * (f4 - f3) * (f2 - f1) + v43 * (fm - f3) * (f2 - f1) + v12 * (fm - f1) * (f4 - f3);
			double tnorm = tangent.norm();

			tripletes.push_back(Trip(i + ninner * 13, lv1, weight_loose * (-d31 * (f4 - f3) - d43 * (fm - f3) - d12 * (f4 - f3)) / tnorm * scale));
			tripletes.push_back(Trip(i + ninner * 13, lv2, weight_loose * (d31 * (f4 - f3) + d43 * (fm - f3)) / tnorm * scale));
			tripletes.push_back(Trip(i + ninner * 13, lv3, weight_loose * (-d31 * (f2 - f1) - d43 * (f2 - f1) - d12 * (fm - f1)) / tnorm * scale));
			tripletes.push_back(Trip(i + ninner * 13, lv4, weight_loose * (d31 * (f2 - f1) + d12 * (fm - f1)) / tnorm * scale));
			tripletes.push_back(Trip(i + ninner * 13, lvm, weight_loose * (d43 * (f2 - f1) + d12 * (f4 - f3)) / tnorm * scale));

			tripletes.push_back(Trip(i + ninner * 13, lrayx, weight_loose * tangent[0] / tnorm * scale));
			tripletes.push_back(Trip(i + ninner * 13, lrayy, weight_loose * tangent[1] / tnorm * scale));
			tripletes.push_back(Trip(i + ninner * 13, lrayz, weight_loose * tangent[2] / tnorm * scale));

			Energy[i + ninner * 13] = weight_loose * ray.dot(tangent) / tnorm * scale;
		}
		if (condition_type == 0) // let the light pass through
		{
			// std::cout<<"pass";

			// ray * np = 0
			tripletes.push_back(Trip(i + ninner * 12, lnpx, ray[0] * scale));
			tripletes.push_back(Trip(i + ninner * 12, lnpy, ray[1] * scale));
			tripletes.push_back(Trip(i + ninner * 12, lnpz, ray[2] * scale));

			tripletes.push_back(Trip(i + ninner * 12, lrayx, np[0] * scale));
			tripletes.push_back(Trip(i + ninner * 12, lrayy, np[1] * scale));
			tripletes.push_back(Trip(i + ninner * 12, lrayz, np[2] * scale));

			Energy[i + ninner * 12] = np.dot(ray) * scale;
		}
		if(condition_type == 3){// reflection
			
			// ray * np = direction_ground * np
			tripletes.push_back(Trip(i + ninner * 12, lrayx, np[0] * scale));
			tripletes.push_back(Trip(i + ninner * 12, lrayy, np[1] * scale));
			tripletes.push_back(Trip(i + ninner * 12, lrayz, np[2] * scale));

			tripletes.push_back(Trip(i + ninner * 12, lnpx, (ray[0] - direction_ground[0]) * scale));
			tripletes.push_back(Trip(i + ninner * 12, lnpy, (ray[1] - direction_ground[1]) * scale));
			tripletes.push_back(Trip(i + ninner * 12, lnpz, (ray[2] - direction_ground[2]) * scale));

			Energy[i + ninner * 12] = (ray.dot(np) - direction_ground.dot(np)) * scale;

			// direction_ground x ray .dot(np)=0
			Eigen::Vector3d gxray = direction_ground.cross(ray);
			tripletes.push_back(Trip(i + ninner * 13, lnpx, gxray[0] * scale));
			tripletes.push_back(Trip(i + ninner * 13, lnpy, gxray[1] * scale));
			tripletes.push_back(Trip(i + ninner * 13, lnpz, gxray[2] * scale));

			double cray0 = np[1] * direction_ground[2] - np[2] * direction_ground[1];
			double cray1 = -np[0] * direction_ground[2] + np[2] * direction_ground[0];
			double cray2 = np[0] * direction_ground[1] - np[1] * direction_ground[0];
			tripletes.push_back(Trip(i + ninner * 13, lrayx, (cray0) * scale));
			tripletes.push_back(Trip(i + ninner * 13, lrayy, (cray1) * scale));
			tripletes.push_back(Trip(i + ninner * 13, lrayz, (cray2) * scale));

			Energy[i + ninner * 13] = gxray.dot(np) * scale;
		}
		if(condition_type == 2){ // if this is a transition point, we do not apply ANY condition, but let the fairness term do its work
			// std::cout<<"skip transition, ";
			continue;
		}
		// r dot (vm+(t1-1)*vf-t1*vt)
		// vf = v1, vt = v2
		tripletes.push_back(Trip(i, lrx, ((V(v1, 0) - V(v2, 0)) * vars[lvm] + (V(v2, 0) - V(vm, 0)) * vars[lv1] + (V(vm, 0) - V(v1, 0)) * vars[lv2]) / dis0 * scale * weight_binormal));
		tripletes.push_back(Trip(i, lry, ((V(v1, 1) - V(v2, 1)) * vars[lvm] + (V(v2, 1) - V(vm, 1)) * vars[lv1] + (V(vm, 1) - V(v1, 1)) * vars[lv2]) / dis0 * scale * weight_binormal));
		tripletes.push_back(Trip(i, lrz, ((V(v1, 2) - V(v2, 2)) * vars[lvm] + (V(v2, 2) - V(vm, 2)) * vars[lv1] + (V(vm, 2) - V(v1, 2)) * vars[lv2]) / dis0 * scale * weight_binormal));

		double r12 = (V.row(v1) - V.row(v2)).dot(r);
		double rm1 = (V.row(vm) - V.row(v1)).dot(r);
		double r2m = (V.row(v2) - V.row(vm)).dot(r);
		tripletes.push_back(Trip(i, lvm, r12 / dis0 * scale * weight_binormal));
		tripletes.push_back(Trip(i, lv1, r2m / dis0 * scale * weight_binormal));
		tripletes.push_back(Trip(i, lv2, rm1 / dis0 * scale * weight_binormal));
		Energy[i] = (r12 * vars[lvm] + rm1 * vars[lv2] + r2m * vars[lv1]) / dis0 * scale * weight_binormal;

		// vf = v3, vt = v4
		tripletes.push_back(Trip(i + ninner, lrx, ((V(v3, 0) - V(v4, 0)) * vars[lvm] + (V(v4, 0) - V(vm, 0)) * vars[lv3] + (V(vm, 0) - V(v3, 0)) * vars[lv4]) / dis1 * scale * weight_binormal));
		tripletes.push_back(Trip(i + ninner, lry, ((V(v3, 1) - V(v4, 1)) * vars[lvm] + (V(v4, 1) - V(vm, 1)) * vars[lv3] + (V(vm, 1) - V(v3, 1)) * vars[lv4]) / dis1 * scale * weight_binormal));
		tripletes.push_back(Trip(i + ninner, lrz, ((V(v3, 2) - V(v4, 2)) * vars[lvm] + (V(v4, 2) - V(vm, 2)) * vars[lv3] + (V(vm, 2) - V(v3, 2)) * vars[lv4]) / dis1 * scale * weight_binormal));

		r12 = (V.row(v3) - V.row(v4)).dot(r);
		rm1 = (V.row(vm) - V.row(v3)).dot(r);
		r2m = (V.row(v4) - V.row(vm)).dot(r);
		tripletes.push_back(Trip(i + ninner, lvm, r12 / dis1 * scale * weight_binormal));
		tripletes.push_back(Trip(i + ninner, lv3, r2m / dis1 * scale * weight_binormal));
		tripletes.push_back(Trip(i + ninner, lv4, rm1 / dis1 * scale * weight_binormal));

		Energy[i + ninner] = (r12 * vars[lvm] + rm1 * vars[lv4] + r2m * vars[lv3]) / dis1 * scale * weight_binormal;

		// r*r=1
		tripletes.push_back(Trip(i + ninner * 2, lrx, 2 * vars(lrx) * scale * weight_binormal));
		tripletes.push_back(Trip(i + ninner * 2, lry, 2 * vars(lry) * scale * weight_binormal));
		tripletes.push_back(Trip(i + ninner * 2, lrz, 2 * vars(lrz) * scale * weight_binormal));

		Energy[i + ninner * 2] = (r.dot(r) - 1) * scale * weight_binormal;
		// std::cout<<"c3 got"<<std::endl;
		// if (analizer.ShadSpecial.size() != ninner)
		// {
		// 	std::cout << "Please check the surface patches that are parallel to the light" << std::endl;
		// }
		// if (analizer.ShadSpecial[i] == 0 && !enable_let_ray_through && !enable_reflection) // if this patch is parallel to the light, skip
		// {
		// 	continue;
		// }
		
		
		// np * r = 0
		tripletes.push_back(Trip(i + ninner * 3, lrx, np[0] * scale));
		tripletes.push_back(Trip(i + ninner * 3, lry, np[1] * scale));
		tripletes.push_back(Trip(i + ninner * 3, lrz, np[2] * scale));

		tripletes.push_back(Trip(i + ninner * 3, lnpx, r[0] * scale));
		tripletes.push_back(Trip(i + ninner * 3, lnpy, r[1] * scale));
		tripletes.push_back(Trip(i + ninner * 3, lnpz, r[2] * scale));

		Energy[i + ninner * 3] = r.dot(np) * scale;

		// tangent * np = 0
		Eigen::Vector3d v12 = V.row(v1) - V.row(v2);
		double d31 = np.dot(v31);
		double d43 = np.dot(v43);
		double d12 = np.dot(v12);
		Eigen::Vector3d tangent = v31 * (f4 - f3) * (f2 - f1) + v43 * (fm - f3) * (f2 - f1) + v12 * (fm - f1) * (f4 - f3);
		double tnorm = tangent.norm();

		tripletes.push_back(Trip(i + ninner * 4, lv1, 1 * (-d31 * (f4 - f3) - d43 * (fm - f3) - d12 * (f4 - f3)) / tnorm * scale));
		tripletes.push_back(Trip(i + ninner * 4, lv2, 1 * (d31 * (f4 - f3) + d43 * (fm - f3)) / tnorm * scale));
		tripletes.push_back(Trip(i + ninner * 4, lv3, 1 * (-d31 * (f2 - f1) - d43 * (f2 - f1) - d12 * (fm - f1)) / tnorm * scale));
		tripletes.push_back(Trip(i + ninner * 4, lv4, 1 * (d31 * (f2 - f1) + d12 * (fm - f1)) / tnorm * scale));
		tripletes.push_back(Trip(i + ninner * 4, lvm, 1 * (d43 * (f2 - f1) + d12 * (f4 - f3)) / tnorm * scale));

		tripletes.push_back(Trip(i + ninner * 4, lnpx, 1 * tangent[0] / tnorm * scale));
		tripletes.push_back(Trip(i + ninner * 4, lnpy, 1 * tangent[1] / tnorm * scale));
		tripletes.push_back(Trip(i + ninner * 4, lnpz, 1 * tangent[2] / tnorm * scale));

		Energy[i + ninner * 4] = 1 * np.dot(tangent) / tnorm * scale;

		// np * np = 1
		tripletes.push_back(Trip(i + ninner * 11, lnpx, 2 * np[0] * scale));
		tripletes.push_back(Trip(i + ninner * 11, lnpy, 2 * np[1] * scale));
		tripletes.push_back(Trip(i + ninner * 11, lnpz, 2 * np[2] * scale));

		Energy[i + ninner * 11] = (np.dot(np) - 1) * scale;

		
		// z - zmin - zl^2 = 0
		tripletes.push_back(Trip(i + ninner * 5, lrayz, 1 * scale));
		tripletes.push_back(Trip(i + ninner * 5, lzl, (-2 * zl) * scale));

		Energy[i + ninner * 5] = (ray[2] - zmin - zl * zl) * scale;

		// zmax - z - zr^2 = 0
		tripletes.push_back(Trip(i + ninner * 6, lrayz, -1 * scale));
		tripletes.push_back(Trip(i + ninner * 6, lzr, (-2 * zr) * scale));

		Energy[i + ninner * 6] = (zmax - ray[2] - zr * zr) * scale;
		
		// x/y >= tanmin -> x - tanmin * y + xl^2 = 0, since y < 0
		tripletes.push_back(Trip(i + ninner * 7, lrayx, 1 * scale));
		tripletes.push_back(Trip(i + ninner * 7, lrayy, -tanmin * scale));
		tripletes.push_back(Trip(i + ninner * 7, lxl, (2 * xl) * scale));

		Energy[i + ninner * 7] = (ray[0] - tanmin * ray[1] + xl * xl) * scale;
		
		// x/y <= tanmax -> x - tanmax * y - xr^2 = 0, since y < 0
		tripletes.push_back(Trip(i + ninner * 8, lrayx, 1 * scale));
		tripletes.push_back(Trip(i + ninner * 8, lrayy, -tanmax * scale));
		tripletes.push_back(Trip(i + ninner * 8, lxr, (-2 * xr) * scale));

		Energy[i + ninner * 8] = (ray[0] - tanmax * ray[1] - xr * xr) * scale;

		// y <= 0 -> y + yr*yr = 0
		tripletes.push_back(Trip(i + ninner * 9, lrayy, 1 * scale));
		tripletes.push_back(Trip(i + ninner * 9, lyr, 2 * yr * scale));

		Energy[i + ninner * 9] = (ray[1] + yr * yr) * scale;

		// ray * ray =1
		tripletes.push_back(Trip(i + ninner * 10, lrayx, 2 * ray[0] * scale)); 
		tripletes.push_back(Trip(i + ninner * 10, lrayy, 2 * ray[1] * scale)); 
		tripletes.push_back(Trip(i + ninner * 10, lrayz, 2 * ray[2] * scale));

		Energy[i + ninner * 10] = (ray.dot(ray) - 1) * scale;

	}
	double max_e = 0;
	int max_loc = -1;
	for (int i = 0; i < Energy.size(); i++)
	{
		if (abs(Energy[i]) > max_e)
		{
			max_e = abs(Energy[i]);
			max_loc = i;
		}
	}
	int howmany = max_loc / ninner;
	std::cout << "me, " << max_e << ", mloc, " << howmany << ",";
}

// the principle normal on a circle.
void lsTools::calculate_shading_init(Eigen::VectorXd& vars,
	const LSAnalizer &analizer, const int vars_start_loc, const int aux_start_loc, std::vector<Trip>& tripletes, Eigen::VectorXd& Energy) {
	int vnbr = V.rows();
	int ninner = analizer.LocalActInner.size();
	double latitude_radian = ShadingLatitude * LSC_PI / 180.;
	// this is to give the correct sun height (in angle) of the architecture.
	double rho = (LSC_PI / 2 - latitude_radian);
	Eigen::Matrix3d rotation;
	rotation << 1, 0, 0,
		0, cos(rho), -sin(rho),
		0, sin(rho), cos(rho);
	// initialize the orthogonal vectors on the faces.
	if (Fslope.size() == 0)
	{

		int whichtype = 0; // shading
		if (analizer.Special.size() > 0)
		{
			whichtype = 2; // partly shading
		}
		if (enable_let_ray_through)
		{
			whichtype = 1; // through
			if (analizer.Special.size() > 0)
			{
				whichtype = 3; // mixing shading and through
			}
		}
		orthogonal_slope_for_different_shading_types(whichtype, InnerV, analizer.Special,
													 V, F, norm_f, Reference_theta, Reference_phi, Theta_tol, Phi_tol,
													 Reference_theta1, Reference_phi1,
													 Fslope, OrthoSlope);
		std::cout<<"####3Shading Type: "<<whichtype<<"\n";
		int count1 = 0, count2 = 0;
		for(int i=0;i<analizer.Special.size();i++){
			if(analizer.Special[i]==1){
				count1++;
			}
			else{
				count2++;
			}
		}
		std::cout<<"Special ver nbr "<<count1<<" out of "<<count2<<" vertices\n";
	}
	tripletes.clear();
	tripletes.reserve(Fslope.size() * 3); //todo;
	Energy = Eigen::VectorXd::Zero(Fslope.size());// todo;
	for (int itr = 0; itr < Fslope.size(); itr++)
	{
		int fid = Fslope[itr];
		Eigen::Vector3d otho = OrthoSlope[itr];
		int v0 = F(fid, 0);
		int v1 = F(fid, 1);
		int v2 = F(fid, 2);
		// locations
		int loc0 = vars_start_loc + v0;
		int loc1 = vars_start_loc + v1;
		int loc2 = vars_start_loc + v2;
		
		double f0 = Glob_lsvars[loc0];
		double f1 = Glob_lsvars[loc1];
		double f2 = Glob_lsvars[loc2];
		Eigen::Vector3d ver0 = V.row(v0);
		Eigen::Vector3d ver1 = V.row(v1);
		Eigen::Vector3d ver2 = V.row(v2);

		Eigen::Vector3d direction = f0 * (ver1 - ver2) + f2 * (ver0 - ver1) + f1 * (ver2 - ver0);
		double scale = direction.norm();
		// direction * otho / scale = 0;
		tripletes.push_back(Trip(itr, loc0, (ver1 - ver2).dot(otho) / scale));
		tripletes.push_back(Trip(itr, loc1, (ver2 - ver0).dot(otho) / scale));
		tripletes.push_back(Trip(itr, loc2, (ver0 - ver1).dot(otho) / scale));

		Energy[itr] = direction.dot(otho) / scale;
	}
}
// deal with asymptotic and geodesic, for using fewer auxiliary variables
void lsTools::calculate_extreme_pseudo_geodesic_values(Eigen::VectorXd &vars, const bool asymptotic,
													   const LSAnalizer &analizer, const int vars_start_loc, std::vector<Trip> &tripletes, Eigen::VectorXd &Energy)
{
	int vnbr = V.rows();
	if(Binormals.rows()!=vnbr){
		Binormals = Eigen::MatrixXd::Zero(vnbr, 3);
	}
	int ninner = analizer.LocalActInner.size();
	tripletes.clear();
	tripletes.reserve(ninner * 15);				// the number of rows is ninner*4, the number of cols is aux_start_loc + ninner * 3 (all the function values and auxiliary vars)
	Energy = Eigen::VectorXd::Zero(ninner * 2); // mesh total energy values
	double max_ele = 0;
	double max_eng = 0;
#ifdef LSC_VARIABLE_CURVE_ANGLES
	std::vector<int> type_record =  PG_Vertex_Types[1];
	for (int itr = 0; itr < type_record.size(); itr++)
	{
		int i = type_record[itr]; // the id in inner vers

#else
	for (int i = 0; i < ninner; i++)
	{
		
		
#endif
		if (analizer.LocalActInner[i] == false)
		{
			std::cout << "singularity" << std::endl;
			continue;
		}
		int vm = IVids[i];
		
		CGMesh::HalfedgeHandle inhd = analizer.heh0[i], outhd = analizer.heh1[i];
		int v1 = lsmesh.from_vertex_handle(inhd).idx();
		int v2 = lsmesh.to_vertex_handle(inhd).idx();
		int v3 = lsmesh.from_vertex_handle(outhd).idx();
		int v4 = lsmesh.to_vertex_handle(outhd).idx();

		double t1 = analizer.t1s[i];
		double t2 = analizer.t2s[i];
		
		Eigen::Vector3d ver0 = V.row(v1) + (V.row(v2) - V.row(v1)) * t1;
		Eigen::Vector3d ver1 = V.row(vm);
		Eigen::Vector3d ver2 = V.row(v3) + (V.row(v4) - V.row(v3)) * t2;
		// the locations
		int lvm = vars_start_loc + vm;
		int lv1 = vars_start_loc + v1;
		int lv2 = vars_start_loc + v2;
		int lv3 = vars_start_loc + v3;
		int lv4 = vars_start_loc + v4;

		// the weights
		double dis0 = ((V.row(v1) - V.row(v2)) * vars[lvm] + (V.row(v2) - V.row(vm)) * vars[lv1] + (V.row(vm) - V.row(v1)) * vars[lv2]).norm();
		double dis1 = ((V.row(v3) - V.row(v4)) * vars[lvm] + (V.row(v4) - V.row(vm)) * vars[lv3] + (V.row(vm) - V.row(v3)) * vars[lv4]).norm();
		double scale = mass_uniform.coeff(vm, vm);
		Eigen::Vector3d norm = norm_v.row(vm);
		if (asymptotic)//asymptotic conditions: n * d1 = 0, n * d2 = 0
		{
			Eigen::Vector3d norm = norm_v.row(vm);
			double r12 = (V.row(v1) - V.row(v2)).dot(norm);
			double rm1 = (V.row(vm) - V.row(v1)).dot(norm);
			double r2m = (V.row(v2) - V.row(vm)).dot(norm);

			tripletes.push_back(Trip(i, lvm, r12 / dis0 * scale));
			tripletes.push_back(Trip(i, lv1, r2m / dis0 * scale));
			tripletes.push_back(Trip(i, lv2, rm1 / dis0 * scale));
			Energy[i] = (r12 * vars[lvm] + rm1 * vars[lv2] + r2m * vars[lv1]) / dis0 * scale;
			if (fabs(r12 / dis0) > max_ele)
			{
				max_ele = fabs(r12 / dis0);
			}
			if (fabs(r2m / dis0) > max_ele)
			{
				max_ele = fabs(r2m / dis0);
			}
			if (fabs(rm1 / dis0) > max_ele)
			{
				max_ele = fabs(rm1 / dis0);
			}
			r12 = (V.row(v3) - V.row(v4)).dot(norm);
			rm1 = (V.row(vm) - V.row(v3)).dot(norm);
			r2m = (V.row(v4) - V.row(vm)).dot(norm);
			tripletes.push_back(Trip(i + ninner, lvm, r12 / dis1 * scale));
			tripletes.push_back(Trip(i + ninner, lv3, r2m / dis1 * scale));
			tripletes.push_back(Trip(i + ninner, lv4, rm1 / dis1 * scale));
			Energy[i + ninner] = (r12 * vars[lvm] + rm1 * vars[lv4] + r2m * vars[lv3]) / dis1 * scale;
			if (fabs(r12 / dis1) > max_ele)
			{
				max_ele = fabs(r12 / dis1);
			}
			if (fabs(r2m / dis1) > max_ele)
			{
				max_ele = fabs(r2m / dis1);
			}
			if (fabs(rm1 / dis1) > max_ele)
			{
				max_ele = fabs(rm1 / dis1);
			}
			
		}
		else{// geodesic condition: norm, d1, d2 coplanar
			
			// if(use_given_direction){
			// 	norm = ray.normalized();
			// 	// std::cout<<"correct given direction "<<norm.transpose()<<std::endl;
			// }
			Eigen::Vector3d r12 = norm.cross(Eigen::Vector3d(V.row(v1) - V.row(v2)));
			Eigen::Vector3d r2m = norm.cross(Eigen::Vector3d(V.row(v2) - V.row(vm)));
			Eigen::Vector3d rm1 = norm.cross(Eigen::Vector3d(V.row(vm) - V.row(v1)));
			Eigen::Vector3d v34 = V.row(v3) - V.row(v4);
			Eigen::Vector3d v4m = V.row(v4) - V.row(vm);
			Eigen::Vector3d vm3 = V.row(vm) - V.row(v3);
			double s1 = dis0 * dis1 / scale;
			double fm = vars[lvm], f1 = vars[lv1], f2 = vars[lv2], f3 = vars[lv3], f4 = vars[lv4];
			Eigen::Vector3d vec_l = r12 * fm + r2m * f1 + rm1 * f2;
			Eigen::Vector3d vec_r = v34 * fm + v4m * f3 + vm3 * f4;
			tripletes.push_back(Trip(i, lvm, (r12.dot(vec_r) + v34.dot(vec_l)) / s1));
			tripletes.push_back(Trip(i, lv1, r2m.dot(vec_r) / s1));
			tripletes.push_back(Trip(i, lv2, rm1.dot(vec_r) / s1));
			tripletes.push_back(Trip(i, lv3, v4m.dot(vec_l) / s1));
			tripletes.push_back(Trip(i, lv4, vm3.dot(vec_l) / s1));
			Energy[i] = vec_l.dot(vec_r) / s1;
			
		}
		Eigen::Vector3d cross = (ver1 - ver0).cross(ver2 - ver1);
		Binormals.row(vm) = cross.normalized();
		double eng = cross.normalized().dot(norm);
		eng = abs(eng);
		if (eng > max_eng)
		{
			max_eng = eng;
		}
	}
	// //std::cout<<"max_ele, "<<max_ele<<", ";
	// if(!asymptotic){
	// 	std::cout<<"max_eng, "<<max_eng<<", ";
	// }
	
}

// this function only need be called after initializing the level set
void lsTools::get_traced_boundary_triangle_direction_derivatives() {
	int size = tracing_start_edges.size();

	spMat jacobian;
	jacobian.resize(size, V.rows());
	std::vector<Trip> triplets;
	triplets.reserve(size * 3);
	// refids.clear();
	// refids.reserve(size * 3);
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
		// refids.push_back(F(fid, 0));
		// refids.push_back(F(fid, 1));
		// refids.push_back(F(fid, 2));
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

void lsTools::calculate_binormal_regulizer(Eigen::VectorXd& vars,
	const LSAnalizer &analizer, const int vars_start_loc, const int aux_start_loc, std::vector<Trip> &tripletes, Eigen::VectorXd& Energy) {
	int vnbr = V.rows();
	int ninner = analizer.LocalActInner.size();

	tripletes.clear();
	tripletes.reserve(ninner * 12);				//
	Energy = Eigen::VectorXd::Zero(ninner * 4); // mesh total energy values
	for (int i = 0; i < ninner; i++)
	{
		if (analizer.LocalActInner[i] == false) {
			std::cout << "singularity" << std::endl;
			continue;
		}
		int vm = IVids[i];
		CGMesh::HalfedgeHandle inhd = analizer.heh0[i], outhd = analizer.heh1[i];
		int v1 = lsmesh.from_vertex_handle(inhd).idx();
		int v2 = lsmesh.to_vertex_handle(inhd).idx();
		int v3 = lsmesh.from_vertex_handle(outhd).idx();
		int v4 = lsmesh.to_vertex_handle(outhd).idx();
		
		double t1 = analizer.t1s[i];
		double t2 = analizer.t2s[i];
		
		Eigen::Vector3d ver0 = V.row(v1) + (V.row(v2) - V.row(v1)) * t1;
		Eigen::Vector3d ver1 = V.row(vm);
		Eigen::Vector3d ver2 = V.row(v3) + (V.row(v4) - V.row(v3)) * t2;
		
		// the locations
		int lvm = vars_start_loc + vm;
		int lv1 = vars_start_loc + v1;
		int lv2 = vars_start_loc + v2;
		int lv3 = vars_start_loc + v3;
		int lv4 = vars_start_loc + v4;

		int lrx = i + aux_start_loc;
		int lry = i + aux_start_loc + ninner;
		int lrz = i + aux_start_loc + ninner * 2;

		Eigen::Vector3d norm = norm_v.row(vm);
		Eigen::Vector3d v31 = V.row(v3) - V.row(v1);
		Eigen::Vector3d v43 = V.row(v4) - V.row(v3);
		Eigen::Vector3d v21 = V.row(v2) - V.row(v1);
		double f1 = vars(lv1);
		double f2 = vars(lv2);
		double f3 = vars(lv3);
		double f4 = vars(lv4);
		double fm = vars(lvm);

		int r1inner = InnerV[v1];
		int r2inner = InnerV[v2];
		int r3inner = InnerV[v3];
		int r4inner = InnerV[v4];
		
		int l1x = r1inner + aux_start_loc;
		int l1y = r1inner + aux_start_loc + ninner;
		int l1z = r1inner + aux_start_loc + ninner * 2;

		int l2x = r2inner + aux_start_loc;
		int l2y = r2inner + aux_start_loc + ninner;
		int l2z = r2inner + aux_start_loc + ninner * 2;

		int l3x = r3inner + aux_start_loc;
		int l3y = r3inner + aux_start_loc + ninner;
		int l3z = r3inner + aux_start_loc + ninner * 2;

		int l4x = r4inner + aux_start_loc;
		int l4y = r4inner + aux_start_loc + ninner;
		int l4z = r4inner + aux_start_loc + ninner * 2;

		Eigen::Vector3d r = Eigen::Vector3d(vars[lrx], vars[lry], vars[lrz]);
		Eigen::Vector3d r1 = Eigen::Vector3d(vars[l1x], vars[l1y], vars[l1z]);
		Eigen::Vector3d r2 = Eigen::Vector3d(vars[l2x], vars[l2y], vars[l2z]);
		Eigen::Vector3d r3 = Eigen::Vector3d(vars[l3x], vars[l3y], vars[l3z]);
		Eigen::Vector3d r4 = Eigen::Vector3d(vars[l4x], vars[l4y], vars[l4z]);
		
		// r * r1 +-1 = 0
		double sign = r.dot(r1) > 0 ? 1 : -1;
		tripletes.push_back(Trip(i, lrx, r1[0]));
		tripletes.push_back(Trip(i, lry, r1[1]));
		tripletes.push_back(Trip(i, lrz, r1[2]));

		tripletes.push_back(Trip(i, l1x, r[0]));
		tripletes.push_back(Trip(i, l1y, r[1]));
		tripletes.push_back(Trip(i, l1z, r[2]));

		Energy[i] = r.dot(r1) - sign;

		// r * r2 +- 1 = 0
		sign = r.dot(r2) > 0 ? 1 : -1;
		tripletes.push_back(Trip(i + ninner, lrx, r2[0]));
		tripletes.push_back(Trip(i + ninner, lry, r2[1]));
		tripletes.push_back(Trip(i + ninner, lrz, r2[2]));

		tripletes.push_back(Trip(i + ninner, l2x, r[0]));
		tripletes.push_back(Trip(i + ninner, l2y, r[1]));
		tripletes.push_back(Trip(i + ninner, l2z, r[2]));

		Energy[i + ninner] = r.dot(r2) - sign;

		// r * r3 +- 1 = 0
		sign = r.dot(r3) > 0 ? 1 : -1;
		tripletes.push_back(Trip(i + ninner * 2, lrx, r3[0]));
		tripletes.push_back(Trip(i + ninner * 2, lry, r3[1]));
		tripletes.push_back(Trip(i + ninner * 2, lrz, r3[2]));

		tripletes.push_back(Trip(i + ninner * 2, l3x, r[0]));
		tripletes.push_back(Trip(i + ninner * 2, l3y, r[1]));
		tripletes.push_back(Trip(i + ninner * 2, l3z, r[2]));

		Energy[i + ninner * 2] = r.dot(r3) - sign;

		// r * r4 +- 1 = 0
		sign = r.dot(r4) > 0 ? 1 : -1;
		tripletes.push_back(Trip(i + ninner * 3, lrx, r4[0]));
		tripletes.push_back(Trip(i + ninner * 3, lry, r4[1]));
		tripletes.push_back(Trip(i + ninner * 3, lrz, r4[2]));

		tripletes.push_back(Trip(i + ninner * 3, l4x, r[0]));
		tripletes.push_back(Trip(i + ninner * 3, l4y, r[1]));
		tripletes.push_back(Trip(i + ninner * 3, l4z, r[2]));

		Energy[i + ninner * 3] = r.dot(r4) - sign;

		/////////////////////////////
	}
}

void lsTools::assemble_solver_binormal_regulizer(Eigen::VectorXd &vars,
												 const LSAnalizer &analizer, const int vars_start_loc,
												 const int aux_start_loc, spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy)
{
	std::vector<Trip> tripletes;
	calculate_binormal_regulizer(vars, analizer, vars_start_loc, aux_start_loc, tripletes, energy);
	int nvars = vars.size();
	int ncondi = energy.size();
	// int ninner = analizer.LocalActInner.size();
	// int rep = 
	spMat J;
	J.resize(ncondi, nvars);
	J.setFromTriplets(tripletes.begin(), tripletes.end());
	H = J.transpose() * J;
	B = -J.transpose() * energy;
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
void lsTools::assemble_solver_boundary_condition_part(const Eigen::VectorXd& func, spMat& H, Eigen::VectorXd& B, Eigen::VectorXd &bcfvalue) {
	spMat bcJacobian;
	assert(trace_vers.size() == trace_hehs.size());
	int size = 0;// the number of constraints
	std::vector<double> vec_elements;
	std::vector<Trip> triplets;
	vec_elements.reserve(lsmesh.n_edges());
	triplets.reserve(lsmesh.n_edges());
	assert(assigned_trace_ls.size() == trace_vers.size());

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
	bcJacobian.setFromTriplets(triplets.begin(), triplets.end());
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
// make some edges have assigned_value
void lsTools::assemble_solver_fixed_values_part(const std::vector<CGMesh::HalfedgeHandle> &hds, const double assigned_value,
												const Eigen::VectorXd &func, spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &bcfvalue)
{
	spMat bcJacobian;
	int size = 0;// the number of constraints
	std::vector<double> vec_elements;
	std::vector<Trip> triplets;
	vec_elements.reserve(hds.size()*2);
	triplets.reserve(hds.size()*2);
	Eigen::VectorXi checked = Eigen::VectorXi::Zero(V.rows());
	for (int i = 0; i < hds.size(); i++)
	{
		CGMesh::HalfedgeHandle edge = hds[i];
		int id_from = lsmesh.from_vertex_handle(edge).idx();
		int id_to = lsmesh.to_vertex_handle(edge).idx();
		if (checked[id_from] == 0)
		{
			triplets.push_back(Trip(size, id_from, 1));
			checked[id_from] = 1;
			vec_elements.push_back(func[id_from] - assigned_value);
			size++;
		}
		if (checked[id_to] == 0)
		{
			triplets.push_back(Trip(size, id_to, 1));
			checked[id_to] = 1;
			vec_elements.push_back(func[id_to] - assigned_value);
			size++;
		}
	}
	// std::cout<<"after loop"<<std::endl;
	// get boundary condition: the jacobian matrix
	bcJacobian.resize(size, V.rows());
	bcJacobian.setFromTriplets(triplets.begin(), triplets.end());
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
void lsTools::assemble_solver_pesudo_geodesic_energy_part_vertex_based(Eigen::VectorXd& vars, const std::vector<double>& angle_degree, 
const LSAnalizer &analizer, const int vars_start_loc, const int aux_start_loc, spMat& H, Eigen::VectorXd& B, Eigen::VectorXd& energy)
{
	std::vector<Trip> tripletes;
	calculate_pseudo_geodesic_opt_expanded_function_values(vars, angle_degree,
		analizer, vars_start_loc, aux_start_loc, tripletes, energy);
	int nvars = vars.size();
	int ncondi = energy.size();
	// int ninner = analizer.LocalActInner.size();
	// int rep = 
	spMat J;
	J.resize(ncondi, nvars);
	J.setFromTriplets(tripletes.begin(), tripletes.end());
	H = J.transpose() * J;
	B = -J.transpose() * energy;
	PGE = energy;
}

// this function gives a matrix which remove the influence of weight_binormal from the energy, thus gives real energy
// size is the size of the pseudo-geodesic energy (for shading),
// 
spMat remove_weight_binormal_influence_in_energy(const int ninner, const int size, const double weight_binormal){
	
	Eigen::VectorXd result = Eigen::VectorXd::Ones(size);
	Eigen::VectorXd weight;
	if(weight_binormal>0){
		weight = Eigen::VectorXd::Ones(ninner * 3) * (1 / weight_binormal);
	}
	else{
		std::cout<<"Warining: You set weight_binormal as 0, meaning you are optimizing nothing!"<<std::endl;
		weight = Eigen::VectorXd::Ones(ninner * 3);
	}
	result.segment(0, ninner * 3) = weight;
	return spMat(result.asDiagonal());
}

// asymptotic, geodesic, shading
void lsTools::assemble_solver_extreme_cases_part_vertex_based(Eigen::VectorXd &vars, const bool asymptotic,
															  const bool use_given_direction,
															  const LSAnalizer &analizer, const int vars_start_loc, spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy)
{
	std::vector<Trip> tripletes;
	int vnbr = V.rows();
	int ninner = analizer.LocalActInner.size();
	if (!asymptotic && use_given_direction)
	{ // shading condition
		if (!enable_shading_init)
		{
			// calculate_shading_condition_auxiliary_vars(vars,
			// 										   analizer, vars_start_loc, vnbr, tripletes, energy);
			calculate_shading_condition_inequivalent(vars,
													   analizer, vars_start_loc, vnbr, tripletes, energy);
			spMat rescale = remove_weight_binormal_influence_in_energy(ninner, energy.size(), weight_binormal);
			energy = rescale * energy;
		}
		else
		{
			calculate_shading_init(vars,
								   analizer, vars_start_loc, vnbr, tripletes, energy);
		}

		//    std::cout<<"trip got"<<std::endl;
	}
	else{// 
		calculate_extreme_pseudo_geodesic_values(vars, asymptotic, analizer, vars_start_loc, tripletes, energy);
	}
	
	int nvars = vars.size();
	int ncondi = energy.size();
	spMat J;
	J.resize(ncondi, nvars);
	J.setFromTriplets(tripletes.begin(), tripletes.end());
	H = J.transpose() * J;
	B = -J.transpose() * energy;
	PGE = energy;
	// if(asymptotic){
	// 	PGE = energy;
	// }
	
}

// the condition that f0+f1+f2=0;
void assemble_AAG_extra_condition(const int mat_size, const int vnbr,
								  const Eigen::VectorXd &func0, const Eigen::VectorXd &func1, const Eigen::VectorXd &func2,
								  spMat &H, Eigen::VectorXd &B, Eigen::VectorXd& energy)
{
	
	
	energy.resize(vnbr); // the nbr of conditions
	std::vector<Trip> tripletes;
	tripletes.reserve(vnbr * 3);
	for(int i=0;i<vnbr;i++){
		tripletes.push_back(Trip(i, i, 1));
		tripletes.push_back(Trip(i, i + vnbr, 1));
		tripletes.push_back(Trip(i, i + vnbr * 2, 1));
		energy[i] = func0[i] + func1[i] + func2[i];
	}
	spMat J;
	J.resize(vnbr, mat_size);
	J.setFromTriplets(tripletes.begin(), tripletes.end());
	H=J.transpose()*J;
	B=-J.transpose()*energy;

}
// the optimized levelset should be othogonal to the given directions on each face id of fids
void lsTools::assemble_solver_othogonal_to_given_face_directions(const Eigen::VectorXd &func, const Eigen::MatrixXd &directions,
														const Eigen::VectorXi &fids, spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy)
{
	std::vector<Trip> tripletes;
	int nvars = V.rows();
	int ncondi = fids.size();
	tripletes.reserve(ncondi*3);
	energy.resize(ncondi);
	for(int i=0;i<ncondi;i++){
		int fid=fids[i];
		int v0=F(fid,0);
		int v1=F(fid,1);
		int v2=F(fid,2);
		// dir0 * V0 + dir1 * V1 + dir2 * V2 is the iso-line direction
		Eigen::Vector3d dir0 = V.row(v2) - V.row(v1);
		Eigen::Vector3d dir1 = V.row(v0) - V.row(v2);
		Eigen::Vector3d dir2 = V.row(v1) - V.row(v0);
		Eigen::Vector3d iso = dir0 * func[v0] + dir1 * func[v1] + dir2 * func[v2];
		double c0 = directions.row(i).dot(dir0);
		double c1 = directions.row(i).dot(dir1);
		double c2 = directions.row(i).dot(dir2);
		double scale = directions.row(i).norm() * iso.norm();
		tripletes.push_back(Trip(i, v0, c0 / scale));
		tripletes.push_back(Trip(i, v1, c1 / scale));
		tripletes.push_back(Trip(i, v2, c2 / scale));
		energy[i] = directions.row(i).dot(iso) / scale;
	}
	spMat J;
	J.resize(ncondi, nvars);
	J.setFromTriplets(tripletes.begin(), tripletes.end());
	H = J.transpose() * J;
	B = -J.transpose() * energy;
}
// the iso-line direction and the gradient direction construct the local frame
// the directions and the grads are the iso-lines and the gradients of the reference level set
// remind that this function does not consider about the location of the variables. If you use this function please
// remember to put H and B into bigger matrices / vectors
void lsTools::assemble_solver_fix_angle_to_given_face_directions(const Eigen::VectorXd &func, const Eigen::MatrixXd &directions,
																 const Eigen::MatrixXd &grads, const double angle_fix,
																 const Eigen::VectorXi &fids, spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy)
{

	std::vector<Trip> tripletes;
	int nvars = V.rows();
	int ncondi = fids.size();
	tripletes.reserve(ncondi * 6);
	energy = Eigen::VectorXd::Zero(ncondi * 2);
	double angle_radian = angle_fix * LSC_PI / 180.;
	for (int i = 0; i < ncondi; i++)
	{
		int fid = fids[i];
		int v0 = F(fid, 0);
		int v1 = F(fid, 1);
		int v2 = F(fid, 2);
		// dir0 * f0 + dir1 * f1 + dir2 * f2 is the iso-line direction

		// direction * iso = cos(theta)
		Eigen::Vector3d dir0 = V.row(v2) - V.row(v1);
		Eigen::Vector3d dir1 = V.row(v0) - V.row(v2);
		Eigen::Vector3d dir2 = V.row(v1) - V.row(v0);
		Eigen::Vector3d iso = dir0 * func[v0] + dir1 * func[v1] + dir2 * func[v2];
		double c0 = directions.row(i).dot(dir0);
		double c1 = directions.row(i).dot(dir1);
		double c2 = directions.row(i).dot(dir2);
		double scale = directions.row(i).norm() * iso.norm();
		tripletes.push_back(Trip(i, v0, c0 / scale));
		tripletes.push_back(Trip(i, v1, c1 / scale));
		tripletes.push_back(Trip(i, v2, c2 / scale));
		energy[i] = (directions.row(i).dot(iso)) / scale - cos(angle_radian);

		// gradient * iso = sin(theta)
		
		c0 = grads.row(i).dot(dir0);
		c1 = grads.row(i).dot(dir1);
		c2 = grads.row(i).dot(dir2);
		scale = grads.row(i).norm() * iso.norm();

		tripletes.push_back(Trip(i + ncondi, v0, c0 / scale));
		tripletes.push_back(Trip(i + ncondi, v1, c1 / scale));
		tripletes.push_back(Trip(i + ncondi, v2, c2 / scale));

		energy[i + ncondi] = (grads.row(i).dot(iso)) / scale - sin(angle_radian);
	}
	spMat J;
	J.resize(ncondi * 2, nvars);
	J.setFromTriplets(tripletes.begin(), tripletes.end());
	H = J.transpose() * J;
	B = -J.transpose() * energy;
}
void lsTools::assemble_solver_fix_two_ls_angle(const Eigen::VectorXd &vars, const double angle_fix,
											   spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy)
{
	std::vector<Trip> tripletes;
	int vnbr = V.rows();
	int fnbr = F.rows();
	tripletes.reserve(fnbr * 12);
	energy = Eigen::VectorXd::Zero(fnbr * 2);
	double angle_radian = angle_fix * LSC_PI / 180.;
	// the directions of the first levelset
	Eigen::MatrixXd directions0 = get_each_face_direction(V, F, vars.segment(0, vnbr));
	// the directions of the second levelset
	Eigen::MatrixXd directions1 = get_each_face_direction(V, F, vars.segment(vnbr, vnbr));

	for (int i = 0; i < fnbr; i++)
	{
		int fid = i;
		int v0 = F(fid, 0);
		int v1 = F(fid, 1);
		int v2 = F(fid, 2);

		int lref0 = v0;
		int lref1 = v1;
		int lref2 = v2;
		int lf0 = v0 + vnbr;
		int lf1 = v1 + vnbr;
		int lf2 = v2 + vnbr;

		double f00 = vars[lref0];
		double f01 = vars[lref1];
		double f02 = vars[lref2];
		double f10 = vars[lf0];
		double f11 = vars[lf1];
		double f12 = vars[lf2];

		// gradient = g0 * fi0 + g1 * fi1 + g2 * fi2
		Eigen::Vector3d g0 =  get_coff_vec_for_gradient(gradVF, fid, v0);
		Eigen::Vector3d g1 =  get_coff_vec_for_gradient(gradVF, fid, v1);
		Eigen::Vector3d g2 =  get_coff_vec_for_gradient(gradVF, fid, v2);
		
		// iso = dir0 * f0 + dir1 * f1 + dir2 * f2 
		Eigen::Vector3d dir0 = V.row(v2) - V.row(v1);
		Eigen::Vector3d dir1 = V.row(v0) - V.row(v2);
		Eigen::Vector3d dir2 = V.row(v1) - V.row(v0);

		Eigen::Vector3d dir = directions1.row(fid);
		Eigen::Vector3d diref = directions0.row(fid);
		Eigen::Vector3d gradref = g0 * f00 + g1 * f01 + g2 * f02;

		// dir * diref = cos(theta)
		double scale = dir.norm() * diref.norm();

		tripletes.push_back(Trip(i, lref0, dir0.dot(dir) / scale));
		tripletes.push_back(Trip(i, lref1, dir1.dot(dir) / scale));
		tripletes.push_back(Trip(i, lref2, dir2.dot(dir) / scale));

		tripletes.push_back(Trip(i, lf0, dir0.dot(diref) / scale));
		tripletes.push_back(Trip(i, lf1, dir1.dot(diref) / scale));
		tripletes.push_back(Trip(i, lf2, dir2.dot(diref) / scale));

		energy[i] = (dir.dot(diref)) / scale - cos(angle_radian);

		// dir * gradref = sin(theta)
		scale = dir.norm() * gradref.norm();

		tripletes.push_back(Trip(i + fnbr, lref0, g0.dot(dir) / scale));
		tripletes.push_back(Trip(i + fnbr, lref1, g1.dot(dir) / scale));
		tripletes.push_back(Trip(i + fnbr, lref2, g2.dot(dir) / scale));

		tripletes.push_back(Trip(i + fnbr, lf0, dir0.dot(gradref) / scale));
		tripletes.push_back(Trip(i + fnbr, lf1, dir1.dot(gradref) / scale));
		tripletes.push_back(Trip(i + fnbr, lf2, dir2.dot(gradref) / scale));

		energy[i + fnbr] = (dir.dot(gradref)) / scale - sin(angle_radian);
	}
	spMat J;
	J.resize(energy.size(), vars.size());
	J.setFromTriplets(tripletes.begin(), tripletes.end());
	H = J.transpose() * J;
	B = -J.transpose() * energy;

}
void lsTools::Run_Level_Set_Opt() {
	
	Eigen::MatrixXd GradValueF, GradValueV;
	Eigen::VectorXd PGEnergy;
	Eigen::VectorXd func = fvalues;
	
	std::vector<double> angle_degree;
	angle_degree.resize(1);
	angle_degree[0] = pseudo_geodesic_target_angle_degree;
	
	int vnbr = V.rows();
	int fnbr = F.rows();
	
	// initialize the level set with some number
	if (func.size() != vnbr)
	{
		if(trace_hehs.size()==0){
			std::cout<<"Please First Set Up The Boundary Condition (Tracing)"<<std::endl;
			return;
		}
		initialize_level_set_accroding_to_parametrization();
		std::cout << "level set get initialized" << std::endl;
		return;
	}
	get_gradient_hessian_values(func, GradValueV, GradValueF);
	if (Glob_lsvars.size() == 0 && trace_hehs.size() == 0)// unit scale when no tracing and before the first iteration
	{
		levelset_unit_scale(func, GradValueF, strip_width);
		std::cout<<"Unit Scale Levelset"<<std::endl;
	}
	// std::cout<<"check "<<func.norm()<<std::endl;
	analysis_pseudo_geodesic_on_vertices(func, analizers[0]);
	int ninner = analizers[0].LocalActInner.size();
	int final_size;
	if (enable_extreme_cases)
	{
		if (Given_Const_Direction)
		{
			final_size = vnbr + ninner * 15; // shading
		}
		else
		{
			final_size = vnbr + ninner * 3; // A or G
		}
	}
	else{
		final_size = ninner * 10 + vnbr;// pseudo - geodesic
	}
	 

	bool need_update_trace_info = (DBdirections.rows() == 0 && trace_hehs.size() > 0) || (Last_Opt_Mesh && trace_hehs.size() > 0); // traced but haven't compute the info
	if (need_update_trace_info)
	{

		get_traced_boundary_triangle_direction_derivatives();
	}
	//  update quantities associated with level set values
	
	
	
	if (Glob_lsvars.size() != final_size) {
		
		if (Glob_lsvars.size() == 0)
		{
			std::cout << "Initializing Global Variable For LevelSet Opt ... " << std::endl;
			Glob_lsvars = Eigen::VectorXd::Zero(final_size); // We change the size if opt more than 1 level set
			Glob_lsvars.segment(0, vnbr) = func;			 
		}
		else{
			if (Glob_lsvars.size() > final_size)
			{ // need to truncate
				Eigen::VectorXd tmp = Glob_lsvars.segment(0, final_size);
				Glob_lsvars = tmp;
				if (Glob_lsvars.size() != final_size)
				{
					std::cout << "size error" << std::endl
							  << std::endl;
				}
			}
			else// need to extend, and re-compute the auxiliaries
			{
				Eigen::VectorXd tmp = Glob_lsvars;
				Glob_lsvars = Eigen::VectorXd::Zero(final_size); // We change the size if opt more than 1 level set
				Glob_lsvars.segment(0, tmp.size()) = tmp;
				Compute_Auxiliaries = true;
			}
		}
	}
	if(Last_Opt_Mesh){
		Compute_Auxiliaries = true;
		std::cout<<"Recomputing Auxiliaries"<<std::endl;
	}
	
	spMat H;
	H.resize(vnbr, vnbr);
	Eigen::VectorXd B=Eigen::VectorXd::Zero(vnbr);
	
	spMat LTL;  // left of laplacian
	Eigen::VectorXd mLTF; // right of laplacian
	assemble_solver_biharmonic_smoothing(func, LTL, mLTF);
	H += weight_laplacian * LTL;
	B += weight_laplacian * mLTF;
	assert(mass.rows() == vnbr);

	// fix inner vers and smooth boundary
	if (enable_inner_vers_fixed)
	{
		// H += weight_mass * InnerV.asDiagonal();
		std::cout << "Please don't use fixing_inner_vertex energy" << std::endl;
	}
	Eigen::VectorXd bcfvalue = Eigen::VectorXd::Zero(vnbr);
	Eigen::VectorXd fbdenergy = Eigen::VectorXd::Zero(vnbr);
	// boundary condition (traced as boundary condition)
	if (trace_hehs.size() > 0)
	{ // if traced, we use boundary condition
		if (!enable_boundary_angles)
		{
			spMat bc_JTJ;
			Eigen::VectorXd bc_mJTF;
			assemble_solver_boundary_condition_part(func, bc_JTJ, bc_mJTF, bcfvalue);
			H += weight_boundary * bc_JTJ;
			B += weight_boundary * bc_mJTF;
		}
		// boundary edge angles
		if (enable_boundary_angles)
		{
			spMat JTJ;
			Eigen::VectorXd mJTF;
			assemble_solver_fixed_boundary_direction_part(GradValueF, tracing_start_edges, func, JTJ, mJTF, fbdenergy);
			H += weight_boundary * JTJ;
			B += weight_boundary * mJTF;
		}
	}
	if (interactive_flist.size() > 0)
    { // if traced, we use boundary condition

        spMat bc_JTJ;
        Eigen::VectorXd bc_mJTF;
        assemble_solver_interactive_boundary_condition_part(func, bc_JTJ, bc_mJTF, bcfvalue);
        H += weight_boundary * bc_JTJ;
        B += weight_boundary * bc_mJTF;
    }

	// strip width condition
	if (enable_strip_width_energy)
	{
		spMat sw_JTJ;
		Eigen::VectorXd sw_mJTF;
		assemble_solver_strip_width_part(GradValueF, sw_JTJ, sw_mJTF);// by default the strip width is 1. Unless tracing info updated the info
		H += weight_strip_width * sw_JTJ;
		B += weight_strip_width * sw_mJTF;
	}

	
	assert(H.rows() == vnbr);
	assert(H.cols() == vnbr);
	// pseudo geodesic
	spMat Hlarge;
	Eigen::VectorXd Blarge;

	Eigen::VectorXd e_smbi;
	Hlarge.resize(final_size, final_size);
	Blarge = Eigen::VectorXd::Zero(final_size);
	
	Hlarge = sum_uneven_spMats(H, Hlarge);
	Blarge = sum_uneven_vectors(B, Blarge);

	if (enable_pseudo_geodesic_energy)
	{

		spMat pg_JTJ;
		Eigen::VectorXd pg_mJTF;
		spMat smbi_H;
		Eigen::VectorXd smbi_B;
		int vars_start_loc = 0;
		int aux_start_loc = vnbr;
		if (!enable_extreme_cases) {

			assemble_solver_pesudo_geodesic_energy_part_vertex_based(Glob_lsvars, angle_degree, analizers[0],
																	vars_start_loc, aux_start_loc, pg_JTJ, pg_mJTF, PGEnergy);
			
		}
		else {
			bool asymptotic = true;
			if (angle_degree[0] != 0)
			{
				asymptotic = false;
			}
			if(Given_Const_Direction){
				std::cout << "Direc, "<<angle_ray_converter(Reference_theta, Reference_phi).transpose();

			}
			// get the vertices associated to the second angle
			analizers[0].Special = shading_condition_info;

			// mark the angles
			Eigen::VectorXd diff;
			// analizers[0].ShadSpecial = shading_detect_parallel_patch(Reference_theta, Reference_phi, diff);

			if (enable_max_energy_check && analizers[0].HighEnergy.size() == 0) // mark the max energy points
			{
				mark_high_energy_vers(PGE, ninner, max_energy_percentage, IVids, analizers[0].HighEnergy, refids);
			}

			// std::cout<<"extreme before"<<std::endl;
			assemble_solver_extreme_cases_part_vertex_based(Glob_lsvars, asymptotic, Given_Const_Direction,
															analizers[0], vars_start_loc, pg_JTJ, pg_mJTF, PGEnergy);
			// std::cout<<"extreme computed"<<std::endl;
			assemble_solver_binormal_regulizer(Glob_lsvars, analizers[0], vars_start_loc, aux_start_loc, smbi_H, smbi_B, e_smbi);
			Hlarge = sum_uneven_spMats(Hlarge, weight_smt_binormal * smbi_H);
			Blarge = sum_uneven_vectors(Blarge, weight_smt_binormal * smbi_B);

			
		}

		Hlarge = sum_uneven_spMats(Hlarge, weight_pseudo_geodesic_energy * pg_JTJ);
		Blarge = sum_uneven_vectors(Blarge, weight_pseudo_geodesic_energy * pg_mJTF);
		
		
		// std::cout<<"extreme computed 1"<<std::endl;
		// if (enable_extreme_cases && Given_Const_Direction)
		// {
		// 	Hlarge += sum_uneven_spMats(Hlarge, weight_shading * sd_H);
		// 	Blarge += sum_uneven_vectors(Blarge, weight_shading * sd_B);
		// }
		Compute_Auxiliaries = false;
	}
	if(vector_contains_NAN(Blarge)){
        std::cout<<"energy contains NAN"<<std::endl;
		return;
    }
	Hlarge += 1e-6 * weight_mass * spMat(Eigen::VectorXd::Ones(final_size).asDiagonal());

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
	fvalues = func;
	Glob_lsvars += dx;
	double energy_biharmonic = (QcH * func).norm();
	double energy_boundary = (bcfvalue).norm();
	std::cout << "energy: harm " << energy_biharmonic << ", bnd " << energy_boundary << ", ";
	if (enable_pseudo_geodesic_energy)
	{
		if (!enable_extreme_cases) {
			double energy_pg = PGEnergy.norm();
			double max_energy_ls = PGEnergy.lpNorm<Eigen::Infinity>();
			std::cout << "pg, " << energy_pg << ", " << "lsmax," << max_energy_ls << ",";
		}
		else {
			double planar_energy = PGEnergy.norm();
			double max_energy_ls = PGEnergy.lpNorm<Eigen::Infinity>();
			std::cout << "extreme_pg, " << planar_energy << " lsmax," << max_energy_ls<<", ";
			if(e_smbi.size()>0){
				std::cout<< "smbi, "<<e_smbi.norm()<<", ";
			}
			 
		}
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
	std::cout << "step " << step_length;
	if (interactive_flist.size() > 0) // if start to solve pseudo-geodeic, we refer to the strokes
	{
		pseudo_geodesic_target_angle_degree = 180. / LSC_PI * get_interactive_angle(func, analizers[0], lsmesh, norm_v, V, F, interactive_flist, interactive_bclist, InnerV);
		std::cout << ", target_angle, " << pseudo_geodesic_target_angle_degree;
	}

	std::cout << std::endl;
	Last_Opt_Mesh = false;
}
// the type starts from -1
void get_ver_nbr_for_each_class(const Eigen::VectorXi &info, std::vector<std::vector<int>>& types)
{
	int type_low = info.minCoeff();
	int type_up = info.maxCoeff();
	int nbr_types = type_up + 2;
	int ninner = info.size();

	types.resize(nbr_types);
	std::cout<<"Types nbr: "<<nbr_types<<", min "<<type_low<<", max "<<type_up<<std::endl;
	for(int i =0;i<ninner;i++){
		int type = info[i] + 1;// the type starts from -1
		types[type].push_back(i);
	}
}
// 4 types: the transition (-1), the asymptotic (0), the angle of 45 degree area (1), the geodesic (2)
void lsTools::Run_Level_Set_Opt_Angle_Variable() {
	
	Eigen::MatrixXd GradValueF, GradValueV;
	Eigen::VectorXd PGEnergy[3];
	Eigen::VectorXd func = fvalues;

	
	int vnbr = V.rows();
	int fnbr = F.rows();
	
	// initialize the level set with some number
	if (func.size() != vnbr)
	{
		if(trace_hehs.size()==0){
			std::cout<<"Please First Set Up The Boundary Condition (Tracing)"<<std::endl;
			return;
		}
		initialize_level_set_accroding_to_parametrization();
		std::cout << "level set get initialized" << std::endl;
		return;
	}
	get_gradient_hessian_values(func, GradValueV, GradValueF);
	if (Glob_lsvars.size() == 0 && trace_hehs.size() == 0)// unit scale when no tracing and before the first iteration
	{
		levelset_unit_scale(func, GradValueF, strip_width);
		std::cout<<" Unit Scale Levelset finished"<<std::endl;
	}
	// std::cout<<"check "<<func.norm()<<std::endl;
	analysis_pseudo_geodesic_on_vertices(func, analizers[0]);
	int ninner = analizers[0].LocalActInner.size();

	if (shading_condition_info.size() == 0)
	{
		std::cout << "Please load the reference meshes and the base mesh" << std::endl;
	}

	if (Glob_lsvars.size() == 0)
	{
		std::cout<<"printing ver nbr for each class..."<<std::endl;
		get_ver_nbr_for_each_class(shading_condition_info, PG_Vertex_Types);
		std::cout << "Types of Vers: " << PG_Vertex_Types.size() << "\n ";
		for (int i = 0; i < PG_Vertex_Types.size(); i++)
		{
			std::cout << "type " << i << ", nbr of vers, " << PG_Vertex_Types[i].size() << std::endl;
		}
		// std::cout<<"Before go inside"<<std::endl;
		get_geodesic_distance(V, F, IVids, PG_Vertex_Types[1], PGGeoDistance);
		assign_pg_angle_based_on_geodesic_distance(PGGeoDistance, pseudo_geodesic_target_angle_degree, PGVariationalAngles);
	}
	
	// pseudo_geodesic points PG_Vertex_Types[0], auxiliary: vnbr + size * 10
	// asymptotic points PG_Vertex_Types[1], auxiliary: vnbr

	int final_size = vnbr + ninner * 10;
	
	if (Glob_lsvars.size() != final_size) {
		
		if (Glob_lsvars.size() == 0)
		{
			std::cout << "Initializing Global Variable For LevelSet Opt ... " << std::endl;
			Glob_lsvars = Eigen::VectorXd::Zero(final_size); // We change the size if opt more than 1 level set
			Glob_lsvars.segment(0, vnbr) = func;			 
		}
		else{
			if (Glob_lsvars.size() > final_size)
			{ // need to truncate
				Eigen::VectorXd tmp = Glob_lsvars.segment(0, final_size);
				Glob_lsvars = tmp;
				if (Glob_lsvars.size() != final_size)
				{
					std::cout << "size error" << std::endl
							  << std::endl;
				}
			}
			else// need to extend, and re-compute the auxiliaries
			{
				Eigen::VectorXd tmp = Glob_lsvars;
				Glob_lsvars = Eigen::VectorXd::Zero(final_size); // We change the size if opt more than 1 level set
				Glob_lsvars.segment(0, tmp.size()) = tmp;
				Compute_Auxiliaries = true;
			}
		}
	}
	
	spMat H;
	H.resize(vnbr, vnbr);
	Eigen::VectorXd B=Eigen::VectorXd::Zero(vnbr);
	
	spMat LTL;  // left of laplacian
	Eigen::VectorXd mLTF; // right of laplacian
	assemble_solver_biharmonic_smoothing(func, LTL, mLTF);
	H += weight_laplacian * LTL;
	B += weight_laplacian * mLTF;
	assert(mass.rows() == vnbr);

	// strip width condition
	if (enable_strip_width_energy)
	{
		spMat sw_JTJ;
		Eigen::VectorXd sw_mJTF;
		assemble_solver_strip_width_part(GradValueF, sw_JTJ, sw_mJTF);// by default the strip width is 1. Unless tracing info updated the info
		H += weight_strip_width * sw_JTJ;
		B += weight_strip_width * sw_mJTF;
	}

	
	assert(H.rows() == vnbr);
	assert(H.cols() == vnbr);
	// pseudo geodesic
	spMat Hlarge;
	Eigen::VectorXd Blarge;

	Eigen::VectorXd e_smbi;
	Hlarge.resize(final_size, final_size);
	Blarge = Eigen::VectorXd::Zero(final_size);
	
	Hlarge = sum_uneven_spMats(H, Hlarge);
	Blarge = sum_uneven_vectors(B, Blarge);

	if (enable_pseudo_geodesic_energy)
	{

		spMat pg_JTJ[3];
		Eigen::VectorXd pg_mJTF[3];
		spMat smbi_H;
		Eigen::VectorXd smbi_B;
		int vars_start_loc = 0;
		int aux_start_loc = vnbr;
		// if (!enable_extreme_cases) {

		// 	assemble_solver_pesudo_geodesic_energy_part_vertex_based(Glob_lsvars, angle_degree, analizers[0],
		// 															vars_start_loc, aux_start_loc, pg_JTJ, pg_mJTF, PGEnergy);
		// }
		// else {
			
		// }
		bool asymptotic = true;
		assemble_solver_extreme_cases_part_vertex_based(Glob_lsvars, asymptotic, false,
														analizers[0], vars_start_loc, pg_JTJ[0], pg_mJTF[0], PGEnergy[0]);
		
		// assemble_solver_binormal_regulizer(Glob_lsvars, analizers[0], vars_start_loc, aux_start_loc, smbi_H, smbi_B, e_smbi);
		// Hlarge = sum_uneven_spMats(Hlarge, weight_smt_binormal * smbi_H);
		// Blarge = sum_uneven_vectors(Blarge, weight_smt_binormal * smbi_B);
		Hlarge = sum_uneven_spMats(Hlarge, weight_pseudo_geodesic_energy * pg_JTJ[0]);
		Blarge = sum_uneven_vectors(Blarge, weight_pseudo_geodesic_energy * pg_mJTF[0]);
		if (enable_extreme_cases)
		{
			std::vector<double> angle_list;
			int vars_start_loc = 0;
			int aux_start_loc = vnbr;
			assemble_solver_pesudo_geodesic_energy_part_vertex_based(Glob_lsvars, angle_list, analizers[0], vars_start_loc, aux_start_loc, pg_JTJ[1],
																	 pg_mJTF[1], PGEnergy[1]);

			Hlarge = sum_uneven_spMats(Hlarge, weight_pseudo_geodesic_energy * pg_JTJ[1]);
			Blarge = sum_uneven_vectors(Blarge, weight_pseudo_geodesic_energy * pg_mJTF[1]);
			Compute_Auxiliaries = false;
		}

		
	}
	if(vector_contains_NAN(Blarge)){
        std::cout<<"energy contains NAN"<<std::endl;
		return;
    }
	Hlarge += 1e-6 * weight_mass * spMat(Eigen::VectorXd::Ones(final_size).asDiagonal());

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
	fvalues = func;
	Glob_lsvars += dx;
	double energy_biharmonic = (QcH * func).norm();
	std::cout << "energy: harm " << energy_biharmonic << ", ";
	if (enable_pseudo_geodesic_energy)
	{
		std::cout<<"pg, "<<PGEnergy[0].norm()<<", ";
		if(enable_extreme_cases){
			std::cout << PGEnergy[1].norm() << ", ";
		}
		
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

	step_length = dx.norm();
	std::cout << "step " << step_length;
	std::cout << std::endl;
	Last_Opt_Mesh = false;
}
void lsTools::Run_AAG(Eigen::VectorXd& func0, Eigen::VectorXd& func1, Eigen::VectorXd& func2){
	Eigen::MatrixXd GradValueF[3], GradValueV[3];
	Eigen::VectorXd PGEnergy[3];

	int vnbr = V.rows();
	int fnbr = F.rows();
	bool first_compute = true; // if we need initialize auxiliary vars
	get_gradient_hessian_values(func0, GradValueV[0], GradValueF[0]);
	get_gradient_hessian_values(func1, GradValueV[1], GradValueF[1]);
	
	// initialize the level set with some number
	if (Glob_lsvars.size() == 0)
	{ // initialize the 3rd levelset only here
		levelset_unit_scale(func0, GradValueF[0], 1);
		levelset_unit_scale(func1, GradValueF[1], 1);
		func2 = -func0 - func1; // func0 + func1 + func2 = 0
	}
	get_gradient_hessian_values(func2, GradValueV[2], GradValueF[2]);

	analysis_pseudo_geodesic_on_vertices(func0, analizers[0]);
	analysis_pseudo_geodesic_on_vertices(func1, analizers[1]);
	analysis_pseudo_geodesic_on_vertices(func2, analizers[2]);
	int ninner = analizers[0].LocalActInner.size();
	int final_size = vnbr * 3; // Change this when using more auxilary vars.

	
	
	
	if (Glob_lsvars.size() == 0) {
		
		first_compute = true;
		std::cout << "Initializing Global Variable For LevelSet Opt ... " << std::endl;
		Glob_lsvars = Eigen::VectorXd::Zero(final_size);// We change the size if opt more than 1 level set
		Glob_lsvars.segment(0, vnbr) = func0;
		Glob_lsvars.segment(vnbr, vnbr) = func1;
		Glob_lsvars.segment(vnbr * 2, vnbr) = func2;
	}
	if(Last_Opt_Mesh){
		Compute_Auxiliaries = true;
		std::cout<<"Recomputing Auxiliaries"<<std::endl;
	}
	
	spMat H;
	H.resize(final_size, final_size);
	Eigen::VectorXd B=Eigen::VectorXd::Zero(final_size);

	spMat LTL0, LTL1, LTL2;				 // left of laplacian
	Eigen::VectorXd mLTF0, mLTF1, mLTF2; // right of laplacian
	assemble_solver_biharmonic_smoothing(func0, LTL0, mLTF0);
	assemble_solver_biharmonic_smoothing(func1, LTL1, mLTF1);
	assemble_solver_biharmonic_smoothing(func2, LTL2, mLTF2);
	H += weight_laplacian * three_spmat_in_diag(LTL0, LTL1, weight_geodesic* LTL2, final_size);
	B += weight_laplacian * three_vec_in_row(mLTF0, mLTF1, weight_geodesic* mLTF2, final_size);

	// strip width condition
	
	if (enable_strip_width_energy)
	{
		spMat sw_JTJ[3];
		Eigen::VectorXd sw_mJTF[3];
		assemble_solver_strip_width_part(GradValueF[0], sw_JTJ[0], sw_mJTF[0]);// by default the strip width is 1. Unless tracing info updated the info
		assemble_solver_strip_width_part(GradValueF[1], sw_JTJ[1], sw_mJTF[1]);// by default the strip width is 1. Unless tracing info updated the info
		assemble_solver_strip_width_part(GradValueF[2], sw_JTJ[2], sw_mJTF[2]);// by default the strip width is 1. Unless tracing info updated the info
		H += weight_strip_width * three_spmat_in_diag(sw_JTJ[0], sw_JTJ[1], weight_geodesic* sw_JTJ[2], final_size);
		B += weight_strip_width * three_vec_in_row(sw_mJTF[0], sw_mJTF[1], weight_geodesic* sw_mJTF[2], final_size);
	}
	spMat extraH;
	Eigen::VectorXd extraB;
	Eigen::VectorXd extra_energy;
	assemble_AAG_extra_condition(final_size, vnbr, func0, func1, func2, extraH, extraB, extra_energy);
	H += weight_boundary * extraH;
	B += weight_boundary * extraB;
	if (enable_pseudo_geodesic_energy )
	{

		spMat pg_JTJ[3];
		Eigen::VectorXd pg_mJTF[3];
		int vars_start_loc = 0;
		
		assemble_solver_extreme_cases_part_vertex_based(Glob_lsvars, true, false, analizers[0], vars_start_loc, pg_JTJ[0], pg_mJTF[0], PGEnergy[0]);
		vars_start_loc = vnbr;
		assemble_solver_extreme_cases_part_vertex_based(Glob_lsvars, true, false, analizers[1], vars_start_loc, pg_JTJ[1], pg_mJTF[1], PGEnergy[1]);
		vars_start_loc = vnbr * 2;
		int aux_start_loc = vnbr * 3;
		std::vector<double> angle_degree(1);
		angle_degree[0]=90;

		assemble_solver_extreme_cases_part_vertex_based(Glob_lsvars, false, false, analizers[2],
														vars_start_loc,
														pg_JTJ[2], pg_mJTF[2], PGEnergy[2]);
		Compute_Auxiliaries = false;
		H += weight_pseudo_geodesic_energy * (pg_JTJ[0] + pg_JTJ[1] + weight_geodesic* pg_JTJ[2]);
		B += weight_pseudo_geodesic_energy * (pg_mJTF[0] + pg_mJTF[1] + weight_geodesic* pg_mJTF[2]);
	}
	
	H += 1e-6 * weight_mass * spMat(Eigen::VectorXd::Ones(final_size).asDiagonal());

	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver(H);

	// assert(solver.info() == Eigen::Success);
	if (solver.info() != Eigen::Success)
	{
		// solving failed
		std::cout << "solver fail" << std::endl;
		return;
	}
	// std::cout<<"solved successfully"<<std::endl;
	Eigen::VectorXd dx = solver.solve(B).eval();
	dx *= 0.75;
	// std::cout << "step length " << dx.norm() << std::endl;
	double level_set_step_length = dx.norm();
	if (level_set_step_length > max_step_length)
	{
		dx *= max_step_length / level_set_step_length;
	}
	

	// energy evaluation
	double energy_biharmonic[3];
	energy_biharmonic[0] = (QcH * func0).norm();
	energy_biharmonic[1] = (QcH * func1).norm();
	energy_biharmonic[2] = (QcH * func2).norm();
	std::cout << "energy: harm " << energy_biharmonic[0] << ", "<<energy_biharmonic[1] << ", "<<energy_biharmonic[2] << ", ";
	if (enable_pseudo_geodesic_energy)
	{
		double energy_pg[3];
		energy_pg[0]=PGEnergy[0].norm();
		energy_pg[1]=PGEnergy[1].norm();
		energy_pg[2]=PGEnergy[2].norm();
		std::cout << "pg, " << energy_pg[0]<<", "<< energy_pg[1]<<", "<< energy_pg[2]<<", pgmax, "<<PGEnergy[0].lpNorm<Eigen::Infinity>()
		<<", "<<PGEnergy[1].lpNorm<Eigen::Infinity>()<<", "<<PGEnergy[2].lpNorm<Eigen::Infinity>()<<", ";
	}

	if (enable_strip_width_energy)
	{

		Eigen::VectorXd ener0, ener1, ener2; 
		ener0= GradValueF[0].rowwise().norm();
		ener1= GradValueF[1].rowwise().norm();
		ener2= GradValueF[2].rowwise().norm();
		ener0 = ener0.asDiagonal() * ener0;
		ener1 = ener1.asDiagonal() * ener1;
		ener2 = ener2.asDiagonal() * ener2;
		Eigen::VectorXd wds = Eigen::VectorXd::Ones(fnbr) * strip_width * strip_width;

		ener0 -= wds;
		ener1 -= wds;
		ener2 -= wds;
		double stp_energy[3];
		stp_energy[0] = Eigen::VectorXd(ener0).dot(ener0);
		stp_energy[1] = Eigen::VectorXd(ener1).dot(ener1);
		stp_energy[2] = Eigen::VectorXd(ener2).dot(ener2);
		std::cout << "strip, " << stp_energy[0] << ", "<< stp_energy[1] << ", "<< stp_energy[2] << ", ";
	}
	std::cout<<"extra, "<<extra_energy.norm()<<", ";
	step_length = dx.norm();
	std::cout << "step " << step_length << std::endl;

	func0 += dx.topRows(vnbr);
	func1 += dx.middleRows(vnbr,vnbr);
	func2 += dx.middleRows(vnbr * 2, vnbr);
	Glob_lsvars += dx;
	Last_Opt_Mesh = false;
}

// in the order of AGG
void lsTools::Run_AGG(Eigen::VectorXd& func0, Eigen::VectorXd& func1, Eigen::VectorXd& func2){
	Eigen::MatrixXd GradValueF[3], GradValueV[3];
	Eigen::VectorXd PGEnergy[3], Eangle;

	int vnbr = V.rows();
	int fnbr = F.rows();
	bool first_compute = true; // if we need initialize auxiliary vars
	get_gradient_hessian_values(func0, GradValueV[0], GradValueF[0]);
	get_gradient_hessian_values(func1, GradValueV[1], GradValueF[1]);
	
	// initialize the level set with some number
	if (Glob_lsvars.size() == 0)
	{ // initialize the 3rd levelset only here
		levelset_unit_scale(func0, GradValueF[0], 1);
		levelset_unit_scale(func1, GradValueF[1], 1);
		func2 = -func0 - func1; // func0 + func1 + func2 = 0
		// return;
	}
	get_gradient_hessian_values(func2, GradValueV[2], GradValueF[2]);

	analysis_pseudo_geodesic_on_vertices(func0, analizers[0]);
	analysis_pseudo_geodesic_on_vertices(func1, analizers[1]);
	analysis_pseudo_geodesic_on_vertices(func2, analizers[2]);
	int ninner = analizers[0].LocalActInner.size();
	int final_size = vnbr * 3; // Change this when using more auxilary vars. 
	
	if (Glob_lsvars.size() == 0) {
		
		first_compute = true;
		std::cout << "Initializing Global Variable For LevelSet Opt ... " << std::endl;
		Glob_lsvars = Eigen::VectorXd::Zero(final_size);// We change the size if opt more than 1 level set
		Glob_lsvars.segment(0, vnbr) = func0;
		Glob_lsvars.segment(vnbr, vnbr) = func1;
		Glob_lsvars.segment(vnbr * 2, vnbr) = func2;
	}
	if(Last_Opt_Mesh){
		Compute_Auxiliaries = true;
		std::cout<<"Recomputing Auxiliaries"<<std::endl;
	}
	
	spMat H;
	H.resize(final_size, final_size);
	Eigen::VectorXd B=Eigen::VectorXd::Zero(final_size);

	spMat LTL0, LTL1, LTL2;				 // left of laplacian
	Eigen::VectorXd mLTF0, mLTF1, mLTF2; // right of laplacian
	assemble_solver_biharmonic_smoothing(func0, LTL0, mLTF0);
	assemble_solver_biharmonic_smoothing(func1, LTL1, mLTF1);
	assemble_solver_biharmonic_smoothing(func2, LTL2, mLTF2);
	H += weight_laplacian * three_spmat_in_diag(LTL0, LTL1, weight_geodesic* LTL2, final_size);
	B += weight_laplacian * three_vec_in_row(mLTF0, mLTF1, weight_geodesic* mLTF2, final_size);

	// strip width condition
	
	if (enable_strip_width_energy)
	{
		spMat sw_JTJ[3];
		Eigen::VectorXd sw_mJTF[3];
		assemble_solver_strip_width_part(GradValueF[0], sw_JTJ[0], sw_mJTF[0]);// by default the strip width is 1. Unless tracing info updated the info
		assemble_solver_strip_width_part(GradValueF[1], sw_JTJ[1], sw_mJTF[1]);// by default the strip width is 1. Unless tracing info updated the info
		assemble_solver_strip_width_part(GradValueF[2], sw_JTJ[2], sw_mJTF[2]);// by default the strip width is 1. Unless tracing info updated the info
		H += weight_strip_width * three_spmat_in_diag(sw_JTJ[0], sw_JTJ[1], weight_geodesic* sw_JTJ[2], final_size);
		B += weight_strip_width * three_vec_in_row(sw_mJTF[0], sw_mJTF[1], weight_geodesic* sw_mJTF[2], final_size);
	}
	spMat extraH;
	Eigen::VectorXd extraB;
	Eigen::VectorXd extra_energy;
	assemble_AAG_extra_condition(final_size, vnbr, func0, func1, func2, extraH, extraB, extra_energy);
	H += weight_boundary * extraH;
	B += weight_boundary * extraB;
	if (enable_pseudo_geodesic_energy)
	{
		spMat pg_JTJ[3], faH;
		Eigen::VectorXd pg_mJTF[3], faB;
		int vars_start_loc = 0;
		// A
		assemble_solver_extreme_cases_part_vertex_based(Glob_lsvars, true, false, analizers[0], vars_start_loc, pg_JTJ[0], pg_mJTF[0], PGEnergy[0]);
		vars_start_loc = vnbr;
		// G
		assemble_solver_extreme_cases_part_vertex_based(Glob_lsvars, false, false, analizers[1], vars_start_loc, pg_JTJ[1], pg_mJTF[1], PGEnergy[1]);
		vars_start_loc = vnbr * 2;
		// G
		assemble_solver_extreme_cases_part_vertex_based(Glob_lsvars, false, false, analizers[2], vars_start_loc, pg_JTJ[2], pg_mJTF[2], PGEnergy[2]);
		Compute_Auxiliaries = false;
		H += weight_pseudo_geodesic_energy * (pg_JTJ[0] + pg_JTJ[1] + weight_geodesic* pg_JTJ[2]);
		B += weight_pseudo_geodesic_energy * (pg_mJTF[0] + pg_mJTF[1] + weight_geodesic * pg_mJTF[2]);
		if (fix_angle_of_two_levelsets)
		{
			// fix angle between the first (A) and the second (G) level set
			assemble_solver_fix_two_ls_angle(Glob_lsvars, angle_between_two_levelsets, faH, faB, Eangle);
			H = sum_uneven_spMats(H, weight_fix_two_ls_angle * faH);
			B = sum_uneven_vectors(B, weight_fix_two_ls_angle * faB);
		}
	}
	
	H += 1e-6 * weight_mass * spMat(Eigen::VectorXd::Ones(final_size).asDiagonal());

	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver(H);

	// assert(solver.info() == Eigen::Success);
	if (solver.info() != Eigen::Success)
	{
		// solving failed
		std::cout << "solver fail" << std::endl;
		return;
	}
	// std::cout<<"solved successfully"<<std::endl;
	Eigen::VectorXd dx = solver.solve(B).eval();
	dx *= 0.75;
	// std::cout << "step length " << dx.norm() << std::endl;
	double level_set_step_length = dx.norm();
	if (level_set_step_length > max_step_length)
	{
		dx *= max_step_length / level_set_step_length;
	}
	

	// energy evaluation
	double energy_biharmonic[3];
	energy_biharmonic[0] = (QcH * func0).norm();
	energy_biharmonic[1] = (QcH * func1).norm();
	energy_biharmonic[2] = (QcH * func2).norm();
	std::cout << "energy: harm " << energy_biharmonic[0] << ", "<<energy_biharmonic[1] << ", "<<energy_biharmonic[2] << ", ";
	if (enable_pseudo_geodesic_energy)
	{
		double energy_pg[3];
		energy_pg[0]=PGEnergy[0].norm();
		energy_pg[1]=PGEnergy[1].norm();
		energy_pg[2]=PGEnergy[2].norm();
		std::cout << "AGGpg, " << energy_pg[0]<<", "<< energy_pg[1]<<", "<< energy_pg[2]<<", pgmax, "<<PGEnergy[0].lpNorm<Eigen::Infinity>()
		<<", "<<PGEnergy[1].lpNorm<Eigen::Infinity>()<<", "<<PGEnergy[2].lpNorm<Eigen::Infinity>()<<", ";
	}

	if (enable_strip_width_energy)
	{

		Eigen::VectorXd ener0, ener1, ener2; 
		ener0= GradValueF[0].rowwise().norm();
		ener1= GradValueF[1].rowwise().norm();
		ener2= GradValueF[2].rowwise().norm();
		ener0 = ener0.asDiagonal() * ener0;
		ener1 = ener1.asDiagonal() * ener1;
		ener2 = ener2.asDiagonal() * ener2;
		Eigen::VectorXd wds = Eigen::VectorXd::Ones(fnbr) * strip_width * strip_width;

		ener0 -= wds;
		ener1 -= wds;
		ener2 -= wds;
		double stp_energy[3];
		stp_energy[0] = Eigen::VectorXd(ener0).dot(ener0);
		stp_energy[1] = Eigen::VectorXd(ener1).dot(ener1);
		stp_energy[2] = Eigen::VectorXd(ener2).dot(ener2);
		std::cout << "strip, " << stp_energy[0] << ", "<< stp_energy[1] << ", "<< stp_energy[2] << ", ";
	}
	std::cout<<"extra, "<<extra_energy.norm()<<", ";
	std::cout<<"fixangle, "<<Eangle.norm()<<", ";
	step_length = dx.norm();
	std::cout << "step " << step_length << std::endl;

	func0 += dx.topRows(vnbr);
	func1 += dx.middleRows(vnbr,vnbr);
	func2 += dx.middleRows(vnbr * 2, vnbr);
	Glob_lsvars += dx;
	Last_Opt_Mesh = false;
}

// in the order of ppg: pseudo-geodesic, pseudo-geodesic and geodesic
void lsTools::Run_PPG(Eigen::VectorXd& func0, Eigen::VectorXd& func1, Eigen::VectorXd& func2){
	Eigen::MatrixXd GradValueF[3], GradValueV[3];
	Eigen::VectorXd PGEnergy[3], Eangle;

	int vnbr = V.rows();
	int fnbr = F.rows();
	bool first_compute = true; // if we need initialize auxiliary vars
	get_gradient_hessian_values(func0, GradValueV[0], GradValueF[0]);
	get_gradient_hessian_values(func1, GradValueV[1], GradValueF[1]);
	
	// initialize the level set with some number
	if (Glob_lsvars.size() == 0)
	{ // initialize the 3rd levelset only here
		levelset_unit_scale(func0, GradValueF[0], 1);
		levelset_unit_scale(func1, GradValueF[1], 1);
		func2 = -func0 - func1; // func0 + func1 + func2 = 0
		// return;
	}
	get_gradient_hessian_values(func2, GradValueV[2], GradValueF[2]);

	analysis_pseudo_geodesic_on_vertices(func0, analizers[0]);
	analysis_pseudo_geodesic_on_vertices(func1, analizers[1]);
	analysis_pseudo_geodesic_on_vertices(func2, analizers[2]);
	int ninner = analizers[0].LocalActInner.size();
	int final_size = (ninner * 10 + vnbr) * 2 + vnbr; // Change this when using more auxilary vars.

	if (Glob_lsvars.size() == 0) {
		
		first_compute = true;
		std::cout << "Initializing Global Variable For LevelSet Opt ... " << std::endl;
		Glob_lsvars = Eigen::VectorXd::Zero(final_size);// We change the size if opt more than 1 level set
		Glob_lsvars.segment(0, vnbr) = func0;
		Glob_lsvars.segment(vnbr, vnbr) = func1;
		Glob_lsvars.segment(vnbr * 2, vnbr) = func2;
	}
	if(Last_Opt_Mesh){
		Compute_Auxiliaries = true;
		std::cout<<"Recomputing Auxiliaries"<<std::endl;
	}
	
	spMat H;
	H.resize(final_size, final_size);
	Eigen::VectorXd B = Eigen::VectorXd::Zero(final_size);

	spMat LTL0, LTL1, LTL2;				 // left of laplacian
	Eigen::VectorXd mLTF0, mLTF1, mLTF2; // right of laplacian
	assemble_solver_biharmonic_smoothing(func0, LTL0, mLTF0);
	assemble_solver_biharmonic_smoothing(func1, LTL1, mLTF1);
	assemble_solver_biharmonic_smoothing(func2, LTL2, mLTF2);
	H += weight_laplacian * three_spmat_in_diag(LTL0, LTL1, weight_geodesic* LTL2, final_size);
	B += weight_laplacian * three_vec_in_row(mLTF0, mLTF1, weight_geodesic* mLTF2, final_size);

	// strip width condition
	
	if (enable_strip_width_energy)
	{
		spMat sw_JTJ[3];
		Eigen::VectorXd sw_mJTF[3];
		assemble_solver_strip_width_part(GradValueF[0], sw_JTJ[0], sw_mJTF[0]);// by default the strip width is 1. Unless tracing info updated the info
		assemble_solver_strip_width_part(GradValueF[1], sw_JTJ[1], sw_mJTF[1]);// by default the strip width is 1. Unless tracing info updated the info
		assemble_solver_strip_width_part(GradValueF[2], sw_JTJ[2], sw_mJTF[2]);// by default the strip width is 1. Unless tracing info updated the info
		H += weight_strip_width * three_spmat_in_diag(sw_JTJ[0], sw_JTJ[1], weight_geodesic* sw_JTJ[2], final_size);
		B += weight_strip_width * three_vec_in_row(sw_mJTF[0], sw_mJTF[1], weight_geodesic* sw_mJTF[2], final_size);
	}
	spMat extraH;
	Eigen::VectorXd extraB;
	Eigen::VectorXd extra_energy;
	assemble_AAG_extra_condition(final_size, vnbr, func0, func1, func2, extraH, extraB, extra_energy);
	H += weight_boundary * extraH;
	B += weight_boundary * extraB;
	spMat pg_JTJ[3], faH;
	Eigen::VectorXd pg_mJTF[3], faB;
	if (enable_pseudo_geodesic_energy)
	{
		// the variables:
		// func0 (vnbr), func1 (vnbr), func2 (vnbr), aux0（ninner * 10), aux1（ninner * 10)
		std::vector<double> target_angles_0(1), target_angles_1(1);
		target_angles_0[0] = pseudo_geodesic_target_angle_degree;
		target_angles_1[0] = pseudo_geodesic_target_angle_degree_2;

		// P1
		int vars_start_loc = 0;
		int aux_start_loc = vnbr * 3;
		assemble_solver_pesudo_geodesic_energy_part_vertex_based(Glob_lsvars, target_angles_0, analizers[0], vars_start_loc, aux_start_loc,
																 pg_JTJ[0], pg_mJTF[0], PGEnergy[0]);
		
		// P2
		vars_start_loc = vnbr;
		aux_start_loc = ninner * 10 + vnbr * 3;
		assemble_solver_pesudo_geodesic_energy_part_vertex_based(Glob_lsvars, target_angles_1, analizers[1], vars_start_loc, aux_start_loc,
																 pg_JTJ[1], pg_mJTF[1], PGEnergy[1]);
	
		vars_start_loc = vnbr * 2;
		// G
		assemble_solver_extreme_cases_part_vertex_based(Glob_lsvars, false, false, analizers[2], vars_start_loc, pg_JTJ[2], pg_mJTF[2], PGEnergy[2]);
		Compute_Auxiliaries = false;
		H += weight_pseudo_geodesic_energy * (pg_JTJ[0] + pg_JTJ[1] + weight_geodesic* pg_JTJ[2]);
		B += weight_pseudo_geodesic_energy * (pg_mJTF[0] + pg_mJTF[1] + weight_geodesic * pg_mJTF[2]);
	}
	if (fix_angle_of_two_levelsets)
	{
		// fix angle between the first (A) and the second (G) level set
		assemble_solver_fix_two_ls_angle(Glob_lsvars, angle_between_two_levelsets, faH, faB, Eangle);
		H = sum_uneven_spMats(H, weight_fix_two_ls_angle * faH);
		B = sum_uneven_vectors(B, weight_fix_two_ls_angle * faB);
	}

	H += 1e-6 * weight_mass * spMat(Eigen::VectorXd::Ones(final_size).asDiagonal());

	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver(H);

	// assert(solver.info() == Eigen::Success);
	if (solver.info() != Eigen::Success)
	{
		// solving failed
		std::cout << "solver fail" << std::endl;
		return;
	}
	// std::cout<<"solved successfully"<<std::endl;
	Eigen::VectorXd dx = solver.solve(B).eval();
	dx *= 0.75;
	// std::cout << "step length " << dx.norm() << std::endl;
	double level_set_step_length = dx.norm();
	if (level_set_step_length > max_step_length)
	{
		dx *= max_step_length / level_set_step_length;
	}
	

	// energy evaluation
	double energy_biharmonic[3];
	energy_biharmonic[0] = (QcH * func0).norm();
	energy_biharmonic[1] = (QcH * func1).norm();
	energy_biharmonic[2] = (QcH * func2).norm();
	std::cout << "energy: harm " << energy_biharmonic[0] << ", "<<energy_biharmonic[1] << ", "<<energy_biharmonic[2] << ", ";
	if (enable_pseudo_geodesic_energy)
	{
		double energy_pg[3];
		energy_pg[0]=PGEnergy[0].norm();
		energy_pg[1]=PGEnergy[1].norm();
		energy_pg[2]=PGEnergy[2].norm();
		std::cout << "PPGpg, " << energy_pg[0]<<", "<< energy_pg[1]<<", "<< energy_pg[2]<<", pgmax, "<<PGEnergy[0].lpNorm<Eigen::Infinity>()
		<<", "<<PGEnergy[1].lpNorm<Eigen::Infinity>()<<", "<<PGEnergy[2].lpNorm<Eigen::Infinity>()<<", ";
	}

	if (enable_strip_width_energy)
	{

		Eigen::VectorXd ener0, ener1, ener2; 
		ener0= GradValueF[0].rowwise().norm();
		ener1= GradValueF[1].rowwise().norm();
		ener2= GradValueF[2].rowwise().norm();
		ener0 = ener0.asDiagonal() * ener0;
		ener1 = ener1.asDiagonal() * ener1;
		ener2 = ener2.asDiagonal() * ener2;
		Eigen::VectorXd wds = Eigen::VectorXd::Ones(fnbr) * strip_width * strip_width;

		ener0 -= wds;
		ener1 -= wds;
		ener2 -= wds;
		double stp_energy[3];
		stp_energy[0] = Eigen::VectorXd(ener0).dot(ener0);
		stp_energy[1] = Eigen::VectorXd(ener1).dot(ener1);
		stp_energy[2] = Eigen::VectorXd(ener2).dot(ener2);
		std::cout << "strip, " << stp_energy[0] << ", "<< stp_energy[1] << ", "<< stp_energy[2] << ", ";
	}
	std::cout<<"extra, "<<extra_energy.norm()<<", ";
	std::cout<<"fixangle, "<<Eangle.norm()<<", ";
	step_length = dx.norm();
	std::cout << "step " << step_length << std::endl;

	func0 += dx.topRows(vnbr);
	func1 += dx.middleRows(vnbr,vnbr);
	func2 += dx.middleRows(vnbr * 2, vnbr);
	Glob_lsvars += dx;
	Last_Opt_Mesh = false;
}


// get a levelset othogonal to the reference one, or has a certain angle with the reference one.
void lsTools::Run_Othogonal_Levelset(const Eigen::VectorXd &func_ref)
{
	Eigen::MatrixXd GradValueF, GradValueV, ref_gf, ref_gv;
	Eigen::VectorXd PGEnergy, Eshading, Eangle;
	Eigen::VectorXd func = fvalues;

	int vnbr = V.rows();
	int fnbr = F.rows();

	// initialize the level set with some number
	if (func.size() != vnbr)
	{
		if (trace_hehs.size() == 0)
		{
			std::cout << "Please First Set Up The Boundary Condition (Tracing)" << std::endl;
			return;
		}
		initialize_level_set_accroding_to_parametrization();
		std::cout << "level set get initialized" << std::endl;
		return;
	}
	analysis_pseudo_geodesic_on_vertices(func, analizers[0]);
	Eigen::MatrixXd directions = get_each_face_direction(V, F, func_ref);
	int final_size = vnbr; // Change this when using more auxilary vars

	//  update quantities associated with level set values

	get_gradient_hessian_values(func, GradValueV, GradValueF);
	get_gradient_hessian_values(func_ref, ref_gv, ref_gf);
	if (Glob_lsvars.size() == 0)
	{
		std::cout << "Initializing Global Variable For LevelSet Opt ... " << std::endl;
		Glob_lsvars = Eigen::VectorXd::Zero(final_size); // We change the size if opt more than 1 level set
		Glob_lsvars.segment(0, vnbr) = func;
	}
	if (Last_Opt_Mesh)
	{
		Compute_Auxiliaries = true;
		std::cout << "Recomputing Auxiliaries" << std::endl;
	}

	spMat H;
	H.resize(vnbr, vnbr);
	Eigen::VectorXd B = Eigen::VectorXd::Zero(vnbr);

	spMat LTL;			  // left of laplacian
	Eigen::VectorXd mLTF; // right of laplacian
	assemble_solver_biharmonic_smoothing(func, LTL, mLTF);
	H += weight_laplacian * LTL;
	B += weight_laplacian * mLTF;
	assert(mass.rows() == vnbr);

	Eigen::VectorXd bcfvalue = Eigen::VectorXd::Zero(vnbr);
	Eigen::VectorXd fbdenergy = Eigen::VectorXd::Zero(vnbr);
	// boundary condition (traced as boundary condition)
	if (trace_hehs.size() > 0)
	{ // if traced, we count the related energies in

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
		assemble_solver_strip_width_part(GradValueF, sw_JTJ, sw_mJTF); // by default the strip width is 1. Unless tracing info updated the info
		H += weight_strip_width * sw_JTJ;
		B += weight_strip_width * sw_mJTF;
	}

	assert(H.rows() == vnbr);
	assert(H.cols() == vnbr);
	// pseudo geodesic
	spMat Hlarge;
	Eigen::VectorXd Blarge;
	Hlarge.resize(final_size, final_size);
	Blarge = Eigen::VectorXd::Zero(final_size);

	Hlarge = sum_uneven_spMats(H, Hlarge);
	Blarge = sum_uneven_vectors(B, Blarge);

	if (enable_pseudo_geodesic_energy)
	{

		spMat pg_JTJ, faH;
		Eigen::VectorXd pg_mJTF, faB;
		int vars_start_loc = 0;
		int aux_start_loc = vnbr;
		Eigen::VectorXi fids(fnbr);
		for (int i = 0; i < fnbr; i++)
		{
			fids[i] = i;
		}
		if(fix_angle_of_two_levelsets){
			// fix angle between the target and the refrence level set
			assemble_solver_fix_angle_to_given_face_directions(func, directions, ref_gf, angle_between_two_levelsets, fids, faH, faB, Eangle);
			Hlarge = sum_uneven_spMats(Hlarge, weight_fix_two_ls_angle * faH);
			Blarge = sum_uneven_vectors(Blarge, weight_fix_two_ls_angle * faB);
			// optimize the current levelset to make it a geodesic level set.
			assemble_solver_extreme_cases_part_vertex_based(Glob_lsvars, false, false, analizers[0], vars_start_loc, pg_JTJ, pg_mJTF, PGEnergy);
		}
		else{
			assemble_solver_othogonal_to_given_face_directions(func, directions, fids, pg_JTJ, pg_mJTF, PGEnergy);
		}
		

		Hlarge = sum_uneven_spMats(Hlarge, weight_pseudo_geodesic_energy * pg_JTJ);
		Blarge = sum_uneven_vectors(Blarge, weight_pseudo_geodesic_energy * pg_mJTF);
	}

	Hlarge += 1e-6 * weight_mass * spMat(Eigen::VectorXd::Ones(final_size).asDiagonal());

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
	if (enable_pseudo_geodesic_energy) // only pg energy and boundary angles requires step length
	{
		if (level_set_step_length > max_step_length)
		{
			dx *= max_step_length / level_set_step_length;
		}
	}

	func += dx.topRows(vnbr);
	fvalues = func;
	Glob_lsvars += dx;
	double energy_biharmonic = (QcH * func).norm();
	double energy_boundary = (bcfvalue).norm();
	std::cout << "energy: harm " << energy_biharmonic << ", bnd " << energy_boundary << ", ";
	if (enable_pseudo_geodesic_energy)
	{
		double energy_pg = PGEnergy.norm();
		double max_energy_ls = PGEnergy.lpNorm<Eigen::Infinity>();
		std::cout << "pg_energy, " << energy_pg << ", "
				  << "lsmax," << max_energy_ls << ", ";
	}
	if(fix_angle_of_two_levelsets){
		double eangle = Eangle.norm();
		std::cout << "fix_angle, " << eangle << ", ";
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
	step_length = dx.norm();
	std::cout << "step " << step_length << std::endl;

	Last_Opt_Mesh = false;
}