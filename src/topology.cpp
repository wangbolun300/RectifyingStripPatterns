#include <lsc/basic.h>
#include <lsc/tools.h>
// This version also consider about seam edges on a mesh.
// directions[0] is the inward pointing to the point, and direction[1] is outward shooting from the point. 
// handles are the two halfedges opposite to the point. 
// l2s shows if the value is from large to small along the halfedge direction. if true, then it is the inward direction.
bool find_active_faces_and_directions_around_ver_bnd(CGMesh& lsmesh, const Eigen::VectorXd& fvalues,
	const Eigen::MatrixXd& V, const int vid,
	std::vector<bool>& l2s,// if the halfedge handle is from large to small
	std::vector<CGMesh::HalfedgeHandle>& handles)
{
    l2s.clear();
    handles.clear();
    l2s.reserve(2);
    handles.reserve(2);
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
	// if (dircs.size() != 2) {// it is a singularity
	// 	return false;
	// }
	// assert(order[0] + order[1] == 1);
    if (order.size() == 2)
    {
        handles.resize(2);
        l2s.resize(2);
        for (int i = 0; i < 2; i++)
        {
            handles[order[i]] = hes[i];
            l2s[order[i]] = large_to_small[i];
        }
    }
    else{
        assert(order.size() == 1);
        handles = hes;
        l2s = large_to_small;
    }

    assert(lsmesh.face_handle(handles[0]).idx() != lsmesh.face_handle(handles[1]).idx() && "the two halfedges correspond to different faces");
	return true;

}

// consider the overlap boundary pairs
// correspondance shows the vertex the inner and outer triangles corresponds to
void lsTools::analysis_pseudo_geodesic_on_vertices(const Eigen::VectorXd &func_values,
                                                   const std::vector<std::array<int, 2>> boundary_pairs,
                                                   Eigen::VectorXi &LocalActInner,
                                                   Eigen::VectorXi &LocalActBpair,
                                                   std::vector<bool> &correspondance,
                                                   std::vector<CGMesh::HalfedgeHandle> &heh0, std::vector<CGMesh::HalfedgeHandle> &heh1,
                                                   std::vector<double> &t1s, std::vector<double> &t2s)
{
    int ninner = IVids.size();
    int npair = boundary_pairs.size();
	int vnbr = V.rows();

    LocalActInner.resize(ninner);
    LocalActBpair.resize(npair);
    correspondance.resize(npair);
    heh0.resize(ninner + npair);
    heh1.resize(ninner + npair);
    t1s.resize(ninner + npair);
    t2s.resize(ninner + npair);

    for (int i = 0; i < ninner + npair; i++)
    {

        if (i < ninner)
        {
            int vid;
            vid = IVids[i];
            std::vector<bool> l2s; // if the halfedge handle is from large to small
            std::vector<CGMesh::HalfedgeHandle> handles;
            bool active = find_active_faces_and_directions_around_ver_bnd(lsmesh, func_values, V, vid,
                                                                          l2s, handles);
            if (handles.size() == 2)
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
        else
        {
            int vid1 = BndPairs[i - ninner][0];
            int vid2 = BndPairs[i - ninner][1];

            std::vector<bool> l2s1, l2s2; // if the halfedge handle is from large to small
            std::vector<CGMesh::HalfedgeHandle> handles1, handles2;
            find_active_faces_and_directions_around_ver_bnd(lsmesh, func_values, V, vid1,
                                                                          l2s1, handles1);
            find_active_faces_and_directions_around_ver_bnd(lsmesh, func_values, V, vid2,
                                                                          l2s2, handles2);
            if (handles1.size() == 1 && handles2.size() == 1 && l2s1[0] + l2s2[0] == 1)
            {
                LocalActBpair[i - ninner] = true;
            }
            else{
                LocalActBpair[i - ninner] = false;
                continue;
            }
            CGMesh::HalfedgeHandle hd1;
            CGMesh::HalfedgeHandle hd2;
            bool normalorder = true;
            if (l2s1[0])
            {
                hd1 = handles1[0]; // the inward edge
                hd2 = handles2[0]; // the outward edge
            }
            else{
                hd1 = handles2[0]; // the inward edge
                hd2 = handles1[0]; // the outward edge
                normalorder = false;
            }

            int v1 = lsmesh.from_vertex_handle(hd1).idx();
            int v2 = lsmesh.to_vertex_handle(hd1).idx();
            int v3 = lsmesh.from_vertex_handle(hd2).idx();
            int v4 = lsmesh.to_vertex_handle(hd2).idx();
            int fid1 = lsmesh.face_handle(hd1).idx();
            int fid2 = lsmesh.face_handle(hd2).idx();

            double t1;
            double t2;
            if(normalorder){
                t1 = get_t_of_value(func_values[vid1], func_values[v1], func_values[v2]);
                t2 = get_t_of_value(func_values[vid2], func_values[v3], func_values[v4]);
            }
            else{
                t1 = get_t_of_value(func_values[vid2], func_values[v1], func_values[v2]);
                t2 = get_t_of_value(func_values[vid1], func_values[v3], func_values[v4]);
            }
            correspondance[i - ninner] = normalorder;

            assert(func_values[v1] > func_values[v2]);
            assert(func_values[v4] > func_values[v3]);
            assert(fid1 != fid2);
            heh0[i] = hd1;
            heh1[i] = hd2;
            t1s[i] = t1;
            t2s[i] = t2;
        }
    }
}
void lsTools::analysis_pseudo_geodesic_on_vertices(const Eigen::VectorXd &func_values, const std::vector<std::array<int, 2>> boundary_pairs, LSAnalizer &analizer)
{
    analysis_pseudo_geodesic_on_vertices(func_values, boundary_pairs, analizer.LocalActInner,
                                         analizer.LocalActBpair, analizer.Correspondance, analizer.heh0, analizer.heh1,
                                         analizer.t1s, analizer.t2s);
}
void lsTools::assemble_solver_pesudo_geodesic_energy_part_stitch_boundary(Eigen::VectorXd &vars, const std::vector<double> &angle_degree,
                                                                          const LSAnalizer &analizer, const int vars_start_loc, const int aux_start_loc, spMat &H, 
                                                                          Eigen::VectorXd &B, Eigen::VectorXd &Energy)
{
    double cos_angle, sin_angle;
	int vnbr = V.rows();

    int ninner = analizer.LocalActInner.size();
    int nbpair = analizer.LocalActBpair.size();
	if (angle_degree.size() == 1) {
		double angle_radian = angle_degree[0] * LSC_PI / 180.; // the angle in radian
		cos_angle = cos(angle_radian);
		sin_angle = sin(angle_radian);
	}
	// if (Analyze_Optimized_LS_Angles)
	// {
	// 	// std::cout<<"Check 1"<<std::endl;
	// 	AnalizedAngelVector.resize(ninner);
	// }

	std::vector<Trip> tripletes;
	tripletes.reserve(nbpair * 65); 
	Energy = Eigen::VectorXd::Zero(nbpair * 12);// mesh total energy values

	assert(angle_degree.size() == 1 || angle_degree.size() == vnbr);
	for (int i = 0; i < nbpair; i++)
	{
        int vm0, vm1;
        // make sure vm0 is on the inward side, vm1 is on the outward side.
        if (analizer.Correspondance[i]) 
        {
            vm0 = BndPairs[i][0];
            vm1 = BndPairs[i][1];
        }
        else{
            vm0 = BndPairs[i][1];
            vm1 = BndPairs[i][0];
        }

		if (angle_degree.size() == vnbr) // in this case, the angle on vm0 should be the same as vm1
		{
			double angle_radian = angle_degree[vm0] * LSC_PI / 180.; // the angle in radian
			cos_angle = cos(angle_radian);
			sin_angle = sin(angle_radian);
		}


		CGMesh::HalfedgeHandle inhd = analizer.heh0[i], outhd = analizer.heh1[i];
		int v1 = lsmesh.from_vertex_handle(inhd).idx();
		int v2 = lsmesh.to_vertex_handle(inhd).idx();
		int v3 = lsmesh.from_vertex_handle(outhd).idx();
		int v4 = lsmesh.to_vertex_handle(outhd).idx();
		
		double t1 = analizer.t1s[i];
		double t2 = analizer.t2s[i];
		
		Eigen::Vector3d ver0 = V.row(v1) + (V.row(v2) - V.row(v1)) * t1;
		Eigen::Vector3d ver1 = V.row(vm0);
		Eigen::Vector3d ver2 = V.row(v3) + (V.row(v4) - V.row(v3)) * t2;
		
		// the locations
		int lvm0 = vars_start_loc + vm0;
        int lvm1 = vars_start_loc + vm1;
		int lv1 = vars_start_loc + v1;
		int lv2 = vars_start_loc + v2;
		int lv3 = vars_start_loc + v3;
		int lv4 = vars_start_loc + v4;

		int lrx = i + aux_start_loc;
		int lry = i + aux_start_loc + nbpair;
		int lrz = i + aux_start_loc + nbpair * 2;
		int lux = i + aux_start_loc + nbpair * 3;
		int luy = i + aux_start_loc + nbpair * 4;
		int luz = i + aux_start_loc + nbpair * 5;
		int lsx = i + aux_start_loc + nbpair * 6;
		int lsy = i + aux_start_loc + nbpair * 7;
		int lsz = i + aux_start_loc + nbpair * 8;
		int lh = i + aux_start_loc + nbpair * 9;
        // the normal vector is regarded as the average of both the boundary normals.
        Eigen::Vector3d norm = (norm_v.row(vm0) + norm_v.row(vm1)).normalized();
        Eigen::Vector3d v31 = V.row(v3) - V.row(v1);
		Eigen::Vector3d v43 = V.row(v4) - V.row(v3);
		Eigen::Vector3d v21 = V.row(v2) - V.row(v1);
		double f1 = vars(lv1);
		double f2 = vars(lv2);
		double f3 = vars(lv3);
		double f4 = vars(lv4);
		double fm0 = vars(lvm0);
        double fm1 = vars(lvm1);
		
		Eigen::Vector3d real_u = norm.cross(ver2 - ver0);
		real_u = real_u.normalized();
		if (Compute_Auxiliaries) {
			// std::cout<<"init auxiliaries"<<std::endl;
			Eigen::Vector3d real_s = v31 * (f4 - f3) * (f2 - f1) + v43 * (fm1 - f3) * (f2 - f1) - v21 * (fm0 - f1) * (f4 - f3);
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
		Binormals.row(vm0) = r.dot(norm) < 0 ? -r : r;
        Binormals.row(vm1) = Binormals.row(vm0);
		// the weights
		double dis0 = ((V.row(v1) - V.row(v2)) * vars[lvm0] + (V.row(v2) - V.row(vm0)) * vars[lv1] + (V.row(vm0) - V.row(v1)) * vars[lv2]).norm();
		double dis1 = ((V.row(v3) - V.row(v4)) * vars[lvm1] + (V.row(v4) - V.row(vm1)) * vars[lv3] + (V.row(vm1) - V.row(v3)) * vars[lv4]).norm();

		assert(dis0 != 0 && dis1 != 0);
		
        double scale0 = mass_uniform.coeff(vm0, vm0);
        double scale1 = mass_uniform.coeff(vm1, vm1);

		double scale = scale0 + scale1;
		double EhanceBinormal = 1;
		// if (Analyze_Optimized_LS_Angles && !Use_Fitting_Angles)
		// {
		// 	EhanceBinormal = 100;
		// }
		assert(scale != 0);
		// r dot (vm+(t1-1)*vf-t1*vt)
		// vf = v1, vt = v2
		tripletes.push_back(Trip(i, lrx, ((V(v1, 0) - V(v2, 0)) * vars[lvm0] + (V(v2, 0) - V(vm0, 0)) * vars[lv1] + (V(vm0, 0) - V(v1, 0)) * vars[lv2]) / dis0 * scale * EhanceBinormal));
		tripletes.push_back(Trip(i, lry, ((V(v1, 1) - V(v2, 1)) * vars[lvm0] + (V(v2, 1) - V(vm0, 1)) * vars[lv1] + (V(vm0, 1) - V(v1, 1)) * vars[lv2]) / dis0 * scale * EhanceBinormal));
		tripletes.push_back(Trip(i, lrz, ((V(v1, 2) - V(v2, 2)) * vars[lvm0] + (V(v2, 2) - V(vm0, 2)) * vars[lv1] + (V(vm0, 2) - V(v1, 2)) * vars[lv2]) / dis0 * scale * EhanceBinormal));

		double r12 = (V.row(v1) - V.row(v2)).dot(r);
		double rm1 = (V.row(vm0) - V.row(v1)).dot(r);
		double r2m = (V.row(v2) - V.row(vm0)).dot(r);
		tripletes.push_back(Trip(i, lvm0, r12 / dis0 * scale * EhanceBinormal));
		tripletes.push_back(Trip(i, lv1, r2m / dis0 * scale * EhanceBinormal));
		tripletes.push_back(Trip(i, lv2, rm1 / dis0 * scale * EhanceBinormal));
		Energy[i] = (r12 * vars[lvm0] + rm1 * vars[lv2] + r2m * vars[lv1]) / dis0 * scale * EhanceBinormal;

		// vf = v3, vt = v4
		tripletes.push_back(Trip(i + nbpair, lrx, ((V(v3, 0) - V(v4, 0)) * vars[lvm1] + (V(v4, 0) - V(vm1, 0)) * vars[lv3] + (V(vm1, 0) - V(v3, 0)) * vars[lv4]) / dis1 * scale * EhanceBinormal));
		tripletes.push_back(Trip(i + nbpair, lry, ((V(v3, 1) - V(v4, 1)) * vars[lvm1] + (V(v4, 1) - V(vm1, 1)) * vars[lv3] + (V(vm1, 1) - V(v3, 1)) * vars[lv4]) / dis1 * scale * EhanceBinormal));
		tripletes.push_back(Trip(i + nbpair, lrz, ((V(v3, 2) - V(v4, 2)) * vars[lvm1] + (V(v4, 2) - V(vm1, 2)) * vars[lv3] + (V(vm1, 2) - V(v3, 2)) * vars[lv4]) / dis1 * scale * EhanceBinormal));

		r12 = (V.row(v3) - V.row(v4)).dot(r);
		rm1 = (V.row(vm1) - V.row(v3)).dot(r);
		r2m = (V.row(v4) - V.row(vm1)).dot(r);
		tripletes.push_back(Trip(i + nbpair, lvm1, r12 / dis1 * scale * EhanceBinormal));
		tripletes.push_back(Trip(i + nbpair, lv3, r2m / dis1 * scale * EhanceBinormal));
		tripletes.push_back(Trip(i + nbpair, lv4, rm1 / dis1 * scale * EhanceBinormal));

		Energy[i + nbpair] = (r12 * vars[lvm1] + rm1 * vars[lv4] + r2m * vars[lv3]) / dis1 * scale * EhanceBinormal;

		// r*r=1
		tripletes.push_back(Trip(i + nbpair * 2, lrx, 2 * vars(lrx) * scale * EhanceBinormal));
		tripletes.push_back(Trip(i + nbpair * 2, lry, 2 * vars(lry) * scale * EhanceBinormal));
		tripletes.push_back(Trip(i + nbpair * 2, lrz, 2 * vars(lrz) * scale * EhanceBinormal));

		Energy[i + nbpair * 2] = (r.dot(r) - 1) * scale * EhanceBinormal;
		// if (Analyze_Optimized_LS_Angles)
		// {
		// 	Eigen::Vector3d nu = u.normalized();
		// 	Eigen::Vector3d nn = norm;
		// 	Eigen::Vector3d nb = r.normalized();
		// 	double cosina = nb.dot(nn);
		// 	double sina = nb.dot(nu);
		// 	if (sina < 0)
        //     {
        //         cosina *= -1;
        //     }
        //     double angle = acos(cosina);
		// 	AnalizedAngelVector[i] = angle;
		// 	if(Use_Opt_Only_BNMS){
		// 		continue;
		// 	}
		// }
		// (r*norm)^2 - cos^2 = 0

		tripletes.push_back(Trip(i + nbpair * 3, lrx,
								 (2 * norm(0) * norm(0) * vars(lrx) + 2 * norm(0) * norm(1) * vars(lry) + 2 * norm(0) * norm(2) * vars(lrz)) * scale));
		tripletes.push_back(Trip(i + nbpair * 3, lry,
								 (2 * norm(1) * norm(1) * vars(lry) + 2 * norm(0) * norm(1) * vars(lrx) + 2 * norm(2) * norm(1) * vars(lrz)) * scale));
		tripletes.push_back(Trip(i + nbpair * 3, lrz,
								 (2 * norm(2) * norm(2) * vars(lrz) + 2 * norm(0) * norm(2) * vars(lrx) + 2 * norm(1) * norm(2) * vars(lry)) * scale));

		double ndr = norm.dot(r);
		Energy[i + nbpair * 3] = (ndr * ndr - cos_angle * cos_angle) * scale;

		// u*u=1
		
		tripletes.push_back(Trip(i + nbpair * 4, lux, 2 * vars(lux) * scale));
		tripletes.push_back(Trip(i + nbpair * 4, luy, 2 * vars(luy) * scale));
		tripletes.push_back(Trip(i + nbpair * 4, luz, 2 * vars(luz) * scale));

		Energy[i + nbpair * 4] = (u.dot(u) - 1) * scale;

		// u * norm = 0
		tripletes.push_back(Trip(i + nbpair * 5, lux, norm(0) * scale));
		tripletes.push_back(Trip(i + nbpair * 5, luy, norm(1) * scale));
		tripletes.push_back(Trip(i + nbpair * 5, luz, norm(2) * scale));

		Energy[i + nbpair * 5] = norm.dot(u) * scale;

		// u*s = 0

		tripletes.push_back(Trip(i + nbpair * 6, lux, s(0) / ns * scale));
		tripletes.push_back(Trip(i + nbpair * 6, luy, s(1) / ns * scale));
		tripletes.push_back(Trip(i + nbpair * 6, luz, s(2) / ns * scale));

		tripletes.push_back(Trip(i + nbpair * 6, lsx, u(0) / ns * scale));
		tripletes.push_back(Trip(i + nbpair * 6, lsy, u(1) / ns * scale));
		tripletes.push_back(Trip(i + nbpair * 6, lsz, u(2) / ns * scale));

		Energy[i + nbpair * 6] = s.dot(u) / ns * scale;
		assert(ns != 0);
		// xxx-s = 0
		
		// x axis
		tripletes.push_back(Trip(i + nbpair * 7, lsx, -1 * scale));
		tripletes.push_back(Trip(i + nbpair * 7, lv1, (-v31(0) * (f4 - f3) - v43(0) * (fm1 - f3) + v21(0) * (f4 - f3)) * scale));
		tripletes.push_back(Trip(i + nbpair * 7, lv2, (v31(0) * (f4 - f3) + v43(0) * (fm1 - f3)) * scale));
		tripletes.push_back(Trip(i + nbpair * 7, lv3, (-v31(0) * (f2 - f1) - v43(0) * (f2 - f1) + v21(0) * (fm0 - f1)) * scale));
		tripletes.push_back(Trip(i + nbpair * 7, lv4, (v31(0) * (f2 - f1) - v21(0) * (fm0 - f1)) * scale));
		tripletes.push_back(Trip(i + nbpair * 7, lvm0, -v21(0) * (f4 - f3) * scale));
        tripletes.push_back(Trip(i + nbpair * 7, lvm1, v43(0) * (f2 - f1) * scale));
		Energy[i + nbpair * 7] = (v31(0) * (f4 - f3) * (f2 - f1) + v43(0) * (fm1 - f3) * (f2 - f1) - v21(0) * (fm0 - f1) * (f4 - f3) - s(0)) * scale;

		// y axis
		tripletes.push_back(Trip(i + nbpair * 8, lsy, -1 * scale));
		tripletes.push_back(Trip(i + nbpair * 8, lv1, (-v31(1) * (f4 - f3) - v43(1) * (fm1 - f3) + v21(1) * (f4 - f3)) * scale));
		tripletes.push_back(Trip(i + nbpair * 8, lv2, (v31(1) * (f4 - f3) + v43(1) * (fm1 - f3)) * scale));
		tripletes.push_back(Trip(i + nbpair * 8, lv3, (-v31(1) * (f2 - f1) - v43(1) * (f2 - f1) + v21(1) * (fm0 - f1)) * scale));
		tripletes.push_back(Trip(i + nbpair * 8, lv4, (v31(1) * (f2 - f1) - v21(1) * (fm0 - f1)) * scale));
		tripletes.push_back(Trip(i + nbpair * 8, lvm0, -v21(1) * (f4 - f3) * scale));
        tripletes.push_back(Trip(i + nbpair * 8, lvm1, v43(1) * (f2 - f1) * scale));
		Energy[i + nbpair * 8] = (v31(1) * (f4 - f3) * (f2 - f1) + v43(1) * (fm1 - f3) * (f2 - f1) - v21(1) * (fm0 - f1) * (f4 - f3) - s(1)) * scale;

		// z axis
		tripletes.push_back(Trip(i + nbpair * 9, lsz, -1 * scale));
		tripletes.push_back(Trip(i + nbpair * 9, lv1, (-v31(2) * (f4 - f3) - v43(2) * (fm1 - f3) + v21(2) * (f4 - f3)) * scale));
		tripletes.push_back(Trip(i + nbpair * 9, lv2, (v31(2) * (f4 - f3) + v43(2) * (fm1 - f3)) * scale));
		tripletes.push_back(Trip(i + nbpair * 9, lv3, (-v31(2) * (f2 - f1) - v43(2) * (f2 - f1) + v21(2) * (fm0 - f1)) * scale));
		tripletes.push_back(Trip(i + nbpair * 9, lv4, (v31(2) * (f2 - f1) - v21(2) * (fm0 - f1)) * scale));
		tripletes.push_back(Trip(i + nbpair * 9, lvm0, -v21(2) * (f4 - f3) * scale));
        tripletes.push_back(Trip(i + nbpair * 9, lvm1, v43(2) * (f2 - f1) * scale));
		Energy[i + nbpair * 9] = (v31(2) * (f4 - f3) * (f2 - f1) + v43(2) * (fm1 - f3) * (f2 - f1) - v21(2) * (fm0 - f1) * (f4 - f3) - s(2)) * scale;

		// r*u - h =0
		tripletes.push_back(Trip(i + nbpair * 10, lrx, u(0) * scale));
		tripletes.push_back(Trip(i + nbpair * 10, lry, u(1) * scale));
		tripletes.push_back(Trip(i + nbpair * 10, lrz, u(2) * scale));

		tripletes.push_back(Trip(i + nbpair * 10, lux, r(0) * scale));
		tripletes.push_back(Trip(i + nbpair * 10, luy, r(1) * scale));
		tripletes.push_back(Trip(i + nbpair * 10, luz, r(2) * scale));

		tripletes.push_back(Trip(i + nbpair * 10, lh, -1 * scale));

		Energy[i + nbpair * 10] = (r.dot(u) - vars[lh]) * scale;

		// h times (r.dot(n)) - sin*cos = 0
		double h = vars[lh];
		tripletes.push_back(Trip(i + nbpair * 11, lh, r.dot(norm) * scale));
		tripletes.push_back(Trip(i + nbpair * 11, lrx, norm(0) * h * scale));
		tripletes.push_back(Trip(i + nbpair * 11, lry, norm(1) * h * scale));
		tripletes.push_back(Trip(i + nbpair * 11, lrz, norm(2) * h * scale));

		Energy[i + nbpair * 11] = (r.dot(norm) * h - sin_angle * cos_angle) * scale;

		//
		
	}

	int nvars = vars.size();
	int ncondi = Energy.size();
	// int ninner = analizer.LocalActInner.size();
	// int rep = 
	spMat J;
	J.resize(ncondi, nvars);
	J.setFromTriplets(tripletes.begin(), tripletes.end());
	H = J.transpose() * J;
	B = -J.transpose() * Energy;
	PGE = Energy;
}