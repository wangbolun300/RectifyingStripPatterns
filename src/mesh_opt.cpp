#include<lsc/mesh_opt.h>
#include <lsc/tools.h>
#include <igl/grad.h>
#include <igl/hessian.h>
#include <igl/curved_hessian_energy.h>
#include <igl/repdiag.h>
#include <igl/Timer.h>


bool vector_contains_NAN(Eigen::VectorXd& B){
    for(int i=0;i<B.size();i++){
        if(isnan(B[i])){
            return true;
        }
    }
    return false;
}
// the variables are sorted as v0x, v1x, ...vnx, v0y, v1y,...,vny, v0z, v1z, ..., vnz
// xyz is 0, 1, or 2
void location_in_sparse_matrix(const int vnbr, const int vderivate, const int dxyz, const int vcoff, const int cxyz, int& row, int &col){
    int dlocation = vderivate + dxyz * vnbr;
    int clocation = vcoff + cxyz * vnbr;
    row=dlocation;
    col=clocation;
}
void location_in_sparse_matrix(const int vnbr, const int vderivate, const int dxyz, int &col){
    col = vderivate + dxyz * vnbr;
}
// the coefficient matrix of the cross pair.
// the results are 3 sparse matrices M0, M1, M2, where alpha* (nx*M0+ny*M1+nz*M2)*fvalues is the jacobian of 
// the term alpha*(V(vid0) X V(vid1)).dot(norm)
void get_derivate_coeff_cross_pair(const int vid0, const int vid1, const int vsize, const double alpha, std::array<spMat, 3> &cmats)
{
    cmats[0].resize(vsize*3, vsize*3);
    cmats[1].resize(vsize*3, vsize*3);
    cmats[2].resize(vsize*3, vsize*3);
    std::vector<Trip> triplets;
    int row, col;
    triplets.resize(4);

    // M0
    location_in_sparse_matrix(vsize,vid0,1,vid1,2,row,col);// y0z1
    triplets[0]=Trip(row,col,alpha);
    location_in_sparse_matrix(vsize,vid1,2,vid0,1,row,col);// y0z1
    triplets[1]=Trip(row,col,alpha);
    location_in_sparse_matrix(vsize,vid1,1,vid0,2,row,col);// -y1z0
    triplets[2]=Trip(row,col,-alpha);
    location_in_sparse_matrix(vsize,vid0,2,vid1,1,row,col);// -y1z0
    triplets[3]=Trip(row,col,-alpha);
    cmats[0].setFromTriplets(triplets.begin(),triplets.end());

     // M1
    location_in_sparse_matrix(vsize,vid1,0,vid0,2,row,col);// x1z0
    triplets[0]=Trip(row,col,alpha);
    location_in_sparse_matrix(vsize,vid0,2,vid1,0,row,col);// x1z0
    triplets[1]=Trip(row,col,alpha);
    location_in_sparse_matrix(vsize,vid0,0,vid1,2,row,col);// -x0z1
    triplets[2]=Trip(row,col,-alpha);
    location_in_sparse_matrix(vsize,vid1,2,vid0,0,row,col);// -x0z1
    triplets[3]=Trip(row,col,-alpha);
    cmats[1].setFromTriplets(triplets.begin(),triplets.end());

     // M2
    location_in_sparse_matrix(vsize,vid0,0,vid1,1,row,col);// x0y1
    triplets[0]=Trip(row,col,alpha);
    location_in_sparse_matrix(vsize,vid1,1,vid0,0,row,col);// x0y1
    triplets[1]=Trip(row,col,alpha);
    location_in_sparse_matrix(vsize,vid1,0,vid0,1,row,col);// -x1y0
    triplets[2]=Trip(row,col,-alpha);
    location_in_sparse_matrix(vsize,vid0,1,vid1,0,row,col);// -x1y0
    triplets[3]=Trip(row,col,-alpha);
    cmats[2].setFromTriplets(triplets.begin(),triplets.end());
}

// for each pair of cross product v1 x v2, alpha * v1 x v2.dot(norm) is the form of energy
void get_alpha_associate_with_cross_pairs(const double t1, const double t2, std::array<double,8> &result)
{
    result[0] = 1 - t1;
    result[1] = t1;
    result[2] = 1 - t2;
    result[3] = t2;
    result[4] = t1 * (t2 - 1);
    result[5] = -t1 * t2;
    result[6] = (t1 - 1) * (1 - t2);
    result[7] = t2 * (t1 - 1);
}

// the size of vars should be nvars, if there is no auxiliary variables, nvars = vnbr*3.
// The auxiliary vars:
// r: the bi-normals. located from aux_start_loc to ninner * 3 + aux_start_loc
// u: the side vector. located from aux_start_loc + ninner * 3 to aux_start_loc + ninner * 6
// h: r.dot(u). located from aux_start_loc + ninner * 6 + aux_start_loc + ninner * 7
void lsTools::calculate_mesh_opt_expanded_function_values(Eigen::VectorXd &vars,
                                                          const LSAnalizer &analizer,
                                                          const std::vector<double> &angle_degree,
                                                          const int aux_start_loc, std::vector<Trip> &tripletes, Eigen::VectorXd &MTenergy)
{

    double cos_angle;
    double sin_angle;
    if (angle_degree.size() == 1) {
        double angle_radian = angle_degree[0] * LSC_PI / 180.; // the angle in radian
        cos_angle = cos(angle_radian);
        sin_angle = sin(angle_radian);
    }
    
    int ninner = analizer.LocalActInner.size();
    int vnbr = V.rows();
    if (Binormals.rows() != vnbr)
    {
        Binormals = Eigen::MatrixXd::Zero(vnbr, 3);
    }
    //PeWeight = Eigen::VectorXd::Zero(ninner);
    //lens = Eigen::VectorXd::Zero(ninner);
    tripletes.clear();
    tripletes.reserve(ninner * 70);// the number of rows is ninner*4, the number of cols is aux_start_loc + ninner * 3 (all the vertices and auxiliary vars)
	MTenergy = Eigen::VectorXd::Zero(ninner * 9); // mesh total energy values

    for (int i = 0; i < ninner; i++)
    {
        if (analizer.LocalActInner[i] == false) {
            continue;
        }
        int vm = IVids[i];
        double scale = mass_uniform.coeff(vm, vm);
        if (angle_degree.size() == vnbr) {
            double angle_radian = angle_degree[vm] * LSC_PI / 180.; // the angle in radian
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
        Eigen::Vector3d ver1 = V.row(vm);
        Eigen::Vector3d ver2 = V.row(v3) + (V.row(v4) - V.row(v3)) * t2;
        // the locations
        int lrx = i + aux_start_loc;
        int lry = i + aux_start_loc + ninner;
        int lrz = i + aux_start_loc + ninner * 2;
        int lux = i + aux_start_loc + ninner * 3;
        int luy = i + aux_start_loc + ninner * 4;
        int luz = i + aux_start_loc + ninner * 5;
        int lh = i + aux_start_loc + ninner * 6;
        int l1x = v1;
        int l1y = v1 + vnbr;
        int l1z = v1 + vnbr * 2;
        int l2x = v2;
        int l2y = v2 + vnbr;
        int l2z = v2 + vnbr * 2;
        int l3x = v3;
        int l3y = v3 + vnbr;
        int l3z = v3 + vnbr * 2;
        int l4x = v4;
        int l4y = v4 + vnbr;
        int l4z = v4 + vnbr * 2;
        int lmx = vm;
        int lmy = vm + vnbr;
        int lmz = vm + vnbr * 2;
        Eigen::Vector3d norm = norm_v.row(vm);
        Eigen::Vector3d real_u = norm.cross(ver2 - ver0);
        real_u = real_u.normalized();
        if (Compute_Auxiliaries_Mesh) {
            Eigen::Vector3d real_r = (ver1 - ver0).cross(ver2 - ver1);
            real_r = real_r.normalized();
            vars[lrx] = real_r[0];
            vars[lry] = real_r[1];
            vars[lrz] = real_r[2];
            
            vars[lux] = real_u[0];
            vars[luy] = real_u[1];
            vars[luz] = real_u[2];

            vars[lh] = real_r.dot(real_u);
        }
        Eigen::Vector3d r = Eigen::Vector3d(vars[lrx], vars[lry], vars[lrz]);
        Eigen::Vector3d u = Eigen::Vector3d(vars[lux], vars[luy], vars[luz]);
        Binormals.row(vm) = r.dot(norm) < 0 ? -r : r;
        double h = vars[lh];
        // r is allowed to flip to the opposite position, but u is not allowed, since it is the side vector.
        if (u.dot(real_u) < 0)
		{
			vars(lux) *= -1;
			vars(luy) *= -1;
			vars(luz) *= -1;
			u *= -1;
		}

        int lfx = v1;
        int lfy = v1 + vnbr;
        int lfz = v1 + vnbr * 2;
        int ltx = v2;
        int lty = v2 + vnbr;
        int ltz = v2 + vnbr * 2;

        double d0 = (ver0 - ver1).norm();
        double d1 = (ver2 - ver1).norm();
        double qp = (ver2 - ver0).norm();
        // r dot (vm+(t1-1)*vf-t1*vt)
        // vf = v1, vt = v2
        tripletes.push_back(Trip(i, lrx, (vars[lmx] + (t1 - 1) * vars[lfx] - t1 * vars[ltx]) / d0 * scale));
        tripletes.push_back(Trip(i, lry, (vars[lmy] + (t1 - 1) * vars[lfy] - t1 * vars[lty]) / d0 * scale));
        tripletes.push_back(Trip(i, lrz, (vars[lmz] + (t1 - 1) * vars[lfz] - t1 * vars[ltz]) / d0 * scale));

        tripletes.push_back(Trip(i, lmx, vars[lrx] / d0 * scale));
        tripletes.push_back(Trip(i, lmy, vars[lry] / d0 * scale));
        tripletes.push_back(Trip(i, lmz, vars[lrz] / d0 * scale));

        tripletes.push_back(Trip(i, lfx, (t1 - 1) * vars[lrx] / d0 * scale));
        tripletes.push_back(Trip(i, lfy, (t1 - 1) * vars[lry] / d0 * scale));
        tripletes.push_back(Trip(i, lfz, (t1 - 1) * vars[lrz] / d0 * scale));

        tripletes.push_back(Trip(i, ltx, -t1 * vars[lrx] / d0 * scale));
        tripletes.push_back(Trip(i, lty, -t1 * vars[lry] / d0 * scale));
        tripletes.push_back(Trip(i, ltz, -t1 * vars[lrz] / d0 * scale));

        MTenergy[i] = r.dot(ver1 - ver0) / d0 * scale;

        // vf = v3, vt = v4
        lfx = v3;
        lfy = v3 + vnbr;
        lfz = v3 + vnbr * 2;
        ltx = v4;
        lty = v4 + vnbr;
        ltz = v4 + vnbr * 2;

        tripletes.push_back(Trip(i + ninner, lrx, (vars[lmx] + (t2 - 1) * vars[lfx] - t2 * vars[ltx]) / d1 * scale));
        tripletes.push_back(Trip(i + ninner, lry, (vars[lmy] + (t2 - 1) * vars[lfy] - t2 * vars[lty]) / d1 * scale));
        tripletes.push_back(Trip(i + ninner, lrz, (vars[lmz] + (t2 - 1) * vars[lfz] - t2 * vars[ltz]) / d1 * scale));

        tripletes.push_back(Trip(i + ninner, lmx, vars[lrx] / d1 * scale));
        tripletes.push_back(Trip(i + ninner, lmy, vars[lry] / d1 * scale));
        tripletes.push_back(Trip(i + ninner, lmz, vars[lrz] / d1 * scale));

        tripletes.push_back(Trip(i + ninner, lfx, (t2 - 1) * vars[lrx] / d1 * scale));
        tripletes.push_back(Trip(i + ninner, lfy, (t2 - 1) * vars[lry] / d1 * scale));
        tripletes.push_back(Trip(i + ninner, lfz, (t2 - 1) * vars[lrz] / d1 * scale));

        tripletes.push_back(Trip(i + ninner, ltx, -t2 * vars[lrx] / d1 * scale));
        tripletes.push_back(Trip(i + ninner, lty, -t2 * vars[lry] / d1 * scale));
        tripletes.push_back(Trip(i + ninner, ltz, -t2 * vars[lrz] / d1 * scale));

        MTenergy[i + ninner] = r.dot(ver1 - ver2) / d1 * scale;

        // r*r=1
        tripletes.push_back(Trip(i + ninner * 2, lrx, 2 * vars(lrx) * scale));
        tripletes.push_back(Trip(i + ninner * 2, lry, 2 * vars(lry) * scale));
        tripletes.push_back(Trip(i + ninner * 2, lrz, 2 * vars(lrz) * scale));

        MTenergy[i + ninner * 2] = (r.dot(r) - 1) * scale;

        // (r*norm)^2 - cos^2 = 0

        tripletes.push_back(Trip(i + ninner * 3, lrx,
                                 (2 * norm(0) * norm(0) * vars(lrx) + 2 * norm(0) * norm(1) * vars(lry) + 2 * norm(0) * norm(2) * vars(lrz)) * scale));
        tripletes.push_back(Trip(i + ninner * 3, lry,
                                 (2 * norm(1) * norm(1) * vars(lry) + 2 * norm(0) * norm(1) * vars(lrx) + 2 * norm(2) * norm(1) * vars(lrz)) * scale));
        tripletes.push_back(Trip(i + ninner * 3, lrz,
                                 (2 * norm(2) * norm(2) * vars(lrz) + 2 * norm(0) * norm(2) * vars(lrx) + 2 * norm(1) * norm(2) * vars(lry)) * scale));
        double nr = norm.dot(r);

        MTenergy[i + ninner * 3] = (nr * nr - cos_angle * cos_angle) * scale;

        // u * u = 1

        tripletes.push_back(Trip(i + ninner * 4, lux, 2 * vars(lux) * scale));
        tripletes.push_back(Trip(i + ninner * 4, luy, 2 * vars(luy) * scale));
        tripletes.push_back(Trip(i + ninner * 4, luz, 2 * vars(luz) * scale));

        MTenergy[i + ninner * 4] = (u.dot(u) - 1) * scale;

        // u * norm = 0
        tripletes.push_back(Trip(i + ninner * 5, lux, norm(0) * scale));
        tripletes.push_back(Trip(i + ninner * 5, luy, norm(1) * scale));
        tripletes.push_back(Trip(i + ninner * 5, luz, norm(2) * scale));

        MTenergy[i + ninner * 5] = norm.dot(u) * scale;

        // u * (ver2-ver0) = 0
        double c1 = t1 - 1;
        double c2 = -t1;
        double c3 = 1 - t2;
        double c4 = t2;

        tripletes.push_back(Trip(i + ninner * 6, lux, (c1 * vars[l1x] + c2 * vars[l2x] + c3 * vars[l3x] + c4 * vars[l4x]) / qp * scale));
        tripletes.push_back(Trip(i + ninner * 6, luy, (c1 * vars[l1y] + c2 * vars[l2y] + c3 * vars[l3y] + c4 * vars[l4y]) / qp * scale));
        tripletes.push_back(Trip(i + ninner * 6, luz, (c1 * vars[l1z] + c2 * vars[l2z] + c3 * vars[l3z] + c4 * vars[l4z]) / qp * scale));

        tripletes.push_back(Trip(i + ninner * 6, l1x, c1 * u[0] / qp * scale));
        tripletes.push_back(Trip(i + ninner * 6, l2x, c2 * u[0] / qp * scale));
        tripletes.push_back(Trip(i + ninner * 6, l3x, c3 * u[0] / qp * scale));
        tripletes.push_back(Trip(i + ninner * 6, l4x, c4 * u[0] / qp * scale));

        tripletes.push_back(Trip(i + ninner * 6, l1y, c1 * u[1] / qp * scale));
        tripletes.push_back(Trip(i + ninner * 6, l2y, c2 * u[1] / qp * scale));
        tripletes.push_back(Trip(i + ninner * 6, l3y, c3 * u[1] / qp * scale));
        tripletes.push_back(Trip(i + ninner * 6, l4y, c4 * u[1] / qp * scale));

        tripletes.push_back(Trip(i + ninner * 6, l1z, c1 * u[2] / qp * scale));
        tripletes.push_back(Trip(i + ninner * 6, l2z, c2 * u[2] / qp * scale));
        tripletes.push_back(Trip(i + ninner * 6, l3z, c3 * u[2] / qp * scale));
        tripletes.push_back(Trip(i + ninner * 6, l4z, c4 * u[2] / qp * scale));

        MTenergy[i + ninner * 6] = u.dot(ver2 - ver0) / qp * scale;

        // r * u - h =0
        tripletes.push_back(Trip(i + ninner * 7, lrx, u[0] * scale));
        tripletes.push_back(Trip(i + ninner * 7, lry, u[1] * scale));
        tripletes.push_back(Trip(i + ninner * 7, lrz, u[2] * scale));

        tripletes.push_back(Trip(i + ninner * 7, lux, r[0] * scale));
        tripletes.push_back(Trip(i + ninner * 7, luy, r[1] * scale));
        tripletes.push_back(Trip(i + ninner * 7, luz, r[2] * scale));

        tripletes.push_back(Trip(i + ninner * 7, lh, -1 * scale));

        MTenergy[i + ninner * 7] = (r.dot(u) - h) * scale;

        // r.dot(n)*h = sin * cos
        tripletes.push_back(Trip(i + ninner * 8, lrx, h * norm[0] * scale));
        tripletes.push_back(Trip(i + ninner * 8, lry, h * norm[1] * scale));
        tripletes.push_back(Trip(i + ninner * 8, lrz, h * norm[2] * scale));

        tripletes.push_back(Trip(i + ninner * 8, lh, r.dot(norm) * scale));

        MTenergy[i + ninner * 8] = (r.dot(norm) * h - sin_angle * cos_angle) * scale;
    }
}
// L.cross(R).dot(n)
// pos shows the var is L or R. cor indicates if it is x y or z
double LxRdN_differential(const int pos, const Eigen::Vector3d& a, const Eigen::Vector3d& b, const Eigen::Vector3d& n, const int cor){
    if (pos == 0)
    { // a is the variable
        if (cor == 0)
        { // to x
            return -b(2) * n(1) + b(1) * n(2);
        }
        if (cor == 1)
        {
            return b(2) * n(0) - b(0) * n(2);
        }
        if (cor == 2)
        {
            return -b(1) * n(0) + b(0) * n(1);
        }
        std::cout << "LxRdN_differential ERROR" << std::endl;
    }
    if (pos == 1)
    { // b is the variable
        if (cor == 0)
        { // to x
            return a(2) * n(1) - a(1) * n(2);
        }
        if (cor == 1)
        {
            return -a(2) * n(0) + a(0) * n(2);
        }
        if (cor == 2)
        {
            return a(1) * n(0) - a(0) * n(1);
        }
        std::cout << "LxRdN_differential ERROR" << std::endl;
    }
    std::cout << "LxRdN_differential ERROR" << std::endl;
}
void lsTools::calculate_mesh_opt_extreme_values(Eigen::VectorXd &vars, const int aux_start_loc, const Eigen::VectorXd &func, const bool asymptotic, const bool use_given_direction, const Eigen::Vector3d &ray,
                                                const LSAnalizer &analizer, std::vector<Trip> &tripletes, Eigen::VectorXd &MTenergy)
{

    int ninner = analizer.LocalActInner.size();
    int vnbr = V.rows();

    tripletes.clear();
    tripletes.reserve(ninner * 30);
    int ncondi = ninner * 2;
    if (!asymptotic)
    {
        ncondi = ninner * 4;
    }
    MTenergy = Eigen::VectorXd::Zero(ncondi); // mesh total energy values
    
    for (int i = 0; i < ninner; i++)
    {
        if (analizer.LocalActInner[i] == false) {
            std::cout<<"Singularity" <<std::endl;
            continue;
        }
        int vm = IVids[i];
        CGMesh::HalfedgeHandle inhd = analizer.heh0[i], outhd = analizer.heh1[i];
        int v1 = lsmesh.from_vertex_handle(inhd).idx();
        int v2 = lsmesh.to_vertex_handle(inhd).idx();
        int v3 = lsmesh.from_vertex_handle(outhd).idx();
        int v4 = lsmesh.to_vertex_handle(outhd).idx();
        double scale = mass_uniform.coeff(vm, vm);
        double dis0 = ((V.row(v1) - V.row(v2)) * func[vm] + (V.row(v2) - V.row(vm)) * func[v1] + (V.row(vm) - V.row(v1)) * func[v2]).norm();
        double dis1 = ((V.row(v3) - V.row(v4)) * func[vm] + (V.row(v4) - V.row(vm)) * func[v3] + (V.row(vm) - V.row(v3)) * func[v4]).norm();
        double t1 = analizer.t1s[i];
        double t2 = analizer.t2s[i];
        Eigen::Vector3d ver0 = V.row(v1) + (V.row(v2) - V.row(v1)) * t1;
        Eigen::Vector3d ver1 = V.row(vm);
        Eigen::Vector3d ver2 = V.row(v3) + (V.row(v4) - V.row(v3)) * t2;
        int lmx = vm;
        int lmy = vm + vnbr;
        int lmz = vm + vnbr * 2;
        int l1x = v1;
        int l1y = v1 + vnbr;
        int l1z = v1 + vnbr * 2;
        int l2x = v2;
        int l2y = v2 + vnbr;
        int l2z = v2 + vnbr * 2;
        int l3x = v3;
        int l3y = v3 + vnbr;
        int l3z = v3 + vnbr * 2;
        int l4x = v4;
        int l4y = v4 + vnbr;
        int l4z = v4 + vnbr * 2;
        Eigen::Vector3d norm = norm_v.row(vm);
        
        if (asymptotic)
        {
            int lfx = v1;
            int lfy = v1 + vnbr;
            int lfz = v1 + vnbr * 2;
            int ltx = v2;
            int lty = v2 + vnbr;
            int ltz = v2 + vnbr * 2;
            

            Eigen::Vector3d r12 = (func[v1] - func[v2]) * norm;
            Eigen::Vector3d rm1 = (func[vm] - func[v1]) * norm;
            Eigen::Vector3d r2m = (func[v2] - func[vm]) * norm;

            tripletes.push_back(Trip(i, lmx, r12[0] / dis0 * scale));
            tripletes.push_back(Trip(i, lmy, r12[1] / dis0 * scale));
            tripletes.push_back(Trip(i, lmz, r12[2] / dis0 * scale));

            tripletes.push_back(Trip(i, lfx, r2m[0] / dis0 * scale));
            tripletes.push_back(Trip(i, lfy, r2m[1] / dis0 * scale));
            tripletes.push_back(Trip(i, lfz, r2m[2] / dis0 * scale));

            tripletes.push_back(Trip(i, ltx, rm1[0] / dis0 * scale));
            tripletes.push_back(Trip(i, lty, rm1[1] / dis0 * scale));
            tripletes.push_back(Trip(i, ltz, rm1[2] / dis0 * scale));

            MTenergy[i] = (V.row(vm).dot(r12) + V.row(v1).dot(r2m) + V.row(v2).dot(rm1)) / dis0 * scale;

            // vf = v3, vt = v4
            lfx = v3;
            lfy = v3 + vnbr;
            lfz = v3 + vnbr * 2;
            ltx = v4;
            lty = v4 + vnbr;
            ltz = v4 + vnbr * 2;

            r12 = (func[v3] - func[v4]) * norm;
            rm1 = (func[vm] - func[v3]) * norm;
            r2m = (func[v4] - func[vm]) * norm;

            tripletes.push_back(Trip(i + ninner, lmx, r12[0] / dis1 * scale));
            tripletes.push_back(Trip(i + ninner, lmy, r12[1] / dis1 * scale));
            tripletes.push_back(Trip(i + ninner, lmz, r12[2] / dis1 * scale));

            tripletes.push_back(Trip(i + ninner, lfx, r2m[0] / dis1 * scale));
            tripletes.push_back(Trip(i + ninner, lfy, r2m[1] / dis1 * scale));
            tripletes.push_back(Trip(i + ninner, lfz, r2m[2] / dis1 * scale));

            tripletes.push_back(Trip(i + ninner, ltx, rm1[0] / dis1 * scale));
            tripletes.push_back(Trip(i + ninner, lty, rm1[1] / dis1 * scale));
            tripletes.push_back(Trip(i + ninner, ltz, rm1[2] / dis1 * scale));

            MTenergy[i + ninner] = (V.row(vm).dot(r12) + V.row(v3).dot(r2m) + V.row(v4).dot(rm1)) / dis1 * scale;
        }
        else{
            // std::cout<<"geodesic"<<std::endl;
            if(use_given_direction){
                // std::cout<<"dbgdbg"<<std::endl;
                norm=ray.normalized();
            }
            // std::cout<<"Dirc "<<norm.transpose()<<", ";
            double c1 = 1 - t1;
            double c2 = t1;
            double c3 = 1 - t2;
            double c4 = t2;
            double cm = -1;
            int lrx = i + aux_start_loc;
            int lry = i + aux_start_loc + ninner;
            int lrz = i + aux_start_loc + ninner * 2;
            // the two directions
            Eigen::Vector3d vec_l = c1 * V.row(v1) + c2 * V.row(v2) + cm * V.row(vm);
            Eigen::Vector3d vec_r = c3 * V.row(v3) + c4 * V.row(v4) + cm * V.row(vm);
            double dl = vec_l.norm();
            double dr = vec_r.norm();
            double s1 = vec_l.norm() * vec_r.norm() / scale;
            /////////////////
            // use the auxiliaries to improve the stability
            if (Compute_Auxiliaries_Mesh) // init the auxiliaries
            {
                Eigen::Vector3d real_r = vec_l.cross(vec_r);
                real_r = real_r.normalized();
                vars[lrx] = real_r[0];
                vars[lry] = real_r[1];
                vars[lrz] = real_r[2];
            }
            Eigen::Vector3d r = Eigen::Vector3d(vars[lrx], vars[lry], vars[lrz]);
            // r.dot(r) = 1
            tripletes.push_back(Trip(i, lrx, 2 * r(0) * scale));
            tripletes.push_back(Trip(i, lry, 2 * r(1) * scale));
            tripletes.push_back(Trip(i, lrz, 2 * r(2) * scale));
            MTenergy[i] = (r.dot(r) - 1) * scale;

            // r.dot(dl) = 0
            tripletes.push_back(Trip(i + ninner, lrx, vec_l(0) / dl * scale));
            tripletes.push_back(Trip(i + ninner, lry, vec_l(1) / dl * scale));
            tripletes.push_back(Trip(i + ninner, lrz, vec_l(2) / dl * scale));

            tripletes.push_back(Trip(i + ninner, l1x, c1 * r(0) / dl * scale));
            tripletes.push_back(Trip(i + ninner, l1y, c1 * r(1) / dl * scale));
            tripletes.push_back(Trip(i + ninner, l1z, c1 * r(2) / dl * scale));

            tripletes.push_back(Trip(i + ninner, l2x, c2 * r(0) / dl * scale));
            tripletes.push_back(Trip(i + ninner, l2y, c2 * r(1) / dl * scale));
            tripletes.push_back(Trip(i + ninner, l2z, c2 * r(2) / dl * scale));

            tripletes.push_back(Trip(i + ninner, lmx, cm * r(0) / dl * scale));
            tripletes.push_back(Trip(i + ninner, lmy, cm * r(1) / dl * scale));
            tripletes.push_back(Trip(i + ninner, lmz, cm * r(2) / dl * scale));

            MTenergy[i + ninner] = r.dot(vec_l) / dl * scale;

            // r.dot(dr) = 0
            tripletes.push_back(Trip(i + ninner * 2, lrx, vec_r(0) / dr * scale));
            tripletes.push_back(Trip(i + ninner * 2, lry, vec_r(1) / dr * scale));
            tripletes.push_back(Trip(i + ninner * 2, lrz, vec_r(2) / dr * scale));

            tripletes.push_back(Trip(i + ninner * 2, l3x, c3 * r(0) / dr * scale));
            tripletes.push_back(Trip(i + ninner * 2, l3y, c3 * r(1) / dr * scale));
            tripletes.push_back(Trip(i + ninner * 2, l3z, c3 * r(2) / dr * scale));

            tripletes.push_back(Trip(i + ninner * 2, l4x, c4 * r(0) / dr * scale));
            tripletes.push_back(Trip(i + ninner * 2, l4y, c4 * r(1) / dr * scale));
            tripletes.push_back(Trip(i + ninner * 2, l4z, c4 * r(2) / dr * scale));

            tripletes.push_back(Trip(i + ninner * 2, lmx, cm * r(0) / dr * scale));
            tripletes.push_back(Trip(i + ninner * 2, lmy, cm * r(1) / dr * scale));
            tripletes.push_back(Trip(i + ninner * 2, lmz, cm * r(2) / dr * scale));

            MTenergy[i + ninner * 2] = r.dot(vec_r) / dr * scale;

            // r.dot(n) = 0
            tripletes.push_back(Trip(i + ninner * 3, lrx, norm(0) * scale));
            tripletes.push_back(Trip(i + ninner * 3, lry, norm(1) * scale));
            tripletes.push_back(Trip(i + ninner * 3, lrz, norm(2) * scale));

            MTenergy[i + ninner * 3] = r.dot(norm) * scale;

            ///////////////////
        }
    }
    PGE = MTenergy;
    // if(asymptotic){
    //     PGE = MTenergy;
    // }
}

void lsTools::calculate_mesh_opt_shading_condition_values(const Eigen::VectorXd &func, const Eigen::Vector3d &ray,
                                                const LSAnalizer &analizer, std::vector<Trip> &tripletes, Eigen::VectorXd &MTenergy)
{

    int ninner = analizer.LocalActInner.size();
    int vnbr = V.rows();

    tripletes.clear();
    tripletes.reserve(ninner * 18);
    MTenergy = Eigen::VectorXd::Zero(ninner * 2); // mesh total energy values
    
    for (int i = 0; i < ninner; i++)
    {
        if (analizer.LocalActInner[i] == false) {
            std::cout<<"Singularity" <<std::endl;
            continue;
        }
        int vm = IVids[i];
        CGMesh::HalfedgeHandle inhd = analizer.heh0[i], outhd = analizer.heh1[i];
        int v1 = lsmesh.from_vertex_handle(inhd).idx();
        int v2 = lsmesh.to_vertex_handle(inhd).idx();
        int v3 = lsmesh.from_vertex_handle(outhd).idx();
        int v4 = lsmesh.to_vertex_handle(outhd).idx();
        double dis0 = ((V.row(v1) - V.row(v2)) * func[vm] + (V.row(v2) - V.row(vm)) * func[v1] + (V.row(vm) - V.row(v1)) * func[v2]).norm();
        double dis1 = ((V.row(v3) - V.row(v4)) * func[vm] + (V.row(v4) - V.row(vm)) * func[v3] + (V.row(vm) - V.row(v3)) * func[v4]).norm();
        double t1 = analizer.t1s[i];
        double t2 = analizer.t2s[i];
        Eigen::Vector3d ver0 = V.row(v1) + (V.row(v2) - V.row(v1)) * t1;
        Eigen::Vector3d ver1 = V.row(vm);
        Eigen::Vector3d ver2 = V.row(v3) + (V.row(v4) - V.row(v3)) * t2;
        int lmx = vm;
        int lmy = vm + vnbr;
        int lmz = vm + vnbr * 2;
        int l1x = v1;
        int l1y = v1 + vnbr;
        int l1z = v1 + vnbr * 2;
        int l2x = v2;
        int l2y = v2 + vnbr;
        int l2z = v2 + vnbr * 2;
        int l3x = v3;
        int l3y = v3 + vnbr;
        int l3z = v3 + vnbr * 2;
        int l4x = v4;
        int l4y = v4 + vnbr;
        int l4z = v4 + vnbr * 2;

        Eigen::Vector3d norm = ray.normalized();

        double c1 = t1 - 1;
        double c2 = -t1;
        double c3 = 1 - t2;
        double c4 = t2;
        // double cm = -1;

        Eigen::Vector3d tangent = c1 * V.row(v1) + c2 * V.row(v2) + c3 * V.row(v3) + c4 * V.row(v4);
        double scale = tangent.norm();
        // to v1
        tripletes.push_back(Trip(i, l1x, c1 * norm(0) / scale));
        tripletes.push_back(Trip(i, l1y, c1 * norm(1) / scale));
        tripletes.push_back(Trip(i, l1z, c1 * norm(2) / scale));
        // to v2
        tripletes.push_back(Trip(i, l2x, c2 * norm(0) / scale));
        tripletes.push_back(Trip(i, l2y, c2 * norm(1) / scale));
        tripletes.push_back(Trip(i, l2z, c2 * norm(2) / scale));
        // to v1
        tripletes.push_back(Trip(i, l3x, c3 * norm(0) / scale));
        tripletes.push_back(Trip(i, l3y, c3 * norm(1) / scale));
        tripletes.push_back(Trip(i, l3z, c3 * norm(2) / scale));
        // to v1
        tripletes.push_back(Trip(i, l4x, c4 * norm(0) / scale));
        tripletes.push_back(Trip(i, l4y, c4 * norm(1) / scale));
        tripletes.push_back(Trip(i, l4z, c4 * norm(2) / scale));

        MTenergy[i] = norm.dot(tangent) / scale;
    }
}
void lsTools::assemble_solver_shading_mesh_opt(const Eigen::VectorXd &func, const Eigen::Vector3d &ray,
                                           const LSAnalizer &analizer,
                                           spMat &JTJ, Eigen::VectorXd &B, Eigen::VectorXd &MTEnergy)
{
    std::vector<Trip> tripletes;
    int vsize = V.rows();
    calculate_mesh_opt_shading_condition_values(func, ray,analizer, tripletes, MTEnergy);
    int nvars = vsize * 3;
	int ncondi = MTEnergy.size();
    spMat J;
    J.resize(ncondi, nvars);
    J.setFromTriplets(tripletes.begin(), tripletes.end());
    JTJ = J.transpose() * J;
    B = -J.transpose() * MTEnergy;
}
spMat Jacobian_transpose_mesh_opt_on_ver(const std::array<spMat, 3>& JC,
	const Eigen::Vector3d& norm, const spMat& SMfvalues) {

	spMat result = (JC[0] * norm[0] + JC[1] * norm[1] + JC[2] * norm[2]) * SMfvalues;

	return result;
}
void lsTools::assemble_solver_mesh_opt_part( Eigen::VectorXd& vars,
    const LSAnalizer &analizer,
    const std::vector<double>& angle_degrees, const int aux_start_loc, spMat& JTJ, Eigen::VectorXd& B, Eigen::VectorXd& MTEnergy){
	std::vector<Trip> tripletes;
    int vsize = V.rows();
    int ninner = analizer.LocalActInner.size();
	calculate_mesh_opt_expanded_function_values(vars, analizer, angle_degrees, aux_start_loc, tripletes, MTEnergy);
        
    int nvars = vars.size();
    int ncondi = MTEnergy.size();
    spMat J;
    J.resize(ncondi, nvars);
    J.setFromTriplets(tripletes.begin(), tripletes.end());
    JTJ = J.transpose() * J;
    B = -J.transpose() * MTEnergy;
}
void lsTools::assemble_solver_mesh_extreme(Eigen::VectorXd& vars, const int aux_start_loc, const Eigen::VectorXd &func, const bool asymptotic, const bool use_given_direction, const Eigen::Vector3d &ray,
                                           const LSAnalizer &analizer,
                                           spMat &JTJ, Eigen::VectorXd &B, Eigen::VectorXd &MTEnergy)
{
    std::vector<Trip> tripletes;
    int vsize = V.rows();
    calculate_mesh_opt_extreme_values(vars, aux_start_loc, func, asymptotic, use_given_direction, ray, analizer, tripletes, MTEnergy);
    int nvars = vars.size();
    // int ninner = analizer.LocalActInner.size();
    
    int ncondi = MTEnergy.size();
    spMat J;
    J.resize(ncondi, nvars);
    J.setFromTriplets(tripletes.begin(), tripletes.end());
    JTJ = J.transpose() * J;
    B = -J.transpose() * MTEnergy;
}
void lsTools::assemble_solver_mesh_smoothing(const Eigen::VectorXd &vars, spMat& H, Eigen::VectorXd &B){
    spMat JTJ=igl::repdiag(QcH,3);// the matrix size nx3 x nx3
    Eigen::VectorXd mJTF=-JTJ*vars;
    H=JTJ;
    B=mJTF;
}
void lsTools::update_mesh_properties(){
    // first update all the vertices in openmesh
    int nv = lsmesh.n_vertices();
	for (CGMesh::VertexIter v_it = lsmesh.vertices_begin(); v_it != lsmesh.vertices_end(); ++v_it)
	{
        int vid=v_it.handle().idx();
		lsmesh.point(*v_it)=CGMesh::Point(V(vid,0),V(vid,1),V(vid,2));
	}
    igl::cotmatrix(V, F, Dlps);
    igl::curved_hessian_energy(V, F, QcH);

    // 1
    get_mesh_normals_per_face();
    // 2
    get_mesh_angles();// only useful when computing gradient ver to F. and the vertex normals 
    // // 3 if use this, please don't use get_I_and_II_locally()
    get_mesh_normals_per_ver();
    // 4
    get_face_rotation_matices();// not useful, just for debug

    // 5
    get_rotated_edges_for_each_face();// not useful
    // 6
    get_function_gradient_vertex();

    // 7
    get_function_hessian_vertex();// not useful

    // get_I_and_II_locally(); // useful when we consider curvatures
    get_vertex_rotation_matices();// not useful
    get_all_the_edge_normals();// necessary. In case some one wants to trace again.
}

void solve_mean_value_laplacian_mat(CGMesh& lsmesh,const std::vector<int>& IVids, spMat& mat) {
    int nbr = lsmesh.n_vertices();
    std::vector<Trip> triplets;
    triplets.reserve(nbr * 7);
    spMat mat2;

    mat2.resize(IVids.size() * 3, nbr * 3);
    for (int i = 0; i < IVids.size(); i++) {
        int vid = IVids[i];
        assert(vid >= 0);
        CGMesh::VertexHandle vh = lsmesh.vertex_handle(vid);
        double valence = lsmesh.valence(vh);
        for (CGMesh::VertexVertexIter vv_it = lsmesh.vv_begin(vh); vv_it != lsmesh.vv_end(vh); ++vv_it)
        {
            int id = vv_it.handle().idx();
            int loc;
            location_in_sparse_matrix(nbr, id, 0, loc);
            triplets.push_back(Trip(i, loc, -1 / valence));//x

            location_in_sparse_matrix(nbr, id, 1, loc);
            triplets.push_back(Trip(i + IVids.size(), loc, -1 / valence));//y

            location_in_sparse_matrix(nbr, id, 2, loc);
            triplets.push_back(Trip(i + IVids.size() * 2, loc, -1 / valence));//z
		}
		int loc;
		location_in_sparse_matrix(nbr, vid, 0, loc);
		triplets.push_back(Trip(i, loc, 1));
		location_in_sparse_matrix(nbr, vid, 1, loc);
		triplets.push_back(Trip(i + IVids.size(), loc, 1));
		location_in_sparse_matrix(nbr, vid, 2, loc);
		triplets.push_back(Trip(i + IVids.size() * 2, loc, 1));
    }
    mat2.setFromTriplets(triplets.begin(), triplets.end());
    mat = mat2;

}
void lsTools::assemble_solver_mean_value_laplacian(const Eigen::VectorXd& vars, spMat& H, Eigen::VectorXd& B) {
    spMat JTJ = MVLap.transpose()*MVLap;
    Eigen::VectorXd mJTF = -JTJ * vars;
    H = JTJ;
    B = mJTF;
}

void lsTools::solve_edge_length_matrix(const Eigen::MatrixXd& V, const Eigen::MatrixXi& E, spMat& mat) {
    int enbr = E.rows();
    int nver = V.rows();
    std::vector<Trip> tripletes;
    tripletes.reserve(enbr*6);
    ElStored.resize(enbr * 3);
    for (int i = 0; i < enbr; i++) {
        int vid0 = E(i, 0);
		int vid1 = E(i, 1);
		tripletes.push_back(Trip(3 * i, vid0, 1));
		tripletes.push_back(Trip(3 * i, vid1, -1));
		tripletes.push_back(Trip(3 * i + 1, nver + vid0, 1));
		tripletes.push_back(Trip(3 * i + 1, nver + vid1, -1));
		tripletes.push_back(Trip(3 * i + 2, nver * 2 + vid0, 1));
		tripletes.push_back(Trip(3 * i + 2, nver * 2 + vid1, -1));
        ElStored[3 * i] = V(vid0, 0) - V(vid1, 0);
        ElStored[3 * i+1] = V(vid0, 1) - V(vid1, 1);
        ElStored[3 * i+2] = V(vid0, 2) - V(vid1, 2);
    }
    mat.resize(enbr * 3, nver * 3);
    mat.setFromTriplets(tripletes.begin(), tripletes.end());
}
void lsTools::assemble_solver_mesh_edge_length_part(const Eigen::VectorXd vars, spMat& H, Eigen::VectorXd& B, 
Eigen::VectorXd& ElEnergy) {
    H = Elmat.transpose() * Elmat;
    int enbr = E.rows();
    int nver = V.rows();
    B.resize(enbr * 3);
    for (int i = 0; i < enbr; i++) {
        int vid0 = E(i, 0);
        int vid1 = E(i, 1);
        B[i * 3] = vars[vid0] - vars[vid1] - ElStored[i * 3];
		B[i * 3 + 1] = vars[nver + vid0] - vars[nver + vid1] - ElStored[i * 3 + 1];
        B[i * 3 + 2] = vars[nver*2 + vid0] - vars[nver*2 + vid1] - ElStored[i * 3 + 2];
    }
    ElEnergy = B;
    B = -Elmat.transpose() * B;
}

// the mat must be a symmetric matrix
void spmat_on_corner(const spMat& mat, spMat& target, int nmat, int ntarget) {
    spMat tmp;
    tmp.resize(nmat, ntarget);
    tmp.leftCols(nmat) = mat;
    target.leftCols(nmat) = tmp.transpose();
}
// they must be symmetric matrices
spMat sum_uneven_spMats(const spMat& mat_small, const spMat& mat_large) {
    int snbr = mat_small.rows();
    int lnbr = mat_large.rows();
    if (snbr == lnbr) {
        return mat_small + mat_large;
    }
    spMat tmp;
    tmp.resize(snbr, lnbr);
    tmp.leftCols(snbr) = mat_small;
    spMat result;
    result.resize(lnbr, lnbr);
    result.leftCols(snbr) = tmp.transpose();
	return result + mat_large;
}
spMat three_spmat_in_diag(const spMat& mat0, const spMat& mat1, const spMat& mat2, const int ntarget){
    int ns=mat0.rows();
    spMat tmp0(ns, ntarget), tmp1(ns, ntarget), tmp2(ns, ntarget);
    tmp0.leftCols(ns)=mat0;
    tmp1.middleCols(ns,ns)=mat1;
    tmp2.middleCols(ns*2, ns)=mat2;
    tmp0=tmp0.transpose();
    tmp1=tmp1.transpose();
    tmp2=tmp2.transpose();
    spMat result(ntarget, ntarget);
    result.leftCols(ns)=tmp0;
    result.middleCols(ns, ns)=tmp1;
    result.middleCols(ns*2, ns)= tmp2;
    return result;
}
Eigen::VectorXd three_vec_in_row(const Eigen::VectorXd& ve0, const Eigen::VectorXd& ve1, const Eigen::VectorXd& ve2, const int ntarget){
    int nsize= ve0.size();
    Eigen::VectorXd result=Eigen::VectorXd::Zero(ntarget);
    result.segment(0,nsize)=ve0;
    result.segment(nsize,nsize)=ve1;
    result.segment(nsize*2, nsize)=ve2;
    return result;
}

Eigen::VectorXd sum_uneven_vectors(const Eigen::VectorXd& vsmall, const Eigen::VectorXd& vlarge) {
    int snbr = vsmall.size();
    int lnbr = vlarge.size();
    if (snbr == lnbr) {
        return vsmall + vlarge;
    }
    Eigen::VectorXd tmp=Eigen::VectorXd::Zero(lnbr);
    tmp.topRows(snbr) = vsmall;
    return tmp + vlarge;

}
void lsTools::Run_Mesh_Opt(){
    
    Eigen::VectorXd func = fvalues;

    std::vector<double> angle_degrees(1);
    angle_degrees[0] = pseudo_geodesic_target_angle_degree;
    Eigen::MatrixXd GradValueF, GradValueV;
    
    get_gradient_hessian_values(func, GradValueV, GradValueF);
    bool first_compute = true; // if we need initialize auxiliary vars
    if(!Last_Opt_Mesh){
        std::cout<<"init mesh opt"<<std::endl;
        Compute_Auxiliaries_Mesh=true;
        first_compute = true;// if last time opt levelset, we re-compute the auxiliary vars
    }
    analysis_pseudo_geodesic_on_vertices(func, anas[0]);
    int ninner = anas[0].LocalActInner.size();
    int vnbr=V.rows();
	int final_size = ninner * 7 + vnbr * 3;// Change this when using more auxilary vars
   
    Eigen::VectorXd vars;// the variables feed to the energies without auxiliary variables
    vars.resize(vnbr * 3);
    if (Glob_Vars.size() == 0) {
        std::cout << "Initializing Global Variable For Mesh Opt ... " << std::endl;
        Glob_Vars=Eigen::VectorXd::Zero(final_size);// We change the size if opt more than 1 level set
        Glob_Vars.segment(0, vnbr) = V.col(0);
        Glob_Vars.segment(vnbr, vnbr) = V.col(1);
		Glob_Vars.segment(vnbr * 2, vnbr) = V.col(2);
        std::cout << "Initialized Global Variable" << std::endl;
    }
    vars.topRows(vnbr) = V.col(0);
    vars.middleRows(vnbr, vnbr) = V.col(1);
    vars.bottomRows(vnbr) = V.col(2);
    
    spMat H;
    Eigen::VectorXd B;
    H.resize(vnbr * 3, vnbr * 3);
    B = Eigen::VectorXd::Zero(vnbr * 3);
    
    
    spMat Hsmooth;
    Eigen::VectorXd Bsmooth;
    assemble_solver_mesh_smoothing(vars, Hsmooth, Bsmooth);
    //assemble_solver_mean_value_laplacian(vars, Hsmooth, Bsmooth);
    H += weight_Mesh_smoothness * Hsmooth;
    B += weight_Mesh_smoothness * Bsmooth;
    
    spMat Hel;
    Eigen::VectorXd Bel;
    Eigen::VectorXd ElEnergy;
    assemble_solver_mesh_edge_length_part(vars, Hel, Bel, ElEnergy);
    H += weight_Mesh_edgelength * Hel;
    B += weight_Mesh_edgelength * Bel;
    
	spMat Hpg;
    Eigen::VectorXd Bpg;
    int aux_start_loc = vnbr * 3;// For mesh opt the auxiliary vars start from vnbr*3
    Eigen::VectorXd MTEnergy;
    if (!enable_extreme_cases) {
        assemble_solver_mesh_opt_part(Glob_Vars,
            anas[0], angle_degrees, aux_start_loc, Hpg, Bpg, MTEnergy);
    }
    else {
        bool asymptotic = true;
        if (angle_degrees[0] != 0)
        {
            asymptotic = false;
        }
        Eigen::Vector3d any_ray;
        // in mesh opt we do not optimize for shading
        assemble_solver_mesh_extreme(Glob_Vars, aux_start_loc, func, asymptotic, false, any_ray, anas[0], Hpg, Bpg, MTEnergy);

    }
    Compute_Auxiliaries_Mesh = false;
    spMat Htotal(final_size, final_size);
    Eigen::VectorXd Btotal=Eigen::VectorXd::Zero(final_size);
    // first make the naive energies the correct size
    Htotal = sum_uneven_spMats(H, Htotal);
    Btotal = sum_uneven_vectors(B, Btotal);

    // add the PG energy
    Htotal = sum_uneven_spMats(weight_Mesh_pesudo_geodesic * Hpg, Htotal);
    Btotal = sum_uneven_vectors(weight_Mesh_pesudo_geodesic * Bpg, Btotal);    
    Htotal += spMat(weight_mass * 1e-6 * Eigen::VectorXd::Ones(final_size).asDiagonal());

    if(vector_contains_NAN(Btotal)){
        std::cout<<"energy contains NAN"<<std::endl;
    }

    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver(Htotal);

    if (solver.info() != Eigen::Success)
    {
        std::cout << "solver fail" << std::endl;
        return;
    }

    Eigen::VectorXd dx = solver.solve(Btotal).eval(); 

    dx *= 0.75;
    double mesh_opt_step_length = dx.norm();
    // double inf_norm=dx.cwiseAbs().maxCoeff();
    if (mesh_opt_step_length > Mesh_opt_max_step_length)
    {
        dx *= Mesh_opt_max_step_length / mesh_opt_step_length;
    }
    vars+=dx.topRows(vnbr*3);
    Glob_Vars += dx;
    V.col(0)=vars.topRows(vnbr);
    V.col(1)=vars.middleRows(vnbr,vnbr);
    V.col(2)=vars.bottomRows(vnbr);
    if (vars != Glob_Vars.topRows(vnbr * 3)) {
        std::cout << "vars diff from glob vars" << std::endl;
    }
    double energy_smooth=(V.transpose()*QcH*V).norm();
    /*double energy_mvl = (MVLap * vars).norm();*/
    std::cout<<"Mesh Opt: smooth, "<< energy_smooth <<", ";
    if (!enable_extreme_cases) {
        double energy_ls = MTEnergy.norm();
        double max_energy_ls = MTEnergy.lpNorm<Eigen::Infinity>();
        std::cout << "pg, " << energy_ls << ", MaxEnergy, " << max_energy_ls << ", ";
        double max_ls_angle_energy = MTEnergy.bottomRows(ninner).norm();
        std::cout << "total angle energy, " << max_ls_angle_energy << ", ";
    }
    else {
        double planar_energy = MTEnergy.norm();
        std::cout << "pg, " << planar_energy << ", max "<<MTEnergy.lpNorm<Eigen::Infinity>()<<", ";
    }
    
    double energy_el = ElEnergy.norm();
    std::cout << "el, " << energy_el << ", maxloc, ";
    step_length=dx.norm();
    std::cout<<"step "<<step_length<<std::endl;
    update_mesh_properties();
    Last_Opt_Mesh=true;

}
void levelset_unit_scale(Eigen::VectorXd& func, Eigen::MatrixXd &GradValueF, const double length){
    int fnbr = GradValueF.rows();
    Eigen::VectorXd gn = GradValueF.rowwise().norm();// the lengths of gradients
    double avg = gn.dot(Eigen::VectorXd::Ones(fnbr)) / fnbr; // the average gradient length
    double ratio = length / avg;
    func *= ratio;
    GradValueF *= ratio;
    std::cout<<"Target strip width, "<<length<<", Current average length, "<<avg<<std::endl;
    std::cout<<"the scale ratio, "<<ratio<<std::endl;
    gn = GradValueF.rowwise().norm();
    avg = gn.dot(Eigen::VectorXd::Ones(fnbr)) / fnbr;
    std::cout<<"gradient avg, "<<avg;
}
void lsTools::Run_AAG_Mesh_Opt(Eigen::VectorXd& func0, Eigen::VectorXd& func1, Eigen::VectorXd& func2){

    // levelset unit scale
    bool first_compute = true; // if we need initialize auxiliary vars
    if (!Last_Opt_Mesh)
    {
        Eigen::MatrixXd GradValueF[3], GradValueV[3];
        get_gradient_hessian_values(func0, GradValueV[0], GradValueF[0]);
        get_gradient_hessian_values(func1, GradValueV[1], GradValueF[1]);
        levelset_unit_scale(func0, GradValueF[0], 1);
        levelset_unit_scale(func1, GradValueF[1], 1);
        func2 = -func0 - func1;
        std::cout << "** Mesh Opt initialization done" << std::endl;
        Compute_Auxiliaries_Mesh=true;
        analysis_pseudo_geodesic_on_vertices(func0, anas[0]);
        analysis_pseudo_geodesic_on_vertices(func1, anas[1]);
        analysis_pseudo_geodesic_on_vertices(func2, anas[2]);
        first_compute = true;// if last time opt levelset, we re-compute the auxiliary vars
    }
    int ninner = anas[0].LocalActInner.size();
    int vnbr=V.rows();
    int final_size = vnbr * 3 + ninner * 3; // Change this when using more auxilary vars. Only G use auxiliaries

    spMat H;
    Eigen::VectorXd B;
    H.resize(vnbr * 3, vnbr * 3);
    B = Eigen::VectorXd::Zero(vnbr * 3);
    if (Glob_Vars.size() == 0) {
        std::cout << "Initializing Global Variable For Mesh Opt ... " << std::endl;
        Glob_Vars=Eigen::VectorXd::Zero(final_size);// We change the size if opt more than 1 level set
        Glob_Vars.segment(0, vnbr) = V.col(0);
        Glob_Vars.segment(vnbr, vnbr) = V.col(1);
		Glob_Vars.segment(vnbr * 2, vnbr) = V.col(2);
    }
    Eigen::VectorXd vars;// the variables feed to the energies without auxiliary variables
    vars.resize(vnbr * 3);
    vars.topRows(vnbr) = V.col(0);
    vars.middleRows(vnbr, vnbr) = V.col(1);
    vars.bottomRows(vnbr) = V.col(2);
    
    spMat Hsmooth;
    Eigen::VectorXd Bsmooth;
    assemble_solver_mesh_smoothing(vars, Hsmooth, Bsmooth);// as smooth as possible
    //assemble_solver_mean_value_laplacian(vars, Hsmooth, Bsmooth);
    H += weight_Mesh_smoothness * Hsmooth;
    B += weight_Mesh_smoothness * Bsmooth;
    
    spMat Hel;
    Eigen::VectorXd Bel;
    Eigen::VectorXd ElEnergy;
    assemble_solver_mesh_edge_length_part(vars, Hel, Bel, ElEnergy);// as rigid as possible
    H += weight_Mesh_edgelength * Hel;
    B += weight_Mesh_edgelength * Bel;

	spMat Hpg[3];
    Eigen::VectorXd Bpg[3];
    int aux_start_loc = vnbr * 3;// For mesh opt the auxiliary vars start from vnbr*3
    Eigen::VectorXd MTEnergy[3];
    Eigen::Vector3d any_ray;
    // A
    assemble_solver_mesh_extreme(Glob_Vars, aux_start_loc, func0, true, false, any_ray, anas[0], Hpg[0], Bpg[0], MTEnergy[0]);
    // A
    assemble_solver_mesh_extreme(Glob_Vars, aux_start_loc, func1, true, false, any_ray, anas[1], Hpg[1], Bpg[1], MTEnergy[1]);
    // G
    assemble_solver_mesh_extreme(Glob_Vars, aux_start_loc, func2, false, false, any_ray, anas[2], Hpg[2], Bpg[2], MTEnergy[2]);

    Compute_Auxiliaries_Mesh = false;
    spMat Htotal(final_size, final_size);
    Eigen::VectorXd Btotal=Eigen::VectorXd::Zero(final_size);
    // first make the naive energies the correct size
    Htotal = sum_uneven_spMats(H, Htotal);
    Btotal = sum_uneven_vectors(B, Btotal);

    // add the PG energy
    Htotal = sum_uneven_spMats(weight_Mesh_pesudo_geodesic * (Hpg[0] + Hpg[1]), Htotal);
    Htotal += weight_Mesh_pesudo_geodesic * weight_geodesic * Hpg[2];
    Btotal = sum_uneven_vectors(weight_Mesh_pesudo_geodesic * (Bpg[0] + Bpg[1]), Btotal);
    Btotal += weight_Mesh_pesudo_geodesic * weight_geodesic * Bpg[2];

    Htotal += spMat(weight_mass * 1e-6 * Eigen::VectorXd::Ones(final_size).asDiagonal());

    if(vector_contains_NAN(Btotal)){
        std::cout<<"energy contains NAN"<<std::endl;
    }

    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver(Htotal);

    if (solver.info() != Eigen::Success)
    {
        std::cout << "solver fail" << std::endl;
        return;
    }

    Eigen::VectorXd dx = solver.solve(Btotal).eval(); 

    dx *= 0.75;
    double mesh_opt_step_length = dx.norm();
    // double inf_norm=dx.cwiseAbs().maxCoeff();
    if (mesh_opt_step_length > Mesh_opt_max_step_length)
    {
        dx *= Mesh_opt_max_step_length / mesh_opt_step_length;
    }
    vars += dx.topRows(vnbr * 3);
    Glob_Vars += dx;
    V.col(0)=vars.topRows(vnbr);
    V.col(1)=vars.middleRows(vnbr,vnbr);
    V.col(2)=vars.bottomRows(vnbr);
    if (vars != Glob_Vars.topRows(vnbr * 3)) {
        std::cout << "vars diff from glob vars" << std::endl;
    }
    double energy_smooth=(V.transpose()*QcH*V).norm();
    std::cout<<"Mesh Opt: smooth, "<< energy_smooth <<", ";
    /*double energy_mvl = (MVLap * vars).norm();*/
    double energy_ls[3];
    energy_ls[0] = MTEnergy[0].norm();
    energy_ls[1] = MTEnergy[1].norm();
    energy_ls[2] = MTEnergy[2].norm();
    std::cout << "pg, " << energy_ls[0]<<", "<< energy_ls[1]<<", "<< energy_ls[2]<<", pgmax, "<<MTEnergy[0].lpNorm<Eigen::Infinity>()
		<<", "<<MTEnergy[1].lpNorm<Eigen::Infinity>()<<", "<<MTEnergy[2].lpNorm<Eigen::Infinity>()<<", ";

    double energy_el = ElEnergy.norm();
    std::cout << "el, " << energy_el << ", ";
    step_length=dx.norm();
    std::cout<<"step "<<step_length<<std::endl;
    update_mesh_properties();
    Last_Opt_Mesh=true;
}

// in the order of agg
void lsTools::Run_AGG_Mesh_Opt(Eigen::VectorXd& func0, Eigen::VectorXd& func1, Eigen::VectorXd& func2){

    // levelset unit scale
    bool first_compute = true; // if we need initialize auxiliary vars
    if (!Last_Opt_Mesh)
    {
        Eigen::MatrixXd GradValueF[3], GradValueV[3];
        get_gradient_hessian_values(func0, GradValueV[0], GradValueF[0]);
        get_gradient_hessian_values(func1, GradValueV[1], GradValueF[1]);
        levelset_unit_scale(func0, GradValueF[0], 1);
        levelset_unit_scale(func1, GradValueF[1], 1);
        func2 = -func0 - func1;
        std::cout << "** Mesh Opt initialization done" << std::endl;
        Compute_Auxiliaries_Mesh=true;
        analysis_pseudo_geodesic_on_vertices(func0, anas[0]);
        analysis_pseudo_geodesic_on_vertices(func1, anas[1]);
        analysis_pseudo_geodesic_on_vertices(func2, anas[2]);
        first_compute = true;// if last time opt levelset, we re-compute the auxiliary vars
    }
    int ninner = anas[0].LocalActInner.size();
    int vnbr = V.rows();
    int final_size = vnbr * 3 + ninner * 6; // Change this when using more auxilary vars. Only G use auxiliaries

    spMat H;
    Eigen::VectorXd B;
    H.resize(vnbr * 3, vnbr * 3);
    B = Eigen::VectorXd::Zero(vnbr * 3);
    if (Glob_Vars.size() == 0) {
        std::cout << "Initializing Global Variable For Mesh Opt ... " << std::endl;
        Glob_Vars=Eigen::VectorXd::Zero(final_size);// We change the size if opt more than 1 level set
        Glob_Vars.segment(0, vnbr) = V.col(0);
        Glob_Vars.segment(vnbr, vnbr) = V.col(1);
		Glob_Vars.segment(vnbr * 2, vnbr) = V.col(2);
    }
    Eigen::VectorXd vars;// the variables feed to the energies without auxiliary variables
    vars.resize(vnbr * 3);
    vars.topRows(vnbr) = V.col(0);
    vars.middleRows(vnbr, vnbr) = V.col(1);
    vars.bottomRows(vnbr) = V.col(2);
    
    spMat Hsmooth;
    Eigen::VectorXd Bsmooth;
    assemble_solver_mesh_smoothing(vars, Hsmooth, Bsmooth);// as smooth as possible
    //assemble_solver_mean_value_laplacian(vars, Hsmooth, Bsmooth);
    H += weight_Mesh_smoothness * Hsmooth;
    B += weight_Mesh_smoothness * Bsmooth;
    
    spMat Hel;
    Eigen::VectorXd Bel;
    Eigen::VectorXd ElEnergy;
    assemble_solver_mesh_edge_length_part(vars, Hel, Bel, ElEnergy);// as rigid as possible
    H += weight_Mesh_edgelength * Hel;
    B += weight_Mesh_edgelength * Bel;

	spMat Hpg[3];
    Eigen::VectorXd Bpg[3];
    
    Eigen::VectorXd MTEnergy[3];
    Eigen::Vector3d any_ray;
    // A
    int aux_start_loc = vnbr * 3;// For mesh opt the auxiliary vars start from vnbr*3
    assemble_solver_mesh_extreme(Glob_Vars, aux_start_loc, func0, true, false, any_ray, anas[0], Hpg[0], Bpg[0], MTEnergy[0]);
    // G
    assemble_solver_mesh_extreme(Glob_Vars, aux_start_loc, func1, false, false, any_ray, anas[1], Hpg[1], Bpg[1], MTEnergy[1]);
    // G
    aux_start_loc = vnbr * 3 + ninner * 3;
    assemble_solver_mesh_extreme(Glob_Vars, aux_start_loc, func2, false, false, any_ray, anas[2], Hpg[2], Bpg[2], MTEnergy[2]);

    Compute_Auxiliaries_Mesh = false;
    spMat Htotal(final_size, final_size);
    Eigen::VectorXd Btotal=Eigen::VectorXd::Zero(final_size);
    // first make the naive energies the correct size
    Htotal = sum_uneven_spMats(H, Htotal);
    Btotal = sum_uneven_vectors(B, Btotal);

    // add the PG energy
    Htotal = sum_uneven_spMats(weight_Mesh_pesudo_geodesic * (Hpg[0] + Hpg[1]), Htotal);
    Htotal += weight_Mesh_pesudo_geodesic * weight_geodesic * Hpg[2];
    Btotal = sum_uneven_vectors(weight_Mesh_pesudo_geodesic * (Bpg[0] + Bpg[1]), Btotal);
    Btotal += weight_Mesh_pesudo_geodesic * weight_geodesic * Bpg[2];

    Htotal += spMat(weight_mass * 1e-6 * Eigen::VectorXd::Ones(final_size).asDiagonal());

    if(vector_contains_NAN(Btotal)){
        std::cout<<"energy contains NAN"<<std::endl;
    }

    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver(Htotal);

    if (solver.info() != Eigen::Success)
    {
        std::cout << "solver fail" << std::endl;
        return;
    }

    Eigen::VectorXd dx = solver.solve(Btotal).eval(); 

    dx *= 0.75;
    double mesh_opt_step_length = dx.norm();
    // double inf_norm=dx.cwiseAbs().maxCoeff();
    if (mesh_opt_step_length > Mesh_opt_max_step_length)
    {
        dx *= Mesh_opt_max_step_length / mesh_opt_step_length;
    }
    vars += dx.topRows(vnbr * 3);
    Glob_Vars += dx;
    V.col(0)=vars.topRows(vnbr);
    V.col(1)=vars.middleRows(vnbr,vnbr);
    V.col(2)=vars.bottomRows(vnbr);
    if (vars != Glob_Vars.topRows(vnbr * 3)) {
        std::cout << "vars diff from glob vars" << std::endl;
    }
    double energy_smooth=(V.transpose()*QcH*V).norm();
    std::cout<<"Mesh Opt: smooth, "<< energy_smooth <<", ";
    /*double energy_mvl = (MVLap * vars).norm();*/
    double energy_ls[3];
    energy_ls[0] = MTEnergy[0].norm();
    energy_ls[1] = MTEnergy[1].norm();
    energy_ls[2] = MTEnergy[2].norm();
    std::cout << "pg, " << energy_ls[0]<<", "<< energy_ls[1]<<", "<< energy_ls[2]<<", pgmax, "<<MTEnergy[0].lpNorm<Eigen::Infinity>()
		<<", "<<MTEnergy[1].lpNorm<Eigen::Infinity>()<<", "<<MTEnergy[2].lpNorm<Eigen::Infinity>()<<", ";

    double energy_el = ElEnergy.norm();
    std::cout << "el, " << energy_el << ", ";
    step_length=dx.norm();
    std::cout<<"step "<<step_length<<std::endl;
    update_mesh_properties();
    Last_Opt_Mesh=true;
}
