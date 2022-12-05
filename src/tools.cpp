#include <lsc/basic.h>
#include <lsc/tools.h>

bool triangles_coplanar(const Eigen::Vector3d &t0, const Eigen::Vector3d &t1, const Eigen::Vector3d &t2,
                        const Eigen::Vector3d &p0, const Eigen::Vector3d &p1, const Eigen::Vector3d &p2)
{
    Eigen::Vector3d norm1 = (t1 - t0).cross(t2 - t1).normalized();
    Eigen::Vector3d norm2 = (p1 - p0).cross(p2 - p1).normalized();
    Eigen::Vector3d cross = norm1.cross(norm2);

    if (cross.norm() < 1e-8)
    {
        return true;
    }
    return false;
}
bool triangles_coplanar(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const int fid1, const int fid2){
    Eigen::Vector3d t0, t1, t2, p0, p1, p2;
    t0 = V.row(F(fid1, 0));
    t1 = V.row(F(fid1, 1));
    t2 = V.row(F(fid1, 2));

    p0 = V.row(F(fid2, 0));
    p1 = V.row(F(fid2, 1));
    p2 = V.row(F(fid2, 2));
    return triangles_coplanar(t0,t1,t2,p0,p1,p2);  
}
bool segments_colinear(const Eigen::Vector3d& s1, const Eigen::Vector3d& s2){
    Eigen::Vector3d ss1=s1.normalized();
    Eigen::Vector3d ss2=s2.normalized();
    double dot=ss1.dot(ss2);
    if(fabs(dot)>1-1e-8){
        return true;
    }
    return false;
}
// -1: boundary edge
// 0 not coplanar
// 1 co-planar
int edge_is_coplanar(const CGMesh &lsmesh, const CGMesh::EdgeHandle &edge_middle, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
{

    OpenMesh::HalfedgeHandle e1 = lsmesh.halfedge_handle(edge_middle, 0);
    OpenMesh::HalfedgeHandle e2 = lsmesh.halfedge_handle(edge_middle, 1);

    int f1 = lsmesh.face_handle(e1).idx();
    int f2 = lsmesh.face_handle(e2).idx();
    if (f1 < 0 || f2 < 0)
    {
        return -1;
    }
    Eigen::Vector3d t0 = V.row(F(f1, 0));
    Eigen::Vector3d t1 = V.row(F(f1, 1));
    Eigen::Vector3d t2 = V.row(F(f1, 2));
    Eigen::Vector3d p0 = V.row(F(f2, 0));
    Eigen::Vector3d p1 = V.row(F(f2, 1));
    Eigen::Vector3d p2 = V.row(F(f2, 2));
    return triangles_coplanar(t0, t1, t2, p0, p1, p2);
}

std::vector<Trip> to_triplets(spMat &M)
{
    std::vector<Trip> v;
    v.reserve(M.rows());
    for (int i = 0; i < M.outerSize(); i++)
    {
        for (spMat::InnerIterator it(M, i); it; ++it)
        {
            v.push_back(Trip(it.row(), it.col(), it.value()));
        }
    }

    return v;
}
spMat sparse_vec_to_sparse_maxrix(Efunc &vec)
{
    spMat mat;
    mat.resize(1, vec.size());
    std::vector<Trip> triplets;
    for (Efunc::InnerIterator it(vec); it; ++it)
    {
        triplets.push_back(Trip(0, it.index(), it.value()));
    }
    mat.setFromTriplets(triplets.begin(), triplets.end());
    // std::cout << "converted" << std::endl;
    return mat;
}
Efunc sparse_mat_col_to_sparse_vec(const spMat &mat, const int col)
{
    Efunc vec;
    vec.resize(mat.rows());
    assert(!mat.IsRowMajor);// it works for colomn major
    for (spMat::InnerIterator it(mat, col); it; ++it)
    {
        vec.coeffRef(it.index()) = it.value();
    }
    return vec;
}
Efunc dense_vec_to_sparse_vec(const Eigen::VectorXd &vec)
{
    Efunc result;
    result.resize(vec.size());
    for (int i = 0; i < vec.size(); i++)
    {
        if (vec[i] != 0)
        {
            result.coeffRef(i) = vec[i];
        }
    }
    return result;
}
void mat_col_to_triplets(const spMat &mat, const int col, const int ref, const bool inverse, std::vector<Trip> &triplets)
{
    assert(!mat.IsRowMajor);

    for (spMat::InnerIterator it(mat, col); it; ++it)
    {
        int id = it.index();
        double value = it.value();
        if (inverse)
        {
            triplets.push_back(Trip(ref, id, value));
        }
        else
        {
            triplets.push_back(Trip(id, ref, value));
        }
    }
}
void lsTools::debug_tool(int id, int id2, double value)
{
    trace_vers.resize(1);
    std::cout << "checking one pseudo geodesic" << std::endl;
    // debug init
    ver_dbg.resize(0, 0);
    ver_dbg1.resize(0, 0);
    flag_dbg = false;
    E0_dbg.resize(0, 0);
    direction_dbg.resize(0, 0);

    pnorm_dbg.resize(0, 0);
    pnorm_list_dbg.clear();
    flag_dbg = true;
    id_dbg = id2;

    // debug init end
    double start_point_para = 0.5;
    double start_angle_degree = 60;
    double target_angle = value;
    std::vector<Eigen::Vector3d> curve;
    std::vector<CGMesh::HalfedgeHandle> handles;
    trace_single_pseudo_geodesic_curve_pseudo_vertex_method(target_angle, Boundary_Edges[id], start_point_para, start_angle_degree,
                                                            curve, handles);
    flag_dbg = false;
    int nbr_midpts = curve.size() - 2;
    E0_dbg.resize(nbr_midpts, 3);
    direction_dbg.resize(nbr_midpts, 3);
    Eigen::MatrixXd binormals_real(nbr_midpts, 3);
    for (int i = 0; i < nbr_midpts; i++)
    {
        E0_dbg.row(i) = curve[i + 1];
        Eigen::Vector3d vec0 = (pseudo_vers_dbg[i + 1] - pseudo_vers_dbg[i]).normalized();
        Eigen::Vector3d vec1 = (curve[i + 2] - pseudo_vers_dbg[i + 1]).normalized();
        direction_dbg.row(i) = vec0.cross(vec1).normalized(); // bi normal of the pseudo-curve
        vec0 = (curve[i + 1] - curve[i]).normalized();
        vec1 = (curve[i + 2] - curve[i + 1]).normalized();
        binormals_real.row(i) = vec0.cross(vec1).normalized();
    }
    if (pnorm_list_dbg.size() != nbr_midpts)
    {
        std::cout << "dangerous size! " << pnorm_list_dbg.size() << ", " << nbr_midpts << std::endl;
    }
    if (target_angle > 90 - ANGLE_TOLERANCE && target_angle < 90 + ANGLE_TOLERANCE)
    {
        std::cout << "Traced geodesic" << std::endl;
    }
    else
    {
        assert(pnorm_list_dbg.size() == nbr_midpts);
        pnorm_dbg = vec_list_to_matrix(pnorm_list_dbg);
        // check the abs(consin)
        for (int i = 0; i < nbr_midpts; i++)
        {
            Eigen::Vector3d tmp_dir1 = direction_dbg.row(i);
            Eigen::Vector3d tmp_dir2 = pnorm_list_dbg[i].normalized();
            double cosin = tmp_dir1.dot(tmp_dir2);
            std::cout << i << "th cosin^2 " << cosin * cosin << std::endl;
            tmp_dir1 = binormals_real.row(i);
            cosin = tmp_dir1.dot(tmp_dir2);
            std::cout << i << "th cosin^2 " << cosin * cosin << std::endl;
        }
    }

    // TODO temporarily checking one trace line
    trace_vers[0] = curve;
    trace_hehs.push_back(handles);
}
void lsTools::debug_tool_v2(const std::vector<int> &ids, const std::vector<double> values)
{
    // cylinder_example(5, 10, 50, 30);
    // exit(0);
    int nbr_itv = ids[0]; // every nbr_itv boundary edges we shoot one curve
    if (nbr_itv < 1)
    {
        std::cout << "Please set up the parameter nbr_itv " << std::endl;
        return;
    }
    double target_angle = values[0];
    double start_angel = values[1];
    OpenMesh::HalfedgeHandle init_edge = Boundary_Edges[0];
    if (lsmesh.face_handle(init_edge).idx() >= 0)
    {
        init_edge = lsmesh.opposite_halfedge_handle(init_edge);
    }
    OpenMesh::HalfedgeHandle checking_edge = init_edge;
    int curve_id = 0;
    while (1)
    {
        ver_dbg.resize(0, 0);
        ver_dbg1.resize(0, 0);
        flag_dbg = false;
        E0_dbg.resize(0, 0);
        direction_dbg.resize(0, 0);
        pnorm_dbg.resize(0, 0);
        pnorm_list_dbg.clear();
        flag_dbg = true;
        id_dbg = ids[1];
        std::vector<Eigen::Vector3d> curve;
        std::vector<CGMesh::HalfedgeHandle> handles;
        trace_single_pseudo_geodesic_curve_pseudo_vertex_method(target_angle, checking_edge, 0.5, start_angel,
                                                                curve, handles);
        flag_dbg = false;
        curve_id++;
        if (curve_id == -1)
        {
            bool stop_flag = false;

            for (int i = 0; i < nbr_itv; i++)
            {
                checking_edge = lsmesh.next_halfedge_handle(checking_edge);
                if (checking_edge == init_edge)
                {
                    stop_flag = true;
                    break;
                }
            }
            continue;
        }
        trace_vers.push_back(curve);
        trace_hehs.push_back(handles);
        bool stop_flag = false;

        if (curve_id == -1)
        {
            break;
        }
        for (int i = 0; i < nbr_itv; i++)
        {
            checking_edge = lsmesh.next_halfedge_handle(checking_edge);
            if (checking_edge == init_edge)
            {
                stop_flag = true;
                break;
            }
        }
        if (stop_flag)
        {
            break;
        }
    }
}

void lsTools::debug_tool_v3(int id, int id2, double value)
{
    trace_vers.resize(1);
    std::cout << "checking one pseudo geodesic" << std::endl;
    // debug init
    ver_dbg.resize(0, 0);
    ver_dbg1.resize(0, 0);
    flag_dbg = false;
    E0_dbg.resize(0, 0);
    direction_dbg.resize(0, 0);
    trace_hehs.clear();
    pnorm_dbg.resize(0, 0);
    pnorm_list_dbg.clear();
    flag_dbg = true;
    id_dbg = id2;

    // debug init end
    double start_point_para = 0.5;
    double start_angle_degree = 60;
    double target_angle = value;
    std::vector<Eigen::Vector3d> curve;
    std::vector<CGMesh::HalfedgeHandle> handles;
    trace_single_pseudo_geodesic_curve(target_angle, Boundary_Edges[id], start_point_para, start_angle_degree,
                                       curve, handles);
    flag_dbg = false;
    int nbr_midpts = curve.size() - 2;
    E0_dbg.resize(nbr_midpts, 3);
    direction_dbg.resize(nbr_midpts, 3);
    Eigen::MatrixXd binormals_real(nbr_midpts, 3);
    for (int i = 0; i < nbr_midpts; i++)
    {
        E0_dbg.row(i) = curve[i + 1];
        Eigen::Vector3d vec0 = (curve[i + 1] - curve[i]).normalized();
        Eigen::Vector3d vec1 = (curve[i + 2] - curve[i + 1]).normalized();
        direction_dbg.row(i) = vec0.cross(vec1).normalized(); // bi normal of the pseudo-curve
    }

    if (target_angle > 90 - ANGLE_TOLERANCE && target_angle < 90 + ANGLE_TOLERANCE)
    {
        std::cout << "Traced geodesic" << std::endl;
    }
    else
    {

        for (int i = 0; i < nbr_midpts; i++)
        {
            Eigen::Vector3d tmp_dir1 = direction_dbg.row(i);
            Eigen::Vector3d tmp_dir2 = pnorm_list_dbg[i].normalized();
            double cosin = tmp_dir1.dot(tmp_dir2);
            std::cout << i << "th cosin^2 " << cosin * cosin << std::endl;
        }
    }

    // TODO temporarily checking one trace line
    trace_vers[0] = curve;
    trace_hehs.push_back(handles);
}
void lsTools::debug_tool_v4(const std::vector<int> &ids, const std::vector<double> values)
{
    // cylinder_example(5, 10, 50, 30);
    // exit(0);
    trace_vers.clear();
    trace_hehs.clear();
    int nbr_itv = ids[0]; // every nbr_itv boundary edges we shoot one curve
    if (nbr_itv < 1)
    {
        std::cout << "Please set up the parameter nbr_itv " << std::endl;
        return;
    }
    double target_angle = values[0];
    double start_angel = values[1];
    OpenMesh::HalfedgeHandle init_edge = Boundary_Edges[0];
    if (lsmesh.face_handle(init_edge).idx() >= 0)
    {
        init_edge = lsmesh.opposite_halfedge_handle(init_edge);
    }
    OpenMesh::HalfedgeHandle checking_edge = init_edge;
    int curve_id = 0;
    while (1)
    {
        ver_dbg.resize(0, 0);
        ver_dbg1.resize(0, 0);
        E0_dbg.resize(0, 0);
        direction_dbg.resize(0, 0);
        pnorm_dbg.resize(0, 0);
        pnorm_list_dbg.clear();
        flag_dbg = true;
        id_dbg = ids[1];
        std::vector<Eigen::Vector3d> curve;
        std::vector<CGMesh::HalfedgeHandle> handles;
        trace_single_pseudo_geodesic_curve(target_angle, checking_edge, 0.5, start_angel,
                                           curve, handles);
        flag_dbg = false;
        curve_id++;
        if (curve_id == -1)
        {
            bool stop_flag = false;

            for (int i = 0; i < nbr_itv; i++)
            {
                checking_edge = lsmesh.next_halfedge_handle(checking_edge);
                if (checking_edge == init_edge)
                {
                    stop_flag = true;
                    break;
                }
            }
            continue;
        }
        trace_vers.push_back(curve);
        trace_hehs.push_back(handles);
        bool stop_flag = false;

        if (curve_id == -1)
        {
            break;
        }
        for (int i = 0; i < nbr_itv; i++)
        {
            checking_edge = lsmesh.next_halfedge_handle(checking_edge);
            if (checking_edge == init_edge)
            {
                stop_flag = true;
                break;
            }
        }
        if (stop_flag)
        {
            break;
        }
    }
}

void extend_triplets_offset(std::vector<Trip> &triplets, const spMat &mat, int offrow, int offcol)
{
    for (int i = 0; i < mat.outerSize(); i++)
    {
        for (spMat::InnerIterator it(mat, i); it; ++it)
        {
            triplets.push_back(Trip(it.row() + offrow, it.col() + offcol, it.value()));
        }
    }
}
spMat dense_mat_list_as_sparse_diagnal(const std::vector<Eigen::MatrixXd> &mlist)
{
    spMat result;
    int outer = mlist.size();
    int inner_rows = mlist[0].rows();
    int inner_cols = mlist[0].cols();
    result.resize(outer * inner_rows, outer * inner_cols);
    std::vector<Trip> triplets;
    triplets.reserve(inner_cols * inner_rows * outer);
    for (int i = 0; i < outer; i++)
    {
        extend_triplets_offset(triplets, mlist[i].sparseView(), i * inner_rows, i * inner_cols);
    }
    result.setFromTriplets(triplets.begin(), triplets.end());
    return result;
}

std::vector<double> polynomial_simplify(const std::vector<double> &poly)
{
    std::vector<double> result = poly;
    int size = poly.size();
    for (int i = 0; i < size - 1; i++)
    {
        if (result[size - i - 1] == 0)
        {
            result.pop_back();
        }
        else
        {
            break;
        }
    }
    return result;
}

std::vector<double> polynomial_add(const std::vector<double> &poly1, const std::vector<double> &poly2)
{
    int size = std::max(poly1.size(), poly2.size());
    std::vector<double> result(size);
    for (int i = 0; i < size; i++)
    {
        bool flag1 = i < poly1.size();
        bool flag2 = i < poly2.size();
        if (flag1 && flag2)
        {
            result[i] = poly1[i] + poly2[i];
        }
        else if (flag1)
        {
            result[i] = poly1[i];
        }
        else
        {
            result[i] = poly2[i];
        }
    }
    return polynomial_simplify(result);
}
std::vector<double> polynomial_times(const std::vector<double> &poly1, const std::vector<double> &poly2)
{
    int size = poly1.size() + poly2.size() - 1;
    std::vector<double> result(size);
    for (int i = 0; i < size; i++)
    { // initialize the result
        result[i] = 0;
    }

    for (int i = 0; i < poly1.size(); i++)
    {
        for (int j = 0; j < poly2.size(); j++)
        {
            result[i + j] += poly1[i] * poly2[j];
        }
    }
    return polynomial_simplify(result);
}
std::vector<double> polynomial_times(const std::vector<double> &poly1, const double &nbr)
{
    std::vector<double> result;
    if (nbr == 0)
    {
        result.resize(1);
        result[0] = 0;
        return result;
    }
    result = poly1;
    for (int i = 0; i < result.size(); i++)
    {
        result[i] *= nbr;
    }

    return polynomial_simplify(result);
}
double polynomial_value(const std::vector<double> &poly, const double para)
{
    double result = 0;
    for (int i = 0; i < poly.size(); i++)
    {
        result += poly[i] * std::pow(para, i);
    }
    return result;
}
std::vector<double> polynomial_integration(const std::vector<double> &poly)
{
    std::vector<double> result(poly.size() + 1);
    result[0] = 0;
    for (int i = 1; i < result.size(); i++)
    {
        result[i] = poly[i - 1] / i;
    }
    return polynomial_simplify(result);
}
double polynomial_integration(const std::vector<double> &poly, const double lower, const double upper)
{
    double up = polynomial_value(polynomial_integration(poly), upper);
    double lw = polynomial_value(polynomial_integration(poly), lower);
    return up - lw;
}

Eigen::MatrixXd vec_list_to_matrix(const std::vector<Eigen::Vector3d> &vec)
{
    Eigen::MatrixXd mat;
    mat.resize(vec.size(), 3);
    for (int i = 0; i < mat.rows(); i++)
    {
        mat.row(i) = vec[i];
    }
    return mat;
}

CGMesh::HalfedgeHandle boundary_halfedge(const CGMesh& lsmesh, const CGMesh::HalfedgeHandle& boundary_edge){
    assert(lsmesh.is_boundary(lsmesh.from_vertex_handle(boundary_edge)));
    assert(lsmesh.is_boundary(lsmesh.to_vertex_handle(boundary_edge)));
    CGMesh::HalfedgeHandle result=boundary_edge;
    if(lsmesh.face_handle(result).idx()>=0){
        result=lsmesh.opposite_halfedge_handle(result);
    }
    return result;
}
void split_mesh_boundary_by_corner_detection(CGMesh &lsmesh, const Eigen::MatrixXd& V, const double threadshold_angel_degree,
                                             const std::vector<CGMesh::HalfedgeHandle> &Boundary_Edges, 
                                             std::vector<std::vector<CGMesh::HalfedgeHandle>> &boundaries)
{
    boundaries.clear();
    int vnbr=V.rows();
    SpVeci bedges;// 0: not a boundary edge; 1: not checked boundary edge;
    SpVeci corner_vers;
    corner_vers.resize(vnbr);
    int enbr=lsmesh.n_edges();
    int bsize=Boundary_Edges.size();
    std::cout<<"boundary edges size "<<Boundary_Edges.size()<<std::endl;
    bedges.resize(enbr);
    for(int i=0;i<Boundary_Edges.size();i++){
        CGMesh::EdgeHandle eh=lsmesh.edge_handle(Boundary_Edges[i]);
        bedges.coeffRef(eh.idx())=1; // mark every boundary edges
    }
    
    
    int nbr_corners=0;
    for (int i = 0; i < bsize; i++)
    { // iterate all the boundary edges
        CGMesh::HalfedgeHandle start_he=boundary_halfedge(lsmesh, Boundary_Edges[i]);
        CGMesh::HalfedgeHandle next_he=lsmesh.next_halfedge_handle(start_he);
        // std::cout<<"this face id "<<lsmesh.face_handle(start_he).idx()<<" opposite "<<lsmesh.face_handle(lsmesh.opposite_halfedge_handle(start_he))<<std::endl;
        // std::cout<<"next face id "<<lsmesh.face_handle(next_he).idx()<<" opposite "<<lsmesh.face_handle(lsmesh.opposite_halfedge_handle(next_he))<<std::endl;
        assert(lsmesh.face_handle(next_he).idx()<0);
        Eigen::Vector3d p0=V.row(lsmesh.from_vertex_handle(start_he).idx());
        Eigen::Vector3d p1=V.row(lsmesh.to_vertex_handle(start_he).idx());
        Eigen::Vector3d p2=V.row(lsmesh.to_vertex_handle(next_he).idx());
        Eigen::Vector3d dir1=(p1-p0).normalized();
        Eigen::Vector3d dir2=(p2-p1).normalized();
        int id_middle_ver=lsmesh.to_vertex_handle(start_he).idx();
        double dot_product=dir1.dot(dir2);
        double real_angle=acos(dot_product) * 180 / LSC_PI;
        // std::cout<<"real angle "<<real_angle<<std::endl;
        if (real_angle < 180 - threadshold_angel_degree)
        {
            continue;
        }
        else{// current start_he and next_he are edges of a corner.
            corner_vers.coeffRef(id_middle_ver)=1;
            std::cout<<"find one corner"<<std::endl;
            nbr_corners++; 
        }
    }
    std::cout<<"nbr of corners "<<nbr_corners<<std::endl;
    if (nbr_corners == 0)
    {
        SpVeci recorded_edges;
        recorded_edges.resize(enbr);
        for (int i = 0; i < enbr; i++)
        {
            if (recorded_edges.coeffRef(i) == 0 && bedges.coeffRef(i) == 1)
            { // if this is a boundary edge, and not checked yet
                recorded_edges.coeffRef(i) = 1;
                CGMesh::HalfedgeHandle start_he = lsmesh.halfedge_handle(lsmesh.edge_handle(i), 0);
                start_he = boundary_halfedge(lsmesh, start_he);
                // CGMesh::HalfedgeHandle start_copy = start_he;
                std::vector<CGMesh::HalfedgeHandle> loop;
                loop.push_back(start_he);
                for (int j = 0; j < bsize; j++)
                {
                    CGMesh::HalfedgeHandle next_he = lsmesh.next_halfedge_handle(start_he);
                    int edge_id=lsmesh.edge_handle(next_he).idx();
                    if (recorded_edges.coeffRef(edge_id)==0)// if the next edge is not checked yet
                    {
                        recorded_edges.coeffRef(edge_id)=1;
                        loop.push_back(next_he);
                        start_he = next_he;
                    }
                    else
                    {
                        break;
                    }
                }
                boundaries.push_back(loop);
            }
        }
        return;
    }
    SpVeci startlist;// this list marks if this corner vertex is already used.
    startlist.resize(vnbr);
    // the number of corner points is the same as (or more than, in some extreme cases) the number of segments

    for(int i=0;i<bsize;i++){ 
        CGMesh::HalfedgeHandle start_he=boundary_halfedge(lsmesh, Boundary_Edges[i]);
        int id_start_ver=lsmesh.from_vertex_handle(start_he).idx();
        

        if(corner_vers.coeffRef(id_start_ver)==1&&startlist.coeffRef(id_start_ver)==0){//if this point is a corner and not used as a start point
            startlist.coeffRef(id_start_ver) = 1;// mark it as used point
            std::vector<CGMesh::HalfedgeHandle> loop;
            loop.push_back(start_he);
            for(int j=0;j<bsize;j++){
                CGMesh::HalfedgeHandle next_he=lsmesh.next_halfedge_handle(start_he);
                int id_end_ver=lsmesh.to_vertex_handle(next_he).idx();
                loop.push_back(next_he);
                
                if(corner_vers.coeffRef(id_end_ver)==1){// the end point is a vertex
                    break;
                }
                start_he=next_he;
            }
            boundaries.push_back(loop);
        }
       
    }
}


Eigen::Vector3d sphere_function(double r, double theta, double phi)
{
    double x = r * sin(theta) * cos(phi);
    double y = r * sin(theta) * sin(phi);
    double z = r * cos(theta);
    return Eigen::Vector3d(x, y, z);
}
// create a triangle mesh sphere.
#include <igl/write_triangle_mesh.h>
void sphere_example(double radius, double theta, double phi, int nt, int np)
{

    Eigen::MatrixXd ver;
    Eigen::MatrixXi faces;
    ver.resize(nt * np, 3);
    faces.resize(2 * (nt - 1) * (np - 1), 3);
    int verline = 0;
    double titv = 2 * theta / (nt - 1);
    double pitv = 2 * phi / (np - 1);
    for (int i = 0; i < nt; i++)
    {
        for (int j = 0; j < np; j++)
        {
            double upara = 0.5 * 3.1415926 - theta + i * titv;
            double vpara = -phi + j * pitv;
            ver.row(verline) = sphere_function(radius, upara, vpara);
            verline++;
        }
    }
    faces.resize(2 * (nt - 1) * (np - 1), 3);
    int fline = 0;
    for (int i = 0; i < nt - 1; i++)
    {
        for (int j = 0; j < np - 1; j++)
        {
            int id0 = np * i + j;
            int id1 = np * (i + 1) + j;
            int id2 = np * (i + 1) + j + 1;
            int id3 = np * i + j + 1;
            faces.row(fline) = Eigen::Vector3i(id0, id1, id2);
            faces.row(fline + 1) = Eigen::Vector3i(id0, id2, id3);
            fline += 2;
        }
    }
    std::string path("/Users/wangb0d/bolun/D/vs/levelset/level-set-curves/data/");
    igl::write_triangle_mesh(path + "sphere_" + std::to_string(radius) + "_" + std::to_string(theta) + "_" +
                                 std::to_string(phi) + "_" + std::to_string(nt) + "_" + std::to_string(np) + ".obj",
                             ver, faces);
    std::cout << "sphere mesh file saved " << std::endl;
}
void cylinder_example(double radius, double height, int nr, int nh)
{
    Eigen::MatrixXd ver;
    Eigen::MatrixXi faces;
    ver.resize(nh * nr, 3);
    faces.resize(2 * (nh - 1) * nr, 3);
    int verline = 0;
    double hitv = height / (nh - 1);
    double ritv = 2 * LSC_PI / nr;
    for (int i = 0; i < nr; i++)
    {
        double angle = ritv * i;
        double x = cos(angle) * radius;
        double y = sin(angle) * radius;
        for (int j = 0; j < nh; j++)
        {
            double z = j * hitv;
            ver.row(verline) << x, y, z;
            verline++;
        }
    }

    int fline = 0;
    for (int i = 0; i < nr; i++)
    {
        if (i < nr - 1)
        {
            for (int j = 0; j < nh - 1; j++)
            {
                int id0 = nh * i + j;
                int id1 = nh * (i + 1) + j;
                int id2 = nh * (i + 1) + j + 1;
                int id3 = nh * i + j + 1;
                faces.row(fline) = Eigen::Vector3i(id0, id1, id2);
                faces.row(fline + 1) = Eigen::Vector3i(id0, id2, id3);
                fline += 2;
            }
        }
        else
        {
            for (int j = 0; j < nh - 1; j++)
            {
                int id0 = nh * i + j;
                int id1 = j;
                int id2 = j + 1;
                int id3 = nh * i + j + 1;
                faces.row(fline) = Eigen::Vector3i(id0, id1, id2);
                faces.row(fline + 1) = Eigen::Vector3i(id0, id2, id3);
                fline += 2;
            }
        }
    }
    std::string path("/Users/wangb0d/bolun/D/vs/levelset/level-set-curves/data/");
    igl::write_triangle_mesh(path + "cylinder_" + std::to_string(radius) + "_" + std::to_string(height) + "_" +
                                 std::to_string(nr) + "_" + std::to_string(nh) + ".obj",
                             ver, faces);
    std::cout << "sphere mesh file saved " << std::endl;
}

void cylinder_open_example(double radius, double height, int nr, int nh)
{
    Eigen::MatrixXd ver;
    Eigen::MatrixXi faces;
    ver.resize(nh * nr, 3);
    faces.resize(2 * (nh - 1) * (nr-1), 3);
    int verline = 0;
    double hitv = height / (nh - 1);
    double ritv =  LSC_PI / (nr-1);// we only draw half a cylinder
    for (int i = 0; i < nr; i++)
    {
        double angle = ritv * i;
        double x = cos(angle) * radius;
        double y = sin(angle) * radius;
        for (int j = 0; j < nh; j++)
        {
            double z = j * hitv;
            ver.row(verline) << x, y, z;
            verline++;
        }
    }

    int fline = 0;
    for (int i = 0; i < nr; i++)
    {
        if (i < nr - 1)
        {
            for (int j = 0; j < nh - 1; j++)
            {
                int id0 = nh * i + j;
                int id1 = nh * (i + 1) + j;
                int id2 = nh * (i + 1) + j + 1;
                int id3 = nh * i + j + 1;
                faces.row(fline) = Eigen::Vector3i(id0, id1, id2);
                faces.row(fline + 1) = Eigen::Vector3i(id0, id2, id3);
                fline += 2;
            }
        }
        // else
        // {
        //     for (int j = 0; j < nh - 1; j++)
        //     {
        //         int id0 = nh * i + j;
        //         int id1 = j;
        //         int id2 = j + 1;
        //         int id3 = nh * i + j + 1;
        //         faces.row(fline) = Eigen::Vector3i(id0, id1, id2);
        //         faces.row(fline + 1) = Eigen::Vector3i(id0, id2, id3);
        //         fline += 2;
        //     }
        // }
    }
    std::string path("/Users/wangb0d/bolun/D/vs/levelset/level-set-curves/data/");
    igl::write_triangle_mesh(path + "oc_" + std::to_string(radius) + "_" + std::to_string(height) + "_" +
                                 std::to_string(nr) + "_" + std::to_string(nh) + ".obj",
                             ver, faces);
    std::cout << "cylinder file saved " << std::endl;
}

#include<igl/file_dialog_open.h>
#include<igl/file_dialog_save.h>
bool save_levelset(const Eigen::VectorXd &ls, const Eigen::MatrixXd& binormals){
    std::cout<<"SAVING LevelSet..."<<std::endl;
    std::string fname = igl::file_dialog_save();
    std::ofstream file;
    file.open(fname);
    for(int i=0;i<ls.size();i++){
        file<<ls[i]<<std::endl;
    }
    file.close();
    std::cout<<"LevelSet SAVED"<<std::endl;
    if (binormals.rows() > 0)
    {
        std::cout << "SAVING BI-NORMALS..." << std::endl;
        fname = igl::file_dialog_save();
        file.open(fname);
        for (int i = 0; i < binormals.rows(); i++)
        {
            file << binormals(i, 0) << "," << binormals(i, 1) << "," << binormals(i, 2) << std::endl;
        }
        file.close();
        std::cout << "BI-NORMALS SAVED, size " << binormals.rows() << std::endl;
    }

    return true;
}
bool save_levelset(const Eigen::VectorXd &ls){
    std::cout<<"SAVING LevelSet..."<<std::endl;
    std::string fname = igl::file_dialog_save();
    std::ofstream file;
    file.open(fname);
    for(int i=0;i<ls.size();i++){
        file<<ls[i]<<std::endl;
    }
    file.close();
    std::cout<<"LevelSet SAVED"<<std::endl;

    return true;
}
bool read_levelset(Eigen::VectorXd &ls){
    std::string fname = igl::file_dialog_open();
    std::cout<<"reading "<<fname<<std::endl;
    if (fname.length() == 0)
        return false;

    std::ifstream infile;
    std::vector<double> results;

    infile.open(fname);
    if (!infile.is_open())
    {
        std::cout << "Path Wrong!!!!" << std::endl;
        std::cout << "path, " << fname << std::endl;
        return false;
    }

    int l = 0;
    while (infile) // there is input overload classfile
    {
        std::string s;
        if (!getline(infile, s))
            break;

        if (s[0] != '#')
        {
            std::istringstream ss(s);
            std::string record;

            while (ss)
            {
                std::string line;
                if (!getline(ss, line, ','))
                    break;
                try
                {

                    record = line;

                }
                catch (const std::invalid_argument e)
                {
                    std::cout << "NaN found in file " << fname
                              <<  std::endl;
                    e.what();
                }
            }

            double x = std::stod(record);
            results.push_back(x);
        }
    }

    ls.resize(results.size());
    for(int i=0;i<results.size();i++){
        ls[i]=results[i];
    }
    if (!infile.eof())
    {
        std::cerr << "Could not read file " << fname << "\n";
    }
    std::cout<<fname<<" get readed"<<std::endl;
    return true;
}


bool read_bi_normals(Eigen::MatrixXd &bn){
    std::string fname = igl::file_dialog_open();
    std::cout<<"reading "<<fname<<std::endl;
    if (fname.length() == 0)
        return false;

    std::ifstream infile;
    std::vector<Eigen::Vector3d> results;

    infile.open(fname);
    if (!infile.is_open())
    {
        std::cout << "Path Wrong!!!!" << std::endl;
        std::cout << "path, " << fname << std::endl;
        return false;
    }

    int l = 0;
    while (infile) // there is input overload classfile
    {
        std::string s;
        if (!getline(infile, s))
            break;

        if (s[0] != '#')
        {
            std::istringstream ss(s);
            std::array<std::string, 3> record;
            int c = 0;
            while (ss)
            {
                std::string line;
                if (!getline(ss, line, ','))
                    break;
                try
                {

                    record[c] = line;
                    c++;

                }
                catch (const std::invalid_argument e)
                {
                    std::cout << "NaN found in file " << fname
                              <<  std::endl;
                    e.what();
                }
            }

            double x = std::stod(record[0]);
            double y = std::stod(record[1]);
            double z = std::stod(record[2]);
            results.push_back(Eigen::Vector3d(x, y, z));
        }
    }

    bn.resize(results.size(), 3);
    for(int i=0;i<results.size();i++){
        bn.row(i)=results[i];
    }
    if (!infile.eof())
    {
        std::cerr << "Could not read file " << fname << "\n";
    }
    std::cout<<fname<<" get readed"<<std::endl;
    return true;
}

// get vid1, vid2, vidA and vidB.
// here we do not assume fvalues[vid1]>fvalues[vid2], we need to consider it outside
std::array<int,4> get_vers_around_edge(CGMesh& lsmesh, int edgeid, int& fid1, int &fid2){
    CGMesh::EdgeHandle edge_middle=lsmesh.edge_handle(edgeid);
    CGMesh::HalfedgeHandle he = lsmesh.halfedge_handle(edge_middle, 0);
    int vid1 = lsmesh.to_vertex_handle(he).idx();
    int vid2 = lsmesh.from_vertex_handle(he).idx();
    
    CGMesh::HalfedgeHandle nexthandle=lsmesh.next_halfedge_handle(he);
    int vidA=lsmesh.to_vertex_handle(nexthandle).idx();
    CGMesh::HalfedgeHandle oppohandle=lsmesh.opposite_halfedge_handle(he);
    fid1=lsmesh.face_handle(he).idx();
    fid2=lsmesh.face_handle(oppohandle).idx();
    CGMesh::HalfedgeHandle opnexthandle=lsmesh.next_halfedge_handle(oppohandle);
    int vidB=lsmesh.to_vertex_handle(opnexthandle).idx();
    std::array<int,4> result;
    result[0]=vid1;
    result[1]=vid2;
    result[2]=vidA;
    result[3]=vidB;
    return result;
}

void get_level_set_sample_values(const Eigen::VectorXd &ls, const int nbr, Eigen::VectorXd &values){
    double vmax=ls.maxCoeff();
    double vmin=ls.minCoeff();
    double itv=(vmax-vmin)/(nbr+1);
    values.resize(nbr);
    for(int i=0;i<nbr;i++){
        values[i]=vmin+(i+1)*itv;
    }
    return;
}
void get_level_set_sample_values_even_pace(const Eigen::VectorXd &ls, const int nbr, const double pace, Eigen::VectorXd &values){
    double vmax=ls.maxCoeff();
    double vmin=ls.minCoeff();
    double remain=(vmax-vmin)-(nbr-1)*pace;
    remain/=2; 
    if(remain<0){
        std::cout<<"computing wrong"<<std::endl;
    }
    values.resize(nbr);
    for(int i=0;i<nbr;i++){
        values[i] = remain + vmin + i * pace;
    }
    return;
}


// 
void get_iso_lines(const Eigen::MatrixXd &V,
                   const Eigen::MatrixXi &F, const Eigen::VectorXd &ls, const double value,
                   std::vector<Eigen::Vector3d> &poly0, std::vector<Eigen::Vector3d> &poly1,
                   std::vector<int> &fids)
{
    for (int i = 0; i < F.rows(); i++)
    {
        Eigen::Vector3d E0, E1;
        bool found = find_one_ls_segment_on_triangle(value, F, V,
                                                     ls, i, E0, E1);
        if(found){
            poly0.push_back(E0);
            poly1.push_back(E1);
            fids.push_back(i);
        }
    }
}
bool halfedge_has_value(const CGMesh &lsmesh, const Eigen::MatrixXd &V,
                         const Eigen::VectorXd &ls, const CGMesh::HalfedgeHandle &hd, const double value, Eigen::Vector3d &pt)
{
    int id_f=lsmesh.from_vertex_handle(hd).idx();
    int id_t=lsmesh.to_vertex_handle(hd).idx();
    double value_f=ls[id_f];
    double value_t=ls[id_t];
    double t=get_t_of_value(value, value_f, value_t);
    if(t<0||t>1){
        return false;
    }
    Eigen::Vector3d ver_f=V.row(id_f);
    Eigen::Vector3d ver_t= V.row(id_t);
    pt=get_3d_ver_from_t(t, ver_f, ver_t);
    // std::cout<<"has value "<<value<<", vf "<<value_f<<" vt "<<value_t<<", idf "<<id_f<<", idt "<<id_t<<std::endl;
    return true;
}

// the mesh provides the connectivity, loop is the boundary loop
// left_large indicates if the from ver of each boundary edge is larger value
void get_iso_lines(const CGMesh &lsmesh, const std::vector<CGMesh::HalfedgeHandle>& loop, const Eigen::MatrixXd &V,
                   const Eigen::MatrixXi &F, const Eigen::VectorXd &ls, double value, std::vector<std::vector<Eigen::Vector3d>> &pts,
                   std::vector<bool>& left_large)
{
    pts.clear();
    left_large.clear();
    std::vector<CGMesh::HalfedgeHandle> has_value; // these boundary edges has value
    std::vector<Eigen::Vector3d> start_pts; // these are the computed start/end points
    for(int i=0;i<loop.size();i++){ // find the points on boundary
        
        CGMesh::HalfedgeHandle bhd;
        if (lsmesh.face_handle(loop[i]).idx()<0)
        {
            bhd=loop[i];
        }
        else{
            bhd=lsmesh.opposite_halfedge_handle(loop[i]);
        }
        
        assert(lsmesh.face_handle(bhd).idx()<0);
        Eigen::Vector3d intersection;
        bool found=halfedge_has_value(lsmesh,V,ls, bhd, value,intersection);
        if(!found){
            continue;
        }
        has_value.push_back(bhd);
        start_pts.push_back(intersection);
        int id0=lsmesh.from_vertex_handle(bhd).idx();
        int id1=lsmesh.to_vertex_handle(bhd).idx();
        if(ls[id0]>ls[id1]){
            left_large.push_back(true);
        }
        else{
            left_large.push_back(false);
        }
    }
    int bsize=has_value.size();
    // std::vector<bool> checked(bsize); // show if this boundary edge is checked
    // for(int i=0;i<bsize;i++){
    //     checked[i]=false;
    // }

    for(int i=0;i<bsize;i++){
        // if(checked[i]){
        //     continue;
        // }
        // checked[i]=true;
        if(left_large[i]==false){ // only count for the shoouting in boundaries
            continue;
        }

        std::vector<Eigen::Vector3d> tmpts; // the vers for each polyline
        CGMesh::HalfedgeHandle thd=has_value[i]; // the start handle
        tmpts.push_back(start_pts[i]);
        while(1){
            CGMesh::HalfedgeHandle ophd=lsmesh.opposite_halfedge_handle(thd);
            int fid = lsmesh.face_handle(ophd).idx();
            if(fid<0){ // if already reached a boundary
                for(int j=0;j<bsize;j++){
                    if(ophd==has_value[j]){
                        // if(checked[j]){
                        //     std::cout<<"TOPOLOGY WRONG"<<std::endl;
                        // }
                        // checked[j]=true;
                        break;
                    }
                }
                pts.push_back(tmpts);
                break;
            }
            int id_f=lsmesh.from_vertex_handle(ophd).idx();
            int id_t=lsmesh.to_vertex_handle(ophd).idx();
            CGMesh::HalfedgeHandle prev=lsmesh.prev_halfedge_handle(ophd);
            CGMesh::HalfedgeHandle next=lsmesh.next_halfedge_handle(ophd);
            int bid=lsmesh.from_vertex_handle(prev).idx();
            if(bid==id_f||bid==id_t){
                std::cout<<"wrong topology"<<std::endl;
            }
            if(ls[bid]==value){// a degenerate case where the value is on a vertex of the mesh
                value+=1e-16;
            }
            Eigen::Vector3d intersect;
            bool found =halfedge_has_value(lsmesh, V, ls, prev,value, intersect);
            if(found){
                tmpts.push_back(intersect);
                thd=prev;
                continue;
            }
            found = halfedge_has_value(lsmesh, V, ls, next, value, intersect);
            if(found){
                tmpts.push_back(intersect);
                thd=next;
                continue;
            }
            std::cout<<"Both two halfedges did not find the correct value, error"<<std::endl;
        }
    }

}
bool two_triangles_connected(const Eigen::MatrixXi& F, const int fid0, const int fid1){
    if(fid0==fid1){
        return true;
    }
    for(int i=0;i<3;i++){
        int vid0=F(fid0,i);
        for(int j=0;j<3;j++){
            int vid1=F(fid1,j);
            if(vid0==vid1){
                return true;
            }
        }
    }
    return false;
}

#include<igl/segment_segment_intersect.h>
// TODO this is quadratic computation, need use tree structure
bool get_polyline_intersection(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                               const std::vector<Eigen::Vector3d> &poly0_e0, const std::vector<Eigen::Vector3d> &poly0_e1,
                               const std::vector<Eigen::Vector3d> &poly1_e0, const std::vector<Eigen::Vector3d> &poly1_e1,
                               const std::vector<int> &fid0, const std::vector<int> &fid1, Eigen::Vector3d &result)
{
    assert(poly0_e0.size() == poly0_e1.size());
    assert(poly1_e0.size() == poly1_e1.size());
    for (int i = 0; i < poly0_e0.size(); i++)
    {
        int f0 = fid0[i];
        Eigen::Vector3d e0 = poly0_e0[i];
        Eigen::Vector3d dir0 = poly0_e1[i] - poly0_e0[i];
        for (int j = 0; j < poly1_e0.size(); j++)
        {
            Eigen::Vector3d e1 = poly1_e0[j];
            Eigen::Vector3d dir1 = poly1_e1[j] - poly1_e0[j];
            int f1 = fid1[j];

            if (!two_triangles_connected(F, f0, f1))
            {
                continue;
            }
            double u, v;
            bool intersect = igl::segment_segment_intersect(e0, dir0, e1, dir1, u, v);
            if (intersect)
            {
                result = e0 + u * dir0;
                return true;
            }
        }
    }
    return false;
}
// diff means the input is different from output
Eigen::VectorXi filter_short_quad_rows(const Eigen::VectorXi &quads, const int threads, bool &diff)
{
    int size = quads.size();
    std::vector<std::vector<int>> deleted;
    std::vector<int> segs;
    diff = true;
    for (int i = 0; i < size; i++)
    {
        if (quads[i] < 0) // it is not a ver
        {
            if (segs.size() <= threads && segs.size() > 0)
            {
                deleted.push_back(segs);
                segs.clear();
            }
        }
        else
        {
            segs.push_back(i);
        }
    }
    if (segs.size() <= threads && segs.size() > 0)
    {
        deleted.push_back(segs);
        segs.clear();
    }

    Eigen::VectorXi result = quads;
    if (deleted.size() == 0)
    {
        diff = false;
    }

    for (int i = 0; i < deleted.size(); i++)
    {
        for (int j = 0; j < deleted[i].size(); j++)
        {
            result[deleted[i][j]] = -1;
        }
    }
    return result;
}

Eigen::MatrixXi remove_short_quads_in_face_mat(const Eigen::MatrixXi& qds, const int threads){
    int row = qds.rows();
    int col = qds.cols();
    Eigen::MatrixXi result = qds;
    bool diff;
    for(int i=0;i<row;i++){// we only delete short rows, since each row means same iso line
        Eigen::VectorXi tmp = qds.row(i);
        result.row(i)=filter_short_quad_rows(tmp, threads, diff);
    }
    return result;

}
void get_quad_left_right_ver_id_from_web_mat(const Eigen::MatrixXi &mat, const bool rowbased, const int row, const int &col_start, const int col_end,
                                             Eigen::VectorXi &left, Eigen::VectorXi &right)
{
    int nbr_elements = col_end - col_start + 1;
    if (rowbased)
    {
        left = mat.row(row).segment(col_start, nbr_elements);
        right = mat.row(row + 1).segment(col_start, nbr_elements);
    }
    else
    {
        left = mat.col(row).segment(col_start, nbr_elements);
        right = mat.col(row + 1).segment(col_start, nbr_elements);
    }
}

void get_row_and_col_ver_ids_from_web_mat(const Eigen::MatrixXi &mat, const Eigen::MatrixXi &quads,
                                          std::vector<Eigen::VectorXi> &vrl, std::vector<Eigen::VectorXi> &vrr, std::vector<Eigen::VectorXi> &vcl,
                                          std::vector<Eigen::VectorXi> &vcr)
{
    std::vector<int> ids;
    int maxsize1;
    int maxsize2;
    maxsize1 = 0;
    maxsize2=0;
    vrl.clear();
    vrr.clear();
    vcl.clear();
    vcr.clear();
    for (int i = 0; i < quads.rows(); i++)
    {
        ids.clear();
        for (int j = 0; j < quads.cols(); j++)
        { // make sure there are only 2 elements, recording the first and last non -1 id
            if (quads(i, j) >= 0 && ids.size() == 0)
            {
                ids.push_back(j);
                ids.push_back(j);
            }
            if (quads(i, j) >= 0)
            {
                ids[1] = j;
            }
        }

        if (ids.size() > 0)
        {
            if (ids[1] - ids[0] + 1 > maxsize1) // ids[1] - ids[0] is the number of quads, meaning nbr of vers is +1
            {
                maxsize1 = ids[1] - ids[0] + 1;
            }
        }
        if (ids.size() > 0)
        {
            Eigen::VectorXi left, right;
            get_quad_left_right_ver_id_from_web_mat(mat, true, i, ids[0], ids[1], left, right);
            vrl.push_back(left);
            vrr.push_back(right);
        }
    }

    for (int i = 0; i < quads.cols(); i++)
    {
        ids.clear();
        for (int j = 0; j < quads.rows(); j++)
        { // make sure there are only 2 elements, recording the first and last non -1 id
            if (quads(j, i) >= 0 && ids.size() == 0)
            {
                ids.push_back(j);
                ids.push_back(j);
            }
            if (quads(j, i) >= 0)
            {
                ids[1] = j;
            }
        }

        if (ids.size() > 0)
        {
            if (ids[1] - ids[0] + 1 > maxsize2) // ids[1] - ids[0] is the number of quads, meaning nbr of vers is +1
            {
                maxsize2 = ids[1] - ids[0] + 1;
            }
        }
        if (ids.size() > 0)
        {
            Eigen::VectorXi left, right;
            get_quad_left_right_ver_id_from_web_mat(mat, false, i, ids[0], ids[1], left, right);
            vcl.push_back(left);
            vcr.push_back(right);
        }
    }
    Eigen::VectorXi vecbase1 = Eigen::VectorXi::Ones(maxsize1) * -1;
    Eigen::VectorXi vecbase2 = Eigen::VectorXi::Ones(maxsize2) * -1;
    for (int i = 0; i < vrl.size(); i++)
    {
        Eigen::VectorXi tmp = vecbase1;
        int size = vrl[i].size();
        tmp.segment(0, size) = vrl[i];
        vrl[i] = tmp;

        tmp = vecbase1;
        tmp.segment(0, size) = vrr[i];
        vrr[i] = tmp;
    }

    for (int i = 0; i < vcl.size(); i++)
    {
        Eigen::VectorXi tmp = vecbase2;
        int size = vcl[i].size();
        tmp.segment(0, size) = vcl[i];
        vcl[i] = tmp;

        tmp = vecbase2;
        tmp.segment(0, size) = vcr[i];
        vcr[i] = tmp;
    }
}
void extract_web_from_index_mat(const Eigen::MatrixXi &mat, Eigen::MatrixXi &F, const int threads, std::vector<Eigen::VectorXi> &vrl,
                                std::vector<Eigen::VectorXi> &vrr, std::vector<Eigen::VectorXi> &vcl,
                                std::vector<Eigen::VectorXi> &vcr)
{
    std::array<int, 4> face;
    std::vector<std::array<int, 4>> tface, filted_face;
    int row = mat.rows();
    int col = mat.cols();
    Eigen::MatrixXi qds=Eigen::MatrixXi::Ones(row - 1, col - 1) * -1;// quad mat initialized as -1
    tface.reserve(mat.rows() * mat.cols());
    filted_face.reserve(mat.rows() * mat.cols());
    for(int i=0;i<mat.rows()-1;i++){
        for(int j=0;j<mat.cols()-1;j++){
            int f0=mat(i,j);
            int f1=mat(i,j+1);
            int f2=mat(i+1,j+1);
            int f3=mat(i+1, j);
            if(f0<0||f1<0||f2<0||f3<0){
                continue;
            }
            face={f0,f1,f2,f3};
            qds(i,j)=tface.size();
            tface.push_back(face);

        }
    }
    if (threads > 0)
    {
        qds = remove_short_quads_in_face_mat(qds, threads);
        for (int i = 0; i < qds.rows(); i++)
        {
            for (int j = 0; j < qds.cols(); j++)
            {
                if (qds(i, j) >= 0)
                {
                    filted_face.push_back(tface[qds(i, j)]);
                }
            }
        }
    }
    else{
        filted_face = tface;
    }

    F.resize(filted_face.size(),4);
    for(int i=0;i<filted_face.size();i++){
        F(i,0)=filted_face[i][0];
        F(i,1)=filted_face[i][1];
        F(i,2)=filted_face[i][2];
        F(i,3)=filted_face[i][3];

    }
    get_row_and_col_ver_ids_from_web_mat(mat, qds, vrl, vrr, vcl, vcr);
}
void extract_web_mxn(const int rows, const int cols, Eigen::MatrixXi& F){
    std::array<int, 4> face;
    std::vector<std::array<int,4>> tface; 
    tface.reserve(rows*cols);
    for(int i=0;i<rows-1;i++){
        for(int j=0;j<cols-1;j++){
            int f0 = cols * i + j;
            int f1 = cols * (i + 1) + j;
            int f2 = cols * (i + 1) + (j + 1);
            int f3 = cols * i + j + 1;
            face={f0,f1,f2,f3};
            tface.push_back(face);
        }
    }
    F.resize(tface.size(),4);
    for(int i=0;i<tface.size();i++){
        F(i,0)=tface[i][0];
        F(i,1)=tface[i][1];
        F(i,2)=tface[i][2];
        F(i,3)=tface[i][3];

    }
}
void write_one_info_file(const std::string &fname, const std::vector<Eigen::VectorXi> &current, std::ofstream &file)
{
    file.open(fname);
    for (int i = 0; i < current.size(); i++)
    {
        for (int j = 0; j < current[i].size(); j++)
        {
            if (j == current[i].size() - 1)
            {
                file << current[i][j] << std::endl;
            }
            else
            {
                file << current[i][j] << ",";
            }
        }
    }
    file.close();
}
void save_quad_left_right_info(const std::string &namebase, std::vector<Eigen::VectorXi> &vrl,
                               std::vector<Eigen::VectorXi> &vrr, std::vector<Eigen::VectorXi> &vcl,
                               std::vector<Eigen::VectorXi> &vcr)
{
    std::ofstream file;
    std::vector<Eigen::VectorXi> current;
    std::string fname = namebase + "_rowleft.csv";
    current=vrl;
    write_one_info_file(fname, current, file);
    
    fname = namebase + "_rowright.csv";
    current=vrr;
    write_one_info_file(fname, current, file);

    fname = namebase + "_colleft.csv";
    current=vcl;
    write_one_info_file(fname, current, file);

    fname = namebase + "_colright.csv";
    current=vcr;
    write_one_info_file(fname, current, file);
}
// threadshold_nbr is the nbr of quads we want to discard when too few quads in a row
void extract_levelset_web(const CGMesh &lsmesh, const Eigen::MatrixXd &V,
                          const Eigen::MatrixXi &F, const Eigen::VectorXd &ls0, const Eigen::VectorXd &ls1,
                          const int expect_nbr_ls0, const int expect_nbr_ls1, const int threadshold_nbr,
                          Eigen::MatrixXd &vers, Eigen::MatrixXi &Faces, bool even_pace=false)
{
    std::vector<std::vector<Eigen::Vector3d>> ivs;                                    // the intersection vertices
    Eigen::MatrixXi gridmat;                                                          // the matrix for vertex
    Eigen::VectorXd lsv0, lsv1;                                                       // the extracted level set values;
    std::vector<std::vector<Eigen::Vector3d>> poly0_e0, poly0_e1, poly1_e0, poly1_e1; // polylines
    std::vector<std::vector<int>> fid0, fid1;                                         // face ids of each polyline segment
    std::vector<Eigen::Vector3d> verlist;
    if(ls0.size()!=ls1.size()){
        std::cout<<"ERROR, Please use the correct level sets"<<std::endl;
    }
    int nbr_ls0, nbr_ls1;
    if(!even_pace){
        nbr_ls0=expect_nbr_ls0;
        nbr_ls1=expect_nbr_ls1;
        get_level_set_sample_values(ls0, nbr_ls0, lsv0);
        get_level_set_sample_values(ls1, nbr_ls1, lsv1);
    }
    else{
        double pace = std::min((ls0.maxCoeff() - ls0.minCoeff()) / (expect_nbr_ls0 + 1), (ls1.maxCoeff() - ls1.minCoeff()) / (expect_nbr_ls1 + 1));
        nbr_ls0 = (ls0.maxCoeff() - ls0.minCoeff()) / pace + 1;
        nbr_ls1 = (ls1.maxCoeff() - ls1.minCoeff()) / pace + 1;
        get_level_set_sample_values_even_pace(ls0, nbr_ls0, pace, lsv0);
        get_level_set_sample_values_even_pace(ls1, nbr_ls1, pace, lsv1);
    }

    ivs.resize(nbr_ls0);
    verlist.reserve(nbr_ls0 * nbr_ls1);
    gridmat = Eigen::MatrixXi::Ones(nbr_ls0, nbr_ls1) * -1; // initially there is no quad patterns
    lsv0.resize(nbr_ls0);
    lsv1.resize(nbr_ls1);
    poly0_e0.resize(nbr_ls0);
    poly0_e1.resize(nbr_ls0);
    poly1_e0.resize(nbr_ls1);
    poly1_e1.resize(nbr_ls1);
    fid0.resize(nbr_ls0);
    fid1.resize(nbr_ls1);
    for (int i = 0; i < nbr_ls0; i++)
    {
        ivs[i].resize(nbr_ls1);
        poly0_e0[i].reserve(nbr_ls1);
        poly0_e1[i].reserve(nbr_ls1);
        fid0[i].reserve(nbr_ls1);
    }
    for (int i = 0; i < nbr_ls1; i++)
    {
        poly1_e0[i].reserve(nbr_ls0);
        poly1_e1[i].reserve(nbr_ls0);
        fid1[i].reserve(nbr_ls0);
    }
    
    // std::cout<<"sp_0 \n"<<lsv0.transpose()<<"\nsp_1\n"<<lsv1.transpose()<<std::endl;
    for (int i = 0; i < nbr_ls0; i++)
    {
        double vl = lsv0[i];
        get_iso_lines(V, F, ls0, vl, poly0_e0[i], poly0_e1[i], fid0[i]);
        // std::cout<<i<<" size "<<poly0_e0[i].size()<<"\n";
    }
    for (int i = 0; i < nbr_ls1; i++)
    {
        double vl = lsv1[i];
        get_iso_lines(V, F, ls1, vl, poly1_e0[i], poly1_e1[i], fid1[i]);
        // std::cout<<i<<" size "<<poly1_e0[i].size()<<"\n";
    }
    int vnbr = 0;
    for (int i = 0; i < nbr_ls0; i++)
    {
        for (int j = 0; j < nbr_ls1; j++)
        {
            Eigen::Vector3d ipoint;
            bool intersect = get_polyline_intersection(V, F, poly0_e0[i], poly0_e1[i], poly1_e0[j], poly1_e1[j], fid0[i], fid1[j], ipoint);
            if (intersect)
            {
                verlist.push_back(ipoint);
                gridmat(i, j) = vnbr;
                vnbr++;
            }
        }
    }
    std::cout << "extracted ver nbr " << verlist.size() << std::endl;
    vers = vec_list_to_matrix(verlist);
    
    if (even_pace)
    {
        std::cout << "The even pace quad, pace, " << std::min((ls0.maxCoeff() - ls0.minCoeff()) / (expect_nbr_ls0 + 1), (ls1.maxCoeff() - ls1.minCoeff()) / (expect_nbr_ls1 + 1))
                  << std::endl;
    }
    std::vector<Eigen::VectorXi> vrl;
    std::vector<Eigen::VectorXi> vrr;
    std::vector<Eigen::VectorXi> vcl;
    std::vector<Eigen::VectorXi> vcr;
    extract_web_from_index_mat(gridmat, Faces, threadshold_nbr, vrl, vrr, vcl, vcr);
    std::cout<<"Saving the info for each row or col of quads"<<std::endl;
    std::string fname = igl::file_dialog_save();
    save_quad_left_right_info(fname, vrl, vrr, vcl, vcr);
    std::cout<<"Each row or col of quads got saved"<<std::endl;
}
double polyline_length(const std::vector<Eigen::Vector3d>& line){
    int size = line.size();
    double total=0;
    for(int i=0;i<size-1;i++){
        Eigen::Vector3d dire=line[i]-line[i+1];
        double dis=dire.norm();
        total+=dis;
    }
    return total;
}
int select_longest_polyline(const std::vector<std::vector<Eigen::Vector3d>>& lines, double & length){
    int size = lines.size();
    assert(size>0);
    double max_length=0;
    int max_id=0;
    for(int i = 0;i<size; i++){
        double length = polyline_length(lines[i]);
        if(length>max_length){
            max_length=length;
            max_id=i;
        }
    }
    length=max_length;
    return max_id;
}
double polyline_total_length(const std::vector<std::vector<Eigen::Vector3d>> &lines)
{
    double total = 0;
    for (int i = 0; i < lines.size(); i++)
    {
        total += polyline_length(lines[i]);
    }
    return total;
}
void find_next_pt_on_polyline(const int start_seg, const std::vector<Eigen::Vector3d> &polyline, const double length,
                              const Eigen::Vector3d &pstart, int &seg, Eigen::Vector3d &pt)
{
    int nbr = polyline.size();
    double ocu_dis = 0;
    double dis_to_start = length + (polyline[start_seg] - pstart).norm();
    for (int i = start_seg; i < nbr - 1; i++)
    {
        double dis = (polyline[i] - polyline[i + 1]).norm();
        ocu_dis += dis;
        if (ocu_dis > dis_to_start)
        {                                         // the point should between i and i+1
            double diff = ocu_dis - dis_to_start; // the distance the point to i+1
            double t = diff / dis;
            assert(t >= 0 && t <= 1);
            Eigen::Vector3d p = get_3d_ver_from_t(t, polyline[i + 1], polyline[i]);
            pt = p;
            seg = i;
            return;
        }
    }
    std::cout<<"ERROR OUT OF SEGMENT"<<std::endl;
}
// nbr - 1 is the nbr of segments
// length is the length of the polyline
void sample_polyline_and_extend_verlist(const std::vector<Eigen::Vector3d>& polyline, const int nbr, const double length, std::vector<Eigen::Vector3d>& verlist){
    assert(nbr >= 3);
    double avg = length / (nbr - 1);
    verlist.push_back(polyline[0]);
    int start=0;
    for (int i = 0; i < nbr - 2; i++)
    {
        int seg;
        Eigen::Vector3d pt;
        find_next_pt_on_polyline(start, polyline, avg, verlist.back(), seg, pt);
        verlist.push_back(pt);
        start=seg;
    }
    verlist.push_back(polyline.back());
}
std::vector<Eigen::Vector3d> sample_one_polyline_based_on_length(const std::vector<Eigen::Vector3d>& polyline, const double avg){
    std::vector<Eigen::Vector3d> pts;
    double length = polyline_length(polyline);

    if (length < 0.5 * length) // if the polyline is short, we discard this one
    {
        return pts;
    }
    int nbr = length / avg + 2; // nbr of vers. make sure each segment is no longer than avg
    if (nbr == 2)
    {
        pts.push_back(polyline[0]);
        pts.push_back((polyline[0] + polyline.back()) / 2);
        pts.push_back(polyline.back());
        return pts;
    }
    sample_polyline_and_extend_verlist(polyline, nbr, length, pts);
    return pts;
}
// void sample_polyline_based_on_length(const std::vector<Eigen::Vector3d>& polyline, const double length)
void sample_polylines(const std::vector<std::vector<Eigen::Vector3d>> &polylines, const int nbr, const double length_total,
                      std::vector<std::vector<Eigen::Vector3d>> &verlists, double& avg)
{
    int lnbr = polylines.size();
    avg = length_total / (nbr - 1);// expected average length;
    verlists.resize(lnbr);
    for(int i=0;i<lnbr;i++){
        double length = polyline_length(polylines[i]);
        verlists[i] = sample_one_polyline_based_on_length(polylines[i], avg);
    }
}
void extract_Origami_web(const CGMesh &lsmesh, const Eigen::MatrixXd &V,const std::vector<CGMesh::HalfedgeHandle>& loop,
                         const Eigen::MatrixXi &F, const Eigen::VectorXd &ls,
                         const int expect_nbr_ls, const int expect_nbr_dis,
                         Eigen::MatrixXd &vers, Eigen::MatrixXi &Faces)
{
    Eigen::VectorXd lsv0;
    int nbr_ls0;

    nbr_ls0 = expect_nbr_ls;
    get_level_set_sample_values(ls, nbr_ls0, lsv0);
    std::vector<Eigen::Vector3d> verlist;
    verlist.reserve(nbr_ls0*expect_nbr_dis);
    for(int i=0;i<nbr_ls0;i++){
        std::vector<std::vector<Eigen::Vector3d>> polylines;
        double value = lsv0[i];
        std::vector<bool> left_large;
        get_iso_lines(lsmesh, loop, V, F, ls, value, polylines, left_large);

        double length ;
        int longest = select_longest_polyline(polylines, length);
        std::cout<<"longest got, size "<<polylines[longest].size()<<" length "<<length<<std::endl;
        // for(int j=0;j<polylines[longest].size();j++){
        //     std::cout<<"v "<<polylines[longest][j].transpose()<<std::endl;
        // }
        // exit(0);
        sample_polyline_and_extend_verlist(polylines[longest], expect_nbr_dis, length, verlist);
        std::cout<<"vers found"<<std::endl;
    }
    vers=vec_list_to_matrix(verlist);
    extract_web_mxn(nbr_ls0, expect_nbr_dis, Faces);

}

void lsTools::show_binormals(const Eigen::VectorXd &func, Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, Eigen::MatrixXd& binormals,  double ratio){
    
    E0 = V - Binormals * ratio;
    E1 = V + Binormals * ratio;
}
#include <igl/point_mesh_squared_distance.h>
#include <fstream>
bool write_quad_mesh_with_binormal(const std::string &fname, const Eigen::MatrixXd &Vt, const Eigen::MatrixXi &Ft, const Eigen::MatrixXd &bi,
                                   const Eigen::MatrixXd &Vq, const Eigen::MatrixXi &Fq)
{
    Eigen::MatrixXd C;
    Eigen::VectorXi I;
    Eigen::VectorXd D;
    igl::point_mesh_squared_distance(Vq, Vt, Ft, D, I, C);
    int nq = Vq.rows();
    Eigen::MatrixXd B;
    B.resize(nq, 3);
    for (int i = 0; i < nq; i++)
    {
        int fid = I(i);
        int v0 = Ft(fid, 0);
        int v1 = Ft(fid, 1);
        int v2 = Ft(fid, 2);
        Eigen::Vector3d dir = bi.row(v0) + bi.row(v1) + bi.row(v2);
        dir = dir.normalized();
        B.row(i) = dir;
    }
    std::ofstream file;
    file.open(fname);
    for (int i = 0; i < nq; i++)
    {
        file << "v " << Vq(i, 0) << " " << Vq(i, 1) << " " << Vq(i, 2) << " " << B(i, 0) << " " << B(i, 1) << " " << B(i, 2) << std::endl;
    }
    for (int i = 0; i < Fq.rows(); i++)
    {
        file << "f " << Fq(i, 0) + 1 << " " << Fq(i, 1) + 1 << " " << Fq(i, 2) + 1 << " " << Fq(i, 3) + 1 << std::endl;
    }
    file.close();
    return true;
}
int max_element_in_vec(const Eigen::VectorXd &vec)
{
    int size = vec.size();
    int mid = 0;
    double value = vec[0];
    for (int i = 0; i < size; i++)
    {
        if (vec[i] > value)
        {
            mid = i;
            value = vec[i];
        }
    }
    return mid;
}
// assign ids to the polylines so that we can get the quads easily
// start is the start ver id of this polyline
// last - 1 is the last id. meaning last is the first of the next lines
void assign_verid_to_sampled_polylines(const std::vector<std::vector<Eigen::Vector3d>>& sampled, const int start, 
std::vector<std::vector<int>>& ids, int& last){
    ids.resize(sampled.size());
    last=start;
    for(int i=0;i<sampled.size();i++){
        int lsize = sampled[i].size();
        ids[i].resize(lsize);
        for(int j=0;j<lsize;j++){
            ids[i][j] = last;
            last++;
        }
    }
}
std::vector<Eigen::Vector3d> invert_one_polyline(const std::vector<Eigen::Vector3d>& line){
    int size = line.size();
    std::vector<Eigen::Vector3d> result(size);
    for(int i=0;i<size;i++){
        result[size - i - 1] = line[i];
    }
    return result;
}
std::vector<std::vector<Eigen::Vector3d>> invert_polyline_set(const std::vector<std::vector<Eigen::Vector3d>>& pl){
    std::vector<std::vector<Eigen::Vector3d>> result(pl.size());
    for(int i =0;i<pl.size();i++){
        result[i] = invert_one_polyline(pl[pl.size() - i - 1]);
    }
    return result;
}
// after means the found one is after the given one
// checked shows if this segment is already included
int find_the_closest_line(const std::vector<std::vector<Eigen::Vector3d>> &pls, std::vector<bool> &checked, Eigen::Vector3d &v0,
                          Eigen::Vector3d &v1, int& which, bool& after)
{
    Eigen::Vector3d tv0, tv1;
    bool found = false;
    int id = 0;
    double dmin = std::numeric_limits<double>::max();
    which = 0; // 0: 0-0. 1: 0-1. 2: 1-0. 3: 1-1
    for (int i = 0; i < pls.size(); i++)
    {
        if (checked[i])
        {
            continue;
        }
        if (pls[i].size() < 2)
        {
            checked[i] = true;
            continue;
        }
        found = true;
        tv0 = pls[i].front();
        tv1 = pls[i].back();
        if ((v0 - tv0).norm() < dmin)
        {
            dmin = (v0 - tv0).norm();
            id = i;
            which = 0;
        }
        if ((v0 - tv1).norm() < dmin)
        {
            dmin = (v0 - tv1).norm();
            id = i;
            which = 1;
        }
        if ((v1 - tv0).norm() < dmin)
        {
            dmin = (v1 - tv0).norm();
            id = i;
            which = 2;
        }
        if ((v1 - tv1).norm() < dmin)
        {
            dmin = (v1 - tv1).norm();
            id = i;
            which = 3;
        }
    }
    if(!found){
        return -1;
    }
    checked[id]=true;
    if(which==0){// invert
        after = false;
        v0=pls[id].back();
        v1=v1;
    }
    if(which==1){
        after = false;
        v0=pls[id].front();
        v1=v1;
    }
    if(which==2){
        after = true;
        v0=v0;
        v1=pls[id].back();
    }
    if(which==3){// invert
        after = true;
        v0=v0;
        v1=pls[id].front();
    }
    return id;
    
}
std::vector<std::vector<Eigen::Vector3d>> extend_sorted_list(const std::vector<std::vector<Eigen::Vector3d>> &current, const bool invert,
                                                             const bool after, const std::vector<Eigen::Vector3d>& candidate)
{
    std::vector<std::vector<Eigen::Vector3d>> result;
    std::vector<Eigen::Vector3d> candi;
    if(invert){
        candi=invert_one_polyline(candidate);
    }
    else{
        candi=candidate;
    }
    if(after){
        result=current;
        result.push_back(candi);
    }
    else{
        result.push_back(candi);
        for(int i=0;i<current.size();i++){
            result.push_back(current[i]);
        }
    }
    return result;
}
void get_diff_polylines_order(const std::vector<std::vector<Eigen::Vector3d>> &pls, std::vector<std::vector<Eigen::Vector3d>> &sorted,
                              const int threadshold)
{
    sorted.clear();
    int nbr0 = pls.size();
    std::vector<bool> checked(nbr0, false);
    Eigen::Vector3d v0, v1;
    for (int i = 0; i < nbr0; i++)
    { // get the first
        if (pls[i].size() < 2)
        {
            checked[i] = true;
        }
        else
        {
            v0 = pls[i].front();
            v1 = pls[i].back();
            checked[i] = true;
            sorted.push_back(pls[i]);
            break;
        }
    }
    if (sorted.size() == 0)
    {
        std::cout << "please input the correct polylines " << std::endl;
        return;
    }
    double closet_dis = std::numeric_limits<double>::max();
    while (1)
    {
        int which;
        bool after;
        int find = find_the_closest_line(pls, checked, v0, v1, which, after);
        if (find == -1)
        {
            break;
        }
        bool invert = which == 0 || which == 3 ? true : false;
        if (pls[find].size() >= threadshold)
        {
            sorted = extend_sorted_list(sorted, invert, after, pls[find]);
        }
    }
}
std::vector<std::vector<Eigen::Vector3d>> invert_the_whole_polylines(const std::vector<std::vector<Eigen::Vector3d>> &ply)
{
    std::vector<std::vector<Eigen::Vector3d>> result(ply.size());
    for (int i = 0; i < ply.size(); i++)
    {
        result[i] = invert_one_polyline(ply[ply.size() - 1 - i]);
    }
    return result;
}
void sample_sort_and_discard_short(const std::vector<std::vector<Eigen::Vector3d>> &raw, const Eigen::Vector3d &direction,
                                   const double avg, const int threadshold, std::vector<std::vector<Eigen::Vector3d>> &result)
{
    std::vector<std::vector<Eigen::Vector3d>> tmp, sorted;
    for (int i = 0; i < raw.size(); i++)
    {
        tmp.push_back(sample_one_polyline_based_on_length(raw[i], avg)); // get the sampled points
    }

    get_diff_polylines_order(tmp, sorted, threadshold); // get them in order and delete short segments

    Eigen::Vector3d d = sorted.back().back() - sorted.front().front(); // the direction from the first point to the last
    if (d.dot(direction) < 0)
    {
        // need to invert all the lines
        result = invert_the_whole_polylines(sorted);
    }
    else
    {
        result = sorted;
    }
}

// // resort the polylines so that they are in the right order
// std::vector<std::vector<Eigen::Vector3d>> resort_polylines(const std::vector<std::vector<Eigen::Vector3d>> &polylines)
// {
//     std::vector<std::vector<Eigen::Vector3d>> tmp;
//     tmp.resize(polylines.size());
//     Eigen::Vector3d ref_dir;
//     bool found = false;
//     for (int i = 0; i < polylines.size(); i++)
//     {
//         if (polylines[i].size() > 2)
//         {
//             ref_dir = polylines[i].back() - polylines[i].front();
//             found = true;
//             break;
//         }
//     }
//     if (found == false)
//     {
//         std::cout << "please use a correct polyline for resorting!!!" << std::endl;
//         return tmp;
//     }
//     for (int i = 0; i < polylines.size(); i++)
//     {
//         if (polylines[i].size() < 2)
//         {
//             continue;
//         }
//         Eigen::Vector3d direction = polylines[i].back() - polylines[i].front();
//         if (direction.dot(ref_dir) < 0)
//         { // need invert
//             tmp[i] = invert_one_polyline(polylines[i]);
//         }
//         else
//         {xx
//             tmp[i] = polylines[i];
//         }
//     }
// }
bool project_point_on_lines(const Eigen::Vector3d &p, const std::vector<std::vector<Eigen::Vector3d>> &second,
                                 bool checkused, const std::vector<std::vector<bool>> &used,
                                 int &indi0, int &indi1)
{
    double dismin = std::numeric_limits<double>::max();
    int i0 = -1, i1 = -1;
    for (int i = 0; i < second.size(); i++)
    {
        for (int j = 0; j < second[i].size(); j++)
        {
            if (checkused)
            {
                if (used[i][j])// we need to skip all the points before the last used one
                {
                    dismin = std::numeric_limits<double>::max();
                    i0 = -1;
                    i1 = -1;
                    continue;
                }
            }

            double distance = (p - second[i][j]).norm();
            if (distance < dismin)
            {
                dismin = distance;
                i0 = i;
                i1 = j;
            }
        }
    }
    if(i0==-1||i1==-1){
        return false;
    }
    else{
        indi0=i0;
        indi1=i1;
        return true;
    }

}
// i0,i1 is the first point. j0, j1 is the second point
void connect_quads(const std::vector<std::vector<Eigen::Vector3d>> &first, const std::vector<std::vector<int>> &id0,
                   const std::vector<std::vector<Eigen::Vector3d>> &second, const std::vector<std::vector<int>> &id1,
                   const int i0, const int i1, const int j0, const int j1,
                   std::vector<Eigen::Vector4i> &quads, std::vector<std::vector<bool>> checked, int& final_i0, int& final_i1)
{
    int remain = std::min(first[i0].size() - i1, second[j0].size() - j1) - 1; // the number of quads
    for (int i = 0; i < remain; i++)
    {
        int v0 = id0[i0][i1 + i];
        int v1 = id0[i0][i1 + i + 1];
        int v2 = id1[j0][j1 + i + 1];
        int v3 = id1[j0][j1 + i];
        checked[j0][j1 + i] = true;
        checked[j0][j1 + i + 1] = true;
        final_i0=i0;
        final_i1=i1 + i + 1;
        quads.push_back(Eigen::Vector4i(v0, v1, v2, v3));
    }
}

// always extend the first polylines as the quads
void project_polylines_and_extend_quads(const std::vector<std::vector<Eigen::Vector3d>> &first, const std::vector<std::vector<int>> &id0,
                                        const std::vector<std::vector<Eigen::Vector3d>> &second, const std::vector<std::vector<int>> &id1, 
                                        std::vector<Eigen::Vector4i> &quads)
{
    int size1 = first.size();
    int size2 = second.size();
    std::vector<std::vector<bool>> used; // if the points in the second line is used for connecting quads, they cannot be found by the new points in the first line
    // init the used as false
    used.resize(size2);
    for (int i = 0; i < size2; i++)
    {
        std::vector<bool> tbool(second[i].size(), false);
        used[i] = tbool;
    }

    // start projection and connect quads
    for(int i=0;i<size1;i++){
        bool need_break = false;
        for (int j = 0; j < first[i].size() - 1; j++)
        {
            int rid0, rid1;
            bool found = project_point_on_lines(first[i][j],second,true,used,rid0,rid1);
            if(!found){// if there are no points left, break the loop
                need_break=true;
                break;
            }
            if(rid1==second[rid0].size()-1){// if the point is projected to the end of second line, we project the next point
                
                continue;
            }
            if(rid1==0){// if it is projected to the first point, we project back 
                // project again
                int tid0, tid1;//
                found = project_point_on_lines(second[rid0][rid1], first, false, used, tid0, tid1);
                if (tid1 == first[tid0].size() - 1)
                { // if project back to the end of some seg, means no quad can be extracted, check the next point
                    continue;
                }
                if(tid0<i){// if projected back to the previous seg, continue
                    continue;
                }
                if (tid0 == i && tid1 < j)
                { // if projected to the previous seg, continue;
                    continue;
                }
                int final_i0, final_i1;
                connect_quads(first, id0, second, id1, tid0, tid1, rid0, rid1, quads, used, final_i0, final_i1);
                i=final_i0;
                j=final_i1;
                // the point is projected back somewhere, 
            }
            else{
                int final_i0, final_i1;
                // connect quads
                connect_quads(first, id0, second, id1, i, j, rid0, rid1, quads, used, final_i0, final_i1);
                i=final_i0;
                j=final_i1;
            }
        }
        if(need_break){
            break;
        }
    }
}

Eigen::MatrixXd extract_vers(const std::vector<std::vector<std::vector<Eigen::Vector3d>>> &vertices)
{
    Eigen::MatrixXd result;
    int counter = 0;
    for (int i = 0; i < vertices.size(); i++)
    {
        for (int j = 0; j < vertices[i].size(); j++)
        {
            for (int k = 0; k < vertices[i][j].size(); k++)
            {
                counter++;
            }
        }
    }
    result.resize(counter,3);
    counter = 0;
    for (int i = 0; i < vertices.size(); i++)
    {
        for (int j = 0; j < vertices[i].size(); j++)
        {
            for (int k = 0; k < vertices[i][j].size(); k++)
            {
                result.row(counter)=vertices[i][j][k];
                counter++;
            }
        }
    }
    return result;
}
Eigen::MatrixXi extract_quads(const std::vector<Eigen::Vector4i> &quads){
    int nbr=quads.size();
    Eigen::MatrixXi result(nbr, 4);
    for(int i=0;i<nbr;i++){
        result.row(i)=quads[i];
    }
    return result;
}
// preserves the boundary with given accuracy decided by the size of the mesh
// the threadshold is the least number of vers in each seg
void extract_Quad_Mesh_Zigzag(const CGMesh &lsmesh,const std::vector<CGMesh::HalfedgeHandle>& loop, const Eigen::MatrixXd &V,
                         const Eigen::MatrixXi &F, const Eigen::VectorXd &ls,
                         const int expect_nbr_ls, const int expect_nbr_dis, const int threadshold_nbr,
                         Eigen::MatrixXd &vers, Eigen::MatrixXi &Faces)
{
    Eigen::VectorXd lsv0;
    int nbr_ls0;
    std::vector<Eigen::Vector4i> quads;
    nbr_ls0 = expect_nbr_ls;
    get_level_set_sample_values(ls, nbr_ls0, lsv0);
    std::vector<Eigen::Vector3d> verlist;
    verlist.reserve(nbr_ls0*expect_nbr_dis);
    std::vector<std::vector<std::vector<Eigen::Vector3d>>> curves_all;
    Eigen::VectorXd clengths(nbr_ls0);
    for(int i=0;i<nbr_ls0;i++){
        std::vector<std::vector<Eigen::Vector3d>> polylines;
        double value = lsv0[i];
        std::vector<bool> left_large;
        std::cout<<"iso..."<<std::endl;
        get_iso_lines(lsmesh, loop, V, F, ls, value, polylines, left_large);
        // std::cout<<"iso got, size "<<polylines.size()<<std::endl;
        // double length ;
        // int longest = select_longest_polyline(polylines, length);
        double total_length = polyline_total_length(polylines);
        curves_all.push_back(polylines);
        clengths[i] = (total_length);
        // std::cout<<"longest got, size "<<polylines[longest].size()<<" length "<<length<<std::endl;
        // // for(int j=0;j<polylines[longest].size();j++){
        // //     std::cout<<"v "<<polylines[longest][j].transpose()<<std::endl;
        // // }
        // // exit(0);
        // sample_polyline_and_extend_verlist(polylines[longest], expect_nbr_dis, length, verlist);
        // std::cout<<"vers found"<<std::endl;
    }
    // vers=vec_list_to_matrix(verlist);
    // extract_web_mxn(nbr_ls0, expect_nbr_dis, Faces);
    int longest = max_element_in_vec(clengths); // get the longest iso-lines
    std::vector<std::vector<Eigen::Vector3d>> first_lines, sorted;
    std::vector<std::vector<int>> first_ids;
    double expect_length;
    int first = 0;
    int last;

    sample_polylines(curves_all[longest], expect_nbr_dis, clengths[longest], first_lines, expect_length); // first sample the longest line
    get_diff_polylines_order(first_lines, sorted, threadshold_nbr);
    assign_verid_to_sampled_polylines(sorted, first, first_ids, last);
    first_lines=sorted;


    std::vector<std::vector<std::vector<Eigen::Vector3d>>> vertices;
    std::vector<std::vector<std::vector<int>>> ids;
    std::vector<Eigen::Vector4i> faces;
    vertices.push_back(first_lines);
    ids.push_back(first_ids);
    std::cout<<"before extracting quads: nbr of lines, "<<curves_all.size()<<std::endl;
    std::cout<<"longetst, "<<longest<<std::endl;
    std::vector<std::vector<Eigen::Vector3d>> this_lines = first_lines; // the first line
    std::vector<std::vector<int>> this_ids = first_ids;                 // the first ids
    for (int i = longest + 1; i < curves_all.size(); i++)// to the right
    {
        std::vector<std::vector<Eigen::Vector3d>> next_raw = curves_all[i];                // the next polylines
        Eigen::Vector3d direction = this_lines.back().back() - this_lines.front().front(); // the reference orientation
        std::vector<std::vector<Eigen::Vector3d>> next_lines;                              // the sampled and sorted lines
        std::vector<std::vector<int>> next_ids;
        first = last; // the ver ids

        sample_sort_and_discard_short(next_raw, direction, expect_length, threadshold_nbr, next_lines);
        assign_verid_to_sampled_polylines(next_lines, first, next_ids, last);
        project_polylines_and_extend_quads(this_lines, this_ids, next_lines, next_ids, quads);
        this_lines=next_lines;
        this_ids=next_ids;
        vertices.push_back(next_lines);
        // std::cout<<"first ,segs "<<first_lines.size()<<", second segs, "<<next_lines.size()<<std::endl;
        // std::cout<<"1 size, "<<first_lines[0].size()<<", 2 size, "<<next_lines[0].size()<<std::endl;
        // break;
        
        // std::cout<<"forward, "<<longest<<std::endl;
    }
    this_lines = first_lines; // the first line
    this_ids = first_ids;                 // the first ids
    for (int i = longest - 1; i >= 0; i--) // to the right
    {
        std::vector<std::vector<Eigen::Vector3d>> next_raw = curves_all[i];                // the next polylines
        Eigen::Vector3d direction = this_lines.back().back() - this_lines.front().front(); // the reference orientation
        std::vector<std::vector<Eigen::Vector3d>> next_lines;                              // the sampled and sorted lines
        std::vector<std::vector<int>> next_ids;
        first = last; // the ver ids

        sample_sort_and_discard_short(next_raw, direction, expect_length, threadshold_nbr, next_lines);
        assign_verid_to_sampled_polylines(next_lines, first, next_ids, last);
        project_polylines_and_extend_quads(next_lines, next_ids, this_lines, this_ids, quads);
        this_lines=next_lines;
        this_ids=next_ids;
        vertices.push_back(next_lines);
    }
    vers = extract_vers(vertices);
    Faces=extract_quads(quads);
}

Eigen::MatrixXd get_each_face_direction(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXd &func)
{
    int vnbr = V.rows();
    int fnbr = F.rows();
    Eigen::MatrixXd result(fnbr, 3);
    for (int i = 0; i < fnbr; i++)
    {
        int v0 = F(i, 0);
        int v1 = F(i, 1);
        int v2 = F(i, 2);
        // dir0 * V0 + dir1 * V1 + dir2 * V2 is the iso-line direction
        Eigen::Vector3d dir0 = V.row(v2) - V.row(v1);
        Eigen::Vector3d dir1 = V.row(v0) - V.row(v2);
        Eigen::Vector3d dir2 = V.row(v1) - V.row(v0);
        Eigen::Vector3d iso = dir0 * func[v0] + dir1 * func[v1] + dir2 * func[v2];
        result.row(i) = iso;
    }
    return result;
}