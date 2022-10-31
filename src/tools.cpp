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
spMat rib_method_arrange_matrices_rows(const std::vector<spMat> &mats)
{
    int nvec = mats.size();
    int inner_rows = mats[0].rows();
    int inner_cols = mats[0].cols();
    spMat result;
    result.resize(nvec * inner_rows, inner_cols); // the rows get more, the column number remain the same.
    std::vector<Trip> triplets;
    triplets.reserve(inner_cols * inner_rows * nvec);
    for (int i = 0; i < nvec; i++)
    {
        // TODO
    }
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
bool save_levelset(const Eigen::VectorXd &ls){
    std::string fname = igl::file_dialog_save();
    std::ofstream file;
    file.open(fname);
    for(int i=0;i<ls.size();i++){
        file<<ls[i]<<std::endl;
    }
    file.close();
    return true;
}
bool read_levelset(Eigen::VectorXd &ls){
    std::string fname = igl::file_dialog_open();

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

    return true;
}

// get vid1, vid2, vidA and vidB.
// here we do not assume fvalues[vid1]>fvalues[vid2], we need to consider it outside
std::array<int,4> get_vers_around_edge(CGMesh& lsmesh, int edgeid, int& fid1, int &fid2){
    CGMesh::EdgeHandle edge_middle=lsmesh.edge_handle(edgeid);
    CGMesh::HalfedgeHandle he = lsmesh.halfedge_handle(edge_middle, 0);
    int vid1 = lsmesh.to_vertex_handle(he).idx();
    int vid2 = lsmesh.from_vertex_handle(he).idx();
    // if(fvalues[vid1]<fvalues[vid2]){
    //     he = lsmesh.halfedge_handle(edge_middle, 1);
    //     vid1 = lsmesh.to_vertex_handle(he).idx();
    //     vid2 = lsmesh.from_vertex_handle(he).idx();
    // }
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
void lsTools::show_binormals(Eigen::MatrixXd& E0, Eigen::MatrixXd& E1, double ratio){
    std::cout<<"show binormals "<<std::endl;
    int ninner=IVids.size();
    if(ninner<=0){
        std::cout<<"please optimize pseudo-geodesic energy first"<<std::endl;
        return;
    }
    std::cout<<"there are "<<ninner<<" bi-normals"<<std::endl;
    E0.resize(ninner,3);
    E1.resize(ninner,3);
    for (int i = 0; i < ninner; i++)
    {
        int vid = IVids[i];
        Eigen::Vector3d ver=V.row(vid);
        if(ActInner[i] == false){// this is a singularity
            E0.row(vid)=ver;
            E1.row(vid)=ver;
            continue;
        }
        E0.row(i)=ver;
        E1.row(i)=ver+ratio*Eigen::Vector3d(vBinormal.row(i));

    }
}