#include <lsc/basic.h>
#include <lsc/tools.h>
#include <igl/readOBJ.h>

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
// void lsTools::debug_tool(int id, int id2, double value)
// {
//     trace_vers.resize(1);
//     std::cout << "checking one pseudo geodesic" << std::endl;
//     // debug init
//     ver_dbg.resize(0, 0);
//     ver_dbg1.resize(0, 0);
//     flag_dbg = false;
//     E0_dbg.resize(0, 0);
//     direction_dbg.resize(0, 0);

//     pnorm_dbg.resize(0, 0);
//     pnorm_list_dbg.clear();
//     flag_dbg = true;
//     id_dbg = id2;

//     // debug init end
//     double start_point_para = 0.5;
//     double start_angle_degree = 60;
//     double target_angle = value;
//     std::vector<Eigen::Vector3d> curve;
//     std::vector<CGMesh::HalfedgeHandle> handles;
//     trace_single_pseudo_geodesic_curve_pseudo_vertex_method(target_angle, Boundary_Edges[id], start_point_para, start_angle_degree,
//                                                             curve, handles);
//     flag_dbg = false;
//     int nbr_midpts = curve.size() - 2;
//     E0_dbg.resize(nbr_midpts, 3);
//     direction_dbg.resize(nbr_midpts, 3);
//     Eigen::MatrixXd binormals_real(nbr_midpts, 3);
//     for (int i = 0; i < nbr_midpts; i++)
//     {
//         E0_dbg.row(i) = curve[i + 1];
//         Eigen::Vector3d vec0 = (pseudo_vers_dbg[i + 1] - pseudo_vers_dbg[i]).normalized();
//         Eigen::Vector3d vec1 = (curve[i + 2] - pseudo_vers_dbg[i + 1]).normalized();
//         direction_dbg.row(i) = vec0.cross(vec1).normalized(); // bi normal of the pseudo-curve
//         vec0 = (curve[i + 1] - curve[i]).normalized();
//         vec1 = (curve[i + 2] - curve[i + 1]).normalized();
//         binormals_real.row(i) = vec0.cross(vec1).normalized();
//     }
//     if (pnorm_list_dbg.size() != nbr_midpts)
//     {
//         std::cout << "dangerous size! " << pnorm_list_dbg.size() << ", " << nbr_midpts << std::endl;
//     }
//     if (target_angle > 90 - ANGLE_TOLERANCE && target_angle < 90 + ANGLE_TOLERANCE)
//     {
//         std::cout << "Traced geodesic" << std::endl;
//     }
//     else
//     {
//         assert(pnorm_list_dbg.size() == nbr_midpts);
//         pnorm_dbg = vec_list_to_matrix(pnorm_list_dbg);
//         // check the abs(consin)
//         for (int i = 0; i < nbr_midpts; i++)
//         {
//             Eigen::Vector3d tmp_dir1 = direction_dbg.row(i);
//             Eigen::Vector3d tmp_dir2 = pnorm_list_dbg[i].normalized();
//             double cosin = tmp_dir1.dot(tmp_dir2);
//             std::cout << i << "th cosin^2 " << cosin * cosin << std::endl;
//             tmp_dir1 = binormals_real.row(i);
//             cosin = tmp_dir1.dot(tmp_dir2);
//             std::cout << i << "th cosin^2 " << cosin * cosin << std::endl;
//         }
//     }

//     // TODO temporarily checking one trace line
//     trace_vers[0] = curve;
//     trace_hehs.push_back(handles);
// }
// void lsTools::debug_tool_v2(const std::vector<int> &ids, const std::vector<double> values)
// {
//     // cylinder_example(5, 10, 50, 30);
//     // exit(0);
//     int nbr_itv = ids[0]; // every nbr_itv boundary edges we shoot one curve
//     if (nbr_itv < 1)
//     {
//         std::cout << "Please set up the parameter nbr_itv " << std::endl;
//         return;
//     }
//     double target_angle = values[0];
//     double start_angel = values[1];
//     OpenMesh::HalfedgeHandle init_edge = Boundary_Edges[0];
//     if (lsmesh.face_handle(init_edge).idx() >= 0)
//     {
//         init_edge = lsmesh.opposite_halfedge_handle(init_edge);
//     }
//     OpenMesh::HalfedgeHandle checking_edge = init_edge;
//     int curve_id = 0;
//     while (1)
//     {
//         ver_dbg.resize(0, 0);
//         ver_dbg1.resize(0, 0);
//         flag_dbg = false;
//         E0_dbg.resize(0, 0);
//         direction_dbg.resize(0, 0);
//         pnorm_dbg.resize(0, 0);
//         pnorm_list_dbg.clear();
//         flag_dbg = true;
//         id_dbg = ids[1];
//         std::vector<Eigen::Vector3d> curve;
//         std::vector<CGMesh::HalfedgeHandle> handles;
//         trace_single_pseudo_geodesic_curve_pseudo_vertex_method(target_angle, checking_edge, 0.5, start_angel,
//                                                                 curve, handles);
//         flag_dbg = false;
//         curve_id++;
//         if (curve_id == -1)
//         {
//             bool stop_flag = false;

//             for (int i = 0; i < nbr_itv; i++)
//             {
//                 checking_edge = lsmesh.next_halfedge_handle(checking_edge);
//                 if (checking_edge == init_edge)
//                 {
//                     stop_flag = true;
//                     break;
//                 }
//             }
//             continue;
//         }
//         trace_vers.push_back(curve);
//         trace_hehs.push_back(handles);
//         bool stop_flag = false;

//         if (curve_id == -1)
//         {
//             break;
//         }
//         for (int i = 0; i < nbr_itv; i++)
//         {
//             checking_edge = lsmesh.next_halfedge_handle(checking_edge);
//             if (checking_edge == init_edge)
//             {
//                 stop_flag = true;
//                 break;
//             }
//         }
//         if (stop_flag)
//         {
//             break;
//         }
//     }
// }

// void lsTools::debug_tool_v3(int id, int id2, double value)
// {
//     trace_vers.resize(1);
//     std::cout << "checking one pseudo geodesic" << std::endl;
//     // debug init
//     ver_dbg.resize(0, 0);
//     ver_dbg1.resize(0, 0);
//     flag_dbg = false;
//     E0_dbg.resize(0, 0);
//     direction_dbg.resize(0, 0);
//     trace_hehs.clear();
//     pnorm_dbg.resize(0, 0);
//     pnorm_list_dbg.clear();
//     flag_dbg = true;
//     id_dbg = id2;

//     // debug init end
//     double start_point_para = 0.5;
//     double start_angle_degree = 60;
//     double target_angle = value;
//     std::vector<Eigen::Vector3d> curve;
//     std::vector<CGMesh::HalfedgeHandle> handles;
//     trace_single_pseudo_geodesic_curve(target_angle, Boundary_Edges[id], start_point_para, start_angle_degree,
//                                        curve, handles);
//     flag_dbg = false;
//     int nbr_midpts = curve.size() - 2;
//     E0_dbg.resize(nbr_midpts, 3);
//     direction_dbg.resize(nbr_midpts, 3);
//     Eigen::MatrixXd binormals_real(nbr_midpts, 3);
//     for (int i = 0; i < nbr_midpts; i++)
//     {
//         E0_dbg.row(i) = curve[i + 1];
//         Eigen::Vector3d vec0 = (curve[i + 1] - curve[i]).normalized();
//         Eigen::Vector3d vec1 = (curve[i + 2] - curve[i + 1]).normalized();
//         direction_dbg.row(i) = vec0.cross(vec1).normalized(); // bi normal of the pseudo-curve
//     }

//     if (target_angle > 90 - ANGLE_TOLERANCE && target_angle < 90 + ANGLE_TOLERANCE)
//     {
//         std::cout << "Traced geodesic" << std::endl;
//     }
//     else
//     {

//         for (int i = 0; i < nbr_midpts; i++)
//         {
//             Eigen::Vector3d tmp_dir1 = direction_dbg.row(i);
//             Eigen::Vector3d tmp_dir2 = pnorm_list_dbg[i].normalized();
//             double cosin = tmp_dir1.dot(tmp_dir2);
//             std::cout << i << "th cosin^2 " << cosin * cosin << std::endl;
//         }
//     }

//     // TODO temporarily checking one trace line
//     trace_vers[0] = curve;
//     trace_hehs.push_back(handles);
// }
// void lsTools::debug_tool_v4(const std::vector<int> &ids, const std::vector<double> values)
// {
//     // cylinder_example(5, 10, 50, 30);
//     // exit(0);
//     trace_vers.clear();
//     trace_hehs.clear();
//     int nbr_itv = ids[0]; // every nbr_itv boundary edges we shoot one curve
//     if (nbr_itv < 1)
//     {
//         std::cout << "Please set up the parameter nbr_itv " << std::endl;
//         return;
//     }
//     double target_angle = values[0];
//     double start_angel = values[1];
//     OpenMesh::HalfedgeHandle init_edge = Boundary_Edges[0];
//     if (lsmesh.face_handle(init_edge).idx() >= 0)
//     {
//         init_edge = lsmesh.opposite_halfedge_handle(init_edge);
//     }
//     OpenMesh::HalfedgeHandle checking_edge = init_edge;
//     int curve_id = 0;
//     while (1)
//     {
//         ver_dbg.resize(0, 0);
//         ver_dbg1.resize(0, 0);
//         E0_dbg.resize(0, 0);
//         direction_dbg.resize(0, 0);
//         pnorm_dbg.resize(0, 0);
//         pnorm_list_dbg.clear();
//         flag_dbg = true;
//         id_dbg = ids[1];
//         std::vector<Eigen::Vector3d> curve;
//         std::vector<CGMesh::HalfedgeHandle> handles;
//         trace_single_pseudo_geodesic_curve(target_angle, checking_edge, 0.5, start_angel,
//                                            curve, handles);
//         flag_dbg = false;
//         curve_id++;
//         if (curve_id == -1)
//         {
//             bool stop_flag = false;

//             for (int i = 0; i < nbr_itv; i++)
//             {
//                 checking_edge = lsmesh.next_halfedge_handle(checking_edge);
//                 if (checking_edge == init_edge)
//                 {
//                     stop_flag = true;
//                     break;
//                 }
//             }
//             continue;
//         }
//         trace_vers.push_back(curve);
//         trace_hehs.push_back(handles);
//         bool stop_flag = false;

//         if (curve_id == -1)
//         {
//             break;
//         }
//         for (int i = 0; i < nbr_itv; i++)
//         {
//             checking_edge = lsmesh.next_halfedge_handle(checking_edge);
//             if (checking_edge == init_edge)
//             {
//                 stop_flag = true;
//                 break;
//             }
//         }
//         if (stop_flag)
//         {
//             break;
//         }
//     }
// }

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
    if (vec.size() > 0)
    {
        mat.resize(vec.size(), 3);
        for (int i = 0; i < mat.rows(); i++)
        {
            mat.row(i) = vec[i];
        }
    }
    else{
        std::cout<<"The vertex vector is empty"<<std::endl;
    }
    return mat;
}
Eigen::VectorXd vec_list_to_vector(const std::vector<double>& vec){
    Eigen::VectorXd result;
    result.resize(vec.size());
    for(int i=0;i<vec.size();i++){
        result[i] = vec[i];
    }
    return result;
}
Eigen::VectorXi vec_list_to_vector(const std::vector<int> &vec)
{
    Eigen::VectorXi result;
    result.resize(vec.size());
    for (int i = 0; i < vec.size(); i++)
    {
        result[i] = vec[i];
    }
    return result;
}
std::vector<Eigen::Vector3d> mat_to_vec_list(const Eigen::MatrixXd &m)
{
    std::vector<Eigen::Vector3d> result(m.rows());
    for (int i = 0; i < m.rows(); i++)
    {
        result[i] = m.row(i);
    }
    return result;
}

CGMesh::HalfedgeHandle boundary_halfedge(const CGMesh &lsmesh, const CGMesh::HalfedgeHandle &boundary_edge)
{
    assert(lsmesh.is_boundary(lsmesh.from_vertex_handle(boundary_edge)));
    assert(lsmesh.is_boundary(lsmesh.to_vertex_handle(boundary_edge)));
    CGMesh::HalfedgeHandle result = boundary_edge;
    if (lsmesh.face_handle(result).idx() >= 0)
    {
        result = lsmesh.opposite_halfedge_handle(result);
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

void select_mesh_boundary_curve_on_boundary_loop(CGMesh &lsmesh, const Eigen::MatrixXd &V, const int loopid,
                                                 const int start_edge, const int nbr_edges,
                                                 const std::vector<CGMesh::HalfedgeHandle> &Boundary_Edges,
                                                 std::vector<CGMesh::HalfedgeHandle> &selected)
{
    std::cout << "\nboundary edges size " << Boundary_Edges.size() << std::endl;
    selected.clear();
    int vnbr = V.rows();
    int bsize = Boundary_Edges.size();
    int enbr = lsmesh.n_edges();
    std::vector<std::vector<CGMesh::HalfedgeHandle>> boundaries;
    //////////////////////////////////////////
    SpVeci bedges; // 0: not a boundary edge; 1: not checked boundary edge;
    SpVeci recorded_edges;
    recorded_edges.resize(enbr);
    bedges.resize(enbr);
    for(int i=0;i<Boundary_Edges.size();i++){
        CGMesh::EdgeHandle eh=lsmesh.edge_handle(Boundary_Edges[i]);
        bedges.coeffRef(eh.idx())=1; // mark every boundary edges
    }
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
                int edge_id = lsmesh.edge_handle(next_he).idx();
                if (recorded_edges.coeffRef(edge_id) == 0) // if the next edge is not checked yet
                {
                    recorded_edges.coeffRef(edge_id) = 1;
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
    //////////////////////////////////////////

    // // first get the boundary loops
    // for (int i = 0; i < bsize; i++)
    // {
    //     if (recorded_edges(i) == 0)
    //     { // if this is a boundary edge, and not checked yet
    //         recorded_edges(i) = 1;
    //         CGMesh::HalfedgeHandle start_he = Boundary_Edges[i];
    //         start_he = boundary_halfedge(lsmesh, start_he);
    //         // CGMesh::HalfedgeHandle start_copy = start_he;
    //         std::vector<CGMesh::HalfedgeHandle> loop;
    //         loop.push_back(start_he);
    //         for (int j = 0; j < bsize; j++)
    //         {
    //             CGMesh::HalfedgeHandle next_he = lsmesh.next_halfedge_handle(start_he);
    //             if (recorded_edges(j) == 0) // if the next edge is not checked yet
    //             {
    //                 recorded_edges(j) = 1;
    //                 loop.push_back(next_he);
    //                 start_he = next_he;
    //             }
    //             else
    //             {
    //                 break;
    //             }
    //         }
    //         boundaries.push_back(loop);
    //     }
    // }
    std::cout<<"Total boundary loops nbr, "<<boundaries.size()<<std::endl;

    if (loopid < 0 || loopid >= boundaries.size())
    {
        std::cout << "Please select a correct boundary loop" << std::endl;
        return;
    }

    std::vector<CGMesh::HalfedgeHandle> loop = boundaries[loopid];
    int cp_nbr_edges = nbr_edges;
    if (cp_nbr_edges > loop.size())
    {
        cp_nbr_edges =loop.size();
        std::cout<<"Please reduce the selected number of edges"<<std::endl;
    }

    for (int i = start_edge;; i++)
    {
        selected.push_back(loop[i]);
        if (selected.size() == cp_nbr_edges)
        {
            break;
        }
        if (i + 1 == loop.size())// if the next iteration exceeds the limit, start over
        {
            i = -1;
        }
    }
    
    return;
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
// make a cylinder mesh
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

std::vector<std::vector<Eigen::Vector3d>> xyz_to_vec_list(const std::vector<std::vector<double>> &x, const std::vector<std::vector<double>> &y, const std::vector<std::vector<double>> &z)
{
    std::vector<std::vector<Eigen::Vector3d>> result;
    result.reserve(x.size());
    for (int i = 0; i < x.size(); i++)
    {
        std::vector<Eigen::Vector3d> line;
        for (int j = 0; j < x[i].size(); j++)
        {
            Eigen::Vector3d ver;
            ver[0] = x[i][j];
            ver[1] = y[i][j];
            ver[2] = z[i][j];
            line.push_back(ver);
        }
        result.push_back(line);
    }
    return result;
}
bool read_polylines(std::vector<std::vector<Eigen::Vector3d>>& ply){
    ply.clear();
    std::cout<<"Reading the binormal files and polyline files. Please type down the prefix"<<std::endl;
    std::string fname = igl::file_dialog_save();
    if (fname.length() == 0)
    {
        std::cout<<"Please type something "<<std::endl;
        return false;
    }
    std::vector<std::vector<double>> vx, vy, vz;
    read_csv_data_lbl(fname+"_x.csv", vx);
    read_csv_data_lbl(fname+"_y.csv", vy);
    read_csv_data_lbl(fname + "_z.csv", vz);
    ply = xyz_to_vec_list(vx, vy, vz);
    if (ply.size() > 0)
        return true;
    return false;
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

// read csv data line by line and save the data by rows
bool read_csv_data_lbl(const std::string fname, std::vector<std::vector<double>> &data){
    std::cout<<"reading "<<fname<<std::endl;
    if (fname.length() == 0)
        return false;

    std::ifstream infile;
    std::vector<std::vector<double>> results;

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
            std::vector<double> record;
            int c = 0;
            while (ss)
            {
                std::string line;
                if (!getline(ss, line, ','))
                    break;
                try
                {

                    record.push_back(std::stod(line));
                    c++;

                }
                catch (const std::invalid_argument e)
                {
                    std::cout << "NaN found in file " << fname
                              <<  std::endl;
                    e.what();
                }
            }

            
            results.push_back(record);
        }
    }
    data = results;
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



bool halfedge_has_value(const CGMesh &lsmesh, const Eigen::MatrixXd &V,
                         const Eigen::VectorXd &ls, const CGMesh::HalfedgeHandle &hd, const double value, Eigen::Vector3d &pt)
{
    int id_f=lsmesh.from_vertex_handle(hd).idx();
    int id_t=lsmesh.to_vertex_handle(hd).idx();
    double value_f=ls[id_f];
    double value_t=ls[id_t];
    if(value_f==value_t){ // if this edge is parallel to the level set, skip
        return false;
    }
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
bool halfedge_has_value(const CGMesh &lsmesh, const Eigen::MatrixXd &V,
                         const Eigen::VectorXd &ls, const CGMesh::HalfedgeHandle &hd, const double value, double &t, Eigen::Vector3d &pt)
{
    int id_f=lsmesh.from_vertex_handle(hd).idx();
    int id_t=lsmesh.to_vertex_handle(hd).idx();
    double value_f=ls[id_f];
    double value_t=ls[id_t];
    if(value_f==value_t){ // if this edge is parallel to the level set, skip
        return false;
    }
    t=get_t_of_value(value, value_f, value_t);
    if(t<0||t>1){
        return false;
    }
    Eigen::Vector3d ver_f=V.row(id_f);
    Eigen::Vector3d ver_t= V.row(id_t);
    pt=get_3d_ver_from_t(t, ver_f, ver_t);
    // std::cout<<"has value "<<value<<", vf "<<value_f<<" vt "<<value_t<<", idf "<<id_f<<", idt "<<id_t<<std::endl;
    return true;
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
// the mesh provides the connectivity, loop is the boundary loop
// left_large indicates if the from ver of each boundary edge is larger value
void get_iso_lines_baricenter_coord(const CGMesh &lsmesh, const std::vector<CGMesh::HalfedgeHandle>& loop, const Eigen::MatrixXd &V,
                   const Eigen::MatrixXi &F, const Eigen::VectorXd &ls, double value, std::vector<std::vector<double>> &paras,
                   std::vector<std::vector<CGMesh::HalfedgeHandle>> &handles,
                   std::vector<bool>& left_large)
{
    paras.clear();
    handles.clear();
    left_large.clear();
    std::vector<CGMesh::HalfedgeHandle> has_value; // these boundary edges has value
    std::vector<Eigen::Vector3d> start_pts; // these are the computed start/end points
    std::vector<double> start_ts;
    for (int i = 0; i < loop.size(); i++)
    { // find the points on boundary

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
        double t;
        bool found = halfedge_has_value(lsmesh, V, ls, bhd, value, t, intersection);
        if(!found){
            continue;
        }
        // find the point on bnd, record the halfedge, the intersection pt and the para.
        has_value.push_back(bhd);
        start_pts.push_back(intersection);
        start_ts.push_back(t);
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
        std::vector<double> tmp_ts;
        std::vector<CGMesh::HalfedgeHandle> tmp_hds;
        CGMesh::HalfedgeHandle thd=has_value[i]; // the start handle
        tmpts.push_back(start_pts[i]);
        tmp_ts.push_back(start_ts[i]);
        tmp_hds.push_back(thd);
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
                paras.push_back(tmp_ts);
                handles.push_back(tmp_hds);
                break;
            }
            // find the next pt
            int id_f=lsmesh.from_vertex_handle(ophd).idx();
            int id_t=lsmesh.to_vertex_handle(ophd).idx();
            CGMesh::HalfedgeHandle prev=lsmesh.prev_halfedge_handle(ophd);
            CGMesh::HalfedgeHandle next=lsmesh.next_halfedge_handle(ophd);
            int bid=lsmesh.from_vertex_handle(prev).idx();
            if (bid == id_f || bid == id_t)
            {
                std::cout << "wrong topology" << std::endl;
            }
            if (ls[bid] == value)
            { // a degenerate case where the value is on a vertex of the mesh
                value += 1e-16;
            }
            Eigen::Vector3d intersect;
            double t;
            bool found = halfedge_has_value(lsmesh, V, ls, prev, value, t, intersect);
            if(found){
                tmpts.push_back(intersect);
                tmp_ts.push_back(t);
                tmp_hds.push_back(prev);
                thd=prev;
                continue;
            }
            found = halfedge_has_value(lsmesh, V, ls, next, value, t, intersect);
            if(found){
                tmpts.push_back(intersect);
                tmp_ts.push_back(t);
                tmp_hds.push_back(next);
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

double point_seg_distance(const Eigen::Vector3d &pt, const Eigen::Vector3d &E0, const Eigen::Vector3d &E1)
{
    double d0 = (E0 - pt).dot(E1 - E0);
    double d1 = (E1 - E0).dot(E1 - E0);
    double t = -d0 / d1;
    if (t >= 0 && t <= 1)
    {
        Eigen::Vector3d proj = E0 + t * (E1 - E0);
        return (pt - proj).norm();
    }
    return std::min((pt - E0).norm(), (pt - E1).norm());
}

#include<igl/segment_segment_intersect.h>
#include<igl/project_to_line.h>
// TODO this is quadratic computation, need use tree structure
bool get_polyline_intersection(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                               const std::vector<Eigen::Vector3d> &poly0_e0, const std::vector<Eigen::Vector3d> &poly0_e1,
                               const std::vector<Eigen::Vector3d> &poly1_e0, const std::vector<Eigen::Vector3d> &poly1_e1,
                               const std::vector<int> &fid0, const std::vector<int> &fid1, Eigen::Vector3d &result, bool print = false)
{
    assert(poly0_e0.size() == poly0_e1.size());
    assert(poly1_e0.size() == poly1_e1.size());
    double dist_tol = 1e-5;
    double distance = 1e10;
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
            // double utol = dist_tol / dir0.norm();
            // double vtol = dist_tol / dir1.norm();
            // if ((u - 1.) > utol || u < -utol)// convert distance tolerance to parameter tolerance
            //     continue;

            // if ((v - 1.) > vtol || v < -vtol)
            //     continue;
            Eigen::Vector3d point = e0 + u * dir0;
            result = point;
            double dis0 = point_seg_distance(point, poly0_e0[i], poly0_e1[i]);
            double dis1 = point_seg_distance(point, poly1_e0[j], poly1_e1[j]);
            if (dis0 > dist_tol || dis1 > dist_tol)
            {
                continue;
            }

            if(print){
                std::cout<<"intersect! the intersection point is ("<<point.transpose()<<"\n";
            }

            return true;

            // if(print){
            //     std::cout<<"e0, "<<e0.transpose()<<" e0e, "<<poly0_e1[i].transpose()<<", e1, "<<e1.transpose()<<" e1e, "<<poly1_e1[j].transpose()
            //     <<", u and v, "<<u<<", "<<v<<std::endl;
            // }
        }
    }
    return false;
}

int get_the_fid_of_segment(const CGMesh &lsmesh, const CGMesh::HalfedgeHandle &h0, const CGMesh::HalfedgeHandle &h1, bool &oppo0, bool &oppo1){
    CGMesh::HalfedgeHandle h0op= lsmesh.opposite_halfedge_handle(h0);
    CGMesh::HalfedgeHandle h1op= lsmesh.opposite_halfedge_handle(h1);
    int f0 = lsmesh.face_handle(h0).idx();
    int f0op = lsmesh.face_handle(h0op).idx();
    int f1 = lsmesh.face_handle(h1).idx();
    int f1op = lsmesh.face_handle(h1op).idx();
    int fid = -1;
    if (f0 == f1)
    {
        fid = f0;
        oppo0 = false;
        oppo1 = false;
    }
    if (f0 == f1op)
    {
        fid = f0;
        oppo0 = false;
        oppo1 = true;
    }
    if (f0op == f1)
    {
        fid = f0op;
        oppo0 = true;
        oppo1 = false;
    }
    if (f0op == f1op)
    {
        fid = f0op;
        oppo0 = true;
        oppo1 = true;
    }
    return fid;
}

// the barycenter coordinates should be (u,v,t). since we assume the point is on heh which corresponds to the same face as href, t = 1-u-v.
// we can just say u and v is enough to represent the coordinates
void get_barycenter_coordinate_of_pts_on_edge(const CGMesh &lsmesh, const CGMesh::HalfedgeHandle &href, const CGMesh::HalfedgeHandle &heh,
                                              const double para, double &u, double &v)
{
    if (heh == href)
    {
        u = 0;
        v = 1 - para;
        return;
    }
    if (heh == lsmesh.next_halfedge_handle(href))
    {
        u = para;
        v = 0;
        return;
    }
    if (heh == lsmesh.prev_halfedge_handle(href))
    {
        u = 1 - para;
        v = para;
        return;
    }
    std::cout << "The href and heh do not belong to the same triangle, ERROR!!!" << std::endl;
}
bool seg_seg_intersection_2d(double u0, double u1, double u2, double u3, double v0, double v1, double v2, double v3, double &u, double &v)
{
    Eigen::Matrix2d mat, inv;
    mat(0, 0) = u1 - u0;
    mat(0, 1) = -(u3 - u2);
    mat(1, 0) = v1 - v0;
    mat(1, 1) = -(v3 - v2);
    bool invertable;

    mat.computeInverseWithCheck(inv, invertable);
    if (!invertable)
    { // very likely to be parallel
    // std::cout<<"!!!!invertable"<<std::endl;
        return false;
    }
    Eigen::Vector2d B, result;
    B[0] = u2 - u0;
    B[1] = v2 - v0;
    result = inv * B;
    double alpha = result[0];
    double beta = result[1];
    u = u0 + (u1 - u0) * alpha;
    v = v0 + (v1 - v0) * alpha;
    double t = 1 - u - v;
    // std::cout<<"computed u v t "<<u<<", "<<v<<", "<<t<<std::endl;
    if (u - 1 > SCALAR_ZERO || u < -SCALAR_ZERO)
    {
        return false;
    }
    if (v - 1 > SCALAR_ZERO || v < -SCALAR_ZERO)
    {
        return false;
    }
    
    if (t - 1 > SCALAR_ZERO || t < -SCALAR_ZERO)
    {
        return false;
    }
    return true;
}

bool seg_seg_intersection_barycenter_coordinate(const CGMesh &lsmesh, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                                                const CGMesh::HalfedgeHandle &h0, const CGMesh::HalfedgeHandle &h1,
                                                const CGMesh::HalfedgeHandle &h2, const CGMesh::HalfedgeHandle &h3,
                                                const double t0, const double t1, const double t2, const double t3, Eigen::Vector3d &pt)
{
    bool oppo0, oppo1, oppo2, oppo3;
    int fid0 = get_the_fid_of_segment(lsmesh, h0, h1, oppo0, oppo1);
    int fid1 = get_the_fid_of_segment(lsmesh, h2, h3, oppo2, oppo3);
    if (fid0 == -1 || fid1 == -1)
    {
        std::cout << "Topology wrong in seg_seg_intersection_barycenter_coordinate() !!!" << std::endl;
        return false;
    }
    if (fid0 != fid1)// the two segments are not in the same triangle.
    {
        return false;
    }
    CGMesh::HalfedgeHandle heh0 = h0, heh1= h1, heh2 = h2, heh3 = h3;// the oriented half edges
    double tc0 = t0, tc1 = t1, tc2 = t2, tc3 = t3;// the paras on the oriented half edges
    if (oppo0)
    {
        heh0 = lsmesh.opposite_halfedge_handle(heh0);
        tc0 = 1 - tc0;
    }
    if (oppo1)
    {
        heh1 = lsmesh.opposite_halfedge_handle(heh1);
        tc1 = 1 - tc1;
    }
    if (oppo2)
    {
        heh2 = lsmesh.opposite_halfedge_handle(heh2);
        tc2 = 1 - tc2;
    }
    if (oppo3)
    {
        heh3 = lsmesh.opposite_halfedge_handle(heh3);
        tc3 = 1 - tc3;
    }
    // now the half edges are correctly oriented
    double u0, v0, u1, v1, u2, v2, u3, v3, u, v;
    CGMesh::HalfedgeHandle refhd = heh0;
    get_barycenter_coordinate_of_pts_on_edge(lsmesh, refhd, heh0, tc0, u0, v0);
    get_barycenter_coordinate_of_pts_on_edge(lsmesh, refhd, heh1, tc1, u1, v1);
    get_barycenter_coordinate_of_pts_on_edge(lsmesh, refhd, heh2, tc2, u2, v2);
    get_barycenter_coordinate_of_pts_on_edge(lsmesh, refhd, heh3, tc3, u3, v3);
    // if (fid0 == 585)
    // {
    //     std::cout << "WWWWWe found this triangle" << std::endl;
    //     std::cout << "vids, " << F.row(fid0) << std::endl;
    //     std::cout << "the ref hd, " << lsmesh.from_vertex_handle(refhd).idx() << ", " << lsmesh.to_vertex_handle(refhd).idx() << std::endl;
    //     std::cout<<"u: "<<u0<<", "<<u1<<", "<<u2<<", "<<u3<<std::endl;
    //     std::cout<<"v: "<<v0<<", "<<v1<<", "<<v2<<", "<<v3<<std::endl;
    //     bool intersect = seg_seg_intersection_2d(u0, u1, u2, u3, v0, v1, v2, v3, u, v);
    // }
    // return false;

    bool intersect = seg_seg_intersection_2d(u0, u1, u2, u3, v0, v1, v2, v3, u, v);
    
    if(!intersect){
        return false;
    }
    double t = 1 - u - v;
    int vid0 = lsmesh.from_vertex_handle(heh0).idx();
    int vid1 = lsmesh.to_vertex_handle(heh0).idx();
    CGMesh::HalfedgeHandle next = lsmesh.next_halfedge_handle(heh0);
    int vid2 = lsmesh.to_vertex_handle(next).idx();

    pt = V.row(vid0) * v + V.row(vid1) * t + V.row(vid2) * u;
    return true;
}

bool get_polyline_intersection_barycenter_coordinate(const CGMesh &lsmesh, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                                                     const std::vector<std::vector<CGMesh::HalfedgeHandle>> &hlist0,
                                                     const std::vector<std::vector<CGMesh::HalfedgeHandle>> &hlist1,
                                                     const std::vector<std::vector<double>> &plist0,
                                                     const std::vector<std::vector<double>> &plist1,
                                                     Eigen::Vector3d &result, bool print = false)
{
    double dist_tol = 1e-5;
    double distance = 1e10;
    for (int i = 0; i < plist0.size(); i++)
    {
        for (int j = 0; j < plist0[i].size() - 1; j++)
        {
            double t0 = plist0[i][j];
            double t1 = plist0[i][j + 1];
            CGMesh::HalfedgeHandle h0 = hlist0[i][j];
            CGMesh::HalfedgeHandle h1 = hlist0[i][j + 1];
            for (int k = 0; k < plist1.size(); k++)
            {
                for (int l = 0; l < plist1[k].size() - 1; l++)
                {
                    double t2 = plist1[k][l];
                    double t3 = plist1[k][l + 1];
                    CGMesh::HalfedgeHandle h2 = hlist1[k][l];
                    CGMesh::HalfedgeHandle h3 = hlist1[k][l + 1];
                    Eigen::Vector3d pt;
                    bool intersect = seg_seg_intersection_barycenter_coordinate(lsmesh, V, F, h0, h1, h2, h3, t0, t1, t2, t3, pt);
                    if (intersect)
                    {
                        result = pt;
                        return true;
                    }
                }
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
    int nbr_elements = col_end - col_start + 2;
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
                                          std::vector<Eigen::VectorXi> &vcr, std::vector<Eigen::VectorXi> &vr, std::vector<Eigen::VectorXi> &vc )
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
    
    vr.clear();
    vc.clear();
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
            if (ids[1] - ids[0] + 2 > maxsize1) // ids[1] - ids[0] is the number of quads, meaning nbr of vers is +2
            {
                maxsize1 = ids[1] - ids[0] + 2;
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
            if (ids[1] - ids[0] + 2 > maxsize2) // ids[1] - ids[0] is the number of quads, meaning nbr of vers is +1
            {
                maxsize2 = ids[1] - ids[0] + 2;
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

    maxsize1 = 0;
    // get the whole polylines in each row
    for (int i = 0; i < mat.rows(); i++)
    {
        ids.clear();
        for (int j = 0; j < mat.cols(); j++)
        { // make sure there are only 2 elements, recording the first and last non -1 id
            if (mat(i, j) >= 0 && ids.size() == 0)
            {
                ids.push_back(j);
                ids.push_back(j);
            }
            if (mat(i, j) >= 0)
            {
                ids[1] = j;
            }
        }

        if (ids.size() > 0)
        {
            if (ids[1] - ids[0] + 1 > maxsize1) // ids[1] - ids[0] is the number of quads, meaning nbr of vers is +2
            {
                maxsize1 = ids[1] - ids[0] + 1;
            }
        }
        if (ids.size() > 0)
        {
            Eigen::VectorXi left = mat.row(i).segment(ids[0], ids[1] - ids[0] + 1);

            vr.push_back(left);
            
        }
    }
    vecbase1 = Eigen::VectorXi::Ones(maxsize1) * -1;
    for (int i = 0; i < vr.size(); i++)
    {
        Eigen::VectorXi tmp = vecbase1;
        int size = vr[i].size();
        tmp.segment(0, size) = vr[i];
        vr[i] = tmp;
    }

    maxsize1 = 0;
    // get the whole polylines in each col
    for (int i = 0; i < mat.cols(); i++)
    {
        ids.clear();
        for (int j = 0; j < mat.rows(); j++)
        { // make sure there are only 2 elements, recording the first and last non -1 id
            if (mat(j, i) >= 0 && ids.size() == 0)
            {
                ids.push_back(j);
                ids.push_back(j);
            }
            if (mat(j, i) >= 0)
            {
                ids[1] = j;
            }
        }

        if (ids.size() > 0)
        {
            if (ids[1] - ids[0] + 1 > maxsize1) // ids[1] - ids[0] is the number of quads, meaning nbr of vers is +2
            {
                maxsize1 = ids[1] - ids[0] + 1;
            }
        }
        if (ids.size() > 0)
        {
            Eigen::VectorXi left = mat.col(i).segment(ids[0], ids[1] - ids[0] + 1);

            vc.push_back(left);
            
        }
    }
    vecbase1 = Eigen::VectorXi::Ones(maxsize1) * -1;
    for (int i = 0; i < vc.size(); i++)
    {
        Eigen::VectorXi tmp = vecbase1;
        int size = vc[i].size();
        tmp.segment(0, size) = vc[i];
        vc[i] = tmp;
    }
}
void construct_duplication_mapping(const int vnbr, const Eigen::MatrixXi& F, Eigen::VectorXi &mapping, int& nbr)
{
    mapping = Eigen::VectorXi::Ones(vnbr) * -1;
    for (int i = 0; i < F.rows(); i++)
    {
        for (int j = 0; j < F.cols(); j++)
        {
            int vid = F(i, j);
            assert(vid >= 0);
            mapping[vid] = 1;
        }
    }
    // Eigen::MatrixXi cpmat = mat;
    
    // assert(qds.rows() + 1 == mat.rows() && qds.cols() + 1 == mat.cols());
    // // first we point out the non -1 elements in mapping
    // for (int i = 0; i < qds.rows(); i++)
    // {
    //     for (int j = 0; j < qds.cols(); j++)
    //     {
    //         if (qds(i, j) > 0)
    //         {
    //             int v0 = cpmat(i, j);
    //             int v1 = cpmat(i, j + 1);
    //             int v2 = cpmat(i + 1, j + 1);
    //             int v3 = cpmat(i + 1, j);
    //             mapping[v0] = 1;
    //             mapping[v1] = 1;
    //             mapping[v2] = 1;
    //             mapping[v3] = 1;
    //         }
    //     }
    // }
    int loc = 0;
    for (int i = 0; i < vnbr; i++)
    {
        if (mapping[i] == 1)
        {
            mapping[i] = loc;
            loc++;
        }
    }
    assert(vnbr >= loc);
    nbr = loc;
}
Eigen::MatrixXd remove_ver_duplicated(const Eigen::MatrixXd& ver, const Eigen::VectorXi& mapping, const int vnbr){
    Eigen::MatrixXd result(vnbr, 3);
    Eigen::VectorXi ck = Eigen::VectorXi::Zero(vnbr);
    int counter = 0;
    for (int i = 0; i < mapping.size(); i++)
    {
        if (mapping(i) >= 0)
        {
            result.row(mapping(i)) = ver.row(i);
            ck[mapping(i)] = 1;
            counter++;
        }
    }
    for (int i = 0; i < vnbr; i++)
    {
        assert(ck[i] == 1); // check if every ver is assigned
    }
    assert(counter == vnbr);
    return result;
}
Eigen::MatrixXi remove_fac_duplicated(const Eigen::MatrixXi &f, const Eigen::MatrixXi &mapping)
{
    Eigen::MatrixXi result(f.rows(), f.cols());
    for (int i = 0; i < f.rows(); i++)
    {
        for (int j = 0; j < f.cols(); j++)
        {
            result(i, j) = mapping(f(i, j));
            assert(result(i, j) >= 0);
        }
    }
    // std::cout<<"F\n"<<result<<std::endl;
    return result;
}
std::vector<Eigen::VectorXi> remove_quad_info_duplicated(const std::vector<Eigen::VectorXi> &vrl, const Eigen::MatrixXi &mapping)
{
    int size = vrl[0].size();
    std::vector<Eigen::VectorXi> result(vrl.size());
    for (int i = 0; i < vrl.size(); i++)
    {
        result[i].resize(vrl[i].size());
        for (int j = 0; j < vrl[i].size(); j++)
        {
            if (vrl[i][j] >= 0)
            {
                result[i][j] = mapping(vrl[i][j]);
            }
            else
            {
                result[i][j] = vrl[i][j];
            }
        }
    }
    // filter out the first few -1s
    std::vector<Eigen::VectorXi> filtered;
    for (int i = 0; i < result.size(); i++)
    {
        int start = -1;
        Eigen::VectorXi tmp = Eigen::VectorXi::Ones(size) * -1;
        for (int j = 0; j < size; j++)
        {
            if (result[i][j] >= 0)
            {
                start = j;
                break;
            }
        }
        if (start == 0)
        {
            // std::cout<<"here 1, size "<<result[i].size()<<std::endl;
            filtered.push_back(result[i]);
            continue;
        }
        if (start == -1)
        {
            continue;
        }
        // std::cout<<"here 2"<<std::endl;
        int nbr_eles = size - start;
        tmp.topRows(nbr_eles) = result[i].segment(start, nbr_eles);
        filtered.push_back(tmp);
    }
    return filtered;
}
void extract_web_from_index_mat(const Eigen::MatrixXi &mat, Eigen::MatrixXi &F, const int threads, std::vector<Eigen::VectorXi> &vrl,
                                std::vector<Eigen::VectorXi> &vrr, std::vector<Eigen::VectorXi> &vcl,
                                std::vector<Eigen::VectorXi> &vcr, std::vector<Eigen::VectorXi> &vr, std::vector<Eigen::VectorXi> &vc, Eigen::MatrixXi &qds)
{
    std::array<int, 4> face;
    std::vector<std::array<int, 4>> tface, filted_face;
    int row = mat.rows();
    int col = mat.cols();
    qds = Eigen::MatrixXi::Ones(row - 1, col - 1) * -1; // quad mat initialized as -1
    tface.reserve(mat.rows() * mat.cols());
    filted_face.reserve(mat.rows() * mat.cols());
    for (int i = 0; i < mat.rows() - 1; i++)
    {
        for (int j = 0; j < mat.cols() - 1; j++)
        {
            int f0 = mat(i, j);
            int f1 = mat(i, j + 1);
            int f2 = mat(i + 1, j + 1);
            int f3 = mat(i + 1, j);
            if (f0 < 0 || f1 < 0 || f2 < 0 || f3 < 0)
            {
                continue;
            }
            face = {f0, f1, f2, f3};
            qds(i, j) = tface.size();
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
    else
    {
        filted_face = tface;
    }

    F.resize(filted_face.size(), 4);
    for (int i = 0; i < filted_face.size(); i++)
    {
        F(i, 0) = filted_face[i][0];
        F(i, 1) = filted_face[i][1];
        F(i, 2) = filted_face[i][2];
        F(i, 3) = filted_face[i][3];
    }
    get_row_and_col_ver_ids_from_web_mat(mat, qds, vrl, vrr, vcl, vcr, vr, vc);
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
                               std::vector<Eigen::VectorXi> &vcr, std::vector<Eigen::VectorXi> &vr, std::vector<Eigen::VectorXi> &vc)
{
    std::ofstream file;
    std::vector<Eigen::VectorXi> current;
    std::string fname = namebase + "_rowleft.csv";
    current = vrl;
    write_one_info_file(fname, current, file);

    fname = namebase + "_rowright.csv";
    current = vrr;
    write_one_info_file(fname, current, file);

    fname = namebase + "_colleft.csv";
    current = vcl;
    write_one_info_file(fname, current, file);

    fname = namebase + "_colright.csv";
    current = vcr;
    write_one_info_file(fname, current, file);

    fname = namebase + "_rows.csv";
    current = vr;
    write_one_info_file(fname, current, file);

    fname = namebase + "_cols.csv";
    current = vc;
    write_one_info_file(fname, current, file);
}
void save_strokes(const std::vector<std::vector<int>> &flist, const std::vector<std::vector<Eigen::Vector3f>> &bclist, 
const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
{
    std::cout << "Saving the info for the strokes, please provide the prefix: " << std::endl;
    std::string fname = igl::file_dialog_save();
    std::ofstream file;
    file.open(fname + "_flist.csv");
    for (int i = 0; i < flist.size(); i++)
    {
        for (int j = 0; j < flist[i].size() - 1; j++)
        {
            file << flist[i][j] << ",";
        }
        file << flist[i].back() << std::endl;
    }
    file.close();
    file.open(fname + "_bc.csv");
    for (int i = 0; i < bclist.size(); i++)
    {
        for (int j = 0; j < bclist[i].size() - 1; j++)
        {
            file << bclist[i][j][0] << ","<<bclist[i][j][1]<<","<<bclist[i][j][2]<<",";
        }
        file << bclist[i].back()[0] << "," << bclist[i].back()[1] << "," << bclist[i].back()[2] << "," << std::endl;
    }
    file.close();
    std::ofstream filex, filey, filez;
    for (int i = 0; i < flist.size(); i++)
    {
        filex.open(fname + "_" + std::to_string(i) + "x.csv");
        filey.open(fname + "_" + std::to_string(i) + "y.csv");
        filez.open(fname + "_" + std::to_string(i) + "z.csv");
        for (int j = 0; j < flist[i].size(); j++)
        {
            Eigen::Vector3f bc = bclist[i][j];
            int fid = flist[i][j];
            Eigen::Vector3d v0 = V.row(F(fid,0));
            Eigen::Vector3d v1 = V.row(F(fid,1));
            Eigen::Vector3d v2 = V.row(F(fid,2));
            Eigen::Vector3d pt = bc[0] * v0 + bc[1] * v1 + bc[2] * v2;
            filex << pt[0];
            filey << pt[1];
            filez << pt[2];
            if (j < flist[i].size() - 1)
            {
                filex << ",";
                filey << ",";
                filez << ",";
            }
        }
        filex << "\n";
        filey << "\n";
        filez << "\n";
        filex.close();
        filey.close();
        filez.close();
    }
    std::cout<<"files saved"<<std::endl;
}
void read_strokes( std::vector<std::vector<int>> &flist,  std::vector<std::vector<Eigen::Vector3f>> &bclist){
    std::cout<<"Reading stroke files, please provide the prefix: "<<std::endl;
    std::string fname = igl::file_dialog_save();
    std::ofstream file;
    std::vector<std::vector<double>> tmp_f, tmp_v;
    read_csv_data_lbl(fname + "_flist.csv",tmp_f);
    read_csv_data_lbl(fname + "_bc.csv",tmp_v);
    flist.resize(tmp_f.size());
    bclist.resize(tmp_f.size());
    for(int i=0;i<tmp_f.size();i++){
        std::vector<int> silist;
        std::vector<Eigen::Vector3f> sivers;
        for(int j=0;j<tmp_f[i].size();j++){
            silist.push_back(tmp_f[i][j]);
            sivers.push_back(Eigen::Vector3f(tmp_v[i][j * 3], tmp_v[i][j * 3 + 1], tmp_v[i][j * 3 + 2]));
        }
        flist[i]=silist;
        bclist[i]=sivers;
    }
}
void read_ver_write_into_xyz_files()
{
    std::string fname = igl::file_dialog_open();
    if(fname.empty()){
        std::cout << "please read a proper ver list in obj format" << std::endl;
    }
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::readOBJ(fname, V, F);
    std::ofstream file;
    file.open(fname + "_x.csv");
    for (int i = 0; i < V.rows(); i++)
    {
        file << V(i, 0);
        if (i < V.rows() - 1)
        {
            file << ",";
        }
    }
    file.close();
    file.open(fname + "_y.csv");
    for (int i = 0; i < V.rows(); i++)
    {
        file << V(i, 1);
        if (i < V.rows() - 1)
        {
            file << ",";
        }
    }
    file.close();
    file.open(fname + "_z.csv");
    for (int i = 0; i < V.rows(); i++)
    {
        file << V(i, 2);
        if (i < V.rows() - 1)
        {
            file << ",";
        }
    }
    file.close();

    std::cout<<"files saved"<<std::endl;
}
std::vector<int> element_ids_for_not_computed_vers(const Eigen::VectorXi &vec)
{
    std::vector<int> tmp, result;
    int found = -1; // record the last one has value
    for (int i = 0; i < vec.size(); i++)
    {
        if (vec[i] >= 0)
        {
            found = i;
        }
        else
        {
            if (found >= 0)
            {
                tmp.push_back(i);
            }
        }
    }
    for (int i = 0; i < tmp.size(); i++)
    {
        if (tmp[i] < found)
        {
            result.push_back(tmp[i]);
        }
    }
    return result;
}

// threadshold_nbr is the nbr of quads we want to discard when too few quads in a row
void visual_extract_levelset_web(const CGMesh &lsmesh, const Eigen::MatrixXd &V,
                                 const Eigen::MatrixXi &F, const Eigen::VectorXd &ls0, const Eigen::VectorXd &ls1,
                                 const int expect_nbr_ls0, const int expect_nbr_ls1, Eigen::MatrixXd &E0, Eigen::MatrixXd &E1,
                                 Eigen::MatrixXd &E2, Eigen::MatrixXd &E3,
                                 bool even_pace, bool debug, int dbg0, int dbg1)
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
    for (int i = 0; i < nbr_ls0; i++)
    {
        if(debug){
            if(i != dbg0){
                continue;
            }
        }
        double vl = lsv0[i];
        get_iso_lines(V, F, ls0, vl, poly0_e0[i], poly0_e1[i], fid0[i]);
        // std::cout<<i<<" size "<<poly0_e0[i].size()<<"\n";
    }
    for (int i = 0; i < nbr_ls1; i++)
    {
        if(debug){
            if(i != dbg1){
                continue;
            }
        }
        double vl = lsv1[i];
        get_iso_lines(V, F, ls1, vl, poly1_e0[i], poly1_e1[i], fid1[i]);
        // std::cout<<i<<" size "<<poly1_e0[i].size()<<"\n";
    }
    int enbr = 0;
    for (auto pl : poly0_e0)
    {
        enbr += pl.size();
    }
    E0.resize(enbr, 3);
    E1.resize(enbr, 3);
    
    enbr = 0;
    for (auto pl : poly1_e0)
    {
        enbr += pl.size();
    }
    E2.resize(enbr, 3);
    E3.resize(enbr, 3);

    enbr = 0;
    for (int i = 0; i < poly0_e0.size(); i++)
    {
        for (int j = 0; j < poly0_e0[i].size(); j++)
        {
            E0.row(enbr) = poly0_e0[i][j];
            E1.row(enbr) = poly0_e1[i][j];
            enbr++;
        }
    }
    enbr = 0;
    for (int i = 0; i < poly1_e0.size(); i++)
    {
        for (int j = 0; j < poly1_e0[i].size(); j++)
        {
            E2.row(enbr) = poly1_e0[i][j];
            E3.row(enbr) = poly1_e1[i][j];
            enbr++;
        }
    }
    if(debug){
        Eigen::Vector3d ipoint;
        bool intersect = get_polyline_intersection(V, F, poly0_e0[dbg0], poly0_e1[dbg0], poly1_e0[dbg1], poly1_e1[dbg1], fid0[dbg0], fid1[dbg1], ipoint, true);
        std::cout<<"intersect? "<<intersect<<"\nintersection point\n"<<ipoint.transpose()<<std::endl;
    }
    
}
// threadshold_nbr is the nbr of quads we want to discard when too few quads in a row
void visual_extract_levelset_web_stable(const CGMesh &lsmesh, const std::vector<CGMesh::HalfedgeHandle>& loop, const Eigen::MatrixXd &V,
                                 const Eigen::MatrixXi &F, const Eigen::VectorXd &ls0, const Eigen::VectorXd &ls1,
                                 const int expect_nbr_ls0, const int expect_nbr_ls1, Eigen::MatrixXd &E0, Eigen::MatrixXd &E1,
                                 Eigen::MatrixXd &E2, Eigen::MatrixXd &E3,
                                 bool even_pace, bool debug, int dbg0, int dbg1)
{
    std::vector<std::vector<Eigen::Vector3d>> ivs;                                    // the intersection vertices
    Eigen::MatrixXi gridmat;                                                          // the matrix for vertex
    Eigen::VectorXd lsv0, lsv1;                                                       // the extracted level set values;
    std::vector<std::vector<Eigen::Vector3d>> poly0_e0, poly0_e1, poly1_e0, poly1_e1; // polylines
    std::vector<std::vector<std::vector<double>>> paras0, paras1;
    std::vector<std::vector<std::vector<CGMesh::HalfedgeHandle>>> hds0, hds1; // face ids of each polyline segment
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
    paras0.resize(nbr_ls0);
    paras1.resize(nbr_ls1);
    hds0.resize(nbr_ls0);
    hds1.resize(nbr_ls1);
    
    // std::cout<<"sp_0 \n"<<lsv0.transpose()<<"\nsp_1\n"<<lsv1.transpose()<<std::endl;
    for (int i = 0; i < nbr_ls0; i++)
    {
        double vl = lsv0[i];
        std::vector<bool> left_large;
        get_iso_lines_baricenter_coord(lsmesh, loop, V, F, ls0, vl, paras0[i], hds0[i], left_large);
        // std::cout<<i<<" size "<<poly0_e0[i].size()<<"\n";
    }
    for (int i = 0; i < nbr_ls1; i++)
    {
        double vl = lsv1[i];
        std::vector<bool> left_large;
        get_iso_lines_baricenter_coord(lsmesh, loop, V, F, ls1, vl, paras1[i], hds1[i], left_large);
        // std::cout<<i<<" size "<<poly1_e0[i].size()<<"\n";
    }
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
    for (int i = 0; i < nbr_ls0; i++)
    {
        if(debug){
            if(i != dbg0){
                continue;
            }
        }
        double vl = lsv0[i];
        get_iso_lines(V, F, ls0, vl, poly0_e0[i], poly0_e1[i], fid0[i]);
        // std::cout<<i<<" size "<<poly0_e0[i].size()<<"\n";
    }
    for (int i = 0; i < nbr_ls1; i++)
    {
        if(debug){
            if(i != dbg1){
                continue;
            }
        }
        double vl = lsv1[i];
        get_iso_lines(V, F, ls1, vl, poly1_e0[i], poly1_e1[i], fid1[i]);
        // std::cout<<i<<" size "<<poly1_e0[i].size()<<"\n";
    }
    int enbr = 0;
    for (auto pl : poly0_e0)
    {
        enbr += pl.size();
    }
    E0.resize(enbr, 3);
    E1.resize(enbr, 3);
    
    enbr = 0;
    for (auto pl : poly1_e0)
    {
        enbr += pl.size();
    }
    E2.resize(enbr, 3);
    E3.resize(enbr, 3);

    enbr = 0;
    for (int i = 0; i < poly0_e0.size(); i++)
    {
        for (int j = 0; j < poly0_e0[i].size(); j++)
        {
            E0.row(enbr) = poly0_e0[i][j];
            E1.row(enbr) = poly0_e1[i][j];
            enbr++;
        }
    }
    enbr = 0;
    for (int i = 0; i < poly1_e0.size(); i++)
    {
        for (int j = 0; j < poly1_e0[i].size(); j++)
        {
            E2.row(enbr) = poly1_e0[i][j];
            E3.row(enbr) = poly1_e1[i][j];
            enbr++;
        }
    }
    if(debug){
        Eigen::Vector3d ipoint;
        bool intersect = get_polyline_intersection_barycenter_coordinate(lsmesh, V, F, hds0[dbg0], hds1[dbg1], paras0[dbg0], paras1[dbg1], ipoint, true);
        std::cout << "intersect? " << intersect << "\nintersection point\n"
                  << ipoint.transpose() << std::endl;
    }
    
}

// void filter_out_isolate_vers(void)

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
                // double dis = 
                verlist.push_back(ipoint);
                gridmat(i, j) = vnbr;
                vnbr++;
            }
        }
        // The following pars are for debugging
        // Eigen::VectorXi cv=gridmat.row(i);
        // std::vector<int> problems = element_ids_for_not_computed_vers(cv);
        // if(problems.size()>0){
        //     std::cout<<i<<" this row has missing points \n"<<cv.transpose()<<" size, "<<problems.size()<<std::endl;
            
        //     for(int j: problems){
        //         Eigen::Vector3d ipoint;
        //         get_polyline_intersection(V, F, poly0_e0[i], poly0_e1[i], poly1_e0[j], poly1_e1[j], fid0[i], fid1[j], ipoint, true);
        //     }
        // }
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
    std::vector<Eigen::VectorXi> vr;
    std::vector<Eigen::VectorXi> vc;
    Eigen::MatrixXi qds;
    extract_web_from_index_mat(gridmat, Faces, threadshold_nbr, vrl, vrr, vcl, vcr, vr, vc, qds);

    // remove duplicated vertices
    Eigen::VectorXi mapping;
    int real_nbr;
    construct_duplication_mapping(vers.rows(), Faces, mapping, real_nbr);
    // check if mapping is correct.
    // for(int i=0;i<mapping.size();i++){

    // }


    std::cout<<"mapping constructed"<<std::endl;
    // std::cout<<mapping<<std::endl;
    vers = remove_ver_duplicated(vers, mapping, real_nbr);
    Faces = remove_fac_duplicated(Faces, mapping);
    vrl = remove_quad_info_duplicated(vrl, mapping);
    vrr = remove_quad_info_duplicated(vrr, mapping);
    vcl = remove_quad_info_duplicated(vcl, mapping);
    vcr = remove_quad_info_duplicated(vcr, mapping);
    vr = remove_quad_info_duplicated(vr, mapping);
    vc = remove_quad_info_duplicated(vc, mapping);
    
    std::cout<<"Saving the info for each row or col of quads"<<std::endl;
    std::string fname = igl::file_dialog_save();
    save_quad_left_right_info(fname, vrl, vrr, vcl, vcr, vr, vc);
    std::cout<<"Each row and col of quads got saved"<<std::endl;
}

void getWebFromMat(const Eigen::MatrixXi &mat, Eigen::VectorXi &l0, Eigen::VectorXi &l1, Eigen::VectorXi &l2, Eigen::VectorXi &l3,
                   Eigen::VectorXi &l4, Eigen::VectorXi &l5)
{
    std::vector<int> e0, e1, e2, e3, e4, e5;
    int rows = mat.rows();
    int cols = mat.cols();
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < cols; j++)
        {
            if (j < cols - 1)
            {
                if (mat(i, j) >= 0 && mat(i, j + 1) >= 0)
                {
                    e0.push_back(mat(i, j));
                    e1.push_back(mat(i, j + 1));
                }
            }
            if (i < rows - 1)
            {
                if (mat(i, j) >= 0 && mat(i + 1, j) >= 0)
                {
                    e2.push_back(mat(i, j));
                    e3.push_back(mat(i + 1, j));
                }
            }

            if (i < rows - 1 && j < cols + 1)
            {
                if (mat(i, j + 1) >= 0 && mat(i + 1, j) >= 0)
                {
                    e4.push_back(mat(i, j + 1));
                    e5.push_back(mat(i + 1, j));
                }
            }
        }
    }
    l0 = vec_list_to_vector(e0);
    l1 = vec_list_to_vector(e1);
    l2 = vec_list_to_vector(e2);
    l3 = vec_list_to_vector(e3);
    l4 = vec_list_to_vector(e4);
    l5 = vec_list_to_vector(e5);
}

void writeEdgeFile(const std::string& name, const Eigen::MatrixXd& V, const Eigen::VectorXi& l0, const Eigen::VectorXi& l1){
    std::ofstream fout;
    fout.open(name);
    for (int i = 0; i < V.rows(); i++)
    {
        fout << "v " << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << "\n";
    }
    for(int i=0;i<l0.size();i++){
        fout<<"l "<<l0[i]+1<<" "<<l1[i]+1<<"\n";
    }
    fout.close();
}

// threadshold_nbr is the nbr of quads we want to discard when too few quads in a row
// this function uses the baricenter coordinates to solve intersection points
// diagSegments is true then write the diagonal segments.
void extract_levelset_web_stable(const CGMesh &lsmesh, const std::vector<CGMesh::HalfedgeHandle>& loop, const Eigen::MatrixXd &V,
                          const Eigen::MatrixXi &F, const Eigen::VectorXd &ls0, const Eigen::VectorXd &ls1,
                          const int expect_nbr_ls0, const int expect_nbr_ls1, const int threadshold_nbr,
                          Eigen::MatrixXd &vers, Eigen::MatrixXi &Faces, bool even_pace=false, bool diagSegments = false)
{
    Eigen::MatrixXi gridmat;    // the matrix for vertex
    Eigen::VectorXd lsv0, lsv1; // the extracted level set values;
    std::vector<std::vector<std::vector<double>>> paras0, paras1;
    std::vector<std::vector<std::vector<CGMesh::HalfedgeHandle>>> hds0, hds1; // face ids of each polyline segment
    std::vector<Eigen::Vector3d> verlist;
    if (ls0.size() != ls1.size())
    {
        std::cout << "ERROR, Please use the correct level sets" << std::endl;
    }
    int nbr_ls0, nbr_ls1;
    if (!even_pace)
    {
        nbr_ls0 = expect_nbr_ls0;
        nbr_ls1 = expect_nbr_ls1;
        get_level_set_sample_values(ls0, nbr_ls0, lsv0);
        get_level_set_sample_values(ls1, nbr_ls1, lsv1);
    }
    else
    {
        double pace = std::min((ls0.maxCoeff() - ls0.minCoeff()) / (expect_nbr_ls0 + 1), (ls1.maxCoeff() - ls1.minCoeff()) / (expect_nbr_ls1 + 1));
        nbr_ls0 = (ls0.maxCoeff() - ls0.minCoeff()) / pace + 1;
        nbr_ls1 = (ls1.maxCoeff() - ls1.minCoeff()) / pace + 1;
        get_level_set_sample_values_even_pace(ls0, nbr_ls0, pace, lsv0);
        get_level_set_sample_values_even_pace(ls1, nbr_ls1, pace, lsv1);
    }

    verlist.reserve(nbr_ls0 * nbr_ls1);
    gridmat = Eigen::MatrixXi::Ones(nbr_ls0, nbr_ls1) * -1; // initially there is no quad patterns
    lsv0.resize(nbr_ls0);
    lsv1.resize(nbr_ls1);

    paras0.resize(nbr_ls0);
    paras1.resize(nbr_ls1);
    hds0.resize(nbr_ls0);
    hds1.resize(nbr_ls1);
    
    // std::cout<<"sp_0 \n"<<lsv0.transpose()<<"\nsp_1\n"<<lsv1.transpose()<<std::endl;
    for (int i = 0; i < nbr_ls0; i++)
    {
        double vl = lsv0[i];
        std::vector<bool> left_large;
        get_iso_lines_baricenter_coord(lsmesh, loop, V, F, ls0, vl, paras0[i], hds0[i], left_large);
        // std::cout<<i<<" size "<<poly0_e0[i].size()<<"\n";
    }
    for (int i = 0; i < nbr_ls1; i++)
    {
        double vl = lsv1[i];
        std::vector<bool> left_large;
        get_iso_lines_baricenter_coord(lsmesh, loop, V, F, ls1, vl, paras1[i], hds1[i], left_large);
        // std::cout<<i<<" size "<<poly1_e0[i].size()<<"\n";
    }
    int vnbr = 0;
    for (int i = 0; i < nbr_ls0; i++)
    {
        for (int j = 0; j < nbr_ls1; j++)
        {
            Eigen::Vector3d ipoint;
            bool intersect = get_polyline_intersection_barycenter_coordinate(lsmesh, V, F, hds0[i], hds1[j], paras0[i], paras1[j], ipoint);
            if (intersect)
            {
                // double dis = 
                verlist.push_back(ipoint);
                gridmat(i, j) = vnbr;
                vnbr++;
            }
        }
        // The following pars are for debugging
        // Eigen::VectorXi cv=gridmat.row(i);
        // std::vector<int> problems = element_ids_for_not_computed_vers(cv);
        // if(problems.size()>0){
        //     std::cout<<i<<" this row has missing points \n"<<cv.transpose()<<" size, "<<problems.size()<<std::endl;
            
        //     for(int j: problems){
        //         Eigen::Vector3d ipoint;
        //         get_polyline_intersection(V, F, poly0_e0[i], poly0_e1[i], poly1_e0[j], poly1_e1[j], fid0[i], fid1[j], ipoint, true);
        //     }
        // }
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
    std::vector<Eigen::VectorXi> vr;
    std::vector<Eigen::VectorXi> vc;
    Eigen::MatrixXi qds;
    extract_web_from_index_mat(gridmat, Faces, threadshold_nbr, vrl, vrr, vcl, vcr, vr, vc, qds);

    // remove duplicated vertices
    Eigen::VectorXi mapping;
    int real_nbr;
    construct_duplication_mapping(vers.rows(), Faces, mapping, real_nbr);
    // check if mapping is correct.
    // for(int i=0;i<mapping.size();i++){

    // }


    std::cout<<"mapping constructed"<<std::endl;
    // std::cout<<mapping<<std::endl;
    vers = remove_ver_duplicated(vers, mapping, real_nbr);
    Faces = remove_fac_duplicated(Faces, mapping);
    vrl = remove_quad_info_duplicated(vrl, mapping);
    vrr = remove_quad_info_duplicated(vrr, mapping);
    vcl = remove_quad_info_duplicated(vcl, mapping);
    vcr = remove_quad_info_duplicated(vcr, mapping);
    vr = remove_quad_info_duplicated(vr, mapping);
    vc = remove_quad_info_duplicated(vc, mapping);
    if (!diagSegments)
    {
        std::cout << "Saving the info for each row or col of quads" << std::endl;
        std::string fname = igl::file_dialog_save();
        save_quad_left_right_info(fname, vrl, vrr, vcl, vcr, vr, vc);
        std::cout << "Each row and col of quads got saved" << std::endl;
    }

    if(diagSegments){
        Eigen::VectorXi l0, l1, l2, l3, l4, l5;
        getWebFromMat(gridmat, l0, l1, l2, l3, l4, l5);
        std::cout<<"Before writting the web edges, prefix: \n";
        Eigen::MatrixXd Vall = vec_list_to_matrix(verlist);
        std::string fname = igl::file_dialog_save();
        writeEdgeFile(fname + "R.obj", Vall, l0, l1);
        writeEdgeFile(fname + "C.obj", Vall, l2, l3);
        writeEdgeFile(fname + "D.obj", Vall, l4, l5);
        std::cout<<"after writting the web edges: \n";
    }
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
// 
void find_next_pt_on_polyline(const int start_seg, const std::vector<Eigen::Vector3d> &polyline, const double length,
                              const Eigen::Vector3d &pstart, int &seg, Eigen::Vector3d &pt, double &tlocal)
{
    int nbr = polyline.size();
    double ocu_dis = 0;
    double dis_to_start = length + (polyline[start_seg] - pstart).norm();
    seg = -1;
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
            tlocal = 1- t;
            return;
        }
    }
    std::cout<<"ERROR OUT OF SEGMENT"<<std::endl;
}
// pstart is the start point where t = 0 on this seg
void find_prev_pt_on_polyline(const int start_seg, const std::vector<Eigen::Vector3d> &polyline, const double length,
                              const Eigen::Vector3d &pstart, int &seg, Eigen::Vector3d &pt)
{
    assert(length >= 0);
    int nbr = polyline.size();
    double ocu_dis = 0;
    double dis_to_start = length - (polyline[start_seg] - pstart).norm();
    if (dis_to_start < 0)
    {
        seg = start_seg;
        pt = polyline[seg] + (polyline[seg + 1] - polyline[seg]) * dis_to_start / (polyline[seg + 1] - polyline[seg]).norm();
        return;
    }
    seg = -1;
    for (int i = start_seg - 1; i >= 0; i--)
    {
        double dis = (polyline[i] - polyline[i + 1]).norm();
        ocu_dis += dis;
        if (ocu_dis > dis_to_start)
        {                                           // the point should between i and i+1
            double diff = (ocu_dis - dis_to_start); // the distance the point to i+1
            double t = diff / dis;
            assert(t >= 0 && t <= 1);
            Eigen::Vector3d p = get_3d_ver_from_t(t, polyline[i], polyline[i + 1]);
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
    if(nbr < 3){
        std::cout<<"nbr is less than 3: "<<nbr<<std::endl;
    }
    assert(nbr >= 3);
    double avg = length / (nbr - 1);
    verlist.push_back(polyline[0]);
    int start=0;
    for (int i = 0; i < nbr - 2; i++)
    {
        int seg;
        Eigen::Vector3d pt;
        double tlocal;
        find_next_pt_on_polyline(start, polyline, avg, verlist.back(), seg, pt, tlocal);
        verlist.push_back(pt);
        start=seg;
    }
    verlist.push_back(polyline.back());
}
void sample_polyline_and_extend_verlist(const std::vector<Eigen::Vector3d> &polyline, const int nbr, const double avg, std::vector<Eigen::Vector3d> &verlist,
                                        std::vector<int> &segid, std::vector<double> &t)
{
    if (nbr < 3)
    {
        std::cout << "nbr is less than 3: " << nbr << std::endl;
    }
    assert(nbr >= 3);
    // double avg = length / (nbr - 1); // this is the relationship between avg and nbr
    verlist.push_back(polyline[0]);
    segid.push_back(0);
    t.push_back(0);
    int start = 0;
    for (int i = 0; i < nbr - 2; i++)
    {
        int seg;
        Eigen::Vector3d pt;
        double tl;
        find_next_pt_on_polyline(start, polyline, avg, verlist.back(), seg, pt, tl);
        verlist.push_back(pt);
        segid.push_back(seg);
        double local_t = get_t_of_segment(pt, polyline[seg], polyline[seg + 1]);
        t.push_back(local_t);
        assert(local_t >= 0 && local_t <= 1);
        start = seg;
    }
    verlist.push_back(polyline.back());
    segid.push_back(polyline.size() - 2);
    t.push_back(1);
    assert(t.size()==segid.size()&& t.size()==verlist.size());
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
std::vector<Eigen::Vector3d> sample_one_polyline_and_binormals_based_on_length(const std::vector<Eigen::Vector3d> &polyline, const int nbr,
                                                                               const std::vector<Eigen::Vector3d> &binormals, std::vector<Eigen::Vector3d> &bn_out)
{
    double length = polyline_length(polyline);
    double avg = length / (nbr - 1);
    std::vector<Eigen::Vector3d> verlist;
    std::vector<int> segid;
    std::vector<double> tlist;
    // std::cout<<"check in"<<std::endl;
    sample_polyline_and_extend_verlist(polyline, nbr, avg, verlist, segid, tlist);
    // std::cout<<"check out"<<std::endl;
    for(int i=0;i<verlist.size();i++){
        int id0 = segid[i];
        double t = tlist[i];
        Eigen::Vector3d b1 = binormals[id0], b2 = binormals[id0 + 1];
        if (b1.dot(b2) < 0)
        {
            b2 *= -1;
        }
        Eigen::Vector3d bn = b1 * (1 - t) + b2 * t;
        bn.normalize();
        if(isnan(bn[0])||isnan(bn[1]) ||isnan(bn[2])){
            std::cout<<"Sampled Binormal has NAN"<<std::endl;
            std::cout<<"b1, "<<b1.transpose()<<std::endl;
            std::cout<<"b2, "<<b2.transpose()<<std::endl;
            std::cout<<"t, "<<t<<std::endl;
            std::cout<<"id0, "<<id0<<std::endl;
            std::cout<<"line size , "<<binormals.size()<<std::endl;

        }

        bn_out.push_back(bn);
    }
    return verlist;
    // std::cout<<"check out out out"<<std::endl;
}
std::vector<Eigen::Vector3d> sample_one_polyline_based_on_length(const std::vector<Eigen::Vector3d> &polyline, const int nbr)
{
    double length = polyline_length(polyline);
    double avg = length / (nbr - 1);
    std::vector<Eigen::Vector3d> verlist;
    std::vector<int> segid;
    std::vector<double> tlist;
    // std::cout<<"check in"<<std::endl;
    sample_polyline_and_extend_verlist(polyline, nbr, avg, verlist, segid, tlist);
    // std::cout<<"check out"<<std::endl;
    // for(int i=0;i<verlist.size();i++){
    //     int id0 = segid[i];
    //     double t = tlist[i];
    //     Eigen::Vector3d b1 = binormals[id0], b2 = binormals[id0 + 1];
    //     if (b1.dot(b2) < 0)
    //     {
    //         b2 *= -1;
    //     }
    //     Eigen::Vector3d bn = b1 * (1 - t) + b2 * t;
    //     bn.normalize();
    //     if(isnan(bn[0])||isnan(bn[1]) ||isnan(bn[2])){
    //         std::cout<<"Sampled Binormal has NAN"<<std::endl;
    //         std::cout<<"b1, "<<b1.transpose()<<std::endl;
    //         std::cout<<"b2, "<<b2.transpose()<<std::endl;
    //         std::cout<<"t, "<<t<<std::endl;
    //         std::cout<<"id0, "<<id0<<std::endl;

    //     }
    // }
    return verlist;
    // std::cout<<"check out out out"<<std::endl;
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

void write_polyline_xyz(const std::vector<std::vector<Eigen::Vector3d>> &lines, const std::string prefix)
{
    std::ofstream file;
    file.open(prefix + "_x.csv");
    for (int i = 0; i < lines.size(); i++)
    {
        for (int j = 0; j < lines[i].size() - 1; j++)
        {
            file << lines[i][j][0] << ",";
        }
        file << lines[i].back()[0] << std::endl;
    }
    file.close();
    file.open(prefix + "_y.csv");
    for (int i = 0; i < lines.size(); i++)
    {
        for (int j = 0; j < lines[i].size() - 1; j++)
        {
            file << lines[i][j][1] << ",";
        }
        file << lines[i].back()[1] << std::endl;
    }
    file.close();
    file.open(prefix + "_z.csv");
    for (int i = 0; i < lines.size(); i++)
    {
        for (int j = 0; j < lines[i].size() - 1; j++)
        {
            file << lines[i][j][2] << ",";
        }
        file << lines[i].back()[2] << std::endl;
    }
    file.close();
}
Eigen::Vector3d orient_vector(const Eigen::Vector3d& base, const Eigen::Vector3d &vec){
    if(vec.dot(base)<0){
        return -vec;
    }
    return vec;
}
std::array<double, 3> barycenter_coordinate(const Eigen::Vector3d &v0, const Eigen::Vector3d &v1, const Eigen::Vector3d &v2, const Eigen::Vector3d &p)
{
    double t0 = (v1 - v2).cross(v2 - p).norm();
    double t1 = (v0 - v2).cross(v2 - p).norm();
    double t2 = (v1 - v0).cross(v1 - p).norm();
    double all = (v1 - v0).cross(v1 - v2).norm();
    if (t0 < -1e-3 || t1 < -1e-3 || t2 < -1e-3)
    {
        std::cout << "COMPUTATION ERROR IN barycenter_coordinate()" << std::endl;
        std::cout << "v0, " << v0.transpose() << std::endl;
        std::cout << "v1, " << v1.transpose() << std::endl;
        std::cout << "v2, " << v2.transpose() << std::endl;
        std::cout << "p, " << p.transpose() << std::endl;
        std::cout<<"coor, "<<t0<<", "<<t1<<", "<<t2<<" all, "<<all<<std::endl;

    }
        // if ((t0 + t1 + t2) / all > 1 + 1e-3 || (t0 + t1 + t2) / all < 1 - 1e-3)
        // {
        //     std::cout << "barycenter coordinate not accurate, " << t0 + t1 + t2 << " " << all << std::endl;
        //     std::cout <<v0.transpose()<< "\n"<<v1.transpose()<< "\n"<<v2.transpose()<< "\n"<<p<< "\n";
        // }
    std::array<double, 3> result;
    all = t0 + t1 + t2;
    result[0] = t0 / all;
    result[1] = t1 / all;
    result[2] = t2 / all;
    return result;
}

void generate_estimate_binormals_out_of_polylines(const std::vector<std::vector<Eigen::Vector3d>>& plys, std::vector<std::vector<Eigen::Vector3d>>& bnms){
    bnms = plys;
    for (int i = 0; i < plys.size(); i++)
    {
        for (int j = 1; j < plys[i].size() - 1; j++)
        {
            Eigen::Vector3d vf = plys[i][j - 1];
            Eigen::Vector3d v = plys[i][j];
            Eigen::Vector3d vb = plys[i][j + 1];

            Eigen::Vector3d bvec = (vf - v).cross(vb - v).normalized();
            bnms[i][j] = bvec;
        }
        int csize = plys[i].size();
        bnms[i][0] = bnms[i][1];
        bnms[i][csize - 1] = bnms[i][csize - 2];
    }
    return;
}

#include <igl/point_mesh_squared_distance.h>
void extract_shading_lines(const CGMesh &lsmesh, const Eigen::MatrixXd &V, const std::vector<CGMesh::HalfedgeHandle> &loop,
                           const Eigen::MatrixXi &F, const Eigen::VectorXd &ls,
                           const int expect_nbr_ls, const bool write_binormals)
{
    Eigen::VectorXd lsv0;
    int nbr_ls0;
    std::vector<std::vector<Eigen::Vector3d>> lines;

    nbr_ls0 = expect_nbr_ls;
    get_level_set_sample_values(ls, nbr_ls0, lsv0);
    std::vector<Eigen::Vector3d> verlist;
    verlist.reserve(nbr_ls0 * F.rows());

    for (int i = 0; i < nbr_ls0; i++)
    {
        std::vector<std::vector<Eigen::Vector3d>> polylines;
        double value = lsv0[i];
        std::vector<bool> left_large;
        get_iso_lines(lsmesh, loop, V, F, ls, value, polylines, left_large);
        for (int j = 0; j < polylines.size(); j++)
        {
            if(polylines[j].size()>3){// remove too short elements
                lines.push_back(polylines[j]);
            }
            
        }
    }
    std::cout << "Writing the Polylines, please provide the file prefix" << std::endl;
    std::string fname = igl::file_dialog_save();
    write_polyline_xyz(lines, fname);
    if (!write_binormals)
    {
        return;
    }
    std::cout << "Reading Bi-normals of triangle mesh. If you don't have reference binormal vector file, esc here" << std::endl;
    for (int i = 0; i < lines.size(); i++)
    {
        for (int j = 0; j < lines[i].size(); j++)
        {
            verlist.push_back(lines[i][j]);
        }
    }
    Eigen::MatrixXd bn;
    bool bireaded = read_bi_normals(bn);
    if(!bireaded)
    {
        std::vector<std::vector<Eigen::Vector3d>> binormals;
        generate_estimate_binormals_out_of_polylines(lines, binormals);
        write_polyline_xyz(binormals, fname + "_b");
        std::cout << "Binormals got written" << std::endl;
        return;
    }
    std::cout << "Binormals readed" << std::endl;
    Eigen::MatrixXd C;
    Eigen::VectorXi I;
    Eigen::VectorXd D;
    int nq = verlist.size();
    Eigen::MatrixXd B;
    B.resize(nq, 3);
    Eigen::MatrixXd vers = vec_list_to_matrix(verlist);
    igl::point_mesh_squared_distance(vers, V, F, D, I, C);
    for (int i = 0; i < nq; i++)
    {
        int fid = I(i);
        int v0 = F(fid, 0);
        int v1 = F(fid, 1);
        int v2 = F(fid, 2);
        Eigen::Vector3d pt = verlist[i];
        Eigen::Vector3d ver0 = V.row(v0);
        Eigen::Vector3d ver1 = V.row(v1);
        Eigen::Vector3d ver2 = V.row(v2);
        std::array<double, 3> coor = barycenter_coordinate(ver0, ver1, ver2, pt);

        Eigen::Vector3d dir = coor[0] * bn.row(v0) + coor[1] * bn.row(v1) + coor[2] * bn.row(v2);
        dir = dir.normalized();
        if(dir == Eigen::Vector3d(0,0,0)){
            dir = Eigen::Vector3d(1,0,0);
        }
        B.row(i) = dir;
    }
    std::vector<std::vector<Eigen::Vector3d>> bilist = lines;
    int counter = 0;
    Eigen::Vector3d bbase = B.row(0);
    for (int i = 0; i < lines.size(); i++)
    {
        Eigen::Vector3d base = B.row(counter);
        if(base.dot(bbase)<0){
            base *= -1;
        }

        for (int j = 0; j < lines[i].size(); j++)
        {
            Eigen::Vector3d direction = orient_vector(base, B.row(counter));
            base = direction;
            bilist[i][j] = direction;
            counter++;
        }
    }
    write_polyline_xyz(bilist, fname+"_b");
    std::cout << "Binormals got written" << std::endl;
}

void lsTools::show_binormals(const Eigen::VectorXd &func, Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, Eigen::MatrixXd& binormals,  double ratio){
    
    E0 = V - Binormals * ratio;
    E1 = V + Binormals * ratio;
}

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

void lsTools::show_traced_binormals(Eigen::MatrixXd &bE0, Eigen::MatrixXd &bE1, Eigen::MatrixXd &nE0, Eigen::MatrixXd &nE1, const double scale){
    int nbr = 0; // the number of bi-normals
    for (int i = 0; i < trace_vers.size(); i++)
    {
        if (trace_vers[i].size() > 2)
            nbr += trace_vers[i].size() - 2;
    }
    bE0.resize(nbr, 3);
    bE1.resize(nbr, 3);
    nE0.resize(nbr, 3);
    nE1.resize(nbr, 3);
    int counter = 0;
    for (int i = 0; i < trace_vers.size(); i++) // vi, vi+1, vi+2
    {
        for (int j = 0; j < trace_vers[i].size() - 2; j++)
        {
            Eigen::Vector3d v0 = trace_vers[i][j];
            Eigen::Vector3d v1 = trace_vers[i][j + 1];
            Eigen::Vector3d v2 = trace_vers[i][j + 2];
            Eigen::Vector3d dir1 = v1 - v0;
            Eigen::Vector3d dir2 = v2 - v1;
            Eigen::Vector3d binormal = dir1.cross(dir2).normalized();
            if (counter > bE0.rows() - 1)
            {
                std::cout<<"size wrong with trace vers"<<std::endl;
            }
            bE0.row(counter) = v1 + binormal * scale;
            bE1.row(counter) = v1 - binormal * scale;
            CGMesh::HalfedgeHandle heh = trace_hehs[i][j + 1];
            int eid = lsmesh.edge_handle(heh).idx();
            Eigen::Vector3d enormal = norm_e.row(eid);
            nE0.row(counter) = v1;
            nE1.row(counter) = v1 + enormal * scale;
            counter++;
        }
    }
}
void lsTools::show_max_pg_energy(Eigen::VectorXd &e)
{
    int ninner = analizers[0].LocalActInner.size();
    if(ninner==0){
        return;
    }
    int vnbr = V.rows();
    int rep = PGE.size() / ninner;
    // Eigen::VectorXd result = Eigen::VectorXd::Zero(vnbr);
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(vnbr, rep);
    for (int i = 0; i < ninner; i++)
    {
        int vm = IVids[i];
        for (int j = 0; j < rep; j++)
        {
            result(vm, j) = PGE[i + j * ninner];
        }
    }
    e = result.rowwise().norm();
}
void lsTools::show_max_pg_energy_all(Eigen::MatrixXd &energy)
{
    int ninner = analizers[0].LocalActInner.size();
    if(ninner==0){
        std::cout<<"analizer is empty"<<std::endl;
        return;
    }
    int vnbr = V.rows();
    int rep = PGE.size() / ninner;
    Eigen::MatrixXd result = Eigen::MatrixXd::Zero(vnbr, rep);
    for (int i = 0; i < ninner; i++)
    {
        int vm = IVids[i];
        for (int j = 0; j < rep; j++)
        {
            result(vm, j) = PGE[i + j * ninner];
        }
    }
    energy = result;
}
// ref is the reference vertices marked in the vertex list of the mesh
void get_one_ring_vertices_simple(const Eigen::VectorXi &ref, CGMesh &lsmesh, Eigen::VectorXi &ring1)
{
    int vsize = ref.size();
    ring1 = ref;
    // std::cout<<"inside 1"<<std::endl;
    for (int i = 0; i < vsize; i++)
    {
        if (ref[i] < 0)
        {
            continue;
        }
        int vid = i;
        std::vector<int> pts;
        // std::cout<<"inside 2"<<std::endl;
        get_one_ring_vertices(lsmesh, vid, pts);
        // std::cout<<"inside 3"<<std::endl;
        for (int p : pts)
        {
            if(p<0||p>=ref.size()){
                std::cout<<"p out of boundary, "<<p<<std::endl;
            }
            ring1[p] = 1;
        }
        // std::cout<<"inside 4"<<std::endl;
    }
}

// rotate the z axis to make it align with the earth axis (pointing to the north)
void rotate_z_axis_to_earth_axis(const Eigen::MatrixXd& Vori, Eigen::MatrixXd& Vnew, const double latitude_degree){
    double latitude_radian = latitude_degree * LSC_PI / 180.;
    double rho = (LSC_PI / 2 - latitude_radian); // this is correct rotation
    Eigen::Matrix3d rotation;
    rotation << 1, 0, 0,
        0, cos(rho), -sin(rho),
        0, sin(rho), cos(rho);
    // rotation * Vold = Vnew
    Eigen::MatrixXd Vt = Vori.transpose(); // 3 x n
    Vt = rotation * Vt; // 3 x n
    Vnew = Vt.transpose();
}
Eigen::Vector3d get_light_rotated_back_from_earth_axis(const double latitude_degree, const double theta, const double phi)
{
    double latitude_radian = latitude_degree * LSC_PI / 180.;
    double rho = -(LSC_PI / 2 - latitude_radian); // this is correct rotation
    Eigen::Matrix3d rotation;
    rotation << 1, 0, 0,
        0, cos(rho), -sin(rho),
        0, sin(rho), cos(rho);
    Eigen::Vector3d light = angle_ray_converter(theta, phi);
    std::cout<<"after rotation, "<<light<<std::endl;
    light = rotation * light;
    return light;
}
void mesh_unit_scale(const Eigen::MatrixXd &V, Eigen::MatrixXd &Vout)
{
    double xmin = V.col(0).minCoeff();
    double xmax = V.col(0).maxCoeff();

    double ymin = V.col(1).minCoeff();
    double ymax = V.col(1).maxCoeff();

    double zmin = V.col(2).minCoeff();
    double zmax = V.col(2).maxCoeff();

    Eigen::Vector3d Vmin = Eigen::Vector3d(xmin, ymin, zmin);
    Eigen::Vector3d Vmax = Eigen::Vector3d(xmax, ymax, zmax);
    Eigen::Vector3d Vcent = (Vmax + Vmin) / 2;
    double diagnal = (Vmax - Vmin).norm();
    double ratio = 1. / diagnal;
    Vout = V;
    for (int i = 0; i < V.rows(); i++)
    {
        Vout.row(i) -= Vcent;
    }
    Vout *= ratio;
    std::cout<<"Unit Scale the Mesh, the scaling ratio, "<<ratio<<", the translation, "<<Vcent.norm()<<std::endl;
}

template <typename Tn>
std::vector<int> sort_indices(const std::vector<Tn> &v)
{
    std::vector<int> idx(v.size());
    std::iota(idx.begin(), idx.end(), 0);
    std::sort(idx.begin(), idx.end(), [&v](int i1, int i2)
              { return v[i1] > v[i2]; });
    return idx;
}
// // mark the high energy vertices as 0 in the "he" vector
// void mark_high_energy_vers(const Eigen::VectorXd &energy, const int ninner, const double percentage,
//                            const std::vector<int> &IVids, Eigen::VectorXi &he, std::vector<int>& refid)
// {
//     he.resize(ninner);
//     refid.clear();
//     // first get all the energy on each ver
//     int rep = energy.size() / ninner;
//     Eigen::MatrixXd emat = Eigen::MatrixXd::Zero(ninner, rep);
//     std::vector<double> evec(ninner);
//     for (int i = 0; i < ninner; i++)
//     {
//         for (int j = 0; j < rep; j++)
//         {
//             emat(i, j) = energy[i + j * ninner];
//         }
//         evec[i] = emat.row(i).norm();
//     }
//     std::vector<int> order = sort_indices(evec);
//     std::cout << "first and last of sorting " << evec[order.front()] << ", " << evec[order.back()] << std::endl;
//     he = Eigen::VectorXi::Ones(ninner);
//     int npick = ninner * percentage / 100;
//     std::cout << "Marking Points\n";
//     for (int i = 0; i < npick; i++)
//     {
//         he[order[i]] = 0;
//         int vm = IVids[order[i]];
//         std::cout << "id, " << vm << ", energy, " << evec[order[i]] << std::endl;
//         refid.push_back(vm);
//     }

// }


// real length means we don't do unit scale.
CGMesh polyline_to_strip_mesh(const std::vector<std::vector<Eigen::Vector3d>> &ply, const std::vector<std::vector<Eigen::Vector3d>> &bi, 
const double ratio, const double ratio_back, const bool &realLength)
{
    CGMesh mesh;
    if (ply.empty())
    {
        std::cout << "Please load the polylines" << std::endl;
        return mesh;
    }

    int vnbr = 0;
    for (auto line : ply)
    {
        for (auto ver : line)
        {
            vnbr++;
        }
    }
    std::vector<std::vector<CGMesh::VertexHandle>> vhd_ori, vhd_off;
    vhd_ori.reserve(ply.size());
    vhd_off.reserve(ply.size());
    for (int i = 0; i < ply.size(); i++)
    {
        std::vector<CGMesh::VertexHandle> vhd0, vhd1;
        for (int j = 0; j < ply[i].size(); j++)
        {
            int idf = std::max(j - 1, 0);
            int idb = std::min(j + 1, int(ply[i].size() - 1));
            Eigen::Vector3d tangent = (ply[i][idf] - ply[i][idb]).normalized();
            Eigen::Vector3d direction = bi[i][j];
            if (!realLength)
            {
                direction.normalize();
                if (abs(direction.dot(tangent)) > 1e-3)
                { // not orthogonal, need to get the correct ratio to make the norm is 1
                    double vt = direction.dot(tangent);
                    double alpha = sqrt(1 / (1 - vt * vt));
                    direction *= alpha;
                }
            }

            Eigen::Vector3d p0 = ply[i][j] - direction * ratio_back;
            Eigen::Vector3d p1 = ply[i][j] + direction * ratio;
            vhd0.push_back(mesh.add_vertex(CGMesh::Point(p0[0], p0[1], p0[2])));
            vhd1.push_back(mesh.add_vertex(CGMesh::Point(p1[0], p1[1], p1[2])));
        }
        vhd_ori.push_back(vhd0);
        vhd_off.push_back(vhd1);
    }
    for (int i = 0; i < ply.size(); i++)
    {
        for (int j = 0; j < ply[i].size() - 1; j++)
        {
            std::vector<CGMesh::VertexHandle> fhds;
            fhds.push_back(vhd_ori[i][j]);
            fhds.push_back(vhd_ori[i][j + 1]);
            fhds.push_back(vhd_off[i][j + 1]);
            fhds.push_back(vhd_off[i][j]);
            mesh.add_face(fhds);
        }
    }
    return mesh;
}


void read_plylines_and_binormals(std::vector<std::vector<Eigen::Vector3d>>& ply, std::vector<std::vector<Eigen::Vector3d>>& bin){
    ply.clear();
    bin.clear();
    std::cout<<"Reading the binormal files and polyline files. Please type down the prefix"<<std::endl;
    std::string fname = igl::file_dialog_save();
    if (fname.length() == 0)
    {
        std::cout<<"Please type something "<<std::endl;
        return;
    }
    std::vector<std::vector<double>> vx, vy, vz, nx, ny, nz;
    read_csv_data_lbl(fname+"_x.csv", vx);
    read_csv_data_lbl(fname+"_y.csv", vy);
    read_csv_data_lbl(fname+"_z.csv", vz);
    read_csv_data_lbl(fname+"_b_x.csv", nx);
    read_csv_data_lbl(fname+"_b_y.csv", ny);
    read_csv_data_lbl(fname+"_b_z.csv", nz);

    ply = xyz_to_vec_list(vx, vy, vz);
    bin = xyz_to_vec_list(nx, ny, nz);

}

void read_origami_and_convert_to_polylines(std::vector<std::vector<Eigen::Vector3d>>& ply, std::vector<std::vector<Eigen::Vector3d>>& bin){
    ply.clear();
    bin.clear();
    // get the quad mesh
    std::cout << "\nreading quad mesh file with attached binormals" << std::endl;
    std::string fname = igl::file_dialog_open();
    if (fname.length() == 0)
    {
        std::cout << "reading failed" << std::endl;
        return;
    }
    Eigen::MatrixXd Vqd;
    Eigen::MatrixXi Fqd;
    igl::readOBJ(fname, Vqd, Fqd);
    std::cout<<"Mesh readed, nbr cols, "<<Vqd.cols()<<std::endl<<std::endl;
    if(Vqd.cols()!=6){
        std::cout<<"Please read the quad mesh attached with binormal vector info"<<std::endl;
        return;
    }
    Eigen::MatrixXd V = Vqd.leftCols(3);
    Eigen::MatrixXd bn = Vqd.rightCols(3);
    // std::cout<<"V"<<V<<std::endl;
    // std::cout<<"B"<<bn<<std::endl;
    // get the binormals 

    std::cout<<"reading the row-info"<<std::endl;
    fname = igl::file_dialog_open();
    if (fname.length() == 0)
    {
        std::cout << "reading failed" << std::endl;
        return;
    }
    std::vector<std::vector<double>> row;
    read_csv_data_lbl(fname, row);
    for(int i = 0;i<row.size();i++){
        std::vector<Eigen::Vector3d> line, binormal;
        for(int j=0;j<row[i].size();j++){
            int id = row[i][j];
            if(id<0){
                break;
            }
            Eigen::Vector3d vertex = V.row(id), vbin = bn.row(id);
            line.push_back(vertex);
            binormal.push_back(vbin);
        }
        ply.push_back(line);
        bin.push_back(binormal);
    }
    std::cout<<"readed polyline: the nbr, "<<ply.size()<<std::endl;
}

void update_qd_mesh_with_plylines(){
    std::vector<std::vector<Eigen::Vector3d>> ply, bin;
    std::cout << "\nreading quad mesh file with attached binormals" << std::endl;
    std::string fname = igl::file_dialog_open();
    if (fname.length() == 0)
    {
        std::cout << "reading failed" << std::endl;
        return;
    }
    Eigen::MatrixXd Vqd;
    Eigen::MatrixXi Fqd;
    igl::readOBJ(fname, Vqd, Fqd);
    std::cout << "Mesh readed, nbr cols, " << Vqd.cols() << std::endl
              << std::endl;
    if (Vqd.cols() != 6)
    {
        std::cout << "Please read the quad mesh attached with binormal vector info" << std::endl;
        return;
    }
    read_plylines_and_binormals(ply, bin);

    int counter = 0;
    for(int i=0;i<ply.size();i++){
        for(int j=0;j<ply[i].size();j++){
            Eigen::Vector3d vertex = ply[i][j];
            Eigen::Vector3d binormal = bin[i][j];
            Vqd.row(counter) << vertex[0], vertex[1], vertex[2], binormal[0], binormal[1], binormal[2];
            counter ++;
        }
    }
    std::cout<<"now saving the updated quad mesh "<<std::endl;
    std::ofstream file;
    fname = igl::file_dialog_save();
    file.open(fname);
    for (int i = 0; i < counter; i++)
    {
        file << "v " << Vqd(i, 0) << " " << Vqd(i, 1) << " " << Vqd(i, 2) << " " << Vqd(i, 3) << " " << Vqd(i, 4) << " " << Vqd(i, 5) << std::endl;
    }
    for (int i = 0; i < Fqd.rows(); i++)
    {
        file << "f " << Fqd(i, 0) + 1 << " " << Fqd(i, 1) + 1 << " " << Fqd(i, 2) + 1 << " " << Fqd(i, 3) + 1 << std::endl;
    }
    file.close();
}
int closest_point_id_in_F(const Eigen::Vector3d& query, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const int fid){
    std::array<int,3> ids;
    ids[0] = F(fid, 0);
    ids[1] = F(fid, 1);
    ids[2] = F(fid, 2);
    double mindis = 1000;
    int minid = -1;
    for(int i=0;i<3;i++){
        Eigen::Vector3d ver = V.row(ids[i]);
        double dis = (ver-query).norm();
        if(dis<mindis){
            mindis = dis;
            minid = ids[i];
        }
    }
    if (minid == -1 || mindis > 1e-6)
    {
        std::cout<<"Error in closest_point_id_in_F, minid, "<<minid<<" mindis, "<<mindis<<std::endl;
    }
    return minid;
}

int the_nbr_of_ones(const Eigen::VectorXi &vec)
{
    int nbr = 0;
    for (int i = 0; i < vec.size(); i++)
    {
        if (vec[i] == 1)
        {
            nbr++;
        }
    }
    return nbr;
}
// classify the vertices of the base into 3 types: the patch which corresponds to ref, the transition part, and the rest of the mesh.
void project_mesh_and_get_shading_info(CGMesh &ref, CGMesh &base, const int nbr_rings, Eigen::VectorXi &info,
                                       Eigen::MatrixXd &P1, Eigen::MatrixXd &P2)
{
    lsTools reftool(ref);
    lsTools bastool(base);
    reftool.initialize_mesh_properties();
    bastool.initialize_mesh_properties();
    int rvnbr = reftool.V.rows();
    int bvnbr = bastool.V.rows();
    int ninner = bastool.IVids.size();
    Eigen::VectorXi InnerV = bastool.InnerV;
    Eigen::MatrixXd C;
    Eigen::VectorXi I;
    Eigen::VectorXd D;
    Eigen::VectorXi mapping = Eigen::VectorXi::Ones(bvnbr) * -1;  // the list (size is the same as the vertex list) shows the type of the conditions
    Eigen::VectorXi neighbours = Eigen::VectorXi::Ones(bvnbr) * -1; // the transition points
    std::vector<int> bnd;// boundary of the first type
    igl::point_mesh_squared_distance(reftool.V, bastool.V, bastool.F, D, I, C);
    std::vector<int> removed; // mark all the vertices of the reference
    std::vector<Eigen::Vector3d> points1, points2;
    int nbr_bnd = 0;
    for (int i = 0; i < rvnbr; i++)
    {
        int fid = I(i);
        int v0 = bastool.F(fid, 0);
        int v1 = bastool.F(fid, 1);
        int v2 = bastool.F(fid, 2);
        // the vertex id on base mesh
        int vid = closest_point_id_in_F(reftool.V.row(i), bastool.V, bastool.F, fid);
        if(vid<0){
            std::cout<<"did not find the closest one!"<<std::endl;
        }
        mapping[vid] = 1; // the point is special
        removed.push_back(vid);
        CGMesh::VertexHandle vh = ref.vertex_handle(i);
        if(ref.is_boundary(vh)){ // this is a boundary vertex.
            mapping[vid] = 2;
            neighbours[vid] = 1;
            nbr_bnd++;
        }
    }
    // std::cout<<"nbr bnd, "<<nbr_bnd<<std::endl;
    info = Eigen::VectorXi::Ones(ninner) * -1;
    for (int i = 0; i < nbr_rings; i++)
    {
        Eigen::VectorXi ring;
        // std::cout<<"check1"<<std::endl;
        get_one_ring_vertices_simple(neighbours, base, ring);
        // std::cout<<"check2"<<std::endl;
        // remove innver vers
        for (int j = 0; j < removed.size(); j++)
        {
            ring[removed[j]] = -1;
        }
        // std::cout<<"nbr elements 0, "<<the_nbr_of_ones(neighbours)<<std::endl;
        // std::cout<<"nbr elements 1, "<<the_nbr_of_ones(ring)<<std::endl;
        neighbours = ring;
    }
    for (int i = 0; i < neighbours.size(); i++)
    {
        int location = -1;
        // if (neighbours[i] < 0 && mapping[i] < 0)
        // {
        //     continue;
        // }
        location = InnerV[i];
        if (location < 0)// boundary point of base mesh, ignor
        {
            continue;
        }
        info[location] = 0;// ordinary points
        if(neighbours[i]>0){// this is a transition point
            info[location] = 2;
            points2.push_back(bastool.V.row(i));
        }
        if(mapping[i]>0){// this is a special condition point, corresponding to the reference
            info[location] = 1;
            points1.push_back(bastool.V.row(i));
        }
        if(info[location]<0){ // this is a ordinary point
            info[location] = 0;
        }
        
    }
    P1 = vec_list_to_matrix(points1);
    P2 = vec_list_to_matrix(points2);
}

// -1 means transition, 0 ... , n are types.
void project_mesh_and_get_vertex_class(std::vector<CGMesh> &ref, CGMesh &base, Eigen::VectorXi &info)
{
    int nbr_types = ref.size();
    std::vector<Eigen::VectorXi> infos;
    infos.resize(nbr_types);
    Eigen::MatrixXd P1, P2;
    for (int i = 0; i < nbr_types; i++)
    {
        project_mesh_and_get_shading_info(ref[i], base, 0, infos[i],
                                          P1, P2);
    }
    // std::cout<<"check 1"<<std::endl;
    // get the infos into one vector
    lsTools bastool(base);
    // bastool.initialize_mesh_properties();
    // std::cout<<"check 2"<<std::endl;
    int ninner = bastool.IVids.size();
    info = Eigen::VectorXi::Ones(ninner) * -1; // set all the points as -1.
    for (int i = 0; i < ninner; i++)
    {
        int type = -1;
        for (int j = 0; j < nbr_types; j++)
        {
            if (infos[j][i] == 1)
            {
                type = j;
                // break;
            }
        }
        info[i] = type;
    }
    // std::cout<<"check 3"<<std::endl;

}

// put the matrix in the middle of a matrix which is 3 times larger.
// mat should be a symmetric matrix
spMat put_mat_in_middle(const spMat &mat, const int sizemat)
{
    spMat tmp1(sizemat, sizemat * 3);
    tmp1.middleCols(sizemat, sizemat) = mat;
    spMat tmp2(sizemat * 3, sizemat * 3);
    tmp2.middleCols(sizemat, sizemat) = tmp1.transpose();
    return tmp2;
}
Eigen::VectorXd put_vec_in_middle(const Eigen::VectorXd &vec)
{
    int sizevec = vec.size();
    Eigen::VectorXd result = Eigen::VectorXd::Zero(sizevec * 3);
    result.segment(sizevec, sizevec) = vec;
    return result;
}
Eigen::VectorXd duplicate_valus(const double value, const int n){
    Eigen::VectorXd result;
    result.resize(n);
    for(int i=0;i<n;i++){
        result[i]=value;
    }
    return result;
}
Eigen::MatrixXd duplicate_vector(const Eigen::Vector2d& vec, const int n){
    Eigen::MatrixXd result;
    result.resize(n,2);
    for(int i=0;i<n;i++){
        result.row(i)=vec;
    }
    return result;
}
void lsTools::save_traced_curves(const std::string filename){
    std::ofstream file;

    for (int i = 0; i < trace_vers.size(); i++)
    {
        file.open(filename + "_" + std::to_string(i) + ".csv");
        for(int j=0;j<trace_vers[i].size();j++){
            // the vertex coordinates and the halfedge handle id
            Eigen::Vector3d ver = trace_vers[i][j];
            CGMesh::HalfedgeHandle hdl = trace_hehs[i][j];
            CGMesh::EdgeHandle edge = lsmesh.edge_handle(hdl);
            int id = hdl.idx();
            file<<ver[0]<<","<<ver[1]<<","<<ver[2]<<","<<id<<std::endl;
        }
        file.close();
    }
    
}
void lsTools::load_one_traced_curve(){
    std::string fname = igl::file_dialog_open();
    std::vector<Eigen::Vector3d> pts;
    std::vector<CGMesh::HalfedgeHandle> hds;
    std::vector<std::vector<double>> data;
    std::vector<Eigen::Vector2d> tparas;
    std::ofstream file;
    bool readed = read_csv_data_lbl(fname, data);
    if(data.empty()){
        std::cout<<"load traced curves failed."<<std::endl;
        return;
    }
    for (int i = 0; i < data.size(); i++)
    {
        Eigen::Vector3d ver(data[i][0], data[i][1], data[i][2]);
        CGMesh::HalfedgeHandle halfedge = lsmesh.halfedge_handle(data[i][3]);
        pts.push_back(ver);
        hds.push_back(halfedge);
        int id0 = lsmesh.from_vertex_handle(halfedge).idx();
        int id1 = lsmesh.to_vertex_handle(halfedge).idx();
        Eigen::Vector3d p0 = V.row(id0);
        Eigen::Vector3d p1 = V.row(id1);
        double t = get_t_of_segment(ver, p0, p1);
        if(t<0||t>1){
            std::cout<<"inaccurate t!!!!, "<<t<<std::endl;
        }
        Eigen::Vector2d vpara = paras.row(id0) * (1 - t) + paras.row(id1) * t;
        tparas.push_back(vpara);
    }

    trace_hehs.push_back(hds);
    trace_vers.push_back(pts);
    traced_paras.push_back(tparas);
    if (assigned_trace_ls.empty())
    {
        assigned_trace_ls.push_back(0);
    }
    else{
        assigned_trace_ls.push_back(assigned_trace_ls.back() + 1);
    }
    if (trace_vers.size() >= 2)
    {
        estimate_strip_width_according_to_tracing();
        std::cout<<"The estimated strip width: "<<strip_width<<std::endl;
    }

    std::cout << "loaded one curve succeed. currently there are " << trace_hehs.size() << " curves" << std::endl;
}
void lsTools::clear_traced_curves(){
    trace_hehs.clear();
    trace_vers.clear();
    assigned_trace_ls.clear();
    traced_paras.clear();
}
CGMesh lsTools::write_parameterization_mesh(){
    CGMesh mesh = lsmesh;
    int nv = mesh.n_vertices();
	for (CGMesh::VertexIter v_it = mesh.vertices_begin(); v_it != mesh.vertices_end(); ++v_it)
	{
        int vid=v_it.handle().idx();
        mesh.point(*v_it) = CGMesh::Point(paras(vid, 0), paras(vid, 1), 0);
    }
    return mesh;
}
void lsTools::show_traced_curve_params(Eigen::MatrixXd &curve){
    int ncurves = traced_paras.size();
    int count = 0;
    std::vector<Eigen::Vector3d> pts;
    for(int i=0;i<ncurves;i++){
        for(int j=0;j<traced_paras[i].size();j++){
            auto pr = traced_paras[i][j];
            pts.push_back(Eigen::Vector3d(pr[0], pr[1], 0));
        }
    }
    curve = vec_list_to_matrix(pts);
}
void read_pts_csv_and_write_xyz_files(){
    std::vector<std::vector<Eigen::Vector3d>> ply, bin;
    std::cout<<"Reading ply file"<<std::endl;
    std::string fname = igl::file_dialog_open();
    if (fname.length() == 0)
    {
        std::cout<<"Read Error "<<std::endl;
        return;
    }
    std::vector<std::vector<double>> data;
    read_csv_data_lbl(fname, data);

    std::vector<Eigen::Vector3d> pts;
    pts.resize(data.size());
    std::vector<Eigen::Vector3d> binormals;
    binormals.resize(data.size());
    for (int i = 0; i < data.size(); i++)
    {
        pts[i][0] = data[i][0];
        pts[i][1] = data[i][1];
        pts[i][2] = data[i][2];
    }
    for (int i = 1; i < data.size() - 1; i++)
    {
        Eigen::Vector3d v0 = pts[i - 1];
        Eigen::Vector3d v1 = pts[i];
        Eigen::Vector3d v2 = pts[i + 1];
        Eigen::Vector3d blocal = (v0 - v1).cross(v1 - v2);
        blocal.normalize();
        if(i>1)
        {
            Eigen::Vector3d biprev = binormals[i - 1];
            if (blocal.dot(biprev) < 0)
            {
                blocal *= -1;
            }
        }
        binormals[i] = blocal;
    }
    binormals[0] = binormals[1];
    binormals.back() = binormals[data.size() - 2];
    ply.push_back(pts);
    bin.push_back(binormals);
    std::cout << "Saving the poly files" << std::endl;
    fname = igl::file_dialog_save();
    if (fname.length() == 0)
    {
        std::cout<<"Please type something "<<std::endl;
        return;
    }

    write_polyline_xyz(ply, fname);
    write_polyline_xyz(bin, fname + "_b");
    std::cout << "files get saved" << std::endl;

}

void recover_end_pt_direction_from_crease(const Eigen::Vector3d &crease, const Eigen::Vector3d &v01, Eigen::Vector3d &bi)
{
    Eigen::Vector3d np = crease.cross(v01).normalized();
    bi = v01.cross(np).normalized();
    if (bi.dot(crease) < 0)
    {
        bi *= -1;
    }
    return;
}

void recover_polyline_endpts(){
    std::vector<std::vector<Eigen::Vector3d>> ply0, ply1, bin0, bin1, ply2, bin2;
    std::cout<<"\nFirst read the plyline with creases"<<std::endl;
    read_plylines_and_binormals(ply0, bin0);
    std::cout<<"second read the plyline with binormals"<<std::endl;
    read_plylines_and_binormals(ply1, bin1);
    if(ply0.empty() || ply1.empty()){
        std::cout<<"Please type down correct prefixes for the files"<<std::endl;
        return;
    }
    int ncurve = ply0.size();
    if(ncurve != ply1.size()){
        std::cout<<"The sizes of the curves not match: "<<ncurve<<", "<<ply1.size()<<std::endl;
    }
    int crease_cid = 0;
    for (int i = 0; i < ncurve; i++)
    {
        int vnbr = ply1[i].size();
        if (vnbr < 4)
        {
            ply2.push_back(ply1[i]);
            bin2.push_back(bin1[i]);
            crease_cid++;
            continue;
        }

        std::vector<Eigen::Vector3d> tmp_ply(vnbr), tmp_bin(vnbr), crs_ply = ply0[i], recover_ply(vnbr);
        int vnbr_crease = ply0[crease_cid].size();
        std::vector<Eigen::Vector3d> current_crease = bin0[crease_cid];
        if (vnbr - vnbr_crease != 2)
        {
            std::cout << "ERROR in recover_polyline_endpts(): polyline size not match: "<< vnbr - vnbr_crease<< std::endl;
            return;
        }
        // recover the first pt
        tmp_ply = ply1[i];

        Eigen::Vector3d v01 = tmp_ply[0]-tmp_ply[1];
        Eigen::Vector3d crease = current_crease[0];
        recover_end_pt_direction_from_crease(crease, v01, tmp_bin[0]);
        recover_ply.front() = v01 + crs_ply.front();

        // recover the last pt
        v01 = tmp_ply[vnbr - 1] - tmp_ply[vnbr - 2];
        crease = current_crease.back();
        recover_end_pt_direction_from_crease(crease, v01, tmp_bin.back());
        recover_ply.back() = v01 + crs_ply.back();

        for (int j = 1; j < vnbr-1; j++)
        {
            tmp_bin[j] = current_crease[j - 1];
            recover_ply[j] = crs_ply[j - 1];
        }
        ply2.push_back(tmp_ply);
        bin2.push_back(tmp_bin);

        crease_cid++;
    }
    std::cout<<"Save the polylines with recovered endpts"<<std::endl;
    std::string fname = igl::file_dialog_save();
    if (fname.length() == 0)
    {
        std::cout<<"Please type something "<<std::endl;
        return;
    }

    write_polyline_xyz(ply2, fname);
    write_polyline_xyz(bin2, fname + "_b");
    std::cout << "files get saved" << std::endl;
}

// the nbr of length = the nbr of angles + 1.
void write_unfold_single_strip(int which_curve)
{
    std::vector<std::vector<Eigen::Vector3d>> polylines, creases;
    std::cout<<"Read the crease plyline"<<std::endl;
    read_plylines_and_binormals(polylines, creases);
    std::vector<Eigen::Vector3d> ply, crs;
    ply = polylines[which_curve];
    crs = creases[which_curve];
    int nver = ply.size();
    std::vector<double> lengths(nver - 1);
    std::vector<double> angles(nver - 2);
    std::vector<double> alefts(nver - 2);
    for(int i = 0; i < nver - 2; i++){
        Eigen::Vector3d c = crs[i + 1].normalized();
        Eigen::Vector3d poT = (ply[i + 2] - ply[i + 1]).normalized();
        Eigen::Vector3d neT = (ply[i] - ply[i + 1]).normalized();
        double ang_left = acos(c.dot(neT));
        double ang_right = acos(c.dot(poT));
        angles[i] = ang_left + ang_right;
        lengths[i] = (ply[i] - ply[i + 1]).norm();
        std::cout<<i<<", angle "<<angles[i]<<", length, "<<lengths[i]<<std::endl;
        alefts[i] = ang_left;
    }
    lengths.back() = (ply[nver - 1] - ply[nver - 2]).norm();
    std::cout<<"length back "<<lengths.back()<<std::endl;
    // set the initial direction as (0, 0) ->(0, 1)
    int vnbr = lengths.size() + 1;
    Eigen::MatrixXd pts;
    double total_length = 0;
    pts.resize(vnbr, 2);
    pts.row(0) = Eigen::Vector2d(0, 0);
    pts.row(1) = Eigen::Vector2d(0, lengths[0]);
    total_length += lengths[0];
    for (int i = 2; i < vnbr; i++)
    {
        // solve pts[i]
        // the reference direction from i-2 to i-1,
        Eigen::Vector2d ref_dirc = (pts.row(i - 1) - pts.row(i - 2)).normalized();
        double ang_diff = LSC_PI - angles[i - 2]; // counter clock rotation
        double len = lengths[i - 1];
        total_length+=len;
        double a = cos(ang_diff);
        double b = sin(ang_diff);
        Eigen::Vector2d ref_otho(-ref_dirc[1], ref_dirc[0]); // the othogonal vector is (-y, x)
        // new direction is a * ref_dirc + b * ref_otho.
        Eigen::Vector2d direction = a * ref_dirc + b * ref_otho;
        Eigen::Vector2d newpt = Eigen::Vector2d(pts.row(i - 1)) + direction * len;
        pts.row(i) = newpt;
    }
    Eigen::VectorXd iden = Eigen::VectorXd::Ones(vnbr);
    Eigen::Vector2d centroid = pts.transpose() * iden;
    centroid /= vnbr;
    Eigen::MatrixXd y = pts.rowwise() - centroid.transpose();// y is a matrix of n x 2
    Eigen::Matrix2d YTY = y.transpose() * y;// 2 x 2
    Eigen::EigenSolver<Eigen::MatrixXd> eg(YTY);
    Eigen::VectorXcd ev = eg.eigenvalues();
    int minevid = ev[0].real() < ev[1].real() ? 0 : 1;
    Eigen::VectorXcd eigenvecC = eg.eigenvectors().col(minevid);
    Eigen::Vector2d eigenvec(eigenvecC[0].real(), eigenvecC[1].real());
    Eigen::Vector2d n = eigenvec.normalized();// this is othogonal to the regression line.

    Eigen::VectorXd proj_dis = (y * n).rowwise().norm();
    std::cout<<"The maximal derivate: "<<proj_dis.maxCoeff()<<", the total length of this polyline: "<<total_length<<std::endl;
    std::vector<Eigen::Vector2d> directions;
    directions.resize(vnbr);
    
    Eigen::Vector2d tangent = (pts.row(1) - pts.row(0)).normalized();
    Eigen::Vector2d orthogonal = Eigen::Vector2d(-tangent[1], tangent[0]);
    directions[0] = orthogonal;

    tangent = (pts.row(vnbr - 1) - pts.row(vnbr - 2)).normalized();
    orthogonal = Eigen::Vector2d(-tangent[1], tangent[0]);
    directions.back() = orthogonal;

    for (int i = 1; i < vnbr - 1; i++)
    {
        // the tangent direction
        tangent = (pts.row(i) - pts.row(i - 1)).normalized();
        orthogonal = Eigen::Vector2d(-tangent[1], tangent[0]);
        double ang = LSC_PI - alefts[i - 1];

        directions[i] = tangent * cos(ang) + orthogonal * sin(ang);
        directions[i] *= 1 / sin(ang); // make sure the direction project on orthogonal is 1
    }
    std::vector<std::vector<Eigen::Vector3d>> ply_out, bin_out;
    ply_out.resize(1);
    bin_out.resize(1);
    ply_out[0].resize(vnbr);
    bin_out[0].resize(vnbr);
    for(int i = 0;i<vnbr;i++){
        ply_out[0][i] = Eigen::Vector3d(pts(i,0), pts(i,1), 0);
        bin_out[0][i] = Eigen::Vector3d(directions[i][0], directions[i][1], 0);
    }
    std::cout<<"Write out the unfolded curve: "<<std::endl;
    std::string fname = igl::file_dialog_save();
    if (fname.length() == 0)
    {
        std::cout<<"Please type something "<<std::endl;
        return;
    }

    write_polyline_xyz(ply_out, fname);
    write_polyline_xyz(bin_out, fname + "_b");
    std::cout << "files get saved" << std::endl;
}


void construct_single_developable_strips_by_intersect_rectifying(const int which)
{

    std::vector<std::vector<Eigen::Vector3d>> vertices_in;
    std::vector<std::vector<Eigen::Vector3d>> binormals_in;
    std::cout
        << "Read the plyline and binormal files" << std::endl;
    
    read_plylines_and_binormals(vertices_in, binormals_in);
    std::vector<Eigen::Vector3d> vertices = vertices_in[which];
    std::vector<Eigen::Vector3d> binormals = binormals_in[which];
    int vnbr = vertices.size();
    // the creases lines are constructed by foot points and creases
    // the foot is the intersection of the line and the right plane. vout and crease construct the strip
    std::vector<Eigen::Vector3d> tangents, creases, foot;
    tangents.resize(vnbr);
    creases.resize(vnbr);
    foot.resize(vnbr);


    tangents.front() = (vertices[1] - vertices[0]).normalized();
    tangents.back() = (vertices[vnbr - 1] - vertices[vnbr - 2]).normalized();
    tangents[1] = ((vertices[1] - vertices[0]).normalized() + (vertices[2] - vertices[1]).normalized()).normalized();
    

    foot.front() = vertices.front();
    foot[1] = vertices[1];
    foot.back() = vertices.back();

    creases[0] = binormals[0];
    creases[1] = binormals[1];
    creases[vnbr - 1] = binormals[vnbr - 1];

    for (int i = 2; i < vnbr - 1; i++)
    {
        Eigen::Vector3d dirl = (vertices[i] - vertices[i-1]).normalized();
        Eigen::Vector3d dirr = (vertices[i+1] - vertices[i]).normalized();

        // tangents
        Eigen::Vector3d tleft = tangents[i - 1];
        Eigen::Vector3d tright = (dirl + dirr).normalized();
        // binormals
        Eigen::Vector3d bleft = binormals[i - 1];
        Eigen::Vector3d bright = binormals[i];
        // normals of the both planes
        Eigen::Vector3d nleft = tleft.cross(bleft).normalized();
        Eigen::Vector3d nright = tright.cross(bright).normalized();
        // this crease direction is the intersection of the both planes. if the planes are parallel, use the binormal.
        Eigen::Vector3d c = nleft.cross(nright);
        bool parallel = false;
        if (c.norm() < 1e-4)
        {
            c = bleft;
            parallel = true;
        }
        c.normalize();
        if (c.dot(bleft) < 0)
        {
            c *= -1;
        }
        double scale = 1 / c.dot(bleft);
        c *= scale;
        // the foot point is the intersection of the left tangent with the right plane
        Eigen::Vector3d f;
        if(parallel){
            f = (vertices[i - 1] + vertices[i]) / 2;
        }
        else{
            // (vleft + alpha * tleft - v).dot(nright) = 0
            double alpha = (vertices[i] - vertices[i - 1]).dot(nright) / tleft.dot(nright);
            f = vertices[i - 1] + alpha * tleft;
        }
        tright = (vertices[i] - f).normalized();
        tangents[i] = tright;
        foot[i] = f;
        creases[i] = c;
    }
    std::vector<std::vector<Eigen::Vector3d>> vout, bout;
    vout.push_back(foot);
    bout.push_back(creases);
    std::cout<<"Write out the space curve: "<<std::endl;
    std::string fname = igl::file_dialog_save();
    if (fname.length() == 0)
    {
        std::cout<<"Please type something "<<std::endl;
        return;
    }

    write_polyline_xyz(vout, fname);
    write_polyline_xyz(bout, fname + "_b");
    std::cout << "files get saved" << std::endl;
}
void construct_single_developable_strips_by_intersect_rectifying(
    const std::vector<Eigen::Vector3d> &vertices, const std::vector<Eigen::Vector3d> &binormals,
    std::vector<Eigen::Vector3d> &vout, std::vector<Eigen::Vector3d> &bout)
{

    int vnbr = vertices.size();
    // the creases lines are constructed by foot points and creases
    // the foot is the intersection of the line and the right plane. vout and crease construct the strip
    std::vector<Eigen::Vector3d> tangents, creases, foot;
    tangents.resize(vnbr);
    creases.resize(vnbr);
    foot.resize(vnbr);


    tangents.front() = (vertices[1] - vertices[0]).normalized();
    tangents.back() = (vertices[vnbr - 1] - vertices[vnbr - 2]).normalized();
    tangents[1] = ((vertices[1] - vertices[0]).normalized() + (vertices[2] - vertices[1]).normalized()).normalized();
    

    foot.front() = vertices.front();
    foot[1] = vertices[1];
    foot.back() = vertices.back();

    creases[0] = binormals[0];
    creases[1] = binormals[1];
    creases[vnbr - 1] = binormals[vnbr - 1];

    for (int i = 2; i < vnbr - 1; i++)
    {
        Eigen::Vector3d dirl = (vertices[i] - vertices[i-1]).normalized();
        Eigen::Vector3d dirr = (vertices[i+1] - vertices[i]).normalized();

        // tangents
        Eigen::Vector3d tleft = tangents[i - 1];
        Eigen::Vector3d tright = (dirl + dirr).normalized();
        // binormals
        Eigen::Vector3d bleft = binormals[i - 1];
        Eigen::Vector3d bright = binormals[i];
        // normals of the both planes
        Eigen::Vector3d nleft = tleft.cross(bleft).normalized();
        Eigen::Vector3d nright = tright.cross(bright).normalized();
        // this crease direction is the intersection of the both planes. if the planes are parallel, use the binormal.
        Eigen::Vector3d c = nleft.cross(nright);
        bool parallel = false;
        if (c.norm() < 1e-4)
        {
            c = bleft;
            parallel = true;
        }
        c.normalize();
        if (c.dot(bleft) < 0)
        {
            c *= -1;
        }
        double scale = 1 / c.dot(bleft);
        c *= scale;
        // the foot point is the intersection of the left tangent with the right plane
        Eigen::Vector3d f;
        if(parallel){
            f = (vertices[i - 1] + vertices[i]) / 2;
        }
        else{
            // (vleft + alpha * tleft - v).dot(nright) = 0
            double alpha = (vertices[i] - vertices[i - 1]).dot(nright) / tleft.dot(nright);
            f = vertices[i - 1] + alpha * tleft;
        }
        tright = (vertices[i] - f).normalized();
        tangents[i] = tright;
        foot[i] = f;
        creases[i] = c;
    }

    vout = foot;
    bout = creases;
}
// this is for initialize AAG, thus the offsets should be uniform, the generated mesh should also be uniform.
void construct_single_developable_strips_by_intersect_rectifying_AAG(
    const std::vector<Eigen::Vector3d> &vertices, const std::vector<Eigen::Vector3d> &binormals,
    std::vector<Eigen::Vector3d> &vout, std::vector<Eigen::Vector3d> &bout)
{

    int vnbr = vertices.size();
    // the creases lines are constructed by foot points and creases
    // the foot is the intersection of the line and the right plane. vout and crease construct the strip
    std::vector<Eigen::Vector3d> tangents, creases, foot;
    tangents.resize(vnbr);
    creases.resize(vnbr);
    foot.resize(vnbr);


    tangents.front() = (vertices[1] - vertices[0]).normalized();
    tangents.back() = (vertices[vnbr - 1] - vertices[vnbr - 2]).normalized();
    tangents[1] = ((vertices[1] - vertices[0]).normalized() + (vertices[2] - vertices[1]).normalized()).normalized();

    foot.front() = (vertices[0] + vertices[1]) / 2;
    // foot[1] = vertices[1];
    foot[vnbr - 2] = (vertices[vnbr - 1] + vertices[vnbr - 2]) / 2;// set up the second last foot for AAG

    creases[0] = binormals[0];
    // creases[1] = binormals[1];
    creases[vnbr - 2] = binormals[vnbr - 2]; // AAG
    int nc = 0; // computed intersections not skipped
    // crease[i] is associated with vi and vi+1
    for (int i = 1; i < vnbr - 2; i++)
    {
        Eigen::Vector3d dirl = (vertices[i + 1] - vertices[i]).normalized();
        Eigen::Vector3d dirr = (vertices[i + 2] - vertices[i + 1]).normalized();

        // tangents
        Eigen::Vector3d tleft = tangents[i];
        Eigen::Vector3d tright = (dirl + dirr).normalized();
        // binormals
        Eigen::Vector3d bleft = binormals[i];
        Eigen::Vector3d bright = binormals[i + 1];
        // normals of the both planes
        Eigen::Vector3d nleft = tleft.cross(bleft).normalized();
        Eigen::Vector3d nright = tright.cross(bright).normalized();
        // this crease direction is the intersection of the both planes. if the planes are parallel, use the binormal.
        Eigen::Vector3d c = nleft.cross(nright);
        bool parallel = false;
        if (c.norm() < 1e-4)
        {
            c = bleft;
            parallel = true;
        }
        c.normalize();
        if (c.dot(bleft) < 0)
        {
            c *= -1;
        }
        double scale = 1 / c.dot(bleft);
        c *= scale;
        // the foot point is the intersection of the left tangent with the right plane
        Eigen::Vector3d f;
        if(parallel){
            f = (vertices[i] + vertices[i + 1]) / 2;
        }
        else{
            nc++;
            // (vleft + alpha * tleft - v).dot(nright) = 0
            double alpha = (vertices[i + 1] - vertices[i]).dot(nright) / tleft.dot(nright);
            f = vertices[i] + alpha * tleft;
        }
        tright = (vertices[i] - f).normalized();
        tangents[i + 1] = tright;
        foot[i] = f;
        creases[i] = c;
    }
    // just got the foot and crease from 1~vnbr-3. for 0 and vnbr-2 is easy to obtain, next deal with vnbr-1
    // use 3 order smoothness condition v3 = v0 - 3v1 + 3v2
    foot[vnbr - 1] = foot[vnbr - 4] - 3 * foot[vnbr - 3] + 3 * foot[vnbr - 2]; // the target foot point
    Eigen::Vector3d cm4 = foot[vnbr - 4] + creases[vnbr - 4];
    Eigen::Vector3d cm3 = foot[vnbr - 3] + creases[vnbr - 3];
    Eigen::Vector3d cm2 = foot[vnbr - 2] + creases[vnbr - 2];
    Eigen::Vector3d cm = cm4 - 3 * cm3 + 3 * cm2; // the target offset point

    creases[vnbr - 1] = cm - foot[vnbr - 1]; // the target crease.

    // unitscale the creases so that the shape can be uniform.
    for (int i = 0; i < vnbr - 1; i++){
        double length = (vertices[i] - vertices[i + 1]).norm();
        creases[i] *= length * sqrt(3) / 2;
    }
    creases[vnbr - 1] *= creases[vnbr - 2].norm() / creases[vnbr - 1].norm();
    vout = foot;
    bout = creases;
    std::cout<<"not skipped intersection, "<<nc<<"\n";
}
void construct_developable_strips_by_intersect_rectifying(){
    std::vector<std::vector<Eigen::Vector3d>> vertices_in, vertices_out;
    std::vector<std::vector<Eigen::Vector3d>> binormals_in, binormals_out;
    std::cout
        << "Read the plyline and binormal files" << std::endl;
    
    read_plylines_and_binormals(vertices_in, binormals_in);
    int cnbr = vertices_in.size();
    for (int i = 0; i < cnbr; i++)
    {
        std::vector<Eigen::Vector3d> vers, bnms;
        construct_single_developable_strips_by_intersect_rectifying(vertices_in[i], binormals_in[i], vers, bnms);
        vertices_out.push_back(vers);
        binormals_out.push_back(bnms);
    }
    std::string fname = igl::file_dialog_save();
    if (fname.length() == 0)
    {
        std::cout<<"Please type something "<<std::endl;
        return;
    }

    write_polyline_xyz(vertices_out, fname);
    write_polyline_xyz(binormals_out, fname + "_b");
    std::cout << "files get saved" << std::endl;
}

void draw_catenaries_on_cylinder(){
    double radius = 5;
    double height = 10;
    double angle_degree = 60;

    double angle_radian = angle_degree * LSC_PI / 180.;
    int vnbr = 133;
    double tmin = -LSC_PI / 2 - 0.288514;
    double tmax = LSC_PI / 2 - 0.288514;
    double titv = (tmax - tmin) / (vnbr - 1);
    
    std::string fname = igl::file_dialog_save();
    if (fname.length() == 0)
    {
        std::cout<<"Please type something "<<std::endl;
        return;
    }
    double clength = 0;
    std::ofstream file;
    file.open(fname);
    for (int i = 0; i < vnbr; i++)
    {
        double t = tmin + i * titv;
        double x = radius * cos(t + LSC_PI / 2);// add pi/2 to have a rotation of 90 degree. 
        double y = radius * sin(t + LSC_PI / 2);// add pi/2 to have a rotation of 90 degree. 
        double tangent = tan(angle_radian);
        double z = radius * tangent * cosh(t / tangent) - 5;//minus 5 to make a transportation.

        file<<"v "<<x<<" "<<y<<" "<<z<<std::endl;
        

    }
    file.close();
}

void get_single_isoline(const CGMesh &lsmesh, const std::vector<CGMesh::HalfedgeHandle>& loop, const Eigen::MatrixXd &V,
                   const Eigen::MatrixXi &F, const Eigen::VectorXd &ls, 
                   Eigen::MatrixXd &pts)
{
    double value = (ls.maxCoeff() + ls.minCoeff()) / 2;
    std::vector<std::vector<Eigen::Vector3d>> pts_list;
    std::vector<bool> left_large;
    get_iso_lines(lsmesh, loop, V, F, ls, value, pts_list, left_large);
    int counter = 0;
    for (auto p1 : pts_list)
    {
        for (auto p2 : p1)
        {
            counter++;
        }
    }
    pts.resize(counter, 3);
    counter = 0;
    for (auto p1 : pts_list)
    {
        for (auto p2 : p1)
        {
            pts.row(counter) = p2;
            counter++;
        }
    }
}

void write_curve_into_xyz_file(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &normals,
                               const Eigen::MatrixXd &pts, const std::string &prefix)
{
    
    Eigen::MatrixXd C;
    Eigen::VectorXi I;
    Eigen::VectorXd D;
    int nq = pts.rows();

    igl::point_mesh_squared_distance(pts, V, F, D, I, C);
    std::ofstream fout;
    fout.open(prefix + ".xyz");
    for (int i = 0; i < nq; i++)
    {
        int f = I(i);
        Eigen::Vector3d p = C.row(i);
        Eigen::Vector3d v0 = V.row(F(f, 0));
        Eigen::Vector3d v1 = V.row(F(f, 1));
        Eigen::Vector3d v2 = V.row(F(f, 2));
        std::array<double, 3> coor = barycenter_coordinate(v0, v1, v2, p);
        Eigen::Vector3d norm = coor[0] * normals.row(F(f, 0)) + coor[1] * normals.row(F(f, 1)) + coor[2] * normals.row(F(f, 2));
        assert(normals.row(F(f, 0)).dot(normals.row(F(f, 1))) > 0 && normals.row(F(f, 0)).dot(normals.row(F(f, 2))) > 0);
        norm.normalize();
        fout << pts(i, 0) << " " << pts(i, 1) << " " << pts(i, 2) << " " << norm(0) << " " << norm(1) << " " << norm(2) << "\n";
    }
    fout.close();
}

// unfinished
// ply is the moving one, curve is the fixed one
bool ply_ply_distance(const std::vector<Eigen::Vector3d> &ply, const std::vector<Eigen::Vector3d> &curve,
                      double &distance)
{

    double dis_max = 0;
    double clsu, clsv;
    Eigen::Vector3d closest;
    // std::cout<<"sizes, "<<ply.size()<<", "<<curve.size()<<std::endl;
    for (int i = 1; i < ply.size() - 1; i++)// skip boundary vertices
    {
        Eigen::Vector3d v = ply[i];
        double local_min = 1e6;
        int idmin = 0;
        for (int j = 0; j < curve.size() - 1; j++) 
        {
            Eigen::Vector3d vs = curve[j];
            Eigen::Vector3d ve = curve[j + 1];
            Eigen::Vector3d vsv = vs - v;
            Eigen::Vector3d ves = ve - vs;
            double t = -vsv.dot(ves) / ves.dot(ves);
            if (t < 0)
            {
                t = 0;
            }
            if (t > 1)
            {
                t = 1;
            }
            Eigen::Vector3d pjp = vs + (ve - vs) * t;
            double dis = (v - pjp).norm();
            if (dis < local_min)
            {
                local_min = dis;
            }
            // std::cout<<"DBG, check here"<<std::endl;
        }
        if (dis_max < local_min)
        {
            dis_max = local_min;
        }
    }
    distance = dis_max;

}


void match_the_two_catenaries(const CGMesh &lsmesh, const std::vector<CGMesh::HalfedgeHandle> &loop, const Eigen::MatrixXd &V,
                              const Eigen::MatrixXi &F, const Eigen::MatrixXd &normals, const Eigen::VectorXd &ls,
                              Eigen::MatrixXd& Vcout, Eigen::MatrixXd& Vlout)
{
    // get the lowest point
    double x = 0;
    double y = 5;
    double tangent = tan(60 * LSC_PI / 180);
    double z = 5 * tangent * cosh(0) - 5; // minus 5 to make a transportation.
    std::cout<<"middle point is actually "<<x<<" "<<y<<" "<<z<<std::endl;
    Eigen::MatrixXd pts;
    get_single_isoline(lsmesh, loop, V, F, ls, pts); // get the isoline from mesh
    std::cout<<"Reading the catenary file "<<std::endl;
    std::string fname = igl::file_dialog_open();
    if (fname.length() == 0)
    {
        std::cout<<"Please open something "<<std::endl;
        return;
    }
    Eigen::MatrixXd Vc;
    Eigen::MatrixXi Fc;
    igl::readOBJ(fname, Vc, Fc);

    int lowest_c;
    double zlow = Vc(0, 2);
    for (int i = 0; i < Vc.rows(); i++)
    {
        double z = Vc(i, 2);
        if (zlow > z)
        {
            zlow = z;
            lowest_c = i;
        }
    }
    int lowest_l;
    zlow = pts(0, 2);
    for (int i = 0; i < pts.rows(); i++)
    {
        double z = pts(i, 2);
        if (zlow > z)
        {
            zlow = z;
            lowest_l = i;
        }
    }
    // int lowest_l1;
    // double d0 = abs(pts(lowest_l, 2) - pts(lowest_l - 1, 2));
    // double d1 = abs(pts(lowest_l, 2) - pts(lowest_l + 1, 2));
    // if (d0 < d1)
    // {
    //     lowest_l1 = lowest_l - 1;
    // }
    // else
    // {
    //     lowest_l1 = lowest_l + 1;
    // }
    std::cout<<"The lowest point of ls: "<<pts.row(lowest_l)<<", caten: "<<Vc.row(lowest_c)<<std::endl;
    // Eigen::Vector3d pl = (pts.row(lowest_l) + pts.row(lowest_l1)) / 2;
     Eigen::Vector3d pl = pts.row(lowest_l);
    Eigen::Vector3d pc = Eigen::Vector3d(x, y, z);

    Eigen::Vector3d nl(pl[0], pl[1], 0), nc(pc[0], pc[1], 0);
    nl.normalize();
    nc.normalize();
    double angle = acos(nl.dot(nc));// this is roughly correct for small angle. For arbitrary angle need to modify.
    if (nc.cross(nl)[2] < 0)
    {
        angle *= -1;
    }
    std::cout << "Roughly rotated by angle " << angle << std::endl;
    // counterclockwise rotate angle
    Eigen::Matrix3d rot;
    rot << cos(angle), -sin(angle), 0,
        sin(angle), cos(angle), 0,
        0, 0, 1;
    double trans = (pl - pc)[2];
    Vcout = (rot * Vc.transpose()).transpose();
    Vcout.rowwise() += Eigen::Vector3d(0, 0, trans).transpose();
    Vlout = pts;
    std::cout<<"Writing curves into xyz format. type down the prefix:"<<std::endl;
    fname = igl::file_dialog_save();
    write_curve_into_xyz_file(V, F, normals, Vcout, fname + "_c");
    write_curve_into_xyz_file(V, F, normals, Vlout, fname + "_l");
    std::cout<<"finished writing xyz files"<<std::endl;

}
// angle is from 0 to pi
void slope_of_given_angle(const Eigen::Vector3d& norm, const double angle_radian, double &slope){
    double x = sin(angle_radian);
    double y = cos(angle_radian);
    double z = -(norm[0] * x + norm[1] * y) / norm[2];
    slope = z;
}
// phi is in degree.
// range is in radian, from small to large
// active if false, skip
// Fselect is either 1 or -1.
void shading_slope_range_compute(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::VectorXi& Fselect,
    const Eigen::MatrixXd& norm_f, const double phi_min,
    const double phi_max, std::vector<std::array<double,2>>& range, std::vector<bool> &active){
    double phi0 = phi_min * LSC_PI / 180 - LSC_PI / 2;
    double phi1 = phi_max * LSC_PI / 180 - LSC_PI / 2;
    int fnbr = F.rows();
    int vnbr = V.rows();
    range.resize(fnbr);
    active.resize(fnbr);
    for (int i = 0; i < fnbr; i++)
    {
        if (Fselect.size() > 0)
        {
            if (Fselect[i] < 0)
            {
                active[i] = false;
                continue;
            }
        }
        active[i] = true;
        Eigen::Vector3d v0 = V.row(F(i, 0));
        Eigen::Vector3d v1 = V.row(F(i, 1));
        Eigen::Vector3d v2 = V.row(F(i, 2));
        Eigen::Vector3d norm = norm_f.row(i);
        if (abs(abs(norm(2)) - 1) < 1e-4)
        { // the triangle is almost horizontal.
            active[i] = false;
            continue;
        }
        Eigen::Vector3d deepest;
        if (abs(norm(2)) < 1e-4) // the triangle is almost verticle, all slopes exist
        {
            range[i][0] = -LSC_PI / 2;
            range[i][1] = LSC_PI / 2;
            continue;
        }

        deepest[0] = norm[0];
        deepest[1] = norm[1];
        deepest[2] = -(deepest[0] * norm[0] + deepest[1] * norm[1]) / norm[2]; // deepest * norm = 0
        if (deepest[0] < 0)
        {
            deepest *= -1;
        }
        // get the corresponding phi of this vector
        double scale = deepest[0] * deepest[0] + deepest[1] * deepest[1];
        scale = sqrt(scale);
        // the sin and cos value of the angle. the range is from 0 to pi
        double sn = deepest[0] / scale;
        double cn = deepest[1] / scale;
        double deep_radian = acos(cn);
        bool including_deep = false;
        if (deep_radian < phi1 && deep_radian > phi0)
        {
            including_deep = true;
        }
        double slope0, slope1;
        slope_of_given_angle(norm, phi0, slope0);
        slope_of_given_angle(norm, phi1, slope1);
        if (!including_deep)
        {
            if (slope0 < slope1)
            {
                range[i][0] = atan(slope0);
                range[i][1] = atan(slope1);
            }
            else
            {
                range[i][0] = atan(slope1);
                range[i][1] = atan(slope0);
            }
        }
        else{
            double slopem;

            slopem = deepest[2] / scale;

            if (slope0 < slopem)
            {
                range[i][0] = atan(std::min(slope0, slope1));
                range[i][1] = atan(slopem);
            }
            else
            {
                range[i][0] = atan(slopem);
                range[i][1] = atan(std::max(slope0, slope1));
            }
        }
        
        // if(i<10){
        //     std::cout<<"i "<<i<<" angle min, "<<range[i][0]<<", max, "<<range[i][1]<<std::endl;
        // }
    }
}

// the angle is the angle of slopvec against y axis, from 0 to 2pi
// all_angle means if the norm is close to the reference so that the angle is random
// slope_vec is the vector that x*x + y*y = 1
// here the angle is determined by the vector = norm X ref. users should be careful when using this function about if the vector takes opposite direction.
void get_angle_of_vector_ortho_to_reference(const Eigen::Vector3d &norm, const Eigen::Vector3d& ref, Eigen::Vector3d& slopevec, double& angle_radian, bool& all_angle){
    all_angle = false;
    // if the norm is parallel with ref, all angles are allowed
    if (abs(abs(norm.dot(ref)) - 1) < 1e-4)
    {
        all_angle = true;
        return;
    }
    Eigen::Vector3d tangent = norm.cross(ref).normalized();
    double scale = tangent[0] * tangent[0] + tangent[1] * tangent[1];
    scale = sqrt(scale);
    slopevec = tangent/scale;
    double sn = slopevec[0];
    double cn = slopevec[1];
    double deep_radian = acos(cn); // if the range is from 0 to pi, then sn need not consider.
    if (sn < 0)
    {
        deep_radian = 2 * LSC_PI - deep_radian;
    }
    angle_radian = deep_radian;
}
// ref0 corresponds to phimin, ref1 corresponds to phimax.
void get_angle_range_of_vector_ortho_to_reference(const Eigen::Vector3d &norm, const Eigen::Vector3d &ref0, const Eigen::Vector3d &ref1,
                                                  double &angle_radian0, double &angle_radian1, bool &all_angle)
{
    all_angle = false;
    Eigen::Vector3d slopevec0, slopevec1;
    double angle0, angle1;
    bool aa0, aa1;
    get_angle_of_vector_ortho_to_reference(norm, ref0, slopevec0, angle0, aa0);
    get_angle_of_vector_ortho_to_reference(norm, ref1, slopevec1, angle1, aa1);
    if (aa0 || aa1)
    {
        all_angle = true;
        return;
    }

    // 2d vectors of them. they are already normalized
    Eigen::Vector3d s2d0(slopevec0[0], slopevec0[1], 0), s2d1(slopevec1[0], slopevec1[1], 0);
    Eigen::Vector3d cross = s2d0.cross(s2d1);
    double adiff = acos(s2d0.dot(s2d1));
    double a_ref; // angle of the slopevec1, it should be the smaller-angle one of the two vectors
    // make sure the cross is positive to the z axis
    bool invert = false;
    if (cross[2] < 0)
    {
        invert = true;
    }
    // if no need to change the angles
    if (!invert)
    {
        a_ref = angle1;
    }
    else{// if need to use the inverted direction
        a_ref = angle1 - LSC_PI;
    }
    double b_ref = a_ref + adiff;
    angle_radian0 = a_ref;
    angle_radian1 = b_ref;
}
void shading_slope_range_compute(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXi &Fselect,
                                 const Eigen::MatrixXd &norm_f, const double theta_degree, const double phi_min, const double phi_max,
                                 std::vector<std::array<double, 2>> &range, std::vector<bool> &active)
{
    Eigen::Vector3d vecref0, vecref1;
    vecref0 = angle_ray_converter(theta_degree, phi_min);
    vecref1 = angle_ray_converter(theta_degree, phi_max);
    // phi0 and phi1 are for the range of tangent vectors!!!
    int fnbr = F.rows();
    int vnbr = V.rows();
    range.resize(fnbr);
    active.resize(fnbr);
    for (int i = 0; i < fnbr; i++)
    {
        if (Fselect.size() > 0)
        {
            if (Fselect[i] < 0)
            {
                active[i] = false;
                continue;
            }
        }
        active[i] = true;
        Eigen::Vector3d v0 = V.row(F(i, 0));
        Eigen::Vector3d v1 = V.row(F(i, 1));
        Eigen::Vector3d v2 = V.row(F(i, 2));
        Eigen::Vector3d norm = norm_f.row(i);
        if (abs(abs(norm(2)) - 1) < 1e-4)
        { // the triangle is almost horizontal, it only has slope = 0.
            active[i] = false;
            continue;
        }
        double ang0, ang1;
        bool all_angle;
        // The ang0 ang1 will be almost in the range of (0, pi)
        get_angle_range_of_vector_ortho_to_reference(norm, vecref0, vecref1, ang0, ang1, all_angle);

        // the triangle is almost verticle, all slopes exist, but it also creates difficulty in getting correct orientation of the tangents.
        if (abs(norm(2)) < 1e-4 || all_angle) 
        {
            // range[i][0] = -LSC_PI / 2;
            // range[i][1] = LSC_PI / 2;
            active[i] = false;
            continue;
        }
        Eigen::Vector3d deepest; // it is a vector that lies in the line projected by the norm to the xy plane, and its x>0, 
        deepest[0] = norm[0];
        deepest[1] = norm[1];
        deepest[2] = -(deepest[0] * norm[0] + deepest[1] * norm[1]) / norm[2]; // deepest * norm = 0
        if (deepest[0] < 0)
        {
            deepest *= -1;
        }
        // get the corresponding phi of this vector
        double scale = deepest[0] * deepest[0] + deepest[1] * deepest[1];
        scale = sqrt(scale);
        // the sin and cos value of the angle against y axis. the range is from 0 to pi
        double sn = deepest[0] / scale;
        double cn = deepest[1] / scale;
        double deep_radian = acos(cn);// if the range is from 0 to pi, then sn need not consider.
        bool including_deep1 = false; // include it or the opposite of it
        if (deep_radian < ang1 && deep_radian > ang0)
        {
            including_deep1 = true;
        }
        // if (deep_radian > 0 && deep_radian < LSC_PI / 2)// if deepest angle is [0, pi/2], then its opposite (angle + pi) could be in the range
        // {
        //     if (deep_radian + LSC_PI < ang1 && deep_radian + LSC_PI > ang0)
        //     {
        //         including_deep2 = true;
        //     }
        // }
        double slope0, slope1;
        slope_of_given_angle(norm, ang0, slope0);
        slope_of_given_angle(norm, ang1, slope1);
        if (!including_deep1)
        {
            if (slope0 < slope1)
            {
                range[i][0] = atan(slope0);
                range[i][1] = atan(slope1);
            }
            else
            {
                range[i][0] = atan(slope1);
                range[i][1] = atan(slope0);
            }
        }
        else
        {
            double slopem;

            slopem = deepest[2] / scale;
          
            if (slope0 < slopem)
            {
                range[i][0] = atan(std::min(slope0, slope1));
                range[i][1] = atan(slopem);
            }
            else
            {
                range[i][0] = atan(slopem);
                range[i][1] = atan(std::max(slope0, slope1));
            }
        }

        // if(i<10){
        //     std::cout<<"i "<<i<<" angle min, "<<range[i][0]<<", max, "<<range[i][1]<<std::endl;
        // }
    }
}
// theta, phi is in degree.
// range is in radian
// active if false, skip
// Fselect is either 1 or -1.
void lighting_slope_range_compute(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::VectorXi& Fselect,
                                  const Eigen::MatrixXd &norm_f,  const double theta, const double phi,
                                  std::vector<std::array<double, 2>> &range, std::vector<bool> &active)
{
    double phi_radian = phi * LSC_PI / 180;
    double theta_radian = theta * LSC_PI / 180;
    Eigen::Vector3d Refer_dir;
    Refer_dir << cos(theta_radian) * sin(phi_radian), cos(theta_radian) * cos(phi_radian), sin(theta_radian);
    std::cout << "The reference axis for lighting, " << Refer_dir.transpose() << "\n";
    int fnbr = F.rows();
    int vnbr = V.rows();
    range.resize(fnbr);
    active.resize(fnbr);
    for (int i = 0; i < fnbr; i++)
    {
         if (Fselect.size() > 0)
        {
            if (Fselect[i] < 0)
            {
                active[i] = false;
                continue;
            }
        }
        active[i] = true;
        Eigen::Vector3d v0 = V.row(F(i, 0));
        Eigen::Vector3d v1 = V.row(F(i, 1));
        Eigen::Vector3d v2 = V.row(F(i, 2));
        Eigen::Vector3d norm = norm_f.row(i);
        if (abs(abs(norm.dot(Refer_dir)) - 1) < 1e-4)
        { // the triangle is almost horizontal.
            active[i] = false;
            continue;
        }
        
        if (abs(norm.dot(Refer_dir)) < 1e-4) // the triangle is almost verticle, all slopes exist
        {
            range[i][0] = -LSC_PI / 2;
            range[i][1] = LSC_PI / 2;
            continue;
        }
        Eigen::Vector3d lxn = Refer_dir.cross(norm).normalized();

        Eigen::Vector3d proj0 = norm.cross(lxn).normalized(); // corresponds to maximal angle
        if (proj0.dot(Refer_dir) < 0)
        {
            proj0 *= -1;
        }
        double radian_max = LSC_PI / 2 - acos(proj0.dot(Refer_dir));
        double radian_min = - radian_max;
        range[i][0] = radian_min;
        range[i][1] = radian_max;
        assert(radian_min > - LSC_PI / 2);
        assert(radian_max < LSC_PI / 2);
    }
}

void vote_for_slope(const std::vector<std::array<double, 2>> &range, const std::vector<bool> &active,
                    std::vector<int> &flist, double &angle_degree)
{
    int fnbr = range.size();
    // only consider about angle in [0, pi]. resolution is 5 degree. totally 37 different values.
    std::array<std::vector<int>,37> candidates;
    for (int i = 0; i < 37; i++)
    {
        candidates[i].reserve(100);
    }
    assert(active.size()==fnbr);
    for (int i = 0; i < fnbr; i++)
    {
        if (active[i] == false)
        {
            continue;
        }
        double lower = range[i][0] * 180 / LSC_PI;
        double upper = range[i][1] * 180 / LSC_PI;
        lower = int(lower / 5) * 5;
        upper = int(upper / 5) * 5;
        if (lower == upper)
        {
            int location = lower / 5 + 18;
            assert(location >= 0);
            candidates[location].push_back(i);
            continue;
        }
        int lowid = lower / 5 + 18;
        int upid = upper / 5 + 18;
        assert(lowid <= upid);
        for (int j = lowid; j <= upid; j++)
        {
            assert(j>=0);
            assert(j<37);
            candidates[j].push_back(i);
        }
    }
    Eigen::VectorXd nvotes(37);
    int max_id = 0;
    int max_nbr = 0;
    for (int i = 0; i < 37; i++)
    {
        nvotes[i] = candidates[i].size();
        if (nvotes[i] > max_nbr)
        {
            max_nbr = nvotes[i];
            max_id = i;
        }
        std::cout << "Voting: angle: " << (i - 18) * 5 << ", nbr, " << nvotes[i] << std::endl;
    }
    std::cout << "Statistics of slopes: the most common angle " << (max_id-18) * 5 << ", the nbr of corresponding faces, " << max_nbr << std::endl;
    flist = candidates[max_id];
    angle_degree = (max_id - 18) * 5;
}

bool validate_slope_root(const Eigen::Vector3d &norm, const double z, const double c)
{
    if (c > 1)
    {
        return false;
    }
    double x = sqrt(1 - c * c);
    double y = c;
    Eigen::Vector3d vec(x, y, z);
    if (abs(norm.dot(vec)) > 1e-3)
    {
        return false;
    }
    return true;
}
// the angle is between direc and the plane orthogonal to the axis
// projection is the axis projected to the plane
bool validate_slope_root_angle(const Eigen::Vector3d &norm, const Eigen::Vector3d &axis,
                               const Eigen::Vector3d &direc, const Eigen::Vector3d &proj, const double angle_radian)
{
    // orthogonal to the normal vector
    if (abs(norm.dot(direc)) > 1e-3)
    {
        return false;
    }
    // the angle_radian is the angle between direc and the plane orthogonal to axis
    double angle_compute = LSC_PI / 2 - acos(axis.dot(direc));
    if (abs(angle_compute - angle_radian) > 0.05)
    {
        return false;
    }
    // if the angel is 0, then only one solution.
    if (angle_radian == 0)
    {
        return true;
    }
    // filter out the other solution
    Eigen::Vector3d cross = direc.cross(proj).normalized();
    if (cross.dot(norm) < 0)
    {
        return false;
    }
    return true;
}

// phi_min, phi_max are in degree.
void get_orthogonal_vector_of_selected_slope_shading(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::VectorXi &Fselect,
                                             const Eigen::MatrixXd &norm_f, const double phi_min,
                                             const double phi_max, std::vector<int> &fout, std::vector<Eigen::Vector3d> &ortho)
{
    fout.clear();
    ortho.clear();
    std::vector<std::array<double, 2>> range;
    std::vector<bool> active;
    shading_slope_range_compute(V, F,Fselect, norm_f, phi_min, phi_max, range, active);
    std::vector<int> flist;
    double angle_degree;
    vote_for_slope(range, active, flist, angle_degree);
    double angle_radian = angle_degree * LSC_PI / 180;
    int case1=0, case2=0, case3 =0;
    for (int i = 0; i < flist.size(); i++)
    {
        int fid = flist[i];
        Eigen::Vector3d norm = norm_f.row(fid);
        double z = tan(angle_radian);
        double a0 = norm[0];
        double a1 = norm[1];
        double a2 = norm[2];
        std::vector<double> coeffs(3);
        coeffs[0] = a2 * a2 * z * z - a0 * a0;
        coeffs[1] = 2 * a1 * a2 * z;
        coeffs[2] = a0 * a0 + a1 * a1;
        std::array<double, 2> roots;
        bool solve = quadratic_solver(coeffs, roots);
        if(!solve){
            case1++;
            continue;
        }
        bool vali0 = validate_slope_root(norm, z, roots[0]);
        bool vali1 = validate_slope_root(norm, z, roots[1]);
        // stratigy: if no root or two valid roots, skip. if only one root, accept it.
        if(vali0 * vali1 == 1){//no roots or two roots
            case2++;
            continue;
        }
        Eigen::Vector3d direction;
        if(vali0){
            double c = roots[0];
            double x = sqrt(1 - c * c);
            double y = c;
            direction << x, y, z;
        }
        if(vali1){
            double c = roots[1];
            double x = sqrt(1 - c * c);
            double y = c;
            direction << x, y, z;
        }
        fout.push_back(fid);
        Eigen::Vector3d oout = norm.cross(direction).normalized();
        ortho.push_back(oout);
    }
    std::cout<<"For Shading: The final face nbr of slopes: "<<fout.size()<<std::endl;
    std::cout<<"No roots: "<<case1<<", 2 or 0 valid roots: "<<case2<<std::endl;
}
// phi_min, phi_max are in degree. This version consider about theta
void get_orthogonal_vector_of_selected_slope_shading(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, Eigen::VectorXi &Fselect,
                                                     const Eigen::MatrixXd &norm_f, const double theta_degree,
                                                     const double phi_min, const double phi_max, std::vector<int> &fout, std::vector<Eigen::Vector3d> &ortho)
{
    fout.clear();
    ortho.clear();
    std::vector<std::array<double, 2>> range;
    std::vector<bool> active;
    shading_slope_range_compute(V, F,Fselect, norm_f, theta_degree, phi_min, phi_max, range, active);
    std::vector<int> flist;
    double angle_degree;
    vote_for_slope(range, active, flist, angle_degree);
    double angle_radian = angle_degree * LSC_PI / 180;
    int case1=0, case2=0, case3 =0;
    for (int i = 0; i < flist.size(); i++)
    {
        int fid = flist[i];
        Eigen::Vector3d norm = norm_f.row(fid);
        double z = tan(angle_radian);
        double a0 = norm[0];
        double a1 = norm[1];
        double a2 = norm[2];
        std::vector<double> coeffs(3);
        coeffs[0] = a2 * a2 * z * z - a0 * a0;
        coeffs[1] = 2 * a1 * a2 * z;
        coeffs[2] = a0 * a0 + a1 * a1;
        std::array<double, 2> roots;
        bool solve = quadratic_solver(coeffs, roots);
        if(!solve){
            case1++;
            continue;
        }
        bool vali0 = validate_slope_root(norm, z, roots[0]);
        bool vali1 = validate_slope_root(norm, z, roots[1]);
        // stratigy: if no root or two valid roots, skip. if only one root, accept it.
        if(vali0 * vali1 == 1){//no roots or two roots
            case2++;
            continue;
        }
        Eigen::Vector3d direction;
        if(vali0){
            double c = roots[0];
            double x = sqrt(1 - c * c);
            double y = c;
            direction << x, y, z;
        }
        if(vali1){
            double c = roots[1];
            double x = sqrt(1 - c * c);
            double y = c;
            direction << x, y, z;
        }
        fout.push_back(fid);
        Eigen::Vector3d oout = norm.cross(direction).normalized();
        ortho.push_back(oout);
    }
    std::cout<<"For Shading: The final face nbr of slopes: "<<fout.size()<<std::endl;
    std::cout<<"No roots: "<<case1<<", 2 or 0 valid roots: "<<case2<<std::endl;
}

void print_validation_vector(const Eigen::Vector3d &norm, const Eigen::Vector3d &axis,
                             const Eigen::Vector3d &direc, const Eigen::Vector3d& proj, const double angle_radian)
{
    double dot = direc.dot(axis);
    double angle_compute = LSC_PI/2 - acos(dot);
    Eigen::Vector3d cross = direc.cross(proj).normalized();
    std::cout << "Target angle " << angle_radian << ", compute angle " << angle_compute << std::endl;
    std::cout << "Orthongl to norm, " << norm.dot(direc) << ", orientation wrt norm "<<cross.dot(norm)<<"\n" ;

}

// phi_min, phi_max are in degree.
void get_orthogonal_vector_of_selected_slope_lighting(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,Eigen::VectorXi &Fselect,
                                             const Eigen::MatrixXd &norm_f, const double theta,
                                             const double phi, std::vector<int> &fout, std::vector<Eigen::Vector3d> &ortho)
{
    fout.clear();
    ortho.clear();
    std::vector<std::array<double, 2>> range;
    std::vector<bool> active;
    lighting_slope_range_compute(V, F, Fselect, norm_f, theta, phi, range, active);
    std::vector<int> flist;
    double angle_degree;
    vote_for_slope(range, active, flist, angle_degree);
    double angle_radian = angle_degree * LSC_PI / 180;
    int case1 = 0, case2 = 0, case3 = 0;
    double phi_radian = phi * LSC_PI / 180;
    double theta_radian = theta * LSC_PI / 180;
    Eigen::Vector3d Refer_dir;
    Refer_dir << cos(theta_radian) * sin(phi_radian), cos(theta_radian) * cos(phi_radian), sin(theta_radian);
    for (int i = 0; i < flist.size(); i++)
    {
        int fid = flist[i];
        Eigen::Vector3d norm = norm_f.row(fid);
        
        Eigen::Vector3d lxn = Refer_dir.cross(norm).normalized();
        Eigen::Vector3d proj0 = norm.cross(lxn).normalized(); // corresponds to maximal angle
        if (angle_radian == 0)
        { // if angle is 0, there is only one solution
            Eigen::Vector3d oout;

            oout = proj0;
            if (validate_slope_root_angle(norm, Refer_dir, lxn, proj0, angle_radian))
            {
                ortho.push_back(oout);
                fout.push_back(fid);
            }
            continue;
        }
        else{
            // the vector is a * lxn + b * proj0
            double t1 = lxn.dot(Refer_dir);
            double t2 = proj0.dot(Refer_dir);
            double c = cos(LSC_PI / 2 - angle_radian); // the angle between the axis and the vector direction
            std::vector<double> coeffs(3);
            coeffs[0] = c * c - t1 * t1;
            coeffs[1] = -2 * c * t2;
            coeffs[2] = t1 * t1 + t2 * t2;

            std::array<double, 2> roots;
            bool solve = quadratic_solver(coeffs, roots);
            if (!solve)
            {
                case1++;
                continue;
            }
            double b1 = roots[0], b2 = roots[0], b3 = roots[1], b4 = roots[1];
            double a1 = sqrt(1 - b1 * b1), a2 = -sqrt(1 - b1 * b1), a3 = sqrt(1 - b3 * b3), a4 = -sqrt(1 - b3 * b3);
            Eigen::Vector3d vec1 = a1 * lxn + b1 * proj0;
            Eigen::Vector3d vec2 = a2 * lxn + b2 * proj0;
            Eigen::Vector3d vec3 = a3 * lxn + b3 * proj0;
            Eigen::Vector3d vec4 = a4 * lxn + b4 * proj0;
            bool vali1 = validate_slope_root_angle(norm, Refer_dir, vec1, proj0, angle_radian);
            bool vali2 = validate_slope_root_angle(norm, Refer_dir, vec2, proj0, angle_radian);
            bool vali3 = validate_slope_root_angle(norm, Refer_dir, vec3, proj0, angle_radian);
            bool vali4 = validate_slope_root_angle(norm, Refer_dir, vec4, proj0, angle_radian);
            if(vali1)
            {
                ortho.push_back(vec1);
                fout.push_back(fid);
                continue;
            }
            if(vali2)
            {
                ortho.push_back(vec2);
                fout.push_back(fid);
                continue;
            }
            if(vali3)
            {
                ortho.push_back(vec3);
                fout.push_back(fid);
                continue;
            }
            if(vali4)
            {
                ortho.push_back(vec4);
                fout.push_back(fid);
                continue;
            }
            case2++;
            std::cout<<"\nPrinting Four Roots"<<std::endl;
            print_validation_vector(norm, Refer_dir, vec1, proj0, angle_radian);
            print_validation_vector(norm, Refer_dir, vec2, proj0, angle_radian);
            print_validation_vector(norm, Refer_dir, vec3, proj0, angle_radian);
            print_validation_vector(norm, Refer_dir, vec4, proj0, angle_radian);
            
        }
    }
    std::cout<<"For Lighting: The final face nbr of slopes: "<<fout.size()<<std::endl;
    std::cout<<"No roots: "<<case1<<", No valid roots: "<<case2<<std::endl;
}

// types: 0: shading. 1: light through. 2: partly shading. 3. mixing shading and light through
// info is marked on each vertex: 2: transition, 1: reference patch, 0: ordinary points
// theta2 and phi2 are for light through when type is 3.
void orthogonal_slope_for_different_shading_types(const int whichtype, Eigen::VectorXi InnerV, const Eigen::VectorXi &info,
                                                  const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                                                  const Eigen::MatrixXd &norm_f, const double theta,
                                                  const double phi, const double theta_tol, const double phi_tol,
                                                  const double theta2, const double phi2,
                                                  std::vector<int> &fout, std::vector<Eigen::Vector3d> &ortho)
{
    fout.clear();
    ortho.clear();
    int fnbr = F.rows();
    if (whichtype == 0) // simply shading
    {
        Eigen::VectorXi Fselect; // this is empty, meaning not selecting
        double phi_min = phi - phi_tol;
        double phi_max = phi + phi_tol;
        get_orthogonal_vector_of_selected_slope_shading(V, F, Fselect, norm_f, theta, phi_min, phi_max, fout, ortho);
    }
    if (whichtype == 1) // simply light throughs
    {
        Eigen::VectorXi Fselect; // this is empty, meaning not selecting
        get_orthogonal_vector_of_selected_slope_lighting(V, F, Fselect, norm_f, theta, phi, fout, ortho);
    }
    if (whichtype == 2)// partly shading
    {
        double phi_min = phi - phi_tol;
        double phi_max = phi + phi_tol;

        Eigen::VectorXi Fselect = -1 * Eigen::VectorXi::Ones(fnbr);
        for (int i = 0; i < fnbr; i++)
        {
            int vid0 = F(i, 0);
            int vid1 = F(i, 1);
            int vid2 = F(i, 2);
            int vi0 = InnerV[vid0];
            int vi1 = InnerV[vid1];
            int vi2 = InnerV[vid2];
            if (vi0 < 0 || vi1 < 0 || vi2 < 0)
            {
                continue;
            }
            bool active_face = info[vi0] == 1 && info[vi1] == 1 && info[vi2] == 1;
            if(active_face){
                Fselect[i] = 1;
            }
        }
        get_orthogonal_vector_of_selected_slope_shading(V, F, Fselect, norm_f,theta, phi_min, phi_max, fout, ortho);

    }
    if (whichtype == 3)
    {
        double phi_min = phi - phi_tol;
        double phi_max = phi + phi_tol;
        // 1 for shading, 2 for lighting
        std::vector<int> fout1, fout2;
        std::vector<Eigen::Vector3d> ortho1, ortho2;
        Eigen::VectorXi Fselect1 = -1 * Eigen::VectorXi::Ones(fnbr), Fselect2 = -1 * Eigen::VectorXi::Ones(fnbr);
        for (int i = 0; i < fnbr; i++)
        {
            int vid0 = F(i, 0);
            int vid1 = F(i, 1);
            int vid2 = F(i, 2);
            int vi0 = InnerV[vid0];
            int vi1 = InnerV[vid1];
            int vi2 = InnerV[vid2];
            if (vi0 < 0 || vi1 < 0 || vi2 < 0)
            {
                continue;
            }
            bool active_face_shading = info[vi0] == 0 && info[vi1] == 0 && info[vi2] == 0;
            bool active_face_through = info[vi0] == 1 && info[vi1] == 1 && info[vi2] == 1; // active for light through
            if (active_face_shading)
            {
                Fselect1[i] = 1;
            }
            if(active_face_through){
                Fselect2[i] = 1;
            }
        }
        get_orthogonal_vector_of_selected_slope_shading(V, F, Fselect1, norm_f, theta, phi_min, phi_max, fout1, ortho1);
        get_orthogonal_vector_of_selected_slope_lighting(V, F, Fselect2, norm_f, theta2, phi2, fout2, ortho2);

        fout = fout1;
        ortho = ortho1;
        for (int i = 0; i < fout2.size(); i++)
        {
            fout.push_back(fout2[i]);
            ortho.push_back(ortho2[i]);
        }
    }
}

void lsTools::show_slopes(const double scaling, Eigen::MatrixXd& E0, Eigen::MatrixXd& E1){
    if(Fslope.size()==0){
        return;
    }
    E0.resize(Fslope.size(), 3);
    E1.resize(Fslope.size(), 3);
    for (int i = 0; i < Fslope.size(); i++)
    {
        int fid = Fslope[i];
        int vid0 = F(fid, 0);
        int vid1 = F(fid, 1);
        int vid2 = F(fid, 2);

        Eigen::Vector3d center = (V.row(vid0) + V.row(vid1) + V.row(vid2)) / 3;
        Eigen::Vector3d norm = norm_f.row(fid);
        Eigen::Vector3d ortho = OrthoSlope[i];
        Eigen::Vector3d direction = norm.cross(ortho).normalized();
        E0.row(i) = center + scaling * direction;
        E1.row(i) = center - scaling * direction;
    }
    std::cout<<"The nbr of edges drawn "<<E0.rows()<<std::endl;
    // std::cout<<"Showing the slopes, printing "<<std::endl;
    // for (int i = 0; i < Fslope.size(); i++)
    // {
    //     int fid = Fslope[i];
    //     int vid0 = F(fid, 0);
    //     int vid1 = F(fid, 1);
    //     int vid2 = F(fid, 2);

    //     Eigen::Vector3d center = (V.row(vid0) + V.row(vid1) + V.row(vid2)) / 3;
    //     Eigen::Vector3d norm = norm_f.row(fid);
    //     Eigen::Vector3d ortho = OrthoSlope[i];
    //     Eigen::Vector3d direction = norm.cross(ortho).normalized();
    //     double x = direction[0];
    //     double y = direction[1];
    //     double z = direction[2];
    //     double scale = x * x + y * y;
    //     scale = sqrt(scale);
    //     z = z/scale;
    //     if(z<0){
    //         z=-z;
    //     }
    //     std::cout<<z<<", ";
    // }
    // std::cout<<"\n";

}
// if positively curved, compute the direction orthogonal to the smallest curvature.
// if negatively curved, compute the cos^2. cos is between curvDir[i][0] and the target direction
void orth_smallest_curvature_and_c2(const std::vector<std::vector<double>> &curv,
                                        const std::vector<std::vector<Eigen::Vector3d>> &curvDir, 
                                        std::vector<int> &idspos, 
                                        std::vector<int> & idsneg, std::vector<Eigen::Vector3d> &ortho, 
                                        std::vector<double>& coscos)
{
    assert(curv.size() == curvDir.size());
    idspos.clear();
    idsneg.clear();
    ortho.clear();
    coscos.clear();
    int fnbr = curv.size();
    for (int i = 0; i < fnbr; i++)
    {
        assert(curvDir[i].size() == 2);
        double c1 = curv[i][0];
        double c2 = curv[i][1];
        idspos.push_back(i);
        if (abs(c1) < abs(c2)) // c1 is smallest, pick the direction of c2
        {
            ortho.push_back(curvDir[i][1].normalized());
        }
        else
        {
            ortho.push_back(curvDir[i][0].normalized());
        }
        // if (c1 * c2 > 0) // same sign, directly pick the smallest absolute curvature
        // {
        //     idspos.push_back(i);
        //     if (abs(c1) < abs(c2)) // c1 is smallest, pick the direction of c2
        //     {
        //         ortho.push_back(curvDir[i][1].normalized());
        //     }
        //     else
        //     {
        //         ortho.push_back(curvDir[i][0].normalized());
        //     }
        // }
        // else
        // {
        //     idsneg.push_back(i);
            
        //     double cc = -c2 / (c1 - c2);
        //     double ss = c1 / (c1 - c2);
        //     coscos.push_back(cc);
        // }
    }
}

void get_orthogonal_direction_minimal_principle_curvature(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                                                          std::vector<int> &idspos,
                                                          std::vector<int> &idsneg, std::vector<Eigen::Vector3d> &ortho,
                                                          std::vector<double> &coscos, std::vector<std::vector<Eigen::Vector3d>>& CurvDir)
{
    int vnbr = V.rows();
    int fnbr = F.rows();
    // Precomputation
    QuadricCalculator cc;
    cc.init(V, F);
    cc.computeCurvature_faces();
    // get the vectors orthogonal to the directions of the smallest absolute curvature value
    orth_smallest_curvature_and_c2(cc.curv, cc.curvDir, idspos, idsneg, ortho, coscos);
    assert(idspos.size() == ortho.size() && idsneg.size() == coscos.size());
    CurvDir = cc.curvDir;
}
void lsTools::show_minimal_curvature_directions(Eigen::MatrixXd& E0, Eigen::MatrixXd& E1, const double scaling){
    std::vector<int> idspos;
    std::vector<int> idsneg;
    std::vector<Eigen::Vector3d> ortho;
    std::vector<double> coscos;
    std::vector<std::vector<Eigen::Vector3d>> CurvDir;
    get_orthogonal_direction_minimal_principle_curvature(V, F, idspos, idsneg, ortho, coscos, CurvDir);
    int fnbr = F.rows();
    E0.resize(fnbr + idsneg.size(), 3);
    E1.resize(fnbr + idsneg.size(), 3);
    for (int i = 0; i < idspos.size(); i++)
    {
        int fid = idspos[i];
        Eigen::Vector3d v0 = V.row(F(fid, 0));
        Eigen::Vector3d v1 = V.row(F(fid, 1));
        Eigen::Vector3d v2 = V.row(F(fid, 2));
        Eigen::Vector3d vm = (v0 + v1 + v2) / 3;
        Eigen::Vector3d direction = ortho[i].cross(Eigen::Vector3d(norm_f.row(fid))).normalized();
        E0.row(fid) = vm + direction * scaling;
        E1.row(fid) = vm - direction * scaling;
    }
    for (int itr = 0; itr < idsneg.size(); itr++)
    {
        int fid = idsneg[itr];
        Eigen::Vector3d v0 = V.row(F(fid, 0));
        Eigen::Vector3d v1 = V.row(F(fid, 1));
        Eigen::Vector3d v2 = V.row(F(fid, 2));
        Eigen::Vector3d vm = (v0 + v1 + v2) / 3;
        double cos_angle = sqrt(coscos[itr]);
        double sin_angle = sqrt(1 - coscos[itr]);
        Eigen::Vector3d axis1 = CurvDir[fid][0].normalized();
        Eigen::Vector3d axis2 = CurvDir[fid][1].normalized();
        Eigen::Vector3d direction1 = axis1 * cos_angle + axis2 * sin_angle;
        Eigen::Vector3d direction2 = -axis1 * cos_angle + axis2 * sin_angle;
        E0.row(fid) = vm + direction1 * scaling;
        E1.row(fid) = vm - direction1 * scaling;
        
        E0.row(fnbr + itr) = vm + direction2 * scaling;
        E1.row(fnbr + itr) = vm - direction2 * scaling;
        
    }

}

void get_overlap_ver_correspondance(CGMesh& mesh, std::vector<std::array<int, 2>>& list)
{   
    lsTools tools;
    tools.init(mesh);
    int vnbr = tools.V.rows();
    std::vector<int> binv; // boundary vertices in vertex list
    binv.resize(vnbr, 0);
    std::vector<int> bvers; // boundary vertices
    bvers.reserve(vnbr);

    for (CGMesh::HalfedgeHandle hd : tools.Boundary_Edges)
    {
        int vfr = tools.lsmesh.from_vertex_handle(hd).idx();
        int vto = tools.lsmesh.to_vertex_handle(hd).idx();
        binv[vfr] = 1;
        binv[vto] = 1;
    }
    for (int i = 0; i < vnbr; i++)
    {
        if(binv[i])
        {
            bvers.push_back(i);
        }
    }
    int bsize = bvers.size();
    std::vector<bool> checked(bsize, false);
    std::vector<std::array<int, 2>> pairs;
    pairs.reserve(100);
    for (int i = 0; i < bsize; i++)
    {
        if (checked[i])
        {
            continue;
        }
        for (int j = 0; j < bsize; j++)
        {
            if (i == j || checked[j])
            {
                continue;
            }
            Eigen::Vector3d p0 = tools.V.row(bvers[i]);
            Eigen::Vector3d p1 = tools.V.row(bvers[j]);
            double distance = (p0 - p1).norm();
            if (distance < 1e-6)
            {
                checked[i] = true;
                checked[j] = true;
                pairs.push_back({bvers[i],bvers[j]});
            }
        }
    }
    list = pairs;
    std::cout << "the nbr of overlap boundary vertex pairx: " << list.size() << std::endl;
}

void obj2csv()
{
    std::cout << "Reading curve obj file" << std::endl;
    std::string fname = igl::file_dialog_open();
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::readOBJ(fname, V, F);
    fname = igl::file_dialog_save();
    std::ofstream fout;
    fout.open(fname);
    for (int i = 0; i < V.rows(); i++)
    {
        Eigen::Vector3d xyz = V.row(i);
        fout << xyz[0] << "," << xyz[1] << "," << xyz[2] << "\n";
    }
    fout.close();
}
void csv2objcurves(){
    std::cout<<"reading the 1 family of curves, please provide the prefixes"<<std::endl;
    std::string fname1 = igl::file_dialog_save();
    // std::string fname2 = igl::file_dialog_save();
    if (fname1.length() == 0)
    {
        std::cout << "Please type down the prefix" << std::endl;
        return;
    }
    std::vector<std::vector<double>> vx1, vy1, vz1;
    read_csv_data_lbl(fname1+"_x.csv", vx1);
    read_csv_data_lbl(fname1+"_y.csv", vy1);
    read_csv_data_lbl(fname1+"_z.csv", vz1);
    std::cout << "saving obj edge file" << std::endl;
    std::string fname = igl::file_dialog_save();
    std::ofstream fout;
    fout.open(fname);
    
    for (int i = 0; i < vx1.size(); i++)
    {
        for (int j = 0; j < vx1[i].size(); j++)
        {
            fout << "v " << vx1[i][j] << " " << vy1[i][j] << " " << vz1[i][j] << "\n";
        }
    }
    int count = 1; // id starts from 1
    for (int i = 0; i < vx1.size(); i++)
    {
        for (int j = 0; j < vx1[i].size(); j++)
        {
            if (j < vx1[i].size() - 1)
            {
                fout << "l " << count << " " << count + 1 << "\n";
            }
            count++;
        }
    }
    fout.close();
    std::cout << "file saved" << std::endl;
}

void lsTools::debug_tool(){
    Eigen::VectorXi info;

    orthogonal_slope_for_different_shading_types(1, InnerV, info, V, F, norm_f, 0, 180, 10,
                                                 30, 0, 180, Fslope, OrthoSlope);
}

bool project_1_plyline_on_1_curve(const std::vector<Eigen::Vector3d>& ply, const std::vector<Eigen::Vector3d>& curve, 
    Eigen::Vector3d& proj){
    
    double dis_min = 1e6;
    double clsu, clsv;
    Eigen::Vector3d closest;
    for (int i = 0; i < ply.size() - 1; i++)
    {
        Eigen::Vector3d v0s = ply[i];
        Eigen::Vector3d v0e = ply[i + 1];
        Eigen::Vector3d dir0 = (v0e - v0s);
        for (int j = 0; j < curve.size() - 1; j++)
        {
            Eigen::Vector3d v1s = curve[j];
            Eigen::Vector3d v1e = curve[j + 1];
            Eigen::Vector3d dir1 = (v1e - v1s);

            double u, v;
            igl::segment_segment_intersect(v0s, dir0, v1s, dir1, u, v);
            if (u > 1)
            {
                u = 1;
            }
            if (u < 0)
            {
                u = 0;
            }
            if (v > 1)
            {
                v = 1;
            }
            if (v < 0)
            {
                v = 0;
            }
            Eigen::Vector3d p0 = v0s + dir0 * u;
            Eigen::Vector3d p1 = v1s + dir1 * v;
            double distance = (p0 - p1).norm();
            if (dis_min > distance)
            {
                dis_min = distance;
                closest = p1;
            }
            // the four end points
            p0 = v0s;
            p1 = v1s;
            distance = (p0 - p1).norm();
            if (dis_min > distance)
            {
                dis_min = distance;
                closest = p1;
            }

            p0 = v0e;
            p1 = v1s;
            distance = (p0 - p1).norm();
            if (dis_min > distance)
            {
                dis_min = distance;
                closest = p1;
            }

            p0 = v0s;
            p1 = v1e;
            distance = (p0 - p1).norm();
            if (dis_min > distance)
            {
                dis_min = distance;
                closest = p1;
            }
            p0 = v0e;
            p1 = v1e;
            distance = (p0 - p1).norm();
            if (dis_min > distance)
            {
                dis_min = distance;
                closest = p1;
            }
        }
    }
    double threads = 1e-3;
    bool c1 = (closest - curve.front()).norm() < threads;
    bool c2 = (closest - curve.back()).norm() < threads;
    bool c3 = (closest - ply.front()).norm() < threads;
    bool c4 = (closest - ply.back()).norm() < threads;
    if (c1 || c2 || c3 || c4)
    {
        // std::cout<<"wrong pt: cloest: "<<closest.transpose()<<"\nv1, "<<curve.front().transpose()
        // <<"\nv2, "<<curve.back().transpose()
        // <<"\nv3, "<<ply.front().transpose()
        // <<"\nv4, "<<ply.back().transpose()<<std::endl;

        return false;
    }
    proj = closest;
    return true;
}

// std::vector<Eigen::Vector3d> filter_projection_pts(const std::vector<std::vector<Eigen::Vector3d>> &plys,
//                                                    const std::vector<Eigen::Vector3d> &pts)
// {
//     Eigen::Vector3d vmin = plys[0][0], vmax = plys[0][0];
//     for (auto ply : plys)
//     {
//         for (auto p : ply)
//         {
//             if (p[0] < vmin[0])
//             {
//                 vmin[0] = p[0];
//             }
//             if (p[1] < vmin[1])
//             {
//                 vmin[1] = p[1];
//             }
//             if (p[2] < vmin[2])
//             {
//                 vmin[2] = p[2];
//             }
//             if (p[0] > vmax[0])
//             {
//                 vmax[0] = p[0];
//             }
//             if (p[1] > vmax[1])
//             {
//                 vmax[1] = p[1];
//             }
//             if (p[2] > vmax[2])
//             {
//                 vmax[2] = p[2];
//             }
//         }
//     }
//     double diag = (vmax - vmin).norm();
//     double threadshold = 0.1;
//     std::vector<Eigen::Vector3d> result;
//     for(int i=0;i<)
// }

void project_polylines_on_shading_curves_and_save_results()
{
    std::vector<std::vector<Eigen::Vector3d>> plys, curs, pros;
    std::cout << "Reading orthogonal curves first, then reading reference shading curves" << std::endl;
    read_polylines(plys);
    read_polylines(curs);
    if (plys.size() == 0 || curs.size() == 0)
    {
        return;
    }

    for (int i = 0; i < plys.size(); i++)
    {
        std::cout<<"Processing ply "<<i<<std::endl;
        std::vector<Eigen::Vector3d> pts;
        for (int j = 0; j < curs.size(); j++)
        {
            Eigen::Vector3d p;
            bool intersect =
                project_1_plyline_on_1_curve(plys[i], curs[j], p);
            if (intersect)
            {
                pts.push_back(p);
                // std::cout << "find one intersection point " << std::endl;
            }
            else
            {
                // std::cout << "ply " << i << " and curve " <<j<<" has no closest point"<<std::endl;
            }
        }
        if (!pts.empty())
        {
            pros.push_back(pts);
        }

        else
        {
            std::cout << "Ply " << i << " has 0 curves " << std::endl;
        }
    }
    // write_polyline_xyz()
    std::cout << "Saving the projections" << std::endl;
    std::string fname = igl::file_dialog_save();
    if (fname.length() == 0)
    {
        std::cout << "Please type down the prefix" << std::endl;
        return;
    }
    for (int i = 0; i < pros.size(); i++)
    {
        if (pros[i].size() > 0)
        {
            Eigen::MatrixXd V = vec_list_to_matrix(pros[i]);
            Eigen::MatrixXi F;
            igl::writeOBJ(fname + std::to_string(i) + ".obj", V, F);
            std::cout << fname + std::to_string(i) + ".obj"
                      << " got written" << std::endl;
        }
    }
    std::cout << "All the files are saved" << std::endl;
}

void read_draw_pts_from_plylines(Eigen::MatrixXd &ver)
{
    std::vector<std::vector<Eigen::Vector3d>> plys;
    read_polylines(plys);
    std::vector<Eigen::Vector3d> vlist;
    for (auto line : plys)
    {
        for (auto pt : line)
        {
            vlist.push_back(pt);
        }
    }
    ver = vec_list_to_matrix(vlist);
}

// take polylines as input, offset to both sides of the creases, and save as mesh
void read_plylines_extract_offset_mesh(const double scale_front, const double scale_back, CGMesh& mesh){
    std::vector<std::vector<Eigen::Vector3d>> plys, bnms;
    read_plylines_and_binormals(plys, bnms);
    mesh = polyline_to_strip_mesh(plys, bnms, scale_front, scale_back);
}

void make_example_comparing_two_plylines_distance(){
    Eigen::MatrixXd ply_move, ply_fix;
    Eigen::MatrixXi F;
    igl::readOBJ("/Users/wangb0d/Desktop/tmp/angle/cylinder_p/test/continuous.obj", ply_move, F);
    igl::readOBJ("/Users/wangb0d/Desktop/tmp/angle/cylinder_p/test/discrete.obj", ply_fix, F);

    std::vector<Eigen::Vector3d> pm = mat_to_vec_list(ply_move), pf = mat_to_vec_list(ply_fix);
    double dis;
    ply_ply_distance(pm, pf, dis);
    std::cout<<"maximal distance is "<<dis<<std::endl;
}

void compare_multiple_correspoding_curves_dis(){
    std::cout<<"reading the two families of curves, please provide the prefixes"<<std::endl;
    std::string fname1 = igl::file_dialog_save();
    std::string fname2 = igl::file_dialog_save();
    if (fname1.length() == 0 || fname2.length() == 0)
    {
        std::cout << "Please type down the prefix" << std::endl;
        return;
    }
    std::vector<std::vector<double>> vx2, vy2, vz2, vx1, vy1, vz1;
    read_csv_data_lbl(fname1+"_x.csv", vx1);
    read_csv_data_lbl(fname1+"_y.csv", vy1);
    read_csv_data_lbl(fname1+"_z.csv", vz1);

    read_csv_data_lbl(fname2+"_x.csv", vx2);
    read_csv_data_lbl(fname2+"_y.csv", vy2);
    read_csv_data_lbl(fname2 + "_z.csv", vz2);
    assert(vx1.size() == vx2.size());
    int csize = vx2.size();
    std::vector<std::vector<Eigen::Vector3d>> c1s(csize), c2s(csize);
    for (int i = 0; i < csize; i++)
    {
        c1s[i].resize(vx1[i].size());
        c2s[i].resize(vx2[i].size());
        for (int j = 0; j < c1s[i].size(); j++)
        {
            c1s[i][j] = Eigen::Vector3d(vx1[i][j], vy1[i][j], vz1[i][j]);
        }
        for (int j = 0; j < c2s[i].size(); j++)
        {
            c2s[i][j] = Eigen::Vector3d(vx2[i][j], vy2[i][j], vz2[i][j]);
        }
    }


    double max_dis = -1;
    for (int i = 0; i < csize; i++)
    {
        double cdis = -1;
        ply_ply_distance(c1s[i], c2s[i], cdis);
        std::cout << i << "th, dis, " << cdis << std::endl;
        if (cdis > max_dis)
            max_dis = cdis;
    }
    std::cout<<"absolute max distance between two families of segs: "<<max_dis<<std::endl;
}


std::vector<int> sort_plylines_based_on_midpts(const std::vector<std::vector<Eigen::Vector3d>> &ply)
{
    int cnbr = ply.size();
    std::vector<std::pair<Eigen::Vector3d, int>> pairs(cnbr);
    double largest_dis = 0;
    double pid1, pid2;
    // get the mid pts to represent the lines
    std::vector<Eigen::Vector3d> midpts(cnbr);
    for (int i = 0; i < cnbr; i++)
    {
        int pnbr = ply[i].size();
        Eigen::Vector3d midpt = ply[i][pnbr / 2];
        midpts[i] = midpt;
    }
    // get the two corner polylines.
    for (int i = 0; i < cnbr; i++)
    {
        for (int j = i; j < cnbr; j++)
        {
            double distance = (midpts[i] - midpts[j]).norm();
            if (distance > largest_dis)
            {
                largest_dis = distance;
                pid1 = i;
                pid2 = j;
            }
        }
    }
    // distances from pid1
    std::vector<double> dislist(cnbr);
    for(int i=0;i<cnbr;i++)
    {
        dislist[i] = (midpts[i] - midpts[pid1]).norm();
    }
    std::vector<int> checked(cnbr, false);
    checked[pid1] = true;
    int pushed = 0;
    std::vector<int> order;
    order.push_back(pid1);
    for (int i = 0; i < cnbr; i++)
    {
        double mindis = 1e3;
        int minid = -1;
        for (int j = 0; j < cnbr; j++)
        {
            if (checked[j])
            {
                continue;
            }
            double dis = dislist[j];
            if (dis < mindis)
            {
                mindis = dis;
                minid = j;
            }
        }
        if(minid<0){
            break;
        }
        order.push_back(minid);
        checked[minid] = true;
    }
    return order;
}

void run_sort_polylines()
{
    std::cout<<"sorting the polylines"<<std::endl;
    std::vector<std::vector<Eigen::Vector3d>> ply, bnm;
    read_plylines_and_binormals(ply, bnm);
    std::cout << "Saving the binormal files, please provide the prefix" << std::endl;
    std::string fname = igl::file_dialog_save();
    if (fname.length() == 0)
    {
        std::cout << "\nLSC: save mesh failed, please type down the correct name" << std::endl;
        return;
    }
    std::vector<int> order = sort_plylines_based_on_midpts(ply);

    std::vector<std::vector<Eigen::Vector3d>> lines, binormals;
    for (int i = 0; i < ply.size(); i++)
    {
        lines.push_back(ply[order[i]]);
        binormals.push_back(bnm[order[i]]);
    }

    write_polyline_xyz(lines, fname);
    write_polyline_xyz(binormals, fname + "_b");
    std::cout << "files get saved" << std::endl;
}

void show_rotback_slopes(const Eigen::MatrixXd &E0, const Eigen::MatrixXd &E1, const double latitude_degree,
                         Eigen::MatrixXd &rotmids, Eigen::MatrixXd &rotdirs)
{
    double latitude_radian = latitude_degree * LSC_PI / 180.;
    double rho = -(LSC_PI / 2 - latitude_radian); // this is correct rotation
    Eigen::MatrixXd midpts = (E0 + E1) / 2;
    Eigen::MatrixXd directions = E1 - E0;
    for (int i = 0; i < directions.rows(); i++)
    {
        Eigen::Vector3d dir = directions.row(i);
        directions.row(i) = dir.normalized();
    }

    Eigen::Matrix3d rotation;
    rotation << 1, 0, 0,
        0, cos(rho), -sin(rho),
        0, sin(rho), cos(rho);
    // assert(ply_extracted.size() > 0);
    
    rotmids = (rotation * midpts.transpose()).transpose();
    rotdirs = (rotation * directions.transpose()).transpose();
}

void save_rotback_slopes(const Eigen::MatrixXd &rotmids, const Eigen::MatrixXd &rotdirs){
    std::cout<<"saving the slopes as csv files, type down the prefix"<<std::endl;
    std::string fname = igl::file_dialog_save();
    if (fname.length() == 0 || rotmids.rows()==0)
    {
        std::cout << "\nLSC: save mesh failed, please type down the correct name, or feed the correct data" << std::endl;
        return;
    }
    std::ofstream file;
    file.open(fname + ".csv");
    for (int i = 0; i < rotmids.rows(); i++)
    {
        file << rotmids(i, 0) << "," << rotmids(i, 1) << "," << rotmids(i, 2) << "," << rotdirs(i, 0) << "," << rotdirs(i, 1) << "," << rotdirs(i, 2) << "\n";

    }
    file.close();
    std::cout<<"file saved"<<std::endl;

}
void decrease_ply_nbr_by_half(){
    std::vector<std::vector<Eigen::Vector3d>> ply, bnm, plyout, bnmout;
    read_plylines_and_binormals(ply, bnm);
    int count = 1;
    for (int i = 0; i < ply.size(); i++)
    {
        if (count < ply.size())
        {
            plyout.push_back(ply[count]);
            bnmout.push_back(bnm[count]);
            count += 2;
            std::cout<<"select nbr "<<count;
        }
        else
        {
            break;
        }
    }
    std::string fname = igl::file_dialog_save();
    if (fname.length() == 0)
    {
        std::cout << "\nLSC: save mesh failed, please type down the correct name" << std::endl;
        return;
    }

    write_polyline_xyz(plyout, fname);
    write_polyline_xyz(bnmout, fname + "_b");
    std::cout << "files get saved" << std::endl;
}

void save_invert_levelset(const Eigen::VectorXd& func){
    save_levelset(-func);
}

void write_polyline_joints(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, const Eigen::MatrixXd &norm_v){
    std::vector<std::vector<Eigen::Vector3d>> ply1, ply2, bin1, bin2;
    read_plylines_and_binormals(ply1, bin1);
    read_plylines_and_binormals(ply2, bin2);
    std::vector<std::vector<std::vector<int>>> correspondance;
    std::vector<Eigen::Vector3d> pts, dir0s, dir1s;
    double threads = 0.5;
    double normlength = 2;// each normal length
    correspondance.resize(ply1.size());

    for (int i = 0; i < ply1.size(); i++)
    {
        correspondance[i].resize(ply1[i].size());
        for (int j = 1; j < ply1[i].size() - 1; j++)
        {
            Eigen::Vector3d p = ply1[i][j];
            double mindis = 1e10;
            int kmin = 0;
            int lmin = 0;
            for (int k = 0; k < ply2.size(); k++)
            {
                for (int l = 1; l < ply2[k].size() - 1; l++)
                {
                    Eigen::Vector3d pc = ply2[k][l];
                    double dis = (pc - p).norm();
                    if (dis < mindis)
                    {
                        mindis = dis;
                        kmin = k;
                        lmin = l;
                    }
                }
            }
            if (lmin == 0 || lmin == ply2[kmin].size() - 1)
            {
                continue;
            }
            if (mindis > threads)
            {
                continue;
            }
            else
            {
                correspondance[i][j].push_back(kmin);
                correspondance[i][j].push_back(lmin);
                pts.push_back(p);
                dir0s.push_back(bin1[i][j]);
                if (bin1[i][j].dot(bin2[kmin][lmin]) > 0)// make sure they are in opposite direction
                {
                    dir1s.push_back(bin2[kmin][lmin]);
                }
                else
                {
                    dir1s.push_back(-bin2[kmin][lmin]);
                }
            }
        }
    }
    Eigen::MatrixXd C;
    Eigen::VectorXi I;
    Eigen::VectorXd D;
    int nq = pts.size();
    Eigen::MatrixXd norms;
    norms.resize(nq, 3);
    Eigen::MatrixXd vers = vec_list_to_matrix(pts);
    igl::point_mesh_squared_distance(vers, V, F, D, I, C);
    for (int i = 0; i < nq; i++)
    {
        int fid = I(i);
        int v0 = F(fid, 0);
        int v1 = F(fid, 1);
        int v2 = F(fid, 2);
        Eigen::Vector3d pt = pts[i];
        Eigen::Vector3d ver0 = V.row(v0);
        Eigen::Vector3d ver1 = V.row(v1);
        Eigen::Vector3d ver2 = V.row(v2);
        std::array<double, 3> coor = barycenter_coordinate(ver0, ver1, ver2, pt);

        Eigen::Vector3d dir = coor[0] * norm_v.row(v0) + coor[1] * norm_v.row(v1) + coor[2] * norm_v.row(v2);
        dir = dir.normalized();
        if(dir == Eigen::Vector3d(0,0,0)){
            dir = Eigen::Vector3d(1,0,0);
        }
        norms.row(i) = dir;
    }
    std::vector<Eigen::Vector3d> pout;
    std::vector<std::vector<int>> lines(nq), tris;
    std::string vwrite;
    for (int i = 0; i < nq; i++)
    {
        // E0 E1 are the end points of the axis
        Eigen::Vector3d E0 = pts[i] + normlength / 2 * Eigen::Vector3d(norms.row(i));
        Eigen::Vector3d E1 = pts[i] - normlength / 2 * Eigen::Vector3d(norms.row(i));
        // E2 is on E0 side.
        Eigen::Vector3d E2, E3, T2, T3;
        Eigen::Vector3d c1 = dir0s[i];
        Eigen::Vector3d c2 = dir1s[i];
        Eigen::Vector3d n = norms.row(i);
        // double alpha1 = normlength / (2 * c1.dot(n));
        // double alpha2 = normlength / (2 * c2.dot(n));
        double alpha1 = normlength / 2 / c1.norm();
        double alpha2 = normlength / 2 / c2.norm();
        if (alpha1 > 0)
        { //c1 is on the positive direction wrt the normal vector
            E2 = pts[i] + alpha1 * c1;
            E3 = pts[i] - alpha2 * c2;
        }
        else{
            // c1 is on the negative direction 
            E2 = pts[i] + alpha2 * c2;
            E3 = pts[i] - alpha1 * c1;
        }
        int le0 = pout.size() + 1;
        int le1 = pout.size() + 2;
        int le2 = pout.size() + 3;
        int le3 = pout.size() + 4;
        int le4 = pout.size() + 5; // this is the middle ver
        // lines
        lines[i].push_back(le0);
        lines[i].push_back(le1);

        // triangles
        std::vector<int> tempt;
        tempt.push_back(le4);
        tempt.push_back(le2);
        tempt.push_back(le0);
        tris.push_back(tempt);
        tempt.clear();
        tempt.push_back(le4);
        tempt.push_back(le3);
        tempt.push_back(le1);
        tris.push_back(tempt);

        // points;
        pout.push_back(E0);
        pout.push_back(E1);
        pout.push_back(E2);
        pout.push_back(E3);
        pout.push_back(pts[i]);
        // check the angles
        double angle1 = acos(c1.normalized().dot(n)) * 180 / LSC_PI;
        double angle2 = acos(c2.normalized().dot(n)) * 180 / LSC_PI;
        std::cout << "Angles, " << angle1 << ", " << angle2 << std::endl;
    }
    for (int i = 0; i < pout.size(); i++)
    {
        Eigen::Vector3d pt = pout[i];
        vwrite = vwrite + "v " + std::to_string(pt[0]) + " " + std::to_string(pt[1]) + " " + std::to_string(pt[2]) + "\n";
    }
    std::cout << "Please write down the prefix" << std::endl;
    std::string fname = igl::file_dialog_save();
    if (fname.length() == 0)
    {
        std::cout << "\nLSC: save mesh failed, please type down the correct name" << std::endl;
        return;
    }
    std::ofstream fout;
    fout.open(fname + "l.obj");
    fout << vwrite;
    for (int i = 0; i < lines.size(); i++)
    {
        fout << "l " << lines[i][0] << " " << lines[i][1] << std::endl;
    }
    fout.close();

    fout.open(fname + "tri.obj");
    fout << vwrite;
    for (int i = 0; i < tris.size(); i++)
    {
        fout << "f " << tris[i][0] << " " << tris[i][1] << " " << tris[i][2] << std::endl;
    }
    fout.close();
    std::cout << "finished" << std::endl;
}

// Eigen::Vector3d eq35(double u, double v){
//     double sru = sqrt(u);
//     double x = exp((1-u)*v)*sru*(sru*sin())
// }

// void write_functional_surface(double radius, double height, int nr, int nh)
// {
//     Eigen::MatrixXd ver;
//     Eigen::MatrixXi faces;
//     ver.resize(nh * nr, 3);
//     faces.resize(2 * (nh - 1) * (nr-1), 3);
//     int verline = 0;
//     double hitv = height / (nh - 1);
//     double ritv =  LSC_PI / (nr-1);// we only draw half a cylinder
//     for (int i = 0; i < nr; i++)
//     {
//         double angle = ritv * i;
//         double x = cos(angle) * radius;
//         double y = sin(angle) * radius;
//         for (int j = 0; j < nh; j++)
//         {
//             double z = j * hitv;
//             ver.row(verline) << x, y, z;
//             verline++;
//         }
//     }

//     int fline = 0;
//     for (int i = 0; i < nr; i++)
//     {
//         if (i < nr - 1)
//         {
//             for (int j = 0; j < nh - 1; j++)
//             {
//                 int id0 = nh * i + j;
//                 int id1 = nh * (i + 1) + j;
//                 int id2 = nh * (i + 1) + j + 1;
//                 int id3 = nh * i + j + 1;
//                 faces.row(fline) = Eigen::Vector3i(id0, id1, id2);
//                 faces.row(fline + 1) = Eigen::Vector3i(id0, id2, id3);
//                 fline += 2;
//             }
//         }
//         // else
//         // {
//         //     for (int j = 0; j < nh - 1; j++)
//         //     {
//         //         int id0 = nh * i + j;
//         //         int id1 = j;
//         //         int id2 = j + 1;
//         //         int id3 = nh * i + j + 1;
//         //         faces.row(fline) = Eigen::Vector3i(id0, id1, id2);
//         //         faces.row(fline + 1) = Eigen::Vector3i(id0, id2, id3);
//         //         fline += 2;
//         //     }
//         // }
//     }
//     std::string path("/Users/wangb0d/bolun/D/vs/levelset/level-set-curves/data/");
//     igl::write_triangle_mesh(path + "oc_" + std::to_string(radius) + "_" + std::to_string(height) + "_" +
//                                  std::to_string(nr) + "_" + std::to_string(nh) + ".obj",
//                              ver, faces);
//     std::cout << "cylinder file saved " << std::endl;
// }



void generate_rulings_outof_plylines()
{
    double scale = 0.5;
    std::vector<std::vector<Eigen::Vector3d>> ply1, bin1;
    read_plylines_and_binormals(ply1, bin1);
    std::vector<Eigen::Vector3d> pts;
    std::vector<std::vector<int>> lines;
    
    for (int i = 0; i < ply1.size(); i++)
    {
        for (int j = 1; j < ply1[i].size() - 2; j++)
        {
            std::vector<int> l;
            l.push_back(pts.size() + 1);
            l.push_back(pts.size() + 2);
            lines.push_back(l);
            Eigen::Vector3d p1 = ply1[i][j] + bin1[i][j] * scale;
            Eigen::Vector3d p2 = ply1[i][j] - bin1[i][j] * scale;

            pts.push_back(p1);
            pts.push_back(p2);
        }
    }
    std::string vwrite;
    for (int i = 0; i < pts.size(); i++)
    {
        Eigen::Vector3d pt = pts[i];
        vwrite = vwrite + "v " + std::to_string(pt[0]) + " " + std::to_string(pt[1]) + " " + std::to_string(pt[2]) + "\n";
    }
    std::cout << "Please write down the prefix" << std::endl;
    std::string fname = igl::file_dialog_save();
    if (fname.length() == 0)
    {
        std::cout << "\nLSC: save mesh failed, please type down the correct name" << std::endl;
        return;
    }
    std::ofstream fout;
    fout.open(fname + "l.obj");
    fout << vwrite;
    for (int i = 0; i < lines.size(); i++)
    {
        fout << "l " << lines[i][0] << " " << lines[i][1] << std::endl;
    }
    fout.close();
}


// threadshold_nbr is the nbr of quads we want to discard when too few quads in a row
void compute_ls_intersections_and_assign_parallel_joints(const CGMesh &lsmesh, const Eigen::Vector3d &joint,
                                                         const double scale, const Eigen::MatrixXd &V,
                                                         const Eigen::MatrixXi &F, const Eigen::VectorXd &ls0, const Eigen::VectorXd &ls1,
                                                         const int expect_nbr_ls0, const int expect_nbr_ls1,
                                                         bool even_pace = false)
{
    std::vector<std::vector<Eigen::Vector3d>> ivs;                                    // the intersection vertices
    Eigen::MatrixXi gridmat;                                                          // the matrix for vertex
    Eigen::VectorXd lsv0, lsv1;                                                       // the extracted level set values;
    std::vector<std::vector<Eigen::Vector3d>> poly0_e0, poly0_e1, poly1_e0, poly1_e1; // polylines
    std::vector<std::vector<int>> fid0, fid1;                                         // face ids of each polyline segment
    std::vector<Eigen::Vector3d> verlist, jlist;
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
    std::vector<std::vector<int>> lines;
    for (int i = 0; i < nbr_ls0; i++)
    {
        for (int j = 0; j < nbr_ls1; j++)
        {
            Eigen::Vector3d ipoint;
            bool intersect = get_polyline_intersection(V, F, poly0_e0[i], poly0_e1[i], poly1_e0[j], poly1_e1[j], fid0[i], fid1[j], ipoint);
            if (intersect)
            {
                // double dis =
                verlist.push_back(ipoint);
                gridmat(i, j) = vnbr;
                std::vector<int> ids;
                ids.push_back(jlist.size() + 1);
                ids.push_back(jlist.size() + 2);
                lines.push_back(ids);
                jlist.push_back(ipoint + scale * joint);
                jlist.push_back(ipoint - scale * joint);
                
                vnbr++;
            }
        }
        // The following pars are for debugging
        // Eigen::VectorXi cv=gridmat.row(i);
        // std::vector<int> problems = element_ids_for_not_computed_vers(cv);
        // if(problems.size()>0){
        //     std::cout<<i<<" this row has missing points \n"<<cv.transpose()<<" size, "<<problems.size()<<std::endl;
            
        //     for(int j: problems){
        //         Eigen::Vector3d ipoint;
        //         get_polyline_intersection(V, F, poly0_e0[i], poly0_e1[i], poly1_e0[j], poly1_e1[j], fid0[i], fid1[j], ipoint, true);
        //     }
        // }
    }
    std::string vwrite;
    for (int i = 0; i < jlist.size(); i++)
    {
        Eigen::Vector3d pt = jlist[i];
        vwrite = vwrite + "v " + std::to_string(pt[0]) + " " + std::to_string(pt[1]) + " " + std::to_string(pt[2]) + "\n";
    }
    std::cout << "write the joints. Please write down the prefix" << std::endl;
    std::string fname = igl::file_dialog_save();
    if (fname.length() == 0)
    {
        std::cout << "\nLSC: save mesh failed, please type down the correct name" << std::endl;
        return;
    }
    std::ofstream fout;
    fout.open(fname + "l.obj");
    fout << vwrite;
    for (int i = 0; i < lines.size(); i++)
    {
        fout << "l " << lines[i][0] << " " << lines[i][1] << std::endl;
    }
    fout.close();
    std::cout<<"writing done. "<<std::endl;
    
}
void assign_const_rulings_to_strips(const Eigen::Vector3d& bnm){
    std::vector<std::vector<Eigen::Vector3d>> ply1, bin1;
    read_polylines(ply1);
    bin1 = ply1;
    // read_plylines_and_binormals(ply2, bin2);
    for (int i = 0; i < ply1.size(); i++)
    {
        for (int j = 0; j < ply1[i].size(); j++)
        {
            bin1[i][j] = bnm;
        }
    }
    std::cout << "write prefix for saving new files" << std::endl;
    std::string fname = igl::file_dialog_save();
    if (fname.length() == 0)
    {
        std::cout << "\nLSC: save mesh failed, please type down the correct name" << std::endl;
        return;
    }

    write_polyline_xyz(ply1, fname);
    write_polyline_xyz(bin1, fname + "_b");
    std::cout << "files get saved" << std::endl;
}

void assign_ls_to_subpatch(CGMesh &ref, CGMesh &base, const Eigen::VectorXd &func, Eigen::VectorXd &funout)
{
    lsTools reftool(ref);
    lsTools bastool(base);
    reftool.initialize_mesh_properties();
    bastool.initialize_mesh_properties();
    int rvnbr = reftool.V.rows();
    int bvnbr = bastool.V.rows();
    int ninner = bastool.IVids.size();
    Eigen::VectorXi InnerV = bastool.InnerV;
    Eigen::MatrixXd C;
    Eigen::VectorXi I;
    Eigen::VectorXd D;

    funout.resize(rvnbr);

    igl::point_mesh_squared_distance(reftool.V, bastool.V, bastool.F, D, I, C);
    for (int i = 0; i < rvnbr; i++)
    {
        int f = I(i);

        Eigen::Vector3d p = C.row(i);
        Eigen::Vector3d v0 = bastool.V.row(bastool.F(f, 0));
        Eigen::Vector3d v1 = bastool.V.row(bastool.F(f, 1));
        Eigen::Vector3d v2 = bastool.V.row(bastool.F(f, 2));
        std::array<double, 3> coor = barycenter_coordinate(v0, v1, v2, p);
        double fun0 = func[bastool.F(f, 0)];
        double fun1 = func[bastool.F(f, 1)];
        double fun2 = func[bastool.F(f, 2)];

        double value = coor[0] * fun0 + coor[1] * fun1 + coor[2] * fun2;

        funout[i] = value;
        
    }
}


// the following code is used to generate arrows aligning with given directions.
// merge the two mesh into one file
void extendMesh(Eigen::MatrixXd& V, Eigen::MatrixXi& F, Eigen::MatrixXd& Va, Eigen::MatrixXi& Fa){
    int vnbr0 = V.rows();
    int vnbr1 = Va.rows();
    int fnbr0 = F.rows();
    int fnbr1 = Fa.rows();
    Eigen::MatrixXd Vnew(vnbr0 + vnbr1, 3);
    Eigen::MatrixXi Fnew(fnbr0 + fnbr1, 3);
    Vnew.topRows(vnbr0) = V;
    Vnew.bottomRows(vnbr1) = Va;
    Fnew.topRows(fnbr0) = F;
    Fnew.bottomRows(fnbr1) = Fa;
    Fnew.bottomRows(fnbr1).array() += vnbr0;
    V = Vnew;
    F = Fnew;
}

// void RotationMatrix(const Eigen::Vector3d& vbefore,const Eigen::Vector3d& vafter, Eigen::Matrix3d& mat){
//     Eigen::Vector3d va = vbefore.normalized();
//     Eigen::Vector3d vb = vafter.normalized();
//     Eigen::Vector3d vs = vb.cross(va);
//     Eigen::Vector3d v = vs.normalized();
//     double ca = vb.dot(va);
//     double scale = 1 - ca;
//     Eigen::Vector3d vt = v * scale;
//     mat(0,0) = vt[0] * v[0] + ca;
//     mat(1,1) = vt[1] * v[1] + ca;
//     mat(2,2) = vt[2] * v[2] + ca;
//     vt[0] *= v[1];
//     vt[2] *= v[0];
//     vt[1] *= v[2];
//     mat(0, 1) = vt[0] - vs[2];
//     mat(0, 2) = vt[2] - vs[1];
//     mat(1, 0) = vt[0] - vs[2];
//     mat(1, 2) = vt[1] - vs[0];
//     mat(2, 0) = vt[2] - vs[1];
//     mat(2, 1) = vt[1] - vs[0];
// }

// Rodrigues' rotation formula1: rotate a vector along a random vector
void getRotationVector(const Eigen::Vector3d &vectorRef, const double angleDegree,
                                  const Eigen::Vector3d &vectorIn, Eigen::Vector3d &vectorOut)
{
    const auto &k = vectorRef;
    const auto &v = vectorIn;
    const double angle_radian = angleDegree * LSC_PI / 180;
    vectorOut = v * cos(angle_radian) + (1 - cos(angle_radian)) * v.dot(k) * k + sin(angle_radian) * v.cross(k);
    return;
}



// Rodrigues' rotation formula2: rotate a vector to the target direction
Eigen::Matrix3d getRotationMatrix(Eigen::Vector3d vectorBefore, Eigen::Vector3d vectorAfter)
        {
            Eigen::Vector3d vb = vectorBefore.normalized();
            Eigen::Vector3d va = vectorAfter.normalized();

            Eigen::Vector3d vs = vb.cross(va);
            Eigen::Vector3d v = vs.normalized();
            double ca = vb.dot(va);

            double scale = 1 - ca;
            Eigen::Vector3d vt(v[0] * scale, v[1] * scale, v[2] * scale);

            Eigen::Matrix3d rotationMatrix;

            rotationMatrix(0, 0) = vt[0] * v[0] + ca;
            rotationMatrix(1, 1) = vt[1] * v[1] + ca;
            rotationMatrix(2, 2) = vt[2] * v[2] + ca;
            vt[0] *= v[1];
            vt[2] *= v[0];
            vt[1] *= v[2];

            rotationMatrix(0, 1) = vt[0] - vs[2];
            rotationMatrix(0, 2) = vt[2] + vs[1];
            rotationMatrix(1, 0) = vt[0] + vs[2];
            rotationMatrix(1, 2) = vt[1] - vs[0];
            rotationMatrix(2, 0) = vt[2] - vs[1];
            rotationMatrix(2, 1) = vt[1] + vs[0];

            return rotationMatrix;
        }



// the arrow is originally in (0,0,0) and pointing to z axis
void orientMeshToDirection(const Eigen::Vector3d& pstart, const Eigen::Vector3d& pend, Eigen::MatrixXd& V){
    // put the mesh to the location
    int vnbr = V.rows();
    
    // orient the mesh
    Eigen::Matrix3d rot;
    Eigen::Vector3d dtarget = (pend - pstart).normalized();
    Eigen::Vector3d dstart = Eigen::Vector3d(0, 0, 1);
    rot = getRotationMatrix(dstart, dtarget);
    Eigen::MatrixXd Vtra = V.transpose();
    Vtra = rot * Vtra;
    V = Vtra.transpose();
    V.col(0) = V.col(0) + Eigen::VectorXd::Ones(vnbr) * double(pstart(0));
    V.col(1) = V.col(1) + Eigen::VectorXd::Ones(vnbr) * double(pstart(1));
    V.col(2) = V.col(2) + Eigen::VectorXd::Ones(vnbr) * double(pstart(2));

}

void orientVectors(Eigen::MatrixXd &d, const double sign = 1)
{
    Eigen::Vector3d ref = d.row(0) * sign;
    for (int i = 0; i < d.rows(); i++)
    {
        Eigen::Vector3d vec = d.row(i);
        if (vec.dot(ref) < 0)
        {
            d.row(i) *= -1;
        }
    }
}

void orientEndpts(Eigen::MatrixXd& ref, const Eigen::MatrixXd &start, Eigen::MatrixXd &end, const double sign = 1)
{
    ref *= sign;
    for (int i = 0; i < start.rows(); i++)
    {
        Eigen::Vector3d vec = end.row(i) - start.row(i);
        if (vec.dot(Eigen::Vector3d(ref.row(i))) < 0)
        {
            end.row(i) = -vec + Eigen::Vector3d(start.row(i));
        }
    }
}

void runRotateArrows(double scaling){
    // for normal vectors
    bool orient = false;
    double orientSign = 1;
    std::string normEndfile = "/Users/wangb0d/bolun/F/文件/个人/postdoc/presentation/materials/arrows/pts_normal_end.obj";
    std::string arrowfile = "/Users/wangb0d/bolun/F/文件/个人/postdoc/presentation/materials/arrows/arrowTri.obj";
    std::string pstartfile = "/Users/wangb0d/bolun/F/文件/个人/postdoc/presentation/materials/arrows/wrong_results/pts_start.obj";
    std::string pendfile = "/Users/wangb0d/bolun/F/文件/个人/postdoc/presentation/materials/arrows/wrong_results/pts_normal_end.obj";
    std::string exportfile = "/Users/wangb0d/bolun/F/文件/个人/postdoc/presentation/materials/arrows/normalUnionWrong.obj";

    ////////////////////////
    // // for binormal vectors
    // pendfile = "/Users/wangb0d/bolun/F/文件/个人/postdoc/presentation/materials/arrows/pts_binormal_end.obj";
    // exportfile = "/Users/wangb0d/bolun/F/文件/个人/postdoc/presentation/materials/arrows/binormalVectors.obj";
    // orient = true;
    // orientSign = -1;
    // ///////////////////////////////

    // ////////////////////////
    // // for messy vectors
    // pstartfile = "/Users/wangb0d/bolun/F/文件/个人/postdoc/presentation/materials/arrows/wrong_results/pts_start.obj";
    // pendfile = "/Users/wangb0d/bolun/F/文件/个人/postdoc/presentation/materials/arrows/wrong_results/pts_binormal_end.obj";
    // exportfile = "/Users/wangb0d/bolun/F/文件/个人/postdoc/presentation/materials/arrows/binormalVectorsWrong.obj";
    // orient = true;
    // orientSign = 1;
    ///////////////////////////////


    Eigen::MatrixXd Varrow, Vstart, Vend, V, VNormEnd;
    Eigen::MatrixXi Farrow, Fstart, Fend, F, Ftmp;
    igl::readOBJ(arrowfile, Varrow, Farrow);
    igl::readOBJ(pstartfile, Vstart, Fstart);
    igl::readOBJ(pendfile, Vend, Fend);
    std::cout << "readed " << std::endl;
    // rescaling and orienting the arrow
    Varrow *= scaling;
    if (orient)
    {
        igl::readOBJ(normEndfile, VNormEnd, Ftmp);
        Eigen::MatrixXd ref = VNormEnd - Vstart;
        orientEndpts(ref, Vstart, Vend, orientSign);
    }

    Eigen::Vector3d pstart = Vstart.row(0);
    Eigen::Vector3d pend = Vend.row(0);
    Eigen::MatrixXd Varrowtmp = Varrow;
    Eigen::MatrixXi Farrowtmp = Farrow;
    std::cout << "check  Varrowtmp " << Varrowtmp.rows() << std::endl;
    orientMeshToDirection(pstart, pend, Varrowtmp);
    std::cout << "check  Varrowtmp " << Varrowtmp.rows() << std::endl;
    std::cout << "check  Farrowtmp " << Farrowtmp.rows() << std::endl;
    V = Varrowtmp;
    F = Farrowtmp;
    int nbr = Vstart.rows();
    std::cout << "nbr of segs " << nbr << std::endl;

    for (int i = 1; i < nbr; i++)
    {
        pstart = Vstart.row(i);
        pend = Vend.row(i);
        Varrowtmp = Varrow;
        orientMeshToDirection(pstart, pend, Varrowtmp);
        extendMesh(V, F, Varrowtmp, Farrowtmp);
    }
    igl::writeOBJ(exportfile, V, F);
}

// quads
void generateFlist(const int nvstart, Eigen::MatrixXi& F){
    F.resize(nvstart - 1, 4);
    for (int i = 0; i < F.rows(); i++)
    {
        Eigen::Vector4i face;
        // i, i+1, i+nvstart+1, i+nvstart
        face << i, i + 1, i + nvstart + 1, i + nvstart;
        F.row(i) = face;
    }
}

// triangles
void generateFlistTri(const int nvstart, Eigen::MatrixXi& F){
    F.resize((nvstart - 1) * 2, 4);
    for (int i = 0; i < F.rows() / 2; i++)
    {
        Eigen::Vector3i face0, face1;
        // i, i+1, i+nvstart+1, i+nvstart
        face0 << i, i + 1, i + nvstart + 1;
        face1 << i, i + nvstart + 1, i + nvstart;
        F.row(2 * i) = face0;
        F.row(2 * i + 1) = face1;
    }
}

// quads
void generateFQuands(const Eigen::MatrixXd &Vstart, const Eigen::MatrixXd &Vend, Eigen::MatrixXd &Vall, Eigen::MatrixXi &F,
                     const double scaling)
{
    Eigen::MatrixXd D = Vend - Vstart;
    int nqd = D.rows() - 1;
    Eigen::MatrixXd Vstart0(nqd, 3), Vstart1(nqd, 3);
    Eigen::MatrixXd Vend0(nqd, 3), Vend1(nqd, 3);

    for (int i = 0; i < D.rows() - 1; i++)
    {
        Eigen::Vector3d direction = D.row(i) * scaling;
        Eigen::Vector3d tangent = (Vstart.row(i + 1) - Vstart.row(i)).normalized();
        Vstart0.row(i) = Eigen::Vector3d(Vstart.row(i)) + tangent * scaling / 2;
        Vstart1.row(i) = Eigen::Vector3d(Vstart.row(i)) - tangent * scaling / 2;
        Vend0.row(i) = direction + Eigen::Vector3d(Vstart0.row(i));
        Vend1.row(i) = direction + Eigen::Vector3d(Vstart1.row(i));
    }

    Vall.resize(nqd * 4, 3); 
    Vall.topRows(nqd) = Vstart0;
    Vall.middleRows(nqd, nqd) = Vstart1;
    Vall.middleRows(nqd * 2, nqd) = Vend0;
    Vall.bottomRows(nqd) = Vend1;
    
    F.resize(nqd, 4);

    for (int i = 0; i < F.rows(); i++)
    {
        Eigen::Vector4i face;
        // i, i+nqd, i+3*nqd, i+2*nqd
        face << i, i + nqd, i + 3 * nqd, i + 2 * nqd;
        F.row(i) = face;
    }
}

// scaling is to rescale the quad length
void writeRectifyingPlanes(double scaling)
{
    // for normal vectors
    bool orient = false;
    double orientSign = 1;
    std::string normEndfile = "/Users/wangb0d/bolun/F/文件/个人/postdoc/presentation/materials/arrows/pts_normal_end.obj";
    std::string arrowfile = "/Users/wangb0d/bolun/F/文件/个人/postdoc/presentation/materials/arrows/arrowTri.obj";
    std::string pstartfile = "/Users/wangb0d/bolun/F/文件/个人/postdoc/presentation/materials/arrows/pts_start.obj";
    std::string pendfile = "/Users/wangb0d/bolun/F/文件/个人/postdoc/presentation/materials/arrows/pts_normal_end.obj";
    std::string exportfile = "/Users/wangb0d/bolun/F/文件/个人/postdoc/presentation/materials/arrows/arrowUnion.obj";

    ////////////////////////
    // for binormal vectors
    pendfile = "/Users/wangb0d/bolun/F/文件/个人/postdoc/presentation/materials/arrows/pts_binormal_end.obj";
    exportfile = "/Users/wangb0d/bolun/F/文件/个人/postdoc/presentation/materials/arrows/quads.obj";
    orient = true;
    orientSign = 1; 
    ///////////////////////////////

    ////////////////////////
    // for messy vectors
    pstartfile = "/Users/wangb0d/bolun/F/文件/个人/postdoc/presentation/materials/arrows/wrong_results/pts_start.obj";
    pendfile = "/Users/wangb0d/bolun/F/文件/个人/postdoc/presentation/materials/arrows/wrong_results/pts_binormal_end.obj";
    exportfile = "/Users/wangb0d/bolun/F/文件/个人/postdoc/presentation/materials/arrows/quadsWrong.obj";
    orient = true;
    orientSign = 1;
    ///////////////////////////////


    Eigen::MatrixXd Varrow, Vstart, Vend, V, VNormEnd;
    Eigen::MatrixXi Farrow, Fstart, Fend, F, Ftmp;
    igl::readOBJ(arrowfile, Varrow, Farrow);
    igl::readOBJ(pstartfile, Vstart, Fstart);
    igl::readOBJ(pendfile, Vend, Fend);
    std::cout << "readed start size and end size "<<Vstart.rows()<<", "<<Vend.rows() << std::endl;
    // rescaling and orienting the arrow
    Varrow *= scaling;
    if (orient)
    {
        igl::readOBJ(normEndfile, VNormEnd, Ftmp);
        Eigen::MatrixXd ref = VNormEnd - Vstart;
        orientEndpts(ref, Vstart, Vend, orientSign);
    }
    int nvstart = Vstart.rows();
    // Eigen::MatrixXd Vall(2 * nvstart, 3);
    // Eigen::MatrixXi Fall;
    // Vall.topRows(nvstart) = Vstart;
    // Vend = Vstart + (Vend - Vstart) * scaling;
    // Vall.bottomRows(nvstart) = Vend;

    // generateFlist(nvstart, Fall);
    // generateFlistTri(nvstart, Fall);
    Eigen::MatrixXd Vout;
    Eigen::MatrixXi Fout;
    generateFQuands(Vstart, Vend, Vout, Fout, scaling);

    igl::writeOBJ(exportfile, Vout, Fout);
}

void writeOBJEdges(const std::string &fileout, const Eigen::MatrixXd &V) {
    std::ofstream fout;
    fout.open(fileout);
    int vnbr = V.rows();
    for (int i = 0; i < vnbr; i++)
    {
        fout << "v " << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << "\n";
    }
    for (int i = 0; i < vnbr - 1; i++)
    {
        fout << "l " << i + 1 << " " << i + 2 << "\n";
    }
    fout.close();
}

void writeMessyPoints(double scaling)
{
    // for normal vectors
    bool orient = false;
    double orientSign = 1;
    std::string normEndfile = "/Users/wangb0d/bolun/F/文件/个人/postdoc/presentation/materials/arrows/pts_normal_end.obj";
    std::string arrowfile = "/Users/wangb0d/bolun/F/文件/个人/postdoc/presentation/materials/arrows/arrowTri.obj";
    std::string pstartfile = "/Users/wangb0d/bolun/F/文件/个人/postdoc/presentation/materials/arrows/pts_start.obj";
    std::string pendfile = "/Users/wangb0d/bolun/F/文件/个人/postdoc/presentation/materials/arrows/pts_normal_end.obj";
    std::string exportfile = "/Users/wangb0d/bolun/F/文件/个人/postdoc/presentation/materials/arrows/arrowUnion.obj";

    ////////////////////////
    // for binormal vectors
    pendfile = "/Users/wangb0d/bolun/F/文件/个人/postdoc/presentation/materials/arrows/pts_binormal_end.obj";
    exportfile = "/Users/wangb0d/bolun/F/文件/个人/postdoc/presentation/materials/arrows/lines.obj";
    orient = true;
    orientSign = 1; 
    ///////////////////////////////

    ////////////////////////
    // for messy vectors
    // pstartfile = "/Users/wangb0d/bolun/F/文件/个人/postdoc/presentation/materials/arrows/wrong_results/pts_start.obj";
    // pendfile = "/Users/wangb0d/bolun/F/文件/个人/postdoc/presentation/materials/arrows/wrong_results/pts_binormal_end.obj";
    // exportfile = "/Users/wangb0d/bolun/F/文件/个人/postdoc/presentation/materials/arrows/linesWrong.obj";
    // orient = true;
    // orientSign = 1;
    ///////////////////////////////


    Eigen::MatrixXd Varrow, Vstart, Vend, V, VNormEnd;
    Eigen::MatrixXi Farrow, Fstart, Fend, F, Ftmp;
    igl::readOBJ(arrowfile, Varrow, Farrow);
    igl::readOBJ(pstartfile, Vstart, Fstart);
    igl::readOBJ(pendfile, Vend, Fend);
    // std::cout << "readed start size and end size "<<Vstart.rows()<<", "<<Vend.rows() << std::endl;
    // // rescaling and orienting the arrow
    // Varrow *= scaling;
    // if (orient)
    // {
    //     igl::readOBJ(normEndfile, VNormEnd, Ftmp);
    //     Eigen::MatrixXd ref = VNormEnd - Vstart;
    //     orientEndpts(ref, Vstart, Vend, orientSign);
    // }
    // int nvstart = Vstart.rows();
    // // Eigen::MatrixXd Vall(2 * nvstart, 3);
    // // Eigen::MatrixXi Fall;
    // // Vall.topRows(nvstart) = Vstart;
    // // Vend = Vstart + (Vend - Vstart) * scaling;
    // // Vall.bottomRows(nvstart) = Vend;

    // // generateFlist(nvstart, Fall);
    // // generateFlistTri(nvstart, Fall);
    // Eigen::MatrixXd Vout;
    // Eigen::MatrixXi Fout;
    // generateFQuands(Vstart, Vend, Vout, Fout, scaling);

    writeOBJEdges(exportfile, Vstart);
}

// whichHe is selected from 0, 1, 2
void getQuadFacesFromTriMesh(CGMesh &trimesh, const int whichHe, std::vector<std::array<int, 4>> &quads)
{
    quads.clear();
    for (CGMesh::VertexIter v_it = trimesh.vertices_begin(); v_it != (trimesh.vertices_end()); ++v_it)
    {
        CGMesh::VertexHandle vh = v_it.handle();
        if (trimesh.is_boundary(vh))
        {
            continue;
        }
        int counter = 0;
        CGMesh::HalfedgeHandle he0;
        for (CGMesh::VertexOHalfedgeIter voh = trimesh.voh_begin(vh); voh != trimesh.voh_end(vh); ++voh)
        {
            if (counter == whichHe)
            {
                he0 = voh.handle();
                break;
            }
            counter++;
        }
        CGMesh::HalfedgeHandle next = trimesh.next_halfedge_handle(he0);
        
        int v0 = trimesh.from_vertex_handle(he0).idx();
        int v1 = trimesh.to_vertex_handle(he0).idx();
        int v3 = trimesh.to_vertex_handle(next).idx();
        CGMesh::HalfedgeHandle oppo = trimesh.opposite_halfedge_handle(next);
        if (trimesh.is_boundary(oppo))
        {
            continue;
        }

        CGMesh::HalfedgeHandle opnx = trimesh.next_halfedge_handle(oppo);
        int v2 = trimesh.to_vertex_handle(opnx).idx();
        // verused[v0] = true;
        // verused[v1] = true;
        // verused[v2] = true;
        // verused[v3] = true;
        quads.push_back({v0, v1, v2, v3});
    }
}

void quadlist2Matrix(const std::vector<std::array<int, 4>> &quads, Eigen::MatrixXi &F)
{
    F.resize(quads.size(), 4);
    for (int i = 0; i < quads.size(); i++)
    {
        F(i, 0) = quads[i][0];
        F(i, 1) = quads[i][1];
        F(i, 2) = quads[i][2];
        F(i, 3) = quads[i][3];
    }
}

// write 2 of 3 families of polylines of the triangle mesh
void readTriMeshConvert2QuadMesh(){
    std::string fname = igl::file_dialog_open();
    if (fname.length() == 0)
    {
        std::cout << "\nLSC: read mesh failed" << std::endl;
        return;
    }
    CGMesh trimesh;
    OpenMesh::IO::read_mesh(trimesh, fname);
    std::cout << "\nMesh Readed" << std::endl;
    Eigen::MatrixXd V;
    Eigen::MatrixXi Ftri, Fqd0, Fqd1, Fqd2;
    std::vector<bool> verused(trimesh.n_vertices(), false);

    std::vector<std::array<int,4>> quads1, quads2, quads0;
    getQuadFacesFromTriMesh(trimesh, 0, quads0);
    getQuadFacesFromTriMesh(trimesh, 1, quads1);
    getQuadFacesFromTriMesh(trimesh, 2, quads2);
    MeshProcessing mp;
    mp.mesh2Matrix(trimesh, V, Ftri);
    quadlist2Matrix(quads0, Fqd0);
    quadlist2Matrix(quads1, Fqd1);
    quadlist2Matrix(quads2, Fqd2);
    std::cout<<"Writing the output quads, the prefix...\n";
    std::string fnameout = igl::file_dialog_save();
    igl::writeOBJ(fnameout + "f0.obj", V, Fqd0);
    igl::writeOBJ(fnameout + "f1.obj", V, Fqd1);
    igl::writeOBJ(fnameout + "f2.obj", V, Fqd2);
}

// whichDiag is chosen from 0 and 1
void QuadMeshToTri(CGMesh &quadmesh, Eigen::MatrixXi &F0, Eigen::MatrixXi& F1)
{
    std::vector<std::array<int, 3>> tris0, tris1;

    for (CGMesh::FaceIter f_it = quadmesh.faces_begin(); f_it != (quadmesh.faces_end()); ++f_it)
    {
        int fid = f_it.handle().idx();
        int ecounter = 0;
        std::array<int, 4> vlocal;
        for (CGMesh::FaceVertexIter fv_it = quadmesh.fv_begin(f_it); fv_it != (quadmesh.fv_end(f_it)); ++fv_it)
        {
            vlocal[ecounter] = fv_it.handle().idx();
            ecounter++;
        }
        tris0.push_back({vlocal[0], vlocal[1], vlocal[2]});
        tris0.push_back({vlocal[0], vlocal[2], vlocal[3]});

        tris1.push_back({vlocal[0], vlocal[1], vlocal[3]});
        tris1.push_back({vlocal[1], vlocal[2], vlocal[3]});
    }
    F0.resize(tris0.size(), 3);
    F1.resize(tris0.size(), 3);
    for (int i = 0; i < F0.rows(); i++)
    {
        for (int j = 0; j < 3; j++)
        {
            F0(i, j) = tris0[i][j];
            F1(i, j) = tris1[i][j];
        }
    }
}

void readQuadMesh2TriMesh()
{
    std::string fname = igl::file_dialog_open();
    if (fname.length() == 0)
    {
        std::cout << "\nLSC: read mesh failed" << std::endl;
        return;
    }
    CGMesh quadmesh;
    OpenMesh::IO::read_mesh(quadmesh, fname);
    std::cout << "\nMesh Readed" << std::endl;
    Eigen::MatrixXd V;
    Eigen::MatrixXi Fqd, Ftr0, Ftr1;

    std::vector<std::array<int, 4>> quads1, quads2, quads0;
    QuadMeshToTri(quadmesh, Ftr0, Ftr1);
    MeshProcessing mp;
    mp.mesh2Matrix(quadmesh, V, Fqd);
    std::cout << "Writing the output quads, the prefix...\n";
    std::string fnameout = igl::file_dialog_save();
    igl::writeOBJ(fnameout + "t0.obj", V, Ftr0);
    igl::writeOBJ(fnameout + "t1.obj", V, Ftr1);
}
// this is to deal with a regular patch
void readQuadMesh2TriMesh(const int vinrow)
{
    std::string fname = igl::file_dialog_open();
    if (fname.length() == 0)
    {
        std::cout << "\nLSC: read mesh failed" << std::endl;
        return;
    }
    CGMesh quadmesh;
    OpenMesh::IO::read_mesh(quadmesh, fname);
    std::cout << "\nMesh Readed" << std::endl;
    Eigen::MatrixXd V;
    Eigen::MatrixXi Fqd, Ftr0, Ftr1;
    MeshProcessing mp;
    mp.mesh2Matrix(quadmesh, V, Fqd);
    Ftr0.resize(Fqd.rows() * 2, 3);
    Ftr1 = Ftr0;
    
    int cnbr = V.rows() / vinrow;
    std::cout<<"there are "<<vinrow<<" vertices in each row, nbr of vers "<<V.rows()<<"\n";
    std::cout<<"there are "<<cnbr<<"columns, nbr of triangles "<<Ftr0.rows()<<"\n";
    int counter = 0;
    for (int i = 0; i < vinrow - 1; i++)
    {
        for (int j = 0; j < cnbr - 1; j++)
        {
            int id0 = i + j * vinrow;
            int id1 = id0 + 1;
            int id2 = i + 1 + (j + 1) * vinrow;
            int id3 = i + (j + 1) * vinrow;
            Ftr0.row(counter) << id0, id1, id2;
            Ftr0.row(counter + 1) << id0, id2, id3;
            Ftr1.row(counter) << id0, id1, id3;
            Ftr1.row(counter + 1) << id1, id2, id3;
            counter +=2;
        }
    }

    std::cout << "Writing the output quads, the prefix...\n";
    std::string fnameout = igl::file_dialog_save();
    igl::writeOBJ(fnameout + "t0.obj", V, Ftr0);
    igl::writeOBJ(fnameout + "t1.obj", V, Ftr1);
}



void evaluateGGGConsineConstraints()
{
    std::string fname = igl::file_dialog_open();
    if (fname.length() == 0)
    {
        std::cout << "\nLSC: read mesh failed" << std::endl;
        return;
    }
    CGMesh trimesh;
    OpenMesh::IO::read_mesh(trimesh, fname);
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    MeshProcessing mp;
    mp.mesh2Matrix(trimesh, V, F);
    std::cout << "\nMesh Readed" << std::endl;
    int vnbr = V.rows();
    Eigen::VectorXd energy = Eigen::VectorXd::Zero(V.rows() * 3);
    for (CGMesh::VertexIter v_it = trimesh.vertices_begin(); v_it != (trimesh.vertices_end()); ++v_it)
    {
        CGMesh::VertexHandle vh = v_it.handle();
        if (trimesh.is_boundary(vh))
        {
            continue;
        }
        if (trimesh.valence(vh) != 6)
        {
            continue;
        }
        int vid = vh.idx();
        Eigen::Vector3d v = V.row(vh.idx());
        std::vector<Eigen::Vector3d> vs;
        for (CGMesh::VertexOHalfedgeIter voh = trimesh.voh_begin(vh); voh != trimesh.voh_end(vh); ++voh)
        {
            CGMesh::HalfedgeHandle he0 = voh.handle();
            vs.push_back(V.row(trimesh.to_vertex_handle(he0).idx()));
        }
        Eigen::Vector3d vv0 = (vs[0] - v).normalized();
        Eigen::Vector3d vv1 = (vs[1] - v).normalized();
        Eigen::Vector3d vv2 = (vs[2] - v).normalized();
        Eigen::Vector3d vv3 = (vs[3] - v).normalized();
        Eigen::Vector3d vv4 = (vs[4] - v).normalized();
        Eigen::Vector3d vv5 = (vs[5] - v).normalized();

        double e0 = vv0.dot(vv1) - vv3.dot(vv4);
        double e1 = vv1.dot(vv2) - vv4.dot(vv5);
        double e2 = vv2.dot(vv3) - vv5.dot(vv0);
        energy[vid] = e0;
        energy[vid + vnbr] = e1;
        energy[vid + vnbr * 2] = e2;
    }
    std::cout << "the total energy is " << energy.norm() << "\n";
}

// project the point onto a curve
bool projectPointOnCurve(const std::vector<Eigen::Vector3d> &curve, const Eigen::Vector3d &pt,
                         int &segid, double &tlocal, Eigen::Vector3d &plocal, Eigen::Vector3d& tangent)
{
    segid = -1;
    double dis_max = std::numeric_limits<double>::max();
    for (int i = 0; i < curve.size() - 1; i++)
    {
        Eigen::Vector3d ps = curve[i], pe = curve[i + 1];
        double t = (pt - ps).dot(pe - ps) / (pe - ps).dot(pe - ps);
        if (t < 1e-6 || t > 1 + 1e-6)
        {
            continue;
        }
        Eigen::Vector3d tangent_left, tangent_right, tangent_this;
        tangent_this = pe - ps;
        if (i > 0)
        {
            tangent_left = curve[i] - curve[i - 1];
        }
        else{
            tangent_left = tangent_this;
        }
        if (i < curve.size() - 2)
        {
            tangent_right = curve[i + 2] - curve[i + 1];
        }
        else{
            tangent_right = tangent_this;
        }
        if (t < 0.333)
        {
            tangent_this = tangent_this + tangent_left;
        }
        if (t > 0.667)
        {
            tangent_this = tangent_this + tangent_right;
        }
        double distance = (ps + (pe - ps) * t - pt).norm(); // the distance between the projection and the point
        if (dis_max > distance)
        {
            dis_max = distance;
            segid = i;
            tlocal = t;
            plocal = ps + (pe - ps) * t;
            tangent = tangent_this.normalized();
        }
    }
    if (segid >= 0)
    {
        return true;
    }

    return false;
}
void upsample_and_smooth_curve()
{
    double ratio = 0.3;
    std::string fname = igl::file_dialog_open();
    if(fname.length()==0){
        std::cout << "please read curve\n";
        return;
    }
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::readOBJ(fname, V, F);
    std::vector<Eigen::Vector3d> vout(V.rows() * 2 - 1);
    for (int i = 0; i < V.rows(); i++)
    {
        vout[i * 2] = V.row(i);
    }
    // now vout[0], vout[2], ... are the original data
    for (int itr = 0; itr < 10; itr++)
    {
        double error = 0;
        for (int i = 0; i < V.rows() - 1; i++)
        {
            int vid = 2 * i + 1; // the vid we are looking at
            if (itr == 0)
            { // init
                vout[vid] = (vout[vid - 1] + vout[vid + 1]) / 2;
                continue;
            }
            Eigen::Vector3d delta = (vout[vid - 1] + vout[vid + 1]) / 2 - vout[vid];
            error += delta.dot(delta);
            vout[vid] += delta * ratio;
            if (i > 0)
            {
                delta = 2 * vout[vid - 1] - vout[vid - 2] - vout[vid];
                vout[vid] += delta * ratio;
                error += delta.dot(delta);
            }
            if (i < V.rows() - 2)
            {
                delta = 2 * vout[vid + 1] - vout[vid + 2] - vout[vid];
                vout[vid] += delta * ratio;
                error += delta.dot(delta);
            }
        }
        std::cout<<"current smoothness error, "<<error<<"\n";
    }
    std::cout << "saving files\n";
    fname = igl::file_dialog_save();
    if(fname.length()==0){
        std::cout << "please write down curve name\n";
        return;
    }
    igl::writeOBJ(fname, vec_list_to_matrix(vout), F);
}


void smooth_curve()
{
    double ratio = 1;
    std::string fname = igl::file_dialog_open();
    if(fname.length()==0){
        std::cout << "please read curve\n";
        return;
    }
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    igl::readOBJ(fname, V, F);
    std::vector<Eigen::Vector3d> vout = mat_to_vec_list(V);
    // now vout[0], vout[2], ... are the original data
    for (int itr = 0; itr < 10; itr++)
    {
        double error = 0;
        for (int i = 1; i < V.rows() - 1; i++)
        {
            int vid = i; // the vid we are looking at
            if (itr == 0)
            { // init
                vout[vid] = (vout[vid - 1] + vout[vid + 1]) / 2;
                continue;
            }
            Eigen::Vector3d delta = (vout[vid - 1] + vout[vid + 1]) / 2 - vout[vid];
            error += delta.dot(delta);
            vout[vid] += delta * ratio;
        }
        std::cout<<"current smoothness error, "<<error<<"\n";
    }
    std::cout << "saving files\n";
    fname = igl::file_dialog_save();
    if(fname.length()==0){
        std::cout << "please write down curve name\n";
        return;
    }
    igl::writeOBJ(fname, vec_list_to_matrix(vout), F);
}

void AggFirstStripWithGuideCurve(const std::vector<Eigen::Vector3d> &plyin, const std::vector<Eigen::Vector3d> &binin,
                                QuadOpt &quad_tool, const Eigen::MatrixXd &RefCurve, std::vector<Eigen::Vector3d> &versOut)
{
    double distance = (plyin[0] - plyin[1]).norm();
    // find the projection on the curve, and use this projected point to find the next one.
    // march forward from plocal0 to plocal1
    int segid0, segid1;
    double tlocal0, tlocal1;
    Eigen::Vector3d plocal0, tangent0, plocal1;
    if (quad_tool.V.rows() == 0)
    {
        segid0 = 0;
        tlocal0 = 0;
        plocal0 = RefCurve.row(0);
    }
    else
        projectPointOnCurve(mat_to_vec_list(RefCurve), plyin[0], segid0, tlocal0, plocal0, tangent0);
    std::cout << "start from point segid, " << segid0 << ", t, " << tlocal0 << ", point, " << plocal0.transpose() << "\n";
    std::cout << "march distance, " << distance << "\n";
    find_next_pt_on_polyline(segid0, mat_to_vec_list(RefCurve), distance,
                             plocal0, segid1, plocal1, tlocal1);
    std::cout << "find point segid, " << segid1 << ", t, " << tlocal1 << ", point, " << plocal1.transpose() << "\n";
    aggFirstStripGuideGeodesic(plyin, binin, versOut, mat_to_vec_list(RefCurve), segid1, tlocal1, false);
}