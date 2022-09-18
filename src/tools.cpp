#include <lsc/basic.h>
#include <lsc/tools.h>

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
    assert(!mat.IsRowMajor);
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
    ver_dbg.resize(0,0);
    ver_dbg1.resize(0,0);
    flag_dbg = false;
    E0_dbg.resize(0,0);
    direction_dbg.resize(0,0);
    pnorm_dbg.resize(0,0);
    pnorm_list_dbg.clear();
    flag_dbg = true;
    id_dbg = id2;

    // debug init end
    double start_point_para = 0.5;
    double start_angle_degree = 60;
    double target_angle = value;
    std::vector<Eigen::Vector3d> curve;
    
    trace_single_pseudo_geodesic_curve(target_angle, Boundary_Edges[id], start_point_para, start_angle_degree,
                                       curve);
    flag_dbg = false;
    int nbr_midpts = curve.size() - 2;
    E0_dbg.resize(nbr_midpts, 3);
    direction_dbg.resize(nbr_midpts, 3);
    for (int i = 0; i < nbr_midpts; i++)
    {
        E0_dbg.row(i) = curve[i + 1];
        Eigen::Vector3d vec0 = (pseudo_vers_dbg[i + 1] - pseudo_vers_dbg[i]).normalized();
        Eigen::Vector3d vec1 = (curve[i + 2] - pseudo_vers_dbg[i + 1]).normalized();
        direction_dbg.row(i) = vec0.cross(vec1).normalized();// pseudo normal
    }
    if(pnorm_list_dbg.size() != nbr_midpts){
        std::cout<<"dangerous size! "<<pnorm_list_dbg.size()<<", "<<nbr_midpts<<std::endl;
    }
    if (target_angle > 90 - ANGLE_TOLERANCE && target_angle < 90 + ANGLE_TOLERANCE){
        std::cout<<"Traced geodesic"<<std::endl;
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
        }
    }

    // TODO temporarily checking one trace line
    trace_vers[0] = curve;
}
void lsTools::debug_tool_v2(const std::vector<int>& ids, const std::vector<double> values){
    // cylinder_example(5, 10, 50, 30);
    // exit(0);
    int nbr_itv= ids[0]; // every nbr_itv boundary edges we shoot one curve
    if(nbr_itv<1){
        std::cout<<"Please set up the parameter nbr_itv "<<std::endl;
        return;
    }
    double target_angle=values[0];
    double start_angel=values[1];
    OpenMesh::HalfedgeHandle init_edge=Boundary_Edges[0];
    if(lsmesh.face_handle(init_edge).idx()>=0){
        init_edge=lsmesh.opposite_halfedge_handle(init_edge);
    }
    OpenMesh::HalfedgeHandle checking_edge=init_edge;
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
        trace_single_pseudo_geodesic_curve(target_angle, checking_edge, 0.5, start_angel,
                                           curve);
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
        bool stop_flag=false;
        
        if(curve_id==-1){
            break;
        }
        for (int i = 0; i < nbr_itv; i++)
        {
            checking_edge = lsmesh.next_halfedge_handle(checking_edge);
            if (checking_edge == init_edge)
            {
                stop_flag=true;
                break;
            }
        }
        if(stop_flag){
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