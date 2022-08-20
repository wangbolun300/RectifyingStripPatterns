#include<lsc/basic.h>

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
void lsTools::debug_tool(int id, double value)
    {
        int vnbr=V.rows();
        refids.clear();
        for(int i=0;i<vnbr;i++){
            if(ORB.coeffRef(i,i)>0){
                refids.push_back(i);
            }
        }
    }
void extend_triplets_offset(std::vector<Trip>& triplets, const spMat& mat, int offrow, int offcol){
    for (int i = 0; i < mat.outerSize(); i++)
    {
        for (spMat::InnerIterator it(mat, i); it; ++it)
        {
            triplets.push_back(Trip(it.row()+offrow, it.col()+offcol, it.value()));
        }
    }
}
spMat dense_mat_list_as_sparse_diagnal(const std::vector<Eigen::MatrixXd>& mlist){
    spMat result;
    int outer=mlist.size();
    int inner_rows=mlist[0].rows();
    int inner_cols=mlist[0].cols();
    result.resize(outer*inner_rows,outer*inner_cols);
    std::vector<Trip> triplets;
    triplets.reserve(inner_cols*inner_rows*outer);
    for(int i=0;i<outer;i++){
        extend_triplets_offset(triplets,mlist[i].sparseView(),i*inner_rows,i*inner_cols);
    }
    result.setFromTriplets(triplets.begin(),triplets.end());
    return result;
}
spMat rib_method_arrange_matrices_rows(const std::vector<spMat>& mats){
    int nvec=mats.size();
    int inner_rows=mats[0].rows();
    int inner_cols=mats[0].cols();
    spMat result;
    result.resize(nvec*inner_rows,inner_cols);// the rows get more, the column number remain the same.
    std::vector<Trip> triplets;
    triplets.reserve(inner_cols*inner_rows*nvec);
    for(int i=0;i<nvec;i++){
        //TODO
    }
}