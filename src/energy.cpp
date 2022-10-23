#include <lsc/basic.h>
#include <lsc/tools.h>

// directions[0] is the inward pointing to the point, and direction[1] is outward shooting from the point. 
// handles are the two halfedges opposite to the point. 
// l2s shows if the value is from large to small along the halfedge direction
bool find_active_faces_and_directions_around_ver(CGMesh &lsmesh, const Eigen::VectorXd &fvalues, 
                                        const Eigen::MatrixXd& V,
                                        const Eigen::Vector3d& normal, const int vid,
                                       std::array<Eigen::Vector3d,2> &directions,
                                       std::array<bool,2> &l2s,// if the halfedge handle is from large to small
                                      std::array<CGMesh::HalfedgeHandle,2>& handles)
{
    CGMesh::VertexHandle vh=lsmesh.vertex_handle(vid);
    std::vector<Eigen::Vector3d> dircs;
    std::vector<int> order;// 0 means the first segment, 1 means the second segment
    std::vector<CGMesh::HalfedgeHandle> hes;
    std::vector<bool> large_to_small;
    double value=fvalues[vid];
    for(CGMesh::VertexOHalfedgeIter voh=lsmesh.voh_begin(vh);voh!=lsmesh.voh_end(vh);++voh){
        CGMesh::HalfedgeHandle next_he=lsmesh.next_halfedge_handle(voh);
        int id_from=lsmesh.from_vertex_handle(next_he).idx();
        int id_to=lsmesh.to_vertex_handle(next_he).idx();
        double value1=fvalues[id_from];
        double value2=fvalues[id_to];
        if(value1==value2) {
            continue;
        }
        double value_min;
        double value_max;
        int idmin,idmax;
        bool tmp_l2s;
        if(value1<value2){
            value_min=value1;
            value_max=value2;
            idmin=id_from;
            idmax=id_to;
            tmp_l2s=false;
        }
        else{
            value_min=value2;
            value_max=value1;
            idmin=id_to;
            idmax=id_from;
            tmp_l2s=true;
        }
        if(value<value_min||value>=value_max){
            continue;
        }
        double t;
        t=  get_t_of_value(value, value1, value2);
        Eigen::Vector3d vedge=get_3d_ver_from_t(t,V.row(id_from),V.row(id_to));
        // Eigen::Vector3d big2small=V.row(idmin)-V.row(idmax);
        Eigen::Vector3d edge2ver=Eigen::Vector3d(V.row(vid))-vedge;
        // double flag=big2small.cross(edge2ver).dot(normal);
        if(tmp_l2s){// the halfedge is from large to small, it is the shoot-in direction
            dircs.push_back(edge2ver);
            order.push_back(0);
        }
        else{// it is the shoot-out direction
            dircs.push_back(-edge2ver);
            order.push_back(1);
        }
        int fid=lsmesh.face_handle(next_he).idx();
        hes.push_back(next_he);
        large_to_small.push_back(tmp_l2s);

    }
    if(dircs.size()!=2){// it is a singularity
        return false;
    }
    assert(order[0]+order[1]==1);
    for(int i=0;i<2;i++){
        handles[order[i]]=hes[i];
        directions[order[i]]=dircs[i];
        l2s[order[i]]=large_to_small[i];
    }
    assert(lsmesh.face_handle(handles[0]).idx()==lsmesh.face_handle(handles[1]).idx());
    return true;
   
}
// TODO this can be pre-computed before iterations: for each triangle we compute 3 vectors
Eigen::Vector3d get_coff_vec_for_gradient( std::array<spMat, 3> &gradVF, const int fid, const int vid)
{
    Eigen::Vector3d result;
    result[0]=gradVF[0].coeffRef(fid,vid);
    result[1]=gradVF[1].coeffRef(fid,vid);
    result[2]=gradVF[2].coeffRef(fid,vid);
    return result;
}
// v1, v2, vM is one triangle. v3, v4, vM is one triangle. 
// This code works only for v1, v2, v3, v4 are not duplicated
std::vector<Trip> get_triplets_for_gradient_partial_parts_vertex_based( std::array<spMat, 3>& gradVF, const int v1, const int v2, const int v3, const int v4, const int vM,
 const int fid1, const int fid2, const Eigen::Vector3d& norm, const Eigen::Vector3d& norm1, 
const Eigen::Vector3d& norm2){
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
    triplets.push_back(Trip(vM, vM, 2*cfmfm));
    return triplets;
}

std::vector<Trip> get_triplets_for_gradient_partial_parts_edge_based( std::array<spMat, 3> &gradVF,
                                                                     const int v1, const int v2, const int vA, const int vB,
                                                                     const int fid1, const int fid2, const Eigen::Vector3d &norm, const Eigen::Vector3d &norm1,
                                                                     const Eigen::Vector3d &norm2)
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
void lsTools::calculate_gradient_partial_parts(){
    
    int act_size=Actid.size();
    PEJC.resize(act_size);
    for(int i=0;i<act_size;i++){
        PEJC[i].resize(V.rows(),V.rows());
    }
    // std::array<Eigen::MatrixXd,8> 
    for(int i=0;i<act_size;i++){
        int fid1, fid2;
        std::array<int,4> vids;
        // get the v1, v2, vA, vB.
        vids=get_vers_around_edge(lsmesh,Actid[i], fid1,fid2);
        int v1=vids[0];
        int v2=vids[1];
        int vA=vids[2];
        int vB=vids[3];
       
        Eigen::Vector3d norm=norm_e.row(Actid[i]);
        Eigen::Vector3d norm1=norm_f.row(fid1);
        Eigen::Vector3d norm2=norm_f.row(fid2);
        std::vector<Trip> triplets = get_triplets_for_gradient_partial_parts_edge_based(gradVF,v1,v2,vA,vB,
        fid1,fid2,norm,norm1,norm2);
       
        PEJC[i].setFromTriplets(triplets.begin(),triplets.end());
    }
}
// calculate vector of matrix A, A[i]*f is the jacobian of (g1xg2).dot(norm).
// here g1, g2 are the rotated gradients
void lsTools::calculate_gradient_partial_parts_ver_based(){
    int ninner=IVids.size();
    bool first_time_compute = false; // if it is first time, we need to compute all; otherwise we just need to update
    if (VPEJC.empty())
    {
        VPEJC.resize(ninner);
        Vheh0.resize(ninner);
        Vheh1.resize(ninner);
        ActInner.resize(ninner);
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
        int fid1=lsmesh.face_handle(hd1).idx();
        int fid2=lsmesh.face_handle(hd2).idx();
        assert(fid1!=fid2);
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
            VPEJC[i].setFromTriplets(triplets.begin(), triplets.end());
            Vheh0[i] = hd1;
            Vheh1[i] = hd2;
        }
    }
}

void lsTools::calculate_pseudo_energy_function_values(const double angle_degree){

    int act_size=Actid.size();
    EPEvalue.resize(act_size);
    LsOrient.resize(act_size);// first make g1 g2 same direction, then check if g1xg2 or g2xg1.
    PeWeight.resize(act_size);
    double angle_radian = angle_degree * LSC_PI / 180.; // the angle in radian
    double cos_angle=cos(angle_radian);
    // std::array<Eigen::MatrixXd,8>
    for (int i = 0; i < act_size; i++)
    {
        int fid1, fid2;
        std::array<int, 4> vids;
        // get the v1, v2, vA, vB.
        vids = get_vers_around_edge(lsmesh, Actid[i], fid1, fid2);
        int v1=vids[0];
        int v2=vids[1];
        int vA=vids[2];
        int vB=vids[3];
        Eigen::Vector3d norm = norm_e.row(Actid[i]);
        Eigen::Vector3d norm1=norm_f.row(fid1);
        Eigen::Vector3d norm2=norm_f.row(fid2);
        Eigen::Vector3d g1=gfvalue.row(fid1);
        Eigen::Vector3d g2=gfvalue.row(fid2);
        // rotate gradients to get iso-curve directions
        g1=norm1.cross(g1);
        g2=norm2.cross(g2);
        Eigen::Vector3d g1xg2=g1.cross(g2);

        
        // deal with edge orientation: shoot from small value to big value
        Eigen::Vector3d edirec = V.row(v1) - V.row(v2);
        double flag1=g1.cross(edirec).dot(norm);
        double flag2=g2.cross(edirec).dot(norm);
        double flag=flag1*flag2;
        if(flag>=0){
            LsOrient[i]=1;
        }
        else{
            LsOrient[i]=-1;
        }
        // check if left triangle is on the left
        Eigen::Vector3d ldirec=V.row(vA)-V.row(v1);
        assert(edirec.cross(ldirec).dot(norm)>=0);

        // deal with PeWeight
        double lg1=g1.norm();
        double lg2=g2.norm();
        if(lg1<1e-16 || lg2<1e-16){
            PeWeight[i]=0;// it means the g1xg2 will not be accurate
        }
        else{
            Eigen::Vector3d g1n = g1.normalized();
            Eigen::Vector3d g2n = LsOrient[i]*g2.normalized();// get the oriented g2, to make g1 g2 goes to the same direction
            double cos_real=g1n.dot(g2n);
            PeWeight[i]=(cos_real+1)/2;// cos = 1, weight is 1; cos = -1, means it is very sharp turn, weight = 0
        }

        if(fvalues[v1]<fvalues[v2]){ // orient g1xg2 or g2xg1
            LsOrient[i]*=-1;
        }
        double value = LsOrient[i] * g1xg2.dot(norm) - g1xg2.norm() * cos_angle;
        EPEvalue.coeffRef(i)=value;
    }
}
void lsTools::calculate_pseudo_energy_function_values_vertex_based(const double angle_degree){
    int ninner=IVids.size();
    VPEvalue.resize(ninner);
    LsOrient.resize(ninner);// first make g1 g2 same direction, then check if g1xg2 or g2xg1.
    PeWeight.resize(ninner);
    double angle_radian = angle_degree * LSC_PI / 180.; // the angle in radian
    double cos_angle=cos(angle_radian);
    TODO
    for (int i = 0; i < act_size; i++)
    {
        int fid1, fid2;
        std::array<int, 4> vids;
        // get the v1, v2, vA, vB.
        vids = get_vers_around_edge(lsmesh, Actid[i], fid1, fid2);
        int v1=vids[0];
        int v2=vids[1];
        int vA=vids[2];
        int vB=vids[3];
        Eigen::Vector3d norm = norm_e.row(Actid[i]);
        Eigen::Vector3d norm1=norm_f.row(fid1);
        Eigen::Vector3d norm2=norm_f.row(fid2);
        Eigen::Vector3d g1=gfvalue.row(fid1);
        Eigen::Vector3d g2=gfvalue.row(fid2);
        // rotate gradients to get iso-curve directions
        g1=norm1.cross(g1);
        g2=norm2.cross(g2);
        Eigen::Vector3d g1xg2=g1.cross(g2);

        
        // deal with edge orientation: shoot from small value to big value
        Eigen::Vector3d edirec = V.row(v1) - V.row(v2);
        double flag1=g1.cross(edirec).dot(norm);
        double flag2=g2.cross(edirec).dot(norm);
        double flag=flag1*flag2;
        if(flag>=0){
            LsOrient[i]=1;
        }
        else{
            LsOrient[i]=-1;
        }
        // check if left triangle is on the left
        Eigen::Vector3d ldirec=V.row(vA)-V.row(v1);
        assert(edirec.cross(ldirec).dot(norm)>=0);

        // deal with PeWeight
        double lg1=g1.norm();
        double lg2=g2.norm();
        if(lg1<1e-16 || lg2<1e-16){
            PeWeight[i]=0;// it means the g1xg2 will not be accurate
        }
        else{
            Eigen::Vector3d g1n = g1.normalized();
            Eigen::Vector3d g2n = LsOrient[i]*g2.normalized();// get the oriented g2, to make g1 g2 goes to the same direction
            double cos_real=g1n.dot(g2n);
            PeWeight[i]=(cos_real+1)/2;// cos = 1, weight is 1; cos = -1, means it is very sharp turn, weight = 0
        }

        if(fvalues[v1]<fvalues[v2]){ // orient g1xg2 or g2xg1
            LsOrient[i]*=-1;
        }
        double value = LsOrient[i] * g1xg2.dot(norm) - g1xg2.norm() * cos_angle;
        EPEvalue.coeffRef(i)=value;
    }
}
void lsTools::assemble_solver_boundary_condition_part(spMat& H, Efunc& B){
    spMat bcJacobian;
    Efunc bcVector;
    assert(trace_vers.size()==trace_hehs.size());
    int size=0;// the number of constraints
    std::vector<double> vec_elements;
    std::vector<Trip> triplets;
    vec_elements.reserve(lsmesh.n_edges());
    triplets.reserve(lsmesh.n_edges());
    assert(assigned_trace_ls.size()==trace_vers.size());
    // assigned_trace_ls.resize(trace_vers.size());
    // for (int i = 0; i < trace_vers.size(); i++)// TODO assign assigned_trace_ls when tracing
    // {
    //     assigned_trace_ls[i]=i;
    // }
    // std::cout<<"before going into for loop"<<std::endl;
    for (int i = 0; i < trace_vers.size(); i++)
    {
        
        assert(trace_vers[i].size()==trace_hehs[i].size());
        for(int j=0;j<trace_vers[i].size();j++){
            
            CGMesh::HalfedgeHandle edge=trace_hehs[i][j];
            Eigen::Vector3d point_middle=trace_vers[i][j];
            int id_from=lsmesh.from_vertex_handle(edge).idx();
            int id_to=lsmesh.to_vertex_handle(edge).idx();
            Eigen::Vector3d point_from=V.row(id_from);
            Eigen::Vector3d point_to=V.row(id_to);
            double d1=(point_middle-point_from).norm();
            double d2=(point_to-point_from).norm();
            triplets.push_back(Trip(size, id_from, d2 - d1));
            triplets.push_back(Trip(size, id_to, d1));
            double right_value=(d2-d1)*fvalues[id_from]+d1*fvalues[id_to]-d2*assigned_trace_ls[i];
            vec_elements.push_back(right_value);
            size++;
        }
    }
    // std::cout<<"after loop"<<std::endl;
    // get boundary condition: the jacobian matrix
    bcJacobian.resize(size,V.rows());
    bcJacobian.setFromTriplets(triplets.begin(),triplets.end());// TODO calculate this before optimization
    bcVector.resize(size);
    // std::cout<<"calculated a and b"<<std::endl;
    for(int i=0;i<size;i++){// get the f(x)
        double value=vec_elements[i];
        bcVector.coeffRef(i) = value;
    }
    bcfvalue = bcVector;
    // get JTJ
    H = bcJacobian.transpose() * bcJacobian;
    // get the -(J^T)*f(x)
    B = -bcJacobian.transpose() * bcVector;
}
void lsTools::assemble_solver_pesudo_geodesic_energy_part(spMat &H, Efunc &B)
{
    calculate_pseudo_energy_function_values(pseudo_geodesic_target_angle_degree);
    // std::cout<<"function values get calculated"<<std::endl;

    int act_size = Actid.size();
    int vnbr = V.rows();
    // std::cout<<"vnbr "<<vnbr<<std::endl;
    spMat JTJ;
    spMat JTf;
    JTJ.resize(vnbr, vnbr);
    JTf.resize(vnbr,1);
    Eigen::MatrixXd Mfvalues(vnbr, 1);
    Mfvalues.col(0) = fvalues;
    spMat SMfvalues = Mfvalues.sparseView();
    
    for (int i = 0; i < act_size; i++)
    {
        // std::cout<<"SMfvalues size "<<SMfvalues.rows()<<" "<<SMfvalues.cols()<<std::endl;
        spMat JT = LsOrient[i] * PEJC[i] * SMfvalues;
        // std::cout<<"J size "<<J.rows()<<" "<<J.cols()<<std::endl;
        spMat tjtj=JT* JT.transpose();
        // std::cout<<"tempJTJ size "<<tjtj.rows()<<" "<<tjtj.cols()<<std::endl;
        JTJ += PeWeight[i] * tjtj;
        // std::cout<<"JTJ size "<<JTJ.rows()<<" "<<JTJ.cols()<<std::endl;
        spMat temp=PeWeight[i] * JT * EPEvalue.coeffRef(i);
        // std::cout<<"tmp size "<<temp.rows()<<" "<<temp.cols()<<std::endl;
        JTf +=temp;
    }
    // std::cout<<"assembled"<<std::endl;
    H = JTJ;
    B = sparse_mat_col_to_sparse_vec(-JTf,0);
}

// check if active faces exist around one vertex
// the information can be used to construct pseudo-geodesic energy


// before calling this function, please get the traced curves and assign values 
// for these curves
void lsTools::optimize_laplacian_with_traced_boundary_condition(){
    int vnbr = V.rows();
    int fnbr = F.rows();
    assert(assigned_trace_ls.size()&& "Please trace the curves first before solve the energy");
    assert(assigned_trace_ls.size()==trace_vers.size());
    // initialize the level set with some number
    if (fvalues.size() != vnbr)
    {
        initialize_level_set_accroding_to_parametrization();
        std::cout<<"level set get initialized for smoothing"<<std::endl;
        return;
    }
    // after initialization we can calculate some quantities associated with level set values
    get_gradient_hessian_values();
    // the matrices
    spMat H;
    H.resize(vnbr,vnbr);
    Efunc B;
    B.resize(vnbr);
    spMat bc_JTJ;
    Efunc bc_mJTF;
    // get the boundary condition bcJacobian (J) and bcVector F, 
    // std::cout<<"before assembling matrices"<<std::endl;
    assemble_solver_boundary_condition_part(bc_JTJ, bc_mJTF);
    // std::cout<<"boundary condition set up"<<std::endl;
    spMat LTL;// left of laplacian
    Efunc mLTF;// right of laplacian
    // assemble_solver_laplacian_part(LTL, mLTF);
    assemble_solver_biharmonic_smoothing(LTL,mLTF);
    // std::cout<<"laplacian condition set up"<<std::endl;
    assert(mass.rows() == vnbr);
    
    // gravity + laplacian + boundary condition
    H = weight_mass * mass + weight_laplacian * LTL + weight_boundary * bc_JTJ;
    B=weight_laplacian*mLTF+weight_boundary*bc_mJTF;
    if(enable_pseudo_geodesic_energy&&weight_pseudo_geodesic_energy==0){
        std::cout<<"Pseudo-geodesic Energy weight is 0"<<std::endl;
    }
    if (enable_pseudo_geodesic_energy && weight_pseudo_geodesic_energy > 0)
    {
        spMat pg_JTJ;
        Efunc pg_mJTF;
        // std::cout<<"before solving pg energy"<<std::endl;
        assemble_solver_pesudo_geodesic_energy_part(pg_JTJ, pg_mJTF);
        H += weight_pseudo_geodesic_energy * pg_JTJ;
        B += weight_pseudo_geodesic_energy * pg_mJTF;
    }
    if(enable_strip_width_energy){
        spMat sw_JTJ;
        Efunc sw_mJTF;
        assemble_solver_strip_width_part(sw_JTJ, sw_mJTF);
        H += weight_strip_width * sw_JTJ;
        B += weight_strip_width * sw_mJTF;
    }
    assert(H.rows() == vnbr);
    assert(H.cols() == vnbr);
    // std::cout<<"before solving"<<std::endl;
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver(H);
    assert(solver.info() == Eigen::Success);
    // std::cout<<"solved successfully"<<std::endl;
    Eigen::VectorXd dx = solver.solve(B).eval();
    dx*=0.75;
    // std::cout << "step length " << dx.norm() << std::endl;
    double level_set_step_length = dx.norm();
    // double inf_norm=dx.cwiseAbs().maxCoeff();
    if (enable_pseudo_geodesic_energy && level_set_step_length > max_step_length)
    {
        dx *= max_step_length / level_set_step_length;
    }
    fvalues += dx;
    double energy_laplacian=(Dlps*fvalues).norm();
    double energy_biharmonic=fvalues.transpose()* QcH*fvalues;
    double energy_boundary=(bcfvalue).norm();
    std::cout << "energy: harm " << energy_biharmonic << ", bnd " << energy_boundary << ", ";
    if (enable_pseudo_geodesic_energy)
    {
        double energy_pg = (PeWeight.asDiagonal() * EPEvalue).norm();
        std::cout << "pg, " << energy_pg << ", ";
        // std::cout<<"PeWeight\n"<<PeWeight<<std::endl;
    }
    if (enable_strip_width_energy)
    {
      
        Eigen::VectorXd ener = gfvalue.rowwise().norm();
        ener=ener.asDiagonal()*ener;
        Eigen::VectorXd wds = Eigen::VectorXd::Ones(fnbr)*strip_width*strip_width;
        
        ener-=wds;
        double stp_energy=Eigen::VectorXd(ener).dot(ener);
        std::cout << "strip, " << stp_energy << ", ";
    }
    step_length=dx.norm();
    std::cout<<"step "<<step_length<<std::endl;


}