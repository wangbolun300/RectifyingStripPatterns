#include <lsc/basic.h>
#include <lsc/tools.h>
#include <igl/lscm.h>
#include <igl/harmonic.h>
#include <igl/boundary_loop.h>
#include <igl/grad.h>
#include <igl/hessian.h>
#include <igl/curved_hessian_energy.h>

lsTools::lsTools(CGMesh &mesh)
{
    init(mesh);
}
void lsTools::init(CGMesh &mesh)
{
    MP.mesh2Matrix(mesh, V, F);
    MP.meshEdges(mesh, E);
    lsmesh = mesh;
    std::cout << "parametrization start, getting boundary loop, F size " << F.rows() << std::endl;
    igl::boundary_loop(F, bnd);
    std::cout << "parametrization start, bnd size " << bnd.size() << std::endl;
    Eigen::VectorXi b(2, 1); // IDs of the two fixed points on the boundary.
    b(0) = bnd(0);
    b(1) = bnd(bnd.size() / 2);
    Eigen::MatrixXd bc(2, 2);
    bc << 0, 0, 1, 0;

    // LSCM parametrization
    std::cout << "parametrization start, using LSCM " << std::endl;
    igl::lscm(V, F, b, bc, paras);
    assert(paras.cols() == 2);
    std::cout << "parametrization finished, parameter matrix size " << paras.rows() << " x " << paras.cols() << std::endl;

    igl::cotmatrix(V, F, Dlps);
    igl::curved_hessian_energy(V, F, QcH);
    initialize_mesh_properties();
    // restore the original mesh info: edge lengths, aabb tree, and the original vertices
    int enbr = E.rows();
    ElStored.resize(enbr);
    for (int i = 0; i < enbr; i++)
    {
        int vid0 = E(i, 0);
        int vid1 = E(i, 1);
        ElStored[i] = (V.row(vid0) - V.row(vid1)).norm();
    }
    Vstored = V;
    aabbtree.init(V, F);
    Nstored = norm_v;
}

void lsTools::prepare_level_set_solving(const EnergyPrepare &Energy_initializer)
{
    weight_mass = Energy_initializer.weight_gravity;
    weight_laplacian = Energy_initializer.weight_lap;
    weight_boundary = Energy_initializer.weight_bnd;
    enable_pseudo_geodesic_energy = Energy_initializer.solve_pseudo_geodesic;
    pseudo_geodesic_target_angle_degree = Energy_initializer.target_angle;
    weight_pseudo_geodesic_energy = Energy_initializer.weight_pg;
    max_step_length=Energy_initializer.max_step_length;
    enable_strip_width_energy=Energy_initializer.solve_strip_width_on_traced;
    weight_strip_width = Energy_initializer.weight_strip_width;
    enable_inner_vers_fixed=Energy_initializer.enable_inner_vers_fixed;
    enable_functional_angles=Energy_initializer.enable_functional_angles;
    if(Energy_initializer.enable_functional_angles){
        
        pseudo_geodesic_target_min_angle_degree=Energy_initializer.target_min_angle;
        pseudo_geodesic_target_max_angle_degree=Energy_initializer.target_max_angle;
    }
    pseudo_geodesic_start_angle_degree=Energy_initializer.start_angle;
    enable_boundary_angles=Energy_initializer.enable_boundary_angles;
    enable_extreme_cases = Energy_initializer.enable_extreme_cases;
     
    Given_Const_Direction= Energy_initializer.Given_Const_Direction;
  
}
void lsTools::prepare_mesh_optimization_solving(const MeshEnergyPrepare& initializer){
    Mesh_opt_max_step_length=initializer.Mesh_opt_max_step_length;
    weight_Mesh_pesudo_geodesic=initializer.weight_Mesh_pesudo_geodesic;
    weight_Mesh_smoothness=initializer.weight_Mesh_smoothness;
    pseudo_geodesic_target_angle_degree=initializer.target_angle;
    weight_Mesh_edgelength = initializer.weight_Mesh_edgelength;
    enable_extreme_cases = initializer.enable_extreme_cases;
    weight_mass=initializer.weight_mass;
    Given_Const_Direction= initializer.Given_Const_Direction;
}
void lsTools::convert_paras_as_meshes(CGMesh &output)
{
    Eigen::MatrixXd paras_visual; // 3d visulazation of the paras, nx3.
    Eigen::MatrixXd zeros = Eigen::MatrixXd::Zero(paras.rows(), 1);
    paras_visual.resize(paras.rows(), 3);
    paras_visual << paras, zeros;
    output = lsmesh;
    int i = 0;
    for (CGMesh::VertexIter viter = output.vertices_begin(); viter != output.vertices_end(); ++viter)
    {
        CGMesh::Point pt;
        pt[0] = paras(i, 0);
        pt[1] = paras(i, 1);
        pt[2] = 0;
        i++;
        output.set_point(viter, pt);
    }
}
// each element is the sqrt of area
spMat get_uniformed_mass(spMat& mass_in){
    spMat mass = mass_in;
    for (int i = 0; i < mass.rows(); i++)
    {
        double area = mass.coeffRef(i, i);
        mass.coeffRef(i, i) = sqrt(area);
    }
    int nbr = mass.rows();
    double avg = 0; // avrage area
    for (int i = 0; i < nbr; i++)
    {
        avg += mass.coeffRef(i, i);
    }
    avg /= nbr;
    return mass / avg;
}

// get the area and normal
void lsTools::get_mesh_normals_per_face()
{
    int fsize = F.rows();
    norm_f.resize(fsize, 3);
    areaF.resize(F.rows());
    areaPF.resize(F.rows());
    for (int i = 0; i < fsize; i++)
    {
        int id0 = F(i, 0);
        int id1 = F(i, 1);
        int id2 = F(i, 2);
        Eigen::Vector3d cross = Eigen::Vector3d(V.row(id0) - V.row(id1)).cross(Eigen::Vector3d(V.row(id0) - V.row(id2)));
        norm_f.row(i) = cross.normalized();
        areaF(i) = cross.norm() / 2;
        Eigen::Vector2d p0 = paras.row(id0), p1 = paras.row(id1), p2 = paras.row(id2);
        double le0 = (p2 - p1).norm();
        double le1 = (p0 - p2).norm();
        double le2 = (p0 - p1).norm();
        double lp = (le0 + le1 + le2) / 2;
        areaPF(i) = sqrt(lp * (lp - le0) * (lp - le1) * (lp - le2));
    }
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_BARYCENTRIC, mass);
    mass_uniform = get_uniformed_mass(mass);
}

void lsTools::get_mesh_angles()
{
    int fsize = F.rows();
    angF.resize(fsize, 3);
    for (int i = 0; i < fsize; i++)
    {
        int id0 = F(i, 0);
        int id1 = F(i, 1);
        int id2 = F(i, 2);
        Eigen::Vector3d v01 = V.row(id0) - V.row(id1);
        Eigen::Vector3d v12 = V.row(id1) - V.row(id2);
        Eigen::Vector3d v20 = V.row(id2) - V.row(id0);
        double l01 = v01.norm();
        double l12 = v12.norm();
        double l20 = v20.norm();
        double cos0 = v01.dot(-v20) / (l01 * l20);
        double cos1 = v12.dot(-v01) / (l01 * l12);
        double cos2 = v20.dot(-v12) / (l12 * l20);
        double angle0 = acos(cos0);
        double angle1 = acos(cos1);
        double angle2 = acos(cos2);
        angF.row(i) << angle0, angle1, angle2;
    }
}

// this method use angle of the triangles to get a better smoothness of the mesh
void lsTools::get_mesh_normals_per_ver()
{
    norm_v.resize(V.rows(), 3);
    for (int i = 0; i < V.rows(); i++)
    {
        norm_v.row(i) << 0, 0, 0;
        CGMesh::VertexHandle vh = lsmesh.vertex_handle(i); // for each vertex, iterate all the faces
        for (CGMesh::VertexFaceIter vf_it = lsmesh.vf_begin(vh); vf_it != lsmesh.vf_end(vh); ++vf_it)
        {
            int fid = vf_it.handle().idx();
            int pinf = -1;
            for (int j = 0; j < 3; j++)
            {
                if (F(fid, j) == i)
                {
                    pinf = j;
                    break;
                }
            }
            if (pinf == -1)
            {
                std::cout << "ERROR: no find correct one ring face. \npid: " << vh.idx() << " or " << i << std::endl;
                std::cout << "fid: " << fid << std::endl;
                std::cout << "the vertices of this face: " << F.row(fid) << std::endl
                          << std::endl;
            }
            assert(pinf > -1);
            norm_v.row(i) += angF(fid, pinf) * norm_f.row(fid);
        }
        norm_v.row(i) = Eigen::Vector3d(norm_v.row(i)).normalized();
    }
}

void lsTools::get_face_rotation_matices()
{
    int fsize = F.rows();
    Rotate.resize(fsize);
    for (int i = 0; i < fsize; i++)
    {
        Rotate[i].resize(3, 3);
        Eigen::Vector3d norm = norm_f.row(i);
        double x = norm(0);
        double y = norm(1);
        double z = norm(2);
        Rotate[i] << x * x, x * y - z, x * z + y,
            x * y + z, y * y, y * z - x,
            x * z - y, y * z + x, z * z;
    }
}
void lsTools::get_vertex_rotation_matices()
{
    int vsize = V.rows();
    RotateV.resize(vsize);
    for (int i = 0; i < vsize; i++)
    {
        RotateV[i].resize(3, 3);
        Eigen::Vector3d norm = norm_v.row(i);
        double x = norm(0);
        double y = norm(1);
        double z = norm(2);
        RotateV[i] << x * x, x * y - z, x * z + y,
            x * y + z, y * y, y * z - x,
            x * z - y, y * z + x, z * z;
    }
}
void lsTools::get_rotated_edges_for_each_face()
{
    Erotate.resize(F.rows());
    int debug_counter = 0;
    // the 3 edges are in the order of the openmesh face-edge iterator provides us
    for (CGMesh::FaceIter f_it = lsmesh.faces_begin(); f_it != (lsmesh.faces_end()); ++f_it)
    {
        int fid = f_it.handle().idx();
        int ecounter = 0;
        for (CGMesh::FaceHalfedgeIter fh_it = lsmesh.fh_begin(f_it); fh_it != (lsmesh.fh_end(f_it)); ++fh_it)
        {

            CGMesh::HalfedgeHandle heh = fh_it.current_halfedge_handle();
            int vid1 = lsmesh.from_vertex_handle(heh).idx();
            int vid2 = lsmesh.to_vertex_handle(heh).idx();
            if (debug_counter < 100)
            { // debug tool
                debug_counter++;
                assert(lsmesh.face_handle(heh).idx() == fid);
                // std::cout<<"the face id is "<<fid<<", the corresponding edge is "<<vid1<<" "<<vid2<<std::endl;
                auto opposite = lsmesh.opposite_face_handle(heh);
                // if (lsmesh.is_boundary(lsmesh.from_vertex_handle(heh)) && lsmesh.is_boundary(lsmesh.to_vertex_handle(heh)))
                // {
                //     // std::cout<<"it is already boundary"<<std::endl;
                // }
                // else
                // {
                //     // std::cout<<"opposite is "<<opposite.idx();
                //     // std::cout<<", the corresponding edges are "<<F.row(opposite.idx())<<std::endl<<std::endl;
                // }
            }
            // if(1)
            // //if(!(vid1==F(fid,0)||vid1==F(fid,1)||vid1==F(fid,2))||!(vid2==F(fid,0)||vid2==F(fid,1)||vid2==F(fid,2)))
            // {
            //     std::cout << "F and edge not match.\nF " << fid << ", vers, " << F.row(fid) << std::endl;
            //     std::cout << "E vers: " << vid1<<", "<<vid2 << std::endl
            //               << std::endl;
            // }

            // the rotated vector of V[1]-V[0] is obtained here:
            // the order of the 3 edges are:
            // V(F(fid,0))-V(F(fid,2)),V(F(fid,1))-V(F(fid,0)),V(F(fid,2))-V(F(fid,1))
            Erotate[fid][ecounter] = Eigen::Vector3d(Rotate[fid] * Eigen::Vector3d(V.row(vid2) - V.row(vid1)));
            assert(fid < Erotate.size());
            assert(vid1 == F(fid, 0) || vid1 == F(fid, 1) || vid1 == F(fid, 2));
            assert(vid2 == F(fid, 0) || vid2 == F(fid, 1) || vid2 == F(fid, 2));
            assert(vid1 != vid2);
            if (ecounter == 0)
            {
                assert(vid1 == F(fid, 2) && vid2 == F(fid, 0));
            }
            if (ecounter == 1)
            {
                assert(vid1 == F(fid, 0) && vid2 == F(fid, 1));
            }
            if (ecounter == 2)
            {
                assert(vid1 == F(fid, 1) && vid2 == F(fid, 2));
            }
            ecounter++;
        }
    }
}
// bool dbg01=true;

// output is the 3 matrices of fx, fy and fz.
void lsTools::gradient_v2f(std::array<spMat, 3> &output)
{

    int vsize = V.rows();
    int fsize = F.rows();
    output[0].resize(fsize, vsize); // each face has a output value, which is the combination of function input of each vertex
    output[1].resize(fsize, vsize); // each face has a output value, which is the combination of function input of each vertex
    output[2].resize(fsize, vsize); // each face has a output value, which is the combination of function input of each vertex

    spMat G; // gradient matrix
    igl::grad(V, F, G);
    // spMat Gtrans=G.transpose();// transport of G to get 3*vnbr cols
    // std::vector<Trip> triplets0, triplets1, triplets2;
    // triplets0.reserve(vsize*fsize);
    // triplets1.reserve(vsize*fsize);
    // triplets2.reserve(vsize*fsize);
    // for (int i = 0; i < fsize; i++)
    // {
    //     int colid0 = 3 * i;
    //     int colid1 = colid0 + 1;
    //     int colid2 = colid0 + 2;
    //     for (spMat::InnerIterator it(Gtrans, colid0); it; ++it)
    //     {
    //         triplets0.push_back(Trip(i,it.row(), it.value()));
    //     }
    //     for (spMat::InnerIterator it(Gtrans, colid1); it; ++it)
    //     {
    //         triplets1.push_back(Trip(i,it.row(), it.value()));
    //     }
    //     for (spMat::InnerIterator it(Gtrans, colid2); it; ++it)
    //     {
    //         triplets2.push_back(Trip(i,it.row(), it.value()));
    //     }
    // }
    // output[0].setFromTriplets(triplets0.begin(), triplets0.end());
    // output[1].setFromTriplets(triplets1.begin(), triplets1.end());
    // output[2].setFromTriplets(triplets2.begin(), triplets2.end());
    output[0] = G.topRows(fsize);
    output[1] = G.middleRows(fsize, fsize);
    output[2] = G.bottomRows(fsize);
}

// the output is a weight matrix acorrding to angles, size is vsize*fsize
void lsTools::gradient_f2v(spMat &output)
{
    int vsize = V.rows();
    int fsize = F.rows();
    std::vector<Trip> triplets, triplets_as;
    triplets.reserve(vsize * fsize);
    triplets_as.reserve(vsize);
    output.resize(vsize, fsize);
    spMat ASI(vsize, vsize); // a diagnal matrix with 1/(sum of the angles) as elements
    for (int i = 0; i < V.rows(); i++)
    {
        CGMesh::VertexHandle vh = lsmesh.vertex_handle(i); // for each vertex, iterate all the faces
        double anglesum = 0;
        for (CGMesh::VertexFaceIter vf_it = lsmesh.vf_begin(vh); vf_it != lsmesh.vf_end(vh); ++vf_it)
        {
            assert(vh.idx() == i);
            int fid = vf_it.handle().idx();
            int vvid = -1;
            for (int vv = 0; vv < 3; vv++)
            {
                if (i == F(fid, vv))
                {
                    vvid = vv;
                }
            }
            assert(vvid > -1);
            double angle = angF(fid, vvid);
            assert(angle > 0);
            anglesum += angle;
            triplets.push_back(Trip(i, fid, angle));
        }
        triplets_as.push_back(Trip(i, i, 1 / anglesum));
    }
    output.setFromTriplets(triplets.begin(), triplets.end());
    ASI.setFromTriplets(triplets_as.begin(), triplets_as.end());
    output = ASI * output;
}
void lsTools::get_bnd_vers_and_handles()
{
    int vnbr = V.rows();
    int fnbr = F.rows();
    std::vector<Trip> triplets;
    triplets.reserve(vnbr);
    IVids.clear();
    IVids.reserve(vnbr);
    Eigen::VectorXi flags = Eigen::VectorXi::Ones(vnbr) * -1; // first by default set all vertices as non-boundary vertices
    for (CGMesh::VertexIter v_it = lsmesh.vertices_begin(); v_it != (lsmesh.vertices_end()); ++v_it)
    {
        int vid = v_it.handle().idx();
        if (lsmesh.is_boundary(v_it.handle()))
        { // find all the boundary vertices

            
        }
        else
        {
            flags(vid) = IVids.size();
            IVids.push_back(vid);
        }
    }
    InnerV = flags;
    for(int i = 0;i<IVids.size();i++){
        int vid = IVids[i];
        if(InnerV[vid]!=i){
            std::cout<<"The inner vertex mapping is wrong!"<<std::endl;
        }
    }
    // get all the boundary halfedges
    std::vector<CGMesh::HalfedgeHandle> hdls;
    for (CGMesh::EdgeIter e_it = lsmesh.edges_begin(); e_it != lsmesh.edges_end(); ++e_it)
    {
        OpenMesh::HalfedgeHandle hh = lsmesh.halfedge_handle(e_it, 0);
        CGMesh::HalfedgeHandle ophe=lsmesh.opposite_halfedge_handle(hh);
    
        if (lsmesh.face_handle(hh).idx()<0)
        {
            // this is a boundary halfedge.
            hdls.push_back(hh);
        }
        else if (lsmesh.face_handle(ophe).idx()<0)
        {
            // this is a boundary halfedge.
            hdls.push_back(hh);
        }
    }
    Boundary_Edges = hdls;
    for (int i = 0; i < Boundary_Edges.size(); i++)// orient the boundary edges
    {
        Boundary_Edges[i] = boundary_halfedge(lsmesh, Boundary_Edges[i]);
    }
    std::cout<<"Boundary Edges Nbr: "<<Boundary_Edges.size()<<std::endl;
    solve_mean_value_laplacian_mat(lsmesh, IVids, MVLap);
}
void lsTools::get_function_gradient_vertex()
{

    gradient_v2f(gradVF);
    spMat f2vmat;
    gradient_f2v(f2vmat);
    gradV[0] = f2vmat * gradVF[0];
    gradV[1] = f2vmat * gradVF[1];
    gradV[2] = f2vmat * gradVF[2];
}

void lsTools::get_function_hessian_vertex()
{

    // the code below is to use derivate of derivate
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            HessianV[i][j] = gradV[j] * gradV[i];
        }
    }
}

void lsTools::get_gradient_hessian_values(const Eigen::VectorXd& func, Eigen::MatrixXd& Vgrad, Eigen::MatrixXd& Fgrad)
{

    int vsize = V.rows();
    int fsize = F.rows();
    Vgrad.resize(vsize, 3); // gradient on v
    Fgrad.resize(fsize, 3); // gradient on f
    
    // gv
    Eigen::VectorXd gvx = gradV[0] * func;
    Eigen::VectorXd gvy = gradV[1] * func;
    Eigen::VectorXd gvz = gradV[2] * func;
    Vgrad.col(0) = gvx;
    Vgrad.col(1) = gvy;
    Vgrad.col(2) = gvz;

    // gf
    Eigen::VectorXd gfx = gradVF[0] * func;
    Eigen::VectorXd gfy = gradVF[1] * func;
    Eigen::VectorXd gfz = gradVF[2] * func;
    Fgrad.col(0) = gfx;
    Fgrad.col(1) = gfy;
    Fgrad.col(2) = gfz;
}
void lsTools::show_level_set(Eigen::VectorXd &val)
{
    val = fvalues;
}


void lsTools::show_current_reference_points(Eigen::MatrixXd &pts)
{
    int size = refids.size();
    pts.resize(size, 3);
    for (int i = 0; i < size; i++)
    {
        pts.row(i) = V.row(refids[i]);
    }
}
void lsTools::show_1_order_derivate(Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, Eigen::MatrixXd &E2, Eigen::MatrixXd &E3, double ratio)
{

    E0 = V + Deriv1[0] * ratio;
    E1 = V - Deriv1[0] * ratio;
    E2 = V + Deriv1[1] * ratio;
    E3 = V - Deriv1[1] * ratio;
}
void lsTools::show_vertex_normal(Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, double ratio)
{
    E0.resize(V.rows(), 3);
    E1 = E0;
    E0 = V;
    E1 = E0 + norm_v * ratio;
}

// this version use libigl
void lsTools::assemble_solver_laplacian_part(spMat &H, Efunc &B)
{
    // laplacian matrix

    spMat JTJ = Dlps.transpose() * Dlps;
    Efunc mJTF = dense_vec_to_sparse_vec(-JTJ * fvalues);

    H = mass * JTJ;
    B = mass * mJTF;
    assert(B.size() == V.rows());
    assert(H.rows() == H.cols() && H.rows() == V.rows());
}
void lsTools::assemble_solver_biharmonic_smoothing(const Eigen::VectorXd& func, spMat &H, Eigen::VectorXd &B){
    H = QcH;
    B = -QcH * func;
    assert(B.size() == V.rows());
    assert(H.rows() == H.cols() && H.rows() == V.rows());
}
// assemble matrices for \sum{area(i)*(||gradient(i)||-strip_width)^2}
void lsTools::assemble_solver_strip_width_part(const Eigen::MatrixXd& GradFValue, spMat &H, Eigen::VectorXd &B)
{
    int vsize = V.rows();
    int fsize = F.rows();
    // std::cout<<"compute J"<<std::endl;
    spMat J = (GradFValue.col(0)).asDiagonal() * gradVF[0] + (GradFValue.col(1)).asDiagonal() * gradVF[1] + (GradFValue.col(2)).asDiagonal() * gradVF[2];
    J *= 2;
    // std::cout<<"compute JTJ"<<std::endl;
    spMat JTJ = J.transpose() *areaF.asDiagonal()* J;
    // std::cout<<"compute f"<<std::endl;
    Eigen::VectorXd f = GradFValue.col(0).asDiagonal() * GradFValue.col(0) + GradFValue.col(1).asDiagonal() * GradFValue.col(1) + GradFValue.col(2).asDiagonal() * GradFValue.col(2);
    f -= Eigen::VectorXd::Ones(fsize) * strip_width * strip_width;

    Eigen::VectorXd mJTF_dense = -J.transpose()*areaF.asDiagonal() * f;
    H = JTJ;
    B = mJTF_dense;
}
void lsTools::assemble_solver_mesh_strip_width(spMat &JTJ, Eigen::VectorXd &B, Eigen::VectorXd &energy){
    
}
void lsTools::get_all_the_edge_normals()
{
    int ne = lsmesh.n_edges();
    norm_e.resize(ne, 3);
    ActE.resize(lsmesh.n_edges());
    for (CGMesh::EdgeIter e_it = lsmesh.edges_begin(); e_it != lsmesh.edges_end(); ++e_it)
    {
        int eid = e_it.handle().idx();
        OpenMesh::HalfedgeHandle hh_ori = lsmesh.halfedge_handle(e_it, 0);
        OpenMesh::HalfedgeHandle hh_opp = lsmesh.halfedge_handle(e_it, 1);
        assert(hh_ori == lsmesh.opposite_halfedge_handle(hh_opp));
        int fid1 = lsmesh.face_handle(hh_ori).idx();
        int fid2 = lsmesh.face_handle(hh_opp).idx();
        int vid1=lsmesh.from_vertex_handle(hh_ori).idx();
        int vid2=lsmesh.to_vertex_handle(hh_ori).idx();
        Eigen::Vector3d n1 = norm_v.row(vid1);
        Eigen::Vector3d n2 = norm_v.row(vid2);
        
        // if (fid1 >= 0)
        // {
        //     n1 = norm_f.row(fid1);
        // }
        // if (fid2 >= 0)
        // {
        //     n2 = norm_f.row(fid2);
        // }
        Eigen::Vector3d direction = n1 + n2;
        direction = direction / direction.norm();
        norm_e.row(eid) = direction;
        int coplanarity=edge_is_coplanar(lsmesh,e_it.handle(),V,F);
        if(coplanarity==0){// not on boundary and not co-planar
            ActE.coeffRef(eid)=1;
        }
    }
    Actid.reserve(lsmesh.n_edges());
    for(int i=0;i<ActE.size();i++){
        if(ActE.coeffRef(i)==1){
            Actid.push_back(i);
        }
    }
    // std::cout<<"The edge normal vectors are solved for tracing method"<<std::endl;
}


double get_t_of_segment(const Eigen::Vector3d &ver, const Eigen::Vector3d &start, const Eigen::Vector3d &end)
{
    double dis1 = (ver - start).norm();
    double dis2 = (end - start).norm();
    assert(start != end);
    double t=dis1/dis2;
    return t;
}
double get_t_of_value(const double &ver, const double &start, const double &end)
{
    double dis1 = ver - start;
    double dis2 = end - start;
    assert(start != end);
    double t=dis1/dis2;
    return t;
}
Eigen::Vector2d get_2d_ver_from_t(const double t, const Eigen::Vector2d& start, const Eigen::Vector2d& end){
    return start+t*(end-start);
}
Eigen::Vector3d get_3d_ver_from_t(const double t, const Eigen::Vector3d& start, const Eigen::Vector3d& end){
    return start+t*(end-start);
}

void lsTools::initialize_level_set_accroding_to_parametrization(){
    int lssize=assigned_trace_ls.size();
    double value0=assigned_trace_ls[0];
    double value1=assigned_trace_ls[lssize-1];
    int ver0id=lsmesh.from_vertex_handle(trace_hehs[0][0]).idx();
    int ver1id=lsmesh.to_vertex_handle(trace_hehs[0][0]).idx();
    int ver2id=lsmesh.from_vertex_handle(trace_hehs[lssize-1][0]).idx();
    int ver3id=lsmesh.to_vertex_handle(trace_hehs[lssize-1][0]).idx();

    Eigen::Vector3d ver0=V.row(ver0id);
    Eigen::Vector3d ver1=V.row(ver1id);
    Eigen::Vector3d ver2=V.row(ver2id);
    Eigen::Vector3d ver3=V.row(ver3id);
    assert(paras.rows()==V.rows());
    Eigen::Vector2d ver2d0=paras.row(ver0id);
    Eigen::Vector2d ver2d1=paras.row(ver1id);
    Eigen::Vector2d ver2d2=paras.row(ver2id);
    Eigen::Vector2d ver2d3=paras.row(ver3id);
    Eigen::Vector3d ver_value0=trace_vers[0][0];
    Eigen::Vector3d ver_value1=trace_vers[lssize-1][0];

    double t1=get_t_of_segment(ver_value0,ver0,ver1);
    double t2=get_t_of_segment(ver_value1,ver2,ver3);
    Eigen::Vector2d para_value0=get_2d_ver_from_t(t1,ver2d0,ver2d1);
    Eigen::Vector2d para_value1=get_2d_ver_from_t(t1,ver2d2,ver2d3);
    Eigen::Vector2d trans=(value1-value0)*(para_value1-para_value0)/(para_value1-para_value0).dot(para_value1-para_value0);
    Eigen::VectorXd dupvalue=duplicate_valus(value0,V.rows());
    Eigen::MatrixXd dupmatrix=duplicate_vector(para_value0,V.rows());
    Eigen::MatrixXd tmp=paras-dupmatrix;// matrix of Vi-V0
    Eigen::VectorXd result=tmp*trans;
    result=result-dupvalue;
    fvalues=result;
}
void lsTools::Trace_One_Guide_Pseudo_Geodesic() {

}


void lsTools::initialize_level_set_by_tracing(const TracingPrepare& Tracing_initializer){
    // cylinder_open_example(5, 10, 50, 30);
    // exit(0);
    trace_vers.clear();
    trace_hehs.clear();
    assigned_trace_ls.clear();
    double target_angle = Tracing_initializer.target_angle;// 
    double start_angel = Tracing_initializer.start_angle;
    trace_start_angle_degree=start_angel;
    double threadshold_angel_degree = Tracing_initializer.threadshold_angel_degree; // threadshold for checking mesh boundary corners
    int nbr_itv = Tracing_initializer.every_n_edges; // every nbr_itv boundary edges we shoot one curve
    int which_segment=Tracing_initializer.which_boundary_segment;
    
    std::vector<std::vector<CGMesh::HalfedgeHandle>> boundaries;
    split_mesh_boundary_by_corner_detection(lsmesh, V, threadshold_angel_degree,Boundary_Edges, boundaries);
    std::cout<<"get the boundary segments, how many: "<<boundaries.size()<<std::endl;
    if (nbr_itv < 1)
    {
        std::cout << "Please set up the parameter nbr_itv " << std::endl;
        return;
    }
    if(which_segment>=boundaries.size()){
        std::cout<<"Please set up correct boundary id, the maximal boundary segment number is "<<boundaries.size()<<std::endl;
        return;
    }
    std::vector<CGMesh::HalfedgeHandle> boundary_segment=boundaries[which_segment];
    tracing_start_edges=boundary_segment;
    std::cout<<"the number of edges on this segment "<<boundary_segment.size()<<std::endl;
    OpenMesh::HalfedgeHandle init_edge = boundary_segment[0];
    
    OpenMesh::HalfedgeHandle checking_edge = init_edge;
    int beid=0;
    double lsvalue=0;
    int curvecount=0;
    while (1){
        checking_edge=boundary_segment[beid];
        std::vector<Eigen::Vector3d> curve;
        std::vector<CGMesh::HalfedgeHandle> handles;
        trace_single_pseudo_geodesic_curve(target_angle, checking_edge, 0.5, start_angel,
                                                                curve, handles);
        curvecount++;
        std::cout<<"one curver traced, size "<<curve.size()<<" the ith "<<curvecount<<std::endl;
        
        assert(curve.size()>0);
        trace_vers.push_back(curve);
        trace_hehs.push_back(handles);
        assigned_trace_ls.push_back(lsvalue);
        lsvalue+=1;
        int nextbeid=beid+nbr_itv;
        if(nextbeid<boundary_segment.size()){
            beid=nextbeid;
        }
        else{
            std::cout<<"run out of this edge, final segment id "<<beid<<std::endl;
            break;
        }
        
    }
    estimate_strip_width_according_to_tracing();
    std::cout<<"check the numbr of curves "<<trace_vers.size()<<std::endl;
   
}

// assign values and angles on one boundary curve splited by the corners
void lsTools::initialize_level_set_by_boundary_assignment(const TracingPrepare& Tracing_initializer){
    // cylinder_open_example(5, 10, 50, 30);
    // exit(0);
    trace_vers.clear();
    trace_hehs.clear();
    assigned_trace_ls.clear();
    double target_angle = Tracing_initializer.target_angle;// 
    double start_angel = Tracing_initializer.start_angle;
    trace_start_angle_degree=start_angel;
    double threadshold_angel_degree = Tracing_initializer.threadshold_angel_degree; // threadshold for checking mesh boundary corners
    int nbr_itv = Tracing_initializer.every_n_edges; // every nbr_itv boundary edges we shoot one curve
    int which_segment=Tracing_initializer.which_boundary_segment;
    
    std::vector<std::vector<CGMesh::HalfedgeHandle>> boundaries;
    split_mesh_boundary_by_corner_detection(lsmesh, V, threadshold_angel_degree,Boundary_Edges, boundaries);
    std::cout<<"get the boundary segments, how many: "<<boundaries.size()<<std::endl;
    if (nbr_itv < 1)
    {
        std::cout << "Please set up the parameter nbr_itv " << std::endl;
        return;
    }
    if(which_segment>=boundaries.size()){
        std::cout<<"Please set up correct boundary id, the maximal boundary segment number is "<<boundaries.size()<<std::endl;
        return;
    }
    std::vector<CGMesh::HalfedgeHandle> boundary_segment=boundaries[which_segment];
    tracing_start_edges=boundary_segment;
    std::cout<<"the number of edges on this segment "<<boundary_segment.size()<<std::endl;
    OpenMesh::HalfedgeHandle init_edge = boundary_segment[0];
    
    OpenMesh::HalfedgeHandle checking_edge = init_edge;
    int beid=0;
    double lsvalue=0;
    int curvecount=0;
    while (1){
        checking_edge = boundary_segment[beid];
        std::vector<Eigen::Vector3d> curve;
        std::vector<CGMesh::HalfedgeHandle> handles;

        curve.clear();
        handles.clear();
        CGMesh::HalfedgeHandle intersected_handle_tmp;
        Eigen::Vector3d intersected_point_tmp;
        double start_point_para=0.5;
        bool found = init_pseudo_geodesic_first_segment(checking_edge, start_point_para, start_angel,
                                                        intersected_handle_tmp, intersected_point_tmp);
       
        if (!found)
        {
            std::cout << "error in initialization the boundary direction" << std::endl;
            return;
        }

        Eigen::Vector3d first_point = get_3d_ver_from_t(start_point_para, V.row(lsmesh.from_vertex_handle(checking_edge).idx()),
                                                        V.row(lsmesh.to_vertex_handle(checking_edge).idx()));

        curve.push_back(first_point);
        handles.push_back(checking_edge);
        curve.push_back(intersected_point_tmp);
        handles.push_back(intersected_handle_tmp);
        curvecount++;
        assert(curve.size()>0);
        trace_vers.push_back(curve);
        trace_hehs.push_back(handles);
        assigned_trace_ls.push_back(lsvalue);
        lsvalue+=1;
        int nextbeid=beid+nbr_itv;
        if(nextbeid<boundary_segment.size()){
            beid=nextbeid;
        }
        else{
            break;
        }
        
    }
    estimate_strip_width_according_to_tracing();
    std::cout<<"the numbr of values assigned to this boundary:  "<<trace_vers.size()<<std::endl;
   
}


// assign values and angles on a boundary curve truncated by two selected boundary edges
void lsTools::initialize_level_set_by_select_boundary_segment(const TracingPrepare& Tracing_initializer, bool trace){
    // cylinder_open_example(5, 10, 50, 30);
    // exit(0);
    trace_vers.clear();
    trace_hehs.clear();
    assigned_trace_ls.clear();
    double target_angle = Tracing_initializer.target_angle;// 
    double start_angel = Tracing_initializer.start_angle;
    trace_start_angle_degree=start_angel;
    double threadshold_angel_degree = Tracing_initializer.threadshold_angel_degree; // threadshold for checking mesh boundary corners
    int nbr_itv = Tracing_initializer.every_n_edges; // every nbr_itv boundary edges we shoot one curve
    int which_segment=Tracing_initializer.which_boundary_segment;
    int start_he = Tracing_initializer.start_bnd_he;
    int nbr_edges = Tracing_initializer.nbr_edges;


    
    if (nbr_itv < 1)
    {
        std::cout << "Please set up the parameter nbr_itv " << std::endl;
        return;
    }

    std::vector<CGMesh::HalfedgeHandle> boundary_segment;
    select_mesh_boundary_curve_on_boundary_loop(lsmesh, V, which_segment,
                                               start_he, nbr_edges, Boundary_Edges, boundary_segment);
    std::cout<<"the number of edges on this segment "<<boundary_segment.size()<<std::endl;
    OpenMesh::HalfedgeHandle init_edge = boundary_segment[0];
    
    OpenMesh::HalfedgeHandle checking_edge = init_edge;
    int beid=0;
    double lsvalue=0;
    int curvecount=0;
    while (1){
        checking_edge = boundary_segment[beid];
        std::vector<Eigen::Vector3d> curve;
        std::vector<CGMesh::HalfedgeHandle> handles;

        curve.clear();
        handles.clear();

        if (!trace)// 
        {
            CGMesh::HalfedgeHandle intersected_handle_tmp;
            Eigen::Vector3d intersected_point_tmp;
            double start_point_para = 0.5;
            bool found = init_pseudo_geodesic_first_segment(checking_edge, start_point_para, start_angel,
                                                            intersected_handle_tmp, intersected_point_tmp);

            if (!found)
            {
                std::cout << "error in initialization the boundary direction" << std::endl;
                return;
            }

            Eigen::Vector3d first_point = get_3d_ver_from_t(start_point_para, V.row(lsmesh.from_vertex_handle(checking_edge).idx()),
                                                            V.row(lsmesh.to_vertex_handle(checking_edge).idx()));

            curve.push_back(first_point);
            handles.push_back(checking_edge);
            curve.push_back(intersected_point_tmp);
            handles.push_back(intersected_handle_tmp);
        }
        else
        {

            trace_single_pseudo_geodesic_curve(target_angle, checking_edge, 0.5, start_angel,
                                               curve, handles);
        }

        curvecount++;
        assert(curve.size()>0);
        trace_vers.push_back(curve);
        trace_hehs.push_back(handles);
        assigned_trace_ls.push_back(lsvalue);
        lsvalue+=1;
        int nextbeid=beid+nbr_itv;
        if(nextbeid<boundary_segment.size()){
            beid=nextbeid;
        }
        else{
            break;
        }
        
    }
    estimate_strip_width_according_to_tracing();
    std::cout<<"the numbr of values assigned to this boundary:  "<<trace_vers.size()<<std::endl;
   
}
bool find_one_ls_segment_on_triangle(const double value, const Eigen::MatrixXi &F, const Eigen::MatrixXd &V,
                                      const Eigen::VectorXd &fvalues, const int fid, Eigen::Vector3d &E0, Eigen::Vector3d &E1)
{
    int vid0 = F(fid, 0);
    int vid1 = F(fid, 1);
    int vid2 = F(fid, 2);

    std::vector<std::array<int, 2>> eids;
    eids.reserve(2);
    double fv0 = fvalues[vid0];
    double fv1 = fvalues[vid1];
    double fv2 = fvalues[vid2];
    std::array<int, 2> tmp;
    // check [0, 1)
    if ((value <= fv0 && value > fv1) || (value < fv1 && value >= fv0))
    {
        tmp[0] = vid0;
        tmp[1] = vid1;
        eids.push_back(tmp);
    }
    // check [1, 2)
    if ((value <= fv1 && value > fv2) || (value < fv2 && value >= fv1))
    {
        tmp[0] = vid1;
        tmp[1] = vid2;
        eids.push_back(tmp);
    }
    // check [2, 0)
    if ((value <= fv2 && value > fv0) || (value < fv0 && value >= fv2))
    {
        tmp[0] = vid2;
        tmp[1] = vid0;
        eids.push_back(tmp);
    }
    if (eids.empty() || eids.size() == 1)
    { // if no value, or the value is on one vertex, no segment.
        int nb = 0;
        if (value == fv0)
        {
            tmp[nb] = vid0;
            nb++;
        }
        if (value == fv1)
        {
            tmp[nb] = vid1;
            nb++;
        }
        if (value == fv2 && nb == 2) // three vertices with same value, singularity
        {
            return false;
        }
        if (value == fv2)
        {
            tmp[nb] = vid2;
            nb++;
        }
        if (nb == 2)
        {
            E0 = V.row(tmp[0]);
            E1 = V.row(tmp[1]);
            return true;
        }
        return false;
    }
    assert(eids.size()!=3);

    double t = get_t_of_value(value, fvalues[eids[0][0]], fvalues[eids[0][1]]);
    assert(t >= 0 && t <= 1);
    E0 = get_3d_ver_from_t(t, V.row(eids[0][0]), V.row(eids[0][1]));
    t = get_t_of_value(value, fvalues[eids[1][0]], fvalues[eids[1][1]]);
    E1 = get_3d_ver_from_t(t, V.row(eids[1][0]), V.row(eids[1][1]));
    return true;
}
// bool find_one_ls_segment_on_edge(const double value, const Eigen::MatrixXi &E, const Eigen::MatrixXd &V,
//                                  const Eigen::VectorXd &fvalues, const int fid, Eigen::Vector3d &E0, Eigen::Vector3d &E1)
// {

// }
void lsTools::extract_one_curve(const double value, Eigen::MatrixXd& E0, Eigen::MatrixXd &E1){
    std::vector<Eigen::Vector3d> e0list;
    std::vector<Eigen::Vector3d> e1list;
    e0list.reserve(F.rows()/2);
    e1list.reserve(F.rows()/2);
    for(int i=0;i<F.rows();i++){
        Eigen::Vector3d e0tmp;
        Eigen::Vector3d e1tmp;
        bool found=find_one_ls_segment_on_triangle(value,F,V,fvalues,i,e0tmp,e1tmp);
        if(found){
            e0list.push_back(e0tmp);
            e1list.push_back(e1tmp);
        }
    }
    E0=vec_list_to_matrix(e0list);
    E1=vec_list_to_matrix(e1list);
}
void lsTools::extract_levelset_curves(const int nbr, std::vector<Eigen::MatrixXd> &E0, std::vector<Eigen::MatrixXd> &E1){
    E0.clear();
    E1.clear();
    E0.reserve(nbr);
    E1.reserve(nbr);
    double lsmin=fvalues.minCoeff();
    double lsmax=fvalues.maxCoeff();
    double itv=(lsmax-lsmin)/(nbr+1);
    for(int i=0;i<nbr;i++){
        double lsvalue=lsmin+(i+1)*itv;
        assert(lsvalue<lsmax);
        Eigen::MatrixXd e0tmp, e1tmp;
        extract_one_curve(lsvalue, e0tmp,e1tmp);
        assert(e0tmp.rows()>=1);
        E0.push_back(e0tmp);
        E1.push_back(e1tmp);
    }
}
void lsTools::estimate_strip_width_according_to_tracing(){
    if(trace_vers.size()<2){
        std::cout<<"Too few traced curves. Please trace more curves or load more curves"<<std::endl;
        return;
    }
    
    double dis=0;// total lengths
    int tnbr=trace_vers.size();

    
    for(int i=0;i<tnbr-1;i++){
        dis += (trace_vers[i + 1][0] - trace_vers[i][0]).norm();
    }
    double angle_radian = trace_start_angle_degree * LSC_PI / 180.;
    double dis_real=abs(dis*sin(angle_radian));
    double function_value_diff=abs(assigned_trace_ls.back()-assigned_trace_ls[0]);
    strip_width=function_value_diff/dis_real;
}
void lsTools::print_info(const int vid)
{
    if (analizers[0].LocalActInner.size() > 0)
    {
        int ninner = analizers[0].LocalActInner.size();
        int vnbr = V.rows();
        int vm = vid;
        int i = InnerV[vm]; // the ith element in inner list is vm
        
        CGMesh::HalfedgeHandle inhd = analizers[0].heh0[i], outhd = analizers[0].heh1[i];
        int v1 = lsmesh.from_vertex_handle(inhd).idx();
        int v2 = lsmesh.to_vertex_handle(inhd).idx();
        int v3 = lsmesh.from_vertex_handle(outhd).idx();
        int v4 = lsmesh.to_vertex_handle(outhd).idx();

        double t1 = analizers[0].t1s[i];
        double t2 = analizers[0].t2s[i];

        Eigen::Vector3d ver0 = V.row(v1) + (V.row(v2) - V.row(v1)) * t1;
        Eigen::Vector3d ver1 = V.row(vm);
        Eigen::Vector3d ver2 = V.row(v3) + (V.row(v4) - V.row(v3)) * t2;
        Eigen::Vector3d norm = norm_v.row(vm);
        Eigen::Vector3d real_u = norm.cross(ver2 - ver0);
        real_u = real_u.normalized();
        Eigen::Vector3d real_r = (ver1 - ver0).normalized().cross((ver2 - ver1).normalized());
        if (real_r.norm() < 1e-6)
        { // almost a straight line
            std::cout<<"This is almost a straight line"<<std::endl;
        }
        real_r = real_r.normalized();
        double cos_angle = real_r.dot(norm);
        double sin_angle = real_r.dot(real_u);
        std::cout<<"cos, sin "<<cos_angle<<", "<<sin_angle<<std::endl;
        if (sin_angle < 0)
        {
            cos_angle *= -1;
        }
         double angle = acos(cos_angle);
         std::cout << "angle " << 180. / LSC_PI * angle << std::endl;
    }
    // int xloc = vm + vnbr + ninner * 8;
    // int yloc = vm + vnbr + ninner * 9;
    // int zloc = vm + vnbr + ninner * 10;
    // int lnpx = vm + vnbr + ninner * 11;
    // int lnpy = vm + vnbr + ninner * 12;
    // int lnpz = vm + vnbr + ninner * 13;
    // int lbnx = vm + vnbr + ninner * 0;
    // int lbny = vm + vnbr + ninner * 1;
    // int lbnz = vm + vnbr + ninner * 2;
    // if (enable_shading_init)
    // {
    //     xloc = vm + vnbr + ninner * 3;
    //     yloc = vm + vnbr + ninner * 4;
    //     zloc = vm + vnbr + ninner * 5;
    // }
    // Eigen::Vector3d bn = Eigen::Vector3d(Glob_lsvars[lbnx], Glob_lsvars[lbny], Glob_lsvars[lbnz]);
    // std::cout << "the light, " << Glob_lsvars[xloc] << ", " << Glob_lsvars[yloc] << ", " << Glob_lsvars[zloc] << "\n";
    // std::cout << "the pn, " << bn[0] << ", " << bn[1] << ", " << bn[2] << "\n";
    // std::cout << "the binormal, " << Glob_lsvars[xloc] << ", " << Glob_lsvars[yloc] << ", " << Glob_lsvars[zloc] << "\n";
    // double latitude_radian = ShadingLatitude * LSC_PI / 180.;
    // double rho = (LSC_PI / 2 - latitude_radian);
    // Eigen::Matrix3d rotation;
    // rotation << 1, 0, 0,
    //     0, cos(rho), -sin(rho),
    //     0, sin(rho), cos(rho);
    // Eigen::Vector3d direction_ground = Eigen::Vector3d(0, 0, -1); // the ground direction
    // direction_ground = rotation * direction_ground;
    // Eigen::Vector3d light = Eigen::Vector3d(Glob_lsvars[xloc], Glob_lsvars[yloc], Glob_lsvars[zloc]);
    // Eigen::Vector3d pn = Eigen::Vector3d(Glob_lsvars[lnpx], Glob_lsvars[lnpy], Glob_lsvars[lnpz]);

    // std::cout << "the ground, " << direction_ground[0] << ", " << direction_ground[1] << ", " << direction_ground[2] << "\n";

    // std::cout << "light * pn = " << Eigen::Vector3d(Glob_lsvars[xloc], Glob_lsvars[yloc], Glob_lsvars[zloc]).dot(Eigen::Vector3d(Glob_lsvars[lnpx], Glob_lsvars[lnpy], Glob_lsvars[lnpz]))
    // <<", ground * pn = "<<direction_ground.dot(Eigen::Vector3d(Glob_lsvars[lnpx], Glob_lsvars[lnpy], Glob_lsvars[lnpz]))<<std::endl;
    // std::cout<<"(light, pn, ground) = "<<light.cross(pn).dot(direction_ground)<<std::endl;
// }

    // LSAnalizer anl;
    // analysis_pseudo_geodesic_on_vertices(fvalues,anl);
    // int ninner = anl.LocalActInner.size();
    // Eigen::MatrixXd tangents = Eigen::MatrixXd::Zero(ninner, 3);
    // for (int i = 0; i < ninner; i++)
	// {
	// 	if (anl.LocalActInner[i] == false)
	// 	{
	// 		std::cout << "singularity" << std::endl;
	// 		continue;
	// 	}
	// 	int vm = IVids[i];

	// 	CGMesh::HalfedgeHandle inhd = anl.heh0[i], outhd = anl.heh1[i];
	// 	int v1 = lsmesh.from_vertex_handle(inhd).idx();
	// 	int v2 = lsmesh.to_vertex_handle(inhd).idx();
	// 	int v3 = lsmesh.from_vertex_handle(outhd).idx();
	// 	int v4 = lsmesh.to_vertex_handle(outhd).idx();

	// 	double t1 = anl.t1s[i];
	// 	double t2 = anl.t2s[i];
		
	// 	Eigen::Vector3d ver0 = V.row(v1) + (V.row(v2) - V.row(v1)) * t1;
	// 	Eigen::Vector3d ver1 = V.row(vm);
	// 	Eigen::Vector3d ver2 = V.row(v3) + (V.row(v4) - V.row(v3)) * t2;
	// 	// the locations
    //     Eigen::Vector3d tg = ver2 - ver0;
    //     tg = tg.normalized();
        
    //     tangents.row(i)=tg;
    // }
    // Eigen::Vector3d ray = Reference_ray.normalized();
    // Eigen::VectorXd dot = tangents * ray;
    // std::cout<<"the angle between ray and tangents, "<<dot.norm()<<", max, "<<dot.maxCoeff()<<std::endl;
    // std::cout<<"the tangent values\n "<<tangents<<std::endl;
}

#include <igl/gaussian_curvature.h>
void lsTools::show_posi_negt_curvature(Eigen::MatrixXd &color)
{
    Eigen::VectorXd K;
    // Compute integral of Gaussian curvature
    igl::gaussian_curvature(V, F, K);
    color.resize(V.rows(), V.cols());
    for (int i = 0; i < K.size(); i++)
    {
        if (K[i] > 1e-3)
        {
            color.row(i) << 1, 0, 0; // positive, red
        }
        else
        {
            if (K[i] < -1e-3)
            {
                color.row(i) << 0, 1, 0; // negative, green
            }
            else{
                color.row(i) << 0, 0, 1; // zero, blue
            }
        }

    }
    std::cout<<"max Gaussian curvature, "<<K.maxCoeff()<<", min, "<<K.minCoeff()<<std::endl;
}
