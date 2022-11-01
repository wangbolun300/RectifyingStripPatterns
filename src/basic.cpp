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
    
}

void lsTools::prepare_level_set_solving(const EnergyPrepare &Energy_initializer)
{
    weight_mass = Energy_initializer.weight_gravity;
    weight_laplacian = Energy_initializer.weight_lap;
    weight_boundary = Energy_initializer.weight_bnd;
    enable_pseudo_geodesic_energy = Energy_initializer.solve_pseudo_geodesic;
    if (Energy_initializer.solve_pseudo_geodesic)
    {
        weight_pseudo_geodesic_energy = Energy_initializer.weight_pg;
        pseudo_geodesic_target_angle_degree = Energy_initializer.target_angle;
    }
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
// void lsTools::get_mesh_normals_per_ver()
// {
//     norm_v.resize(V.rows(), 3);
//     for (int i = 0; i < V.rows(); i++)
//     {
//         norm_v.row(i) << 0, 0, 0;
//         CGMesh::VertexHandle vh = lsmesh.vertex_handle(i); // for each vertex, iterate all the faces
//         for (CGMesh::VertexFaceIter vf_it = lsmesh.vf_begin(vh); vf_it != lsmesh.vf_end(vh); ++vf_it)
//         {
//             int fid = vf_it.handle().idx();
//             int pinf = -1;
//             for (int j = 0; j < 3; j++)
//             {
//                 if (F(fid, j) == i)
//                 {
//                     pinf = j;
//                     break;
//                 }
//             }
//             if (pinf == -1)
//             {
//                 std::cout << "ERROR: no find correct one ring face. \npid: " << vh.idx() << " or " << i << std::endl;
//                 std::cout << "fid: " << fid << std::endl;
//                 std::cout << "the vertices of this face: " << F.row(fid) << std::endl
//                           << std::endl;
//             }
//             assert(pinf > -1);
//             norm_v.row(i) += angF(fid, pinf) * norm_f.row(fid);
//         }
//         norm_v.row(i) = Eigen::Vector3d(norm_v.row(i)).normalized();
//     }
// }

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
    IVids.reserve(vnbr);
    Eigen::VectorXd flags = Eigen::VectorXd::Ones(vnbr); // first by default set all vertices as non-boundary vertices
    for (CGMesh::VertexIter v_it = lsmesh.vertices_begin(); v_it != (lsmesh.vertices_end()); ++v_it)
    {
        int vid = v_it.handle().idx();
        if (lsmesh.is_boundary(v_it.handle()))
        { // find all the boundary vertices
            
            flags(vid) = 0;
        }
        else{
            IVids.push_back(vid);
        }
    }
    InnerV = flags;
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

void lsTools::get_gradient_hessian_values()
{

    int vsize = V.rows();
    int fsize = F.rows();
    gvvalue.resize(vsize, 3); // gradient on v
    gfvalue.resize(fsize, 3); // gradient on f
    hfvalue.resize(vsize);    // hessian
    // gv
    Eigen::VectorXd gvx = gradV[0] * fvalues;
    Eigen::VectorXd gvy = gradV[1] * fvalues;
    Eigen::VectorXd gvz = gradV[2] * fvalues;
    gvvalue.col(0) = gvx;
    gvvalue.col(1) = gvy;
    gvvalue.col(2) = gvz;

    // gf
    Eigen::VectorXd gfx = gradVF[0] * fvalues;
    Eigen::VectorXd gfy = gradVF[1] * fvalues;
    Eigen::VectorXd gfz = gradVF[2] * fvalues;
    gfvalue.col(0) = gfx;
    gfvalue.col(1) = gfy;
    gfvalue.col(2) = gfz;

    // hessian
    std::array<std::array<Eigen::VectorXd, 3>, 3> fii;
    for (int i = 0; i < 3; i++)
    {
        for (int j = 0; j < 3; j++)
        {
            fii[i][j] = HessianV[i][j] * fvalues;
        }
    }
    for (int k = 0; k < vsize; k++)
    {
        for (int i = 0; i < 3; i++)
        {
            for (int j = 0; j < 3; j++)
            {
                hfvalue[k](i, j) = fii[i][j](k);
            }
        }
    }
    get_grad_norm_jacobian();
}

void lsTools::make_sphere_ls_example(int rowid)
{
    int vsize = V.rows();
    int fsize = F.rows();
    int rnbr = 20; // # vertices in each row
    fvalues.resize(vsize);

    for (int i = 0; i < vsize; i++)
    {
        fvalues(i) = V(i, 2); // make function value = z
    }
    get_gradient_hessian_values();
    refids.clear();
    bool writefile = false;
    std::ofstream fout;
    std::string filename = std::string("/Users/wangb0d/bolun/D/vs/levelset/level-set-curves/build/") + std::to_string(rowid) + "_info.csv";
    if (writefile)
    {
        fout.open(filename);
        fout << "angle,kn,kg\n";
    }

    for (int i = 0; i < rnbr; i++)
    {
        int vid = rowid * rnbr + i;

        refids.push_back(vid);
        Eigen::Matrix2d LmMN;
        LmMN << II_L[vid], -II_M[vid],
            -II_M[vid], II_N[vid];
        Eigen::Matrix<double, 3, 2> rurv, rvru;
        rurv.col(0) = Deriv1[0].row(vid);
        rurv.col(1) = Deriv1[1].row(vid);
        rvru.col(0) = Deriv1[1].row(vid);
        rvru.col(1) = Deriv1[0].row(vid);
        Eigen::Vector3d gradient = gvvalue.row(vid);
        Eigen::Matrix3d hessian = hfvalue[vid];
        Eigen::Matrix<double, 3, 2> mrvru; // [-rv, ru]
        mrvru.col(0) = -1 * Deriv1[1].row(vid);
        mrvru.col(1) = Deriv1[0].row(vid);
        double normgradf = gradient.norm();
        Eigen::MatrixXd inner_left = normgradf * rvru * LmMN * rvru.transpose();
        Eigen::MatrixXd inner_right = mrvru * rurv.transpose() * hessian * rurv * mrvru.transpose();
        double left = gradient.transpose() * inner_left * gradient;
        double right = gradient.transpose() * inner_right * gradient;
        assert(right != 0);
        double b = right / left;
        std::cout << "vid " << vid << "\n";
        std::cout << "the current row nbr " << rowid << " with pid " << i << "\n**b is " << b << "\ngradient, " << gradient.transpose() << "\n";
        Eigen::Vector3d rotgrad = RotateV[vid] * gradient;
        Eigen::Vector3d s1 = hessian * rotgrad;
        double s2 = rotgrad.dot(s1);
        double kg = -(s2) / (normgradf * normgradf * normgradf);
        std::cout << "Kg is " << kg << std::endl;
        std::cout << "Kn (estimate) is " << kg / b << "\n";
        Eigen::Vector3d ru = Deriv1[0].row(vid);
        Eigen::Vector3d rv = Deriv1[1].row(vid);
        double gfrv = gradient.dot(rv);
        double gfru = gradient.dot(ru);
        // Eigen::Matrix2d LMN;
        // LMN << II_L[vid], II_M[vid],
        //     II_M[vid], II_N[vid];
        Eigen::Vector2d IIl = Eigen::Vector2d(gfrv, gfru);
        double up = IIl.transpose() * LmMN * IIl;
        double cosuv = ru.dot(rv) / (ru.norm() * rv.norm());
        double sin2 = 1 - cosuv * cosuv;
        assert(sin2 >= 0);
        double down = ru.norm() * ru.norm() * rv.norm() * rv.norm() * normgradf * normgradf * sin2;
        double kn = up / down;
        std::cout << "Kn (calculate) is " << up / down << std::endl;
        std::cout << "b recalculate is " << kg / kn << std::endl;

        Eigen::Vector3d q = -(gradient.dot(rv)) * ru + gradient.dot(ru) * rv;
        Eigen::Vector3d rf_test = q.normalized() * gradient.norm();
        // std::cout<<"gradf size is "<<gradient.norm()<<"\ngrad test angle diff "<<rf_test.dot(gradient)<<std::endl;
        std::cout << "qnorm calculated and norm " << gradient.norm() * ru.norm() * rv.norm() * sqrt(sin2) << " " << q.norm() << std::endl;
        double ang = atan(kg / kn) * 180 / 3.1415926;
        std::cout << "angle is " << ang << std::endl;
        std::cout << "angle calculated is " << atan(b) * 180 / 3.1415926 << std::endl;
        if (writefile)
        {
            fout << ang << "," << up / down << "," << kg << "\n";
        }
        // std::cout<<"__\nthe"
        std::cout << "\n";
    }
    if (writefile)
    {
        fout.close();
    }
    std::cout << "__________________________________________________" << std::endl;
}
void lsTools::show_level_set(Eigen::VectorXd &val)
{
    val = fvalues;
}

void lsTools::show_gradients(Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, double ratio)
{

    E0 = V + gvvalue * ratio;
    E1 = V - gvvalue * ratio;
    assert(V.rows() == gvvalue.rows());
}

void lsTools::show_hessian(Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, double ratio, int which)
{
    Eigen::MatrixXd direction(V.rows(), 3);
    for (int i = 0; i < V.rows(); i++)
    {
        direction.row(i) = hfvalue[i].row(which);
    }
    // spMat hess_igl;
    // igl::hessian(V,F,hess_igl);
    // std::cout<<"hessian size "<<hess_igl.rows()<<" "<<hess_igl.cols()<<std::endl;

    E0 = V + direction * ratio;
    E1 = V;
}
void lsTools::show_gradient_scalar(Eigen::VectorXd &values, int which)
{
    values = gvvalue.col(which);
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
void lsTools::show_face_gradients(Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, double ratio)
{
    Eigen::MatrixXd fcent(F.rows(), 3);
    for (int i = 0; i < F.rows(); i++)
    {
        fcent.row(i) = (V.row(F(i, 0)) + V.row(F(i, 1)) + V.row(F(i, 2))) / 3;
    }
    Eigen::MatrixXd dirc(F.rows(), 3);
    // spMat G;
    // igl::grad(V,F,G);
    // Eigen::VectorXd total;
    // total=G*fvalues;
    // int counter=0;
    // for(int j=0;j<3;j++){
    // for(int i=0;i<F.rows();i++){

    //         dirc(i,j)=total(counter);
    //         counter++;
    //     }
    // }
    dirc = gfvalue;

    E0 = fcent + dirc * ratio;
    E1 = fcent - dirc * ratio;
    assert(fcent.rows() == dirc.rows());
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

void lsTools::get_grad_norm_jacobian()
{

    int vsize = V.rows();
    
    spMat fx_diag, fy_diag, fz_diag;     // diagnal matrices with gradeint values as matrices
    Eigen::VectorXd fxvec, fyvec, fzvec; // gradient values, fx, fy, fz
    fxvec = gvvalue.col(0);
    fyvec = gvvalue.col(1);
    fzvec = gvvalue.col(2);

    fx_diag = fxvec.asDiagonal();
    fy_diag = fyvec.asDiagonal();
    fz_diag = fzvec.asDiagonal();

    Eigen::VectorXd grad_norm = gvvalue.rowwise().norm(); // take 1/||gradV(i)|| as ith diagnal element
    // invert grad_norm
    for (int i = 0; i < vsize; i++)
    {
        if (grad_norm(i) > 0)
        {
            grad_norm(i) = 1. / grad_norm(i);
        }
    }
    Eigen::MatrixXd grad_norm_mat(vsize, 1);
    grad_norm_mat.col(0) = grad_norm;
    Eigen::DiagonalMatrix<double, Eigen::Dynamic, Eigen::Dynamic> grad_norm_diag = grad_norm_mat.asDiagonal();
    // each row of Dgrad_norm is the jacobian of ||gradV(i)||
    Dgrad_norm = grad_norm_diag * (fx_diag * gradV[0] + fy_diag * gradV[1] + fz_diag * gradV[2]);
    derivates_calculated = true;
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
void lsTools::assemble_solver_biharmonic_smoothing(spMat &H, Efunc &B){
    spMat JTJ=QcH; 
    Efunc mJTF=dense_vec_to_sparse_vec(-JTJ * fvalues);
    H = JTJ;
    B = mJTF;// the mass matrix is already integrated in QcH.
    assert(B.size() == V.rows());
    assert(H.rows() == H.cols() && H.rows() == V.rows());
}
// assemble matrices for \sum{area(i)*(||gradient(i)||-strip_width)^2}
// TODO modify it to be a triangle face based method
void lsTools::assemble_solver_strip_width_part(spMat &H, Efunc &B)
{
    int vsize = V.rows();
    int fsize = F.rows();
    // std::cout<<"compute J"<<std::endl;
    spMat J = (gfvalue.col(0)).asDiagonal() * gradVF[0] + (gfvalue.col(1)).asDiagonal() * gradVF[1] + (gfvalue.col(2)).asDiagonal() * gradVF[2];
    J *= 2;
    // std::cout<<"compute JTJ"<<std::endl;
    spMat JTJ = J.transpose() *areaF.asDiagonal()* J;
    // std::cout<<"compute f"<<std::endl;
    Eigen::VectorXd f = gfvalue.col(0).asDiagonal() * gfvalue.col(0) + gfvalue.col(1).asDiagonal() * gfvalue.col(1) + gfvalue.col(2).asDiagonal() * gfvalue.col(2);
    f -= Eigen::VectorXd::Ones(fsize) * strip_width * strip_width;

    Eigen::VectorXd mJTF_dense = -J.transpose()*areaF.asDiagonal() * f;
    H = JTJ;
    B = dense_vec_to_sparse_vec(mJTF_dense);
}

// // we summarize the problem as:
// // grad(F).transpose()*grad(F).norm()*R*grad(F)-grad(F).transpose()*D.transpose()*hess(F)*D*grad(F).
// void lsTools::assemble_solver_pseudo_geodesic_part(spMat &H, Efunc& B){

//     int vnbr=V.rows();
//     int fnbr=F.rows();
//     std::vector<Eigen::MatrixXd> Rlist(vnbr);// contains pseudo_geodesic_ratio
//     std::vector<Eigen::MatrixXd> Dlist(vnbr);
//     std::vector<Eigen::MatrixXd> Glist(vnbr);
//     for (int i = 0; i < vnbr; i++)
//     {
//         Eigen::MatrixXd rvru(3, 2);
//         rvru.col(0) = Deriv1[1].row(i);
//         rvru.col(1) = Deriv1[0].row(i);
//         Eigen::Matrix<double, 3, 2> rurv;
//         rurv.col(0) = Deriv1[0].row(i);
//         rurv.col(1) = Deriv1[1].row(i);
//         Eigen::Matrix<double, 3, 2> mrvru; // [-rv, ru]
//         mrvru.col(0) = -1 * Deriv1[1].row(i);
//         mrvru.col(1) = Deriv1[0].row(i);
//         Eigen::Matrix2d LmMN;
//         LmMN << II_L[i], -II_M[i],
//             -II_M[i], II_N[i];
//         Eigen::MatrixXd R=pseudo_geodesic_ratio*rvru*LmMN*rvru.transpose();
//         Rlist[i]=R;
//         Eigen::MatrixXd D=rurv*mrvru.transpose();
//         Dlist[i]=D;
//         Eigen::MatirxXd gf(3,1);
//         gf.col(0)=gvvalue.row(i);
//         Glist[i]=gf;
//     }
//     // solve the first part: derivate(grad(F)).transpose()*grad(F).norm()*R*grad(F)
//     Eigen::VectorXd grad_norm = gvvalue.rowwise().norm(); // take 1/||gradV(i)|| as ith diagnal element
//     spMat diag_grad_norm_trip=grad_norm.replicate(3,1).asDiagonal();
//     spMat diag_Rs, diag_Grad;
//     spMat Jgrad;
//     diag_Rs.resize(3*vnbr,3*vnbr);
//     diag_Grad.resize(3*vnbr,vnbr);
//     Jgrad.resize(3*vnbr,vnbr);
//     assert(diag_grad_norm_trip.rows()==3*vnbr);
//     diag_Rs= dense_mat_list_as_sparse_diagnal(Rlist);
//     diag_Grad= dense_mat_list_as_sparse_diagnal(Glist);
//     for(int i=0;i<vnbr;i++){

//         Jgrad.block(3 * i, 0, 1, vnbr) = gradV[0].block(i, 0, 1, vnbr);
//         Jgrad.block(3 * i + 1, 0, 1, vnbr) = gradV[1].block(i, 0, 1, vnbr);
//         Jgrad.block(3 * i + 2, 0, 1, vnbr) = gradV[2].block(i, 0, 1, vnbr);
//     }
//     spMat JFirst=Jgrad.transpose()*diag_grad_norm_trip*diag_Rs*diag_Grad;

// }
// min()
void lsTools::initialize_and_smooth_level_set_by_laplacian()
{
    int vnbr = V.rows();
    int fnbr = F.rows();

    refids.clear();
    refids.push_back(F(assign_face_id, 0));
    refids.push_back(F(assign_face_id, 1));
    refids.push_back(F(assign_face_id, 2));
    if (fvalues.size() != vnbr)
    {
        fvalues.setZero(vnbr);
        for (int i = 0; i < 3; i++)
        {
            fvalues(F(assign_face_id, i)) = assign_value[i];
        }
    }
    get_gradient_hessian_values();

    // the matrices
    spMat H;
    Efunc B;
    // H*dx=B, H: vnbr*vnbr, B:vnbr.
    assemble_solver_laplacian_part(H, B);

    for (int i = 0; i < 3; i++)
    {
        H.coeffRef(F(assign_face_id, i), F(assign_face_id, i)) += weight_assign_face_value;
        B.coeffRef(F(assign_face_id, i)) += (assign_value[i] - fvalues(F(assign_face_id, i))) * weight_assign_face_value;
    }
    // std::cout << "finish assemble solver with assigned values" << std::endl;
    // now add mass matrix to make the matrix full rank
    assert(mass.rows() == vnbr);
    spMat mass_inverse = mass;
    for (int i = 0; i < vnbr; i++)
    {
        mass_inverse.coeffRef(i, i) = 1 / mass_inverse.coeffRef(i, i);
    }
    H += weight_mass * mass_inverse;

    assert(H.rows() == vnbr);
    assert(H.cols() == vnbr);
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver(H);
    assert(solver.info() == Eigen::Success);

    Eigen::VectorXd dx = solver.solve(B).eval();
    // std::cout << "step length " << dx.norm() << std::endl;
    // level_set_step_length = dx.norm();
    fvalues += dx;
    // for (int i = 0; i < 3; i++)
    // {
    //     fvalues(F(assign_face_id, i)) = assign_value[i];
    // }
}
void lsTools::initialize_and_optimize_strip_width()
{
    int vnbr = V.rows();
    int fnbr = F.rows();

    refids.clear();
    refids.push_back(F(assign_face_id, 0));
    refids.push_back(F(assign_face_id, 1));
    refids.push_back(F(assign_face_id, 2));
    if (fvalues.size() != vnbr)
    {
        fvalues.setZero(vnbr);
        for (int i = 0; i < 3; i++)
        {
            fvalues(F(assign_face_id, i)) = assign_value[i];
        }
    }
    get_gradient_hessian_values();


    // the matrices
    spMat H;
    Efunc B;
    assemble_solver_strip_width_part(H, B);
    H += weight_mass * mass;
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver(H);
    assert(solver.info() == Eigen::Success);

    Eigen::VectorXd dx = solver.solve(B).eval();
    // level_set_step_length = dx.norm();
    fvalues += dx;
    // next check if the energy is getting lower
    Eigen::VectorXd gnorms = gvvalue.rowwise().norm();
    double energy = 0;
    for (int i = 0; i < vnbr; i++)
    {
        gnorms(i) -= strip_width;
    }
    gnorms = gnorms; // check only the interior points
    for (int i = 0; i < vnbr; i++)
    {
        energy += mass.coeffRef(i, i) * gnorms(i) * gnorms(i);
    }
    std::cout << "energy:: " << energy << std::endl;
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
    double dis1 = abs(ver - start);
    double dis2 = abs(end - start);
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
    std::cout<<"test 1 "<<std::endl;
    Eigen::Vector3d ver_value0=trace_vers[0][0];
    Eigen::Vector3d ver_value1=trace_vers[lssize-1][0];

    double t1=get_t_of_segment(ver_value0,ver0,ver1);
    double t2=get_t_of_segment(ver_value1,ver2,ver3);
    Eigen::Vector2d para_value0=get_2d_ver_from_t(t1,ver2d0,ver2d1);
    Eigen::Vector2d para_value1=get_2d_ver_from_t(t1,ver2d2,ver2d3);
    std::cout<<"test 2 "<<std::endl;
    Eigen::Vector2d trans=(value1-value0)*(para_value1-para_value0)/(para_value1-para_value0).dot(para_value1-para_value0);
    std::cout<<"test 4 "<<std::endl;
    Eigen::VectorXd dupvalue=duplicate_valus(value0,V.rows());
    std::cout<<"test 5 "<<std::endl;
    Eigen::MatrixXd dupmatrix=duplicate_vector(para_value0,V.rows());
    std::cout<<"test 3 "<<std::endl;
    Eigen::MatrixXd tmp=paras-dupmatrix;// matrix of Vi-V0
    Eigen::VectorXd result=tmp*trans;
    result=result-dupvalue;
    fvalues=result;
}


void lsTools::initialize_level_set_by_tracing(const std::vector<int> &ids, const std::vector<double> values){
    // cylinder_open_example(5, 10, 50, 30);
    // exit(0);
    trace_vers.clear();
    trace_hehs.clear();
    assigned_trace_ls.clear();
    double target_angle = values[0];// 
    double start_angel = values[1];
    trace_start_angle_degree=start_angel;
    double threadshold_angel_degree = values[2]; // threadshold for checking mesh boundary corners
    int nbr_itv = ids[0]; // every nbr_itv boundary edges we shoot one curve
    id_dbg = ids[1];
    int which_segment=ids[2];
    std::vector<std::vector<CGMesh::HalfedgeHandle>> boundaries;
    split_mesh_boundary_by_corner_detection(lsmesh, V, threadshold_angel_degree,Boundary_Edges, boundaries);
    std::cout<<"get the boundary segments, how many: "<<boundaries.size()<<std::endl;
    if (nbr_itv < 1)
    {
        std::cout << "Please set up the parameter nbr_itv " << std::endl;
        return;
    }
    std::vector<CGMesh::HalfedgeHandle> boundary_segment=boundaries[which_segment];
    std::cout<<"the number of edges on this segment "<<boundary_segment.size()<<std::endl;
    OpenMesh::HalfedgeHandle init_edge = boundary_segment[0];
    
    OpenMesh::HalfedgeHandle checking_edge = init_edge;
    int beid=0;
    double lsvalue=0;
    int curvecount=0;
    while (1){
        ver_dbg.resize(0, 0);
        ver_dbg1.resize(0, 0);
        flag_dbg = false;
        E0_dbg.resize(0, 0);
        direction_dbg.resize(0, 0);
        pnorm_dbg.resize(0, 0);
        pnorm_list_dbg.clear();
        flag_dbg = true;
        checking_edge=boundary_segment[beid];
        std::vector<Eigen::Vector3d> curve;
        std::vector<CGMesh::HalfedgeHandle> handles;
        trace_single_pseudo_geodesic_curve(target_angle, checking_edge, 0.5, start_angel,
                                                                curve, handles);
        curvecount++;
        std::cout<<"one curver traced, size "<<curve.size()<<" the ith "<<curvecount<<std::endl;
        
        flag_dbg = false;
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
    if(eids.empty()||eids.size()==1){// if no value, or the value is on one vertex, no segment.
        return false;
    }
    assert(eids.size()!=3);
    
    double t = get_t_of_value(value,fvalues[eids[0][0]],fvalues[eids[0][1]]);
    E0 = get_3d_ver_from_t(t, V.row(eids[0][0]), V.row(eids[0][1]));
    t=get_t_of_value(value,fvalues[eids[1][0]],fvalues[eids[1][1]]);
    E1 = get_3d_ver_from_t(t, V.row(eids[1][0]), V.row(eids[1][1]));
    return true;
}
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
    assert(trace_vers.size()!=0); // must have at least one curve
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
void lsTools::print_info(const int vid){
    int i=-1;
    int ninner=IVids.size();
    if(ActInner.size()<=0){
        std::cout<<"please opt pseudo-geodesic energy"<<std::endl;
        return;
    }
    for (int j = 0; j < ninner; j++)
    {
        if(IVids[j]==vid){
            i=j;
            break;
        }
    }
    if(i==-1){
        std::cout<<"please check inner vers"<<std::endl;
        return;
    }
    if (ActInner[i] == false)
    { // this is a singularity
        std::cout<<"singularity"<<std::endl;
        return;
    }

    CGMesh::HalfedgeHandle heh1 = Vheh0[i];
    CGMesh::HalfedgeHandle heh2 = Vheh1[i];
    int fid1, fid2;
    fid1 = lsmesh.face_handle(heh1).idx();
    fid2 = lsmesh.face_handle(heh2).idx();
    assert(fid1 != fid2);
    Eigen::Vector3d norm = norm_v.row(vid);
    Eigen::Vector3d norm1 = norm_f.row(fid1);
    Eigen::Vector3d norm2 = norm_f.row(fid2);
    Eigen::Vector3d g1 = gfvalue.row(fid1);
    Eigen::Vector3d g2 = gfvalue.row(fid2);
    std::cout<<"gradients "<<g1.transpose()<<"\n"<<g2.transpose()<<std::endl;
    // rotate gradients to get iso-curve directions
    g1 = norm1.cross(g1);
    g2 = norm2.cross(g2);
    std::cout<<"iso-line \n"<<g1.transpose()<<"\n"<<g2.transpose()<<std::endl;
    Eigen::Vector3d g1xg2 = g1.cross(g2);
    
    std::cout<<"g1xg2 "<<g1xg2.transpose()<<", norm, "<<g1xg2.norm()<<"\n";
    // std::cout<<"energy, "<<g1xg2.dot(norm) / g1xg2.norm() - cos_angle
    Eigen::Vector3d binormal=LsOrient[i]*g1xg2.normalized();
    std::cout<<"binormal "<<binormal.transpose()<<"\n";
    std::cout<<"vBinormal "<<vBinormal.row(i)<<"\n";
    double cos_real=norm.dot(binormal);
    std::cout<<"Vnorm "<<norm.transpose()<<"\ncos_real "<<cos_real<<std::endl;
    std::cout<<"energy value "<<VPEvalue.coeffRef(i)<<std::endl;
    
    if(enable_functional_angles){
        std::cout<<"target angle "<<pseudo_geodesic_angles_per_ver[vid]<<std::endl;
        double angle_radian = pseudo_geodesic_angles_per_ver[vid] * LSC_PI / 180.; // the angle in radian
        double cos_angle = cos(angle_radian);
        std::cout<<"energy recomput "<<LsOrient[i] * g1xg2.dot(norm) / g1xg2.norm() - cos_angle<<std::endl;

    }
    else{
        std::cout<<"target angle "<<pseudo_geodesic_target_angle_degree<<std::endl;
    }
    // std::cout<<"pew\n"<<PeWeight.transpose()<<std::endl;
    // std::cout<<"peV\n"<<VPEvalue<<std::endl;
    std::cout<<"n(peV) "<<VPEvalue.norm()<<std::endl;
    
}