#include <lsc/basic.h>
#include <igl/lscm.h>
#include <igl/harmonic.h>
#include <igl/boundary_loop.h>
#include <igl/grad.h>
#include <igl/hessian.h>

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
void lsTools::get_bnd_and_bnd_one_ring()
{
    int vnbr=V.rows();
    int fnbr=F.rows();
    std::vector<Trip> triplets;
    triplets.reserve(vnbr);
    Eigen::VectorXd flags=Eigen::VectorXd::Ones(vnbr);// first by default set all vertices as non-boundary vertices
    for (CGMesh::VertexIter v_it = lsmesh.vertices_begin(); v_it != (lsmesh.vertices_end()); ++v_it)
    {
        if(lsmesh.is_boundary(v_it.handle())){// find all the boundary vertices
            int vid=v_it.handle().idx();
            flags(vid)=0;
        }
        else// find all the vertices whose one-ring vertices contains boundary.
        {
            for (CGMesh::VertexVertexIter vv_it = lsmesh.vv_begin(v_it.handle()); vv_it != lsmesh.vv_end(v_it.handle()); ++vv_it)
            {
                if (lsmesh.is_boundary(vv_it))
                {
                    int vid=v_it.handle().idx();
                    flags(vid)=0;
                    break;
                }
            }
        }
    }
    ORB=flags.asDiagonal();
    // for(int i=0;i<vnbr;i++){
    //     if(flags(i)>0){
    //         triplets.push_back(Trip(i,i,1.));
    //     }
    // }

    // ORB.setFromTriplets(triplets.begin(),triplets.end());
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
    for(int i=0;i<3;i++){
        for(int j=0;j<3;j++){
            HessianV[i][j]=gradV[j]*gradV[i];
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
        std::cout<<"angle calculated is "<<atan(b) * 180 / 3.1415926<<std::endl;
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


void lsTools::get_all_the_derivate_matrices()
{

    int vsize = V.rows();
    igl::cotmatrix(V, F, Dlps);
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
    Eigen::MatrixXd grad_norm_mat(vsize,1);
    grad_norm_mat.col(0)=grad_norm;
    Eigen::DiagonalMatrix<double, Eigen::Dynamic, Eigen::Dynamic> grad_norm_diag=grad_norm_mat.asDiagonal();
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

// assemble matrices for \sum{area(i)*(||gradient(i)||-strip_width)^2}
void lsTools::assemble_solver_strip_width_part(spMat &H, Efunc &B)
{
    int vsize = V.rows();
    Eigen::VectorXd f = gvvalue.rowwise().norm(); // f=||gradient{i}||-strip_width
    for (int i = 0; i < vsize; i++)
    {
        f(i) -= strip_width;
    }
    spMat J=Dgrad_norm;
    // now only consider about the interior vertices but not the boundaries and one ring from them
    J =  J;
    f =  f;
    // end of interior part
    spMat JTJ = J.transpose() * mass * J;
    Eigen::VectorXd mJTF_dense = -J.transpose() * mass * f;
    H = JTJ;
    B = dense_vec_to_sparse_vec(mJTF_dense);
}

// we summarize the problem as:
// grad(F).transpose()*grad(F).norm()*R*grad(F)-grad(F).transpose()*D.transpose()*hess(F)*D*grad(F).
void lsTools::assemble_solver_pseudo_geodesic_part(spMat &H, Efunc& B){
    
    int vnbr=V.rows();
    int fnbr=F.rows();
    std::vector<Eigen::MatrixXd> Rlist(vnbr);
    std::vector<Eigen::MatrixXd> Dlist(vnbr);
    std::vector<Eigen::MatrixXd> Glist(vnbr);
    for (int i = 0; i < vnbr; i++)
    {
        Eigen::MatrixXd rvru(3, 2);
        rvru.col(0) = Deriv1[1].row(i);
        rvru.col(1) = Deriv1[0].row(i);
        Eigen::Matrix<double, 3, 2> rurv;
        rurv.col(0) = Deriv1[0].row(vid);
        rurv.col(1) = Deriv1[1].row(vid);
        Eigen::Matrix<double, 3, 2> mrvru; // [-rv, ru]
        mrvru.col(0) = -1 * Deriv1[1].row(vid);
        mrvru.col(1) = Deriv1[0].row(vid);
        Eigen::Matrix2d LmMN;
        LmMN << II_L[i], -II_M[i],
            -II_M[i], II_N[i];
        Eigen::MatrixXd R=pseudo_geodesic_ratio*rvru*LmMN*rvru.transpose();
        Rlist[i]=R;
        Eigen::MatrixXd D=rurv*mrvru.transpose();
        Dlist[i]=D;
        Eigen::MatirxXd gf(3,1);
        gf.col(0)=gvvalue.row(i);
        Glist[i]=gf;
    }
    // solve the first part: derivate(grad(F)).transpose()*grad(F).norm()*R*grad(F)
    Eigen::VectorXd grad_norm = gvvalue.rowwise().norm(); // take 1/||gradV(i)|| as ith diagnal element
    spMat diag_grad_norm_trip=grad_norm.replicate(3,1).asDiagonal();
    spMat diag_Rs, diag_Grad;
    spMat Jgrad;
    diag_Rs.resize(3*vnbr,3*vnbr);
    diag_Rs.array()*=0;
    diag_Grad.resize(3*vnbr,vnbr);
    Jgrad.resize(3*vnbr,vnbr);
    assert(diag_grad_norm_trip.rows()==3*vnbr);
    for(int i=0;i<vnbr;i++){
        diag_Rs.block(i*3,i*3,3,3)=Rlist[i].sparseView();
        diag_Grad.block(i*3,i,3,1)=Glist[i].sparseView();
        Jgrad.block(3*i,0,1,vnbr)=gradV[0].block(i,0,1,vnbr);
        Jgrad.block(3*i+1,0,1,vnbr)=
        Jgrad.block(3*i+2,0,1,vnbr)=
    }
    

}
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

    // std::cout << "finished calculate gradient and hessian" << std::endl;
    get_all_the_derivate_matrices();
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
    level_set_step_length = dx.norm();
    fvalues += dx;
    // for (int i = 0; i < 3; i++)
    // {
    //     fvalues(F(assign_face_id, i)) = assign_value[i];
    // }
}
void lsTools::initialize_and_optimize_strip_width(){
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

    // std::cout << "finished calculate gradient and hessian" << std::endl;
    get_all_the_derivate_matrices();
    // the matrices
    spMat H;
    Efunc B;
    assemble_solver_strip_width_part(H,B);
    H += weight_mass * mass;
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver(H);
    assert(solver.info() == Eigen::Success);

    Eigen::VectorXd dx = solver.solve(B).eval();
    level_set_step_length = dx.norm();
    fvalues += dx;
    // next check if the energy is getting lower
    Eigen::VectorXd gnorms =  gvvalue.rowwise().norm();
    double energy = 0;
    for (int i = 0; i < vnbr; i++)
    {
        gnorms(i) -= strip_width;
    }
    gnorms=gnorms;// check only the interior points
    for (int i = 0; i < vnbr; i++)
    {
        energy += mass.coeffRef(i, i) * gnorms(i) * gnorms(i);
    }
    std::cout<<"energy:: "<<energy<<std::endl;
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