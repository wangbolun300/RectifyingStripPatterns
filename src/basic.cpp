#include <lsc/basic.h>
#include <igl/lscm.h>
#include <igl/harmonic.h>
#include <igl/boundary_loop.h>

lsTools::lsTools(CGMesh &mesh)
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
        double le0=(p2-p1).norm();
        double le1=(p0-p2).norm();
        double le2=(p0-p1).norm();
        double lp=(le0+le1+le2)/2;
        areaPF(i) = sqrt(lp*(lp-le0)*(lp-le1)*(lp-le2));
        
    }
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
// Vectorlf is a class wich restore the value of each vertex.
// it can be a Vectorlf, or a Eigen::Vector of double values
// the size of Vectorlf is already assigned .
// the 3 dimentions of output is for x, y and z
void lsTools::gradient_v2f(Vectorlf &input, std::array<Vectorlf, 3> &output)
{
    int vsize = V.rows();
    int fsize = F.rows();
    output[0].resize(fsize, vsize); // each face has a output value, which is the combination of function input of each vertex
    output[1].resize(fsize, vsize); // each face has a output value, which is the combination of function input of each vertex
    output[2].resize(fsize, vsize); // each face has a output value, which is the combination of function input of each vertex

    assert(input.size() == vsize); // each vertex has a input value
    for (int fid = 0; fid < F.rows(); fid++)
    {
        int id0 = F(fid, 0);
        int id1 = F(fid, 1);
        int id2 = F(fid, 2);
        Eigen::Vector3d rot20 = Erotate[fid][0]; // v0-v2
        Eigen::Vector3d rot01 = Erotate[fid][1]; // v1-v0
        Eigen::Vector3d rot12 = Erotate[fid][2]; // v2-v1
        double area = areaF(fid);
        assert(area > 0);
        // assert(dimin == 1 || dimin == 3);

        auto value0 = input(id0); // value on the vertex
        auto value1 = input(id1);
        auto value2 = input(id2);
        // if(dbg01){
        //     std::cout<<"the values on v0, v1 and v2, "<<value0<<", "<<value1<<", "<<value2<<std::endl;
        //     std::cout<<"the rotated edges \n"<<rot20.transpose()<<"\n"<<rot01.transpose()<<"\n"<<rot12.transpose()<<"\n";
        //     dbg01=false;
        // }
        for (int itr = 0; itr < 3; itr++)
        {

            output[itr](fid) = (value0 * rot12[itr] + value1 * rot20[itr] + value2 * rot01[itr]) / (2 * area);
        }
    }
}

void lsTools::gradient_f2v(std::array<Vectorlf, 3> &input, std::array<Vectorlf, 3> &output)
{
    int vsize = V.rows();
    int fsize = F.rows();
    assert(input[0].size() == fsize);
    output[0].resize(vsize, vsize);
    output[1].resize(vsize, vsize);
    output[2].resize(vsize, vsize);
    for (int i = 0; i < V.rows(); i++)
    {
        CGMesh::VertexHandle vh = lsmesh.vertex_handle(i); // for each vertex, iterate all the faces

        double anglesum = 0;
        for (CGMesh::VertexFaceIter vf_it = lsmesh.vf_begin(vh); vf_it != lsmesh.vf_end(vh); ++vf_it)
        {
            assert(vh.idx() == i);
            int fid = vf_it.handle().idx();
            int vvid=-1;
            for(int vv=0;vv<3;vv++){
                if(i==F(fid,vv)){
                    vvid=vv;
                }
            }
            assert(vvid>-1);
            double angle = angF(fid,vvid);
            assert(angle>0);
            anglesum += angle;
            assert(i < output[0].size());
            output[0](i) += angle * input[0](fid);
            output[1](i) += angle * input[1](fid);
            output[2](i) += angle * input[2](fid);
        }
        output[0](i) = output[0](i) / anglesum;
        output[1](i) = output[1](i) / anglesum;
        output[2](i) = output[2](i) / anglesum;
    }
}
void lsTools::initialize_function_coefficient(Vectorlf &input)
{
    int vsize = V.rows();
    int fsize = F.rows();
    // for each vertex there is a function value
    // #input# is in the form of an identity matrix
    input.resize(vsize, vsize);
    for (int i = 0; i < vsize; i++)
    {
        input(i).coeffRef(i) = 1;
    }
}
void lsTools::gradient_easy_interface(Vectorlf &func, std::array<Vectorlf, 3> &gonf, std::array<Vectorlf, 3> &gonv)
{
    gradient_v2f(func, gonf); // calculate the gradients on faces;
    gradient_f2v(gonf, gonv); // calculate the gradients on vertices by area-weighted averaging of face gradients
}

void lsTools::get_function_gradient_vertex()
{
    Vectorlf iden; // this is the initial function coefficients, in the form of an identity matrix
    initialize_function_coefficient(iden);
    gradient_easy_interface(iden, gradVF, gradV);
}

void lsTools::get_function_hessian_vertex()
{
    std::array<Vectorlf, 3> temp_f; // temprary value calculated on the faces
    Vectorlf fx = gradV[0];         // partial(f)(x)
    Vectorlf fy = gradV[1];         // partial(f)(y)
    Vectorlf fz = gradV[2];         // partial(f)(z)
    temp_f[0].clear();
    temp_f[1].clear();
    temp_f[2].clear();
    gradient_easy_interface(fx, temp_f, HessianV[0]); // fxx, fxy, fxz
    temp_f[0].clear();
    temp_f[1].clear();
    temp_f[2].clear();
    gradient_easy_interface(fy, temp_f, HessianV[1]); // fyx, fyy, fyz
    temp_f[0].clear();
    temp_f[1].clear();
    temp_f[2].clear();
    gradient_easy_interface(fz, temp_f, HessianV[2]); // fzx, fzy, fzz
}

void lsTools::get_lf_value(const Vectorlf &coff, Eigen::VectorXd &res)
{
    int length = coff.size();
    
    res.resize(length);
    for (int i = 0; i < length; i++)
    {
        res(i) = coff(i).dot(fvalues);
    }
}
void lsTools::get_gradient_hessian_values()
{
    int vsize = V.rows();
    int fsize = F.rows();
    gvvalue.resize(vsize, 3);
    hfvalue.resize(vsize);

    Eigen::VectorXd gx, gy, gz;
    get_lf_value(gradV[0], gx);
    get_lf_value(gradV[1], gy);
    get_lf_value(gradV[2], gz);
    gvvalue.col(0) = gx;
    gvvalue.col(1) = gy;
    gvvalue.col(2) = gz;
    gfvalue.resize(fsize,3);
    get_lf_value(gradVF[0],gx);
    get_lf_value(gradVF[1],gy);
    get_lf_value(gradVF[2],gz);
    gfvalue.col(0) = gx;
    gfvalue.col(1) = gy;
    gfvalue.col(2) = gz;
    Eigen::VectorXd hxx, hxy, hxz, hyx, hyy, hyz, hzx, hzy, hzz;
    get_lf_value(HessianV[0][0], hxx);
    get_lf_value(HessianV[0][1], hxy);
    get_lf_value(HessianV[0][2], hxz);

    get_lf_value(HessianV[1][0], hyx);
    get_lf_value(HessianV[1][1], hyy);
    get_lf_value(HessianV[1][2], hyz);

    get_lf_value(HessianV[2][0], hzx);
    get_lf_value(HessianV[2][1], hzy);
    get_lf_value(HessianV[2][2], hzz);

    for (int i = 0; i < vsize; i++)
    {
        hfvalue[i] << hxx[i], hyx[i], hzx[i],
            hxy[i], hyy[i], hzy[i],
            hxz[i], hyz[i], hzz[i];
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
            double upara = 0.5*3.1415926-theta + i * titv;
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
            int id0 = nt * j + i;
            int id1 = nt * (j + 1) + i;
            int id2 = nt * (j + 1) + i + 1;
            int id3 = nt * j + i + 1;
            faces.row(fline) = Eigen::Vector3i(id0, id1, id2);
            faces.row(fline + 1) = Eigen::Vector3i(id0, id2, id3);
            fline += 2;
        }
    }
    std::string path("/Users/wangb0d/bolun/D/vs/levelset/level-set-curves/data/");
    igl::write_triangle_mesh(path+"sphere_"+std::to_string(radius)+"_"+std::to_string(theta)+"_"+
    std::to_string(phi)+"_"+std::to_string(nt)+"_"+std::to_string(np)+".obj", ver, faces);
    std::cout<<"sphere mesh file saved "<<std::endl;
}

void lsTools::make_sphere_ls_example(int rowid){
    int vsize=V.rows();
    int fsize=F.rows();
    int rnbr=20;// # vertices in each row
    fvalues.resize(vsize);
    
    for(int i=0;i<vsize;i++){
        fvalues(i)=V(i,2);// make function value = z
    }
    get_gradient_hessian_values();
    refids.clear();

    for(int i=0;i<rnbr;i++)
    {
        int vid = rowid * rnbr + i;

        refids.push_back(vid);
        Eigen::Matrix2d LmMN;
        LmMN << II_L[vid], -II_M[vid],
            -II_M[vid], II_N[vid];
        Eigen::Matrix<double, 3, 2> rurv, rvru;
        rurv.col(0)=Deriv1[0].row(vid);
        rurv.col(1)=Deriv1[1].row(vid);
        rvru.col(0)=Deriv1[1].row(vid);
        rvru.col(1)=Deriv1[0].row(vid);
        Eigen::Vector3d gradient=gvvalue.row(vid);
        Eigen::Matrix3d hessian=hfvalue[vid];
        Eigen::Matrix<double, 3,2> mrvru; // [-rv, ru]
        mrvru.col(0)=-1*Deriv1[1].row(vid);
        mrvru.col(1)=Deriv1[0].row(vid);
        double normgradf=gradient.norm();
        Eigen::MatrixXd inner_left=normgradf*rvru*LmMN*rvru.transpose();
        Eigen::MatrixXd inner_right=mrvru*rurv.transpose()*hessian*rurv*mrvru.transpose();
        double left=gradient.transpose()*inner_left*gradient;
        double right=gradient.transpose()*inner_right*gradient;
        assert(right!=0);
        double b=right/left;
        std::cout<<"vid "<<vid<<"\n";
        std::cout<<"the current row nbr "<<rowid<<" with pid "<<i<<"\n**b is "<<b<<"\ngradient, "<<gradient.transpose()<<"\n";
        Eigen::Vector3d rotgrad=RotateV[vid]*gradient;
        Eigen::Vector3d s1=hessian*rotgrad;
        double s2=rotgrad.dot(s1);
        double kg=-(s2)/(normgradf*normgradf*normgradf);
        std::cout<<"Kg is "<<kg<<std::endl;
        std::cout<<"Kn (estimate) is "<<kg/b<<"\n";
        Eigen::Vector3d ru=Deriv1[0].row(vid);
        Eigen::Vector3d rv=Deriv1[1].row(vid);
        double gfrv=gradient.dot(rv);
        double gfru=gradient.dot(ru);
        // Eigen::Matrix2d LMN;
        // LMN << II_L[vid], II_M[vid],
        //     II_M[vid], II_N[vid];
        Eigen::Vector2d IIl = Eigen::Vector2d(gfrv, gfru);
        double up=IIl.transpose()*LmMN*IIl;
        double cosuv=ru.dot(rv)/(ru.norm()*rv.norm());
        double sin2=1-cosuv*cosuv;assert(sin2>=0);
        double down=ru.norm()*ru.norm()*rv.norm()*rv.norm()*normgradf*normgradf*sin2;
        double kn=up/down;
        std::cout<<"Kn (calculate) is "<<up/down<<std::endl;
        std::cout<<"b recalculate is "<<kg/kn<<std::endl;
  
        Eigen::Vector3d q=-(gradient.dot(rv))*ru+gradient.dot(ru)*rv;
        Eigen::Vector3d rf_test=q.normalized()*gradient.norm();
        //std::cout<<"gradf size is "<<gradient.norm()<<"\ngrad test angle diff "<<rf_test.dot(gradient)<<std::endl;
        std::cout<<"qnorm calculated and norm "<<gradient.norm()*ru.norm()*rv.norm()*sqrt(sin2)<<" "<<q.norm()<<std::endl;
        // std::cout<<"__\nthe"
        std::cout<<"\n";
    }

    std::cout<<"__________________________________________________"<<std::endl;
}
void lsTools::show_level_set(Eigen::VectorXd &val){
   val= fvalues;
}

void lsTools::show_gradients(Eigen::MatrixXd& E0, Eigen::MatrixXd &E1, double ratio){
    
    E0=V+gvvalue*ratio;
    E1=V-gvvalue*ratio;
    assert(V.rows()==gvvalue.rows());
}
void lsTools::show_current_reference_points(Eigen::MatrixXd& pts){
    int size= refids.size();
    pts.resize(size,3);
    for(int i=0;i<size;i++){
        pts.row(i)=V.row(refids[i]);
    }
}
void lsTools::show_face_gradients(Eigen::MatrixXd& E0, Eigen::MatrixXd &E1, double ratio){
    Eigen::MatrixXd fcent(F.rows(),3);
    for(int i=0;i<F.rows();i++){
        fcent.row(i)=(V.row(F(i,0))+V.row(F(i,1))+V.row(F(i,2)))/3;
    }
    Eigen::MatrixXd dirc(F.rows(),3);
    Eigen::VectorXd gx, gy, gz;
    get_lf_value(gradVF[0], gx);
    get_lf_value(gradVF[1], gy);
    get_lf_value(gradVF[2], gz);
    dirc.col(0) = gx;
    dirc.col(1) = gy;
    dirc.col(2) = gz;


    E0=fcent+dirc*ratio;
    E1=fcent-dirc*ratio;
    assert(fcent.rows()==dirc.rows());
}
void lsTools::show_1_order_derivate(Eigen::MatrixXd& E0, Eigen::MatrixXd &E1,Eigen::MatrixXd &E2,Eigen::MatrixXd &E3, double ratio){
    
    E0=V+Deriv1[0]*ratio;
    E1=V-Deriv1[0]*ratio;
    E2=V+Deriv1[1]*ratio;
    E3=V-Deriv1[1]*ratio;
}
void lsTools::show_vertex_normal(Eigen::MatrixXd& E0, Eigen::MatrixXd& E1, double ratio){
    E0.resize(V.rows(),3);
    E1=E0;
    E0=V;
    E1=E0+norm_v*ratio;
}
void lsTools::show_face_grad_max_angle(Eigen::MatrixXd& points){
    
}