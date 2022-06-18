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
        Eigen::Vector2d p1 = paras.row(id0), p2 = paras.row(id1), p3 = paras.row(id2);
        areaPF(i) = fabs(0.5 * (p1[0] * p2[1] - p2[0] * p1[1]) + (p2[0] * p3[1] - p3[0] * p2[1]) + (p3[0] * p1[1] - p1[0] * p3[1]));
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

void lsTools::get_rotated_parameter_edges()
{
    Eigen::Matrix2d rotate2d; // 2d rotation matrix
    rotate2d << 0, -1,
        1, 0;

    Erotate2d.resize(F.rows());

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
            Erotate2d[fid][ecounter] = Eigen::Vector2d(rotate2d * Eigen::Vector2d(paras.row(vid2) - paras.row(vid1)));
            assert(fid < Erotate2d.size());
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
        assert(ecounter == 3);
    }
}
// the input should be the values on each vertex, nx3, (x, y, z) or (x_u, y_u, z_u)
// the output is the derivate on each facet, for u and v separately.
void lsTools::surface_derivate_v2f(const Eigen::MatrixXd &vvalues, std::array<Eigen::MatrixXd, 2> &DonF)
{
    Eigen::Matrix<double, 3, 3> pos; // the 3 vertices of the face
    Eigen::Matrix<double, 3, 2> par; // the parameters of the 3 vertices
    assert(vvalues.rows() == V.rows());
    int vsize = V.rows();
    int fsize = F.rows();
    DonF[0].resize(fsize, 3); // partial derivate to u: x_u, y_u, z_u
    DonF[1].resize(fsize, 3); // partial derivate to v

    for (int fid = 0; fid < F.rows(); fid++)
    {
        int id0 = F(fid, 0);
        int id1 = F(fid, 1);
        int id2 = F(fid, 2);
        Eigen::Vector2d rot20 = Erotate2d[fid][0]; // v0-v2
        Eigen::Vector2d rot01 = Erotate2d[fid][1]; // v1-v0
        Eigen::Vector2d rot12 = Erotate2d[fid][2]; // v2-v1
        par.row(0) = rot12;
        par.row(1) = rot20;
        par.row(2) = rot01;
        pos.col(0) = vvalues.row(id0);
        pos.col(1) = vvalues.row(id1);
        pos.col(2) = vvalues.row(id2);

        double areap = areaPF(fid);
        assert(areap > 0);
        Eigen::MatrixXd res = pos * par / (2 * areap); // a 3x2 matrix, [xu,xv;yu,yv;zu,zv]
        DonF[0].row(fid) = res.col(0);
        DonF[1].row(fid) = res.col(1);
    }
}

void lsTools::surface_derivate_f2v(const std::array<Eigen::MatrixXd, 2> &DonF, std::array<Eigen::MatrixXd, 2> &DonV)
{
    int vsize = V.rows();
    int fsize = F.rows();
    assert(DonF[0].rows() == fsize);
    DonV[0] = Eigen::MatrixXd::Zero(vsize, 3);
    DonV[1] = Eigen::MatrixXd::Zero(vsize, 3);
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
            DonV[0].row(i) += angle * DonF[0].row(fid);
            DonV[1].row(i) += angle * DonF[1].row(fid);
        }
        assert(anglesum>0);
        DonV[0].row(i) = DonV[0].row(i) / anglesum;
        DonV[1].row(i) = DonV[1].row(i) / anglesum;
    }
}
void lsTools::surface_derivate_easy_interface(const Eigen::MatrixXd &vvalues, std::array<Eigen::MatrixXd, 2> &DonV)
{
    assert(vvalues.rows() == V.rows() && vvalues.cols() == 3);
    std::array<Eigen::MatrixXd, 2> DonF;
    surface_derivate_v2f(vvalues, DonF);
    surface_derivate_f2v(DonF, DonV);
}
void lsTools::get_surface_derivate()
{
    surface_derivate_easy_interface(V, Deriv1); // get (x_u, y_u, z_u) and (x_v, y_v, z_v)
    std::array<Eigen::MatrixXd, 2> tmp;
    surface_derivate_easy_interface(Deriv1[0], tmp); // get (x_uu, y_uu, z_uu) and (x_uv, y_uv, z_uv)
    Deriv2[0] = tmp[0];
    Deriv2[1] = tmp[1];
    surface_derivate_easy_interface(Deriv1[1], tmp); // get (x_vu, y_vu, z_vu) and (x_vv, y_vv, z_vv)
    Deriv2[2] = tmp[0];
    Deriv2[3] = tmp[1];
}
void lsTools::get_surface_II_each_ver()
{
    int vsize = V.rows();
    II_L.resize(vsize);
    II_M.resize(vsize);
    II_N.resize(vsize);
    int dbgcount = 0;
    for (int i = 0; i < vsize; i++)
    {
        II_L[i] = Deriv2[0].row(i).dot(norm_v.row(i));     // r_uu*n
        II_N[i] = Deriv2[3].row(i).dot(norm_v.row(i));     // r_vv*n
        double r_vu = Deriv2[2].row(i).dot(norm_v.row(i)); // r_vu*n
        double r_uv = Deriv2[1].row(i).dot(norm_v.row(i));
        II_M[i] = (r_uv + r_vu) / 2; // r_uv*n
        // if(fabs(II_M[i]-r_vu)>SCALAR_ZERO){
        //     dbgcount++;
        //     std::cout<<"calculation of II is not accurate:\n"<<II_M[i]<<", "<<r_vu<<", diff: "<<fabs(II_M[i]-r_vu)<<std::endl
        //     <<"nbr, "<<dbgcount<<std::endl;

        // }
    }
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
        Eigen::Matrix2d LMN;
        LMN << II_L[vid], -II_M[vid],
            -II_M[vid], II_N[vid];
        Eigen::Matrix<double, 3, 2> rvru;
        rvru.col(0)=Deriv1[1].row(vid);// x_v, y_v, z_v
        rvru.col(1)=Deriv1[0].row(vid);// x_u, y_u, z_u
        Eigen::Vector3d gradient=gvvalue.row(vid);
        Eigen::Matrix3d hessian=hfvalue[vid];
        Eigen::Matrix<double, 3,2> mrurv; // [-ru, rv]
        mrurv.col(0)=-1*Deriv1[0].row(vid);
        mrurv.col(1)=Deriv1[1].row(vid);
        double normgradf=gradient.norm();
        Eigen::MatrixXd inner_left=normgradf*rvru*LMN*rvru.transpose();
        Eigen::MatrixXd inner_right=rvru*mrurv.transpose()*hessian*mrurv*rvru.transpose();
        double left=gradient.transpose()*inner_left*gradient;
        double right=gradient.transpose()*inner_right*gradient;
        assert(right!=0);
        double b=left/right;
        std::cout<<"the current row nbr "<<rowid<<" with pid "<<i<<"\n**b is "<<b<<"\ngradient, "<<gradient.transpose()<<"\n\n";
        Eigen::Vector3d rotgrad=RotateV[vid]*gradient;
        Eigen::Vector3d s1=hessian*rotgrad;
        double s2=rotgrad.dot(s1);
        double kg=(s2)/(normgradf*normgradf*normgradf);
        std::cout<<"Kg is "<<kg<<std::endl;
    }
    // double zav=ztotal/rnbr;
    // double zdiffmax=0;
    // for(int i=0;i<rnbr;i++)
    // {
    //     int vid = rowid * rnbr + i;
    //     double zdiff=fabs(zav-V(vid,2));
    //     if(zdiffmax<zdiff){
    //         zdiffmax=zdiff;
    //     }
    // }
    // double max_edge_error=0;
    // double max_rot_err=0;
    // int maxrid=-1, maxeid=-1;
    // for(int i=0;i<F.rows();i++){
    //     int v0id=F(i,0);
    //     int v1id=F(i,1);
    //     int v2id=F(i,2);
    //     double error0=Erotate[i][0].dot(Eigen::Vector3d(V.row(v0id)-V.row(v2id)));
    //     double error1=Erotate[i][1].dot(Eigen::Vector3d(V.row(v1id)-V.row(v0id)));
    //     double error2=Erotate[i][2].dot(Eigen::Vector3d(V.row(v2id)-V.row(v1id)));
    //     // std::cout<<"edge dot, fid "<<i<<", "<<error0<<" "<<error1
    //     // <<" "<<error2<<"\n";
    //     if(max_edge_error<fabs(error0)){
    //         max_edge_error=fabs(error0);
    //         maxeid=i;
    //     }
    //     if(max_edge_error<fabs(error0)){
    //         max_edge_error=fabs(error0);
    //         maxeid=i;
    //     }
    //     if(max_edge_error<fabs(error1)){
    //         max_edge_error=fabs(error1);
    //         maxeid=i;
    //     }
    //     if(max_edge_error<fabs(error2)){
    //         max_edge_error=fabs(error2);
    //         maxeid=i;
    //     }
    //     double rdt=Rotate[i].determinant();
    //     rdt=fabs(rdt-1);
    //     if(max_rot_err<rdt){
    //         max_rot_err=rdt;
    //         maxrid=i;
    //     }
    // }
    
    // std::cout<<"max edge dot error "<<max_edge_error<<" in fid "<<maxeid<<std::endl;
    // std::cout<<"max rotate deter error "<<max_rot_err<<" in fid "<<maxrid<<std::endl;
    // std::cout<<"max z value diff "<<zdiffmax<<std::endl;

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
void lsTools::show_face_grad_max_angle(Eigen::MatrixXd& points){
    // double max_diff=0;
    // int diff_id=0;
    // int fsize=F.rows();
    // for(int i=0;i<fsize;i++){
    //     Eigen::Vector3d v0=Eigen::Vector3d(V.row(F(i,0)));
    //     Eigen::Vector3d v1=Eigen::Vector3d(V.row(F(i,1)));
    //     Eigen::Vector3d v2=Eigen::Vector3d(V.row(F(i,2)));
    //     Eigen::Vector3d e01=v0-v1;
    //     Eigen::Vector3d e12=v1-v2;
    //     Eigen::Vector3d e20=v2-v0;
    //     Eigen::Vector3d grad=Eigen::Vector3d(gfvalue.row(i));
    //     double diff01=fabs(e01.normalized().dot(grad.normalized()));
    //     double diff12=fabs(e12.normalized().dot(grad.normalized()));
    //     double diff20=fabs(e20.normalized().dot(grad.normalized()));
    //     double diff=std::min(std::min(diff01,diff12),diff20);
    //     if(diff>max_diff){
    //         max_diff=diff;
    //         diff_id=i;
    //     }
    // }
    // std::cout<<"max grad-edge dot diff is "<<max_diff<<", the angle is "<<acos(max_diff)*180/3.1415926<<", fid "<<diff_id<<std::endl;
    // points.resize(3,3);
    // points.row(0)=V.row(F(diff_id,0));
    // points.row(1)=V.row(F(diff_id,1));
    // points.row(2)=V.row(F(diff_id,2));
    // std::cout<<"the three vertices\n "<<points.row(0)<<"\n"<<points.row(1)<<"\n"<<points.row(2)<<"\n";
    // std::cout<<"the fvalues "<<fvalues(F(diff_id,0))<<", "<<fvalues(F(diff_id,1))<<", "<<fvalues(F(diff_id,2))<<"\n";
    // std::cout<<"rotation matrix\n "<<Rotate[diff_id]<<"\n";
    // std::cout<<"normal "<<norm_f.row(diff_id)<<"\n";
    // std::cout<<"gradient "<<gfvalue.row(diff_id)<<std::endl;
    // std::cout<<"gradient*norm, "<<Eigen::Vector3d(gfvalue.row(diff_id)).dot(Eigen::Vector3d(norm_f.row(diff_id)))<<"\n"<<std::endl;
    
    // Eigen::Vector3d re02 = Erotate[diff_id][0];
    // Eigen::Vector3d re10 = Erotate[diff_id][1];
    // Eigen::Vector3d re21 = Erotate[diff_id][2];
    // Eigen::Vector3d norm = norm_f.row(diff_id);
    // std::cout<<"rotated edges * norm, "<<re02.dot(norm)<<", "<<re10.dot(norm)<<", "<<re21.dot(norm)<<std::endl;
    // Eigen::Vector3d grad=Eigen::Vector3d(gfvalue.row(diff_id));
    // std::cout<<"sum of rotated "<<re02+re10+re21<<std::endl;
    // Eigen::Vector3d v0 = Eigen::Vector3d(V.row(F(diff_id, 0)));
    // Eigen::Vector3d v1 = Eigen::Vector3d(V.row(F(diff_id, 1)));
    // Eigen::Vector3d v2 = Eigen::Vector3d(V.row(F(diff_id, 2)));
    // Eigen::Vector3d e01 = v0 - v1;
    // Eigen::Vector3d e12 = v1 - v2;
    // Eigen::Vector3d e20 = v2 - v0;
    // std::cout << "the three edges dot gradient " <<e01.dot(grad)<<", "<<e12.dot(grad)<<", "<<e20.dot(grad)<<std::endl;
    // Eigen::Vector3d fake_grad=fvalues(F(diff_id,0))*re21+fvalues(F(diff_id,1))*re02+fvalues(F(diff_id,2))*re10;
    // double area=areaF(diff_id);
    // fake_grad/=2*area;
    // std::cout<<"recalculated grad, "<<fake_grad<<std::endl;
    // std::cout<<"recalculated gradient dot edges "<<e01.dot(fake_grad)<<", "<<e12.dot(fake_grad)<<", "<<e20.dot(fake_grad)<<std::endl;
}