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

        double areasum = 0;
        for (CGMesh::VertexFaceIter vf_it = lsmesh.vf_begin(vh); vf_it != lsmesh.vf_end(vh); ++vf_it)
        {
            assert(vh.idx() == i);
            int fid = vf_it.handle().idx();
            double area = areaF(fid);
            areasum += area;
            assert(i < output[0].size());
            output[0](i) += area * input[0](fid);
            output[1](i) += area * input[1](fid);
            output[2](i) += area * input[2](fid);
        }
        output[0](i) = output[0](i) / areasum;
        output[1](i) = output[1](i) / areasum;
        output[2](i) = output[2](i) / areasum;
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
        input(i).coeffRef(i) = i;
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
        double areasum = 0;
        for (CGMesh::VertexFaceIter vf_it = lsmesh.vf_begin(vh); vf_it != lsmesh.vf_end(vh); ++vf_it)
        {
            assert(vh.idx() == i);
            int fid = vf_it.handle().idx();
            double area = areaPF(fid);
            areasum += area;
            DonV[0].row(i) += area * DonF[0].row(fid);
            DonV[1].row(i) += area * DonF[1].row(fid);
        }
        DonV[0].row(i) = DonV[0].row(i) / areasum;
        DonV[1].row(i) = DonV[1].row(i) / areasum;
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
    assert(length = V.rows());
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
    gfvalue.resize(vsize, 3);
    hfvalue.resize(vsize);

    Eigen::VectorXd gx, gy, gz;
    get_lf_value(gradV[0], gx);
    get_lf_value(gradV[1], gy);
    get_lf_value(gradV[2], gz);
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