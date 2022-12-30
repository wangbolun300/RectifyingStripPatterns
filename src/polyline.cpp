#include <lsc/basic.h>
#include <lsc/tools.h>
#include<igl/file_dialog_open.h>
#include<igl/file_dialog_save.h>
void sample_polylines_and_binormals_evenly(const int nbr_segs, const std::vector<std::vector<Eigen::Vector3d>> &ply_in, const std::vector<std::vector<Eigen::Vector3d>> &bi_in,
                            std::vector<std::vector<Eigen::Vector3d>> &ply, std::vector<std::vector<Eigen::Vector3d>> &bi)
{
    double max_length = 0;
    ply.clear();
    bi.clear();
    std::vector<double> plylengths;
    for (int i = 0; i < ply_in.size(); i++)
    {
        double length = polyline_length(ply_in[i]);
        plylengths.push_back(length);
        if (length > max_length)
        {
            max_length = length;
        }
    }
    // std::cout<<"check 1"<<std::endl;
    double avg = max_length / nbr_segs; 
    for (int i = 0; i < ply_in.size(); i++)
    {
        std::vector<Eigen::Vector3d> ply_current = ply_in[i], bin_current = bi_in[i];
        std::vector<Eigen::Vector3d> ply_sampled, bin_sampled;
        int nbr = plylengths[i] / avg + 1;

        ply_sampled = sample_one_polyline_and_binormals_based_on_length(ply_current, nbr, bin_current, bin_sampled);
        ply.push_back(ply_sampled);
        bi.push_back(bin_sampled);
        // std::cout<<"pushed"<<std::endl;
    }
    // std::cout<<"check 2"<<std::endl;


}

void PolyOpt::init(const std::vector<std::vector<Eigen::Vector3d>> &ply_in, const std::vector<std::vector<Eigen::Vector3d>> &bi_in, const int sample_nbr)
{
    std::vector<std::vector<Eigen::Vector3d>> ply;
    std::vector<std::vector<Eigen::Vector3d>> bi;
    if (sample_nbr == -1) // do not sample the polylines
    {
        ply = ply_in;
        bi = bi_in;
    }
    else{
        sample_polylines_and_binormals_evenly(sample_nbr, ply_in, bi_in, ply, bi);
        std::cout<<"polyline get sampled"<<std::endl;
    }
    int vnbr = 0;
    for (auto line : ply)
    {
        for (auto ver : line)
        {
            vnbr++;
        }
    }
    // x0, x1, ..., xn, y0, ..., zn, bx0, ... bxn, ..., bzn
    PlyVars.resize(vnbr * 2 * 3);
    Front.resize(vnbr);
    Back.resize(vnbr);
    int counter = 0;
    for (int i = 0; i < ply.size(); i++)
    {
        for (int j = 0; j < ply[i].size(); j++)
        {
            int lpx = counter;
            int lpy = counter + vnbr;
            int lpz = counter + vnbr * 2;
            int lnx = counter + vnbr * 3;
            int lny = counter + vnbr * 4;
            int lnz = counter + vnbr * 5;
            Eigen::Vector3d ver = ply[i][j];
            Eigen::Vector3d bin = bi[i][j];
            PlyVars[lpx] = ver[0];
            PlyVars[lpy] = ver[1];
            PlyVars[lpz] = ver[2];
            PlyVars[lnx] = bin[0];
            PlyVars[lny] = bin[1];
            PlyVars[lnz] = bin[2];
            if (j == 0)
            {
                Front[counter] = -1;
            }
            else
            {
                Front[counter] = counter - 1;
            }
            if (j == ply[i].size() - 1)
            {
                Back[counter] = -1;
            }
            else
            {
                Back[counter] = counter + 1;
            }

            counter++;
        }
    }
    OriVars = PlyVars;
    RecMesh = polyline_to_strip_mesh(ply, bi, strip_scale);
    ply_extracted = ply;
    bin_extracted = bi;
}

void PolyOpt::assemble_gravity(spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy)
{
    if (Front.size() == 0)
    {
        std::cout << "please use assemble_gravity after init the PolyOpt" << std::endl;
        return;
    }

    // std::vector<Trip> tripletes;
    int vnbr = Front.size();
    H = Eigen::VectorXd::Ones(vnbr * 6).asDiagonal();
    energy = Eigen::VectorXd::Zero(vnbr * 6);
    // tripletes.reserve(vnbr * 6);
    // for (int i = 0; i < PlyVars.size(); i++)
    // {
    //     int rep = i / vnbr;
    //     int vid = i - rep * vnbr;
    //     energy[i] =
    //     // tripletes.push_back
    // }
    energy.segment(0, vnbr * 3) = (PlyVars - OriVars).segment(0, vnbr * 3);

    B = -energy;
}
void PolyOpt::assemble_polyline_smooth(spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy)
{
    if (Front.size() == 0)
    {
        std::cout << "please use assemble_polyline_smooth after init the PolyOpt" << std::endl;
        return;
    }

    std::vector<Trip> tripletes;
    int vnbr = Front.size();
    energy = Eigen::VectorXd::Zero(vnbr * 5);
    tripletes.reserve(vnbr * 25);
    // std::cout<<"check 1"<<std::endl;
    for (int i = 0; i < vnbr; i++)
    {

        int vid = i;
        int fr = Front[vid];
        int bk = Back[vid];
        bool compute_line_smth = true;
        bool compute_front = true;
        bool compute_back = true;
        if (fr == -1 || bk == -1)
        {
            compute_line_smth = false;
        }
        
        if (fr == -1)
        {
            fr = 0;
            compute_front = false;
        }
        if (bk == -1)
        {
            bk = 0;
            compute_back = false;
        }
        // locations
        int lx = vid;
        int ly = vid + vnbr;
        int lz = vid + vnbr * 2;
        int lnx = vid + vnbr * 3;
        int lny = vid + vnbr * 4;
        int lnz = vid + vnbr * 5;
        int lfx = fr;
        int lfy = fr + vnbr;
        int lfz = fr + vnbr * 2;
        int lbx = bk;
        int lby = bk + vnbr;
        int lbz = bk + vnbr * 2;
        int lfnx = fr + vnbr * 3;
        int lfny = fr + vnbr * 4;
        int lfnz = fr + vnbr * 5;
        int lbnx = bk + vnbr * 3;
        int lbny = bk + vnbr * 4;
        int lbnz = bk + vnbr * 5;
        Eigen::Vector3d ver(PlyVars[lx], PlyVars[ly], PlyVars[lz]);
        Eigen::Vector3d vfr(PlyVars[lfx], PlyVars[lfy], PlyVars[lfz]);
        Eigen::Vector3d vbk(PlyVars[lbx], PlyVars[lby], PlyVars[lbz]);
        Eigen::Vector3d nbi(PlyVars[lnx], PlyVars[lny], PlyVars[lnz]);
        Eigen::Vector3d nfr(PlyVars[lfnx], PlyVars[lfny], PlyVars[lfnz]);
        Eigen::Vector3d nbk(PlyVars[lbnx], PlyVars[lbny], PlyVars[lbnz]);

        if (compute_line_smth)
        {
            Eigen::Vector3d Over(OriVars[lx], OriVars[ly], OriVars[lz]);
            Eigen::Vector3d Ovfr(OriVars[lfx], OriVars[lfy], OriVars[lfz]);
            Eigen::Vector3d Ovbk(OriVars[lbx], OriVars[lby], OriVars[lbz]);

            double dis0 = (Over - Ovfr).norm();
            double dis1 = (Over - Ovbk).norm();
            double t = dis0 / (dis1 + dis0);
            // x = (1-t) * x_front + t * x_back
            tripletes.push_back(Trip(i, lx, 1));
            tripletes.push_back(Trip(i, lfx, t - 1));
            tripletes.push_back(Trip(i, lbx, -t));

            energy[i] = ver[0] - (1 - t) * vfr[0] - t * vbk[0];

            // y = (1-t) * y_front + t * y_back
            tripletes.push_back(Trip(i + vnbr, ly, 1));
            tripletes.push_back(Trip(i + vnbr, lfy, t - 1));
            tripletes.push_back(Trip(i + vnbr, lby, -t));

            energy[i + vnbr] = ver[1] - (1 - t) * vfr[1] - t * vbk[1];

            // z = (1-t) * z_front + t * z_back
            tripletes.push_back(Trip(i + vnbr * 2, lz, 1));
            tripletes.push_back(Trip(i + vnbr * 2, lfz, t - 1));
            tripletes.push_back(Trip(i + vnbr * 2, lbz, -t));

            energy[i + vnbr * 2] = ver[2] - (1 - t) * vfr[2] - t * vbk[2];
        }

        double sign;
        if (compute_front)
        {
            sign = nbi.dot(nfr) > 0 ? 1 : -1;
            // nbi * nfr +- 1 = 0
            tripletes.push_back(Trip(i + vnbr * 3, lnx, binormal_ratio * nfr[0]));
            tripletes.push_back(Trip(i + vnbr * 3, lny, binormal_ratio * nfr[1]));
            tripletes.push_back(Trip(i + vnbr * 3, lnz, binormal_ratio * nfr[2]));

            tripletes.push_back(Trip(i + vnbr * 3, lfnx, binormal_ratio * nbi[0]));
            tripletes.push_back(Trip(i + vnbr * 3, lfny, binormal_ratio * nbi[1]));
            tripletes.push_back(Trip(i + vnbr * 3, lfnz, binormal_ratio * nbi[2]));

            energy[i + vnbr * 3] = binormal_ratio * (nbi.dot(nfr) - sign);
        }

        if (compute_back)
        {
            // nbi * nbk +- 1 = 0
            sign = nbi.dot(nbk) > 0 ? 1 : -1;

            tripletes.push_back(Trip(i + vnbr * 4, lnx, binormal_ratio * nbk[0]));
            tripletes.push_back(Trip(i + vnbr * 4, lny, binormal_ratio * nbk[1]));
            tripletes.push_back(Trip(i + vnbr * 4, lnz, binormal_ratio * nbk[2]));

            tripletes.push_back(Trip(i + vnbr * 4, lbnx, binormal_ratio * nbi[0]));
            tripletes.push_back(Trip(i + vnbr * 4, lbny, binormal_ratio * nbi[1]));
            tripletes.push_back(Trip(i + vnbr * 4, lbnz, binormal_ratio * nbi[2]));

            energy[i + vnbr * 4] = binormal_ratio * (nbi.dot(nbk) - sign);
        }
    }

    // std::cout<<"check 2"<<std::endl;
    spMat J;
    J.resize(energy.size(), vnbr * 6);
    J.setFromTriplets(tripletes.begin(), tripletes.end());
    H = J.transpose() * J;
    B = -J.transpose() * energy;
    // std::cout<<"check 3"<<std::endl;
}
void PolyOpt::assemble_binormal_condition(spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy)
{
    if (Front.size() == 0)
    {
        std::cout << "please use assemble_binormal_condition after init the PolyOpt" << std::endl;
        return;
    }

    std::vector<Trip> tripletes;
    int vnbr = Front.size();
    energy = Eigen::VectorXd::Zero(vnbr * 3);
    tripletes.reserve(vnbr * 25);
    for (int i = 0; i < vnbr; i++)
    {
        int vid = i;
        int fr = Front[vid];
        int bk = Back[vid];
        bool compute_front = true;
        bool compute_back = true;
        
        if (fr == -1)
        {
            fr = 0;
            compute_front = false;
        }
        if (bk == -1)
        {
            bk = 0;
            compute_back = false;
        }
        // locations
        int lx = vid;
        int ly = vid + vnbr;
        int lz = vid + vnbr * 2;
        int lnx = vid + vnbr * 3;
        int lny = vid + vnbr * 4;
        int lnz = vid + vnbr * 5;
        int lfx = fr;
        int lfy = fr + vnbr;
        int lfz = fr + vnbr * 2;
        int lbx = bk;
        int lby = bk + vnbr;
        int lbz = bk + vnbr * 2;
        int lfnx = fr + vnbr * 3;
        int lfny = fr + vnbr * 4;
        int lfnz = fr + vnbr * 5;
        int lbnx = bk + vnbr * 3;
        int lbny = bk + vnbr * 4;
        int lbnz = bk + vnbr * 5;
        Eigen::Vector3d ver(PlyVars[lx], PlyVars[ly], PlyVars[lz]);
        Eigen::Vector3d vfr(PlyVars[lfx], PlyVars[lfy], PlyVars[lfz]);
        Eigen::Vector3d vbk(PlyVars[lbx], PlyVars[lby], PlyVars[lbz]);
        Eigen::Vector3d nbi(PlyVars[lnx], PlyVars[lny], PlyVars[lnz]);
        Eigen::Vector3d nfr(PlyVars[lfnx], PlyVars[lfny], PlyVars[lfnz]);
        Eigen::Vector3d nbk(PlyVars[lbnx], PlyVars[lbny], PlyVars[lbnz]);

        double dis0 = (ver - vfr).norm();
        double dis1 = (ver - vbk).norm();

        Eigen::Vector3d seg;
        if (compute_front)
        {
            // (ver-vfr).dot(nbi) = 0
            tripletes.push_back(Trip(i, lx, nbi[0] / dis0));
            tripletes.push_back(Trip(i, ly, nbi[1] / dis0));
            tripletes.push_back(Trip(i, lz, nbi[2] / dis0));

            tripletes.push_back(Trip(i, lfx, -nbi[0] / dis0));
            tripletes.push_back(Trip(i, lfy, -nbi[1] / dis0));
            tripletes.push_back(Trip(i, lfz, -nbi[2] / dis0));

            seg = ver - vfr;
            tripletes.push_back(Trip(i, lnx, seg[0] / dis0));
            tripletes.push_back(Trip(i, lny, seg[1] / dis0));
            tripletes.push_back(Trip(i, lnz, seg[2] / dis0));

            energy[i] = seg.dot(nbi) / dis0;
        }
        if (compute_back)
        {
            // (ver-vbk).dot(nbi) = 0
            tripletes.push_back(Trip(i + vnbr, lx, nbi[0] / dis1));
            tripletes.push_back(Trip(i + vnbr, ly, nbi[1] / dis1));
            tripletes.push_back(Trip(i + vnbr, lz, nbi[2] / dis1));

            tripletes.push_back(Trip(i + vnbr, lbx, -nbi[0] / dis1));
            tripletes.push_back(Trip(i + vnbr, lby, -nbi[1] / dis1));
            tripletes.push_back(Trip(i + vnbr, lbz, -nbi[2] / dis1));

            seg = ver - vbk;
            tripletes.push_back(Trip(i + vnbr, lnx, seg[0] / dis1));
            tripletes.push_back(Trip(i + vnbr, lny, seg[1] / dis1));
            tripletes.push_back(Trip(i + vnbr, lnz, seg[2] / dis1));

            energy[i + vnbr] = seg.dot(nbi) / dis1;
        }

        // nbi * nbi = 1
        tripletes.push_back(Trip(i + vnbr * 2, lnx, 2 * nbi[0]));
        tripletes.push_back(Trip(i + vnbr * 2, lny, 2 * nbi[1]));
        tripletes.push_back(Trip(i + vnbr * 2, lnz, 2 * nbi[2]));

        energy[i + vnbr * 2] = nbi.dot(nbi) - 1;
    }
    spMat J;
    J.resize(energy.size(), vnbr * 6);
    J.setFromTriplets(tripletes.begin(), tripletes.end());
    H = J.transpose() * J;
    B = -J.transpose() * energy;
}

void PolyOpt::extract_rectifying_plane_mesh()
{
    int nv = RecMesh.n_vertices() / 2;
    Eigen::MatrixXd v0list(nv, 3), v1list(nv, 3);
    Eigen::MatrixXd binlist(nv, 3);
    for (int i = 0; i < nv; i++)
    {
        v0list(i, 0) = PlyVars(i);
        v0list(i, 1) = PlyVars(i + nv);
        v0list(i, 2) = PlyVars(i + nv * 2);
        binlist(i, 0) = PlyVars(i + nv * 3);
        binlist(i, 1) = PlyVars(i + nv * 4);
        binlist(i, 2) = PlyVars(i + nv * 5);
    }

    v1list = v0list + strip_scale * binlist;

    for (CGMesh::VertexIter v_it = RecMesh.vertices_begin(); v_it != RecMesh.vertices_end(); ++v_it)
    {
        int vid = v_it.handle().idx();
        int v0_id = vid / 2;
        Eigen::Vector3d pt;
        if (vid % 2 == 0) // this is a vertex of the polyline
        {
            pt = v0list.row(v0_id);
        }
        else
        { // this is a offset vertex
            pt = v1list.row(v0_id);
        }
        RecMesh.point(*v_it) = CGMesh::Point(pt(0), pt(1), pt(2));
    }
}
void PolyOpt::polyline_to_matrix(const std::vector<std::vector<Eigen::Vector3d>> &ply, Eigen::MatrixXd &V)
{
    std::vector<Eigen::Vector3d> verlist;
    for (auto line : ply)
    {
        for (auto pt : line)
        {
            verlist.push_back(pt);
        }
    }
    V = vec_list_to_matrix(verlist);

}
void PolyOpt::vertex_matrix_to_polyline(std::vector<std::vector<Eigen::Vector3d>> &ply, Eigen::MatrixXd& V){
    int counter = 0;
    for (int i =0;i<ply.size();i++)
    {
        for (int j =0;j<ply[i].size();j++)
        {
            Eigen::Vector3d pt = V.row(counter);
            ply[i][j] = pt;
            counter++;
        }
    }
}
void PolyOpt::rotate_back_the_model_to_horizontal_coordinates(const double latitude_degree)
{
    double latitude_radian = latitude_degree * LSC_PI / 180.;
    double rho = -(LSC_PI / 2 - latitude_radian); // this is correct rotation
    Eigen::Matrix3d rotation;
    rotation << 1, 0, 0,
        0, cos(rho), -sin(rho),
        0, sin(rho), cos(rho);
    ply_rotated = ply_extracted;
    bin_rotated = bin_extracted;
    for (int i = 0; i < ply_rotated.size(); i++)
    {
        for (int j = 0; j < ply_rotated[i].size(); j++)
        {
            ply_rotated[i][j] = rotation * ply_rotated[i][j];
            bin_rotated[i][j] = rotation * bin_rotated[i][j];
        }
    }
}
void PolyOpt::extract_polylines_and_binormals()
{
    int vnbr = PlyVars.size() / 6;
    int counter = 0;
    for (int i = 0; i < ply_extracted.size(); i++)
    {
        for (int j = 0; j < ply_extracted[i].size(); j++)
        {
            ply_extracted[i][j][0] = PlyVars[counter];
            ply_extracted[i][j][1] = PlyVars[counter + vnbr];
            ply_extracted[i][j][2] = PlyVars[counter + vnbr * 2];
            bin_extracted[i][j][0] = PlyVars[counter + vnbr * 3];
            bin_extracted[i][j][1] = PlyVars[counter + vnbr * 4];
            bin_extracted[i][j][2] = PlyVars[counter + vnbr * 5];
            counter ++;
        }
    }
    assert(counter == vnbr);
}
// if rotated, we save the rotated mesh; otherwise we save the mesh extracted from optimizer.
void PolyOpt::save_polyline_and_binormals_as_files(const bool rotated)
{
    std::cout << "Saving the binormal files, please provide the prefix" << std::endl;
    std::string fname = igl::file_dialog_save();
    std::vector<std::vector<Eigen::Vector3d>> lines, binormals;
    if (rotated)
    {
        lines = ply_rotated;
        binormals = bin_rotated;
    }
    else
    {
        lines = ply_extracted;
        binormals = bin_extracted;
    }
    write_polyline_xyz(lines, fname);
    write_polyline_xyz(binormals, fname + "_b");
    std::cout << "files get saved" << std::endl;
}
void PolyOpt::save_polyline_and_binormals_as_files(const std::vector<std::vector<Eigen::Vector3d>> &ply, const std::vector<std::vector<Eigen::Vector3d>> &bi)
{
    std::string fname = igl::file_dialog_save();
    write_polyline_xyz(ply, fname);
    write_polyline_xyz(bi, fname + "_b");
    std::cout << "files get saved" << std::endl;
}
void PolyOpt::save_polyline_and_binormals_as_files(const std::string &fname,
                                                   const std::vector<std::vector<Eigen::Vector3d>> &ply, const std::vector<std::vector<Eigen::Vector3d>> &bi)
{
    write_polyline_xyz(ply, fname);
    write_polyline_xyz(bi, fname + "_b");
    std::cout << "files get saved" << std::endl;
}
void PolyOpt::orient_binormals_of_plyline(const std::vector<std::vector<Eigen::Vector3d>> &bi, std::vector<std::vector<Eigen::Vector3d>> &bout)
{
    std::vector<std::vector<Eigen::Vector3d>> bilist = bi;
    // int counter = 0;
    Eigen::Vector3d bbase = bi[0][0];
    for (int i = 0; i < bi.size(); i++)
    {
        Eigen::Vector3d base = bi[i][0];
        if (base.dot(bbase) < 0)
        {
            base *= -1;
        }

        for (int j = 0; j < bi[i].size(); j++)
        {
            Eigen::Vector3d direction = orient_vector(base, bi[i][j]);
            base = direction;
            bilist[i][j] = direction;
        }
    }
    bout = bilist;
}
void PolyOpt::opt()
{
    spMat H;
    Eigen::VectorXd B;

    spMat Hmass;
    Eigen::VectorXd Bmass, emass;

    assemble_gravity(Hmass, Bmass, emass);

    H = weight_mass * 1e-6 * Hmass;
    B = weight_mass * 1e-6 * Bmass;

    spMat Hsmt;
    Eigen::VectorXd Bsmt, esmt;

    assemble_polyline_smooth(Hsmt, Bsmt, esmt);

    H += weight_smooth * Hsmt;
    B += weight_smooth * Bsmt;

    spMat Hbin;
    Eigen::VectorXd Bbin, ebin;

    assemble_binormal_condition(Hbin, Bbin, ebin);

    H += weight_binormal * Hbin;
    B += weight_binormal * Bbin;

    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver(H);

    // assert(solver.info() == Eigen::Success);
    if (solver.info() != Eigen::Success)
    {
        // solving failed
        std::cout << "solver fail" << std::endl;
        return;
    }
    // std::cout<<"solved successfully"<<std::endl;
    Eigen::VectorXd dx = solver.solve(B).eval();
    double level_set_step_length = dx.norm();
    if (level_set_step_length > max_step)
    {
        dx *= max_step / level_set_step_length;
    }
    PlyVars += dx;

    double energy_gravity = emass.norm();
    double energy_smoothness = esmt.norm();
    double energy_binormal = ebin.norm();
    double stplength = dx.norm();

    std::cout << "Ply mass, " << energy_gravity << ", smth, " << energy_smoothness << ", binormal, " << energy_binormal << ", step, " << stplength << std::endl;
    extract_rectifying_plane_mesh();
    extract_polylines_and_binormals();
    //
}