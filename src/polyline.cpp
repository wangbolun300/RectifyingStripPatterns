#include <lsc/basic.h>
#include <lsc/tools.h>
#include<igl/file_dialog_open.h>
#include<igl/file_dialog_save.h>
bool print_info = false;
void solve_project_vector_on_plane(const Eigen::Vector3d &vec, const Eigen::Vector3d &a1, const Eigen::Vector3d &a2,
                                   Eigen::Vector3d &v2d)
{

    Eigen::Matrix2d mat, matinv;
    mat << a1.dot(a1), a1.dot(a2),
        a1.dot(a2), a2.dot(a2);
    Eigen::Vector2d B(vec.dot(a1), vec.dot(a2));
    bool invertable;

    mat.computeInverseWithCheck(matinv, invertable);
    if(!invertable){
        std::cout<<"ERROR: the two input vectors are linearly dependent"<<std::endl;
        assert(0);
        return;
    }
    Eigen::Vector2d ab = matinv * B;
    v2d = ab[0] * a1 + ab[1] * a2;
    if(print_info){
        std::cout<<"\nmat,\n"<<mat<<"\nB,"<<B.transpose()<<"\nmatinv\n"<<matinv<<"\na,b, "<<ab[0]<<", "<<ab[1]<<"\nv2d, "
        <<v2d.transpose()<<"\n";
    }
}

void convert_polyline_to_developable_strips(const std::vector<Eigen::Vector3d> &ply,
                                            const std::vector<Eigen::Vector3d> &bnm, std::vector<Eigen::Vector3d> &creases)
{
    // creases[0] = bnm[0];
    creases = bnm;
    int vnbr = ply.size();
    print_info = false;
    //ply.size() > 3;
    // double cosbeta1 = 1; // the angle between the left plane and the x-z plane. The first is the
    // double cosbeta2 = 0; // between the left plane and the x-y plane.
    Eigen::Vector3d last_crease = creases[0];
    for (int i = 1; i < vnbr - 1; i++)
    {
        Eigen::Vector3d vf = ply[i - 1];
        Eigen::Vector3d vb = ply[i + 1];
        Eigen::Vector3d pt = ply[i];

        Eigen::Vector3d d1 = (pt - vf).normalized();
        Eigen::Vector3d d2 = (vb - pt).normalized();
        Eigen::Vector3d bi = bnm[i];
        bi = bi.normalized();
        Eigen::Vector3d yy = (bi.cross(-d1)).normalized();
        // project -d1 to x axis, yy to y axis, bi to z axis.
        Eigen::Matrix3d Rinv, R;
        Rinv.col(0) = -d1;
        Rinv.col(1) = yy;
        Rinv.col(2) = bi;
        bool invertable;
        Rinv.computeInverseWithCheck(R, invertable);
        if (!invertable)
        {
            std::cout << "Please check why it is not invertable here" << std::endl;
            return;
        }
        // project the crease to y-z plane
        Eigen::Vector3d proj_crease;
        solve_project_vector_on_plane(last_crease, bi, yy, proj_crease);
        proj_crease = proj_crease.normalized();
        double cosbeta1 = proj_crease.dot(bi);
        double cosbeta2 = proj_crease.dot(yy);
        // R is the projection matrix, project the 3 vectors to the standard axis.
        // the following calculation is in the new coordinate system

        Eigen::Vector3d D = R * d2;
        Eigen::Vector3d C(0, cosbeta2, cosbeta1);
        double D0 = D(0), D1 = D(1), D2 = D(2);
        double a, b;                   // the crease is defined as a*C+b*(1,0,0)
        if (abs(1 + D0) < SCALAR_ZERO) // unlike to happen
        {
            a = 0;
            b = 1;
        }
        else
        {
            a = 1;
            b = -(D1 * cosbeta2 + D2 * cosbeta1) / (1 + D0);
        }
        Eigen::Vector3d A = Eigen::Vector3d(b, a * cosbeta2, a * cosbeta1);
        if (A(2) < 0) // to make sure A is the same orient with the binormal vector
        {
            A *= -1;
        }
        double ratio = sqrt(1 / (A.dot(A) - A(0) * A(0)));
        A *= ratio;
        creases[i] = Rinv * A;
        last_crease = creases[i];
        if (print_info)
        {
            std::cout << "R\n"
                      << R << "\nbi, " << bi.transpose() << "\nyy, " << yy.transpose() << "\nproj_crease, " << proj_crease.transpose()
                      << "\ncrease, " << last_crease.transpose() << "\n";
            exit(0);
        }
    }
    // solve for the last vector
    Eigen::Vector3d d1 = (ply[vnbr - 1] - ply[vnbr - 2]).normalized();
    Eigen::Vector3d bi = bnm[vnbr - 1];
    bi = bi.normalized();
    Eigen::Vector3d yy = (bi.cross(-d1)).normalized();
    // the last crease
    solve_project_vector_on_plane(last_crease, bi, yy, creases[vnbr - 1]);
}

// rul_left is the ruling of the left plane,
void reflect_wrt_rectifying(const Eigen::Vector3d &ver0, const Eigen::Vector3d &ver1, const Eigen::Vector3d &ver2, const Eigen::Vector3d &rul_left,
                            const Eigen::Vector3d &bi_middle, Eigen::Vector3d &crease, Eigen::Vector3d rul_right)
{
    Eigen::Vector3d dir0 = (ver1 - ver0).normalized();
    Eigen::Vector3d dir1 = (ver2 - ver1).normalized();

    double straight_tol = 1e-2; // 0.01
    if (dir0.dot(dir1) > 1 - straight_tol)
    { // the curve does not bend much, meaning that it is a flat one
        crease = rul_left;
        rul_right = rul_left;
        return;
    }
    if (rul_left.dot(bi_middle) > 1 - straight_tol)
    { // the ruling of the left surface is the same as the right one, is the binormal
        rul_right = bi_middle;
        crease = bi_middle;
        return;
    }

    // P0, P1 define the mirror
    Eigen::Vector3d P0 = (dir0 + dir1).normalized();
    Eigen::Vector3d P1 = bi_middle.normalized();

    // the other edge of the left strip
    Eigen::Vector3d prl_start = ver0 + rul_left;
    // the direction is dir0. (prl_start + t * dir0 - ver1, p0, p1) = 0
    // solve the intersection point
    Eigen::Matrix3d c0, c1;
    c0.col(0) = prl_start - ver1;
    c0.col(1) = P0;
    c0.col(2) = P1;
    c1.col(0) = dir0;
    c1.col(1) = P0;
    c1.col(2) = P1;
    double coff0 = c0.determinant();
    double coff1 = c1.determinant();
    if (abs(coff1) < 1e-6)
    {
        std::cout << "ERROR: the vector is parallel to the plane" << std::endl;
    }
    double t = -coff0 / coff1;
    Eigen::Vector3d proj = prl_start + t * dir0; // the intersection point between the other side of the strip and the mirror
    crease = proj - ver1;
    // the ruling is proj - (v1 + t * dir1)
    t = (proj - ver1).dot(dir1);
    rul_right = proj - (ver1 + t * dir1);
    double scale = rul_right.norm();
    if (abs(scale - 1) > SCALAR_ZERO)
    {
        std::cout << "The ruling is not accurate, the length of it is " << scale << std::endl;
    }
    rul_right.normalize();
}

void convert_polyline_to_developable_strips_reflection_method(const std::vector<Eigen::Vector3d> &ply,
                                                              const std::vector<Eigen::Vector3d> &bnm, std::vector<Eigen::Vector3d> &creases)
{
    creases = bnm;
    int vnbr = ply.size();
    Eigen::Vector3d ruling;
    ruling = (bnm[0] + bnm[1]).normalized();
    for (int i = 1; i < vnbr - 1; i++)
    {
        Eigen::Vector3d rul_right, crs;
        reflect_wrt_rectifying(ply[i - 1], ply[i], ply[i + 1], ruling, bnm[i], crs, rul_right);
        creases[i] = crs;
        ruling = rul_right;
    }
    creases[vnbr - 1] = ruling;
}



void evaluate_strip_developability(const std::vector<Eigen::Vector3d> &ply,
                                   const std::vector<Eigen::Vector3d> &crease, double &planarity, double &angle_diff, double &max_angle_error)
{
    // use the distances between diagonals to evaluate planarity
    double planar_total = 0;
    double angle_total = 0;
    double angle_max = 0;
    int vnbr = ply.size();
    for (int i = 0; i < vnbr - 1; i++)
    {
        Eigen::Vector3d v0 = ply[i];
        Eigen::Vector3d v1 = ply[i + 1];
        Eigen::Vector3d v2 = ply[i + 1] + crease[i + 1];
        Eigen::Vector3d v3 = ply[i] + crease[i];
        // dirction of diagonals
        Eigen::Vector3d dir0 = (v2 - v0).normalized();
        Eigen::Vector3d dir1 = (v3 - v1).normalized();
        double distance = dir0.cross(dir1).dot(v0 - v1);
        distance = abs(distance) / dir0.cross(dir1).norm();
        planar_total += distance;

        // double dis_left = (v0 - v1).norm();
        // double dis_right = (v2 - v3).norm();
        // dis_total_left += dis_left;
        // dis_total_right += dis_right;
        if (i == 0)
        {
            continue;
        }
        Eigen::Vector3d dir_middle = crease[i].normalized();
        Eigen::Vector3d dir_left = (ply[i - 1] - ply[i]).normalized();
        Eigen::Vector3d dir_right = (ply[i + 1] - ply[i]).normalized();
        double angle_left = acos(dir_middle.dot(dir_left));
        double angle_right = acos(dir_middle.dot(dir_right));
        double angle_err = fabs(angle_left + angle_right - LSC_PI);
        if (angle_err > angle_max)
        {
            angle_max = angle_err;
            if (print_info)
            {
                std::cout << "find max: i " << i << ", angle diff, " << angle_err << std::endl;
            }
        }

        angle_total += angle_err;
    }
    planarity = planar_total;
    angle_diff = angle_total;
    max_angle_error = angle_max;
    // length_diff = abs(dis_total_left - dis_total_right);
}

void evaluate_and_print_strip_developability(const std::vector<std::vector<Eigen::Vector3d>> &ply,
                                             const std::vector<std::vector<Eigen::Vector3d>> &crease)
{
    double planar_total = 0;
    double anglediff_total = 0;
    double max_angle_error = 0;
    for (int i = 0; i < ply.size(); i++)
    {
        double planar, anglediff, angle_max;
        evaluate_strip_developability(ply[i], crease[i], planar, anglediff, angle_max);
        planar_total += planar;
        anglediff_total += anglediff;
        if (angle_max > max_angle_error)
        {
            max_angle_error = angle_max;
        }
    }
    std::cout << "Evaluate current polylines: planarity: " << planar_total << ", total angle changes of the strip, " << anglediff_total << ", maximal angle error, "
              << max_angle_error << std::endl;
}

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
    Eigen::VectorXd endpts;
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
    VerNbr = vnbr;
    // x0, x1, ..., xn, y0, ..., zn, bx0, ... bxn, ..., bzn
    PlyVars.resize(vnbr * 2 * 3);
    endpts = Eigen::VectorXd::Zero(vnbr * 6);
    VinPly = Eigen::VectorXd::Ones(vnbr * 6) * -1;
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
            if (j == 0) // first point of the seg
            {
                endpts[lpx] = 1;
                endpts[lpy] = 1;
                endpts[lpz] = 1;
                // endpts[lnx] = 1;
                // endpts[lny] = 1;
                // endpts[lnz] = 1;
                Front[counter] = -1;
            }
            else
            {
                Front[counter] = counter - 1;
            }
            if (j == ply[i].size() - 1) // last point of the seg
            {
                endpts[lpx] = 1;
                endpts[lpy] = 1;
                endpts[lpz] = 1;
                // endpts[lnx] = 1;
                // endpts[lny] = 1;
                // endpts[lnz] = 1;
                Back[counter] = -1;
            }
            else
            {
                Back[counter] = counter + 1;
            }
            VinPly[lpx] = i;
            VinPly[lpy] = i;
            VinPly[lpz] = i;

            VinPly[lnx] = i;
            VinPly[lny] = i;
            VinPly[lnz] = i;

            counter++;
        }
    }
    OriVars = PlyVars;
    RecMesh = polyline_to_strip_mesh(ply, bi, strip_scale);
    ply_extracted = ply;
    bin_extracted = bi;
    endpts_signs = endpts.asDiagonal();
}
// use this when avoiding flipping of the strips.
void PolyOpt::force_smoothing_binormals(){
    int counter = 0;
    for (int i = 0; i < ply_extracted.size(); i++)
    {
        int lnx = counter + VerNbr * 3;
        int lny = counter + VerNbr * 4;
        int lnz = counter + VerNbr * 5;
        
        Eigen::Vector3d bin = Eigen::Vector3d(PlyVars[lnx], PlyVars[lny], PlyVars[lnz]);
        for (int j = 0; j < ply_extracted[i].size(); j++)
        {
            lnx = counter + VerNbr * 3;
            lny = counter + VerNbr * 4;
            lnz = counter + VerNbr * 5;
            if (pick_single_line)
            {
                if (pick_line_id == i)
                {
                    PlyVars[lnx] = bin[0];
                    PlyVars[lny] = bin[1];
                    PlyVars[lnz] = bin[2];
                }
            }
            else
            {
                PlyVars[lnx] = bin[0];
                PlyVars[lny] = bin[1];
                PlyVars[lnz] = bin[2];
            }

            counter++;
        }
    }
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
    H = spMat(Eigen::VectorXd::Ones(vnbr * 6).asDiagonal()) + endpts_signs * (ratio_endpts - 1);
    energy = Eigen::VectorXd::Zero(vnbr * 6);
    energy.segment(0, vnbr * 3) = (PlyVars - OriVars).segment(0, vnbr * 3);

    // pick the selected one.
    if (pick_single_line)
    {
        for (int i = 0; i < VinPly.size(); i++)
        {
            if (VinPly[i] != pick_line_id)
            {
                H.coeffRef(i, i) = 0;
                energy[i] = 0;
            }
        }
    }

    B = -H * energy;
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
        if (pick_single_line)
        {
            if (pick_line_id != VinPly[vid])
            {
                continue;
            }
        }
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
void PolyOpt::assemble_smooth_neighbouring_polyline_binormal(spMat& H, Eigen::VectorXd& B, Eigen::VectorXd &energy){
    
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
        if (pick_single_line)
        {
            if (pick_line_id != VinPly[vid])
            {
                continue;
            }
        }
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
// the binormals has an angle theta with surface normals
void PolyOpt::assemble_angle_condition(spMat& H, Eigen::VectorXd& B, Eigen::VectorXd &energy){
    
    if (Front.size() == 0)
    {
        std::cout << "please use assemble_angle_condition after init the PolyOpt" << std::endl;
        return;
    }

    std::vector<Trip> tripletes;
    int vnbr = Front.size();
    energy = Eigen::VectorXd::Zero(vnbr * 2);
    tripletes.reserve(vnbr * 6);
    double angle_radian = target_angle * LSC_PI / 180.;
    double cos_angle = cos(angle_radian);
    double sin_angle = sin(angle_radian);
    for (int i = 0; i < vnbr; i++)
    {
        int vid = i;
        if (pick_single_line)
        {
            if (pick_line_id != VinPly[vid])
            {
                continue;
            }
        }
        int fr = Front[vid];
        int bk = Back[vid];
        bool compute_front = true;
        bool compute_back = true;
        
        if (fr == -1)
        {
            fr = 0;
            compute_front = false;
            continue;
        }
        if (bk == -1)
        {
            bk = 0;
            compute_back = false;
            continue;
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

        // int lux = vid + vnbr * 6;
        // int luy = vid + vnbr * 7;
        // int luz = vid + vnbr * 8;
        // int lh = vid + vnbr * 9;
        Eigen::Vector3d ver(PlyVars[lx], PlyVars[ly], PlyVars[lz]);
        Eigen::Vector3d vfr(PlyVars[lfx], PlyVars[lfy], PlyVars[lfz]);
        Eigen::Vector3d vbk(PlyVars[lbx], PlyVars[lby], PlyVars[lbz]);
        Eigen::Vector3d nbi(PlyVars[lnx], PlyVars[lny], PlyVars[lnz]);
        // Eigen::Vector3d nfr(PlyVars[lfnx], PlyVars[lfny], PlyVars[lfnz]);
        // Eigen::Vector3d nbk(PlyVars[lbnx], PlyVars[lbny], PlyVars[lbnz]);
        Eigen::Vector3d tangent = -(vbk - vfr);// This sign is changable depends on what user defines of the curve tangent directions
        Eigen::Vector3d norm = norm_v.row(vid);
        // if(first_compute){
        //     Eigen::Vector3d real_u = norm.cross(tangent).normalized();
        //     Eigen::Vector3d real_h = nbi.dot(real_u);
        //     PlyVars[lux] = real_u[0];
        //     PlyVars[luy] = real_u[1];
        //     PlyVars[luz] = real_u[2];
        //     PlyVars[lh] = real_h;
        //      first_compute = false;
        // }
        Eigen::Vector3d real_u = (norm.cross(tangent)).normalized();

        double dis0 = (ver - vfr).norm();
        double dis1 = (ver - vbk).norm();

        // (binormal * n)^2 = cos^2
        double ndb = norm.dot(nbi);
        tripletes.push_back(Trip(i, lnx, 2 * ndb * norm(0)));
        tripletes.push_back(Trip(i, lny, 2 * ndb * norm(1)));
        tripletes.push_back(Trip(i, lnz, 2 * ndb * norm(2)));

        energy[i] = ndb * ndb - cos_angle * cos_angle;

        // (binormal * n) * (binormal * u) = sin * cos
        // regard n and u as constant values in each iteration, since u is related to the tangent vector, which will not change much.
        double udb = real_u.dot(nbi);

        tripletes.push_back(Trip(i + vnbr, lnx, ndb * real_u[0] + udb * norm[0]));
        tripletes.push_back(Trip(i + vnbr, lny, ndb * real_u[1] + udb * norm[1]));
        tripletes.push_back(Trip(i + vnbr, lnz, ndb * real_u[2] + udb * norm[2]));

        energy[i + vnbr] = ndb * udb - sin_angle * cos_angle;
    }
    spMat J;
    J.resize(energy.size(), vnbr * 6);
    J.setFromTriplets(tripletes.begin(), tripletes.end());
    H = J.transpose() * J;
    B = -J.transpose() * energy;
}
#include <igl/point_mesh_squared_distance.h>
// please use me after initializing the polylines
// Vr, Fr and nr are for the reference surfaces
void PolyOpt::get_normal_vector_from_reference(const Eigen::MatrixXd &Vr, const Eigen::MatrixXi &Fr, const Eigen::MatrixXd &nr)
{
    Eigen::MatrixXd C;
    Eigen::VectorXi I;
    Eigen::VectorXd D;
    Eigen::MatrixXd vers;
    int vnbr = VerNbr;
    vers.resize(vnbr,3);
    norm_v.resize(vnbr,3);
    for (int i = 0; i < vnbr; i++)
    {
        vers(i, 0) = PlyVars[i];
        vers(i, 1) = PlyVars[i + vnbr];
        vers(i, 2) = PlyVars[i + vnbr * 2];
    }
    igl::point_mesh_squared_distance(vers, Vr, Fr, D, I, C);
    for (int i = 0; i < vnbr; i++)
    {
        int fid = I(i);
        int v0 = Fr(fid, 0);
        int v1 = Fr(fid, 1);
        int v2 = Fr(fid, 2);
        Eigen::Vector3d pt = vers.row(i);
        Eigen::Vector3d ver0 = Vr.row(v0);
        Eigen::Vector3d ver1 = Vr.row(v1);
        Eigen::Vector3d ver2 = Vr.row(v2);
        std::array<double, 3> coor = barycenter_coordinate(ver0, ver1, ver2, pt);

        Eigen::Vector3d dir = coor[0] * nr.row(v0) + coor[1] * nr.row(v1) + coor[2] * nr.row(v2);
        dir = dir.normalized();
        norm_v.row(i) = dir;
    }

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
    std::cout<<"saving polyline files, please provide the prefix"<<std::endl;
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
        bbase = base;

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

    spMat Hangle;
    Eigen::VectorXd Bangle, Eangle;
    if(weight_angle > 0){
        assemble_angle_condition(Hangle, Bangle,Eangle);
        H += weight_angle * Hangle;
        B += weight_angle * Bangle;
    }
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
    
    std::cout << "Ply mass, " << energy_gravity << ", smth, " << energy_smoothness << ", binormal, " << energy_binormal << ", step, " << stplength;
    if(weight_angle > 0){
        double energy_angle = Eangle.norm();
        std::cout<<", Energy angle, "<<energy_angle<<", Eangle max, "<< Eangle.lpNorm<Eigen::Infinity>();
    }
    std::cout<<"\n";
    extract_rectifying_plane_mesh();
    extract_polylines_and_binormals();
    //
}
std::vector<std::vector<int>> load_row_col_info(const std::vector<std::vector<double>> &rowstmp)
{
    std::vector<std::vector<double>> tmp = rowstmp;
    std::vector<std::vector<int>> result;
    result.resize(rowstmp.size());
    for (int i = 0; i < result.size(); i++)
    {
        std::vector<int> curve;
        curve.resize(rowstmp[i].size());
        for (int j = 0; j < curve.size(); j++)
        {
            curve[j] = rowstmp[i][j];
        }
        result[i]=curve;
    }
    return result;
}
void load_row_col_info(const std::vector<std::vector<double>> &rowstmp, Eigen::VectorXi &row_f, Eigen::VectorXi &row_b)
{
    for (int i = 0; i < rowstmp.size(); i++)
    {
        for (int j = 0; j < rowstmp[i].size(); j++)
        {
            int vid = rowstmp[i][j];
            if (vid < 0)
            {
                continue;
            }
            int front = -1;
            int back = -1;
            if (j != 0)
            {
                front = rowstmp[i][j - 1];
            }
            if (j != rowstmp[i].size() - 1)
            {
                back = rowstmp[i][j + 1];
            }
            row_f[vid] = front;
            row_b[vid] = back;
        }
    }
    return;
}

void load_diagonal_info(const Eigen::VectorXi &row_front, const Eigen::VectorXi &row_back, const Eigen::VectorXi &col_front,
                        const Eigen::VectorXi &col_back, Eigen::VectorXi &d0_front, Eigen::VectorXi &d1_front,
                        Eigen::VectorXi &d0_back, Eigen::VectorXi &d1_back)
{
    int vnbr = row_front.size();
    d0_front = Eigen::VectorXi::Ones(vnbr) * -1;
    d0_back = d0_front;
    d1_front = d0_front;
    d1_back = d0_front;
    for (int i = 0; i < vnbr; i++)
    {
        int id = i;
        int rfid = row_front(id);
        int rbid = row_back(id);
        int cfid = col_front(id);
        int cbid = col_back(id);
        if (rfid < 0 || rbid < 0 || cfid < 0 || cbid < 0) // if it is a boundary vertex, skip
        {
            continue;
        }
        int lu = row_front(cfid); // upleft
        int ld = row_back(cfid);
        int ru = row_front(cbid);
        int rd = row_back(cbid);
        d0_front[id] = lu;
        d0_back[id] = rd;
        d1_front[id] = ld;
        d1_back[id] = ru;
    }
}

void QuadOpt::init(CGMesh &mesh_in, const std::string &prefix)
{
    MP.mesh2Matrix(mesh_in, V, F);
    std::vector<std::vector<double>> rowstmp;
    std::vector<std::vector<double>> colstmp;
    read_csv_data_lbl(prefix + "_rows.csv", rowstmp);
    read_csv_data_lbl(prefix + "_cols.csv", colstmp);
    rowinfo = rowstmp;
    colinfo = colstmp;

    int vnbr = V.rows();
    if (OptType == 0)// AAG
    {
        // vertices vnbr * 3, normals vnbr * 3, binormals for G vnbr * 3  
        varsize = vnbr * 9;
    }
    if (OptType == 1)// AGG
    {
        // vertices vnbr * 3, normals vnbr * 3, binormals for G0 vnbr * 3, binormals for G1 vnbr * 3.    
        varsize = vnbr * 12;
    }
    if (OptType == 2)// PP
    {
        // vertices vnbr * 3, normals vnbr * 3, binormals vnbr * 6
        varsize = vnbr * 12;
   
    }
    if (OptType == 3)// PPG
    {
        // vertices vnbr * 3, normals vnbr * 3, binormals vnbr * 9
        varsize = vnbr * 15;

    }
    
    row_front = Eigen::VectorXi::Ones(vnbr) * -1;
    row_back = Eigen::VectorXi::Ones(vnbr) * -1;
    col_front = Eigen::VectorXi::Ones(vnbr) * -1;
    col_back = Eigen::VectorXi::Ones(vnbr) * -1;


    load_row_col_info(rowstmp, row_front,row_back);
    load_row_col_info(colstmp, col_front, col_back);

    load_diagonal_info(row_front, row_back, col_front, col_back, d0_front, d1_front, d0_back, d1_back);
    OrigVars = Eigen::VectorXd::Zero(varsize);
    OrigVars.segment(0, vnbr) = V.col(0);
    OrigVars.segment(vnbr, vnbr) = V.col(1);
    OrigVars.segment(vnbr * 2, vnbr) = V.col(2);
    mesh_original = mesh_in;
    mesh_update = mesh_original;
}

void QuadOpt::reset()
{
    // re-choose a diagonal
    if (OptType == 0) // AAG
    {                                  // set up G
        std::cout << "Set Up AAG Diagonals" << std::endl;
        if (WhichDiagonal == 0)
        {
            d0_type = 2;
            d1_type = 0;
        }
        if (WhichDiagonal == 1)
        { // set up G
            d0_type = 0;
            d1_type = 2;
        }
    }
    if (OptType == 1) // GGA
    {                                  // set up A
        std::cout << "Set Up GGA Diagonals" << std::endl;
        if (WhichDiagonal == 0)
        {
            d0_type = 1;
            d1_type = 0;
        }
        if (WhichDiagonal == 1)
        { // set up A
            d0_type = 0;
            d1_type = 1;
        }
    }
    MP.mesh2Matrix(mesh_original, V, F);
    int vnbr = V.rows();
    // reset the var size
    if (OptType == 0)// AAG
    {
        // vertices vnbr * 3, normals vnbr * 3, binormals for G vnbr * 3  
        varsize = vnbr * 9;
    }
    if (OptType == 1)// AGG
    {
        // vertices vnbr * 3, normals vnbr * 3, binormals for G0 vnbr * 3, binormals for G1 vnbr * 3.    
        varsize = vnbr * 12;
    }
    if (OptType == 2)// PP
    {
        // vertices vnbr * 3, normals vnbr * 3, binormals vnbr * 6
        varsize = vnbr * 12;
   
    }
    if (OptType == 3)// PPG
    {
        // vertices vnbr * 3, normals vnbr * 3, binormals vnbr * 9
        varsize = vnbr * 15;

    }

    OrigVars = Eigen::VectorXd::Zero(varsize);
    OrigVars.segment(0, vnbr) = V.col(0);
    OrigVars.segment(vnbr, vnbr) = V.col(1);
    OrigVars.segment(vnbr * 2, vnbr) = V.col(2);
    GlobVars.resize(0);
    ComputeAuxiliaries = true;
    mesh_update = mesh_original;
}

void QuadOpt::assemble_gravity(spMat& H, Eigen::VectorXd& B, Eigen::VectorXd &energy){

    // std::vector<Trip> tripletes;
    int vnbr = V.rows();
    H = spMat(Eigen::VectorXd::Ones(GlobVars.size()).asDiagonal());
    energy = Eigen::VectorXd::Zero(GlobVars.size());
    energy.segment(0, vnbr * 3) = (GlobVars - OrigVars).segment(0, vnbr * 3);

    B = -H * energy;
}
// lv is ver location, lf is front location, lb is back location, cid is the id of the condition
// ratiof is the ratio of lengths.
// v - (1-ratio) f - ratio * b = 0
void push_fairness_conditions(std::vector<Trip> &tripletes, Eigen::VectorXd &energy,
                              const int lv, const int lf, const int lb,
                              const double vval, const double fval, const double bval,
                              const int cid, const double ratiof)
{
    tripletes.push_back(Trip(cid, lv, 1));
    tripletes.push_back(Trip(cid, lf, ratiof - 1));
    tripletes.push_back(Trip(cid, lb, -ratiof));
    energy[cid] = vval - (1 - ratiof) * fval - ratiof * bval;
}

void QuadOpt::assemble_fairness(spMat& H, Eigen::VectorXd& B, Eigen::VectorXd &energy){
    

    std::vector<Trip> tripletes;
    int vnbr = V.rows();
    energy = Eigen::VectorXd::Zero(vnbr * 12);//todo
    tripletes.reserve(vnbr * 12 * 3);
    int counter = 0;
    for (int i = 0; i < vnbr; i++)
    {
        // the ids
        int vid = i;
        int rf = row_front[vid];
        int rb = row_back[vid];
        int cf = col_front[vid];
        int cb = col_back[vid];
        int d0f = d0_front[vid];
        int d0b = d0_back[vid];
        int d1f = d1_front[vid];
        int d1b = d1_back[vid];
        // the locations 
        int lvx = vid;
        int lvy = vid + vnbr;
        int lvz = vid + 2 * vnbr;

        int lrfx = rf;
        int lrfy = rf + vnbr;
        int lrfz = rf + vnbr * 2;

        int lrbx = rb;
        int lrby = rb + vnbr;
        int lrbz = rb + vnbr * 2;

        int lcfx = cf;
        int lcfy = cf + vnbr;
        int lcfz = cf + vnbr * 2;

        int lcbx = cb;
        int lcby = cb + vnbr;
        int lcbz = cb + vnbr * 2;

        int ld0fx = d0f;
        int ld0fy = d0f + vnbr;
        int ld0fz = d0f + vnbr * 2;

        int ld0bx = d0b;
        int ld0by = d0b + vnbr;
        int ld0bz = d0b + vnbr * 2;

        int ld1fx = d1f;
        int ld1fy = d1f + vnbr;
        int ld1fz = d1f + vnbr * 2;

        int ld1bx = d1b;
        int ld1by = d1b + vnbr;
        int ld1bz = d1b + vnbr * 2;
        
        bool row_smt = true;
        bool col_smt = true;
        bool d0_smt = true;
        bool d1_smt = true;

        if (rf < 0 || rb < 0)
        {
            row_smt = false;
        }
        if (cf < 0 || cb < 0)
        {
            col_smt = false;
        }

        if (d0_type == 0 || d0f < 0 || d0b < 0)
        {
            d0_smt = false;
        }
        if (d1_type == 0 || d1f < 0 || d1b < 0)
        {
            d1_smt = false;
        }
        if (row_smt)
        {
            assert(rf != rb);
            Eigen::Vector3d Ver(GlobVars[lvx], GlobVars[lvy], GlobVars[lvz]);
            Eigen::Vector3d Vfr(GlobVars[lrfx], GlobVars[lrfy], GlobVars[lrfz]);
            Eigen::Vector3d Vbk(GlobVars[lrbx], GlobVars[lrby], GlobVars[lrbz]);
            double ratiof = (Ver - Vfr).norm() / (Vbk - Vfr).norm();
            assert((Vbk - Vfr).norm() > 0);
            counter = i;
            push_fairness_conditions(tripletes, energy, lvx, lrfx, lrbx, Ver[0], Vfr[0], Vbk[0], counter, ratiof);
            counter = i + vnbr;
            push_fairness_conditions(tripletes, energy, lvy, lrfy, lrby, Ver[1], Vfr[1], Vbk[1], counter, ratiof);
            counter = i + vnbr * 2;
            push_fairness_conditions(tripletes, energy, lvz, lrfz, lrbz, Ver[2], Vfr[2], Vbk[2], counter, ratiof);
            
        }

        if (col_smt)
        {
            Eigen::Vector3d Ver(GlobVars[lvx], GlobVars[lvy], GlobVars[lvz]);
            Eigen::Vector3d Vfr(GlobVars[lcfx], GlobVars[lcfy], GlobVars[lcfz]);
            Eigen::Vector3d Vbk(GlobVars[lcbx], GlobVars[lcby], GlobVars[lcbz]);
            double ratiof = (Ver - Vfr).norm() / (Vbk - Vfr).norm();
            assert((Vbk - Vfr).norm() > 0);
            counter = i + vnbr * 3;
            push_fairness_conditions(tripletes, energy, lvx, lcfx, lcbx, Ver[0], Vfr[0], Vbk[0], counter, ratiof);
            counter = i + vnbr * 4;
            push_fairness_conditions(tripletes, energy, lvy, lcfy, lcby, Ver[1], Vfr[1], Vbk[1], counter, ratiof);
            counter = i + vnbr * 5;
            push_fairness_conditions(tripletes, energy, lvz, lcfz, lcbz, Ver[2], Vfr[2], Vbk[2], counter, ratiof);
            
        }

        if (d0_smt)
        {
            Eigen::Vector3d Ver(GlobVars[lvx], GlobVars[lvy], GlobVars[lvz]);
            Eigen::Vector3d Vfr(GlobVars[ld0fx], GlobVars[ld0fy], GlobVars[ld0fz]);
            Eigen::Vector3d Vbk(GlobVars[ld0bx], GlobVars[ld0by], GlobVars[ld0bz]);
            double ratiof = (Ver - Vfr).norm() / (Vbk - Vfr).norm();
            assert((Vbk - Vfr).norm() > 0);
            counter = i + vnbr * 6;
            push_fairness_conditions(tripletes, energy, lvx, ld0fx, ld0bx, Ver[0], Vfr[0], Vbk[0], counter, ratiof);
            counter = i + vnbr * 7;
            push_fairness_conditions(tripletes, energy, lvy, ld0fy, ld0by, Ver[1], Vfr[1], Vbk[1], counter, ratiof);
            counter = i + vnbr * 8;
            push_fairness_conditions(tripletes, energy, lvz, ld0fz, ld0bz, Ver[2], Vfr[2], Vbk[2], counter, ratiof);
        }
        if (d1_smt)
        {
            Eigen::Vector3d Ver(GlobVars[lvx], GlobVars[lvy], GlobVars[lvz]);
            Eigen::Vector3d Vfr(GlobVars[ld1fx], GlobVars[ld1fy], GlobVars[ld1fz]);
            Eigen::Vector3d Vbk(GlobVars[ld1bx], GlobVars[ld1by], GlobVars[ld1bz]);
            double ratiof = (Ver - Vfr).norm() / (Vbk - Vfr).norm();
            assert((Vbk - Vfr).norm() > 0);
            counter = i + vnbr * 9;
            push_fairness_conditions(tripletes, energy, lvx, ld1fx, ld1bx, Ver[0], Vfr[0], Vbk[0], counter, ratiof);
            counter = i + vnbr * 10;
            push_fairness_conditions(tripletes, energy, lvy, ld1fy, ld1by, Ver[1], Vfr[1], Vbk[1], counter, ratiof);
            counter = i + vnbr * 11;
            push_fairness_conditions(tripletes, energy, lvz, ld1fz, ld1bz, Ver[2], Vfr[2], Vbk[2], counter, ratiof);
        }
    }

    // std::cout<<"check 2"<<std::endl;
    spMat J;
    J.resize(energy.size(), GlobVars.size());
    J.setFromTriplets(tripletes.begin(), tripletes.end());
    H = J.transpose() * J;
    B = -J.transpose() * energy;
}

// bnm_start is the id where the binormal starts in the variable matrix
// from bnm_start to bnm_start + vnbr * 3 are the binormal variables of this family of curves
// family: 0: row. 1: col. 2: d0. 3: d1.
void QuadOpt::assemble_binormal_conditions(spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy, int family, int bnm_start)
{
    // std::cout<<"check 1, family "<<family<<std::endl;
    std::vector<Trip> tripletes;
    int vnbr = V.rows();
    energy = Eigen::VectorXd::Zero(vnbr * 3); 
    tripletes.reserve(vnbr * 25); 
    int counter = 0;
    for (int i = 0; i < vnbr; i++)
    {
        // the ids
        int vid = i;
        int rf = row_front[vid];
        int rb = row_back[vid];
        int cf = col_front[vid];
        int cb = col_back[vid];
        int d0f = d0_front[vid];
        int d0b = d0_back[vid];
        int d1f = d1_front[vid];
        int d1b = d1_back[vid];
        // the locations
        // the vertex
        int lvx = vid;
        int lvy = vid + vnbr;
        int lvz = vid + 2 * vnbr;
        // the front and back vertices
        int lfx, lfy, lfz, lbx, lby, lbz;
        // the binormal vector locations
        int lrx = bnm_start + vid;
        int lry = bnm_start + vid + vnbr;
        int lrz = bnm_start + vid + vnbr * 2;
        bool compute = true;
        if (family == 0) // rows
        {
            lfx = rf;
            lfy = rf + vnbr;
            lfz = rf + vnbr * 2;

            lbx = rb;
            lby = rb + vnbr;
            lbz = rb + vnbr * 2;
        }
        if (family == 1) // cols
        {
            lfx = cf;
            lfy = cf + vnbr;
            lfz = cf + vnbr * 2;

            lbx = cb;
            lby = cb + vnbr;
            lbz = cb + vnbr * 2;
        }
        if (family == 2) // d0
        {
            lfx = d0f;
            lfy = d0f + vnbr;
            lfz = d0f + vnbr * 2;

            lbx = d0b;
            lby = d0b + vnbr;
            lbz = d0b + vnbr * 2;
        }

        if (family == 3) // d1
        {

            lfx = d1f;
            lfy = d1f + vnbr;
            lfz = d1f + vnbr * 2;

            lbx = d1b;
            lby = d1b + vnbr;
            lbz = d1b + vnbr * 2;
        }
        if (lfx < 0 || lbx < 0)
        {
            compute = false;
        }
        
        if(compute == false){ // the vertex is on the boundary
            if(!ComputeAuxiliaries){
                // here we choose doing nothing and care only about inner vertices
            }
            continue;
        }
        Eigen::Vector3d Ver(GlobVars[lvx], GlobVars[lvy], GlobVars[lvz]);
        Eigen::Vector3d Vf(GlobVars[lfx], GlobVars[lfy], GlobVars[lfz]);
        Eigen::Vector3d Vb(GlobVars[lbx], GlobVars[lby], GlobVars[lbz]);
        if(ComputeAuxiliaries){
            //init the binormals
            Eigen::Vector3d real_r = ((Ver - Vf).cross(Vb - Ver)).normalized();
            if (real_r.norm() == 0)
            {
                real_r = Eigen::Vector3d(0, 0, 1);
            }
            GlobVars[lrx] = real_r[0];
            GlobVars[lry] = real_r[1];
            GlobVars[lrz] = real_r[2];
        }
        
        Eigen::Vector3d r(GlobVars[lrx], GlobVars[lry], GlobVars[lrz]);
        double dis0 = (Vf - Ver).norm();
        double dis1 = (Vb - Ver).norm();
        double dis = dis0 * dis1;
        Eigen::Vector3d fc = Vf - Ver;
        Eigen::Vector3d bc = Vb - Ver;
        // r * (Vf - Ver) = 0
        tripletes.push_back(Trip(i, lrx, (Vf(0) - Ver(0)) / dis0));
        tripletes.push_back(Trip(i, lry, (Vf(1) - Ver(1)) / dis0));
        tripletes.push_back(Trip(i, lrz, (Vf(2) - Ver(2)) / dis0));

        tripletes.push_back(Trip(i, lfx, (r(0)) / dis0));
        tripletes.push_back(Trip(i, lfy, (r(1)) / dis0));
        tripletes.push_back(Trip(i, lfz, (r(2)) / dis0));

        tripletes.push_back(Trip(i, lvx, -(r(0)) / dis0));
        tripletes.push_back(Trip(i, lvy, -(r(1)) / dis0));
        tripletes.push_back(Trip(i, lvz, -(r(2)) / dis0));

        energy[i] = r.dot(Vf - Ver) / dis0;
        // r * (Vb - Ver) = 0
        tripletes.push_back(Trip(i + vnbr, lrx, (Vb(0) - Ver(0)) / dis1));
        tripletes.push_back(Trip(i + vnbr, lry, (Vb(1) - Ver(1)) / dis1));
        tripletes.push_back(Trip(i + vnbr, lrz, (Vb(2) - Ver(2)) / dis1));

        tripletes.push_back(Trip(i + vnbr, lbx, (r(0)) / dis1));
        tripletes.push_back(Trip(i + vnbr, lby, (r(1)) / dis1));
        tripletes.push_back(Trip(i + vnbr, lbz, (r(2)) / dis1));

        tripletes.push_back(Trip(i + vnbr, lvx, -(r(0)) / dis1));
        tripletes.push_back(Trip(i + vnbr, lvy, -(r(1)) / dis1));
        tripletes.push_back(Trip(i + vnbr, lvz, -(r(2)) / dis1));

        energy[i + vnbr] = r.dot(Vb - Ver) / dis1;
        // r * r = 1
        tripletes.push_back(Trip(i + vnbr * 2, lrx, 2 * r(0)));
        tripletes.push_back(Trip(i + vnbr * 2, lry, 2 * r(1)));
        tripletes.push_back(Trip(i + vnbr * 2, lrz, 2 * r(2)));
        energy[i + vnbr * 2] = r.dot(r) - 1;
    }
    spMat J;
    J.resize(energy.size(), GlobVars.size());
    J.setFromTriplets(tripletes.begin(), tripletes.end());
    H = J.transpose() * J;
    B = -J.transpose() * energy;
}
// lv is ver location, lf is front location, lb is back location, cid is the id of the condition
void push_asymptotic_condition(std::vector<Trip> &tripletes, Eigen::VectorXd &energy, const int cid,
                        const int lvx, const int lvy, const int lvz,
                        const int lfx, const int lfy, const int lfz,
                        const int lbx, const int lby, const int lbz,
                        const int lnx, const int lny, const int lnz,
                        const Eigen::Vector3d &Vf, const Eigen::Vector3d &Ver, const Eigen::Vector3d &Vb, 
                        const Eigen::Vector3d &norm)
{

    double dis0 = (Vf - Ver).norm();
    double dis1 = (Vb - Ver).norm();
    // norm * (Vf-Ver) = 0
    tripletes.push_back(Trip(cid, lnx, (Vf(0) - Ver(0)) / dis0));
    tripletes.push_back(Trip(cid, lny, (Vf(1) - Ver(1)) / dis0));
    tripletes.push_back(Trip(cid, lnz, (Vf(2) - Ver(2)) / dis0));

    tripletes.push_back(Trip(cid, lfx, (norm(0)) / dis0));
    tripletes.push_back(Trip(cid, lfy, (norm(1)) / dis0));
    tripletes.push_back(Trip(cid, lfz, (norm(2)) / dis0));

    tripletes.push_back(Trip(cid, lvx, -(norm(0)) / dis0));
    tripletes.push_back(Trip(cid, lvy, -(norm(1)) / dis0));
    tripletes.push_back(Trip(cid, lvz, -(norm(2)) / dis0));

    energy[cid] = norm.dot(Vf - Ver) / dis0;

    // norm * (Vb-Ver) = 0
    tripletes.push_back(Trip(cid + 1, lnx, (Vb(0) - Ver(0)) / dis1));
    tripletes.push_back(Trip(cid + 1, lny, (Vb(1) - Ver(1)) / dis1));
    tripletes.push_back(Trip(cid + 1, lnz, (Vb(2) - Ver(2)) / dis1));

    tripletes.push_back(Trip(cid + 1, lbx, (norm(0)) / dis1));
    tripletes.push_back(Trip(cid + 1, lby, (norm(1)) / dis1));
    tripletes.push_back(Trip(cid + 1, lbz, (norm(2)) / dis1));

    tripletes.push_back(Trip(cid + 1, lvx, -(norm(0)) / dis1));
    tripletes.push_back(Trip(cid + 1, lvy, -(norm(1)) / dis1));
    tripletes.push_back(Trip(cid + 1, lvz, -(norm(2)) / dis1));

    energy[cid + 1] = norm.dot(Vb - Ver) / dis1;
    
}

void push_geodesic_condition(std::vector<Trip> &tripletes, Eigen::VectorXd &energy, const int cid,
                             const int lnx, const int lny, const int lnz,
                             const int lrx, const int lry, const int lrz,
                             const Eigen::Vector3d &norm, const Eigen::Vector3d &r)
{
    // norm * r = 0
    tripletes.push_back(Trip(cid, lnx, r(0)));
    tripletes.push_back(Trip(cid, lny, r(1)));
    tripletes.push_back(Trip(cid, lnz, r(2)));

    tripletes.push_back(Trip(cid, lrx, norm(0)));
    tripletes.push_back(Trip(cid, lry, norm(1)));
    tripletes.push_back(Trip(cid, lrz, norm(2)));

    energy[cid] = norm.dot(r);
}

// type: 0 disabled. 1 means asymptotic, 2 means geodesic, 3 pseudo-geodesic
// aux_start_location is the start location of the auxiliaries of this family
void QuadOpt::assemble_pg_extreme_cases(spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy,
                          const int type, const int family, const int bnm_start_location)
{

    std::vector<Trip> tripletes;
    int vnbr = V.rows();
    int energy_size = 0;
    if(type == 1){ // asymptotic
        energy_size = vnbr * 2;
    }
    if(type == 2){
        energy_size = vnbr;
    }
    if(type == 3){
        std::cout<<"ERROR: Please do not use this function for general pseudo-geodesic solving"<<std::endl;
        return;
    }
    energy = Eigen::VectorXd::Zero(energy_size);
    tripletes.reserve(vnbr * 18);
    int counter = 0;
    for (int i = 0; i < vnbr; i++)
    {
        // the ids
        int vid = i;
        int rf = row_front[vid];
        int rb = row_back[vid];
        int cf = col_front[vid];
        int cb = col_back[vid];
        int d0f = d0_front[vid];
        int d0b = d0_back[vid];
        int d1f = d1_front[vid];
        int d1b = d1_back[vid];

         // the locations
        // the vertex
        int lvx = vid;
        int lvy = vid + vnbr;
        int lvz = vid + 2 * vnbr;
        // the front and back vertices
        int lfx, lfy, lfz, lbx, lby, lbz;
        int lnx = vnbr * 3 + i;
        int lny = vnbr * 4 + i;
        int lnz = vnbr * 5 + i;
        bool compute = true;
        if (family == 0) // rows
        {
            lfx = rf;
            lfy = rf + vnbr;
            lfz = rf + vnbr * 2;

            lbx = rb;
            lby = rb + vnbr;
            lbz = rb + vnbr * 2;
        }
        if (family == 1) // cols
        {
            lfx = cf;
            lfy = cf + vnbr;
            lfz = cf + vnbr * 2;

            lbx = cb;
            lby = cb + vnbr;
            lbz = cb + vnbr * 2;
        }
        if (family == 2) // d0
        {
            lfx = d0f;
            lfy = d0f + vnbr;
            lfz = d0f + vnbr * 2;

            lbx = d0b;
            lby = d0b + vnbr;
            lbz = d0b + vnbr * 2;
        }

        if (family == 3) // d1
        {

            lfx = d1f;
            lfy = d1f + vnbr;
            lfz = d1f + vnbr * 2;

            lbx = d1b;
            lby = d1b + vnbr;
            lbz = d1b + vnbr * 2;
        }
        if (lfx < 0 || lbx < 0)
        {
            compute = false;
        }
        if (!compute)
        {
            continue;
        }
        Eigen::Vector3d Ver(GlobVars[lvx], GlobVars[lvy], GlobVars[lvz]);
        Eigen::Vector3d Vfr(GlobVars[lfx], GlobVars[lfy], GlobVars[lfz]);
        Eigen::Vector3d Vbk(GlobVars[lbx], GlobVars[lby], GlobVars[lbz]);
        Eigen::Vector3d norm(GlobVars[lnx], GlobVars[lny], GlobVars[lnz]);
        if (type == 1)
        { // asymptotic
            counter = i * 2;
            push_asymptotic_condition(tripletes, energy, counter, lvx, lvy, lvz, lfx, lfy, lfz, lbx, lby, lbz, lnx, lny, lnz,
                                      Vfr, Ver, Vbk, norm);
        }
        if(type == 2){// geodesic
            counter = i;
            int lrx = bnm_start_location + i;
            int lry = bnm_start_location + i + vnbr;
            int lrz = bnm_start_location + i + vnbr * 2;
            Eigen::Vector3d r(GlobVars[lrx], GlobVars[lry], GlobVars[lrz]);
            push_geodesic_condition(tripletes, energy, counter, lnx, lny, lnz, lrx, lry, lrz, norm, r);
        }
    }

    // std::cout<<"check 2"<<std::endl;
    spMat J;
    J.resize(energy.size(), GlobVars.size());
    J.setFromTriplets(tripletes.begin(), tripletes.end());
    H = J.transpose() * J;
    B = -J.transpose() * energy;
}
void QuadOpt::assemble_pg_cases(const double angle_radian, spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy,
                                const int family, const int bnm_start_location)
{
    std::vector<Trip> tripletes;
    int vnbr = V.rows();
    int energy_size = 0;

    energy_size = vnbr * 2;

    energy = Eigen::VectorXd::Zero(energy_size);
    tripletes.reserve(vnbr * 6);
    for (int i = 0; i < vnbr; i++)
    {
        // the ids
        int vid = i;
        int rf = row_front[vid];
        int rb = row_back[vid];
        int cf = col_front[vid];
        int cb = col_back[vid];
        int d0f = d0_front[vid];
        int d0b = d0_back[vid];
        int d1f = d1_front[vid];
        int d1b = d1_back[vid];

         // the locations
        // the vertex
        int lvx = vid;
        int lvy = vid + vnbr;
        int lvz = vid + 2 * vnbr;
        // the front and back vertices
        int lfx, lfy, lfz, lbx, lby, lbz;
        int lnx = vnbr * 3 + i;
        int lny = vnbr * 4 + i;
        int lnz = vnbr * 5 + i;
        int lrx = bnm_start_location + i;
        int lry = bnm_start_location + i + vnbr;
        int lrz = bnm_start_location + i + vnbr * 2;
        
        bool compute = true;
        if (family == 0) // rows
        {
            lfx = rf;
            lfy = rf + vnbr;
            lfz = rf + vnbr * 2;

            lbx = rb;
            lby = rb + vnbr;
            lbz = rb + vnbr * 2;
        }
        if (family == 1) // cols
        {
            lfx = cf;
            lfy = cf + vnbr;
            lfz = cf + vnbr * 2;

            lbx = cb;
            lby = cb + vnbr;
            lbz = cb + vnbr * 2;
        }
        if (family == 2) // d0
        {
            lfx = d0f;
            lfy = d0f + vnbr;
            lfz = d0f + vnbr * 2;

            lbx = d0b;
            lby = d0b + vnbr;
            lbz = d0b + vnbr * 2;
        }

        if (family == 3) // d1
        {

            lfx = d1f;
            lfy = d1f + vnbr;
            lfz = d1f + vnbr * 2;

            lbx = d1b;
            lby = d1b + vnbr;
            lbz = d1b + vnbr * 2;
        }
        if (lfx < 0 || lbx < 0)
        {
            compute = false;
        }
        if (!compute)
        {
            continue;
        }
        Eigen::Vector3d Ver(GlobVars[lvx], GlobVars[lvy], GlobVars[lvz]);
        Eigen::Vector3d Vfr(GlobVars[lfx], GlobVars[lfy], GlobVars[lfz]);
        Eigen::Vector3d Vbk(GlobVars[lbx], GlobVars[lby], GlobVars[lbz]);
        Eigen::Vector3d norm(GlobVars[lnx], GlobVars[lny], GlobVars[lnz]);
        Eigen::Vector3d r(GlobVars[lrx], GlobVars[lry], GlobVars[lrz]);

        Eigen::Vector3d tangent = -(Vbk -Vfr);// This sign is changable depends on what user defines of the curve tangent directions
        
        Eigen::Vector3d real_u = (norm.cross(tangent)).normalized();

        double dis0 = (Ver - Vfr).norm();
        double dis1 = (Ver - Vbk).norm();
        double cos_angle = cos(angle_radian);
        double sin_angle = sin(angle_radian);

        // (binormal * n)^2 = cos^2
        double ndb = norm.dot(r);
        tripletes.push_back(Trip(i, lnx, 2 * ndb * norm(0)));
        tripletes.push_back(Trip(i, lny, 2 * ndb * norm(1)));
        tripletes.push_back(Trip(i, lnz, 2 * ndb * norm(2)));

        energy[i] = ndb * ndb - cos_angle * cos_angle;

        // (binormal * n) * (binormal * u) = sin * cos
        // regard n and u as constant values in each iteration, since u is related to the tangent vector, which will not change much.
        double udb = real_u.dot(r);

        tripletes.push_back(Trip(i + vnbr, lnx, ndb * real_u[0] + udb * norm[0]));
        tripletes.push_back(Trip(i + vnbr, lny, ndb * real_u[1] + udb * norm[1]));
        tripletes.push_back(Trip(i + vnbr, lnz, ndb * real_u[2] + udb * norm[2]));

        energy[i + vnbr] = ndb * udb - sin_angle * cos_angle;
        
    }

    // std::cout<<"check 2"<<std::endl;
    spMat J;
    J.resize(energy.size(), GlobVars.size());
    J.setFromTriplets(tripletes.begin(), tripletes.end());
    H = J.transpose() * J;
    B = -J.transpose() * energy;
}

void QuadOpt::assemble_normal_conditions(spMat& H, Eigen::VectorXd& B, Eigen::VectorXd &energy){
    std::vector<Trip> tripletes;
    int vnbr = V.rows();
    energy = Eigen::VectorXd::Zero(vnbr * 4); // todo
    tripletes.reserve(vnbr * 12 * 3); // todo
    int counter = 0;
    for (int i = 0; i < vnbr; i++)
    {
        // the ids
        int vid = i;
        int rf = row_front[vid];
        int rb = row_back[vid];
        int cf = col_front[vid];
        int cb = col_back[vid];
        int d0f = d0_front[vid];
        int d0b = d0_back[vid];
        int d1f = d1_front[vid];
        int d1b = d1_back[vid];
        // the locations
        // the vertex
        int lvx = vid;
        int lvy = vid + vnbr;
        int lvz = vid + 2 * vnbr;
        
        int lrfx = rf;
        int lrfy = rf + vnbr;
        int lrfz = rf + vnbr * 2;

        int lrbx = rb;
        int lrby = rb + vnbr;
        int lrbz = rb + vnbr * 2;

        int lcfx = cf;
        int lcfy = cf + vnbr;
        int lcfz = cf + vnbr * 2;

        int lcbx = cb;
        int lcby = cb + vnbr;
        int lcbz = cb + vnbr * 2;

        int lnx = i + vnbr * 3;
        int lny = i + vnbr * 4;
        int lnz = i + vnbr * 5;

        bool row_smt = true;
        bool col_smt = true;
        
        if (rf < 0 || rb < 0)
        {
            row_smt = false;
        }
        if (cf < 0 || cb < 0)
        {
            col_smt = false;
        }
        if (row_smt == false || col_smt == false)
        { // the vertex is on the boundary
            if(!ComputeAuxiliaries){
                // here we choose doing nothing and care only about inner vertices
            }
            continue;
        }
        
        Eigen::Vector3d Ver(GlobVars[lvx], GlobVars[lvy], GlobVars[lvz]);
        Eigen::Vector3d Vrf(GlobVars[lrfx], GlobVars[lrfy], GlobVars[lrfz]);
        Eigen::Vector3d Vrb(GlobVars[lrbx], GlobVars[lrby], GlobVars[lrbz]);
        Eigen::Vector3d Vcf(GlobVars[lcfx], GlobVars[lcfy], GlobVars[lcfz]);
        Eigen::Vector3d Vcb(GlobVars[lcbx], GlobVars[lcby], GlobVars[lcbz]);
        if(ComputeAuxiliaries){
            //init the binormals
            Eigen::Vector3d real_n = ((Vrf - Vrb).cross(Vcf - Vcb)).normalized();
            GlobVars[lnx] = real_n[0];
            GlobVars[lny] = real_n[1];
            GlobVars[lnz] = real_n[2];
        }
        
        
        Eigen::Vector3d norm(GlobVars[lnx], GlobVars[lny], GlobVars[lnz]);
        double dis0 = (Vrf - Ver).norm();
        double dis1 = (Vrb - Ver).norm();
        double dis2 = (Vcf - Ver).norm();
        double dis3 = (Vcb - Ver).norm();

        // n * (Vrf - Ver) = 0
        tripletes.push_back(Trip(i, lnx, (Vrf(0) - Ver(0)) / dis0));
        tripletes.push_back(Trip(i, lny, (Vrf(1) - Ver(1)) / dis0));
        tripletes.push_back(Trip(i, lnz, (Vrf(2) - Ver(2)) / dis0));

        tripletes.push_back(Trip(i, lrfx, (norm(0)) / dis0));
        tripletes.push_back(Trip(i, lrfy, (norm(1)) / dis0));
        tripletes.push_back(Trip(i, lrfz, (norm(2)) / dis0));

        tripletes.push_back(Trip(i, lvx, -(norm(0)) / dis0));
        tripletes.push_back(Trip(i, lvy, -(norm(1)) / dis0));
        tripletes.push_back(Trip(i, lvz, -(norm(2)) / dis0));

        energy[i] = norm.dot(Vrf - Ver) / dis0;

        // norm * (Vrb - Ver) = 0
        tripletes.push_back(Trip(i + vnbr, lnx, (Vrb(0) - Ver(0)) / dis1));
        tripletes.push_back(Trip(i + vnbr, lny, (Vrb(1) - Ver(1)) / dis1));
        tripletes.push_back(Trip(i + vnbr, lnz, (Vrb(2) - Ver(2)) / dis1));

        tripletes.push_back(Trip(i + vnbr, lrbx, (norm(0)) / dis1));
        tripletes.push_back(Trip(i + vnbr, lrby, (norm(1)) / dis1));
        tripletes.push_back(Trip(i + vnbr, lrbz, (norm(2)) / dis1));

        tripletes.push_back(Trip(i + vnbr, lvx, -(norm(0)) / dis1));
        tripletes.push_back(Trip(i + vnbr, lvy, -(norm(1)) / dis1));
        tripletes.push_back(Trip(i + vnbr, lvz, -(norm(2)) / dis1));

        energy[i + vnbr] = norm.dot(Vrb - Ver) / dis1;

        // norm * (Vcf - Ver) = 0
        tripletes.push_back(Trip(i + vnbr * 2, lnx, (Vcf(0) - Ver(0)) / dis2));
        tripletes.push_back(Trip(i + vnbr * 2, lny, (Vcf(1) - Ver(1)) / dis2));
        tripletes.push_back(Trip(i + vnbr * 2, lnz, (Vcf(2) - Ver(2)) / dis2));

        tripletes.push_back(Trip(i + vnbr * 2, lcfx, (norm(0)) / dis2));
        tripletes.push_back(Trip(i + vnbr * 2, lcfy, (norm(1)) / dis2));
        tripletes.push_back(Trip(i + vnbr * 2, lcfz, (norm(2)) / dis2));

        tripletes.push_back(Trip(i + vnbr * 2, lvx, -(norm(0)) / dis2));
        tripletes.push_back(Trip(i + vnbr * 2, lvy, -(norm(1)) / dis2));
        tripletes.push_back(Trip(i + vnbr * 2, lvz, -(norm(2)) / dis2));

        energy[i + vnbr * 2] = norm.dot(Vcf - Ver) / dis2;

        // norm * (Vcb - Ver) = 0
        tripletes.push_back(Trip(i + vnbr * 3, lnx, (Vcb(0) - Ver(0)) / dis3));
        tripletes.push_back(Trip(i + vnbr * 3, lny, (Vcb(1) - Ver(1)) / dis3));
        tripletes.push_back(Trip(i + vnbr * 3, lnz, (Vcb(2) - Ver(2)) / dis3));

        tripletes.push_back(Trip(i + vnbr * 3, lcbx, (norm(0)) / dis3));
        tripletes.push_back(Trip(i + vnbr * 3, lcby, (norm(1)) / dis3));
        tripletes.push_back(Trip(i + vnbr * 3, lcbz, (norm(2)) / dis3));

        tripletes.push_back(Trip(i + vnbr * 3, lvx, -(norm(0)) / dis3));
        tripletes.push_back(Trip(i + vnbr * 3, lvy, -(norm(1)) / dis3));
        tripletes.push_back(Trip(i + vnbr * 3, lvz, -(norm(2)) / dis3));

        energy[i + vnbr * 3] = norm.dot(Vcb - Ver) / dis3;
    }
    spMat J;
    J.resize(energy.size(), GlobVars.size());
    J.setFromTriplets(tripletes.begin(), tripletes.end());
    H = J.transpose() * J;
    B = -J.transpose() * energy;
}

void QuadOpt::opt(){
    int vnbr = V.rows();
    if (GlobVars.size() != varsize)
    {
        ComputeAuxiliaries = true;
        GlobVars = Eigen::VectorXd::Zero(varsize);
        GlobVars.segment(0, vnbr) = V.col(0);
        GlobVars.segment(vnbr, vnbr) = V.col(1);
        GlobVars.segment(vnbr * 2, vnbr) = V.col(2);
    }
    spMat H;
    H.resize(varsize, varsize);
    Eigen::VectorXd B = Eigen::VectorXd::Zero(varsize);
    Eigen::VectorXd Egravity, Esmth, Enorm, Ebnm0, Ebnm1, Ebnm2, Epg[3];

    spMat Hgravity;  // approximation
	Eigen::VectorXd Bgravity; // right of laplacian
	assemble_gravity(Hgravity, Bgravity, Egravity);
    H += weight_gravity * Hgravity;
	B += weight_gravity * Bgravity;

    spMat Hsmth;
    Eigen::VectorXd Bsmth;
	assemble_fairness(Hsmth, Bsmth, Esmth);
    H += weight_fairness * Hsmth;
	B += weight_fairness * Bsmth;

    spMat Hnorm;
    Eigen::VectorXd Bnorm;
    assemble_normal_conditions(Hnorm, Bnorm, Enorm);
    H += weight_pg * Hnorm;
	B += weight_pg * Bnorm;

    if (OptType == 0)
    { // AAG, the diagonal is geodesic that requires binormals
        spMat Hbnm;
        Eigen::VectorXd Bbnm;
        int family = -1;
        if (d0_type > 0)
        {
            family = 2;
        }
        if (d1_type > 0)
        {
            family = 3;
        }
        // vertices vnbr * 3, normals vnbr * 3, binormals for G vnbr * 3  
        int bnm_start = vnbr * 6;
        assemble_binormal_conditions(Hbnm, Bbnm, Ebnm0, family, bnm_start);
        H += weight_pg * Hbnm;
        B += weight_pg * Bbnm;


        spMat Hpg[3];
        Eigen::VectorXd Bpg[3];
        // G
        int type = 2; // G
        assemble_pg_extreme_cases(Hpg[0], Bpg[0], Epg[0], type, family, bnm_start);

        // A
        type = 1; // A0
        family = 0;// row
        bnm_start = -1;
        assemble_pg_extreme_cases(Hpg[1], Bpg[1], Epg[1], type, family, bnm_start);

        type = 1; // A1
        family = 1;// col
        bnm_start = -1;
        assemble_pg_extreme_cases(Hpg[2], Bpg[2], Epg[2], type, family, bnm_start);

        H += weight_pg * (Hpg[0] + Hpg[1] + Hpg[2]);
        B += weight_pg * (Bpg[0] + Bpg[1] + Bpg[2]);
    }

    if (OptType == 1) 
    { // GGA, the diagonal is asymptotic
        spMat Hbnm[2];
        Eigen::VectorXd Bbnm[2];
        int family = -1;
        family = 0;// row
        // vertices vnbr * 3, normals vnbr * 3, binormals for G0 vnbr * 3, binormals for G1 vnbr * 3    
        int bnm_start = vnbr * 6;
        assemble_binormal_conditions(Hbnm[0], Bbnm[0], Ebnm0, family, bnm_start);

        family = 1; // column
        bnm_start = vnbr * 9;
        assemble_binormal_conditions(Hbnm[1], Bbnm[1], Ebnm1, family, bnm_start);

        H += weight_pg * (Hbnm[0] + Hbnm[1]);
        B += weight_pg * (Bbnm[0] + Bbnm[1]);
        
        if (d0_type > 0)
        {
            family = 2;
        }
        if (d1_type > 0)
        {
            family = 3;
        }

        spMat Hpg[3];
        Eigen::VectorXd Bpg[3];
        // A
        int type = 1; // A
        bnm_start = -1;
        assemble_pg_extreme_cases(Hpg[0], Bpg[0], Epg[0], type, family, bnm_start);

        // G
        type = 2; // G0
        family = 0;// row
        bnm_start = vnbr * 6;
        assemble_pg_extreme_cases(Hpg[1], Bpg[1], Epg[1], type, family, bnm_start);

        type = 2; // G1
        family = 1;// col
        bnm_start = vnbr * 9;
        assemble_pg_extreme_cases(Hpg[2], Bpg[2], Epg[2], type, family, bnm_start);

        H += weight_pg * (Hpg[0] + Hpg[1] + Hpg[2]);
        B += weight_pg * (Bpg[0] + Bpg[1] + Bpg[2]);
    }
    if (OptType == 2 || OptType == 3) 
    { // PP or PPG
        spMat Hbnm[3];
        Eigen::VectorXd Bbnm[3];
        int family = -1;
        family = 0;// row
        // vertices vnbr * 3, normals vnbr * 3, binormals for P0 vnbr * 3, binormals for P1 vnbr * 3    
        int bnm_start = vnbr * 6;
        assemble_binormal_conditions(Hbnm[0], Bbnm[0], Ebnm0, family, bnm_start);

        family = 1; // column
        bnm_start = vnbr * 9;
        assemble_binormal_conditions(Hbnm[1], Bbnm[1], Ebnm1, family, bnm_start);

        H += weight_pg * (Hbnm[0] + Hbnm[1]);
        B += weight_pg * (Bbnm[0] + Bbnm[1]);


        spMat Hpg[3];
        Eigen::VectorXd Bpg[3];
        // P
        bnm_start = vnbr * 6;
        family = 0;
        double angle_radian = angle_degree0 * LSC_PI / 180.;
        assemble_pg_cases(angle_radian, Hpg[0], Bpg[0], Epg[0], family, bnm_start);

        // P
        bnm_start = vnbr * 9;
        family = 1;
        angle_radian = angle_degree1 * LSC_PI / 180.;
        assemble_pg_cases(angle_radian, Hpg[0], Bpg[0], Epg[0], family, bnm_start);
        
        H += weight_pg * (Hpg[0] + Hpg[1]);
        B += weight_pg * (Bpg[0] + Bpg[1]);
        // type = 2; // G1
        // family = 1;// col
        // bnm_start = vnbr * 9;
        // assemble_pg_extreme_cases(Hpg[2], Bpg[2], Epg[2], type, family, bnm_start);


        
        if (OptType == 3)// G as diagonal 
        {
            if (WhichDiagonal == 0)
            {
                family = 2;
            }
            if (WhichDiagonal == 1)
            {
                family = 3;
            }

            bnm_start = vnbr * 12;
            assemble_binormal_conditions(Hbnm[2], Bbnm[2], Ebnm2, family, bnm_start);
            H += weight_pg * Hbnm[2];
            B += weight_pg * Bbnm[2];

            int type = 2; 
            assemble_pg_extreme_cases(Hpg[2], Bpg[2], Epg[2], type, family, bnm_start);
            H += weight_pg * Hpg[2];
            B += weight_pg * Bpg[2];
        }
        
    }

    // assemble together
    H += 1e-6 * spMat(Eigen::VectorXd::Ones(varsize).asDiagonal());
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
    dx *= 0.75;
    GlobVars += dx;
    if(OptType == 0){
        std::cout<<"AAG, ";
    }
    if(OptType == 1){
        std::cout<<"GGA, ";
    }
    if(OptType == 2){
        std::cout<<"PP, ";
    }
    if(OptType == 3){
        std::cout<<"PPG, ";
    }

    double ev_appro = Egravity.norm();
    double ev_smt = Esmth.norm();
    double ev_norm = Enorm.norm();
    std::cout << "Eclose, " << ev_appro << ", lap, " << ev_smt << ", normal vector, " << ev_norm << ", ";
    if (OptType == 0)
    {
        double ebi = Ebnm0.norm();

        std::cout << "bnm, " << ebi << ", GAA, " << Epg[0].norm() << ", " << Epg[1].norm() << ", " << Epg[2].norm() << ", ";
    }
    if (OptType == 1)
    {
        double ebi0 = Ebnm0.norm();
        double ebi1 = Ebnm1.norm();
        std::cout << "bnm, " << ebi0 << ", " << ebi1 << ", AGG, " << Epg[0].norm() << ", " << Epg[1].norm() << ", " << Epg[2].norm() << ", ";
    }
    if (OptType == 2 || OptType == 3)
    {
        double ebi0 = Ebnm0.norm();
        double ebi1 = Ebnm1.norm();

        std::cout << "bnm, " << ebi0 << ", " << ebi1;
        if (OptType == 3)
        {
            double ebi2 = Ebnm2.norm();
            std::cout << ", " << ebi2 << ", PPG, ";
        }
        else
        {
            std::cout << ", PP, ";
        }
        std::cout << Epg[0].norm() << ", " << Epg[1].norm();
        if (OptType == 3)
        {
            std::cout << ", " <<Epg[2].norm();
        }
    }
    real_step_length = dx.norm();
    std::cout << ", stp, " << dx.norm() << ", ";
    std::cout << "\n";

    // convert the data to the mesh format
    ComputeAuxiliaries = false;
    V.col(0) = GlobVars.segment(0, vnbr);
    V.col(1) = GlobVars.segment(vnbr, vnbr);
    V.col(2) = GlobVars.segment(vnbr * 2, vnbr);
    for (CGMesh::VertexIter v_it = mesh_update.vertices_begin(); v_it != mesh_update.vertices_end(); ++v_it)
    {
        int vid = v_it.handle().idx();
        mesh_update.point(*v_it) = CGMesh::Point(V(vid, 0), V(vid, 1), V(vid, 2));
    }
}

void QuadOpt::show_curve_families(std::array<Eigen::MatrixXd, 3>& edges){
    int row_m = rowinfo.size() / 2;
    int col_m = rowinfo[row_m].size() / 2;
    int vid = rowinfo[row_m][col_m]; // choose a point in the middle

    std::vector<Eigen::Vector3d> vtp;
    int vcurrent = vid;
    for (int i = 0;; i++)// rows
    {
        if (row_back[vcurrent] < 0)
        {
            break;
        }
        int vnext = row_back[vcurrent];
        vtp.push_back(V.row(vnext));
        vcurrent = vnext;
    }
    edges[0] = vec_list_to_matrix(vtp);

    vtp.clear();
    vcurrent = vid;
    for (int i = 0;; i++)// cols
    {
        if (col_back[vcurrent] < 0)
        {
            break;
        }
        int vnext = col_back[vcurrent];
        vtp.push_back(V.row(vnext));
        vcurrent = vnext;
    }
    edges[1] = vec_list_to_matrix(vtp);

    vtp.clear();
    vcurrent = vid;
    Eigen::VectorXi ajc;
    if (d0_type > 0)
    {
        ajc = d0_back;
    }
    if (d1_type > 0)
    {
        ajc = d1_back;
    }
    
    for (int i = 0;; i++) // diagonals
    {
        if (ajc[vcurrent] < 0)
        {
            break;
        }
        int vnext = ajc[vcurrent];
        vtp.push_back(V.row(vnext));
        vcurrent = vnext;
    }
    edges[2] = vec_list_to_matrix(vtp);
    

}