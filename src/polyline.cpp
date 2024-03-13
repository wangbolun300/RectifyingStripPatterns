#include <lsc/basic.h>
#include <lsc/tools.h>
#include<igl/file_dialog_open.h>
#include<igl/file_dialog_save.h>
#include <igl/point_mesh_squared_distance.h>
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

// the nbr of vers - 1 = the nbr of length = the nbr of angles + 1.
void check_strip_straightness(const std::vector<double>& lengths, const std::vector<double> &angles, double& devratio, bool print){
    // set the initial direction as (0, 0) ->(0, 1)
    int vnbr = lengths.size() + 1;
    Eigen::MatrixXd pts;
    double total_length = 0;
    pts.resize(vnbr, 2);
    pts.row(0) = Eigen::Vector2d(0, 0);
    pts.row(1) = Eigen::Vector2d(0, lengths[0]);
    total_length+=lengths[0];
    for (int i = 2; i < vnbr; i++)
    {
        // solve pts[i]
        // the reference direction from i-2 to i-1,
        Eigen::Vector2d ref_dirc = (pts.row(i - 1) - pts.row(i - 2)).normalized();
        double ang_diff = LSC_PI - angles[i - 2]; // counter clock rotation
        double len = lengths[i - 1];
        total_length+=len;
        double a = cos(ang_diff);
        double b = sin(ang_diff);
        Eigen::Vector2d ref_otho(-ref_dirc[1], ref_dirc[0]); // the othogonal vector is (-y, x)
        // new direction is a * ref_dirc + b * ref_otho.
        Eigen::Vector2d direction = a * ref_dirc + b * ref_otho;
        Eigen::Vector2d newpt = Eigen::Vector2d(pts.row(i - 1)) + direction * len;
        pts.row(i) = newpt;
    }
    Eigen::VectorXd iden = Eigen::VectorXd::Ones(vnbr);
    Eigen::Vector2d centroid = pts.transpose() * iden;
    centroid /= vnbr;
    Eigen::MatrixXd y = pts.rowwise() - centroid.transpose();// y is a matrix of n x 2
    Eigen::Matrix2d YTY = y.transpose() * y;// 2 x 2
    Eigen::EigenSolver<Eigen::MatrixXd> eg(YTY);
    Eigen::VectorXcd ev = eg.eigenvalues();
    int minevid = ev[0].real() < ev[1].real() ? 0 : 1;
    Eigen::VectorXcd eigenvecC = eg.eigenvectors().col(minevid);
    Eigen::Vector2d eigenvec(eigenvecC[0].real(), eigenvecC[1].real());
    Eigen::Vector2d n = eigenvec.normalized();// this is othogonal to the regression line.

    Eigen::VectorXd proj_dis = (y * n).rowwise().norm();
    if (print)
    {
        std::cout << "The maximal derivate: " << proj_dis.maxCoeff() << ", the total length of this polyline: " << total_length << std::endl;
    }

    devratio = proj_dis.maxCoeff() / total_length;
}


// the tangents are on each vertex. 
// evaluate the rectifying planes intersects with each other
void evaluate_intersect_rectifying(const std::vector<Eigen::Vector3d> &vertices, const std::vector<Eigen::Vector3d> &tangents,
                                   const std::vector<Eigen::Vector3d> &binormals, const std::vector<double>& lengths)
{
    int vnbr = vertices.size();
    double angle_diff = 0;
    std::vector<double> angles;
    angles.resize(vnbr - 1);
    
    for (int i = 0; i < vnbr - 1; i++)
    {
        Eigen::Vector3d normal1 = (tangents[i].cross(binormals[i])).normalized();
        Eigen::Vector3d normal2 = (tangents[i + 1].cross(binormals[i + 1])).normalized();
        Eigen::Vector3d crease;
        if (normal1.cross(normal2).norm() < 1e-6)
        { // if the two planes are already coplanar, then the crease is the binormal
            crease = binormals[i];
        }
        else
        {
            crease = normal1.cross(normal2).normalized();
        }
        if (crease.dot(binormals[i]) < 0)
        { // make sure the crease is the same orientation as the binormals.
            crease *= -1;
        }
        double ang_left = acos(crease.dot(-tangents[i]));
        double ang_right = acos(crease.dot(tangents[i + 1]));
        angles[i] = ang_left + ang_right;
    }
    double ratio;
    check_strip_straightness(lengths, angles, ratio, true);

}
void evaluate_strip_straightness(const std::vector<std::vector<Eigen::Vector3d>>& plys,const std::vector<std::vector<Eigen::Vector3d>>& bnms){
    double max_ratio = 0;
    for(int i=0;i<plys.size();i++){
        int pnbr = plys[i].size();
        std::vector<double> lengths(pnbr - 1);
        std::vector<double> angles(pnbr - 2);
        for (int j = 1; j < pnbr - 1; j++)
        {
            Eigen::Vector3d vs = plys[i][j-1];
            Eigen::Vector3d vm = plys[i][j];
            Eigen::Vector3d ve = plys[i][j + 1];
            lengths[j - 1] = (vs - vm).norm();
            Eigen::Vector3d crease = bnms[i][j];
            Eigen::Vector3d mtangent = (vs - vm).normalized();
            Eigen::Vector3d tangent = (ve - vm).normalized();
            

            double ang_left = acos(crease.dot(mtangent));
            double ang_right = acos(crease.dot(tangent));
            angles[j - 1] = ang_left + ang_right;
        }
        lengths.back() = (plys[i][pnbr - 1] - plys[i][pnbr - 2]).norm();
        double ratio;
        check_strip_straightness(lengths, angles, ratio, false);
        if (ratio > max_ratio)
        {
            max_ratio = ratio;
        }
    }
    std::cout << "Evaluating the straightness: maximal deviation is " << max_ratio << "\n";
}

void get_polyline_rectifying_plane_on_inner_vers(const std::vector<Eigen::Vector3d> &ply_in, const std::vector<Eigen::Vector3d> &bnm_in,
                                                 std::vector<Eigen::Vector3d> &vertices, std::vector<Eigen::Vector3d> &tangents, std::vector<Eigen::Vector3d> &binormals)
{
    int vnbr = ply_in.size()-2;
    vertices.resize(vnbr);
    tangents.resize(vnbr);
    binormals.resize(vnbr);
    std::vector<double> lengths;// the lengths of each segment estimated on each inner vertex.
    lengths.resize(vnbr);
    for (int i = 0; i < vnbr; i++)
    {
        // Eigen::Vector3d dir0 = (ply_in[i + 1] - ply_in[i]).normalized();
        // Eigen::Vector3d dir1 = (ply_in[i + 2] - ply_in[i + 1]).normalized();
        // Eigen::Vector3d bisector = (dir0 + dir1).normalized();
        Eigen::Vector3d bisector = (ply_in[i + 2] - ply_in[i]).normalized();
        vertices[i] = ply_in[i + 1];
        binormals[i] = bnm_in[i + 1];
        tangents[i] = bisector;
        lengths[i] = (ply_in[i + 2] - ply_in[i]).norm() / 2;
    }

    evaluate_intersect_rectifying(vertices, tangents,binormals,lengths);
    
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
                                             const std::vector<std::vector<Eigen::Vector3d>> &bnm)
{
    // double planar_total = 0;
    // double anglediff_total = 0;
    // double max_angle_error = 0;
    for (int i = 0; i < ply.size(); i++)
    {
        std::vector<Eigen::Vector3d> vertices;
        std::vector<Eigen::Vector3d> tangents;
        std::vector<Eigen::Vector3d> binormals;
        // double planar, anglediff, angle_max;
        if (ply[i].size() < 4)// 4 vertices, resulting in 2 rectifying planes, which has 1 crease
        {
            continue;
        }
        get_polyline_rectifying_plane_on_inner_vers(ply[i], bnm[i],
                                                 vertices, tangents, binormals);                              
        // evaluate_strip_developability(ply[i], crease[i], planar, anglediff, angle_max);
        // planar_total += planar;
        // anglediff_total += anglediff;
        // if (angle_max > max_angle_error)
        // {
        //     max_angle_error = angle_max;
        // }
    }
    // std::cout << "Evaluate current polylines: planarity: " << planar_total << ", total angle changes of the strip, " << anglediff_total << ", maximal angle error, "
    //           << max_angle_error << std::endl;
}
// compute the rectifying planes on each inner vertex.
void get_polyline_rectifying_planes(const std::vector<std::vector<Eigen::Vector3d>> &ply,
                                    const std::vector<std::vector<Eigen::Vector3d>> &bnm,
                                    std::vector<std::vector<Eigen::Vector3d>>& vertices,
                                    std::vector<std::vector<Eigen::Vector3d>>& tangents,
                                    std::vector<std::vector<Eigen::Vector3d>>& binormals)
{
    vertices.clear();
    tangents.clear();
    binormals.clear();
     for (int i = 0; i < ply.size(); i++)
    {
        std::vector<Eigen::Vector3d> ver;
        std::vector<Eigen::Vector3d> tan;
        std::vector<Eigen::Vector3d> bin;
        if (ply[i].size() < 4)// 4 vertices, resulting in 2 rectifying planes, which has 1 crease
        {
            continue;
        }
        get_polyline_rectifying_plane_on_inner_vers(ply[i], bnm[i],
                                                 ver, tan, bin);
        vertices.push_back(ver);
        tangents.push_back(tan);
        binormals.push_back(bin);                 
    }

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
    bool hasnan = false;
    for (int i = 0; i < ply_in.size(); i++)
    {
        std::vector<Eigen::Vector3d> ply_current = ply_in[i], bin_current = bi_in[i];
        std::vector<Eigen::Vector3d> ply_sampled, bin_sampled;
        int nbr = plylengths[i] / avg + 1;
        if (nbr < 3)
        {
            ply_sampled.push_back(ply_current.front());
            ply_sampled.push_back(ply_current.back());
            bin_sampled.push_back(bin_current.front());
            bin_sampled.push_back(bin_current.back());
            Eigen::Vector3d bi0 = bin_current.front();
            Eigen::Vector3d bi1 = bin_current.back();
            if (isnan(bi0[0]) || isnan(bi0[1]) || isnan(bi0[2]) || isnan(bi1[0]) || isnan(bi1[1]) || isnan(bi1[2]))
            {
                std::cout << "NAN in SHORT" << std::endl;
            }
        }
        else{
            ply_sampled = sample_one_polyline_and_binormals_based_on_length(ply_current, nbr, bin_current, bin_sampled);
        }
        
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
    // Eigen::VectorXd endpts;
    if (sample_nbr == -1) // do not sample the polylines
    {
        ply = ply_in;
        bi = bi_in;
    }
    else{
        sample_polylines_and_binormals_evenly(sample_nbr, ply_in, bi_in, ply, bi);
        std::cout<<"polyline get sampled"<<std::endl;
    }
    bool hasnan = false;
    for (int i = 0; i < bi.size(); i++)
    {
        for (int j = 0; j < bi[i].size(); j++)
        {
            if (isnan(bi[i][j][0]) || isnan(bi[i][j][1]) || isnan(bi[i][j][2]))
            {
                hasnan = true;
                break;
            }
        }
        if (hasnan)
        {
            break;
        }
    }
    if (hasnan)
    {
        std::cout << "INIT ERROR: SAMPLED BINORMALS HAS NAN" << std::endl;
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
    PlyVars = Eigen::VectorXd::Zero(vnbr * 2 * 3);
    // endpts = Eigen::VectorXd::Zero(vnbr * 6);
    VinPly = Eigen::VectorXd::Ones(vnbr * 6) * -1;
    Front.resize(vnbr);
    Back.resize(vnbr);
    int counter = 0;
    // int csize = 0;
    for (int i = 0; i < ply.size(); i++)
    {
        // if (csize < ply[i].size())
        // {
        //     csize = ply[i].size();
        // }
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
                // endpts[lpx] = 1;
                // endpts[lpy] = 1;
                // endpts[lpz] = 1;
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
                // endpts[lpx] = 1;
                // endpts[lpy] = 1;
                // endpts[lpz] = 1;
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
    ply_extracted = ply;// to record the topology of the vers
    bin_extracted = bi;
    // endpts_signs = endpts.asDiagonal();
    opt_for_polyline = true;
    opt_for_crease = false;
    // MaxNbrPinC = csize;
    if(vector_contains_NAN(PlyVars)){
        std::cout<<"INIT ERROR: variables contain NAN"<<std::endl;
    }
    
}

void PolyOpt::init(const Eigen::MatrixXd& ply_in)
{
    int vnbr = ply_in.rows();
    VerNbr = vnbr;
    // x0, x1, ..., xn, y0, ..., zn, bx0, ... bxn, ..., bzn
    PlyVars = Eigen::VectorXd::Zero(vnbr * 2 * 3);
    // endpts = Eigen::VectorXd::Zero(vnbr * 6);
    Front.resize(vnbr);
    Back.resize(vnbr);
    int counter = 0;
    std::vector<Eigen::Vector3d> binormals;
    // int csize = ply_in.size();
    binormals.resize(ply_in.rows());

    for (int j = 0; j < ply_in.rows(); j++)
    {
        int lpx = counter;
        int lpy = counter + vnbr;
        int lpz = counter + vnbr * 2;
        int lnx = counter + vnbr * 3;
        int lny = counter + vnbr * 4;
        int lnz = counter + vnbr * 5;
        Eigen::Vector3d ver = ply_in.row(j);

        PlyVars[lpx] = ver[0];
        PlyVars[lpy] = ver[1];
        PlyVars[lpz] = ver[2];
        if (j > 0 && j < ply_in.rows() - 1)
        {
            Eigen::Vector3d ver0 = ply_in.row(j - 1);
            Eigen::Vector3d ver1 = ply_in.row(j + 1);
            Eigen::Vector3d bin = (ver0 - ver).cross(ver1 - ver);
            if (bin.norm() == 0)
            {
                bin = Eigen::Vector3d(0, 0, 1);
            }
            else
            {
                bin.normalize();
            }
            PlyVars[lnx] = bin[0];
            PlyVars[lny] = bin[1];
            PlyVars[lnz] = bin[2];
            binormals[j] = bin;
        }
        if (j == 1) // the first binormal will take from the second
        {
            PlyVars[vnbr * 3] = binormals[1][0];
            PlyVars[vnbr * 4] = binormals[1][1];
            PlyVars[vnbr * 5] = binormals[1][2];
            binormals[0] = binormals[1];
        }
        if (j == 0) // first point of the seg
        {
            Front[counter] = -1;
        }
        else
        {
            Front[counter] = counter - 1;
        }
        if (j == ply_in.rows() - 1) // last point of the seg
        {
            PlyVars[lnx] = binormals[j - 1][0]; // binormals are taken from the one
            PlyVars[lny] = binormals[j - 1][1];
            PlyVars[lnz] = binormals[j - 1][2];
            binormals[j] = binormals[j - 1];
            Back[counter] = -1;
        }
        else
        {
            Back[counter] = counter + 1;
        }

        counter++;
    }
    OriVars = PlyVars;
    std::vector<std::vector<Eigen::Vector3d>> ply(1), bi(1);
    ply[0] = mat_to_vec_list(ply_in);
    bi[0] = binormals;
    RecMesh = polyline_to_strip_mesh(ply, bi, strip_scale);

    ply_extracted = ply; // to record the topology of the vers
    bin_extracted = bi;
    // endpts_signs = endpts.asDiagonal();
    opt_for_polyline = true;
    opt_for_crease = false;
    // MaxNbrPinC = csize;
    if(vector_contains_NAN(PlyVars)){
        std::cout<<"INIT ERROR: variables contain NAN"<<std::endl;
    }
    singleLineProcessing = true;
    
}
// this function is to adjust the offset for the propagation of AAG 
// it goes alternatively to avoid auxiliary variables
// verFix is the fixed vertices, vers is the foot points around verFix, creases + vers is the offset points
void adjustAagOffset(const std::vector<Eigen::Vector3d> &verFix,
                     const std::vector<Eigen::Vector3d> &vers, std::vector<Eigen::Vector3d> &creases)
{
    const double step = 0.8;
    int vnbr = vers.size();
    spMat H;
    Eigen::VectorXd B;
    std::vector<Eigen::Vector3d> versCp = vers;
    std::vector<Eigen::Vector3d> creasesCp = creases;
    for (int i = 1; i < vnbr - 1; i++)
    {
        Eigen::Vector3d cver = vers[i] + creases[i]; // offset point
        Eigen::Vector3d pf = vers[i - 1] + creases[i - 1]; // offset front
        Eigen::Vector3d pb = vers[i + 1] + creases[i + 1]; // offset back
        Eigen::Vector3d norm = (cver - verFix[i - 1]).cross(cver - verFix[i]).normalized();
        Eigen::Vector3d f = vers[i]; // the foot point.
        Eigen::Vector3d c = creases[i];// the direction
        Eigen::Vector3d A = c.cross(f-pb)+(f-pf).cross(c);
        Eigen::Vector3d B = (f - pf).cross(f - pb);
        // b = alpha * A + B
        // solve for b*n = 0
        double alpha = -B.dot(norm) / (A.dot(norm));
        double delta = (alpha - 1) * step; // rescale to be more accurate.
        creases[i] *= delta + 1;
    }
    creases[0] = creases[1];
    creases[vnbr - 1] = creases[vnbr - 2];

    // evaluation
    double errorAll = 0;
    for (int i = 1; i < vnbr - 1; i++)
    {
        Eigen::Vector3d cver = vers[i] + creases[i]; // offset point
        Eigen::Vector3d pf = vers[i - 1] + creases[i - 1]; // offset front
        Eigen::Vector3d pb = vers[i + 1] + creases[i + 1]; // offset back
        Eigen::Vector3d norm = (cver - verFix[i - 1]).cross(cver - verFix[i]).normalized();
        double error = (cver - pf).normalized().cross(cver - pb).normalized().dot(norm);
        error = error * error;
        errorAll += error;
    }
    std::cout << "error is " << errorAll << "\n";
}

// this function is designed to get the first strip for AGG evolution
// para1 is the parameter for the foot point, by default 0.5. para2 is the ratio of largest and smallest offset distances, default 1
// para3 and para4 are the start and end angles (in degree) of offset direction rotating along the tangent vector.
void aggFirstStrip(const std::vector<Eigen::Vector3d> &pts, const std::vector<Eigen::Vector3d> &bnm,
                   std::vector<Eigen::Vector3d> &pout, const bool invertDirection, const double para1, const double para2, 
                   const double para3, const double para4)
{
    int vnbr = pts.size();
    std::vector<Eigen::Vector3d> result(vnbr);
    assert(pts.size() == bnm.size());
    double angle_radian0 = para3 * LSC_PI / 180;
    double angle_radian1 = para4 * LSC_PI / 180;
    double angle_delta = (angle_radian1 - angle_radian0) / (vnbr - 1);
    // std::cout<<"the offset length along n, ";
    // the m-p is orthogonal to the plane spanned by tangent and n. m is the mid point of the edge.
    for (int i = 0; i < vnbr - 1; i++)
    {
        double angle_radian = angle_radian0 + i * angle_delta;
        double offRatio = i * (para2 - 1) / (vnbr - 1) + 1;
        Eigen::Vector3d normals = ((bnm[i] + bnm[i + 1]) / 2).normalized();
        double lengths = (pts[i] - pts[i + 1]).norm();
        Eigen::Vector3d p, direction, tangent;
        tangent = pts[i + 1] - pts[i];
        tangent.normalize();
        direction = normals.cross(tangent).normalized();
        if(invertDirection){
            direction *= -1;
        }
        // rotate against tangent
        double l = tan(angle_radian);
        // std::cout<<l<<",";
        direction = direction + l * normals;
        direction.normalize();
        p = direction * lengths * offRatio * 1.732 / 2 + (pts[i + 1] - pts[i]) * para1 + pts[i];
        result[i] = p;
    }
    Eigen::Vector3d tangent = pts[vnbr - 1] - pts[vnbr - 2];
    // Eigen::Vector3d m = (pts[vnbr - 2] + pts[vnbr - 1]) / 2;
    // auto foot = tangent + m;
    result[vnbr - 1] = result[vnbr - 2] + tangent;
    pout = result;
}
// this version uses a more robust version to get normal vectors
// ptslast is the curve next to pts.
void aggFirstStrip_following(const std::vector<Eigen::Vector3d> &pts, const std::vector<Eigen::Vector3d> &ptslast,
                   std::vector<Eigen::Vector3d> &pout, const bool invertDirection, const double para1, const double para2, 
                   const double para3, const double para4)
{
    int vnbr = pts.size();
    std::vector<Eigen::Vector3d> result(vnbr);
    double angle_radian0 = para3 * LSC_PI / 180;
    double angle_radian1 = para4 * LSC_PI / 180;
    double angle_delta = (angle_radian1 - angle_radian0) / (vnbr - 1);
    // std::cout<<"the offset length along n, ";
    // the m-p is orthogonal to the plane spanned by tangent and n. m is the mid point of the edge.
    std::vector<Eigen::Vector3d> ns(vnbr);
    for (int i = 0; i < vnbr; i++)
    {
        int fr = i + 1;
        int bk = i - 1;
        if (i == 0)
        {
            bk = i;
        }
        if (i == vnbr - 1)
        {
            fr = i;
        }
        ns[i] = (ptslast[i] - pts[i]).cross(pts[fr] - pts[bk]).normalized();
    }
    for (int i = 0; i < vnbr - 1; i++)
    {
        double angle_radian = angle_radian0 + i * angle_delta;
        double offRatio = i * (para2 - 1) / (vnbr - 1) + 1;
        Eigen::Vector3d normals= ((ns[i] + ns[i + 1]) / 2).normalized();
        double lengths = (pts[i] - pts[i + 1]).norm();
        Eigen::Vector3d p, direction, tangent;
        tangent = pts[i + 1] - pts[i];
        tangent.normalize();
        direction = normals.cross(tangent).normalized();
        if(invertDirection){
            direction *= -1;
        }
        // rotate against tangent
        double l = tan(angle_radian);
        // std::cout<<l<<",";
        direction = direction + l * normals;
        direction.normalize();
        p = direction * lengths * offRatio * 1.732 / 2 + (pts[i + 1] - pts[i]) * para1 + pts[i];
        result[i] = p;
    }
    Eigen::Vector3d tangent = pts[vnbr - 1] - pts[vnbr - 2];
    // Eigen::Vector3d m = (pts[vnbr - 2] + pts[vnbr - 1]) / 2;
    // auto foot = tangent + m;
    result[vnbr - 1] = result[vnbr - 2] + tangent;
    pout = result;
}
// take a guide curve as the first geodesic 
// extendVid is the vertex id of the guide curve that is new extended to the mesh
void aggFirstStripGuideGeodesic(const std::vector<Eigen::Vector3d> &pts, const std::vector<Eigen::Vector3d> &bnm,
                   std::vector<Eigen::Vector3d> &pout, const std::vector<Eigen::Vector3d>& guide,
                   const int segid, const double tlocal,
                    const bool invertDirection)
{
    int vnbr = pts.size();
    std::vector<Eigen::Vector3d> result(vnbr);
    assert(pts.size() == bnm.size());
    // the m-p is orthogonal to the plane spanned by tangent and n. m is the mid point of the edge.

    result[0] = tlocal * (guide[segid + 1] - guide[segid]) + guide[segid];// the first point comes from the guide curve
    double l, t;
    t = (result[0] - pts[0]).dot(pts[1] - pts[0]) / (pts[1] - pts[0]).dot(pts[1] - pts[0]);
    l = (pts[0] - result[0] + (pts[1] - pts[0]) * t).norm();

    for (int i = 1; i < vnbr - 1; i++)
    {
        
        Eigen::Vector3d normals = ((bnm[i] + bnm[i + 1]) / 2).normalized();
        // double lengths = (pts[i] - pts[i + 1]).norm();
        Eigen::Vector3d p, direction, tangent;
        tangent = pts[i + 1] - pts[i];
        tangent.normalize();
        direction = normals.cross(tangent).normalized();
        if(invertDirection){
            direction *= -1;
        }
        // double offRatio = i * (para2 - 1) / (vnbr - 1) + 1;
        Eigen::Vector3d foot = pts[i] + (pts[i + 1] - pts[i]) * t;
        p = direction * l  + foot;
        result[i] = p;
    }
    Eigen::Vector3d tangent = pts[vnbr - 1] - pts[vnbr - 2];
    // Eigen::Vector3d m = (pts[vnbr - 2] + pts[vnbr - 1]) / 2;
    // auto foot = tangent + m;
    result[vnbr - 1] = result[vnbr - 2] + tangent;
    pout = result;
}




// adjust the last boundary of AGG after propagating.
// the bnms are the binormals of the second last boundary
void adjustAggOffset(const std::vector<Eigen::Vector3d> &pts_all,
                     const int vinrow, std::vector<Eigen::Vector3d> &pout)
{
    pout.clear();
    pout.resize(vinrow);
    int vnbr = pts_all.size();
    int ncol = pts_all.size() / vinrow;
    assert(ncol > 2 && "we don't apply this method on the first strip"); 
    int counter = 0;
    std::vector<Eigen::Vector3d> pout_all;
    pout_all = pts_all;
    double ratio = 0.3;
    double error = 0;
    std::vector<Eigen::Vector3d> ns(vinrow);
    for (int i = 0; i < vinrow; i++)
    {
        int fr = (ncol - 1) * vinrow + i + 1;
        int bk = (ncol - 1) * vinrow + i - 1;
        if (i == 0)
        {
            bk = (ncol - 1) * vinrow + i;
        }
        if (i == vinrow - 1)
        {
            fr = (ncol - 1) * vinrow + i;
        }
        ns[i] = (pts_all[(ncol - 2) * vinrow + i] - pts_all[(ncol - 1) * vinrow + i]).cross(pts_all[fr] - pts_all[bk]).normalized();
    }
    for (int i = 0; i < vinrow - 2; i++)
    {
        int bid = (ncol - 1) * vinrow + i; // boundary point id
        Eigen::Vector3d p = pts_all[bid];
        Eigen::Vector3d p1 = pts_all[bid - vinrow];
        Eigen::Vector3d p2 = pts_all[bid - vinrow + 1];
        Eigen::Vector3d p3 = pts_all[bid - vinrow * 2]; // there is a useless point between p3 and p4
        Eigen::Vector3d p4 = pts_all[bid - vinrow * 2 + 2];
        Eigen::Vector3d v1 = p - p1;
        Eigen::Vector3d v2 = p1 - p3;
        Eigen::Vector3d v3 = p - p2;
        Eigen::Vector3d v4 = p2 - p4;
        Eigen::Vector3d n1 = ns[i];
        Eigen::Vector3d n2 = ns[i + 1];

        double r11 = v1.cross(n1).dot(v2);
        double r12 = v3.cross(n1).dot(v2);
        double r21 = v1.cross(n2).dot(v4);
        double r22 = v3.cross(n2).dot(v4);
        double d1 = (p1 - p).cross(n1).dot(v2);
        double d2 = (p2 - p).cross(n2).dot(v4);
        error += pow(((p - p1).normalized()).cross(n1).dot(v2.normalized()), 2);
        error += pow(((p - p2).normalized()).cross(n2).dot(v4.normalized()), 2);
        Eigen::Matrix2d mat, inv;
        mat << r11, r12, r21, r22;
        Eigen::Vector2d d;
        d << d1, d2;
        bool invertable;
        mat.computeInverseWithCheck(inv, invertable);
        Eigen::Vector2d ab;
        
        if (invertable)
        {
            ab = inv * d * 0.3;
        }
        else
        {
            ab << 0, 0;
            if (abs(r11) > 1e-4)
            {
                ab[0] = (d1 - ab[1] * r12) / r11 * ratio;
            }
            if (abs(r12) > 1e-4)
            {
                ab[1] = (d1 - ab[0] * r11) / r12 * ratio;
            }
            if (abs(r21) > 1e-4)
            {
                ab[0] = (d2 - ab[1] * r22) / r21 * ratio;
            }
            if (abs(r22) > 1e-4)
            {
                ab[1] = (d2 - ab[0] * r21) / r22 * ratio;
            }
        }
        pout_all[bid] += ab[0] * v1 + ab[1] * v3;
        pout[i] = pout_all[bid];
    }
    // deal with the second last point: only requires p1, p3, n1
    int bid = vnbr - 2; // boundary point id
    Eigen::Vector3d p = pts_all[bid];
    Eigen::Vector3d p1 = pts_all[bid - vinrow];
    Eigen::Vector3d p2 = pts_all[bid - vinrow + 1];
    Eigen::Vector3d p3 = pts_all[bid - vinrow * 2]; // there is a useless point between p3 and p4
    Eigen::Vector3d v1 = p - p1;
    Eigen::Vector3d v2 = p1 - p3;
    Eigen::Vector3d v3 = p - p2;
    Eigen::Vector3d n1 = ns[vinrow - 2];
    double r11 = v1.cross(n1).dot(v2);
    double r12 = v3.cross(n1).dot(v2);
    double d1 = (p1 - p).cross(n1).dot(v2);
    Eigen::Vector2d ab(0, 0);
    if (abs(r11) > 1e-4)
    {
        ab[0] = (d1 - ab[1] * r12) / r11 * ratio;
    }
    if (abs(r12) > 1e-4)
    {
        ab[1] = (d1 - ab[0] * r11) / r12 * ratio;
    }
    // pout_all[bid] += ab[0] * v1 + ab[1] * v3; 


    // the second last point smoothness
    Eigen::Vector3d delta = (pout_all[bid + 1] + pout_all[bid - 1]) / 2 - pout_all[bid];
    pout_all[bid] += delta * ratio; // smoothness 1: row

    double alpha = 1/(pout_all[bid] - pout_all[bid - vinrow]).norm();
    double beta = 1/(pout_all[bid - vinrow] - pout_all[bid - 2 * vinrow]).norm();
    delta = (alpha + beta) / alpha * pout_all[bid - vinrow] - beta / alpha * pout_all[bid - 2 * vinrow] - pout_all[bid];
    // delta = 2 * pout_all[bid - vinrow] - pout_all[bid - 2 * vinrow] - pout_all[bid]; // one order smoothness
    pout_all[bid] += delta * ratio; // smoothness 2: col
    pout[vinrow - 2] = pout_all[bid];
    error += pow(((p - p1).normalized()).cross(n1).dot(v2.normalized()), 2);

    // the last point
    bid = vnbr - 1; // boundary point id
    // the first smooth
    delta = 2 * pout_all[bid - 1] - pout_all[bid - 2] - pout_all[bid]; // one order smoothness
    pout_all[bid] += delta ;                                    // smoothness
    // the second smooth
    delta = 2 * pout_all[bid - vinrow] - pout_all[bid - 2 * vinrow] - pout_all[bid]; // one order smoothness
    pout_all[bid] += delta ;                                                  // smoothness

    pout[vinrow - 1] = pout_all[bid];
    std::cout << "squared error before opt, " << error << "\n";
    
}

// adjust the last boundary of AGG after propagating.
// the bnms are the binormals of the second last boundary
// this version takes a input curve as the first geodesic
void adjustAggOffsetGuideGeodesic(const std::vector<Eigen::Vector3d> &pts_all, const std::vector<Eigen::Vector3d> &bnms,
                     const int vinrow, const std::vector<Eigen::Vector3d> &guide, int& segid, double& tlocal,
                     std::vector<Eigen::Vector3d> &pout)
{
    pout.clear();
    pout.resize(vinrow);
    int vnbr = pts_all.size();
    int ncol = pts_all.size() / vinrow;
    assert(ncol > 2 && "we don't apply this method on the first strip"); 
    assert(vinrow == bnms.size());
    int counter = 0;
    std::vector<Eigen::Vector3d> pout_all;
    pout_all = pts_all;
    double ratio = 0.3;
    double error = 0;
    for (int i = 0; i < vinrow - 2; i++)
    {
        int bid = (ncol - 1) * vinrow + i; // boundary point id
        Eigen::Vector3d p = pts_all[bid];
        Eigen::Vector3d p1 = pts_all[bid - vinrow];
        Eigen::Vector3d p2 = pts_all[bid - vinrow + 1];
        Eigen::Vector3d p3 = pts_all[bid - vinrow * 2]; // there is a useless point between p3 and p4
        Eigen::Vector3d p4 = pts_all[bid - vinrow * 2 + 2];
        Eigen::Vector3d v1 = p - p1;
        Eigen::Vector3d v2 = p1 - p3;
        Eigen::Vector3d v3 = p - p2;
        Eigen::Vector3d v4 = p2 - p4;
        Eigen::Vector3d n1 = bnms[i];
        Eigen::Vector3d n2 = bnms[i + 1];

        double r11 = v1.cross(n1).dot(v2);
        double r12 = v3.cross(n1).dot(v2);
        double r21 = v1.cross(n2).dot(v4);
        double r22 = v3.cross(n2).dot(v4);
        double d1 = (p1 - p).cross(n1).dot(v2);
        double d2 = (p2 - p).cross(n2).dot(v4);
        error += pow(((p - p1).normalized()).cross(n1).dot(v2.normalized()), 2);
        error += pow(((p - p2).normalized()).cross(n2).dot(v4.normalized()), 2);

        if (i == 0) // the first point is on the curve: 
        {
            Eigen::Vector3d pt0 = guide[segid], pt1 = guide[segid + 1]; // the current position should be pt0 + t * (pt1 - pt0)
            double left = (pt1 - pt0).cross(n2).dot(v4);
            double right = (-pt0).cross(n2).dot(v4);
            double t_target = right / left; // the target position if move a full step. 
            t_target = (t_target - tlocal) * ratio + tlocal; // this is the target position
            double target_length = (t_target - tlocal) * (pt0 - pt1).norm(); // the distance need to march along the curve.
            Eigen::Vector3d pstart = tlocal * (pt1 - pt0) + pt0;
            int segid_out;
            Eigen::Vector3d pt;
            if (target_length > 0)
            {
                double tl;
                find_next_pt_on_polyline(segid, guide, target_length, pstart, segid_out, pt, tl);   
            }
            else{
                find_prev_pt_on_polyline(segid, guide, target_length, pstart, segid_out, pt);
            }
            segid = segid_out;
            tlocal = (pt - guide[segid]).norm() / (guide[segid + 1] - guide[segid]).norm();

            pout_all[bid] = pt;
            pout[i] = pout_all[bid];
            continue;
        }

        Eigen::Matrix2d mat, inv;
        mat << r11, r12, r21, r22;
        Eigen::Vector2d d;
        d << d1, d2;
        bool invertable;
        mat.computeInverseWithCheck(inv, invertable);
        Eigen::Vector2d ab;
        
        if (invertable)
        {
            ab = inv * d * 0.3;
        }
        else
        {
            ab << 0, 0;
            if (abs(r11) > 1e-4)
            {
                ab[0] = (d1 - ab[1] * r12) / r11 * ratio;
            }
            if (abs(r12) > 1e-4)
            {
                ab[1] = (d1 - ab[0] * r11) / r12 * ratio;
            }
            if (abs(r21) > 1e-4)
            {
                ab[0] = (d2 - ab[1] * r22) / r21 * ratio;
            }
            if (abs(r22) > 1e-4)
            {
                ab[1] = (d2 - ab[0] * r21) / r22 * ratio;
            }
        }
        pout_all[bid] += ab[0] * v1 + ab[1] * v3;
        pout[i] = pout_all[bid];
    }
    // deal with the second last point: only requires p1, p3, n1
    int bid = vnbr - 2; // boundary point id
    Eigen::Vector3d p = pts_all[bid];
    Eigen::Vector3d p1 = pts_all[bid - vinrow];
    Eigen::Vector3d p2 = pts_all[bid - vinrow + 1];
    Eigen::Vector3d p3 = pts_all[bid - vinrow * 2]; // there is a useless point between p3 and p4
    Eigen::Vector3d v1 = p - p1;
    Eigen::Vector3d v2 = p1 - p3;
    Eigen::Vector3d v3 = p - p2;
    Eigen::Vector3d n1 = bnms[vinrow - 2];
    double r11 = v1.cross(n1).dot(v2);
    double r12 = v3.cross(n1).dot(v2);
    double d1 = (p1 - p).cross(n1).dot(v2);
    Eigen::Vector2d ab(0, 0);
    if (abs(r11) > 1e-4)
    {
        ab[0] = (d1 - ab[1] * r12) / r11 * ratio;
    }
    if (abs(r12) > 1e-4)
    {
        ab[1] = (d1 - ab[0] * r11) / r12 * ratio;
    }
    // pout_all[bid] += ab[0] * v1 + ab[1] * v3; 


    // the second last point smoothness
    Eigen::Vector3d delta = (pout_all[bid + 1] + pout_all[bid - 1]) / 2 - pout_all[bid];
    pout_all[bid] += delta * ratio; // smoothness 1: row

    double alpha = 1/(pout_all[bid] - pout_all[bid - vinrow]).norm();
    double beta = 1/(pout_all[bid - vinrow] - pout_all[bid - 2 * vinrow]).norm();
    delta = (alpha + beta) / alpha * pout_all[bid - vinrow] - beta / alpha * pout_all[bid - 2 * vinrow] - pout_all[bid];
    // delta = 2 * pout_all[bid - vinrow] - pout_all[bid - 2 * vinrow] - pout_all[bid]; // one order smoothness
    pout_all[bid] += delta * ratio; // smoothness 2: col
    pout[vinrow - 2] = pout_all[bid];
    error += pow(((p - p1).normalized()).cross(n1).dot(v2.normalized()), 2);

    // the last point
    bid = vnbr - 1; // boundary point id
    // the first smooth
    delta = 2 * pout_all[bid - 1] - pout_all[bid - 2] - pout_all[bid]; // one order smoothness
    pout_all[bid] += delta ;                                    // smoothness
    // the second smooth
    delta = 2 * pout_all[bid - vinrow] - pout_all[bid - 2 * vinrow] - pout_all[bid]; // one order smoothness
    pout_all[bid] += delta ;                                                  // smoothness

    pout[vinrow - 1] = pout_all[bid];
    std::cout << "squared error before opt, " << error << "\n";
    
}



// flat points: where franet frame is not properly defined. including: low curvature points, inflection points
void statics_for_curvatures(const Eigen::VectorXd &cvec, std::vector<bool>& flat, double &cthreadshold)
{
    Eigen::VectorXd curvatures = cvec;
    int nend = 0;
    int nlv0 = 0, nlv1 = 0, nlv2 = 0, nlv3 = 0;
    double cmax = curvatures.maxCoeff();
    flat.resize(curvatures.size(), false);
    for (int i = 0; i < cvec.size(); i++)
    {
        if (cvec[i] < 0)
        {
            curvatures[i] = 0;
            nend++;
        }
        else{
            if(curvatures[i]<0.01*cmax){
                nlv0++;
                flat[i] = true;
            }
            if(curvatures[i]<0.005*cmax){
                nlv1++;
            }
            if(curvatures[i]<0.001*cmax){
                nlv2++;
            }
            
        }
    }
    cthreadshold = 0.01 * cmax;
    std::cout << "Curvature Statics: Max curvature, " << cmax << ", below 1%, " << nlv0 << ", 5e-3, " << nlv1 << ", 1e-4, " << nlv2 << ", out of, "
              << curvatures.size() << "\n";
}

void PolyOpt::init_crease_opt(const std::vector<std::vector<Eigen::Vector3d>> &vertices,
                              const std::vector<std::vector<Eigen::Vector3d>> &tangents,
                              const std::vector<std::vector<Eigen::Vector3d>> &binormals)
{

    int vnbr = 0;
    for (auto line : vertices)
    {
        for (auto ver : line)
        {
            vnbr++;
        }
    }
    VerNbr = vnbr;
    
    Front.resize(vnbr);
    Back.resize(vnbr);
    vertices_cp.resize(vnbr, 3);
    tangents_cp.resize(vnbr, 3);
    binormal_cp.resize(vnbr, 3);
    VinPly = Eigen::VectorXd::Ones(vnbr) * -1;
    int counter = 0;
    for (int i = 0; i < vertices.size(); i++)
    {
        for (int j = 0; j < vertices[i].size(); j++)
        {
            int location = counter;
            if (j == 0) // first point of the seg
            {
                Front[counter] = -1;
            }
            else
            {
                Front[counter] = counter - 1;
            }
            if (j == vertices[i].size() - 1) // last point of the seg
            {
                Back[counter] = -1;
            }
            else
            {
                Back[counter] = counter + 1;
            }
            VinPly[location] = i;
            vertices_cp.row(counter) = vertices[i][j];
            tangents_cp.row(counter) = tangents[i][j];
            binormal_cp.row(counter) = binormals[i][j];

            counter++;
        }
    }
    // the variables are the vertices, creases, normals of the faces, binormals, principle normals.
    PlyVars = Eigen::VectorXd::Zero(vnbr * 15);
    // vertices
    PlyVars.segment(0, vnbr) = vertices_cp.col(0);
    PlyVars.segment(vnbr, vnbr) = vertices_cp.col(1);
    PlyVars.segment(vnbr * 2, vnbr) = vertices_cp.col(2);
    // creases
    PlyVars.segment(vnbr * 3, vnbr) = binormal_cp.col(0);
    PlyVars.segment(vnbr * 4, vnbr) = binormal_cp.col(1);
    PlyVars.segment(vnbr * 5, vnbr) = binormal_cp.col(2);
    // binormals
    PlyVars.segment(vnbr * 9, vnbr) = binormal_cp.col(0);
    PlyVars.segment(vnbr * 10, vnbr) = binormal_cp.col(1);
    PlyVars.segment(vnbr * 11, vnbr) = binormal_cp.col(2);
    
    int ncompu = 0;

    // to record the curvature on each point.
    AngCollector = Eigen::VectorXd::Ones(vnbr) * -1;
    Inflecs = Eigen::VectorXd::Ones(vnbr) * -1;
    // std::vector<bool> skips(vnbr, false);
    for (int i = 0; i < vnbr; i++)
    {
        int bk = Back[i];
        if (bk == -1)
        {
            continue;
        }
        int lnx = i + vnbr * 6;
        int lny = i + vnbr * 7;
        int lnz = i + vnbr * 8;
        int lpx = i + vnbr * 12;
        int lpy = i + vnbr * 13;
        int lpz = i + vnbr * 14;
        Eigen::Vector3d Ver = vertices_cp.row(i);
        Eigen::Vector3d Vbk = vertices_cp.row(bk);
        Eigen::Vector3d Bver = binormal_cp.row(i);
        if (abs((Ver - Vbk).normalized().dot(Bver)) > 0.1)
        {
            std::cout<<"WARNING: The Binormal is inaccurate!!!"<<std::endl;
        }
        Eigen::Vector3d norm = (Ver-Vbk).cross(Bver).normalized();
        PlyVars[lnx] = norm[0];
        PlyVars[lny] = norm[1];
        PlyVars[lnz] = norm[2];

        Eigen::Vector3d tangent = tangents_cp.row(i);
        Eigen::Vector3d pnormal = tangent.cross(Bver).normalized();
        PlyVars[lpx] = tangent[0];
        PlyVars[lpy] = tangent[1];
        PlyVars[lpz] = tangent[2];

        ncompu++;
        // // compute the curvature, and predicate inflation points
        // int fr = Front[i];
        // if(fr == -1){
        //     continue;
        // }
        // // curvature
        // Eigen::Vector3d Vfr = vertices_cp.row(fr);
        // Eigen::Vector3d tleft = (Ver - Vfr).normalized();
        // Eigen::Vector3d tright = (Vbk - Ver).normalized();
        // double lleft = (Ver - Vfr).norm();
        // double lright = (Vbk - Ver).norm();
        // double length = (lleft + lright) / 2;
        // Eigen::Vector3d cvector = (tright - tleft) / length;
        // double curvature = cvector.norm();
        // AngCollector[i] = curvature;
        // // inflection point
        // if (Back[bk] == -1)
        // {
        //     continue;
        // }
        // Eigen::Vector3d Vbb = vertices_cp.row(Back[bk]);
        // Eigen::Vector3d Bfr = (Ver - Vfr).cross(Vbk - Ver);
        // Eigen::Vector3d Bbk = (Vbk - Ver).cross(Vbb - Vbk);
        // if (Bfr.dot(Bbk) < 0) // check inflection points using the binormals
        // {
        //     // Inflecs[i] = 1;
        //     // Inflecs[bk] = 1;
        //     skips[i] = true;
        //     skips[bk] = true;
        //     continue;
            
        // }
        

        // // Eigen::Vector3d Bfr = binormal_cp.row(fr);
        // // Eigen::Vector3d Bbk = binormal_cp.row(bk);
        // // Eigen::Vector3d bcross0 = Bver.cross(Bfr);
        // // Eigen::Vector3d bcross1 = Bbk.cross(Bver);
        // // if (bcross0.dot(bcross1) < 0) // check inflection points using the binormals
        // // {
        // //     Inflecs[i] = 1;
        // //     ninflc++;
        // // }
    }
    ply_extracted = vertices;// to record the topology of the vers
    bin_extracted = binormals;
    OriVars = PlyVars;
    RecMesh = polyline_to_strip_mesh(vertices, binormals, strip_scale);
    opt_for_crease = true;
    opt_for_polyline = false;
    double cthreadshold;
    FlatPts.resize(vnbr, false);
    // statics_for_curvatures(AngCollector, FlatPts, cthreadshold);
    int ninflc = 0;
    // for (int i = 0; i < vnbr; i++)
    // {
    //     if (Inflecs[i] == 1)
    //     {
    //         ninflc++;
    //     }
    //     if(skips[i]){
    //         continue;
    //     }
    //     int bk = Back[i];
    //     if (bk == -1)
    //     {
    //         continue;
    //     }
    //     int fr = Front[i];
    //     if(fr == -1){
    //         continue;
    //     }
    //     if(AngCollector[i]<cthreadshold){
    //         continue;
    //     }
    //     Eigen::Vector3d Ver = vertices_cp.row(i);
    //     Eigen::Vector3d Vbk = vertices_cp.row(bk);
    //     Eigen::Vector3d Bver = binormal_cp.row(i);
    //     Eigen::Vector3d Vfr = vertices_cp.row(fr);
    //     Eigen::Vector3d tleft = (Ver - Vfr).normalized();
    //     Eigen::Vector3d tright = (Vbk - Ver).normalized();
    //     double lleft = (Ver - Vfr).norm();
    //     double lright = (Vbk - Ver).norm();
    //     double length = (lleft + lright) / 2;
    //     // use the tortion to 
    //     Eigen::Vector3d Bbk = binormal_cp.row(bk); 
    //     Eigen::Vector3d Bfr = binormal_cp.row(fr);

    //     Eigen::Vector3d Bbm = (Bbk + Bver).normalized();
    //     Eigen::Vector3d Bfm = (Bfr + Bver).normalized();

    //     Eigen::Vector3d Bdiff = (Bbm - Bfm) / length;
    //     double tortion = Bdiff.norm();
    //     // normalize tortion and curvature
    //     double tocu = sqrt(tortion * tortion + AngCollector[i] * AngCollector[i]);
    //     PlyVars[i] = tortion / tocu;
    //     PlyVars[i + vnbr] = AngCollector[i] / tocu;
    // }
    // std::cout<<"Nbr of Inflection Points: "<<ninflc<<std::endl;
   
    
    std::cout<<"Initialization done. Ver nbr: "<<vertices_cp.rows()<<", nbr of curves, "<<vertices.size()<<
    ", nbr of faces, "<<ncompu<<std::endl;

}
void PolyOpt::draw_inflections(Eigen::MatrixXd& pts){
    std::vector<Eigen::Vector3d> tmp;
    for(int i=0;i<vertices_cp.rows();i++){
        if(Inflecs[i]>0){
            tmp.push_back(vertices_cp.row(i));
        }
    }
    pts = vec_list_to_matrix(tmp);

}
// use this when avoiding flipping of the strips.
void PolyOpt::force_smoothing_binormals(){
    int counter = 0;
    for (int i = 0; i < ply_extracted.size(); i++)
    {
        
        // counter is the first point. need to use the middle point as reference
        int mid = counter + ply_extracted[i].size() / 2;
        int lnx = mid + VerNbr * 3;
        int lny = mid + VerNbr * 4;
        int lnz = mid + VerNbr * 5;

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
// use this when avoiding flipping of the strips.
void PolyOpt::force_invert_binormals(){
    int counter = 0;
    for (int i = 0; i < ply_extracted.size(); i++)
    {
        for (int j = 0; j < ply_extracted[i].size(); j++)
        {
            int lnx = counter + VerNbr * 3;
            int lny = counter + VerNbr * 4;
            int lnz = counter + VerNbr * 5;
            PlyVars[lnx] *= -1;
            PlyVars[lny] *= -1;
            PlyVars[lnz] *= -1;
            counter++;
        }
    }
}
void PolyOpt::related_error(){
    int vnbr = VerNbr;
    Eigen::MatrixXd C;
    Eigen::VectorXi I;
    Eigen::VectorXd D;
    
    Eigen::MatrixXd vers(vnbr, 3);
    vers.col(0) = PlyVars.segment(0, vnbr);
    vers.col(1) = PlyVars.segment(vnbr, vnbr);
    vers.col(2) = PlyVars.segment(vnbr * 2, vnbr);
    if(Vref.cols() != vers.cols()){
        std::cout<<"In PolyOpt::assemble_gravity: Vref.cols() "<<Vref.cols()<<std::endl;
        return;
    }
    assert(Vref.cols() == vers.cols());
    igl::point_mesh_squared_distance(vers, Vref, Fref, D, I, C);
    double maxdis = sqrt(D.maxCoeff());
    Eigen::Vector3d vmin(Vref.col(0).minCoeff(), Vref.col(1).minCoeff(), Vref.col(2).minCoeff());
    Eigen::Vector3d vmax(Vref.col(0).maxCoeff(), Vref.col(1).maxCoeff(), Vref.col(2).maxCoeff());
    double bbd = (vmin - vmax).norm();
    std::cout << "the ratio of approximate distance / bounding box diagonal: " << maxdis / bbd << std::endl;
}
void PolyOpt::assemble_gravity(spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy)
{
    int vnbr = VerNbr;
    std::vector<Eigen::Vector3d> vprojs(vnbr);
    std::vector<Eigen::Vector3d> Nlocal(vnbr);
    Eigen::MatrixXd C;
    Eigen::VectorXi I;
    Eigen::VectorXd D;
    int nq = vnbr;
    Eigen::MatrixXd vers(vnbr, 3);
    vers.col(0) = PlyVars.segment(0, vnbr);
    vers.col(1) = PlyVars.segment(vnbr, vnbr);
    vers.col(2) = PlyVars.segment(vnbr * 2, vnbr);
    if(Vref.cols() != vers.cols()){
        std::cout<<"In PolyOpt::assemble_gravity: Vref.cols() "<<Vref.cols()<<std::endl;
    }
    assert(Vref.cols() == vers.cols());
    igl::point_mesh_squared_distance(vers, Vref, Fref, D, I, C);

    for (int i = 0; i < vnbr; i++)
    {
        int f = I(i);

        Eigen::Vector3d p = C.row(i);
        Eigen::Vector3d v0 = Vref.row(Fref(f, 0));
        Eigen::Vector3d v1 = Vref.row(Fref(f, 1));
        Eigen::Vector3d v2 = Vref.row(Fref(f, 2));
        std::array<double, 3> coor = barycenter_coordinate(v0, v1, v2, p);
        Eigen::Vector3d norm = coor[0] * Nref.row(Fref(f, 0)) + coor[1] * Nref.row(Fref(f, 1)) + coor[2] * Nref.row(Fref(f, 2));
        assert(Nref.row(Fref(f, 0)).dot(Nref.row(Fref(f, 1))) > 0 && Nref.row(Fref(f, 0)).dot(Nref.row(Fref(f, 2))) > 0);
        norm.normalize();
        vprojs[i] = p;
        Nlocal[i] = norm;

        // debug
        Eigen::Vector3d rec_pro = coor[0] * v0 + coor[1] * v1 + coor[2] * v2;
        if ((rec_pro - p).norm() > 0.5)
        {
            std::cout<<"INACCURATE SOLVING: coor: "<<coor[0]<<", "<<coor[1]<<", "<<coor[2]<<std::endl;
            std::cout<<"v0, "<<v0.transpose()<<std::endl;
            std::cout<<"v1, "<<v1.transpose()<<std::endl;
            std::cout<<"v2, "<<v2.transpose()<<std::endl;
            std::cout<<"p "<<p.transpose()<<std::endl;
            std::cout<<"p_compute "<<rec_pro.transpose()<<std::endl;
        }
    }
    // Ppro0 = vec_list_to_matrix(vprojs);
    // Npro0 = vec_list_to_matrix(Nlocal);
    std::vector<Trip> tripletes;
    tripletes.reserve(vnbr * 9);
    energy = Eigen::VectorXd::Zero(vnbr * 3);
    for (int i = 0; i < vnbr; i++)
    {
        int vid = i;
        int lx = vid;
        int ly = vid + vnbr;
        int lz = vid + vnbr * 2;
        int lnx = vid + vnbr * 3;
        int lny = vid + vnbr * 4;
        int lnz = vid + vnbr * 5;

        Eigen::Vector3d normal = Nlocal[vid]; // the normal vector
        Eigen::Vector3d closest = vprojs[vid]; // the closest point
        Eigen::Vector3d ver(PlyVars[lx], PlyVars[ly], PlyVars[lz]);

        // ver.dot(n^*) - ver^*.dot(n^*) = 0
        tripletes.push_back(Trip(i, lx, normal[0]));
        tripletes.push_back(Trip(i, ly, normal[1]));
        tripletes.push_back(Trip(i, lz, normal[2]));

        energy[i] = ver.dot(normal) - closest.dot(normal);

        // vertices not far away from the original ones
        // (ver - ver^*)^2 = 0
        double scale = 1;
        Eigen::Vector3d verori;
        verori = closest;
        // Eigen::Vector3d normal_local = normal;
        // verori[0] = OriVars[lx];
        // verori[1] = OriVars[ly];
        // verori[2] = OriVars[lz];
        Eigen::Vector3d vdiff = ver - verori;
        // if (vdiff.dot(normal_local) < 0)
        // {
        //     normal_local *= -1;
        // }
        // double cos_angle = vdiff.normalized().dot(normal_local);
        // double angle_radian = acos(cos_angle);
        // if (angle_radian > 30)
        // {
        //     continue;
        // }
        // std::cout<<"Through here"<<std::endl;
        tripletes.push_back(Trip(i + vnbr, lx, 2 * vdiff[0] * scale));
        tripletes.push_back(Trip(i + vnbr, ly, 2 * vdiff[1] * scale));
        tripletes.push_back(Trip(i + vnbr, lz, 2 * vdiff[2] * scale));

        energy[i + vnbr] = vdiff.dot(vdiff) * scale;

        // close to the original position
        verori[0] = OriVars[lx];
        verori[1] = OriVars[ly];
        verori[2] = OriVars[lz];
        vdiff = ver - verori;
        scale = 1;

        tripletes.push_back(Trip(i + vnbr * 2, lx, 2 * vdiff[0] * scale));
        tripletes.push_back(Trip(i + vnbr * 2, ly, 2 * vdiff[1] * scale));
        tripletes.push_back(Trip(i + vnbr * 2, lz, 2 * vdiff[2] * scale));

        energy[i + vnbr * 2] = vdiff.dot(vdiff) * scale;
    }
    int nvars = PlyVars.size();
    int ncondi = energy.size();
    spMat J;
    J.resize(ncondi, nvars);
    J.setFromTriplets(tripletes.begin(), tripletes.end());
    H = J.transpose() * J;
    B = -J.transpose() * energy;
}

void PolyOpt::assemble_gravity_single_ply(spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy)
{
    int vnbr = VerNbr;
    
    std::vector<Trip> tripletes;
    tripletes.reserve(vnbr * 3);
    energy = Eigen::VectorXd::Zero(vnbr);
    for (int i = 0; i < vnbr; i++)
    {
        int vid = i;
        int lx = vid;
        int ly = vid + vnbr;
        int lz = vid + vnbr * 2;
        int lnx = vid + vnbr * 3;
        int lny = vid + vnbr * 4;
        int lnz = vid + vnbr * 5;

        Eigen::Vector3d ver(PlyVars[lx], PlyVars[ly], PlyVars[lz]);
        // vertices not far away from the original ones
        // (ver - ver^*)^2 = 0

        Eigen::Vector3d verori;

        // close to the original position
        verori[0] = OriVars[lx];
        verori[1] = OriVars[ly];
        verori[2] = OriVars[lz];
        auto vdiff = ver - verori;
        double scale = 1;

        tripletes.push_back(Trip(i, lx, 2 * vdiff[0] * scale));
        tripletes.push_back(Trip(i, ly, 2 * vdiff[1] * scale));
        tripletes.push_back(Trip(i, lz, 2 * vdiff[2] * scale));

        energy[i] = vdiff.dot(vdiff) * scale;
    }
    int nvars = PlyVars.size();
    int ncondi = energy.size();
    spMat J;
    J.resize(ncondi, nvars);
    J.setFromTriplets(tripletes.begin(), tripletes.end());
    H = J.transpose() * J;
    B = -J.transpose() * energy;
}
void PolyOpt::assemble_polyline_smooth(const bool crease, spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy)
{
    if (Front.size() == 0)
    {
        std::cout << "please use assemble_polyline_smooth after init the PolyOpt" << std::endl;
        return;
    }

    std::vector<Trip> tripletes;
    int vnbr = Front.size();
    energy = Eigen::VectorXd::Zero(vnbr * 6);
    tripletes.reserve(vnbr * 25);
    // std::cout<<"check 1"<<std::endl;
    for (int i = 0; i < vnbr; i++)
    {

        int vid = i;
        if (pick_single_line && VinPly.size() > 0)
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

        if(crease){ // if opt for crease, the binormal is in a different place
            lnx = vid + vnbr * 9;
            lny = vid + vnbr * 10;
            lnz = vid + vnbr * 11;

            lfnx = fr + vnbr * 9;
            lfny = fr + vnbr * 10;
            lfnz = fr + vnbr * 11;

            lbnx = bk + vnbr * 9;
            lbny = bk + vnbr * 10;
            lbnz = bk + vnbr * 11;
        }
        Eigen::Vector3d ver(PlyVars[lx], PlyVars[ly], PlyVars[lz]);
        Eigen::Vector3d vfr(PlyVars[lfx], PlyVars[lfy], PlyVars[lfz]);
        Eigen::Vector3d vbk(PlyVars[lbx], PlyVars[lby], PlyVars[lbz]);
        Eigen::Vector3d nbi(PlyVars[lnx], PlyVars[lny], PlyVars[lnz]);
        Eigen::Vector3d nfr(PlyVars[lfnx], PlyVars[lfny], PlyVars[lfnz]);
        Eigen::Vector3d nbk(PlyVars[lbnx], PlyVars[lbny], PlyVars[lbnz]);

        if (compute_line_smth)
        {
            int bbk = Back[bk];
            if (bbk >= 0)// the last point exist
            {
                int lbbx = bbk;
                int lbby = bbk + vnbr;
                int lbbz = bbk + vnbr * 2;
                // vfr - 3 * ver + 3 * vbk - bbk = 0
                Eigen::Vector3d vbb(PlyVars[lbbx], PlyVars[lbby], PlyVars[lbbz]);
                Eigen::Vector3d err = vfr - 3 * ver + 3 * vbk - vbb;
                // x
                tripletes.push_back(Trip(i, lfx, 1));
                tripletes.push_back(Trip(i, lx, -3));
                tripletes.push_back(Trip(i, lbx, 3));
                tripletes.push_back(Trip(i, lbbx, -1));
                energy[i] = err[0];

                // y
                tripletes.push_back(Trip(i + vnbr, lfy, 1));
                tripletes.push_back(Trip(i + vnbr, ly, -3));
                tripletes.push_back(Trip(i + vnbr, lby, 3));
                tripletes.push_back(Trip(i + vnbr, lbby, -1));
                energy[i + vnbr] = err[1];

                // z
                tripletes.push_back(Trip(i + vnbr * 2, lfz, 1));
                tripletes.push_back(Trip(i + vnbr * 2, lz, -3));
                tripletes.push_back(Trip(i + vnbr * 2, lbz, 3));
                tripletes.push_back(Trip(i + vnbr * 2, lbbz, -1));
                energy[i + vnbr * 2] = err[2];

                // double dis0 = (ver - vfr).norm();
                // double dis1 = (ver - vbk).norm();
                // Eigen::Vector3d error = (ver - vfr) / dis0 + (ver - vbk) / dis1;

                // // (v - fr)/||v-fr|| + (v-bk)/||v-bk|| = 0

                // // x
                // tripletes.push_back(Trip(i, lx, 1 / dis0 + 1 / dis1));
                // tripletes.push_back(Trip(i, lfx, -1 / dis0));
                // tripletes.push_back(Trip(i, lbx, -1 / dis1));
                // energy[i] = error[0];

                // // y
                // tripletes.push_back(Trip(i + vnbr, ly, 1 / dis0 + 1 / dis1));
                // tripletes.push_back(Trip(i + vnbr, lfy, -1 / dis0));
                // tripletes.push_back(Trip(i + vnbr, lby, -1 / dis1));

                // energy[i + vnbr] = error[1];

                // // z
                // tripletes.push_back(Trip(i + vnbr * 2, lz, 1 / dis0 + 1 / dis1));
                // tripletes.push_back(Trip(i + vnbr * 2, lfz, -1 / dis0));
                // tripletes.push_back(Trip(i + vnbr * 2, lbz, -1 / dis1));

                // energy[i + vnbr * 2] = error[2];
            }
        }
        // the 2nd order smoothness of the binormals
        if (compute_line_smth){
            double sign1, sign2;
            sign1 = nbi.dot(nfr) > 0 ? 1 : -1;
            sign2 = nbi.dot(nbk) > 0 ? 1 : -1;
            // sign1 * nfr + sign2 * nbk - 2 * nbi = 0
            Eigen::Vector3d err = sign1 * nfr + sign2 * nbk - 2 * nbi;
            tripletes.push_back(Trip(i + vnbr * 3, lfnx, binormal_ratio * sign1));
            tripletes.push_back(Trip(i + vnbr * 3, lbnx, binormal_ratio * sign2));
            tripletes.push_back(Trip(i + vnbr * 3, lnx, -2 * binormal_ratio));
            energy[i + vnbr * 3] = binormal_ratio * err[0];

            tripletes.push_back(Trip(i + vnbr * 4, lfny, binormal_ratio * sign1));
            tripletes.push_back(Trip(i + vnbr * 4, lbny, binormal_ratio * sign2));
            tripletes.push_back(Trip(i + vnbr * 4, lny, -2 * binormal_ratio));
            energy[i + vnbr * 4] = binormal_ratio * err[1];

            tripletes.push_back(Trip(i + vnbr * 5, lfnz, binormal_ratio * sign1));
            tripletes.push_back(Trip(i + vnbr * 5, lbnz, binormal_ratio * sign2));
            tripletes.push_back(Trip(i + vnbr * 5, lnz, -2 * binormal_ratio));
            energy[i + vnbr * 5] = binormal_ratio * err[2];
            if(isnan(binormal_ratio)){
                std::cout<<"binormal_ratio is nan"<<std::endl;
            }
            if(isnan(sign1)){
                std::cout<<"sign1 is nan"<<std::endl;
            }
            if(isnan(sign2)){
                std::cout<<"sign2 is nan"<<std::endl;
            }
            if(isnan(err[0]) || isnan(err[1]) || isnan(err[2])){
                std::cout<<"err[0] is nan"<<std::endl;
                // std::cout<<"error, "<<err<<std::endl;
                // std::cout << "sign1, " << sign1 << std::endl;
                // std::cout << "sign2, " << sign2 << std::endl;
                // std::cout << "nfr, " << nfr << std::endl;
                // std::cout << "nbk, " << nbk << std::endl;
                // std::cout << "nbi, " << nbi << std::endl;
                // std::cout << "original nfr, " << OriVars[lfnx]<<","<< OriVars[lfny]<<", "<< OriVars[lfnz] << std::endl;
                // std::cout << "original nbk, " << OriVars[lbnx]<<","<< OriVars[lbny]<<", "<< OriVars[lbnz] << std::endl;
                // std::cout<<"line id "<<VinPly[lfnx]<<", "<<VinPly[lbnx]<<std::endl;
                
                
            }
            // if (vector_contains_NAN(PlyVars))
            // {
            //     std::cout << "Polyvars has nan" << std::endl;
            // }
        }
        
    }
    // assert(!vector_contains_NAN(energy));
    int nanposition;
    if(vector_contains_NAN(energy, nanposition)){
        std::cout<<"NAN in Ply smooth, location: "<<nanposition<<", with vnbr "<<vnbr<<std::endl;
    }
    // std::cout<<"check 2"<<std::endl;
    spMat J;
    J.resize(energy.size(), PlyVars.size());
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
        if (pick_single_line && VinPly.size() > 0)
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

// k is from 1.
bool within_k_pts_to_bnd(const Eigen::VectorXi &Front, const Eigen::VectorXi &Back, const int vid, const int k)
{
    int tid = vid;
    for (int i = 0; i < k; i++)
    {
        tid = Front[tid];
        if (tid == -1)
        {
            return true;
        }
    }
    tid = vid;
    for (int i = 0; i < k; i++)
    {
        tid = Back[tid];
        if (tid == -1)
        {
            return true;
        }
    }
    return false;
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
    double percent_filter = 5;//filter 5% of the points close to the end points
    // int nfilter = std::max(MaxNbrPinC * percent_filter / 100 / 2, 1.0); // filter n points on each side.
    // WARNING I removed this nfilter but maybe lower the robustness
    for (int i = 0; i < vnbr; i++)
    {
        int vid = i;
        // if (within_k_pts_to_bnd(Front, Back, vid, nfilter))
        // {
        //     continue;
        // }

        if (pick_single_line && VinPly.size() > 0)
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

// please use me after initializing the polylines
// Vr, Fr and nr are for the reference surfaces
void PolyOpt::get_normal_vector_from_reference(const Eigen::MatrixXd &Vr, const Eigen::MatrixXi &Fr, const Eigen::MatrixXd &nr)
{
    if(Vref.rows() != Vr.rows()){
        Vref = Vr;
        Fref = Fr;
        Nref = nr;
    }
    
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
void PolyOpt::extract_rectifying_plane_mesh_from_crease(){
    int nv = Front.size();
    Eigen::MatrixXd v0list(nv, 3), v1list(nv, 3);
    Eigen::MatrixXd binlist(nv, 3);
    
    for (int i = 0; i < nv; i++)
    {
        // the crease is defined as alpha * tangent + beta * binormal. This is to 
        // make sure that the strip width is 1.
        int lvx = i;
        int lvy = i + nv;
        int lvz = i + nv * 2;
        int lcx = i + nv * 3;
        int lcy = i + nv * 4;
        int lcz = i + nv * 5;
        int lbx = i + nv * 9;
        int lby = i + nv * 10;
        int lbz = i + nv * 11;

        Eigen::Vector3d bnm(PlyVars[lbx],PlyVars[lby],PlyVars[lbz]);
        Eigen::Vector3d crs(PlyVars[lcx],PlyVars[lcy],PlyVars[lcz]);
        Eigen::Vector3d ver(PlyVars[lvx],PlyVars[lvy],PlyVars[lvz]);
        double scale = 1 / (crs.dot(bnm));
        v0list.row(i) = ver;
        binlist.row(i) = crs * scale;
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
    assert(ply_extracted.size() > 0);
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
    assert(ply_extracted.size() > 0);
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

void PolyOpt::extract_polylines_and_binormals_from_creases(){
    int nv = Front.size();
    Eigen::MatrixXd v0list(nv, 3);
    Eigen::MatrixXd binlist(nv, 3);
    for (int i = 0; i < nv; i++)
    {
        int lvx = i;
        int lvy = i + nv;
        int lvz = i + nv * 2;
        int lcx = i + nv * 3;
        int lcy = i + nv * 4;
        int lcz = i + nv * 5;
        int lbx = i + nv * 9;
        int lby = i + nv * 10;
        int lbz = i + nv * 11;
        Eigen::Vector3d ver(PlyVars[lvx],PlyVars[lvy],PlyVars[lvz]);
        Eigen::Vector3d bnm(PlyVars[lbx],PlyVars[lby],PlyVars[lbz]);
        Eigen::Vector3d crs(PlyVars[lcx],PlyVars[lcy],PlyVars[lcz]);
        double scale = 1 / (crs.dot(bnm));
        binlist.row(i) = crs * scale;
        v0list.row(i) = ver;
        
    }
    int counter = 0;
    assert(ply_extracted.size() > 0);
    for (int i = 0; i < ply_extracted.size(); i++)
    {
        for (int j = 0; j < ply_extracted[i].size(); j++)
        {
            ply_extracted[i][j] = v0list.row(counter);
            bin_extracted[i][j] = binlist.row(counter);
            counter ++;
        }
    }
}
// if rotated, we save the rotated mesh; otherwise we save the mesh extracted from optimizer.
void PolyOpt::save_polyline_and_binormals_as_files(const bool rotated)
{
    std::cout << "Saving the binormal files, please provide the prefix" << std::endl;
    std::string fname = igl::file_dialog_save();
    if (fname.length() == 0)
    {
        std::cout << "\nLSC: save mesh failed, please type down the correct name" << std::endl;
        return;
    }
    std::vector<std::vector<Eigen::Vector3d>> lines, binormals;
    if (rotated)
    {
        lines = ply_rotated;
        binormals = bin_rotated;
    }
    else
    {
        assert(ply_extracted.size() > 0);
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
    Eigen::Vector3d bbase = OrientEndPts ? bi[0].back() : bi[0][0];
    if (!OrientEndPts)
    {
        for (int i = 0; i < bi.size(); i++)
        {
            Eigen::Vector3d base = OrientEndPts ? bi[i].back() : bi[i][0];
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
    }
    else{
        for (int i = 0; i < bi.size(); i++)
        {
            Eigen::Vector3d base = OrientEndPts ? bi[i].back() : bi[i][0];
            if (base.dot(bbase) < 0)
            {
                base *= -1;
            }
            bbase = base;

            for (int j = bi[i].size() - 1; j > -1; j--)
            {
                Eigen::Vector3d direction = orient_vector(base, bi[i][j]);
                base = direction;
                bilist[i][j] = direction;
            }
        }
    }

    bout = bilist;
    OrientEndPts = !OrientEndPts;
}
void PolyOpt::opt()
{
    if (opt_for_crease)
    {
        std::cout << "The environment is polluted, Please re-start the program" << std::endl;
        return;
    }
    spMat H;
    Eigen::VectorXd B;

    spMat Hmass;
    Eigen::VectorXd Bmass, emass;

    if (singleLineProcessing)
    {
        assemble_gravity_single_ply(Hmass, Bmass, emass);
    }
    else
    {
        assemble_gravity(Hmass, Bmass, emass);
    }

    H = weight_mass * Hmass;
    B = weight_mass * Bmass;

    spMat Hsmt;
    Eigen::VectorXd Bsmt, esmt;

    assemble_polyline_smooth(false, Hsmt, Bsmt, esmt);

    H += weight_smooth * Hsmt;
    B += weight_smooth * Bsmt;

    spMat Hbin;
    Eigen::VectorXd Bbin, ebin;

    assemble_binormal_condition(Hbin, Bbin, ebin);

    H += weight_binormal * Hbin;
    B += weight_binormal * Bbin;

    spMat Hangle;
    Eigen::VectorXd Bangle, Eangle;
    if(weight_angle > 0 && singleLineProcessing == false){
        assemble_angle_condition(Hangle, Bangle,Eangle);
        H += weight_angle * Hangle;
        B += weight_angle * Bangle;
    }
    H += 1e-6 * spMat(Eigen::VectorXd::Ones(H.rows()).asDiagonal());
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
    if(weight_angle > 0 && singleLineProcessing == false){
        double energy_angle = Eangle.norm();
        std::cout<<", Energy angle, "<<energy_angle<<", Eangle max, "<< Eangle.lpNorm<Eigen::Infinity>();
    }
    std::cout<<"\n";
    extract_rectifying_plane_mesh();
    extract_polylines_and_binormals();
    opt_for_polyline = true;
    //
}

// try to make the variable to 0
void PolyOpt::assemble_gravity_crease(spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy)
{
    std::vector<Trip> tripletes;
    int vnbr = Front.size();
    tripletes.reserve(vnbr * 15);
    energy = Eigen::VectorXd::Zero(vnbr * 5);
    
    for(int i = 0;i<vnbr;i++){
        int vid = i;
        int bk = Back[vid];
        int fr = Front[vid];
        bool compute_front = true;
        bool compute_back = true;

        if(fr == -1){
            compute_front = false;
            fr = vid;
        }

        if (bk == -1)
        {
            compute_back = false;
            bk = vid;
        }
        // locations: vertices, creases, normals of the faces, binormals, principle normals.
        int lvx = vid;
        int lvy = vid + vnbr;
        int lvz = vid + vnbr * 2;
        int lfrx = fr;
        int lfry = fr + vnbr;
        int lfrz = fr + vnbr * 2;
        int lbkx = bk;
        int lbky = bk + vnbr;
        int lbkz = bk + vnbr * 2;

        int lbx = vid + vnbr * 9;
        int lby = vid + vnbr * 10;
        int lbz = vid + vnbr * 11;

        double xo = OriVars[lvx];
        double yo = OriVars[lvy];
        double zo = OriVars[lvz];
        Eigen::Vector3d bio(OriVars[lbx],OriVars[lby],OriVars[lbz]);

        Eigen::Vector3d ver(PlyVars[lvx],PlyVars[lvy],PlyVars[lvz]);
        Eigen::Vector3d vfr(PlyVars[lfrx],PlyVars[lfry],PlyVars[lfrz]);
        Eigen::Vector3d vbk(PlyVars[lbkx],PlyVars[lbky],PlyVars[lbkz]);
        
        // ver[0] - xo = 0
        tripletes.push_back(Trip(i, lvx, 1));
        energy[i] = ver[0] - xo;
        
        // ver[1] - yo = 0
        tripletes.push_back(Trip(i + vnbr, lvy, 1));
        energy[i + vnbr] = ver[1] - yo;

        // ver[2] - zo = 0
        tripletes.push_back(Trip(i + vnbr * 2, lvz, 1));
        energy[i + vnbr * 2] = ver[2] - zo;

        // if (compute_front && compute_back)
        if (1)
        {
            int lcx = vid + vnbr * 3;
            int lcy = vid + vnbr * 4;
            int lcz = vid + vnbr * 5;
            Eigen::Vector3d crease(PlyVars[lcx],PlyVars[lcy],PlyVars[lcz]);

            // (bio x tangent) * crease = 0
            Eigen::Vector3d tangent = (vbk - vfr).normalized();
            Eigen::Vector3d bxt = bio.cross(tangent);
            double err = bxt.dot(crease);
            tripletes.push_back(Trip(i + vnbr * 3, lcx, bxt[0]));
            tripletes.push_back(Trip(i + vnbr * 3, lcy, bxt[1]));
            tripletes.push_back(Trip(i + vnbr * 3, lcz, bxt[2]));
            energy[i + vnbr * 3] = err;

            // // bio * (ver - vfr) / scale = 0
            // double scale = (ver - vfr).norm();
            // tripletes.push_back(Trip(i + vnbr * 3, lvx, binormal_ratio * bio[0] / scale));
            // tripletes.push_back(Trip(i + vnbr * 3, lvy, binormal_ratio * bio[1] / scale));
            // tripletes.push_back(Trip(i + vnbr * 3, lvz, binormal_ratio * bio[2] / scale));

            // tripletes.push_back(Trip(i + vnbr * 3, lfrx, - binormal_ratio * bio[0] / scale));
            // tripletes.push_back(Trip(i + vnbr * 3, lfry, - binormal_ratio * bio[1] / scale));
            // tripletes.push_back(Trip(i + vnbr * 3, lfrz, - binormal_ratio * bio[2] / scale));

            // energy[i + vnbr * 3] = binormal_ratio * bio.dot(ver - vfr) / scale;

            // // bio * (ver - vbk) / scale = 0;
            // scale = (ver - vbk).norm();
            // tripletes.push_back(Trip(i + vnbr * 4, lvx, binormal_ratio * bio[0] / scale));
            // tripletes.push_back(Trip(i + vnbr * 4, lvy, binormal_ratio * bio[1] / scale));
            // tripletes.push_back(Trip(i + vnbr * 4, lvz, binormal_ratio * bio[2] / scale));

            // tripletes.push_back(Trip(i + vnbr * 4, lbkx, - binormal_ratio * bio[0] / scale));
            // tripletes.push_back(Trip(i + vnbr * 4, lbky, - binormal_ratio * bio[1] / scale));
            // tripletes.push_back(Trip(i + vnbr * 4, lbkz, - binormal_ratio * bio[2] / scale));

            // energy[i + vnbr * 4] = binormal_ratio * bio.dot(ver - vbk) / scale;
        }
    }

    spMat J;
    J.resize(energy.size(), PlyVars.size());
    J.setFromTriplets(tripletes.begin(), tripletes.end());
    H = J.transpose() * J;
    B = -J.transpose() * energy;
}
void PolyOpt::assemble_crease_planarity(spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy)
{
    if (Front.size() == 0)
    {
        std::cout << "Please correctly initialize the creases" << std::endl;
        return; 
    }
    
    std::vector<Trip> tripletes;
    int vnbr = Front.size();
    // Eigen::VectorXi Markers = Eigen::VectorXi::Zero(vnbr);
    energy = Eigen::VectorXd::Zero(vnbr * 9);
    tripletes.reserve(vnbr * 60);
    FlatDeterm = Eigen::VectorXd::Zero(vnbr);;
    for (int i = 0; i < vnbr; i++)
    {
        int vid = i;
        int bk = Back[vid];
        int fr = Front[vid];
        bool compute_front = true;
        bool compute_back = true;

        if(fr == -1){
            compute_front = false;
            fr = vid;
        }

        if (bk == -1)
        {
            compute_back = false;
            bk = vid;
        }
        // locations: vertices, creases, normals of the faces, binormals, principle normals.
        int lvx = vid;
        int lvy = vid + vnbr;
        int lvz = vid + vnbr * 2;
        int lfrx = fr;
        int lfry = fr + vnbr;
        int lfrz = fr + vnbr * 2;
        int lbkx = bk;
        int lbky = bk + vnbr;
        int lbkz = bk + vnbr * 2;

        
        int lcx = vid + vnbr * 3;
        int lcy = vid + vnbr * 4;
        int lcz = vid + vnbr * 5;

        int lcbx = bk + vnbr * 3;
        int lcby = bk + vnbr * 4;
        int lcbz = bk + vnbr * 5;
        
        int lnx = vid + vnbr * 6;
        int lny = vid + vnbr * 7;
        int lnz = vid + vnbr * 8;

        int lbx = vid + vnbr * 9;
        int lby = vid + vnbr * 10;
        int lbz = vid + vnbr * 11;

        int lpx = vid + vnbr * 12;
        int lpy = vid + vnbr * 13;
        int lpz = vid + vnbr * 14;

        Eigen::Vector3d bnm(PlyVars[lbx],PlyVars[lby],PlyVars[lbz]);
        Eigen::Vector3d pnm(PlyVars[lpx],PlyVars[lpy],PlyVars[lpz]);
        Eigen::Vector3d crs(PlyVars[lcx],PlyVars[lcy],PlyVars[lcz]);
        Eigen::Vector3d bcrs(PlyVars[lcbx],PlyVars[lcby],PlyVars[lcbz]);
        
        Eigen::Vector3d ver(PlyVars[lvx],PlyVars[lvy],PlyVars[lvz]);
        Eigen::Vector3d norm(PlyVars[lnx],PlyVars[lny],PlyVars[lnz]);
        Eigen::Vector3d vfr(PlyVars[lfrx],PlyVars[lfry],PlyVars[lfrz]);
        Eigen::Vector3d vbk(PlyVars[lbkx],PlyVars[lbky],PlyVars[lbkz]);

        if(compute_back){
            Eigen::Vector3d cl = crs.normalized();
            Eigen::Vector3d cr = bcrs.normalized();
            Eigen::Vector3d tangent = (ver - vbk).normalized();

            FlatDeterm[i] = cl.cross(cr).dot(tangent);
        }
        if(compute_front){
            // binormal * (ver - vfr) = 0
            double err = bnm.dot((ver - vfr).normalized());
            double scale = (ver - vfr).norm();
            tripletes.push_back(Trip(i, lbx, (ver - vfr)[0] / scale));
            tripletes.push_back(Trip(i, lby, (ver - vfr)[1] / scale));
            tripletes.push_back(Trip(i, lbz, (ver - vfr)[2] / scale));

            tripletes.push_back(Trip(i, lvx, bnm[0] / scale));
            tripletes.push_back(Trip(i, lvy, bnm[1] / scale));
            tripletes.push_back(Trip(i, lvz, bnm[2] / scale));

            tripletes.push_back(Trip(i, lfrx, -bnm[0] / scale));
            tripletes.push_back(Trip(i, lfry, -bnm[1] / scale));
            tripletes.push_back(Trip(i, lfrz, -bnm[2] / scale));

            energy[i] = err;

        }
        if(compute_back){
            // binormal * (ver - vbk) = 0
            double err = bnm.dot((ver - vbk).normalized());
            double scale = (ver - vbk).norm();
            tripletes.push_back(Trip(i + vnbr, lbx, (ver - vbk)[0] / scale));
            tripletes.push_back(Trip(i + vnbr, lby, (ver - vbk)[1] / scale));
            tripletes.push_back(Trip(i + vnbr, lbz, (ver - vbk)[2] / scale));

            tripletes.push_back(Trip(i + vnbr, lvx, bnm[0] / scale));
            tripletes.push_back(Trip(i + vnbr, lvy, bnm[1] / scale));
            tripletes.push_back(Trip(i + vnbr, lvz, bnm[2] / scale));

            tripletes.push_back(Trip(i + vnbr, lbkx, -bnm[0] / scale));
            tripletes.push_back(Trip(i + vnbr, lbky, -bnm[1] / scale));
            tripletes.push_back(Trip(i + vnbr, lbkz, -bnm[2] / scale));

            energy[i + vnbr] = err;
        }
        // binormal * binormal = 1
        tripletes.push_back(Trip(i + vnbr * 2, lbx, 2 * bnm[0]));
        tripletes.push_back(Trip(i + vnbr * 2, lby, 2 * bnm[1]));
        tripletes.push_back(Trip(i + vnbr * 2, lbz, 2 * bnm[2]));

        energy[i + vnbr * 2] = bnm.dot(bnm) - 1;

        if (compute_back)
        {
            // norm * (ver - vbk) = 0
            double scale = (ver - vbk).norm();
            double err = norm.dot(ver - vbk) / scale;
            tripletes.push_back(Trip(i + vnbr * 3, lnx, (ver - vbk)[0] / scale));
            tripletes.push_back(Trip(i + vnbr * 3, lny, (ver - vbk)[1] / scale));
            tripletes.push_back(Trip(i + vnbr * 3, lnz, (ver - vbk)[2] / scale));

            tripletes.push_back(Trip(i + vnbr * 3, lvx, norm[0] / scale));
            tripletes.push_back(Trip(i + vnbr * 3, lvy, norm[1] / scale));
            tripletes.push_back(Trip(i + vnbr * 3, lvz, norm[2] / scale));

            tripletes.push_back(Trip(i + vnbr * 3, lbkx, -norm[0] / scale));
            tripletes.push_back(Trip(i + vnbr * 3, lbky, -norm[1] / scale));
            tripletes.push_back(Trip(i + vnbr * 3, lbkz, -norm[2] / scale));
            energy[i + vnbr * 3] = err;
            // norm * crs = 0
            tripletes.push_back(Trip(i + vnbr * 4, lnx, crs[0]));
            tripletes.push_back(Trip(i + vnbr * 4, lny, crs[1]));
            tripletes.push_back(Trip(i + vnbr * 4, lnz, crs[2]));

            tripletes.push_back(Trip(i + vnbr * 4, lcx, norm[0]));
            tripletes.push_back(Trip(i + vnbr * 4, lcy, norm[1]));
            tripletes.push_back(Trip(i + vnbr * 4, lcz, norm[2]));

            energy[i + vnbr * 4] = norm.dot(crs);

            // norm * bcrs = 0

            tripletes.push_back(Trip(i + vnbr * 5, lnx, bcrs[0]));
            tripletes.push_back(Trip(i + vnbr * 5, lny, bcrs[1]));
            tripletes.push_back(Trip(i + vnbr * 5, lnz, bcrs[2]));

            tripletes.push_back(Trip(i + vnbr * 5, lcbx, norm[0]));
            tripletes.push_back(Trip(i + vnbr * 5, lcby, norm[1]));
            tripletes.push_back(Trip(i + vnbr * 5, lcbz, norm[2]));

            energy[i + vnbr * 5] = norm.dot(bcrs);

            // norm * norm = 1
            tripletes.push_back(Trip(i + vnbr * 6, lnx, 2 * norm[0]));
            tripletes.push_back(Trip(i + vnbr * 6, lny, 2 * norm[1]));
            tripletes.push_back(Trip(i + vnbr * 6, lnz, 2 * norm[2]));

            energy[i + vnbr * 6] = norm.dot(norm) - 1;
        }
        if(compute_front && compute_back){
            double scale0 = (vfr - ver).norm();
            double scale1 = (vbk - ver).norm();
            // crease * (vfr - ver) / scale0 + crease * (vbk - ver) / scale1 = 0
            tripletes.push_back(Trip(i + vnbr * 7, lcx, ((vfr - ver) / scale0 + (vbk - ver) / scale1)[0]));
            tripletes.push_back(Trip(i + vnbr * 7, lcy, ((vfr - ver) / scale0 + (vbk - ver) / scale1)[1]));
            tripletes.push_back(Trip(i + vnbr * 7, lcz, ((vfr - ver) / scale0 + (vbk - ver) / scale1)[2]));

            tripletes.push_back(Trip(i + vnbr * 7, lfrx, crs[0] / scale0));
            tripletes.push_back(Trip(i + vnbr * 7, lfry, crs[1] / scale0));
            tripletes.push_back(Trip(i + vnbr * 7, lfrz, crs[2] / scale0));

            tripletes.push_back(Trip(i + vnbr * 7, lbkx, crs[0] / scale1));
            tripletes.push_back(Trip(i + vnbr * 7, lbky, crs[1] / scale1));
            tripletes.push_back(Trip(i + vnbr * 7, lbkz, crs[2] / scale1));

            tripletes.push_back(Trip(i + vnbr * 7, lvx, -crs[0] / scale0 - crs[0] / scale1));
            tripletes.push_back(Trip(i + vnbr * 7, lvy, -crs[1] / scale0 - crs[1] / scale1));
            tripletes.push_back(Trip(i + vnbr * 7, lvz, -crs[2] / scale0 - crs[2] / scale1));

            energy[i + vnbr * 7] = crs.dot(vfr - ver) / scale0 + crs.dot(vbk - ver) / scale1;
        }
        //crease * crease = 1
        tripletes.push_back(Trip(i + vnbr * 8, lcx, 2 * crs[0]));
        tripletes.push_back(Trip(i + vnbr * 8, lcy, 2 * crs[1]));
        tripletes.push_back(Trip(i + vnbr * 8, lcz, 2 * crs[2]));

        energy[i + vnbr * 8] = crs.dot(crs) - 1;



        // pnm * bnm = 0?
        

        // pnm * (vfr - vbk) = 0?

        // pnm * pnm = 1?

        // crease * pnm = 0?

    }
    spMat J;
    J.resize(energy.size(), PlyVars.size());
    J.setFromTriplets(tripletes.begin(), tripletes.end());
    H = J.transpose() * J;
    B = -J.transpose() * energy;
}

void PolyOpt::opt_planarity(){
    if(opt_for_polyline){
        std::cout<<"The environment is polluted, Please re-start the program"<<std::endl;
        return;
    }
    int vnbr = Front.size();
    spMat H;
    Eigen::VectorXd B;
    spMat Hmass;
    Eigen::VectorXd Bmass, emass;
     assemble_gravity_crease(Hmass, Bmass, emass);
    H = weight_mass * Hmass;
    B = weight_mass * Bmass;

    spMat Hsmt;
    Eigen::VectorXd Bsmt, esmt;
    assemble_polyline_smooth(true, Hsmt, Bsmt, esmt);
    H += weight_smooth * Hsmt;
    B += weight_smooth * Bsmt;

    spMat Hangle;
    Eigen::VectorXd Bangle, Eangle;

    if(weight_angle > 0){
        assemble_crease_planarity(Hangle, Bangle,Eangle);
        H += weight_angle * Hangle;
        B += weight_angle * Bangle;
    }
    H += 1e-6 * spMat(Eigen::VectorXd::Ones(H.rows()).asDiagonal());
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

    double stplength = dx.norm();
    double energy_gravity = emass.norm();

    std::cout << "Crease, gravity, " << energy_gravity<<", gravity max, "<<emass.lpNorm<Eigen::Infinity>();
    if (binormal_ratio != 0)
    {
        // energy of approximating binormals
        Eigen::VectorXd eab = emass.segment(vnbr * 3, 2 * vnbr) / binormal_ratio;
        std::cout << ", bnm appro max, " << eab.lpNorm<Eigen::Infinity>();
    }
    std::cout << ", smt, " << esmt.norm() << ", step, " << stplength;
    if(weight_angle > 0){
        double energy_angle = Eangle.norm();
        std::cout<<", Energy Planarity, "<<energy_angle<<", Planar max, "<< Eangle.lpNorm<Eigen::Infinity>()
        <<", planar determinate, "<<FlatDeterm.norm()<<", max determinate, "<<FlatDeterm.lpNorm<Eigen::Infinity>();
    }
    std::cout<<"\n";
    
    extract_rectifying_plane_mesh_from_crease();
    extract_polylines_and_binormals_from_creases();
    opt_for_crease = true;
}
