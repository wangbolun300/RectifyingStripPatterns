#include <lsc/basic.h>
#include <lsc/tools.h>
// takes edge point and pseudo-normal as inputs, output the
IGL_INLINE void QuadricCalculator::get_pseudo_vertex_on_edge(
    const int center_id,
    const Eigen::Vector3d &epoint, const Eigen::Vector3d &pnormal, Eigen::Vector3d &pvertex)
{
    // CHECK che esista la mesh
    const size_t vertices_count = vertices.rows();

    if (vertices_count == 0)
        return;

    std::vector<int> vv;
    std::vector<int> vvtmp;
    vv.clear();
    vvtmp.clear();
    Eigen::Vector3d me = vertices.row(center_id);
    getKRing(center_id, kRing, vv);
    if (vv.size() < 6)
    {
        std::cerr << "Error in calculating pseudo-vertices, because too few neighbouring vertices " << std::endl;
        return;
    }
    if (projectionPlaneCheck)
    {
        vvtmp.reserve(vv.size());
        applyProjOnPlane(pnormal, vv, vvtmp);
        if (vvtmp.size() >= 6 && vvtmp.size() < vv.size())
            vv = vvtmp;
    }
    if (vv.size() < 6)
    {
        std::cerr << "Error in calculating pseudo-vertices, because too few neighbouring vertices " << std::endl;
        return;
    }

    std::vector<Eigen::Vector3d> ref(3);
    computeReferenceFrame(center_id, pnormal, ref);
    Quadric q;
    fitQuadric(me, ref, vv, &q);
    const double a = q.a();
    const double b = q.b();
    const double c = q.c();
    const double d = q.d();
    const double e = q.e();
    // the surface is a*x*x+b*x*y+c*y*y+d*x*e*y
    // (ep_x, ep_y, ep_z) is the local position of the pseudo-vertex
    double ep_x;
    double ep_y;
    double ep_z;
    Eigen::Vector3d relative = epoint - me;
    ep_x = relative.dot(ref[0]);
    ep_y = relative.dot(ref[1]);
    ep_z = a * ep_x * ep_x + b * ep_x * ep_y + c * ep_y * ep_y + d * ep_x + e * ep_y;
    Eigen::Vector3d global_vec = ref[0] * ep_x + ref[1] * ep_y + ref[2] * ep_z;
    pvertex = me + global_vec;
}

bool quadratic_solver(const std::vector<double> &func, std::array<double, 2> &roots)
{
    assert(func.size() <= 3);
    if (func.size() < 2)
    {
        // std::cout<<"failed in quadratic solver: no roots"<<std::endl;
        return false;
    }
    if (func.size() == 2)
    {
        // the function is f1*x+f0=0
        double f1 = func[1];
        double f0 = func[0];
        double result = -1. * f0 / f1;
        roots[0] = result;
        roots[1] = result;
        return true;
    }
    double a = func[2], b = func[1], c = func[0];
    double lambda = b * b - 4 * a * c;
    if (lambda < 0)
    {
        return false;
    }
    roots[0] = (-1. * b + sqrt(lambda)) / (2 * a);
    roots[1] = (-1. * b - sqrt(lambda)) / (2 * a);
    return true;
}

// double

// this function is designed to handle geodesic cases, which does not need to solve a quadratic function
// if p1_is_ver, then search the edge vs-ve; if p1_is_ver is false, then only search the triangle which corresponds to edge_middle
bool find_geodesic_intersection_p1_is_ver(const Eigen::Vector3d &p0, const Eigen::Vector3d &p1,
                                          const Eigen::Vector3d &vs, const Eigen::Vector3d &ve, const Eigen::Vector3d &pnorm,
                                          Eigen::Vector3d &p_end)
{
    Eigen::Vector3d ves = ve - vs;
    Eigen::Vector3d vs1 = vs - p1;
    Eigen::Vector3d v01 = p0 - p1;

    // normal director of the plane is norm_degree1*t+norm_degree0
    Eigen::Vector3d norm_degree1 = ves.cross(v01);
    Eigen::Vector3d norm_degree0 = vs1.cross(v01);

    std::vector<double> poly_norm_x(2), poly_norm_y(2), poly_norm_z(2);
    poly_norm_x[0] = norm_degree0[0];
    poly_norm_y[0] = norm_degree0[1];
    poly_norm_z[0] = norm_degree0[2];
    poly_norm_x[1] = norm_degree1[0];
    poly_norm_y[1] = norm_degree1[1];
    poly_norm_z[1] = norm_degree1[2];
    std::vector<double> left;
    left = polynomial_times(poly_norm_x, pnorm[0]);
    left = polynomial_add(left, polynomial_times(poly_norm_y, pnorm[1]));
    left = polynomial_add(left, polynomial_times(poly_norm_z, pnorm[2]));
    assert(left.size() < 3);
    if (left.size() < 2)
    {
        return false;
    }
    double t = -left[0] / left[1];
    if (t < -MERGE_VERTEX_RATIO || t > 1 + MERGE_VERTEX_RATIO)
    {
        return false;
    }
    if (t >= -MERGE_VERTEX_RATIO && t < 0)
    {
        t = 0;
    }
    if (t > 1 && t <= 1 + MERGE_VERTEX_RATIO)
    {
        t = 1;
    }
    if (t >= 0 && t <= 1)
    {
        p_end = vs + t * ves;
        return true;
    }
    return false;
}

// p0 p1 is the first segment, vs ve is the middle edge, vf, vt is the searching edge
bool solve_next_geodesic_point(const Eigen::Vector3d &p0, const Eigen::Vector3d &p1,
                               const Eigen::Vector3d &vf, const Eigen::Vector3d &vt,
                               const Eigen::Vector3d &pnorm, Eigen::Vector3d &p_end)
{
    Eigen::Vector3d vtf = vt - vf;
    Eigen::Vector3d p01 = p0 - p1;
    Eigen::Vector3d vf1 = vf - p1;
    double degree1 = vtf.cross(p01).dot(pnorm);
    double degree0 = vf1.cross(p01).dot(pnorm);
    if (degree1 == 0)
    {
        return false;
    }
    double t = -degree0 / degree1;
    if (t < -MERGE_VERTEX_RATIO || t > 1 + MERGE_VERTEX_RATIO)
    {
        return false;
    }
    if (t >= -MERGE_VERTEX_RATIO && t < 0)
    {
        t = 0;
    }
    if (t > 1 && t <= 1 + MERGE_VERTEX_RATIO)
    {
        t = 1;
    }
    if (t >= 0 && t <= 1)
    {
        p_end = vf + t * vtf;
        return true;
    }
    return false;
}

// caution: this function does not check if the middle edge is the boundary. this should be done before calling this function
bool lsTools::find_geodesic_intersection_p1_is_NOT_ver(const Eigen::Vector3d &p0, const Eigen::Vector3d &p1,
                                                       const CGMesh::HalfedgeHandle &edge_middle, const Eigen::Vector3d &pnorm,
                                                       CGMesh::HalfedgeHandle &edge_out, Eigen::Vector3d &p_end)
{
    CGMesh::HalfedgeHandle ophe = lsmesh.opposite_halfedge_handle(edge_middle);
    // vs ve is the middle edge
    // Eigen::Vector3d vs = V.row(lsmesh.from_vertex_handle(edge_middle).idx());
    // Eigen::Vector3d ve = V.row(lsmesh.to_vertex_handle(edge_middle).idx());
    // check next halfedge
    CGMesh::HalfedgeHandle checking_he = lsmesh.next_halfedge_handle(ophe);
    Eigen::Vector3d vf;
    Eigen::Vector3d vt;
    assert(lsmesh.face_handle(edge_middle).idx() != lsmesh.face_handle(ophe).idx());
    assert(lsmesh.face_handle(ophe).idx() == lsmesh.face_handle(checking_he).idx());
    vf = V.row(lsmesh.from_vertex_handle(checking_he).idx());
    vt = V.row(lsmesh.to_vertex_handle(checking_he).idx());
    bool found = solve_next_geodesic_point(p0, p1, vf, vt, pnorm, p_end);
    if (found)
    {
        edge_out = checking_he;
        return true;
    }
    // check the previous halfedge
    checking_he = lsmesh.prev_halfedge_handle(ophe);
    assert(lsmesh.face_handle(edge_middle).idx() != lsmesh.face_handle(ophe).idx());
    assert(lsmesh.face_handle(ophe).idx() == lsmesh.face_handle(checking_he).idx());
    vf = V.row(lsmesh.from_vertex_handle(checking_he).idx());
    vt = V.row(lsmesh.to_vertex_handle(checking_he).idx());
    found = solve_next_geodesic_point(p0, p1, vf, vt, pnorm, p_end);
    if (found)
    {
        edge_out = checking_he;
        return true;
    }
    return false;
}
// the first segment is p0, p1. the edge we search on is vs, ve.
// to find all the vertices on all edges that form the osculating plane with a specific angle
// when angle is 0 degree, the curve will be a asymptotic, when angle is 90 degree, the curve will be a geodesic,
// but this function does not solve geodesic cases.
// the input angle should be in radian, angle is degree*LSC_PI/180
bool find_osculating_plane_intersection_not_geodesic_p1_is_ver(const Eigen::Vector3d &p0, const Eigen::Vector3d &p1,
                                                               const Eigen::Vector3d &vs, const Eigen::Vector3d &ve,
                                                               const Eigen::Vector3d &pnorm,
                                                               const double angle,
                                                               std::vector<Eigen::Vector3d> &p_end)
{
    p_end.clear();
    Eigen::Vector3d ves = ve - vs;
    Eigen::Vector3d vs1 = vs - p1;
    Eigen::Vector3d v10 = p1 - p0;
    Eigen::Vector3d normal0 = v10.cross(vs1);
    Eigen::Vector3d normal1 = v10.cross(ves);
    double up0 = normal0.dot(pnorm);
    double up1 = normal1.dot(pnorm);
    double left0 = up0 * up0;
    double left1 = 2 * up0 * up1;
    double left2 = up1 * up1;
    double pnorm_length = pnorm.norm();
    double pnormsqr = pnorm_length * pnorm_length;
    double cos_angle = cos(angle);
    double right_coff = pnormsqr * cos_angle * cos_angle;
    double right0 = normal0.dot(normal0);
    double right1 = 2 * normal0.dot(normal1);
    double right2 = normal1.dot(normal1);
    right0 *= right_coff;
    right1 *= right_coff;
    right2 *= right_coff;
    double poly0 = left0 - right0;
    double poly1 = left1 - right1;
    double poly2 = left2 - right2;

    // left = norm.dot(pnorm), right=|norm|*|pnorm|*cos_angle.
    // left^2= right^2

    std::vector<double> equation(3);
    equation[0] = poly0;
    equation[1] = poly1;
    equation[2] = poly2;
    equation = polynomial_simplify(equation);
    assert(equation.size() <= 3); // at most a quadratic function.
    std::array<double, 2> roots;
    bool solved = quadratic_solver(equation, roots);
    if (!solved)
    {
        return false; // did not find such a plane
    }
    for (int i = 0; i < 2; i++)
    {
        double t = roots[i];
        if (t < -MERGE_VERTEX_RATIO || t > 1 + MERGE_VERTEX_RATIO)
        {
            continue;
        }
        std::cout << "real t " << t << std::endl;
        if (t >= -MERGE_VERTEX_RATIO && t < 0)
        {
            t = 0;
        }
        if (t > 1 && t <= 1 + MERGE_VERTEX_RATIO)
        {
            t = 1;
        }
        if (t >= 0 && t <= 1)
        {
            p_end.push_back(vs + t * ves);
            // this is just for debug
            Eigen::Vector3d direction0 = (p1 - p0).normalized();
            Eigen::Vector3d direction1 = (vs + t * ves - p1).normalized();
            Eigen::Vector3d pseudo_norm = direction0.cross(direction1).normalized();
            double dot_product = pseudo_norm.dot(pnorm);
            double cos_sqr = dot_product * dot_product;
            double real_angle = acos(dot_product) * 180 / LSC_PI;

            std::cout << "**INSIDE, cosin^2 = " << cos_sqr << std::endl;
            // std::cout<<"the pseudo_norm "<<pseudo_norm.transpose()<<std::endl;
            // std::cout<<"the pnorm "<<pnorm.transpose()<<std::endl;
            // std::cout<<"dire0 = "<<direction0.transpose()<<std::endl;
            // std::cout<<"dire1 = "<<direction1.transpose()<<std::endl;
            std::cout << "t = " << t << std::endl;
            std::cout << "poly value, " << polynomial_value(equation, t) << std::endl;
            std::cout << "real angle is " << real_angle << std::endl;
            std::cout << "p0 " << p0.transpose() << std::endl;
            std::cout << "p1 " << p1.transpose() << std::endl;
            std::cout << "ve " << ve.transpose() << std::endl;
            std::cout << "vs " << vs.transpose() << std::endl;
            // std::cout<<"computed "<<(vs + t * ves-p1).transpose()<<std::endl;
            std::cout << "p1-p0 " << (p1 - p0).transpose() << std::endl;
            std::cout << "pn-p1 " << (vs + t * ves - p1).transpose() << std::endl;
            std::cout << "ve-vs " << ves.transpose() << std::endl;
            std::cout << "p1p0, vevs relation cos " << (p1 - p0).normalized().dot(ves.normalized()) << std::endl;
            std::cout << "** p1p0, pnp1 cross " << (p1 - p0).cross(vs + t * ves - p1).transpose() << std::endl;

            std::cout << "poly " << poly0 << " " << poly1 << " " << poly2 << std::endl;
            std::cout << "b^2-4ac " << poly1 * poly1 - 4 * poly0 * poly2 << std::endl;

            std::cout << std::endl;
        }
    }

    // debug end
    return p_end.size();
}
// not geodesic
void find_intersection_on_halfedge(const CGMesh& lsmesh, const Eigen::MatrixXd& V, CGMesh::HalfedgeHandle& checking_he, 
const Eigen::Vector3d &p0, const Eigen::Vector3d &p1, const Eigen::Vector3d &pnorm, const double angle,
std::vector<CGMesh::HalfedgeHandle> &edge_out, std::vector<Eigen::Vector3d> &p_end){
    Eigen::Vector3d vs;
    Eigen::Vector3d ve;
    vs = V.row(lsmesh.from_vertex_handle(checking_he).idx());
    ve = V.row(lsmesh.to_vertex_handle(checking_he).idx());
    std::vector<Eigen::Vector3d> temp_pend;
    bool found = find_osculating_plane_intersection_not_geodesic_p1_is_ver(p0, p1, vs, ve, pnorm, angle, temp_pend);
    if (found)
    {
        for (int i = 0; i < temp_pend.size(); i++)
        {
            edge_out.push_back(checking_he);
            p_end.push_back(temp_pend[i]);
        }
    }
}

bool lsTools::find_osculating_plane_intersection_not_geodesic_p1_is_not_ver(const Eigen::Vector3d &p0, const Eigen::Vector3d &p1,
                                                                            const CGMesh::HalfedgeHandle &edge_middle, const Eigen::Vector3d &pnorm, const double angle,
                                                                            std::vector<CGMesh::HalfedgeHandle> &edge_out, std::vector<Eigen::Vector3d> &p_end)
{
    edge_out.clear();
    p_end.clear();
    CGMesh::HalfedgeHandle ophe = lsmesh.opposite_halfedge_handle(edge_middle);
    // check next halfedge
    CGMesh::HalfedgeHandle checking_he = lsmesh.next_halfedge_handle(ophe);
    assert(lsmesh.face_handle(edge_middle).idx() != lsmesh.face_handle(ophe).idx());
    assert(lsmesh.face_handle(ophe).idx() == lsmesh.face_handle(checking_he).idx());
    find_intersection_on_halfedge(lsmesh,V, checking_he,p0,p1,pnorm,angle,edge_out,p_end);
    // check the previous halfedge
    checking_he = lsmesh.prev_halfedge_handle(ophe);
    assert(lsmesh.face_handle(edge_middle).idx() != lsmesh.face_handle(ophe).idx());
    assert(lsmesh.face_handle(ophe).idx() == lsmesh.face_handle(checking_he).idx());
    find_intersection_on_halfedge(lsmesh,V, checking_he,p0,p1,pnorm,angle,edge_out,p_end);
    // the next parts are for the cases the edge_middle is tangent to the curve
    checking_he = lsmesh.next_halfedge_handle(edge_middle);
    find_intersection_on_halfedge(lsmesh,V, checking_he,p0,p1,pnorm,angle,edge_out,p_end);
    checking_he = lsmesh.prev_halfedge_handle(edge_middle);
    find_intersection_on_halfedge(lsmesh,V, checking_he,p0,p1,pnorm,angle,edge_out,p_end);

    std::cout << "lower level check" << std::endl;
    for (int i = 0; i < p_end.size(); i++)
    {
        Eigen::Vector3d direc0 = (p1 - p0).normalized();
        Eigen::Vector3d direc1 = (p_end[i] - p1).normalized();
        Eigen::Vector3d dddnorm = direc0.cross(direc1).normalized();
        double cosangle = pnorm.dot(dddnorm);
        std::cout << "lower, cos^2 " << cosangle * cosangle << std::endl;
    }
}

void extend_vector(std::vector<Eigen::Vector3d> &first, std::vector<Eigen::Vector3d> &second)
{
    first.insert(first.end(), second.begin(), second.end());
}
void extend_vector(std::vector<CGMesh::HalfedgeHandle> &first, std::vector<CGMesh::HalfedgeHandle> &second)
{
    first.insert(first.end(), second.begin(), second.end());
}
template <typename Tp>
void vector_duplicate(const Tp &input, const int size, std::vector<Tp> &result)
{
    result.resize(size);
    for (int i = 0; i < size; i++)
    {
        result[i] = input;
    }
}
bool angles_match(const double angle_degree1, const double angle_degree2)
{
    double uniformed = angle_degree1;
    if (angle_degree1 > 180)
    {
        uniformed -= 180;
    }
    if (angle_degree1 < 0)
    {
        uniformed += 180;
    }
    if (fabs(uniformed - angle_degree2) < ANGLE_TOLERANCE)
    {
        return true;
    }
    return false;
}

// caution: this is not an initialization to find the first and second point
// TODO here is the bug: We should not take curve, but take pseudo-curve
// 1. the bi-normal should not change much. PARAMETER [correct_normal]
// 2. the angle should not be too far away from the target angle (due to the numerical error) PARAMETER [correct_angle]
// 3. the curve should be smooth. No sharp turns
void pseudo_geodesic_intersection_filter_by_closeness(
    const std::vector<Eigen::Vector3d> &curve, const std::array<Eigen::Vector3d,3>& pcurve_local, const double angle_degree,
    const std::vector<Eigen::Vector3d> &pnorm_list, const Eigen::Vector3d &pnorm,
    const std::vector<Eigen::Vector3d> &candi_points, int &id)
{
    assert(curve.size() >= 2); // take previous points to find the next one
    // assert(candi_points.size() > 0);
    if (candi_points.size() == 0)
    {
        id = -1;
    }
    if (candi_points.size() == 1)
    {
        id = 0;
        return;
    }
    int size = curve.size();
    Eigen::Vector3d dire1 = (pcurve_local[2] - pcurve_local[1]).normalized();
    std::vector<bool> flag_right;
    int candi_size = candi_points.size();
    flag_right.resize(candi_size);
    for (int i = 0; i < candi_size; i++)
    {
        flag_right[i] = true;
    }
    if (size >= 3)
    { // when there is an osculating plane in the previous point
        Eigen::Vector3d dire0 = (pcurve_local[1] - pcurve_local[0]).normalized();
        Eigen::Vector3d osnormal = dire0.cross(dire1);
        double dot_product0 = pnorm_list[size - 3].dot(osnormal);
        for (int i = 0; i < candi_size; i++)
        {
            Eigen::Vector3d dire2 = (candi_points[i] - pcurve_local[2]).normalized();
            Eigen::Vector3d osnormal1 = dire1.cross(dire2).normalized();
            double dot_product1 = pnorm.dot(osnormal1);
            bool correct_normal = dot_product0 * dot_product1 > 0; // remove negative normals
            double real_angle = acos(dot_product1) * 180 / LSC_PI;
            bool correct_angle = angles_match(real_angle, angle_degree); // remove planar segments
            flag_right[i] = correct_normal && correct_angle;             // check if the signs are the same. if so, set as true;
            std::cout << "correct normal " << correct_normal << std::endl;
            std::cout << "correct angle " << correct_angle << std::endl;
            std::cout << "real angle, " << real_angle << ", angle_degree, " << angle_degree << std::endl;
            std::cout << "dot_product1, " << dot_product1 << std::endl;
        }
    }

    int closest_id = -1;
    // TODO should we consider the angle consistancy? to avoid cases like 5 degree flip to 355 degree?
    // int closest_id_consider_angle=-1;
    double closest_angle_diff_radian = -5;
    // double closest_angle_diff_radian_consider_angle = 370. * LSC_PI / 180.;
    // int quadrant; // show which quadrant the normal is in respect to the pnormal
    // if(angle_degree>=0&& angle_degree<=90){
    //   quadrant=1;
    // }
    // if(angle_degree>=90&& angle_degree<=180){
    //   quadrant=2;
    // }
    // if(angle_degree>=180&& angle_degree<=270){
    //   quadrant=3;
    // }
    // if(angle_degree>=270&& angle_degree<=360){
    //   quadrant=4;
    // }
    double angle_radian = angle_degree * LSC_PI / 180;
    for (int i = 0; i < candi_points.size(); i++)
    {
        if (!flag_right[i])
            continue;
        Eigen::Vector3d dire2 = (candi_points[i] - pcurve_local[2]).normalized();
        
        double dot_product = dire1.dot(dire2);
        if(dot_product<0){// avoid sharp turn
            continue;
        }
        if (dot_product > closest_angle_diff_radian) // select the most smooth one
        {
            closest_angle_diff_radian = dot_product;
            closest_id = i;
        }
    }
    // already found the vertex that forms the angle most close to the given one
    id = closest_id;
    return;
}
void initial_segment_intersection_filter_by_closeness(
    const double angle_degree,
    const Eigen::Vector3d &start_point,
    const Eigen::Vector3d &reference_direction,
    const std::vector<Eigen::Vector3d> &candi_points, int &id)
{
    double angle_radian = angle_degree * LSC_PI / 180;
    assert(candi_points.size() > 0);
    if (candi_points.size() == 1)
    {
        id = 0;
        return;
    }
    int closest_id = -1;
    double closest_angle_diff_radian = 370. * LSC_PI / 180.;

    for (int i = 0; i < candi_points.size(); i++)
    {
        Eigen::Vector3d direction = (candi_points[i] - start_point).normalized();
        double inner_product = direction.dot(reference_direction);
        double angle_tmp_radian = acos(inner_product);
        double angle_diff = std::min(fabs(angle_tmp_radian - angle_radian), fabs(2 * LSC_PI - angle_tmp_radian - angle_radian));
        if (angle_diff < closest_angle_diff_radian)
        {
            closest_angle_diff_radian = angle_diff;
            closest_id = i;
        }
    }
    // already found the vertex that forms the angle most close to the given one
    id = closest_id;
    return;
}

// please initialize quadricCalculator before tracing a single step
// the idea is to first check if point_middle is a vertex.
// if yes, then search in the one-ring edges to find the candidate intersections;
// otherwise, search the other two edges of the neighbouring face.
// after finding all the candidates, find the best point, which satisfies:
// 1. form an obtuse angle with the previous segment; 2. choose the one that form the angle most similar with the previous segment
// TODO should also consider about the inner product between the binormal and normal
// TODO maybe we consider the start point of each step as the pseudo-vertex of the last step?
// TODO since the intersections cannot be found only when tracing asymptotics, we can set a tolerance for enable pseudo-vertices.
// angle is in degree, when degree is 0, is a geodesic. when degree is 90, it is an asymptotic
bool lsTools::get_pseudo_vertex_and_trace_forward(
    QuadricCalculator &cc,
    const std::vector<Eigen::Vector3d> &curve, std::array<Eigen::Vector3d,3>& pcurve_local, const double angle_degree,
    const CGMesh::HalfedgeHandle &edge_middle,
    const Eigen::Vector3d &point_in, const Eigen::Vector3d &point_middle,
    const bool calculate_pseudo_vertex, CGMesh::HalfedgeHandle &edge_out,
    Eigen::Vector3d &point_out, bool &generate_pseudo_vertex, Eigen::Vector3d &pseudo_vertex_out)
{
    std::cout << "inside trace forward, curve size = " << curve.size() << std::endl;
    unsigned radius = 2;
    bool is_geodesic = false;
    pseudo_vertex_out = point_middle; // by default the pseudo vertex is the point middle
    generate_pseudo_vertex = false;
    double angle_radian = angle_degree * LSC_PI / 180.; // the angle in radian
    // std::vector<CGMesh::halfedge_handle>
    std::vector<Eigen::Vector3d> candidate_pts;
    std::vector<CGMesh::HalfedgeHandle> candidate_handles;
    // Precomputation
    // cc.init(V, F);// TODO move init out of this function

    cc.kRing = radius;
    cc.st = K_RING_SEARCH;

    // start to calculate the pseudo-vertex
    int fid1 = lsmesh.face_handle(edge_middle).idx();
    if (lsmesh.is_boundary(lsmesh.from_vertex_handle(edge_middle)) && lsmesh.is_boundary(lsmesh.to_vertex_handle(edge_middle)))
    {
        // this is a boundary edge on the boundary loop.
        return false;
    }
    // deal with the case where the point_middle is a mesh vertex
    int ver_from_id = lsmesh.from_vertex_handle(edge_middle).idx();
    int ver_to_id = lsmesh.to_vertex_handle(edge_middle).idx();
    Eigen::Vector3d ver_from = V.row(ver_from_id);
    Eigen::Vector3d ver_to = V.row(ver_to_id);

    double dist_from = (point_middle - ver_from).norm();
    double dist_to = (point_middle - ver_to).norm();
    double dist_total = dist_from + dist_to;
    double from_ratio = dist_from / dist_total;
    double to_ratio = dist_to / dist_total;
    if (angle_degree > 90 - ANGLE_TOLERANCE && angle_degree < 90 + ANGLE_TOLERANCE)
    { // it means it is a geodesic
        is_geodesic = true;
        // std::cout << "**tracing a geodesic" << std::endl;
    }
    if (from_ratio <= MERGE_VERTEX_RATIO || to_ratio <= MERGE_VERTEX_RATIO) // the point_middle is a vertex
    {
        std::cout << "\nnext middle is ver" << std::endl;
        CGMesh::VertexHandle center_handle;
        if (from_ratio <= MERGE_VERTEX_RATIO)
        {
            center_handle = lsmesh.from_vertex_handle(edge_middle);
        }
        else
        {
            center_handle = lsmesh.to_vertex_handle(edge_middle);
        }
        Eigen::Vector3d pnorm = norm_v.row(center_handle.idx());

        for (CGMesh::VertexOHalfedgeIter voh_itr = lsmesh.voh_begin(center_handle);
             voh_itr != lsmesh.voh_end(center_handle); ++voh_itr)
        {
            auto heh = voh_itr.handle();
            CGMesh::HalfedgeHandle edge_to_check = lsmesh.next_halfedge_handle(heh);
            assert(lsmesh.from_vertex_handle(edge_to_check).idx() != center_handle.idx());
            assert(lsmesh.to_vertex_handle(edge_to_check).idx() != center_handle.idx());
            Eigen::Vector3d vs = V.row(lsmesh.from_vertex_handle(edge_to_check).idx());
            Eigen::Vector3d ve = V.row(lsmesh.to_vertex_handle(edge_to_check).idx());

            if (is_geodesic)
            { // it means it is a geodesic
                bool found = find_geodesic_intersection_p1_is_ver(point_in, point_middle, vs, ve, pnorm, point_out);
                if (found)
                {
                    edge_out = edge_to_check;
                    return true;
                }
            } // geodesic end
            else
            {
                std::vector<Eigen::Vector3d> tmp_candidates;
                std::vector<CGMesh::HalfedgeHandle> tmp_hehs;
                bool found = find_osculating_plane_intersection_not_geodesic_p1_is_ver(point_in, point_middle, vs, ve, pnorm, angle_radian, tmp_candidates);
                if (found)
                {
                    extend_vector(candidate_pts, tmp_candidates);
                    vector_duplicate(edge_to_check, tmp_candidates.size(), tmp_hehs);
                    extend_vector(candidate_handles, tmp_hehs);
                }
            }
        } // end of searching all the one-ring edges
        if (is_geodesic)
        {
            return false;
        }
        else
        {

            // std::cout<<"before picking the points "<<std::endl;
            // pick the point
            int id;
            if (flag_dbg)
            {
                // std::cout << "we are debugging" << std::endl;
                if (curve.size() - 1 == id_dbg)
                {
                    std::cout << "+++++++++++++++++++++++++++++++++++\nCHECKING START vertex" << std::endl;
                    // std::cout << "output debug vertices" << std::endl;
                    ver_dbg1.resize(1, 3);
                    ver_dbg1.row(0) = point_middle;
                    ver_dbg = vec_list_to_matrix(candidate_pts);
                    flag_dbg = false;
                    std::cout << "+++++++++++++++++++++++++++++++++++\nCHECKING END" << std::endl;
                }
            }
            std::cout << "before filtering, size " << candidate_pts.size() << std::endl;
            pseudo_geodesic_intersection_filter_by_closeness(curve,pcurve_local, angle_degree, pnorm_list_dbg, pnorm, candidate_pts, id);
            if (id < 0)
            {
                // TODO we may not need to re-compute pver, since we are on vertices
                if (calculate_pseudo_vertex)
                {
                    std::cout << "TTTTT after computing pseudo vertex, still no result " << std::endl;
                    return false;
                }
                else
                {
                    // no need to update pcurve_local[2] since we don't compute pver on vertices.
                    bool recalculate_result = get_pseudo_vertex_and_trace_forward(
                        cc,
                        curve, pcurve_local,angle_degree,
                        edge_middle,
                        point_in, point_middle,
                        true, edge_out,
                        point_out, generate_pseudo_vertex, pseudo_vertex_out);
                    return recalculate_result;
                }
            }

            pnorm_list_dbg.push_back(pnorm);
            edge_out = candidate_handles[id];
            point_out = candidate_pts[id];
            std::cout << "**Higher level checking" << std::endl;
            Eigen::Vector3d direc0 = (point_middle - point_in).normalized();
            Eigen::Vector3d direc1 = (point_out - point_middle).normalized();
            Eigen::Vector3d dddnorm = (direc0.cross(direc1)).normalized();
            double tmpcos = pnorm.dot(dddnorm);
            std::cout << "cos^2 = " << tmpcos * tmpcos << std::endl;
            std::cout << std::endl;
            return true;
        }
    }

    // if the middle_point is on the edge
    int fid2 = lsmesh.opposite_face_handle(edge_middle).idx();
    assert(fid1 != fid2);
    int cent_id = lsmesh.from_vertex_handle(edge_middle).idx();
    Eigen::Vector3d norm1 = norm_f.row(fid1);
    Eigen::Vector3d norm2 = norm_f.row(fid2);
    Eigen::Vector3d pnorm = (norm1 + norm2).normalized(); // the normal of the pseudo-vertex

    Eigen::Vector3d pver;
    if (calculate_pseudo_vertex)
    {
        std::cout << "TTTTT getting pseudo vertex " << std::endl;
        cc.get_pseudo_vertex_on_edge(cent_id, point_middle, pnorm, pver);
        generate_pseudo_vertex = true; // TODO add tolerance
    }
    else
    {
        pver = point_middle;
    }

    std::cout << "    pseudo vertex diff " << (point_middle - pver).norm() << std::endl;
    
    pseudo_vertex_out = pver;
    pcurve_local[2]=pver;
    if (is_geodesic)
    {
        bool found = find_geodesic_intersection_p1_is_NOT_ver(point_in, pver, edge_middle, pnorm, edge_out, point_out);
        return found;
    }
    else
    {
        bool found = find_osculating_plane_intersection_not_geodesic_p1_is_not_ver(point_in, pver, edge_middle, pnorm,
                                                                                   angle_radian, candidate_handles, candidate_pts);

        // pick the point
        int id;
        if (flag_dbg)
        {
            // std::cout << "we are debugging" << std::endl;
            if (curve.size() - 1 == id_dbg)
            {
                std::cout << "+++++++++++++++++++++++++++++++++++\nCHECKING START" << std::endl;
                ver_dbg1.resize(1, 3);
                ver_dbg1.row(0) = point_middle;
                ver_dbg = vec_list_to_matrix(candidate_pts);
                flag_dbg = false;
                std::cout << "+++++++++++++++++++++++++++++++++++\nCHECKING END" << std::endl;
            }
        }
        std::cout << "before filtering, size " << candidate_pts.size() << std::endl;
        pseudo_geodesic_intersection_filter_by_closeness(curve, pcurve_local,angle_degree, pnorm_list_dbg, pnorm, candidate_pts, id);
        if (id < 0)
        {
            if (calculate_pseudo_vertex)
            {
                std::cout << "TTTTT after computing pseudo vertex, still no result " << std::endl;
                return false;
            }
            else
            {
                
                bool recalculate_result = get_pseudo_vertex_and_trace_forward(
                    cc,
                    curve,pcurve_local, angle_degree,
                    edge_middle,
                    point_in, point_middle,
                    true, edge_out,
                    point_out, generate_pseudo_vertex, pseudo_vertex_out);
                return recalculate_result;
            }
        }
        pnorm_list_dbg.push_back(pnorm);
        edge_out = candidate_handles[id];
        point_out = candidate_pts[id];
        std::cout << "**Higher level checking" << std::endl;
        Eigen::Vector3d direc0 = (pver - point_in).normalized();
        Eigen::Vector3d direc1 = (point_out - pver).normalized();
        Eigen::Vector3d dddnorm = (direc0.cross(direc1)).normalized();
        double tmpcos = pnorm.dot(dddnorm);
        std::cout << "cos^2 = " << tmpcos * tmpcos << std::endl;
        std::cout << std::endl;
        return true;
    }
    return false;
}
int conservative_sign(const double input)
{
    if (input > QUADRANT_TOLERANCE)
    {
        return 1;
    }
    if (input < -QUADRANT_TOLERANCE)
    {
        return -1;
    }
    return 0;
}

bool angle_satisfy_quadrant(const double angle_degree, const Eigen::Vector3d &reference_direction,
                            const Eigen::Vector3d &normal, const Eigen::Vector3d &query_direction)
{
    double angle_radian = angle_degree * LSC_PI / 180.;
    Eigen::Vector3d yframe = normal.cross(reference_direction).normalized();
    double x = query_direction.dot(reference_direction);
    double y = query_direction.dot(yframe);
    int x_sign = conservative_sign(x);
    int y_sign = conservative_sign(y);
    bool quadrant[4];
    quadrant[0] = 0;
    quadrant[1] = 0;
    quadrant[2] = 0;
    quadrant[3] = 0;
    if (x_sign >= 0 && y_sign >= 0)
    {
        quadrant[0] = 1;
    }
    if (x_sign <= 0 && y_sign >= 0)
    {
        quadrant[1] = 1;
    }
    if (x_sign <= 0 && y_sign <= 0)
    {
        quadrant[2] = 1;
    }
    if (x_sign >= 0 && y_sign <= 0)
    {
        quadrant[3] = 1;
    }
    int target_quadrant;
    if (angle_degree >= 0 && angle_degree <= 90)
    {
        target_quadrant = 0;
    }
    if (angle_degree >= 90 && angle_degree <= 180)
    {
        target_quadrant = 1;
    }
    if (angle_degree >= 180 && angle_degree <= 270)
    {
        target_quadrant = 2;
    }
    if (angle_degree >= 270 && angle_degree <= 360)
    {
        target_quadrant = 3;
    }
    if (quadrant[target_quadrant])
    {
        return true;
    }
    return false;
}
void pseudo_geodesic_intersection_satisfy_angle_quadrant(
    const double angle_degree, const Eigen::Vector3d &reference_direction,
    const Eigen::Vector3d &start_point,
    const Eigen::Vector3d &normal, const std::vector<CGMesh::HalfedgeHandle> &handle_in,
    const std::vector<Eigen::Vector3d> &point_in, std::vector<CGMesh::HalfedgeHandle> &handle_out,
    std::vector<Eigen::Vector3d> &point_out)
{
    for (int i = 0; i < point_in.size(); i++)
    {
        Eigen::Vector3d direction = (point_in[i] - start_point).normalized();
        if (angle_satisfy_quadrant(angle_degree, reference_direction, normal, direction))
        {
            handle_out.push_back(handle_in[i]);
            point_out.push_back(point_in[i]);
        }
    }
}
// this code extend the candidate point list.
// after calling this code, filter the points according to the quadrant and closeness.
bool find_initial_direction_intersection_on_edge(const Eigen::Vector3d &start_point, const double angle_degree,
                                                 const Eigen::Vector3d &reference_direction, const Eigen::Vector3d &vs, const Eigen::Vector3d &ve,
                                                 const Eigen::Vector3d &normal, const CGMesh::HalfedgeHandle &handle_in, std::vector<CGMesh::HalfedgeHandle> &handle_out,
                                                 std::vector<Eigen::Vector3d> &point_out)
{
    double angle_radian = angle_degree * LSC_PI / 180.;
    Eigen::Vector3d ves = ve - vs;
    Eigen::Vector3d vssp = vs - start_point;
    double degree1 = ves.dot(reference_direction);
    double degree0 = vssp.dot(reference_direction);
    double left_square_0 = degree0 * degree0;
    double left_square_1 = 2 * degree0 * degree1;
    double left_square_2 = degree1 * degree1;
    double right_square_0 = vssp.dot(vssp);
    double right_square_1 = 2 * vssp.dot(ves);
    double right_square_2 = ves.dot(ves);
    double right_cons = cos(angle_degree) * reference_direction.dot(reference_direction);
    right_square_0 *= right_cons;
    right_square_1 *= right_cons;
    right_square_2 *= right_cons;
    double eq0 = left_square_0 - right_square_0;
    double eq1 = left_square_1 - right_square_1;
    double eq2 = left_square_2 - right_square_2;
    std::vector<double> equation(3);
    equation[0] = eq0;
    equation[1] = eq1;
    equation[2] = eq2;
    equation = polynomial_simplify(equation);
    std::array<double, 2> roots;
    bool found = quadratic_solver(equation, roots);
    std::vector<Eigen::Vector3d> p_end;
    for (int i = 0; i < 2; i++)
    {
        double t = roots[i];
        if (t < -MERGE_VERTEX_RATIO || t > 1 + MERGE_VERTEX_RATIO)
        {
            continue;
        }
        if (t >= -MERGE_VERTEX_RATIO && t < 0)
        {
            t = 0;
        }
        if (t > 1 && t <= 1 + MERGE_VERTEX_RATIO)
        {
            t = 1;
        }
        if (t >= 0 && t <= 1)
        {
            p_end.push_back(vs + t * ves);
        }
    }
    if (p_end.size() == 0)
    {
        return false;
    }
    for (int i = 0; i < p_end.size(); i++)
    {
        handle_out.push_back(handle_in);
        point_out.push_back(p_end[i]);
    }
    return true;
}
bool lsTools::init_pseudo_geodesic_first_segment(
    const CGMesh::HalfedgeHandle &start_boundary_edge_pre, const double &start_point_para,
    const double start_boundary_angle_degree,
    CGMesh::HalfedgeHandle &intersected_handle,
    Eigen::Vector3d &intersected_point)
{
    CGMesh::HalfedgeHandle start_boundary_edge;
    if (lsmesh.face_handle(start_boundary_edge_pre).idx() < 0)
    {
        start_boundary_edge = lsmesh.opposite_halfedge_handle(start_boundary_edge_pre);
    }
    else
    {
        start_boundary_edge = start_boundary_edge_pre;
    }
    Eigen::Vector3d vf = V.row(lsmesh.from_vertex_handle(start_boundary_edge).idx());
    Eigen::Vector3d vt = V.row(lsmesh.to_vertex_handle(start_boundary_edge).idx());
    Eigen::Vector3d start_point = vf + start_point_para * (vt - vf);
    Eigen::Vector3d reference_direction = (vt - vf).normalized();
    std::vector<CGMesh::HalfedgeHandle> handle_out, handle_out2;
    std::vector<Eigen::Vector3d> point_out, point_out2;

    if (start_point_para == 0 || start_point_para == 1)
    {
        CGMesh::VertexHandle center_handle;
        if (start_point_para == 0)
        {
            center_handle = lsmesh.from_vertex_handle(start_boundary_edge);
        }
        else
        {
            center_handle = lsmesh.to_vertex_handle(start_boundary_edge);
        }
        Eigen::Vector3d normal = V.row(center_handle.idx());
        for (CGMesh::VertexOHalfedgeIter voh_itr = lsmesh.voh_begin(center_handle);
             voh_itr != lsmesh.voh_end(center_handle); ++voh_itr)
        {

            auto heh = voh_itr.handle();
            CGMesh::HalfedgeHandle edge_to_check = lsmesh.next_halfedge_handle(heh);
            assert(lsmesh.from_vertex_handle(edge_to_check).idx() != center_handle.idx());
            assert(lsmesh.to_vertex_handle(edge_to_check).idx() != center_handle.idx());
            Eigen::Vector3d vs = V.row(lsmesh.from_vertex_handle(edge_to_check).idx());
            Eigen::Vector3d ve = V.row(lsmesh.to_vertex_handle(edge_to_check).idx());
            find_initial_direction_intersection_on_edge(start_point, start_boundary_angle_degree, reference_direction,
                                                        vs, ve, normal, edge_to_check, handle_out, point_out);
        }
    }

    // check next halfedge
    int fid = lsmesh.face_handle(start_boundary_edge).idx();
    Eigen::Vector3d normal = norm_f.row(fid);
    CGMesh::HalfedgeHandle checking_he = lsmesh.next_halfedge_handle(start_boundary_edge);
    Eigen::Vector3d vs;
    Eigen::Vector3d ve;
    assert(lsmesh.face_handle(start_boundary_edge).idx() == lsmesh.face_handle(checking_he).idx());
    vs = V.row(lsmesh.from_vertex_handle(checking_he).idx());
    ve = V.row(lsmesh.to_vertex_handle(checking_he).idx());
    find_initial_direction_intersection_on_edge(start_point, start_boundary_angle_degree, reference_direction,
                                                vs, ve, normal, checking_he, handle_out, point_out);
    // check the previous halfedge
    checking_he = lsmesh.prev_halfedge_handle(start_boundary_edge);
    assert(lsmesh.face_handle(start_boundary_edge).idx() == lsmesh.face_handle(checking_he).idx());
    vs = V.row(lsmesh.from_vertex_handle(checking_he).idx());
    ve = V.row(lsmesh.to_vertex_handle(checking_he).idx());
    find_initial_direction_intersection_on_edge(start_point, start_boundary_angle_degree, reference_direction,
                                                vs, ve, normal, checking_he, handle_out, point_out);
    pseudo_geodesic_intersection_satisfy_angle_quadrant(start_boundary_angle_degree, reference_direction, start_point,
                                                        normal, handle_out, point_out, handle_out2, point_out2);
    if (point_out2.size() == 0)
    {
        std::cout << "ERROR: no point satisfies the angle quadrant in initialization" << std::endl;
        return false;
    }
    int result_id = -1;
    initial_segment_intersection_filter_by_closeness(start_boundary_angle_degree, start_point,
                                                     reference_direction, point_out2, result_id);
    intersected_handle = handle_out2[result_id];
    intersected_point = point_out2[result_id];
    return true;
}
void update_tracing_list(const Eigen::Vector3d &ver, const Eigen::Vector3d &pver, const CGMesh::HalfedgeHandle &handle,
                         std::vector<Eigen::Vector3d> &verlist, std::vector<Eigen::Vector3d> &pverlist, std::vector<CGMesh::HalfedgeHandle> &handle_list)
{
    verlist.push_back(ver);
    pverlist.push_back(pver);
    handle_list.push_back(handle);
}
Eigen::Vector3d the_first_point(const Eigen::Vector3d &vs, const Eigen::Vector3d &vt, const double para)
{
    return vs + (vt - vs) * para;
}
bool lsTools::trace_single_pseudo_geodesic_curve(const double target_angle_degree,
                                                 const CGMesh::HalfedgeHandle &start_boundary_edge, const double &start_point_para,
                                                 const double start_boundary_angle_degree,
                                                 std::vector<Eigen::Vector3d> &curve)
{
    curve.clear();
    CGMesh::HalfedgeHandle intersected_handle_tmp;
    Eigen::Vector3d intersected_point_tmp;
    bool found = init_pseudo_geodesic_first_segment(start_boundary_edge, start_point_para, start_boundary_angle_degree,
                                                    intersected_handle_tmp, intersected_point_tmp);
    std::cout << "INITIALIZARION DONE!!!" << std::endl;
    if (!found)
    {
        std::cout << "error in initialization the boundary direction" << std::endl;
        return false;
    }
    QuadricCalculator cc;
    std::vector<Eigen::Vector3d> pcurve;
    std::vector<CGMesh::HalfedgeHandle> handles;
    Eigen::Vector3d first_point = the_first_point(
        V.row(lsmesh.from_vertex_handle(start_boundary_edge).idx()),
        V.row(lsmesh.to_vertex_handle(start_boundary_edge).idx()),
        start_point_para);

    update_tracing_list(first_point, first_point, start_boundary_edge, curve, pcurve, handles);
    curve.push_back(intersected_point_tmp);
    handles.push_back(intersected_handle_tmp);
    cc.init(V, F); // TODO move this out of this function
    std::cout << "target angle " << target_angle_degree << std::endl;
    for (int i = 0;; i++)
    {
        std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
        CGMesh::HalfedgeHandle edge_out;
        Eigen::Vector3d point_out, pseudo_vertex_out;
        bool generate_pseudo_vertex;
        std::array<Eigen::Vector3d,3> pcurve_local; // the local pseudo-curve
        pcurve_local[2]=curve[curve.size()-1];// temporarily set it as the end point of the curve
        pcurve_local[1]=pcurve[curve.size()-2];// this is the end point of the pcurve.
        if(curve.size()>=3){
            pcurve_local[0]=pcurve[curve.size()-3];
        }
        bool calculate_pseudo_vertex_default=true;
        found = get_pseudo_vertex_and_trace_forward(cc, curve, pcurve_local, target_angle_degree, intersected_handle_tmp,
                                                    first_point,
                                                    intersected_point_tmp, calculate_pseudo_vertex_default, edge_out, point_out,
                                                    generate_pseudo_vertex, pseudo_vertex_out);
        if (found)
        {

            update_tracing_list(point_out, pseudo_vertex_out, edge_out, curve, pcurve, handles);

            intersected_handle_tmp = edge_out; // middle edge = edge out
            first_point = pseudo_vertex_out;   // first point = pseudo middle point
            intersected_point_tmp = point_out; // point middle = point out. but may eventually be converted to a pseudo point
        }
        else
        {
            // the last pseudo point is the middle point of this step, aka the out point of last step
            pcurve.push_back(intersected_point_tmp);
            break;
        }
    }
    pseudo_vers_dbg = pcurve;
    assert(pcurve.size() == curve.size());
    return true;
}
void lsTools::show_pseudo_geodesic_curve(std::vector<Eigen::MatrixXd> &E0, std::vector<Eigen::MatrixXd> &E1, Eigen::MatrixXd &vers)
{
    E0.resize(trace_vers.size());
    E1.resize(trace_vers.size());
    int pnbr=0;
    for (int i = 0; i < trace_vers.size(); i++)
    {
        int vsize = trace_vers[i].size();
        E0[i].resize(vsize - 1, 3);
        E1[i].resize(vsize - 1, 3);
        for (int j = 0; j < vsize - 1; j++)
        {
            E0[i].row(j) = trace_vers[i][j];
            E1[i].row(j) = trace_vers[i][j + 1];
            pnbr++;
        }
        pnbr++;
        
    }
    vers.resize(pnbr,3);
    pnbr=0;
    for(int i=0;i<trace_vers.size();i++){
        for(int j=0;j<trace_vers[i].size();j++){
            vers.row(pnbr)=trace_vers[i][j];
            pnbr++;
        }
    }
}