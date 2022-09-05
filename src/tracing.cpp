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
                                                               const Eigen::Vector3d &vs, const Eigen::Vector3d &ve, const Eigen::Vector3d &pnorm,
                                                               const double angle, std::vector<Eigen::Vector3d> &p_end)
{
    p_end.clear();
    Eigen::Vector3d ves = ve - vs;
    Eigen::Vector3d vs1 = vs - p1;
    Eigen::Vector3d v01 = p0 - p1;

    // normal director of the plane is norm_degree1*t+norm_degree0
    Eigen::Vector3d norm_degree1 = ves.cross(v01);
    Eigen::Vector3d norm_degree0 = vs1.cross(v01);

    double cos_angle = cos(angle);
    std::vector<double> poly_norm_x(2), poly_norm_y(2), poly_norm_z(2);
    poly_norm_x[0] = norm_degree0[0];
    poly_norm_y[0] = norm_degree0[1];
    poly_norm_z[0] = norm_degree0[2];
    poly_norm_x[1] = norm_degree1[0];
    poly_norm_y[1] = norm_degree1[1];
    poly_norm_z[1] = norm_degree1[2];
    // TODO we can get the function left easily by cross products and dot products. we leave this to fucture.
    // left = norm.dot(pnorm), right=|norm|*|pnorm|*cos_angle.
    // left^2= right^2
    std::vector<double> left;
    left = polynomial_times(poly_norm_x, pnorm[0]);
    left = polynomial_add(left, polynomial_times(poly_norm_y, pnorm[1]));
    left = polynomial_add(left, polynomial_times(poly_norm_z, pnorm[2]));
    left = polynomial_times(left, left);
    std::vector<double> right;
    right = polynomial_times(poly_norm_x, poly_norm_x);
    right = polynomial_add(right, polynomial_times(poly_norm_y, poly_norm_y));
    right = polynomial_add(right, polynomial_times(poly_norm_z, poly_norm_z));
    double pnorm_square = pnorm.norm();
    pnorm_square = pnorm_square * pnorm_square;
    right = polynomial_times(right, pnorm_square * cos_angle * cos_angle);
    std::vector<double> right_inverse = polynomial_times(right, -1.);
    std::vector<double> equation = polynomial_add(left, right_inverse);
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
            Eigen::Vector3d direction0=(p1-p0).normalized();
            Eigen::Vector3d direction1=(vs + t * ves-p1).normalized();
            Eigen::Vector3d pseudo_norm=direction0.cross(direction1).normalized();
            double dot_product=pseudo_norm.dot(pnorm);
            double cos_sqr=dot_product*dot_product;

            
            std::cout << "**INSIDE, cosin^2 = "<<cos_sqr<<std::endl;
            std::cout<<"the pseudo_norm "<<pseudo_norm.transpose()<<std::endl;
            std::cout<<"the pnorm "<<pnorm.transpose()<<std::endl;
        }
    }

    // debug end
    return p_end.size();
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
    Eigen::Vector3d vs;
    Eigen::Vector3d ve;
    assert(lsmesh.face_handle(edge_middle).idx() != lsmesh.face_handle(ophe).idx());
    assert(lsmesh.face_handle(ophe).idx() == lsmesh.face_handle(checking_he).idx());
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
    // check the previous halfedge
    checking_he = lsmesh.prev_halfedge_handle(ophe);
    assert(lsmesh.face_handle(edge_middle).idx() != lsmesh.face_handle(ophe).idx());
    assert(lsmesh.face_handle(ophe).idx() == lsmesh.face_handle(checking_he).idx());
    vs = V.row(lsmesh.from_vertex_handle(checking_he).idx());
    ve = V.row(lsmesh.to_vertex_handle(checking_he).idx());
    found = find_osculating_plane_intersection_not_geodesic_p1_is_ver(p0, p1, vs, ve, pnorm, angle, temp_pend);
    if (found)
    {
        for (int i = 0; i < temp_pend.size(); i++)
        {
            edge_out.push_back(checking_he);
            p_end.push_back(temp_pend[i]);
        }
    }
    return p_end.size();
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

// caution: this is not an initialization to find the first and second point
void pseudo_geodesic_intersection_filter_by_closeness(
    const std::vector<Eigen::Vector3d> &curve, const double angle_degree,
    const Eigen::Vector3d &pnorm,
    const std::vector<Eigen::Vector3d> &candi_points, int &id)
{
    assert(curve.size() >= 2); // take previous points to find the next one
    assert(candi_points.size() > 0);
    if (candi_points.size() == 1)
    {
        id = 0;
        return;
    }
    int size = curve.size();
    Eigen::Vector3d dire1 = (curve[size - 1] - curve[size - 2]).normalized();
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
        Eigen::Vector3d dire2 = (candi_points[i] - curve[size - 1]).normalized();
        // Eigen::Vector3d normal = dire1.cross(dire2).normalized();
        // double normal_pnormal = normal.dot(pnorm);
        // double dire1_normal = dire1.dot(normal);
        // double angle_tmp_radian = acos(normal_pnormal);
        // double angle_diff = std::min(fabs(angle_tmp_radian - angle_radian), fabs(2 * LSC_PI - angle_tmp_radian - angle_radian));
        double dot_product=dire1.dot(dire2);
        if (dot_product > closest_angle_diff_radian)// select the most smooth one
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
    const std::vector<Eigen::Vector3d> &curve, const double angle_degree,
    const CGMesh::HalfedgeHandle &edge_middle,
    const Eigen::Vector3d &point_in, const Eigen::Vector3d &point_middle, CGMesh::HalfedgeHandle &edge_out,
    Eigen::Vector3d &point_out, bool &generate_pseudo_vertex, Eigen::Vector3d &pseudo_vertex_out)
{
    std::cout<<"inside trace forward, curve size = "<<curve.size()<<std::endl;
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
        pnorm_list_dbg.push_back(pnorm);
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
            if (candidate_pts.size() == 0)
            {
                std::cout<<"ERROR did not find any candidate points"<<std::endl;
                return false;
            }
            // std::cout<<"before picking the points "<<std::endl;
            // pick the point
            int id;
            if (flag_dbg)
            {
                // std::cout << "we are debugging" << std::endl;
                if (curve.size() == id_dbg)
                {
                    // std::cout << "output debug vertices" << std::endl;
                    ver_dbg1.resize(1, 3);
                    ver_dbg1.row(0) = point_middle;
                    ver_dbg = vec_list_to_matrix(candidate_pts);
                    flag_dbg = false;
                }
            }
            
            pseudo_geodesic_intersection_filter_by_closeness(curve, angle_degree, pnorm, candidate_pts, id);
            edge_out = candidate_handles[id];
            point_out = candidate_pts[id];
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
    pnorm_list_dbg.push_back(pnorm);
    Eigen::Vector3d pver;
    cc.get_pseudo_vertex_on_edge(cent_id, point_middle, pnorm, pver);
    generate_pseudo_vertex = true; // TODO add tolerance
    pseudo_vertex_out = pver;
    if (is_geodesic)
    {
        bool found = find_geodesic_intersection_p1_is_NOT_ver(point_in, pver, edge_middle, pnorm, edge_out, point_out);
        return found;
    }
    else
    {
        bool found = find_osculating_plane_intersection_not_geodesic_p1_is_not_ver(point_in, pver, edge_middle, pnorm,
                                                                                   angle_radian, candidate_handles, candidate_pts);
        if (candidate_pts.size() == 0)
        {
            return false;
        }
        // pick the point
        int id;
        if (flag_dbg)
        {
            // std::cout << "we are debugging" << std::endl;
            if (curve.size() == id_dbg)
            {
                // std::cout << "output debug vertices" << std::endl;
                ver_dbg1.resize(1, 3);
                ver_dbg1.row(0) = point_middle;
                ver_dbg = vec_list_to_matrix(candidate_pts);
                flag_dbg = false;
            }
        }
        pseudo_geodesic_intersection_filter_by_closeness(curve, angle_degree, pnorm, candidate_pts, id);
        edge_out = candidate_handles[id];
        point_out = candidate_pts[id];
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
        CGMesh::HalfedgeHandle edge_out;
        Eigen::Vector3d point_out, pseudo_vertex_out;
        bool generate_pseudo_vertex;
        found = get_pseudo_vertex_and_trace_forward(cc, curve, target_angle_degree, intersected_handle_tmp, first_point,

                                                    intersected_point_tmp, edge_out, point_out, generate_pseudo_vertex, pseudo_vertex_out);
        if (found)
        {

            update_tracing_list(intersected_point_tmp, pseudo_vertex_out, edge_out, curve, pcurve, handles);

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
    assert(pcurve.size() == curve.size());
    return true;
}
void lsTools::show_pseudo_geodesic_curve(Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, Eigen::MatrixXd &vers)
{

    int vsize = trace_vers.size();
    E0.resize(vsize - 1, 3);
    E1.resize(vsize - 1, 3);
    for (int i = 0; i < vsize - 1; i++)
    {
        E0.row(i) = trace_vers[i];
        E1.row(i) = trace_vers[i + 1];
    }
    vers.resize(vsize, 3);
    for (int i = 0; i < vsize; i++)
    {
        vers.row(i) = trace_vers[i];
    }
}