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
        if ((p1 - p0).dot(p_end - p1) < 0) // to avoid trace back
        {
            return false;
        }
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
            // Eigen::Vector3d direction0 = (p1 - p0).normalized();
            // Eigen::Vector3d direction1 = (vs + t * ves - p1).normalized();
            // Eigen::Vector3d pseudo_norm = direction0.cross(direction1).normalized();
            // double dot_product = pseudo_norm.dot(pnorm);
            // double cos_sqr = dot_product * dot_product;
            // double real_angle = acos(dot_product) * 180 / LSC_PI;
        }
    }

    // debug end
    return p_end.size();
}
// not geodesic
void find_intersection_on_halfedge(const CGMesh &lsmesh, const Eigen::MatrixXd &V, CGMesh::HalfedgeHandle &checking_he,
                                   const Eigen::Vector3d &p0, const Eigen::Vector3d &p1, const Eigen::Vector3d &pnorm, const double angle,
                                   std::vector<CGMesh::HalfedgeHandle> &edge_out, std::vector<Eigen::Vector3d> &p_end)
{
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
    find_intersection_on_halfedge(lsmesh, V, checking_he, p0, p1, pnorm, angle, edge_out, p_end);
    // check the previous halfedge
    checking_he = lsmesh.prev_halfedge_handle(ophe);
    assert(lsmesh.face_handle(edge_middle).idx() != lsmesh.face_handle(ophe).idx());
    assert(lsmesh.face_handle(ophe).idx() == lsmesh.face_handle(checking_he).idx());
    find_intersection_on_halfedge(lsmesh, V, checking_he, p0, p1, pnorm, angle, edge_out, p_end);
    // the next parts are for the cases the edge_middle is tangent to the curve
    checking_he = lsmesh.next_halfedge_handle(edge_middle);
    find_intersection_on_halfedge(lsmesh, V, checking_he, p0, p1, pnorm, angle, edge_out, p_end);
    checking_he = lsmesh.prev_halfedge_handle(edge_middle);
    find_intersection_on_halfedge(lsmesh, V, checking_he, p0, p1, pnorm, angle, edge_out, p_end);

    /*std::cout << "lower level check" << std::endl;
    for (int i = 0; i < p_end.size(); i++)
    {
        Eigen::Vector3d direc0 = (p1 - p0).normalized();
        Eigen::Vector3d direc1 = (p_end[i] - p1).normalized();
        Eigen::Vector3d dddnorm = direc0.cross(direc1).normalized();
        double cosangle = pnorm.dot(dddnorm);
        std::cout << "lower, cos^2 " << cosangle * cosangle << std::endl;
    }*/
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

// angle_degree1 can be larger than 180 or smaller than 0.
bool angles_match(const double angle_degree1, const double angle_degree2)
{
    double uniformed = angle_degree1;
    if (angle_degree1 > 180)// we allow inverse directions
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

bool binormal_correct_angle(const Eigen::Vector3d &seg_in, const Eigen::Vector3d &seg_out, const Eigen::Vector3d &pnorm, const double angle_degree)
{
    Eigen::Vector3d dirin = seg_in.normalized();
    Eigen::Vector3d dirout = seg_out.normalized();
    if (dirin.dot(dirout) < 0)
    { // avoid sharp turn.
        return false;
    }
    Eigen::Vector3d norm = dirin.cross(dirout).normalized();
    double dot_product1 = norm.dot(pnorm);
    double real_angle = acos(dot_product1) * 180 / LSC_PI;// 0 ~ 180
    if (isnan(real_angle))
    {
        return false;
    }
    // if degree is 30, then 30 or 150 or 210 or 330 will be allowed.
    // but still need to filer out 150 and 330
    bool match = angles_match(angle_degree, real_angle) || angles_match(180 - angle_degree, real_angle);
    if (!match)
    {
        return false;
    }

    // if we are checking asymptotic, then don't need further check
    if (angles_match(angle_degree, 0) || angles_match(angle_degree, 180))
    {
        return true;
    }
    // tangent is othogonal to surface normal and the curve bi-normal
    Eigen::Vector3d tangent = norm.cross(pnorm).normalized();
    // if the target is the same as real, then target is in II and III, real is in I and IV. cross should be same direction as tangent
    if (angles_match(angle_degree, real_angle))
    {
        if (tangent.dot(seg_in) > 0)
        {
            return true;
        }
        return false;
    }
    // if the target is -real, then target is in II and III, and real is in IV and I. cross should be inverse direction as tangent
    if (angles_match(180 - angle_degree, real_angle))
    {
        if (tangent.dot(seg_in) < 0)
        {
            return true;
        }
        return false;
    }

    return false;
}


void pseudo_geodesic_intersection_filter_by_closeness(
    const std::vector<Eigen::Vector3d> &curve, const double angle_degree, const Eigen::Vector3d &pnorm,
    const std::vector<Eigen::Vector3d> &candi_points, int &id)
{
    assert(curve.size() >= 2); // take previous points to find the next one
    // assert(candi_points.size() > 0);
    if (candi_points.size() == 0)
    {
        id = -1;
    }
    // if (candi_points.size() == 1)
    // {
    //     id = 0;
    //     return;
    // }
    int size = curve.size();
    Eigen::Vector3d dire1 = (curve[size - 1] - curve[size - 2]).normalized();
    int closest_id = -1;

    // int closest_id_consider_angle=-1;
    double closest_distance = std::numeric_limits<double>::max();

    for (int i = 0; i < candi_points.size(); i++)
    {

        Eigen::Vector3d dire2 = (candi_points[i] - curve[size - 1]).normalized();

        if (!binormal_correct_angle(curve[size - 1] - curve[size - 2], candi_points[i] - curve[size - 1], pnorm, angle_degree))
        {
            continue;
        }
        if (dire1.dot(dire2) < 0)
        {
            continue;
        }
        double distance = (candi_points[i] - curve[size - 1]).norm();
        if (distance < closest_distance) // select the cloest one
        {
            closest_distance = distance;
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
    const Eigen::Vector3d &normal,
    const std::vector<Eigen::Vector3d> &candi_points, int &id)
{
    //std::cout<<"before filtering size "<<candi_points.size()<<std::endl;
    int closest_id = -1;
    double closest_distance = std::numeric_limits<double>::max();
    for(int i=0;i<candi_points.size(); i++){
        Eigen::Vector3d direction=candi_points[i]-start_point;
        direction=direction.normalized();
        double dot_product=direction.dot(reference_direction.normalized());
        double real_angle = acos(dot_product) * 180 / LSC_PI;
        if(!angles_match(real_angle,angle_degree)){
            //std::cout<<"angle not match, id "<<i<<" angle "<<real_angle <<std::endl;
            continue;
        }
        if (isnan(real_angle))
        {
            std::cout<<"angle nan, id "<<i<<std::endl;
            continue;
        }
        Eigen::Vector3d tangent = reference_direction.cross(direction).normalized();
        if (angle_degree >= 0 && angle_degree <= 180)
        {
            if (tangent.dot(normal) >= 0)
            {
            }
            else{
                //std::cout<<"angle cross wrong, id "<<i<<std::endl;
                continue;
            }
        }
        else
        { // 180~360 degree
            if (tangent.dot(normal) <= 0)
            {
            }
            else{
                //std::cout<<"angle cross wrong, id "<<i<<std::endl;
                continue;
            }
        }
        double dis=(candi_points[i]-start_point).norm();
        if(dis<closest_distance){
            closest_distance=dis;
            closest_id=i;
        }
    }

    // already found the vertex that forms the angle most close to the given one
    id = closest_id;
    return;
}
void get_one_ring_vertices(CGMesh &lsmesh, const int id, std::vector<int> &pts)
{
    pts.clear();
    CGMesh::VertexHandle vh = lsmesh.vertex_handle(id);
    for (CGMesh::VertexVertexIter vvi = lsmesh.vv_begin(vh); vvi != lsmesh.vv_end(vh); ++vvi)
    {
        CGMesh::VertexHandle ver = vvi.handle();
        assert(ver.idx()>=0);
        pts.push_back(ver.idx());
        //std::cout << ver.idx() << ", ";
    }
    //std::cout << std::endl;
}
bool is_boundary_edge(CGMesh &lsmesh, const CGMesh::HalfedgeHandle &he)
{
    CGMesh::HalfedgeHandle ophe=lsmesh.opposite_halfedge_handle(he);
    bool result=(lsmesh.face_handle(he).idx()<0)||(lsmesh.face_handle(ophe).idx()<0);
    return result;
}

// in the first iteration, start_point_ids is empty so we can search the adjecent trianlge or the one-ring
// edges (depends on if it is a vertex or not)
// start_point_ids may contain duplicated points
bool lsTools::get_checking_edges(const std::vector<int> &start_point_ids, const CGMesh::HalfedgeHandle &edge_middle, const Eigen::Vector3d &point_middle,
                                 Efunc &edges_checked, Efunc &points_checked, NeighbourInfo &ninfo, std::vector<int> &point_to_check)
{
    //std::cout << "--------------getting checking edges, itr " << ninfo.round << std::endl;
    ninfo.edges.clear();
    point_to_check.clear();
    int eid=lsmesh.edge_handle(edge_middle).idx();
    if (ninfo.round > 2)
    { // search only two rings
        return false;
    }

    if (start_point_ids.size() == 0)
    { // we initialize the searching
        ninfo.round += 1;
        ninfo.is_vertex = false;
        int ver_from_id = lsmesh.from_vertex_handle(edge_middle).idx();
        int ver_to_id = lsmesh.to_vertex_handle(edge_middle).idx();
        Eigen::Vector3d ver_from = V.row(ver_from_id);
        Eigen::Vector3d ver_to = V.row(ver_to_id);

        double dist_from = (point_middle - ver_from).norm();
        double dist_to = (point_middle - ver_to).norm();
        double dist_total = dist_from + dist_to;
        double from_ratio = dist_from / dist_total;
        double to_ratio = dist_to / dist_total;
        if (from_ratio <= MERGE_VERTEX_RATIO)
        {

            ninfo.is_vertex = true;
            ninfo.center_handle = lsmesh.from_vertex_handle(edge_middle);
        }
        if (to_ratio <= MERGE_VERTEX_RATIO)
        {

            ninfo.is_vertex = true;
            ninfo.center_handle = lsmesh.to_vertex_handle(edge_middle);
        }
        if (ninfo.is_vertex) // if it is a vertex, we search for one-ring opposite edges.
        {
            //std::cout << "fall in vertex point " << std::endl;
            points_checked.coeffRef(ninfo.center_handle.idx()) = 1;
            ninfo.pnorm = norm_v.row(ninfo.center_handle.idx());
            for (CGMesh::VertexOHalfedgeIter voh_itr = lsmesh.voh_begin(ninfo.center_handle);
                 voh_itr != lsmesh.voh_end(ninfo.center_handle); ++voh_itr)
            {
                auto heh = voh_itr.handle();
                CGMesh::HalfedgeHandle edge_to_check = lsmesh.next_halfedge_handle(heh);

                point_to_check.push_back(lsmesh.to_vertex_handle(edge_to_check).idx());
                assert(lsmesh.from_vertex_handle(edge_to_check).idx() != ninfo.center_handle.idx());
                assert(lsmesh.to_vertex_handle(edge_to_check).idx() != ninfo.center_handle.idx());
                ninfo.edges.push_back(edge_to_check);
                /*std::cout << "edge, \n"
                          << V.row(lsmesh.from_vertex_handle(edge_to_check).idx()) << "\n"
                          << V.row(lsmesh.to_vertex_handle(edge_to_check).idx()) << std::endl;*/
            }
        }
        else // it is not a vertex and we search on edges of the opposite triangle
        {
            /*std::cout << "fall in middle point " << std::endl;
            std::cout << "from id " << ver_from_id << std::endl;
            std::cout << "to id " << ver_to_id << std::endl;*/
            int fid1 = lsmesh.face_handle(edge_middle).idx();
            int fid2 = lsmesh.opposite_face_handle(edge_middle).idx();
            assert(fid1 != fid2);
            ninfo.pnorm = norm_e.row(eid); // the normal of the pseudo-vertex

            CGMesh::HalfedgeHandle ophe = lsmesh.opposite_halfedge_handle(edge_middle);
            CGMesh::HalfedgeHandle ophe_next = lsmesh.next_halfedge_handle(ophe);
            ninfo.edges.push_back(lsmesh.next_halfedge_handle(ophe));
            ninfo.edges.push_back(lsmesh.prev_halfedge_handle(ophe));
            std::vector<int> temp_v;
            get_one_ring_vertices(lsmesh, ver_from_id, temp_v);
            for (int vid : temp_v)
            {
                point_to_check.push_back(vid);
            }
            get_one_ring_vertices(lsmesh, ver_to_id, temp_v);
            for (int vid : temp_v)
            {
                point_to_check.push_back(vid);
            }
            if (!is_boundary_edge(lsmesh, edge_middle))
            { // if the opposite vertex exists, we push it into the list
                int ver_oppo_id = lsmesh.to_vertex_handle(ophe_next).idx();
                assert(ver_oppo_id != ver_from_id && ver_oppo_id != ver_to_id);
                get_one_ring_vertices(lsmesh, ver_oppo_id, temp_v);
                for (int vid : temp_v)
                {
                    point_to_check.push_back(vid);
                }
            }

            // ninfo.edges.push_back(lsmesh.next_halfedge_handle(edge_middle));
            // ninfo.edges.push_back(lsmesh.prev_halfedge_handle(edge_middle));
        }
        // updated checked lists
        for (int i = 0; i < ninfo.edges.size(); i++)
        {
            int id = lsmesh.edge_handle(ninfo.edges[i]).idx();
            edges_checked.coeffRef(id) = 1;
        }
    }
    else
    { // it is not the first iteration
        ninfo.round++;
        for (int id : start_point_ids)
        {
            CGMesh::VertexHandle vh = lsmesh.vertex_handle(id);
            if (points_checked.coeffRef(id) == 1)
            { // if this point is already checked, skip. otherwise, mark it as checked.
                continue;
            }
            points_checked.coeffRef(id) = 1;
            std::vector<CGMesh::HalfedgeHandle> tmp_hds;
            for (CGMesh::VertexOHalfedgeIter voh_itr = lsmesh.voh_begin(vh);
                 voh_itr != lsmesh.voh_end(vh); ++voh_itr) // for each edge shooting out from this edge
            {

                CGMesh::HalfedgeHandle heh = voh_itr.handle();                      // the neibouring edges
                CGMesh::HalfedgeHandle heh1 = lsmesh.next_halfedge_handle(voh_itr); // the opposite
                CGMesh::VertexHandle ver = lsmesh.to_vertex_handle(heh);            // the one-ring vertices
                if (points_checked.coeffRef(ver.idx()) != 1)                        // if this point is not checked yet, we add this to be checked
                {
                    point_to_check.push_back(ver.idx());
                }
                int neighbour_eid = lsmesh.edge_handle(heh).idx();
                if (edges_checked.coeffRef(neighbour_eid) != 1)
                { // if edge is not checked yet, add into list, and mark it as checked
                    ninfo.edges.push_back(heh);
                    edges_checked.coeffRef(neighbour_eid) = 1;
                }
                int opposite_eid = lsmesh.edge_handle(heh1).idx();
                if (edges_checked.coeffRef(opposite_eid) != 1)
                {
                    ninfo.edges.push_back(heh);
                    edges_checked.coeffRef(opposite_eid) = 1;
                }
            }
        }
    }
    //std::cout << "--------------round " << ninfo.round << std::endl;
    return point_to_check.size();
    return true;
}
bool lsTools::trace_pseudo_geodesic_forward(
    const NeighbourInfo &ninfo,
    const std::vector<Eigen::Vector3d> &curve, const double angle_degree,
    const CGMesh::HalfedgeHandle &edge_middle,
    const Eigen::Vector3d &point_in, const Eigen::Vector3d &point_middle,
    CGMesh::HalfedgeHandle &edge_out,
    Eigen::Vector3d &point_out)
{
    //std::cout << "inside trace forward, curve size = " << curve.size() << std::endl;

    bool is_geodesic = false;
    // int coplanar=edge_is_coplanar(lsmesh, lsmesh.edge_handle(edge_middle),  V, F);
    // if(coplanar==1){// TODO calculate coplanarity outside
    //     std::cout<<"COPLANARITY HAPPENED"<<std::endl;
    //     is_geodesic=true;
    // }

    double angle_radian = angle_degree * LSC_PI / 180.; // the angle in radian
    // std::vector<CGMesh::halfedge_handle>
    std::vector<Eigen::Vector3d> candidate_pts;
    std::vector<CGMesh::HalfedgeHandle> candidate_handles;

    // start to calculate the pseudo-vertex

    if (is_boundary_edge(lsmesh, edge_middle))
    {
        // this is a boundary edge on the boundary loop.
        return false;
    }

    if (angle_degree > 90 - ANGLE_TOLERANCE && angle_degree < 90 + ANGLE_TOLERANCE)
    { // it means it is a geodesic
        is_geodesic = true;
        // std::cout << "**tracing a geodesic" << std::endl;
    }
    for (CGMesh::HalfedgeHandle edge_to_check : ninfo.edges)
    {
        Eigen::Vector3d vs = V.row(lsmesh.from_vertex_handle(edge_to_check).idx());
        Eigen::Vector3d ve = V.row(lsmesh.to_vertex_handle(edge_to_check).idx());

        if (is_geodesic)
        { // it means it is a geodesic
            bool found;

            found = find_geodesic_intersection_p1_is_ver(point_in, point_middle, vs, ve, ninfo.pnorm, point_out);
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
            bool found;

            found = find_osculating_plane_intersection_not_geodesic_p1_is_ver(point_in, point_middle, vs, ve, ninfo.pnorm, angle_radian, tmp_candidates);

            if (found)
            {
                extend_vector(candidate_pts, tmp_candidates);
                vector_duplicate(edge_to_check, tmp_candidates.size(), tmp_hehs);
                extend_vector(candidate_handles, tmp_hehs);
            }
        }
    }
    if (is_geodesic)
    {
        return false;
    }
    else // filter the points
    {
        int id;
        //std::cout << "before filtering, size " << candidate_pts.size() << std::endl;
        pseudo_geodesic_intersection_filter_by_closeness(curve, angle_degree, ninfo.pnorm, candidate_pts, id);
        // // DEBUG
        // if(curve.size()==15){
        //     id=0;
        // }
        if (id < 0)
        {
            //std::cout << "after filtering, no results" << std::endl;
            return false;
        }

        edge_out = candidate_handles[id];
        point_out = candidate_pts[id];
        return true;
    }
    ////////////////////////

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

// this code extend the candidate point list.
// after calling this code, filter the points according to the quadrant and closeness.
bool find_initial_direction_intersection_on_edge(const Eigen::Vector3d &start_point, const double angle_degree,
                                                 const Eigen::Vector3d &reference_direction, const Eigen::Vector3d &vs, const Eigen::Vector3d &ve,
                                                 const Eigen::Vector3d &normal, const CGMesh::HalfedgeHandle &handle_in, std::vector<CGMesh::HalfedgeHandle> &handle_out,
                                                 std::vector<Eigen::Vector3d> &point_out)
{
    //////////////////
    // Eigen::Vector3d d1=vs-start_point;
    // Eigen::Vector3d d2=ve-start_point;
    // d1=d1.normalized();
    // d2=d2.normalized();

    std::array<double, 2> roots;
    std::vector<double> equation(3);
    // std::cout<<"vs angle "<<
    double angle_radian = angle_degree * LSC_PI / 180.;
    Eigen::Vector3d ves = ve - vs;
    Eigen::Vector3d vssp = vs - start_point;
    double degree1 = ves.dot(reference_direction);
    double degree0 = vssp.dot(reference_direction);
    if (!angles_match(angle_degree, 90))
    {
        double left_square_0 = degree0 * degree0;
        double left_square_1 = 2 * degree0 * degree1;
        double left_square_2 = degree1 * degree1;
        double right_square_0 = vssp.dot(vssp);
        double right_square_1 = 2 * vssp.dot(ves);
        double right_square_2 = ves.dot(ves);
        double right_cons = cos(angle_radian) * cos(angle_radian) * reference_direction.dot(reference_direction);
        // std::cout<<"## right cos "<<right_cons<<std::endl;
        right_square_0 *= right_cons;
        right_square_1 *= right_cons;
        right_square_2 *= right_cons;
        double eq0 = left_square_0 - right_square_0;
        double eq1 = left_square_1 - right_square_1;
        double eq2 = left_square_2 - right_square_2;
        
        equation[0] = eq0;
        equation[1] = eq1;
        equation[2] = eq2;
        
        /*std::cout << "equation " << eq0 << " " << eq1 << " " << eq2 << std::endl;
        std::cout << "b^2-4ac= " << eq1 * eq1 - 4 * eq0 * eq2 << std::endl;*/
        equation = polynomial_simplify(equation);

        bool found = quadratic_solver(equation, roots);
        if (!found)
        {
            return false;
        }
    }
    else
    {
        if (degree1 == 0)
        {
            return false;
        }
        double value = -degree0 / degree1;
        roots[0] = value;
        roots[1] = -1;
    }

    std::vector<Eigen::Vector3d> p_end;
    for (int i = 0; i < 2; i++)
    {
        double t = roots[i];
        if (!angles_match(angle_degree, 90))
        {
            //std::cout << "t " << t << " equation value " << polynomial_value(equation, t) << std::endl;
        }

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
        Eigen::Vector3d direction = (p_end[i] - start_point).normalized();
        double dot_product = direction.dot(reference_direction.normalized());
        double real_angle = acos(dot_product) * 180 / LSC_PI;
        //std::cout << "recheck angle " << real_angle << std::endl;
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
    assert((vt-vf).norm()>1e-8);
    Eigen::Vector3d reference_direction = (vt - vf).normalized();
    std::vector<CGMesh::HalfedgeHandle> handle_out;
    std::vector<Eigen::Vector3d> point_out;

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
        Eigen::Vector3d normal = norm_v.row(center_handle.idx());
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
    assert(lsmesh.face_handle(start_boundary_edge).idx()>=0);
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
    // pseudo_geodesic_intersection_satisfy_angle_quadrant(start_boundary_angle_degree, reference_direction, start_point,
    //                                                     normal, handle_out, point_out, handle_out2, point_out2);
    if (point_out.size() == 0)
    {
        std::cout << "ERROR: no point satisfies the angle quadrant in initialization" << std::endl;
        return false;
    }
    int result_id = -1;
    initial_segment_intersection_filter_by_closeness(start_boundary_angle_degree, start_point,
                                                     reference_direction,normal, point_out, result_id);
    if(result_id==-1){
        std::cout<<"in init_pseudo_geodesic_first_segment() not finding any point"<<std::endl;
        return false;
    }
    intersected_handle = handle_out[result_id];
    intersected_point = point_out[result_id];
    return true;
}
void update_tracing_list(const Eigen::Vector3d &ver, const Eigen::Vector3d &pver, const CGMesh::HalfedgeHandle &handle,
                         std::vector<Eigen::Vector3d> &verlist, std::vector<Eigen::Vector3d> &pverlist, std::vector<CGMesh::HalfedgeHandle> &handle_list)
{
    verlist.push_back(ver);
    pverlist.push_back(pver);
    handle_list.push_back(handle);
}


bool lsTools::trace_single_pseudo_geodesic_curve(const double target_angle_degree,
                                                 const CGMesh::HalfedgeHandle &start_boundary_edge, const double &start_point_para,
                                                 const double start_boundary_angle_degree,
                                                 std::vector<Eigen::Vector3d> &curve,
                                                 std::vector<CGMesh::HalfedgeHandle> &handles)
{
    curve.clear();
    handles.clear();
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

    Eigen::Vector3d first_point = get_3d_ver_from_t(start_point_para, V.row(lsmesh.from_vertex_handle(start_boundary_edge).idx()),
        V.row(lsmesh.to_vertex_handle(start_boundary_edge).idx()));

    curve.push_back(first_point);
    handles.push_back(start_boundary_edge);
    curve.push_back(intersected_point_tmp);
    handles.push_back(intersected_handle_tmp);
    //std::cout << "target angle " << target_angle_degree << std::endl;
    for (int i = 0;; i++)
    {
        //std::cout << "XXXXXXXXXXXXXXXXXXXXXXXXXXXXX" << std::endl;
        CGMesh::HalfedgeHandle edge_out;
        Eigen::Vector3d point_out;
        NeighbourInfo ninfo;
        std::vector<int> start_point_ids; // initially this is empty
        Efunc edges_checked;
        edges_checked.resize(E.rows());
        Efunc points_checked;
        points_checked.resize(V.rows());
        bool edges_found = true;
        for (int j = 0;; j++)
        {
            //
            std::vector<int> point_to_check;
            edges_found = get_checking_edges(start_point_ids, intersected_handle_tmp, intersected_point_tmp,
                                             edges_checked, points_checked, ninfo, point_to_check);
            //std::cout << "get checked edges, edge size " << ninfo.edges.size() << std::endl;
            if (!edges_found)
            {
                //std::cout<<"edges not found"<<std::endl;
                break;
            }
            found = trace_pseudo_geodesic_forward(ninfo, curve, target_angle_degree, intersected_handle_tmp, // edge middle
                                                  first_point,                                               // point in
                                                  intersected_point_tmp,                                     // point middle
                                                  edge_out, point_out);
            if (found)
            {
                curve.push_back(point_out);
                handles.push_back(edge_out);
                intersected_handle_tmp = edge_out;   // middle edge = edge out
                first_point = intersected_point_tmp; // first point = point middle
                intersected_point_tmp = point_out;   // point middle = point out. but may eventually be converted to a pseudo point
                break;
            }
            else
            {
                start_point_ids = point_to_check;
            }
        }
        if (!edges_found)
        {
            break;
        }
    }
    return true;
}
void lsTools::show_pseudo_geodesic_curve(std::vector<Eigen::MatrixXd> &E0, std::vector<Eigen::MatrixXd> &E1, Eigen::MatrixXd &vers)
{
    //std::cout<<"check 1"<<std::endl;
    // for debug
    

    // for debug
    // if(trace_vers.size()==0){
    //     std::cout<<"Please Trace The Curves First"<<std::endl;
    // }
    E0.resize(trace_vers.size());
    E1.resize(trace_vers.size());
    int pnbr = 0;
    for (int i = 0; i < trace_vers.size(); i++)
    {
        int vsize = trace_vers[i].size();
        assert(vsize>=1);
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
    vers.resize(pnbr, 3);
    pnbr = 0;
    for (int i = 0; i < trace_vers.size(); i++)
    {
        for (int j = 0; j < trace_vers[i].size(); j++)
        {
            vers.row(pnbr) = trace_vers[i][j];
            pnbr++;
        }
    }
}