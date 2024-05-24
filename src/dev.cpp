#include <lsc/tools.h>
#include <lsc/dev.h>
#include <lsc/basic.h>
#include <igl/readOBJ.h>
#include <igl/segment_segment_intersect.h>
#include <igl/predicates/predicates.h>
#include <igl/writeOBJ.h>
inline int orient2d(
    const Eigen::Vector2d& p, const Eigen::Vector2d& q, const Eigen::Vector2d& r)
{
    igl::predicates::exactinit();
    return (int)igl::predicates::orient2d(p, q, r);
}

inline bool point_on_segment2d(const Eigen::Vector2d &p, const Eigen::Vector2d &a, const Eigen::Vector2d &b, bool know_colinear)
{
    if (!know_colinear)
    {
        if (orient2d(p, a, b) != 0)
        {
            return false;
        }
    }
    // colinear

    // on one of the point.
    if (p == a || p == b)
    {
        return true;
    }
    // degeneration
    if (a == b)
    {
        return false;
    }
    // find a point not on the line
    Eigen::Vector2d r = Eigen::Vector2d::Random();
    if (orient2d(r, a, b) == 0)
    {
        while (1)
        {
            r = Eigen::Vector2d::Random();
            if (orient2d(r, a, b) != 0)
            {
                break;
            }
        }
    }
    int o1 = orient2d(p, r, a);
    int o2 = orient2d(p, b, r);
    if (o1 == o2) // p is on the same side of the two edges
    {
        return true;
    }
    return false;
}

inline bool segment_segment_intersection2d(const Eigen::Vector2d &p0, const Eigen::Vector2d &p1,
                                           const Eigen::Vector2d &t0, const Eigen::Vector2d &t1)
{
    if (p0 == t0 || p0 == t1 || p1 == t0 || p1 == t1)
    {
        return true;
    }
    std::array<int, 4> os;
    os[0] = orient2d(p0, t0, t1);
    os[1] = orient2d(p1, t0, t1);
    os[2] = orient2d(t0, p0, p1);
    os[3] = orient2d(t1, p0, p1);

    if (os[0] == 0)
    {
        if (point_on_segment2d(p0, t0, t1, true))
        {
            return true;
        }
    }
    if (os[1] == 0)
    {
        if (point_on_segment2d(p1, t0, t1, true))
        {
            return true;
        }
    }
    if (os[2] == 0)
    {
        if (point_on_segment2d(t0, p0, p1, true))
        {
            return true;
        }
    }
    if (os[3] == 0)
    {
        if (point_on_segment2d(t1, p0, p1, true))
        {
            return true;
        }
    }
    if (os[0] * os[1] * os[2] * os[3] == 0)
    {
        return false;
    }
    // check if the two segments crosses
    if (os[0] * os[1] == -1 && os[2] * os[3] == -1)
    {
        return true;
    }
    return false;
}

// 0: no intersectin, 1: intersect interior. 2. intersect q0 or q1. 3. the point p0 is on the segment
// the code is only for the purpose of our code, not exactly solving the generic problems
// the input p0 p1 defines the ray, p1 is far away enough
int ray_segment_intersection2d(const Eigen::Vector2d &p0, const Eigen::Vector2d &p1,
                               const Eigen::Vector2d &q0, const Eigen::Vector2d &q1)
{
    if (p0 == q0 || p0 == q1)
    {
        return 3;
    }
    int o1 = orient2d(p0, p1, q0);
    int o2 = orient2d(p0, p1, q1);
    if (o1 * o2 > 0) // q0, q2 on the same side of the ray. no intersection
    {
        return 0;
    }
    if(o1 * o2 == 0)
    {
        return 2;
    }
    if(point_on_segment2d(p0,q0,q1,false))
    {
        return 3;
    }
    return segment_segment_intersection2d(p0, p1, q0, q1);
}

// p0, p1 defines a ray (actually it is a segment, p1 is far enough from p0), 
// -1: need to shoot another ray.
// 0: not intersect
// 1: intersect
int ray_source_in_loop(const Eigen::Vector3d& p0in, const Eigen::Vector3d& p1in, const Eigen::MatrixXd &Vl)
{
    Eigen::Vector2d p0(p0in[0], p0in[1]);
    Eigen::Vector2d p1(p1in[0], p1in[1]);
    int ni = 0;
    int vnbr = Vl.rows();
    for (int i = 0; i < vnbr; i++)
    {
        int id0 = i;
        int id1 = (i + 1) % vnbr;
        Eigen::Vector2d q0(Vl(id0, 0), Vl(id0, 1));
        Eigen::Vector2d q1(Vl(id1, 0), Vl(id1, 1));
        int check = ray_segment_intersection2d(p0, p1, q0, q1);
        if(!check)
        {
            continue;
        }
        if(check == 1)
        {
            ni++;
        }
        if(check == 2)
        {
            return -1; // shoot on the vertices of the segment. need another shoot.
        }
        if(check == 3) // the sorce point p0 is on the segment, the point is inside the loop
        {
            return 1;
        }

    }
    return ni % 2;
}
Eigen::Vector2d edge_loop_intersectionPoint(const Eigen::Vector2d &p0, const Eigen::Vector2d &p1,
                                            const Eigen::MatrixXd &Vl, bool& found, int& round)
{
    int vnbr = Vl.rows();
    int eid = -1;
    found = false;
    Eigen::Vector2d result(-1, -1);
    round = -1;
    for (int i = 0; i < vnbr; i++)
    {
        int id0 = i;
        int id1 = (i + 1) % vnbr;
        Eigen::Vector2d q0(Vl(id0, 0), Vl(id0, 1));
        Eigen::Vector2d q1(Vl(id1, 0), Vl(id1, 1));
        bool check = segment_segment_intersection2d(p0, p1, q0, q1);
        if(check)
        {
            if (eid >= 0)
            {
                std::cout << "this edge is intersecting the loop with more than 1 point!\n";
                exit(0);
            }
            eid = i;
            double t0, t1;
            Eigen::Vector3d p03(p0[0], p0[1], 0);
            Eigen::Vector3d q03(q0[0], q0[1], 0);
            Eigen::Vector3d p13(p1[0], p1[1], 0);
            Eigen::Vector3d q13(q1[0], q1[1], 0);
            bool intersect = igl::segment_segment_intersect(p03, p13 - p03, q03, q13 - q03, t0, t1);
            if (!intersect)
            {
                std::cout << "numerical problem happens, we better implement this function by ourselves\n";
                exit(0);
            }
            if (t0 < 1e-6)
            {
                t0 = 0;
                round = 0;
            }
            if (t0 > 1 - 1e-6)
            {
                t0 = 1;
                round = 1;
            }
            Eigen::Vector3d i3 = p03 + t0 * (p13 - p03);
            result[0] = i3[0];
            result[1] = i3[1];
            found = true;

        }
    }
    return result;
}

Eigen::Vector2d findRandomPointOutCircle(const Eigen::Vector3d& bcIn, const double r)
{
    double radius = (Eigen::Vector3d::Random()[0] + 1) / 2 * 3.1415926;
    Eigen::Vector2d result;
    Eigen::Vector2d bc(bcIn[0], bcIn[1]);
    return Eigen::Vector2d(r * (1.001) * cos(radius), r * (1.001) * sin(radius)) + bc;
}

// the vertices are 3d, but we take a projection onto the xy plane
std::vector<int> pointsInLoop(const Eigen::MatrixXd &V, const Eigen::MatrixXd &Vl)
{
    Eigen::Vector3d vmin(Vl.col(0).minCoeff(), Vl.col(1).minCoeff(), Vl.col(2).minCoeff());
    Eigen::Vector3d vmax(Vl.col(0).maxCoeff(), Vl.col(1).maxCoeff(), Vl.col(2).maxCoeff());
    std::cout<<"min and max of the scene: \n"<<vmin.transpose()<<", \n"<<vmax.transpose()<<"\n";
    Eigen::Vector3d centroid = (vmin + vmax) / 2;
    double r = (vmax - vmin).norm() / 2 * (1.001); // radius
    std::vector<int> inLoop;
    for (int i = 0; i < V.rows(); i++)
    {
        Eigen::Vector3d p = V.row(i);
        if (p[0] < vmin[0] || p[1] < vmin[1] || p[0] > vmax[0] || p[1] > vmax[1])
        {
            continue;
        }
        int check = 0;
        while(1)
        {
            Eigen::Vector2d p12d = findRandomPointOutCircle(centroid, r);
            check = ray_source_in_loop(p, Eigen::Vector3d(p12d[0], p12d[1], 0), Vl);
            if (check == -1)
            {
                continue;
            }
            break;
        }
        if(check)
        {
            inLoop.push_back(i);
        }
    }
    return inLoop;
}

// classify the original quads into preserved full quads or quads need to be cutted
// fvin records if the four vertices of a quad are inside or outside.
// flags records the inner vertices as true.
void classifyOriginalQuads(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const std::vector<int> &pIn,
                           std::vector<int> &quadFull, std::vector<int> &quadCut, std::vector<std::array<bool, 4>> &fvin, 
                           std::vector<bool> &flags)
{
    int fnbr = F.rows();
    int vnbr = V.rows();
    flags = std::vector<bool> (vnbr, false);
    for (int i = 0; i < pIn.size(); i++)
    {
        flags[pIn[i]] = true;
    }
    for (int i = 0; i < fnbr; i++)
    {
        int counter = 0;
        std::array<bool,4> inside;
        for (int j = 0; j < 4; j++)
        {
            if (flags[F(i, j)])
            {
                counter++;
                inside[j] = true;
            }
            else
            {
                inside[j] = false;
            }
        }
        if (counter == 4)
        {
            quadFull.push_back(i);
            continue;
        }
        if (counter > 0 && counter < 4)
        {
            quadCut.push_back(i);
            fvin.push_back(inside);
        }
    }
}

void findALlIntersections(const CGMesh &mesh, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                          const Eigen::MatrixXd &Vl, std::vector<CGMesh::EdgeHandle> &ehs,
                          std::vector<Eigen::Vector2d> &intersections)
{
    for (CGMesh::EdgeIter e_it = mesh.edges_begin(); e_it != mesh.edges_end(); ++e_it)
    {   
        OpenMesh::EdgeHandle eh = e_it.handle();
        OpenMesh::HalfedgeHandle hh = mesh.halfedge_handle(eh, 0);
        int vid0 = mesh.from_vertex_handle(hh).idx();
        int vid1 = mesh.to_vertex_handle(hh).idx();
        Eigen::Vector2d p0(V(vid0, 0), V(vid0, 1)), p1(V(vid1, 0), V(vid1, 1));
        bool found = false;
        int round;
        Eigen::Vector2d intersection = edge_loop_intersectionPoint(p0, p1, Vl, found, round);
        if(found)
        {
            ehs.push_back(eh);
            intersections.push_back(intersection);
        }
    }
}

// -1: the edge is not cutted.
// other: the intersection point id (in the intersection list) on the edge
int edgeCutted(const int vf, const int vt, const Eigen::Vector4i &quad, CGMesh &mesh, const std::vector<CGMesh::EdgeHandle> &ehs,
               const std::array<bool, 4> &fvin)
{
    int result = -1;
    bool in0 = false, in1 = false;
    bool f0 = false, f1 = false;
    for (int i = 0; i < 4; i++)
    {
        if (vf == quad[i])
        {
            in0 = fvin[i];
            f0 = true;
        }
        if (vt == quad[i])
        {
            in1 = fvin[i];
            f1 = true;
        }
    }
    if (!f0 || !f1)
    {
        std::cout<<"the edge is not in the quad!!!\n";
        exit(0);
    }
    if (in0 == in1)
    {
        return -1;
    }
    // this edge is cutted, find the corresponding halfedge.
    bool found = false;
    int counter = 0;
    for (auto eh : ehs)
    {
        CGMesh::HalfedgeHandle hh = mesh.halfedge_handle(eh, 0);
        int v0 = mesh.from_vertex_handle(hh).idx();
        int v1 = mesh.to_vertex_handle(hh).idx();
        bool c1 = vf == v0 && vt == v1;
        bool c2 = vf == v1 && vt == v0;
        if (c1 || c2)
        {
            found = true;
            return counter;
        }
        counter++;
    }
    if(!found)
    {
        std::cout<<"ERROR: the cutted edge is not found in the list\n";
        exit(0);
    }
}
// from a cut point to another cut point
std::vector<int> fromCutPtToCutPt(const std::vector<int> &vloop, const int start, const int vnbr)
{
    std::vector<int> result;
    if (vloop[start] < vnbr)
    {
        std::cout << "please take a proper start point!\n";
        exit(0);
    }
    result.push_back(vloop[start]);
    for (int i = 1; i < vloop.size(); i++)
    {
        int id = (i + start) % vloop.size();
        result.push_back(vloop[id]);
        if (vloop[id] >= vnbr)
        {
            break;
        }
    }
    return result;
}

// pick a inner polygon from the two polygons of the cutted quad
// Vflags: size is vnbr, showing which is inside. 
// be: boundary edge
void classifyVloop(const std::vector<bool> &Vflags, const std::vector<int> &vloop, const int vnbr,
                   std::vector<int> &poly, std::array<int,2>& be)
{
    std::array<std::vector<int>, 2> ps;
    int counter = 0;
    for (int i = 0; i < vloop.size(); i++)
    {
        if (vloop[i] >= vnbr) // this is a cut point
        {
            ps[counter] = fromCutPtToCutPt(vloop, i, vnbr);
            be[counter] = vloop[i];
            counter++;
        }
    }
    if(Vflags[ps[0][1]])
    {
        poly = ps[0];
    }
    else{
        poly = ps[1];
        if(!Vflags[ps[1][1]])
        {
            std::cout<<"Both the two polygons are outside!\n";
            exit(0);
        }
    }
    return;
}


// bLoop is the boundary loop
void splitIntersectedQuads(CGMesh &mesh, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                           const Eigen::MatrixXd &Vl, const std::vector<CGMesh::EdgeHandle> &ehs,
                           const std::vector<int> &quadCut, const std::vector<std::array<bool, 4>> &fvin,
                           const std::vector<bool> &Vflags, std::vector<std::array<int, 4>> &rQuads,
                           std::vector<std::array<int, 3>> &rTris, std::vector<std::array<int, 5>> &rPenta,
                           std::vector<int> &bLoop)
{
    int vnbr = V.rows();
    // the generated quads, triangles and pentagons.
    
    std::vector<std::array<int,2>> bEdges; // boundary edges
    for (int i = 0; i < quadCut.size(); i++)
    {
        CGMesh::FaceHandle fh = mesh.face_handle(quadCut[i]);
        std::vector<int> vloop; // the loop for the face
        // iterate over the 4 edges, if there is a cut, add the vertex into the loop.
        for (CGMesh::FaceHalfedgeIter fh_it = mesh.fh_begin(fh); fh_it != mesh.fh_end(fh); ++fh_it)
        {
            OpenMesh::HalfedgeHandle hh = mesh.halfedge_handle(fh_it);
            int vfrom = mesh.from_vertex_handle(hh).idx();
            int vto = mesh.to_vertex_handle(hh).idx();
            vloop.push_back(vfrom);
            int iid = edgeCutted(vfrom, vto, F.row(quadCut[i]), mesh, ehs, fvin[i]);
            if (iid < 0)
            { // not cutted, continue;
                continue;
            }
            else
            {
                vloop.push_back(iid + vnbr);
            }
        }
        if (vloop.size() != 6)
        {
            std::cout << "The vloop has " << vloop.size() << " vertices!\n";
            exit(0);
        }
        // classify the cutted quad into two polygons. keep the inner one.
        std::vector<int> poly;
        std::array<int, 2> be;
        classifyVloop(Vflags, vloop, vnbr, poly, be);
        bEdges.push_back(be);
        int psize = poly.size();
        if(psize<3||psize>5)
        {
            std::cout<<"psize wrong!!!\n";
            exit(0);
        }
        if (psize == 3)
        {
            rTris.push_back({poly[0], poly[1], poly[2]});
        }
        if (psize == 4)
        {
            rQuads.push_back({poly[0], poly[1], poly[2], poly[3]});
        }
        if (psize == 5)
        {
            rPenta.push_back({poly[0], poly[1], poly[2], poly[3], poly[4]});
        }
    }
    // sort the boundary edges into a loop
}

// void findConnectivityForInnerVertices()


void cutBoundaryGenerateTopology()
{
    Eigen::MatrixXd Vquad, Vcurve;
    Eigen::MatrixXi Fquad, Fcurve;

    igl::readOBJ("/Users/wangb0d/Desktop/tmp/Khusrav/CRPC/send/0.obj", Vquad, Fquad);
    igl::readOBJ("/Users/wangb0d/Desktop/tmp/Khusrav/CRPC/send/curve2d.obj", Vcurve, Fcurve);
    int vnbr = Vquad.rows();
    int fnbr = Fquad.rows();
    Eigen::MatrixXd Vproj = Vquad; // prject the vertices onto 2d
    Vproj.col(2) = Eigen::VectorXd::Zero(vnbr);
    std::cout<<"The quad has "<<Vquad.rows()<<" vertices\n";
    // mark the points that are inside the loop.
    std::vector<int> pin = pointsInLoop(Vproj, Vcurve);
    std::cout<<"there are "<<pin.size()<<" vertices in the circle\n";

    Eigen::MatrixXd Vin(pin.size(), 3);
    for (int i = 0; i < pin.size(); i++)
    {
        Vin.row(i) = Vquad.row(pin[i]);
    }
    std::vector<int> quadFull;
    std::vector<int> quadCut;
    std::vector<std::array<bool, 4>> fvin;
    std::vector<bool> vFlags;
    // record the quads that are cutted.
    classifyOriginalQuads(Vquad, Fquad, pin, quadFull, quadCut, fvin, vFlags);
    std::cout<<"there are "<<quadFull.size()<<" quads are full, and "<<quadCut.size()<<" are cutted\n";
    MeshProcessing mp;
    CGMesh mesh;
    mp.matrix2Mesh(mesh, Vquad, Fquad);
    std::vector<CGMesh::EdgeHandle> ehs;
    std::vector<Eigen::Vector2d> intersections;
    // find all the intersections on the edges.
    findALlIntersections(mesh, Vquad, Fquad, Vcurve, ehs, intersections);
    std::cout << "find intersection point nbr, " << ehs.size() << "\n";
    // get the cutted polygons.
    // the triangles, the quads and the pentagons
    std::vector<std::array<int,3>> rT;
    std::vector<std::array<int,4>> rQ;
    std::vector<std::array<int,5>> rP;
    splitIntersectedQuads(mesh, Vquad, Fquad, Vcurve, ehs, quadCut, fvin, vFlags, rQ, rT, rP);

    // get the connection for all the inner vertices.

    // write the points inside
    igl::writeOBJ("/Users/wangb0d/Desktop/tmp/Khusrav/CRPC/send/pIn.obj", Vin, Eigen::MatrixXi());
}
