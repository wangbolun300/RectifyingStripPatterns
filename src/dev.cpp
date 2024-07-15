#include <lsc/tools.h>
#include <lsc/dev.h>
#include <lsc/basic.h>
#include <igl/readOBJ.h>
#include <igl/segment_segment_intersect.h>
#include <igl/predicates/predicates.h>
#include <igl/writeOBJ.h>
#include <fstream>
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
Eigen::Vector3d edge_loop_intersectionPoint(const Eigen::Vector2d &p0, const Eigen::Vector2d &p1,
                                            const Eigen::Vector3d &pp0, const Eigen::Vector3d &pp1,
                                            const Eigen::MatrixXd &Vl, bool &found, int &round)
{
    int vnbr = Vl.rows();
    int eid = -1;
    found = false;
    Eigen::Vector3d result(-1, -1, -1);
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
            Eigen::Vector3d i3 = pp0 + t0 * (pp1 - pp0);
            result = i3;
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
    // std::cout<<"cutted fids\n";
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
            // std::cout<<i<<", ";
            // if(i==304)
            // {
            //     std::cout << "checking " << i << ", the counter, " << counter << " verin, " << inside[0] << ", " << inside[1] << ", " << inside[2]
            //               << ", " << inside[3] << "\n";
            //     std::cout<<"the face, "<< F.row(i)<<"\n";      
            // }
        }
        
    }
    std::cout<<"\n";
}

void findALlIntersections(const CGMesh &mesh, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                          const Eigen::MatrixXd &Vl, std::vector<CGMesh::EdgeHandle> &ehs,
                          std::vector<Eigen::Vector3d> &intersections)
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
        Eigen::Vector3d intersection = edge_loop_intersectionPoint(p0, p1, V.row(vid0), V.row(vid1), Vl, found, round);
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
        std::cout<<"the edge: "<<vf<<", "<<vt<<", the quad, "<<quad.transpose()<<"\n";
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
    bool found = false;
    for (int i = 1; i < vloop.size(); i++)
    {
        int id = (i + start) % vloop.size();
        result.push_back(vloop[id]);
        if (vloop[id] >= vnbr)
        {
            found = true;
            break;
        }
    }
    if(!found)
    {
        std::cout<<"the loop end point is not found!\n";
        exit(0);
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
            std::cout<<"the vnbr "<<vnbr<<"\n";
            for(auto r : ps[0])
            {
                std::cout<<r<<", ("<<Vflags[r]<<"), ";
            }
            std::cout<<"\n";
            for(auto r : ps[1])
            {
                std::cout<<r<<", ("<<Vflags[r]<<"), ";
            }
            std::cout<<"\n";
            exit(0);
        }
    }
    return;
}

void sortBoundaryEdgesIntoLoop(const std::vector<std::array<int, 2>> &bes, std::vector<int> &loop)
{
    int enbr = bes.size();
    std::vector<bool> echecked(enbr, false); 
    loop.push_back(bes[0][0]);
    loop.push_back(bes[0][1]);
    echecked[0] = true;
    int nbrChecked = 1;
    while(1)
    {
        if (nbrChecked == enbr - 1)
        {
            break;
        }
        int id = loop.back();
        for (int i = 0; i < enbr; i++)
        {
            if(echecked[i])
            {
                continue;
            }
            if (id == bes[i][0] || id == bes[i][1])
            {
                if (id == bes[i][0])
                {
                    loop.push_back(bes[i][1]);
                }
                else{
                    loop.push_back(bes[i][0]);
                }
                echecked[i] = true;
                nbrChecked++;
            }
        }
    }


}

bool verInQuad(const int vid, const Eigen::Vector4i& face)
{
    if(vid==face[0]||vid==face[1]||vid==face[2]||vid==face[3])
    {
        return true;
    }
    return false;
}
// bLoop is the boundary loop
// split the quads into 
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
        std::vector<int> vloop;                // the loop for the face
        // iterate over the 4 edges, if there is a cut, add the vertex into the loop.
        CGMesh::FaceHalfedgeIter fh_it = mesh.fh_begin(fh);
        OpenMesh::HalfedgeHandle hh = fh_it.handle();

        for (int j = 0; j < 4; j++)
        {
            hh = mesh.next_halfedge_handle(hh);
            int vfrom = mesh.from_vertex_handle(hh).idx();
            int vto = mesh.to_vertex_handle(hh).idx();
            for(int c : vloop)
            {
                if(c==vfrom)
                {
                    std::cout<<"the vertex is already in \nwhich edge "<<j<<"\n";
                    exit(0);
                }
            }

            vloop.push_back(vfrom);
            bool notin = false;
            if(!verInQuad(vfrom, F.row(quadCut[i])))
            {
                std::cout<<"ver is not in!\n"<<vfrom<<"\n";
                notin = true;
            }
            if(!verInQuad(vto, F.row(quadCut[i])))
            {
                std::cout<<"ver is not in!\n"<<vfrom<<"\n";
                notin = true;
            }
            if(notin)
            {
                std::cout<<"face , "<<F.row(quadCut[i])<<"\nwhich edge, "<<j<<"\n";
                std::cout<<"the face vertices itr:\n";
                for(CGMesh::FaceVertexIter fvitr = mesh.fv_begin(fh);fvitr != mesh.fv_end(fh); ++fvitr)
                {
                    std::cout<<fvitr.handle().idx()<<", ";
                }
                std::cout<<"\n";
            }
            
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

        // for (CGMesh::FaceEdgeIter fh_it = mesh.fe_begin(fh); fh_it != mesh.fe_end(fh); ++fh_it)
        // {
        //     OpenMesh::HalfedgeHandle hh = mesh.halfedge_handle(fh_it, 0);
        //     int vfrom = mesh.from_vertex_handle(hh).idx();
        //     int vto = mesh.to_vertex_handle(hh).idx();
        //     for(int c : vloop)
        //     {
        //         if(c==vfrom)
        //         {
        //             std::cout<<"the vertex is already in \n";
        //             exit(0);
        //         }
        //     }

        //     vloop.push_back(vfrom);

        //     if(!verInQuad(vfrom, F.row(quadCut[i])))
        //     {
        //         std::cout<<"ver is not in!\n";
        //         exit(0);
        //     }
        //     if(!verInQuad(vto, F.row(quadCut[i])))
        //     {
        //         std::cout<<"ver is not in!\n";
        //         exit(0);
        //     }
            
        //     int iid = edgeCutted(vfrom, vto, F.row(quadCut[i]), mesh, ehs, fvin[i]);
        //     if (iid < 0)
        //     { // not cutted, continue;
        //         continue;
        //     }
        //     else
        //     {
        //         vloop.push_back(iid + vnbr);
        //     }
        // }
        if (vloop.size() != 6)
        {
            std::cout << "The vloop has " << vloop.size() << " vertices!\n";
            exit(0);
        }
        // classify the cutted quad into two polygons. keep the inner one.
        std::vector<int> poly;
        std::array<int, 2> be;
        bool foundin = false;
        for(int vid : vloop)
        {
            if (vid < vnbr)
            {
                if(Vflags[vid])
                {
                    foundin = true;
                }
            }
        }
        if (!foundin)
        {
            std::cout << "error: the cutted quad doesn't have any inside ver, fid, " << quadCut[i] << "\n";
            std::cout << "face, " << F.row(quadCut[i]) << "\n";
            std::cout << "loop: \n";
            for (int vid : vloop)
            {
                std::cout << vid << ", ";
            }
            std::cout << "\n";

            // exit(0);
        }
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
    sortBoundaryEdgesIntoLoop(bEdges, bLoop);
}

// returns the halfedge id of the edge v0-v1
int whichHalfedge(const int v0, const int v1, CGMesh &mesh, const std::vector<CGMesh::EdgeHandle> &ehs)
{

    for (int i = 0; i < ehs.size(); i++)
    {
        CGMesh::HalfedgeHandle hh = mesh.halfedge_handle(ehs[i], 0);
        int vf = mesh.from_vertex_handle(hh).idx();
        int vt = mesh.to_vertex_handle(hh).idx();
        if ((v0 == vf && v1 == vt) || (v0 == vt && v1 == vf))
        {
            return i;
        }
    }
    std::cout << "We cannot find the halfedge!\n";
    exit(0);
}

void findConnectivityForInnerVertices(const std::vector<bool> &Vflags, const std::vector<int> &pin, const int vnbr,
                                      CGMesh &mesh, const std::vector<CGMesh::EdgeHandle> &ehs,
                                      std::vector<std::array<int,5>>& connect)
{
    for (int i = 0; i < pin.size(); i++)
    {
        CGMesh::VertexHandle vh = mesh.vertex_handle(pin[i]);
        std::vector<int> nbs; // neighbours
        for (CGMesh::VertexOHalfedgeIter voh_it = mesh.voh_begin(vh); voh_it != mesh.voh_end(vh); ++voh_it)
        {
            CGMesh::HalfedgeHandle hh = voh_it.handle();
            int idout = mesh.to_vertex_handle(hh).idx();
            if (idout == pin[i])
            {
                std::cout << "from and to are the same!\n";
                exit(0);
            }
            if(Vflags[idout]) // this vertex is inside the boundary loop, add it into the neighbours
            {
                nbs.push_back(idout);
            }
            else // this vertex is outside, then check which halfedge it is.
            {
                int vid = whichHalfedge(pin[i], idout, mesh, ehs) + vnbr;
                nbs.push_back(vid);
            }
        }
        if(nbs.size() == 4) // record only the valence 4 vertices
            connect.push_back({pin[i], nbs[0], nbs[1], nbs[2], nbs[3]});
    }
}

Eigen::MatrixXd connectVertexList(const Eigen::MatrixXd&V0, const std::vector<Eigen::Vector3d>& V1)
{
    int vnbr = V0.rows()+V1.size();
    Eigen::MatrixXd V(vnbr, 3);
    V.topRows(V0.rows()) = V0;
    V.bottomRows(V1.size()) = vec_list_to_matrix(V1);
    return V;
}

void remapVertices(const int vnbr, const int vnbrAll, const Eigen::MatrixXd &Vall,
                   const std::vector<bool> &flags, std::vector<int> &map, Eigen::MatrixXd &Vclean)
{
    int counter = 0;
    map = std::vector<int>(vnbrAll, -1);
    for (int i = 0; i < flags.size(); i++)
    {
        if (flags[i])
        {
            map[i] = counter;
            counter++;
        }
    }
    for (int i = flags.size(); i < vnbrAll; i++)
    {
        map[i] = counter;
        counter++;
    }
    Vclean.resize(counter, 3);
    for(int i=0;i<map.size();i++)
    {
        if(map[i]>=0)
        {
            Vclean.row(map[i]) = Vall.row(i);
        }
    }
}

void cleanFaceList(const std::vector<int> &map, const std::vector<int> &quadFull, const Eigen::MatrixXi &F,
                   std::vector<std::array<int,4>> &Fclean)
{
    Fclean.resize(quadFull.size());
    for (int i = 0; i < quadFull.size(); i++)
    {
        Eigen::Vector4i qd = F.row(quadFull[i]);
        Fclean[i] ={ map[qd[0]], map[qd[1]],
            map[qd[2]], map[qd[3]]};
        if (map[qd[0]] < 0 || map[qd[1]] < 0 ||
            map[qd[2]] < 0 || map[qd[3]] < 0)
        {
            std::cout<<"wrong in map!\n";
            exit(0);
        }
    }

}
template <typename T>
void cleanFaceList(const std::vector<int> &map, std::vector<T>& faces, const int order)
{
    std::vector<T> fcp = faces;
    for(int i=0;i<faces.size();i++)
    {
        for(int j=0;j<order;j++)
        {
            fcp[i][j] = map[faces[i][j]];
            if(map[faces[i][j]]<0)
            {
                std::cout<<"error in templated face cleaner!\n";
                exit(0);
            }
        }
    }
    faces = fcp;
}
void cleanLoopList(const std::vector<int> &map, std::vector<int>& loop)
{
    std::vector<int> fcp = loop;
    for (int i = 0; i < loop.size(); i++)
    {
        fcp[i] = map[loop[i]];
    }
    loop = fcp;
}

template <typename T>
void writeObjPly(const std::string filename, const Eigen::MatrixXd &V, const std::vector<int> &vmap,
                 const std::vector<T> &F, const int order)
{
    std::ofstream fout;
    fout.open(filename);
    for (int i = 0; i < V.rows(); i++)
    {
        if (vmap[i] >= 0)
            fout << "v " << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << "\n";
    }
    for (int i = 0; i < F.size(); i++)
    {
        fout << "f";
       
            for (int j = 0; j < order; j++)
            {
                fout << " " << F[i][j] + 1;
            }
        fout << "\n";
    }
    fout.close();
}
void writeObjPly(const std::string filename, const Eigen::MatrixXd &V, const std::vector<int> &vmap,
                 const std::vector<int> &F)
{
    std::ofstream fout;
    fout.open(filename);
    for (int i = 0; i < V.rows(); i++)
    {
        if (vmap[i] >= 0)
            fout << "v " << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << "\n";
    }
    for (int i = 0; i < F.size(); i++)
    {
        fout << "f";

        fout << "," << F[i] + 1;

        fout << "\n";
    }
    fout.close();
}

void writeObjLoop(const std::string filename, const Eigen::MatrixXd &V, const std::vector<int> &vmap,
                  const std::vector<int> &Loop)
{
    std::ofstream fout;
    fout.open(filename);
    for (int i = 0; i < V.rows(); i++)
    {
        if (vmap[i] >= 0)
            fout << "v " << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << "\n";
    }
    for (int i = 0; i < Loop.size() - 1; i++)
    {
        fout << "l " << Loop[i] + 1 << " " << Loop[i + 1] + 1;

        fout << "\n";
    }
    fout << "l " << Loop.back() + 1 << " " << Loop.front() + 1;
    fout.close();
}

void writeSomePoints(const std::string filename, const Eigen::MatrixXd &V, const std::vector<int> &vmap,
                  const std::vector<std::array<int, 5>> &connect, int cid)
{
    std::ofstream fout;
    fout.open(filename);
    for (int i = 0; i < V.rows(); i++)
    {
        if (vmap[i] >= 0)
            fout << "v " << V(i, 0) << " " << V(i, 1) << " " << V(i, 2) << "\n";
    }
    fout<<"f "<<connect[cid][0]+1<<" "<<connect[cid][1]+1<<" "<<connect[cid][2]+1<<"\n";
    fout<<"f "<<connect[cid][0]+1<<" "<<connect[cid][3]+1<<" "<<connect[cid][4]+1<<"\n";
}

 
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
    std::vector<Eigen::Vector3d> intersections;
    // find all the intersections on the edges.
    findALlIntersections(mesh, Vquad, Fquad, Vcurve, ehs, intersections);
    std::cout << "find intersection point nbr, " << ehs.size() << "\n";
    // get the cutted polygons and the boundary loop
    // the triangles, the quads and the pentagons
    std::vector<std::array<int,3>> rT;
    std::vector<std::array<int,4>> rQ;
    std::vector<std::array<int,5>> rP;
    std::vector<int> bLoop;
    splitIntersectedQuads(mesh, Vquad, Fquad, Vcurve, ehs, quadCut, fvin, vFlags, rQ, rT, rP, bLoop);
    std::cout<<"new generated quads, "<<rQ.size()<<", trians, "<<rT.size()<<", pentas, "<<rP.size()<<"\n";
    std::vector<std::array<int,5>> connect;
    // get the connection for all the inner vertices.
    findConnectivityForInnerVertices(vFlags, pin, vnbr, mesh, ehs, connect);
    Eigen::MatrixXd Vall = connectVertexList(Vquad, intersections);

    // clean the data
    std::vector<int> map;
    Eigen::MatrixXd Vclean;
    remapVertices(vnbr, Vall.rows(), Vall, vFlags, map, Vclean);
    std::vector<std::array<int,4>> quadsClean;
    // clean the face lists
    cleanFaceList(map, quadFull, Fquad, quadsClean);
    cleanFaceList(map, rQ, 4);
    cleanFaceList(map, rT, 3);
    cleanFaceList(map, rP, 5);
    // clean the connectivity
    cleanFaceList(map, connect, 5);
    // clean the boundary loop
    cleanLoopList(map, bLoop);
    // write them
    writeObjPly("/Users/wangb0d/Desktop/tmp/Khusrav/CRPC/send/rQF.obj", Vall, map, quadsClean, 4);
    writeObjPly("/Users/wangb0d/Desktop/tmp/Khusrav/CRPC/send/rQ.obj", Vall, map, rQ, 4);
    writeObjPly("/Users/wangb0d/Desktop/tmp/Khusrav/CRPC/send/rT.obj", Vall, map, rT, 3);
    writeObjPly("/Users/wangb0d/Desktop/tmp/Khusrav/CRPC/send/rP.obj", Vall, map, rP, 5);
    writeObjPly("/Users/wangb0d/Desktop/tmp/Khusrav/CRPC/send/connect.txt", Vall, map, connect, 5);
    writeObjLoop("/Users/wangb0d/Desktop/tmp/Khusrav/CRPC/send/loop.obj", Vall, map, bLoop);
    writeObjPly("/Users/wangb0d/Desktop/tmp/Khusrav/CRPC/send/loop.txt", Vall, map, bLoop);
    std::cout<<"connection size "<<connect.size()<<"\n";
    writeSomePoints("/Users/wangb0d/Desktop/tmp/Khusrav/CRPC/send/test.obj", Vall, map, connect, 766);

    // test and plot the connectivity



        // write the points inside
    // igl::writeOBJ("/Users/wangb0d/Desktop/tmp/Khusrav/CRPC/send/pIn.obj", Vin, Eigen::MatrixXi());
}
