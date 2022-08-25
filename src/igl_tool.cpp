// This file is mostly from libigl, principal_curvature.cpp. We borrow some functions to get I and II fundamental
// forms of the triangle mesh.

// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include <lsc/igl_tool.h>
#include <lsc/basic.h>
#include <lsc/tools.h>

class comparer
{
public:
  IGL_INLINE bool operator()(const std::pair<int, double> &lhs, const std::pair<int, double> &rhs) const
  {
    return lhs.second > rhs.second;
  }
};

IGL_INLINE QuadricCalculator::QuadricCalculator()
{
  this->localMode = true;
  this->projectionPlaneCheck = true;
  this->sphereRadius = 5;
  this->st = SPHERE_SEARCH;
  this->nt = AVERAGE;
  this->montecarlo = false;
  this->montecarloN = 0;
  this->kRing = 3;
  this->curvatureComputed = false;
  this->expStep = true;
}

IGL_INLINE void QuadricCalculator::init(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F)
{
  // Normalize vertices
  vertices = V;

  //  vertices = vertices.array() - vertices.minCoeff();
  //  vertices = vertices.array() / vertices.maxCoeff();
  //  vertices = vertices.array() * (1.0/igl::avg_edge_length(V,F));

  faces = F;
  igl::adjacency_list(F, vertex_to_vertices);
  igl::vertex_triangle_adjacency(V, F, vertex_to_faces, vertex_to_faces_index);
  igl::per_face_normals(V, F, face_normals);
  igl::per_vertex_normals(V, F, face_normals, vertex_normals);
}

IGL_INLINE void QuadricCalculator::fitQuadric(const Eigen::Vector3d &v, const std::vector<Eigen::Vector3d> &ref, const std::vector<int> &vv, Quadric *q)
{
  std::vector<Eigen::Vector3d> points;
  points.reserve(vv.size());

  for (unsigned int i = 0; i < vv.size(); ++i)
  {

    Eigen::Vector3d cp = vertices.row(vv[i]);

    // vtang non e` il v tangente!!!
    Eigen::Vector3d vTang = cp - v;

    double x = vTang.dot(ref[0]);
    double y = vTang.dot(ref[1]);
    double z = vTang.dot(ref[2]);
    points.push_back(Eigen::Vector3d(x, y, z));
  }
  if (points.size() < 5)
  {
    std::cerr << "ASSERT FAILED! fit function requires at least 5 points: Only " << points.size() << " were given." << std::endl;
    *q = Quadric(0, 0, 0, 0, 0);
  }
  else
  {
    *q = Quadric::fit(points);
  }
}

IGL_INLINE void QuadricCalculator::finalEigenStuff(int i, const std::vector<Eigen::Vector3d> &ref, Quadric &q)
{

  const double a = q.a();
  const double b = q.b();
  const double c = q.c();
  const double d = q.d();
  const double e = q.e();

  //  if (fabs(a) < 10e-8 || fabs(b) < 10e-8)
  //  {
  //    std::cout << "Degenerate quadric: " << i << std::endl;
  //  }

  double E = 1.0 + d * d;
  double F = d * e;
  double G = 1.0 + e * e;

  Eigen::Vector3d n = Eigen::Vector3d(-d, -e, 1.0).normalized();

  double L = 2.0 * a * n[2];
  double M = b * n[2];
  double N = 2 * c * n[2];

  // ----------------- Eigen stuff
  Eigen::Matrix2d m;
  m << L * G - M * F, M * E - L * F, M * E - L * F, N * E - M * F;
  m = m / (E * G - F * F);
  Eigen::SelfAdjointEigenSolver<Eigen::Matrix2d> eig(m);

  Eigen::Vector2d c_val = eig.eigenvalues();
  Eigen::Matrix2d c_vec = eig.eigenvectors();

  // std::cerr << "c_val:" << c_val << std::endl;
  // std::cerr << "c_vec:" << c_vec << std::endl;

  // std::cerr << "c_vec:" << c_vec(0) << " "  << c_vec(1) << std::endl;

  c_val = -c_val;

  Eigen::Vector3d v1, v2;
  v1[0] = c_vec(0);
  v1[1] = c_vec(1);
  v1[2] = 0; // d * v1[0] + e * v1[1];

  v2[0] = c_vec(2);
  v2[1] = c_vec(3);
  v2[2] = 0; // d * v2[0] + e * v2[1];

  // v1 = v1.normalized();
  // v2 = v2.normalized();

  Eigen::Vector3d v1global = ref[0] * v1[0] + ref[1] * v1[1] + ref[2] * v1[2];
  Eigen::Vector3d v2global = ref[0] * v2[0] + ref[1] * v2[1] + ref[2] * v2[2];

  v1global.normalize();
  v2global.normalize();

  v1global *= c_val(0);
  v2global *= c_val(1);

  if (c_val[0] > c_val[1])
  {
    curv[i] = std::vector<double>(2);
    curv[i][0] = c_val(1);
    curv[i][1] = c_val(0);
    curvDir[i] = std::vector<Eigen::Vector3d>(2);
    curvDir[i][0] = v2global;
    curvDir[i][1] = v1global;
  }
  else
  {
    curv[i] = std::vector<double>(2);
    curv[i][0] = c_val(0);
    curv[i][1] = c_val(1);
    curvDir[i] = std::vector<Eigen::Vector3d>(2);
    curvDir[i][0] = v1global;
    curvDir[i][1] = v2global;
  }
  // ---- end Eigen stuff
}

IGL_INLINE void QuadricCalculator::getKRing(const int start, const double r, std::vector<int> &vv)
{
  int bufsize = vertices.rows();
  vv.reserve(bufsize);
  std::list<std::pair<int, int>> queue;
  std::vector<bool> visited(bufsize, false);
  queue.push_back(std::pair<int, int>(start, 0));
  visited[start] = true;
  while (!queue.empty())
  {
    int toVisit = queue.front().first;
    int distance = queue.front().second;
    queue.pop_front();
    vv.push_back(toVisit);
    if (distance < (int)r)
    {
      for (unsigned int i = 0; i < vertex_to_vertices[toVisit].size(); ++i)
      {
        int neighbor = vertex_to_vertices[toVisit][i];
        if (!visited[neighbor])
        {
          queue.push_back(std::pair<int, int>(neighbor, distance + 1));
          visited[neighbor] = true;
        }
      }
    }
  }
}

IGL_INLINE void QuadricCalculator::getSphere(const int start, const double r, std::vector<int> &vv, int min)
{
  int bufsize = vertices.rows();
  vv.reserve(bufsize);
  std::list<int> queue;
  std::vector<bool> visited(bufsize, false);
  queue.push_back(start);
  visited[start] = true;
  Eigen::Vector3d me = vertices.row(start);
  std::priority_queue<std::pair<int, double>, std::vector<std::pair<int, double>>, comparer> extra_candidates;
  while (!queue.empty())
  {
    int toVisit = queue.front();
    queue.pop_front();
    vv.push_back(toVisit);
    for (unsigned int i = 0; i < vertex_to_vertices[toVisit].size(); ++i)
    {
      int neighbor = vertex_to_vertices[toVisit][i];
      if (!visited[neighbor])
      {
        Eigen::Vector3d neigh = vertices.row(neighbor);
        double distance = (me - neigh).norm();
        if (distance < r)
          queue.push_back(neighbor);
        else if ((int)vv.size() < min)
          extra_candidates.push(std::pair<int, double>(neighbor, distance));
        visited[neighbor] = true;
      }
    }
  }
  while (!extra_candidates.empty() && (int)vv.size() < min)
  {
    std::pair<int, double> cand = extra_candidates.top();
    extra_candidates.pop();
    vv.push_back(cand.first);
    for (unsigned int i = 0; i < vertex_to_vertices[cand.first].size(); ++i)
    {
      int neighbor = vertex_to_vertices[cand.first][i];
      if (!visited[neighbor])
      {
        Eigen::Vector3d neigh = vertices.row(neighbor);
        double distance = (me - neigh).norm();
        extra_candidates.push(std::pair<int, double>(neighbor, distance));
        visited[neighbor] = true;
      }
    }
  }
}

IGL_INLINE Eigen::Vector3d QuadricCalculator::project(const Eigen::Vector3d &v, const Eigen::Vector3d &vp, const Eigen::Vector3d &ppn)
{
  return (vp - (ppn * ((vp - v).dot(ppn))));
}

IGL_INLINE void QuadricCalculator::computeReferenceFrame(int i, const Eigen::Vector3d &normal, std::vector<Eigen::Vector3d> &ref)
{

  Eigen::Vector3d longest_v = Eigen::Vector3d(vertices.row(vertex_to_vertices[i][0]));

  longest_v = (project(vertices.row(i), longest_v, normal) - Eigen::Vector3d(vertices.row(i))).normalized();

  /* L'ultimo asse si ottiene come prodotto vettoriale tra i due
   * calcolati */
  Eigen::Vector3d y_axis = (normal.cross(longest_v)).normalized();
  ref[0] = longest_v;
  ref[1] = y_axis;
  ref[2] = normal;
}

IGL_INLINE void QuadricCalculator::getAverageNormal(int j, const std::vector<int> &vv, Eigen::Vector3d &normal)
{
  normal = (vertex_normals.row(j)).normalized();
  if (localMode)
    return;

  for (unsigned int i = 0; i < vv.size(); ++i)
  {
    normal += vertex_normals.row(vv[i]).normalized();
  }
  normal.normalize();
}

IGL_INLINE void QuadricCalculator::getProjPlane(int j, const std::vector<int> &vv, Eigen::Vector3d &ppn)
{
  int nr;
  double a, b, c;
  double nx, ny, nz;
  double abcq;

  a = b = c = 0;

  if (localMode)
  {
    for (unsigned int i = 0; i < vertex_to_faces.at(j).size(); ++i)
    {
      Eigen::Vector3d faceNormal = face_normals.row(vertex_to_faces.at(j).at(i));
      a += faceNormal[0];
      b += faceNormal[1];
      c += faceNormal[2];
    }
  }
  else
  {
    for (unsigned int i = 0; i < vv.size(); ++i)
    {
      a += vertex_normals.row(vv[i])[0];
      b += vertex_normals.row(vv[i])[1];
      c += vertex_normals.row(vv[i])[2];
    }
  }
  nr = rotateForward(&a, &b, &c);
  abcq = a * a + b * b + c * c;
  nx = sqrt(a * a / abcq);
  ny = sqrt(b * b / abcq);
  nz = sqrt(1 - nx * nx - ny * ny);
  rotateBackward(nr, &a, &b, &c);
  rotateBackward(nr, &nx, &ny, &nz);

  ppn = chooseMax(Eigen::Vector3d(nx, ny, nz), Eigen::Vector3d(a, b, c), a * b);
  ppn.normalize();
}

IGL_INLINE double QuadricCalculator::getAverageEdge()
{
  double sum = 0;
  int count = 0;

  for (int i = 0; i < faces.rows(); ++i)
  {
    for (short unsigned j = 0; j < 3; ++j)
    {
      Eigen::Vector3d p1 = vertices.row(faces.row(i)[j]);
      Eigen::Vector3d p2 = vertices.row(faces.row(i)[(j + 1) % 3]);

      double l = (p1 - p2).norm();

      sum += l;
      ++count;
    }
  }

  return (sum / (double)count);
}

IGL_INLINE void QuadricCalculator::applyProjOnPlane(const Eigen::Vector3d &ppn, const std::vector<int> &vin, std::vector<int> &vout)
{
  for (std::vector<int>::const_iterator vpi = vin.begin(); vpi != vin.end(); ++vpi)
    if (vertex_normals.row(*vpi) * ppn > 0.0)
      vout.push_back(*vpi);
}

IGL_INLINE void QuadricCalculator::applyMontecarlo(const std::vector<int> &vin, std::vector<int> *vout)
{
  if (montecarloN >= vin.size())
  {
    *vout = vin;
    return;
  }

  float p = ((float)montecarloN) / (float)vin.size();
  for (std::vector<int>::const_iterator vpi = vin.begin(); vpi != vin.end(); ++vpi)
  {
    float r;
    if ((r = ((float)rand() / RAND_MAX)) < p)
    {
      vout->push_back(*vpi);
    }
  }
}

IGL_INLINE void QuadricCalculator::computeCurvature()
{
  // CHECK che esista la mesh
  const size_t vertices_count = vertices.rows();

  if (vertices_count == 0)
    return;

  curvDir = std::vector<std::vector<Eigen::Vector3d>>(vertices_count);
  curv = std::vector<std::vector<double>>(vertices_count);

  scaledRadius = getAverageEdge() * sphereRadius;

  std::vector<int> vv;
  std::vector<int> vvtmp;
  Eigen::Vector3d normal;

  // double time_spent;
  // double searchtime=0, ref_time=0, fit_time=0, final_time=0;

  for (size_t i = 0; i < vertices_count; ++i)
  {
    vv.clear();
    vvtmp.clear();
    Eigen::Vector3d me = vertices.row(i);
    switch (st)
    {
    case SPHERE_SEARCH:
      getSphere(i, scaledRadius, vv, 6);
      break;
    case K_RING_SEARCH:
      getKRing(i, kRing, vv);
      break;
    default:
      fprintf(stderr, "Error: search type not recognized");
      return;
    }

    if (vv.size() < 6)
    {
      // std::cerr << "Could not compute curvature of radius " << scaledRadius << std::endl;
      continue;
    }

    if (projectionPlaneCheck)
    {
      vvtmp.reserve(vv.size());
      applyProjOnPlane(vertex_normals.row(i), vv, vvtmp);
      if (vvtmp.size() >= 6 && vvtmp.size() < vv.size())
        vv = vvtmp;
    }

    switch (nt)
    {
    case AVERAGE:
      getAverageNormal(i, vv, normal);
      break;
    case PROJ_PLANE:
      getProjPlane(i, vv, normal);
      break;
    default:
      fprintf(stderr, "Error: normal type not recognized");
      return;
    }
    if (vv.size() < 6)
    {
      // std::cerr << "Could not compute curvature of radius " << scaledRadius << std::endl;
      continue;
    }
    if (montecarlo)
    {
      if (montecarloN < 6)
        break;
      vvtmp.reserve(vv.size());
      applyMontecarlo(vv, &vvtmp);
      vv = vvtmp;
    }

    if (vv.size() < 6)
      return;
    std::vector<Eigen::Vector3d> ref(3);
    computeReferenceFrame(i, normal, ref);

    Quadric q;
    fitQuadric(me, ref, vv, &q);
    finalEigenStuff(i, ref, q);
  }

  lastRadius = sphereRadius;
  curvatureComputed = true;
}

IGL_INLINE void QuadricCalculator::get_I_and_II(Eigen::MatrixXd &ru, Eigen::MatrixXd &rv, Eigen::MatrixXd &ruu, Eigen::MatrixXd &ruv,
                                                Eigen::MatrixXd &rvv, std::vector<double> &L_list, std::vector<double> &M_list, std::vector<double> &N_list, Eigen::MatrixXd &normals)
{
  // CHECK che esista la mesh
  const size_t vertices_count = vertices.rows();

  if (vertices_count == 0)
    return;

  // Eigen::MatrixXd ru, rv,ruu,rvv,ruv;
  // std::vector<double> L_list, M_list, N_list;
  // Eigen::Vector3d local_ru, local_rv, local_ruu, local_rvv, local_ruv;
  ru.resize(vertices_count, 3);
  rv.resize(vertices_count, 3);
  ruv.resize(vertices_count, 3);
  ruu.resize(vertices_count, 3);
  rvv.resize(vertices_count, 3);
  L_list.resize(vertices_count);
  M_list.resize(vertices_count);
  N_list.resize(vertices_count);
  normals.resize(vertices_count, 3);
  std::vector<int> vv;
  std::vector<int> vvtmp;
  Eigen::Vector3d normal;

  for (size_t i = 0; i < vertices_count; ++i)
  {
    vv.clear();
    vvtmp.clear();
    Eigen::Vector3d me = vertices.row(i);
    getKRing(i, kRing, vv);
    if (vv.size() < 6)
    {
      std::cerr << "Could not compute I and II using libigl, because too few neighbouring vertices " << std::endl;
      continue;
    }

    if (projectionPlaneCheck)
    {
      vvtmp.reserve(vv.size());
      applyProjOnPlane(vertex_normals.row(i), vv, vvtmp);
      if (vvtmp.size() >= 6 && vvtmp.size() < vv.size())
        vv = vvtmp;
    }

    switch (nt)
    {
    case AVERAGE:
      getAverageNormal(i, vv, normal); // this normal is just to give the z axis of the frame, thus allows some error
      break;
    case PROJ_PLANE:
      getProjPlane(i, vv, normal);
      break;
    default:
      fprintf(stderr, "Error: normal type not recognized");
      return;
    }

    if (vv.size() < 6)
    {
      std::cerr << "Could not compute I and II using libigl, because too few neighbouring vertices " << std::endl;
      return;
    }

    std::vector<Eigen::Vector3d> ref(3);
    computeReferenceFrame(i, normal, ref);

    Quadric q;
    fitQuadric(me, ref, vv, &q);
    const double a = q.a();
    const double b = q.b();
    const double c = q.c();
    const double d = q.d();
    const double e = q.e();

    //  if (fabs(a) < 10e-8 || fabs(b) < 10e-8)
    //  {
    //    std::cout << "Degenerate quadric: " << i << std::endl;
    //  }
    //   ru.row(i) << 1, 0, d;
    //   rv.row(i) << 0, 1, e;
    //   ruu.row(i) << 0, 0, 2 * a;
    //   ruv.row(i) << 0, 0, b;
    //   rvv.row(i) << 0, 0, 2 * c;
    ru.row(i) = ref[0] + ref[2] * d;
    rv.row(i) = ref[1] + ref[2] * e;
    ruu.row(i) = 2 * a * ref[2];
    ruv.row(i) = b * ref[2];
    rvv.row(i) = 2 * c * ref[2];

    Eigen::Vector3d n = Eigen::Vector3d(-d, -e, 1.0).normalized();
    normals.row(i) = ref[0] * n[0] + ref[1] * n[1] + ref[2] * n[2];
    double L = 2.0 * a * n[2];
    double M = b * n[2];
    double N = 2 * c * n[2];
    L_list[i] = L;
    M_list[i] = M;
    N_list[i] = N;
  }

  lastRadius = sphereRadius;
}

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

IGL_INLINE void QuadricCalculator::printCurvature(const std::string &outpath)
{
  using namespace std;
  if (!curvatureComputed)
    return;

  std::ofstream of;
  of.open(outpath.c_str());

  if (!of)
  {
    fprintf(stderr, "Error: could not open output file %s\n", outpath.c_str());
    return;
  }

  int vertices_count = vertices.rows();
  of << vertices_count << endl;
  for (int i = 0; i < vertices_count; ++i)
  {
    of << curv[i][0] << " " << curv[i][1] << " " << curvDir[i][0][0] << " " << curvDir[i][0][1] << " " << curvDir[i][0][2] << " " << curvDir[i][1][0] << " " << curvDir[i][1][1] << " " << curvDir[i][1][2] << endl;
  }

  of.close();
}

void lsTools::get_I_and_II_locally()
{
  unsigned radius = 2;

  // Precomputation
  QuadricCalculator cc;
  cc.init(V, F);

  cc.kRing = radius;
  cc.st = K_RING_SEARCH;

  // Compute
  cc.get_I_and_II(Deriv1[0], Deriv1[1], Deriv2[0], Deriv2[1], Deriv2[3], II_L, II_M, II_N, norm_v);
  Deriv2[2] = Deriv2[1]; // now r_uv = r_vu
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
    edge_out=checking_he;
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
    edge_out=checking_he;
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
    }
  }
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
    for(int i=0;i<temp_pend.size();i++){
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
    for(int i=0;i<temp_pend.size();i++){
      edge_out.push_back(checking_he);
      p_end.push_back(temp_pend[i]);
    }
  }
  return p_end.size();
}

void extend_vector(const std::vector<Eigen::Vector3d>& first, const std::vector<Eigen::Vector3d>& second, std::vector<Eigen::Vector3d>& result){
  result=first.insert(first.end(),second.begin(),second.end());
}
void extend_vector(const std::vector<CGMesh::HalfedgeHandle>& first, const std::vector<CGMesh::HalfedgeHandle>& second, std::vector<CGMesh::HalfedgeHandle>& result){
  result=first.insert(first.end(),second.begin(),second.end());
}
template<typename Tp>
void vector_duplicate(const Tp& input, const int size, std::vector<Tp>& result){
  result.resize(size);
  for(int i=0;i<size;i++){
    result[i]=input;
  }
}

// caution: this is not an initialization to find the first and second point
void pseudo_geodesic_intersection_filter(
    const std::vector<Eigen::Vector3d> &curve, const double angle_degree,
    const Eigen::Vector3d &pnorm,
    const std::vector<Eigen::Vector3d> &candi_points, int &id)
{
  assert(curve.size()>=2);// take previous points to find the next one
  assert(candi_points.size()>0);
  if(candi_points.size()==1){
    id=0;
    return;
  }
  int size=curve.size();
  Eigen::Vector3d dire1 = (curve[size - 1] - curve[size - 2]).normalized();
  int closest_id=-1;
  // TODO should we consider the angle consistancy? to avoid cases like 5 degree flip to 355 degree?
  // int closest_id_consider_angle=-1;
  double closest_angle_diff_radian = 370. * LSC_PI / 180.;
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
  for(int i=0;i<candi_points.size();i++){
    Eigen::Vector3d dire2= (candi_points[i]-curve[size - 1]).normalized();
    Eigen::Vector3d normal=dire1.cross(dire2).normalized();
    double normal_pnormal=normal.dot(pnorm);
    double dire1_normal=dire1.dot(normal);
    double angle_tmp_radian = acos(normal_pnormal);
    double angle_diff = std::min(fabs(angle_tmp_radian - angle_radian), fabs(2 * LSC_PI - angle_tmp_radian - angle_radian));
    if (angle_diff < closest_angle_diff_radian)
    {
      closest_angle_diff_radian = angle_diff;
      closest_id = i;
    }
  }
  // already found the vertex that forms the angle most close to the given one
  id=closest_id;
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
// angle is in degree, when degree is 90, is a geodesic.
bool lsTools::get_pseudo_vertex_and_trace_forward(
    QuadricCalculator &cc,
    const std::vector<Eigen::Vector3d> &curve, const double angle_degree,
    const CGMesh::HalfedgeHandle &edge_middle,
    const Eigen::Vector3d &point_in, const Eigen::Vector3d &point_middle, CGMesh::HalfedgeHandle &edge_out,
    Eigen::Vector3d &point_out, bool &generate_pseudo_vertex, Eigen::Vector3d &pseudo_vertex_out)
{
  unsigned radius = 2;
  bool is_geodesic=false;
  generate_pseudo_vertex=false;
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
  if(angle_degree>90-ANGLE_TOLERANCE||angle_degree<90+ANGLE_TOLERANCE){// it means it is a geodesic
        is_geodesic=true;
  }
  if (from_ratio <= MERGE_VERTEX_RATIO || to_ratio <= MERGE_VERTEX_RATIO)// the point_middle is a vertex
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
    Eigen::Vector3d pnorm=norm_v.row(center_handle.idx());
    for (CGMesh::VertexOHalfedgeIter voh_itr = lsmesh.voh_begin(center_handle);
         voh_itr != lsmesh.voh_end(center_handle); ++voh_itr)
    {
      auto heh = voh_itr.handle();
      CGMesh::HalfedgeHandle edge_to_check = lsmesh.next_halfedge_handle(heh);
      assert(lsmesh.from_vertex_handle(edge_to_check).idx() != center_handle.idx());
      assert(lsmesh.to_vertex_handle(edge_to_check).idx() != center_handle.idx());
      Eigen::Vector3d vs=V.row(lsmesh.from_vertex_handle(edge_to_check).idx());
      Eigen::Vector3d ve=V.row(lsmesh.to_vertex_handle(edge_to_check).idx());
      
      if(is_geodesic){// it means it is a geodesic
        bool found=find_geodesic_intersection_p1_is_ver(point_in,point_middle,vs,ve,pnorm,point_out);
        if(found){
          edge_out=edge_to_check;
          return true;
        }
      }// geodesic end
      else{
        std::vector<Eigen::Vector3d> tmp_candidates;
        std::vector<CGMesh::HalfedgeHandle> tmp_hehs;
        bool found = find_osculating_plane_intersection_not_geodesic_p1_is_ver(point_in,point_middle,vs,ve,pnorm,angle_radian,tmp_candidates);
        if(found){
          extend_vector(candidate_pts,tmp_candidates,candidate_pts);
          vector_duplicate(edge_to_check,tmp_candidates.size(),tmp_hehs);
          extend_vector(candidate_handles,tmp_hehs,candidate_handles);
        }
      }
    }// end of searching all the one-ring edges
    if(is_geodesic){
    return false;
    }
    else{
      if(candidate_pts.size()==0){
        return false;
      }
      // pick the point
      int id;
      pseudo_geodesic_intersection_filter(curve,angle_degree,pnorm,candidate_pts,id);
      edge_out=candidate_handles[id];
      point_out=candidate_pts[id];
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
  cc.get_pseudo_vertex_on_edge(cent_id, point_middle, pnorm, pver);
  generate_pseudo_vertex = true; // TODO add tolerance
  pseudo_vertex_out = pver;
  if (is_geodesic)
  {
    bool found = find_geodesic_intersection_p1_is_NOT_ver(point_in, pver, edge_middle, pnorm, edge_out, point_out);
    return found;
  }
  else{
    bool found = find_osculating_plane_intersection_not_geodesic_p1_is_not_ver(point_in, pver, edge_middle, pnorm,
                                                                               angle_radian, candidate_handles, candidate_pts);
    if(candidate_pts.size()==0){
        return false;
      }
      // pick the point
      int id;
      pseudo_geodesic_intersection_filter(curve,angle_degree,pnorm,candidate_pts,id);
      edge_out=candidate_handles[id];
      point_out=candidate_pts[id];
      return true;
  }
  return false;
}