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


class comparer
{
public:
  IGL_INLINE bool operator()(const std::pair<int, double> &lhs, const std::pair<int, double> &rhs) const
  {
    return lhs.second > rhs.second;
  }
};

QuadricCalculator::QuadricCalculator()
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

void QuadricCalculator::computeReferenceFrame(int i, const Eigen::Vector3d &normal, std::vector<Eigen::Vector3d> &ref)
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
