
#pragma once
#include <iostream>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <queue>
#include <list>
#include <cmath>
#include <limits>

#include <Eigen/SparseCholesky>

// Lib IGL includes
#include <igl/adjacency_list.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/avg_edge_length.h>
#include <igl/vertex_triangle_adjacency.h>


typedef enum
{
  SPHERE_SEARCH,
  K_RING_SEARCH
} searchType;

typedef enum
{
  AVERAGE,
  PROJ_PLANE
} normalType;
class QuadricCalculator
{
public:
  /* Row number i represents the i-th vertex, whose columns are:
   curv[i][0] : K2
   curv[i][1] : K1
   curvDir[i][0] : PD1
   curvDir[i][1] : PD2
   */
  std::vector< std::vector<double> > curv;
  std::vector< std::vector<Eigen::Vector3d> > curvDir;
  bool curvatureComputed;
  class Quadric
  {
  public:

    IGL_INLINE Quadric ()
    {
      a() = b() = c() = d() = e() = 1.0;
    }

    IGL_INLINE Quadric(double av, double bv, double cv, double dv, double ev)
    {
      a() = av;
      b() = bv;
      c() = cv;
      d() = dv;
      e() = ev;
    }

    IGL_INLINE double& a() { return data[0];}
    IGL_INLINE double& b() { return data[1];}
    IGL_INLINE double& c() { return data[2];}
    IGL_INLINE double& d() { return data[3];}
    IGL_INLINE double& e() { return data[4];}

    double data[5];

    IGL_INLINE double evaluate(double u, double v)
    {
      return a()*u*u + b()*u*v + c()*v*v + d()*u + e()*v;
    }

    IGL_INLINE double du(double u, double v)
    {
      return 2.0*a()*u + b()*v + d();
    }

    IGL_INLINE double dv(double u, double v)
    {
      return 2.0*c()*v + b()*u + e();
    }

    IGL_INLINE double duv(double u, double v)
    {
      return b();
    }

    IGL_INLINE double duu(double u, double v)
    {
      return 2.0*a();
    }

    IGL_INLINE double dvv(double u, double v)
    {
      return 2.0*c();
    }


    IGL_INLINE static Quadric fit(const std::vector<Eigen::Vector3d> &VV)
    {
      assert(VV.size() >= 5);
      if (VV.size() < 5)
      {
        std::cerr << "ASSERT FAILED! fit function requires at least 5 points: Only " << VV.size() << " were given." << std::endl;
        exit(0);
      }

      Eigen::MatrixXd A(VV.size(),5);
      Eigen::MatrixXd b(VV.size(),1);
      Eigen::MatrixXd sol(5,1);

      for(unsigned int c=0; c < VV.size(); ++c)
      {
        double u = VV[c][0];
        double v = VV[c][1];
        double n = VV[c][2];

        A(c,0) = u*u;
        A(c,1) = u*v;
        A(c,2) = v*v;
        A(c,3) = u;
        A(c,4) = v;

        b(c) = n;
      }

      sol=A.jacobiSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(b);

      return Quadric(sol(0),sol(1),sol(2),sol(3),sol(4));
    }
  };

public:

  Eigen::MatrixXd vertices;
  // Face list of current mesh    (#F x 3) or (#F x 4)
  // The i-th row contains the indices of the vertices that forms the i-th face in ccw order
  Eigen::MatrixXi faces;

  std::vector<std::vector<int> > vertex_to_vertices;
  std::vector<std::vector<int> > vertex_to_faces;
  std::vector<std::vector<int> > vertex_to_faces_index;
  Eigen::MatrixXd face_normals;
  Eigen::MatrixXd vertex_normals;

  /* Size of the neighborhood */
  double sphereRadius;
  int kRing;

  bool localMode; /* Use local mode */
  bool projectionPlaneCheck; /* Check collected vertices on tangent plane */
  bool montecarlo;
  unsigned int montecarloN;

  searchType st; /* Use either a sphere search or a k-ring search */
  normalType nt;

  double lastRadius;
  double scaledRadius;
  std::string lastMeshName;

  /* Benchmark related variables */
  bool expStep; /* True if we want the radius to increase exponentially */
  int step;  /* If expStep==false, by how much rhe radius increases on every step */
  int maxSize; /* The maximum limit of the radius in the benchmark */

  QuadricCalculator();
  IGL_INLINE void init(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F);

  IGL_INLINE void finalEigenStuff(int, const std::vector<Eigen::Vector3d>&, Quadric&);
  IGL_INLINE void fitQuadric(const Eigen::Vector3d&, const std::vector<Eigen::Vector3d>& ref, const std::vector<int>& , Quadric *);
  IGL_INLINE void applyProjOnPlane(const Eigen::Vector3d&, const std::vector<int>&, std::vector<int>&);
  IGL_INLINE void getSphere(const int, const double, std::vector<int>&, int min);
  IGL_INLINE void getKRing(const int, const double,std::vector<int>&);
  IGL_INLINE Eigen::Vector3d project(const Eigen::Vector3d&, const Eigen::Vector3d&, const Eigen::Vector3d&);
  void computeReferenceFrame(int i, const Eigen::Vector3d& normal, std::vector<Eigen::Vector3d>& ref);
  IGL_INLINE void getAverageNormal(int, const std::vector<int>&, Eigen::Vector3d&);
  IGL_INLINE void getProjPlane(int, const std::vector<int>&, Eigen::Vector3d&);
  IGL_INLINE void applyMontecarlo(const std::vector<int>&,std::vector<int>*);
  IGL_INLINE void computeCurvature();
  // get the derivates, get the second fundamental form, and get the corresponding normals
  IGL_INLINE void get_I_and_II(Eigen::MatrixXd& ru, Eigen::MatrixXd& rv, Eigen::MatrixXd& ruu, Eigen::MatrixXd& ruv,
Eigen::MatrixXd& rvv, std::vector<double>& L_list,std::vector<double>& M_list, std::vector<double> &N_list,Eigen::MatrixXd &normals);
  IGL_INLINE void get_pseudo_vertex_on_edge(
    const int center_id,
    const Eigen::Vector3d &epoint, const Eigen::Vector3d &pnormal,Eigen::Vector3d &pvertex);
  IGL_INLINE void printCurvature(const std::string& outpath);
  IGL_INLINE double getAverageEdge();

  IGL_INLINE static int rotateForward (double *v0, double *v1, double *v2)
  {
    double t;

    if (std::abs(*v2) >= std::abs(*v1) && std::abs(*v2) >= std::abs(*v0))
      return 0;

    t = *v0;
    *v0 = *v2;
    *v2 = *v1;
    *v1 = t;

    return 1 + rotateForward (v0, v1, v2);
  }

  IGL_INLINE static void rotateBackward (int nr, double *v0, double *v1, double *v2)
  {
    double t;

    if (nr == 0)
      return;

    t = *v2;
    *v2 = *v0;
    *v0 = *v1;
    *v1 = t;

    rotateBackward (nr - 1, v0, v1, v2);
  }

  IGL_INLINE static Eigen::Vector3d chooseMax (Eigen::Vector3d n, Eigen::Vector3d abc, double ab)
  {
    int max_i;
    double max_sp;
    Eigen::Vector3d nt[8];

    n.normalize ();
    abc.normalize ();

    max_sp = - std::numeric_limits<double>::max();

    for (int i = 0; i < 4; ++i)
    {
      nt[i] = n;
      if (ab > 0)
      {
        switch (i)
        {
          case 0:
            break;

          case 1:
            nt[i][2] = -n[2];
            break;

          case 2:
            nt[i][0] = -n[0];
            nt[i][1] = -n[1];
            break;

          case 3:
            nt[i][0] = -n[0];
            nt[i][1] = -n[1];
            nt[i][2] = -n[2];
            break;
        }
      }
      else
      {
        switch (i)
        {
          case 0:
            nt[i][0] = -n[0];
            break;

          case 1:
            nt[i][1] = -n[1];
            break;

          case 2:
            nt[i][0] = -n[0];
            nt[i][2] = -n[2];
            break;

          case 3:
            nt[i][1] = -n[1];
            nt[i][2] = -n[2];
            break;
        }
      }

      if (nt[i].dot(abc) > max_sp)
      {
        max_sp = nt[i].dot(abc);
        max_i = i;
      }
    }
    return nt[max_i];
  }

};


