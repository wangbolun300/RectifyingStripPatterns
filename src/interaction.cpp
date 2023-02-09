#include <lsc/interaction.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/heat_geodesics.h>


void collect_2d_positions(igl::opengl::glfw::Viewer &viewer, const double time_limit,
                          std::vector<double> &xs, std::vector<double> &ys)
{
    double time_interval = 0.2; // time interval in seconds
    int nbr_itr = time_limit / time_interval;
    xs.reserve(nbr_itr);
    ys.reserve(nbr_itr);
    for (int i = 0;; i++)
    {
        if (i >= nbr_itr)
        {
            break;
        }
        double x = viewer.current_mouse_x;
        double y = viewer.core().viewport(3) - viewer.current_mouse_y;
        xs.push_back(x);
        ys.push_back(y);
        std::cout<<"before getting one pt"<<std::endl;
        sleep(time_interval);
        std::cout<<" one pt got, ("<<x<<", "<<y<<")"<<std::endl;
    }
}

// this version only leave 2 vertices for each face. We give up this version since the density of sampling is often low.
// void get_stroke_intersection_on_edges(const CGMesh &mesh, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
//                                       const std::vector<int> &flist,
//                                       const std::vector<Eigen::Vector3f> &bclist)
// {
//     int fprev = -1;
//     int fsize = flist.size();
//     std::vector<int> ffilter;
//     std::vector<Eigen::Vector3f> bfil1, bfil2;
//     ffilter.reserve(fsize);
//     bfil1.reserve(fsize);
//     bfil2.reserve(fsize);
//     std::vector<std::vector<int>> fsort;
//     std::vector<int> ftmp; // record ids of flist
//     // get the segments intersect with edges
//     for (int i = 0; i < fsize; i++)
//     {
//         int fcurrent = flist[i];
//         if (fcurrent != fprev)
//         {
//             if (!ftmp.empty())
//             {
//                 ftmp.push_back(i - 1);
//                 fsort.push_back(ftmp);
//                 ftmp.clear();
//             }
//             ftmp.push_back(i);
//             fprev = fcurrent;
//         }
//     }
//     fsize = fsort.size();
//     for (int i = 0; i < fsize; i++)
//     {
//         for (int j = 0; j < fsort[i].size(); j++)
//         {
//             std::cout << fsort[i][j] << ", ";
//         }
//         std::cout << "\n";
//     }
// }

// filter to leave at most 1 point in each triangle.
void get_stroke_intersection_on_edges(const CGMesh &mesh, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                                      const std::vector<int> &flist,
                                      const std::vector<Eigen::Vector3f> &bclist, std::vector<int> &fout, 
                                      std::vector<Eigen::Vector3d>& bcout)
{
    int fprev = -1;
    int fsize = flist.size();

    std::vector<Eigen::Vector3f> bfil1, bfil2;
    bfil1.reserve(fsize);
    bfil2.reserve(fsize);
    std::vector<std::vector<int>> fsort;
    std::vector<std::vector<Eigen::Vector3d>> bcsort;
    std::vector<int> ftmp; // record ids
    std::vector<Eigen::Vector3d> bctmp;
    // get the segments intersect with edges
    for (int i = 0; i < fsize; i++)
    {
        int fcurrent = flist[i];
        Eigen::Vector3d bccurrent(bclist[i][0], bclist[i][1], bclist[i][2]);
        if (fcurrent != fprev)
        {
            if (!ftmp.empty())
            {
                fsort.push_back(ftmp);
                bcsort.push_back(bctmp);
                ftmp.clear();
                bctmp.clear();
            }
            fprev = fcurrent;
        }
        ftmp.push_back(fcurrent);
        bctmp.push_back(bccurrent);
    }
    fsort.push_back(ftmp);
    bcsort.push_back(bctmp);
    fsize = fsort.size();
    // std::cout << "the selected fids:" << std::endl;
    for (int i = 0; i < fsize; i++)
    {
        int middle;
        int tsize = fsort[i].size();
        if (tsize == 0)
        {
            std::cout << "error: no point is collected here" << std::endl;
            exit(0);
        }
        middle = tsize / 2;
        assert(middle >= 0 && middle < tsize);
        int select_id = fsort[i][middle];
        Eigen::Vector3d select_bc = bcsort[i][middle];
        fout.push_back(select_id);
        bcout.push_back(select_bc);
    }
}
void draw_stroke_on_mesh(const CGMesh &mesh, igl::opengl::glfw::Viewer &viewer,
                         const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const std::vector<double> &xs,
                         const std::vector<double> &ys, Eigen::MatrixXd& Vers)
{
    std::cout<<"drawing one curve, projection size "<<xs.size()<<std::endl;
    int ptsize = xs.size();
    assert(xs.size() == ys.size());
    std::vector<int> flist;
    std::vector<Eigen::Vector3f> bclist;
    flist.reserve(ptsize);
    bclist.reserve(ptsize);
    // get projection 3d points on the mesh
    for (int i = 0; i < ptsize; i++)
    {
        int fid;
        Eigen::Vector3f bc;
        double x = xs[i];
        double y = ys[i];
        bool proj = igl::unproject_onto_mesh(
            Eigen::Vector2f(x, y),
            viewer.core().view,
            viewer.core().proj,
            viewer.core().viewport,
            viewer.data().V,
            viewer.data().F,
            fid,
            bc);
        // std::cout<<"projected fid "<<fid<<std::endl;
        if (proj)
        {
            flist.push_back(fid);
            bclist.push_back(bc);
        }
    }
    std::vector<int> fout; 
    std::vector<Eigen::Vector3d> bcout;
    get_stroke_intersection_on_edges(mesh, V, F,
                                     flist,
                                     bclist, fout, bcout);
    // get 3d points according to the barycenter coordinates
    int vnbr = fout.size();
    Vers.resize(vnbr, 3);
    for(int i=0;i<vnbr;i++){
        int fid = fout[i];
        int vid0= F(fid,0);
        int vid1= F(fid,1);
        int vid2= F(fid,2);
        double u = bcout[i][0];
        double v = bcout[i][1];
        double t = bcout[i][2];
        Eigen::Vector3d ver0 = V.row(vid0);
        Eigen::Vector3d ver1 = V.row(vid1);
        Eigen::Vector3d ver2 = V.row(vid2);
        Vers.row(i) = ver0 * u + ver1 * v + ver2 * t;
    }
    std::cout<<"one curve drawn, the filtered size "<<Vers.rows()<<std::endl;
}

// void get_geodesic_dis_and_sort_pts(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, 
// const igl::HeatGeodesicsData<double> &data, const std::vector<int> &flist, std::vector<Eigen::Vector3f>& bclist)


// // this function takes the strokes as input, check if there are self-intersections, 
// bool lsTools::receive_interactive_strokes(const std::vector<int> &flist, const std::vector<Eigen::Vector3f>& bclist){
//     std::vector<int> tmp_flist;
//     std::vector<Eigen::Vector3f> tmp_bclist;
//     // check the distance between the first points of each curve. And get the sorted order of the curves
//     double t = std::pow(igl::avg_edge_length(V, F), 2);
//     igl::HeatGeodesicsData<double> data;
//     bool computed = igl::heat_geodesics_precompute(V,F,t,data);
//     if(!computed){
//         std::cout<<"Computing HeatGeodesic failed"<<std::endl;
//         return false;
//     }
//     int ncurves = flist.size();
//     std::vector<int> order;
//     for (int i = 0; i < ncurves; i++)
//     {
//         order.push_back(i);
//     }


//     Eigen::VectorXd D;
//     D = Eigen::VectorXd::Zero(V.rows());
//         for(int cid = 0;cid<3;cid++)
//         {
//           const int vid = F(fid,cid);
//           Eigen::VectorXd Dc;
//           igl::heat_geodesics_solve(data,(Eigen::VectorXi(1,1)<<vid).finished(),Dc);
//           D += Dc*bc(cid);


//     return false;
// }