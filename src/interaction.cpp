#include <lsc/interaction.h>
#include <lsc/tools.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/heat_geodesics.h>
// assign default arguments for interactive design.
void AssignAutoRunDefaultArgs(AutoRunArgs &args, const bool compute_pg)
{
    if (!compute_pg)
    { // just for a smooth field following the strokes
        args.compute_pg = false;
        args.stop_step_length = 1e-3; // when the step length lower than this, stop
        
        args.parts.resize(3);

        // first small boundary condition
        args.parts[0].iterations = 30;
        args.parts[0].weight_gravity = 1;
        args.parts[0].weight_lap = 10;
        args.parts[0].weight_bnd = 0.1;
        args.parts[0].weight_pg = 0;
        args.parts[0].weight_strip_width = 10;

        // second use large boundary condition
        args.parts[1] = args.parts[0];
        args.parts[1].weight_bnd = 1;

        // second use large boundary condition
        args.parts[2] = args.parts[1];
        args.parts[2].iterations = 100;
        args.parts[2].weight_bnd = 30;
    }
    else{//paras for pseodu-geodesic
        args.compute_pg = true;
        args.stop_step_length = 5e-2; // when the step length lower than this, stop
        args.stop_energy_sqrt = 1e-3;
        args.parts.resize(5);
        // 
        args.parts[0].iterations = 20;
        args.parts[0].weight_gravity = 100;
        args.parts[0].weight_lap = 0.001;
        args.parts[0].weight_bnd = 0.0001;
        args.parts[0].weight_pg = 0.001;
        args.parts[0].weight_strip_width = 0.0001;

        args.parts[1] = args.parts[0];
        args.parts[1].iterations = 20;
        args.parts[0].weight_lap = 0.0001;
        args.parts[1].weight_pg = 0.1;
        args.parts[1].weight_bnd = 0;
        args.parts[0].weight_strip_width = 0;

        args.parts[2] = args.parts[1];
        args.parts[2].iterations = 20;
        args.parts[2].weight_pg = 3;
        

        args.parts[3] = args.parts[2];
        args.parts[3].iterations = 100;
        args.parts[3].weight_pg = 30;

        // args.parts[4] = args.parts[3];
        // args.parts[4].iterations = 100;
        // args.parts[4].weight_pg = 10;
    }
}

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
#ifdef WIN32
		Sleep(time_interval);
#else
		sleep(time_interval);
#endif
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
void get_stroke_intersection_on_triangles(const CGMesh &mesh, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                                      const std::vector<int> &flist,
                                      const std::vector<Eigen::Vector3f> &bclist, std::vector<int> &fout, 
                                      std::vector<Eigen::Vector3f>& bcout)
{
    int fprev = -1;
    int fsize = flist.size();

    std::vector<Eigen::Vector3f> bfil1, bfil2;
    bfil1.reserve(fsize);
    bfil2.reserve(fsize);
    std::vector<std::vector<int>> fsort;
    std::vector<std::vector<Eigen::Vector3f>> bcsort;
    std::vector<int> ftmp; // record ids
    std::vector<Eigen::Vector3f> bctmp;
    // get the segments intersect with edges
    for (int i = 0; i < fsize; i++)
    {
        int fcurrent = flist[i];
        Eigen::Vector3f bccurrent(bclist[i][0], bclist[i][1], bclist[i][2]);
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
    Eigen::VectorXi flags = Eigen::VectorXi::Zero(F.rows());
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
        Eigen::Vector3f select_bc = bcsort[i][middle];
        if (flags[select_id] == 1)//this face is already appeared before
        {
            continue;
        }
        flags[select_id] = 1;
        fout.push_back(select_id);
        bcout.push_back(select_bc);
    }
}
void draw_stroke_on_mesh(const CGMesh &mesh, igl::opengl::glfw::Viewer &viewer,
                         const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const std::vector<double> &xs,
                         const std::vector<double> &ys, Eigen::MatrixXd& Vers, std::vector<int> &fout, 
                         std::vector<Eigen::Vector3f> &bcout)
{
    fout.clear();
    bcout.clear();
    Vers.resize(0,0);
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
    
    get_stroke_intersection_on_triangles(mesh, V, F,
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

//
bool get_furthest_end_pts_and_get_distance(const Eigen::MatrixXd &paras, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
                                           const std::vector<int> &flist, const std::vector<Eigen::Vector3f> &bclist,
                                           int &id0, int &id1, double &geodis)
{
    int ncandi = flist.size();
    assert(bclist.size()==ncandi);
    if (ncandi < 2)
    {
        std::cout<<"Too few Curves Drawn, Please retry with more curves"<<std::endl;
        return false;
    }
    std::vector<Eigen::Vector2d> pparas(ncandi);
    for (int i = 0; i < ncandi; i++)
    {
        int fid = flist[i];
        int vid0 = F(fid, 0);
        int vid1 = F(fid, 1);
        int vid2 = F(fid, 2);
        double u = bclist[i](0);
        double v = bclist[i](1);
        double t = bclist[i](2);
        pparas[i] = paras.row(vid0) * u + paras.row(vid1) * v + paras.row(vid2) * t;
    }
    double max_dis = 0;
    int tmpid0 = -1;
    int tmpid1 = -1;
    for (int i = 0; i < ncandi; i++)
    {
        // std::cout<<"pparas "<<i<<", "<<pparas[i].transpose()<<std::endl;
        for (int j = 0; j < ncandi; j++)
        {
            double dis = (pparas[i] - pparas[j]).norm();
            if (dis > max_dis)
            {
                max_dis = dis;
                tmpid0 = i;
                tmpid1 = j;
            }
        }
    }
    id0 = tmpid0;
    id1 = tmpid1;
    // std::cout<<"in here ids, "<<id0<<", "<<id1<<std::endl;
    double t = std::pow(igl::avg_edge_length(V, F), 2);
    igl::HeatGeodesicsData<double> data;
    bool computed = igl::heat_geodesics_precompute(V,F,t,data);
    if(!computed){
        std::cout<<"Computing HeatGeodesic failed"<<std::endl;
        return false;
    }
    // pick id0 as reference pt, and compute the geodesic distance between id0 and id1
    int fid = flist[id0];
    Eigen::Vector3f bc = bclist[id0];
    Eigen::VectorXd D;
    D = Eigen::VectorXd::Zero(V.rows());
    for (int cid = 0; cid < 3; cid++)
    {
        const int vid = F(fid, cid);
        Eigen::VectorXd Dc;
        igl::heat_geodesics_solve(data, (Eigen::VectorXi(1, 1) << vid).finished(), Dc);
        D += Dc * bc(cid);
    }
    int vid0 = F(fid, 0);
    int vid1 = F(fid, 1);
    int vid2 = F(fid, 2);
    geodis = D(vid0) * bc(0) + D(vid1) * bc(1) + D(vid2) * bc(2);
    return true;
}
void print_facelist(const std::vector<std::vector<int>> &flist)
{
    std::cout << "printing curve face list" << std::endl;
    for (int i = 0; i < flist.size(); i++)
    {
        for (int j = 0; j < flist[i].size(); j++)
        {
            std::cout << flist[i][j] << ", ";
        }
        std::cout << std::endl;
    }
}

// ep means end point
bool check_stroke_topology(const int fnbr, const std::vector<std::vector<int>> &flist,
                           const std::vector<std::vector<Eigen::Vector3f>> &bclist,
                           std::vector<int> &fep, std::vector<Eigen::Vector3f> &bep)
{
    int nc = flist.size();
    fep.resize(nc);
    bep.resize(nc);
    Eigen::VectorXi flags = Eigen::VectorXi::Zero(fnbr);
    Eigen::VectorXi ilist = Eigen::VectorXi::Ones(fnbr) * -1;
    Eigen::VectorXi jlist = Eigen::VectorXi::Ones(fnbr) * -1;
    for (int i = 0; i < flist.size(); i++)
    {
        fep[i] = flist[i][0];
        bep[i] = bclist[i][0];
        for (int j = 0; j < flist[i].size(); j++)
        {
            int fid = flist[i][j];
            if (flags[fid] == 1)
            {
                std::cout<<"Self intersected curve id: "<<i<<" pt id, "<<j<<
                "\nIt has intersection with "<<ilist[fid]<<", "<<jlist[fid]<<std::endl;
                print_facelist(flist);
                return false;
            }
            flags[fid] = 1;
            ilist[fid] = i;
            jlist[fid] = j;
        }
    }
    return true;
}

// this function takes the strokes as input, check if there are self-intersections, and use it as an initialization 
// of the levelset.
bool lsTools::receive_interactive_strokes_and_init_ls(const std::vector<std::vector<int>> &flist,
                                          const std::vector<std::vector<Eigen::Vector3f>> &bclist)
{
    int fnbr = F.rows();
    std::vector<int> fep;
    std::vector<Eigen::Vector3f> bep;
    bool topology = check_stroke_topology(fnbr, flist, bclist, fep, bep);
    if (!topology)
    {
        std::cout << "Curves self-intersected! Please re-draw the curves" << std::endl;
        return false;
    }
    int id0, id1;
    double geodis;
    topology = get_furthest_end_pts_and_get_distance(paras, V, F,
                                                     fep, bep, id0, id1, geodis);                                            
    if (!topology)
    {
        return false;
    }
    assert(id0 < fep.size() && id1 < fep.size());
    int fid0 = flist[id0][0];
    int fid1 = flist[id1][0];
    Eigen::Vector3f bc0 = bclist[id0][0];
    Eigen::Vector3f bc1 = bclist[id1][0];
    Eigen::Vector2d para_value0 = paras.row(F(fid0, 0)) * bc0(0) + paras.row(F(fid0, 1)) * bc0(1) + paras.row(F(fid0, 2)) * bc0(2);
    Eigen::Vector2d para_value1 = paras.row(F(fid1, 0)) * bc1(0) + paras.row(F(fid1, 1)) * bc1(1) + paras.row(F(fid1, 2)) * bc1(2);
    Eigen::Vector2d trans = (geodis) * (para_value1 - para_value0) / (para_value1 - para_value0).dot(para_value1 - para_value0);
    Eigen::VectorXd dupvalue = duplicate_valus(0, V.rows());
    Eigen::MatrixXd dupmatrix = duplicate_vector(para_value0, V.rows());
    Eigen::MatrixXd tmp = paras - dupmatrix; // matrix of Vi-V0
    Eigen::VectorXd result = tmp * trans;
    result = result - dupvalue;
    fvalues = result;
    interactive_flist = flist;
    interactive_bclist = bclist;
    Compute_Auxiliaries = true;// every time init the ls, we compute auxiliary variables
    return true;
}

double compute_stroke_avg_ls_value(const Eigen::VectorXd &func, const Eigen::MatrixXi& F, const std::vector<int> &flist,
                                   const std::vector<Eigen::Vector3f> &bclist)
{
    double result = 0;
    int npt = flist.size();
    for (int i = 0; i < npt; i++)
    {
        int fid = flist[i];
        double u = bclist[i][0];
        double v = bclist[i][1];
        double t = bclist[i][2];
        int v0 = F(fid, 0);
        int v1 = F(fid, 1);
        int v2 = F(fid, 2);

        result += func[v0] * u + func[v1] * v + func[v2] * t;
    }
    return result / npt;
}
void lsTools::assemble_solver_interactive_boundary_condition_part(const Eigen::VectorXd &func, spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &bcfvalue){
    spMat bcJacobian;
	int size = 0;// the number of constraints
	std::vector<double> vec_elements;
	std::vector<Trip> triplets;
    int ncurves = interactive_flist.size();
    assert(ncurves == interactive_bclist.size());
    std::vector<double> avg_vls(ncurves); // the average ls values
    for(int i = 0; i<ncurves; i++){
        avg_vls[i] = compute_stroke_avg_ls_value(func, F, interactive_flist[i], interactive_bclist[i]);
    }
	vec_elements.reserve(lsmesh.n_edges());
	triplets.reserve(lsmesh.n_edges());

	for (int i = 0; i < ncurves; i++)
	{
        double avg_func = avg_vls[i];
		for (int j = 0; j < interactive_flist[i].size(); j++) {

            int fid = interactive_flist[i][j];
            double u = interactive_bclist[i][j][0];
            double v = interactive_bclist[i][j][1];
            double t = interactive_bclist[i][j][2];
            int v0 = F(fid, 0);
            int v1 = F(fid, 1);
            int v2 = F(fid, 2);

            triplets.push_back(Trip(size, v0, u));
            triplets.push_back(Trip(size, v1, v));
            triplets.push_back(Trip(size, v2, t));
            double right_value = u * func[v0] + v * func[v1] + t * func[v2] - avg_func;
            vec_elements.push_back(right_value);
			size++;
		}
	}
	// get boundary condition: the jacobian matrix
	bcJacobian.resize(size, V.rows());
	bcJacobian.setFromTriplets(triplets.begin(), triplets.end());
	bcfvalue.resize(size);
	// std::cout<<"calculated a and b"<<std::endl;
	for (int i = 0; i < size; i++) {// get the f(x)
		double value = vec_elements[i];
		bcfvalue(i) = value;
	}

    // get JTJ
    H = bcJacobian.transpose() * bcJacobian;
    // get the -(J^T)*f(x)
    B = -bcJacobian.transpose() * bcfvalue;
}
template<typename Q>
void print_queue(Q q)
{
    // NB: q is passed by value because there is no way to traverse
    // priority_queue's content without erasing the queue.
    std::cout << '\n';
    for (; !q.empty(); q.pop())
        std::cout << q.top() << ", ";
    std::cout << '\n';
}


bool debug_flag = false;
double get_interactive_angle(const Eigen::VectorXd &func, const LSAnalizer &analizer, 
                             const CGMesh &lsmesh, const Eigen::MatrixXd &norm_v,
                             const Eigen::MatrixXd &V,
                             const Eigen::MatrixXi &F, const std::vector<std::vector<int>> interactive_flist,
                             const std::vector<std::vector<Eigen::Vector3f>> &interactive_bclist, const Eigen::VectorXi &InnerV)
{
    assert(func.size()>0);
    if(debug_flag){
        std::cout<<"check in 1, ";
    }
    int vnbr = V.rows();
    int ninner = analizer.LocalActInner.size();
    int ncurves = interactive_flist.size();
    double angle_avg = 0;
    int counter = 0;
    std::priority_queue<double> dbg_angles;
    // std::vector<bool> dbg_lesszero;
    if(debug_flag){
        std::cout<<"check in 2, ";
    }
    for (int i = 0; i < ncurves; i++)
    {
        int nbrpts = interactive_bclist[i].size();
        int start = 0.25 * nbrpts;
        int end = 0.75 * nbrpts;
        for (int j = start; j < end; j++)
        {
            assert(i < interactive_bclist.size() && j < interactive_bclist[i].size());
            assert(start <= end);
            double u = interactive_bclist[i][j][0];
            double v = interactive_bclist[i][j][1];
            double t = interactive_bclist[i][j][2];
            int k = -1;
            double closest_value = 1;
            int iterator;
            for (iterator = 0; iterator < 3; iterator++)
            {
                if (interactive_bclist[i][j][iterator] < closest_value)
                {
                    closest_value = interactive_bclist[i][j][iterator];
                    k = iterator;
                }
            }
            assert(k >= 0);
            int fid = interactive_flist[i][j];
            int vid[3];
            vid[0] = F(fid, 0);
            vid[1] = F(fid, 1);
            vid[2] = F(fid, 2);

            int vm = vid[k];
            int ith = InnerV[vm];
            if (ith < 0)// meaning the point is too close to the boundary, which can not be used to analyse the angle
            {
                continue;
            }
            assert(ith < analizer.LocalActInner.size() && ith >= 0);
            if (analizer.LocalActInner[ith] == false)
            {
                std::cout << "singularity, " << vm << std::endl;
                continue;
            }
            CGMesh::HalfedgeHandle inhd = analizer.heh0[ith], outhd = analizer.heh1[ith];
            int v1 = lsmesh.from_vertex_handle(inhd).idx();
            int v2 = lsmesh.to_vertex_handle(inhd).idx();
            int v3 = lsmesh.from_vertex_handle(outhd).idx();
            int v4 = lsmesh.to_vertex_handle(outhd).idx();

            double t1 = analizer.t1s[ith];
            double t2 = analizer.t2s[ith];

            Eigen::Vector3d ver0 = V.row(v1) + (V.row(v2) - V.row(v1)) * t1;
            Eigen::Vector3d ver1 = V.row(vm);
            Eigen::Vector3d ver2 = V.row(v3) + (V.row(v4) - V.row(v3)) * t2;
            Eigen::Vector3d norm = norm_v.row(vm);
            Eigen::Vector3d real_u = norm.cross(ver2 - ver0);
            real_u = real_u.normalized();
            Eigen::Vector3d real_r = (ver1 - ver0).normalized().cross((ver2 - ver1).normalized());
            if (real_r.norm() < 1e-6)
            { // almost a straight line
                continue;
            }
            real_r = real_r.normalized();
            double cos_angle = real_r.dot(norm);
            double sin_angle = real_r.dot(real_u);
            if (sin_angle < 0)
            {
                cos_angle *= -1;
            }
            double angle = acos(cos_angle);
            dbg_angles.push(180 / LSC_PI * angle);

            // converting the angle to [0~PI] is resonable
            // if sin_angle>0, just take the theta;
            // if sin_angle<0, take PI - theta
            if (sin_angle < 0)
            {
                // angle = LSC_PI - angle;
                // dbg_lesszero.push_back(1);
            }
            else
            {
                // dbg_lesszero.push_back(0);
            }
            angle_avg += angle;
            counter++;
        }
    }
    if(debug_flag){
        std::cout<<"check in 3, ";
    }
    // print_queue(dbg_angles);
    return angle_avg / counter;
}

void lsTools::Run_Level_Set_Opt_interactive(const bool compute_pg)
{

    Eigen::MatrixXd GradValueF, GradValueV;
    Eigen::VectorXd PGEnergy;
    // initial fvalues already computed using the strokes
    Eigen::VectorXd func = fvalues;

    int vnbr = V.rows();
    int fnbr = F.rows();
    if(func.size()==0){
        std::cout<<"Please init the level set using the interactive curves before solving optimizations"<<std::endl;
        return;
    }
    get_gradient_hessian_values(func, GradValueV, GradValueF);
    if(debug_flag){
        std::cout<<"check 1, ";
    }
    // std::cout<<"check "<<func.norm()<<std::endl;
    analysis_pseudo_geodesic_on_vertices(func, analizers[0]);
    int ninner = analizers[0].LocalActInner.size();
    int final_size;
    if (!compute_pg) // only the boundary condition
    {
        final_size = vnbr;
    }
    else
    {
        final_size = ninner * 10 + vnbr; // pseudo - geodesic
    }
    

	if (Glob_lsvars.size() != final_size) {
        Compute_Auxiliaries = true;
		if (Glob_lsvars.size() == 0)
		{
			std::cout << "Initializing Global Variable For LevelSet Opt ... " << std::endl;
			Glob_lsvars = Eigen::VectorXd::Zero(final_size); // We change the size if opt more than 1 level set
			Glob_lsvars.segment(0, vnbr) = func;
            std::cout << "Global Variable For LevelSet is initialized " << std::endl;
		}
		else{
			if (Glob_lsvars.size() > final_size)
			{ // need to truncate
				Eigen::VectorXd tmp = Glob_lsvars.segment(0, final_size);
				Glob_lsvars = tmp;
				if (Glob_lsvars.size() != final_size)
				{
					std::cout << "size error" << std::endl
							  << std::endl;
				}
			}
			else// need to extend, and re-compute the auxiliaries
			{
				Eigen::VectorXd tmp = Glob_lsvars;
				Glob_lsvars = Eigen::VectorXd::Zero(final_size); // We change the size if opt more than 1 level set
				Glob_lsvars.segment(0, tmp.size()) = tmp;
				
			}
		}
	}
    if(debug_flag){
        std::cout<<"check 2, ";
    }
	// std::cout<<"check 1"<<std::endl;
	spMat H;
	H.resize(vnbr, vnbr);
    Eigen::VectorXd B = Eigen::VectorXd::Zero(vnbr);

    spMat LTL;  // left of laplacian
	Eigen::VectorXd mLTF; // right of laplacian
	assemble_solver_biharmonic_smoothing(func, LTL, mLTF);
	H += weight_laplacian * LTL;
	B += weight_laplacian * mLTF;
	assert(mass.rows() == vnbr);
    // std::cout<<"check 2"<<std::endl;
	// fix inner vers and smooth boundary
	Eigen::VectorXd bcfvalue = Eigen::VectorXd::Zero(vnbr);
	// boundary condition
    if (interactive_flist.size() > 0)
    { // if traced, we use boundary condition

        spMat bc_JTJ;
        Eigen::VectorXd bc_mJTF;
        assemble_solver_interactive_boundary_condition_part(func, bc_JTJ, bc_mJTF, bcfvalue);
        H += weight_boundary * bc_JTJ;
        B += weight_boundary * bc_mJTF;
    }
    // std::cout<<"check 3"<<std::endl;
    // strip width condition
    if(debug_flag){
        std::cout<<"check 3, ";
    }
    spMat sw_JTJ;
    Eigen::VectorXd sw_mJTF;
    assemble_solver_strip_width_part(GradValueF, sw_JTJ, sw_mJTF); // by default the strip width is 1. Unless tracing info updated the info
    H += weight_strip_width * sw_JTJ;
    B += weight_strip_width * sw_mJTF;
    if(debug_flag){
        std::cout<<"check 4, ";
    }
    assert(H.rows() == vnbr);
	assert(H.cols() == vnbr);
	// pseudo geodesic
	spMat Hlarge;
	Eigen::VectorXd Blarge;

	Eigen::VectorXd e_smbi;
	Hlarge.resize(final_size, final_size);
	Blarge = Eigen::VectorXd::Zero(final_size);
	
	Hlarge = sum_uneven_spMats(H, Hlarge);
	Blarge = sum_uneven_vectors(B, Blarge);
    // std::cout<<"check 4"<<std::endl;
	if (compute_pg)
	{

		spMat pg_JTJ;
		Eigen::VectorXd pg_mJTF;
		int vars_start_loc = 0;
		int aux_start_loc = vnbr;
        std::vector<double> angle_degree(1);
        angle_degree[0] = pseudo_geodesic_target_angle_degree;
        if(Compute_Auxiliaries){
            std::cout<<"Computing Auxiliaries"<<std::endl;
        }
        assemble_solver_pesudo_geodesic_energy_part_vertex_based(Glob_lsvars, angle_degree, analizers[0],
                                                                 vars_start_loc, aux_start_loc, pg_JTJ, pg_mJTF, PGEnergy);

		Hlarge = sum_uneven_spMats(Hlarge, weight_pseudo_geodesic_energy * pg_JTJ);
		Blarge = sum_uneven_vectors(Blarge, weight_pseudo_geodesic_energy * pg_mJTF);
        // std::cout<<"weight_pseudo_geodesic_energy, "<<weight_pseudo_geodesic_energy<<std::endl;
		Compute_Auxiliaries = false;
	}
    if(debug_flag){
        std::cout<<"check 5, ";
    }
	if(vector_contains_NAN(Blarge)){
        std::cout<<"energy contains NAN"<<std::endl;
		return;
    }
    // std::cout<<"check 5"<<std::endl;
	Hlarge += 1e-6 * weight_mass * spMat(Eigen::VectorXd::Ones(final_size).asDiagonal());

	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver(Hlarge);

	// assert(solver.info() == Eigen::Success);
	if (solver.info() != Eigen::Success)
	{
		// solving failed
		std::cout << "solver fail" << std::endl;
		return;
	}
	// std::cout<<"solved successfully"<<std::endl;
	Eigen::VectorXd dx = solver.solve(Blarge).eval();
	dx *= 0.75;
	// std::cout << "step length " << dx.norm() << std::endl;
	double level_set_step_length = dx.norm();


	func += dx.topRows(vnbr);
	fvalues = func;
	Glob_lsvars += dx;
	double energy_biharmonic = (QcH * func).norm();
	double energy_boundary = (bcfvalue).norm();
	std::cout << "energy: harm " << energy_biharmonic << ", bnd " << energy_boundary << ", ";
    // std::cout<<"check 5.5"<<std::endl;
    if (compute_pg)
    {
        double energy_pg = PGEnergy.norm();
        pgerror = energy_pg * energy_pg;
        double max_energy_ls = PGEnergy.lpNorm<Eigen::Infinity>();
        std::cout << "pg, " << energy_pg << ", "
                  << "lsmax," << max_energy_ls << ",";
    }
    // std::cout<<"check 6"<<std::endl;
    Eigen::VectorXd ener = GradValueF.rowwise().norm();
    ener = ener.asDiagonal() * ener;
    Eigen::VectorXd wds = Eigen::VectorXd::Ones(fnbr) * strip_width * strip_width;

    ener -= wds;
    double stp_energy = Eigen::VectorXd(ener).dot(ener);
    std::cout << "strip, " << stp_energy << ", ";

    step_length = dx.norm();
	std::cout << "step " << step_length;
    if(debug_flag){
        std::cout<<"check 6, ";
    }
    // get the average angle around the strokes
    if (!compute_pg)
    {
        pseudo_geodesic_target_angle_degree = 180. / LSC_PI * get_interactive_angle(func, analizers[0], lsmesh, norm_v,
         V, F, interactive_flist, interactive_bclist, InnerV);
    }
    if(debug_flag){
        std::cout<<"check 7, ";
    }
    std::cout<<", target_angle, "<<pseudo_geodesic_target_angle_degree;
    std::cout<<"\n";
    Last_Opt_Mesh = false;
}
