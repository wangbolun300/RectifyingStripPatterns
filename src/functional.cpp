// This file deals with functional target angle values
#include <lsc/basic.h>
#include <lsc/tools.h>
#include <igl/exact_geodesic.h>
void assign_angles_based_on_funtion_values(const Eigen::VectorXd &fvalues, const double angle_large_function, const double angle_small_function,
                                           std::vector<double>& angle_degree)
{
    int vnbr=fvalues.size();
    angle_degree.resize(vnbr);
    
    double fmax=fvalues.maxCoeff();
    double fmin=fvalues.minCoeff();
    double favg=(fmax+fmin)/2;
    double itv=(angle_large_function-angle_small_function) /(fmax-fmin);
    for(int i=0;i<vnbr;i++){
        double value=fvalues[i];
        double ang=angle_small_function+(value-fmin)*itv;
        angle_degree[i]=ang;
    }

}

// IVids is a vector of size ninner, mapping inner ver ids to ver ids.
void get_geodesic_distance(const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const std::vector<int> &IVids,
                           const std::vector<int> &idinner, Eigen::VectorXd &D)
{
    
    Eigen::VectorXi VS, FS, VT, FT;
    VS.resize(idinner.size());
    for (int i = 0; i < idinner.size(); i++)
    {
        int id = idinner[i]; // the id in inner vers
        int vm = IVids[id];
        VS[i] = vm;
    }
    // All vertices are the targets
    VT.setLinSpaced(V.rows(), 0, V.rows() - 1);
    igl::exact_geodesic(V, F, VS, FS, VT, FT, D);
}

// the strategy is: within the percentage % of the distance, the angle increase linearly.
void assign_pg_angle_based_on_geodesic_distance(const Eigen::VectorXd &D, const double percentage, Eigen::VectorXd &angle)
{
    double dismax = D.maxCoeff();
    double dthreads = dismax * percentage / 100.;
    angle.resize(D.size());
    // std::cout<<"dis max, "<<dismax<<", threadshold "<<dthreads<<", percentage, "<<percentage<<std::endl;
    for (int i = 0; i < angle.size(); i++)
    {
        double dis = D[i];
        double a;
        if (dis < dthreads)
        {
            a = dis / dthreads * 90;
        }
        else
        {
            a = 90;
        }
        angle[i] = a;
    }
}

// void get_lengths_of_plys()
// idinner is the inner vers of asymptotic area.
// IVids[idinner] is the vid of the inner vers
// void assign_pg_angle_based_on_ls_length(const CGMesh &lsmesh, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
//                                         const Eigen::VectorXd &ls,
//                                         const std::vector<CGMesh::HalfedgeHandle> &loop,
//                                         const std::vector<int> &IVids, const std::vector<int> &idinner,
//                                         const double percentage,  
//                                         std::vector<double>& ts)
// {
//     vidangle.clear();
//     vidangle.reserve(V.rows() / 2);
//     double lsmin = ls.minCoeff();
//     double lsmax = ls.maxCoeff();
//     double value = (lsmin + lsmax) / 2;
//     std::vector<std::vector<double>> paras;
//     std::vector<std::vector<CGMesh::HalfedgeHandle>> handles;
//     std::vector<bool> left_large;
//     std::vector<Eigen::Vector3d> pts;
//     get_iso_lines_baricenter_coord(lsmesh, loop, V, ls, value, paras, handles, left_large);
    
//     if()
//     int maxsize = 0;
//     int maxid
    

// }

// angles are in radian, size is ninner, target are in degree, size is vnbr, raw_angles size is ninner
// increase the angle, and fit it into a polynomial to get the target angles.

void assign_angle_get_closer_angle_to_asymptotic(const double percent, const Eigen::VectorXd &angles, const int vnbr,
                                    const Eigen::VectorXd &fvalues,
                                    const std::vector<int> &IVids, std::vector<double> &target, Eigen::VectorXd& raw_angles)
{
    int ninner = angles.size();
    double aasym = 0; // target angle is 0 or 180. 180 is included in 0.
    Eigen::VectorXd Xs(ninner), Ys(ninner);
    Eigen::MatrixXd A(angles.size(), 4);
    raw_angles.resize(angles.size());
    for (int i = 0; i < angles.size(); i++)
    {
        double a = angles[i];
        a += (aasym - a) * percent / 100.;
        raw_angles[i] = a;

        int vm = IVids[i];
        double x = fvalues[vm];
        Xs[i] = x;
        Ys[i] = a;
        A.row(i) << x * x * x, x * x, x, 1;
        // double angle_degree = a * 180 / LSC_PI;
        // target[vm] = angle_degree;
    }
    Eigen::Matrix4d ATA = A.transpose() * A, inv;
    bool invertable;
    ATA.computeInverseWithCheck(inv, invertable);
    if(!invertable){
        std::cout<<"MATRIX NOT INVERTABLE IN void assign_angle_get_closer_angle_to_asymptotic()"<<std::endl;
        return;
    }
    Eigen::Vector4d paras = inv*A.transpose()*Ys;
    // output the fitting results
    target.resize(vnbr);
    double a = paras[0];
    double b = paras[1];
    double c = paras[2];
    double d = paras[3];
    for (int i = 0; i < ninner; i++)
    {
        // i is the id in inner, vm is the id in global
        int vm = IVids[i];
        double f = fvalues[vm];
        double ang = a * f * f * f + b * f * f + c * f + d; // this is in radian
        ang = ang * 180. / LSC_PI; // this is in degree.
        if (ang > 90)// don't exceed 90 degree.
        {
            ang = 90;
        }
        if (ang < 0)
        {
            ang = 0;
        }
        target[vm] = ang;
    }
    std::cout<<"TargetFitting, "<<a<<", "<<b<<", "<<c<<", "<<d<<", ";
}

//
void assign_angle_by_poly_fitting(const Eigen::VectorXd &func, const std::vector<double> &angle_degree,
                                  const Eigen::VectorXi &InnerV, const int ninner,
                                  const double perc_conservative, std::vector<double> &angle_out)
{
    int vnbr = func.size();
    Eigen::VectorXd Xs(ninner), Ys(ninner);
    for (int i = 0; i < vnbr; i++) // first make the angle conservative.
    {
        int loc = InnerV[i];// location in inner vers
        if (loc < 0)
        {
            continue;
        }
        double ang_ref = angle_degree[i];
        double ang_conserv = ang_ref * (1 + perc_conservative);
        if (ang_conserv > 90)
        {
            ang_conserv = ang_ref;
        }
        Xs[loc] = func[i];
        Ys[loc] = ang_conserv;
    }
    // fit it to a polynomial of degree 3
    Eigen::MatrixXd A(ninner, 4);
    for(int i=0;i<ninner;i++){
        double x = Xs[i];
        A(i, 0) = x * x * x;
        A(i, 1) = x * x;
        A(i, 2) = x;
        A(i, 3) = 1;
    }
    Eigen::Matrix4d ATA = A.transpose() * A, inv;
    bool invertable;
    ATA.computeInverseWithCheck(inv, invertable);
    if(!invertable){
        std::cout<<"MATRIX NOT INVERTABLE IN void assign_angle_by_poly_fitting()"<<std::endl;
        return;
    }
    Eigen::Vector4d paras = inv*A.transpose()*Ys;
    // output the fitting results
    angle_out.resize(vnbr);
    double a = paras[0];
    double b = paras[1];
    double c = paras[2];
    double d = paras[3];
    for (int i = 0; i < vnbr; i++)
    {
        int loc = InnerV[i];// location in inner vers
        if (loc < 0)
        {
            continue;
        }
        double f = func[i];
        
        // angle = a*f^3 + b*f^2 + c*f + d
        double ang = a * f * f * f + b * f * f + c * f + d;
        if (ang > 90)
        {
            ang = 90;
        }
        if (ang < 0)
        {
            ang = 0;
        }
        angle_out[i] = ang;
    }
    std::cout<<"Fitting, "<<a<<", "<<b<<", "<<c<<", "<<d<<", ";
}
void lsTools::assemble_solver_follow_min_abs_curvature(spMat &H, Eigen::VectorXd &B, Eigen::VectorXd& Eng){
    if(CurvatureDirections.empty()){
        get_orthogonal_direction_minimal_principle_curvature(V, F, PosFids, NegFids, OrthMinK, CosSqr, CurvatureDirections);
    }
    int fnbr = F.rows();
    std::vector<Trip> tripletes;
    Eigen::VectorXd energy;
    tripletes.reserve(3 * (PosFids.size() + NegFids.size()));
    energy.resize(PosFids.size() + NegFids.size());
    for (int itr = 0; itr < PosFids.size(); itr++) // positive areas, orthogonal to the direction
	{
        
        int fid = PosFids[itr];
		Eigen::Vector3d otho = OrthMinK[itr];
		int v0 = F(fid, 0);
		int v1 = F(fid, 1);
		int v2 = F(fid, 2);
		// locations
		int loc0 = v0;
		int loc1 = v1;
		int loc2 = v2;
		
		double f0 = Glob_lsvars[loc0];
		double f1 = Glob_lsvars[loc1];
		double f2 = Glob_lsvars[loc2];
		Eigen::Vector3d ver0 = V.row(v0);
		Eigen::Vector3d ver1 = V.row(v1);
		Eigen::Vector3d ver2 = V.row(v2);

		Eigen::Vector3d direction = f0 * (ver1 - ver2) + f2 * (ver0 - ver1) + f1 * (ver2 - ver0);
		double scale = direction.norm();
		// direction * otho / scale = 0;
		tripletes.push_back(Trip(itr, loc0, (ver1 - ver2).dot(otho) / scale));
		tripletes.push_back(Trip(itr, loc1, (ver2 - ver0).dot(otho) / scale));
		tripletes.push_back(Trip(itr, loc2, (ver0 - ver1).dot(otho) / scale));

		energy[itr] = direction.dot(otho) / scale;

    }
    int condistart = PosFids.size();
    for (int itr = 0; itr < NegFids.size(); itr++) // negative areas, the angle between the tangent direction and the 
	{
        
        int fid = NegFids[itr];
		Eigen::Vector3d reference = CurvatureDirections[fid][0].normalized();
        double cc = CosSqr[itr];

		int v0 = F(fid, 0);
		int v1 = F(fid, 1);
		int v2 = F(fid, 2);
		// locations
		int loc0 = v0;
		int loc1 = v1;
		int loc2 = v2;
		
		double f0 = Glob_lsvars[loc0];
		double f1 = Glob_lsvars[loc1];
		double f2 = Glob_lsvars[loc2];
		Eigen::Vector3d ver0 = V.row(v0);
		Eigen::Vector3d ver1 = V.row(v1);
		Eigen::Vector3d ver2 = V.row(v2);

		Eigen::Vector3d direction = f0 * (ver1 - ver2) + f2 * (ver0 - ver1) + f1 * (ver2 - ver0);
		double scale = direction.norm();
        double dot = direction.dot(reference);
        double v12r = (ver1 - ver2).dot(reference);
        double v20r = (ver2 - ver0).dot(reference);
        double v01r = (ver0 - ver1).dot(reference);
        // ((direction / scale) * reference)^2 = cc
        tripletes.push_back(Trip(itr + condistart, loc0, 2 * v12r * dot));
        tripletes.push_back(Trip(itr + condistart, loc1, 2 * v20r * dot));
        tripletes.push_back(Trip(itr + condistart, loc2, 2 * v01r * dot));

        energy[itr + condistart] = dot * dot - scale * scale * cc;
    }
    spMat J;
    J.resize(energy.size(), Glob_lsvars.size());
    J.setFromTriplets(tripletes.begin(), tripletes.end());
    H = J.transpose() * J;
    B = -J.transpose() * energy;
    Eng = energy;
}
void  lsTools::Run_AsOrthAsPossible_LS(){
    Analyze_Optimized_LS_Angles = true;
    Eigen::MatrixXd GradValueF, GradValueV;
	Eigen::VectorXd PGEnergy;
	Eigen::VectorXd func = fvalues;
	
	std::vector<double> angle_degree;

    int vnbr = V.rows();
	int fnbr = F.rows();
    bool early_terminate = false;

    
    // initialize the level set with some number
	if (func.size() != vnbr)
	{
		if(trace_hehs.size()==0){
			std::cout<<"Please First Set Up The Boundary Condition (Tracing)"<<std::endl;
			return;
		}
		initialize_level_set_accroding_to_parametrization();
		std::cout << "level set get initialized" << std::endl;
		return;
	}
	get_gradient_hessian_values(func, GradValueV, GradValueF);

	analysis_pseudo_geodesic_on_vertices(func, analizers[0]);
	int ninner = analizers[0].LocalActInner.size();
    // the logic: 1. For the first iteration, we compute the real angles, and terminate.
    // 2. Then, we optimize with increasing angles. Sometimes we just stop increasing, 
    // but use the fixed target angles to optimize until convergence.
    // 3. in the last phase, we optimize with the conservatively fitted angles.
    if (AnalizedAngelVector.size() == 0)
    {
        early_terminate = true;
        angle_degree.push_back(0); // put any angle as target, since we terminate early anyway.
    }
    else
    {
        if(Disable_Changed_Angles){ // use the unchanged angels until converges.
            angle_degree = changed_angles;
        }
        else{
            double percentage = 5;
            assign_angle_get_closer_angle_to_asymptotic(percentage, AnalizedAngelVector, vnbr, fvalues, IVids, angle_degree, RawTargetAngles);
            changed_angles = angle_degree;
        }

        if (Use_Fitting_Angles) // use the fitting result.
        {
            double percentage = 10;
            assign_angle_by_poly_fitting(func, changed_angles, InnerV, ninner, percentage, angle_degree);
        }
    }

	int final_size;
    final_size = ninner * 10 + vnbr;// pseudo - geodesic 
	
	
	
	if (Glob_lsvars.size() != final_size) {
		
		if (Glob_lsvars.size() == 0)
		{
			std::cout << "Initializing Global Variable For LevelSet Opt ... " << std::endl;
			Glob_lsvars = Eigen::VectorXd::Zero(final_size); // We change the size if opt more than 1 level set
			Glob_lsvars.segment(0, vnbr) = func;			 
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
				Compute_Auxiliaries = true;
			}
		}
	}
	if(Last_Opt_Mesh){
		Compute_Auxiliaries = true;
		std::cout<<"Recomputing Auxiliaries"<<std::endl;
	}
	
	spMat H;
	H.resize(vnbr, vnbr);
    Eigen::VectorXd B = Eigen::VectorXd::Zero(vnbr);

    spMat LTL;  // left of laplacian
	Eigen::VectorXd mLTF; // right of laplacian
	assemble_solver_biharmonic_smoothing(func, LTL, mLTF);
	H += weight_laplacian * LTL;
	B += weight_laplacian * mLTF;
	assert(mass.rows() == vnbr);

	// fix inner vers and smooth boundary
	if (enable_inner_vers_fixed)
	{
		// H += weight_mass * InnerV.asDiagonal();
		std::cout << "Please don't use fixing_inner_vertex energy" << std::endl;
	}
	Eigen::VectorXd bcfvalue = Eigen::VectorXd::Zero(vnbr);
	
	if (interactive_flist.size() > 0)
    { // if traced, we use boundary condition

        spMat bc_JTJ;
        Eigen::VectorXd bc_mJTF;
        assemble_solver_interactive_boundary_condition_part(func, bc_JTJ, bc_mJTF, bcfvalue);
        H += weight_boundary * bc_JTJ;
        B += weight_boundary * bc_mJTF;
    }

	// strip width condition
	if (enable_strip_width_energy)
	{
		spMat sw_JTJ;
		Eigen::VectorXd sw_mJTF;
		assemble_solver_strip_width_part(GradValueF, sw_JTJ, sw_mJTF);// by default the strip width is 1. Unless tracing info updated the info
		H += weight_strip_width * sw_JTJ;
		B += weight_strip_width * sw_mJTF;
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
    
	if (enable_pseudo_geodesic_energy)
	{

		spMat pg_JTJ;
		Eigen::VectorXd pg_mJTF;
		spMat smbi_H;
		Eigen::VectorXd smbi_B;
		int vars_start_loc = 0;
		int aux_start_loc = vnbr;

        assemble_solver_pesudo_geodesic_energy_part_vertex_based(Glob_lsvars, angle_degree, analizers[0],
                                                                 vars_start_loc, aux_start_loc, pg_JTJ, pg_mJTF, PGEnergy);

		Hlarge = sum_uneven_spMats(Hlarge, weight_pseudo_geodesic_energy * pg_JTJ);
		Blarge = sum_uneven_vectors(Blarge, weight_pseudo_geodesic_energy * pg_mJTF);

		Compute_Auxiliaries = false;
	}
    if(early_terminate){
        std::cout<<"Early Terminate for the first iteration"<<std::endl;
        return;
    }
	if(vector_contains_NAN(Blarge)){
        std::cout<<"energy contains NAN"<<std::endl;
		return;
    }
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
	// double inf_norm=dx.cwiseAbs().maxCoeff();
	if (enable_boundary_angles || enable_pseudo_geodesic_energy)// only pg energy and boundary angles requires step length
	{
		if (level_set_step_length > max_step_length)
		{
			dx *= max_step_length / level_set_step_length;
		}
	}

	func += dx.topRows(vnbr);
	fvalues = func;
	Glob_lsvars += dx;
	double energy_biharmonic = (QcH * func).norm();
	double energy_boundary = (bcfvalue).norm();
    std::cout << "energy: harm " << energy_biharmonic << ", bnd " << energy_boundary << ", ";
    if (enable_pseudo_geodesic_energy)
    {
        double energy_pg = PGEnergy.norm();
        double max_energy_ls = PGEnergy.lpNorm<Eigen::Infinity>();
        std::cout << "pg, " << energy_pg << ", "
                  << "lsmax," << max_energy_ls << ",";
    }

    if (enable_strip_width_energy)
	{

		Eigen::VectorXd ener = GradValueF.rowwise().norm();
		ener = ener.asDiagonal() * ener;
		Eigen::VectorXd wds = Eigen::VectorXd::Ones(fnbr) * strip_width * strip_width;

		ener -= wds;
		double stp_energy = Eigen::VectorXd(ener).dot(ener);
		std::cout << "strip, " << stp_energy << ", ";
	}
	step_length = dx.norm();
	std::cout << "step " << step_length;
	if (interactive_flist.size() > 0) // if start to solve pseudo-geodeic, we refer to the strokes
	{
		pseudo_geodesic_target_angle_degree = 180. / LSC_PI * get_interactive_angle(func, analizers[0], lsmesh, norm_v, V, F, interactive_flist, interactive_bclist, InnerV);
		std::cout << ", target_angle, " << pseudo_geodesic_target_angle_degree;
	}
    assert(AnalizedAngelVector.size() > 0);
    std::cout<<", current angle: max, "<<AnalizedAngelVector.maxCoeff()<<", min, "<<AnalizedAngelVector.minCoeff()
    <<", ang_average, "<<AnalizedAngelVector.mean()<<", ";

	std::cout << std::endl;
	Last_Opt_Mesh = false;
}
void lsTools::RunInitAOAP(){
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

    // std::cout<<"check "<<func.norm()<<std::endl;
    analysis_pseudo_geodesic_on_vertices(func, analizers[0]);
    int ninner = analizers[0].LocalActInner.size();
    int final_size;
    final_size = vnbr;
    

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

	spMat H;
	H.resize(vnbr, vnbr);
    Eigen::VectorXd B = Eigen::VectorXd::Zero(vnbr);

    spMat LTL;  // left of laplacian
	Eigen::VectorXd mLTF; // right of laplacian
	assemble_solver_biharmonic_smoothing(func, LTL, mLTF);
	H += weight_laplacian * LTL;
	B += weight_laplacian * mLTF;
	assert(mass.rows() == vnbr);


    // strip width condition
    spMat sw_JTJ;
    Eigen::VectorXd sw_mJTF;
    assemble_solver_strip_width_part(GradValueF, sw_JTJ, sw_mJTF); // by default the strip width is 1. Unless tracing info updated the info
    H += weight_strip_width * sw_JTJ;
    B += weight_strip_width * sw_mJTF;

    assert(H.rows() == vnbr);
	assert(H.cols() == vnbr);
	// pseudo geodesic
	spMat Hlarge;
	Eigen::VectorXd Blarge;
	Hlarge.resize(final_size, final_size);
	Blarge = Eigen::VectorXd::Zero(final_size);
	
	Hlarge = sum_uneven_spMats(H, Hlarge);
    Blarge = sum_uneven_vectors(B, Blarge);
    // std::cout<<"check 4"<<std::endl;

    spMat pg_JTJ;
    Eigen::VectorXd pg_mJTF;
    
    assemble_solver_follow_min_abs_curvature(pg_JTJ, pg_mJTF, PGEnergy);

    Hlarge = sum_uneven_spMats(Hlarge, weight_pseudo_geodesic_energy * pg_JTJ);
    Blarge = sum_uneven_vectors(Blarge, weight_pseudo_geodesic_energy * pg_mJTF);
    // std::cout<<"weight_pseudo_geodesic_energy, "<<weight_pseudo_geodesic_energy<<std::endl;
    Compute_Auxiliaries = false;

    if (vector_contains_NAN(Blarge))
    {
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
	// double inf_norm=dx.cwiseAbs().maxCoeff();
	// if (enable_boundary_angles || enable_pseudo_geodesic_energy)// only pg energy and boundary angles requires step length
	// {
	// 	if (level_set_step_length > max_step_length)
	// 	{
	// 		dx *= max_step_length / level_set_step_length;
	// 	}
	// }

	func += dx.topRows(vnbr);
	fvalues = func;
	Glob_lsvars += dx;
    double energy_biharmonic = (QcH * func).norm();
    std::cout << "INIT AOAP: energy: harm " << energy_biharmonic << ", ";
    // std::cout<<"check 5.5"<<std::endl;

    double energy_pg = PGEnergy.norm();
    double max_energy_ls = PGEnergy.lpNorm<Eigen::Infinity>();
    std::cout << "pg, " << energy_pg << ", "
              << "lsmax," << max_energy_ls << ",";

    // std::cout<<"check 6"<<std::endl;
    Eigen::VectorXd ener = GradValueF.rowwise().norm();
    ener = ener.asDiagonal() * ener;
    Eigen::VectorXd wds = Eigen::VectorXd::Ones(fnbr) * strip_width * strip_width;

    ener -= wds;
    double stp_energy = Eigen::VectorXd(ener).dot(ener);
    std::cout << "strip, " << stp_energy << ", ";

    step_length = dx.norm();
	std::cout << "step " << step_length;
    std::cout<<"\n";
    Last_Opt_Mesh = false;
}
#include<igl/file_dialog_save.h>
void lsTools::write_fitting_data(){
    // write: 1. real angles: AnalizedAngelVector (size is ninner), 2. the target angles: changed_angles (size is vnbr), 
    // 3. fvalues: size is vnbr. 4. raw target angles RawTargetAngles: (size is ninner)
    // we only write the inner ver values
    int ninner = analizers[0].LocalActInner.size();
    std::string fname = igl::file_dialog_save();
    std::ofstream file;
    file.open(fname);
    for (int i = 0; i < ninner; i++)
    {
        int vloc = IVids[i];
        double realangle = AnalizedAngelVector[i];
        double targangle = changed_angles[vloc] * LSC_PI / 180.;
        double realfunc = fvalues[vloc];
        double tar_raw = RawTargetAngles[i];
        file<<realangle<<","<<targangle<<","<<realfunc<<","<<tar_raw<<"\n";

    }
    file.close();
    std::cout<<"Necessary info saved"<<std::endl;
}

// the auxiliaries:, the axis, the target cos value, the target sin value ,
// 0 ~ vnbr-1 : function values, 
// vnbr, vnbr + 1, vnbr + 2: the axis
// vnbr + 3, cos
// vnbr + 4, sin

void lsTools::assemble_solver_constant_slope(const Eigen::VectorXd &vars,
                                             const bool fix_axis, const Eigen::Vector3d &axis_in, const bool fix_angle, const double angle_degree,
                                             const LSAnalizer &analizer,
                                             spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy)
{
    std::vector<Trip> tripletes;
    int vnbr = V.rows();
    int fnbr = F.rows();
    double cos_angle_fix, sin_angle_fix;
    double angle_radian = angle_degree * LSC_PI / 180.; // the angle in radian
    cos_angle_fix = cos(angle_radian);
    sin_angle_fix = sin(angle_radian);
    Eigen::Vector3d axis_fix = axis_in.normalized();
    tripletes.reserve(10 + fnbr * 11);              //
    energy = Eigen::VectorXd::Zero(2 * fnbr + 5); // mesh total energy values

    int lax = vnbr;
    int lay = vnbr + 1;
    int laz = vnbr + 2;
    int lcos = vnbr + 3;
    int lsin = vnbr + 4;
    if (Compute_Auxiliaries)
    {
        if (fix_axis)
        {
            Glob_lsvars[lax] = axis_fix[0];
            Glob_lsvars[lay] = axis_fix[1];
            Glob_lsvars[laz] = axis_fix[2];
        }
        else
        {
            Glob_lsvars[lax] = 0;
            Glob_lsvars[lay] = 0;
            Glob_lsvars[laz] = 1;
        }
        if (fix_angle)
        {
            Glob_lsvars[lcos] = cos_angle_fix;
            Glob_lsvars[lsin] = sin_angle_fix;
        }
        else
        {
            Glob_lsvars[lcos] = 1 / 1.414;
            Glob_lsvars[lsin] = 1 / 1.414;
        }
    }
    double cos_angle = Glob_lsvars[lcos];
    double sin_angle = Glob_lsvars[lsin];
    Eigen::Vector3d axis(Glob_lsvars[lax], Glob_lsvars[lay], Glob_lsvars[laz]);
    for (int i = 0; i < fnbr; i++)
    {
        int fid = i;
        int v0 = F(fid, 0);
        int v1 = F(fid, 1);
        int v2 = F(fid, 2);
        // dir0 * V0 + dir1 * V1 + dir2 * V2 is the iso-line direction
        Eigen::Vector3d dir0 = V.row(v2) - V.row(v1);
        Eigen::Vector3d dir1 = V.row(v0) - V.row(v2);
        Eigen::Vector3d dir2 = V.row(v1) - V.row(v0);
        Eigen::Vector3d iso = dir0 * vars[v0] + dir1 * vars[v1] + dir2 * vars[v2];
        Eigen::Vector3d norm = norm_f.row(fid);
        
        Eigen::Vector3d bxis;
        if (axis.cross(norm).norm() > 1e-4)
        {
            bxis = axis.cross(norm).normalized();
        }
        else
        {
            continue;
        }
        double scale = iso.norm() * axis.norm();
        // iso * axis / scale - cos = 0;
        tripletes.push_back(Trip(i, v0, dir0.dot(axis) / scale));
        tripletes.push_back(Trip(i, v1, dir1.dot(axis) / scale));
        tripletes.push_back(Trip(i, v2, dir2.dot(axis) / scale));

        tripletes.push_back(Trip(i, lax, iso[0] / scale));
        tripletes.push_back(Trip(i, lay, iso[1] / scale));
        tripletes.push_back(Trip(i, laz, iso[2] / scale));

        tripletes.push_back(Trip(i, lcos, -1));

        energy[i] = iso.dot(axis) / scale - cos_angle;

        // iso * bxis / scale - sin = 0
        scale = iso.norm();
        tripletes.push_back(Trip(i + fnbr, v0, dir0.dot(bxis) / scale));
        tripletes.push_back(Trip(i + fnbr, v1, dir1.dot(bxis) / scale));
        tripletes.push_back(Trip(i + fnbr, v2, dir2.dot(bxis) / scale));

        tripletes.push_back(Trip(i + fnbr, lsin, -1));

        energy[i + fnbr] = iso.dot(bxis) / scale - sin_angle;
    }
    // sin^2 + cos^2 = 1
    tripletes.push_back(Trip(2 * fnbr, lcos, 2 * cos_angle));
    tripletes.push_back(Trip(2 * fnbr, lsin, 2 * sin_angle));

    energy[2 * fnbr] = sin_angle * sin_angle + cos_angle * cos_angle - 1;

    // axis * axis = 1
    tripletes.push_back(Trip(2 * fnbr + 1, lax, 2 * axis[0]));
    tripletes.push_back(Trip(2 * fnbr + 1, lay, 2 * axis[1]));
    tripletes.push_back(Trip(2 * fnbr + 1, laz, 2 * axis[2]));

    energy[2 * fnbr + 1] = axis.dot(axis) - 1;

    if (fix_axis)
    {
        double scale = axis.norm();
        // axis * axis_fix / scale = 1;
        tripletes.push_back(Trip(2 * fnbr + 2, lax, axis_fix[0] / scale));
        tripletes.push_back(Trip(2 * fnbr + 2, lay, axis_fix[1] / scale));
        tripletes.push_back(Trip(2 * fnbr + 2, laz, axis_fix[2] / scale));

        energy[2 * fnbr + 2] = axis.dot(axis_fix) / scale - 1;
    }
    if (fix_angle)
    {
        // sin_angle - sin_angle_fix = 0
        tripletes.push_back(Trip(2 * fnbr + 3, lsin, 1));
        energy[2 * fnbr + 3] = sin_angle - sin_angle_fix;

        // cos_angle - cos_angle_fix = 0
        tripletes.push_back(Trip(2 * fnbr + 4, lcos, 1));
        energy[2 * fnbr + 4] = cos_angle - cos_angle_fix;
    }
    spMat J;
	J.resize(energy.size(), Glob_lsvars.size());
	J.setFromTriplets(tripletes.begin(), tripletes.end());
	H = J.transpose() * J;
	B = -J.transpose() * energy;
}

void lsTools::Run_ConstSlopeOpt(){
    Eigen::MatrixXd GradValueF, GradValueV;
	Eigen::VectorXd PGEnergy;
	Eigen::VectorXd func = fvalues;

	int vnbr = V.rows();
	int fnbr = F.rows();

	// initialize the level set with some number
	if (func.size() != vnbr)
	{
		if (trace_hehs.size() == 0)
		{
			std::cout << "Please First Set Up The Boundary Condition (Tracing)" << std::endl;
			return;
		}
		initialize_level_set_accroding_to_parametrization();
		std::cout << "level set get initialized" << std::endl;
		return;
	}
	analysis_pseudo_geodesic_on_vertices(func, analizers[0]);

	int final_size = vnbr + 5; // Change this when using more auxilary vars

	//  update quantities associated with level set values

	get_gradient_hessian_values(func, GradValueV, GradValueF);

	if (Glob_lsvars.size() == 0)
	{
		std::cout << "Initializing Global Variable For LevelSet Opt ... " << std::endl;
		Glob_lsvars = Eigen::VectorXd::Zero(final_size); // We change the size if opt more than 1 level set
		Glob_lsvars.segment(0, vnbr) = func;
	}
	if (Last_Opt_Mesh)
	{
		Compute_Auxiliaries = true;
		std::cout << "Recomputing Auxiliaries" << std::endl;
	}

	spMat H;
	H.resize(vnbr, vnbr);
	Eigen::VectorXd B = Eigen::VectorXd::Zero(vnbr);

	spMat LTL;			  // left of laplacian
	Eigen::VectorXd mLTF; // right of laplacian
	assemble_solver_biharmonic_smoothing(func, LTL, mLTF);
	H += weight_laplacian * LTL;
	B += weight_laplacian * mLTF;
	assert(mass.rows() == vnbr);

	Eigen::VectorXd bcfvalue = Eigen::VectorXd::Zero(vnbr);
	Eigen::VectorXd fbdenergy = Eigen::VectorXd::Zero(vnbr);
	// boundary condition (traced as boundary condition)
	if (trace_hehs.size() > 0)
	{ // if traced, we count the related energies in

		spMat bc_JTJ;
		Eigen::VectorXd bc_mJTF;
		assemble_solver_boundary_condition_part(func, bc_JTJ, bc_mJTF, bcfvalue);
		H += weight_boundary * bc_JTJ;
		B += weight_boundary * bc_mJTF;
	}

	// strip width condition
	if (enable_strip_width_energy)
	{
		spMat sw_JTJ;
		Eigen::VectorXd sw_mJTF;
		assemble_solver_strip_width_part(GradValueF, sw_JTJ, sw_mJTF); // by default the strip width is 1. Unless tracing info updated the info
		H += weight_strip_width * sw_JTJ;
		B += weight_strip_width * sw_mJTF;
	}

	assert(H.rows() == vnbr);
	assert(H.cols() == vnbr);
	// pseudo geodesic
	spMat Hlarge;
	Eigen::VectorXd Blarge;
	Hlarge.resize(final_size, final_size);
	Blarge = Eigen::VectorXd::Zero(final_size);

	Hlarge = sum_uneven_spMats(H, Hlarge);
	Blarge = sum_uneven_vectors(B, Blarge);

	if (enable_pseudo_geodesic_energy)
	{

		spMat pg_JTJ;
		Eigen::VectorXd pg_mJTF;
        assemble_solver_constant_slope(func, AxisFixedForSlopes, AxisFixIn, AnglesFixedForSlopes, AngleFixIn, analizers[0],
        pg_JTJ, pg_mJTF, PGEnergy);

		Hlarge = sum_uneven_spMats(Hlarge, weight_pseudo_geodesic_energy * pg_JTJ);
		Blarge = sum_uneven_vectors(Blarge, weight_pseudo_geodesic_energy * pg_mJTF);
	}

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
	// double inf_norm=dx.cwiseAbs().maxCoeff();
	if (enable_pseudo_geodesic_energy) // only pg energy and boundary angles requires step length
	{
		if (level_set_step_length > max_step_length)
		{
			dx *= max_step_length / level_set_step_length;
		}
	}

	func += dx.topRows(vnbr);
	fvalues = func;
	Glob_lsvars += dx;
	double energy_biharmonic = (QcH * func).norm();
	double energy_boundary = (bcfvalue).norm();
	std::cout << "energy: harm " << energy_biharmonic << ", bnd " << energy_boundary << ", ";
	if (enable_pseudo_geodesic_energy)
	{
		double energy_pg = PGEnergy.norm();
		double max_energy_ls = PGEnergy.lpNorm<Eigen::Infinity>();
		std::cout << "pg_energy, " << energy_pg << ", "
				  << "lsmax," << max_energy_ls << ", ";
	}

	if (enable_strip_width_energy)
	{

		Eigen::VectorXd ener = GradValueF.rowwise().norm();
		ener = ener.asDiagonal() * ener;
		Eigen::VectorXd wds = Eigen::VectorXd::Ones(fnbr) * strip_width * strip_width;

		ener -= wds;
		double stp_energy = Eigen::VectorXd(ener).dot(ener);
		std::cout << "strip, " << stp_energy << ", ";
	}
	step_length = dx.norm();
    Eigen::Vector3d axis(Glob_lsvars[vnbr], Glob_lsvars[vnbr + 1], Glob_lsvars[vnbr + 2]);
    double cos_angle = Glob_lsvars[vnbr + 3];
    double sin_angle = Glob_lsvars[vnbr + 4];
    std::cout<<"Axis ("<<axis[0]<<", "<<axis[1]<<", "<<axis[2]<<"), cos_angle "<<cos_angle;
    std::cout << ", step " << step_length << std::endl;

	Last_Opt_Mesh = false;
}