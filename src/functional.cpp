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

// angles are in radian
void assign_angle_get_closer_angle_to_asymptotic(const double percent, const Eigen::VectorXd &angles, const int vnbr,
                                    const std::vector<int> &IVids, std::vector<double> &target)
{
    target.resize(vnbr);
    double aasym = 0; // target angle is 0 or 180. 180 is included in 0.
    for (int i = 0; i < angles.size(); i++)
    {
        double a = angles[i];
        a += (aasym - a) * percent / 100.;
        int vm = IVids[i];
        target[vm] = a;
    }
}

void  lsTools::Run_AsOrthAsPossible_LS(){
    Analyze_Optimized_LS_Angles = true;
    Eigen::MatrixXd GradValueF, GradValueV;
	Eigen::VectorXd PGEnergy;
	Eigen::VectorXd func = fvalues;
	
	std::vector<double> angle_degree; todo;
	
	int vnbr = V.rows();
	int fnbr = F.rows();
	
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
	// if (Glob_lsvars.size() == 0 && trace_hehs.size() == 0)// In this case we don't unit scale the level set
	// {
	// 	levelset_unit_scale(func, GradValueF, strip_width);
	// 	std::cout<<"Unit Scale Levelset"<<std::endl;
	// }
	// std::cout<<"check "<<func.norm()<<std::endl;
	analysis_pseudo_geodesic_on_vertices(func, analizers[0]);
	int ninner = analizers[0].LocalActInner.size();
	int final_size;
    todo: add more auxiliary variables;
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
    if(AnalizedAngelVector.size()!=0){
        
    }
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
        double max_ls_angle_energy = PGEnergy.bottomRows(ninner).norm();
        std::cout << "total angle energy, " << max_ls_angle_energy << ", ";
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

	std::cout << std::endl;
	Last_Opt_Mesh = false;
}