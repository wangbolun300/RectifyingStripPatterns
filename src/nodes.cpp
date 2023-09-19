#include <lsc/basic.h>
#include <lsc/tools.h>
#include <igl/file_dialog_save.h>

// ver and ray defines a line, the norm is the normal vector on the point ver
// please normalize ray and norm before using it
double line_point_signed_distance(const Eigen::Vector3d &pt, const Eigen::Vector3d &ver, const Eigen::Vector3d &ray,
								  const Eigen::Vector3d &norm)
{
	 Eigen::Vector3d bxis = norm.cross(ray).normalized();
	 // norm, ray and bxis is a frame
	 // vp*bxis is the signed distance
	 Eigen::Vector3d vp = pt - ver;
	 double b = bxis.dot(vp);
	 return b;

}

// larger or smaller: 1: larger, -1, smaller
// the poly: a*x*x + b*y + c
void quadratic_function_compare_0(const double a, const double b, const double c, const double xmin, const double xmax,
								  std::vector<int> &largerOrSmaller, std::vector<std::array<double, 2>> &regions)
{
	largerOrSmaller.clear();
	regions.clear();
	std::vector<double> paras(3);
	std::array<double, 2> roots;
	paras[0] = c;
	paras[1] = b;
	paras[2] = a;
	paras = polynomial_simplify(paras);
	bool solved = quadratic_solver(paras, roots);
	if (!solved)
	{ // always larger or less than 0
		if (a > 0)
		{
			largerOrSmaller.push_back(1);
			regions.push_back({xmin, xmax});
			return;
		}
		if (a < 0)
		{
			largerOrSmaller.push_back(-1);
			regions.push_back({xmin, xmax});
			return;
		}
		if(a==0){
			// normally this cannot happen
			std::cout << "ERROR: unexpected in quadratic_function_compare_0, polynomial\n"
					  << a << "," << b << "," << c << std::endl;
		}
		return;
	}
	if (roots[0] == roots[1])
	{
		if (paras.size() == 3)
		{
			if (a > 0)
			{
				largerOrSmaller.push_back(1);
				regions.push_back({xmin, xmax});
				return;
			}
			if (a < 0)
			{
				largerOrSmaller.push_back(-1);
				regions.push_back({xmin, xmax});
				return;
			}
		}
		else{// the function is a straight line
			if (roots[0] >= xmax || roots[0] <= xmin)
			{
				double testvalue = polynomial_value(paras, (xmin + xmax) / 2);

				if (testvalue > 0)
				{
					largerOrSmaller.push_back(1);
					regions.push_back({xmin, xmax});
					return;
				}
				if (testvalue < 0)
				{
					largerOrSmaller.push_back(-1);
					regions.push_back({xmin, xmax});
					return;
				}
			}
			else
			{
				double testvalue0 = polynomial_value(paras, xmin);
				double testvalue1 = polynomial_value(paras, xmax);
				if (testvalue0 > testvalue1)// first large then less
				{ 
					largerOrSmaller.push_back(1);
					regions.push_back({xmin, roots[0]});

					largerOrSmaller.push_back(-1);
					regions.push_back({roots[0], xmax});
					return;
				}
				else// first less then larger
				{
					largerOrSmaller.push_back(-1);
					regions.push_back({xmin, roots[0]});

					largerOrSmaller.push_back(1);
					regions.push_back({roots[0], xmax});
					return;
				}
			}
		}
	}
	else// two different roots
	{
		std::vector<double> ordered_roots;
		if (roots[0] < roots[1])
		{
			ordered_roots.push_back(roots[0]);
			ordered_roots.push_back(roots[1]);
		}
		else
		{
			ordered_roots.push_back(roots[1]);
			ordered_roots.push_back(roots[0]);
		}
		std::vector<double> rinr;// roots in the region
		if (ordered_roots[0] < xmax && ordered_roots[0] > xmin)
		{
			rinr.push_back(ordered_roots[0]);
		}
		if (ordered_roots[1] < xmax && ordered_roots[1] > xmin)
		{
			rinr.push_back(ordered_roots[1]);
		}
		if (rinr.size() == 0)// no root in the region
		{
			double testvalue = polynomial_value(paras, (xmin + xmax) / 2);
			if (testvalue > 0)
			{
				largerOrSmaller.push_back(1);
				regions.push_back({xmin, xmax});
				return;
			}
			else{
				largerOrSmaller.push_back(-1);
				regions.push_back({xmin, xmax});
				return;
			}
		}
		if (rinr.size() == 1)// one root in the region
		{
			double slope = 2 * a * rinr[0] + b; // slope is 2*a*x + b
			if (slope > 0)
			{
				largerOrSmaller.push_back(-1);
				regions.push_back({xmin, rinr[0]});

				largerOrSmaller.push_back(1);
				regions.push_back({rinr[0], xmax});
				return;
			}
			else
			{
				largerOrSmaller.push_back(1);
				regions.push_back({xmin, rinr[0]});

				largerOrSmaller.push_back(-1);
				regions.push_back({rinr[0], xmax});
				return;
			}
		}
		if (rinr.size() == 2)// two roots in the region
		{
			double testvalue0 = polynomial_value(paras, xmin);
			double testvalue1 = polynomial_value(paras, (rinr[0] + rinr[1]) / 2);
			if (testvalue0 < testvalue1)
			{ //<,>,<
				largerOrSmaller.push_back(-1);
				regions.push_back({xmin, rinr[0]});

				largerOrSmaller.push_back(1);
				regions.push_back({rinr[0], rinr[1]});

				largerOrSmaller.push_back(-1);
				regions.push_back({rinr[1], xmax});

				return;
			}
			if (testvalue0 > testvalue1)
			{ //>,<,>
				largerOrSmaller.push_back(1);
				regions.push_back({xmin, rinr[0]});

				largerOrSmaller.push_back(-1);
				regions.push_back({rinr[0], rinr[1]});

				largerOrSmaller.push_back(1);
				regions.push_back({rinr[1], xmax});

				return;
			}
		}

	}
	return;

}

void lsTools::test_quadratic_compare_0(const double a, const double b, const double c, const double xmin, const double xmax)
{
	std::cout << "Testing quadratic polynomials" << std::endl;
	std::vector<int> signs;
	std::vector<std::array<double, 2>> regions;
	quadratic_function_compare_0(a,b,c,xmin,xmax, signs, regions);
	int size = signs.size();
	for(int i=0;i<size;i++){
		std::cout << "sign, " << signs[i] << ", region, [" << regions[i][0] << ", " << regions[i][1] << "]\n";
	}


}
double get_minimum_of_poly3(const double a, const double b, const double c, const double d, const double xmin, const double xmax)
{
	double qa = 3 * a;
	double qb = 2 * b;
	double qc = c;
	std::vector<int> largerOrSmaller;
	std::vector<std::array<double, 2>> regions;
	quadratic_function_compare_0(qa, qb, qc, xmin, xmax, largerOrSmaller, regions);
	int nseg = largerOrSmaller.size();
	assert(nseg == regions.size());
	double result;
	std::vector<double> ply = {d, c, b, a};
	if (nseg == 1)
	{
		if (largerOrSmaller[0] > 0)// increase
		{
			return polynomial_value(ply, xmin);
		}
		else// decrease
		{
			return polynomial_value(ply, xmax);
		}
	}
	if (nseg == 2)
	{
		if(largerOrSmaller[0] > 0){//increase and decrease
			return std::min(polynomial_value(ply, xmin), polynomial_value(ply, xmax));
		}
		else{ // decrease and increase
			double value1 = regions[0][1];
			return polynomial_value(ply, value1);
		}
	}
	if (nseg == 3)
	{
		std::array<double, 4> values;
		values[0] = xmin;
		values[1] = regions[0][1];
		values[2] = regions[1][1];
		values[3] = xmax;
		double minvalue = polynomial_value(ply, xmin);
		for (int i = 1; i < 4; i++)
		{
			double tmp = polynomial_value(ply, values[i]);
			if (tmp < minvalue)
			{
				minvalue = tmp;
			}
		}
		return minvalue;
	}

	std::cout<<"ERROR in get_minimum_of_poly3."<<std::endl;
	exit(0);
}

// the size of radius_in is ninner, the size of radius_out is vnbr
// this function fits X (func) and Y (radius_in^2) into a polynomial
void assign_distances_by_poly_fitting(const Eigen::VectorXd &fvalues, const std::vector<double> &radius_in,
									  const std::vector<int> &IVids, const int ninner, std::vector<double> &radius_out,
									  std::vector<double> &plyout, Eigen::VectorXd &Xori, Eigen::VectorXd &Yori)
{
	// std::cout<<"In 1"<<std::endl;
    Eigen::VectorXd Xs(ninner), Ys(ninner);
    Eigen::MatrixXd A(ninner, 4);
    for (int i = 0; i < ninner; i++)
    {
        double a = radius_in[i] * radius_in[i]; // here we take the squared value, since the radius is actually |radius|
        int vm = IVids[i];
        double x = fvalues[vm];
        Xs[i] = x;
        Ys[i] = a;
        A.row(i) << x * x * x, x * x, x, 1;
        // double angle_degree = a * 180 / LSC_PI;
        // target[vm] = angle_degree;
    }
	Xori = Xs;
	Yori = Ys;
	// std::cout<<"In 2"<<std::endl;
    Eigen::Matrix4d ATA = A.transpose() * A, inv;
    bool invertable;
    ATA.computeInverseWithCheck(inv, invertable);
    if(!invertable){
        std::cout<<"MATRIX NOT INVERTABLE IN void assign_angle_get_closer_angle_to_asymptotic()"<<std::endl;
        return;
    }
    Eigen::Vector4d paras = inv*A.transpose()*Ys;
    // output the fitting results
    double a = paras[0];
    double b = paras[1];
	double c = paras[2];
	double d = paras[3];
	double xmin = fvalues.minCoeff();
	double xmax = fvalues.maxCoeff();
	// make the polynomial always larger than 0
	double minvalue = get_minimum_of_poly3(a, b, c, d, xmin, xmax);
	if (minvalue <= 0)
	{
		d = d - minvalue + 1e-6;
		std::cout<<"lowest: "<<minvalue<<",";
	}
	radius_out.resize(fvalues.size());
	for (int i = 0; i < ninner; i++)
    {
        // i is the id in inner, vm is the id in global
        int vm = IVids[i];
        double f = fvalues[vm];
		double value = a * f * f * f + b * f * f + c * f + d; // this is the squared distance
		if (value < 0)
		{
			std::cout << "ERROR in assign_distances_by_poly_fitting: the polynomial should be larger than 0" << std::endl;
			exit(0);
		}
		value = sqrt(value);
		radius_out[vm] = value;
	}
	std::cout << "TargetFitting, " << a << ", " << b << ", " << c << ", " << d << ", ";
	plyout = {d, c, b, a};
}

void lsTools::assemble_solver_co_point_ruling_node(const Eigen::VectorXd &vars,
                                                   const bool fix_centroid, const Eigen::Vector3d &centroid_in, 
												   const bool fix_radius, const std::vector<double>& radius_input,
                                                   const LSAnalizer &analizer,
                                                   spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy)
{
	// std::cout << "check 1" << std::endl;

	int vnbr = V.rows();
	int ninner = analizer.LocalActInner.size();
	RadiusCollected.resize(ninner);
    std::vector<Trip> tripletes;
	tripletes.reserve(ninner * 15);				// the number of rows is ninner*4, the number of cols is aux_start_loc + ninner * 3 (all the function values and auxiliary vars)
	energy = Eigen::VectorXd::Zero(ninner * 3); // mesh total energy values
    // bxis0 = Eigen::MatrixXd::Zero(ninner, 3);
    // bxis1 = bxis0;

	Eigen::Vector3d centroid_fix = centroid_in;

	// the sphere centroid
	int lox = vnbr + ninner * 5;
	int loy = vnbr + ninner * 5 + 1;
	int loz = vnbr + ninner * 5 + 2;
	


	if (Compute_Auxiliaries)
    {
		if(fix_centroid){
			Glob_lsvars[lox] = centroid_fix[0];
			Glob_lsvars[loy] = centroid_fix[1];
			Glob_lsvars[loz] = centroid_fix[2];
		}
		else{
			Glob_lsvars[lox] = 0;
			Glob_lsvars[loy] = 0;
			Glob_lsvars[loz] = 1;
		}

	}
	Eigen::Vector3d centroid(Glob_lsvars[lox], Glob_lsvars[loy], Glob_lsvars[loz]);
	if(fix_centroid){
		// std::cout << "fix centroid, " << centroid_fix.transpose() << "\n";
		centroid = centroid_fix;
		Glob_lsvars[lox] = centroid_fix[0];
		Glob_lsvars[loy] = centroid_fix[1];
		Glob_lsvars[loz] = centroid_fix[2];
	} 

    for (int i = 0; i < ninner; i++)
	{
		if (analizer.LocalActInner[i] == false)
		{
			std::cout << "singularity" << std::endl;
			continue;
		}
		
		
		int vm = IVids[i];

		CGMesh::HalfedgeHandle inhd = analizer.heh0[i], outhd = analizer.heh1[i];
		int v1 = lsmesh.from_vertex_handle(inhd).idx();
		int v2 = lsmesh.to_vertex_handle(inhd).idx();
		int v3 = lsmesh.from_vertex_handle(outhd).idx();
		int v4 = lsmesh.to_vertex_handle(outhd).idx();

		double t1 = analizer.t1s[i];
		double t2 = analizer.t2s[i];
		
		Eigen::Vector3d ver0 = V.row(v1) + (V.row(v2) - V.row(v1)) * t1;
		Eigen::Vector3d ver1 = V.row(vm);
		Eigen::Vector3d ver2 = V.row(v3) + (V.row(v4) - V.row(v3)) * t2;
		// the locations
		int lvm = vm;
		int lv1 = v1;
		int lv2 = v2;
		int lv3 = v3;
		int lv4 = v4;
		// ruling
		int lrx = vnbr + i;
		int lry = vnbr + i + ninner;
		int lrz = vnbr + i + ninner * 2;
		// the auxiliary variable (centroid - pm) * tangent / tannorm
		int lr = vnbr + i + ninner * 3;
		// the radius
		int lra = vnbr + i + ninner * 4;
		
		double fm = Glob_lsvars[lvm];
        double f1 = Glob_lsvars[lv1];
        double f2 = Glob_lsvars[lv2];
        double f3 = Glob_lsvars[lv3];
        double f4 = Glob_lsvars[lv4];
        Eigen::Vector3d pm = V.row(vm);
        Eigen::Vector3d p1 = V.row(v1);
        Eigen::Vector3d p2 = V.row(v2);
        Eigen::Vector3d p3 = V.row(v3);
        Eigen::Vector3d p4 = V.row(v4);

		// the contact point on the sphere
		
		// the tangent is 
        // (f2-f1) * (f4-fm) * v3 + (f2-f1) * (fm-f3) * v4 - (f4-f3) * (f2-fm) * v1 - (f4-f3) * (fm-f1) * v2
		Eigen::Vector3d tangent = (f2 - f1) * (f4 - fm) * p3 + (f2 - f1) * (fm - f3) * p4 - (f4 - f3) * (f2 - fm) * p1 - (f4 - f3) * (fm - f1) * p2;
		double tannorm = tangent.norm();
		Eigen::Vector3d norm = norm_v.row(vm);
		Eigen::Vector3d vo = centroid - pm;
		// ruling
		Eigen::Vector3d real_r = (ver1 - ver0).cross(ver2 - ver1).normalized();
		if (Compute_Auxiliaries)
		{
			
			Glob_lsvars[lrx] = real_r[0];
			Glob_lsvars[lry] = real_r[1];
			Glob_lsvars[lrz] = real_r[2];

			double real_projection = (centroid - pm).dot(tangent) / tannorm;
			Glob_lsvars[lr] = real_projection;// don't care about orientation anymore
			Glob_lsvars[lra] = sqrt((centroid - pm).dot(centroid - pm) - real_projection * real_projection);
			RadiusCollected[i] = Glob_lsvars[lra];
			continue;
		}
		RadiusCollected[i] = Glob_lsvars[lra];
		Eigen::Vector3d ruling(Glob_lsvars[lrx], Glob_lsvars[lry], Glob_lsvars[lrz]);
		double r = Glob_lsvars[lr];
		double radius = Glob_lsvars[lra];

		// the weights
		// double dis0 = ((V.row(v1) - V.row(v2)) * vars[lvm] + (V.row(v2) - V.row(vm)) * vars[lv1] + (V.row(vm) - V.row(v1)) * vars[lv2]).norm();
		// double dis1 = ((V.row(v3) - V.row(v4)) * vars[lvm] + (V.row(v4) - V.row(vm)) * vars[lv3] + (V.row(vm) - V.row(v3)) * vars[lv4]).norm();

		double dis10 = (ver1 - ver0).norm();
		double dis21 = (ver2 - ver1).norm();

		// the 1st method: express the rulings and make it passing through the centroid

		// the 2nd method: the tangent vector is tangent to a sphere
		
		// condition 1:  (centroid - pm) * tangent / tnorm = r

		if(!fix_centroid){
			tripletes.push_back(Trip(i, lox, tangent[0] / tannorm));
			tripletes.push_back(Trip(i, loy, tangent[1] / tannorm));
			tripletes.push_back(Trip(i, loz, tangent[2] / tannorm));
		}
		tripletes.push_back(Trip(i, lr, -1));

		if (!RulingOnlyFindingSphere)
		{
			double v1b = p1.dot(centroid - pm) / tannorm;
			double v2b = p2.dot(centroid - pm) / tannorm;
			double v3b = p3.dot(centroid - pm) / tannorm;
			double v4b = p4.dot(centroid - pm) / tannorm;

			tripletes.push_back(Trip(i, lvm, -(f2 - f1) * v3b + (f2 - f1) * v4b + (f4 - f3) * v1b - (f4 - f3) * v2b));
			tripletes.push_back(Trip(i, lv1, -(f4 - fm) * v3b - (fm - f3) * v4b + (f4 - f3) * v2b));
			tripletes.push_back(Trip(i, lv2, (f4 - fm) * v3b + (fm - f3) * v4b - (f4 - f3) * v1b));
			tripletes.push_back(Trip(i, lv3, -(f2 - f1) * v4b + (f2 - fm) * v1b + (fm - f1) * v2b));
			tripletes.push_back(Trip(i, lv4, (f2 - f1) * v3b - (f2 - fm) * v1b - (fm - f1) * v2b));
		}
		else{
			// std::cout << "Finding only sphere" << std::endl;
		}

		energy[i] = (centroid - pm).dot(tangent) / tannorm - r;

		// condition 2:  (centroid - pm) ^2 - r^2 = radius^2
		if (!fix_centroid)
		{
			tripletes.push_back(Trip(i + ninner, lox, 2 * centroid[0] - 2 * pm[0]));
			tripletes.push_back(Trip(i + ninner, loy, 2 * centroid[1] - 2 * pm[1]));
			tripletes.push_back(Trip(i + ninner, loz, 2 * centroid[2] - 2 * pm[2]));
		}

		tripletes.push_back(Trip(i + ninner, lr, -2 * r));
		tripletes.push_back(Trip(i + ninner, lra, -2 * radius));

		energy[i + ninner] = (centroid - pm).dot(centroid - pm) - r * r - radius * radius;

		// condition 3:  radius ^ 2 - radius_input[vm] ^ 2 = 0
		tripletes.push_back(Trip(i + ninner * 2, lra, 2 * radius));
		energy[i + ninner * 2] = radius * radius - radius_input[vm] * radius_input[vm];
		// if(i==5){
		// 	std::cout<<"The 5th target, "<<radius_input[vm]<<", ";
		// }

        
    }
	if (Compute_Auxiliaries)
	{
		H = spMat(Eigen::VectorXd::Zero(Glob_lsvars.size()).asDiagonal());
		B = Eigen::VectorXd::Zero(Glob_lsvars.size());
		return;
	}

	std::cout << "Ruling 1: " << energy.segment(0, ninner).norm() << ", " << energy.segment(ninner, ninner).norm() << ", "
			  << energy.segment(ninner * 2, ninner).norm() << ", 3th avg: " << energy.segment(ninner * 2, ninner).norm() / ninner;
	spMat J;
	J.resize(energy.size(), Glob_lsvars.size());
	J.setFromTriplets(tripletes.begin(), tripletes.end());
	H = J.transpose() * J;
	B = -J.transpose() * energy; 
}

void lsTools::statistic_tangent_vector_to_pt_distance(){
	int vnbr = V.rows();
	analysis_pseudo_geodesic_on_vertices(fvalues, analizers[0]);
	LSAnalizer analizer = analizers[0];
	int ninner = analizer.LocalActInner.size();
	Eigen::Vector3d centroid = RulingCentroid;
	TangentE0 = Eigen::MatrixXd(ninner, 3);
	TangentE1 = TangentE0;
	NodeProjectPts = TangentE0;
	
	std::ofstream file;
	if (write_ls_radius_file)
	{
		std::string fname;
		fname = igl::file_dialog_save();
		if (fname.length() == 0)
		{
			std::cout << "Please type down the prefix" << std::endl;
			return;
		}
		
		file.open(fname + ".csv");
	}
	for (int i = 0; i < ninner; i++)
	{
		if (analizer.LocalActInner[i] == false)
		{
			std::cout << "singularity" << std::endl;
			continue;
		}
		int vm = IVids[i];

		CGMesh::HalfedgeHandle inhd = analizer.heh0[i], outhd = analizer.heh1[i];
		int v1 = lsmesh.from_vertex_handle(inhd).idx();
		int v2 = lsmesh.to_vertex_handle(inhd).idx();
		int v3 = lsmesh.from_vertex_handle(outhd).idx();
		int v4 = lsmesh.to_vertex_handle(outhd).idx();

		double t1 = analizer.t1s[i];
		double t2 = analizer.t2s[i];
		
		Eigen::Vector3d ver0 = V.row(v1) + (V.row(v2) - V.row(v1)) * t1;
		Eigen::Vector3d ver1 = V.row(vm);
		Eigen::Vector3d ver2 = V.row(v3) + (V.row(v4) - V.row(v3)) * t2;
		Eigen::Vector3d tangent = (ver2 - ver0).normalized();
		double dis1 = (centroid - ver1).dot(tangent);
		dis1 = dis1 * dis1;
		double dis2 = (centroid - ver1).dot(centroid - ver1);
		double radius = sqrt(dis2 - dis1);
		std::cout<<"the "<<i<<"th, length: "<<radius<<std::endl;
		TangentE0.row(i) = ver1;
		TangentE1.row(i) = ver1 + tangent;
		NodeProjectPts.row(i) = ver1 + tangent * (centroid - ver1).dot(tangent);
		if (write_ls_radius_file)
		{
			file << fvalues[vm] << "," << radius << std::endl;
		}
	}
	if (write_ls_radius_file)
	{
		file.close();
	}
}

std::vector<double> polynormial_sqrt_values(const std::vector<double> &poly, const Eigen::VectorXd &xvalues)
{
	int vnbr = xvalues.size();
	std::vector<double> result(vnbr);
	for (int i = 0; i < vnbr; i++)
	{
		double value = polynomial_value(poly, xvalues[i]);
		if (value < 0)
		{
			std::cout<<"ERROR in polynormial_sqrt_values: the value should be larger than 0"<<std::endl;
			exit(0);
		}
		result[i] = sqrt(value);
	}
	return result;
}
void lsTools::Run_ruling_opt()
{
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
	int ninner = analizers[0].LocalActInner.size();
	// the variables: scalar field, the rulings, the auxiliary projection value, the radius of the sphere, the centroid of the sphere.
	int final_size = vnbr + ninner * 5 + 3; // Change this when using more auxilary vars

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
		std::vector<double> radius_fitting;
		if (!Compute_Auxiliaries)
		{
			if(RulingFixRadius){
				// radius_fitting = polynormial_sqrt_values(NodeFitPly, func);
				// std::cout<<"Fixed Ply: "<<NodeFitPly[3]<<", "<<NodeFitPly[2]<<", "<<NodeFitPly[1]<<", "<<NodeFitPly[0]<<", ";
				radius_fitting = NodeFitTargetRadius;
			}
			else{
				assign_distances_by_poly_fitting(func, RadiusCollected, IVids, ninner, radius_fitting, NodeFitPly, XsNode, YsNode);
				NodeFitTargetRadius = radius_fitting;
			}
			
		}
		assemble_solver_co_point_ruling_node(func, RulingFixCentroid, RulingCentroid, RulingFixRadius, radius_fitting, analizers[0],
        pg_JTJ, pg_mJTF, PGEnergy);

		Hlarge = sum_uneven_spMats(Hlarge, weight_pseudo_geodesic_energy * pg_JTJ);
		Blarge = sum_uneven_vectors(Blarge, weight_pseudo_geodesic_energy * pg_mJTF);
        Compute_Auxiliaries = false;
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
	Eigen::Vector3d centroid(Glob_lsvars[vnbr + ninner * 5], Glob_lsvars[vnbr + ninner * 5 + 1], Glob_lsvars[vnbr + ninner * 5 + 2]);
    std::cout<<"Centroid ("<<centroid[0]<<", "<<centroid[1]<<", "<<centroid[2]<<")";
    std::cout << ", step " << step_length << std::endl;

	Last_Opt_Mesh = false;
}