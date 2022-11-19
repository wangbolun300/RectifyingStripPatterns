// This file deals with functional target angle values
#include <lsc/basic.h>
#include <lsc/tools.h>
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
