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

void lsTools::calculate_pseudo_energy_function_values_vertex_based(const std::vector<double>& angle_degree, Eigen::VectorXd &lens)
{
    int ninner=IVids.size();
    VPEvalue.resize(ninner);
    LsOrient.resize(ninner);// first make g1 g2 same direction, then check if g1xg2 or g2xg1.
    PeWeight.resize(ninner);
    lens.resize(ninner);
    vBinormal.resize(ninner,3);
    
    for (int i = 0; i < ninner; i++)
    {
        
        if(ActInner[i] == false){// this is a singularity
            continue;
        }
        int vid = IVids[i];
        double angle_radian = angle_degree[i] * LSC_PI / 180.; // the angle in radian
        double cos_angle = cos(angle_radian);
        CGMesh::HalfedgeHandle heh1=Vheh0[i];
        CGMesh::HalfedgeHandle heh2=Vheh1[i];
        int fid1, fid2;
        fid1=lsmesh.face_handle(heh1).idx();
        fid2=lsmesh.face_handle(heh2).idx();
        assert(fid1!=fid2);
        Eigen::Vector3d norm = norm_v.row(vid);
        Eigen::Vector3d norm1=norm_f.row(fid1);
        Eigen::Vector3d norm2=norm_f.row(fid2);
        Eigen::Vector3d g1=gfvalue.row(fid1);
        Eigen::Vector3d g2=gfvalue.row(fid2);
        // rotate gradients to get iso-curve directions
        g1=norm1.cross(g1);
        g2=norm2.cross(g2);
        Eigen::Vector3d g1xg2=g1.cross(g2);
        // deal with the orientation of the level set: the left of the curve should be bigger value
        double flag1=g1.dot(Vdire0[i]);
        double flag2=g2.dot(Vdire1[i]);
        double flag=flag1*flag2;
        if(flag>=0){
            LsOrient[i]=1;
        }
        else{
            LsOrient[i]=-1;
        }

        // deal with PeWeight
        double lg1=g1.norm();
        double lg2=g2.norm();
        if(lg1<1e-16 || lg2<1e-16){
            PeWeight[i]=0;// it means the g1xg2 will not be accurate
        }
        else{
            Eigen::Vector3d g1n = g1.normalized();
            Eigen::Vector3d g2n = LsOrient[i]*g2.normalized();// get the oriented g2, to make g1 g2 goes to the same direction
            vBinormal.row(i)=g1n.cross(g2n).normalized();
            double cos_real=g1n.cross(g2n).normalized().dot(norm);// angle between the binormal and the surface normal
            if(cos_real<-1||cos_real>1){
                std::cout<<"angle is wrong"<<std::endl;
            }
            if(cos_angle<-1||cos_angle>1){
                std::cout<<"angle is wrong1"<<std::endl;
            }
            double cos_diff = fabs(cos_real-cos_angle)/2;
            double angle_real=acos(g1n.dot(g2n)) * 180 / LSC_PI;// the angle between the two iso-lines
            // cos = 1, weight is 1; cos = -1, means it is very sharp turn, weight = 0
            // it is still resonable for vertex-based method. since there can be multiple intersections between the
            // osculating plane and the one-ring of vertex
            //std::cout<<"angle is "<<angle_real<<std::endl;
            // if (angles_match(angle_real, 0))// if it is close to straight line, we skip.
            // {
            //     PeWeight[i] = 0;
            // }
            // else
                PeWeight[i] = cos_diff;
        }

        // if(fvalues[v1]<fvalues[v2]){ // orient g1xg2 or g2xg1
        //     LsOrient[i]*=-1;
        // }// we don't need it for vertex based method, since it is implicitly given by Vheh0 and Vheh1
        double length=g1xg2.norm();
        double value;
        if(length<1e-16){
            value=1e6;// if it is really small, we set the value very big
        }
        else// divided by the length to avoid straight lines
        {
            value = LsOrient[i] * g1xg2.dot(norm) / length - cos_angle;
        }
        VPEvalue.coeffRef(i)=value;
        lens[i]=length;
    }
}