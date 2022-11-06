#include<lsc/mesh_opt.h>
#include <lsc/igl_tool.h>
void lsTools::compute_derivate_pairs()
{
    int ninner=ActInner.size();
    std::vector<std::array<spMat,8>> MJCmats;// mesh Jocobian coffecient matrices
    std::vector<double> t1,t2;
    t1.resize(ninner);
    t2.resize(ninner);
    
    for (int i = 0; i < ninner; i++)
    {
        if(ActInner[i]==false){
            continue;
        }
        int vm = IVids[i];
        CGMesh::HalfedgeHandle inhd = Vheh0[i], outhd=Vheh1[i];
        int in_from=lsmesh.from_vertex_handle(inhd).idx();
        int in_to=lsmesh.to_vertex_handle(inhd).idx();
        int out_from=lsmesh.from_vertex_handle(outhd).idx();
        int out_to=lsmesh.to_vertex_handle(outhd).idx();
        // int in_large=in_from;
        // int in_small=in_to;
        // int out_large=out_from;
        // int out_small=out_to;
        // if(fvalues[in_from]<fvalues[in_to]){
        //     in_large=in_to;
        //     in_small=in_from;
        // }
        // if(fvalues[out_from]<fvalues[out_to]){
        //     out_large=out_to;
        //     out_small=out_from;
        // }
        t1[i]=get_t_of_value(fvalues[vm], fvalues[in_from],fvalues[in_to]);
        t2[i]=get_t_of_value(fvalues[vm], fvalues[out_from],fvalues[out_to]);
        if(t1[i]<-SCALAR_ZERO||t1[i]>1+SCALAR_ZERO||t2[i]<-SCALAR_ZERO||t2[i]>1+SCALAR_ZERO){
            std::cout<<"ERROR in Mesh Opt: Using Wrong Triangle For Finding Level Set"<<std::endl;
            assert(false);
        }
        if(in_to==out_from){
            //TODO
        }
        if(in_from==out_to){
            //TODO
        }
        

    }
}
