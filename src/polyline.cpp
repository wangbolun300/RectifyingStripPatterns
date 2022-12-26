#include <lsc/basic.h>
#include<lsc/tools.h>

class PolyOpt
{
public:
    PolyOpt(){};
    void init(const std::vector<std::vector<Eigen::Vector3d>> &ply, const std::vector<std::vector<Eigen::Vector3d>> &bi);
    Eigen::VectorXd PlyVars; // vars for the polylines
    Eigen::VectorXd OriVars; // original vars
    Eigen::VectorXi Front;   // the vertex id of the front point
    Eigen::VectorXi Back;    // the vertex id of the back point

    double weight_binsmooth;
    double weight_plysmooth;
    double weight_mass;
    double weight_binormal;

    void opt();
    // make the values not far from original data
    void assemble_gravity(spMat& H, Eigen::VectorXd& B, Eigen::VectorXd &energy);
    void assemble_polyline_smooth(spMat& H, Eigen::VectorXd& B, Eigen::VectorXd &energy);
};

void PolyOpt::init(const std::vector<std::vector<Eigen::Vector3d>> &ply, const std::vector<std::vector<Eigen::Vector3d>> &bi){
    int vnbr = 0;
    for (auto line : ply)
    {
        for (auto ver : line)
        {
            vnbr++;
        }
    }
    // x0, x1, ..., xn, y0, ..., zn, bx0, ... bxn, ..., bzn
    PlyVars.resize(vnbr * 2 * 3);
    Front.resize(vnbr);
    Back.resize(vnbr);
    int counter = 0;
    for (int i = 0; i < ply.size(); i++)
    {
        for (int j = 0; j < ply[i].size(); j++)
        {
            int lpx = counter;
            int lpy = counter + vnbr;
            int lpz = counter + vnbr * 2;
            int lnx = counter + vnbr * 3;
            int lny = counter + vnbr * 4;
            int lnz = counter + vnbr * 5;
            Eigen::Vector3d ver = ply[i][j];
            Eigen::Vector3d bin = bi[i][j];
            PlyVars[lpx] = ver[0];
            PlyVars[lpy] = ver[1];
            PlyVars[lpz] = ver[2];
            PlyVars[lnx] = bin[0];
            PlyVars[lny] = bin[1];
            PlyVars[lnz] = bin[2];
            if (j == 0)
            {
                Front[counter] = -1;
            }
            else
            {
                Front[counter] = counter - 1;
            }
            if (j == ply[i].size() - 1)
            {
                Back[counter] = -1;
            }
            else
            {
                Back[counter] = counter + 1;
            }

            counter++;

        }
    }
    OriVars = PlyVars;
}

void PolyOpt::assemble_gravity(spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy)
{
    if (Front.size() == 0)
    {
        std::cout << "please use assemble_gravity after init the PolyOpt" << std::endl;
        return;
    }

    // std::vector<Trip> tripletes;
    int vnbr = Front.size();
    H = Eigen::VectorXd::Ones(vnbr * 6).asDiagonal();
    energy = Eigen::VectorXd::Zero(vnbr * 6);
    // tripletes.reserve(vnbr * 6);
    // for (int i = 0; i < PlyVars.size(); i++)
    // {
    //     int rep = i / vnbr;
    //     int vid = i - rep * vnbr;
    //     energy[i] =
    //     // tripletes.push_back
    // }
    energy = PlyVars - OriVars;

    B = -energy;
}
void PolyOpt::assemble_polyline_smooth(spMat& H, Eigen::VectorXd& B, Eigen::VectorXd &energy){
    if (Front.size() == 0)
    {
        std::cout << "please use assemble_polyline_smooth after init the PolyOpt" << std::endl;
        return;
    }

    std::vector<Trip> tripletes;
    int vnbr = Front.size();
    energy = Eigen::VectorXd::Zero(vnbr * 6);
    tripletes.reserve(vnbr * 6);
    for (int i = 0; i < PlyVars.size(); i++)
    {
        int rep = i / vnbr;
        int vid = i - rep * vnbr;
        int fr = Front[vid];
        int bk = Back[vid];
        if(fr == -1 || bk == -1){
            continue;
        }
        // locations
        int lx = vid;
        int ly = vid + vnbr;
        int lz = vid + vnbr * 2;
        int lnx = vid + vnbr * 3;
        int lny = vid + vnbr * 4;
        int lnz = vid + vnbr * 5;
        int lfx = fr;
        int lfy = fr + vnbr;
        int lfz = fr + vnbr * 2;
        int lbx = bk;
        int lby = bk + vnbr;
        int lbz = bk + vnbr * 2;
        int lfnx = fr + vnbr * 3;
        int lfny = fr + vnbr * 4;
        int lfnz = fr + vnbr * 5;
        int lbnx = bk + vnbr * 3;
        int lbny = bk + vnbr * 4;
        int lbnz = bk + vnbr * 5;
        Eigen::Vector3d ver(PlyVars[lx], PlyVars[ly], PlyVars[lz]);
        Eigen::Vector3d vfr(PlyVars[lfx], PlyVars[lfy], PlyVars[lfz]);
        Eigen::Vector3d vbk(PlyVars[lbx], PlyVars[lby], PlyVars[lbz]);
        Eigen::Vector3d nbi(PlyVars[lnx], PlyVars[lny], PlyVars[lnz]);
        Eigen::Vector3d nfr(PlyVars[lfnx], PlyVars[lfny], PlyVars[lfnz]);
        Eigen::Vector3d nbk(PlyVars[lbnx], PlyVars[lbny], PlyVars[lbnz]);
        // 2*x = x_front + x_back
        tripletes.push_back(Trip(i, lx, 2));
        tripletes.push_back(Trip(i, lfx, -1));
        tripletes.push_back(Trip(i, lbx, -1));

        energy[i] = 2 * ver[0] - vfr[0] - vbk[0];

        // 2*y = y_front + y_back
        tripletes.push_back(Trip(i + vnbr, ly, 2));
        tripletes.push_back(Trip(i + vnbr, lfy, -1));
        tripletes.push_back(Trip(i + vnbr, lby, -1));

        energy[i + vnbr] = 2 * ver[1] - vfr[1] - vbk[1];

        // 2*z = z_front + z_back
        tripletes.push_back(Trip(i + vnbr * 2, lz, 2));
        tripletes.push_back(Trip(i + vnbr * 2, lfz, -1));
        tripletes.push_back(Trip(i + vnbr * 2, lbz, -1));

        energy[i + vnbr * 2] = 2 * ver[2] - vfr[2] - vbk[2];
        // tripletes.push_back
    }
}
void PolyOpt::opt()
{
}