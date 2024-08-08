#include <lsc/basic.h>
#include <lsc/tools.h>
#include<igl/file_dialog_open.h>
#include<igl/file_dialog_save.h>
#include <igl/point_mesh_squared_distance.h>

std::vector<std::vector<int>> load_row_col_info(const std::vector<std::vector<double>> &rowstmp)
{
    std::vector<std::vector<double>> tmp = rowstmp;
    std::vector<std::vector<int>> result;
    result.resize(rowstmp.size());
    for (int i = 0; i < result.size(); i++)
    {
        std::vector<int> curve;
        curve.resize(rowstmp[i].size());
        for (int j = 0; j < curve.size(); j++)
        {
            curve[j] = rowstmp[i][j];
        }
        result[i]=curve;
    }
    return result;
}
void load_row_col_info(const std::vector<std::vector<double>> &rowstmp, Eigen::VectorXi &row_f, Eigen::VectorXi &row_b)
{
    for (int i = 0; i < rowstmp.size(); i++)
    {
        for (int j = 0; j < rowstmp[i].size(); j++)
        {
            int vid = rowstmp[i][j];
            if (vid < 0)
            {
                continue;
            }
            int front = -1;
            int back = -1;
            if (j != 0)
            {
                front = rowstmp[i][j - 1];
            }
            if (j != rowstmp[i].size() - 1)
            {
                back = rowstmp[i][j + 1];
            }
            row_f[vid] = front;
            row_b[vid] = back;
        }
    }
    return;
}
// d0 is in between of rf and cf.
// d1 is in between of rb and cb.
void load_diagonal_info(const Eigen::VectorXi &row_front, const Eigen::VectorXi &row_back, const Eigen::VectorXi &col_front,
                        const Eigen::VectorXi &col_back, Eigen::VectorXi &d0_front, Eigen::VectorXi &d1_front,
                        Eigen::VectorXi &d0_back, Eigen::VectorXi &d1_back)
{
    int vnbr = row_front.size();
    d0_front = Eigen::VectorXi::Ones(vnbr) * -1;
    d0_back = d0_front;
    d1_front = d0_front;
    d1_back = d0_front;
    for (int i = 0; i < vnbr; i++)
    {
        int id = i;
        int rfid = row_front(id);
        int rbid = row_back(id);
        int cfid = col_front(id);
        int cbid = col_back(id);
        if (rfid < 0 || rbid < 0 || cfid < 0 || cbid < 0) // if it is a boundary vertex, skip
        {
            continue;
        }
        if (cfid >= 0)
        {
            int lu = row_front(cfid); // upleft
            d0_front[id] = lu;
            int ld = row_back(cfid);
            d1_front[id] = ld;
        }
        if (cbid >= 0)
        {
            int ru = row_front(cbid);
            int rd = row_back(cbid);
            d0_back[id] = rd;
            d1_back[id] = ru;
        }
    }
}

// get the row and col info directly from openmesh.
void getQuadRowCols(CGMesh &mesh_in, Eigen::VectorXi &rf, Eigen::VectorXi &rb, Eigen::VectorXi &cf, Eigen::VectorXi &cb)
{
    const int vnbr = mesh_in.n_vertices();
    const int fnbr = mesh_in.n_faces();
    rf = -Eigen::VectorXi::Ones(vnbr);
    rb = rf;
    cf = rf;
    cb = rf;
    for (int i = 0; i < fnbr; i++)
    {
        CGMesh::FaceHandle fh = mesh_in.face_handle(i); // get the face handle of face i
        int counter = 0;
        for (CGMesh::FaceHalfedgeIter voh = mesh_in.fh_begin(fh); voh != mesh_in.fh_end(fh); ++voh)
        {
            CGMesh::HalfedgeHandle he = voh.current_halfedge_handle();
            int vstart = mesh_in.from_vertex_handle(he).idx();
            int vend = mesh_in.to_vertex_handle(he).idx();
            if (counter == 0)
            {
                rb[vstart] = vend;
                rf[vend] = vstart;
            }
            if (counter == 2)
            {
                rb[vend] = vstart;
                rf[vstart] = vend;
            }
            if(counter == 1)
            {
                cb[vstart] = vend;
                cf[vend] = vstart;
            }
            if(counter == 3)
            {
                cb[vend] = vstart;
                cf[vstart] = vend;
            }
            
            counter++;
        }
    }
    // validation
    for (int i = 0; i < vnbr; i++)
    {
        std::vector<int> ids;
        ids.reserve(5);
        ids.push_back(i);
        ids.push_back(cf[i]);
        ids.push_back(cb[i]);
        ids.push_back(rb[i]);
        ids.push_back(rf[i]);
        for (int j = 0; j < ids.size(); j++)
        {
            int id0 = ids[j];
            if (id0 < 0)
            {
                continue;
            }
            for (int k = 0; k < ids.size(); k++)
            {
                int id1 = ids[k];
                if (j == k || id1 < 0)
                {
                    continue;
                }
                if (id0 == id1)
                {
                    std::cout << "ERROR in getQuadRowCols, wrong neighboring info!\n";
                    std::cout << ids[0] << ", " << ids[1] << ", " << ids[2] << ", " << ids[3] << ", " << ids[4] << ", \n";
                    std::cout<<"rf, "<<rf.transpose()<<"\n";
                    std::cout<<"rb, "<<rb.transpose()<<"\n";
                    std::cout<<"cf, "<<cf.transpose()<<"\n";
                    std::cout<<"cb, "<<cb.transpose()<<"\n";
                    exit(0);
                }
            }
        }
    }
}

// get the row and col info from a regular mesh, prescribed vnbr in one row: rnbr.
void getQuadRowColsRegular(const int vnbr, const int rnbr, Eigen::VectorXi &rf, Eigen::VectorXi &rb, Eigen::VectorXi &cf, Eigen::VectorXi &cb)
{
    int cnbr = vnbr / rnbr;
    rf = -Eigen::VectorXi::Ones(vnbr);
    rb = rf;
    cf = rf;
    cb = rf;
    for (int i = 0; i < cnbr; i++)
    {
        for (int j = 0; j < rnbr; j++)
        {
            int vid = i * rnbr + j;
            assert(vid>=0);
            assert(vid<vnbr);
            int vfr, vbk, vlf, vrt;
            // front
            if (j == 0)
            {
                vfr = -1;
            }
            else
            {
                vfr = vid - 1;
            }
            // back
            if (j == rnbr - 1)
            {
                vbk = -1;
            }
            else
            {
                vbk = vid + 1;
            }
            // left
            if (i == 0)
            {
                vlf = -1;
            }
            else
            {
                vlf = vid - rnbr;
            }
            // right
            if (i == cnbr - 1)
            {
                vrt = -1;
            }
            else
            {
                vrt = vid + rnbr;
            }

            // assign the values.
            rf[vid] = vfr;
            rb[vid] = vbk;
            cf[vid] = vlf;
            cb[vid] = vrt;
            assert(vfr<vnbr);
            assert(vbk<vnbr);
            assert(vlf<vnbr);
            assert(vrt<vnbr);
        }
    }
}
void QuadOpt::init(CGMesh &mesh_in, const std::string &prefix)
{
    MP.mesh2Matrix(mesh_in, V, F);
    std::vector<std::vector<double>> rowstmp;
    std::vector<std::vector<double>> colstmp;
    read_csv_data_lbl(prefix + "_rows.csv", rowstmp);
    read_csv_data_lbl(prefix + "_cols.csv", colstmp);
    rowinfo = rowstmp;
    colinfo = colstmp;
    // for (int i = 0; i < rowinfo.size(); i++)
    // {
    //     std::cout << "init, row size, " << i << "th, " << rowinfo[i].size() << std::endl;
    // }
    // for (int i = 0; i < colinfo.size(); i++)
    // {
    //     std::cout << "init, col size, " << i << "th, " << colinfo[i].size() << std::endl;
    // }

    int vnbr = V.rows();
    if (OptType == 0)// AAG
    {
        // vertices vnbr * 3, normals vnbr * 3, binormals for G vnbr * 3  
        varsize = vnbr * 9;
    }
    if (OptType == 1)// AGG
    {
        // vertices vnbr * 3, normals vnbr * 3, binormals for G0 vnbr * 3, binormals for G1 vnbr * 3.    
        varsize = vnbr * 12;
    }
    if (OptType == 2)// PP
    {
        // vertices vnbr * 3, normals vnbr * 3, binormals + side vectors:  (vnbr * 6) * 2
        varsize = vnbr * 18;
   
    }
    if (OptType == 3)// PPG
    {
        // vertices vnbr * 3, normals vnbr * 3, binormals + side vectors:  (vnbr * 6) * 2 + vnbr * 3
        varsize = vnbr * 21;
    }
    if (OptType == 4)// AAGG
    {
        // vertices vnbr * 3, normals vnbr * 3, binormals for G vnbr * 3 * 2  
        varsize = vnbr * 12;

    }
    
    row_front = Eigen::VectorXi::Ones(vnbr) * -1;
    row_back = Eigen::VectorXi::Ones(vnbr) * -1;
    col_front = Eigen::VectorXi::Ones(vnbr) * -1;
    col_back = Eigen::VectorXi::Ones(vnbr) * -1;


    load_row_col_info(rowstmp, row_front,row_back);
    load_row_col_info(colstmp, col_front, col_back);

    load_diagonal_info(row_front, row_back, col_front, col_back, d0_front, d1_front, d0_back, d1_back);
    OrigVars = Eigen::VectorXd::Zero(varsize);
    OrigVars.segment(0, vnbr) = V.col(0);
    OrigVars.segment(vnbr, vnbr) = V.col(1);
    OrigVars.segment(vnbr * 2, vnbr) = V.col(2);
    mesh_original = mesh_in;
    mesh_update = mesh_original;
    ComputeAuxiliaries = true;
    Estimate_PG_Angles = true;
    Eigen::VectorXd grav_vec = Eigen::VectorXd::Ones(vnbr * 3);
    Eigen::VectorXd var_vec = Eigen::VectorXd::Zero(varsize);
    var_vec.segment(0, vnbr * 3) = grav_vec;
    gravity_matrix = var_vec.asDiagonal();
    Eigen::Vector3d vmin(V.col(0).minCoeff(), V.col(1).minCoeff(), V.col(2).minCoeff());
    Eigen::Vector3d vmax(V.col(0).maxCoeff(), V.col(1).maxCoeff(), V.col(2).maxCoeff());
    std::cout<<"Initialized quad mesh: Vnbr, "<<V.rows()<<", BBD, "<<(vmin - vmax).norm()<<std::endl;;
}

void QuadOpt::init(CGMesh &mesh_in)
{
    MP.mesh2Matrix(mesh_in, V, F);
    // for (int i = 0; i < rowinfo.size(); i++)
    // {
    //     std::cout << "init, row size, " << i << "th, " << rowinfo[i].size() << std::endl;
    // }
    // for (int i = 0; i < colinfo.size(); i++)
    // {
    //     std::cout << "init, col size, " << i << "th, " << colinfo[i].size() << std::endl;
    // }

    int vnbr = V.rows();
    if (OptType == 0)// AAG
    {
        // vertices vnbr * 3, normals vnbr * 3, binormals for G vnbr * 3  
        varsize = vnbr * 9;
    }
    if (OptType == 1)// AGG
    {
        // vertices vnbr * 3, normals vnbr * 3, binormals for G0 vnbr * 3, binormals for G1 vnbr * 3.    
        varsize = vnbr * 12;
    }
    if (OptType == 2)// PP
    {
        // vertices vnbr * 3, normals vnbr * 3, binormals + side vectors:  (vnbr * 6) * 2
        varsize = vnbr * 18;
   
    }
    if (OptType == 3)// PPG
    {
        // vertices vnbr * 3, normals vnbr * 3, binormals + side vectors:  (vnbr * 6) * 2 + vnbr * 3
        varsize = vnbr * 21;

    }
    if (OptType == 4)// AAGG
    {
        // vertices vnbr * 3, normals vnbr * 3, binormals for G vnbr * 3 * 2  
        varsize = vnbr * 12;

    }

    getQuadRowColsRegular(vnbr, vNbrInRow, row_front, row_back, col_front, col_back);
    std::cout << "Row Col info computed:\n"<<vnbr<<" vertices, each row has "<<vNbrInRow<<" vertices\n";
    load_diagonal_info(row_front, row_back, col_front, col_back, d0_front, d1_front, d0_back, d1_back);
    OrigVars = Eigen::VectorXd::Zero(varsize);
    OrigVars.segment(0, vnbr) = V.col(0);
    OrigVars.segment(vnbr, vnbr) = V.col(1);
    OrigVars.segment(vnbr * 2, vnbr) = V.col(2);
    mesh_original = mesh_in;
    mesh_update = mesh_original;
    ComputeAuxiliaries = true;
    Estimate_PG_Angles = true;
    Eigen::VectorXd grav_vec = Eigen::VectorXd::Ones(vnbr * 3);
    Eigen::VectorXd var_vec = Eigen::VectorXd::Zero(varsize);
    var_vec.segment(0, vnbr * 3) = grav_vec;
    gravity_matrix = var_vec.asDiagonal();
    Eigen::Vector3d vmin(V.col(0).minCoeff(), V.col(1).minCoeff(), V.col(2).minCoeff());
    Eigen::Vector3d vmax(V.col(0).maxCoeff(), V.col(1).maxCoeff(), V.col(2).maxCoeff());
    std::cout<<"Initialized quad mesh: Vnbr, "<<V.rows()<<", BBD, "<<(vmin - vmax).norm()<<std::endl;;
}
void QuadOpt::initISO(const std::vector<Eigen::Vector3d> &Vlist, const int rnbr)
{
	int vnbr = Vlist.size();
	V = vec_list_to_matrix(Vlist);
	// std::cout<<"check V, \n"<<V<<"\n";
	constructRegularF(vnbr, rnbr, F);
	MeshProcessing mp;
	mp.matrix2Mesh(mesh_original, V, F);
	vNbrInRow = rnbr;
	if (OriginalCurve.size() == 0)
	{
		for (int i = 0; i < rnbr * 2; i++) // GGG we preserve 2 strips
		{
			OriginalCurve.push_back(Vlist[i]);
		}
	}
	int fnbr = F.rows();
	varsize = vnbr * 3;
	d0Lengths.resize(fnbr);
	d1Lengths.resize(fnbr);
	d0d1dots.resize(fnbr);
	for (int i = 0; i < fnbr; i++)
	{
		Eigen::Vector3d v0 = V.row(F(i, 0));
		Eigen::Vector3d v1 = V.row(F(i, 1));
		Eigen::Vector3d v2 = V.row(F(i, 2));
		Eigen::Vector3d v3 = V.row(F(i, 3));
		d0Lengths[i] = (v0 - v2).norm();
		d1Lengths[i] = (v1 - v3).norm();
		d0d1dots[i] = (v0 - v2).dot(v1 - v3);
	}
	if (vNbrInRow <= 0)
	{
		std::cout << "Please set the parameter vNbrInRow!\n";
		return;
	}
	std::cout << "check 1\n";
	getQuadRowColsRegular(vnbr, vNbrInRow, row_front, row_back, col_front, col_back);
	std::cout << "Row Col info computed:\n" << vnbr << " vertices, each row has " << vNbrInRow << " vertices\n";
	load_diagonal_info(row_front, row_back, col_front, col_back, d0_front, d1_front, d0_back, d1_back);
	OrigVars = Eigen::VectorXd::Zero(varsize);
	for (int i = 0; i < vnbr; i++)
	{
		OrigVars[3 * i] = V(i, 0);
		OrigVars[3 * i + 1] = V(i, 1);
		OrigVars[3 * i + 2] = V(i, 2);
	}
	std::cout << "check 2\n";
	mesh_update = mesh_original;
	ComputeAuxiliaries = true;
	Eigen::VectorXd grav_vec = Eigen::VectorXd::Ones(vnbr * 3);
	Eigen::VectorXd var_vec = Eigen::VectorXd::Zero(varsize);
	var_vec.segment(0, vnbr * 3) = grav_vec;
	gravity_matrix = var_vec.asDiagonal();
	Eigen::Vector3d vmin(V.col(0).minCoeff(), V.col(1).minCoeff(), V.col(2).minCoeff());
	Eigen::Vector3d vmax(V.col(0).maxCoeff(), V.col(1).maxCoeff(), V.col(2).maxCoeff());
	std::cout << "Initialized quad mesh: Vnbr, " << V.rows() << ", BBD, " << (vmin - vmax).norm() << std::endl;;
}
// the gravity vector show where the vertices are
// type == 0: aag.
// type == 1: agg.
// type == 2: ggg
void assignAagAggVertices(Eigen::VectorXd &vars,Eigen::VectorXd &graVec, const int varsize, const std::vector<Eigen::Vector3d> &vers, int type)
{
    int varForVer = 9;
    if (type == 1)
    {
        varForVer = 12;
    }
	if (type == 2)
	{
		varForVer = 15;
	}
    vars = Eigen::VectorXd::Zero(varsize);
    assert(varsize == vers.size() * varForVer);
    graVec = vars;
    for (int i = 0; i < vers.size(); i++)
    {
        vars[i * varForVer] = vers[i][0];
        vars[i * varForVer + 1] = vers[i][1];
        vars[i * varForVer + 2] = vers[i][2];
        graVec[i * varForVer] = 1;
        graVec[i * varForVer + 1] = 1;
        graVec[i * varForVer + 2] = 1;
    }
}

void constructRegularF(const int vnbr, const int rnbr, Eigen::MatrixXi &F)
{
    
    int cnbr = vnbr / rnbr;
    F.resize((cnbr - 1) * (rnbr - 1), 4);
    int fnbr = 0;
    for (int i = 0; i < cnbr - 1; i++)
    {
        for (int j = 0; j < rnbr - 1; j++)
        {
            int v0 = i * rnbr + j;
            int v1 = i * rnbr + j + 1;
            int v2 = (i + 1) * rnbr + j + 1;
            int v3 = (i + 1) * rnbr + j;
            F(fnbr, 0) = v0;
            F(fnbr, 1) = v1;
            F(fnbr, 2) = v2;
            F(fnbr, 3) = v3;
            fnbr++;
        }
    }
    // std::cout<<"check F, \n"<<F<<"\n";
}

void QuadOpt::initAAG(const std::vector<Eigen::Vector3d> &Vlist, const int rnbr)
{

    int vnbr = Vlist.size();
    V = vec_list_to_matrix(Vlist);
    // std::cout<<"check V, \n"<<V<<"\n";
    constructRegularF(vnbr, rnbr, F);
    MeshProcessing mp;
    mp.matrix2Mesh(mesh_original, V, F);
    vNbrInRow = rnbr;
    if (OriginalCurve.size() == 0)
    {
        for (int i = 0; i < rnbr * 2; i++)
        {
            OriginalCurve.push_back(Vlist[i]);
        }
    }

    // vertices vnbr * 3, normals vnbr * 3, binormals for G vnbr * 3
    varsize = vnbr * 9;
    // row_front and row_back is the geodesic direction.
    getQuadRowColsRegular(vnbr, rnbr, row_front, row_back, col_front, col_back);
    std::cout << "Row Col info computed\n";
    load_diagonal_info(row_front, row_back, col_front, col_back, d0_front, d1_front, d0_back, d1_back);
    std::cout << "Topo obtained\n";
    Eigen::VectorXd var_vec;
    assignAagAggVertices(OrigVars, var_vec, varsize, Vlist, 0);
    std::cout << "variables assigned\n";
    // mesh_original = mesh_in;
    mesh_update = mesh_original;
    verOriginal = Vlist;
    verUpdate = verOriginal;
    ComputeAuxiliaries = true;
    if (WhichDiagonal == 0)
    {
        d0_type = 1;
    }
    if (WhichDiagonal == 1)
    {
        d1_type = 1;
    }
    Estimate_PG_Angles = false;
    gravity_matrix = var_vec.asDiagonal();
    Eigen::Vector3d vmin(V.col(0).minCoeff(), V.col(1).minCoeff(), V.col(2).minCoeff());
    Eigen::Vector3d vmax(V.col(0).maxCoeff(), V.col(1).maxCoeff(), V.col(2).maxCoeff());
    std::cout<<"Initialized quad mesh: Vnbr, "<<V.rows()<<", BBD, "<<(vmin - vmax).norm()<<std::endl;;
}

void QuadOpt::initAGG(const std::vector<Eigen::Vector3d> &Vlist, const int rnbr)
{

    int vnbr = Vlist.size();
    V = vec_list_to_matrix(Vlist);
    // std::cout<<"check V, \n"<<V<<"\n";
    constructRegularF(vnbr, rnbr, F);
    MeshProcessing mp;
    mp.matrix2Mesh(mesh_original, V, F);
    vNbrInRow = rnbr;
    if (OriginalCurve.size() == 0)
    {
        for (int i = 0; i < rnbr * 2; i++) // AGG we preserve 2 strips
        {
            OriginalCurve.push_back(Vlist[i]);
        }
    }

    // vertices vnbr * 3, normals vnbr * 3, binormals for G1 vnbr * 3, binormals for G2 vnbr * 3
    varsize = vnbr * 12;
    // row_front and row_back is the geodesic direction.
    getQuadRowColsRegular(vnbr, rnbr, row_front, row_back, col_front, col_back);
    std::cout << "Row Col info computed\n";
    load_diagonal_info(row_front, row_back, col_front, col_back, d0_front, d1_front, d0_back, d1_back);
    std::cout << "Topo obtained\n";
    Eigen::VectorXd var_vec;
    assignAagAggVertices(OrigVars, var_vec, varsize, Vlist, 1);
    std::cout << "variables assigned\n";
    // mesh_original = mesh_in;
    mesh_update = mesh_original;
    verOriginal = Vlist;
    verUpdate = verOriginal;
    ComputeAuxiliaries = true;
    if (WhichDiagonal == 0)
    {
        d0_type = 1;
    }
    if (WhichDiagonal == 1)
    {
        d1_type = 1;
    }
    Estimate_PG_Angles = false;
    gravity_matrix = var_vec.asDiagonal();
    Eigen::Vector3d vmin(V.col(0).minCoeff(), V.col(1).minCoeff(), V.col(2).minCoeff());
    Eigen::Vector3d vmax(V.col(0).maxCoeff(), V.col(1).maxCoeff(), V.col(2).maxCoeff());
    std::cout<<"Initialized quad mesh: Vnbr, "<<V.rows()<<", BBD, "<<(vmin - vmax).norm()<<std::endl;;
}

void QuadOpt::initGGG(const std::vector<Eigen::Vector3d> &Vlist, const int rnbr)
{

	int vnbr = Vlist.size();
	V = vec_list_to_matrix(Vlist);
	// std::cout<<"check V, \n"<<V<<"\n";
	constructRegularF(vnbr, rnbr, F);
	MeshProcessing mp;
	mp.matrix2Mesh(mesh_original, V, F);
	vNbrInRow = rnbr;
	if (OriginalCurve.size() == 0)
	{
		for (int i = 0; i < rnbr * 2; i++) // GGG we preserve 2 strips
		{
			OriginalCurve.push_back(Vlist[i]);
		}
	}

	// vertices vnbr * 3, normals vnbr * 3, binormals for G1 vnbr * 3, binormals for G2 vnbr * 3, binormals for G3 vnbr * 3
	varsize = vnbr * 15;
	// row_front and row_back is the geodesic direction.
	getQuadRowColsRegular(vnbr, rnbr, row_front, row_back, col_front, col_back);
	std::cout << "Row Col info computed\n";
	load_diagonal_info(row_front, row_back, col_front, col_back, d0_front, d1_front, d0_back, d1_back);
	std::cout << "Topo obtained\n";
	Eigen::VectorXd var_vec;
	assignAagAggVertices(OrigVars, var_vec, varsize, Vlist, 2);
	std::cout << "variables assigned\n";
	// mesh_original = mesh_in;
	mesh_update = mesh_original;
	verOriginal = Vlist;
	verUpdate = verOriginal;
	ComputeAuxiliaries = true;
	if (WhichDiagonal == 0)
	{
		d0_type = 1;
	}
	if (WhichDiagonal == 1)
	{
		d1_type = 1;
	}
	Estimate_PG_Angles = false;
	gravity_matrix = var_vec.asDiagonal();
	Eigen::Vector3d vmin(V.col(0).minCoeff(), V.col(1).minCoeff(), V.col(2).minCoeff());
	Eigen::Vector3d vmax(V.col(0).maxCoeff(), V.col(1).maxCoeff(), V.col(2).maxCoeff());
	std::cout << "Initialized quad mesh: Vnbr, " << V.rows() << ", BBD, " << (vmin - vmax).norm() << std::endl;;
}
void QuadOpt::reset()
{
    // re-choose a diagonal
    if (OptType == 0) // AAG
    {                                  // set up G
        std::cout << "Set Up AAG Diagonals" << std::endl;
        if (WhichDiagonal == 0)
        {
            d0_type = 2;
            d1_type = 0;
        }
        if (WhichDiagonal == 1)
        { // set up G
            d0_type = 0;
            d1_type = 2;
        }
    }
    if (OptType == 1) // GGA
    {                                  // set up A
        std::cout << "Set Up GGA Diagonals" << std::endl;
        if (WhichDiagonal == 0)
        {
            d0_type = 1;
            d1_type = 0;
        }
        if (WhichDiagonal == 1)
        { // set up A
            d0_type = 0;
            d1_type = 1;
        }
    }
    if (OptType == 2) // PPG
    {

        d0_type = 0;
        d1_type = 0;
    }
    if (OptType == 3) // PPG
    {                         
        std::cout << "Set Up GGA Diagonals" << std::endl;
        if (WhichDiagonal == 0)
        {
            d0_type = 2;
            d1_type = 0;
        }
        if (WhichDiagonal == 1)
        { 
            d0_type = 0;
            d1_type = 2;
        }
    }
    if (OptType == 4) // AAGG
    {                         
        std::cout << "Set Up AAGG Diagonals" << std::endl;
        if (WhichDiagonal == 0)
        {
            d0_type = 2;
            d1_type = 2;
        }
        if (WhichDiagonal == 1)
        { 
            d0_type = 2;
            d1_type = 2;
        }
    }
    MP.mesh2Matrix(mesh_original, V, F);
    std::cout<<"Get V and F of the quads: "<<V.rows()<<", "<<F.rows()<<std::endl;
    int vnbr = V.rows();
    // reset the var size
    if (OptType == 0)// AAG
    {
        // vertices vnbr * 3, normals vnbr * 3, binormals for G vnbr * 3  
        varsize = vnbr * 9;
    }
    if (OptType == 1)// AGG
    {
        // vertices vnbr * 3, normals vnbr * 3, binormals for G0 vnbr * 3, binormals for G1 vnbr * 3.    
        varsize = vnbr * 12;
    }
    if (OptType == 2)// PP
    {
        // vertices vnbr * 3, normals vnbr * 3, binormals + side vectors:  (vnbr * 6) * 2
        varsize = vnbr * 18;
   
    }
    if (OptType == 3)// PPG
    {
        // vertices vnbr * 3, normals vnbr * 3, binormals + side vectors:  (vnbr * 6) * 2 + vnbr * 3
        varsize = vnbr * 21;

    }
    if (OptType == 4)// AAGG
    {
        // vertices vnbr * 3, normals vnbr * 3, binormals for G0 vnbr * 3, binormals for G1 vnbr * 3.    
        varsize = vnbr * 12;
    }

    OrigVars = Eigen::VectorXd::Zero(varsize);
    OrigVars.segment(0, vnbr) = V.col(0);
    OrigVars.segment(vnbr, vnbr) = V.col(1);
    OrigVars.segment(vnbr * 2, vnbr) = V.col(2);
    GlobVars.resize(0);
    ComputeAuxiliaries = true;
    Estimate_PG_Angles = true;
    mesh_update = mesh_original;
    Eigen::VectorXd grav_vec = Eigen::VectorXd::Ones(vnbr * 3);
    Eigen::VectorXd var_vec = Eigen::VectorXd::Zero(varsize);
    var_vec.segment(0, vnbr * 3) = grav_vec;
    gravity_matrix = var_vec.asDiagonal();
}

void QuadOpt::assemble_gravity(spMat& H, Eigen::VectorXd& B, Eigen::VectorXd &energy){

    int vnbr = V.rows();
    std::vector<Eigen::Vector3d> vprojs(vnbr);
    std::vector<Eigen::Vector3d> Nlocal(vnbr);
    for (int i = 0; i < vnbr; i++)
    {
        Eigen::Vector3d query = V.row(i);
        int f;
        Eigen::RowVector3d proj;
        aabbtree.squared_distance(Vtri, Ftri, query, f, proj);
        Eigen::Vector3d p = proj;
        Eigen::Vector3d v0 = Vtri.row(Ftri(f, 0));
        Eigen::Vector3d v1 = Vtri.row(Ftri(f, 1));
        Eigen::Vector3d v2 = Vtri.row(Ftri(f, 2));
        std::array<double, 3> coor = 
        barycenter_coordinate(v0, v1, v2, p);
        Eigen::Vector3d norm = coor[0] * Ntri.row(Ftri(f, 0)) + coor[1] * Ntri.row(Ftri(f, 1)) + coor[2] * Ntri.row(Ftri(f, 2));
        assert(Ntri.row(Ftri(f, 0)).dot(Ntri.row(Ftri(f, 1))) > 0 && Ntri.row(Ftri(f, 0)).dot(Ntri.row(Ftri(f, 2))) > 0);
        norm.normalize();
        vprojs[i] = p;
        Nlocal[i] = norm;

        // debug
        Eigen::Vector3d rec_pro = coor[0] * v0 + coor[1] * v1 + coor[2] * v2;
        if ((rec_pro - p).norm() > 0.5)
        {
            std::cout<<"INACCURATE SOLVING: coor: "<<coor[0]<<", "<<coor[1]<<", "<<coor[2]<<std::endl;
            std::cout<<"v0, "<<v0.transpose()<<std::endl;
            std::cout<<"v1, "<<v1.transpose()<<std::endl;
            std::cout<<"v2, "<<v2.transpose()<<std::endl;
            std::cout<<"p "<<p.transpose()<<std::endl;
            std::cout<<"p_compute "<<rec_pro.transpose()<<std::endl;
        }
    }
    // Ppro0 = vec_list_to_matrix(vprojs);
    // Npro0 = vec_list_to_matrix(Nlocal);
    std::vector<Trip> tripletes;
    tripletes.reserve(vnbr * 9);
    energy = Eigen::VectorXd::Zero(vnbr * 3);
    for (int i = 0; i < vnbr; i++)
    {
        int vid = i;
        int lx = vid;
        int ly = vid + vnbr;
        int lz = vid + vnbr * 2;
        Eigen::Vector3d normal = Nlocal[vid]; // the normal vector
        Eigen::Vector3d closest = vprojs[vid]; // the closest point
        Eigen::Vector3d ver = V.row(vid);

        // ver.dot(n^*) - ver^*.dot(n^*) = 0
        tripletes.push_back(Trip(i, lx, normal[0]));
        tripletes.push_back(Trip(i, ly, normal[1]));
        tripletes.push_back(Trip(i, lz, normal[2]));

        energy[i] = ver.dot(normal) - closest.dot(normal);

        // vertices not far away from the original ones
        // (ver - ver^*)^2 = 0
        double scale = 1;
        Eigen::Vector3d verori;
        verori = closest;
        // Eigen::Vector3d normal_local = normal;
        // verori[0] = OrigVars[lx];
        // verori[1] = OrigVars[ly];
        // verori[2] = OrigVars[lz];
        Eigen::Vector3d vdiff = Eigen::Vector3d(V.row(vid)) - verori;
        // if (vdiff.dot(normal_local) < 0)
        // {
        //     normal_local *= -1;
        // }
        // double cos_angle = vdiff.normalized().dot(normal_local);
        // double angle_radian = acos(cos_angle);
        // if (angle_radian > 30)
        // {
        //     continue;
        // }
        // std::cout<<"Through here"<<std::endl;
        tripletes.push_back(Trip(i + vnbr, lx, 2 * vdiff[0] * scale));
        tripletes.push_back(Trip(i + vnbr, ly, 2 * vdiff[1] * scale));
        tripletes.push_back(Trip(i + vnbr, lz, 2 * vdiff[2] * scale));

        energy[i + vnbr] = vdiff.dot(vdiff) * scale;

        verori[0] = OrigVars[lx];
        verori[1] = OrigVars[ly];
        verori[2] = OrigVars[lz];
        scale = 0.1;
        vdiff = Eigen::Vector3d(V.row(vid)) - verori;

        tripletes.push_back(Trip(i + vnbr * 2, lx, 2 * vdiff[0] * scale));
        tripletes.push_back(Trip(i + vnbr * 2, ly, 2 * vdiff[1] * scale));
        tripletes.push_back(Trip(i + vnbr * 2, lz, 2 * vdiff[2] * scale));

        energy[i + vnbr * 2] = vdiff.dot(vdiff) * scale;
    }
    int nvars = GlobVars.size();
    int ncondi = energy.size();
    spMat J;
    J.resize(ncondi, nvars);
    J.setFromTriplets(tripletes.begin(), tripletes.end());
    H = J.transpose() * J;
    B = -J.transpose() * energy;
}
// this function approximates the first point in each row to the given curve.
// this is only used for AAG and AGG
// if the reference curve exists, enable the approximation
void QuadOpt::assemble_approximate_curve_conditions(spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy, const int order)
{
    int varOffset = 9;
    if(order == 2){
        varOffset = 12;
    }
	if (order == 3)
	{
		varOffset = 15;
	}
    int vnbr = V.rows();
    int rnbr = vnbr / vNbrInRow; // the number of rows
    std::vector<Trip> tripletes;
    tripletes.reserve(rnbr * 6);
    // std::cout<<"in the assemble_approximate_curve_conditions\n";
    energy = Eigen::VectorXd::Zero(rnbr * 4);
    if (curveRef.size() > 0)
        for (int i = 0; i < rnbr; i++)
        {
            int vid = i * vNbrInRow; // the first point in each row.
            int lx = vid * varOffset;
            int ly = vid * varOffset + 1;
            int lz = vid * varOffset + 2;
            // std::cout<<"vid "<<vid<<"\n";
            Eigen::Vector3d ver = V.row(vid);
            int segid;
            double tlocal;
            Eigen::Vector3d plocal;
            Eigen::Vector3d tangent;
            // find the projection point of this point on the curve.
            bool project = true;
            if (i == 0)
            {
                plocal = curveRef[0];
            }
            else
                project = projectPointOnCurve(curveRef, ver, segid, tlocal, plocal, tangent);

            // std::cout<<"tangent "<<tangent.transpose()<<<<""<<"\n";
            if (!project)
            {
                continue;
            }

            // vertices not far away from the original ones
            // (ver - ver^*)^2 = 0
            double scale = 1;
            Eigen::Vector3d verori = plocal;
            Eigen::Vector3d vdiff;
            vdiff = Eigen::Vector3d(V.row(vid)) - verori;

            tripletes.push_back(Trip(i, lx, 2 * vdiff[0] * scale));
            tripletes.push_back(Trip(i, ly, 2 * vdiff[1] * scale));
            tripletes.push_back(Trip(i, lz, 2 * vdiff[2] * scale));

            energy[i] = vdiff.dot(vdiff) * scale;
            if (i == 0)
            {
                continue;
            }

            // ver - plocal - tangent * ||ver - plocal|| = 0: the vertex is moving along the tangent directions
            scale = (ver - plocal).norm();
            int sign = (ver - plocal).dot(tangent) > 0 ? -1 : 1;

            tripletes.push_back(Trip(i + rnbr, lx, 1));
            energy[i + rnbr] = ((ver - plocal) + sign * tangent * scale)[0];

            tripletes.push_back(Trip(i + rnbr * 2, ly, 1));
            energy[i + rnbr * 2] = ((ver - plocal) + sign * tangent * scale)[1];

            tripletes.push_back(Trip(i + rnbr * 3, lz, 1));
            energy[i + rnbr * 3] = ((ver - plocal) + sign * tangent * scale)[2];
        }
    int nvars = GlobVars.size();
    int ncondi = energy.size();
    spMat J;
    J.resize(ncondi, nvars);
    J.setFromTriplets(tripletes.begin(), tripletes.end());
    H = J.transpose() * J;
    B = -J.transpose() * energy;
    // std::cout<<"out of the assemble_approximate_curve_conditions\n";
}
void QuadOpt::assemble_isometric_condition(spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy)
{
	int vnbr = V.rows();
	int fnbr = F.rows();
	std::vector<Trip> tripletes;
	tripletes.reserve(24 * fnbr);
	energy = Eigen::VectorXd::Zero(fnbr * 3);
	for (int i = 0; i < fnbr; i++)
	{
		int v0 = F(i, 0);
		int v1 = F(i, 1);
		int v2 = F(i, 2);
		int v3 = F(i, 3);

		int l0x = v0 * 3;
		int l0y = v0 * 3 + 1;
		int l0z = v0 * 3 + 2;

		int l1x = v1 * 3;
		int l1y = v1 * 3 + 1;
		int l1z = v1 * 3 + 2;

		int l2x = v2 * 3;
		int l2y = v2 * 3 + 1;
		int l2z = v2 * 3 + 2;

		int l3x = v3 * 3;
		int l3y = v3 * 3 + 1;
		int l3z = v3 * 3 + 2;
		
		// (v0 - v2)^2 - l^2 = 0
		tripletes.push_back(Trip(3 * i, l0x, 2 * V(v0, 0) - 2 * V(v2, 0)));
		tripletes.push_back(Trip(3 * i, l0y, 2 * V(v0, 1) - 2 * V(v2, 1)));
		tripletes.push_back(Trip(3 * i, l0z, 2 * V(v0, 2) - 2 * V(v2, 2)));

		tripletes.push_back(Trip(3 * i, l2x, 2 * V(v2, 0) - 2 * V(v0, 0)));
		tripletes.push_back(Trip(3 * i, l2y, 2 * V(v2, 1) - 2 * V(v0, 1)));
		tripletes.push_back(Trip(3 * i, l2z, 2 * V(v2, 2) - 2 * V(v0, 2)));
		
		double length = Eigen::Vector3d(V.row(v0) - V.row(v2)).norm();
		energy[3 * i] = length * length - d0Lengths[i] * d0Lengths[i];

		// (v1 - v3)^2 - l^2 = 0
		tripletes.push_back(Trip(3 * i + 1, l1x, 2 * V(v1, 0) - 2 * V(v3, 0)));
		tripletes.push_back(Trip(3 * i + 1, l1y, 2 * V(v1, 1) - 2 * V(v3, 1)));
		tripletes.push_back(Trip(3 * i + 1, l1z, 2 * V(v1, 2) - 2 * V(v3, 2)));

		tripletes.push_back(Trip(3 * i + 1, l3x, 2 * V(v3, 0) - 2 * V(v1, 0)));
		tripletes.push_back(Trip(3 * i + 1, l3y, 2 * V(v3, 1) - 2 * V(v1, 1)));
		tripletes.push_back(Trip(3 * i + 1, l3z, 2 * V(v3, 2) - 2 * V(v1, 2)));

		length = Eigen::Vector3d(V.row(v1) - V.row(v3)).norm();
		energy[3 * i + 1] = length * length - d1Lengths[i] * d1Lengths[i];

		// (v0 - v2) * (v1 - v3) - cosin(angle) = 0
		Eigen::Vector3d V02 = V.row(v0) - V.row(v2);
		Eigen::Vector3d V13 = V.row(v1) - V.row(v3);

		tripletes.push_back(Trip(3 * i + 2, l0x, V13[0]));
		tripletes.push_back(Trip(3 * i + 2, l0y, V13[1]));
		tripletes.push_back(Trip(3 * i + 2, l0z, V13[2]));

		tripletes.push_back(Trip(3 * i + 2, l2x, -V13[0]));
		tripletes.push_back(Trip(3 * i + 2, l2y, -V13[1]));
		tripletes.push_back(Trip(3 * i + 2, l2z, -V13[2]));

		tripletes.push_back(Trip(3 * i + 2, l1x, V02[0]));
		tripletes.push_back(Trip(3 * i + 2, l1y, V02[1]));
		tripletes.push_back(Trip(3 * i + 2, l1z, V02[2]));

		tripletes.push_back(Trip(3 * i + 2, l3x, -V02[0]));
		tripletes.push_back(Trip(3 * i + 2, l3y, -V02[1]));
		tripletes.push_back(Trip(3 * i + 2, l3z, -V02[2]));

		energy[3 * i + 2] = V02.dot(V13) - d0d1dots[0];
	}
	int nvars = GlobVars.size();
	int ncondi = energy.size();
	spMat J;
	J.resize(ncondi, nvars);
	J.setFromTriplets(tripletes.begin(), tripletes.end());
	H = J.transpose() * J;
	B = -J.transpose() * energy;
}

void QuadOpt::assemble_gravity_AAG_AGG(spMat& H, Eigen::VectorXd& B, Eigen::VectorXd &energy, int order){
    int varOffset = 9;
    if(order == 2){
        varOffset = 12;
    }
	if (order == 3)
	{
		varOffset = 15;
	}
	if (order == 4)
	{
		varOffset = 3;
	}
    std::vector<Trip> tripletes;
    int vnbr = V.rows();
    tripletes.reserve(vnbr * 3);
    energy = Eigen::VectorXd::Zero(vnbr);
    for (int i = 0; i < vnbr; i++)
    {
        int vid = i;
        int lx = vid * varOffset;
        int ly = vid * varOffset + 1;
        int lz = vid * varOffset + 2;
        Eigen::Vector3d ver = V.row(vid);

        // vertices not far away from the original ones
        // (ver - ver^*)^2 = 0
        double scale = 1;
        // if (i >= vNbrInRow)
        // {
        //     scale = 0.1; // weak approximate to the second curve
        // }
        Eigen::Vector3d verori;
        Eigen::Vector3d vdiff;
        int nbrOfAppro = AGGAAGApproInitStrip? vNbrInRow * 2: vNbrInRow; // if approximate to the first strip, need the first 2 rows of mesh.
        if (i < nbrOfAppro)
        {
            verori = OriginalCurve[i];
        }
        else
        {
			if (GggEnableEditing) {
				if (i == GggEditingVid)
				{
					verori = GggTargetPosition;
				}
				else
				{
					continue;
				}
			}
			else
				continue;
		}

        vdiff = Eigen::Vector3d(V.row(vid)) - verori;

        tripletes.push_back(Trip(i, lx, 2 * vdiff[0] * scale));
        tripletes.push_back(Trip(i, ly, 2 * vdiff[1] * scale));
        tripletes.push_back(Trip(i, lz, 2 * vdiff[2] * scale));

        energy[i] = vdiff.dot(vdiff) * scale;
    }
    int nvars = GlobVars.size();
    int ncondi = energy.size();
    spMat J;
    J.resize(ncondi, nvars);
    J.setFromTriplets(tripletes.begin(), tripletes.end());
    H = J.transpose() * J;
    B = -J.transpose() * energy;
}
void QuadOpt::assemble_gravity_ApproOriginal(spMat& H, Eigen::VectorXd& B, Eigen::VectorXd &energy)
{
    
    std::vector<Trip> tripletes;
    int vnbr = V.rows();
    tripletes.reserve(vnbr * 3);
    energy = Eigen::VectorXd::Zero(vnbr);
    for (int i = 0; i < vnbr; i++)
    {
        int vid = i;
        int lx = vid;
        int ly = vid + vnbr;
        int lz = vid + vnbr * 2;
        Eigen::Vector3d ver = V.row(vid);

        // vertices not far away from the original ones
        // (ver - ver^*)^2 = 0
        double scale = 1;
        // if (i >= vNbrInRow)
        // {
        //     scale = 0.1; // weak approximate to the second curve
        // }
        Eigen::Vector3d verori;
        verori[0] = OrigVars[lx];
        verori[1] = OrigVars[ly];
        verori[2] = OrigVars[lz];
        Eigen::Vector3d vdiff;

        vdiff = Eigen::Vector3d(V.row(vid)) - verori;

        tripletes.push_back(Trip(i, lx, 2 * vdiff[0] * scale));
        tripletes.push_back(Trip(i, ly, 2 * vdiff[1] * scale));
        tripletes.push_back(Trip(i, lz, 2 * vdiff[2] * scale));

        energy[i] = vdiff.dot(vdiff) * scale;
    }
    int nvars = GlobVars.size();
    int ncondi = energy.size();
    spMat J;
    J.resize(ncondi, nvars);
    J.setFromTriplets(tripletes.begin(), tripletes.end());
    H = J.transpose() * J;
    B = -J.transpose() * energy;
}
void QuadOpt::load_triangle_mesh_tree(const igl::AABB<Eigen::MatrixXd, 3> &tree, const Eigen::MatrixXd &Vt,
                                      const Eigen::MatrixXi &Ft, const Eigen::MatrixXd &Nt)
{
    aabbtree = tree;
    Vtri = Vt;
    Ftri = Ft;
    Ntri = Nt;
    std::cout<<"Tree got loaded from triangle mesh."<<std::endl;
}

// given a vertex id, return the x y z locations of a vertex in variables
void locatorRules(int vid, int vnbr, int order, int &x, int &y, int &z)
{
    if (order == 0)
    {
        x = vid;
        y = vid + vnbr;
        z = vid + vnbr * 2;
		return;
    }
    if (order == 1)
    {
        x = 9 * vid;
        y = 9 * vid + 1;
        z = 9 * vid + 2;
    }
    if (order == 2)
    {
        x = 12 * vid;
        y = 12 * vid + 1;
        z = 12 * vid + 2;
    }
	if (order == 3)
	{
		x = 15 * vid;
		y = 15 * vid + 1;
		z = 15 * vid + 2;
	}
	if (order == 4)
	{
		x = 3 * vid;
		y = 3 * vid + 1;
		z = 3 * vid + 2;
	}
}
// lv is ver location, lf is front location, lb is back location, cid is the id of the condition
// (v-f)/scale0 + (v-b)/scale1 = 0
void push_fairness_conditions(std::vector<Trip> &tripletes, Eigen::VectorXd &energy,
                              const int lv, const int lf, const int lb,
                              const double vval, const double fval, const double bval,
                              const int cid, const double scale0, const double scale1, const int order = 0)
{
    // tripletes.push_back(Trip(cid, lv, 1 / scale0 + 1 / scale1));
    // tripletes.push_back(Trip(cid, lf, -1 / scale0));
    // tripletes.push_back(Trip(cid, lb, -1 / scale1));
    // energy[cid] = (vval - fval) / scale0 + (vval - bval) / scale1;
    if (order == 0)
    {
        tripletes.push_back(Trip(cid, lv, 1 / scale0 + 1 / scale1));
        tripletes.push_back(Trip(cid, lf, -1 / scale0));
        tripletes.push_back(Trip(cid, lb, -1 / scale1));
        energy[cid] = (vval - fval) / scale0 + (vval - bval) / scale1;
    }
    else
    {
        tripletes.push_back(Trip(cid, lv, 2));
        tripletes.push_back(Trip(cid, lf, -1));
        tripletes.push_back(Trip(cid, lb, -1));
        energy[cid] = (vval - fval) + (vval - bval);
    }
}

bool validateOrder1(const int vnbr, const int location)
{
    int mn = location % 9;
    if (mn == 0 || mn == 1 || mn == 2)
    {
        return true;
    }
    return false;
}
void QuadOpt::assemble_fairness(spMat& H, Eigen::VectorXd& B, Eigen::VectorXd &energy, const int order){
    

    std::vector<Trip> tripletes;
    int vnbr = V.rows();
    energy = Eigen::VectorXd::Zero(vnbr * 12);//each family takes vnbr * 3 constraints, there are at most 4 families
    tripletes.reserve(vnbr * 12 * 3);
    int counter = 0;
    for (int i = 0; i < vnbr; i++)
    {
        // the ids
        int vid = i;
        int rf = row_front[vid];
        int rb = row_back[vid];
        int cf = col_front[vid];
        int cb = col_back[vid];
        int d0f = d0_front[vid];
        int d0b = d0_back[vid];
        int d1f = d1_front[vid];
        int d1b = d1_back[vid];
        // the locations 
        int lvx;
        int lvy;
        int lvz;
        locatorRules(vid, vnbr, order, lvx, lvy, lvz);
        int lrfx;
        int lrfy;
        int lrfz;
        locatorRules(rf, vnbr, order, lrfx, lrfy, lrfz);

        int lrbx;
        int lrby;
        int lrbz;
        locatorRules(rb, vnbr, order, lrbx, lrby, lrbz);

        int lcfx;
        int lcfy;
        int lcfz;
        locatorRules(cf, vnbr, order, lcfx, lcfy, lcfz);

        int lcbx;
        int lcby;
        int lcbz;
        locatorRules(cb, vnbr, order, lcbx, lcby, lcbz);

        int ld0fx;
        int ld0fy;
        int ld0fz;
        locatorRules(d0f, vnbr, order, ld0fx, ld0fy, ld0fz);

        int ld0bx;
        int ld0by;
        int ld0bz;
        locatorRules(d0b, vnbr, order, ld0bx, ld0by, ld0bz);

        int ld1fx;
        int ld1fy;
        int ld1fz;
        locatorRules(d1f, vnbr, order, ld1fx, ld1fy, ld1fz);

        int ld1bx;
        int ld1by;
        int ld1bz;
        locatorRules(d1b, vnbr, order, ld1bx, ld1by, ld1bz);
        
        bool row_smt = true;
        bool col_smt = true;
        bool d0_smt = true;
        bool d1_smt = true;
        // if (order == 2 && i < vNbrInRow * 2) // cancel the smoothness for the first 2 rows
        // {
        //     row_smt = false;
        // }

        if (rf < 0 || rb < 0)
        {
            row_smt = false;
        }
        if (cf < 0 || cb < 0)
        {
            col_smt = false;
        }

        if (d0_type == 0 || d0f < 0 || d0b < 0 || OptType == 2)
        {
            d0_smt = false;
        }
        if (d1_type == 0 || d1f < 0 || d1b < 0 || OptType == 2)
        {
            d1_smt = false;
        }
        if (order == 1)
        {
            // std::cout<<"Smoothness checking ... \n";
            if (row_smt)
            {
                assert(validateOrder1(vnbr, lvx) && validateOrder1(vnbr, lvy) && validateOrder1(vnbr, lvz));
                assert(validateOrder1(vnbr, lrfx) && validateOrder1(vnbr, lrfy) && validateOrder1(vnbr, lrfz));
                assert(validateOrder1(vnbr, lrbx) && validateOrder1(vnbr, lrby) && validateOrder1(vnbr, lrbz));
            }
            if (col_smt)
            {
                assert(validateOrder1(vnbr, lvx) && validateOrder1(vnbr, lvy) && validateOrder1(vnbr, lvz));
                assert(validateOrder1(vnbr, lcfx) && validateOrder1(vnbr, lcfy) && validateOrder1(vnbr, lcfz));
                assert(validateOrder1(vnbr, lcbx) && validateOrder1(vnbr, lcby) && validateOrder1(vnbr, lcbz));
            }
            if (d0_smt)
            {
                assert(validateOrder1(vnbr, lvx) && validateOrder1(vnbr, lvy) && validateOrder1(vnbr, lvz));
                assert(validateOrder1(vnbr, ld0fx) && validateOrder1(vnbr, ld0fy) && validateOrder1(vnbr, ld0fz));
                assert(validateOrder1(vnbr, ld0bx) && validateOrder1(vnbr, ld0by) && validateOrder1(vnbr, ld0bz));
            }
            if (d1_smt)
            {
                assert(validateOrder1(vnbr, lvx) && validateOrder1(vnbr, lvy) && validateOrder1(vnbr, lvz));
                assert(validateOrder1(vnbr, ld1fx) && validateOrder1(vnbr, ld1fy) && validateOrder1(vnbr, ld1fz));
                assert(validateOrder1(vnbr, ld1bx) && validateOrder1(vnbr, ld1by) && validateOrder1(vnbr, ld1bz));
            }
        }

        if (row_smt)
        {
            assert(rf != rb);
            Eigen::Vector3d Ver(GlobVars[lvx], GlobVars[lvy], GlobVars[lvz]);
            Eigen::Vector3d Vfr(GlobVars[lrfx], GlobVars[lrfy], GlobVars[lrfz]);
            Eigen::Vector3d Vbk(GlobVars[lrbx], GlobVars[lrby], GlobVars[lrbz]);
            double scale0 = (Ver - Vfr).norm();
            double scale1 = (Ver - Vbk).norm();
            if ((Vbk - Vfr).norm() == 0)
            {
                std::cout << "Error: rf, " << rf << ", rb " << rb << "\n"<<"V\n"<<V<<"\n";
                std::cout<<"Vbk, "<<Vbk.transpose()<<"\nVfr, "<<Vfr.transpose()<<"\n";
                
            }
            assert((Vbk - Vfr).norm() > 0);
            
            counter = i;
            push_fairness_conditions(tripletes, energy, lvx, lrfx, lrbx, Ver[0], Vfr[0], Vbk[0], counter, scale0, scale1, order);
            counter = i + vnbr;
            push_fairness_conditions(tripletes, energy, lvy, lrfy, lrby, Ver[1], Vfr[1], Vbk[1], counter, scale0, scale1, order);
            counter = i + vnbr * 2;
            push_fairness_conditions(tripletes, energy, lvz, lrfz, lrbz, Ver[2], Vfr[2], Vbk[2], counter, scale0, scale1, order);
            
        }

        if (col_smt)
        {
            Eigen::Vector3d Ver(GlobVars[lvx], GlobVars[lvy], GlobVars[lvz]);
            Eigen::Vector3d Vfr(GlobVars[lcfx], GlobVars[lcfy], GlobVars[lcfz]);
            Eigen::Vector3d Vbk(GlobVars[lcbx], GlobVars[lcby], GlobVars[lcbz]);
            double scale0 = (Ver - Vfr).norm();
            double scale1 = (Ver - Vbk).norm();
            assert((Vbk - Vfr).norm() > 0);
            counter = i + vnbr * 3;
            push_fairness_conditions(tripletes, energy, lvx, lcfx, lcbx, Ver[0], Vfr[0], Vbk[0], counter, scale0, scale1, order);
            counter = i + vnbr * 4;
            push_fairness_conditions(tripletes, energy, lvy, lcfy, lcby, Ver[1], Vfr[1], Vbk[1], counter, scale0, scale1, order);
            counter = i + vnbr * 5;
            push_fairness_conditions(tripletes, energy, lvz, lcfz, lcbz, Ver[2], Vfr[2], Vbk[2], counter, scale0, scale1, order);
            
        }

        if (d0_smt)
        {
            Eigen::Vector3d Ver(GlobVars[lvx], GlobVars[lvy], GlobVars[lvz]);
            Eigen::Vector3d Vfr(GlobVars[ld0fx], GlobVars[ld0fy], GlobVars[ld0fz]);
            Eigen::Vector3d Vbk(GlobVars[ld0bx], GlobVars[ld0by], GlobVars[ld0bz]);
            double scale0 = (Ver - Vfr).norm();
            double scale1 = (Ver - Vbk).norm();
            assert((Vbk - Vfr).norm() > 0);
            counter = i + vnbr * 6;
            push_fairness_conditions(tripletes, energy, lvx, ld0fx, ld0bx, Ver[0], Vfr[0], Vbk[0], counter, scale0, scale1, order);
            counter = i + vnbr * 7;
            push_fairness_conditions(tripletes, energy, lvy, ld0fy, ld0by, Ver[1], Vfr[1], Vbk[1], counter, scale0, scale1, order);
            counter = i + vnbr * 8;
            push_fairness_conditions(tripletes, energy, lvz, ld0fz, ld0bz, Ver[2], Vfr[2], Vbk[2], counter, scale0, scale1, order);
        }
        if (d1_smt)
        {
            Eigen::Vector3d Ver(GlobVars[lvx], GlobVars[lvy], GlobVars[lvz]);
            Eigen::Vector3d Vfr(GlobVars[ld1fx], GlobVars[ld1fy], GlobVars[ld1fz]);
            Eigen::Vector3d Vbk(GlobVars[ld1bx], GlobVars[ld1by], GlobVars[ld1bz]);
            double scale0 = (Ver - Vfr).norm();
            double scale1 = (Ver - Vbk).norm();
            assert((Vbk - Vfr).norm() > 0);
            counter = i + vnbr * 9;
            push_fairness_conditions(tripletes, energy, lvx, ld1fx, ld1bx, Ver[0], Vfr[0], Vbk[0], counter, scale0, scale1, order);
            counter = i + vnbr * 10;
            push_fairness_conditions(tripletes, energy, lvy, ld1fy, ld1by, Ver[1], Vfr[1], Vbk[1], counter, scale0, scale1, order);
            counter = i + vnbr * 11;
            push_fairness_conditions(tripletes, energy, lvz, ld1fz, ld1bz, Ver[2], Vfr[2], Vbk[2], counter, scale0, scale1, order);
        }
    }

    // std::cout<<"check 2"<<std::endl;
    spMat J;
    J.resize(energy.size(), GlobVars.size());
    J.setFromTriplets(tripletes.begin(), tripletes.end());
    H = J.transpose() * J;
    B = -J.transpose() * energy;
}

class variableLocator
{
public:
    variableLocator(){};
    void getLocations(int order, int vnbr, int vid, int rf, int rb, int cf, int cb,
                      int d0f, int d0b, int d1f, int d1b, int bnm_start_location, int family);
    // this is to choose which polylines to use for normal vectors.
    void getLocationsNormalCondition(int order, int vnbr, int vid, int &rf, int &rb, int &cf, int &cb,
                                                  int d0f, int d0b, int d1f, int d1b, int bnm_start_location, int whichDiagonal);
    int lvx;
    int lvy;
    int lvz;
    int lfx;
    int lfy;
    int lfz;
    int lbx;
    int lby;
    int lbz;
    int lnx;
    int lny;
    int lnz;
    int lrx;
    int lry;
    int lrz;
    int lrx1;
    int lry1;
    int lrz1;
	int lrx2;
	int lry2;
	int lrz2;

    int lrfx;
    int lrfy;
    int lrfz;

    int lrbx;
    int lrby;
    int lrbz;

    int lcfx;
    int lcfy;
    int lcfz;

    int lcbx;
    int lcby;
    int lcbz;
};



void variableLocator::getLocations(int order, int vnbr, int vid, int rf, int rb, int cf, int cb,
                                   int d0f, int d0b, int d1f, int d1b, int bnm_start_location, int family)
{
    int nbrOffset = 9;
    if (order == 2)
    {
        nbrOffset = 12;
    }
	if (order == 3)
	{
		nbrOffset = 15;
	}
    if (order == 0)
    {
        lvx = vid;
        lvy = vid + vnbr;
        lvz = vid + 2 * vnbr;
        lnx = vnbr * 3 + vid;
        lny = vnbr * 4 + vid;
        lnz = vnbr * 5 + vid;
        lrx = bnm_start_location + vid;
        lry = bnm_start_location + vid + vnbr;
        lrz = bnm_start_location + vid + vnbr * 2;
    }
    else
    {
        lvx = nbrOffset * vid;
        lvy = nbrOffset * vid + 1;
        lvz = nbrOffset * vid + 2;
        lnx = nbrOffset * vid + 3;
        lny = nbrOffset * vid + 4;
        lnz = nbrOffset * vid + 5;
        lrx = nbrOffset * vid + 6;
        lry = nbrOffset * vid + 7;
        lrz = nbrOffset * vid + 8;
    }
    if(order == 2 || order == 3) // AGG, GGG needs another binormal vector
    {
        lrx1 = nbrOffset * vid + 9;
        lry1 = nbrOffset * vid + 10;
        lrz1 = nbrOffset * vid + 11;
		if (order == 3)
		{
			lrx2 = nbrOffset * vid + 12;
			lry2 = nbrOffset * vid + 13;
			lrz2 = nbrOffset * vid + 14;
		}
    }

    if (family == 0) // rows
    {
        locatorRules(rf, vnbr, order, lfx, lfy, lfz);
        locatorRules(rb, vnbr, order, lbx, lby, lbz);
    }
    if (family == 1) // cols
    {
        locatorRules(cf, vnbr, order, lfx, lfy, lfz);
        locatorRules(cb, vnbr, order, lbx, lby, lbz);
    }
    if (family == 2) // d0
    {
        locatorRules(d0f, vnbr, order, lfx, lfy, lfz);
        locatorRules(d0b, vnbr, order, lbx, lby, lbz);
    }

    if (family == 3) // d1
    {
        locatorRules(d1f, vnbr, order, lfx, lfy, lfz);
        locatorRules(d1b, vnbr, order, lbx, lby, lbz);
    }
}

// bnm_start is the id where the binormal starts in the variable matrix
// from bnm_start to bnm_start + vnbr * 3 are the binormal variables of this family of curves
// family: 0: row. 1: col. 2: d0. 3: d1.
// whichBnm: this variable is for the propagating AGG, GGG method (when order == 2 or 3) to point out which binormal is used
void QuadOpt::assemble_binormal_conditions(spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy, int family,
                                           int bnm_start, const int order, const int whichBnm)
{
    // std::cout<<"check 1, family "<<family<<std::endl;
    std::vector<Trip> tripletes;
    int vnbr = V.rows();
    energy = Eigen::VectorXd::Zero(vnbr * 3); 
    tripletes.reserve(vnbr * 25); 
    int counter = 0;
    for (int i = 0; i < vnbr; i++)
    {
        // the ids
        int vid = i;
        int rf = row_front[vid];
        int rb = row_back[vid];
        int cf = col_front[vid];
        int cb = col_back[vid];
        int d0f = d0_front[vid];
        int d0b = d0_back[vid];
        int d1f = d1_front[vid];
        int d1b = d1_back[vid];
        // the locations

        variableLocator locator;
        locator.getLocations(order, vnbr, vid, rf, rb,
                             cf, cb, d0f, d0b, d1f, d1b, bnm_start, family);
            // the vertex
        int lvx = locator.lvx;
        int lvy = locator.lvy;
        int lvz = locator.lvz;
        // the front and back vertices
        int lfx = locator.lfx;
        int lfy = locator.lfy;
        int lfz = locator.lfz;
        int lbx = locator.lbx;
        int lby = locator.lby;
        int lbz = locator.lbz;
        // int lnx = locator.lnx;
        // int lny = locator.lny;
        // int lnz = locator.lnz;
        int lrx = locator.lrx;
        int lry = locator.lry;
        int lrz = locator.lrz;
        if ((order == 2 || order == 3) && whichBnm == 1) // AGG or GGG and use the second binormal attached on this vertex
        {
            // the binormal is the second one.
            lrx = locator.lrx1;
            lry = locator.lry1;
            lrz = locator.lrz1;
        }
		if (order == 3 && whichBnm == 2)// GGG, and use the third binormal on this vertex
		{
			lrx = locator.lrx2;
			lry = locator.lry2;
			lrz = locator.lrz2;
		}
        bool compute = true;
        
        if (lfx < 0 || lbx < 0)
        {
            compute = false;
        }
        
        if(compute == false){ // the vertex is on the boundary
            if (!ComputeAuxiliaries)
            {
                // here we choose doing nothing and care only about inner vertices
            }
            continue;
        }
        Eigen::Vector3d Ver(GlobVars[lvx], GlobVars[lvy], GlobVars[lvz]);
        Eigen::Vector3d Vf(GlobVars[lfx], GlobVars[lfy], GlobVars[lfz]);
        Eigen::Vector3d Vb(GlobVars[lbx], GlobVars[lby], GlobVars[lbz]);
        if(ComputeAuxiliaries){
            //init the binormals
            Eigen::Vector3d real_r = ((Ver - Vf).cross(Vb - Ver)).normalized();
            if (real_r.norm() == 0)
            {
                real_r = Eigen::Vector3d(0, 0, 1);
            }
            GlobVars[lrx] = real_r[0];
            GlobVars[lry] = real_r[1];
            GlobVars[lrz] = real_r[2];
        }
        
        Eigen::Vector3d r(GlobVars[lrx], GlobVars[lry], GlobVars[lrz]);
        double dis0 = (Vf - Ver).norm();
        double dis1 = (Vb - Ver).norm();
        double dis = dis0 * dis1;
        Eigen::Vector3d fc = Vf - Ver;
        Eigen::Vector3d bc = Vb - Ver;
        // r * (Vf - Ver) = 0
        tripletes.push_back(Trip(i, lrx, (Vf(0) - Ver(0)) / dis0));
        tripletes.push_back(Trip(i, lry, (Vf(1) - Ver(1)) / dis0));
        tripletes.push_back(Trip(i, lrz, (Vf(2) - Ver(2)) / dis0));

        tripletes.push_back(Trip(i, lfx, (r(0)) / dis0));
        tripletes.push_back(Trip(i, lfy, (r(1)) / dis0));
        tripletes.push_back(Trip(i, lfz, (r(2)) / dis0));

        tripletes.push_back(Trip(i, lvx, -(r(0)) / dis0));
        tripletes.push_back(Trip(i, lvy, -(r(1)) / dis0));
        tripletes.push_back(Trip(i, lvz, -(r(2)) / dis0));

        energy[i] = r.dot(Vf - Ver) / dis0;
        // r * (Vb - Ver) = 0
        tripletes.push_back(Trip(i + vnbr, lrx, (Vb(0) - Ver(0)) / dis1));
        tripletes.push_back(Trip(i + vnbr, lry, (Vb(1) - Ver(1)) / dis1));
        tripletes.push_back(Trip(i + vnbr, lrz, (Vb(2) - Ver(2)) / dis1));

        tripletes.push_back(Trip(i + vnbr, lbx, (r(0)) / dis1));
        tripletes.push_back(Trip(i + vnbr, lby, (r(1)) / dis1));
        tripletes.push_back(Trip(i + vnbr, lbz, (r(2)) / dis1));

        tripletes.push_back(Trip(i + vnbr, lvx, -(r(0)) / dis1));
        tripletes.push_back(Trip(i + vnbr, lvy, -(r(1)) / dis1));
        tripletes.push_back(Trip(i + vnbr, lvz, -(r(2)) / dis1));

        energy[i + vnbr] = r.dot(Vb - Ver) / dis1;
        // r * r = 1
        tripletes.push_back(Trip(i + vnbr * 2, lrx, 2 * r(0)));
        tripletes.push_back(Trip(i + vnbr * 2, lry, 2 * r(1)));
        tripletes.push_back(Trip(i + vnbr * 2, lrz, 2 * r(2)));
        energy[i + vnbr * 2] = r.dot(r) - 1;
    }
    spMat J;
    J.resize(energy.size(), GlobVars.size());
    J.setFromTriplets(tripletes.begin(), tripletes.end());
    H = J.transpose() * J;
    B = -J.transpose() * energy;
}
// lv is ver location, lf is front location, lb is back location, cid is the id of the condition
void push_asymptotic_condition(std::vector<Trip> &tripletes, Eigen::VectorXd &energy, const int cid,
                        const int lvx, const int lvy, const int lvz,
                        const int lfx, const int lfy, const int lfz,
                        const int lbx, const int lby, const int lbz,
                        const int lnx, const int lny, const int lnz,
                        const Eigen::Vector3d &Vf, const Eigen::Vector3d &Ver, const Eigen::Vector3d &Vb, 
                        const Eigen::Vector3d &norm)
{

    double dis0 = (Vf - Ver).norm();
    double dis1 = (Vb - Ver).norm();
    // norm * (Vf-Ver) = 0
    tripletes.push_back(Trip(cid, lnx, (Vf(0) - Ver(0)) / dis0));
    tripletes.push_back(Trip(cid, lny, (Vf(1) - Ver(1)) / dis0));
    tripletes.push_back(Trip(cid, lnz, (Vf(2) - Ver(2)) / dis0));

    tripletes.push_back(Trip(cid, lfx, (norm(0)) / dis0));
    tripletes.push_back(Trip(cid, lfy, (norm(1)) / dis0));
    tripletes.push_back(Trip(cid, lfz, (norm(2)) / dis0));

    tripletes.push_back(Trip(cid, lvx, -(norm(0)) / dis0));
    tripletes.push_back(Trip(cid, lvy, -(norm(1)) / dis0));
    tripletes.push_back(Trip(cid, lvz, -(norm(2)) / dis0));

    energy[cid] = norm.dot(Vf - Ver) / dis0;

    // norm * (Vb-Ver) = 0
    tripletes.push_back(Trip(cid + 1, lnx, (Vb(0) - Ver(0)) / dis1));
    tripletes.push_back(Trip(cid + 1, lny, (Vb(1) - Ver(1)) / dis1));
    tripletes.push_back(Trip(cid + 1, lnz, (Vb(2) - Ver(2)) / dis1));

    tripletes.push_back(Trip(cid + 1, lbx, (norm(0)) / dis1));
    tripletes.push_back(Trip(cid + 1, lby, (norm(1)) / dis1));
    tripletes.push_back(Trip(cid + 1, lbz, (norm(2)) / dis1));

    tripletes.push_back(Trip(cid + 1, lvx, -(norm(0)) / dis1));
    tripletes.push_back(Trip(cid + 1, lvy, -(norm(1)) / dis1));
    tripletes.push_back(Trip(cid + 1, lvz, -(norm(2)) / dis1));

    energy[cid + 1] = norm.dot(Vb - Ver) / dis1;
    
}

void push_geodesic_condition(std::vector<Trip> &tripletes, Eigen::VectorXd &energy, const int cid,
                             const int lnx, const int lny, const int lnz,
                             const int lrx, const int lry, const int lrz,
                             const Eigen::Vector3d &norm, const Eigen::Vector3d &r)
{
    // norm * r = 0
    tripletes.push_back(Trip(cid, lnx, r(0)));
    tripletes.push_back(Trip(cid, lny, r(1)));
    tripletes.push_back(Trip(cid, lnz, r(2)));

    tripletes.push_back(Trip(cid, lrx, norm(0)));
    tripletes.push_back(Trip(cid, lry, norm(1)));
    tripletes.push_back(Trip(cid, lrz, norm(2)));

    energy[cid] = norm.dot(r);
}

// type: 0 disabled. 1 means asymptotic, 2 means geodesic, 3 pseudo-geodesic
// aux_start_location is the start location of the auxiliaries of this family
// when order = 0, then the variables are: all the vertices, all the normals, all others...
// when order = 1 (AAG), the variables are: ver+normal+binormal, ... ,
// when order = 2, or 3, it means AGG or GGG
void
QuadOpt::assemble_pg_extreme_cases(spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy,
                                   const int type, const int family, const int bnm_start_location, const int order,
                                   const int whichBnm)
{

    std::vector<Trip> tripletes;
    int vnbr = V.rows();
    int energy_size = 0;
    if(type == 1){ // asymptotic
        energy_size = vnbr * 2;
    }
    if(type == 2){ // geodesic
        energy_size = vnbr;
    }
    if(type == 3){
        std::cout<<"ERROR: Please do not use this function for general pseudo-geodesic solving"<<std::endl;
        return;
    }
    energy = Eigen::VectorXd::Zero(energy_size);
    tripletes.reserve(vnbr * 18);
    int counter = 0;
    for (int i = 0; i < vnbr; i++)
    {
        // the ids
        int vid = i;
        int rf = row_front[vid];
        int rb = row_back[vid];
        int cf = col_front[vid];
        int cb = col_back[vid];
        int d0f = d0_front[vid];
        int d0b = d0_back[vid];
        int d1f = d1_front[vid];
        int d1b = d1_back[vid];

         // the locations
        variableLocator locator;
        locator.getLocations(order, vnbr, vid, rf, rb,
                             cf, cb, d0f, d0b, d1f, d1b, bnm_start_location, family);
            // the vertex
        int lvx = locator.lvx;
        int lvy = locator.lvy;
        int lvz = locator.lvz;
        // the front and back vertices
        int lfx = locator.lfx;
        int lfy = locator.lfy;
        int lfz = locator.lfz;
        int lbx = locator.lbx;
        int lby = locator.lby;
        int lbz = locator.lbz;
        int lnx = locator.lnx;
        int lny = locator.lny;
        int lnz = locator.lnz;
        bool compute = true;
        
        if (rf < 0 || rb < 0 || cf < 0 || cb < 0)// on the boundary it does not have a normal vector properly defined
        {
            compute = false;
        }
        if (lfx < 0 || lbx < 0)
        {
            compute = false;
        }
        if (!compute)
        {
            continue;
        }
        if(!(lvx >= 0 && lfx >= 0 && lbx >= 0 && lnx >= 0)){
            std::cout << lvx << ", " << lfx << ", " << lbx << ", " << lnx << std::endl;
        }
        assert(lvx >= 0 && lfx >= 0 && lbx >= 0 && lnx >= 0);
        assert(lvy >= 0 && lfy >= 0 && lby >= 0 && lny >= 0);
        assert(lvz >= 0 && lfz >= 0 && lbz >= 0 && lnz >= 0);

        assert(lvx <= GlobVars.size() && lfx <= GlobVars.size() && lbx <= GlobVars.size() && lnx <= GlobVars.size());
        assert(lvy <= GlobVars.size() && lfy <= GlobVars.size() && lby <= GlobVars.size() && lny <= GlobVars.size());
        assert(lvz <= GlobVars.size() && lfz <= GlobVars.size() && lbz <= GlobVars.size() && lnz <= GlobVars.size());
        Eigen::Vector3d Ver(GlobVars[lvx], GlobVars[lvy], GlobVars[lvz]);
        Eigen::Vector3d Vfr(GlobVars[lfx], GlobVars[lfy], GlobVars[lfz]);
        Eigen::Vector3d Vbk(GlobVars[lbx], GlobVars[lby], GlobVars[lbz]);
        Eigen::Vector3d norm(GlobVars[lnx], GlobVars[lny], GlobVars[lnz]);
        if (type == 1)
        { // asymptotic
            counter = i * 2;
            push_asymptotic_condition(tripletes, energy, counter, lvx, lvy, lvz, lfx, lfy, lfz, lbx, lby, lbz, lnx, lny, lnz,
                                      Vfr, Ver, Vbk, norm);
        }
        if(type == 2){// geodesic
            counter = i;
            int lrx = locator.lrx;
            int lry = locator.lry;
            int lrz = locator.lrz;
			if ((order == 2 || order == 3) && whichBnm == 1) // AGG or GGG and use the second binormal attached on this vertex
			{
				// the binormal is the second one.
				lrx = locator.lrx1;
				lry = locator.lry1;
				lrz = locator.lrz1;
			}
			if (order == 3 && whichBnm == 2)// GGG, and use the third binormal on this vertex
			{
				lrx = locator.lrx2;
				lry = locator.lry2;
				lrz = locator.lrz2;
			}
            Eigen::Vector3d r(GlobVars[lrx], GlobVars[lry], GlobVars[lrz]);
            push_geodesic_condition(tripletes, energy, counter, lnx, lny, lnz, lrx, lry, lrz, norm, r);
        }
    }
    // std::cout<<"check 2"<<std::endl;
    spMat J;
    J.resize(energy.size(), GlobVars.size());
    J.setFromTriplets(tripletes.begin(), tripletes.end());
    H = J.transpose() * J;
    B = -J.transpose() * energy;
}
void QuadOpt::assemble_pg_cases(const double angle_radian, spMat &H, Eigen::VectorXd &B, Eigen::VectorXd &energy,
                                const int family, const int bnm_start_location)
{
    if(angle_radian == 90 * LSC_PI / 180.)
    {
        std::cout<<"G...\n";
        assemble_pg_extreme_cases(H,B, energy,2,family,bnm_start_location);
        return;
    }
    std::vector<Trip> tripletes;
    int vnbr = V.rows();
    int energy_size = 0;
    double angle_avg = 0;
    int counter = 0;

    energy_size = vnbr * 5;

    energy = Eigen::VectorXd::Zero(energy_size);
    tripletes.reserve(vnbr * 40);
    for (int i = 0; i < vnbr; i++)
    {
        // the ids
        int vid = i;
        int rf = row_front[vid];
        int rb = row_back[vid];
        int cf = col_front[vid];
        int cb = col_back[vid];
        int d0f = d0_front[vid];
        int d0b = d0_back[vid];
        int d1f = d1_front[vid];
        int d1b = d1_back[vid];

         // the locations
        // the vertex
        int lvx = vid;
        int lvy = vid + vnbr;
        int lvz = vid + 2 * vnbr;
        // the front and back vertices
        int lfx, lfy, lfz, lbx, lby, lbz;
        // the normal vectors
        int lnx = vnbr * 3 + i;
        int lny = vnbr * 4 + i;
        int lnz = vnbr * 5 + i;
        // the binormals
        int lrx = bnm_start_location + i;
        int lry = bnm_start_location + i + vnbr;
        int lrz = bnm_start_location + i + vnbr * 2;
        // the side vectors
        int lux = bnm_start_location + i + vnbr * 3;
        int luy = bnm_start_location + i + vnbr * 4;
        int luz = bnm_start_location + i + vnbr * 5;
        
        bool compute = true;
        if (family == 0) // rows
        {
            lfx = rf;
            lfy = rf + vnbr;
            lfz = rf + vnbr * 2;

            lbx = rb;
            lby = rb + vnbr;
            lbz = rb + vnbr * 2;
        }
        if (family == 1) // cols
        {
            lfx = cf;
            lfy = cf + vnbr;
            lfz = cf + vnbr * 2;

            lbx = cb;
            lby = cb + vnbr;
            lbz = cb + vnbr * 2;
        }
        if (family == 2) // d0
        {
            lfx = d0f;
            lfy = d0f + vnbr;
            lfz = d0f + vnbr * 2;

            lbx = d0b;
            lby = d0b + vnbr;
            lbz = d0b + vnbr * 2;
        }

        if (family == 3) // d1
        {

            lfx = d1f;
            lfy = d1f + vnbr;
            lfz = d1f + vnbr * 2;

            lbx = d1b;
            lby = d1b + vnbr;
            lbz = d1b + vnbr * 2;
        }
        if (rf < 0 || rb < 0 || cf < 0 || cb < 0)
        {
            compute = false;
        }
        if (lfx < 0 || lbx < 0)
        {
            compute = false;
        }
        if (!compute)
        {
            continue;
        }
        Eigen::Vector3d Ver(GlobVars[lvx], GlobVars[lvy], GlobVars[lvz]);
        Eigen::Vector3d Vfr(GlobVars[lfx], GlobVars[lfy], GlobVars[lfz]);
        Eigen::Vector3d Vbk(GlobVars[lbx], GlobVars[lby], GlobVars[lbz]);
        Eigen::Vector3d norm(GlobVars[lnx], GlobVars[lny], GlobVars[lnz]);
        Eigen::Vector3d r(GlobVars[lrx], GlobVars[lry], GlobVars[lrz]);

        Eigen::Vector3d tangent = Vbk - Vfr; // This sign is changable depends on what user defines of the curve tangent directions
        Eigen::Vector3d real_u = (norm.cross(tangent)).normalized();
        if(ComputeAuxiliaries){
            // std::cout<<"Init Auxiliary Variables ..."<<std::endl;
            GlobVars[lux] = real_u[0];
            GlobVars[luy] = real_u[1];
            GlobVars[luz] = real_u[2];
        }
        Eigen::Vector3d u(GlobVars[lux], GlobVars[luy], GlobVars[luz]);
        if(real_u.dot(u)<0){
            std::cout<<"flip u"<<std::endl;
            u * -1;
            GlobVars[lux] *= -1;
            GlobVars[luy] *= -1;
            GlobVars[luz] *= -1;
        }
        

        double dis0 = (Ver - Vfr).norm();
        double dis1 = (Ver - Vbk).norm();
        double cos_angle = cos(angle_radian);
        double sin_angle = sin(angle_radian);

        // (binormal * n)^2 = cos^2
        double ndb = norm.dot(r);
        tripletes.push_back(Trip(i, lrx, 2 * ndb * norm(0)));
        tripletes.push_back(Trip(i, lry, 2 * ndb * norm(1)));
        tripletes.push_back(Trip(i, lrz, 2 * ndb * norm(2)));

        tripletes.push_back(Trip(i, lnx, 2 * ndb * r(0)));
        tripletes.push_back(Trip(i, lny, 2 * ndb * r(1)));
        tripletes.push_back(Trip(i, lnz, 2 * ndb * r(2)));

        energy[i] = ndb * ndb - cos_angle * cos_angle;

        // (binormal * n) * (binormal * u) = sin * cos
        // regard n and u as constant values in each iteration, since u is related to the tangent vector, which will not change much.
        double udb = u.dot(r);

        tripletes.push_back(Trip(i + vnbr, lnx, udb * r[0]));
        tripletes.push_back(Trip(i + vnbr, lny, udb * r[1]));
        tripletes.push_back(Trip(i + vnbr, lnz, udb * r[2]));

        tripletes.push_back(Trip(i + vnbr, lrx, ndb * u[0] + udb * norm[0]));
        tripletes.push_back(Trip(i + vnbr, lry, ndb * u[1] + udb * norm[1]));
        tripletes.push_back(Trip(i + vnbr, lrz, ndb * u[2] + udb * norm[2]));

        tripletes.push_back(Trip(i + vnbr, lux, ndb * r[0]));
        tripletes.push_back(Trip(i + vnbr, luy, ndb * r[1]));
        tripletes.push_back(Trip(i + vnbr, luz, ndb * r[2]));

        energy[i + vnbr] = ndb * udb - sin_angle * cos_angle;

        // u * u = 1
        tripletes.push_back(Trip(i + vnbr * 2, lux, 2 * u[0]));
        tripletes.push_back(Trip(i + vnbr * 2, luy, 2 * u[1]));
        tripletes.push_back(Trip(i + vnbr * 2, luz, 2 * u[2]));

        energy[i + vnbr * 2] = u.dot(u) - 1;

        // u * n = 0
        tripletes.push_back(Trip(i + vnbr * 3, lux, norm[0]));
        tripletes.push_back(Trip(i + vnbr * 3, luy, norm[1]));
        tripletes.push_back(Trip(i + vnbr * 3, luz, norm[2]));

        tripletes.push_back(Trip(i + vnbr * 3, lnx, u[0]));
        tripletes.push_back(Trip(i + vnbr * 3, lny, u[1]));
        tripletes.push_back(Trip(i + vnbr * 3, lnz, u[2]));

        energy[i + vnbr * 3] = u.dot(norm);

        // u * (Vfr - Vbk) = 0
        double disFB = (Vfr - Vbk).norm();
        tripletes.push_back(Trip(i + vnbr * 4, lux, (Vfr - Vbk)[0] / disFB));
        tripletes.push_back(Trip(i + vnbr * 4, luy, (Vfr - Vbk)[1] / disFB));
        tripletes.push_back(Trip(i + vnbr * 4, luz, (Vfr - Vbk)[2] / disFB));

        tripletes.push_back(Trip(i + vnbr * 4, lfx, u[0] / disFB));
        tripletes.push_back(Trip(i + vnbr * 4, lfy, u[1] / disFB));
        tripletes.push_back(Trip(i + vnbr * 4, lfz, u[2] / disFB));

        tripletes.push_back(Trip(i + vnbr * 4, lbx, -u[0] / disFB));
        tripletes.push_back(Trip(i + vnbr * 4, lby, -u[1] / disFB));
        tripletes.push_back(Trip(i + vnbr * 4, lbz, -u[2] / disFB));

        energy[i + vnbr * 4] = u.dot(Vfr - Vbk) / disFB;

        // get the average pseudo-geodesic angle
        if (Estimate_PG_Angles)
        {
            double cosang = r.dot(norm);
            double sinang = r.dot(u);
            if (sinang < 0)
            {
                cosang *= -1;
            }
            double angle = acos(cosang);// get angle in radian
            angle_avg += 180 / LSC_PI * angle; // get angle in degree.
            counter++;
        }
    }
    if(Estimate_PG_Angles){
        angle_avg /= counter;
        std::cout<<"Estimating PG Angles: "<<angle_avg<<std::endl;
    }
    // std::cout<<"detail, "<<energy.segment(0,vnbr).norm()<<", "<<energy.segment(vnbr, vnbr).norm()<<", tgt, "<<angle_radian<<", ";
    // std::cout<<"check 2"<<std::endl;
    spMat J;
    J.resize(energy.size(), GlobVars.size());
    J.setFromTriplets(tripletes.begin(), tripletes.end());
    H = J.transpose() * J;
    B = -J.transpose() * energy;
}

void variableLocator::getLocationsNormalCondition(int order, int vnbr, int vid, int &rf, int &rb, int &cf, int &cb,
                                                  int d0f, int d0b, int d1f, int d1b, int bnm_start_location, int whichDiagonal)
{
    if (cf < 0)
    {
        cf = vid;
    }
    if (cb < 0)
    {
        cb = vid;
    }
    int varOffset = 9;
    if (order == 2)
    {
        varOffset = 12;
    }
	if (order == 3)
	{
		varOffset = 15;
	}
    if (order == 0)
    {
        lvx = vid;
        lvy = vid + vnbr;
        lvz = vid + 2 * vnbr;
        lnx = vnbr * 3 + vid;
        lny = vnbr * 4 + vid;
        lnz = vnbr * 5 + vid;
        // lrx = bnm_start_location + vid;
        // lry = bnm_start_location + vid + vnbr;
        // lrz = bnm_start_location + vid + vnbr * 2;
    }

    else
    {
        lvx = varOffset * vid;
        lvy = varOffset * vid + 1;
        lvz = varOffset * vid + 2;
        lnx = varOffset * vid + 3;
        lny = varOffset * vid + 4;
        lnz = varOffset * vid + 5;
        // lrx = varOffset * vid + 6;
        // lry = varOffset * vid + 7;
        // lrz = varOffset * vid + 8;
    }
    if (order == 0)
    {
        if (rf < 0)
        {
            rf = vid;
        }
        if (rb < 0)
        {
            rb = vid;
        }

        lrfx = rf;
        lrfy = rf + vnbr;
        lrfz = rf + vnbr * 2;

        lrbx = rb;
        lrby = rb + vnbr;
        lrbz = rb + vnbr * 2;

        lcfx = cf;
        lcfy = cf + vnbr;
        lcfz = cf + vnbr * 2;

        lcbx = cb;
        lcby = cb + vnbr;
        lcbz = cb + vnbr * 2;
    }
    else
    {
        lcfx = varOffset * cf;
        lcfy = varOffset * cf + 1;
        lcfz = varOffset * cf + 2;

        lcbx = varOffset * cb;
        lcby = varOffset * cb + 1;
        lcbz = varOffset * cb + 2;
        // int df = d0f;
        // int db = d0b;
        // if (whichDiagonal == 1) // choose d1.
        // {
        //     df = d1f;
        //     db = d1b;
        // }

        // if (df < 0)
        // {
        //     df = vid;
        // }
        // if (db < 0)
        // {
        //     db = vid;
        // }
        // lrfx = varOffset * df;
        // lrfy = varOffset * df + 1;
        // lrfz = varOffset * df + 2;
        // lrbx = varOffset * db;
        // lrby = varOffset * db + 1;
        // lrbz = varOffset * db + 2;
        if (rf < 0)
        {
            rf = vid;
        }
        if (rb < 0)
        {
            rb = vid;
        }

        lrfx = varOffset * rf;
        lrfy = varOffset * rf + 1;
        lrfz = varOffset * rf + 2;

        lrbx = varOffset * rb;
        lrby = varOffset * rb + 1;
        lrbz = varOffset * rb + 2;
    }
}

void QuadOpt::assemble_normal_conditions(spMat& H, Eigen::VectorXd& B, Eigen::VectorXd &energy, const int order){
    std::vector<Trip> tripletes;
    int vnbr = V.rows();
    energy = Eigen::VectorXd::Zero(vnbr * 3); // todo
    tripletes.reserve(vnbr * 21); // todo
    int counter = 0;
    std::vector<Eigen::Vector3d> edge_nm0, edge_nm1; 
    edge_nm0.reserve(vnbr);
    edge_nm1.reserve(vnbr);
    for (int i = 0; i < vnbr; i++)
    {
        int vid = i;
        int rf = row_front[vid];
        int rb = row_back[vid];
        int cf = col_front[vid];
        int cb = col_back[vid];
        int d0f = d0_front[vid];
        int d0b = d0_back[vid];
        int d1f = d1_front[vid];
        int d1b = d1_back[vid];
 
        // if on one direction, there is no vertex, we use the current vertex to replace it for normal vector 
        // calculation. In such a way, the orientaion of the normal vector will not change.
        variableLocator locator;
        locator.getLocationsNormalCondition(order, vnbr, vid, rf, rb,
                             cf, cb, d0f, d0b, d1f, d1b, -1, WhichDiagonal);
            // the vertex
        int lvx = locator.lvx;
        int lvy = locator.lvy;
        int lvz = locator.lvz;
        // the front and back vertices

        int lnx = locator.lnx;
        int lny = locator.lny;
        int lnz = locator.lnz;

        int lrfx = locator.lrfx;
        int lrfy = locator.lrfy;
        int lrfz = locator.lrfz;

        int lrbx = locator.lrbx;
        int lrby = locator.lrby;
        int lrbz = locator.lrbz;

        int lcfx = locator.lcfx;
        int lcfy = locator.lcfy;
        int lcfz = locator.lcfz;

        int lcbx = locator.lcbx;
        int lcby = locator.lcby;
        int lcbz = locator.lcbz;
        // if (order != 0 && lrbx == lrfx) // when use diagonal to compute normal vectors, this means the diagonal does not exist
        // {
        //     continue;
        // }

        Eigen::Vector3d Ver(GlobVars[lvx], GlobVars[lvy], GlobVars[lvz]);
        Eigen::Vector3d Vrf(GlobVars[lrfx], GlobVars[lrfy], GlobVars[lrfz]);
        Eigen::Vector3d Vrb(GlobVars[lrbx], GlobVars[lrby], GlobVars[lrbz]);
        Eigen::Vector3d Vcf(GlobVars[lcfx], GlobVars[lcfy], GlobVars[lcfz]);
        Eigen::Vector3d Vcb(GlobVars[lcbx], GlobVars[lcby], GlobVars[lcbz]);
        Eigen::Vector3d real_n = ((Vrf - Vrb).cross(Vcf - Vcb)).normalized();
        if(ComputeAuxiliaries){
            //init the binormals
            
            GlobVars[lnx] = real_n[0];
            GlobVars[lny] = real_n[1];
            GlobVars[lnz] = real_n[2];
        }

        
        
        Eigen::Vector3d norm(GlobVars[lnx], GlobVars[lny], GlobVars[lnz]);
        if (real_n.dot(norm) < 0)
        {
            std::cout<<"flip n"<<std::endl;
            norm *= -1;
            GlobVars[lnx] *= -1;
            GlobVars[lny] *= -1;
            GlobVars[lnz] *= -1;
        }
        edge_nm0.push_back(Ver);
        edge_nm1.push_back(Ver + norm);
        double dis0 = (Vrf - Vrb).norm();
        double dis1 = (Vcf - Vcb).norm();

        // n * (Vrf - Vrb) = 0
        tripletes.push_back(Trip(i, lnx, (Vrf(0) - Vrb(0)) / dis0));
        tripletes.push_back(Trip(i, lny, (Vrf(1) - Vrb(1)) / dis0));
        tripletes.push_back(Trip(i, lnz, (Vrf(2) - Vrb(2)) / dis0));

        tripletes.push_back(Trip(i, lrfx, (norm(0)) / dis0));
        tripletes.push_back(Trip(i, lrfy, (norm(1)) / dis0));
        tripletes.push_back(Trip(i, lrfz, (norm(2)) / dis0));

        tripletes.push_back(Trip(i, lrbx, -(norm(0)) / dis0));
        tripletes.push_back(Trip(i, lrby, -(norm(1)) / dis0));
        tripletes.push_back(Trip(i, lrbz, -(norm(2)) / dis0));

        energy[i] = norm.dot(Vrf - Vrb) / dis0;
        if (isnan(energy[i]))
        {
            std::cout << "NAN 1\nnorm, " << norm.transpose() << ", values, " << Vrf.transpose()
             << ", " << Vrb.transpose() << ", dis0, " << dis0 << "\n";
        }

        // norm * (Vcf - Vcb) = 0
        tripletes.push_back(Trip(i + vnbr, lnx, (Vcf(0) - Vcb(0)) / dis1));
        tripletes.push_back(Trip(i + vnbr, lny, (Vcf(1) - Vcb(1)) / dis1));
        tripletes.push_back(Trip(i + vnbr, lnz, (Vcf(2) - Vcb(2)) / dis1));

        tripletes.push_back(Trip(i + vnbr, lcfx, (norm(0)) / dis1));
        tripletes.push_back(Trip(i + vnbr, lcfy, (norm(1)) / dis1));
        tripletes.push_back(Trip(i + vnbr, lcfz, (norm(2)) / dis1));

        tripletes.push_back(Trip(i + vnbr, lcbx, -(norm(0)) / dis1));
        tripletes.push_back(Trip(i + vnbr, lcby, -(norm(1)) / dis1));
        tripletes.push_back(Trip(i + vnbr, lcbz, -(norm(2)) / dis1));

        energy[i + vnbr] = norm.dot(Vcf - Vcb) / dis1;
        // if (isnan(energy[i + vnbr]))
        // {
        //     std::cout<<"NAN 2\n ";
        // }

        // norm * norm = 1
        tripletes.push_back(Trip(i + vnbr * 2, lnx, 2 * norm(0)));
        tripletes.push_back(Trip(i + vnbr * 2, lny, 2 * norm(1)));
        tripletes.push_back(Trip(i + vnbr * 2, lnz, 2 * norm(2)));

        energy[i + vnbr * 2] = norm.dot(norm) - 1;
        // if (isnan(energy[i + vnbr * 2]))
        // {
        //     std::cout<<"NAN 3\n ";
        // }

        
    }
    Enm0 = vec_list_to_matrix(edge_nm0);
    Enm1 = vec_list_to_matrix(edge_nm1);
    spMat J;
    J.resize(energy.size(), GlobVars.size());
    J.setFromTriplets(tripletes.begin(), tripletes.end());
    H = J.transpose() * J;
    B = -J.transpose() * energy;
}

void QuadOpt::opt(){
    if(Ftri.rows()==0){
        std::cout<<"The AABB Tree is NOT Loaded. Please Load the Tree"<<std::endl;
        return;
    }
    int vnbr = V.rows();
    if (GlobVars.size() != varsize)
    {
        ComputeAuxiliaries = true;
        GlobVars = Eigen::VectorXd::Zero(varsize);
        GlobVars.segment(0, vnbr) = V.col(0);
        GlobVars.segment(vnbr, vnbr) = V.col(1);
        GlobVars.segment(vnbr * 2, vnbr) = V.col(2);
    }
    spMat H;
    H.resize(varsize, varsize);
    Eigen::VectorXd B = Eigen::VectorXd::Zero(varsize);
    Eigen::VectorXd Egravity, Esmth, Enorm, Ebnm0, Ebnm1, Ebnm2, Epg[4];

    spMat Hgravity;  // approximation
	Eigen::VectorXd Bgravity; // right of laplacian
	assemble_gravity_ApproOriginal(Hgravity, Bgravity, Egravity);
    H += weight_gravity * Hgravity;
	B += weight_gravity * Bgravity;

    spMat Hsmth;
    Eigen::VectorXd Bsmth;
	assemble_fairness(Hsmth, Bsmth, Esmth);
    H += weight_fairness * Hsmth;
	B += weight_fairness * Bsmth;

    spMat Hnorm;
    Eigen::VectorXd Bnorm;
    assemble_normal_conditions(Hnorm, Bnorm, Enorm);
    H += weight_pg * Hnorm;
	B += weight_pg * Bnorm;

    if (OptType == 0)
    { // AAG, the diagonal is geodesic that requires binormals
        spMat Hbnm;
        Eigen::VectorXd Bbnm;
        int family = -1;
        if (d0_type > 0)
        {
            family = 2;
        }
        if (d1_type > 0)
        {
            family = 3;
        }
        // vertices vnbr * 3, normals vnbr * 3, binormals for G vnbr * 3  
        int bnm_start = vnbr * 6;
        assemble_binormal_conditions(Hbnm, Bbnm, Ebnm0, family, bnm_start);
        H += weight_pg * Hbnm;
        B += weight_pg * Bbnm;


        spMat Hpg[3];
        Eigen::VectorXd Bpg[3];
        // G
        int type = 2; // G
        assemble_pg_extreme_cases(Hpg[0], Bpg[0], Epg[0], type, family, bnm_start);

        // A
        type = 1; // A0
        family = 0;// row
        bnm_start = -1;
        assemble_pg_extreme_cases(Hpg[1], Bpg[1], Epg[1], type, family, bnm_start);

        type = 1; // A1
        family = 1;// col
        bnm_start = -1;
        assemble_pg_extreme_cases(Hpg[2], Bpg[2], Epg[2], type, family, bnm_start);

        H += weight_pg * (pg_ratio * Hpg[0] + Hpg[1] + Hpg[2]);
        B += weight_pg * (pg_ratio * Bpg[0] + Bpg[1] + Bpg[2]);
    }

    if (OptType == 1) 
    { // GGA, the diagonal is asymptotic
        spMat Hbnm[2];
        Eigen::VectorXd Bbnm[2];
        int family = -1;
        family = 0;// row
        // vertices vnbr * 3, normals vnbr * 3, binormals for G0 vnbr * 3, binormals for G1 vnbr * 3    
        int bnm_start = vnbr * 6;
        assemble_binormal_conditions(Hbnm[0], Bbnm[0], Ebnm0, family, bnm_start);

        family = 1; // column
        bnm_start = vnbr * 9;
        assemble_binormal_conditions(Hbnm[1], Bbnm[1], Ebnm1, family, bnm_start);

        H += weight_pg * (Hbnm[0] + Hbnm[1]);
        B += weight_pg * (Bbnm[0] + Bbnm[1]);
        
        if (d0_type > 0)
        {
            family = 2;
        }
        if (d1_type > 0)
        {
            family = 3;
        }

        spMat Hpg[3];
        Eigen::VectorXd Bpg[3];
        // A
        int type = 1; // A
        bnm_start = -1;
        assemble_pg_extreme_cases(Hpg[0], Bpg[0], Epg[0], type, family, bnm_start);

        // G
        type = 2; // G0
        family = 0;// row
        bnm_start = vnbr * 6;
        assemble_pg_extreme_cases(Hpg[1], Bpg[1], Epg[1], type, family, bnm_start);

        type = 2; // G1
        family = 1;// col
        bnm_start = vnbr * 9;
        assemble_pg_extreme_cases(Hpg[2], Bpg[2], Epg[2], type, family, bnm_start);

        H += weight_pg * (Hpg[0] + Hpg[1] + Hpg[2]);
        B += weight_pg * (Bpg[0] + Bpg[1] + Bpg[2]);
    }
    if (OptType == 2 || OptType == 3) 
    { // PP or PPG
        spMat Hbnm[3];
        Eigen::VectorXd Bbnm[3];
        int family = -1;
        family = 0;// row
        // vertices vnbr * 3, normals vnbr * 3, binormals for P0 vnbr * 3, u for P0 vnbr * 3, binormals for P1 vnbr * 3, u for P1 vnbr * 3    
        int bnm_start = vnbr * 6;
        assemble_binormal_conditions(Hbnm[0], Bbnm[0], Ebnm0, family, bnm_start);
        family = 1; // column
        bnm_start = vnbr * 12;
        assemble_binormal_conditions(Hbnm[1], Bbnm[1], Ebnm1, family, bnm_start);
        H += weight_pg * pg_ratio * (Hbnm[0] + Hbnm[1]);
        B += weight_pg * pg_ratio * (Bbnm[0] + Bbnm[1]);


        spMat Hpg[3];
        Eigen::VectorXd Bpg[3];
        // P0
        bnm_start = vnbr * 6;
        family = 0;
        double angle_radian = angle_degree0 * LSC_PI / 180.;
        assemble_pg_cases(angle_radian, Hpg[0], Bpg[0], Epg[0], family, bnm_start);
        // P1
        bnm_start = vnbr * 12;
        family = 1;
        angle_radian = angle_degree1 * LSC_PI / 180.;
        assemble_pg_cases(angle_radian, Hpg[1], Bpg[1], Epg[1], family, bnm_start);
        // H += weight_pg * pg_ratio * (Hpg[1]);
        // B += weight_pg * pg_ratio * (Bpg[1]);
        H += weight_pg * (Hpg[0] + Hpg[1]);
        B += weight_pg * (Bpg[0] + Bpg[1]);

        Estimate_PG_Angles = false;

        
        if (OptType == 3)// G as diagonal 
        {
            if (WhichDiagonal == 0)
            {
                family = 2;
            }
            if (WhichDiagonal == 1)
            {
                family = 3;
            }

            bnm_start = vnbr * 18;
            assemble_binormal_conditions(Hbnm[2], Bbnm[2], Ebnm2, family, bnm_start);
            H += weight_pg * pg_ratio * Hbnm[2];
            B += weight_pg * pg_ratio * Bbnm[2];
            // std::cout<<"check 5"<<std::endl;

            int type = 2; 
            assemble_pg_extreme_cases(Hpg[2], Bpg[2], Epg[2], type, family, bnm_start);
            // std::cout<<"check 6"<<std::endl;
            H += weight_pg * Hpg[2];
            B += weight_pg * Bpg[2];
        }
        
    }
    if (OptType == 4)
    { // AAGG, the diagonal is geodesic that requires binormals
        spMat Hbnm, Hbnm1;
        Eigen::VectorXd Bbnm, Bbnm1;
        // int family = -1;
        // if (d0_type > 0)
        // {
        //     family = 2;
        // }
        // if (d1_type > 0)
        // {
        //     family = 3;
        // }
        // vertices vnbr * 3, normals vnbr * 3, binormals for G0 vnbr * 3, binormals for G1 vnbr * 3  
        int bnm_start0 = vnbr * 6;
        assemble_binormal_conditions(Hbnm, Bbnm, Ebnm0, 2, bnm_start0);
        H += weight_pg * Hbnm;
        B += weight_pg * Bbnm;

        int bnm_start1 = vnbr * 9;
        assemble_binormal_conditions(Hbnm1, Bbnm1, Ebnm1, 3, bnm_start1);
        H += weight_pg * Hbnm1;
        B += weight_pg * Bbnm1;



        spMat Hpg[4];
        Eigen::VectorXd Bpg[4];
        // G
        int type = 2; // G
        assemble_pg_extreme_cases(Hpg[0], Bpg[0], Epg[0], type, 2, bnm_start0);
        // G another
        assemble_pg_extreme_cases(Hpg[1], Bpg[1], Epg[1], type, 3, bnm_start1);

        // A
        type = 1; // A0
        int family = 0;// row
        int bnm_start = -1;
        assemble_pg_extreme_cases(Hpg[2], Bpg[2], Epg[2], type, family, bnm_start);

        type = 1; // A1
        family = 1;// col
        bnm_start = -1;
        assemble_pg_extreme_cases(Hpg[3], Bpg[3], Epg[3], type, family, bnm_start);

        H += weight_pg * (pg_ratio * Hpg[0] + Hpg[1] + Hpg[2] + Hpg[3]);
        B += weight_pg * (pg_ratio * Bpg[0] + Bpg[1] + Bpg[2] + Bpg[3]);
    }

    // assemble together

    H += 1e-6 * (weight_mass * gravity_matrix + spMat(Eigen::VectorXd::Ones(varsize).asDiagonal()));
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver(H);

	// assert(solver.info() == Eigen::Success);
	if (solver.info() != Eigen::Success)
	{
		// solving failed
		std::cout << "solver fail" << std::endl;
		return;
	}
	// std::cout<<"solved successfully"<<std::endl;
    Eigen::VectorXd dx = solver.solve(B).eval();
    dx *= 0.75;
    double step_length = dx.norm();
    if (step_length > max_step)
    {
        dx *= max_step / step_length;
    }
    GlobVars += dx;
    if(OptType == 0){
        std::cout<<"AAG, ";
    }
    if(OptType == 1){
        std::cout<<"GGA, ";
    }
    if(OptType == 2){
        std::cout<<"PP, ";
    }
    if(OptType == 3){
        std::cout<<"PPG, ";
    }
    if(OptType == 4){
        std::cout<<"AAGG, ";
    }

    double ev_appro = Egravity.norm();
    double ev_smt = Esmth.norm();
    double ev_norm = Enorm.norm();
    std::cout << "Eclose, " << ev_appro << ", lap, " << ev_smt << ", normal vector, " << ev_norm << ", ";
    if (OptType == 0)
    {
        double ebi = Ebnm0.norm();

        std::cout << "bnm, " << ebi << ", GAA, " << Epg[0].norm() << ", " << Epg[1].norm() << ", " << Epg[2].norm() << ", ";
    }
    if (OptType == 1)
    {
        double ebi0 = Ebnm0.norm();
        double ebi1 = Ebnm1.norm();
        std::cout << "bnm, " << ebi0 << ", " << ebi1 << ", AGG, " << Epg[0].norm() << ", " << Epg[1].norm() << ", " << Epg[2].norm() << ", ";
    }
    if (OptType == 2 || OptType == 3)
    {
        double ebi0 = Ebnm0.norm();
        double ebi1 = Ebnm1.norm();

        std::cout << "bnm, " << ebi0 << ", " << ebi1;
        if (OptType == 3)
        {
            double ebi2 = Ebnm2.norm();
            std::cout << ", " << ebi2 << ", PPG, ";
        }
        else
        {
            std::cout << ", PP, ";
        }
        std::cout << Epg[0].norm() << ", " << Epg[1].norm();
        if (OptType == 3)
        {
            std::cout << ", " <<Epg[2].norm();
        }
    }
    if (OptType == 4)
    {
        std::cout << "bnm, " << Ebnm0.norm()<<", "<< Ebnm1.norm()<< ", GGAA, " << Epg[0].norm() << ", " << Epg[1].norm() << ", " << Epg[2].norm() << ", "<<Epg[3].norm() ;
    }
    real_step_length = dx.norm();
    std::cout << ", stp, " << dx.norm() << ", diagonal types, "<<d0_type<<", "<<d1_type<<", ";
    std::cout << "\n";

    // convert the data to the mesh format
    ComputeAuxiliaries = false;
    V.col(0) = GlobVars.segment(0, vnbr);
    V.col(1) = GlobVars.segment(vnbr, vnbr);
    V.col(2) = GlobVars.segment(vnbr * 2, vnbr);
    for (CGMesh::VertexIter v_it = mesh_update.vertices_begin(); v_it != mesh_update.vertices_end(); ++v_it)
    {
        int vid = v_it.handle().idx();
        mesh_update.point(*v_it) = CGMesh::Point(V(vid, 0), V(vid, 1), V(vid, 2));
    }
}

void QuadOpt::optAAG()
{
    if (WhichDiagonal == 0)
    {
        d0_type = 1;
    }
    if (WhichDiagonal == 1)
    {
        d1_type = 1;
    }
    int vnbr = V.rows();
    int order = 1; // the order = 1 means AAG
    if (GlobVars.size() != varsize)
    {
        int sizediff = (varsize - GlobVars.size()) / 9;
        std::cout << "Assigning variables based on previous ones, there are " << sizediff << " vertices more\n";
        ComputeAuxiliaries = true;
        Eigen::VectorXd tmpVar = Eigen::VectorXd::Zero(varsize);
        if (GlobVars.size() > 0)
            tmpVar.segment(0, GlobVars.size()) = GlobVars;
        GlobVars = tmpVar;
        for (int i = vnbr - sizediff; i < vnbr; i++)
        {
            GlobVars[i * 9] = V(i, 0);
            GlobVars[i * 9 + 1] = V(i, 1);
            GlobVars[i * 9 + 2] = V(i, 2);
        }
    }
    spMat H;
    H.resize(varsize, varsize);
    Eigen::VectorXd B = Eigen::VectorXd::Zero(varsize);
    Eigen::VectorXd Egravity, Esmth, Enorm, Ebnm0, Ebnm1, Ebnm2, Epg[4];

    spMat Hgravity;  // approximation
	Eigen::VectorXd Bgravity; // right of laplacian
	assemble_gravity_AAG_AGG(Hgravity, Bgravity, Egravity);
    H += weight_gravity * Hgravity;
	B += weight_gravity * Bgravity;

    spMat HCurve;  // approximation to curve
	Eigen::VectorXd BCurve, ECurve; // right of approximation to curve
    assemble_approximate_curve_conditions(HCurve, BCurve, ECurve, order);
    H += weight_curve * HCurve;
	B += weight_curve * BCurve;
    
    spMat Hsmth;
    Eigen::VectorXd Bsmth;
	assemble_fairness(Hsmth, Bsmth, Esmth, order);
    H += weight_fairness * Hsmth;
	B += weight_fairness * Bsmth;
    spMat Hnorm;
    Eigen::VectorXd Bnorm;
    assemble_normal_conditions(Hnorm, Bnorm, Enorm, order); // the 
    H += weight_pg * Hnorm;
	B += weight_pg * Bnorm;
    // AAG, the row is the geodesic, the col is the A, and the chosen diag is another A
    spMat Hbnm;
    Eigen::VectorXd Bbnm;
    int family = -1;
    family = 0; // row is the geodesic
    // vertices vnbr * 3, normals vnbr * 3, binormals for G vnbr * 3
    assemble_binormal_conditions(Hbnm, Bbnm, Ebnm0, family, -1, order);
    H += weight_pg * Hbnm;
    B += weight_pg * Bbnm;
    // G
    int type = 2; // G
    spMat Hpg[3];
    Eigen::VectorXd Bpg[3];
    assemble_pg_extreme_cases(Hpg[0], Bpg[0], Epg[0], type, family, -1, order);

    // A
    type = 1;   // A0
    family = 1; // col
    assemble_pg_extreme_cases(Hpg[1], Bpg[1], Epg[1], type, family, -1, order);

    // A, the diagonal is another A
    if (d0_type > 0)
    {
        family = 2;
    }
    if (d1_type > 0)
    {
        family = 3;
    }
    type = 1;   // A1
    assemble_pg_extreme_cases(Hpg[2], Bpg[2], Epg[2], type, family, -1, order);

    H += weight_pg * pg_ratio * (Hpg[0] + Hpg[1] + Hpg[2]);
    B += weight_pg * pg_ratio * (Bpg[0] + Bpg[1] + Bpg[2]);

    // assemble together

    H += 1e-6 * (weight_mass * gravity_matrix + spMat(Eigen::VectorXd::Ones(varsize).asDiagonal()));
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver(H);

	// assert(solver.info() == Eigen::Success);
	if (solver.info() != Eigen::Success)
	{
		// solving failed
		std::cout << "solver fail" << std::endl;
		return;
	}
	// std::cout<<"solved successfully"<<std::endl;
    Eigen::VectorXd dx = solver.solve(B).eval();
    dx *= 0.75;
    double step_length = dx.norm();
    if (step_length > max_step)
    {
        dx *= max_step / step_length;
    }
    GlobVars += dx;

    std::cout << "AAG, ";

    double ev_appro = Egravity.norm();
    double ev_smt = Esmth.norm();
    double ev_norm = Enorm.norm();
    std::cout << "Eclose, " << ev_appro<<", Ecurve, "<<ECurve.norm() << ", lap, " << ev_smt << ", normal vector, " << ev_norm << ", ";

    double ebi = Ebnm0.norm();
    std::cout << "bnm, " << ebi << ", GAA, " << Epg[0].norm() << ", " << Epg[1].norm() << ", " << Epg[2].norm() << ", ";

    real_step_length = dx.norm();
    std::cout << ", stp, " << dx.norm() << ", diagonal types, "<<d0_type<<", "<<d1_type<<", ";
    std::cout << "\n";
    // std::cout<<"here 1\n";
    // convert the data to the mesh format
    ComputeAuxiliaries = false;
    for (int i = 0; i < vnbr; i++)
    {
        V(i, 0) = GlobVars[i * 9];
        V(i, 1) = GlobVars[i * 9 + 1];
        V(i, 2) = GlobVars[i * 9 + 2];
    }
    // std::cout<<"here 2\nver nbrs "<<mesh_update.n_vertices()<<"\n";
    for (CGMesh::VertexIter v_it = mesh_update.vertices_begin(); v_it != mesh_update.vertices_end(); ++v_it)
    {
        int vid = v_it.handle().idx();
        if (vid >= V.rows() || vid < 0)
        {
            std::cout << "Error in vid " << vid << std::endl;
        }
        assert(vid < V.rows());
        assert(vid >= 0);
        mesh_update.point(*v_it) = CGMesh::Point(V(vid, 0), V(vid, 1), V(vid, 2));
    }
    // std::cout<<"opt finished\n";
}

void collectAggVectors(const Eigen::VectorXd &vars, Eigen::MatrixXd &n,
                       Eigen::MatrixXd &b0,
                       Eigen::MatrixXd &b1)
{
    int vnbr = vars.size() / 12;
    n.resize(vnbr,3);
    b0 = b1 = n;
    for (int i = 0; i < vnbr; i++)
    {
        Eigen::Vector3d p(vars[i * 12],vars[i * 12 + 1],vars[i * 12 + 2]);
        Eigen::Vector3d nvec(vars[i * 12 + 3],vars[i * 12 + 4],vars[i * 12 + 5]);
        Eigen::Vector3d b0vec(vars[i * 12 + 6],vars[i * 12 + 7],vars[i * 12 + 8]);
        Eigen::Vector3d b1vec(vars[i * 12 + 9],vars[i * 12 + 10],vars[i * 12 + 11]);
        n.row(i) = nvec;
        b0.row(i) = b0vec;
        b1.row(i) = b1vec;
    }
}

void QuadOpt::optAGG()
{
    if (WhichDiagonal == 0)
    {
        d0_type = 1;
    }
    if (WhichDiagonal == 1)
    {
        d1_type = 1;
    }
    int vnbr = V.rows();
    int order = 2; // the order = 1 means AAG, order = 2 is for AGG
    int varOffset = 12; // AGG has 12 * n variables
    if (GlobVars.size() != varsize)
    {
        int sizediff = (varsize - GlobVars.size()) / varOffset;
        std::cout << "Assigning variables based on previous ones, there are " << sizediff << " vertices more\n";
        ComputeAuxiliaries = true;
        Eigen::VectorXd tmpVar = Eigen::VectorXd::Zero(varsize);
        if (GlobVars.size() > 0)
            tmpVar.segment(0, GlobVars.size()) = GlobVars;
        GlobVars = tmpVar;
        for (int i = vnbr - sizediff; i < vnbr; i++)
        {
            GlobVars[i * varOffset] = V(i, 0);
            GlobVars[i * varOffset + 1] = V(i, 1);
            GlobVars[i * varOffset + 2] = V(i, 2);
        }
    }
    spMat H;
    H.resize(varsize, varsize);
    Eigen::VectorXd B = Eigen::VectorXd::Zero(varsize);
    Eigen::VectorXd Egravity, Esmth, Enorm, Ebnm0, Ebnm1, Ebnm2, Epg[4];

    spMat Hgravity;  // approximation
	Eigen::VectorXd Bgravity; // right of approximation
	assemble_gravity_AAG_AGG(Hgravity, Bgravity, Egravity, order);
    H += weight_gravity * Hgravity;
	B += weight_gravity * Bgravity;
    spMat HCurve;  // approximation to curve
	Eigen::VectorXd BCurve, ECurve; // right of approximation to curve
    assemble_approximate_curve_conditions(HCurve, BCurve, ECurve, order);
    H += weight_curve * HCurve;
	B += weight_curve * BCurve;

    spMat Hsmth;
    Eigen::VectorXd Bsmth;
	assemble_fairness(Hsmth, Bsmth, Esmth, order);
    H += weight_fairness * Hsmth;
	B += weight_fairness * Bsmth;
    spMat Hnorm;
    Eigen::VectorXd Bnorm;
    assemble_normal_conditions(Hnorm, Bnorm, Enorm, order); // the normal vectors
    H += weight_pg * Hnorm;
	B += weight_pg * Bnorm;
    // AAG, the row is the geodesic, the col is the A, and the chosen diag is another A
    spMat Hbnm0, Hbnm1;
    Eigen::VectorXd Bbnm0, Bbnm1;
    int family = -1;
    
    family = 1; // column is the geodesic
    assemble_binormal_conditions(Hbnm0, Bbnm0, Ebnm0, family, -1, order, 0); // 1 means the second binormal vector
    H += weight_pg * Hbnm0;
    B += weight_pg * Bbnm0;

    // G
    int type = 2; // G1
    family = 1; // column is the geodesic
    spMat Hpg[3];
    Eigen::VectorXd Bpg[3];
    assemble_pg_extreme_cases(Hpg[0], Bpg[0], Epg[0], type, family, -1, order, 0);

    // A
    type = 1;   // A
    family = 0; // row is the asymptotic
    assemble_pg_extreme_cases(Hpg[1], Bpg[1], Epg[1], type, family, -1, order);

    // A, the diagonal is another A
    if (d0_type > 0)
    {
        family = 2;
    }
    if (d1_type > 0)
    {
        family = 3;
    }
    type = 2;   // G2
    assemble_pg_extreme_cases(Hpg[2], Bpg[2], Epg[2], type, family, -1, order, 1);

    // another binormal
    assemble_binormal_conditions(Hbnm1, Bbnm1, Ebnm1, family, -1, order, 1);// 0 means the first binormal vector
    H += weight_pg * Hbnm1;
    B += weight_pg * Bbnm1;

    H += weight_pg * pg_ratio * (Hpg[0] + Hpg[1] + Hpg[2]);
    B += weight_pg * pg_ratio * (Bpg[0] + Bpg[1] + Bpg[2]);

    // assemble together

    H += 1e-6 * (weight_mass * gravity_matrix + spMat(Eigen::VectorXd::Ones(varsize).asDiagonal()));
    Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver(H);

	// assert(solver.info() == Eigen::Success);
	if (solver.info() != Eigen::Success)
	{
		// solving failed
		std::cout << "solver fail" << std::endl;
		return;
	}
	// std::cout<<"solved successfully"<<std::endl;
    Eigen::VectorXd dx = solver.solve(B).eval();
    dx *= 0.75;
    double step_length = dx.norm();
    if (step_length > max_step)
    {
        dx *= max_step / step_length;
    }
    GlobVars += dx;

    std::cout << "AGG, ";

    double ev_appro = Egravity.norm();
    double ev_smt = Esmth.norm();
    double ev_norm = Enorm.norm();
    std::cout << "Eclose, " << ev_appro<<", Ecurve, "<<ECurve.norm() << ", lap, " << ev_smt << ", normal vector, " << ev_norm << ", ";

    double ebi0 = Ebnm0.norm();
    double ebi1 = Ebnm1.norm();
    ebi0 *= ebi0;
    ebi1 *= ebi1;
    double ebi = sqrt(ebi0 + ebi1);
    std::cout << "bnm, " << ebi << ", GAG, " << Epg[0].norm() << ", " << Epg[1].norm() << ", " << Epg[2].norm() << ", ";

    real_step_length = dx.norm();
    std::cout << ", stp, " << dx.norm() << ", diagonal types, "<<d0_type<<", "<<d1_type<<", ";
    std::cout << "\n";
    // std::cout<<"here 1\n";
    // convert the data to the mesh format
    ComputeAuxiliaries = false;
    for (int i = 0; i < vnbr; i++)
    {
        V(i, 0) = GlobVars[i * varOffset];
        V(i, 1) = GlobVars[i * varOffset + 1];
        V(i, 2) = GlobVars[i * varOffset + 2];
    }
    // std::cout<<"here 2\nver nbrs "<<mesh_update.n_vertices()<<"\n";
    for (CGMesh::VertexIter v_it = mesh_update.vertices_begin(); v_it != mesh_update.vertices_end(); ++v_it)
    {
        int vid = v_it.handle().idx();
        if (vid >= V.rows() || vid < 0)
        {
            std::cout << "Error in vid " << vid << std::endl;
        }
        assert(vid < V.rows());
        assert(vid >= 0);
        mesh_update.point(*v_it) = CGMesh::Point(V(vid, 0), V(vid, 1), V(vid, 2));
    }
    collectAggVectors(GlobVars, N, B0, B1);
    // std::cout<<"opt finished\n";
}
void QuadOpt::optGGG()
{
	if (WhichDiagonal == 0)
	{
		d0_type = 1;
	}
	if (WhichDiagonal == 1)
	{
		d1_type = 1;
	}
	int vnbr = V.rows();
	int order = 3; // the order = 1 means AAG, order = 2 is for AGG, order = 3 is for GGG.
	int varOffset = 15; // GGG has 15 * n variables
	if (GlobVars.size() != varsize)
	{
		int sizediff = (varsize - GlobVars.size()) / varOffset;
		std::cout << "Assigning variables based on previous ones, there are " << sizediff << " vertices more\n";
		ComputeAuxiliaries = true;
		Eigen::VectorXd tmpVar = Eigen::VectorXd::Zero(varsize);
		if (GlobVars.size() > 0)
			tmpVar.segment(0, GlobVars.size()) = GlobVars;
		GlobVars = tmpVar;
		for (int i = vnbr - sizediff; i < vnbr; i++)
		{
			GlobVars[i * varOffset] = V(i, 0);
			GlobVars[i * varOffset + 1] = V(i, 1);
			GlobVars[i * varOffset + 2] = V(i, 2);
		}
	}
	spMat H;
	H.resize(varsize, varsize);
	Eigen::VectorXd B = Eigen::VectorXd::Zero(varsize);
	Eigen::VectorXd Egravity, Esmth, Enorm, Ebnm0, Ebnm1, Ebnm2, Epg[4];

	spMat Hgravity;  // approximation
	Eigen::VectorXd Bgravity; // right of approximation
	assemble_gravity_AAG_AGG(Hgravity, Bgravity, Egravity, order);
	H += weight_gravity * Hgravity;
	B += weight_gravity * Bgravity;
	spMat HCurve;  // approximation to curve
	Eigen::VectorXd BCurve, ECurve; // right of approximation to curve
	assemble_approximate_curve_conditions(HCurve, BCurve, ECurve, order);
	H += weight_curve * HCurve;
	B += weight_curve * BCurve;

	spMat Hsmth;
	Eigen::VectorXd Bsmth;
	assemble_fairness(Hsmth, Bsmth, Esmth, order);
	H += weight_fairness * Hsmth;
	B += weight_fairness * Bsmth;
	spMat Hnorm;
	Eigen::VectorXd Bnorm;
	assemble_normal_conditions(Hnorm, Bnorm, Enorm, order); // the normal vectors
	H += weight_pg * Hnorm;
	B += weight_pg * Bnorm;
	// GGG constraints
	spMat Hbnm0, Hbnm1, Hbnm2;
	Eigen::VectorXd Bbnm0, Bbnm1, Bbnm2;
	spMat Hpg[3];
	Eigen::VectorXd Bpg[3];
	int family = -1;
	int type = 2; // Geodesic

	family = 0; // row is the geodesic
	assemble_binormal_conditions(Hbnm0, Bbnm0, Ebnm0, family, -1, order, 0); 
	H += weight_pg * Hbnm0;
	B += weight_pg * Bbnm0;
	// G
	assemble_pg_extreme_cases(Hpg[0], Bpg[0], Epg[0], type, family, -1, order, 0);

	// G
	family = 1; // colomn is the geodesic
	assemble_binormal_conditions(Hbnm1, Bbnm1, Ebnm1, family, -1, order, 1); // 1 means the second binormal vector
	H += weight_pg * Hbnm1;
	B += weight_pg * Bbnm1;
	assemble_pg_extreme_cases(Hpg[1], Bpg[1], Epg[1], type, family, -1, order, 1);

	// A, the diagonal is another A
	if (d0_type > 0)
	{
		family = 2;
	}
	if (d1_type > 0)
	{
		family = 3;
	}
	// another binormal
	assemble_binormal_conditions(Hbnm2, Bbnm2, Ebnm2, family, -1, order, 2);// 0 means the first binormal vector
	H += weight_pg * Hbnm2;
	B += weight_pg * Bbnm2;
	assemble_pg_extreme_cases(Hpg[2], Bpg[2], Epg[2], type, family, -1, order, 2);

	

	H += weight_pg * pg_ratio * (Hpg[0] + Hpg[1] + Hpg[2]);
	B += weight_pg * pg_ratio * (Bpg[0] + Bpg[1] + Bpg[2]);

	// assemble together

	H += 1e-6 * (weight_mass * gravity_matrix + spMat(Eigen::VectorXd::Ones(varsize).asDiagonal()));
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver(H);

	// assert(solver.info() == Eigen::Success);
	if (solver.info() != Eigen::Success)
	{
		// solving failed
		std::cout << "solver fail" << std::endl;
		return;
	}
	// std::cout<<"solved successfully"<<std::endl;
	Eigen::VectorXd dx = solver.solve(B).eval();
	dx *= 0.75;
	double step_length = dx.norm();
	if (step_length > max_step)
	{
		dx *= max_step / step_length;
	}
	GlobVars += dx;

	std::cout << "GGG, ";

	double ev_appro = Egravity.norm();
	double ev_smt = Esmth.norm();
	double ev_norm = Enorm.norm();
	std::cout << "Eclose, " << ev_appro << ", Ecurve, " << ECurve.norm() << ", lap, " << ev_smt << ", normal vector, " << ev_norm << ", ";

	double ebi0 = Ebnm0.norm();
	double ebi1 = Ebnm1.norm();
	double ebi2 = Ebnm2.norm();
	ebi0 *= ebi0;
	ebi1 *= ebi1;
	ebi2 *= ebi2;
	double ebi = sqrt(ebi0 + ebi1 + ebi2);
	std::cout << "bnm, " << ebi << ", GGG, " << Epg[0].norm() << ", " << Epg[1].norm() << ", " << Epg[2].norm() << ", ";

	real_step_length = dx.norm();
	std::cout << ", stp, " << dx.norm() << ", diagonal types, " << d0_type << ", " << d1_type << ", ";
	std::cout << "\n";
	// std::cout<<"here 1\n";
	// convert the data to the mesh format
	ComputeAuxiliaries = false;
	for (int i = 0; i < vnbr; i++)
	{
		V(i, 0) = GlobVars[i * varOffset];
		V(i, 1) = GlobVars[i * varOffset + 1];
		V(i, 2) = GlobVars[i * varOffset + 2];
	}
	// std::cout<<"here 2\nver nbrs "<<mesh_update.n_vertices()<<"\n";
	for (CGMesh::VertexIter v_it = mesh_update.vertices_begin(); v_it != mesh_update.vertices_end(); ++v_it)
	{
		int vid = v_it.handle().idx();
		if (vid >= V.rows() || vid < 0)
		{
			std::cout << "Error in vid " << vid << std::endl;
		}
		assert(vid < V.rows());
		assert(vid >= 0);
		mesh_update.point(*v_it) = CGMesh::Point(V(vid, 0), V(vid, 1), V(vid, 2));
	}
}
void QuadOpt::optISO()
{
	int vnbr = V.rows();
	int order = 4; // the order = 1 means AAG, order = 2 is for AGG, order = 3 is for GGG, order = 4 is for isometric deformation.
	int varOffset = 3; // isometric deformation has 3 * n variables
	if (GlobVars.size() != varsize)
	{
		int sizediff = (varsize - GlobVars.size()) / varOffset;
		std::cout << "Assigning variables based on previous ones, there are " << sizediff << " vertices more\n";
		ComputeAuxiliaries = true;
		Eigen::VectorXd tmpVar = Eigen::VectorXd::Zero(varsize);
		if (GlobVars.size() > 0)
			tmpVar.segment(0, GlobVars.size()) = GlobVars;
		GlobVars = tmpVar;
		for (int i = vnbr - sizediff; i < vnbr; i++)
		{
			GlobVars[i * varOffset] = V(i, 0);
			GlobVars[i * varOffset + 1] = V(i, 1);
			GlobVars[i * varOffset + 2] = V(i, 2);
		}
	}
	spMat H, Hiso;
	H.resize(varsize, varsize);
	Eigen::VectorXd B = Eigen::VectorXd::Zero(varsize), Biso;
	Eigen::VectorXd Egravity, Esmth, Eiso;

	spMat Hgravity;  // approximation
	Eigen::VectorXd Bgravity; // right of approximation
	assemble_gravity_AAG_AGG(Hgravity, Bgravity, Egravity, order);
	H += weight_gravity * Hgravity;
	B += weight_gravity * Bgravity;

	spMat Hsmth;
	Eigen::VectorXd Bsmth;
	assemble_fairness(Hsmth, Bsmth, Esmth, order);
	H += weight_fairness * Hsmth;
	B += weight_fairness * Bsmth;
	
	assemble_isometric_condition(Hiso, Biso, Eiso);

	H += weight_pg * Hiso;
	B += weight_pg * Biso;

	// assemble together

	H += 1e-6 * (weight_mass * gravity_matrix + spMat(Eigen::VectorXd::Ones(varsize).asDiagonal()));
	Eigen::SimplicialLLT<Eigen::SparseMatrix<double>> solver(H);

	// assert(solver.info() == Eigen::Success);
	if (solver.info() != Eigen::Success)
	{
		// solving failed
		std::cout << "solver fail" << std::endl;
		return;
	}
	// std::cout<<"solved successfully"<<std::endl;
	Eigen::VectorXd dx = solver.solve(B).eval();
	dx *= 0.75;
	double step_length = dx.norm();
	if (step_length > max_step)
	{
		dx *= max_step / step_length;
	}
	GlobVars += dx;

	std::cout << "ISO, ";

	double ev_appro = Egravity.norm();
	double ev_smt = Esmth.norm();
	double ev_iso = Eiso.norm();
	std::cout << "Eclose, " << ev_appro << ", lap, " << ev_smt << ", isometric, " << ev_iso;

	real_step_length = dx.norm();
	std::cout << ", stp, " << real_step_length;
	std::cout << "\n";
	// std::cout<<"here 1\n";
	// convert the data to the mesh format
	ComputeAuxiliaries = false;
	for (int i = 0; i < vnbr; i++)
	{
		V(i, 0) = GlobVars[i * varOffset];
		V(i, 1) = GlobVars[i * varOffset + 1];
		V(i, 2) = GlobVars[i * varOffset + 2];
	}
	// std::cout<<"here 2\nver nbrs "<<mesh_update.n_vertices()<<"\n";
	for (CGMesh::VertexIter v_it = mesh_update.vertices_begin(); v_it != mesh_update.vertices_end(); ++v_it)
	{
		int vid = v_it.handle().idx();
		if (vid >= V.rows() || vid < 0)
		{
			std::cout << "Error in vid " << vid << std::endl;
		}
		assert(vid < V.rows());
		assert(vid >= 0);
		mesh_update.point(*v_it) = CGMesh::Point(V(vid, 0), V(vid, 1), V(vid, 2));
	}
}

Eigen::MatrixXd QuadOpt::propagateBoundary()
{
    Eigen::MatrixXd result(vNbrInRow, 3);
    int vnbr = V.rows();
    int counter = 0;
    for (int i = vnbr - vNbrInRow; i < vnbr; i++)
    {
        result.row(counter) = V.row(i);
        counter++;
    }
    return result;
}

void QuadOpt::show_curve_families(std::array<Eigen::MatrixXd, 3>& edges){
    // int row_m = rowinfo.size() / 2;
    // int col_m = rowinfo[row_m].size() / 2;

    // int vid = rowinfo[row_m][col_m]; // choose a point in the middle
    int vid = mesh_original.n_vertices() / 3;
    std::cout<<"rowback\n"<<row_back<<"\n";
    std::cout<<"colback\n"<<col_back<<"\n";
    std::vector<Eigen::Vector3d> vtp;
    int vcurrent = vid;
    for (int i = 0;; i++)// rows
    {
        if (row_back[vcurrent] < 0)
        {
            break;
        }
        int vnext = row_back[vcurrent];
        vtp.push_back(V.row(vnext));
        vcurrent = vnext;
    }

    edges[0] = vec_list_to_matrix(vtp);

    vtp.clear();
    vcurrent = vid;
    for (int i = 0;; i++)// cols
    {
        if (col_back[vcurrent] < 0)
        {
            break;
        }
        int vnext = col_back[vcurrent];
        vtp.push_back(V.row(vnext));
        vcurrent = vnext;
    }
    edges[1] = vec_list_to_matrix(vtp);

    vtp.clear();
    vcurrent = vid;
    Eigen::VectorXi ajc;
    if (WhichDiagonal == 0)
    {
        ajc = d0_back;
    }
    if (WhichDiagonal == 1)
    {
        ajc = d1_back;
    }
    std::cout<<"diag_back\n"<<ajc<<"\n";
    for (int i = 0;; i++) // diagonals
    {
        if (ajc[vcurrent] < 0)
        {
            break;
        }
        int vnext = ajc[vcurrent];
        vtp.push_back(V.row(vnext));
        vcurrent = vnext;
    }
    edges[2] = vec_list_to_matrix(vtp);
    
}

void QuadOpt::show_diagonals(Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, Eigen::MatrixXd &E2, Eigen::MatrixXd &E3)
{
    // int row_m = rowinfo.size() / 2;
    // int col_m = rowinfo[row_m].size() / 2;

    // int vid = rowinfo[row_m][col_m]; // choose a point in the middle
    int vnbr = V.rows();
    Eigen::VectorXi fr = WhichDiagonal == 0 ? d0_front : d1_front;
    std::vector<Eigen::Vector3d> v0, v1, v2, v3;
    // Eigen::VectorXi bk = WhichDiagonal == 0 ? d0_back : d1_back;
    for (int i = 0; i < vnbr; i++)
    {
        if (fr[i] >= 0 && row_front(i) >= 0 && col_front(i) >= 0)
        {
            v0.push_back(V.row(i));
            v1.push_back(V.row(fr[i]));
            v2.push_back(V.row(row_front(i)));
            v3.push_back(V.row(col_front(i)));
        }
    }
    E0 = vec_list_to_matrix(v0);
    E1 = vec_list_to_matrix(v1);
    E2 = vec_list_to_matrix(v2);
    E3 = vec_list_to_matrix(v3);
    std::cout<<"rowfront "<<row_front.transpose()<<"\ncolfront "<<col_front.transpose()<<"\nd0front"<<d0_front.transpose()<<"\n";
}
void QuadOpt::get_Bnd(Eigen::VectorXi& Bnd){
    int vnbr = V.rows();
    Bnd = Eigen::VectorXi::Zero(vnbr);
    for(int i=0;i<vnbr;i++){
        int vid = i;
        int rf = row_front[vid];
        int rb = row_back[vid];
        int cf = col_front[vid];
        int cb = col_back[vid];

        if (rf == -1 || rb == -1 || cf == -1 || cb == -1)
        {
            Bnd[i] = 1;
        }
    }
}
void QuadOpt::extract_binormals(const int family, const int bnm_start, const int vid, Eigen::Vector3d &bi)
{
    int vnbr = V.rows();
    int rf = row_front[vid];
    int rb = row_back[vid];
    int cf = col_front[vid];
    int cb = col_back[vid];
    int d0f = d0_front[vid];
    int d0b = d0_back[vid];
    int d1f = d1_front[vid];
    int d1b = d1_back[vid];
    // the locations
    // the vertex
    int lvx = vid;
    int lvy = vid + vnbr;
    int lvz = vid + 2 * vnbr;
    // the front and back vertices
    int lfx, lfy, lfz, lbx, lby, lbz;
    // the binormal vector locations
    int lrx = bnm_start + vid;
    int lry = bnm_start + vid + vnbr;
    int lrz = bnm_start + vid + vnbr * 2;

    int lnx = vnbr * 3 + vid;
    int lny = vnbr * 4 + vid;
    int lnz = vnbr * 5 + vid;
    bool compute = true;
    if (family == 0) // rows
    {
        lfx = rf;
        lfy = rf + vnbr;
        lfz = rf + vnbr * 2;

        lbx = rb;
        lby = rb + vnbr;
        lbz = rb + vnbr * 2;
    }
    if (family == 1) // cols
    {
        lfx = cf;
        lfy = cf + vnbr;
        lfz = cf + vnbr * 2;

        lbx = cb;
        lby = cb + vnbr;
        lbz = cb + vnbr * 2;
    }
    if (family == 2) // d0
    {
        lfx = d0f;
        lfy = d0f + vnbr;
        lfz = d0f + vnbr * 2;

        lbx = d0b;
        lby = d0b + vnbr;
        lbz = d0b + vnbr * 2;
    }

    if (family == 3) // d1
    {

        lfx = d1f;
        lfy = d1f + vnbr;
        lfz = d1f + vnbr * 2;

        lbx = d1b;
        lby = d1b + vnbr;
        lbz = d1b + vnbr * 2;
    }
    if (lfx < 0 || lbx < 0)// this is an end point, we use the neighbouring binormal vector
    {
        compute = false;
    }
    if (OptType == 2 || OptType == 3) // if it is pseudo-geodesic, then use the optimized values
    {
        if (compute == false)
        { // the vertex is on the boundary
            bi = Eigen::Vector3d(0, 0, 0);
            return;
        }
    }
    if (lfx < 0)
    {
        lfx = lvx;
        lfy = lvx + vnbr;
        lfz = lvx + vnbr * 2;
    }
    if (lbx < 0)
    {
        lbx = lvx;
        lby = lvx + vnbr;
        lbz = lvx + vnbr * 2;
    }
    Eigen::Vector3d Ver(GlobVars[lvx], GlobVars[lvy], GlobVars[lvz]);
    Eigen::Vector3d Vf(GlobVars[lfx], GlobVars[lfy], GlobVars[lfz]);
    Eigen::Vector3d Vb(GlobVars[lbx], GlobVars[lby], GlobVars[lbz]);
    Eigen::Vector3d r(GlobVars[lrx], GlobVars[lry], GlobVars[lrz]);
    Eigen::Vector3d n(GlobVars[lnx], GlobVars[lny], GlobVars[lnz]);
    Eigen::Vector3d tangent = (Vf - Vb).normalized();
    if (OptType == 2 || OptType == 3)
    { // PP or PPG
        if (family == 0 || family == 1)// not the diagonal
        {
            bi = r;
        }
        else // the diagonal is G
        {
            bi = n.cross(tangent).normalized();
        }
    }
    if (OptType == 0)
    { // AAG
        if (family == 0 || family == 1) // A
        {
            bi = n;
        }
        else // G
        {
            bi = n.cross(tangent).normalized();
        }
    }
    if (OptType == 1) // GGA
    {
        if (family == 0 || family == 1) // G
        {
            bi = n.cross(tangent).normalized();
        }
        else // A
        {
            bi = n;
        }
    }
}
void QuadOpt::extract_binormals_propagation(const int family, const int vid, Eigen::Vector3d& bi, int varsOff)
{
    int vnbr = V.rows();
    int rf = row_front[vid];
    int rb = row_back[vid];
    int cf = col_front[vid];
    int cb = col_back[vid];
    int d0f = d0_front[vid];
    int d0b = d0_back[vid];
    int d1f = d1_front[vid];
    int d1b = d1_back[vid];
    // the locations
    // the vertex
    int lvx = vid * varsOff;
    int lvy = vid * varsOff + 1;
    int lvz = vid * varsOff + 2;
    // the front and back vertices
    int lfx, lfy, lfz, lbx, lby, lbz;
    // the binormal vector locations
    if (varsOff == 9)
    { // AAG
        if (family == 0) // g
        {
            bi << GlobVars[vid * varsOff + 6], GlobVars[vid * varsOff + 7], GlobVars[vid * varsOff + 8];
            return;
        }
        if(family == 1 || family == 2 || family == 3) // the column A or the diagonal A
        {
            bi << GlobVars[vid * varsOff + 3], GlobVars[vid * varsOff + 4], GlobVars[vid * varsOff + 5];
            return;
        }
    }
    if (varsOff == 12)
    { // AGG
        if (family == 0) // A
        {
            bi << GlobVars[vid * varsOff + 3], GlobVars[vid * varsOff + 4], GlobVars[vid * varsOff + 5];
            return;
        }
        if(family == 1) // column G
        {
            bi << GlobVars[vid * varsOff + 6], GlobVars[vid * varsOff + 7], GlobVars[vid * varsOff + 8];
            return;
        }
        if (family == 2 || family == 3) // the diagonal G
        {
            bi << GlobVars[vid * varsOff + 9], GlobVars[vid * varsOff + 10], GlobVars[vid * varsOff + 11];
            return;
        }
    }
}
void QuadOpt::extract_diagonals(const int family, std::vector<std::vector<int>> &digs){
    std::vector<std::vector<int>> curves;
    int vnbr = V.rows();
    // Eigen::VectorXi checked = Eigen::VectorXi::Zero(vnbr);
    // for(int i=0;i<vnbr; i++){
    //     if(checked[i]){
    //         continue;
    //     }
    //     for()
    // }
    Eigen::VectorXi Front, Back;
    if (family == 2)
    {
        Front = d0_front;
        Back = d0_back;
    }
    if (family == 3)
    {
        Front = d1_front;
        Back = d1_back;
    }

    Eigen::VectorXi starts = Eigen::VectorXi::Zero(vnbr);
    assert(family == 2 || family == 3);
    int nbr_start = 0;
    int nbr_wrongs = 0;
    std::vector<Eigen::Vector3d> wrongs;
    for (int i = 0; i < vnbr; i++) // get all the start points
    {

        int fr, bk;
        fr = Front[i];
        bk = Back[i];
        if (fr >= 0)
        {
            int ff;
            if (Front[fr] == -1)
            {
                starts[i] = 1;
                nbr_start++;
            }
        }
    }
    std::cout<<"Nbr of start points of the diagonals: "<<nbr_start<<std::endl;

    for (int itr = 0; itr < vnbr; itr++)
    {
        if (starts[itr] == 0)
        {
            continue;
        }
        std::vector<int> tmpc;
        tmpc.push_back(itr);
        int vcurrent = itr;
        for (int j = 0; j < vnbr; j++)
        {
            int d0b = d0_back[vcurrent];
            int d1b = d1_back[vcurrent];
            int bk;
            if (family == 2) // d0
            {
                bk = d0b;
            }

            if (family == 3) // d1
            {
                bk = d1b;
            }
            if (bk >= 0)
            {
                tmpc.push_back(bk);
                vcurrent = bk;
            }
            else{
                break;
            }
        }
        if (tmpc.size() > 0)
        {
            curves.push_back(tmpc);
            if (tmpc.size() == 1) // only have 1 start point but not end point
            {
                wrongs.push_back(V.row(tmpc[0]));
            }
        }
    }
    Debugtool = vec_list_to_matrix(wrongs);
    digs = curves;
}

// void assign_bnd_ver_bnms(const Eigen::VectorXi &Bnd, const std::vector<std::vector<int>> ids, std::vector<std::vector<Eigen::Vector3d>> &bnm)
// {
//     assert(ids.size() == bnm.size());
//     for (int i = 0; i < ids.size(); i++)
//     {
//         assert(ids[i].size() == bnm[i].size());
//         for (int j = 0; j < ids[i].size(); i++)
//         {
            
//         }
//     }
// }
void computeColRowInfo(const int vnbr, const int vinrow, std::vector<std::vector<double>> &rowinfo, std::vector<std::vector<double>> &colinfo)
{
    rowinfo.clear();
    colinfo.clear();
    int cnbr = vnbr / vinrow;
    for (int i = 0; i < cnbr; i++)
    {
        int start = i * vinrow;
        std::vector<double> cvec;
        for (int j = start; j < start + vinrow; j++)
        {
            cvec.push_back(j);
        }
        rowinfo.push_back(cvec);
    }
    for (int i = 0; i < vinrow; i++)
    {
        std::vector<double> rvec;
        for (int j = 0; j < cnbr; j++)
        {
            rvec.push_back(i + vinrow * j);
        }
        colinfo.push_back(rvec);
    }
}

void QuadOpt::write_polyline_info(){
    // the iterations might look weird, because in our data the info are in the form of 0, 1, 2, -1, -1, 3, 4, ...
    if (rowinfo.size() == 0 || colinfo.size() == 0)
    {
        if (vNbrInRow <= 0)
        {
            std::cout<<"please input the vinrow variable to tell me how many vertices in each row\n";
            return;
        }
        computeColRowInfo(V.rows(), vNbrInRow, rowinfo, colinfo);

    }
    PolyOpt plytool;
    Eigen::VectorXi Bnd;
    get_Bnd(Bnd);
    std::string fname = igl::file_dialog_save();
    if (fname.length() == 0)
    {
        std::cout << "\nLSC: save mesh failed, please type down the correct name" << std::endl;
        return;
    }
    std::vector<std::vector<Eigen::Vector3d>> lines_rows, lines_cols,
        bi_rows, bi_cols;
    std::vector<std::vector<int>> idrows, idcols, iddiags;
    int vnbr = V.rows();
    // extract the family 0
    int family = 0;
    int bnm_start = vnbr * 6;
    std::cout<<"Row: ";
    std::vector<Eigen::Vector3d> line, binormal;
    std::vector<int> ids;
    for (int i = 0; i < rowinfo.size(); i++)
    {
        line.clear();
        binormal.clear();
        ids.clear();

        for (int j = 0; j < rowinfo[i].size(); j++)
        {
            int vid = rowinfo[i][j];
            if (vid < 0 || Bnd[vid] == 1)
            {
                if (line.size() > 0)
                {
                    int real_size = line.size();
                    if (real_size > 2)
                    {
                        lines_rows.push_back(line);
                        bi_rows.push_back(binormal);
                        idrows.push_back(ids);
                    }
                    line.clear();
                    binormal.clear();
                    ids.clear();
                }
                continue;
            }
            if (Bnd[vid] != 1)
            {
                Eigen::Vector3d b_local;
                extract_binormals(family, bnm_start, vid, b_local);
                line.push_back(V.row(vid));
                binormal.push_back(b_local);
                ids.push_back(vid);
            }
        }
        if (line.size() > 0)
        {
            int real_size = line.size();
            if (real_size > 2)
            {
                lines_rows.push_back(line);
                bi_rows.push_back(binormal);
                idrows.push_back(ids);
            }
        }
    }
    family = 1;
    bnm_start =  vnbr * 12;
    std::cout<<"Col: ";
    for (int i = 0; i < colinfo.size(); i++)
    {
        line.clear();
        binormal.clear();
        ids.clear();

        for (int j = 0; j < colinfo[i].size(); j++)
        {
            int vid = colinfo[i][j];
            if (vid < 0 || Bnd[vid] == 1)
            {
                if (line.size() > 0)
                {
                    int real_size = line.size();
                    if (real_size > 2)
                    {
                        lines_cols.push_back(line);
                        bi_cols.push_back(binormal);
                        idcols.push_back(ids);
                    }

                    line.clear();
                    binormal.clear();
                    ids.clear();
                }
                continue;
            }
            if (Bnd[vid] != 1)
            {
                Eigen::Vector3d b_local;
                extract_binormals(family, bnm_start, vid, b_local);
                line.push_back(V.row(vid));
                binormal.push_back(b_local);
                ids.push_back(vid);
            }
        }
        if (line.size() > 0)
        {
            int real_size = line.size();
            if (real_size > 2)
            {
                lines_cols.push_back(line);
                bi_cols.push_back(binormal);
                idcols.push_back(ids);
            }
        }
    }

    
    std::vector<std::vector<Eigen::Vector3d>> bout0, bout1, bout2;
    if (OptType == 0 || OptType == 1 || OptType == 3) // AAG, GGA, PPG diagonals
    {
        std::vector<std::vector<int>> diagonals;
        int family = WhichDiagonal == 0 ? 2 : 3;
        extract_diagonals(family, diagonals);
        std::cout<<"Nbr of Diagonals: "<<diagonals.size()<<std::endl;;
        std::vector<std::vector<Eigen::Vector3d>> cv, bn;
        for (int i = 0; i < diagonals.size(); i++)
        {
            // for (int j = 0; j < diagonals[i].size(); j++)
            // {
            //     std::cout << diagonals[i][j] << ", ";
            // }
            // std::cout << "\n";
            std::vector<Eigen::Vector3d> tmpc, tmpb;
            std::vector<int> tmpi;
            for (int j = 0; j < diagonals[i].size(); j++)
            {
                int vid = diagonals[i][j];
                if (Bnd[vid] == 1)
                {
                    if (!tmpc.empty())
                    {
                        cv.push_back(tmpc);
                        bn.push_back(tmpb);
                        iddiags.push_back(tmpi);
                    }
                    tmpc.clear();
                    tmpb.clear();
                    tmpi.clear();
                    continue;
                }
                tmpc.push_back(V.row(diagonals[i][j]));
                Eigen::Vector3d localb;
                extract_binormals(family, 0, diagonals[i][j], localb);
                tmpb.push_back(localb);
                tmpi.push_back(diagonals[i][j]);
            }
            if (!tmpc.empty())
            {
                cv.push_back(tmpc);
                bn.push_back(tmpb);
                iddiags.push_back(tmpi);
            }
        }
        plytool.orient_binormals_of_plyline(bn, bout0);
        write_polyline_xyz(cv, fname + "_d");
        write_polyline_xyz(bout0, fname + "_d_b");
    }
    plytool.orient_binormals_of_plyline(bi_rows, bout1);
    write_polyline_xyz(lines_rows, fname+"_r");
    write_polyline_xyz(bout1, fname + "_r_b");
    plytool.orient_binormals_of_plyline(bi_cols, bout2);
    write_polyline_xyz(lines_cols, fname+"_c");
    write_polyline_xyz(bout2, fname + "_c_b");
    std::cout << "files get saved" << std::endl;
}

void QuadOpt::write_polyline_info_propagation(){
    // the iterations might look weird, because in our data the info are in the form of 0, 1, 2, -1, -1, 3, 4, ...
    int WhichType = varsize == V.rows() * 9 ? 0 : 1; // 0 is AAG, 1 is AGG
    int varsOff = WhichType == 0 ? 9 : 12;
    if (rowinfo.size() == 0 || colinfo.size() == 0)
    {
        std::cout << "Please use the other initialization method to give rowinfo and col info\n";
        return;
    }
    PolyOpt plytool;
    Eigen::VectorXi Bnd;
    get_Bnd(Bnd);
    std::string fname = igl::file_dialog_save();
    if (fname.length() == 0)
    {
        std::cout << "\nLSC: save mesh failed, please type down the correct name" << std::endl;
        return;
    }
    std::vector<std::vector<Eigen::Vector3d>> lines_rows, lines_cols,
        bi_rows, bi_cols;
    std::vector<std::vector<int>> idrows, idcols, iddiags;
    int vnbr = V.rows();
    // extract the family 0
    int family = 0;
    int bnm_start = vnbr * 6;
    std::cout<<"Row: ";
    std::vector<Eigen::Vector3d> line, binormal;
    std::vector<int> ids;
    for (int i = 0; i < rowinfo.size(); i++)
    {
        line.clear();
        binormal.clear();
        ids.clear();

        for (int j = 0; j < rowinfo[i].size(); j++)
        {
            int vid = rowinfo[i][j];
            if (vid < 0 || Bnd[vid] == 1)
            {
                if (line.size() > 0)
                {
                    int real_size = line.size();
                    if (real_size > 2)
                    {
                        lines_rows.push_back(line);
                        bi_rows.push_back(binormal);
                        idrows.push_back(ids);
                    }
                    line.clear();
                    binormal.clear();
                    ids.clear();
                }
                continue;
            }
            if (Bnd[vid] != 1)
            {
                Eigen::Vector3d b_local;
                extract_binormals_propagation(family, vid, b_local, varsOff);
                line.push_back(V.row(vid));
                binormal.push_back(b_local);
                ids.push_back(vid);
            }
        }
        if (line.size() > 0)
        {
            int real_size = line.size();
            if (real_size > 2)
            {
                lines_rows.push_back(line);
                bi_rows.push_back(binormal);
                idrows.push_back(ids);
            }
        }
    }
    family = 1;
    bnm_start =  vnbr * 12;
    std::cout<<"Col: ";
    for (int i = 0; i < colinfo.size(); i++)
    {
        line.clear();
        binormal.clear();
        ids.clear();

        for (int j = 0; j < colinfo[i].size(); j++)
        {
            int vid = colinfo[i][j];
            if (vid < 0 || Bnd[vid] == 1)
            {
                if (line.size() > 0)
                {
                    int real_size = line.size();
                    if (real_size > 2)
                    {
                        lines_cols.push_back(line);
                        bi_cols.push_back(binormal);
                        idcols.push_back(ids);
                    }

                    line.clear();
                    binormal.clear();
                    ids.clear();
                }
                continue;
            }
            if (Bnd[vid] != 1)
            {
                Eigen::Vector3d b_local;
                extract_binormals_propagation(family, vid, b_local, varsOff);
                line.push_back(V.row(vid));
                binormal.push_back(b_local);
                ids.push_back(vid);
            }
        }
        if (line.size() > 0)
        {
            int real_size = line.size();
            if (real_size > 2)
            {
                lines_cols.push_back(line);
                bi_cols.push_back(binormal);
                idcols.push_back(ids);
            }
        }
    }

    
    std::vector<std::vector<Eigen::Vector3d>> bout0, bout1, bout2;
    if (1) // use diagonals
    {
        std::vector<std::vector<int>> diagonals;
        int family = WhichDiagonal == 0 ? 2 : 3;
        extract_diagonals(family, diagonals);
        std::cout<<"Nbr of Diagonals: "<<diagonals.size()<<std::endl;;
        std::vector<std::vector<Eigen::Vector3d>> cv, bn;
        for (int i = 0; i < diagonals.size(); i++)
        {
            // for (int j = 0; j < diagonals[i].size(); j++)
            // {
            //     std::cout << diagonals[i][j] << ", ";
            // }
            // std::cout << "\n";
            std::vector<Eigen::Vector3d> tmpc, tmpb;
            std::vector<int> tmpi;
            for (int j = 0; j < diagonals[i].size(); j++)
            {
                int vid = diagonals[i][j];
                if (Bnd[vid] == 1)
                {
                    if (!tmpc.empty())
                    {
                        cv.push_back(tmpc);
                        bn.push_back(tmpb);
                        iddiags.push_back(tmpi);
                    }
                    tmpc.clear();
                    tmpb.clear();
                    tmpi.clear();
                    continue;
                }
                tmpc.push_back(V.row(diagonals[i][j]));
                Eigen::Vector3d localb;
                // extract_binormals(family, 0, diagonals[i][j], localb);
                extract_binormals_propagation(family, diagonals[i][j], localb, varsOff);
                tmpb.push_back(localb);
                tmpi.push_back(diagonals[i][j]);
            }
            if (!tmpc.empty())
            {
                cv.push_back(tmpc);
                bn.push_back(tmpb);
                iddiags.push_back(tmpi);
            }
        }
        plytool.orient_binormals_of_plyline(bn, bout0);
        write_polyline_xyz(cv, fname + "_d");
        write_polyline_xyz(bout0, fname + "_d_b");
    }
    plytool.orient_binormals_of_plyline(bi_rows, bout1);
    write_polyline_xyz(lines_rows, fname+"_r");
    write_polyline_xyz(bout1, fname + "_r_b");
    plytool.orient_binormals_of_plyline(bi_cols, bout2);
    write_polyline_xyz(lines_cols, fname+"_c");
    write_polyline_xyz(bout2, fname + "_c_b");
    std::cout << "files get saved" << std::endl;
}