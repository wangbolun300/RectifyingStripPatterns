#include <igl/readOFF.h>
#include <igl/slice.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/arap.h>
#include <igl/harmonic.h>
#include <igl/readOBJ.h>
#include <iostream>
#include <igl/unproject_ray.h>
#include <lsc/MeshProcessing.h>
#include <lsc/basic.h>
#include <lsc/tools.h>
#include <lsc/interaction.h>
// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <igl/parula.h>
#include <igl/isolines_map.h>

typedef OpenMesh::PolyMesh_ArrayKernelT<> CGMesh;

// the interface operators
namespace lscif
{

	const Eigen::RowVector3d sea_green(70. / 255., 252. / 255., 167. / 255.);
	const Eigen::RowVector3d hot_red(255. / 255., 12. / 255., 17. / 255.);
	const Eigen::RowVector3d pure_blk(0,0,0);
	MeshProcessing MP;

	// This V and F are only for visulization of the default mesh
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;

	bool keyPress_1 = false;
	bool keyPress_2 = false;
	bool keyPress_d = false;

	bool enable_functional_degrees = false;
	bool enable_boundary_angles = false;
	// bool show_

	// Optimization Parameters
	double weigth_closeness = 0.0;
	double refer_AveEL = 1;
	double weight_mass = 100.;	  // weight of the mass function to
	double weight_boundary = 100; // weight of the boundary condition
	double weight_laplacian = 10;
	double weight_pseudo_geodesic = 10;
	double weight_strip_width = 1;
	double weight_geodesic= 3;
    // double weight_shading = 0.01; // For shading: the light is othogonal to the tangent direction 

	double maximal_step_length = 0.5;
	bool enable_inner_vers_fixed = false;
	bool enable_extreme_cases;
	bool Given_Const_Direction = false;
	int which_levelset = 0;
	bool fix_angle_of_two_levelsets = false;
	double angle_between_two_levelsets = 60;
	double weight_fix_two_ls_angle = 0;

	// Tracing Parameters
	int which_seg_id = 0; // the face id of which we will assign value to
	int nbr_of_intervals = 3;
	int start_bnd_he = 0;// start boundary halfedge, for init the values by assign values on some edges
    int nbr_edges = 20;
	int id_debug_global = 5;
	double target_angle = 60;
	double target_angle_2 = 60;
	double target_min_angle = 90; // target angle for min function value;
	double target_max_angle = 90; // target angle for max function value;
	double start_angle = 60;
	double threadshold_angel_degree = 150; // the threadshold for detecting boundary corners
	int dbg_int;
	int dbg_int2 = 5;
	double dbg_dbl = 30;
	double dbg_dbl2 = 60;
	double vector_scaling = 1;
	int extracted_nbr = 5;
	bool debug_flag = false;
	std::vector<int> fixedVertices;
	double ptInx = 0;
	double ptIny = 0;
	double ptInz = 0;
	int TracingType = 0;
	// int QuadOptType = 0;


	int OpIter = 10;
	CGMesh mesh;
	CGMesh RefMesh;
	static bool enable_pg_energy_checkbox = false;
	static bool enable_strip_width_checkbox = false;
	static bool fixBoundary_checkbox = false;

	// mesh processing for level set
	int nbr_lines_first_ls = 30;
	int nbr_lines_second_ls = 30;
	Eigen::VectorXd readed_LS1;
	Eigen::VectorXd readed_LS2;
	Eigen::VectorXd readed_LS3;
	std::vector<CGMesh> readed_mesh1;
	CGMesh readed_mesh2;
	int Nbr_Iterations_Mesh_Opt=10;
	double weight_Mesh_mass = 100;
	double weight_Mesh_smoothness=0.001;
	double weight_Mesh_edgelength = 0.01;
    double weight_Mesh_pesudo_geodesic=30;
    double Mesh_opt_max_step_length=0.1;
	double weight_Mesh_approximation = 0.01;


	int update_ver_id;
	int ring_nbr;
	std::vector<int> update_verlist;

	static bool selectFEV[30] = {true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true};
	std::vector<std::string> meshFileName;
	std::vector<CGMesh> Meshes;
	lsTools tools;
	PolyOpt poly_tool;
	QuadOpt quad_tool;
	double weight_angle = 0;
	double InputPx = 0; // theta
	double InputPy = 180; // phi

	double InputPx1 = 0;
	double InputPy1 = 180;
	
	double InputThetaTol = 23.5; // give 10 degrees of tolerance
	double InputThetaTol1 = 23.5; // give 10 degrees of tolerance
	double InputPhiTol = 45; // give 10 degrees of tolerance
	double InputPhiTol1 = 45; // give 10 degrees of tolerance

	double Shading_Latitude = 34;
	bool enable_max_energy_check = false;
	bool shading_init = false;
	bool let_ray_through = false;
	bool let_ray_reflect = false;
	bool recompute_auxiliaries = false;
	double max_e_percentage = 10;
	double weight_binormal = 1;
	double weight_smt_binormal = 0;
	double weight_endpoint_ratio = 1;
	bool pick_single_ply = false;
	int pick_line_id = 0;
	std::vector<std::vector<Eigen::Vector3d>> Poly_readed;
	std::vector<std::vector<Eigen::Vector3d>> Bino_readed;
	std::vector<std::vector<Eigen::Vector3d>> Poly_opt;
	std::vector<std::vector<Eigen::Vector3d>> Bino_opt;

	std::vector<int> VertexType;

	// interaction stuff
	bool draw_strokes = false;
	bool left_button_down = false;
	std::vector<std::vector<double>> project_2dx;
	std::vector<std::vector<double>> project_2dy;
	std::vector<double> project_x_tmp;
	std::vector<double> project_y_tmp;
	std::vector<std::vector<int>> flist;
	std::vector<std::vector<Eigen::Vector3f>> bclist;
	std::vector<int> flist_tmp;
	std::vector<Eigen::Vector3f> bclist_tmp;
	std::vector<std::vector<Eigen::Vector3d>> project_pts;
	AutoRunArgs autorunner[2];

	// functional angles stuff.
	bool Disable_Changed_Angles = false;
    bool Use_Fitting_Angles = false; 
	bool Use_Opt_Only_BNMS = false;




	// add a new mesh into the mesh lists, and show the new mesh along with previous showed meshes
	void updateMeshViewer(igl::opengl::glfw::Viewer &viewer, CGMesh &mesh)
	{
		Eigen::MatrixXd V;
		Eigen::MatrixXi F;

		MP.mesh2Matrix(mesh, V, F);
		// Create new data slot and set to selected
		if (!(viewer.data().F.rows() == 0 && viewer.data().V.rows() == 0))
		{
			viewer.append_mesh();
		}
		//		viewer.data().clear();
		viewer.data().set_mesh(V, F);
		viewer.data().is_visible = true;

		// show polygon edges
		Eigen::MatrixXi E;
		MP.meshEdges(mesh, E);
		Eigen::MatrixXd C = Eigen::MatrixXd::Zero(E.rows(), 3);
		viewer.data().set_edges(V, E, C);

		viewer.data().show_lines = false;
		viewer.data().line_width = 2.0;
		viewer.core().align_camera_center(viewer.data().V, viewer.data().F);

		VertexType.resize(V.size(), 0);
	}

	bool pre_draw(igl::opengl::glfw::Viewer &viewer)
	{
		return false;
	}
	bool key_up(igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods)
	{
		switch (key)
		{
		case '1':
		{
			keyPress_1 = false;
			return true;
		}

		case '2':
		{
			keyPress_2 = false;
			return true;
		}
		case '3':
		{
			draw_strokes = false;
			return true;
		}
		case GLFW_KEY_X:
		{
			keyPress_d = false;
			return true;
		}
		}
		return false;
	}
	bool key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods)
	{
		switch (key)
		{

		case '1':
		{
			keyPress_1 = true;
			// std::cout << viewer.current_mouse_x << "___" << viewer.current_mouse_y << std::endl;
			return true;
		}

		case '2':
		{
			keyPress_2 = true;
			return true;
		}
		case '3':
		{
			draw_strokes = true;
			return true;
		}

		case GLFW_KEY_X:
		{
			keyPress_d = true;
			if (keyPress_d)
			{
				// keyPress_d = false;
				int id = viewer.selected_data_index;
				// if (id > 1)
				//	viewer.selected_data_index = id - 2;
				// std::cout << "The current Mesh ID is" << viewer.selected_data_index << std::endl;
				if (lscif::Meshes.size() == 1)
				{
					std::cout << "Do not delete the last mesh" << std::endl;
					// ImGui::End();
				}
				else
				{
					viewer.erase_mesh(id);
					if (id > -1)
					{
						lscif::meshFileName.erase(lscif::meshFileName.begin() + id);
						lscif::Meshes.erase(lscif::Meshes.begin() + id);
					}
					if (lscif::Meshes.size() > id)
					{
						viewer.selected_data_index = id;
					}
					else
					{

						viewer.selected_data_index = id - 1;
					}
				}
			}

			return true;
		}
		}

		return false;
	}

	// when holding "1" or "2", with mouse left button down, points can be choosen
	// TODO fix the bug when select one mesh but click on the other one, and when swich to
	// the other mesh, all the histories are recorded and reproduce points on the second mesh
	bool mouse_down(igl::opengl::glfw::Viewer &viewer, int button, int modifier)
	{
		if (keyPress_1)
		{
			int fid;
			Eigen::Vector3f bc;
			// Cast a ray in the view direction starting from the mouse position
			double x = viewer.current_mouse_x;
			double y = viewer.core().viewport(3) - viewer.current_mouse_y;
			if (igl::unproject_onto_mesh(
					Eigen::Vector2f(x, y),
					viewer.core().view,
					viewer.core().proj,
					viewer.core().viewport,
					viewer.data().V,
					viewer.data().F,
					fid,
					bc))
			{
				int max;
				bc.maxCoeff(&max);
				int vid = viewer.data().F(fid, max);

				if (VertexType[vid] != 1)
					VertexType[vid] = 1;
				else
					VertexType[vid] = 0;

				std::vector<int> fixedVid;
				std::vector<int> controlVid;
				for (int i = 0; i < viewer.data().V.size(); i++)
				{
					if (VertexType[i] == 1)
						fixedVid.push_back(i);
					if (VertexType[i] == 2)
						controlVid.push_back(i);
				}

				Eigen::VectorXi vids = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(fixedVid.data(), fixedVid.size());
				Eigen::VectorXi vids2 = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(controlVid.data(), controlVid.size());

				std::cout << "Selected point:\nposition related to the viewer: " << viewer.current_mouse_x << "___" << viewer.current_mouse_y << std::endl;
				std::cout << "vid, " << vid << ", fid, " << fid << ", mass_unfiform, " << lscif::tools.mass_uniform.coeffRef(vid, vid) << std::endl;
				if (lscif::tools.fvalues.rows() > 0)
				{
					std::cout << "level set value " << lscif::tools.fvalues[vid] << std::endl;
					
				}
				Eigen::VectorXd energy;
				lscif::tools.show_max_pg_energy(energy);
				if (energy.size() > 0)
				{
					std::cout << "local energy, " << energy[vid] << std::endl;
				}
				Eigen::MatrixXd energy_all;
				lscif::tools.show_max_pg_energy_all(energy_all);
				if (energy_all.rows() == viewer.data().V.rows())
				{
					std::cout << energy_all.row(vid) << std::endl;
				}
				else{
					std::cout<<"the size of energy_all, "<<energy_all.rows()<<std::endl;
				}

				std::cout << "point position: (" << viewer.data().V(vid, 0) << ", " << viewer.data().V(vid, 1) << ", " << viewer.data().V(vid, 2) << ")\n\n";
				// lscif::tools.print_info(vid);
				if (lscif::Given_Const_Direction)
				{
					Eigen::VectorXd diff;
					lscif::tools.shading_detect_parallel_patch(lscif::InputPx, lscif::InputPy, diff);

					std::cout << "light surface parallel " << diff[vid] << std::endl;
					double latitude_radian = lscif::Shading_Latitude * LSC_PI / 180.;
					double rho = (LSC_PI / 2 - latitude_radian);
					Eigen::Matrix3d rotation;
					rotation << 1, 0, 0,
						0, cos(rho), -sin(rho),
						0, sin(rho), cos(rho);
					Eigen::Vector3d direction_ground = Eigen::Vector3d(0, 0, -1); // the ground direction
					direction_ground = rotation * direction_ground;
					Eigen::MatrixXd ground = Eigen::MatrixXd::Ones(lscif::tools.V.rows(), lscif::tools.V.cols());
					ground.col(0) *= direction_ground[0];
					ground.col(1) *= direction_ground[1];
					ground.col(2) *= direction_ground[2];
					viewer.data().add_edges(lscif::tools.V - lscif::vector_scaling * ground, lscif::tools.V + lscif::vector_scaling * ground, lscif::pure_blk);
				}

				// viewer.data().add_points(igl::slice(viewer.data().V, vids, 1), hot_red);
				// viewer.data().set_points(igl::slice(viewer.data().V, vids, 1), lscif::hot_red);
				viewer.data().add_points(igl::slice(viewer.data().V, vids, 1), lscif::hot_red);
				Eigen::MatrixXd pts;
				lscif::tools.show_current_reference_points(pts);
				viewer.data().add_points(pts, lscif::sea_green);
				return true;
			}
		}
		if (keyPress_2)
		{
			int fid;
			Eigen::Vector3f bc;
			// Cast a ray in the view direction starting from the mouse position
			double x = viewer.current_mouse_x;
			double y = viewer.core().viewport(3) - viewer.current_mouse_y;
			if (igl::unproject_onto_mesh(
					Eigen::Vector2f(x, y),
					viewer.core().view,
					viewer.core().proj,
					viewer.core().viewport,
					viewer.data().V,
					viewer.data().F,
					fid,
					bc))
			{
				int max;
				bc.maxCoeff(&max);
				int vid = viewer.data().F(fid, max);

				if (VertexType[vid] != 2)
					VertexType[vid] = 2;
				else
					VertexType[vid] = 0;

				std::vector<int> fixedVid;
				std::vector<int> controlVid;
				for (int i = 0; i < viewer.data().V.size(); i++)
				{
					if (VertexType[i] == 1)
						fixedVid.push_back(i);
					if (VertexType[i] == 2)
						controlVid.push_back(i);
				}

				Eigen::VectorXi vids = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(fixedVid.data(), fixedVid.size());
				Eigen::VectorXi vids2 = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(controlVid.data(), controlVid.size());

				std::cout << "Selected point:\nposition related to the viewer: " << viewer.current_mouse_x << "___" << viewer.current_mouse_y << std::endl;
				std::cout << "vid, " << vid << ", fid, " << fid << std::endl;
				if (lscif::tools.fvalues.rows() > 0)
				{
					std::cout << "level set value " << lscif::tools.fvalues[vid] << std::endl;
				}
				std::cout << "point position: (" << viewer.data().V(vid, 0) << ", " << viewer.data().V(vid, 1) << ", " << viewer.data().V(vid, 2) << ")\n\n";
				viewer.data().set_points(igl::slice(viewer.data().V, vids, 1), lscif::hot_red);
				viewer.data().add_points(igl::slice(viewer.data().V, vids2, 1), lscif::sea_green);
				return true;
			}
		}
		
		left_button_down = true;
		return false;
	}
	bool mouse_up(igl::opengl::glfw::Viewer &viewer, int button, int modifier)
	{
		left_button_down = false;
		if (!project_x_tmp.empty()) // if just draw a curve
		{
			project_2dx.push_back(project_x_tmp);
			project_2dy.push_back(project_y_tmp);
			Eigen::MatrixXd Vers;
			draw_stroke_on_mesh(lscif::tools.lsmesh, viewer,
								lscif::tools.V, lscif::tools.F, project_x_tmp,
								project_y_tmp, Vers, flist_tmp, bclist_tmp);
			flist.push_back(flist_tmp);
			bclist.push_back(bclist_tmp);
			project_x_tmp.clear();
			project_y_tmp.clear();
			std::cout << "There are " << project_2dx.size() << " curves" << std::endl;
			viewer.data().add_points(Vers, lscif::hot_red);
		}

		return true;
		return false;
	}
	bool mouse_move(igl::opengl::glfw::Viewer & viewer, int mouse_x, int mouse_y){
		if (draw_strokes && left_button_down)
		{
			// Eigen::Vector3d posEnd = igl::unproject(
            //      Eigen::Vector3f(viewer.current_mouse_x,
            //                      viewer.core().viewport[3] -
            //                        static_cast<float>(viewer.current_mouse_y),
            //                      viewer.down_mouse_z),
            //      viewer.core().view,
            //      viewer.core().proj,
            //      viewer.core().viewport)
            //      .template cast<double>();
			// 	 std::cout<<"posend, "<<posEnd.transpose()<<std::endl;
				

				double x = viewer.current_mouse_x;
				double y = viewer.core().viewport(3) - viewer.current_mouse_y;
				lscif::project_x_tmp.push_back(x);
				lscif::project_y_tmp.push_back(y);
				
				// sleep(0.05);
				

				 return true;
		}
		return false;
	}

	void show_axis(Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, double ratio)
	{
		E0 = Eigen::MatrixXd::Zero(3, 3);
		E1 = ratio * Eigen::MatrixXd::Identity(3, 3);
	}

};

// ----------------------------------------------------------------------------

int main(int argc, char *argv[])
{

	std::map<int, Eigen::RowVector3d> colors;
	int last_selected = -1;
	std::string example_root_path(CHECKER_BOARD_DATA_DIR);
#ifdef _WIN32
	std::string spliter = "\\";
#else
	std::string spliter = "/";
#endif
	std::string modelname = "oc_5.000000_10.000000_50_30.obj";
	std::string fname = example_root_path + std::string(modelname);
	OpenMesh::IO::read_mesh(lscif::mesh, fname);
	lscif::tools.init(lscif::mesh);
	lscif::MP.mesh2Matrix(lscif::mesh, lscif::V, lscif::F);
	lscif::meshFileName.push_back(modelname);
	lscif::Meshes.push_back(lscif::mesh);

	// Init the viewer
	igl::opengl::glfw::Viewer viewer;
	igl::opengl::glfw::imgui::ImGuiPlugin plugin;
	viewer.plugins.push_back(&plugin);
	// Attach a menu plugin
	igl::opengl::glfw::imgui::ImGuiMenu menu;
	plugin.widgets.push_back(&menu);

	// Add content to the default menu window
	menu.callback_draw_viewer_menu = [&]()
	{
		// Draw parent menu content
		menu.draw_viewer_menu();
	};

	// Draw additional windows
	menu.callback_draw_custom_window = [&]()
	{
		// Define next window position + size
		ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(300, 800), ImGuiCond_FirstUseEver);
		ImGui::Begin(
			"Levelset Curves", nullptr,
			ImGuiWindowFlags_NoSavedSettings);

		// Add a button
		if (ImGui::Button("Import Mesh", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
		{

			std::string fname = igl::file_dialog_open();
			if (fname.length() == 0)
			{
				std::cout << "\nLSC: read mesh failed" << std::endl;
				ImGui::End();
				return;
			}

			OpenMesh::IO::read_mesh(lscif::mesh, fname);
			lscif::tools.init(lscif::mesh);
			size_t last_dot = fname.rfind('.');
			size_t last_slash = fname.rfind(spliter); // TODO on linux it should be '/'
			std::string fnameNew = fname.substr(last_slash + 1, (last_dot - last_slash - 1));
			lscif::meshFileName.push_back(fnameNew);
			lscif::Meshes.push_back(lscif::mesh);

			int nbe = 0;
			for (CGMesh::EdgeIter eit = lscif::mesh.edges_begin(); eit != lscif::mesh.edges_end(); ++eit)
			{
				if (lscif::mesh.is_boundary(eit))
					nbe++;
			}

			int nbhe = 0;
			for (CGMesh::HalfedgeIter heit = lscif::mesh.halfedges_begin(); heit != lscif::mesh.halfedges_end(); ++heit)
			{
				if (lscif::mesh.is_boundary(heit))
					nbhe++;
			}

			std::cout << "Mesh is: " << fname << std::endl;

			std::cout << "Vertices/Edges/Faces/HalfEdges/boundaryHalfEdge/boundaryEdge: " << lscif::mesh.n_vertices() << "/" << lscif::mesh.n_edges() << "/" << lscif::mesh.n_faces() << "/" << lscif::mesh.n_halfedges() << "/" << nbhe << "/" << nbe << std::endl;

			lscif::updateMeshViewer(viewer, lscif::mesh);
		}

		ImGui::SameLine();

		// Add a button
		if (ImGui::Button("Save Mesh", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
		{
			std::string fname = igl::file_dialog_save();

			if (fname.length() == 0)
			{
				std::cout << "\nLSC: save mesh failed" << std::endl;
				ImGui::End();
				return;
			}
			fname = fname + ".obj";

			int id = viewer.selected_data_index;
			OpenMesh::IO::write_mesh(lscif::Meshes[id], fname);
		}

		ImGui::SameLine();

		if (ImGui::Button("Import levelset", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
		{

			bool read = read_levelset(lscif::tools.fvalues);
			if (read)
			{
				std::cout << "\nlevel set readed" << std::endl;
			}
		}

		ImGui::SameLine();

		if (ImGui::Button("Save levelset", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
		{
			bool saved = save_levelset(lscif::tools.fvalues, lscif::tools.Binormals);
			if (saved)
			{
				std::cout << "\nlevel set saved" << std::endl;
			}
			
		}
		if (ImGui::Button("Read Binormals", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
		{
			std::cout<<"Select a Bi-normal file for visulization"<<std::endl;
			read_bi_normals(lscif::tools.Binormals);
			std::cout<<"bi-normals got readed"<<std::endl;
		}
		ImGui::SameLine();
		if (ImGui::Button("Debug", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
		{
			
			//  Eigen::Vector3d light =  get_light_rotated_back_from_earth_axis(lscif::Shading_Latitude, lscif::InputPx, lscif::InputPy);
			//  std::cout<<"The light direction rotated back is\n "<<light<<std::endl;
			// std::cout<<"checking the angle between the light and the target light"<<std::endl;
			// Eigen::Vector3d target_light = angle_ray_converter(lscif::InputPx, lscif::InputPy);
			// std::cout<<"target light, "<<target_light.transpose()<<std::endl;
			// auto light = lscif::tools.Lights;
			// for(int i = 0;i<light.rows();i++)
			// {
			// 	std::cout<<light.row(i).dot(target_light)<<std::endl;
			// }
			// Eigen::VectorXd diff;
			// std::cout<<"ploting the angle between the light and the surface"<<std::endl;
			// lscif::tools.shading_detect_parallel_patch(lscif::InputPx, lscif::InputPy, diff);
			// viewer.data().set_colors(diff);

		}
		ImGui::SameLine();
		if (ImGui::Button("UseUnitScaleMesh", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
		{
			std::cout<<"Mesh Unit Scale According to Edge Length"<<std::endl;
			int id = viewer.selected_data_index;
			CGMesh updatedMesh;
			CGMesh inputMesh = lscif::Meshes[id];
			lscif::MP.MeshUnitScale(inputMesh, updatedMesh);
			lscif::updateMeshViewer(viewer, updatedMesh);
			lscif::meshFileName.push_back("Uni_" + lscif::meshFileName[id]);
			lscif::Meshes.push_back(updatedMesh);
			viewer.selected_data_index = id;
			lscif::mesh = updatedMesh;
			lscif::tools.init(updatedMesh);


		}

		// Add new group
		if (ImGui::CollapsingHeader("Pseudo-Geodesic Angle", ImGuiTreeNodeFlags_DefaultOpen))
		{
			ImGui::PushItemWidth(50);
			ImGui::InputDouble("target angle", &lscif::target_angle, 0, 0, "%.4f");
			ImGui::SameLine();
			ImGui::InputDouble("second angle", &lscif::target_angle_2, 0, 0, "%.4f");
			ImGui::InputInt("Which Levelset", &lscif::which_levelset, 0, 0);
			ImGui::InputDouble("weight_geodesic", &lscif::weight_geodesic, 0, 0, "%.4f");
			
		}

		if (ImGui::CollapsingHeader("Energy Solving", ImGuiTreeNodeFlags_CollapsingHeader))
		{
			// Expose variable directly ...
			// ImGui::InputDouble("Closeness", &lscif::weigth_closeness, 0, 0, "%.4f");
			ImGui::InputInt("Iteration", &lscif::OpIter, 0, 0);
			ImGui::InputDouble("weight ls mass(big)", &lscif::weight_mass, 0, 0, "%.4f");
			ImGui::InputDouble("weight boundary (big)", &lscif::weight_boundary, 0, 0, "%.4f");
			ImGui::InputDouble("weight laplacian (small)", &lscif::weight_laplacian, 0, 0, "%.4f");
			ImGui::InputDouble("weight pseudo-geodesic", &lscif::weight_pseudo_geodesic, 0, 0, "%.4f");
			ImGui::InputDouble("weight Strip width(Need SW Enabled)", &lscif::weight_strip_width, 0, 0, "%.4f");
			ImGui::InputDouble("MaxStpLenght(Need PG Enabled)", &lscif::maximal_step_length, 0, 0, "%.4f");
			ImGui::InputDouble("WeightAngle", &lscif::weight_angle, 0, 0, "%.4f");
			ImGui::InputDouble("RatioEndPts", &lscif::weight_endpoint_ratio, 0, 0, "%.4f");

			ImGui::Checkbox("Enable functional Energy", &lscif::enable_functional_degrees);
			ImGui::SameLine();
			ImGui::Checkbox("Enable Extreme Cases", &lscif::enable_extreme_cases);

			if (lscif::enable_functional_degrees)
			{
				ImGui::InputDouble("angle (min func)", &lscif::target_min_angle, 0, 0, "%.4f");
				ImGui::InputDouble("angle (max func)", &lscif::target_max_angle, 0, 0, "%.4f");
			}

			// ImGui::InputDouble("Strip Width Ratio", &lscif::strip_width_ratio, 0, 0, "%.4f");

			ImGui::Checkbox("Enable PG Energy", &lscif::enable_pg_energy_checkbox);
			ImGui::SameLine();
			ImGui::Checkbox("Enable Strip Width", &lscif::enable_strip_width_checkbox);
			ImGui::Checkbox("Enable Smooth Boundary", &lscif::enable_inner_vers_fixed);
			ImGui::SameLine();
			ImGui::Checkbox("Enable Boundary Angles", &lscif::enable_boundary_angles);
			
		}

		// ImGui::Checkbox("Fix Boundary", &lscif::fixBoundary_checkbox);
	

	if (ImGui::CollapsingHeader("Tracing", ImGuiTreeNodeFlags_DefaultOpen))
	{
		ImGui::Combo("TracingType", &lscif::TracingType,
                 "SegTrace\0OneEdgeTrace\0SegDrctn\0\0"), 

		ImGui::InputInt("Select 1 boundary", &lscif::which_seg_id, 0, 0);
		ImGui::InputInt("every i segments", &lscif::nbr_of_intervals, 0, 0);
		ImGui::InputInt("start bnd edge ", &lscif::start_bnd_he, 0, 0);
		ImGui::InputInt("nbr bnd edges ", &lscif::nbr_edges, 0, 0);
		ImGui::InputDouble("start angle", &lscif::start_angle, 0, 0, "%.4f");
		ImGui::InputDouble("corner angle", &lscif::threadshold_angel_degree, 0, 0, "%.4f");
		if (ImGui::Button("Trace Curves", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				int id = viewer.selected_data_index;
				CGMesh updatedMesh;
				CGMesh inputMesh = lscif::Meshes[id];
				TracingPrepare Tracing_initializer;
				Tracing_initializer.every_n_edges = lscif::nbr_of_intervals;
				Tracing_initializer.start_angle = lscif::start_angle;
				Tracing_initializer.target_angle = lscif::target_angle;
				Tracing_initializer.threadshold_angel_degree = lscif::threadshold_angel_degree;
				Tracing_initializer.which_boundary_segment = lscif::which_seg_id;
				Tracing_initializer.start_bnd_he = lscif::start_bnd_he;
				Tracing_initializer.nbr_edges = lscif::nbr_edges;
				if(lscif::TracingType==0){
					lscif::tools.initialize_level_set_by_select_boundary_segment(Tracing_initializer, true);
				}
				if(lscif::TracingType==1){
					lscif::tools.initialize_level_set_by_tracing(Tracing_initializer);
				}	
				if(lscif::TracingType==2){
					lscif::tools.initialize_level_set_by_select_boundary_segment(Tracing_initializer, false);
				}
				
				
				
				std::cout << "finish tracing" << std::endl;
				lscif::updateMeshViewer(viewer, inputMesh);
				lscif::meshFileName.push_back("dbg_" + lscif::meshFileName[id]);
				lscif::Meshes.push_back(inputMesh);
				Eigen::MatrixXd E0, E1;
				Eigen::MatrixXd E2, E3, Ea0, Ea1;
				// lscif::tools.show_gradients(E0,E1, lscif::vector_scaling);
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);
				const Eigen::RowVector3d black(0, 0, 0);
				const Eigen::RowVector3d green(0.2, 0.8, 0.2);
				Eigen::MatrixXd RGB = Eigen::MatrixXd::Identity(3, 3);

				Eigen::MatrixXd pts;
				std::vector<Eigen::MatrixXd> E0list, E1list;
				std::cout << "before ploting the curve" << std::endl;
				lscif::tools.show_pseudo_geodesic_curve(E0list, E1list, pts);
				for (int i = 0; i < E0list.size(); i++) // plot the curves
				{
					E0 = E0list[i];
					E1 = E1list[i];
					std::cout << "edge sizes " << E0.rows() << ", " << E1.size() << std::endl;
					viewer.data().add_edges(E0, E1, red);
				}
				std::cout << "the number of curves " << E0list.size() << std::endl;

				if (1) // plot the vertices of the curves
				{
					viewer.data().add_points(pts, red);
				}
				viewer.selected_data_index = id;
				// std::cout<<"acos(1) "<<acos(1.)<<std::endl;
			}

			ImGui::SameLine();
			if (ImGui::Button("TracedBinormals", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				std::cout << "Showing the bi-normals of the traced curves. The green lines are the bi-normals and the blue lines are the surface normals" << std::endl;
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);
				const Eigen::RowVector3d black(0, 0, 0);
				const Eigen::RowVector3d green(0.2, 0.8, 0.2);
				Eigen::MatrixXd bE0, bE1, nE0, nE1;
				lscif::tools.show_traced_binormals(bE0,bE1,nE0,nE1, lscif::vector_scaling);

				viewer.data().add_edges(bE0, bE1, green);
				viewer.data().add_edges(nE0, nE1, blue);
			}
			ImGui::SameLine();
			if (ImGui::Button("SaveTracedCurves", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				int id = viewer.selected_data_index;
				CGMesh updatedMesh;
				CGMesh inputMesh = lscif::Meshes[id];
				TracingPrepare Tracing_initializer;
				Tracing_initializer.every_n_edges = lscif::nbr_of_intervals;
				Tracing_initializer.start_angle = lscif::start_angle;
				Tracing_initializer.target_angle = lscif::target_angle;
				Tracing_initializer.threadshold_angel_degree = lscif::threadshold_angel_degree;
				Tracing_initializer.which_boundary_segment = lscif::which_seg_id;
				Tracing_initializer.start_bnd_he = lscif::start_bnd_he;
				Tracing_initializer.nbr_edges = lscif::nbr_edges;
				if(lscif::TracingType==0){
					lscif::tools.initialize_level_set_by_select_boundary_segment(Tracing_initializer, true);
				}
				if(lscif::TracingType==1){
					lscif::tools.initialize_level_set_by_tracing(Tracing_initializer);
				}	
				if(lscif::TracingType==2){
					lscif::tools.initialize_level_set_by_select_boundary_segment(Tracing_initializer, false);
				}
				
				
				std::cout << "finish tracing" << std::endl;
				lscif::updateMeshViewer(viewer, inputMesh);
				lscif::meshFileName.push_back("dbg_" + lscif::meshFileName[id]);
				lscif::Meshes.push_back(inputMesh);
				Eigen::MatrixXd E0, E1;
				Eigen::MatrixXd E2, E3, Ea0, Ea1;
				// lscif::tools.show_gradients(E0,E1, lscif::vector_scaling);
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);
				const Eigen::RowVector3d black(0, 0, 0);
				const Eigen::RowVector3d green(0.2, 0.8, 0.2);
				Eigen::MatrixXd RGB = Eigen::MatrixXd::Identity(3, 3);

				Eigen::MatrixXd pts;
				std::vector<Eigen::MatrixXd> E0list, E1list;
				std::cout << "before ploting the curve" << std::endl;
				lscif::tools.show_pseudo_geodesic_curve(E0list, E1list, pts);
				for (int i = 0; i < E0list.size(); i++) // plot the curves
				{
					E0 = E0list[i];
					E1 = E1list[i];
					std::cout << "edge sizes " << E0.rows() << ", " << E1.size() << std::endl;
					viewer.data().add_edges(E0, E1, red);
				}
				std::cout << "the number of curves " << E0list.size() << std::endl;

				if (1) // plot the vertices of the curves
				{
					viewer.data().add_points(pts, red);
				}
				viewer.selected_data_index = id;
				std::cout<<"SAVING The Traced Vertices, please provide the prefix..."<<std::endl;
				std::string fname = igl::file_dialog_save();
				// std::ofstream file;
				// file.open(fname);
				lscif::tools.save_traced_curves(fname);
				
				// file.close();

				std::cout << "Traced Data SAVED" << std::endl;

				// std::cout<<"acos(1) "<<acos(1.)<<std::endl;
			}
			
			if (ImGui::Button("LoadTracedCurves", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				lscif::tools.load_one_traced_curve();
			}
			ImGui::SameLine();
			if (ImGui::Button("ClearTracedCurves", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				lscif::tools.clear_traced_curves();
			}
			ImGui::SameLine();
			if (ImGui::Button("SaveOBJCurve", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				int id = viewer.selected_data_index;
				CGMesh updatedMesh;
				CGMesh inputMesh = lscif::Meshes[id];
				TracingPrepare Tracing_initializer;
				Tracing_initializer.every_n_edges = lscif::nbr_of_intervals;
				Tracing_initializer.start_angle = lscif::start_angle;
				Tracing_initializer.target_angle = lscif::target_angle;
				Tracing_initializer.threadshold_angel_degree = lscif::threadshold_angel_degree;
				Tracing_initializer.which_boundary_segment = lscif::which_seg_id;
				Tracing_initializer.start_bnd_he = lscif::start_bnd_he;
				Tracing_initializer.nbr_edges = lscif::nbr_edges;
				if(lscif::TracingType==0){
					lscif::tools.initialize_level_set_by_select_boundary_segment(Tracing_initializer, true);
				}
				if(lscif::TracingType==1){
					lscif::tools.initialize_level_set_by_tracing(Tracing_initializer);
				}	
				if(lscif::TracingType==2){
					lscif::tools.initialize_level_set_by_select_boundary_segment(Tracing_initializer, false);
				}
				
				
				std::cout << "finish tracing" << std::endl;
				lscif::updateMeshViewer(viewer, inputMesh);
				lscif::meshFileName.push_back("dbg_" + lscif::meshFileName[id]);
				lscif::Meshes.push_back(inputMesh);
				Eigen::MatrixXd E0, E1;
				Eigen::MatrixXd E2, E3, Ea0, Ea1;
				// lscif::tools.show_gradients(E0,E1, lscif::vector_scaling);
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);
				const Eigen::RowVector3d black(0, 0, 0);
				const Eigen::RowVector3d green(0.2, 0.8, 0.2);
				Eigen::MatrixXd RGB = Eigen::MatrixXd::Identity(3, 3);

				Eigen::MatrixXd pts;
				std::vector<Eigen::MatrixXd> E0list, E1list;
				std::cout << "before ploting the curve" << std::endl;
				lscif::tools.show_pseudo_geodesic_curve(E0list, E1list, pts);
				for (int i = 0; i < E0list.size(); i++) // plot the curves
				{
					E0 = E0list[i];
					E1 = E1list[i];
					std::cout << "edge sizes " << E0.rows() << ", " << E1.size() << std::endl;
					viewer.data().add_edges(E0, E1, red);
				}
				std::cout << "the number of curves " << E0list.size() << std::endl;

				if (1) // plot the vertices of the curves
				{
					viewer.data().add_points(pts, red);
				}
				viewer.selected_data_index = id;
				std::cout<<"SAVING The Traced Vertices, Binormals and Normals into obj format, please provide prefix"<<std::endl;
				std::string fname = igl::file_dialog_save();
				std::ofstream file;
				file.open(fname+".obj");
				for (int i = 0; i < pts.rows(); i++)
				{
					file << "v " << pts(i, 0) << " " << pts(i, 1) << " " << pts(i, 2) << std::endl;
				}

				file.close();

				Eigen::MatrixXd bE0, bE1, nE0, nE1;
				lscif::tools.show_traced_binormals(bE0, bE1, nE0, nE1, lscif::vector_scaling);

				
				file.open(fname+"_start.obj");
				assert(bE0.rows() == pts.rows() - 2 && "Please Trace only one curve so we can save the correct data");
				for (int i = 0; i < bE0.rows(); i++)
				{
					file << "v " << pts(i + 1, 0) << " " << pts(i + 1, 1) << " " << pts(i + 1, 2) << std::endl;
				}

				file.close();

				file.open(fname+"_binormal_end.obj");
				for (int i = 0; i < bE1.rows(); i++)
				{
					file << "v " << bE1(i, 0) << " " << bE1(i, 1) << " " << bE1(i, 2) << std::endl;
				}

				file.close();


				file.open(fname+"_normal_end.obj");
				for (int i = 0; i < nE1.rows(); i++)
				{
					file << "v " << nE1(i, 0) << " " << nE1(i, 1) << " " << nE1(i, 2) << std::endl;
				}

				file.close();

				

				std::cout << "Traced Data SAVED" << std::endl;

				// std::cout<<"acos(1) "<<acos(1.)<<std::endl;
			}
			ImGui::SameLine();
			if (ImGui::Button("SaveCurveParam", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{// this is to save the parametrized tracing results.

				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);

				const Eigen::RowVector3d black(0, 0, 0);
				const Eigen::RowVector3d green(0.2, 0.8, 0.2);

				Eigen::MatrixXd pts;
				std::cout << "Load one curve, convert it into 2d, then save it as an obj file" << std::endl;

				lscif::tools.show_traced_curve_params(pts);
				if(pts.rows()==0){
					std::cout<<"Please follow the instructions before using this button"<<std::endl;
					ImGui::End();
					return;
				}
				viewer.data().add_points(pts, red);
				std::cout<<"Saving the curve points in parametrization domain..."<<std::endl;
				std::string fname = igl::file_dialog_save();
				std::ofstream file;
				file.open(fname);
				for (int i = 0; i < pts.rows(); i++)
				{
					file << "v " << pts(i, 0) << " " << pts(i, 1) << " " << pts(i, 2) << std::endl;
				}

				file.close();

				
			}
		// ImGui::SameLine();
		// ImGui::Checkbox("Fix Boundary", &lscif::fixBoundary_checkbox);
		}
		ImGui::Separator();

		if (ImGui::CollapsingHeader("Debug input", ImGuiTreeNodeFlags_CollapsingHeader))
		{
			ImGui::Checkbox("debugFlag", &lscif::debug_flag);
			ImGui::InputInt("int", &lscif::dbg_int, 0, 0);
			ImGui::InputInt("int2", &lscif::dbg_int2, 0, 0);
			ImGui::InputDouble("double", &lscif::dbg_dbl, 0, 0, "%.4f");
			ImGui::InputDouble("double2", &lscif::dbg_dbl2, 0, 0, "%.4f");
			ImGui::InputDouble("vector_scaling", &lscif::vector_scaling, 0, 0, "%.4f");
			ImGui::InputInt("extracted_nbr", &lscif::extracted_nbr, 0, 0);
			ImGui::PushItemWidth(50);
			ImGui::InputDouble("x", &lscif::ptInx, 0, 0, "%.4f");
			ImGui::SameLine();
			ImGui::InputDouble("y", &lscif::ptIny, 0, 0, "%.4f");
			ImGui::SameLine();
			ImGui::InputDouble("z", &lscif::ptInz, 0, 0, "%.4f");
			ImGui::SameLine();
			if (ImGui::Button("drawPt", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				Eigen::MatrixXd pts(1,3);
				pts.row(0)<<lscif::ptInx,lscif::ptIny,lscif::ptInz;
				viewer.data().add_points(pts, lscif::hot_red);
			}
			if (ImGui::Button("ReadVers&Write", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				read_ver_write_into_xyz_files();
			}
			ImGui::SameLine();
			if (ImGui::Button("TMP", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				// read_pts_csv_and_write_xyz_files();
				// recover_polyline_endpts();
				// using intersections of the rectifying planes to get the developable.
				// construct_single_developable_strips_by_intersect_rectifying(0); 
				// write_unfold_single_strip();    
				// draw_catenaries_on_cylinder();
				Eigen::MatrixXd Vcout; 
				Eigen::MatrixXd Vlout;
				match_the_two_catenaries(lscif::tools.lsmesh, lscif::tools.Boundary_Edges, lscif::tools.V,
										 lscif::tools.F, lscif::tools.fvalues, Vcout, Vlout);
				viewer.data().add_points(Vcout, lscif::hot_red);
				viewer.data().add_points(Vlout, lscif::sea_green);
			}
		}
		if (ImGui::CollapsingHeader("LS Processing", ImGuiTreeNodeFlags_DefaultOpen))
		{
			


			if (ImGui::Button("LvSet Opt", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{

				
				EnergyPrepare einit;
				einit.weight_gravity = lscif::weight_mass;
				einit.weight_lap = lscif::weight_laplacian;
				einit.weight_bnd = lscif::weight_boundary;
				einit.weight_pg = lscif::weight_pseudo_geodesic;
				einit.weight_strip_width = lscif::weight_strip_width;
				einit.solve_pseudo_geodesic = lscif::enable_pg_energy_checkbox;
				einit.target_angle = lscif::target_angle;
				einit.max_step_length = lscif::maximal_step_length;
				einit.solve_strip_width_on_traced = lscif::enable_strip_width_checkbox;
				einit.enable_inner_vers_fixed = lscif::enable_inner_vers_fixed;
				einit.enable_functional_angles = lscif::enable_functional_degrees;
				einit.enable_extreme_cases = lscif::enable_extreme_cases;
				einit.target_min_angle = lscif::target_min_angle;
				einit.target_max_angle = lscif::target_max_angle;
				einit.start_angle = lscif::start_angle;
				einit.enable_boundary_angles = lscif::enable_boundary_angles;
				einit.Given_Const_Direction = lscif::Given_Const_Direction;
				// lscif::tools.weight_shading = lscif::weight_shading;
				lscif::tools.prepare_level_set_solving(einit);
				lscif::tools.Second_Ray_vers = lscif::update_verlist;
				lscif::tools.Second_Ray_nbr_rings = lscif::ring_nbr;
				lscif::tools.Reference_theta = lscif::InputPx;
				lscif::tools.Reference_phi = lscif::InputPy;
				lscif::tools.Reference_theta1 = lscif::InputPx1;
				lscif::tools.Reference_phi1 = lscif::InputPy1;
				lscif::tools.Theta_tol = lscif::InputThetaTol;
				lscif::tools.Phi_tol = lscif::InputPhiTol;
				lscif::tools.Theta_tol1 = lscif::InputThetaTol1;
				lscif::tools.Phi_tol1 = lscif::InputPhiTol1;
				lscif::tools.weight_geodesic=lscif::weight_geodesic;
				lscif::tools.enable_max_energy_check = lscif::enable_max_energy_check;
				lscif::tools.max_energy_percentage = lscif::max_e_percentage;
				lscif::tools.enable_shading_init = lscif::shading_init;
				lscif::tools.weight_binormal = lscif::weight_binormal;
				lscif::tools.weight_smt_binormal = lscif::weight_smt_binormal;
				lscif::tools.enable_let_ray_through = lscif::let_ray_through;
				lscif::tools.ShadingLatitude = lscif::Shading_Latitude;
				lscif::tools.enable_reflection = lscif::let_ray_reflect;
				lscif::tools.recompute_auxiliaries = lscif::recompute_auxiliaries;
				

				for (int i = 0; i < lscif::OpIter; i++)
				{
					lscif::tools.Run_Level_Set_Opt();
					if (lscif::tools.step_length < 1e-16 && i != 0)
					{ // step length actually is the value for the last step
						std::cout << "optimization converges " << std::endl;
						break;
					}
				}
				std::cout << "waiting for instruction..." << std::endl;
				// lscif::MP.MeshUnitScale(inputMesh, updatedMesh);
				int id = viewer.selected_data_index;
				CGMesh inputMesh = lscif::tools.lsmesh;
				lscif::updateMeshViewer(viewer, inputMesh);
				lscif::meshFileName.push_back("lso_" + lscif::meshFileName[id]);
				lscif::Meshes.push_back(inputMesh);

				Eigen::VectorXd level_set_values;
				lscif::tools.show_level_set(level_set_values);
				if(level_set_values.size()==0){
					ImGui::End();
					return;
				}
				Eigen::MatrixXd CM;
				// std::cout<<"before compute colormap"<<std::endl;
				igl::parula(Eigen::VectorXd::LinSpaced(21, 0, 1).eval(), false, CM);
				igl::isolines_map(Eigen::MatrixXd(CM), CM);
				// std::cout<<"before set colormap"<<std::endl;
				viewer.data().set_colormap(CM);
				// std::cout<<"before set color"<<std::endl;
				viewer.data().set_data(level_set_values);
				// Eigen::MatrixXd E0, E1;
				// // lscif::tools.show_gradients(E0,E1, lscif::vector_scaling);
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);

				viewer.selected_data_index = id;
			}
			ImGui::SameLine();
			if (ImGui::Button("Mesh Opt", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				if(lscif::tools.fvalues.size()==0){
					std::cout<<"Please load or initialize the level-set before opt the mesh"<<std::endl;
					ImGui::End();
					return;
				}
				int id = viewer.selected_data_index;
				CGMesh updatedMesh;
				CGMesh inputMesh = lscif::Meshes[id];
				MeshEnergyPrepare initializer;
				initializer.Mesh_opt_max_step_length=lscif::Mesh_opt_max_step_length;
				initializer.weight_Mesh_pesudo_geodesic=lscif::weight_Mesh_pesudo_geodesic;
				initializer.weight_Mesh_smoothness=lscif::weight_Mesh_smoothness;
				initializer.weight_mass=lscif::weight_mass;
				initializer.target_angle=lscif::target_angle;
				initializer.weight_Mesh_edgelength = lscif::weight_Mesh_edgelength;
				initializer.enable_extreme_cases=lscif::enable_extreme_cases;
				initializer.Given_Const_Direction = lscif::Given_Const_Direction;
				lscif::tools.weight_Mesh_approximation = lscif::weight_Mesh_approximation;
				lscif::tools.weight_Mesh_mass = lscif::weight_Mesh_mass;
				// lscif::tools.weight_shading = lscif::weight_shading;
				lscif::tools.Reference_theta = lscif::InputPx;
				lscif::tools.Reference_phi = lscif::InputPy;
				// lscif::tools.Reference_theta2 = lscif::InputPx1;
				// lscif::tools.Reference_phi2 = lscif::InputPy1;
				lscif::tools.Theta_tol = lscif::InputThetaTol;
				lscif::tools.Phi_tol = lscif::InputPhiTol;
				// lscif::tools.Theta_tol2 = lscif::InputPhiTol;
				// lscif::tools.Phi_tol2 = lscif::InputPhiTol1;
				lscif::tools.prepare_mesh_optimization_solving(initializer);
				for(int i=0;i<lscif::Nbr_Iterations_Mesh_Opt;i++){
					lscif::tools.Run_Mesh_Opt();
				}
				updatedMesh=lscif::tools.lsmesh;
				lscif::updateMeshViewer(viewer, updatedMesh);
				lscif::meshFileName.push_back("mso_" + lscif::meshFileName[id]);
				lscif::Meshes.push_back(updatedMesh);
				
				viewer.selected_data_index = id;
				std::cout<<"waiting for instructions"<<std::endl;
				
			}
			ImGui::SameLine();
			if (ImGui::Button("Orthogonal LS", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				int id = viewer.selected_data_index;
				CGMesh inputMesh = lscif::tools.lsmesh;
				EnergyPrepare einit;
				einit.weight_gravity = lscif::weight_mass;
				einit.weight_lap = lscif::weight_laplacian;
				einit.weight_bnd = lscif::weight_boundary;
				einit.weight_pg = lscif::weight_pseudo_geodesic;
				einit.weight_strip_width = lscif::weight_strip_width;
				einit.solve_pseudo_geodesic = lscif::enable_pg_energy_checkbox;
				einit.target_angle = lscif::target_angle;
				einit.max_step_length = lscif::maximal_step_length;
				einit.solve_strip_width_on_traced = lscif::enable_strip_width_checkbox;
				einit.enable_inner_vers_fixed = lscif::enable_inner_vers_fixed;
				einit.enable_functional_angles = lscif::enable_functional_degrees;
				einit.enable_extreme_cases = lscif::enable_extreme_cases;
				einit.target_min_angle = lscif::target_min_angle;
				einit.target_max_angle = lscif::target_max_angle;
				einit.start_angle = lscif::start_angle;
				einit.enable_boundary_angles = lscif::enable_boundary_angles;
				lscif::tools.fix_angle_of_two_levelsets = lscif::fix_angle_of_two_levelsets;
				lscif::tools.angle_between_two_levelsets = lscif::angle_between_two_levelsets;
				lscif::tools.weight_fix_two_ls_angle = lscif::weight_fix_two_ls_angle;
				// lscif::tools.weight_shading = lscif::weight_shading;
				lscif::tools.prepare_level_set_solving(einit);
				if(lscif::readed_LS1.size()==0){
					std::cout<<"Please load LevelSet 1 as a reference"<<std::endl;
					ImGui::End();
					return;
				}
				for (int i = 0; i < lscif::OpIter; i++)
				{
					lscif::tools.Run_Othogonal_Levelset(lscif::readed_LS1);
					if (lscif::tools.step_length < 1e-16 && i != 0)
					{ // step length actually is the value for the last step
						std::cout << "optimization converges " << std::endl;
						break;
					}
				}
				std::cout << "waiting for instruction..." << std::endl;
				// lscif::MP.MeshUnitScale(inputMesh, updatedMesh);
				lscif::updateMeshViewer(viewer, inputMesh);
				lscif::meshFileName.push_back("lso_" + lscif::meshFileName[id]);
				lscif::Meshes.push_back(inputMesh);

				Eigen::VectorXd level_set_values;
				lscif::tools.show_level_set(level_set_values);
				if(level_set_values.size()==0){
					ImGui::End();
					return;
				}
				Eigen::MatrixXd CM;
				// std::cout<<"before compute colormap"<<std::endl;
				igl::parula(Eigen::VectorXd::LinSpaced(21, 0, 1).eval(), false, CM);
				igl::isolines_map(Eigen::MatrixXd(CM), CM);
				// std::cout<<"before set colormap"<<std::endl;
				viewer.data().set_colormap(CM);
				// std::cout<<"before set color"<<std::endl;
				viewer.data().set_data(level_set_values);
				// Eigen::MatrixXd E0, E1;
				// // lscif::tools.show_gradients(E0,E1, lscif::vector_scaling);
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);

				viewer.selected_data_index = id;
				
			}
			ImGui::SameLine();
			if (ImGui::Button("AngleVariable", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{

				
				EnergyPrepare einit;
				einit.weight_gravity = lscif::weight_mass;
				einit.weight_lap = lscif::weight_laplacian;
				einit.weight_bnd = lscif::weight_boundary;
				einit.weight_pg = lscif::weight_pseudo_geodesic;
				einit.weight_strip_width = lscif::weight_strip_width;
				einit.solve_pseudo_geodesic = lscif::enable_pg_energy_checkbox;
				einit.target_angle = lscif::target_angle;
				// std::cout<<"Target angle, "<<lscif::target_angle<<std::endl; exit(0);
				einit.max_step_length = lscif::maximal_step_length;
				einit.solve_strip_width_on_traced = lscif::enable_strip_width_checkbox;
				einit.enable_inner_vers_fixed = lscif::enable_inner_vers_fixed;
				einit.enable_functional_angles = lscif::enable_functional_degrees;
				einit.enable_extreme_cases = lscif::enable_extreme_cases;
				einit.target_min_angle = lscif::target_min_angle;
				einit.target_max_angle = lscif::target_max_angle;
				einit.start_angle = lscif::start_angle;
				einit.enable_boundary_angles = lscif::enable_boundary_angles;
				einit.Given_Const_Direction = lscif::Given_Const_Direction;
				// lscif::tools.weight_shading = lscif::weight_shading;
				lscif::tools.prepare_level_set_solving(einit);
				lscif::tools.Second_Ray_vers = lscif::update_verlist;
				lscif::tools.Second_Ray_nbr_rings = lscif::ring_nbr;
				lscif::tools.Reference_theta = lscif::InputPx;
				lscif::tools.Reference_phi = lscif::InputPy;
				lscif::tools.Reference_theta1 = lscif::InputPx1;
				lscif::tools.Reference_phi1 = lscif::InputPy1;
				lscif::tools.Theta_tol = lscif::InputThetaTol;
				lscif::tools.Phi_tol = lscif::InputPhiTol;
				lscif::tools.Theta_tol1 = lscif::InputThetaTol1;
				lscif::tools.Phi_tol1 = lscif::InputPhiTol1;
				lscif::tools.weight_geodesic=lscif::weight_geodesic;
				lscif::tools.enable_max_energy_check = lscif::enable_max_energy_check;
				lscif::tools.max_energy_percentage = lscif::max_e_percentage;
				lscif::tools.enable_shading_init = lscif::shading_init;
				lscif::tools.weight_binormal = lscif::weight_binormal;
				lscif::tools.weight_smt_binormal = lscif::weight_smt_binormal;
				lscif::tools.enable_let_ray_through = lscif::let_ray_through;
				lscif::tools.ShadingLatitude = lscif::Shading_Latitude;
				lscif::tools.enable_reflection = lscif::let_ray_reflect;
				lscif::tools.recompute_auxiliaries = lscif::recompute_auxiliaries;
				

				for (int i = 0; i < lscif::OpIter; i++)
				{
					lscif::tools.Run_Level_Set_Opt_Angle_Variable();
					if (lscif::tools.step_length < 1e-16 && i != 0)
					{ // step length actually is the value for the last step
						std::cout << "optimization converges " << std::endl;
						break;
					}
				}
				std::cout << "waiting for instruction..." << std::endl;
				// lscif::MP.MeshUnitScale(inputMesh, updatedMesh);
				int id = viewer.selected_data_index;
				CGMesh inputMesh = lscif::tools.lsmesh;
				lscif::updateMeshViewer(viewer, inputMesh);
				lscif::meshFileName.push_back("lso_" + lscif::meshFileName[id]);
				lscif::Meshes.push_back(inputMesh);

				Eigen::VectorXd level_set_values;
				lscif::tools.show_level_set(level_set_values);
				if(level_set_values.size()==0){
					ImGui::End();
					return;
				}
				Eigen::MatrixXd CM;
				// std::cout<<"before compute colormap"<<std::endl;
				igl::parula(Eigen::VectorXd::LinSpaced(21, 0, 1).eval(), false, CM);
				igl::isolines_map(Eigen::MatrixXd(CM), CM);
				// std::cout<<"before set colormap"<<std::endl;
				viewer.data().set_colormap(CM);
				// std::cout<<"before set color"<<std::endl;
				viewer.data().set_data(level_set_values);
				// Eigen::MatrixXd E0, E1;
				// // lscif::tools.show_gradients(E0,E1, lscif::vector_scaling);
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);

				viewer.selected_data_index = id;
			}
			if (ImGui::Button("AAG LvSet", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				if(lscif::readed_LS1.size()==0|| lscif::readed_LS2.size()==0){
					std::cout<<"Please Read the correct levelset 1 and levelset 2"<<std::endl;
					ImGui::End();
					return;
				}
				int id = viewer.selected_data_index;
				CGMesh inputMesh = lscif::Meshes[id];
				EnergyPrepare einit;
				einit.weight_gravity = lscif::weight_mass;
				einit.weight_lap = lscif::weight_laplacian;
				einit.weight_bnd = lscif::weight_boundary;
				einit.weight_pg = lscif::weight_pseudo_geodesic;
				einit.weight_strip_width = lscif::weight_strip_width;
				einit.solve_pseudo_geodesic = lscif::enable_pg_energy_checkbox;
				einit.target_angle = lscif::target_angle;
				einit.max_step_length = lscif::maximal_step_length;
				einit.solve_strip_width_on_traced = lscif::enable_strip_width_checkbox;
				einit.enable_inner_vers_fixed = lscif::enable_inner_vers_fixed;
				einit.enable_functional_angles = lscif::enable_functional_degrees;
				einit.enable_extreme_cases = lscif::enable_extreme_cases;
				einit.target_min_angle = lscif::target_min_angle;
				einit.target_max_angle = lscif::target_max_angle;
				einit.start_angle = lscif::start_angle;
				einit.enable_boundary_angles = lscif::enable_boundary_angles;
				
				lscif::tools.prepare_level_set_solving(einit);
				lscif::tools.weight_geodesic=lscif::weight_geodesic;
				
				for (int i = 0; i < lscif::OpIter; i++)
				{
					lscif::tools.Run_AAG(lscif::readed_LS1, lscif::readed_LS2, lscif::readed_LS3);
					if (lscif::tools.step_length < 1e-16 && i != 0)
					{ // step length actually is the value for the last step
						std::cout << "optimization converges " << std::endl;
						break;
					}
				}
				std::cout << "waiting for instruction..." << std::endl;
				// lscif::MP.MeshUnitScale(inputMesh, updatedMesh);
				lscif::updateMeshViewer(viewer, inputMesh);
				lscif::meshFileName.push_back("aag_" + lscif::meshFileName[id]);
				lscif::Meshes.push_back(inputMesh);

				viewer.selected_data_index = id;
				
			}
			ImGui::SameLine();

			if (ImGui::Button("AAG MshOpt", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				if(lscif::readed_LS1.size()==0|| lscif::readed_LS2.size()==0){
					std::cout<<"Please Read the correct levelset 1 and levelset 2"<<std::endl;
					ImGui::End();
					return;
				}
				int id = viewer.selected_data_index;
				CGMesh updatedMesh;
				CGMesh inputMesh = lscif::Meshes[id];
				MeshEnergyPrepare initializer;
				initializer.Mesh_opt_max_step_length=lscif::Mesh_opt_max_step_length;
				initializer.weight_Mesh_pesudo_geodesic=lscif::weight_Mesh_pesudo_geodesic;
				initializer.weight_Mesh_smoothness=lscif::weight_Mesh_smoothness;
				initializer.weight_mass=lscif::weight_mass;
				initializer.target_angle = lscif::target_angle;
				initializer.weight_Mesh_edgelength = lscif::weight_Mesh_edgelength;
				initializer.enable_extreme_cases = lscif::enable_extreme_cases;
				lscif::tools.weight_geodesic=lscif::weight_geodesic;
				lscif::tools.weight_Mesh_mass = lscif::weight_Mesh_mass;
				lscif::tools.weight_Mesh_approximation = lscif::weight_Mesh_approximation;
				lscif::tools.prepare_mesh_optimization_solving(initializer);
				for (int i = 0; i < lscif::Nbr_Iterations_Mesh_Opt; i++)
				{
					lscif::tools.Run_AAG_Mesh_Opt(lscif::readed_LS1, lscif::readed_LS2, lscif::readed_LS3);
				}
				updatedMesh = lscif::tools.lsmesh;
				lscif::updateMeshViewer(viewer, updatedMesh);
				lscif::meshFileName.push_back("aagm_" + lscif::meshFileName[id]);
				lscif::Meshes.push_back(updatedMesh);

				viewer.selected_data_index = id;
				std::cout << "waiting for instructions" << std::endl;
			}
			ImGui::SameLine();
			if (ImGui::Button("AGG LvSet", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				if(lscif::readed_LS1.size()==0|| lscif::readed_LS2.size()==0){
					std::cout<<"Please Read the correct levelset 1 and levelset 2"<<std::endl;
					ImGui::End();
					return;
				}
				int id = viewer.selected_data_index;
				CGMesh inputMesh = lscif::Meshes[id];
				EnergyPrepare einit;
				einit.weight_gravity = lscif::weight_mass;
				einit.weight_lap = lscif::weight_laplacian;
				einit.weight_bnd = lscif::weight_boundary;
				einit.weight_pg = lscif::weight_pseudo_geodesic;
				einit.weight_strip_width = lscif::weight_strip_width;
				einit.solve_pseudo_geodesic = lscif::enable_pg_energy_checkbox;
				einit.target_angle = lscif::target_angle;
				einit.max_step_length = lscif::maximal_step_length;
				einit.solve_strip_width_on_traced = lscif::enable_strip_width_checkbox;
				einit.enable_inner_vers_fixed = lscif::enable_inner_vers_fixed;
				einit.enable_functional_angles = lscif::enable_functional_degrees;
				einit.enable_extreme_cases = lscif::enable_extreme_cases;
				einit.target_min_angle = lscif::target_min_angle;
				einit.target_max_angle = lscif::target_max_angle;
				einit.start_angle = lscif::start_angle;
				einit.enable_boundary_angles = lscif::enable_boundary_angles;
				lscif::tools.fix_angle_of_two_levelsets = lscif::fix_angle_of_two_levelsets;
				lscif::tools.angle_between_two_levelsets = lscif::angle_between_two_levelsets;
				lscif::tools.weight_fix_two_ls_angle = lscif::weight_fix_two_ls_angle;
				
				lscif::tools.prepare_level_set_solving(einit);
				lscif::tools.weight_geodesic=lscif::weight_geodesic;
				
				for (int i = 0; i < lscif::OpIter; i++)
				{
					lscif::tools.Run_AGG(lscif::readed_LS1, lscif::readed_LS2, lscif::readed_LS3);
					if (lscif::tools.step_length < 1e-16 && i != 0)
					{ // step length actually is the value for the last step
						std::cout << "optimization converges " << std::endl;
						break;
					}
				}
				std::cout << "waiting for instruction..." << std::endl;
				// lscif::MP.MeshUnitScale(inputMesh, updatedMesh);
				lscif::updateMeshViewer(viewer, inputMesh);
				lscif::meshFileName.push_back("agg_" + lscif::meshFileName[id]);
				lscif::Meshes.push_back(inputMesh);

				viewer.selected_data_index = id;
			}
			ImGui::SameLine();
			if (ImGui::Button("AGG MshOpt", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				if(lscif::readed_LS1.size()==0|| lscif::readed_LS2.size()==0){
					std::cout<<"Please Read the correct levelset 1 and levelset 2"<<std::endl;
					ImGui::End();
					return;
				}
				int id = viewer.selected_data_index;
				CGMesh updatedMesh;
				CGMesh inputMesh = lscif::Meshes[id];
				MeshEnergyPrepare initializer;
				initializer.Mesh_opt_max_step_length=lscif::Mesh_opt_max_step_length;
				initializer.weight_Mesh_pesudo_geodesic=lscif::weight_Mesh_pesudo_geodesic;
				initializer.weight_Mesh_smoothness=lscif::weight_Mesh_smoothness;
				initializer.weight_mass=lscif::weight_mass;
				initializer.target_angle = lscif::target_angle;
				initializer.weight_Mesh_edgelength = lscif::weight_Mesh_edgelength;
				initializer.enable_extreme_cases = lscif::enable_extreme_cases;
				lscif::tools.weight_geodesic=lscif::weight_geodesic;
				lscif::tools.weight_Mesh_mass = lscif::weight_Mesh_mass;
				lscif::tools.weight_Mesh_approximation = lscif::weight_Mesh_approximation;
				lscif::tools.prepare_mesh_optimization_solving(initializer);
				for (int i = 0; i < lscif::Nbr_Iterations_Mesh_Opt; i++)
				{
					lscif::tools.Run_AGG_Mesh_Opt(lscif::readed_LS1, lscif::readed_LS2, lscif::readed_LS3);
				}
				updatedMesh = lscif::tools.lsmesh;
				lscif::updateMeshViewer(viewer, updatedMesh);
				lscif::meshFileName.push_back("aggm_" + lscif::meshFileName[id]);
				lscif::Meshes.push_back(updatedMesh);

				viewer.selected_data_index = id;
				std::cout << "waiting for instructions" << std::endl;
			}
			
			if (ImGui::Button("PPG LvSet", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				if(lscif::readed_LS1.size()==0|| lscif::readed_LS2.size()==0){
					std::cout<<"Please Read the correct levelset 1 and levelset 2"<<std::endl;
					ImGui::End();
					return;
				}
				int id = viewer.selected_data_index;
				CGMesh inputMesh = lscif::Meshes[id];
				EnergyPrepare einit;
				einit.weight_gravity = lscif::weight_mass;
				einit.weight_lap = lscif::weight_laplacian;
				einit.weight_bnd = lscif::weight_boundary;
				einit.weight_pg = lscif::weight_pseudo_geodesic;
				einit.weight_strip_width = lscif::weight_strip_width;
				einit.solve_pseudo_geodesic = lscif::enable_pg_energy_checkbox;
				einit.target_angle = lscif::target_angle;
				einit.max_step_length = lscif::maximal_step_length;
				einit.solve_strip_width_on_traced = lscif::enable_strip_width_checkbox;
				einit.enable_inner_vers_fixed = lscif::enable_inner_vers_fixed;
				einit.enable_functional_angles = lscif::enable_functional_degrees;
				einit.enable_extreme_cases = lscif::enable_extreme_cases;
				einit.target_min_angle = lscif::target_min_angle;
				einit.target_max_angle = lscif::target_max_angle;
				einit.start_angle = lscif::start_angle;
				einit.enable_boundary_angles = lscif::enable_boundary_angles;
				lscif::tools.fix_angle_of_two_levelsets = lscif::fix_angle_of_two_levelsets;
				lscif::tools.angle_between_two_levelsets = lscif::angle_between_two_levelsets;
				lscif::tools.weight_fix_two_ls_angle = lscif::weight_fix_two_ls_angle;
				lscif::tools.weight_geodesic=lscif::weight_geodesic;
				lscif::tools.pseudo_geodesic_target_angle_degree = lscif::target_angle;
				lscif::tools.pseudo_geodesic_target_angle_degree_2 = lscif::target_angle_2;

				lscif::tools.prepare_level_set_solving(einit);
				
				
				for (int i = 0; i < lscif::OpIter; i++)
				{
					lscif::tools.Run_PPG(lscif::readed_LS1, lscif::readed_LS2, lscif::readed_LS3);
					if (lscif::tools.step_length < 1e-16 && i != 0)
					{ // step length actually is the value for the last step
						std::cout << "optimization converges " << std::endl;
						break;
					}
				}
				std::cout << "waiting for instruction..." << std::endl;
				// lscif::MP.MeshUnitScale(inputMesh, updatedMesh);
				lscif::updateMeshViewer(viewer, inputMesh);
				lscif::meshFileName.push_back("agg_" + lscif::meshFileName[id]);
				lscif::Meshes.push_back(inputMesh);

				viewer.selected_data_index = id;
			}
			ImGui::SameLine();
			if (ImGui::Button("draw AAG", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				if(lscif::readed_LS1.size()==0){
					std::cout<<"Please Read the correct levelset 1 and levelset 2"<<std::endl;
					ImGui::End();
					return;
				}
				Eigen::MatrixXd CM;
				igl::parula(Eigen::VectorXd::LinSpaced(21, 0, 1).eval(), false, CM);
				igl::isolines_map(Eigen::MatrixXd(CM), CM);
				viewer.data().set_colormap(CM);
				if (lscif::which_levelset == 0)
				{
					viewer.data().set_data(lscif::readed_LS1);
				}
				if (lscif::which_levelset == 1)
				{
					viewer.data().set_data(lscif::readed_LS2);
				}
				if (lscif::which_levelset == 2)
				{
					viewer.data().set_data(lscif::readed_LS3);
				}
			}
			ImGui::SameLine();
			if (ImGui::Button("Save AAG", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				bool saved = save_levelset(lscif::readed_LS1);
				if (saved)
				{
					std::cout << "\nlevel 1  get saved" << std::endl;
				}
				saved = save_levelset(lscif::readed_LS2);
				if (saved)
				{
					std::cout << "\nlevel 2  get saved" << std::endl;
				}
				saved = save_levelset(lscif::readed_LS3);
				if (saved)
				{
					std::cout << "\nlevel 3  get saved" << std::endl;
				}
			}

			if (ImGui::Button("draw extracted", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				if(lscif::tools.fvalues.size()==0){
					std::cout<<"Please init the levelset or load a correct levelset"<<std::endl;
					ImGui::End();
					return;
				}

				int id = viewer.selected_data_index;
				CGMesh inputMesh = lscif::Meshes[id];

				// lscif::MP.MeshUnitScale(inputMesh, updatedMesh);
				lscif::updateMeshViewer(viewer, inputMesh);
				lscif::meshFileName.push_back("extracted_" + lscif::meshFileName[id]);
				lscif::Meshes.push_back(inputMesh);

				// viewer.data().set_colors(level_set_values);
				Eigen::MatrixXd E0, E1;
				// // lscif::tools.show_gradients(E0,E1, lscif::vector_scaling);
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);

				const Eigen::RowVector3d black(0, 0, 0);
				const Eigen::RowVector3d green(0.2, 0.8, 0.2);
				std::vector<Eigen::MatrixXd> E0list, E1list;
				lscif::tools.extract_levelset_curves(lscif::extracted_nbr, E0list, E1list);
				for (int i = 0; i < E0list.size(); i++) // plot the curves
				{
					E0 = E0list[i];
					E1 = E1list[i];
					std::cout << "traced seg sizes " << E0.rows() << ", " << E1.size() << std::endl;
					viewer.data().add_edges(E0, E1, red);
				}
				viewer.selected_data_index = id;
			}
			ImGui::SameLine();
			if (ImGui::Button("draw traced", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);

				const Eigen::RowVector3d black(0, 0, 0);
				const Eigen::RowVector3d green(0.2, 0.8, 0.2);

				std::vector<Eigen::MatrixXd> E0list, E1list;
				Eigen::MatrixXd pts, E0, E1;
				std::cout << "before ploting the traced curves" << std::endl;

				lscif::tools.show_pseudo_geodesic_curve(E0list, E1list, pts);
				if (E0list.size() == 0)
				{
					std::cout << "please trace the curves first" << std::endl;
				}

				for (int i = 0; i < E0list.size(); i++) // plot the curves
				{
					E0 = E0list[i];
					E1 = E1list[i];
					std::cout << "edge sizes " << E0.rows() << ", " << E1.size() << std::endl;
					viewer.data().add_edges(E0, E1, red);
				}
				std::cout << "the number of curves " << E0list.size() << std::endl;

				if (1) // plot the vertices of the curves
				{
					viewer.data().add_points(pts, red);
				}
			}
			ImGui::SameLine();
			if (ImGui::Button("draw levelset", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				Eigen::VectorXd level_set_values;
				lscif::tools.show_level_set(level_set_values);
				if(level_set_values.size()==0){
					std::cout<<"Please init the levelset or load a correct levelset"<<std::endl;
					ImGui::End();
					return;
				}
				Eigen::MatrixXd CM;
				igl::parula(Eigen::VectorXd::LinSpaced(21, 0, 1).eval(), false, CM);
				igl::isolines_map(Eigen::MatrixXd(CM), CM);
				viewer.data().set_colormap(CM);
				viewer.data().set_data(level_set_values);
			}
			ImGui::SameLine();
			if (ImGui::Button("draw bi-normals", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);
				const Eigen::RowVector3d yellow(241./255, 196./255, 15./255);

				const Eigen::RowVector3d black(0, 0, 0);
				const Eigen::RowVector3d green(0.2, 0.8, 0.2);

				Eigen::MatrixXd pts, E0, E1;
				std::cout << "before ploting the bi-normals" << std::endl;
				Eigen::MatrixXd binormals;
				lscif::tools.show_binormals(lscif::tools.fvalues, E0, E1, binormals, lscif::vector_scaling);
				viewer.data().add_edges(E0, E1, green);
				viewer.data().add_edges(lscif::tools.V, lscif::tools.V + lscif::vector_scaling * lscif::tools.norm_v, black);
				if (lscif::enable_extreme_cases && lscif::Given_Const_Direction)
				{
					std::cout << "Ploting the light directions" << std::endl;
					Eigen::MatrixXd light(lscif::tools.V.rows(), 3);
					if (lscif::tools.Lights.rows() == 0)
					{
						Eigen::Vector3d direction;
						direction = angle_ray_converter(lscif::InputPx, lscif::InputPy);
						
						for (int i = 0; i < light.rows(); i++)
						{
							light.row(i) = direction;
						}
						
					}
					else{
						light = lscif::tools.Lights;

					}
					Eigen::MatrixXd E2, E3;
					E2 = lscif::tools.V + lscif::vector_scaling * light;
					E3 = lscif::tools.V - lscif::vector_scaling * light;
					viewer.data().add_edges(E2, E3, yellow);
				}
			}

			if (ImGui::Button("draw error", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				Eigen::VectorXd error;
				lscif::tools.show_max_pg_energy(error);
				if(error.size()==0){
					std::cout<<"Please run Levelset Opt to get the pseudo-geodesic error"<<std::endl;
					ImGui::End();
					return;
				}
				
				viewer.data().set_data(error);
			}
			ImGui::SameLine();
			if (ImGui::Button("draw RefPts", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				std::cout<<"drawing the reference points in yellow"<<std::endl;
				Eigen::MatrixXd pts;
				lscif::tools.show_current_reference_points(pts);
				const Eigen::RowVector3d yellow(241./255, 196./255, 15./255);
				viewer.data().add_points(pts, yellow);
				
			}
			ImGui::SameLine();
			if (ImGui::Button("draw Gauss", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				std::cout<<"drawing the sign of Gaussian curvature"<<std::endl;
				
				Eigen::MatrixXd pts, color;
				lscif::tools.show_posi_negt_curvature(color);
				viewer.data().add_points(lscif::tools.V, color);
				
			}
			ImGui::SameLine();
			if (ImGui::Button("draw Paramtrztn", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				int id = viewer.selected_data_index;
				CGMesh inputMesh = lscif::tools.write_parameterization_mesh();

				lscif::updateMeshViewer(viewer, inputMesh);
				lscif::meshFileName.push_back("lso_" + lscif::meshFileName[id]);
				lscif::Meshes.push_back(inputMesh);
				viewer.selected_data_index = id;

				Eigen::VectorXd level_set_values;
				lscif::tools.show_level_set(level_set_values);
				if(level_set_values.size()==0){
					std::cout<<"Please init the levelset or load a correct levelset"<<std::endl;
					ImGui::End();
					return;
				}
				Eigen::MatrixXd CM;
				igl::parula(Eigen::VectorXd::LinSpaced(21, 0, 1).eval(), false, CM);
				igl::isolines_map(Eigen::MatrixXd(CM), CM);
				viewer.data().set_colormap(CM);
				viewer.data().set_data(level_set_values);
				
			}
			if (ImGui::Button("SaveUnitScaleMesh", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				std::cout << "This button is to unit scale the target mesh and save the result mesh\n ";
				Eigen::MatrixXd vout;
				mesh_unit_scale(lscif::tools.V, vout);
				std::string fname = igl::file_dialog_save();
				if (fname.length() == 0)
				{
					ImGui::End();
				}
				else
				{
					igl::writeOBJ(fname, vout, lscif::tools.F);
					std::cout << "mesh saved" << std::endl;
				}
			}
			ImGui::SameLine();
			if (ImGui::Button("Focus Selected Mesh", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				int id = viewer.selected_data_index;
				viewer.core().align_camera_center(viewer.data_list[id].V);
			}
			ImGui::SameLine();

			if (ImGui::Button("Delete Selected Mesh", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				int id = viewer.selected_data_index;
				// if (id > 1)
				//	viewer.selected_data_index = id - 2;
				// std::cout << "The current Mesh ID is" << viewer.selected_data_index << std::endl;
				if (lscif::Meshes.size() == 1)
				{
					std::cout << "Do not delete the last mesh" << std::endl;
					ImGui::End();
					return;
				}
				else
				{
					viewer.erase_mesh(id);
					if (id > -1)
					{
						lscif::meshFileName.erase(lscif::meshFileName.begin() + id);
						lscif::Meshes.erase(lscif::Meshes.begin() + id);
					}
					if (lscif::Meshes.size() > id)
					{
						viewer.selected_data_index = id;
					}
					else
					{

						viewer.selected_data_index = id - 1;
					}
				}
			}
			ImGui::SameLine();

			if (ImGui::Button("Show boundary Ver", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				int id = viewer.selected_data_index;

				CGMesh inputMesh = lscif::Meshes[id];

				std::vector<int> Bvids;
				for (CGMesh::VertexIter v_it = inputMesh.vertices_begin(); v_it != (inputMesh.vertices_end()); ++v_it)
				{

					if (inputMesh.is_boundary(v_it.handle()))
						Bvids.push_back(v_it.handle().idx());
				}
				Eigen::VectorXi vids = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(Bvids.data(), Bvids.size());

				viewer.data().set_points(igl::slice(viewer.data().V, vids, 1), lscif::hot_red);
			}

			if (ImGui::Button("ReadPlylines", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f))){
				int id = viewer.selected_data_index;
				read_plylines_and_binormals(lscif::Poly_readed, lscif::Bino_readed);
				if(lscif::Poly_readed.empty()){
					std::cout << "ERROR, Please load the polylines first" << std::endl;
					ImGui::End();
					return;
				}
				CGMesh updatemesh = polyline_to_strip_mesh(lscif::Poly_readed, lscif::Bino_readed, lscif::vector_scaling);
				lscif::updateMeshViewer(viewer, updatemesh);
				lscif::meshFileName.push_back("ply_" + lscif::meshFileName[id]);
				lscif::Meshes.push_back(updatemesh);
				lscif::poly_tool.init(lscif::Poly_readed, lscif::Bino_readed, lscif::nbr_lines_first_ls);
				std::cout<<"Sample nbr: "<<lscif::nbr_lines_first_ls<<std::endl;
				viewer.selected_data_index = id;
			}
			ImGui::SameLine();
			if (ImGui::Button("OptiPolylines", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f))){
				if(lscif::Poly_readed.empty()){
					std::cout << "ERROR, Please load the polylines first" << std::endl;
					ImGui::End();
					return;
				}
				lscif::poly_tool.weight_smooth= lscif::weight_laplacian;
				lscif::poly_tool.weight_mass = lscif::weight_mass;
				lscif::poly_tool.weight_binormal = lscif::weight_pseudo_geodesic;
				lscif::poly_tool.max_step = lscif::maximal_step_length;
				lscif::poly_tool.strip_scale = lscif::vector_scaling;
				lscif::poly_tool.binormal_ratio = lscif::weight_geodesic;
				lscif::poly_tool.weight_angle = lscif::weight_angle;
				lscif::poly_tool.target_angle = lscif::target_angle;
				lscif::poly_tool.ratio_endpts = lscif::weight_endpoint_ratio;
				lscif::poly_tool.pick_single_line = lscif::pick_single_ply;
				lscif::poly_tool.pick_line_id = lscif::pick_line_id;


				for (int i = 0; i < lscif::OpIter; i++)
				{
					if (lscif::weight_angle > 0)
					{
						lscif::poly_tool.get_normal_vector_from_reference(lscif::tools.V, lscif::tools.F, lscif::tools.norm_v);
					}
					lscif::poly_tool.opt();
				}
				CGMesh updatemesh = lscif::poly_tool.RecMesh;
				int id = viewer.selected_data_index;
				lscif::updateMeshViewer(viewer, updatemesh);
				lscif::meshFileName.push_back("ply_" + lscif::meshFileName[id]);
				lscif::Meshes.push_back(updatemesh);
				viewer.selected_data_index = id;
				std::cout<<"waiting for instructions"<<std::endl;

			}
		}
		ImGui::SameLine();
		if (ImGui::Button("SavePolylines", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
		{
			if (lscif::poly_tool.ply_extracted.empty())
			{
				std::cout << "ERROR, Please opt the polyline first" << std::endl;
				ImGui::End();
				return;
			}
			lscif::poly_tool.save_polyline_and_binormals_as_files(false); // true means save the rotated polylines
		}
		ImGui::SameLine();
		if (ImGui::Button("Orient&SavePlys", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))

		{
			std::cout<<"this function is to orient binormals if there are inverted ones"<<std::endl;
			int id = viewer.selected_data_index;
			read_plylines_and_binormals(lscif::Poly_readed, lscif::Bino_readed);
			std::vector<std::vector<Eigen::Vector3d>> biout;
			lscif::poly_tool.orient_binormals_of_plyline(lscif::Bino_readed, biout);
			CGMesh updatemesh = polyline_to_strip_mesh(lscif::Poly_readed, biout, lscif::vector_scaling);
			lscif::updateMeshViewer(viewer, updatemesh);
			lscif::meshFileName.push_back("ply_" + lscif::meshFileName[id]);
			lscif::Meshes.push_back(updatemesh);
			viewer.selected_data_index = id;
			lscif::poly_tool.save_polyline_and_binormals_as_files(lscif::Poly_readed, biout); // true means save the rotated polylines
		}
		if (ImGui::Button("RotateBackPly", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
		{
			if (lscif::poly_tool.ply_extracted.empty())
			{
				std::cout << "ERROR, Please opt the polyline first" << std::endl;
				ImGui::End();
				return;
			}
			std::cout<<"Rotating the mesh to the coordinate system computed by latitude"<<std::endl;
			lscif::poly_tool.rotate_back_the_model_to_horizontal_coordinates(lscif::Shading_Latitude);
			lscif::poly_tool.save_polyline_and_binormals_as_files(true); // true means save the rotated polylines
		}
		ImGui::SameLine();
		if (ImGui::Button("ForcePlyBinm", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
		{
			if (lscif::poly_tool.ply_extracted.empty())
			{
				std::cout << "ERROR, Please opt the polyline first" << std::endl;
				ImGui::End();
				return;
			}
			lscif::poly_tool.force_smoothing_binormals();
			std::cout<<"The binormals are forced to be smoothed"<<std::endl;

			CGMesh updatemesh = lscif::poly_tool.RecMesh;
			int id = viewer.selected_data_index;
			lscif::updateMeshViewer(viewer, updatemesh);
			lscif::meshFileName.push_back("Smt_" + lscif::meshFileName[id]);
			lscif::Meshes.push_back(updatemesh);
			viewer.selected_data_index = id;
			std::cout << "waiting for instructions" << std::endl;
		}
		ImGui::SameLine();
		if (ImGui::Button("DevelopableStrips", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
		{
			if (lscif::Poly_readed.empty())
			{
				std::cout << "ERROR, Please read the polyline first" << std::endl;
				ImGui::End();
				return;
			}
			// lscif::Poly_readed, lscif::Bino_readed
			std::vector<std::vector<Eigen::Vector3d>> dev_crease;
			for (int i = 0; i < lscif::Poly_readed.size(); i++)
			{
				std::vector<Eigen::Vector3d> creases;
				convert_polyline_to_developable_strips_reflection_method(lscif::Poly_readed[i], lscif::Bino_readed[i], creases);
				dev_crease.push_back(creases);
			}
			int id = viewer.selected_data_index;
			CGMesh updatemesh = polyline_to_strip_mesh(lscif::Poly_readed, dev_crease, lscif::vector_scaling);
			lscif::updateMeshViewer(viewer, updatemesh);
			lscif::meshFileName.push_back("Dev_" + lscif::meshFileName[id]);
			lscif::Meshes.push_back(updatemesh);
			viewer.selected_data_index = id;
			lscif::poly_tool.save_polyline_and_binormals_as_files(lscif::Poly_readed, dev_crease);
			std::cout << "waiting for instructions" << std::endl;
		}
		ImGui::SameLine();
		if (ImGui::Button("EvltDvlp", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
		{
			std::cout<<"evaluate developable of the loaded strips ..."<<std::endl;
			// get_polyline_rectifying_plane_on_inner_vers(const std::vector<Eigen::Vector3d> &ply_in, const std::vector<Eigen::Vector3d> &bnm_in,
            //                                      std::vector<Eigen::Vector3d> &vertices, std::vector<Eigen::Vector3d> &tangents, std::vector<Eigen::Vector3d> &binormals);
			evaluate_and_print_strip_developability(lscif::Poly_readed, lscif::Bino_readed);

		}
		if (ImGui::Button("ReadCreases", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
		{
			int id = viewer.selected_data_index;
			read_plylines_and_binormals(lscif::Poly_readed, lscif::Bino_readed);
			if (lscif::Poly_readed.empty())
			{
				std::cout << "ERROR, Please load the polylines first" << std::endl;
				ImGui::End();
				return;
			}
			CGMesh updatemesh = polyline_to_strip_mesh(lscif::Poly_readed, lscif::Bino_readed, lscif::vector_scaling);
			lscif::updateMeshViewer(viewer, updatemesh);
			lscif::meshFileName.push_back("ply_" + lscif::meshFileName[id]);
			lscif::Meshes.push_back(updatemesh);

			std::vector<std::vector<Eigen::Vector3d>> vertices;
			std::vector<std::vector<Eigen::Vector3d>> tangents;
			std::vector<std::vector<Eigen::Vector3d>> binormals;
			get_polyline_rectifying_planes(lscif::Poly_readed, lscif::Bino_readed, vertices, tangents, binormals);
			lscif::poly_tool.init_crease_opt(vertices, tangents, binormals);

			viewer.selected_data_index = id;
		}
		ImGui::SameLine();
		if (ImGui::Button("OptCrease", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
		{
			if (lscif::Poly_readed.empty())
			{
				std::cout << "ERROR, Please load the polylines first" << std::endl;
				ImGui::End();
				return;
			}
			lscif::poly_tool.weight_smooth = lscif::weight_laplacian;
			lscif::poly_tool.weight_mass = lscif::weight_mass;
			lscif::poly_tool.weight_binormal = lscif::weight_pseudo_geodesic;
			lscif::poly_tool.max_step = lscif::maximal_step_length;
			lscif::poly_tool.strip_scale = lscif::vector_scaling;
			lscif::poly_tool.binormal_ratio = lscif::weight_geodesic;
			lscif::poly_tool.weight_angle = lscif::weight_angle;
			lscif::poly_tool.target_angle = lscif::target_angle;
			lscif::poly_tool.ratio_endpts = lscif::weight_endpoint_ratio;
			lscif::poly_tool.pick_single_line = lscif::pick_single_ply;
			lscif::poly_tool.pick_line_id = lscif::pick_line_id;

			for (int i = 0; i < lscif::OpIter; i++)
			{
				// if (lscif::weight_angle > 0)
				// {
				// 	lscif::poly_tool.get_normal_vector_from_reference(lscif::tools.V, lscif::tools.F, lscif::tools.norm_v);
				// }
				lscif::poly_tool.opt_planarity();
			}
			CGMesh updatemesh = lscif::poly_tool.RecMesh;
			int id = viewer.selected_data_index;
			lscif::updateMeshViewer(viewer, updatemesh);
			lscif::meshFileName.push_back("ply_" + lscif::meshFileName[id]);
			lscif::Meshes.push_back(updatemesh);
			viewer.selected_data_index = id;
			std::cout << "waiting for instructions" << std::endl;
		}
		ImGui::SetNextItemOpen(true);
		if (ImGui::TreeNode("Mesh Management"))
		{

			// HelpMarker("This is a more typical looking tree with selectable nodes.\nClick to select, CTRL+Click to toggle, click on arrows or double-click to open.");

			static int selection_mask = (1 << 2); // Dumb representation of what may be user-side selection state. You may carry selection state inside or outside your objects in whatever format you see fit.
			int node_clicked = -1;				  // Temporary storage of what node we have clicked to process selection at the end of the loop. May be a pointer to your own node type, etc.
			// ImGui::PushStyleVar(ImGuiStyleVar_IndentSpacing, ImGui::GetFontSize() * 3); // Increase spacing to differentiate leaves from expanded contents.
			int id = viewer.data_list.size();

			for (int i = 0; i < id; i++)
			{
				//// Disable the default open on single-click behavior and pass in Selected flag according to our selection state.
				ImGuiTreeNodeFlags node_flags = ImGuiTreeNodeFlags_OpenOnArrow | ImGuiTreeNodeFlags_OpenOnDoubleClick;
				if (selection_mask & (1 << i))
					node_flags |= ImGuiTreeNodeFlags_Selected;
				// static bool check = true;
				// ImGui::Checkbox("checkbox", &check);
				// ImGui::SameLine();
				// ImGui::Checkbox("F", &selectFEV[3 * i]);
				// ImGui::SameLine();
				// ImGui::Checkbox("E", &selectFEV[3 * i + 1]);
				// ImGui::SameLine();
				// ImGui::Checkbox("V", &selectFEV[3 * i + 2]);
				// ImGui::SameLine();
				ImGui::PushID(i);
				ImGui::Checkbox("Visible", &lscif::selectFEV[i]);
				ImGui::PopID();
				ImGui::SameLine();
				if (lscif::selectFEV[i])
					viewer.data_list[i].is_visible = true;
				else
					viewer.data_list[i].is_visible = false;

				node_flags |= ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen; // ImGuiTreeNodeFlags_Bullet
				// node_flags |=  ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_DefaultOpen; // ImGuiTreeNodeFlags_Bullet
				// ImGui::TreeNodeEx((void*)(intptr_t)i, node_flags, "Object %d", i);
				ImGui::TreeNodeEx((void *)(intptr_t)i, node_flags, lscif::meshFileName[i].c_str());

				if (ImGui::IsItemClicked())
					node_clicked = i;
			}
			if (node_clicked != -1)
			{
				// Update selection state. Process outside of tree loop to avoid visual inconsistencies during the clicking-frame.
				if (ImGui::GetIO().KeyCtrl)
					selection_mask ^= (1 << node_clicked); // CTRL+click to toggle
				else									   // if (!(selection_mask & (1 << node_clicked))) // Depending on selection behavior you want, this commented bit preserve selection when clicking on item that is part of the selection
					selection_mask = (1 << node_clicked);  // Click to single-select

				viewer.selected_data_index = node_clicked;
			}

			ImGui::TreePop();
		}
		ImGui::PushStyleColor(ImGuiCol_Button, (ImVec4)ImColor::HSV(0.0f, 0.6f, 0.6f));
		ImGui::PushStyleColor(ImGuiCol_ButtonHovered, (ImVec4)ImColor::HSV(0.0f, 0.7f, 0.7f));
		ImGui::PushStyleColor(ImGuiCol_ButtonActive, (ImVec4)ImColor::HSV(0.0f, 0.8f, 0.8f));
		ImGui::PopStyleColor(3);
		ImGui::End();
	};
	igl::opengl::glfw::imgui::ImGuiMenu menu_mp;
	plugin.widgets.push_back(&menu_mp);

	// Add content to the default menu window
	menu_mp.callback_draw_viewer_menu = [&]()
	{
		// Draw parent menu content
		menu_mp.draw_viewer_menu();
	};

	// Draw additional windows
	menu_mp.callback_draw_custom_window = [&]()
	{
		// Define next window position + size
		ImGui::SetNextWindowPos(ImVec2(1100.f * menu_mp.menu_scaling(), 10), ImGuiCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(300, 600), ImGuiCond_FirstUseEver);
		ImGui::Begin(
			"Mesh Optimization", nullptr,
			ImGuiWindowFlags_NoSavedSettings);
		if (ImGui::CollapsingHeader("Quad Mesh Extraction", ImGuiTreeNodeFlags_DefaultOpen))
		{
			if (ImGui::Button("Load First Level Set", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{

				read_levelset(lscif::readed_LS1);
				std::cout<<"Levelset 1 get readed "<<std::endl;
			}
			ImGui::SameLine();
			if (ImGui::Button("Load Second Level Set", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{

				read_levelset(lscif::readed_LS2);
				std::cout<<"Levelset 2 get readed "<<std::endl;
			}
			ImGui::SameLine();
			if (ImGui::Button("Clear Readed Level Sets", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				lscif::readed_LS1.resize(0);
				lscif::readed_LS2.resize(0);
				std::cout << "The readed Level Sets Got Cleared" << std::endl;
			}
			if (ImGui::Button("Load First mesh", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				std::string fname = igl::file_dialog_open();
				if (fname.length() == 0)
				{
					std::cout << "\nLSC: read mesh failed" << std::endl;
					ImGui::End();
					return;
				}
				CGMesh readed_mesh;
				OpenMesh::IO::read_mesh(readed_mesh, fname);
				lscif::readed_mesh1.push_back(readed_mesh);
				std::cout << "mesh 1 get readed, current mesh nbr: "<< lscif::readed_mesh1.size() << std::endl;
			}
			ImGui::SameLine();
			if (ImGui::Button("Load Second mesh", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				std::string fname = igl::file_dialog_open();
				if (fname.length() == 0)
				{
					std::cout << "\nLSC: read mesh failed" << std::endl;
					ImGui::End();
					return;
				}

				OpenMesh::IO::read_mesh(lscif::readed_mesh2, fname);
				std::cout << "mesh 2 get readed " << std::endl;
			}
			ImGui::SameLine();
			if (ImGui::Button("ShadingTypes", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				Eigen::VectorXi info;
				Eigen::MatrixXd P1;
				Eigen::MatrixXd P2;
				if (lscif::readed_mesh1.empty())
				{
					std::cout << "\nLSC: Please read meshes" << std::endl;
					ImGui::End();
					return;
				}

				project_mesh_and_get_shading_info(lscif::readed_mesh1[0], lscif::readed_mesh2, lscif::ring_nbr, info, P1, P2);
				std::cout << "the first type is green, the second is red" << std::endl;

				viewer.data().add_points(P2, lscif::hot_red);
				viewer.data().add_points(P1, lscif::sea_green);
				std::cout << "Nbr first, " << P1.rows() << std::endl;
				std::cout << "Nbr second, " << P2.rows() << std::endl;

				std::cout<<"Passing the info to the level-set system"<<std::endl;
				lscif::tools.shading_condition_info = info;
			}
			ImGui::SameLine();
			if (ImGui::Button("PGTypes", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				Eigen::VectorXi info;
				Eigen::MatrixXd P1;
				Eigen::MatrixXd P2;
				if (lscif::readed_mesh1.empty())
				{
					std::cout << "\nLSC: Please read meshes" << std::endl;
					ImGui::End();
					return;
				}

				project_mesh_and_get_vertex_class(lscif::readed_mesh1, lscif::readed_mesh2, info);

				std::cout << "Passing the info to the level-set system" << std::endl;
				lscif::tools.shading_condition_info = info;
			}

			if (ImGui::CollapsingHeader(" Quad Mesh Size Paras", ImGuiTreeNodeFlags_DefaultOpen))
			{
				// Expose variable directly ...
				// ImGui::InputDouble("Closeness", &lscif::weigth_closeness, 0, 0, "%.4f");
				ImGui::InputInt("Nbr of U lines", &lscif::nbr_lines_first_ls, 0, 0);
				ImGui::InputInt("Nbr of V lines", &lscif::nbr_lines_second_ls, 0, 0);
				// ImGui::InputInt("Iteration", &lscif::OpIter, 0, 0);
				// ImGui::InputDouble("weight ls mass(big)", &lscif::weight_mass, 0, 0, "%.4f");

				// ImGui::Checkbox("Fix Boundary", &lscif::fixBoundary_checkbox);
			}

			
			// Add a button
			if (ImGui::Button("Otho Fields Quads", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				int id = viewer.selected_data_index;
				CGMesh updatedMesh;

				Eigen::MatrixXd VER;
				Eigen::MatrixXi FAC;
				if (lscif::readed_LS1.size() > 0 && lscif::readed_LS2.size() > 0)
				{
					extract_levelset_web(lscif::tools.lsmesh, lscif::tools.V, lscif::tools.F, lscif::readed_LS1, lscif::readed_LS2, lscif::nbr_lines_first_ls,
										 lscif::nbr_lines_second_ls, 3, VER, FAC, false);
				}
				else
				{
					std::cout << "ERROR, Please load both two level sets" << std::endl;
					ImGui::End();
					return;
				}
				std::cout<<"Saving the quad mesh"<<std::endl;
				std::string fname = igl::file_dialog_save();
				if (fname.length() == 0)
				{
					ImGui::End();
				}
				else
				{
					igl::writeOBJ(fname, VER, FAC);
					std::cout << "Quad mesh saved" << std::endl;
				}
				Eigen::MatrixXd E0, E1, E2, E3;
				visual_extract_levelset_web(lscif::tools.lsmesh, lscif::tools.V, lscif::tools.F, lscif::readed_LS1, lscif::readed_LS2, lscif::nbr_lines_first_ls,
										 lscif::nbr_lines_second_ls,E0, E1, E2, E3,false);
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);
				const Eigen::RowVector3d black(0, 0, 0);
				const Eigen::RowVector3d green(0.2, 0.8, 0.2);
				viewer.data().add_edges(E0, E1, red);
				viewer.data().add_edges(E2, E3, green);

			}
			ImGui::SameLine();
			if (ImGui::Button("Extr Pace Quad", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				std::cout<<"This code is to extract AAG quads"<<std::endl;
				int id = viewer.selected_data_index;
				CGMesh inputMesh = lscif::Meshes[id];
				Eigen::MatrixXd VER;
				Eigen::MatrixXi FAC;
				if (!lscif::debug_flag)
				{
					if (lscif::readed_LS1.size() > 0 && lscif::readed_LS2.size() > 0)
					{
						int filter_nbr = 3;
						extract_levelset_web_stable(lscif::tools.lsmesh, lscif::tools.Boundary_Edges, lscif::tools.V, lscif::tools.F, lscif::readed_LS1, lscif::readed_LS2, lscif::nbr_lines_first_ls,
											 lscif::nbr_lines_second_ls, filter_nbr, VER, FAC, true);
					}
					else
					{
						std::cout << "ERROR, Please load level sets" << std::endl;
					}
					std::string fname = igl::file_dialog_save();
					if (fname.length() == 0)
					{
						std::cout << "Quad mesh saving failed" << std::endl;
					}
					else
					{
						igl::writeOBJ(fname, VER, FAC);
						std::cout << "Quad mesh saved" << std::endl;
					}
				}

				Eigen::MatrixXd E0, E1, E2, E3;
				visual_extract_levelset_web_stable(lscif::tools.lsmesh,lscif::tools.Boundary_Edges, lscif::tools.V, lscif::tools.F, lscif::readed_LS1, lscif::readed_LS2, lscif::nbr_lines_first_ls,
											lscif::nbr_lines_second_ls, E0, E1, E2, E3, true, lscif::debug_flag, lscif::dbg_int, lscif::dbg_int2);
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);
				const Eigen::RowVector3d black(0, 0, 0);
				const Eigen::RowVector3d green(0.2, 0.8, 0.2);

				// lscif::MP.MeshUnitScale(inputMesh, updatedMesh);
				lscif::updateMeshViewer(viewer, inputMesh);
				lscif::meshFileName.push_back("paceqd_" + lscif::meshFileName[id]);
				lscif::Meshes.push_back(inputMesh);
				viewer.data().add_edges(E0, E1, red);
				viewer.data().add_edges(E2, E3, green);
				viewer.selected_data_index = id;
				
			}
			ImGui::SameLine();
			if (ImGui::Button("Extr Shading Lines", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				int id = viewer.selected_data_index;
				
				if(lscif::readed_LS1.size()!=lscif::tools.V.rows()){
					std::cout<<"Please load a correct Level Set 1"<<std::endl;
					ImGui::End();
					return;
				}
				extract_shading_lines(lscif::tools.lsmesh, lscif::tools.V, lscif::tools.Boundary_Edges,
									  lscif::tools.F, lscif::readed_LS1, lscif::nbr_lines_first_ls);
				// std::string fname = igl::file_dialog_save();
				// if (fname.length() == 0)
				// {
				// 	ImGui::End();
				// }
				// else
				// {
				// 	igl::writeOBJ(fname, VER, FAC);
				// 	std::cout << "Quad mesh saved" << std::endl;
				// }
			}
			ImGui::SameLine();
			if (ImGui::Button("Binormal Origami", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				int id = viewer.selected_data_index;
				CGMesh updatedMesh;
				
				Eigen::MatrixXd VER;
				Eigen::MatrixXi FAC;
				std::cout<<"reading Quad obj"<<std::endl;
				std::string fname = igl::file_dialog_open();
				if (fname.length() == 0)
				{
					std::cout << "\nLSC: read mesh failed" << std::endl;
					ImGui::End();
					return;
				}
				igl::readOBJ(fname, VER, FAC);
				std::cout<<"OBJ readed, ver in face "<<FAC.cols()<<std::endl;
				std::cout<<"Reading Bi-normals of triangle mesh"<<std::endl;
				Eigen::MatrixXd bn;
				read_bi_normals(bn);
				std::cout<<"Binormals readed"<<std::endl;
				std::cout<<"generating mesh with bi-normals"<<std::endl;
				fname = igl::file_dialog_save();
				
				if (fname.length() == 0)
				{
					ImGui::End();
					
				}
				
				write_quad_mesh_with_binormal(fname,lscif::tools.V, lscif::tools.F, bn, VER, FAC );
				
			}
			if (ImGui::Button("ExtrtQuad&IsoLines", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				int id = viewer.selected_data_index;
				
				Eigen::MatrixXd VER;
				Eigen::MatrixXi FAC;
				if(lscif::readed_LS1.size()!=lscif::tools.V.rows()){
					std::cout<<"Please load a correct Level Set 1"<<std::endl;
					ImGui::End();
					return;
				}
				extract_Quad_Mesh_Zigzag(lscif::tools.lsmesh, lscif::tools.Boundary_Edges, lscif::tools.V, lscif::tools.F, lscif::readed_LS1,
										 lscif::nbr_lines_first_ls, lscif::nbr_lines_second_ls, 3, VER, FAC);
				std::string fname = igl::file_dialog_save();
				if (fname.length() == 0)
				{
					ImGui::End();
				}
				else
				{
					igl::writeOBJ(fname, VER, FAC);
					std::cout << "Quad mesh saved" << std::endl;
				}
			}
			ImGui::SameLine();
			if (ImGui::Button("ShadingUnitScale", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				std::cout << "This button is to unit scale the target mesh and polylines and save the results\n ";
				std::cout << "Reading triangle mesh\n ";
				std::string fname = igl::file_dialog_open();
				if (fname.length() == 0)
				{
					std::cout << "\nLSC: read mesh failed" << std::endl;
					ImGui::End();
					return;
				}
				Eigen::MatrixXd Vtri, Vply;
				Eigen::MatrixXi Ftri;
				igl::readOBJ(fname, Vtri, Ftri);
				read_plylines_and_binormals(lscif::Poly_readed, lscif::Bino_readed);
				lscif::poly_tool.polyline_to_matrix(lscif::Poly_readed, Vply);

				int size1 = Vtri.rows();
				int size2 = Vply.rows();
				Eigen::MatrixXd Vall(size1 + size2, 3);
				Vall.topRows(size1) = Vtri;
				Vall.bottomRows(size2) = Vply;

				Eigen::MatrixXd vout;
				mesh_unit_scale(Vall, vout);
				Vtri = vout.topRows(size1);
				Vply = vout.bottomRows(size2);
				lscif::poly_tool.vertex_matrix_to_polyline(lscif::Poly_readed, Vply);

				std::cout << "Saving the result mesh and polylines, please provide the prifix" << std::endl;
				fname = igl::file_dialog_save();
				if (fname.length() == 0)
				{
					ImGui::End();
				}
				else
				{
					igl::writeOBJ(fname + ".obj", Vtri, lscif::tools.F);
					lscif::poly_tool.save_polyline_and_binormals_as_files(fname, lscif::Poly_readed, lscif::Bino_readed); // true means save the rotated polylines
					std::cout << "mesh saved" << std::endl;
				}
			}
		}

		ImGui::SameLine();
		if (ImGui::Button("ReadQuadPlyLines", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
		{
			int id = viewer.selected_data_index;
			read_origami_and_convert_to_polylines(lscif::Poly_readed, lscif::Bino_readed);
			CGMesh updatemesh = polyline_to_strip_mesh(lscif::Poly_readed, lscif::Bino_readed, lscif::vector_scaling);
			// std::cout<<"get the strip mesh"<<std::endl;
			lscif::updateMeshViewer(viewer, updatemesh);
			
			// std::cout<<"get the strip mesh"<<std::endl;
			lscif::meshFileName.push_back("ply_" + lscif::meshFileName[id]);
			lscif::Meshes.push_back(updatemesh);
			// std::cout<<"before init the poly opt"<<std::endl;
			lscif::poly_tool.init(lscif::Poly_readed, lscif::Bino_readed, -1); // we do not sample the polylines
		}
		ImGui::SameLine();
		if (ImGui::Button("UpdateQuadWithPlylines", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))

		{
			update_qd_mesh_with_plylines();
		}
		if (ImGui::CollapsingHeader("Mesh Optimization", ImGuiTreeNodeFlags_CollapsingHeader))
		{
			// Expose variable directly ...
			ImGui::InputInt("Iteration Mesh Opt", &lscif::Nbr_Iterations_Mesh_Opt, 0, 0);
			ImGui::InputDouble("Weight Mass", &lscif::weight_Mesh_mass, 0, 0, "%.4f");
			ImGui::InputDouble("Weight Smoothness", &lscif::weight_Mesh_smoothness, 0, 0, "%.4f");
			ImGui::InputDouble("Weight Mesh Pseudo Geodesic", &lscif::weight_Mesh_pesudo_geodesic, 0, 0, "%.4f");
			ImGui::InputDouble("Weight Mesh Egdge Length", &lscif::weight_Mesh_edgelength, 0, 0, "%.4f");
			ImGui::InputDouble("Weight Approx", &lscif::weight_Mesh_approximation, 0, 0, "%.4f");
			ImGui::InputDouble("Max Mesh Step Length", &lscif::Mesh_opt_max_step_length, 0, 0, "%.4f");
			if (ImGui::Button("DrawProjection", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);
				const Eigen::RowVector3d yellow(241. / 255, 196. / 255, 15. / 255);

				const Eigen::RowVector3d black(0, 0, 0);
				const Eigen::RowVector3d green(0.2, 0.8, 0.2);

				if (lscif::tools.Ppro0.size() != 0)
				{
					viewer.data().add_edges(lscif::tools.Ppro0, lscif::tools.V, green);
					// viewer.data().add_edges(lscif::tools.Ppro0,
					// 						lscif::tools.Ppro0 + lscif::vector_scaling * lscif::tools.Npro0, red);
					viewer.data().add_edges(lscif::tools.Ppro0,
											lscif::tools.Ppro0 + lscif::vector_scaling * lscif::tools.Npro0, red);
				}
			}
			// ImGui::InputInt("Iteration", &lscif::OpIter, 0, 0);
			// ImGui::InputDouble("weight ls mass(big)", &lscif::weight_mass, 0, 0, "%.4f");
			// ImGui::Checkbox("Fix Boundary", &lscif::fixBoundary_checkbox);
		}
		if (ImGui::CollapsingHeader("Shading Paras", ImGuiTreeNodeFlags_CollapsingHeader))
		{
			ImGui::InputInt("Ver_id", &lscif::update_ver_id, 0, 0);
			ImGui::InputInt("Ring nbr", &lscif::ring_nbr, 0, 0);
			ImGui::InputDouble("Latitude", &lscif::Shading_Latitude, 0, 0, "%.4f");
			if (ImGui::Button("UpdateVertexID", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				lscif::update_verlist.push_back(lscif::update_ver_id);
				std::cout << "Current ver list" << std::endl;
				for (int i = 0; i < lscif::update_verlist.size(); i++)
				{
					std::cout << lscif::update_verlist[i] << ",";
				}
				std::cout << "\n";
			}
			ImGui::SameLine();
			if (ImGui::Button("LoadVertexID", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				std::cout<<"Read Vertex ID file row by row"<<std::endl;
				std::string fname = igl::file_dialog_open();
				if (fname.length() == 0)
				{
					std::cout << "\nLSC: read file failed" << std::endl;
					ImGui::End();
					return;
				}
				std::vector<std::vector<double>> readed;
				read_csv_data_lbl(fname, readed);
				std::vector<int> list(readed[0].size());
				for (int i = 0; i < list.size(); i++)
				{
					list[i] = readed[0][i];
				}
				lscif::update_verlist = list;
				for (int i = 0; i < lscif::update_verlist.size(); i++)
				{
					std::cout << lscif::update_verlist[i] << ",";
				}
				std::cout << "\n";

			}
			ImGui::SameLine();
			if (ImGui::Button("ClearVertexList", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				lscif::update_verlist.clear();
				std::cout << "ver list got cleared" << std::endl;
			}
			ImGui::SameLine();
			if (ImGui::Button("SaveRotatedMesh", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				std::cout << "Rotate the mesh in the earth axis according to the input latitude.\nNow saving the mesh" << std::endl;
				Eigen::MatrixXd Vnew;
				rotate_z_axis_to_earth_axis(lscif::tools.V, Vnew, lscif::Shading_Latitude);

				std::string fname = igl::file_dialog_save();
				if (fname.length() == 0)
				{
					std::cout << "\nLSC: save mesh failed" << std::endl;
					ImGui::End();
					return;
				}
				igl::writeOBJ(fname, Vnew, lscif::tools.F);
				std::cout << "mesh saved" << std::endl;
			}
			ImGui::Checkbox("Shading", &lscif::Given_Const_Direction);
			ImGui::SameLine();
			ImGui::Checkbox("Light Through", &lscif::let_ray_through);
			ImGui::SameLine();
			ImGui::Checkbox("Light Reflect", &lscif::let_ray_reflect);
			
			if (lscif::Given_Const_Direction)
			{
				ImGui::PushItemWidth(50);
				ImGui::InputDouble("T0", &lscif::InputPx, 0, 0, "%.6f");
				ImGui::SameLine();
				ImGui::InputDouble("P0", &lscif::InputPy, 0, 0, "%.6f");
				ImGui::SameLine();
				ImGui::InputDouble("Ttol0", &lscif::InputThetaTol, 0, 0, "%.6f");
				ImGui::SameLine();
				ImGui::InputDouble("Ptol0", &lscif::InputPhiTol, 0, 0, "%.6f");

				ImGui::InputDouble("T1", &lscif::InputPx1, 0, 0, "%.6f");
				ImGui::SameLine();
				ImGui::InputDouble("P1", &lscif::InputPy1, 0, 0, "%.6f");
				ImGui::SameLine();
				ImGui::InputDouble("Ttol1", &lscif::InputThetaTol1, 0, 0, "%.6f");
				ImGui::SameLine();
				ImGui::InputDouble("Ptol1", &lscif::InputPhiTol1, 0, 0, "%.6f");
				ImGui::Checkbox("MarkMaxEnergy", &lscif::enable_max_energy_check);
				ImGui::SameLine();
				if (ImGui::Button("ClearMaxEnergy", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
				{
					lscif::tools.clear_high_energy_markers_in_analizer();
				}
				ImGui::SameLine();
				ImGui::InputDouble("MaxEnergyPercent", &lscif::max_e_percentage, 0, 0, "%.6f");
				ImGui::Checkbox("ShadingInit", &lscif::shading_init);
				ImGui::SameLine();
				ImGui::InputDouble("weight binormal", &lscif::weight_binormal, 0, 0, "%.4f");
				ImGui::InputDouble("weight smtBinormal", &lscif::weight_smt_binormal, 0, 0, "%.4f");
				ImGui::Checkbox("RecomptAuxiliaries", &lscif::recompute_auxiliaries);
			}
			
		}
		ImGui::PushItemWidth(50);
		ImGui::Checkbox("pickPly", &lscif::pick_single_ply);
		ImGui::SameLine();
		ImGui::InputInt("pickID", &lscif::pick_line_id, 0, 0);
		ImGui::Checkbox("2LvStAngle", &lscif::fix_angle_of_two_levelsets);
		ImGui::SameLine();
		ImGui::InputDouble("AngleOfLSs", &lscif::angle_between_two_levelsets, 0, 0, "%.4f");
		ImGui::SameLine();
		ImGui::InputDouble("WeightAngleOfLSs", &lscif::weight_fix_two_ls_angle, 0, 0, "%.4f");
		
		ImGui::Checkbox("drawStrokes", &lscif::draw_strokes);
		ImGui::SameLine();
		if (ImGui::Button("ClearStrokes", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))

		{
			lscif::project_2dx.clear();
			lscif::project_2dy.clear();
			lscif::project_x_tmp.clear();
			lscif::project_y_tmp.clear();
			lscif::project_pts.clear();
			lscif::flist.clear();
			lscif::bclist.clear();
			CGMesh updatemesh = lscif::tools.lsmesh;
			int id = viewer.selected_data_index;
			lscif::updateMeshViewer(viewer, updatemesh);
			lscif::meshFileName.push_back(lscif::meshFileName[id]);
			lscif::Meshes.push_back(updatemesh);
			std::cout<<"The Strokes Are Cleared"<<std::endl;
		}
		ImGui::SameLine();
		if (ImGui::Button("Strokes2LS", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))

		{
			std::cout<<"Using Strokes to Initialize the Level Set"<<std::endl;
			lscif::tools.receive_interactive_strokes_and_init_ls(lscif::flist, lscif::bclist);
			bool compute_pg = false;
			AssignAutoRunDefaultArgs(lscif::autorunner[0], compute_pg);
			compute_pg = true;
			AssignAutoRunDefaultArgs(lscif::autorunner[1], compute_pg);
			std::cout<<"Using Default Parameters"<<std::endl;
		}
		
		if (ImGui::Button("AutoRunSmoothLS", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))

		{
			AutoRunArgs runner = lscif::autorunner[0];
			for (int i = 0; i < runner.parts.size(); i++)
			{
				lscif::tools.weight_mass = runner.parts[i].weight_gravity;
				lscif::tools.weight_laplacian = runner.parts[i].weight_lap;
				lscif::tools.weight_boundary = runner.parts[i].weight_bnd;
				lscif::tools.weight_pseudo_geodesic_energy = runner.parts[i].weight_pg;
				lscif::tools.weight_strip_width = runner.parts[i].weight_strip_width;
				for (int j = 0; j < runner.parts[i].iterations; j++)
				{
					lscif::tools.Run_Level_Set_Opt_interactive(runner.compute_pg);
					if (lscif::tools.step_length < runner.stop_step_length && j != 0)
					{ 
						std::cout << "optimization converges " << std::endl;
						break;
					}
				}
			}

			std::cout << "A Smooth Scalar Function Following Input Strokes Is Generated ..." << std::endl;
			int id = viewer.selected_data_index;
			CGMesh inputMesh = lscif::tools.lsmesh;
			lscif::updateMeshViewer(viewer, inputMesh);
			lscif::meshFileName.push_back("lso_" + lscif::meshFileName[id]);
			lscif::Meshes.push_back(inputMesh);

			Eigen::VectorXd level_set_values;
			lscif::tools.show_level_set(level_set_values);
			if (level_set_values.size() == 0)
			{
				std::cout<<"Level Set is not initialized"<<std::endl;
				ImGui::End();
				return;
			}
			Eigen::MatrixXd CM;
			igl::parula(Eigen::VectorXd::LinSpaced(21, 0, 1).eval(), false, CM);
			igl::isolines_map(Eigen::MatrixXd(CM), CM);
			viewer.data().set_colormap(CM);
			viewer.data().set_data(level_set_values);
			const Eigen::RowVector3d red(0.8, 0.2, 0.2);
			const Eigen::RowVector3d blue(0.2, 0.2, 0.8);

			viewer.selected_data_index = id;
		}
		ImGui::SameLine();
		if (ImGui::Button("AutoRunPd-Gdsic", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))

		{
			AutoRunArgs runner = lscif::autorunner[1];
			for (int i = 0; i < runner.parts.size(); i++)
			{
				lscif::tools.weight_mass = runner.parts[i].weight_gravity;
				lscif::tools.weight_laplacian = runner.parts[i].weight_lap;
				lscif::tools.weight_boundary = runner.parts[i].weight_bnd;
				lscif::tools.weight_pseudo_geodesic_energy = runner.parts[i].weight_pg;
				lscif::tools.weight_strip_width = runner.parts[i].weight_strip_width;
				for (int j = 0; j < runner.parts[i].iterations; j++)
				{
					lscif::tools.Run_Level_Set_Opt_interactive(runner.compute_pg);
					if (lscif::tools.step_length < runner.stop_step_length && j != 0)
					{ 
						std::cout << "optimization converges " << std::endl;
						break;
					}
				}
				std::cout<<"\n";
			}

			std::cout << "Pseudo-Geodesic Curves in level sets are Generated ..." << std::endl;
			int id = viewer.selected_data_index;
			CGMesh inputMesh = lscif::tools.lsmesh;
			lscif::updateMeshViewer(viewer, inputMesh);
			lscif::meshFileName.push_back("lso_" + lscif::meshFileName[id]);
			lscif::Meshes.push_back(inputMesh);

			Eigen::VectorXd level_set_values;
			lscif::tools.show_level_set(level_set_values);
			if (level_set_values.size() == 0)
			{
				std::cout<<"Level Set is not initialized"<<std::endl;
				ImGui::End();
				return;
			}
			Eigen::MatrixXd CM;
			igl::parula(Eigen::VectorXd::LinSpaced(21, 0, 1).eval(), false, CM);
			igl::isolines_map(Eigen::MatrixXd(CM), CM);
			viewer.data().set_colormap(CM);
			viewer.data().set_data(level_set_values);
			const Eigen::RowVector3d red(0.8, 0.2, 0.2);
			const Eigen::RowVector3d blue(0.2, 0.2, 0.8);

			viewer.selected_data_index = id;
		}
		ImGui::SameLine();
		if (ImGui::Button("SaveStrokes", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))

		{
			save_strokes(lscif::flist, lscif::bclist);
		}
		ImGui::SameLine();
		if (ImGui::Button("ReadStrokes", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))

		{
			Eigen::MatrixXd V = lscif::tools.V;
			Eigen::MatrixXi F = lscif::tools.F;
			read_strokes(lscif::flist, lscif::bclist);
			Eigen::MatrixXd vlist;
			for (int i = 0; i < lscif::flist.size(); i++)
			{
				vlist.resize(lscif::flist[i].size(), 3);
				for (int j = 0; j < lscif::flist[i].size(); j++)
				{
					int fid = lscif::flist[i][j];
					Eigen::Vector3f bc = lscif::bclist[i][j];
					Eigen::Vector3d v0 = V.row(F(fid, 0));
					Eigen::Vector3d v1 = V.row(F(fid, 1));
					Eigen::Vector3d v2 = V.row(F(fid, 2));

					Eigen::Vector3d pt = bc[0] * v0 + bc[1] * v1 + bc[2] * v2;
					vlist.row(j) = pt;
				}
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				viewer.data().add_points(vlist, red); // show rows
			}
		}
		if (ImGui::CollapsingHeader("QuadMeshOpt", ImGuiTreeNodeFlags_DefaultOpen))
		{
			// ImGui::InputInt("Ver_id", &lscif::update_ver_id, 0, 0);
			// ImGui::InputInt("Ring nbr", &lscif::ring_nbr, 0, 0);
			// ImGui::InputDouble("Latitude", &lscif::Shading_Latitude, 0, 0, "%.4f");
			ImGui::Combo("QuadType", &lscif::quad_tool.OptType,
                 "AAG\0GGA\0PP\0PPG\0\0");
			ImGui::SameLine();
			
			ImGui::Combo("FamilyOfDiag", &lscif::quad_tool.WhichDiagonal,
                 "D0\0D1\0\0");
			

			if (ImGui::Button("LoadQuads", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				std::string fname = igl::file_dialog_open();
				if (fname.length() == 0)
				{
					std::cout << "\nLSC: read mesh failed" << std::endl;
					ImGui::End();
					return;
				}
				CGMesh quadmesh;
				OpenMesh::IO::read_mesh(quadmesh, fname);
				std::cout << "\nMesh Readed. Now Please Type Down the Prefix for Row and Column Info" << std::endl;
				std::string savefile = igl::file_dialog_save();
				if (savefile.length() == 0)
				{
					std::cout << "\nPlease Type Down Something For The Row And Col Info" << std::endl;
					ImGui::End();
					return;
				}
				lscif::quad_tool.init(quadmesh, savefile);

				size_t last_dot = fname.rfind('.');
				size_t last_slash = fname.rfind(spliter); // TODO on linux it should be '/'
				std::string fnameNew = fname.substr(last_slash + 1, (last_dot - last_slash - 1));
				lscif::meshFileName.push_back(fnameNew);
				lscif::Meshes.push_back(quadmesh);

				int nbe = 0;
				for (CGMesh::EdgeIter eit = quadmesh.edges_begin(); eit != quadmesh.edges_end(); ++eit)
				{
					if (quadmesh.is_boundary(eit))
						nbe++;
				}

				int nbhe = 0;
				for (CGMesh::HalfedgeIter heit = quadmesh.halfedges_begin(); heit != quadmesh.halfedges_end(); ++heit)
				{
					if (quadmesh.is_boundary(heit))
						nbhe++;
				}

				std::cout << "Mesh is: " << fname << std::endl;

				std::cout << "Vertices/Edges/Faces/HalfEdges/boundaryHalfEdge/boundaryEdge: " << quadmesh.n_vertices() << "/" << quadmesh.n_edges() << "/" << quadmesh.n_faces() << "/" << quadmesh.n_halfedges() << "/" << nbhe << "/" << nbe << std::endl;

				lscif::updateMeshViewer(viewer, quadmesh);
			}
			ImGui::SameLine();
			
			if (ImGui::Button("Set&Reset", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				lscif::quad_tool.reset();
				std::array<Eigen::MatrixXd, 3> Edges;
				lscif::quad_tool.show_curve_families(Edges);
				int id = viewer.selected_data_index;
				lscif::updateMeshViewer(viewer, lscif::quad_tool.mesh_original);
				lscif::meshFileName.push_back("Rst_" + lscif::meshFileName[id]);
				lscif::Meshes.push_back(lscif::quad_tool.mesh_original);
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);
				const Eigen::RowVector3d black(0, 0, 0);
				const Eigen::RowVector3d green(0.2, 0.8, 0.2);

				viewer.data().add_points(Edges[0], red);// show rows 
				viewer.data().add_points(Edges[1], green);// cols
				viewer.data().add_points(Edges[2], blue);// diagonals
				// std::cout<<"Visualing Rows Cols and Diagonals : "<<Edges[0]<<",\n\n "<<Edges[1]<<",\n\n "<<Edges[2]<<"\n";
				std::cout<<"The rows, columns and the diagonals are marked in red, green and blue"<<std::endl;;
				viewer.selected_data_index = id;
			}
			ImGui::SameLine();
			
			if (ImGui::Button("Opt", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				lscif::quad_tool.weight_fairness = lscif::weight_laplacian;
				lscif::quad_tool.weight_gravity = lscif::weight_boundary;
				lscif::quad_tool.weight_pg = lscif::weight_pseudo_geodesic;
				lscif::quad_tool.angle_degree0 = lscif::target_angle;
				lscif::quad_tool.angle_degree1 = lscif::target_angle_2;
				lscif::quad_tool.pg_ratio = lscif::weight_geodesic;
				lscif::quad_tool.weight_mass = lscif::weight_mass;
				lscif::quad_tool.max_step = lscif::maximal_step_length;
				if (lscif::quad_tool.V.rows()==0)
				{
					std::cout << "\nEmpty quad, please load a quad mesh first" << std::endl;
					ImGui::End();
					return;
				}
				for (int i = 0; i < lscif::OpIter; i++)
				{
					lscif::quad_tool.opt();
					if (lscif::quad_tool.real_step_length < 1e-16 && i != 0)
					{ // step length actually is the value for the last step
						std::cout << "optimization converges " << std::endl;
						break;
					}
				}
				std::cout << "waiting for instruction..." << std::endl;
				// lscif::MP.MeshUnitScale(inputMesh, updatedMesh);
				int id = viewer.selected_data_index;
				CGMesh updateMesh = lscif::quad_tool.mesh_update;
				lscif::updateMeshViewer(viewer, updateMesh);
				lscif::meshFileName.push_back("lso_" + lscif::meshFileName[id]);
				lscif::Meshes.push_back(updateMesh);
				viewer.selected_data_index = id;
			}
			ImGui::SameLine();
			
			if (ImGui::Button("DrawNormal", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);
				const Eigen::RowVector3d yellow(241./255, 196./255, 15./255);

				const Eigen::RowVector3d black(0, 0, 0);
				const Eigen::RowVector3d green(0.2, 0.8, 0.2);

				Eigen::MatrixXd pts, E0, E1;

				std::cout << "ploting the normal vectors" << std::endl;
				Eigen::MatrixXd binormals;
				if (lscif::quad_tool.Enm0.rows() == 0)
				{
					std::cout << "\nWrong Usage ..." << std::endl;
					ImGui::End();
					return;
				}

				viewer.data().add_edges(lscif::quad_tool.Enm0,
										lscif::quad_tool.Enm0 + (lscif::quad_tool.Enm1 - lscif::quad_tool.Enm0) * lscif::vector_scaling, green);
			}
			if (ImGui::Button("writePlyInfo", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				lscif::quad_tool.write_polyline_info();

			}
			ImGui::SameLine();
			if (ImGui::Button("LoadTriangTree", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				lscif::quad_tool.load_triangle_mesh_tree(lscif::tools.aabbtree, lscif::tools.Vstored,
														 lscif::tools.F, lscif::tools.Nstored);
			}
		}
		if (ImGui::CollapsingHeader("FunctionalAngles", ImGuiTreeNodeFlags_DefaultOpen))
		{
			// ImGui::InputInt("Ver_id", &lscif::update_ver_id, 0, 0);
			// ImGui::InputInt("Ring nbr", &lscif::ring_nbr, 0, 0);
			// ImGui::InputDouble("Latitude", &lscif::Shading_Latitude, 0, 0, "%.4f");
			ImGui::Checkbox("DisableIncreasing", &lscif::Disable_Changed_Angles);
			ImGui::SameLine();
			ImGui::Checkbox("UseFittingAngles", &lscif::Use_Fitting_Angles);
			ImGui::SameLine();
			ImGui::Checkbox("Opt_Only_BNMS", &lscif::Use_Opt_Only_BNMS);
			
			

			if (ImGui::Button("OrthoAsPsblOpt", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				//binormals as orthogonal as possible.
				EnergyPrepare einit;
				einit.weight_gravity = lscif::weight_mass;
				einit.weight_lap = lscif::weight_laplacian;
				einit.weight_bnd = lscif::weight_boundary;
				einit.weight_pg = lscif::weight_pseudo_geodesic;
				einit.weight_strip_width = lscif::weight_strip_width;
				einit.solve_pseudo_geodesic = lscif::enable_pg_energy_checkbox;
				einit.target_angle = lscif::target_angle;
				einit.max_step_length = lscif::maximal_step_length;
				einit.solve_strip_width_on_traced = lscif::enable_strip_width_checkbox;
				einit.enable_inner_vers_fixed = lscif::enable_inner_vers_fixed;
				einit.enable_functional_angles = lscif::enable_functional_degrees;
				einit.enable_extreme_cases = lscif::enable_extreme_cases;
				einit.target_min_angle = lscif::target_min_angle;
				einit.target_max_angle = lscif::target_max_angle;
				einit.start_angle = lscif::start_angle;
				einit.enable_boundary_angles = lscif::enable_boundary_angles;
				einit.Given_Const_Direction = lscif::Given_Const_Direction;
				// lscif::tools.weight_shading = lscif::weight_shading;
				lscif::tools.prepare_level_set_solving(einit);
				lscif::tools.Second_Ray_vers = lscif::update_verlist;
				lscif::tools.Second_Ray_nbr_rings = lscif::ring_nbr;
				lscif::tools.Reference_theta = lscif::InputPx;
				lscif::tools.Reference_phi = lscif::InputPy;
				lscif::tools.Reference_theta1 = lscif::InputPx1;
				lscif::tools.Reference_phi1 = lscif::InputPy1;
				lscif::tools.Theta_tol = lscif::InputThetaTol;
				lscif::tools.Phi_tol = lscif::InputPhiTol;
				lscif::tools.Theta_tol1 = lscif::InputThetaTol1;
				lscif::tools.Phi_tol1 = lscif::InputPhiTol1;
				lscif::tools.weight_geodesic=lscif::weight_geodesic;
				lscif::tools.enable_max_energy_check = lscif::enable_max_energy_check;
				lscif::tools.max_energy_percentage = lscif::max_e_percentage;
				lscif::tools.enable_shading_init = lscif::shading_init;
				lscif::tools.weight_binormal = lscif::weight_binormal;
				lscif::tools.weight_smt_binormal = lscif::weight_smt_binormal;
				lscif::tools.enable_let_ray_through = lscif::let_ray_through;
				lscif::tools.ShadingLatitude = lscif::Shading_Latitude;
				lscif::tools.enable_reflection = lscif::let_ray_reflect;
				lscif::tools.recompute_auxiliaries = lscif::recompute_auxiliaries;
				lscif::tools.Disable_Changed_Angles = lscif::Disable_Changed_Angles;
				lscif::tools.Use_Fitting_Angles = lscif::Use_Fitting_Angles;
				lscif::tools.Use_Opt_Only_BNMS = lscif::Use_Opt_Only_BNMS;

				for (int i = 0; i < lscif::OpIter; i++)
				{
					lscif::tools.Run_AsOrthAsPossible_LS();
					if (lscif::tools.step_length < 1e-16 && i != 0)
					{ // step length actually is the value for the last step
						std::cout << "optimization converges " << std::endl;
						break;
					}
				}
				std::cout << "waiting for instruction..." << std::endl;
				// lscif::MP.MeshUnitScale(inputMesh, updatedMesh);
				int id = viewer.selected_data_index;
				CGMesh inputMesh = lscif::tools.lsmesh;
				lscif::updateMeshViewer(viewer, inputMesh);
				lscif::meshFileName.push_back("AOAP_" + lscif::meshFileName[id]);
				lscif::Meshes.push_back(inputMesh);

				Eigen::VectorXd level_set_values;
				lscif::tools.show_level_set(level_set_values);
				if(level_set_values.size()==0){
					ImGui::End();
					return;
				}
				Eigen::MatrixXd CM;
				// std::cout<<"before compute colormap"<<std::endl;
				igl::parula(Eigen::VectorXd::LinSpaced(21, 0, 1).eval(), false, CM);
				igl::isolines_map(Eigen::MatrixXd(CM), CM);
				// std::cout<<"before set colormap"<<std::endl;
				viewer.data().set_colormap(CM);
				// std::cout<<"before set color"<<std::endl;
				viewer.data().set_data(level_set_values);
				// Eigen::MatrixXd E0, E1;
				// // lscif::tools.show_gradients(E0,E1, lscif::vector_scaling);
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);

				viewer.selected_data_index = id;
			}
			ImGui::SameLine();
			if (ImGui::Button("WriteFitting", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				lscif::tools.write_fitting_data(); 
			}
		}
		ImGui::PushStyleColor(ImGuiCol_Button, (ImVec4)ImColor::HSV(0.0f, 0.6f, 0.6f));
		ImGui::PushStyleColor(ImGuiCol_ButtonHovered, (ImVec4)ImColor::HSV(0.0f, 0.7f, 0.7f));
		ImGui::PushStyleColor(ImGuiCol_ButtonActive, (ImVec4)ImColor::HSV(0.0f, 0.8f, 0.8f));
		ImGui::PopStyleColor(3);
		ImGui::End();
	};

	// Refresh selected mesh colors
	viewer.callback_pre_draw =
		[&](igl::opengl::glfw::Viewer &)
	{
		if (last_selected != viewer.selected_data_index)
		{
			for (auto &data : viewer.data_list)
			{
				data.set_colors(colors[data.id]);
			}
			viewer.data_list[viewer.selected_data_index].set_colors(Eigen::RowVector3d(0.9, 0.9, 0.1));
			last_selected = viewer.selected_data_index;
		}
		return false;
	};

	// Plot the mesh
	viewer.data().set_mesh(lscif::V, lscif::F);
	viewer.data().show_lines = false;
	viewer.data().is_visible = true;
	// show polygon edges
	Eigen::MatrixXi E;
	lscif::MP.meshEdges(lscif::mesh, E);
	Eigen::MatrixXd C = Eigen::MatrixXd::Zero(E.rows(), 3);
	viewer.data().set_edges(lscif::V, E, C);
	viewer.data().line_width = 2.0;
	viewer.data().invert_normals = true;
	lscif::VertexType.resize(lscif::V.size(), 0);

	viewer.callback_pre_draw = &lscif::pre_draw;
	viewer.callback_key_down = &lscif::key_down;
	viewer.callback_mouse_down = &lscif::mouse_down;
	viewer.callback_key_up = &lscif::key_up;
	viewer.callback_mouse_move = &lscif::mouse_move;
	viewer.callback_mouse_up = &lscif::mouse_up;
	viewer.core().is_animating = false;
	viewer.core().animation_max_fps = 30.;
	viewer.core().background_color = Eigen::Vector4f(1, 1, 1, 100);
	viewer.launch();
}
