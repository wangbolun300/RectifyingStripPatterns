#include <igl/readOFF.h>
#include <igl/slice.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <igl/unproject_onto_mesh.h>
#include <igl/arap.h>
#include <igl/harmonic.h>
#include <imgui/imgui.h>
#include <imgui/imgui.h>
#include <iostream>

#include <igl/unproject_ray.h>

#include <lsc/MeshProcessing.h>
#include <lsc/OptimizerCheckboard.h>
#include <lsc/basic.h>
// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
typedef OpenMesh::PolyMesh_ArrayKernelT<> CGMesh;

// the interface operators
namespace lscif
{

	const Eigen::RowVector3d sea_green(70. / 255., 252. / 255., 167. / 255.);
	const Eigen::RowVector3d hot_red(255. / 255., 12. / 255., 17. / 255.);
	MeshProcessing MP;

	// This V and F are only for visulization of the default mesh
	Eigen::MatrixXd V;
	Eigen::MatrixXi F;

	bool keyPress_1 = false;
	bool keyPress_2 = false;

	// Optimization Parameters
	double weigth_closeness = 0.0;
	double refer_AveEL = 1;
	double weight_mass = 100.;			  // weight of the mass function to
	double weight_boundary=100;	// weight of the boundary condition
	double weight_laplacian=1;
	int which_seg_id = 0;			  // the face id of which we will assign value to
	int nbr_of_intervals = 3;
	int id_debug_global=5;
	double target_angle=60;
	double start_angle=60;
	double threadshold_angel_degree=120;// the threadshold for detecting boundary corners
	int dbg_int;
	int dbg_int2 = 5;
	double dbg_dbl = 30;
	double dbg_dbl2 = 60;
	double vector_scaling = 1;
	double strip_width_ratio = 1.;

	std::vector<int> fixedVertices;

	int OpIter = 10;
	CGMesh mesh;
	CGMesh RefMesh;
	static bool fixCorner_checkbox = false;
	static bool fixBoundary_checkbox = false;
	static bool selectFEV[30] = {true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true, true};
	std::vector<std::string> meshFileName;
	std::vector<CGMesh> Meshes;
	lsTools tools;
	double InputPx = 0.0;
	double InputPy = 0.0;
	double InputPz = 0.0;
	std::vector<int> VertexType;

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
				std::cout << "vid, " << vid << ", fid, " << fid << std::endl;
				if (lscif::tools.fvalues.rows() > 0)
				{
					std::cout << "level set value " << lscif::tools.fvalues[vid] << std::endl;
				}
				std::cout << "point position: (" << viewer.data().V(vid, 0) << ", " << viewer.data().V(vid, 1) << ", " << viewer.data().V(vid, 2) << ")\n\n";
				// viewer.data().add_points(igl::slice(viewer.data().V, vids, 1), hot_red);
				viewer.data().set_points(igl::slice(viewer.data().V, vids, 1), lscif::hot_red);
				viewer.data().add_points(igl::slice(viewer.data().V, vids2, 1), lscif::sea_green);
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
	std::string fname = example_root_path + std::string("oc_5.000000_10.000000_50_30.obj");
	OpenMesh::IO::read_mesh(lscif::mesh, fname);

	lscif::MP.mesh2Matrix(lscif::mesh, lscif::V, lscif::F);
	lscif::meshFileName.push_back("quad_pq");
	lscif::Meshes.push_back(lscif::mesh);

	// Init the viewer
	igl::opengl::glfw::Viewer viewer;

	// Attach a menu plugin
	igl::opengl::glfw::imgui::ImGuiMenu menu;
	viewer.plugins.push_back(&menu);

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
		ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(300, 800), ImGuiSetCond_FirstUseEver);
		ImGui::Begin(
			"Levelset Curves", nullptr,
			ImGuiWindowFlags_NoSavedSettings);

		// Add a button
		if (ImGui::Button("Import Mesh", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
		{

			std::string fname = igl::file_dialog_open();
			if (fname.length() == 0)
				return;

			OpenMesh::IO::read_mesh(lscif::mesh, fname);

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
				return;
			fname = fname + ".obj";

			int id = viewer.selected_data_index;
			OpenMesh::IO::write_mesh(lscif::Meshes[id], fname);
		}

		ImGui::SameLine();

		if (ImGui::Button("Import Para", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
		{
			std::string fname = igl::file_dialog_open();

			if (fname.length() == 0)
				return;

			std::ifstream settings(fname);
			settings >> lscif::weigth_closeness;
		}

		ImGui::SameLine();

		if (ImGui::Button("Save Para", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
		{
			std::string fname = igl::file_dialog_save();

			std::ofstream file;
			file.open(fname);

			file << lscif::weigth_closeness << std::endl;

			file.close();
		}

		static bool updateMeshFlag = false;
		static int updateMeshId = 0;
		// when clicking this button, it will print the info of this mesh
		if (ImGui::Button("Update Reference Mesh", ImVec2(ImGui::GetWindowSize().x * 0.35f, 0.0f)))
		{
			updateMeshFlag = true;
			updateMeshId = viewer.selected_data_index;
			lscif::RefMesh = lscif::Meshes[updateMeshId];
			double refAveEL = 0;
			double totalL, minL, maxL;
			lscif::MP.getEdgeLengthStatistics(lscif::RefMesh, totalL, minL, maxL, refAveEL);
			bool printInfo = true;
			if (printInfo)
			{
				//	std::cout << "Mesh Informatin Vertices/Edges/Faces/HalfEdges:::::" << RefMesh.n_vertices() << "/" << RefMesh.n_edges() << "/" << RefMesh.n_faces() << "/" << RefMesh.n_halfedges() << std::endl;
				std::cout << "Edge lengths Total/minL/maxL/averageL:::::" << totalL << "/" << minL << "/" << maxL << "/" << refAveEL << std::endl;
			}
			lscif::refer_AveEL = refAveEL;
		}
		// show the name of the mesh when printing the info of this mesh
		if (updateMeshFlag)
		{
			ImGui::SameLine();
			std::string str = lscif::meshFileName[updateMeshId];
			const char *meshnameChar = str.c_str();
			ImGui::Text(meshnameChar);
		}

		// Add new group
		if (ImGui::CollapsingHeader("Weights", ImGuiTreeNodeFlags_DefaultOpen))
		{
			// Expose variable directly ...
			// ImGui::InputDouble("Closeness", &lscif::weigth_closeness, 0, 0, "%.4f");
			ImGui::InputInt("Iteration", &lscif::OpIter, 0, 0);
			ImGui::InputDouble("weight ls mass(big)", &lscif::weight_mass, 0, 0, "%.4f");
			ImGui::InputDouble("weight boundary (big)", &lscif::weight_boundary, 0, 0, "%.4f");
			ImGui::InputDouble("weight laplacian (small)", &lscif::weight_laplacian, 0, 0, "%.4f");
			ImGui::InputInt("which boundary segment id", &lscif::which_seg_id, 0, 0);
			ImGui::InputInt("every i segments", &lscif::nbr_of_intervals, 0, 0);
			ImGui::InputInt("debug id", &lscif::id_debug_global, 0, 0);
			ImGui::InputDouble("target angle", &lscif::target_angle, 0, 0, "%.4f");
			ImGui::InputDouble("start angle", &lscif::start_angle, 0, 0, "%.4f");
			ImGui::InputDouble("threadshold angle", &lscif::threadshold_angel_degree, 0, 0, "%.4f");
			// ImGui::InputDouble("Average Edge length", &lscif::refer_AveEL, 0, 0, "%.4f");
			ImGui::InputDouble("Strip Width Ratio", &lscif::strip_width_ratio, 0, 0, "%.4f");

			// ImGui::Checkbox("Fix Corners", &lscif::fixCorner_checkbox);
			// ImGui::SameLine();
			// ImGui::Checkbox("Fix Boundary", &lscif::fixBoundary_checkbox);
		}

		ImGui::Separator();

		if (ImGui::CollapsingHeader("Mesh Processing", ImGuiTreeNodeFlags_DefaultOpen))
		{
			if (ImGui::CollapsingHeader("Debug input", ImGuiTreeNodeFlags_DefaultOpen))
			{
				// Expose variable directly ...
				ImGui::InputInt("int", &lscif::dbg_int, 0, 0);
				ImGui::InputInt("int2", &lscif::dbg_int2, 0, 0);
				ImGui::InputDouble("double", &lscif::dbg_dbl, 0, 0, "%.4f");
				ImGui::InputDouble("double2", &lscif::dbg_dbl2, 0, 0, "%.4f");
				ImGui::InputDouble("vector_scaling", &lscif::vector_scaling, 0, 0, "%.4f");
			}
			if (ImGui::Button("debug", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{

				int id = viewer.selected_data_index;
				CGMesh updatedMesh;
				CGMesh inputMesh = lscif::Meshes[id];
				lscif::tools.init(inputMesh);
				lscif::tools.initialize_mesh_properties();
				// calculate the pseudo-geodesic
				lscif::tools.debug_tool(lscif::dbg_int, lscif::dbg_int2, lscif::dbg_dbl);
				lscif::updateMeshViewer(viewer, inputMesh);
				lscif::meshFileName.push_back("dbg_" + lscif::meshFileName[id]);
				lscif::Meshes.push_back(inputMesh);
				Eigen::MatrixXd E0, E1;
				// lscif::tools.show_gradients(E0,E1, lscif::vector_scaling);
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);
				const Eigen::RowVector3d black(0, 0, 0);
				const Eigen::RowVector3d green(0.2, 0.8, 0.2);
				Eigen::MatrixXd RGB = Eigen::MatrixXd::Identity(3, 3);
				Eigen::MatrixXd E2, E3, Ea0, Ea1;
				Eigen::MatrixXd pts;
				std::vector<Eigen::MatrixXd> E0list, E1list;
				lscif::tools.show_pseudo_geodesic_curve(E0list, E1list, pts);
				E0 = E0list[0];
				E1 = E1list[0];
				viewer.data().add_edges(E0, E1, red);
				if (0)
				{
					lscif::tools.show_current_reference_points(pts);
					viewer.data().add_points(pts, red);
				}
				if (1)
				{
					viewer.data().add_points(pts, red);
				}
				viewer.data().add_points(lscif::tools.ver_dbg1, black);
				viewer.data().add_points(lscif::tools.ver_dbg, blue);
				for (int i = 0; i < lscif::tools.ver_dbg.rows(); i++)
				{
					Eigen::MatrixXd tmp_ver(1, 3);
					tmp_ver.row(0) = lscif::tools.ver_dbg.row(i);
					viewer.data().add_edges(lscif::tools.ver_dbg1, tmp_ver, blue);
				}
				std::cout << "checking query size " << lscif::tools.ver_dbg1.rows() << std::endl;
				std::cout << "checking result size " << lscif::tools.ver_dbg.rows() << std::endl;
				std::cout << "checking results " << lscif::tools.ver_dbg << std::endl;
				Eigen::MatrixXd Edge1_dbg = lscif::tools.E0_dbg + lscif::tools.direction_dbg * lscif::vector_scaling;
				viewer.data().add_edges(lscif::tools.E0_dbg, Edge1_dbg, green);
				if (lscif::tools.pnorm_dbg.rows() == lscif::tools.E0_dbg.rows())
				{
					Eigen::MatrixXd Edge2_dbg = lscif::tools.E0_dbg + lscif::tools.pnorm_dbg * lscif::vector_scaling;
					viewer.data().add_edges(lscif::tools.E0_dbg, Edge2_dbg, black);
				}

				viewer.selected_data_index = id;
			}
			ImGui::SameLine();
			if (ImGui::Button("debug cylinder", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{

				int id = viewer.selected_data_index;
				CGMesh updatedMesh;
				CGMesh inputMesh = lscif::Meshes[id];
				lscif::tools.init(inputMesh);
				lscif::tools.initialize_mesh_properties();
				// calculate the pseudo-geodesic
				std::vector<int> input_int;
				std::vector<double> input_dbl;
				input_int.push_back(lscif::dbg_int);
				input_int.push_back(lscif::dbg_int2);
				input_dbl.push_back(lscif::dbg_dbl);
				input_dbl.push_back(lscif::dbg_dbl2);
				lscif::tools.debug_tool_v2(input_int, input_dbl);
				lscif::updateMeshViewer(viewer, inputMesh);
				lscif::meshFileName.push_back("dbg_" + lscif::meshFileName[id]);
				lscif::Meshes.push_back(inputMesh);
				Eigen::MatrixXd E0, E1;
				// lscif::tools.show_gradients(E0,E1, lscif::vector_scaling);
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);
				const Eigen::RowVector3d black(0, 0, 0);
				const Eigen::RowVector3d green(0.2, 0.8, 0.2);
				Eigen::MatrixXd RGB = Eigen::MatrixXd::Identity(3, 3);
				Eigen::MatrixXd E2, E3, Ea0, Ea1;
				Eigen::MatrixXd pts;
				std::vector<Eigen::MatrixXd> E0list, E1list;
				lscif::tools.show_pseudo_geodesic_curve(E0list, E1list, pts);
				for (int i = 0; i < E0list.size(); i++)
				{
					E0 = E0list[i];
					E1 = E1list[i];
					viewer.data().add_edges(E0, E1, red);
					std::cout << "edge sizes " << E0.rows() << std::endl;
				}

				if (0)
				{
					lscif::tools.show_current_reference_points(pts);
					viewer.data().add_points(pts, red);
				}
				if (1)
				{
					viewer.data().add_points(pts, red);
				}

				// viewer.data().add_points(lscif::tools.ver_dbg1, black);
				// viewer.data().add_points(lscif::tools.ver_dbg, blue);
				// for(int i=0;i<lscif::tools.ver_dbg.rows();i++){
				// 	Eigen::MatrixXd tmp_ver(1,3);
				// 	tmp_ver.row(0)=lscif::tools.ver_dbg.row(i);
				// 	viewer.data().add_edges(lscif::tools.ver_dbg1,tmp_ver,blue);
				// }
				// std::cout << "checking query size " << lscif::tools.ver_dbg1.rows() << std::endl;
				// std::cout << "checking result size " << lscif::tools.ver_dbg.rows() << std::endl;
				// std::cout << "checking results " << lscif::tools.ver_dbg<<std::endl;
				// Eigen::MatrixXd Edge1_dbg=lscif::tools.E0_dbg+lscif::tools.direction_dbg*lscif::vector_scaling;
				// viewer.data().add_edges(lscif::tools.E0_dbg,Edge1_dbg,green);
				// Eigen::MatrixXd Edge2_dbg=lscif::tools.E0_dbg+lscif::tools.pnorm_dbg*lscif::vector_scaling;
				// viewer.data().add_edges(lscif::tools.E0_dbg,Edge2_dbg,black);
				viewer.selected_data_index = id;
			}
			ImGui::SameLine();
			if (ImGui::Button("debug neighbour search", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{

				int id = viewer.selected_data_index;
				CGMesh updatedMesh;
				CGMesh inputMesh = lscif::Meshes[id];
				lscif::tools.init(inputMesh);
				lscif::tools.initialize_mesh_properties();
				// calculate the pseudo-geodesic
				lscif::tools.debug_tool_v3(lscif::dbg_int, lscif::dbg_int2, lscif::dbg_dbl);
				lscif::updateMeshViewer(viewer, inputMesh);
				lscif::meshFileName.push_back("dbg_" + lscif::meshFileName[id]);
				lscif::Meshes.push_back(inputMesh);
				Eigen::MatrixXd E0, E1;
				// lscif::tools.show_gradients(E0,E1, lscif::vector_scaling);
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);
				const Eigen::RowVector3d black(0, 0, 0);
				const Eigen::RowVector3d green(0.2, 0.8, 0.2);
				Eigen::MatrixXd RGB = Eigen::MatrixXd::Identity(3, 3);
				Eigen::MatrixXd E2, E3, Ea0, Ea1;
				Eigen::MatrixXd pts;
				std::vector<Eigen::MatrixXd> E0list, E1list;
				lscif::tools.show_pseudo_geodesic_curve(E0list, E1list, pts);
				E0 = E0list[0];
				E1 = E1list[0];
				viewer.data().add_edges(E0, E1, red);
				if (1)
				{
					// lscif::tools.show_current_reference_points(pts);
					viewer.data().add_points(lscif::tools.visual_pts_dbg, black);
					std::cout << "the nbr of checking points " << lscif::tools.visual_pts_dbg.rows() << std::endl;
					std::cout << lscif::tools.visual_pts_dbg << std::endl;
				}
				if (1)
				{
					viewer.data().add_points(pts, red);
				}
				viewer.data().add_points(lscif::tools.ver_dbg1, black);
				viewer.data().add_points(lscif::tools.ver_dbg, blue);
				for (int i = 0; i < lscif::tools.ver_dbg.rows(); i++)
				{
					Eigen::MatrixXd tmp_ver(1, 3);
					tmp_ver.row(0) = lscif::tools.ver_dbg.row(i);
					viewer.data().add_edges(lscif::tools.ver_dbg1, tmp_ver, blue);
				}
				std::cout << "checking query size " << lscif::tools.ver_dbg1.rows() << std::endl;
				std::cout << "checking result size " << lscif::tools.ver_dbg.rows() << std::endl;
				std::cout << "checking results " << lscif::tools.ver_dbg << std::endl;
				Eigen::MatrixXd Edge1_dbg = lscif::tools.E0_dbg + lscif::tools.direction_dbg * lscif::vector_scaling;
				viewer.data().add_edges(lscif::tools.E0_dbg, Edge1_dbg, green);
				if (lscif::tools.pnorm_dbg.rows() == lscif::tools.E0_dbg.rows())
				{
					Eigen::MatrixXd Edge2_dbg = lscif::tools.E0_dbg + lscif::tools.pnorm_dbg * lscif::vector_scaling;
					viewer.data().add_edges(lscif::tools.E0_dbg, Edge2_dbg, black);
				}

				viewer.selected_data_index = id;
				// std::cout<<"acos(1) "<<acos(1.)<<std::endl;
			}
			
			if (ImGui::Button("levelset tracing", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{

				int id = viewer.selected_data_index;
				CGMesh updatedMesh;
				CGMesh inputMesh = lscif::Meshes[id];
				lscif::tools.init(inputMesh);
				lscif::tools.weight_mass=lscif::weight_mass;
				lscif::tools.weight_laplacian=lscif::weight_laplacian;
				lscif::tools.weight_boundary=lscif::weight_boundary;
				std::vector<int> input_int;
				std::vector<double> input_dbl;
				input_int.push_back(lscif::nbr_of_intervals);
				input_int.push_back(lscif::id_debug_global);
				input_int.push_back(lscif::which_seg_id);
				input_dbl.push_back(lscif::target_angle);
				input_dbl.push_back(lscif::start_angle);
				input_dbl.push_back(lscif::threadshold_angel_degree);
				lscif::tools.initialize_mesh_properties();
				lscif::tools.initialize_level_set_by_tracing(input_int, input_dbl);
				std::cout<<"finish tracing"<<std::endl;
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
				std::cout<<"before ploting the curve"<<std::endl;
				lscif::tools.show_pseudo_geodesic_curve(E0list, E1list, pts);
				for (int i = 0; i < E0list.size(); i++)// plot the curves
				{
					E0 = E0list[i];
					E1 = E1list[i];
					std::cout << "edge sizes " << E0.rows()<<", "<<E1.size() << std::endl;
					viewer.data().add_edges(E0, E1, red);
					
				}
				std::cout<<"the number of curves "<<E0list.size()<<std::endl;
				
				if (1)// plot the vertices of the curves
				{
					viewer.data().add_points(pts, red);
				}
				viewer.selected_data_index = id;
				// std::cout<<"acos(1) "<<acos(1.)<<std::endl;
			}
			ImGui::SameLine();
			if (ImGui::Button("smooth+boundary", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{

				int id = viewer.selected_data_index;
				CGMesh inputMesh = lscif::Meshes[id];
				lscif::tools.weight_mass = lscif::weight_mass;
				lscif::tools.weight_boundary = lscif::weight_boundary;
				lscif::tools.weight_laplacian = lscif::weight_laplacian;
				

				for (int i = 0; i < lscif::OpIter; i++)
				{
					lscif::tools.optimize_laplacian_with_traced_boundary_condition();
					std::cout << "step length " << lscif::tools.level_set_step_length << std::endl;
				}
				std::cout << "waiting for instruction..." << std::endl;
				// lscif::MP.MeshUnitScale(inputMesh, updatedMesh);
				lscif::updateMeshViewer(viewer, inputMesh);
				lscif::meshFileName.push_back("smt_" + lscif::meshFileName[id]);
				lscif::Meshes.push_back(inputMesh);

				Eigen::VectorXd level_set_values;
				lscif::tools.show_level_set(level_set_values);
				viewer.data().set_colors(level_set_values);
				// Eigen::MatrixXd E0, E1;
				// // lscif::tools.show_gradients(E0,E1, lscif::vector_scaling);
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);

				viewer.selected_data_index = id;
			}
			if (ImGui::Button("sphere example", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{

				int id = viewer.selected_data_index;
				CGMesh inputMesh = lscif::Meshes[id];
				lscif::tools.init(inputMesh);
				lscif::tools.initialize_mesh_properties();
				lscif::tools.make_sphere_ls_example(lscif::dbg_int);
				// lscif::MP.MeshUnitScale(inputMesh, updatedMesh);
				lscif::updateMeshViewer(viewer, inputMesh);
				lscif::meshFileName.push_back("dbg_" + lscif::meshFileName[id]);
				lscif::Meshes.push_back(inputMesh);

				Eigen::VectorXd level_set_values;
				lscif::tools.show_level_set(level_set_values);
				viewer.data().set_colors(level_set_values);
				Eigen::MatrixXd E0, E1;
				// lscif::tools.show_gradients(E0,E1, lscif::vector_scaling);
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);
				Eigen::MatrixXd RGB = Eigen::MatrixXd::Identity(3, 3);
				Eigen::MatrixXd E2, E3, Ea0, Ea1;
				Eigen::MatrixXd pts;
				lscif::tools.show_current_reference_points(pts);
				viewer.data().add_points(pts, red);
				// lscif::tools.show_gradients(E2, E3, lscif::vector_scaling);
				int which = 1;
				lscif::tools.show_hessian(E2, E3, lscif::vector_scaling, which);
				// Eigen::VectorXd gradient_xyz_values;
				// lscif::tools.show_gradient_scalar(gradient_xyz_values,which);
				// //  lscif::tools.show_1_order_derivate(E0, E1, E2, E3, lscif::vector_scaling);
				// //  lscif::tools.show_vertex_normal(E0,E1,lscif::vector_scaling);
				// //  viewer.data().add_edges(E0,E1,red);
				viewer.data().add_edges(E2, E3, blue);
				// lscif::show_axis(Ea0,Ea1,lscif::vector_scaling);
				// viewer.data().add_edges(Ea0,Ea1,RGB);
				// viewer.data().set_colors(gradient_xyz_values);
				viewer.selected_data_index = id;
			}
			if (ImGui::Button("init level set", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				// int id = viewer.selected_data_index;
				// std::cout<<"viewer.data().V.rows "<<viewer.data().V.rows()<<"\n";
				// std::cout<<"viewer.data_list[id].V.rows() "<<viewer.data_list[id].V.rows()<<"\n";
				// std::cout<<"If the above two numbers, then viewer.data() refers to the selected data."<<std::endl;
				int id = viewer.selected_data_index;
				CGMesh inputMesh = lscif::Meshes[id];
				lscif::tools.init(inputMesh);
				lscif::tools.initialize_mesh_properties();
				viewer.selected_data_index = id;
			}
			// if (ImGui::Button("smooth level set", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			// {

			// 	int id = viewer.selected_data_index;
			// 	CGMesh inputMesh = lscif::Meshes[id];
			// 	lscif::tools.weight_assign_face_value = lscif::weight_assign_face_value;
			// 	lscif::tools.weight_mass = lscif::weight_mass;
			// 	lscif::tools.assign_face_id = lscif::assign_face_id;
			// 	lscif::tools.assign_value[0] = lscif::assign_value0;
			// 	lscif::tools.assign_value[1] = lscif::assign_value1;
			// 	lscif::tools.assign_value[2] = lscif::assign_value2;

			// 	for (int i = 0; i < lscif::OpIter; i++)
			// 	{
			// 		lscif::tools.initialize_and_smooth_level_set_by_laplacian();
			// 		std::cout << "step length " << lscif::tools.level_set_step_length << std::endl;
			// 	}
			// 	std::cout << "waiting for instruction..." << std::endl;
			// 	// lscif::MP.MeshUnitScale(inputMesh, updatedMesh);
			// 	lscif::updateMeshViewer(viewer, inputMesh);
			// 	lscif::meshFileName.push_back("dbg_" + lscif::meshFileName[id]);
			// 	lscif::Meshes.push_back(inputMesh);

			// 	Eigen::VectorXd level_set_values;
			// 	lscif::tools.show_level_set(level_set_values);
			// 	viewer.data().set_colors(level_set_values);
			// 	// Eigen::MatrixXd E0, E1;
			// 	// // lscif::tools.show_gradients(E0,E1, lscif::vector_scaling);
			// 	const Eigen::RowVector3d red(0.8, 0.2, 0.2);
			// 	const Eigen::RowVector3d blue(0.2, 0.2, 0.8);
			// 	// Eigen::MatrixXd E2, E3;
			// 	// //lscif::tools.show_face_gradients(E2, E3, lscif::vector_scaling);
			// 	// // lscif::tools.show_1_order_derivate(E0, E1, E2, E3, lscif::vector_scaling);
			// 	// // lscif::tools.show_vertex_normal(E0,E1,lscif::vector_scaling);
			// 	// // viewer.data().add_edges(E0,E1,red);
			// 	// // viewer.data().add_edges(E2,E3,blue);
			// 	Eigen::MatrixXd pts;
			// 	lscif::tools.show_current_reference_points(pts);
			// 	viewer.data().add_points(pts, red);

			// 	viewer.selected_data_index = id;
			// }
			// ImGui::SameLine();
			// if (ImGui::Button("opti band width", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			// {

			// 	int id = viewer.selected_data_index;
			// 	CGMesh inputMesh = lscif::Meshes[id];
			// 	lscif::tools.weight_assign_face_value = lscif::weight_assign_face_value;
			// 	lscif::tools.weight_mass = lscif::weight_mass;
			// 	lscif::tools.assign_face_id = lscif::assign_face_id;
			// 	lscif::tools.assign_value[0] = lscif::assign_value0;
			// 	lscif::tools.assign_value[1] = lscif::assign_value1;
			// 	lscif::tools.assign_value[2] = lscif::assign_value2;
			// 	lscif::tools.strip_width = lscif::strip_width_ratio;
			// 	for (int i = 0; i < lscif::OpIter; i++)
			// 	{
			// 		lscif::tools.initialize_and_optimize_strip_width();
			// 		std::cout << "step length " << lscif::tools.level_set_step_length << std::endl;
			// 	}
			// 	std::cout << "waiting for instruction..." << std::endl;
			// 	// lscif::MP.MeshUnitScale(inputMesh, updatedMesh);
			// 	lscif::updateMeshViewer(viewer, inputMesh);
			// 	lscif::meshFileName.push_back("dbg_" + lscif::meshFileName[id]);
			// 	lscif::Meshes.push_back(inputMesh);

			// 	Eigen::VectorXd level_set_values;
			// 	lscif::tools.show_level_set(level_set_values);
			// 	viewer.data().set_colors(level_set_values);
			// 	// Eigen::MatrixXd E0, E1;
			// 	// // lscif::tools.show_gradients(E0,E1, lscif::vector_scaling);
			// 	const Eigen::RowVector3d red(0.8, 0.2, 0.2);
			// 	const Eigen::RowVector3d blue(0.2, 0.2, 0.8);
			// 	// Eigen::MatrixXd E2, E3;
			// 	// //lscif::tools.show_face_gradients(E2, E3, lscif::vector_scaling);
			// 	// // lscif::tools.show_1_order_derivate(E0, E1, E2, E3, lscif::vector_scaling);
			// 	// // lscif::tools.show_vertex_normal(E0,E1,lscif::vector_scaling);
			// 	// // viewer.data().add_edges(E0,E1,red);
			// 	// // viewer.data().add_edges(E2,E3,blue);
			// 	Eigen::MatrixXd pts;
			// 	lscif::tools.show_current_reference_points(pts);
			// 	viewer.data().add_points(pts, red);

			// 	viewer.selected_data_index = id;
			// }
			ImGui::SameLine();
			if (ImGui::Button("MeshUnitScale", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				int id = viewer.selected_data_index;
				CGMesh updatedMesh;
				CGMesh inputMesh = lscif::Meshes[id];

				lscif::MP.MeshUnitScale(inputMesh, updatedMesh);
				lscif::updateMeshViewer(viewer, updatedMesh);
				lscif::meshFileName.push_back("unit_" + lscif::meshFileName[id]);
				lscif::Meshes.push_back(updatedMesh);
				viewer.selected_data_index = id;
			}

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

				viewer.erase_mesh(id);
				if (id > -1)
				{
					lscif::meshFileName.erase(lscif::meshFileName.begin() + id);
					lscif::Meshes.erase(lscif::Meshes.begin() + id);
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
		}

		if (ImGui::CollapsingHeader("Closest point projection", ImGuiTreeNodeFlags_DefaultOpen))
		{
			ImGui::PushItemWidth(100);
			ImGui::InputDouble("Px", &lscif::InputPx, 0, 0, "%.6f");
			ImGui::SameLine();
			ImGui::InputDouble("Py", &lscif::InputPy, 0, 0, "%.6f");
			ImGui::SameLine();
			ImGui::InputDouble("Pz", &lscif::InputPz, 0, 0, "%.6f");
			ImGui::SameLine();
			if (ImGui::Button("Draw Point", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				int id = viewer.selected_data_index;

				CGMesh inputMesh = lscif::Meshes[id];

				OptimizerCheckboard Optimizer;
				Optimizer.referenceMesh = inputMesh;
				CGMesh::Point pt = CGMesh::Point(lscif::InputPx, lscif::InputPy, lscif::InputPz);
				Eigen::MatrixXd ptM(1, 3);
				ptM << pt[0], pt[1], pt[2];
				viewer.data().add_points(ptM, Eigen::RowVector3d(0, 0, 1)); // draw input point
			}
			ImGui::SameLine();
			if (ImGui::Button("Closest Point", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				int id = viewer.selected_data_index;

				CGMesh inputMesh = lscif::Meshes[id];

				OptimizerCheckboard Optimizer;
				Optimizer.referenceMesh = inputMesh;
				CGMesh::Point pt = CGMesh::Point(lscif::InputPx, lscif::InputPy, lscif::InputPz);
				CGMesh::Point Closestpt;
				double d = Optimizer.closest_point_projection(pt, Closestpt);
				std::cout << "Input point is (" << pt[0] << "," << pt[1] << "," << pt[2] << ")." << std::endl;
				std::cout << "Its closest point is (" << Closestpt[0] << "," << Closestpt[1] << "," << Closestpt[2] << ")." << std::endl;
				std::cout << "The closest distance is " << d << std::endl;

				Eigen::MatrixXd ptM(1, 3);
				ptM << pt[0], pt[1], pt[2];
				Eigen::MatrixXd ClosetptM(1, 3);
				ClosetptM << Closestpt[0], Closestpt[1], Closestpt[2];
				viewer.data().add_points(ptM, Eigen::RowVector3d(0, 0, 1));			  // draw input point
				viewer.data().add_points(ClosetptM, Eigen::RowVector3d(1, 0, 0));	  // draw its closest point
				viewer.data().add_edges(ptM, ClosetptM, Eigen::RowVector3d(1, 0, 0)); // draw connection line
			}
			ImGui::PopItemWidth();
		}

		ImGui::SetNextTreeNodeOpen(true);
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

	viewer.callback_key_down =
		[&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)
	{
		if (key == GLFW_KEY_BACKSPACE)
		{
			int old_id = viewer.data().id;
			if (viewer.erase_mesh(viewer.selected_data_index))
			{
				colors.erase(old_id);
				last_selected = -1;
			}
			return true;
		}
		return false;
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
	viewer.core().is_animating = false;
	viewer.core().animation_max_fps = 30.;
	viewer.core().background_color = Eigen::Vector4f(1, 1, 1, 100);
	viewer.launch();
}
