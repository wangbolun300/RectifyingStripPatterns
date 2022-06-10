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

// -------------------- OpenMesh
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
// ----------------------------------------------------------------------------
typedef OpenMesh::PolyMesh_ArrayKernelT<>  CGMesh;
typedef Eigen::SparseMatrix<double> spMat;


MeshProcessing MP;
const Eigen::RowVector3d sea_green(70. / 255., 252. / 255., 167. / 255.);
const Eigen::RowVector3d hot_red(255. / 255., 12. / 255., 17. / 255.);
Eigen::MatrixXd V, U;
Eigen::MatrixXi F, E;
Eigen::VectorXi  b;

bool keyPress_1 = false;
bool keyPress_2 = false;


//Optimization Parameters
double weigth_closeness = 0.0;
double weight_Planarity = 0.0;
double weight_EqualDiagornals = 1;
double weight_polyhedron = 0;
double weight_OrthDiagornals = 1;
double weight_UniformDiagornals = 0;
double weight_Fairness = 0;
double weight_FairnessDiag = 0;
double weight_polyhedron_vertex = 0;
double weight_Normalize_vnormal = 0;
double refer_AveEL = 1;
static int fairnessChangeType = 0;
std::vector<int> fixedVertices;

int OpIter = 10;
CGMesh mesh;
CGMesh RefMesh;
static bool fixCorner_checkbox = false;
static bool fixBoundary_checkbox = false;
static bool relativeFairness_checkbox = false;
static bool selectFEV[30] = { true, true, true, true, true, true, true, true, true, true, true, true, true, true, true,true, true, true, true, true, true, true, true, true, true,true, true, true, true, true };
std::vector<std::string> meshFileName;
std::vector<CGMesh> Meshes;

double InputPx = 0.0;
double InputPy = 0.0;
double InputPz = 0.0;
std::vector<int> VertexType;
void updateMeshViewer(igl::opengl::glfw::Viewer& viewer, CGMesh& mesh)
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

	/*viewer.data().compute_normals();
	viewer.data().uniform_colors(Eigen::Vector3d(51.0 / 255.0, 43.0 / 255.0, 33.3 / 255.0),
		Eigen::Vector3d(255.0 / 255.0, 228.0 / 255.0, 58.0 / 255.0),
		Eigen::Vector3d(255.0 / 255.0, 235.0 / 255.0, 80.0 / 255.0));*/
}
bool pre_draw(igl::opengl::glfw::Viewer& viewer)
{
	return false;
}
bool key_up(igl::opengl::glfw::Viewer& viewer, unsigned char key, int mods)
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
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int mods)
{
	switch (key)
	{

	case '1':
	{
		keyPress_1 = true;
		//std::cout << viewer.current_mouse_x << "___" << viewer.current_mouse_y << std::endl;
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
bool mouse_down(igl::opengl::glfw::Viewer& viewer, int button, int modifier)
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

			std::cout << viewer.current_mouse_x << "___" << viewer.current_mouse_y << std::endl;
			std::cout << vid << "_xxx_" << fid << std::endl;
			//viewer.data().add_points(igl::slice(viewer.data().V, vids, 1), hot_red);
			viewer.data().set_points(igl::slice(viewer.data().V, vids, 1), hot_red);
			viewer.data().add_points(igl::slice(viewer.data().V, vids2, 1), sea_green);
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


			std::cout << viewer.current_mouse_x << "___" << viewer.current_mouse_y << std::endl;
			std::cout << vid << "_xxx_" << fid << std::endl;		
			viewer.data().set_points(igl::slice(viewer.data().V, vids, 1), hot_red);
			viewer.data().add_points(igl::slice(viewer.data().V, vids2, 1), sea_green);
			return true;
		}
	}


	return false;
}

int main(int argc, char* argv[])
{
	std::map<int, Eigen::RowVector3d> colors;
	int last_selected = -1;
	std::string example_root_path(CHECKER_BOARD_DATA_DIR);
#ifdef _WIN32
			std::string spliter="\\";
#else
			std::string spliter="/";
#endif
	std::string fname = example_root_path+std::string("fertility_input.obj");
	OpenMesh::IO::read_mesh(mesh, fname);

	MP.mesh2Matrix(mesh, V, F);
	meshFileName.push_back("fertility");
	Meshes.push_back(mesh);

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
		ImGui::SetNextWindowSize(ImVec2(800, 800), ImGuiSetCond_FirstUseEver);
		ImGui::Begin(
			"Checkerboard Patterns", nullptr,
			ImGuiWindowFlags_NoSavedSettings
		);

		// Add a button
		if (ImGui::Button("Import Mesh", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
		{

			std::string fname = igl::file_dialog_open();
			if (fname.length() == 0)
				return;

			OpenMesh::IO::read_mesh(mesh, fname);

			size_t last_dot = fname.rfind('.');
			size_t last_slash = fname.rfind(spliter);// TODO on linux it should be '/'
			std::string fnameNew = fname.substr(last_slash + 1, (last_dot - last_slash - 1));
			meshFileName.push_back(fnameNew);
			Meshes.push_back(mesh);

			int nbe = 0;
			for (CGMesh::EdgeIter eit = mesh.edges_begin(); eit != mesh.edges_end(); ++eit)
			{
				if (mesh.is_boundary(eit))
					nbe++;
			}

			int nbhe = 0;
			for (CGMesh::HalfedgeIter heit = mesh.halfedges_begin(); heit != mesh.halfedges_end(); ++heit)
			{
				if (mesh.is_boundary(heit))
					nbhe++;
			}

			std::cout << "Mesh is: " << fname << std::endl;

			std::cout << "Vertices/Edges/Faces/HalfEdges/boundaryHalfEdge/boundaryEdge: " << mesh.n_vertices() << "/" << mesh.n_edges() << "/" << mesh.n_faces() << "/" << mesh.n_halfedges() << "/" << nbhe << "/" << nbe << std::endl;

			updateMeshViewer(viewer, mesh);
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
			OpenMesh::IO::write_mesh(Meshes[id], fname);
		}

		ImGui::SameLine();

		if (ImGui::Button("Import Para", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
		{
			std::string fname = igl::file_dialog_open();

			if (fname.length() == 0)
				return;

			std::ifstream settings(fname);
			settings >> weigth_closeness;
			settings >> weight_Planarity;
			settings >> weight_polyhedron_vertex;
			settings >> weight_Normalize_vnormal;
			settings >> weight_EqualDiagornals;
			settings >> weight_OrthDiagornals;
			settings >> weight_UniformDiagornals;
			settings >> weight_Fairness;
			settings >> weight_FairnessDiag;
			settings >> OpIter;
		}

		ImGui::SameLine();

		if (ImGui::Button("Save Para", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
		{
			std::string fname = igl::file_dialog_save();

			std::ofstream file;
			file.open(fname);

			file << weigth_closeness << std::endl;
			file << weight_Planarity << std::endl;
			file << weight_polyhedron_vertex << std::endl;
			file << weight_Normalize_vnormal << std::endl;
			file << weight_EqualDiagornals << std::endl;
			file << weight_OrthDiagornals << std::endl;
			file << weight_UniformDiagornals << std::endl;
			file << weight_Fairness << std::endl;
			file << weight_FairnessDiag << std::endl;
			file << OpIter << std::endl;

			file.close();
		}


		static bool updateMeshFlag = false;
		static int updateMeshId = 0;
		if (ImGui::Button("Update Reference Mesh", ImVec2(ImGui::GetWindowSize().x * 0.35f, 0.0f)))
		{
			updateMeshFlag = true;
			updateMeshId = viewer.selected_data_index;
			RefMesh = Meshes[updateMeshId];
			double refAveEL = 0; double totalL, minL, maxL;
			MP.getEdgeLengthStatistics(RefMesh, totalL, minL, maxL, refAveEL);
			bool printInfo = true;
			if (printInfo)
			{
				//	std::cout << "Mesh Informatin Vertices/Edges/Faces/HalfEdges:::::" << RefMesh.n_vertices() << "/" << RefMesh.n_edges() << "/" << RefMesh.n_faces() << "/" << RefMesh.n_halfedges() << std::endl;
				std::cout << "Edge lengths Total/minL/maxL/averageL:::::" << totalL << "/" << minL << "/" << maxL << "/" << refAveEL << std::endl;
			}
			refer_AveEL = refAveEL;
		}

		if (updateMeshFlag)
		{
			ImGui::SameLine();
			std::string str = meshFileName[updateMeshId];
			const char* meshnameChar = str.c_str();
			ImGui::Text(meshnameChar);
		}


		// Add new group
		if (ImGui::CollapsingHeader("Weights", ImGuiTreeNodeFlags_DefaultOpen))
		{
			// Expose variable directly ...			
			ImGui::InputDouble("Closeness", &weigth_closeness, 0, 0, "%.4f");
			ImGui::InputDouble("Planarity", &weight_Planarity, 0, 0, "%.4f");
			ImGui::InputDouble("Planarity_whiteQuad", &weight_polyhedron_vertex, 0, 0, "%.4f");
			ImGui::InputDouble("NormalizeVnormal", &weight_Normalize_vnormal, 0, 0, "%.4f");
			ImGui::InputDouble("EqualDiagonals", &weight_EqualDiagornals, 0, 0, "%.4f");
			ImGui::InputDouble("OrthDiagonals", &weight_OrthDiagornals, 0, 0, "%.4f");
			ImGui::InputDouble("UniformDiagonals", &weight_UniformDiagornals, 0, 0, "%.4f");
			ImGui::InputDouble("Fairness", &weight_Fairness, 0, 0, "%.4f");

			ImGui::InputInt("Iteration", &OpIter, 0, 0);

			ImGui::InputDouble("Average Edge length", &refer_AveEL, 0, 0, "%.4f");
			ImGui::Checkbox("Fix Corners", &fixCorner_checkbox);
			ImGui::SameLine();
			ImGui::Checkbox("Fix Boundary", &fixBoundary_checkbox);

			ImGui::SameLine();
			ImGui::Checkbox("related fairness", &relativeFairness_checkbox);



		}

		ImGui::Separator();
		ImGui::PushStyleColor(ImGuiCol_Button, (ImVec4)ImColor::HSV(0.0f, 0.6f, 0.6f));
		ImGui::PushStyleColor(ImGuiCol_ButtonHovered, (ImVec4)ImColor::HSV(0.0f, 0.7f, 0.7f));
		ImGui::PushStyleColor(ImGuiCol_ButtonActive, (ImVec4)ImColor::HSV(0.0f, 0.8f, 0.8f));
		if (ImGui::Button("Optimization"))
		{
				int id = viewer.selected_data_index;
				CGMesh inputMesh = Meshes[id];
				mesh = inputMesh;

				if (inputMesh.n_vertices() > 0)
				{
					OptimizerCheckboard Optimizer;
					Optimizer.OpIter = OpIter;
					Optimizer.weight_closeness = weigth_closeness;
					Optimizer.weight_EqualDiagornals = weight_EqualDiagornals;
					Optimizer.weight_polyhedron = weight_polyhedron;
					Optimizer.weight_OrthDiagornals = weight_OrthDiagornals;
					Optimizer.weight_UniformDiagornals = weight_UniformDiagornals;
					Optimizer.weight_Fairness = weight_Fairness;
					Optimizer.weight_FairnessDiag = weight_FairnessDiag;
					Optimizer.weight_polyhedron_vertex = weight_polyhedron_vertex;
					Optimizer.weight_Nornalization_vertexNormal = weight_Normalize_vnormal;
					Optimizer.fairnessChangeType = fairnessChangeType;
					Optimizer.relativeFairness = relativeFairness_checkbox;

					Optimizer.coarseMesh = mesh;
					if (RefMesh.n_vertices() == 0)
						RefMesh = mesh;

					Optimizer.referenceMesh = RefMesh;
					double refAveEL = 0; double totalL, minL, maxL;
					if (RefMesh.n_vertices() > 0)
						MP.getEdgeLengthStatistics(RefMesh, totalL, minL, maxL, refAveEL);

					Optimizer.refAveEL = refAveEL;

					MP.getEdgeLengthStatistics(mesh, totalL, minL, maxL, refAveEL);
					Optimizer.Lsqure = refer_AveEL * refer_AveEL * 2;

					if (fixCorner_checkbox)
					{
						fixedVertices.clear();
						fixedVertices = MP.meshCorners(mesh);
						std::vector<int> corners = MP.meshCorners(mesh);
						int nFix = corners.size();
						b.resize(nFix);
						b = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(corners.data(), corners.size());
					}

					if (fixBoundary_checkbox)
					{
						fixedVertices.clear();
						std::vector<int> bvs;
						for (CGMesh::VertexIter vit = mesh.vertices_begin(); vit != mesh.vertices_end(); ++vit)
						{
							if (mesh.is_boundary(vit))
								bvs.push_back(vit.handle().idx());
						}
						fixedVertices = bvs;
					}
					for (int i = 0; i < mesh.n_vertices(); i++)
						if (VertexType[i] == 1)
							fixedVertices.push_back(i);

					Optimizer.FixedVertices = fixedVertices;


					std::vector<double> sol = Optimizer.fillMatandSolve_CheckBoard();
					CGMesh updatedMesh;
					MP.MeshUpdate(mesh, updatedMesh, sol);


					meshFileName.push_back("OP_" + meshFileName[id]);
					Meshes.push_back(updatedMesh);

					updateMeshViewer(viewer, updatedMesh);
					viewer.selected_data_index = id;
				}
		}
		ImGui::PopStyleColor(3);


		if (ImGui::CollapsingHeader("Mesh Processing", ImGuiTreeNodeFlags_DefaultOpen))
		{
			
			if (ImGui::Button("Black Rectangles", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				int id = viewer.selected_data_index;
				CGMesh updatedMesh;
				CGMesh inputMesh = Meshes[id];

				std::vector<int> Vtype = MP.MeshVertexClassificationBlackWhite(inputMesh);
				MP.MeshEdgeDual(inputMesh, updatedMesh, Vtype, 0.5, 0);
				
				updateMeshViewer(viewer, updatedMesh);
				viewer.data().face_based = true;
				viewer.data().uniform_colors(Eigen::Vector3d(51.0 / 255.0, 43.0 / 255.0, 33.3 / 255.0),
					Eigen::Vector3d(51.0 / 255.0, 43.0 / 255.0, 33.3 / 255.0),
					Eigen::Vector3d(51.0 / 255.0, 43.0 / 255.0, 33.3 / 255.0));
				meshFileName.push_back("BLACK_" + meshFileName[id]);
				Meshes.push_back(updatedMesh);

				// white faces
				CGMesh updatedMesh1;
				MP.MeshEdgeDual(inputMesh, updatedMesh1, Vtype, 0.5, 1);

				updateMeshViewer(viewer, updatedMesh1);
				viewer.data().face_based = true;
				viewer.data().uniform_colors(Eigen::Vector3d(151.0 / 255.0, 151.0 / 255.0, 151.3 / 255.0),
					Eigen::Vector3d(151.0 / 255.0, 151.0 / 255.0, 151.3 / 255.0),
					Eigen::Vector3d(151.0 / 255.0, 151.0 / 255.0, 151.3 / 255.0));
				meshFileName.push_back("white_" + meshFileName[id]);
				Meshes.push_back(updatedMesh1);
				viewer.selected_data_index = id;
			}
			ImGui::SameLine();
			if (ImGui::Button("CaltmullClack", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				int id = viewer.selected_data_index;
				CGMesh updatedMesh;
				CGMesh inputMesh = Meshes[id];


				MP.SubdivCoarseMesh_CaltmullClack(inputMesh, updatedMesh);

				updateMeshViewer(viewer, updatedMesh);
				meshFileName.push_back("Subdiv_" + meshFileName[id]);
				Meshes.push_back(updatedMesh);
				viewer.selected_data_index = id;

			}
			ImGui::SameLine();

			if (ImGui::Button("MeshUnitScale", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				int id = viewer.selected_data_index;
				CGMesh updatedMesh;
				CGMesh inputMesh = Meshes[id];

				MP.MeshUnitScale(inputMesh, updatedMesh);
				updateMeshViewer(viewer, updatedMesh);
				meshFileName.push_back("unit_" + meshFileName[id]);
				Meshes.push_back(updatedMesh);
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
				//if (id > 1)
				//	viewer.selected_data_index = id - 2;
				//std::cout << "The current Mesh ID is" << viewer.selected_data_index << std::endl;

				viewer.erase_mesh(id);
				if (id > -1)
				{
					meshFileName.erase(meshFileName.begin() + id);
					Meshes.erase(Meshes.begin() + id);
				}


			}
			ImGui::SameLine();

			if (ImGui::Button("Show boundary Ver", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				int id = viewer.selected_data_index;

				CGMesh inputMesh = Meshes[id];

				std::vector<int> Bvids;
				for (CGMesh::VertexIter v_it = inputMesh.vertices_begin(); v_it != (inputMesh.vertices_end()); ++v_it)
				{

					if (inputMesh.is_boundary(v_it.handle()))
						Bvids.push_back(v_it.handle().idx());

				}
				Eigen::VectorXi vids = Eigen::Map<Eigen::VectorXi, Eigen::Unaligned>(Bvids.data(), Bvids.size());

				viewer.data().set_points(igl::slice(viewer.data().V, vids, 1), hot_red);

			}

		}

		if (ImGui::CollapsingHeader("Closest point projection", ImGuiTreeNodeFlags_DefaultOpen))
		{
			ImGui::PushItemWidth(100);
			ImGui::InputDouble("Px", &InputPx, 0, 0, "%.6f");
			ImGui::SameLine();
			ImGui::InputDouble("Py", &InputPy, 0, 0, "%.6f");
			ImGui::SameLine();
			ImGui::InputDouble("Pz", &InputPz, 0, 0, "%.6f");
			ImGui::SameLine();
			if (ImGui::Button("Draw Point", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				int id = viewer.selected_data_index;

				CGMesh inputMesh = Meshes[id];

				OptimizerCheckboard Optimizer;
				Optimizer.referenceMesh = inputMesh;
				CGMesh::Point pt = CGMesh::Point(InputPx, InputPy, InputPz);
				Eigen::MatrixXd ptM(1, 3);
				ptM << pt[0], pt[1], pt[2];			
				viewer.data().add_points(ptM, Eigen::RowVector3d(0, 0, 1)); //draw input point				
			}
			ImGui::SameLine();
			if (ImGui::Button("Closest Point", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				int id = viewer.selected_data_index;

				CGMesh inputMesh = Meshes[id];

				OptimizerCheckboard Optimizer;
				Optimizer.referenceMesh = inputMesh;
				CGMesh::Point pt = CGMesh::Point(InputPx, InputPy, InputPz);
				CGMesh::Point Closestpt;
				double d=Optimizer.closest_point_projection(pt, Closestpt);
				std::cout << "Input point is (" << pt[0] << "," << pt[1] << "," << pt[2] << ")." << std::endl;
				std::cout<<"Its closest point is (" << Closestpt[0] << "," << Closestpt[1] << "," << Closestpt[2] << ")." << std::endl;
				std::cout << "The closest distance is " << d << std::endl;

				Eigen::MatrixXd ptM(1, 3);
				ptM << pt[0], pt[1], pt[2];
				Eigen::MatrixXd ClosetptM(1, 3);
				ClosetptM << Closestpt[0], Closestpt[1], Closestpt[2];
				viewer.data().add_points(ptM, Eigen::RowVector3d(0, 0, 1)); //draw input point
				viewer.data().add_points(ClosetptM, Eigen::RowVector3d(1, 0, 0)); // draw its closest point
				viewer.data().add_edges(ptM, ClosetptM, Eigen::RowVector3d(1, 0, 0)); // draw connection line

			}
			ImGui::PopItemWidth();
		}

		ImGui::SetNextTreeNodeOpen(true);
		if (ImGui::TreeNode("Mesh Management"))
		{


			//HelpMarker("This is a more typical looking tree with selectable nodes.\nClick to select, CTRL+Click to toggle, click on arrows or double-click to open.");


			static int selection_mask = (1 << 2); // Dumb representation of what may be user-side selection state. You may carry selection state inside or outside your objects in whatever format you see fit.
			int node_clicked = -1;                // Temporary storage of what node we have clicked to process selection at the end of the loop. May be a pointer to your own node type, etc.
			//ImGui::PushStyleVar(ImGuiStyleVar_IndentSpacing, ImGui::GetFontSize() * 3); // Increase spacing to differentiate leaves from expanded contents.
			int id = viewer.data_list.size();

			for (int i = 0; i < id; i++)
			{
				//// Disable the default open on single-click behavior and pass in Selected flag according to our selection state.
				ImGuiTreeNodeFlags node_flags = ImGuiTreeNodeFlags_OpenOnArrow | ImGuiTreeNodeFlags_OpenOnDoubleClick;
				if (selection_mask & (1 << i))
					node_flags |= ImGuiTreeNodeFlags_Selected;
				//static bool check = true;
				//ImGui::Checkbox("checkbox", &check);
				//ImGui::SameLine();
				//ImGui::Checkbox("F", &selectFEV[3 * i]);
				//ImGui::SameLine();
				//ImGui::Checkbox("E", &selectFEV[3 * i + 1]);
				//ImGui::SameLine();
				//ImGui::Checkbox("V", &selectFEV[3 * i + 2]);
				//ImGui::SameLine();
				ImGui::PushID(i);
				ImGui::Checkbox("Visible", &selectFEV[i]);
				ImGui::PopID();
				ImGui::SameLine();
				if (selectFEV[i])
					viewer.data_list[i].is_visible = true;
				else
					viewer.data_list[i].is_visible = false;





				node_flags |= ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen; // ImGuiTreeNodeFlags_Bullet
				//node_flags |=  ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_DefaultOpen; // ImGuiTreeNodeFlags_Bullet
				//ImGui::TreeNodeEx((void*)(intptr_t)i, node_flags, "Object %d", i);
				ImGui::TreeNodeEx((void*)(intptr_t)i, node_flags, meshFileName[i].c_str());

				if (ImGui::IsItemClicked())
					node_clicked = i;
			}
			if (node_clicked != -1)
			{
				// Update selection state. Process outside of tree loop to avoid visual inconsistencies during the clicking-frame.
				if (ImGui::GetIO().KeyCtrl)
					selection_mask ^= (1 << node_clicked);          // CTRL+click to toggle
				else //if (!(selection_mask & (1 << node_clicked))) // Depending on selection behavior you want, this commented bit preserve selection when clicking on item that is part of the selection
					selection_mask = (1 << node_clicked);           // Click to single-select

				viewer.selected_data_index = node_clicked;
			}

			ImGui::TreePop();
		}

		ImGui::End();
	};

	viewer.callback_key_down =
		[&](igl::opengl::glfw::Viewer&, unsigned int key, int mod)
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
		[&](igl::opengl::glfw::Viewer&)
	{
		if (last_selected != viewer.selected_data_index)
		{
			for (auto& data : viewer.data_list)
			{
				data.set_colors(colors[data.id]);
			}
			viewer.data_list[viewer.selected_data_index].set_colors(Eigen::RowVector3d(0.9, 0.9, 0.1));
			last_selected = viewer.selected_data_index;
		}
		return false;
	};




	// Plot the mesh
	viewer.data().set_mesh(V, F);
	viewer.data().show_lines = false;
	viewer.data().is_visible = true;
	// show polygon edges
	Eigen::MatrixXi E;
	MP.meshEdges(mesh, E);
	Eigen::MatrixXd C = Eigen::MatrixXd::Zero(E.rows(), 3);
	viewer.data().set_edges(V, E, C);
	viewer.data().line_width = 2.0;
	viewer.data().invert_normals = true;
	VertexType.resize(V.size(), 0);



	viewer.callback_pre_draw = &pre_draw;
	viewer.callback_key_down = &key_down;
	viewer.callback_mouse_down = &mouse_down;
	viewer.callback_key_up = &key_up;
	viewer.core().is_animating = false;
	viewer.core().animation_max_fps = 30.;
	viewer.core().background_color = Eigen::Vector4f(1, 1, 1, 100);
	viewer.launch();
}
