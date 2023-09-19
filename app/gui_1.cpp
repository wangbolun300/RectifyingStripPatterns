#include <gui.h>
void lscif::draw_menu1(igl::opengl::glfw::Viewer &viewer, igl::opengl::glfw::imgui::ImGuiPlugin &plugin,
					   igl::opengl::glfw::imgui::ImGuiMenu &menu)
{
	// Draw additional windows
	menu.callback_draw_custom_window = [&]()
	{
		// Define next window position + size
		ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(300, 1000), ImGuiCond_FirstUseEver);
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

			OpenMesh::IO::read_mesh(mesh, fname);
			tools.init(mesh);
			size_t last_dot = fname.rfind('.');
			size_t last_slash = fname.rfind(spliter); // TODO on linux it should be '/'
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
			Eigen::Vector3d vmin(tools.V.col(0).minCoeff(), tools.V.col(1).minCoeff(), tools.V.col(2).minCoeff());
			Eigen::Vector3d vmax(tools.V.col(0).maxCoeff(), tools.V.col(1).maxCoeff(), tools.V.col(2).maxCoeff());

			std::cout << "Vertices/Edges/Faces/HalfEdges/boundaryHalfEdge/boundaryEdge/BBD: "
					  << mesh.n_vertices() << "/" << mesh.n_edges() << "/" << mesh.n_faces() << "/" << mesh.n_halfedges() << "/" << nbhe << "/" << nbe << "/" << (vmin - vmax).norm() << std::endl;

			updateMeshViewer(viewer, mesh);
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
			OpenMesh::IO::write_mesh(Meshes[id], fname);
		}

		ImGui::SameLine();

		if (ImGui::Button("Import levelset", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
		{

			bool read = read_levelset(tools.fvalues);
			if (read)
			{
				std::cout << "\nlevel set readed" << std::endl;
			}
		}

		ImGui::SameLine();

		if (ImGui::Button("Save levelset", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
		{
			bool saved = save_levelset(tools.fvalues, tools.Binormals);
			if (saved)
			{
				std::cout << "\nlevel set saved" << std::endl;
			}
		}
		if (ImGui::Button("Read Binormals", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
		{
			std::cout << "Select a Bi-normal file for visulization" << std::endl;
			read_bi_normals(tools.Binormals);
			std::cout << "bi-normals got readed" << std::endl;
		}
		ImGui::SameLine();
		if (ImGui::Button("Debug", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
		{

			//  Eigen::Vector3d light =  get_light_rotated_back_from_earth_axis(Shading_Latitude, InputPx, InputPy);
			//  std::cout<<"The light direction rotated back is\n "<<light<<std::endl;
			// std::cout<<"checking the angle between the light and the target light"<<std::endl;
			// Eigen::Vector3d target_light = angle_ray_converter(InputPx, InputPy);
			// std::cout<<"target light, "<<target_light.transpose()<<std::endl;
			// auto light = tools.Lights;
			// for(int i = 0;i<light.rows();i++)
			// {
			// 	std::cout<<light.row(i).dot(target_light)<<std::endl;
			// }
		}
		ImGui::SameLine();
		if (ImGui::Button("UseUnitScaleMesh", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
		{
			std::cout << "Mesh Unit Scale According to Edge Length" << std::endl;
			int id = viewer.selected_data_index;
			CGMesh updatedMesh;
			CGMesh inputMesh = Meshes[id];
			MP.MeshUnitScale(inputMesh, updatedMesh);
			updateMeshViewer(viewer, updatedMesh);
			meshFileName.push_back("Uni_" + meshFileName[id]);
			Meshes.push_back(updatedMesh);
			viewer.selected_data_index = id;
			mesh = updatedMesh;
			tools.init(updatedMesh);
		}

		// Add new group
		if (ImGui::CollapsingHeader("Pseudo-Geodesic Angle", ImGuiTreeNodeFlags_DefaultOpen))
		{
			ImGui::PushItemWidth(50);
			ImGui::InputDouble("target angle", &target_angle, 0, 0, "%.4f");
			ImGui::SameLine();
			ImGui::InputDouble("second angle", &target_angle_2, 0, 0, "%.4f");
			ImGui::InputInt("Which Levelset", &which_levelset, 0, 0);
			ImGui::InputDouble("weight_geodesic", &weight_geodesic, 0, 0, "%.4f");
		}

		if (ImGui::CollapsingHeader("Energy Solving", ImGuiTreeNodeFlags_CollapsingHeader))
		{
			// Expose variable directly ...
			// ImGui::InputDouble("Closeness", &weigth_closeness, 0, 0, "%.4f");
			ImGui::InputInt("Iteration", &OpIter, 0, 0);
			ImGui::InputDouble("weight ls mass(big)", &weight_mass, 0, 0, "%.4f");
			ImGui::InputDouble("weight boundary (big)", &weight_boundary, 0, 0, "%.4f");
			ImGui::InputDouble("weight laplacian (small)", &weight_laplacian, 0, 0, "%.4f");
			ImGui::InputDouble("weight pseudo-geodesic", &weight_pseudo_geodesic, 0, 0, "%.4f");
			ImGui::InputDouble("weight Strip width(Need SW Enabled)", &weight_strip_width, 0, 0, "%.4f");
			ImGui::InputDouble("MaxStpLenght(Need PG Enabled)", &maximal_step_length, 0, 0, "%.4f");
			ImGui::InputDouble("WeightAngle", &weight_angle, 0, 0, "%.4f");
			ImGui::InputDouble("RatioEndPts", &weight_endpoint_ratio, 0, 0, "%.4f");

			ImGui::Checkbox("Enable Extreme Cases", &enable_extreme_cases);

			ImGui::Checkbox("Enable PG Energy", &enable_pg_energy_checkbox);
			ImGui::SameLine();
			ImGui::Checkbox("Enable Strip Width", &enable_strip_width_checkbox);
		}

		// ImGui::Checkbox("Fix Boundary", &fixBoundary_checkbox);

		if (ImGui::CollapsingHeader("Tracing", ImGuiTreeNodeFlags_DefaultOpen))
		{
			ImGui::Combo("TracingType", &TracingType,
						 "SegTrace\0OneEdgeTrace\0SegDrctn\0\0"),

				ImGui::InputInt("Select 1 boundary", &which_seg_id, 0, 0);
			ImGui::InputInt("every i segments", &nbr_of_intervals, 0, 0);
			ImGui::InputInt("start bnd edge ", &start_bnd_he, 0, 0);
			ImGui::InputInt("nbr bnd edges ", &nbr_edges, 0, 0);
			ImGui::InputDouble("start angle", &start_angle, 0, 0, "%.4f");
			ImGui::InputDouble("corner angle", &threadshold_angel_degree, 0, 0, "%.4f");
			if (ImGui::Button("Trace Curves", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				int id = viewer.selected_data_index;
				CGMesh updatedMesh;
				CGMesh inputMesh = Meshes[id];
				TracingPrepare Tracing_initializer;
				Tracing_initializer.every_n_edges = nbr_of_intervals;
				Tracing_initializer.start_angle = start_angle;
				Tracing_initializer.target_angle = target_angle;
				Tracing_initializer.threadshold_angel_degree = threadshold_angel_degree;
				Tracing_initializer.which_boundary_segment = which_seg_id;
				Tracing_initializer.start_bnd_he = start_bnd_he;
				Tracing_initializer.nbr_edges = nbr_edges;
				if (TracingType == 0)
				{
					tools.initialize_level_set_by_select_boundary_segment(Tracing_initializer, true);
				}
				if (TracingType == 1)
				{
					tools.initialize_level_set_by_tracing(Tracing_initializer);
				}
				if (TracingType == 2)
				{
					tools.initialize_level_set_by_select_boundary_segment(Tracing_initializer, false);
				}

				std::cout << "finish tracing" << std::endl;
				updateMeshViewer(viewer, inputMesh);
				meshFileName.push_back("dbg_" + meshFileName[id]);
				Meshes.push_back(inputMesh);
				Eigen::MatrixXd E0, E1;
				Eigen::MatrixXd E2, E3, Ea0, Ea1;
				// tools.show_gradients(E0,E1, vector_scaling);
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);
				const Eigen::RowVector3d black(0, 0, 0);
				const Eigen::RowVector3d green(0.2, 0.8, 0.2);
				Eigen::MatrixXd RGB = Eigen::MatrixXd::Identity(3, 3);

				Eigen::MatrixXd pts;
				std::vector<Eigen::MatrixXd> E0list, E1list;
				std::cout << "before ploting the curve" << std::endl;
				tools.show_pseudo_geodesic_curve(E0list, E1list, pts);
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
				tools.show_traced_binormals(bE0, bE1, nE0, nE1, vector_scaling);

				viewer.data().add_edges(bE0, bE1, green);
				viewer.data().add_edges(nE0, nE1, blue);
			}
			ImGui::SameLine();
			if (ImGui::Button("SaveTracedCurves", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				int id = viewer.selected_data_index;
				CGMesh updatedMesh;
				CGMesh inputMesh = Meshes[id];
				TracingPrepare Tracing_initializer;
				Tracing_initializer.every_n_edges = nbr_of_intervals;
				Tracing_initializer.start_angle = start_angle;
				Tracing_initializer.target_angle = target_angle;
				Tracing_initializer.threadshold_angel_degree = threadshold_angel_degree;
				Tracing_initializer.which_boundary_segment = which_seg_id;
				Tracing_initializer.start_bnd_he = start_bnd_he;
				Tracing_initializer.nbr_edges = nbr_edges;
				if (TracingType == 0)
				{
					tools.initialize_level_set_by_select_boundary_segment(Tracing_initializer, true);
				}
				if (TracingType == 1)
				{
					tools.initialize_level_set_by_tracing(Tracing_initializer);
				}
				if (TracingType == 2)
				{
					tools.initialize_level_set_by_select_boundary_segment(Tracing_initializer, false);
				}

				std::cout << "finish tracing" << std::endl;
				updateMeshViewer(viewer, inputMesh);
				meshFileName.push_back("dbg_" + meshFileName[id]);
				Meshes.push_back(inputMesh);
				Eigen::MatrixXd E0, E1;
				Eigen::MatrixXd E2, E3, Ea0, Ea1;
				// tools.show_gradients(E0,E1, vector_scaling);
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);
				const Eigen::RowVector3d black(0, 0, 0);
				const Eigen::RowVector3d green(0.2, 0.8, 0.2);
				Eigen::MatrixXd RGB = Eigen::MatrixXd::Identity(3, 3);

				Eigen::MatrixXd pts;
				std::vector<Eigen::MatrixXd> E0list, E1list;
				std::cout << "before ploting the curve" << std::endl;
				tools.show_pseudo_geodesic_curve(E0list, E1list, pts);
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
				std::cout << "SAVING The Traced Vertices, please provide the prefix..." << std::endl;
				std::string fname = igl::file_dialog_save();
				// std::ofstream file;
				// file.open(fname);
				tools.save_traced_curves(fname);

				// file.close();

				std::cout << "Traced Data SAVED" << std::endl;

				// std::cout<<"acos(1) "<<acos(1.)<<std::endl;
			}

			if (ImGui::Button("LoadTracedCurves", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				tools.load_one_traced_curve();
			}
			ImGui::SameLine();
			if (ImGui::Button("ClearTracedCurves", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				tools.clear_traced_curves();
			}
			ImGui::SameLine();
			if (ImGui::Button("SaveOBJCurve", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				int id = viewer.selected_data_index;
				CGMesh updatedMesh;
				CGMesh inputMesh = Meshes[id];
				TracingPrepare Tracing_initializer;
				Tracing_initializer.every_n_edges = nbr_of_intervals;
				Tracing_initializer.start_angle = start_angle;
				Tracing_initializer.target_angle = target_angle;
				Tracing_initializer.threadshold_angel_degree = threadshold_angel_degree;
				Tracing_initializer.which_boundary_segment = which_seg_id;
				Tracing_initializer.start_bnd_he = start_bnd_he;
				Tracing_initializer.nbr_edges = nbr_edges;
				if (TracingType == 0)
				{
					tools.initialize_level_set_by_select_boundary_segment(Tracing_initializer, true);
				}
				if (TracingType == 1)
				{
					tools.initialize_level_set_by_tracing(Tracing_initializer);
				}
				if (TracingType == 2)
				{
					tools.initialize_level_set_by_select_boundary_segment(Tracing_initializer, false);
				}

				std::cout << "finish tracing" << std::endl;
				updateMeshViewer(viewer, inputMesh);
				meshFileName.push_back("dbg_" + meshFileName[id]);
				Meshes.push_back(inputMesh);
				Eigen::MatrixXd E0, E1;
				Eigen::MatrixXd E2, E3, Ea0, Ea1;
				// tools.show_gradients(E0,E1, vector_scaling);
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);
				const Eigen::RowVector3d black(0, 0, 0);
				const Eigen::RowVector3d green(0.2, 0.8, 0.2);
				Eigen::MatrixXd RGB = Eigen::MatrixXd::Identity(3, 3);

				Eigen::MatrixXd pts;
				std::vector<Eigen::MatrixXd> E0list, E1list;
				std::cout << "before ploting the curve" << std::endl;
				tools.show_pseudo_geodesic_curve(E0list, E1list, pts);
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
				std::cout << "SAVING The Traced Vertices, Binormals and Normals into obj format, please provide prefix" << std::endl;
				std::string fname = igl::file_dialog_save();
				std::ofstream file;
				file.open(fname + ".obj");
				for (int i = 0; i < pts.rows(); i++)
				{
					file << "v " << pts(i, 0) << " " << pts(i, 1) << " " << pts(i, 2) << std::endl;
				}

				file.close();

				Eigen::MatrixXd bE0, bE1, nE0, nE1;
				tools.show_traced_binormals(bE0, bE1, nE0, nE1, vector_scaling);

				file.open(fname + "_start.obj");
				assert(bE0.rows() == pts.rows() - 2 && "Please Trace only one curve so we can save the correct data");
				for (int i = 0; i < bE0.rows(); i++)
				{
					file << "v " << pts(i + 1, 0) << " " << pts(i + 1, 1) << " " << pts(i + 1, 2) << std::endl;
				}

				file.close();

				file.open(fname + "_binormal_end.obj");
				for (int i = 0; i < bE1.rows(); i++)
				{
					file << "v " << bE1(i, 0) << " " << bE1(i, 1) << " " << bE1(i, 2) << std::endl;
				}

				file.close();

				file.open(fname + "_normal_end.obj");
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
			{ // this is to save the parametrized tracing results.

				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);

				const Eigen::RowVector3d black(0, 0, 0);
				const Eigen::RowVector3d green(0.2, 0.8, 0.2);

				Eigen::MatrixXd pts;
				std::cout << "Load one curve, convert it into 2d, then save it as an obj file" << std::endl;

				tools.show_traced_curve_params(pts);
				if (pts.rows() == 0)
				{
					std::cout << "Please follow the instructions before using this button" << std::endl;
					ImGui::End();
					return;
				}
				viewer.data().add_points(pts, red);
				std::cout << "Saving the curve points in parametrization domain..." << std::endl;
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
			// ImGui::Checkbox("Fix Boundary", &fixBoundary_checkbox);
		}
		ImGui::Separator();

		if (ImGui::CollapsingHeader("Debug input", ImGuiTreeNodeFlags_CollapsingHeader))
		{
			ImGui::Checkbox("debugFlag", &debug_flag);
			ImGui::InputInt("int", &dbg_int, 0, 0);
			ImGui::InputInt("int2", &dbg_int2, 0, 0);
			ImGui::InputDouble("double", &dbg_dbl, 0, 0, "%.4f");
			ImGui::InputDouble("double2", &dbg_dbl2, 0, 0, "%.4f");
			ImGui::InputDouble("vector_scaling", &vector_scaling, 0, 0, "%.4f");
			ImGui::SameLine();
			ImGui::InputDouble("vector_scaling2", &vector_scaling2, 0, 0, "%.4f");
			ImGui::InputInt("extracted_nbr", &extracted_nbr, 0, 0);
			ImGui::PushItemWidth(50);
			ImGui::InputDouble("x", &ptInx, 0, 0, "%.4f");
			ImGui::SameLine();
			ImGui::InputDouble("y", &ptIny, 0, 0, "%.4f");
			ImGui::SameLine();
			ImGui::InputDouble("z", &ptInz, 0, 0, "%.4f");
			ImGui::SameLine();
			if (ImGui::Button("drawPt", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				Eigen::MatrixXd pts(1, 3);
				pts.row(0) << ptInx, ptIny, ptInz;
				viewer.data().add_points(pts, hot_red);
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
				// Eigen::MatrixXd Vcout;
				// Eigen::MatrixXd Vlout;
				// match_the_two_catenaries(tools.lsmesh, tools.Boundary_Edges, tools.V,
				// 						 tools.F, tools.norm_v, tools.fvalues, Vcout, Vlout);
				// viewer.data().add_points(Vcout, hot_red);
				// viewer.data().add_points(Vlout, sea_green);
				// Eigen::MatrixXd Vout;
				// poly_tool.draw_inflections(Vout);
				// viewer.data().add_points(Vout, hot_red);
				std::string fname = igl::file_dialog_open();
				if (fname.length() == 0)
				{
					std::cout << "\nLSC: read mesh failed" << std::endl;
					ImGui::End();
					return;
				}
				Eigen::MatrixXd Vtri;
				Eigen::MatrixXi Ftri;
				igl::readOBJ(fname, Vtri, Ftri);
				double length = 0;
				for (int i = 0; i < Vtri.rows() - 1; i++)
				{
					length += Eigen::Vector3d(Vtri.row(i) - Vtri.row(i + 1)).norm();
				}
				std::cout << "length of the curve: " << length << std::endl;
			}
			ImGui::SameLine();
			if (ImGui::Button("DBG", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				// read_pts_csv_and_write_xyz_files();
				// recover_polyline_endpts();
				// using intersections of the rectifying planes to get the developable.
				// construct_single_developable_strips_by_intersect_rectifying(0);
				// write_unfold_single_strip();
				// draw_catenaries_on_cylinder();
				// tools.debug_tool();
				make_example_comparing_two_plylines_distance();
			}
			ImGui::SameLine();
			if (ImGui::Button("OBJ2CSV", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				obj2csv();
			}
			if (ImGui::Button("DrawQdOptDbg", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{

				viewer.data().add_points(quad_tool.Debugtool, hot_red);
			}
			ImGui::SameLine();
			if (ImGui::Button("ResortPlys", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				run_sort_polylines();
			}
			ImGui::SameLine();
			if (ImGui::Button("SaveSlope", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				Eigen::MatrixXd E0, E1;
				tools.show_slopes(vector_scaling, E0, E1);
				Eigen::MatrixXd rotmids, rotdirs;
				// viewer.data().add_edges(E0, E1, sea_green);
				show_rotback_slopes(E0, E1, Shading_Latitude,
									rotmids, rotdirs);
				viewer.data().add_edges(
					rotmids - vector_scaling * rotdirs, rotmids + vector_scaling * rotdirs, hot_red);
				viewer.data().add_points(
					rotmids, hot_red);
				save_rotback_slopes(rotmids, rotdirs);
			}

			ImGui::SameLine();
			if (ImGui::Button("HalfNbrPly", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				decrease_ply_nbr_by_half();
			}

			if (ImGui::Button("WriteTime", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				std::cout << "Total Time: " << time_total << std::endl;
				std::cout << "Total Iterations: " << iteration_total << std::endl;
				std::cout << "AVG Time: " << time_total / iteration_total << std::endl;
			}
			ImGui::SameLine();
			if (ImGui::Button("TimerErrorReset", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				time_total = 0;
				iteration_total = 0;
				timer_global.stop();
				timelist.clear();
				errorlist.clear();
				std::cout << "Timer and error list get initialized " << std::endl;
			}
			ImGui::SameLine();
			if (ImGui::Button("SaveInvertLS", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				if (readed_LS1.size() == 0)
				{
					std::cout << "please read the first levelset" << std::endl;
				}
				else
				{
					save_invert_levelset(readed_LS1);
					std::cout << "finish writing inverted ls" << std::endl;
				}
			}
			ImGui::SameLine();
			if (ImGui::Button("AssignConstRuling", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				// write_polyline_joints(tools.V, tools.F, tools.norm_v);
				assign_const_rulings_to_strips(Eigen::Vector3d(ptInx, ptIny, ptInz));
			}
			if (ImGui::Button("OthoLs2Joints", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{

				compute_ls_intersections_and_assign_parallel_joints(tools.lsmesh, Eigen::Vector3d(ptInx, ptIny, ptInz),
																	vector_scaling, tools.V, tools.F,
																	readed_LS1, readed_LS2,
																	nbr_lines_first_ls,
																	nbr_lines_second_ls,
																	false);
			}
			ImGui::SameLine();
			if (ImGui::Button("AssignLsToPatch", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{

				Eigen::VectorXd funout;
				assign_ls_to_subpatch(readed_mesh1[0], readed_mesh2, tools.fvalues, funout);
				std::cout << "write sub-patch scalar field" << std::endl;
				save_levelset(funout);
			}
			ImGui::SameLine();
			if (ImGui::Button("CmprCrvsDis", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				std::cout << "Compare two faimilies of curves max distance" << std::endl;
				compare_multiple_correspoding_curves_dis();
			}
			if (ImGui::Button("WriteErrorCSV", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				Eigen::VectorXd error_print;
				tools.show_max_pg_energy(error_print);
				if (error_print.size() == 0)
				{
					std::cout<<"the error list is empty"<<std::endl;
					ImGui::End();
				}
				std::cout << "write error on vertices, please provide the prefix" << std::endl;
				std::string fname = igl::file_dialog_save();
				if (fname.length() == 0)
				{
					ImGui::End();
				}
				else
				{
					std::ofstream file;
					file.open(fname+".csv");
					for (int i = 0; i < error_print.size(); i++)
					{
						file<<error_print[i]<<"\n";
					}
					file.close();
					
					std::cout << "csv file saved" << std::endl;
				}
				
			}
			ImGui::SameLine();
			if (ImGui::Button("TimeErrorItr", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				if(timelist.size()!=errorlist.size() || timelist.empty()){
					std::cout<<"the time list and error list need to be reset"<<std::endl;
					ImGui::End();
				}
				std::cout << "write Runtime, error v.s. itr. please provide the prefix" << std::endl;
				std::string fname = igl::file_dialog_save();
				if (fname.length() == 0)
				{
					ImGui::End();
				}
				else
				{
					std::ofstream file;
					file.open(fname+".csv");
					for (int i = 0; i < timelist.size(); i++)
					{
						file<<i<<","<<timelist[i]<<","<<errorlist[i]<<"\n";
					}
					file.close();
					
					std::cout << "csv file saved" << std::endl;
				}
				
			}
			ImGui::SameLine();
			if (ImGui::Button("printApproError", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{

				std::cout << "print the approximation error of the strips to the original surface" << std::endl;
				poly_tool.related_error();
			}
		}
		if (ImGui::CollapsingHeader("LS Processing", ImGuiTreeNodeFlags_DefaultOpen))
		{

			if (ImGui::Button("LvSet Opt", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{

				EnergyPrepare einit;
				einit.weight_gravity = weight_mass;
				einit.weight_lap = weight_laplacian;
				einit.weight_bnd = weight_boundary;
				einit.weight_pg = weight_pseudo_geodesic;
				einit.weight_strip_width = weight_strip_width;
				einit.solve_pseudo_geodesic = enable_pg_energy_checkbox;
				einit.target_angle = target_angle;
				einit.max_step_length = maximal_step_length;
				einit.solve_strip_width_on_traced = enable_strip_width_checkbox;
				einit.enable_extreme_cases = enable_extreme_cases;
				einit.Given_Const_Direction = Given_Const_Direction;
				// tools.weight_shading = weight_shading;
				tools.prepare_level_set_solving(einit);
				tools.Reference_theta = InputPx;
				tools.Reference_phi = InputPy;
				tools.Reference_theta1 = InputPx1;
				tools.Reference_phi1 = InputPy1;
				tools.Theta_tol = InputThetaTol;
				tools.Phi_tol = InputPhiTol;
				tools.Theta_tol1 = InputThetaTol1;
				tools.Phi_tol1 = InputPhiTol1;
				tools.weight_geodesic = weight_geodesic;
				tools.enable_shading_init = shading_init;
				tools.weight_binormal = weight_binormal;
				tools.weight_smt_binormal = weight_smt_binormal;
				tools.enable_let_ray_through = let_ray_through;
				tools.ShadingLatitude = Shading_Latitude;
				tools.enable_reflection = let_ray_reflect;
				tools.recompute_auxiliaries = recompute_auxiliaries;
				igl::Timer lctimer;
				timer_global.start();
				for (int i = 0; i < OpIter; i++)
				{
					lctimer.start();
					tools.Run_Level_Set_Opt();
					lctimer.stop();
					errorlist.push_back(tools.pgerror);
					if (timelist.empty())
					{
						timelist.push_back(lctimer.getElapsedTimeInSec());
					}
					else
					{
						timelist.push_back(timelist.back() + lctimer.getElapsedTimeInSec());
					}
					iteration_total++;
					if (tools.step_length < 1e-16 && i != 0)
					{ // step length actually is the value for the last step
						std::cout << "optimization converges " << std::endl;
						break;
					}
				}
				timer_global.stop();
				time_total += timer_global.getElapsedTimeInSec();
				std::cout << "waiting for instruction..." << std::endl;
				// MP.MeshUnitScale(inputMesh, updatedMesh);
				int id = viewer.selected_data_index;
				CGMesh inputMesh = tools.lsmesh;
				updateMeshViewer(viewer, inputMesh);
				meshFileName.push_back("lso_" + meshFileName[id]);
				Meshes.push_back(inputMesh);

				Eigen::VectorXd level_set_values;
				tools.show_level_set(level_set_values);
				if (level_set_values.size() == 0)
				{
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
				// // tools.show_gradients(E0,E1, vector_scaling);
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);

				viewer.selected_data_index = id;
			}
			ImGui::SameLine();
			if (ImGui::Button("Mesh Opt", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				if (tools.fvalues.size() == 0)
				{
					std::cout << "Please load or initialize the level-set before opt the mesh" << std::endl;
					ImGui::End();
					return;
				}
				int id = viewer.selected_data_index;
				CGMesh updatedMesh;
				CGMesh inputMesh = Meshes[id];
				MeshEnergyPrepare initializer;
				initializer.Mesh_opt_max_step_length = Mesh_opt_max_step_length;
				initializer.weight_Mesh_pesudo_geodesic = weight_Mesh_pesudo_geodesic;
				initializer.weight_Mesh_smoothness = weight_Mesh_smoothness;
				initializer.weight_mass = weight_mass;
				initializer.target_angle = target_angle;
				initializer.weight_Mesh_edgelength = weight_Mesh_edgelength;
				initializer.enable_extreme_cases = enable_extreme_cases;
				initializer.Given_Const_Direction = Given_Const_Direction;
				tools.weight_Mesh_approximation = weight_Mesh_approximation;
				tools.weight_Mesh_mass = weight_Mesh_mass;
				// tools.weight_shading = weight_shading;
				tools.Reference_theta = InputPx;
				tools.Reference_phi = InputPy;
				// tools.Reference_theta2 = InputPx1;
				// tools.Reference_phi2 = InputPy1;
				tools.Theta_tol = InputThetaTol;
				tools.Phi_tol = InputPhiTol;
				// tools.Theta_tol2 = InputPhiTol;
				// tools.Phi_tol2 = InputPhiTol1;
				tools.prepare_mesh_optimization_solving(initializer);
				timer_global.start();
				for (int i = 0; i < Nbr_Iterations_Mesh_Opt; i++)
				{
					tools.Run_Mesh_Opt();
					iteration_total++;
				}
				timer_global.stop();
				time_total += timer_global.getElapsedTimeInSec();
				updatedMesh = tools.lsmesh;
				updateMeshViewer(viewer, updatedMesh);
				meshFileName.push_back("mso_" + meshFileName[id]);
				Meshes.push_back(updatedMesh);

				viewer.selected_data_index = id;
				std::cout << "waiting for instructions" << std::endl;
			}
			ImGui::SameLine();
			if (ImGui::Button("Orthogonal LS", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				int id = viewer.selected_data_index;
				CGMesh inputMesh = tools.lsmesh;
				EnergyPrepare einit;
				einit.weight_gravity = weight_mass;
				einit.weight_lap = weight_laplacian;
				einit.weight_bnd = weight_boundary;
				einit.weight_pg = weight_pseudo_geodesic;
				einit.weight_strip_width = weight_strip_width;
				einit.solve_pseudo_geodesic = enable_pg_energy_checkbox;
				einit.target_angle = target_angle;
				einit.max_step_length = maximal_step_length;
				einit.solve_strip_width_on_traced = enable_strip_width_checkbox;
				einit.enable_extreme_cases = enable_extreme_cases;
				tools.fix_angle_of_two_levelsets = fix_angle_of_two_levelsets;
				tools.angle_between_two_levelsets = angle_between_two_levelsets;
				tools.weight_fix_two_ls_angle = weight_fix_two_ls_angle;
				// tools.weight_shading = weight_shading;
				tools.prepare_level_set_solving(einit);
				if (readed_LS1.size() == 0)
				{
					std::cout << "Please load LevelSet 1 as a reference" << std::endl;
					ImGui::End();
					return;
				}
				for (int i = 0; i < OpIter; i++)
				{
					tools.Run_Othogonal_Levelset(readed_LS1);
					if (tools.step_length < 1e-16 && i != 0)
					{ // step length actually is the value for the last step
						std::cout << "optimization converges " << std::endl;
						break;
					}
				}
				std::cout << "waiting for instruction..." << std::endl;
				// MP.MeshUnitScale(inputMesh, updatedMesh);
				updateMeshViewer(viewer, inputMesh);
				meshFileName.push_back("lso_" + meshFileName[id]);
				Meshes.push_back(inputMesh);

				Eigen::VectorXd level_set_values;
				tools.show_level_set(level_set_values);
				if (level_set_values.size() == 0)
				{
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
				// // tools.show_gradients(E0,E1, vector_scaling);
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);

				viewer.selected_data_index = id;
			}
			ImGui::SameLine();
			if (ImGui::Button("AngleVariable", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{

				EnergyPrepare einit;
				einit.weight_gravity = weight_mass;
				einit.weight_lap = weight_laplacian;
				einit.weight_bnd = weight_boundary;
				einit.weight_pg = weight_pseudo_geodesic;
				einit.weight_strip_width = weight_strip_width;
				einit.solve_pseudo_geodesic = enable_pg_energy_checkbox;
				einit.target_angle = target_angle;
				// std::cout<<"Target angle, "<<target_angle<<std::endl; exit(0);
				einit.max_step_length = maximal_step_length;
				einit.solve_strip_width_on_traced = enable_strip_width_checkbox;
				einit.enable_extreme_cases = enable_extreme_cases;
				einit.Given_Const_Direction = Given_Const_Direction;
				// tools.weight_shading = weight_shading;
				tools.prepare_level_set_solving(einit);
				tools.Reference_theta = InputPx;
				tools.Reference_phi = InputPy;
				tools.Reference_theta1 = InputPx1;
				tools.Reference_phi1 = InputPy1;
				tools.Theta_tol = InputThetaTol;
				tools.Phi_tol = InputPhiTol;
				tools.Theta_tol1 = InputThetaTol1;
				tools.Phi_tol1 = InputPhiTol1;
				tools.weight_geodesic = weight_geodesic;
				tools.enable_shading_init = shading_init;
				tools.weight_binormal = weight_binormal;
				tools.weight_smt_binormal = weight_smt_binormal;
				tools.enable_let_ray_through = let_ray_through;
				tools.ShadingLatitude = Shading_Latitude;
				tools.enable_reflection = let_ray_reflect;
				tools.recompute_auxiliaries = recompute_auxiliaries;

				for (int i = 0; i < OpIter; i++)
				{
					tools.Run_Level_Set_Opt_Angle_Variable();
					if (tools.step_length < 1e-16 && i != 0)
					{ // step length actually is the value for the last step
						std::cout << "optimization converges " << std::endl;
						break;
					}
				}
				std::cout << "waiting for instruction..." << std::endl;
				// MP.MeshUnitScale(inputMesh, updatedMesh);
				int id = viewer.selected_data_index;
				CGMesh inputMesh = tools.lsmesh;
				updateMeshViewer(viewer, inputMesh);
				meshFileName.push_back("lso_" + meshFileName[id]);
				Meshes.push_back(inputMesh);

				Eigen::VectorXd level_set_values;
				tools.show_level_set(level_set_values);
				if (level_set_values.size() == 0)
				{
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
				// // tools.show_gradients(E0,E1, vector_scaling);
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);

				viewer.selected_data_index = id;
			}
			if (ImGui::Button("AAG LvSet", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				if (readed_LS1.size() == 0 || readed_LS2.size() == 0)
				{
					std::cout << "Please Read the correct levelset 1 and levelset 2" << std::endl;
					ImGui::End();
					return;
				}
				int id = viewer.selected_data_index;
				CGMesh inputMesh = Meshes[id];
				EnergyPrepare einit;
				einit.weight_gravity = weight_mass;
				einit.weight_lap = weight_laplacian;
				einit.weight_bnd = weight_boundary;
				einit.weight_pg = weight_pseudo_geodesic;
				einit.weight_strip_width = weight_strip_width;
				einit.solve_pseudo_geodesic = enable_pg_energy_checkbox;
				einit.target_angle = target_angle;
				einit.max_step_length = maximal_step_length;
				einit.solve_strip_width_on_traced = enable_strip_width_checkbox;
				einit.enable_extreme_cases = enable_extreme_cases;

				tools.prepare_level_set_solving(einit);
				tools.weight_geodesic = weight_geodesic;
				timer_global.start();
				for (int i = 0; i < OpIter; i++)
				{
					tools.Run_AAG(readed_LS1, readed_LS2, readed_LS3);
					iteration_total++;
					if (tools.step_length < 1e-16 && i != 0)
					{ // step length actually is the value for the last step
						std::cout << "optimization converges " << std::endl;
						break;
					}
				}
				timer_global.stop();
				time_total += timer_global.getElapsedTimeInSec();

				std::cout << "waiting for instruction..." << std::endl;
				// MP.MeshUnitScale(inputMesh, updatedMesh);
				updateMeshViewer(viewer, inputMesh);
				meshFileName.push_back("aag_" + meshFileName[id]);
				Meshes.push_back(inputMesh);

				viewer.selected_data_index = id;
			}
			ImGui::SameLine();

			if (ImGui::Button("AAG MshOpt", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				if (readed_LS1.size() == 0 || readed_LS2.size() == 0)
				{
					std::cout << "Please Read the correct levelset 1 and levelset 2" << std::endl;
					ImGui::End();
					return;
				}
				int id = viewer.selected_data_index;
				CGMesh updatedMesh;
				CGMesh inputMesh = Meshes[id];
				MeshEnergyPrepare initializer;
				initializer.Mesh_opt_max_step_length = Mesh_opt_max_step_length;
				initializer.weight_Mesh_pesudo_geodesic = weight_Mesh_pesudo_geodesic;
				initializer.weight_Mesh_smoothness = weight_Mesh_smoothness;
				initializer.weight_mass = weight_mass;
				initializer.target_angle = target_angle;
				initializer.weight_Mesh_edgelength = weight_Mesh_edgelength;
				initializer.enable_extreme_cases = enable_extreme_cases;
				tools.weight_geodesic = weight_geodesic;
				tools.weight_Mesh_mass = weight_Mesh_mass;
				tools.weight_Mesh_approximation = weight_Mesh_approximation;
				tools.prepare_mesh_optimization_solving(initializer);
				for (int i = 0; i < Nbr_Iterations_Mesh_Opt; i++)
				{
					tools.Run_AAG_Mesh_Opt(readed_LS1, readed_LS2, readed_LS3);
				}
				updatedMesh = tools.lsmesh;
				updateMeshViewer(viewer, updatedMesh);
				meshFileName.push_back("aagm_" + meshFileName[id]);
				Meshes.push_back(updatedMesh);

				viewer.selected_data_index = id;
				std::cout << "waiting for instructions" << std::endl;
			}
			ImGui::SameLine();
			if (ImGui::Button("AGG LvSet", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				if (readed_LS1.size() == 0 || readed_LS2.size() == 0)
				{
					std::cout << "Please Read the correct levelset 1 and levelset 2" << std::endl;
					ImGui::End();
					return;
				}
				int id = viewer.selected_data_index;
				CGMesh inputMesh = Meshes[id];
				EnergyPrepare einit;
				einit.weight_gravity = weight_mass;
				einit.weight_lap = weight_laplacian;
				einit.weight_bnd = weight_boundary;
				einit.weight_pg = weight_pseudo_geodesic;
				einit.weight_strip_width = weight_strip_width;
				einit.solve_pseudo_geodesic = enable_pg_energy_checkbox;
				einit.target_angle = target_angle;
				einit.max_step_length = maximal_step_length;
				einit.solve_strip_width_on_traced = enable_strip_width_checkbox;
				einit.enable_extreme_cases = enable_extreme_cases;
				tools.fix_angle_of_two_levelsets = fix_angle_of_two_levelsets;
				tools.angle_between_two_levelsets = angle_between_two_levelsets;
				tools.weight_fix_two_ls_angle = weight_fix_two_ls_angle;

				tools.prepare_level_set_solving(einit);
				tools.weight_geodesic = weight_geodesic;
				timer_global.start();
				for (int i = 0; i < OpIter; i++)
				{
					tools.Run_AGG(readed_LS1, readed_LS2, readed_LS3);
					iteration_total++;
					if (tools.step_length < 1e-16 && i != 0)
					{ // step length actually is the value for the last step
						std::cout << "optimization converges " << std::endl;
						break;
					}
				}
				timer_global.stop();
				time_total += timer_global.getElapsedTimeInSec();
				std::cout << "waiting for instruction..." << std::endl;
				// MP.MeshUnitScale(inputMesh, updatedMesh);
				updateMeshViewer(viewer, inputMesh);
				meshFileName.push_back("agg_" + meshFileName[id]);
				Meshes.push_back(inputMesh);

				viewer.selected_data_index = id;
			}
			ImGui::SameLine();
			if (ImGui::Button("AGG MshOpt", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				if (readed_LS1.size() == 0 || readed_LS2.size() == 0)
				{
					std::cout << "Please Read the correct levelset 1 and levelset 2" << std::endl;
					ImGui::End();
					return;
				}
				int id = viewer.selected_data_index;
				CGMesh updatedMesh;
				CGMesh inputMesh = Meshes[id];
				MeshEnergyPrepare initializer;
				initializer.Mesh_opt_max_step_length = Mesh_opt_max_step_length;
				initializer.weight_Mesh_pesudo_geodesic = weight_Mesh_pesudo_geodesic;
				initializer.weight_Mesh_smoothness = weight_Mesh_smoothness;
				initializer.weight_mass = weight_mass;
				initializer.target_angle = target_angle;
				initializer.weight_Mesh_edgelength = weight_Mesh_edgelength;
				initializer.enable_extreme_cases = enable_extreme_cases;
				tools.weight_geodesic = weight_geodesic;
				tools.weight_Mesh_mass = weight_Mesh_mass;
				tools.weight_Mesh_approximation = weight_Mesh_approximation;
				tools.prepare_mesh_optimization_solving(initializer);
				for (int i = 0; i < Nbr_Iterations_Mesh_Opt; i++)
				{
					tools.Run_AGG_Mesh_Opt(readed_LS1, readed_LS2, readed_LS3);
				}
				updatedMesh = tools.lsmesh;
				updateMeshViewer(viewer, updatedMesh);
				meshFileName.push_back("aggm_" + meshFileName[id]);
				Meshes.push_back(updatedMesh);

				viewer.selected_data_index = id;
				std::cout << "waiting for instructions" << std::endl;
			}

			if (ImGui::Button("PPG LvSet", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				if (readed_LS1.size() == 0 || readed_LS2.size() == 0)
				{
					std::cout << "Please Read the correct levelset 1 and levelset 2" << std::endl;
					ImGui::End();
					return;
				}
				int id = viewer.selected_data_index;
				CGMesh inputMesh = Meshes[id];
				EnergyPrepare einit;
				einit.weight_gravity = weight_mass;
				einit.weight_lap = weight_laplacian;
				einit.weight_bnd = weight_boundary;
				einit.weight_pg = weight_pseudo_geodesic;
				einit.weight_strip_width = weight_strip_width;
				einit.solve_pseudo_geodesic = enable_pg_energy_checkbox;
				einit.target_angle = target_angle;
				einit.max_step_length = maximal_step_length;
				einit.solve_strip_width_on_traced = enable_strip_width_checkbox;
				einit.enable_extreme_cases = enable_extreme_cases;
				tools.fix_angle_of_two_levelsets = fix_angle_of_two_levelsets;
				tools.angle_between_two_levelsets = angle_between_two_levelsets;
				tools.weight_fix_two_ls_angle = weight_fix_two_ls_angle;
				tools.weight_geodesic = weight_geodesic;
				tools.pseudo_geodesic_target_angle_degree = target_angle;
				tools.pseudo_geodesic_target_angle_degree_2 = target_angle_2;

				tools.prepare_level_set_solving(einit);
				timer_global.start();
				for (int i = 0; i < OpIter; i++)
				{
					tools.Run_PPG(readed_LS1, readed_LS2, readed_LS3);
					iteration_total++;
					if (tools.step_length < 1e-16 && i != 0)
					{ // step length actually is the value for the last step
						std::cout << "optimization converges " << std::endl;
						break;
					}
				}
				timer_global.stop();
				time_total += timer_global.getElapsedTimeInSec();
				std::cout << "waiting for instruction..." << std::endl;
				// MP.MeshUnitScale(inputMesh, updatedMesh);
				updateMeshViewer(viewer, inputMesh);
				meshFileName.push_back("agg_" + meshFileName[id]);
				Meshes.push_back(inputMesh);

				viewer.selected_data_index = id;
			}
			ImGui::SameLine();
			if (ImGui::Button("draw AAG", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				if (readed_LS1.size() == 0)
				{
					std::cout << "Please Read the correct levelset 1 and levelset 2" << std::endl;
					ImGui::End();
					return;
				}
				Eigen::MatrixXd CM;
				igl::parula(Eigen::VectorXd::LinSpaced(21, 0, 1).eval(), false, CM);
				igl::isolines_map(Eigen::MatrixXd(CM), CM);
				viewer.data().set_colormap(CM);
				if (which_levelset == 0)
				{
					viewer.data().set_data(readed_LS1);
				}
				if (which_levelset == 1)
				{
					viewer.data().set_data(readed_LS2);
				}
				if (which_levelset == 2)
				{
					viewer.data().set_data(readed_LS3);
				}
			}
			ImGui::SameLine();
			if (ImGui::Button("Save AAG", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				bool saved = save_levelset(readed_LS1);
				if (saved)
				{
					std::cout << "\nlevel 1  get saved" << std::endl;
				}
				saved = save_levelset(readed_LS2);
				if (saved)
				{
					std::cout << "\nlevel 2  get saved" << std::endl;
				}
				saved = save_levelset(readed_LS3);
				if (saved)
				{
					std::cout << "\nlevel 3  get saved" << std::endl;
				}
			}

			if (ImGui::Button("draw extracted", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				if (tools.fvalues.size() == 0)
				{
					std::cout << "Please init the levelset or load a correct levelset" << std::endl;
					ImGui::End();
					return;
				}

				int id = viewer.selected_data_index;
				CGMesh inputMesh = Meshes[id];

				// MP.MeshUnitScale(inputMesh, updatedMesh);
				updateMeshViewer(viewer, inputMesh);
				meshFileName.push_back("extracted_" + meshFileName[id]);
				Meshes.push_back(inputMesh);

				// viewer.data().set_colors(level_set_values);
				Eigen::MatrixXd E0, E1;
				// // tools.show_gradients(E0,E1, vector_scaling);
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);

				const Eigen::RowVector3d black(0, 0, 0);
				const Eigen::RowVector3d green(0.2, 0.8, 0.2);
				std::vector<Eigen::MatrixXd> E0list, E1list;
				tools.extract_levelset_curves(extracted_nbr, E0list, E1list);
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

				tools.show_pseudo_geodesic_curve(E0list, E1list, pts);
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
				tools.show_level_set(level_set_values);
				if (level_set_values.size() == 0)
				{
					std::cout << "Please init the levelset or load a correct levelset" << std::endl;
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
				const Eigen::RowVector3d yellow(241. / 255, 196. / 255, 15. / 255);

				const Eigen::RowVector3d black(0, 0, 0);
				const Eigen::RowVector3d green(0.2, 0.8, 0.2);

				Eigen::MatrixXd pts, E0, E1;
				std::cout << "before ploting the bi-normals" << std::endl;
				Eigen::MatrixXd binormals;
				tools.show_binormals(tools.fvalues, E0, E1, binormals, vector_scaling);
				viewer.data().add_edges(E0, E1, green);
				viewer.data().add_edges(tools.V, tools.V + vector_scaling * tools.norm_v, black);
				if (enable_extreme_cases && Given_Const_Direction)
				{
					std::cout << "Ploting the light directions" << std::endl;
					Eigen::MatrixXd light(tools.V.rows(), 3);
					if (tools.Lights.rows() == 0)
					{
						Eigen::Vector3d direction;
						direction = angle_ray_converter(InputPx, InputPy);

						for (int i = 0; i < light.rows(); i++)
						{
							light.row(i) = direction;
						}
					}
					else
					{
						light = tools.Lights;
					}
					Eigen::MatrixXd E2, E3;
					E2 = tools.V + vector_scaling * light;
					E3 = tools.V - vector_scaling * light;
					viewer.data().add_edges(E2, E3, yellow);
				}
			}

			if (ImGui::Button("draw error", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				Eigen::VectorXd error;
				tools.show_max_pg_energy(error);
				if (error.size() == 0)
				{
					std::cout << "Please run Levelset Opt to get the pseudo-geodesic error" << std::endl;
					ImGui::End();
					return;
				}

				viewer.data().set_data(error);
			}
			ImGui::SameLine();
			if (ImGui::Button("draw RefPts", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				std::cout << "drawing the reference points in yellow" << std::endl;
				Eigen::MatrixXd pts;
				tools.show_current_reference_points(pts);
				const Eigen::RowVector3d yellow(241. / 255, 196. / 255, 15. / 255);
				viewer.data().add_points(pts, yellow);
			}
			ImGui::SameLine();
			if (ImGui::Button("draw Gauss", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				std::cout << "drawing the sign of Gaussian curvature" << std::endl;

				Eigen::MatrixXd pts, color;
				tools.show_posi_negt_curvature(color);
				viewer.data().add_points(tools.V, color);
			}
			ImGui::SameLine();
			if (ImGui::Button("draw Paramtrztn", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				int id = viewer.selected_data_index;
				CGMesh inputMesh = tools.write_parameterization_mesh();

				updateMeshViewer(viewer, inputMesh);
				meshFileName.push_back("lso_" + meshFileName[id]);
				Meshes.push_back(inputMesh);
				viewer.selected_data_index = id;

				Eigen::VectorXd level_set_values;
				tools.show_level_set(level_set_values);
				if (level_set_values.size() == 0)
				{
					std::cout << "Please init the levelset or load a correct levelset" << std::endl;
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
				mesh_unit_scale(tools.V, vout);
				std::string fname = igl::file_dialog_save();
				if (fname.length() == 0)
				{
					ImGui::End();
				}
				else
				{
					igl::writeOBJ(fname, vout, tools.F);
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
				if (Meshes.size() == 1)
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
						meshFileName.erase(meshFileName.begin() + id);
						Meshes.erase(Meshes.begin() + id);
					}
					if (Meshes.size() > id)
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

			if (ImGui::Button("ReadPlylines", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				int id = viewer.selected_data_index;
				read_plylines_and_binormals(Poly_readed, Bino_readed);
				if (Poly_readed.empty())
				{
					std::cout << "ERROR, Please load the polylines first" << std::endl;
					ImGui::End();
					return;
				}
				CGMesh updatemesh = polyline_to_strip_mesh(Poly_readed, Bino_readed, vector_scaling);
				updateMeshViewer(viewer, updatemesh);
				meshFileName.push_back("ply_" + meshFileName[id]);
				Meshes.push_back(updatemesh);
				poly_tool.init(Poly_readed, Bino_readed, nbr_lines_first_ls);
				poly_tool.get_normal_vector_from_reference(tools.V, tools.F, tools.norm_v);
				std::cout << "Sample nbr: " << nbr_lines_first_ls << std::endl;
				viewer.selected_data_index = id;
			}
			ImGui::SameLine();
			if (ImGui::Button("OptiPolylines", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				if (Poly_readed.empty())
				{
					std::cout << "ERROR, Please load the polylines first" << std::endl;
					ImGui::End();
					return;
				}
				poly_tool.weight_smooth = weight_laplacian;
				poly_tool.weight_mass = weight_mass;
				poly_tool.weight_binormal = weight_pseudo_geodesic;
				poly_tool.max_step = maximal_step_length;
				poly_tool.strip_scale = vector_scaling;
				poly_tool.binormal_ratio = weight_geodesic;
				poly_tool.weight_angle = weight_angle;
				poly_tool.target_angle = target_angle;
				poly_tool.ratio_endpts = weight_endpoint_ratio;
				poly_tool.pick_single_line = pick_single_ply;
				poly_tool.pick_line_id = pick_line_id;
				if (poly_tool.Vref.rows() == 0)
				{
					poly_tool.get_normal_vector_from_reference(tools.V, tools.F, tools.norm_v);
					std::cout << "Reading Reference Triangle Mesh ... " << std::endl;
				}

				for (int i = 0; i < OpIter; i++)
				{

					poly_tool.opt();
				}
				CGMesh updatemesh = poly_tool.RecMesh;
				int id = viewer.selected_data_index;
				updateMeshViewer(viewer, updatemesh);
				meshFileName.push_back("ply_" + meshFileName[id]);
				Meshes.push_back(updatemesh);
				viewer.selected_data_index = id;
				std::cout << "waiting for instructions" << std::endl;
			}

			ImGui::SameLine();
			if (ImGui::Button("SavePolylines", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				if (poly_tool.ply_extracted.empty())
				{
					std::cout << "ERROR, Please opt the polyline first" << std::endl;
					ImGui::End();
					return;
				}
				poly_tool.save_polyline_and_binormals_as_files(false); // true means save the rotated polylines
			}
			ImGui::SameLine();
			if (ImGui::Button("Orient&SavePlys", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))

			{
				std::cout << "this function is to orient binormals if there are inverted ones" << std::endl;
				int id = viewer.selected_data_index;
				read_plylines_and_binormals(Poly_readed, Bino_readed);
				std::vector<std::vector<Eigen::Vector3d>> biout;
				poly_tool.orient_binormals_of_plyline(Bino_readed, biout);
				CGMesh updatemesh = polyline_to_strip_mesh(Poly_readed, biout, vector_scaling);
				updateMeshViewer(viewer, updatemesh);
				meshFileName.push_back("ply_" + meshFileName[id]);
				Meshes.push_back(updatemesh);
				viewer.selected_data_index = id;
				poly_tool.save_polyline_and_binormals_as_files(Poly_readed, biout); // true means save the rotated polylines
			}
			if (ImGui::Button("RotateBackPly", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				if (poly_tool.ply_extracted.empty())
				{
					std::cout << "ERROR, Please opt the polyline first" << std::endl;
					ImGui::End();
					return;
				}
				std::cout << "Rotating the mesh to the coordinate system computed by latitude" << std::endl;
				poly_tool.rotate_back_the_model_to_horizontal_coordinates(Shading_Latitude);
				poly_tool.save_polyline_and_binormals_as_files(true); // true means save the rotated polylines
			}
			ImGui::SameLine();
			if (ImGui::Button("ForcePlyBinm", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				if (poly_tool.ply_extracted.empty())
				{
					std::cout << "ERROR, Please opt the polyline first" << std::endl;
					ImGui::End();
					return;
				}
				poly_tool.force_smoothing_binormals();
				std::cout << "The binormals are forced to be smoothed" << std::endl;

				CGMesh updatemesh = poly_tool.RecMesh;
				int id = viewer.selected_data_index;
				updateMeshViewer(viewer, updatemesh);
				meshFileName.push_back("Smt_" + meshFileName[id]);
				Meshes.push_back(updatemesh);
				viewer.selected_data_index = id;
				std::cout << "waiting for instructions" << std::endl;
			}
			ImGui::SameLine();
			if (ImGui::Button("DevelopableStrips", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				if (Poly_readed.empty())
				{
					std::cout << "ERROR, Please read the polyline first" << std::endl;
					ImGui::End();
					return;
				}
				// Poly_readed, Bino_readed
				std::vector<std::vector<Eigen::Vector3d>> dev_crease;
				for (int i = 0; i < Poly_readed.size(); i++)
				{
					std::vector<Eigen::Vector3d> creases;
					convert_polyline_to_developable_strips_reflection_method(Poly_readed[i], Bino_readed[i], creases);
					dev_crease.push_back(creases);
				}
				int id = viewer.selected_data_index;
				CGMesh updatemesh = polyline_to_strip_mesh(Poly_readed, dev_crease, vector_scaling);
				updateMeshViewer(viewer, updatemesh);
				meshFileName.push_back("Dev_" + meshFileName[id]);
				Meshes.push_back(updatemesh);
				viewer.selected_data_index = id;
				poly_tool.save_polyline_and_binormals_as_files(Poly_readed, dev_crease);
				std::cout << "waiting for instructions" << std::endl;
			}
			ImGui::SameLine();
			if (ImGui::Button("EvltDvlp", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				std::cout << "evaluate developable of the loaded strips ..." << std::endl;
				// get_polyline_rectifying_plane_on_inner_vers(const std::vector<Eigen::Vector3d> &ply_in, const std::vector<Eigen::Vector3d> &bnm_in,
				//                                      std::vector<Eigen::Vector3d> &vertices, std::vector<Eigen::Vector3d> &tangents, std::vector<Eigen::Vector3d> &binormals);
				evaluate_and_print_strip_developability(Poly_readed, Bino_readed);
			}
			if (ImGui::Button("ReadCreases", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				int id = viewer.selected_data_index;
				read_plylines_and_binormals(Poly_readed, Bino_readed);
				if (Poly_readed.empty())
				{
					std::cout << "ERROR, Please load the polylines first" << std::endl;
					ImGui::End();
					return;
				}
				CGMesh updatemesh = polyline_to_strip_mesh(Poly_readed, Bino_readed, vector_scaling);
				updateMeshViewer(viewer, updatemesh);
				meshFileName.push_back("ply_" + meshFileName[id]);
				Meshes.push_back(updatemesh);

				std::vector<std::vector<Eigen::Vector3d>> vertices;
				std::vector<std::vector<Eigen::Vector3d>> tangents;
				std::vector<std::vector<Eigen::Vector3d>> binormals;
				get_polyline_rectifying_planes(Poly_readed, Bino_readed, vertices, tangents, binormals);
				poly_tool.init_crease_opt(vertices, tangents, binormals);

				viewer.selected_data_index = id;
				evaluate_strip_straightness(Poly_readed, Bino_readed);
			}
			ImGui::SameLine();
			if (ImGui::Button("OptCrease", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				if (Poly_readed.empty())
				{
					std::cout << "ERROR, Please load the polylines first" << std::endl;
					ImGui::End();
					return;
				}
				poly_tool.weight_smooth = weight_laplacian;
				poly_tool.weight_mass = weight_mass;
				poly_tool.weight_binormal = weight_pseudo_geodesic;
				poly_tool.max_step = maximal_step_length;
				poly_tool.strip_scale = vector_scaling;
				poly_tool.binormal_ratio = weight_geodesic;
				poly_tool.weight_angle = weight_angle;
				poly_tool.target_angle = target_angle;
				poly_tool.ratio_endpts = weight_endpoint_ratio;
				poly_tool.pick_single_line = pick_single_ply;
				poly_tool.pick_line_id = pick_line_id;

				for (int i = 0; i < OpIter; i++)
				{
					// if (weight_angle > 0)
					// {
					// 	poly_tool.get_normal_vector_from_reference(tools.V, tools.F, tools.norm_v);
					// }
					poly_tool.opt_planarity();
				}
				CGMesh updatemesh = poly_tool.RecMesh;
				int id = viewer.selected_data_index;
				updateMeshViewer(viewer, updatemesh);
				meshFileName.push_back("ply_" + meshFileName[id]);
				Meshes.push_back(updatemesh);
				viewer.selected_data_index = id;
				std::cout << "waiting for instructions" << std::endl;
			}
			ImGui::SameLine();
			if (ImGui::Button("DevV1", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				construct_developable_strips_by_intersect_rectifying();
			}
			ImGui::SameLine();
			if (ImGui::Button("RecoverEndPts", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				recover_polyline_endpts();
			}
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
				ImGui::Checkbox("Visible", &selectFEV[i]);
				ImGui::PopID();
				ImGui::SameLine();
				if (selectFEV[i])
					viewer.data_list[i].is_visible = true;
				else
					viewer.data_list[i].is_visible = false;

				node_flags |= ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_NoTreePushOnOpen; // ImGuiTreeNodeFlags_Bullet
				// node_flags |=  ImGuiTreeNodeFlags_Leaf | ImGuiTreeNodeFlags_DefaultOpen; // ImGuiTreeNodeFlags_Bullet
				// ImGui::TreeNodeEx((void*)(intptr_t)i, node_flags, "Object %d", i);
				ImGui::TreeNodeEx((void *)(intptr_t)i, node_flags, meshFileName[i].c_str());

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
}

void lscif::draw_menu2(igl::opengl::glfw::Viewer &viewer, igl::opengl::glfw::imgui::ImGuiPlugin &plugin,
					   igl::opengl::glfw::imgui::ImGuiMenu &menu_mp)
{
	menu_mp.callback_draw_custom_window = [&]()
	{
		// Define next window position + size
		ImGui::SetNextWindowPos(ImVec2(1100.f * menu_mp.menu_scaling(), 10), ImGuiCond_FirstUseEver);
		ImGui::SetNextWindowSize(ImVec2(300, 800), ImGuiCond_FirstUseEver);
		ImGui::Begin(
			"Panel 1", nullptr,
			ImGuiWindowFlags_NoSavedSettings);
		if (ImGui::CollapsingHeader("Quad Mesh Extraction", ImGuiTreeNodeFlags_DefaultOpen))
		{
			if (ImGui::Button("Load First Level Set", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{

				read_levelset(readed_LS1);
				std::cout << "Levelset 1 get readed " << std::endl;
			}
			ImGui::SameLine();
			if (ImGui::Button("Load Second Level Set", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{

				read_levelset(readed_LS2);
				std::cout << "Levelset 2 get readed " << std::endl;
			}
			ImGui::SameLine();
			if (ImGui::Button("Clear Readed Level Sets", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				readed_LS1.resize(0);
				readed_LS2.resize(0);
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
				readed_mesh1.push_back(readed_mesh);
				std::cout << "mesh 1 get readed, current mesh nbr: " << readed_mesh1.size() << std::endl;
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

				OpenMesh::IO::read_mesh(readed_mesh2, fname);
				std::cout << "mesh 2 get readed " << std::endl;
			}
			ImGui::SameLine();
			if (ImGui::Button("ShadingTypes", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				Eigen::VectorXi info;
				Eigen::MatrixXd P1;
				Eigen::MatrixXd P2;
				if (readed_mesh1.empty())
				{
					std::cout << "\nLSC: Please read meshes" << std::endl;
					ImGui::End();
					return;
				}

				project_mesh_and_get_shading_info(readed_mesh1[0], readed_mesh2, ring_nbr, info, P1, P2);
				std::cout << "the first type is green, the second is red" << std::endl;

				viewer.data().add_points(P2, hot_red);
				viewer.data().add_points(P1, sea_green);
				std::cout << "Nbr first, " << P1.rows() << std::endl;
				std::cout << "Nbr second, " << P2.rows() << std::endl;

				std::cout << "Passing the info to the level-set system" << std::endl;
				tools.shading_condition_info = info;
			}
			ImGui::SameLine();
			if (ImGui::Button("PGTypes", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				Eigen::VectorXi info;
				Eigen::MatrixXd P1;
				Eigen::MatrixXd P2;
				if (readed_mesh1.empty())
				{
					std::cout << "\nLSC: Please read meshes" << std::endl;
					ImGui::End();
					return;
				}

				project_mesh_and_get_vertex_class(readed_mesh1, readed_mesh2, info);

				std::cout << "Passing the info to the level-set system" << std::endl;
				tools.shading_condition_info = info;
			}

			if (ImGui::CollapsingHeader("Quad Mesh Size Paras", ImGuiTreeNodeFlags_DefaultOpen))
			{
				// Expose variable directly ...
				// ImGui::InputDouble("Closeness", &weigth_closeness, 0, 0, "%.4f");
				ImGui::InputInt("Nbr of U lines", &nbr_lines_first_ls, 0, 0);
				ImGui::InputInt("Nbr of V lines", &nbr_lines_second_ls, 0, 0);
				// ImGui::InputInt("Iteration", &OpIter, 0, 0);
				// ImGui::InputDouble("weight ls mass(big)", &weight_mass, 0, 0, "%.4f");

				// ImGui::Checkbox("Fix Boundary", &fixBoundary_checkbox);
			}

			// Add a button
			if (ImGui::Button("Otho Fields Quads", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				int id = viewer.selected_data_index;
				CGMesh updatedMesh;

				Eigen::MatrixXd VER;
				Eigen::MatrixXi FAC;
				if (readed_LS1.size() > 0 && readed_LS2.size() > 0)
				{
					extract_levelset_web(tools.lsmesh, tools.V, tools.F, readed_LS1, readed_LS2, nbr_lines_first_ls,
										 nbr_lines_second_ls, 3, VER, FAC, false);
				}
				else
				{
					std::cout << "ERROR, Please load both two level sets" << std::endl;
					ImGui::End();
					return;
				}
				std::cout << "Saving the quad mesh" << std::endl;
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
				visual_extract_levelset_web(tools.lsmesh, tools.V, tools.F, readed_LS1, readed_LS2, nbr_lines_first_ls,
											nbr_lines_second_ls, E0, E1, E2, E3, false);
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
				std::cout << "This code is to extract AAG quads" << std::endl;
				int id = viewer.selected_data_index;
				CGMesh inputMesh = Meshes[id];
				Eigen::MatrixXd VER;
				Eigen::MatrixXi FAC;
				if (!debug_flag)
				{
					if (readed_LS1.size() > 0 && readed_LS2.size() > 0)
					{
						int filter_nbr = 3;
						extract_levelset_web_stable(tools.lsmesh, tools.Boundary_Edges, tools.V, tools.F, readed_LS1, readed_LS2, nbr_lines_first_ls,
													nbr_lines_second_ls, filter_nbr, VER, FAC, true);
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
				visual_extract_levelset_web_stable(tools.lsmesh, tools.Boundary_Edges, tools.V, tools.F, readed_LS1, readed_LS2, nbr_lines_first_ls,
												   nbr_lines_second_ls, E0, E1, E2, E3, true, debug_flag, dbg_int, dbg_int2);
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);
				const Eigen::RowVector3d black(0, 0, 0);
				const Eigen::RowVector3d green(0.2, 0.8, 0.2);

				// MP.MeshUnitScale(inputMesh, updatedMesh);
				updateMeshViewer(viewer, inputMesh);
				meshFileName.push_back("paceqd_" + meshFileName[id]);
				Meshes.push_back(inputMesh);
				viewer.data().add_edges(E0, E1, red);
				viewer.data().add_edges(E2, E3, green);
				viewer.selected_data_index = id;
			}
			ImGui::SameLine();
			if (ImGui::Button("ExtrLines&Binormals", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				int id = viewer.selected_data_index;

				if (readed_LS1.size() != tools.V.rows())
				{
					std::cout << "Please load a correct Level Set 1" << std::endl;
					ImGui::End();
					return;
				}
				extract_shading_lines(tools.lsmesh, tools.V, tools.Boundary_Edges,
									  tools.F, readed_LS1, nbr_lines_first_ls, true);
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
				std::cout << "reading Quad obj" << std::endl;
				std::string fname = igl::file_dialog_open();
				if (fname.length() == 0)
				{
					std::cout << "\nLSC: read mesh failed" << std::endl;
					ImGui::End();
					return;
				}
				igl::readOBJ(fname, VER, FAC);
				std::cout << "OBJ readed, ver in face " << FAC.cols() << std::endl;
				std::cout << "Reading Bi-normals of triangle mesh" << std::endl;
				Eigen::MatrixXd bn;
				read_bi_normals(bn);
				std::cout << "Binormals readed" << std::endl;
				std::cout << "generating mesh with bi-normals" << std::endl;
				fname = igl::file_dialog_save();

				if (fname.length() == 0)
				{
					ImGui::End();
				}

				write_quad_mesh_with_binormal(fname, tools.V, tools.F, bn, VER, FAC);
			}
			if (ImGui::Button("ExtractLines", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				int id = viewer.selected_data_index;

				Eigen::MatrixXd VER;
				Eigen::MatrixXi FAC;
				if (readed_LS1.size() != tools.V.rows())
				{
					std::cout << "Please load a correct Level Set 1" << std::endl;
					ImGui::End();
					return;
				}
				extract_shading_lines(tools.lsmesh, tools.V, tools.Boundary_Edges,
									  tools.F, readed_LS1, nbr_lines_first_ls, false);
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
				read_plylines_and_binormals(Poly_readed, Bino_readed);
				poly_tool.polyline_to_matrix(Poly_readed, Vply);

				int size1 = Vtri.rows();
				int size2 = Vply.rows();
				Eigen::MatrixXd Vall(size1 + size2, 3);
				Vall.topRows(size1) = Vtri;
				Vall.bottomRows(size2) = Vply;

				Eigen::MatrixXd vout;
				mesh_unit_scale(Vall, vout);
				Vtri = vout.topRows(size1);
				Vply = vout.bottomRows(size2);
				poly_tool.vertex_matrix_to_polyline(Poly_readed, Vply);

				std::cout << "Saving the result mesh and polylines, please provide the prifix" << std::endl;
				fname = igl::file_dialog_save();
				if (fname.length() == 0)
				{
					ImGui::End();
				}
				else
				{
					igl::writeOBJ(fname + ".obj", Vtri, tools.F);
					poly_tool.save_polyline_and_binormals_as_files(fname, Poly_readed, Bino_readed); // true means save the rotated polylines
					std::cout << "mesh saved" << std::endl;
				}
			}
		}

		ImGui::SameLine();
		if (ImGui::Button("ReadQuadPlyLines", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
		{
			int id = viewer.selected_data_index;
			read_origami_and_convert_to_polylines(Poly_readed, Bino_readed);
			CGMesh updatemesh = polyline_to_strip_mesh(Poly_readed, Bino_readed, vector_scaling);
			// std::cout<<"get the strip mesh"<<std::endl;
			updateMeshViewer(viewer, updatemesh);

			// std::cout<<"get the strip mesh"<<std::endl;
			meshFileName.push_back("ply_" + meshFileName[id]);
			Meshes.push_back(updatemesh);
			// std::cout<<"before init the poly opt"<<std::endl;
			poly_tool.init(Poly_readed, Bino_readed, -1); // we do not sample the polylines
		}
		ImGui::SameLine();
		if (ImGui::Button("UpdateQuadWithPlylines", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))

		{
			update_qd_mesh_with_plylines();
		}
		if (ImGui::Button("ShadingProjLines", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))

		{
			project_polylines_on_shading_curves_and_save_results();
		}
		ImGui::SameLine();
		if (ImGui::Button("ReadDrawPlyPts", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))

		{
			Eigen::MatrixXd ver;

			read_draw_pts_from_plylines(ver);
			const Eigen::RowVector3d red(0.8, 0.2, 0.2);
			viewer.data().add_points(ver, red); // show rows
		}
		ImGui::SameLine();
		if (ImGui::Button("ShadingMesh", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))

		{
			std::cout << "Please Set Up both vector_scaling and vector_scaling2 for generating the mesh" << std::endl;
			int id = viewer.selected_data_index;

			CGMesh updatemesh;
			read_plylines_extract_offset_mesh(vector_scaling, vector_scaling2, updatemesh);
			updateMeshViewer(viewer, updatemesh);
			meshFileName.push_back("Shading_" + meshFileName[id]);
			Meshes.push_back(updatemesh);

			viewer.selected_data_index = id;
		}

		if (ImGui::CollapsingHeader("Mesh Optimization", ImGuiTreeNodeFlags_CollapsingHeader))
		{
			// Expose variable directly ...
			ImGui::InputInt("Iteration Mesh Opt", &Nbr_Iterations_Mesh_Opt, 0, 0);
			ImGui::InputDouble("Weight Mass", &weight_Mesh_mass, 0, 0, "%.4f");
			ImGui::InputDouble("Weight Smoothness", &weight_Mesh_smoothness, 0, 0, "%.4f");
			ImGui::InputDouble("Weight Mesh Pseudo Geodesic", &weight_Mesh_pesudo_geodesic, 0, 0, "%.4f");
			ImGui::InputDouble("Weight Mesh Egdge Length", &weight_Mesh_edgelength, 0, 0, "%.4f");
			ImGui::InputDouble("Weight Approx", &weight_Mesh_approximation, 0, 0, "%.4f");
			ImGui::InputDouble("Max Mesh Step Length", &Mesh_opt_max_step_length, 0, 0, "%.4f");
			if (ImGui::Button("DrawProjection", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);
				const Eigen::RowVector3d yellow(241. / 255, 196. / 255, 15. / 255);

				const Eigen::RowVector3d black(0, 0, 0);
				const Eigen::RowVector3d green(0.2, 0.8, 0.2);

				if (tools.Ppro0.size() != 0)
				{
					viewer.data().add_edges(tools.Ppro0, tools.V, green);
					// viewer.data().add_edges(tools.Ppro0,
					// 						tools.Ppro0 + vector_scaling * tools.Npro0, red);
					viewer.data().add_edges(tools.Ppro0,
											tools.Ppro0 + vector_scaling * tools.Npro0, red);
				}
			}
			ImGui::SameLine();
			if (ImGui::Button("SmtMesh", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				int id = viewer.selected_data_index;
				CGMesh updatedMesh;
				CGMesh inputMesh = Meshes[id];
				MeshEnergyPrepare initializer;
				initializer.Mesh_opt_max_step_length = Mesh_opt_max_step_length;
				initializer.weight_Mesh_smoothness = weight_Mesh_smoothness;
				initializer.weight_mass = weight_mass;
				initializer.weight_Mesh_edgelength = weight_Mesh_edgelength;
				tools.weight_Mesh_approximation = weight_Mesh_approximation;
				tools.weight_Mesh_mass = weight_Mesh_mass;

				tools.prepare_mesh_optimization_solving(initializer);
				for (int i = 0; i < Nbr_Iterations_Mesh_Opt; i++)
				{
					tools.Run_Mesh_Smoothness();
				}
				updatedMesh = tools.lsmesh;
				updateMeshViewer(viewer, updatedMesh);
				meshFileName.push_back("Smt_" + meshFileName[id]);
				Meshes.push_back(updatedMesh);

				viewer.selected_data_index = id;
				std::cout << "waiting for instructions" << std::endl;
			}
			// ImGui::InputInt("Iteration", &OpIter, 0, 0);
			// ImGui::InputDouble("weight ls mass(big)", &weight_mass, 0, 0, "%.4f");
			// ImGui::Checkbox("Fix Boundary", &fixBoundary_checkbox);
		}
		if (ImGui::CollapsingHeader("Shading Paras", ImGuiTreeNodeFlags_CollapsingHeader))
		{
			ImGui::InputInt("Ver_id", &update_ver_id, 0, 0);
			ImGui::InputInt("Ring nbr", &ring_nbr, 0, 0);
			ImGui::InputDouble("Latitude", &Shading_Latitude, 0, 0, "%.4f");
			if (ImGui::Button("UpdateVertexID", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				update_verlist.push_back(update_ver_id);
				std::cout << "Current ver list" << std::endl;
				for (int i = 0; i < update_verlist.size(); i++)
				{
					std::cout << update_verlist[i] << ",";
				}
				std::cout << "\n";
			}
			ImGui::SameLine();
			if (ImGui::Button("LoadVertexID", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				std::cout << "Read Vertex ID file row by row" << std::endl;
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
				update_verlist = list;
				for (int i = 0; i < update_verlist.size(); i++)
				{
					std::cout << update_verlist[i] << ",";
				}
				std::cout << "\n";
			}
			ImGui::SameLine();
			if (ImGui::Button("ClearVertexList", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				update_verlist.clear();
				std::cout << "ver list got cleared" << std::endl;
			}
			ImGui::SameLine();
			if (ImGui::Button("SaveRotatedMesh", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				std::cout << "Rotate the mesh in the earth axis according to the input latitude.\nNow saving the mesh" << std::endl;
				Eigen::MatrixXd Vnew;
				rotate_z_axis_to_earth_axis(tools.V, Vnew, Shading_Latitude);

				std::string fname = igl::file_dialog_save();
				if (fname.length() == 0)
				{
					std::cout << "\nLSC: save mesh failed" << std::endl;
					ImGui::End();
					return;
				}
				igl::writeOBJ(fname, Vnew, tools.F);
				std::cout << "mesh saved" << std::endl;
			}

			ImGui::Checkbox("Shading", &Given_Const_Direction);
			ImGui::SameLine();
			ImGui::Checkbox("Light Through", &let_ray_through);
			ImGui::SameLine();
			ImGui::Checkbox("Light Reflect", &let_ray_reflect);

			ImGui::PushItemWidth(50);
			ImGui::InputDouble("T0", &InputPx, 0, 0, "%.6f");
			ImGui::SameLine();
			ImGui::InputDouble("P0", &InputPy, 0, 0, "%.6f");
			ImGui::SameLine();
			ImGui::InputDouble("Ttol0", &InputThetaTol, 0, 0, "%.6f");
			ImGui::SameLine();
			ImGui::InputDouble("Ptol0", &InputPhiTol, 0, 0, "%.6f");

			ImGui::InputDouble("T1", &InputPx1, 0, 0, "%.6f");
			ImGui::SameLine();
			ImGui::InputDouble("P1", &InputPy1, 0, 0, "%.6f");
			ImGui::SameLine();
			ImGui::InputDouble("Ttol1", &InputThetaTol1, 0, 0, "%.6f");
			ImGui::SameLine();
			ImGui::InputDouble("Ptol1", &InputPhiTol1, 0, 0, "%.6f");

			ImGui::Checkbox("ShadingInit", &shading_init);
			ImGui::SameLine();
			ImGui::InputDouble("weight binormal", &weight_binormal, 0, 0, "%.4f");
			ImGui::InputDouble("weight smtBinormal", &weight_smt_binormal, 0, 0, "%.4f");
			ImGui::Checkbox("RecomptAuxiliaries", &recompute_auxiliaries);
			if (ImGui::Button("draw slopes", ImVec2(ImGui::GetWindowSize().x * 0.23f, 0.0f)))
			{
				Eigen::MatrixXd E0, E1;
				tools.show_slopes(vector_scaling, E0, E1);
				viewer.data().add_edges(E0, E1, sea_green);
				std::cout << "Drawing the slopes" << std::endl;
			}
		}

		ImGui::PushItemWidth(50);
		if (ImGui::CollapsingHeader("PickOnePolyline", ImGuiTreeNodeFlags_CollapsingHeader))
		{
			ImGui::Checkbox("pickPly", &pick_single_ply);
			ImGui::SameLine();
			ImGui::InputInt("pickID", &pick_line_id, 0, 0);
		}
		if (ImGui::CollapsingHeader("TwoLvlStAngle", ImGuiTreeNodeFlags_CollapsingHeader))
		{
			ImGui::Checkbox("2LvStAngle", &fix_angle_of_two_levelsets);
			ImGui::SameLine();
			ImGui::InputDouble("AngleOfLSs", &angle_between_two_levelsets, 0, 0, "%.4f");
			ImGui::SameLine();
			ImGui::InputDouble("WeightAngleOfLSs", &weight_fix_two_ls_angle, 0, 0, "%.4f");
		}
		if (ImGui::CollapsingHeader("InterActiveDesign", ImGuiTreeNodeFlags_DefaultOpen))
		{
			ImGui::Checkbox("drawStrokes", &draw_strokes);
			ImGui::SameLine();
			if (ImGui::Button("ClearStrokes", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))

			{
				project_2dx.clear();
				project_2dy.clear();
				project_x_tmp.clear();
				project_y_tmp.clear();
				project_pts.clear();
				flist.clear();
				bclist.clear();
				CGMesh updatemesh = tools.lsmesh;
				int id = viewer.selected_data_index;
				updateMeshViewer(viewer, updatemesh);
				meshFileName.push_back(meshFileName[id]);
				Meshes.push_back(updatemesh);
				std::cout << "The Strokes Are Cleared" << std::endl;
			}
			ImGui::SameLine();
			if (ImGui::Button("Strokes2LS", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))

			{
				std::cout << "Using Strokes to Initialize the Level Set" << std::endl;
				tools.receive_interactive_strokes_and_init_ls(flist, bclist);
				bool compute_pg = false;
				AssignAutoRunDefaultArgs(autorunner[0], compute_pg);
				compute_pg = true;
				AssignAutoRunDefaultArgs(autorunner[1], compute_pg);
				std::cout << "Using Default Parameters" << std::endl;
			}

			if (ImGui::Button("AutoRunSmoothLS", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))

			{
				AutoRunArgs runner = autorunner[0];
				for (int i = 0; i < runner.parts.size(); i++)
				{
					tools.weight_mass = runner.parts[i].weight_gravity;
					tools.weight_laplacian = runner.parts[i].weight_lap;
					tools.weight_boundary = runner.parts[i].weight_bnd;
					tools.weight_pseudo_geodesic_energy = runner.parts[i].weight_pg;
					tools.weight_strip_width = runner.parts[i].weight_strip_width;
					for (int j = 0; j < runner.parts[i].iterations; j++)
					{
						tools.Run_Level_Set_Opt_interactive(runner.compute_pg);
						if (tools.step_length < runner.stop_step_length && j != 0)
						{
							std::cout << "optimization converges " << std::endl;
							break;
						}
					}
				}

				std::cout << "A Smooth Scalar Function Following Input Strokes Is Generated ..." << std::endl;
				int id = viewer.selected_data_index;
				CGMesh inputMesh = tools.lsmesh;
				updateMeshViewer(viewer, inputMesh);
				meshFileName.push_back("lso_" + meshFileName[id]);
				Meshes.push_back(inputMesh);

				Eigen::VectorXd level_set_values;
				tools.show_level_set(level_set_values);
				if (level_set_values.size() == 0)
				{
					std::cout << "Level Set is not initialized" << std::endl;
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
			if (ImGui::Button("AutoRunPG", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))

			{
				timer_global.start();
				AutoRunArgs runner = autorunner[1];
				igl::Timer lctimer;
				for (int i = 0; i < runner.parts.size(); i++)
				{
					tools.weight_mass = runner.parts[i].weight_gravity;
					tools.weight_laplacian = runner.parts[i].weight_lap;
					tools.weight_boundary = runner.parts[i].weight_bnd;
					tools.weight_pseudo_geodesic_energy = runner.parts[i].weight_pg;
					tools.weight_strip_width = runner.parts[i].weight_strip_width;
					for (int j = 0; j < runner.parts[i].iterations; j++)
					{
						lctimer.start();
						tools.Run_Level_Set_Opt_interactive(runner.compute_pg);
						lctimer.stop();
						errorlist.push_back(tools.pgerror);
						if(timelist.empty()){
							timelist.push_back(lctimer.getElapsedTimeInSec());
						}
						else{
							timelist.push_back(timelist.back() + lctimer.getElapsedTimeInSec());
						}

						iteration_total++;
						if (
							(tools.step_length < runner.stop_step_length || tools.pgerror < runner.stop_energy_sqrt * runner.stop_energy_sqrt) 
							&& j != 0)
						{
							std::cout << "optimization converges " << std::endl;
							break;
						}
					}
					std::cout << "\n";
				}
				timer_global.stop();
				time_total += timer_global.getElapsedTimeInSec();

				std::cout << "Pseudo-Geodesic Curves in level sets are Generated ..." << std::endl;
				int id = viewer.selected_data_index;
				CGMesh inputMesh = tools.lsmesh;
				updateMeshViewer(viewer, inputMesh);
				meshFileName.push_back("lso_" + meshFileName[id]);
				Meshes.push_back(inputMesh);

				Eigen::VectorXd level_set_values;
				tools.show_level_set(level_set_values);
				if (level_set_values.size() == 0)
				{
					std::cout << "Level Set is not initialized" << std::endl;
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
				save_strokes(flist, bclist, tools.V, tools.F);
			}
			ImGui::SameLine();
			if (ImGui::Button("ReadStrokes", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))

			{
				Eigen::MatrixXd V = tools.V;
				Eigen::MatrixXi F = tools.F;
				read_strokes(flist, bclist);
				Eigen::MatrixXd vlist;
				for (int i = 0; i < flist.size(); i++)
				{
					vlist.resize(flist[i].size(), 3);
					for (int j = 0; j < flist[i].size(); j++)
					{
						int fid = flist[i][j];
						Eigen::Vector3f bc = bclist[i][j];
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
		}
		if (ImGui::CollapsingHeader("QuadMeshOpt", ImGuiTreeNodeFlags_CollapsingHeader))
		{
			// ImGui::InputInt("Ver_id", &update_ver_id, 0, 0);
			// ImGui::InputInt("Ring nbr", &ring_nbr, 0, 0);
			// ImGui::InputDouble("Latitude", &Shading_Latitude, 0, 0, "%.4f");
			ImGui::Combo("QuadType", &quad_tool.OptType,
						 "AAG\0GGA\0PP\0PPG\0\0");
			ImGui::SameLine();

			ImGui::Combo("FamilyOfDiag", &quad_tool.WhichDiagonal,
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
				quad_tool.init(quadmesh, savefile);

				size_t last_dot = fname.rfind('.');
				size_t last_slash = fname.rfind(spliter); // TODO on linux it should be '/'
				std::string fnameNew = fname.substr(last_slash + 1, (last_dot - last_slash - 1));
				meshFileName.push_back(fnameNew);
				Meshes.push_back(quadmesh);

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

				updateMeshViewer(viewer, quadmesh);
			}
			ImGui::SameLine();

			if (ImGui::Button("Set&Reset", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				quad_tool.reset();
				std::array<Eigen::MatrixXd, 3> Edges;
				quad_tool.show_curve_families(Edges);
				int id = viewer.selected_data_index;
				updateMeshViewer(viewer, quad_tool.mesh_original);
				meshFileName.push_back("Rst_" + meshFileName[id]);
				Meshes.push_back(quad_tool.mesh_original);
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);
				const Eigen::RowVector3d black(0, 0, 0);
				const Eigen::RowVector3d green(0.2, 0.8, 0.2);

				viewer.data().add_points(Edges[0], red);   // show rows
				viewer.data().add_points(Edges[1], green); // cols
				viewer.data().add_points(Edges[2], blue);  // diagonals
				// std::cout<<"Visualing Rows Cols and Diagonals : "<<Edges[0]<<",\n\n "<<Edges[1]<<",\n\n "<<Edges[2]<<"\n";
				std::cout << "The rows, columns and the diagonals are marked in red, green and blue" << std::endl;
				;
				viewer.selected_data_index = id;
			}
			ImGui::SameLine();

			if (ImGui::Button("Opt", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				quad_tool.weight_fairness = weight_laplacian;
				quad_tool.weight_gravity = weight_boundary;
				quad_tool.weight_pg = weight_pseudo_geodesic;
				quad_tool.angle_degree0 = target_angle;
				quad_tool.angle_degree1 = target_angle_2;
				quad_tool.pg_ratio = weight_geodesic;
				quad_tool.weight_mass = weight_mass;
				quad_tool.max_step = maximal_step_length;
				if (quad_tool.V.rows() == 0)
				{
					std::cout << "\nEmpty quad, please load a quad mesh first" << std::endl;
					ImGui::End();
					return;
				}
				timer_global.start();
				for (int i = 0; i < OpIter; i++)
				{
					quad_tool.opt();
					iteration_total++;
					if (quad_tool.real_step_length < 1e-16 && i != 0)
					{ // step length actually is the value for the last step
						std::cout << "optimization converges " << std::endl;
						break;
					}
				}
				timer_global.stop();
				time_total += timer_global.getElapsedTimeInSec();

				std::cout << "waiting for instruction..." << std::endl;
				// MP.MeshUnitScale(inputMesh, updatedMesh);
				int id = viewer.selected_data_index;
				CGMesh updateMesh = quad_tool.mesh_update;
				updateMeshViewer(viewer, updateMesh);
				meshFileName.push_back("lso_" + meshFileName[id]);
				Meshes.push_back(updateMesh);
				viewer.selected_data_index = id;
			}
			ImGui::SameLine();

			if (ImGui::Button("DrawNormal", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);
				const Eigen::RowVector3d yellow(241. / 255, 196. / 255, 15. / 255);

				const Eigen::RowVector3d black(0, 0, 0);
				const Eigen::RowVector3d green(0.2, 0.8, 0.2);

				Eigen::MatrixXd pts, E0, E1;

				std::cout << "ploting the normal vectors" << std::endl;
				Eigen::MatrixXd binormals;
				if (quad_tool.Enm0.rows() == 0)
				{
					std::cout << "\nWrong Usage ..." << std::endl;
					ImGui::End();
					return;
				}

				viewer.data().add_edges(quad_tool.Enm0,
										quad_tool.Enm0 + (quad_tool.Enm1 - quad_tool.Enm0) * vector_scaling, green);
			}
			if (ImGui::Button("writePlyInfo", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				quad_tool.write_polyline_info();
			}
			ImGui::SameLine();
			if (ImGui::Button("LoadTriangTree", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				quad_tool.load_triangle_mesh_tree(tools.aabbtree, tools.Vstored,
												  tools.F, tools.Nstored);
			}
		}
		if (ImGui::CollapsingHeader("FunctionalAngles", ImGuiTreeNodeFlags_CollapsingHeader))
		{
			// ImGui::InputInt("Ver_id", &update_ver_id, 0, 0);
			// ImGui::InputInt("Ring nbr", &ring_nbr, 0, 0);
			// ImGui::InputDouble("Latitude", &Shading_Latitude, 0, 0, "%.4f");
			ImGui::Checkbox("DisableIncreasing", &Disable_Changed_Angles);
			ImGui::SameLine();
			ImGui::Checkbox("UseFittingAngles", &Use_Fitting_Angles);
			ImGui::SameLine();
			ImGui::Checkbox("Opt_Only_BNMS", &Use_Opt_Only_BNMS);

			if (ImGui::Button("InitAOAP", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				// binormals as orthogonal as possible.
				EnergyPrepare einit;
				einit.weight_gravity = weight_mass;
				einit.weight_lap = weight_laplacian;
				einit.weight_bnd = weight_boundary;
				einit.weight_pg = weight_pseudo_geodesic;
				einit.weight_strip_width = weight_strip_width;
				einit.solve_pseudo_geodesic = enable_pg_energy_checkbox;
				einit.target_angle = target_angle;
				einit.max_step_length = maximal_step_length;
				einit.solve_strip_width_on_traced = enable_strip_width_checkbox;
				einit.enable_extreme_cases = enable_extreme_cases;
				einit.Given_Const_Direction = Given_Const_Direction;
				// tools.weight_shading = weight_shading;
				tools.prepare_level_set_solving(einit);
				tools.Reference_theta = InputPx;
				tools.Reference_phi = InputPy;
				tools.Reference_theta1 = InputPx1;
				tools.Reference_phi1 = InputPy1;
				tools.Theta_tol = InputThetaTol;
				tools.Phi_tol = InputPhiTol;
				tools.Theta_tol1 = InputThetaTol1;
				tools.Phi_tol1 = InputPhiTol1;
				tools.weight_geodesic = weight_geodesic;
				tools.enable_shading_init = shading_init;
				tools.weight_binormal = weight_binormal;
				tools.weight_smt_binormal = weight_smt_binormal;
				tools.enable_let_ray_through = let_ray_through;
				tools.ShadingLatitude = Shading_Latitude;
				tools.enable_reflection = let_ray_reflect;
				tools.recompute_auxiliaries = recompute_auxiliaries;
				tools.Disable_Changed_Angles = Disable_Changed_Angles;
				tools.Use_Fitting_Angles = Use_Fitting_Angles;
				tools.Use_Opt_Only_BNMS = Use_Opt_Only_BNMS;

				for (int i = 0; i < OpIter; i++)
				{
					tools.RunInitAOAP();
					if (tools.step_length < 1e-16 && i != 0)
					{ // step length actually is the value for the last step
						std::cout << "optimization converges " << std::endl;
						break;
					}
				}
				std::cout << "waiting for instruction..." << std::endl;
				// MP.MeshUnitScale(inputMesh, updatedMesh);
				int id = viewer.selected_data_index;
				CGMesh inputMesh = tools.lsmesh;
				updateMeshViewer(viewer, inputMesh);
				meshFileName.push_back("AOAP_" + meshFileName[id]);
				Meshes.push_back(inputMesh);

				Eigen::VectorXd level_set_values;
				tools.show_level_set(level_set_values);
				if (level_set_values.size() == 0)
				{
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
				// // tools.show_gradients(E0,E1, vector_scaling);
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);

				viewer.selected_data_index = id;
			}
			ImGui::SameLine();
			if (ImGui::Button("OrthoAsPsblOpt", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				// binormals as orthogonal as possible.
				EnergyPrepare einit;
				einit.weight_gravity = weight_mass;
				einit.weight_lap = weight_laplacian;
				einit.weight_bnd = weight_boundary;
				einit.weight_pg = weight_pseudo_geodesic;
				einit.weight_strip_width = weight_strip_width;
				einit.solve_pseudo_geodesic = enable_pg_energy_checkbox;
				einit.target_angle = target_angle;
				einit.max_step_length = maximal_step_length;
				einit.solve_strip_width_on_traced = enable_strip_width_checkbox;
				einit.enable_extreme_cases = enable_extreme_cases;
				einit.Given_Const_Direction = Given_Const_Direction;
				// tools.weight_shading = weight_shading;
				tools.prepare_level_set_solving(einit);
				tools.Reference_theta = InputPx;
				tools.Reference_phi = InputPy;
				tools.Reference_theta1 = InputPx1;
				tools.Reference_phi1 = InputPy1;
				tools.Theta_tol = InputThetaTol;
				tools.Phi_tol = InputPhiTol;
				tools.Theta_tol1 = InputThetaTol1;
				tools.Phi_tol1 = InputPhiTol1;
				tools.weight_geodesic = weight_geodesic;
				tools.enable_shading_init = shading_init;
				tools.weight_binormal = weight_binormal;
				tools.weight_smt_binormal = weight_smt_binormal;
				tools.enable_let_ray_through = let_ray_through;
				tools.ShadingLatitude = Shading_Latitude;
				tools.enable_reflection = let_ray_reflect;
				tools.recompute_auxiliaries = recompute_auxiliaries;
				tools.Disable_Changed_Angles = Disable_Changed_Angles;
				tools.Use_Fitting_Angles = Use_Fitting_Angles;
				tools.Use_Opt_Only_BNMS = Use_Opt_Only_BNMS;

				for (int i = 0; i < OpIter; i++)
				{
					tools.Run_AsOrthAsPossible_LS();
					if (tools.step_length < 1e-16 && i != 0)
					{ // step length actually is the value for the last step
						std::cout << "optimization converges " << std::endl;
						break;
					}
				}
				std::cout << "waiting for instruction..." << std::endl;
				// MP.MeshUnitScale(inputMesh, updatedMesh);
				int id = viewer.selected_data_index;
				CGMesh inputMesh = tools.lsmesh;
				updateMeshViewer(viewer, inputMesh);
				meshFileName.push_back("AOAP_" + meshFileName[id]);
				Meshes.push_back(inputMesh);

				Eigen::VectorXd level_set_values;
				tools.show_level_set(level_set_values);
				if (level_set_values.size() == 0)
				{
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
				// // tools.show_gradients(E0,E1, vector_scaling);
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);

				viewer.selected_data_index = id;
			}
			ImGui::SameLine();
			if (ImGui::Button("WriteFitting", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				tools.write_fitting_data();
			}
			ImGui::SameLine();
			if (ImGui::Button("DrawMinK", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				Eigen::MatrixXd E0, E1;
				tools.show_minimal_curvature_directions(E0, E1, vector_scaling);
				viewer.data().add_edges(E0, E1, hot_red);
			}
		}

		if (ImGui::CollapsingHeader("ConstSlope", ImGuiTreeNodeFlags_CollapsingHeader))
		{
			// ImGui::InputInt("Ver_id", &update_ver_id, 0, 0);
			// ImGui::InputInt("Ring nbr", &ring_nbr, 0, 0);
			// ImGui::InputDouble("Latitude", &Shading_Latitude, 0, 0, "%.4f");
			ImGui::Checkbox("AxisFixed", &AxisFixedForSlopes);
			ImGui::SameLine();
			ImGui::Checkbox("AnglesFixed", &AnglesFixedForSlopes);
			ImGui::PushItemWidth(50);
			ImGui::InputDouble("x", &AxisFixInX, 0, 0, "%.4f");
			ImGui::SameLine();
			ImGui::InputDouble("y", &AxisFixInY, 0, 0, "%.4f");
			ImGui::SameLine();
			ImGui::InputDouble("z", &AxisFixInZ, 0, 0, "%.4f");
			ImGui::SameLine();
			ImGui::InputDouble("Angle", &AngleFixIn, 0, 0, "%.4f");

			if (ImGui::Button("Run", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				int id = viewer.selected_data_index;
				CGMesh inputMesh = tools.lsmesh;
				EnergyPrepare einit;
				einit.weight_gravity = weight_mass;
				einit.weight_lap = weight_laplacian;
				einit.weight_bnd = weight_boundary;
				einit.weight_pg = weight_pseudo_geodesic;
				einit.weight_strip_width = weight_strip_width;
				einit.solve_pseudo_geodesic = enable_pg_energy_checkbox;
				einit.target_angle = target_angle;
				einit.max_step_length = maximal_step_length;
				einit.solve_strip_width_on_traced = enable_strip_width_checkbox;
				einit.enable_extreme_cases = enable_extreme_cases;
				tools.fix_angle_of_two_levelsets = fix_angle_of_two_levelsets;
				tools.angle_between_two_levelsets = angle_between_two_levelsets;
				tools.weight_fix_two_ls_angle = weight_fix_two_ls_angle;
				tools.AxisFixedForSlopes = AxisFixedForSlopes;
				tools.AnglesFixedForSlopes = AnglesFixedForSlopes;
				;
				tools.AxisFixIn = Eigen::Vector3d(AxisFixInX, AxisFixInY, AxisFixInZ);
				tools.AngleFixIn = AngleFixIn;
				// tools.weight_shading = weight_shading;
				//
				tools.prepare_level_set_solving(einit);
				timer_global.start();
				for (int i = 0; i < OpIter; i++)
				{
					tools.Run_ConstSlopeOpt();
					iteration_total++;
					if (tools.step_length < 1e-16 && i != 0)
					{ // step length actually is the value for the last step
						std::cout << "optimization converges " << std::endl;
						break;
					}
				}
				timer_global.stop();
				time_total += timer_global.getElapsedTimeInSec();
				std::cout << "waiting for instruction..." << std::endl;
				updateMeshViewer(viewer, inputMesh);
				meshFileName.push_back("lso_" + meshFileName[id]);
				Meshes.push_back(inputMesh);

				Eigen::VectorXd level_set_values;
				tools.show_level_set(level_set_values);
				if (level_set_values.size() == 0)
				{
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
			if (ImGui::Button("DrawBxis", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);
				const Eigen::RowVector3d yellow(241. / 255, 196. / 255, 15. / 255);

				const Eigen::RowVector3d black(0, 0, 0);
				const Eigen::RowVector3d green(0.2, 0.8, 0.2);

				Eigen::MatrixXd E0, E1;
				std::cout << "ploting the bxis" << std::endl;
				Eigen::MatrixXd binormals;
				tools.show_bxis(E0, E1);
				viewer.data().add_edges(E0, E1, green);
			}
		}

		if (ImGui::CollapsingHeader("TorsionFreeNodes", ImGuiTreeNodeFlags_DefaultOpen))
		{
			// ImGui::InputInt("Ver_id", &update_ver_id, 0, 0);
			// ImGui::InputInt("Ring nbr", &ring_nbr, 0, 0);
			// ImGui::InputDouble("Latitude", &Shading_Latitude, 0, 0, "%.4f");
			ImGui::Checkbox("FixCentroid", &CentroidFixedForNodes);
			ImGui::SameLine();
			ImGui::Checkbox("fixPly", &RadiusFixedForNodes);
			ImGui::SameLine();
			ImGui::Checkbox("FindC&R", &FindCandRForNodes);
			ImGui::SameLine();
			ImGui::Checkbox("WriteLs&Radius", &tools.write_ls_radius_file);
			ImGui::PushItemWidth(50);
			ImGui::InputDouble("Cx", &CentroidFixInX, 0, 0, "%.4f");
			ImGui::SameLine();
			ImGui::InputDouble("Cy", &CentroidFixInY, 0, 0, "%.4f");
			ImGui::SameLine();
			ImGui::InputDouble("Cz", &CentroidFixInZ, 0, 0, "%.4f");
			ImGui::SameLine();
			ImGui::InputDouble("Radius", &RadiusFixIn, 0, 0, "%.4f");
			if (ImGui::Button("RunRuling", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				// std::cout<<"in run"<<std::endl;
				int id = viewer.selected_data_index;
				CGMesh inputMesh = tools.lsmesh;
				EnergyPrepare einit;
				einit.weight_gravity = weight_mass;
				einit.weight_lap = weight_laplacian;
				einit.weight_bnd = weight_boundary;
				einit.weight_pg = weight_pseudo_geodesic;
				einit.weight_strip_width = weight_strip_width;
				einit.solve_pseudo_geodesic = enable_pg_energy_checkbox;
				einit.target_angle = target_angle;
				einit.max_step_length = maximal_step_length;
				einit.solve_strip_width_on_traced = enable_strip_width_checkbox;
				einit.enable_extreme_cases = enable_extreme_cases;
				tools.fix_angle_of_two_levelsets = fix_angle_of_two_levelsets;
				tools.angle_between_two_levelsets = angle_between_two_levelsets;
				tools.weight_fix_two_ls_angle = weight_fix_two_ls_angle;

				tools.RulingFixCentroid = CentroidFixedForNodes;
				tools.RulingFixRadius = RadiusFixedForNodes;
				tools.RulingOnlyFindingSphere = FindCandRForNodes;
				tools.RulingCentroid = Eigen::Vector3d(CentroidFixInX, CentroidFixInY, CentroidFixInZ);
				// tools.RulingRadius = RadiusFixIn;
				// tools.weight_shading = weight_shading;
				//
				tools.prepare_level_set_solving(einit);
				timer_global.start();
				for (int i = 0; i < OpIter; i++)
				{
					tools.Run_ruling_opt();
					iteration_total++;
					if (tools.step_length < 1e-16 && i != 0)
					{ // step length actually is the value for the last step
						std::cout << "optimization converges " << std::endl;
						break;
					}
				}
				timer_global.stop();
				time_total += timer_global.getElapsedTimeInSec();
				std::cout << "waiting for instruction..." << std::endl;
				updateMeshViewer(viewer, inputMesh);
				meshFileName.push_back("lso_" + meshFileName[id]);
				Meshes.push_back(inputMesh);

				Eigen::VectorXd level_set_values;
				tools.show_level_set(level_set_values);
				if (level_set_values.size() == 0)
				{
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
			
			if (ImGui::Button("RulingCheckRadius", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				tools.RulingCentroid = Eigen::Vector3d(CentroidFixInX, CentroidFixInY, CentroidFixInZ);
				if(tools.fvalues.size()>0){
					tools.statistic_tangent_vector_to_pt_distance();
				}
				else{
					std::cout<<"please load level sets"<<std::endl;
				}
				
			}
			ImGui::SameLine();
			
			if (ImGui::Button("DrawTangents", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				
				const Eigen::RowVector3d red(0.8, 0.2, 0.2);
				const Eigen::RowVector3d blue(0.2, 0.2, 0.8);
				const Eigen::RowVector3d yellow(241. / 255, 196. / 255, 15. / 255);

				const Eigen::RowVector3d black(0, 0, 0);
				const Eigen::RowVector3d green(0.2, 0.8, 0.2);

				viewer.data().add_edges(tools.TangentE0, tools.TangentE0 + (tools.TangentE1 - tools.TangentE0) * vector_scaling,
										green);
			}
			ImGui::SameLine();
			
			if (ImGui::Button("TestQuadratics", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				tools.test_quadratic_compare_0(CentroidFixInX, CentroidFixInY, CentroidFixInZ, dbg_dbl, dbg_dbl2);
			}
			if (ImGui::Button("Plot2dFitting", ImVec2(ImGui::GetWindowSize().x * 0.25f, 0.0f)))
			{
				//tools.test_quadratic_compare_0(CentroidFixInX, CentroidFixInY, CentroidFixInZ, dbg_dbl, dbg_dbl2);
				if (tools.NodeFitPly.size() == 0)
				{
					ImGui::End();
					return;
				}
				int id = viewer.selected_data_index;
				plot2dFittingResults(viewer, tools.XsNode, tools.YsNode, tools.NodeFitPly);
				meshFileName.push_back("fitting_" + meshFileName[id]);
				CGMesh meshtmp;
				Meshes.push_back(meshtmp);
				viewer.selected_data_index = id;
			}


		}
		// binormals as orthogonal as possible.
		ImGui::PushStyleColor(ImGuiCol_Button, (ImVec4)ImColor::HSV(0.0f, 0.6f, 0.6f));
		ImGui::PushStyleColor(ImGuiCol_ButtonHovered, (ImVec4)ImColor::HSV(0.0f, 0.7f, 0.7f));
		ImGui::PushStyleColor(ImGuiCol_ButtonActive, (ImVec4)ImColor::HSV(0.0f, 0.8f, 0.8f));
		ImGui::PopStyleColor(3);
		ImGui::End();
	};
}
