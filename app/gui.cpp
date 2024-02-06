#include <gui.h>

lscif::lscif()
{
	sea_green = Eigen::RowVector3d(70. / 255., 252. / 255., 167. / 255.);
	hot_red = Eigen::RowVector3d(255. / 255., 12. / 255., 17. / 255.);
	pure_blk = Eigen::RowVector3d(0, 0, 0);
}

// add a new mesh into the mesh lists, and show the new mesh along with previous showed meshes
void lscif::updateMeshViewer(igl::opengl::glfw::Viewer &viewer, CGMesh &mesh)
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
void lscif::updatePointsViewer(igl::opengl::glfw::Viewer &viewer, const Eigen::MatrixXd& V)
{
	Eigen::MatrixXi F;
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
	E.resize(V.rows() - 1, 2);
	for (int i = 0; i < E.rows(); i++)
	{
		E(i, 0) = i;
		E(i, 1) = i + 1;
	}
	Eigen::MatrixXd C = Eigen::MatrixXd::Zero(E.rows(), 3);
	viewer.data().set_edges(V, E, C);

	viewer.data().show_lines = false;
	viewer.data().line_width = 2.0;
	viewer.core().align_camera_center(viewer.data().V, viewer.data().F);

	VertexType.resize(V.size(), 0);
}

void lscif::plot2dFittingResults(igl::opengl::glfw::Viewer &viewer, const Eigen::VectorXd &xs, const Eigen::VectorXd &ys,
								 const std::vector<double> &ply)
{
	// make one triangle as the mesh
	Eigen::MatrixXd Vt(3,3);
	Eigen::MatrixXi Ft(1,3);
	double xmin = xs.minCoeff();
	double xmax = xs.maxCoeff();
	double ymin = ys.minCoeff();
	double ymax = ys.maxCoeff();
	Vt << xmin, ymin, 0,
		xmin, ymax, 0,
		xmax, ymin, 0;
	Ft << 0, 1, 2;
	viewer.append_mesh();
	viewer.data().set_mesh(Vt, Ft);
	viewer.data().is_visible = true;

	viewer.data().show_lines = false;
	viewer.data().line_width = 2.0;
	viewer.data().point_size = 5;
	viewer.core().align_camera_center(viewer.data().V, viewer.data().F);

	int vnbr = xs.size();
	Eigen::MatrixXd pts(vnbr, 3), cur(100, 3);
	pts.col(0) = xs;
	pts.col(1) = ys;
	pts.col(2) = Eigen::VectorXd::Zero(vnbr);
	double itv = (xmax - xmin) / 100;
	for (int i = 0; i < 100; i++)
	{
		double x = xmin + i * itv;
		double y = polynomial_value(ply, x);
		// y = ply[3] * x * x *x + ply[]
		cur.row(i) << x, y, 0;
	}
	Eigen::MatrixXd E0(99, 3);
	E0 = cur.block(0, 0, 99, 3);
	Eigen::MatrixXd E1(99, 3);
	E1 = cur.block(1, 0, 99, 3);

	viewer.data().add_points(pts, hot_red);
	viewer.data().add_edges(E0, E1, sea_green);
}

bool lscif::pre_draw(igl::opengl::glfw::Viewer &viewer)
{
	return false;
}
bool lscif::key_up(igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods)
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
bool lscif::key_down(igl::opengl::glfw::Viewer &viewer, unsigned char key, int mods)
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
			if (Meshes.size() == 1)
			{
				std::cout << "Do not delete the last mesh" << std::endl;
				// ImGui::End();
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

		return true;
	}
	}

	return false;
}

// when holding "1" or "2", with mouse left button down, points can be choosen
// TODO fix the bug when select one mesh but click on the other one, and when swich to
// the other mesh, all the histories are recorded and reproduce points on the second mesh
bool lscif::mouse_down(igl::opengl::glfw::Viewer &viewer, int button, int modifier)
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
			std::cout << "vid, " << vid << ", fid, " << fid << ", mass_unfiform, " << tools.mass_uniform.coeffRef(vid, vid) << std::endl;
			if (tools.fvalues.rows() > 0)
			{
				std::cout << "level set value " << tools.fvalues[vid] << std::endl;
			}
			Eigen::VectorXd energy;
			tools.show_max_pg_energy(energy);
			if (energy.size() > 0)
			{
				std::cout << "local energy, " << energy[vid] << std::endl;
			}
			Eigen::MatrixXd energy_all;
			tools.show_max_pg_energy_all(energy_all);
			if (energy_all.rows() == viewer.data().V.rows())
			{
				std::cout << energy_all.row(vid) << std::endl;
			}
			else
			{
				std::cout << "the size of energy_all, " << energy_all.rows() << std::endl;
			}

			std::cout << "point position: (" << viewer.data().V(vid, 0) << ", " << viewer.data().V(vid, 1) << ", " << viewer.data().V(vid, 2) << ")\n\n";
			// tools.print_info(vid);
			if (Given_Const_Direction)
			{
				double latitude_radian = Shading_Latitude * LSC_PI / 180.;
				double rho = (LSC_PI / 2 - latitude_radian);
				Eigen::Matrix3d rotation;
				rotation << 1, 0, 0,
					0, cos(rho), -sin(rho),
					0, sin(rho), cos(rho);
				Eigen::Vector3d direction_ground = Eigen::Vector3d(0, 0, -1); // the ground direction
				direction_ground = rotation * direction_ground;
				Eigen::MatrixXd ground = Eigen::MatrixXd::Ones(tools.V.rows(), tools.V.cols());
				ground.col(0) *= direction_ground[0];
				ground.col(1) *= direction_ground[1];
				ground.col(2) *= direction_ground[2];
				viewer.data().add_edges(tools.V - vector_scaling * ground, tools.V + vector_scaling * ground, pure_blk);
			}

			// viewer.data().add_points(igl::slice(viewer.data().V, vids, 1), hot_red);
			// viewer.data().set_points(igl::slice(viewer.data().V, vids, 1), hot_red);
			viewer.data().add_points(igl::slice(viewer.data().V, vids, 1), hot_red);
			Eigen::MatrixXd pts;
			tools.show_current_reference_points(pts);
			viewer.data().add_points(pts, sea_green);
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
			if (tools.fvalues.rows() > 0)
			{
				std::cout << "level set value " << tools.fvalues[vid] << std::endl;
			}
			std::cout << "point position: (" << viewer.data().V(vid, 0) << ", " << viewer.data().V(vid, 1) << ", " << viewer.data().V(vid, 2) << ")\n\n";
			viewer.data().set_points(igl::slice(viewer.data().V, vids, 1), hot_red);
			viewer.data().add_points(igl::slice(viewer.data().V, vids2, 1), sea_green);
			return true;
		}
	}

	left_button_down = true;
	return false;
}
bool lscif::mouse_up(igl::opengl::glfw::Viewer &viewer, int button, int modifier)
{
	left_button_down = false;
	if (!project_x_tmp.empty()) // if just draw a curve
	{
		project_2dx.push_back(project_x_tmp);
		project_2dy.push_back(project_y_tmp);
		Eigen::MatrixXd Vers;
		draw_stroke_on_mesh(tools.lsmesh, viewer,
							tools.V, tools.F, project_x_tmp,
							project_y_tmp, Vers, flist_tmp, bclist_tmp);
		flist.push_back(flist_tmp);
		bclist.push_back(bclist_tmp);
		project_x_tmp.clear();
		project_y_tmp.clear();
		std::cout << "There are " << project_2dx.size() << " curves" << std::endl;
		viewer.data().add_points(Vers, hot_red);
	}

	return true;
	return false;
}
bool lscif::mouse_move(igl::opengl::glfw::Viewer &viewer, int mouse_x, int mouse_y)
{
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
		project_x_tmp.push_back(x);
		project_y_tmp.push_back(y);

		// sleep(0.05);

		return true;
	}
	return false;
}

void lscif::show_axis(Eigen::MatrixXd &E0, Eigen::MatrixXd &E1, double ratio)
{
	E0 = Eigen::MatrixXd::Zero(3, 3);
	E1 = ratio * Eigen::MatrixXd::Identity(3, 3);
}

void lscif::viewer_launch_functions(igl::opengl::glfw::Viewer &viewer, igl::opengl::glfw::imgui::ImGuiPlugin &plugin,
									igl::opengl::glfw::imgui::ImGuiMenu &menu, igl::opengl::glfw::imgui::ImGuiMenu &menu_mp)
{
	std::map<int, Eigen::RowVector3d> colors;
	int last_selected = -1;
	std::string example_root_path(CHECKER_BOARD_DATA_DIR);

	std::string modelname = "Hall.obj";
	std::string fname = example_root_path + std::string(modelname);
	OpenMesh::IO::read_mesh(mesh, fname);
	tools.init(mesh);
	MP.mesh2Matrix(mesh, V, F);
	meshFileName.push_back(modelname);
	Meshes.push_back(mesh);

	// Init the viewer

	viewer.plugins.push_back(&plugin);
	// Attach a menu plugin
	plugin.widgets.push_back(&menu);

	// Add content to the default menu window
	menu.callback_draw_viewer_menu = [&]()
	{
		// Draw parent menu content
		menu.draw_viewer_menu();
	};

	draw_menu1(viewer, plugin, menu);
	plugin.widgets.push_back(&menu_mp);

	// Add content to the default menu window
	menu_mp.callback_draw_viewer_menu = [&]()
	{
		// Draw parent menu content
		menu_mp.draw_viewer_menu();
	};

	// Draw additional windows
	draw_menu2(viewer, plugin, menu_mp);
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
}
void launch_viewer()
{
	igl::opengl::glfw::Viewer viewer;
	igl::opengl::glfw::imgui::ImGuiPlugin plugin;
	igl::opengl::glfw::imgui::ImGuiMenu menu, menu1;
	lscif viewerTools;
	viewerTools.viewer_launch_functions(viewer, plugin, menu, menu1);

	using namespace std::placeholders;
	auto func_predraw = std::bind(&lscif::pre_draw, &viewerTools, _1);
	auto func_key_down = std::bind(&lscif::key_down, &viewerTools, _1, _2, _3);
	auto func_mouse_down = std::bind(&lscif::mouse_down, &viewerTools, _1, _2, _3);
	auto func_key_up = std::bind(&lscif::key_up, &viewerTools, _1, _2, _3);
	auto func_mouse_move = std::bind(&lscif::mouse_move, &viewerTools, _1, _2, _3);
	auto func_mouse_up = std::bind(&lscif::mouse_up, &viewerTools, _1, _2, _3);

	viewer.callback_pre_draw = func_predraw; //&pre_draw;
	viewer.callback_key_down = func_key_down;
	viewer.callback_mouse_down = func_mouse_down;
	viewer.callback_key_up = func_key_up;
	viewer.callback_mouse_move = func_mouse_move;
	viewer.callback_mouse_up = func_mouse_up;

	viewer.core().is_animating = false;
	viewer.core().animation_max_fps = 30.;
	viewer.core().background_color = Eigen::Vector4f(1, 1, 1, 100);
	viewer.launch();

}