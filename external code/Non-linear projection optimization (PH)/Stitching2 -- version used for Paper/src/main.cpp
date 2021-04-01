#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/writeOFF.h>
#include <igl/facet_components.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/remove_unreferenced.h>
#include <igl/edges.h>
#include <igl/upsample.h>
#include <igl/PI.h>

#include <igl/gaussian_curvature.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/boundary_loop.h>
#include <igl/hausdorff.h>

#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include "DogSet.hpp"
#include "ProcessSeams.hpp"
#include "objectives.hpp"

#include <fstream>

// Main parameters to set:
//-------------------------


bool upsample       = false;    // upsample?
double gcSmooth     = 0.2;      // smoothness for graph cut

// weights for optimization:
double wSeam        = 10;   // seam smoothness
double wSmooth      = 5;    // overall smoothness
double wPair        = 5;    // patch boundaries should match
double wProx        = 5;    // keep close to the dogs
double wDev         = 100;  // minimize angle defect (progessively increased 1000 * iter^3)
double wLevMa       = 5.0;  // Levenberg-Marquardt regularization (sometimes prevents getting stuck)
double wSeamProx    = .5;   // the mesh boundary should stay close to its initial projection on the dogs

int iterations_elapsed = 0;
double max_defect = 0.001;
bool progressiveDevelopabilityIncrease = false;

double DogSet::weight_normals = 0.0;
bool DogSet::use_distance_smoothness = true;

Eigen::MatrixXd V;
Eigen::MatrixXi F;
Eigen::MatrixXd V0;
Eigen::MatrixXi F0;

Eigen::MatrixXd X;
Eigen::MatrixXd PV0;
OptimizationObjectives* optimization = NULL;
ProcessSeams* seams = NULL;


//double average_mass = 0;
//double max_defect = 0;

std::vector<int> boundary_loop_flat(const Eigen::MatrixXi& F)
{
    std::vector<std::vector<int>> boundary_indices;
    igl::boundary_loop(F, boundary_indices);

    std::vector<int> flat_boundary_indices;
    for (std::vector<int> loop : boundary_indices)
        flat_boundary_indices.insert(flat_boundary_indices.end(), loop.begin(), loop.end());

    //for (vector<int> loop : global_boundary_indices)
    //	log_list(5, loop, "loop: ", false);

    return flat_boundary_indices;
}

void get_gaussian_curvature(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F, Eigen::VectorXd& out_K, bool do_exclude_boundary)
{
    //get gaussian curvature
    igl::gaussian_curvature(V, F, out_K); // Compute integral of Gaussian curvature

    // Compute mass matrix
    Eigen::SparseMatrix<double> M, Minv;
    igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);
    igl::invert_diag(M, Minv);
    
    out_K = (Minv * out_K).eval(); // Divide by area to get integral average
    if (!do_exclude_boundary)
        return;

    //exclude boundary
    std::vector<int> boundary = boundary_loop_flat(F);
    for (int vi : boundary)
        out_K(vi) = 0;
}

double get_BBDiagonal(const Eigen::MatrixXd& V)
{
    const Eigen::RowVectorXd min_point = V.colwise().minCoeff();
    const Eigen::RowVectorXd max_point = V.colwise().maxCoeff();
    const Eigen::RowVectorXd bounding_box = max_point - min_point;

    return bounding_box.norm();
}

double get_hausdorff_distance(const Eigen::MatrixXd& V0, const Eigen::MatrixXi& F0, const Eigen::MatrixXd& V1, const Eigen::MatrixXi& F1)
{
    double d0, d1;
    igl::hausdorff(V0, F0, V1, F1, d0);
    igl::hausdorff(V1, F1, V0, F0, d1);

    double bb0 = get_BBDiagonal(V0);
    double bb1 = get_BBDiagonal(V1);

    double d0_percent = d0 / bb0;
    double d1_percent = d1 / bb1;

    return std::max(d0_percent, d1_percent);
}


void initialize_optimization(Eigen::MatrixXd& V, Eigen::MatrixXi& F, DogSet& dogset)
{
	if (seams != NULL)
		delete seams;
	seams = new ProcessSeams(V, F, dogset.I);

	// find smooth seam lines and split parts accordingly
	seams->splitMeshes();

	// optimization
	auto pairs = seams->getPairs(); //stores pairs of corresponding vertices
	//Eigen::MatrixXd X = seams.V;
	X = seams->V;

	//Eigen::MatrixXd PV0;
	dogset.projectPoints(X, seams->L, PV0);

	//output number of patches (info)
	Eigen::VectorXi pieces;
	igl::facet_components(seams->F, pieces);
	std::cout << "(info) labelling: found " << pieces.maxCoeff() + 1 << " patches" << std::endl;


	if (optimization != NULL)
		delete optimization;
	optimization = new OptimizationObjectives(seams->V, seams->L, seams->F, seams->I, dogset, pairs);
}

void run_optimization(Eigen::MatrixXd& V, Eigen::MatrixXi& F, DogSet& dogset)
{
	optimization->wSeamProx = wSeamProx;
	optimization->wSeamSmooth = wSeam;
	optimization->wSmooth = wSmooth;
	optimization->wPair = wPair;
	optimization->wProx = wProx;
	optimization->wDev = wDev;


	for (int k = 0; k < 100; ++k)
	{
		if (k > 2)
		{
			optimization->wDev = wDev = 10 * k * k * k;
		}

		const int steps = optimization->optimize(X, 1.0, wLevMa);

		Eigen::VectorXd ad;
		optimization->angleDefect(X, ad);

		// take the mean only over non zero values to exclude boundary vertices
		double mean = .0;
		int cnt2 = 0;
		for (int i = 0; i < ad.size(); ++i)
		{
			if (ad(i))
			{
				mean += std::abs(ad(i));
				++cnt2;
			}
		}

		mean /= cnt2;

		std::cout << " gauss: " << ad.norm() << " mean: " << mean << " max: " << ad.cwiseAbs().maxCoeff() << " steps: " << steps << "\n";
		// igl::writeOBJ("../out.obj", X, seams.F);
		if (ad.cwiseAbs().maxCoeff() < 0.005) break;
	}

	igl::writeOBJ("out.obj", X, seams->F);
}

void update_label_view(igl::opengl::glfw::Viewer& viewer, int& idProjection, int& idMinProjection, Eigen::MatrixXd& projColors)
{
    viewer.data(idProjection).clear();
    viewer.data(idProjection).set_mesh(X, seams->F);
    igl::colormap(igl::COLOR_MAP_TYPE_JET, seams->L, 0.0, seams->L.maxCoeff(), projColors);
    viewer.data(idProjection).set_colors(projColors);

    viewer.data(idMinProjection).clear();
    viewer.data(idMinProjection).set_mesh(PV0, seams->F);
    viewer.data(idMinProjection).set_colors(projColors);

    //output number of patches (info)
    Eigen::VectorXi C;
    igl::facet_components(seams->F, C);
    std::cout << "(info) labelling: found " << C.maxCoeff() + 1 << " patches" << std::endl;
}

void export_mesh_data(std::string suffix)
{
	// --- Stitch mesh to show curvature on seams ---
	Eigen::MatrixXd stitched_V;
	Eigen::MatrixXi stitched_F;
	seams->stitchMesh(X, stitched_V, stitched_F);

	igl::writeOBJ("../out_stitched"+suffix+".obj", stitched_V, stitched_F);
	igl::writeOBJ("../out"+suffix+".obj", X, seams->F);


	// --- Compute result metrics and write to file ---
	Eigen::VectorXd K;
	get_gaussian_curvature(X, seams->F, K, true);
	std::cout << "\n(info) EXPORTING mesh with normalized |K|max = " << K.cwiseAbs().maxCoeff() << ", Kmean = " << K.mean() << "\n" << std::endl;

	//hausdorff distance: original -> developable
	double d0;
	igl::hausdorff(V0, F0, X, seams->F, d0);
	double bb0 = get_BBDiagonal(V0);

	//hausdorff distance: developable -> original
	double d1;
	igl::hausdorff(X, seams->F, V0, F0, d1);
	double bb1 = get_BBDiagonal(X);

	std::ofstream file_metrics("../out_metrics"+suffix+".txt");
		file_metrics << "|Kmax | = " << K.cwiseAbs().maxCoeff() << "\n";
		file_metrics << " Kmean  = " << K.mean() << "\n";
		file_metrics << " Kmin   = " << K.minCoeff() << "\n";
		file_metrics << " Kmax   = " << K.maxCoeff() << "\n";

		file_metrics << "\nhausdorff distance = " << std::max(d0 / bb0, d1 / bb1) << "\n";
		file_metrics << "  original -> developable " << "\n";
		file_metrics << "    absolute = " << d0 << "\n";
		file_metrics << "    relative = " << d0 / bb0 << "(wrt diagonal = " << bb0 << ")\n";
		file_metrics << "  developable -> original " << "\n";
		file_metrics << "    absolute = " << d1 << "\n";
		file_metrics << "    relative = " << d1 / bb1 << "(wrt diagonal = " << bb1 << ")\n";
	file_metrics.close();

}

int main(const int argc, const char* argv[])
{
    //std::string dogFile = "../../data/new/side11_40/dog.obj";
    //std::string inputFile = "../../data/new/side11_40/input.obj";
    std::string dogFile =   "../data/dogs.obj";
    std::string inputFile = "../data/input.obj";
    

    if(argc == 3)
    {
        dogFile = argv[1];
        inputFile = argv[2];
    }
    
    //Eigen::MatrixXd V;
    //Eigen::MatrixXi F;
    //Eigen::MatrixXd V0;
    //Eigen::MatrixXi F0;
     
    igl::readOBJ(inputFile, V0, F0);
    V = V0;
    F = F0;

 //   if(upsample)
 //   {
 //       igl::upsample(V0, F0, V, F, 1);
 //   } 
	//else
 //   {
 //       V = V0;
 //       F = F0;
 //   }
        

	//Eigen::SparseMatrix<double> M;
	//igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_VORONOI, M);

	//Eigen::VectorXd ones(V.rows());
	//ones.setOnes();
	//Eigen::VectorXd V_mass = M * ones;
	//average_mass = V_mass.mean();
	//max_defect = (0.5 / 180 * igl::PI) * average_mass;
	//std::cout << "Mass mean = " << average_mass << " --> threshold = " << max_defect << std::endl;



    Eigen::MatrixXd Vdog;
    Eigen::MatrixXi Fdog;
  
    igl::readOBJ(dogFile, Vdog, Fdog);
    DogSet dogset(Vdog, Fdog);
   
    // assign dogs
    dogset.graphCutAssignement(V, F, gcSmooth);
    
    {
        // how many different labels are used
        std::vector<int> labels(dogset.I.data(), dogset.I.data() + dogset.I.size());
        std::sort(labels.begin(), labels.end());
        labels.erase(std::unique(labels.begin(), labels.end()), labels.end());
        std::cout << "#labels: " << labels.size() << std::endl;
    }

    
	initialize_optimization(V, F, dogset);
	//run_optimization(V, F, dogset);

	/*
    ProcessSeams seams(V, F, dogset.I);
   
    // find smooth seam lines and split parts accordingly
    seams.splitMeshes();
    
    // optimization
    auto pairs = seams.getPairs(); //stores pairs of corresponding vertices
    Eigen::MatrixXd X = seams.V;

    Eigen::MatrixXd PV0;
    dogset.projectPoints(X, seams.L, PV0);
    
    OptimizationObjectives optimization(seams.V, seams.L, seams.F, seams.I, dogset, pairs);
    
    
    optimization.wSeamProx = wSeamProx;
    optimization.wSeamSmooth = wSeam;
    optimization.wSmooth = wSmooth;
    optimization.wPair = wPair;
    optimization.wProx = wProx;
    optimization.wDev = wDev;
    
    
    for(int k = 0; k < 100; ++k)
    {
        if(k > 2)
        {
            optimization.wDev = wDev = 10 * k * k * k;
        }
        
        const int steps = optimization.optimize(X, 1.0, wLevMa);
    
        Eigen::VectorXd ad;
        optimization.angleDefect(X, ad);
                
        // take the mean only over non zero values to exclude boundary vertices
        double mean = .0;
        int cnt2 = 0;
        for(int i = 0; i < ad.size(); ++i)
        {
            if(ad(i))
            {
                mean += std::abs(ad(i));
                ++cnt2;
            }
        }
        
        mean /= cnt2;
           
        std::cout << " gauss: " << ad.norm() << " mean: " << mean << " max: " <<  ad.cwiseAbs().maxCoeff() << " steps: " << steps<< "\n";
         // igl::writeOBJ("../out.obj", X, seams.F);
        if(ad.cwiseAbs().maxCoeff() < 0.005) break;
    }
    
    igl::writeOBJ("out.obj", X, seams.F);
    */

#if 1 // 0 for disable GUI

    Eigen::MatrixXd dogColors;
    
    Eigen::VectorXi C;
    igl::facet_components(dogset.Fall, C);
    igl::colormap(igl::COLOR_MAP_TYPE_JET, C, 0.0, C.maxCoeff(), dogColors);
    

    Eigen::MatrixXd origColors;
    igl::colormap(igl::COLOR_MAP_TYPE_JET, seams->L, 0.0, seams->L.maxCoeff(), origColors);
        
    Eigen::MatrixXd PV;
    dogset.projectPoints(seams->V, seams->L, PV);
    
    Eigen::MatrixXd projColors;
    igl::colormap(igl::COLOR_MAP_TYPE_JET, seams->L, 0.0, seams->L.maxCoeff(), projColors);


    igl::opengl::glfw::Viewer viewer;

    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);

    
    bool showOriginal = false;
    bool showDogs = false;
    bool showProjection = true;
    bool showMinProjection = false;
    
    int idOriginal, idDogs, idProjection, idMinProjection;

 
    double t0 = 1.;
    int iterations = 1;
	double smallest_maxK = 1e6;

    // Add content to the default menu window
    menu.callback_draw_viewer_menu = [&]()
    {
        // Draw parent menu content
        //   menu.draw_viewer_menu();
        
        // Add new group
        if (ImGui::CollapsingHeader("View", ImGuiTreeNodeFlags_DefaultOpen))
        {
            
            if(ImGui::Checkbox("Show original", &showOriginal))
            {
                viewer.data(idOriginal).set_visible(showOriginal);
                viewer.selected_data_index = idOriginal;
            }
            
            if(ImGui::Checkbox("Show DOG", &showDogs))
            {
                viewer.data(idDogs).set_visible(showDogs);
                viewer.selected_data_index = idDogs;
            }
            
            if(ImGui::Checkbox("Show projection", &showProjection))
            {
                viewer.data(idProjection).set_visible(showProjection);
                viewer.selected_data_index = idProjection;
            }
            
            if(ImGui::Checkbox("Show naive projection", &showMinProjection))
            {
                viewer.data(idMinProjection).set_visible(showMinProjection);
                viewer.selected_data_index = idMinProjection;
            }
            //if (ImGui::Button("upsample"))
            //{
            //    igl::upsample(V, F, V, F, 1);

            //    dogset.graphCutAssignement(V, F, gcSmooth);
            //    initialize_optimization(V, F, dogset);

            //    update_label_view(viewer, idProjection, idMinProjection, projColors);
            //}
            //if(ImGui::Button("smooth labels"))
            //{
            //    X = PV;
            //    viewer.data(idProjection).set_mesh(X, seams->F);
            //}
        }
        
        if (ImGui::CollapsingHeader("Optimization Weights", ImGuiTreeNodeFlags_DefaultOpen))
        {
			ImGui::InputDouble("Label Smooth", &gcSmooth);
			ImGui::Checkbox("use distance weighted smoothness", &DogSet::use_distance_smoothness);
			ImGui::InputDouble("weight normals", &DogSet::weight_normals);
			if (ImGui::Button("update labels"))
			{
				dogset.graphCutAssignement(V, F, gcSmooth);
				initialize_optimization(V, F, dogset);
                update_label_view(viewer, idProjection, idMinProjection, projColors);
			}

            ImGui::InputDouble("Seam", &wSeam);
            ImGui::InputDouble("Smooth", &wSmooth);
            ImGui::InputDouble("Pair", &wPair);
            ImGui::InputDouble("Proximity", &wProx);
            ImGui::InputDouble("Angles", &wDev);
            ImGui::InputDouble("SeamProx", &wSeamProx);
            
            if (ImGui::CollapsingHeader("Optimization Settings", ImGuiTreeNodeFlags_DefaultOpen))
            {
                ImGui::InputDouble("Levenberg-Marquardt", &wLevMa);
                ImGui::InputDouble("initial timestep", &t0);
				ImGui::Checkbox("increase Developability progressively", &progressiveDevelopabilityIncrease);
                ImGui::InputInt("Iterations", &iterations);
            }
            
            if(ImGui::Button("reset"))
            {
                X = seams->V;
                viewer.data(idProjection).set_mesh(X, seams->F);
                
                Eigen::MatrixXd N;
                igl::per_vertex_normals(X, seams->F, N);
                viewer.data(idProjection).set_normals(N);
            
                viewer.data(idProjection).set_colors(projColors);

                iterations_elapsed = 0;
            }
            
            if(ImGui::Button("recompute"))
            {
              //  OptimizationObjectives optimization(seams.V, seams.L, seams.F, dogset, pairs);
                optimization->wSeamSmooth = wSeam;
                optimization->wSmooth = wSmooth;
                optimization->wPair = wPair;
                optimization->wProx = wProx;
                optimization->wDev = wDev;
                optimization->wSeamProx = wSeamProx;
                
                
                std::cout << "start optimization" << std::endl;
                std::cout << "initial energy: " << optimization->objective(X) << "\n";
                //std::cout << "target defect:  " << max_defect << "\n";
                
                for(int k = 0; k < iterations; ++k)
                {
                    const int steps = optimization->optimize(X, t0, wLevMa);
                    igl::writeOFF("../curr.off", X, seams->F);
                    std::cout << "steps: " << steps<< "\n";
                    

					iterations_elapsed++;
					if (progressiveDevelopabilityIncrease && (iterations_elapsed % iterations) == 0)
					{
						wDev *= 2;
						//wDev = std::pow(iterations_elapsed, 3);
					}

                    const double maxK = optimization->currentDevelopabilityObjective.cwiseAbs().maxCoeff();
                    if (maxK < max_defect)
                    {
                        std::cout << "                                --> reached threshold" << std::endl;
						export_mesh_data("_thres");
                        break;
                    }

					if (maxK < 5*max_defect && maxK < smallest_maxK)
					{
						export_mesh_data("_automin");
						smallest_maxK = maxK;
					}

                    /*
                    Eigen::VectorXd ad;
                    optimization->angleDefect(X, ad);
                    
                    // take the mean only over non zero values to exclude boundary vertices
                    double mean = .0;
                    int cnt2 = 0;
                    for(int i = 0; i < ad.size(); ++i)
                    {
                        if(ad(i))
                        {
                            mean += std::abs(ad(i));
                            ++cnt2;
                        }
                    }
                    
                    mean /= cnt2;
                    
                    std::cout << " gauss: " << ad.norm() << " mean: " << mean << " max: " <<  ad.cwiseAbs().maxCoeff() << " steps: " << steps<< "\n";
					if (ad.cwiseAbs().maxCoeff() <= max_defect)
						std::cout << "    --> reached threshold" << std::endl;
                    
                    std::ofstream file("../ad");
                    file << ad;
                    file.close();
                    */
                }
                

                
                std::cout << "done optimization" << std::endl;
                
                //---
                Eigen::VectorXd K;
                get_gaussian_curvature(X, seams->F, K, true);
                std::cout << "        igl::gaussian_curvature |K|max = " << K.cwiseAbs().maxCoeff() << ", Kmean = " << K.mean() << "\n" << std::endl;
                //---

                viewer.data(idProjection).set_mesh(X, seams->F);
               
                Eigen::MatrixXd N;
                igl::per_vertex_normals(X, seams->F, N);
                viewer.data(idProjection).set_normals(N);
                viewer.data(idProjection).set_colors(projColors);
            }
        }
        
        if (ImGui::CollapsingHeader("Output", ImGuiTreeNodeFlags_DefaultOpen))
        {
            if (ImGui::Button("stitch and export"))
            {
				export_mesh_data("_manual");
				/*
                // --- Stitch mesh to show curvature on seams ---
                Eigen::MatrixXd stitched_V;
                Eigen::MatrixXi stitched_F;
                seams->stitchMesh(X, stitched_V, stitched_F);

                igl::writeOBJ("../out_stitched.obj", stitched_V, stitched_F);
                igl::writeOBJ("../out.obj", X, seams->F);


                // --- Compute result metrics and write to file ---
                Eigen::VectorXd K;
                get_gaussian_curvature(X, seams->F, K, true);
                std::cout << "\n(info) EXPORTING mesh with normalized |K|max = " << K.cwiseAbs().maxCoeff() << ", Kmean = " << K.mean() << "\n" << std::endl;

                //hausdorff distance: original -> developable
                double d0;
                igl::hausdorff(V0, F0, X, seams->F, d0);
                double bb0 = get_BBDiagonal(V0);

                //hausdorff distance: developable -> original
                double d1;
                igl::hausdorff(X, seams->F, V0, F0, d1);
                double bb1 = get_BBDiagonal(X);

                std::ofstream file_metrics("../out_metrics.txt");
                	file_metrics << "|Kmax | = " << K.cwiseAbs().maxCoeff() << "\n";
                	file_metrics << " Kmean  = " << K.mean() << "\n";
                	file_metrics << " Kmin   = " << K.minCoeff() << "\n";
                	file_metrics << " Kmax   = " << K.maxCoeff() << "\n";

                	file_metrics << "\nhausdorff distance = " << std::max(d0/bb0, d1/bb1) << "\n";
                	file_metrics << "  original -> developable " << "\n";
                	file_metrics << "    absolute = " << d0 << "\n";
                	file_metrics << "    relative = " << d0/bb0 << "(wrt diagonal = " << bb0 << ")\n";
                	file_metrics << "  developable -> original " << "\n";
                	file_metrics << "    absolute = " << d1 << "\n";
                	file_metrics << "    relative = " << d1/bb1 << "(wrt diagonal = " << bb1 << ")\n";
                file_metrics.close();
				*/
            }
        }
    };
    
    
    viewer.core().background_color.setOnes();
    
 
    idOriginal = viewer.selected_data_index;
    
    idDogs = viewer.append_mesh();
    idProjection = viewer.append_mesh();
    idMinProjection = viewer.append_mesh();
    
    viewer.data(idOriginal).set_mesh(V, F);
    viewer.data(idOriginal).show_lines = 0;
   
    viewer.data(idDogs).set_mesh(dogset.Vall, dogset.Fall);
    viewer.data(idProjection).set_mesh(X, seams->F);
    viewer.data(idMinProjection).set_mesh(PV0, seams->F);
    
    
    
    viewer.data(idOriginal).set_colors(origColors);
    viewer.data(idDogs).set_colors(dogColors);
    viewer.data(idProjection).set_colors(projColors);
    viewer.data(idMinProjection).set_colors(projColors);
      
    
    viewer.data(idOriginal).set_visible(showOriginal);
    viewer.data(idDogs).set_visible(showDogs);
    viewer.data(idProjection).set_visible(showProjection);
    viewer.data(idMinProjection).set_visible(showMinProjection);
    
    viewer.launch();
#endif
	return 0;
}
