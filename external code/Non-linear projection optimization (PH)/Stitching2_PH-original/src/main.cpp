#include <igl/readOBJ.h>
#include <igl/writeOBJ.h>
#include <igl/writeOFF.h>
#include <igl/facet_components.h>
#include <igl/point_mesh_squared_distance.h>
#include <igl/remove_unreferenced.h>
#include <igl/edges.h>
#include <igl/upsample.h>


#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>

#include "DogSet.hpp"
#include "ProcessSeams.hpp"
#include "objectives.hpp"
#include <fstream>

// Main parameters to set:
//-------------------------


bool upsample = false; // upsample?
double gcSmooth = 0.2; // smoothness for graph cut

// weights for optimization:
double wSeam = 10;          // seam smoothness
double wSmooth = 5;         // overall smoothness
double wPair = 5;           // patch boundaries should match
double wProx = 5;           // keep close to the dogs
double wDev = 1000 * 8;     // minimize angle defect (progessively increased 1000 * iter^3
double wLevMa = 5.0;        // Levenberg-Marquardt regularization (sometimes prevents getting stuck)
double wSeamProx = 0;       // the mesh boundary should stay close to its initial projection on the dogs



int main(const int argc, const char* argv[])
{
    std::string dogFile = "../../data/new/side11_40/dog.obj";
    std::string inputFile = "../../data/new/side11_40/input.obj";
    

    if(argc == 3)
    {
        dogFile = argv[1];
        inputFile = argv[2];
    }
    
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    Eigen::MatrixXd V0;
    Eigen::MatrixXi F0;
     
    igl::readOBJ(inputFile, V0, F0);
  
    
    if(upsample)
    {
        igl::upsample(V0, F0, V, F, 1);
    } else
    {
        V = V0;
        F = F0;
    }
        
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
    
#if 0 // disable GUI

    Eigen::MatrixXd dogColors;
    
    Eigen::VectorXi C;
    igl::facet_components(dogset.Fall, C);
    igl::colormap(igl::COLOR_MAP_TYPE_JET, C, 0.0, C.maxCoeff(), dogColors);
    

    Eigen::MatrixXd origColors;
    igl::colormap(igl::COLOR_MAP_TYPE_JET, seams.L, 0.0, seams.L.maxCoeff(), origColors);
        
    Eigen::MatrixXd PV;
    dogset.projectPoints(seams.V, seams.L, PV);
    
    Eigen::MatrixXd projColors;
    igl::colormap(igl::COLOR_MAP_TYPE_JET, seams.L, 0.0, seams.L.maxCoeff(), projColors);


    igl::opengl::glfw::Viewer viewer;

    // Attach a menu plugin
    igl::opengl::glfw::imgui::ImGuiMenu menu;
    viewer.plugins.push_back(&menu);

    
    bool showOriginal = true;
    bool showDogs = false;
    bool showProjection = false;
    bool showMinProjection = false;
    
    int idOriginal, idDogs, idProjection, idMinProjection;

 
    double t0 = 1.;
    int iterations = 1;
    
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
            
            if(ImGui::Button("smooth labels"))
            {
                X = PV;
                viewer.data(idProjection).set_mesh(X, seams.F);
            }
        }
        
        if (ImGui::CollapsingHeader("Optimization Weights", ImGuiTreeNodeFlags_DefaultOpen))
        {
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
                ImGui::InputInt("Iterations", &iterations);
            }
            
            if(ImGui::Button("reset"))
            {
                X = seams.V;
                viewer.data(idProjection).set_mesh(X, seams.F);
                
                Eigen::MatrixXd N;
                igl::per_vertex_normals(X, seams.F, N);
                viewer.data(idProjection).set_normals(N);
            
                viewer.data(idProjection).set_colors(projColors);
            }
            
            if(ImGui::Button("recompute"))
            {
              //  OptimizationObjectives optimization(seams.V, seams.L, seams.F, dogset, pairs);
                optimization.wSeamSmooth = wSeam;
                optimization.wSmooth = wSmooth;
                optimization.wPair = wPair;
                optimization.wProx = wProx;
                optimization.wDev = wDev;
                optimization.wSeamProx = wSeamProx;
                
                
                std::cout << "start optimization" << std::endl;
                std::cout << "initial energy: " << optimization.objective(X) << "\n";
                
                for(int k = 0; k < iterations; ++k)
                {
                    const int steps = optimization.optimize(X, t0, wLevMa);
                    
                    
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
                
                    igl::writeOFF("../curr.off", X, seams.F);
                    
                    std::ofstream file("../ad");
                    file << ad;
                    file.close();
                }
                

                
                std::cout << "done optimization" << std::endl;
          
                viewer.data(idProjection).set_mesh(X, seams.F);
               
                Eigen::MatrixXd N;
                igl::per_vertex_normals(X, seams.F, N);
                viewer.data(idProjection).set_normals(N);
                viewer.data(idProjection).set_colors(projColors);
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
    viewer.data(idProjection).set_mesh(X, seams.F);
    viewer.data(idMinProjection).set_mesh(PV0, seams.F);
    
    
    
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
