#include <igl/cotmatrix.h>
#include <igl/viewer/Viewer.h>
#include <igl/per_vertex_normals.h>_
#include <igl/readOBJ.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>
#include "nanoflann.hpp"
#include "tinyply.h"
#include "calc_barycentre.h"
#include "icp.h"
#include "ICP_p2plane.h"
#include "point_to_plane.h"

using namespace tinyply;
using namespace Eigen;
using namespace std;
using namespace nanoflann;

Eigen::MatrixXd VA;
Eigen::MatrixXi FA;
Eigen::MatrixXd VB;
Eigen::MatrixXi FB;
Eigen::MatrixXd NB;


bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier);

void add_noise(Eigen::MatrixXd &V, float sigma);



int main()
{
    // Read mesh A

    // Declare matrices to store vertices and faces

    igl::readOBJ("mesh1.obj", VA, FA);

    uint32_t vertexCountA = VA.rows();
    uint32_t faceCountA = FA.rows();
    NB.resize(NB.rows(), 3);

    // Read mesh B

    igl::readOBJ("mesh2.obj", VB, FB);


    // Visualize the two meshes together
    // Concatenate (VA,FA) and (VB,FB) into (V,F)
    Eigen::MatrixXd V;
    V.resize(VA.rows()+VB.rows(),VA.cols());
    V<<VA,VB;
    Eigen::MatrixXi F;
    F.resize(FA.rows()+FB.rows(),FA.cols());
    F<<FA,(FB.array()+VA.rows());

    // blue color for faces of first mesh, orange for second
    Eigen::MatrixXd C;
    C.resize(F.rows(),3);
    C<< Eigen::RowVector3d(0.2,0.3,0.8).replicate(FA.rows(),1),
            Eigen::RowVector3d(1.0,0.7,0.2).replicate(FB.rows(),1);

    //VIEW two meshes
    igl::viewer::Viewer viewer;
    viewer.callback_key_down = &key_down;
    viewer.data.set_mesh(V,F);
    viewer.data.set_colors(C);
    viewer.data.set_face_based(true);
    viewer.launch();

}

bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier) {
    std::cout<<"Key: "<<key<<" "<<(unsigned int)key<<std::endl;
    if (key == '1') {
        // Call ICP
        ICP(VA, VB);


    }
    else if(key == '2'){
        viewer.data.clear();
        //Define transformation and define mesh B, which will be mesh A transformed
        //Define rotation
        Eigen::Matrix3d rot;
        rot = AngleAxisd(0.5, Vector3d(1, 0, 0));

        Eigen::MatrixXd VB_transposed;
        VB.resize(VA.rows(), 3);
        VB_transposed.resize(3, VA.rows());
        FB.resize(FA.rows(), 3);
        Eigen::MatrixXd VA_transposed;
        VA_transposed.resize(3, VA.rows());
        VA_transposed = VA.transpose();
        VB_transposed= rot*VA_transposed;
        VB = VB_transposed.transpose();
        FB = FA;
        Eigen::MatrixXd V;
        V.resize(VA.rows()+VB.rows(),VA.cols());
        V<<VA,VB;
        Eigen::MatrixXi F;
        F.resize(FA.rows()+FB.rows(),FA.cols());
        F<<FA,(FB.array()+VA.rows());

        // blue color for faces of first mesh, orange for second
        Eigen::MatrixXd C;
        C.resize(F.rows(),3);
        C<< Eigen::RowVector3d(0.2,0.3,0.8).replicate(FA.rows(),1),
                Eigen::RowVector3d(1.0,0.7,0.2).replicate(FB.rows(),1);

        //VIEW two meshes
        viewer.data.set_mesh(V,F);
        viewer.data.set_colors(C);
        viewer.data.set_face_based(true);

        // Call ICP
        ICP(VA, VB);
        /*
        V.resize(VA.rows()+VB.rows(),VA.cols());
        V<<VA,VB;
        F.resize(FA.rows()+FB.rows(),FA.cols());
        F<<FA,(FB.array()+VA.rows());


        //VIEW two meshes
        viewer.data.set_mesh(V,F);
        viewer.data.set_colors(C);
        viewer.data.set_face_based(true);
        */

    }
    else if(key == '3'){
        //Add noise
        float sigma = 0.001;
        add_noise(VA, sigma);
        Eigen::MatrixXd V;
        V.resize(VA.rows()+VB.rows(),VA.cols());
        V<<VA,VB;
        Eigen::MatrixXi F;
        F.resize(FA.rows()+FB.rows(),FA.cols());
        F<<FA,(FB.array()+VA.rows());

        // blue color for faces of first mesh, orange for second
        Eigen::MatrixXd C;
        C.resize(F.rows(),3);
        C<< Eigen::RowVector3d(0.2,0.3,0.8).replicate(FA.rows(),1),
                Eigen::RowVector3d(1.0,0.7,0.2).replicate(FB.rows(),1);

        //VIEW two meshes
        viewer.data.set_mesh(V,F);
        viewer.data.set_colors(C);
        viewer.data.set_face_based(true);

        //Call ICP
        //ICP(VA, VB);
        /*
        V.resize(VA.rows()+VB.rows(),VA.cols());
        V<<VA,VB;
        F.resize(FA.rows()+FB.rows(),FA.cols());
        F<<FA,(FB.array()+VA.rows());


        //VIEW two meshes
        viewer.data.set_mesh(V,F);
        viewer.data.set_colors(C);
        viewer.data.set_face_based(true);
         */
    }
    else if (key == '6'){
        //Calculate the normals in the destination mesh
        igl::per_vertex_normals(VB, FB, NB);

        ICP_p2plane(VA, VB, NB);
    }
    viewer.data.clear();
    Eigen::MatrixXd V;
    V.resize(VA.rows()+VB.rows(),VA.cols());
    V<<VA,VB;
    Eigen::MatrixXi F;
    F.resize(FA.rows()+FB.rows(),FA.cols());
    F<<FA,(FB.array()+VA.rows());

    // blue color for faces of first mesh, orange for second
    Eigen::MatrixXd C;
    C.resize(F.rows(),3);
    C<< Eigen::RowVector3d(0.2,0.3,0.8).replicate(FA.rows(),1),
            Eigen::RowVector3d(1.0,0.7,0.2).replicate(FB.rows(),1);

    //VIEW two meshes
    viewer.data.set_mesh(V,F);
    viewer.data.set_colors(C);
    viewer.data.set_face_based(true);
    return false;
}



void add_noise(Eigen::MatrixXd &V, float sigma){
    V = V + sigma*(Eigen::MatrixXd::Random(V.rows(),V.cols()));
}


