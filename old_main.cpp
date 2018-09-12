#include <igl/cotmatrix.h>
#include <igl/viewer/Viewer.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>

#include <vector>
#include <sstream>
#include <fstream>
#include "nanoflann.hpp"
#include "tinyply.h"

using namespace tinyply;
using namespace Eigen;
using namespace std;
using namespace nanoflann;

void read_ply_file(const std::string & filename, std::vector<float> &verts, std::vector<uint32_t> &faces,
                   uint32_t* vertexCount, uint32_t* faceCount);
void ply_to_matrix(std::vector<float> &verts, std::vector<uint32_t> faces, Eigen::MatrixXd &V, Eigen::MatrixXi &F,
                   const uint32_t vertexCount, const uint32_t faceCount);

Eigen::Vector3d calc_barycentre(Eigen::MatrixXd &VA, const uint32_t vertexCountA);

void ICP(Eigen::MatrixXd &VA, Eigen::MatrixXd &VB, const uint32_t vertexCountA);

int main()
{

    // Read mesh A
    // Declare variables to hold vertices and faces
    uint32_t vertexCountA;
    uint32_t faceCountA;
    std::vector<float> vertsA;
    std::vector<uint32_t> facesA;

    // Parse ply file
    read_ply_file("bun_zipper_res2.ply", vertsA, facesA, &vertexCountA, &faceCountA);

    // Declare matrices to store vertices and faces
    Eigen::MatrixXd VA;
    Eigen::MatrixXi FA;
    VA.resize(vertexCountA, 3);
    FA.resize(faceCountA, 3);

    // Store values in the matrices
    ply_to_matrix(vertsA, facesA, VA, FA, vertexCountA, faceCountA);

    //Define transformation and define mesh B, which will be mesh A transformed
    Eigen::Vector3d transl(0.03,0.03,0.03);
    Eigen::MatrixXd VB;
    Eigen::MatrixXi FB;
    VB.resize(vertexCountA, 3);
    FB.resize(faceCountA, 3);

    // Apply transformation to mesh A, obtaining mesh B
    size_t size = VA.rows();
    for(int row=0; row < size; row++){
        Eigen::Vector3d rowA(VA(row, 0), VA(row, 1), VA(row, 2));
        Eigen::Vector3d newRow;
        newRow(0) = rowA(0) + transl(0);
        newRow(1) = rowA(1) + transl(1);
        newRow(2) = rowA(2) + transl(2);
        VB(row, 0) = newRow(0);
        VB(row, 1) = newRow(1);
        VB(row, 2) = newRow(2);
    }
    FB = FA;
    uint32_t vertexCountB = vertexCountA;
    uint32_t faceCountB = faceCountA;

    /*
    //Define rotation
    Eigen::AngleAxis<float> rot(1.5, Vector3f(0,1,0));
    //Translation<float, 3> trans(1, 1, 1);
    Eigen::MatrixXd VB;
    Eigen::MatrixXi FB;


    VB = rot*VA;
    FB = rot*FA;

    uint32_t vertexCountB = vertexCountA;
    uint32_t faceCountB = faceCountA;
    */
    /*
    // Read mesh 2
    uint32_t vertexCountB;
    uint32_t faceCountB;
    std::vector<float> vertsB;
    std::vector<uint32_t> facesB;
    read_ply_file("bun_zipper_res2.ply", vertsB, facesB, &vertexCountB, &faceCountB);

    //Debugging
    std::cout << "\tRead " << vertsB.size() << " total vertices (" << vertexCountB << " properties)." << std::endl;
    std::cout << "\tRead " << facesB.size() << " total faces (triangles) (" << faceCountB << " properties)." << std::endl;

    Eigen::MatrixXd VB;
    Eigen::MatrixXi FB;
    VB.resize(vertexCountB, 3);
    FB.resize(faceCountB, 3);

    ply_to_matrix(vertsB, facesB, VB, FB, vertexCountB, faceCountB);
    */
    /*
     For debugging:

    for(int i=0; i< verts.size(); i++){
        std::cout << "VERTICES:" << verts[i] << std::endl;
    }
    for(int i=0; i< faces.size(); i++){
        std::cout << "FACES:" << faces[i] << std::endl;
    }
    */

    /*For debugging:
    std::string sep = "\n----------------------------------------\n";
    std::cout << "VERTICES AFTER" << std::endl;
    std::cout << V << sep;
    std::cout << "FACES AFTER " << std::endl;
    std::cout << F << sep;
    */

    // Plot the mesh

    /*
    // Visualize one mesh
    igl::viewer::Viewer viewer;
    //viewer.data.set_mesh(VA, FA);
    viewer.data.add_points(VA, Eigen::RowVector3d(0, 0, 0)); //Plots only the points, for debugging
    viewer.launch();
    */

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
    viewer.data.set_mesh(V,F);
    viewer.data.set_colors(C);
    viewer.data.set_face_based(true);
    viewer.launch();

    //Call ICP

    ICP(VA, VB, vertexCountA);

    // Visualize the two meshes after alignment
    // Concatenate (VA,FA) and (VB,FB) into (V,F)
    V = Eigen::MatrixXd::Zero(VA.rows() + VB.rows(), VA.cols());
    //V.resize(VA.rows()+VB.rows(),VA.cols());
    V<<VA,VB;
    //Eigen::MatrixXi F;
    F = Eigen::MatrixXi::Zero(FA.rows()+FB.rows(), FA.cols());
    //F.resize(FA.rows()+FB.rows(),FA.cols());
    F<<FA,(FB.array()+VA.rows());


    // blue color for faces of first mesh, orange for second
    /*Eigen::MatrixXd C;
    //C.resize(F.rows(),3);
    C<< Eigen::RowVector3d(0.2,0.3,0.8).replicate(FA.rows(),1),
            Eigen::RowVector3d(1.0,0.7,0.2).replicate(FB.rows(),1);
    */
    //VIEW two meshes
    //igl::viewer::Viewer viewer;
    viewer.data.set_mesh(V,F);
    viewer.data.set_colors(C);
    viewer.data.set_face_based(true);
    viewer.launch();

}


void read_ply_file(const std::string & filename, std::vector<float> &verts, std::vector<uint32_t> &faces,
                   uint32_t* vertexCount, uint32_t* faceCount)
{
    // Tinyply can and will throw exceptions at you!
    try
    {

        // Read the file and create a std::istringstream suitable
        // for the lib -- tinyply does not perform any file i/o.
        std::ifstream ss(filename, std::ios::binary);

        // Parse the ASCII header fields
        PlyFile file(ss);

        for (auto e : file.get_elements())
        {
            std::cout << "element - " << e.name << " (" << e.size << ")" << std::endl;
            for (auto p : e.properties)
            {
                std::cout << "\tproperty - " << p.name << " (" << PropertyTable[p.propertyType].str << ")" << std::endl;
            }
        }
        std::cout << std::endl;

        for (auto c : file.comments)
        {
            std::cout << "Comment: " << c << std::endl;
        }

        // Define containers to hold the extracted data. The type must match
        // the property type given in the header. Tinyply will interally allocate the
        // the appropriate amount of memory.
        //std::vector<float> verts;
        std::vector<float> norms;
        std::vector<uint8_t> colors;

        //std::vector<uint32_t> faces;
        std::vector<float> uvCoords;

        uint32_t  normalCount, colorCount, faceTexcoordCount, faceColorCount;
        normalCount = colorCount = faceTexcoordCount = faceColorCount = 0;

        // The count returns the number of instances of the property group. The vectors
        // above will be resized into a multiple of the property group size as
        // they are "flattened"... i.e. verts = {x, y, z, x, y, z, ...}
        *vertexCount = file.request_properties_from_element("vertex", { "x", "y", "z" }, verts);
        normalCount = file.request_properties_from_element("vertex", { "nx", "ny", "nz" }, norms);
        colorCount = file.request_properties_from_element("vertex", { "red", "green", "blue", "alpha" }, colors);

        // For properties that are list types, it is possibly to specify the expected count (ideal if a
        // consumer of this library knows the layout of their format a-priori). Otherwise, tinyply
        // defers allocation of memory until the first instance of the property has been found
        // as implemented in file.read(ss)
        *faceCount = file.request_properties_from_element("face", { "vertex_indices" }, faces, 3);
        faceTexcoordCount = file.request_properties_from_element("face", { "texcoord" }, uvCoords, 6);


        file.read(ss);

        std::cout << "Face1:" << std::endl;
        // Good place to put a breakpoint!

        //std::cout << "\tRead " << verts.size() << " total vertices (" << vertexCount << " properties)." << std::endl;
        std::cout << "\tRead " << norms.size() << " total normals (" << normalCount << " properties)." << std::endl;
        std::cout << "\tRead " << colors.size() << " total vertex colors (" << colorCount << " properties)." << std::endl;
        std::cout << "\tRead " << uvCoords.size() << " total texcoords (" << faceTexcoordCount << " properties)." << std::endl;
        std::cout << "\tRead " << verts.size() << " total vertices (" << vertexCount << " properties)." << std::endl;
        std::cout << "\tRead " << faces.size() << " total faces (triangles) (" << faceCount << " properties)." << std::endl;
        /*
        // Fixme - tinyply isn't particularly sensitive to mismatched properties and prefers to crash instead of throw. Use
        // actual data from parsed headers instead of manually defining properties added to a new file like below:

        std::filebuf fb;
        fb.open("converted.ply", std::ios::out | std::ios::binary);
        std::ostream outputStream(&fb);

        PlyFile myFile;

        myFile.add_properties_to_element("vertex", { "x", "y", "z" }, verts);
        myFile.add_properties_to_element("vertex", { "red", "green", "blue" }, colors);
        myFile.add_properties_to_element("face", { "vertex_indices" }, faces, 3, PlyProperty::Type::UINT8);

        myFile.comments.push_back("generated by tinyply");
        myFile.write(outputStream, true);

        fb.close();
        */
    }

    catch (const std::exception & e)
    {
        std::cerr << "Caught exception: " << e.what() << std::endl;
    }
}

void ply_to_matrix(std::vector<float> &verts, std::vector<uint32_t> faces, Eigen::MatrixXd &V, Eigen::MatrixXi &F,
                   const uint32_t vertexCount, const uint32_t faceCount){
    std::vector<double> verts_double(verts.begin(), verts.end());
    std::vector<int> faces_int(faces.begin(), faces.end());



    int index = 0;
    for(int i = 0 ; i < vertexCount; i++) {
        for(int c = 0 ; c < 3 ; c++) {
            //std::cout << "\tindex:" << index<< std::endl; //For debugging
            V(i, c) = verts_double[index];
            index = index + 1;
        }
    }
    index = 0;
    for(int i = 0 ; i < faceCount ; i++) {
        for(int c = 0 ; c < 3 ; c++) {

            //std::cout << "\tindex:" << index << std::endl;// For debugging
            F(i, c) = faces_int[index];
            index = index +1;

        }
    }

}


Eigen::Vector3d calc_barycentre(Eigen::MatrixXd &V, const uint32_t vertexCountA){
    //Calculate barycentre for points A and B
    Eigen::Vector3d sums = Eigen::Vector3d::Zero();
    for(int row = 0; row < V.rows(); row++){
        for(int col = 0; col < 3; col++){
            sums(col) = sums(col) + V(row, col);
        }
    }
    Eigen::Vector3d bary = sums/vertexCountA;
    // For debugging purposes
    // std::cout<< "Sums" << sums(0) << " " << sums(1) << " " << sums(2) << std::endl;
    // std::cout<< "Barycentre" <<  bary << std::endl;
    return bary;
}

void ICP(Eigen::MatrixXd &VA, Eigen::MatrixXd &VB, const uint32_t vertexCountA){
    //
    // START of the ICP algorithm. We want to align the vertices of mesh A (source) to mesh B (model).
    //
    // Create data structure to hold vertices of mesh B, K-D tree.
    size_t nSamples = VB.rows();
    Eigen::Matrix<double, Dynamic, Dynamic>  mat(nSamples, 3);
    for (size_t i = 0; i < VB.rows(); i++)
        for (size_t d = 0; d < 3; d++)
            mat(i,d) = VB(i,d);

    typedef KDTreeEigenMatrixAdaptor <Eigen::Matrix<double, Dynamic, Dynamic>> my_kd_tree_t;

    my_kd_tree_t  mat_index(mat, 10);
    mat_index.index->buildIndex();

    // Matrix that holds the points in mesh B which are closest to mesh A
    Eigen::MatrixXd closest_VB;
    closest_VB.resize(vertexCountA, 3);

    // Error variables
    float oldError = 1;
    float error = 1;
    // Iterate
    for(int iter = 0; iter < 10; iter++) {
        // Store the error of the previous iteration
        oldError = error;
        // MATCHING
        // Iterate over all the points in mesh A
        for (int row = 0; row < VA.rows(); row++) {
            // Define the current query point
            std::vector<double> query_pt(3);
            for (size_t d = 0; d < 3; d++)
                query_pt[d] = VA(row, d);

            // Returned index of the closest point
            size_t ret_index;
            // Distance squared between the two points
            double out_dist_sqr;

            nanoflann::KNNResultSet<double> resultSet(1);

            resultSet.init(&ret_index, &out_dist_sqr);

            mat_index.index->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));

            //For debugging purposes:
            //std::cout << "Index: " << ret_index << std::endl;
            //std::cout << "Out Dist Sqr" << out_dist_sqr << std::endl;
            //std::cout << "Vertex: " << VB(ret_index, 0) << VB(ret_index, 1) << VB(ret_index, 2) << std::endl;
            closest_VB(row, 0) = VB(ret_index, 0);
            closest_VB(row, 1) = VB(ret_index, 1);
            closest_VB(row, 2) = VB(ret_index, 2);

        }

        // COMPUTING THE REGISTRATION

        Eigen::Vector3d baryA = calc_barycentre(VA, vertexCountA);
        Eigen::Vector3d baryB = calc_barycentre(closest_VB, vertexCountA);

        Eigen::MatrixXd VA_hat;
        Eigen::MatrixXd VB_hat;
        VA_hat.resize(vertexCountA, 3);
        VB_hat.resize(vertexCountA, 3);

        for (int row = 0; row < VA.rows(); row++) {
            for (int col = 0; col < 3; col++) {
                VA_hat(row, col) = VA(row, col) - baryA(col);
                VB_hat(row, col) = closest_VB(row, col) - baryB(col);
            }
        }

        Eigen::Matrix3d A = Eigen::Matrix3d::Zero();

        for (int row = 0; row < VA_hat.rows(); row++) {
            Eigen::Vector3d rowA_hat;
            Eigen::Vector3d rowB_hat;

            rowA_hat(0) = VA_hat(row, 0);
            rowA_hat(1) = VA_hat(row, 1);
            rowA_hat(2) = VA_hat(row, 2);

            rowB_hat(0) = VB_hat(row, 0);
            rowB_hat(1) = VB_hat(row, 1);
            rowB_hat(2) = VB_hat(row, 2);

            A = A + rowA_hat * (rowB_hat.transpose());
        }

        std::cout << "matrix: " << A << std::endl;

        // Compute singular value decomposition
        Eigen::JacobiSVD<MatrixXd> svd(A, ComputeThinU | ComputeThinV);
        /* // For debugging purposes
        std::cout << "Its singular values are:" << std::endl << svd.singularValues() << std::endl;
        std::cout << "Its left singular vectors are the columns of the thin U matrix:" << std::endl << svd.matrixU() << std::endl;
        std::cout << "Its right singular vectors are the columns of the thin V matrix:" << std::endl << svd.matrixV() << std::endl;
        */

        // Compute rotation matrix
        Eigen::Matrix3d U_matrix = svd.matrixU();
        Eigen::Matrix3d V_matrix = svd.matrixV();
        Eigen::Matrix3d R = V_matrix * (U_matrix.transpose());

        //std::cout << "Rotation matrix: " << std::endl << R << std::endl;

        // Compute translation vector:
        Eigen::Vector3d t = baryB - R * baryA;

        //std::cout << "Translation vector: " << std::endl << t << std::endl;

        // UPDATING SCAN ALIGNMENT
        float sum = 0;
        for (int row = 0; row < VA.rows(); row++) {
            Eigen::Vector3d rowVA;
            Eigen::Vector3d row_closest_VB;
            rowVA(0) = VA(row, 0);
            rowVA(1) = VA(row, 1);
            rowVA(2) = VA(row, 2);

            row_closest_VB(0) = closest_VB(row, 0);
            row_closest_VB(1) = closest_VB(row, 1);
            row_closest_VB(2) = closest_VB(row, 2);

            rowVA = R * rowVA + t;

            VA(row, 0) = rowVA(0);
            VA(row, 1) = rowVA(1);
            VA(row, 2) = rowVA(2);

            Eigen::Vector3d diff;
            diff = row_closest_VB - rowVA;
            sum = sum + (diff.transpose()) * diff;

        }
        error = sum/VA.rows();
        std::cout << "Iteration: " << iter << std::endl;
        std::cout << "The previous error was: " << oldError << std::endl;
        std::cout << "The current error is: " << error << std::endl;
        std::cout << "The difference (previous - current error) is: " << oldError - error << std::endl;

    }

}