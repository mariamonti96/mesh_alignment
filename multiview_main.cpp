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

using namespace tinyply;
using namespace Eigen;
using namespace std;
using namespace nanoflann;

Eigen::MatrixXd VA;
Eigen::MatrixXi FA;
Eigen::MatrixXd VB;
Eigen::MatrixXi FB;
Eigen::MatrixXd NB;
Eigen::MatrixXd VC;
Eigen::MatrixXi FC;
Eigen::MatrixXd VD;
Eigen::MatrixXi FD;
Eigen::MatrixXd VE;
Eigen::MatrixXi FE;

bool key_down(igl::viewer::Viewer& viewer, unsigned char key, int modifier);

void add_noise(Eigen::MatrixXd &V, float sigma);
void pseudo_inverse(Eigen::MatrixXd &A, Eigen::MatrixXd &A_pinv);

void point_to_plane(Eigen::MatrixXd &VA, Eigen::MatrixXd &VB,Eigen::MatrixXd &NB,
                    Eigen::Matrix3d &R, Eigen::Vector3d &t);

Eigen::Vector3d calc_barycentre(Eigen::MatrixXd &V);

void ICP(Eigen::MatrixXd &VA, Eigen::MatrixXd &VB);
void ICP_p2plane(Eigen::MatrixXd &VA, Eigen::MatrixXd &VB, Eigen::MatrixXd &NB);
void ICP_multiview(Eigen::MatrixXd &VA, Eigen::MatrixXd &VB, Eigen::MatrixXd &VC);
int main()
{

    // Read mesh A

    // Declare matrices to store vertices and faces
    //write_ply_example("bla.ply");

    igl::readOBJ("mesh1.obj", VA, FA);

    NB.resize(NB.rows(), 3);

    // Read mesh 2

    igl::readOBJ("mesh2.obj", VB, FB);

    igl::readOBJ("mesh3.obj", VC, FC);

    igl::readOBJ("mesh4.obj", VD, FD);

    igl::readOBJ("mesh5.obj", VE, FE);



    // Visualize the five meshes together
    // Concatenate (VA,FA) and (VB,FB) into (V,F)
    Eigen::MatrixXd V;
    V.resize(VA.rows()+VB.rows()+VC.rows(),VA.cols());
    V<<VA,VB,VC;
    Eigen::MatrixXi F;
    F.resize(FA.rows()+FB.rows()+FC.rows(),FA.cols());
    F<<FA,(FB.array()+VA.rows()), (FC.array() + VA.rows() + VB.rows());

    // blue color for faces of first mesh, orange for second
    Eigen::MatrixXd C;
    C.resize(F.rows(),3);
    C<< Eigen::RowVector3d(0.2,0.3,0.8).replicate(FA.rows(),1),
            Eigen::RowVector3d(1.0,0.7,0.2).replicate(FB.rows(),1),
            Eigen::RowVector3d(1.0, 0.0, 0.0).replicate(FC.rows(),1);

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

    }
    else if (key == '6'){
        //Calculate the normals in the destination mesh
        igl::per_vertex_normals(VB, FB, NB);

        ICP_p2plane(VA, VB, NB);
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
    }
    else if(key == '5'){
        ICP_multiview(VA, VB, VC);

        // Visualize the five meshes together
        // Concatenate (VA,FA) and (VB,FB) into (V,F)
        viewer.data.clear();
        Eigen::MatrixXd V;
        V.resize(VA.rows()+VB.rows()+VC.rows(),VA.cols());
        V<<VA,VB,VC;
        Eigen::MatrixXi F;
        F.resize(FA.rows()+FB.rows()+FC.rows(),FA.cols());
        F<<FA,(FB.array()+VA.rows()), (FC.array() + VA.rows() + VB.rows());

        // blue color for faces of first mesh, orange for second
        Eigen::MatrixXd C;
        C.resize(F.rows(),3);
        C<< Eigen::RowVector3d(0.2,0.3,0.8).replicate(FA.rows(),1),
                Eigen::RowVector3d(1.0,0.7,0.2).replicate(FB.rows(),1),
                Eigen::RowVector3d(1.0, 0.0, 0.0).replicate(FC.rows(),1);

    }

    return false;
}

Eigen::Vector3d calc_barycentre(Eigen::MatrixXd &V){
    //Calculate barycentre for points A and B
    Eigen::Vector3d sums = Eigen::Vector3d::Zero();
    for(int row = 0; row < V.rows(); row++){
        for(int col = 0; col < 3; col++){
            sums(col) = sums(col) + V(row, col);
        }
    }
    Eigen::Vector3d bary = sums/V.rows();
    // For debugging purposes
    // std::cout<< "Sums" << sums(0) << " " << sums(1) << " " << sums(2) << std::endl;
    // std::cout<< "Barycentre" <<  bary << std::endl;
    return bary;
}

void ICP(Eigen::MatrixXd &VA, Eigen::MatrixXd &VB){
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
    closest_VB.resize(VA.rows(), 3);

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

        Eigen::Vector3d baryA = calc_barycentre(VA);
        Eigen::Vector3d baryB = calc_barycentre(closest_VB);

        Eigen::MatrixXd VA_hat;
        Eigen::MatrixXd VB_hat;
        VA_hat.resize(VA.rows(), 3);
        VB_hat.resize(VA.rows(), 3);

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

void pseudo_inverse(Eigen::MatrixXd &A, Eigen::MatrixXd &A_pinv){
    A_pinv = ((A.transpose()*A).inverse()).transpose();
}

void point_to_plane(Eigen::MatrixXd &VA, Eigen::MatrixXd &VB, Eigen::MatrixXd &NB,
                    Eigen::Matrix3d &R, Eigen::Vector3d &t){

    //std::cout << "Normals: " << N << std::endl; //For debugging purposes
    Eigen::Vector3d rowNB;
    Eigen::Vector3d rowVA;
    Eigen::Vector3d rowVB;
    Eigen::MatrixXd A;
    A.resize(VA.rows(), 6);
    Eigen::VectorXd b;
    b.resize(VB.rows(), 1);


    for(int row= 0; row < VA.rows(); row++){
        for(int col = 0; col < 3; col++){

            rowNB(col) = NB(row, col);
            rowVA(col) = VA(row, col);
            rowVB(col) = VB(row, col);
        }
        Eigen::VectorXd rowA;
        rowA.resize(6, 1);

        rowA << (rowVA.cross(rowNB)), rowNB;

        A.row(row) = rowA.transpose();

        b(row) = rowNB.dot(rowVB) - rowNB.dot(rowVA);
    }

    //std::cout << "A: "<< A << std::endl;
    //std::cout << "b: "<< b << std::endl;
    VectorXd x = A.jacobiSvd(ComputeThinU | ComputeThinV).solve(b);

    //std::cout << "x: " << x << std::endl;
    R(0, 0) = 1;
    R(0, 1) = x(0)*x(1) - x(2);
    R(0, 2) = x(0)*x(2) + x(1);
    R(1, 0) = x(2);
    R(1, 1) = x(0)*x(1)*x(2) + 1;
    R(1, 2) = 1;
    R(2, 0) = -x(1);
    R(2, 1) = x(0);
    R(2, 2) = 1;

    t(0) = x(3);
    t(1) = x(4);
    t(2) = x(5);
}


void ICP_p2plane(Eigen::MatrixXd &VA, Eigen::MatrixXd &VB, Eigen::MatrixXd &NB){
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
    closest_VB.resize(VA.rows(), 3);
    Eigen::MatrixXd closest_NB;
    closest_NB.resize(VA.rows(), 3);

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

            closest_NB(row, 0) = NB(ret_index, 0);
            closest_NB(row, 1) = NB(ret_index, 1);
            closest_NB(row, 2) = NB(ret_index, 2);

        }

        // COMPUTING THE REGISTRATION

        Eigen::Matrix3d R;
        Eigen::Vector3d t;

        point_to_plane(VA, closest_VB, closest_NB, R, t);

        // UPDATING SCAN ALIGNMENT
        float sum = 0;
        for (int row = 0; row < VA.rows(); row++) {
            //std::cout << "iter : " << row << std::endl;
            //std::cout << "R: " << R << std::endl;
            //std::cout << "t: " << t << std::endl;
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

void add_noise(Eigen::MatrixXd &V, float sigma){
    V = V + sigma*(Eigen::MatrixXd::Random(V.rows(),V.cols()));
}


void ICP_multiview(Eigen::MatrixXd &VA, Eigen::MatrixXd &VB, Eigen::MatrixXd &VC){
    //
    // START of the ICP algorithm. We want to align the vertices of mesh A (source) to mesh B (model).
    //
    // Create data structure to hold vertices of mesh B, K-D tree.



    // Error variables
    float oldError = 1;
    float error = 1;
    // Iterate

    for(int iter = 0; iter < 10; iter++) {
        // Store the error of the previous iteration
        std::cout << "Iteration: " << iter << std::endl;
        oldError = error;

        size_t nSamplesA = VA.rows();
        Eigen::Matrix<double, Dynamic, Dynamic>  matA(nSamplesA, 3);
        for (size_t i = 0; i < VA.rows(); i++)
            for (size_t d = 0; d < 3; d++)
                matA(i,d) = VA(i,d);

        size_t nSamplesB = VB.rows();
        Eigen::Matrix<double, Dynamic, Dynamic>  matB(nSamplesB, 3);
        for (size_t i = 0; i < VB.rows(); i++)
            for (size_t d = 0; d < 3; d++)
                matB(i,d) = VB(i,d);

        size_t nSamplesC = VC.rows();
        Eigen::Matrix<double, Dynamic, Dynamic>  matC(nSamplesC, 3);
        for (size_t i = 0; i < VC.rows(); i++)
            for (size_t d = 0; d < 3; d++)
                matC(i,d) = VC(i,d);

        typedef KDTreeEigenMatrixAdaptor <Eigen::Matrix<double, Dynamic, Dynamic>> my_kd_tree_t;

        my_kd_tree_t  matA_index(matA, 10);
        matA_index.index->buildIndex();

        my_kd_tree_t  matB_index(matB, 10);
        matB_index.index->buildIndex();

        my_kd_tree_t  matC_index(matC, 10);
        matC_index.index->buildIndex();

        // Matrix that holds the points which are closest to mesh A
        Eigen::MatrixXd closest_to_A;
        closest_to_A.resize(VA.rows(), 3);

        // Matrix that holds the points which are closest to mesh B
        Eigen::MatrixXd closest_to_B;
        closest_to_B.resize(VB.rows(), 3);

        // Matrix that holds the points which are closest to mesh B
        Eigen::MatrixXd closest_to_C;
        closest_to_C.resize(VC.rows(), 3);

        ///
        ///SCAN A
        ///
        // MATCHING
        std::cout << "Computing scan A" << std::endl;
        for (int row = 0; row < VA.rows(); row++){
            std::vector<double> query_pt(3);
            for (size_t d = 0; d < 3; d++)
                query_pt[d] = VA(row, d);

            //Closest pts in the 2 other meshes
            Eigen::Matrix<double, Dynamic, Dynamic>  closest_to_query_pt(2, 3);
            //Returned indices of the closest pt in the 2 other meshes
            vector<size_t> ret_indices(2);
            //Returned index of the closest pt in the current mesh
            size_t current_ret_index;
            // Distance squared between the two points
            double out_dist_sqr;
            // Find index of closest point in mesh B
            nanoflann::KNNResultSet<double> resultSet(1);
            resultSet.init(&current_ret_index, &out_dist_sqr);
            matB_index.index->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));
            ret_indices[0] = current_ret_index;

            // Find values of closest point in mesh B
            closest_to_query_pt(0, 0) = VB(ret_indices[0], 0);
            closest_to_query_pt(0, 1) = VB(ret_indices[0], 1);
            closest_to_query_pt(0, 2) = VB(ret_indices[0], 2);


            // Find index of closest point in mesh C
            resultSet.init(&current_ret_index, &out_dist_sqr);
            matC_index.index->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));
            ret_indices[1] = current_ret_index;

            // Find values of closest point in mesh C
            closest_to_query_pt(1, 0) = VC(ret_indices[1], 0);
            closest_to_query_pt(1, 1) = VC(ret_indices[1], 1);
            closest_to_query_pt(1, 2) = VC(ret_indices[1], 2);

            //Build K-D tree to hold the closest pts in the other 2 meshes
            my_kd_tree_t  closest_to_query_pt_index(closest_to_query_pt, 10);
            closest_to_query_pt_index.index->buildIndex();
            // Find index of closest point between the 2 other closest points
            resultSet.init(&current_ret_index, &out_dist_sqr);
            closest_to_query_pt_index.index->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));

            // Find values of closest point in mesh C
            closest_to_A(row, 0)= closest_to_query_pt(current_ret_index, 0);
            closest_to_A(row, 1) = closest_to_query_pt(current_ret_index, 1);
            closest_to_A(row, 2) = closest_to_query_pt(current_ret_index, 2);

        }

        std::cout<< "First pt: " << VA(0, 0) << " " << VA(0, 1) << " " << VA(0, 2) << std::endl;
        std::cout<< "First match: " << closest_to_A(0, 0) << " " << closest_to_A(0, 1) << " " << closest_to_A(0, 2) << std::endl;
        // COMPUTING THE REGISTRATION

        Eigen::Vector3d baryA = calc_barycentre(VA);
        Eigen::Vector3d baryClosest_to_A = calc_barycentre(closest_to_A);

        Eigen::MatrixXd VA_hat;
        Eigen::MatrixXd closest_to_A_hat;
        VA_hat.resize(VA.rows(), 3);
        closest_to_A_hat.resize(VA.rows(), 3);

        for (int row = 0; row < VA.rows(); row++) {
            for (int col = 0; col < 3; col++) {
                VA_hat(row, col) = VA(row, col) - baryA(col);
                closest_to_A(row, col) = closest_to_A(row, col) - baryClosest_to_A(col);
            }
        }

        Eigen::Matrix3d A = Eigen::Matrix3d::Zero();

        for (int row = 0; row < VA_hat.rows(); row++) {
            Eigen::Vector3d rowA_hat;
            Eigen::Vector3d rowClosest_to_A_hat;

            rowA_hat(0) = VA_hat(row, 0);
            rowA_hat(1) = VA_hat(row, 1);
            rowA_hat(2) = VA_hat(row, 2);

            rowClosest_to_A_hat(0) = closest_to_A_hat(row, 0);
            rowClosest_to_A_hat(1) = closest_to_A_hat(row, 1);
            rowClosest_to_A_hat(2) = closest_to_A_hat(row, 2);

            A = A + rowA_hat * (rowClosest_to_A_hat.transpose());
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
        Eigen::Vector3d t = baryClosest_to_A - R * baryA;

        //std::cout << "Translation vector: " << std::endl << t << std::endl;

        // UPDATING SCAN ALIGNMENT
        float sum = 0;
        for (int row = 0; row < VA.rows(); row++) {
            Eigen::Vector3d rowVA;
            Eigen::Vector3d row_closest_to_A;
            rowVA(0) = VA(row, 0);
            rowVA(1) = VA(row, 1);
            rowVA(2) = VA(row, 2);

            row_closest_to_A(0) = closest_to_A(row, 0);
            row_closest_to_A(1) = closest_to_A(row, 1);
            row_closest_to_A(2) = closest_to_A(row, 2);

            rowVA = R * rowVA + t;

            VA(row, 0) = rowVA(0);
            VA(row, 1) = rowVA(1);
            VA(row, 2) = rowVA(2);

            Eigen::Vector3d diff;
            diff = row_closest_to_A - rowVA;
            sum = sum + (diff.transpose()) * diff;

        }
        error = sum/VA.rows();



        //std::cout << "The previous error was: " << oldError << std::endl;
        //std::cout << "The current error is: " << error << std::endl;
        //std::cout << "The difference (previous - current error) is: " << oldError - error << std::endl;

        ///
        ///SCAN B
        ///
        // MATCHING
        // SCAN B
        std::cout << "Computing scan B" << std::endl;
        for (int row = 0; row < VB.rows(); row++){
            std::vector<double> query_pt(3);
            for (size_t d = 0; d < 3; d++)
                query_pt[d] = VB(row, d);

            //Closest pts in the 2 other meshes
            Eigen::Matrix<double, Dynamic, Dynamic>  closest_to_query_pt(2, 3);
            //Returned indices of the closest pt in the 2 other meshes
            vector<size_t> ret_indices(2);
            //Returned index of the closest pt in the current mesh
            size_t current_ret_index;
            // Distance squared between the two points
            double out_dist_sqr;

            // Find index of closest point in mesh A
            nanoflann::KNNResultSet<double> resultSet(1);
            resultSet.init(&current_ret_index, &out_dist_sqr);
            matA_index.index->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));
            ret_indices[0] = current_ret_index;

            // Find values of closest point in mesh A
            closest_to_query_pt(0, 0) = VA(ret_indices[0], 0);
            closest_to_query_pt(0, 1) = VA(ret_indices[0], 1);
            closest_to_query_pt(0, 2) = VA(ret_indices[0], 2);


            // Find index of closest point in mesh C
            resultSet.init(&current_ret_index, &out_dist_sqr);
            matC_index.index->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));
            ret_indices[1] = current_ret_index;

            // Find values of closest point in mesh C
            closest_to_query_pt(1, 0) = VC(ret_indices[1], 0);
            closest_to_query_pt(1, 1) = VC(ret_indices[1], 1);
            closest_to_query_pt(1, 2) = VC(ret_indices[1], 2);

            //Build K-D tree to hold the closest pts in the other 2 meshes
            my_kd_tree_t  closest_to_query_pt_index(closest_to_query_pt, 10);
            closest_to_query_pt_index.index->buildIndex();
            // Find index of closest point between the 2 other closest points
            resultSet.init(&current_ret_index, &out_dist_sqr);
            closest_to_query_pt_index.index->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));

            // Find values of closest point
            closest_to_B(row, 0)= closest_to_query_pt(current_ret_index, 0);
            closest_to_B(row, 1) = closest_to_query_pt(current_ret_index, 1);
            closest_to_B(row, 2) = closest_to_query_pt(current_ret_index, 2);

        }


        // COMPUTING THE REGISTRATION

        Eigen::Vector3d baryB = calc_barycentre(VB);
        Eigen::Vector3d baryClosest_to_B = calc_barycentre(closest_to_B);

        Eigen::MatrixXd VB_hat;
        Eigen::MatrixXd closest_to_B_hat;
        VB_hat.resize(VB.rows(), 3);
        closest_to_B_hat.resize(VB.rows(), 3);

        for (int row = 0; row < VB.rows(); row++) {
            for (int col = 0; col < 3; col++) {
                VB_hat(row, col) = VB(row, col) - baryB(col);
                closest_to_B(row, col) = closest_to_B(row, col) - baryClosest_to_B(col);
            }
        }

        A = Eigen::Matrix3d::Zero();

        for (int row = 0; row < VB_hat.rows(); row++) {
            Eigen::Vector3d rowB_hat;
            Eigen::Vector3d rowClosest_to_B_hat;

            rowB_hat(0) = VB_hat(row, 0);
            rowB_hat(1) = VB_hat(row, 1);
            rowB_hat(2) = VB_hat(row, 2);

            rowClosest_to_B_hat(0) = closest_to_B_hat(row, 0);
            rowClosest_to_B_hat(1) = closest_to_B_hat(row, 1);
            rowClosest_to_B_hat(2) = closest_to_B_hat(row, 2);

            A = A + rowB_hat * (rowClosest_to_B_hat.transpose());
        }

        std::cout << "matrix: " << A << std::endl;

        // Compute singular value decomposition
        Eigen::JacobiSVD<MatrixXd> svd_B(A, ComputeThinU | ComputeThinV);

        /* // For debugging purposes
        std::cout << "Its singular values are:" << std::endl << svd.singularValues() << std::endl;
        std::cout << "Its left singular vectors are the columns of the thin U matrix:" << std::endl << svd.matrixU() << std::endl;
        std::cout << "Its right singular vectors are the columns of the thin V matrix:" << std::endl << svd.matrixV() << std::endl;
        */

        // Compute rotation matrix
        U_matrix = svd_B.matrixU();
        V_matrix = svd_B.matrixV();
        R = V_matrix * (U_matrix.transpose());

        //std::cout << "Rotation matrix: " << std::endl << R << std::endl;

        // Compute translation vector:
        t = baryClosest_to_B - R * baryB;

        //std::cout << "Translation vector: " << std::endl << t << std::endl;

        // UPDATING SCAN ALIGNMENT
        sum = 0;
        for (int row = 0; row < VB.rows(); row++) {
            Eigen::Vector3d rowVB;
            Eigen::Vector3d row_closest_to_B;
            rowVB(0) = VB(row, 0);
            rowVB(1) = VB(row, 1);
            rowVB(2) = VB(row, 2);

            row_closest_to_B(0) = closest_to_B(row, 0);
            row_closest_to_B(1) = closest_to_B(row, 1);
            row_closest_to_B(2) = closest_to_B(row, 2);

            rowVB = R * rowVB + t;

            VB(row, 0) = rowVB(0);
            VB(row, 1) = rowVB(1);
            VB(row, 2) = rowVB(2);

            Eigen::Vector3d diff;
            diff = row_closest_to_B - rowVB;
            sum = sum + (diff.transpose()) * diff;

        }
        error = sum/VB.rows();

        ///
        ///SCAN C
        ///
        // MATCHING
        // SCAN C
        std::cout << "Computing scan C" << std::endl;
        for (int row = 0; row < VC.rows(); row++){
            std::vector<double> query_pt(3);
            for (size_t d = 0; d < 3; d++)
                query_pt[d] = VC(row, d);

            //Closest pts in the 2 other meshes
            Eigen::Matrix<double, Dynamic, Dynamic>  closest_to_query_pt(2, 3);
            //Returned indices of the closest pt in the 2 other meshes
            vector<size_t> ret_indices(2);
            //Returned index of the closest pt in the current mesh
            size_t current_ret_index;
            // Distance squared between the two points
            double out_dist_sqr;

            // Find index of closest point in mesh A
            nanoflann::KNNResultSet<double> resultSet(1);
            resultSet.init(&current_ret_index, &out_dist_sqr);
            matA_index.index->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));
            ret_indices[0] = current_ret_index;

            // Find values of closest point in mesh A
            closest_to_query_pt(0, 0) = VA(ret_indices[0], 0);
            closest_to_query_pt(0, 1) = VA(ret_indices[0], 1);
            closest_to_query_pt(0, 2) = VA(ret_indices[0], 2);


            // Find index of closest point in mesh B
            resultSet.init(&current_ret_index, &out_dist_sqr);
            matB_index.index->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));
            ret_indices[1] = current_ret_index;

            // Find values of closest point in mesh B
            closest_to_query_pt(1, 0) = VB(ret_indices[1], 0);
            closest_to_query_pt(1, 1) = VB(ret_indices[1], 1);
            closest_to_query_pt(1, 2) = VB(ret_indices[1], 2);

            //Build K-D tree to hold the closest pts in the other 2 meshes
            my_kd_tree_t  closest_to_query_pt_index(closest_to_query_pt, 10);
            closest_to_query_pt_index.index->buildIndex();
            // Find index of closest point between the 2 other closest points
            resultSet.init(&current_ret_index, &out_dist_sqr);
            closest_to_query_pt_index.index->findNeighbors(resultSet, &query_pt[0], nanoflann::SearchParams(10));

            // Find values of closest point in mesh C
            closest_to_C(row, 0)= closest_to_query_pt(current_ret_index, 0);
            closest_to_C(row, 1) = closest_to_query_pt(current_ret_index, 1);
            closest_to_C(row, 2) = closest_to_query_pt(current_ret_index, 2);


        }

        // COMPUTING THE REGISTRATION

        Eigen::Vector3d baryC = calc_barycentre(VC);
        Eigen::Vector3d baryClosest_to_C = calc_barycentre(closest_to_C);

        Eigen::MatrixXd VC_hat;
        Eigen::MatrixXd closest_to_C_hat;
        VC_hat.resize(VC.rows(), 3);
        closest_to_C_hat.resize(VC.rows(), 3);

        for (int row = 0; row < VC.rows(); row++) {
            for (int col = 0; col < 3; col++) {
                VC_hat(row, col) = VC(row, col) - baryC(col);
                closest_to_C(row, col) = closest_to_C(row, col) - baryClosest_to_C(col);
            }
        }

        A = Eigen::Matrix3d::Zero();

        for (int row = 0; row < VC_hat.rows(); row++) {
            Eigen::Vector3d rowC_hat;
            Eigen::Vector3d rowClosest_to_C_hat;

            rowC_hat(0) = VC_hat(row, 0);
            rowC_hat(1) = VC_hat(row, 1);
            rowC_hat(2) = VC_hat(row, 2);

            rowClosest_to_C_hat(0) = closest_to_C_hat(row, 0);
            rowClosest_to_C_hat(1) = closest_to_C_hat(row, 1);
            rowClosest_to_C_hat(2) = closest_to_C_hat(row, 2);

            A = A + rowC_hat * (rowClosest_to_C_hat.transpose());
        }

        std::cout << "matrix: " << A << std::endl;

        // Compute singular value decomposition
        Eigen::JacobiSVD<MatrixXd> svd_C(A, ComputeThinU | ComputeThinV);

        /* // For debugging purposes
        std::cout << "Its singular values are:" << std::endl << svd.singularValues() << std::endl;
        std::cout << "Its left singular vectors are the columns of the thin U matrix:" << std::endl << svd.matrixU() << std::endl;
        std::cout << "Its right singular vectors are the columns of the thin V matrix:" << std::endl << svd.matrixV() << std::endl;
        */

        // Compute rotation matrix
        U_matrix = svd_C.matrixU();
        V_matrix = svd_C.matrixV();
        R = V_matrix * (U_matrix.transpose());

        //std::cout << "Rotation matrix: " << std::endl << R << std::endl;

        // Compute translation vector:
        t = baryClosest_to_C - R * baryC;

        //std::cout << "Translation vector: " << std::endl << t << std::endl;

        // UPDATING SCAN ALIGNMENT
        sum = 0;
        for (int row = 0; row < VC.rows(); row++) {
            Eigen::Vector3d rowVC;
            Eigen::Vector3d row_closest_to_C;
            rowVC(0) = VC(row, 0);
            rowVC(1) = VC(row, 1);
            rowVC(2) = VC(row, 2);

            row_closest_to_C(0) = closest_to_C(row, 0);
            row_closest_to_C(1) = closest_to_C(row, 1);
            row_closest_to_C(2) = closest_to_C(row, 2);

            rowVC = R * rowVC + t;

            VC(row, 0) = rowVC(0);
            VC(row, 1) = rowVC(1);
            VC(row, 2) = rowVC(2);

            Eigen::Vector3d diff;
            diff = row_closest_to_C - rowVC;
            sum = sum + (diff.transpose()) * diff;

        }
        error = sum/VC.rows();
    }

}

