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
#include "calc_barycentre.h"
#include "icp.h"

using namespace Eigen;
using namespace std;
using namespace nanoflann;

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
