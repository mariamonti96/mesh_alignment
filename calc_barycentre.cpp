#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
#include "calc_barycentre.h"

using namespace Eigen;
using namespace std;

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
