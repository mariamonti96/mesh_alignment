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
#include "calc_barycentre.h"
#include "point_to_plane.h"

using namespace Eigen;
using namespace std;


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


