# Mesh Alignment
Mesh Alignment using the Iterative Closest Points (ICP) Algorithm

The aim of this project is to implement the Iterative Closest Points (ICP) algorithm for the alignment of meshes, which was originally developed by Chen and Medioni [2] and Besl and McKay [1].

It provides a method for accurate and computationally efficient registration of 3D shapes and, in particuar, it always converges monotonically to the nearest local minimum of a mean-square distance metric. 

Since the original development, many variants have been proposed, which affect different stages of ICP. In this work, we examine some of the variants and evaluate the performance of the algorithm in different cases.

The libraries used are the following:

- nanoflann

- libigl

- Eigen

In order to be able to run the code, these libraries should be installed.

When running the code, the libigl viewer will open. To run any of the tasks, the corresponding number should be pressed on the keyboard:

- '1': Visualisation of two models (M1 and M2) and alignment of M2 with M1 using point-to-point ICP. 

- '2': Progressive rotation of model M1, setting M2= R(M1), alignment of M2 with M1 and evaluation of the convergence behaviour. 

- '3': Progressive perturbation of model M2 by adding white Gaussian noise and evaluation of the performance of the ICP algorithm.

- '6': Point-to-plane ICP. Currently not working.

To know more about the ICP algorithm:

[1] P. J. Besl and N. D. McKay. A method for registration of 3-d shapes. IEEE Transac-
tions on Pattern Analysis and Machine Intelligence, 14(2):239{256, Feb 1992.

[2] Y. Chen and G. Medioni. Object modeling by registration of multiple range images. In
Proceedings. 1991 IEEE International Conference on Robotics and Automation, pages
2724{2729 vol.3, Apr 1991.

[3] Kok-Lim Low. Linear least-squares optimization for point-to-plane icp surface regis-
tration. 01 2004.

[4] S. Rusinkiewicz and M. Levoy. Ecient variants of the icp algorithm. In Proceedings
Third International Conference on 3-D Digital Imaging and Modeling, pages 145{152,
2001.
