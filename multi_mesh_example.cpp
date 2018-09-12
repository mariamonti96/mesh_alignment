#include <igl/read_triangle_mesh.h>
#include <igl/viewer/Viewer.h>

int main(int argc, char * argv[])
{
    Eigen::MatrixXd VA,VB;
    Eigen::MatrixXi FA,FB;
    igl::read_triangle_mesh(argv[1],VA,FA);
    igl::read_triangle_mesh(argv[2],VB,FB);

    // Concatenate (VA,FA) and (VB,FB) into (V,F)
    Eigen::MatrixXd V(VA.rows()+VB.rows(),VA.cols());
    V<<VA,VB;
    Eigen::MatrixXi F(FA.rows()+FB.rows(),FA.cols());
    F<<FA,(FB.array()+VA.rows());

    // blue color for faces of first mesh, orange for second
    Eigen::MatrixXd C(F.rows(),3);
    C<<
     Eigen::RowVector3d(0.2,0.3,0.8).replicate(FA.rows(),1),
            Eigen::RowVector3d(1.0,0.7,0.2).replicate(FB.rows(),1);

    igl::viewer::Viewer viewer;
    viewer.data.set_mesh(V,F);
    viewer.data.set_colors(C);
    viewer.data.set_face_based(true);
    viewer.launch();
}
