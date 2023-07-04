#include <dphi_cloth_triangle_dX.h>
//Input:
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the 1x3 vertex indices for this tetrahedron
//  X - the 3D position in the underformed space at which to compute the gradient
//Output:
//  dphi - the 3x3 gradient of the the basis functions wrt to X. The i'th row stores d phi_i/dX

//compute 3x3 deformation gradient 
void dphi_cloth_triangle_dX(Eigen::Matrix3d &dphi, Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, Eigen::Ref<const Eigen::Vector3d> X) {
    dphi.setZero();
    Eigen::Matrix32d T;
    T.col(0) = V.row(element(1)) - V.row(element(0));
    T.col(1) = V.row(element(2)) - V.row(element(0));
    dphi.block(1,0,2,3) = (T.transpose()*T).inverse()*T.transpose();
    dphi.block(0,0,1,3) = -dphi.block(2,0,1,3) - dphi.block(1,0,1,3);
    

}

