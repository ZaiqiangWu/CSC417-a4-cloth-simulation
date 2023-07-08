#include <Eigen/Dense>
#include <EigenTypes.h>

//  q - generalized coordinates for the FEM system
//  dX - the 3x3 matrix containing dphi/dX
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the vertex indices of this triangle
//  area - the area of this triangle
//  mu,lambda - material parameters for the cloth material model
//Output:
//  f - the per-triangle gradient of the membrane potential energy (the linear model described in the README).
void dV_membrane_corotational_dq(Eigen::Vector9d &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda);

Eigen::Matrix3d get_skew(Eigen::Vector3d t);

void ComputedFdq(Eigen::Matrix99d &dFdq, Eigen::Vector3d n, Eigen::Vector3d x0,Eigen::Vector3d x1,Eigen::Vector3d x2, Eigen::Vector3d N, Eigen::Ref<const Eigen::Matrix3d> dX);