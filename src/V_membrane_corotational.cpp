#include <V_membrane_corotational.h>
//  q - generalized coordinates for the FEM system
//  dX - the 3x3 matrix containing dphi/dX
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the vertex indices of this triangle
//  area - the area of this triangle
//  mu,lambda - material parameters for the cloth material model
//Output:
//  energy- the per-triangle potential energy (the linear model described in the README).

//Allowed to use libigl SVD or Eigen SVD for this part
void V_membrane_corotational(double &energy, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda) {
    Eigen::Vector3d x0=q.block(3*element(0),0,3,1);
    Eigen::Vector3d x1=q.block(3*element(1),0,3,1);
    Eigen::Vector3d x2=q.block(3*element(2),0,3,1);
    Eigen::Vector3d n = ((x2-x1).cross(x0-x2)).normalized();
    Eigen::Vector3d X0=V.row(element(0));
    Eigen::Vector3d X1=V.row(element(1));
    Eigen::Vector3d X2=V.row(element(2));
    Eigen::Vector3d N =  ((X2-X1).cross(X0-X2)).normalized();
    Eigen::Matrix3d x;
    x.col(0)=x0;
    x.col(1)=x1;
    x.col(2)=x2;
    Eigen::Matrix3d F = x*dX + n*N.transpose();
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(F,Eigen::ComputeFullU|Eigen::ComputeFullV);
    Eigen::Vector3d s = svd.singularValues();
    energy =0.0;
    for(int i=0;i<3;++i)
    {
        energy += area*mu*(s(i)-1)*(s(i)-1);
    }
    energy += 0.5*area*lambda*(s(0)+s(1)+s(2)-3)*(s(0)+s(1)+s(2)-3);


    

}
