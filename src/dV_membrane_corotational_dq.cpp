#include <dV_membrane_corotational_dq.h>
#include <iostream>

Eigen::Matrix3d get_skew(Eigen::Vector3d t)
{
    Eigen::Matrix3d t_hat;
    t_hat << 0, -t(2), t(1),
            t(2), 0, -t(0),
            -t(1), t(0), 0;
    return t_hat;
}

//  q - generalized coordinates for the FEM system
//  dX - the 3x3 matrix containing dphi/dX
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the vertex indices of this triangle
//  area - the area of this triangle
//  mu,lambda - material parameters for the cloth material model
//Output:
//  f - the per-triangle gradient of the membrane potential energy (the linear model described in the README).
void dV_membrane_corotational_dq(Eigen::Vector9d &dV, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda) {

    //Deformation Gradient
    Eigen::Matrix3d dx; //deformed tangent matrix 
    Eigen::Matrix3d U;
    Eigen::Vector3d S; 
    Eigen::Matrix3d W; 

    //TODO: SVD Here
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
    U = svd.matrixU();
    S = svd.singularValues();
    W = svd.matrixV();

    //Fix for inverted elements (thanks to Danny Kaufman)
    double det = S[0]*S[1];
    
     if(det <= -1e-10)
    {
        if(S[0] < 0) S[0] *= -1;
        if(S[1] < 0) S[1] *= -1;
        if(S[2] < 0) S[2] *= -1;
    }
    
    if(U.determinant() <= 0)
    {
        U(0, 2) *= -1;
        U(1, 2) *= -1;
        U(2, 2) *= -1;
    }
    
    if(W.determinant() <= 0)
    {
        W(0, 2) *= -1;
        W(1, 2) *= -1;
        W(2, 2) *= -1;
    }
    
    //TODO: energy model gradient
    double s0=S(0);
    double s1=S(1);
    double s2=S(2);
    double dphids0 = mu*2*(s0-1) + lambda*(s0+s1+s2-3.0);
    double dphids1 = mu*2*(s1-1) + lambda*(s0+s1+s2-3.0);
    double dphids2 = mu*2*(s2-1) + lambda*(s0+s1+s2-3.0);
    Eigen::Matrix3d dphidF = U*Eigen::Vector3d(dphids0,dphids1,dphids2).asDiagonal()*W.transpose();

    Eigen::Matrix3d dndx0=(Eigen::Matrix3d::Identity()/ sqrt(2*area)-4*area*area*n*n.transpose())* get_skew(x2-x1);
    Eigen::Matrix3d dndx1=(Eigen::Matrix3d::Identity()/ sqrt(2*area)-4*area*area*n*n.transpose())* get_skew(x0-x2);
    Eigen::Matrix3d dndx2=(Eigen::Matrix3d::Identity()/ sqrt(2*area)-4*area*area*n*n.transpose())* get_skew(x1-x0);
    Eigen::Matrix3d dndx[3] = {dndx0,dndx1,dndx2};

    Eigen::Matrix99d dFdq=Eigen::Matrix99d::Zero();
    for(int i=0;i<3;++i)
    {
        for(int j=0;j<3;++j)
        {
            for(int a=0;a<3;++a)
            {
                for(int b=0;b<3;++b)
                {
                    dFdq(3*j+i,3*a+b)=dndx[a](i,b)*N(j);
                    if(b==i)
                    {
                        dFdq(3*j+i,3*a+b)+=dX(a,j);
                    }
                }
            }
        }
    }
    Eigen::VectorXd dphidF_vec(Eigen::Map<Eigen::VectorXd>(dphidF.data(),9));
    dV = area* (dphidF_vec.transpose()*dFdq).transpose();
}
