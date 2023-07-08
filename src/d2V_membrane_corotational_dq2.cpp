#include <d2V_membrane_corotational_dq2.h>
#include <iostream>
#include <dV_membrane_corotational_dq.h>
//  q - generalized coordinates for the FEM system
//  dX - the 3x3 matrix containing dphi/dX
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  element - the vertex indices of this triangle
//  area - the area of this triangle
//  mu,lambda - material parameters for the cloth material model
//Output:
//  H - the per-triangle Hessian of the potential energy (the linear model described in the README).
void d2V_membrane_corotational_dq2(Eigen::Matrix99d &H, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Matrix3d> dX, 
                          Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::RowVectorXi> element, double area, 
                          double mu, double lambda) {
    

    //SVD = USW^T
    Eigen::Matrix3d U;
    Eigen::Vector3d S; 
    Eigen::Matrix3d W; 
    Eigen::Matrix3d F; //deformation gradient
    
    double tol = 1e-5;
    
    //Compute SVD of F here
    Eigen::Vector3d x0=q.block(3*element(0),0,3,1);
    Eigen::Vector3d x1=q.block(3*element(1),0,3,1);
    Eigen::Vector3d x2=q.block(3*element(2),0,3,1);
    Eigen::Vector3d n = ((x1-x0).cross(x2-x0));
    n=n/n.stableNorm();
    Eigen::Vector3d X0=V.row(element(0));
    Eigen::Vector3d X1=V.row(element(1));
    Eigen::Vector3d X2=V.row(element(2));
    Eigen::Vector3d N =  ((X1-X0).cross(X2-X0));
    N = N/N.stableNorm();
    Eigen::Matrix3d x;
    x.col(0)=x0;
    x.col(1)=x1;
    x.col(2)=x2;
    F = x*dX + n*N.transpose();
    Eigen::JacobiSVD<Eigen::Matrix3d> svd(F,Eigen::ComputeFullU|Eigen::ComputeFullV);
    U = svd.matrixU();
    S = svd.singularValues();
    W = svd.matrixV();
    
    //deal with singularity in the svd gradient
    if(std::fabs(S[0] - S[1]) < tol || std::fabs(S[1] - S[2]) < tol || std::fabs(S[0] - S[2]) < tol) {
        F += Eigen::Matrix3d::Random()*tol;
        Eigen::JacobiSVD<Eigen::Matrix3d> svd2(F, Eigen::ComputeFullU | Eigen::ComputeFullV);
        U = svd2.matrixU();
        W = svd2.matrixV();
        S = svd2.singularValues();
    }
    
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

    //TODO: compute H, the hessian of the corotational energy
    Eigen::Tensor3333d dUdF;
    Eigen::Tensor333d dsdF;
    Eigen::Tensor3333d dWdF;
    dsvd(dUdF,dsdF,dWdF,F);
    Eigen::Matrix99d dFdq;
    ComputedFdq(dFdq,n,x0,x1,x2,N,dX);
    Eigen::Matrix99d d2phidF2;
    d2phidF2.setZero();
    double s0=S(0);
    double s1=S(1);
    double s2=S(2);
    double dphids0 = mu*2*(s0-1) + lambda*(s0+s1+s2-3.0);
    double dphids1 = mu*2*(s1-1) + lambda*(s0+s1+s2-3.0);
    double dphids2 = mu*2*(s2-1) + lambda*(s0+s1+s2-3.0);
    Eigen::Matrix3d dS=Eigen::Vector3d(dphids0,dphids1,dphids2).asDiagonal();
    Eigen::Matrix3d d2phids2=2*mu*Eigen::Matrix3d::Identity() + lambda*Eigen::Matrix3d::Ones();
    for(int i=0;i<3;++i)
    {
        for(int j=0;j<3;++j)
        {
            Eigen::Matrix3d dFij;
            dFij.setZero();
            Eigen::Vector3d dsij=d2phids2*dsdF[i][j];

            dFij = dUdF[i][j]*dS*W.transpose()+U*dsij.asDiagonal()*W.transpose()+U*dS*dWdF[i][j].transpose();
            for(int r=0;r<3;++r)
            {
                for(int c=0;c<3;++c)
                {
                    d2phidF2(3*j+i,3*c+r) = dFij(r,c);
                }
            }
        }
    }
    H = area * dFdq.transpose()*d2phidF2*dFdq;

    

    //fix errant eigenvalues
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix99d> es(H);
    
    Eigen::MatrixXd DiagEval = es.eigenvalues().real().asDiagonal();
    Eigen::MatrixXd Evec = es.eigenvectors().real();
    
    for (int i = 0; i < 9; ++i) {
        if (es.eigenvalues()[i]<1e-6) {
            DiagEval(i,i) = 1e-3;
        }
    }
    
    H = Evec * DiagEval * Evec.transpose();
    
}
