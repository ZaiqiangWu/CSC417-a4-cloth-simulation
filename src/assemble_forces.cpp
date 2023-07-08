#include <assemble_forces.h>
#include <iostream>
#include<cassert>
//Input:
//  q - generalized coordinates for the FEM system
//  qdot - generalized velocity for the FEM system
//  dX - an mx9 matrix which stores the flattened dphi/dX matrices for each tetrahedron.
//       Convert this values back to 3x3 matrices using the following code (NOTE YOU MUST USE THE TEMPORARY VARIABLE tmp_row):
//       Eigen::Matrix<double, 1,9> tmp_row
//       tmp_row = dX.row(ei); //ei is the triangle index.
//       Eigen::Map<const Eigen::Matrix3d>(tmp_row.data())
//  V - the nx3 matrix of undeformed vertex positions. Each row is a single undeformed vertex position.
//  F - the mx3 triangle connectivity matrix. Each row contains to indices into V that indicate a spring between those vertices.
//  a0 - the mx1 vector of undeformed triangle areas
//  mu,lambda - material parameters for the cloth material model
//Output:
//  f - the vector 3xn vector of forces acting on each node of the mass-spring system
void assemble_forces(Eigen::VectorXd &f, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::MatrixXd> qdot, Eigen::Ref<const Eigen::MatrixXd> dX,
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, Eigen::Ref<const Eigen::VectorXd> a0,
                     double mu, double lambda) {

    int n_triangles = F.rows();
    int q_size = q.rows();

    f.resize(q_size);
    f.setZero();
    for(int i=0;i<n_triangles;++i)
    {
        Eigen::Vector9d single_force;
        Eigen::RowVectorXi element = F.row(i);
        Eigen::Matrix<double, 1,9> tmp_row;
        tmp_row = dX.row(i); //ei is the triangle index.
        Eigen::Map<const Eigen::Matrix3d> new_dX(tmp_row.data());
        //std::cout<<tmp_row.rows()<<std::endl;
        dV_membrane_corotational_dq(single_force,q,new_dX,V,element,a0(i),mu,lambda);
        if (isnan(single_force.norm())||isinf(single_force.norm()) )
            std::cout<<"single force norm:"<<single_force.norm()<<std::endl;
        assert(!single_force.hasNaN());

        single_force*=-1.0;
        for(int vid=0;vid<3;++vid)
        {
            for(int dim=0;dim<3;++dim)
            {
                assert(!single_force.hasNaN());
                f(element(vid)*3+dim)+=single_force(vid*3+dim);
                assert(!f.hasNaN());
            }
        }
        for(int j=0;j<f.rows();++j)
        {
            if(isnan(f(j)))
                std::cout<<"found nan:"<<j<<f(j)<<std::endl;

        }
        if(f.hasNaN())
        {
            for(int vid=0;vid<3;++vid)
            {
                std::cout<<element(vid)*3<<std::endl;

            }

        }

        assert(!f.hasNaN());
    }


        
       
    };
