#include <assemble_stiffness.h>
#include <iostream>
void assemble_stiffness(Eigen::SparseMatrixd &K, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::VectorXd> qdot, Eigen::Ref<const Eigen::MatrixXd> dX,
                     Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F, Eigen::Ref<const Eigen::VectorXd> a0, 
                     double mu, double lambda) {
    int n_triangles = F.rows();
    int q_size = q.rows();
    //K=-d2Vdq2
    K.setZero();
    K.resize(q_size,q_size);
    std::vector<Eigen::Triplet<double>> TripletList;
    for(int i=0;i<n_triangles;++i)
    {
        Eigen::VectorXi element = F.row(i);
        Eigen::Matrix99d single_stiffness=Eigen::Matrix99d::Zero();
        Eigen::Matrix<double, 1,9> tmp_row;
        tmp_row = dX.row(i); //ei is the triangle index.
        Eigen::Map<const Eigen::Matrix3d> new_dX(tmp_row.data());
        d2V_membrane_corotational_dq2(single_stiffness,q,new_dX,V,element,a0(i),mu,lambda);
        single_stiffness*=-1.0;//K=-H
        for(int vid=0;vid<3;++vid)
        {
            for(int dim=0;dim<3;++dim)
            {
                for(int vid2=0;vid2<3;++vid2)
                {
                    for(int dim2=0;dim2<3;++dim2)
                    {
                        TripletList.emplace_back(element(vid)*3+dim, element(vid2)*3+dim2, single_stiffness(vid*3+dim,vid2*3+dim2));
                    }
                }

            }
        }
    }
    K.setFromTriplets(TripletList.begin(),TripletList.end());
    K.makeCompressed();
       
        
    };
