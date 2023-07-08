#include <mass_matrix_mesh.h>
#include <vector>
#include<iostream>
//Input:
//  q - generalized coordinates for the FEM system
//  V - the nx3 matrix of undeformed vertex positions
//  F - the mx3 matrix of triangle-vertex indices
//  density - the density of the cloth material
//  areas - the mx1 vector of undeformed triangle areas
//Output:
//  M - sparse mass matrix for the entire mesh
void mass_matrix_mesh(Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::VectorXd> q, 
                         Eigen::Ref<const Eigen::MatrixXd> V, Eigen::Ref<const Eigen::MatrixXi> F,
                         double density, Eigen::Ref<const Eigen::VectorXd> areas) {
    int n_triangles = F.rows();
    std::vector<Eigen::Triplet<double>> TripletList;

    int q_size = q.rows();
    M.resize(q_size,q_size);
    M.setZero();
    for(int i=0;i<n_triangles;++i)
    {
        Eigen::RowVectorXi element = F.row(i);
        //std::cout<<"element:"<<element<<std::endl;
        double area = areas(i);
        for(int vid=0;vid<3;++vid)
        {
            for(int vid2=0;vid2<3;++vid2)
            {
                double item =density*area/12.0;
                if(vid==vid2)
                {
                    item = density*area/6.0;
                }
                for(int dim=0;dim<3;++dim)
                {
                    TripletList.emplace_back(3*element(vid)+dim,3*element(vid2)+dim,item);
                }
            }
        }

    }
    /*for(int i=0;i<q_size;++i)
    {
        TripletList.emplace_back(i,i,0.1);
    }*/
    M.setFromTriplets(TripletList.begin(),TripletList.end());
    M.makeCompressed();

    //std::cout<<"det(M):"<<M.toDense().determinant()<<std::endl;//wrong!!
   
}
 