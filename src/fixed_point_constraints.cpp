#include <fixed_point_constraints.h>
#include <algorithm>
void fixed_point_constraints(Eigen::SparseMatrixd &P, unsigned int q_size, const std::vector<unsigned int> indices) {
    std::vector<Eigen::Triplet<double>> tripletlist;
    unsigned int n_fixed = indices.size();
    P.setZero();
    P.resize(q_size-n_fixed*3, q_size);
    std::vector<unsigned int> free_indices;
    for(int i=0;i<q_size/3;++i)
    {
        if(std::find(indices.begin(),indices.end(),i)==indices.end())
        {
            free_indices.push_back(i);
        }
    }

    for(int r=0;r<q_size/3-n_fixed;++r)
    {
        int c = free_indices[r];
        //std::cout<<"r:"<<r<<"c:"<<c<<std::endl;
        tripletlist.emplace_back(3*r+0,3*c+0,1.0);
        tripletlist.emplace_back(3*r+1,3*c+1,1.0);
        tripletlist.emplace_back(3*r+2,3*c+2,1.0);
    }

    P.setFromTriplets(tripletlist.begin(),tripletlist.end());
    P.makeCompressed();
    //std::cout<<"fixed_point_constraint finished"<<std::endl;


    //double deter = (P*P.transpose().toDense()).determinant();
    //std::cout<<"P deter: "<<deter<<std::endl;


}