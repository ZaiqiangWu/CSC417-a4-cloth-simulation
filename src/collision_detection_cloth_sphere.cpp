#include <collision_detection_cloth_sphere.h>
#include <iostream>
//  q - generalized coordinates for the FEM system
//  center - the position of the sphere center in the world space
//  radius - the radius of the sphere in the world space
//Output:
//  cloth_index - the indices of all vertices currently in contact with the sphere
//  normals - the outward facing contact normals for each contacting vertex.
void collision_detection_cloth_sphere(std::vector<unsigned int> &cloth_index, std::vector<Eigen::Vector3d> &normals, Eigen::Ref<const Eigen::VectorXd> q, Eigen::Ref<const Eigen::Vector3d> center, double radius) {

    cloth_index.clear();
    normals.clear();
    int q_size=q.rows();
    for(int i=0;i<q_size/3;++i)
    {
        Eigen::Vector3d vertex = q.block(3*i,0,3,1);
        if((vertex-center).norm()<radius)
        {
            cloth_index.push_back(i);
            Eigen::Vector3d normal=(vertex-center).normalized();
            normals.push_back(normal);
        }
    }

    
}