#include <velocity_filter_cloth_sphere.h>
//Input:
//  qdot - the 3nx1 generalized velocities of the cloth mesh
//  index - a list of collision vertex indices from the collision detector
//  normals - a list of collision normals from the collision detector
//Output:
//  qdot- the filtered 3nx1 generalized velocities
void velocity_filter_cloth_sphere(Eigen::VectorXd &qdot, const std::vector<unsigned int> &indices, 
                                  const std::vector<Eigen::Vector3d> &normals) {
    for(int i=0;i<indices.size();++i)
    {
        unsigned int vid=indices[i];
        Eigen::Vector3d normal=normals[i];
        Eigen::Vector3d velocity = qdot.block(3*vid,0,3,1);
        Eigen::Vector3d v_n,v_t,v_new;
        v_n = velocity.dot(normal)*normal;
        v_t = velocity-v_n;
        v_new = v_t;
        qdot.block(3*vid,0,3,1)=v_new;
    }
    

}