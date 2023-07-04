#include <dV_cloth_gravity_dq.h>
//  M - sparse mass matrix for the entire mesh
//  g - the acceleration due to gravity
//Output:
//  fg - the gradient of the gravitational potential for the entire mesh
void dV_cloth_gravity_dq(Eigen::VectorXd &fg, Eigen::SparseMatrixd &M, Eigen::Ref<const Eigen::Vector3d> g) {
    int q_size = M.rows();
    Eigen::VectorXd g_g;
    g_g.resize(q_size);
    for(int i=0;i<q_size;++i)
    {
        g_g(i) = g(i%3);
    }
    fg = -g_g.transpose()*M;
    
}
