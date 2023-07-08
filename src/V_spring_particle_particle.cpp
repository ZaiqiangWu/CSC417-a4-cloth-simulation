#include <V_spring_particle_particle.h>

//the potential energy of a spring with 3D end points q0 and qd and undeformed length l0

//Input:
//  q0 - the generalized coordinates of the first node of the spring
//  q1 - the generalized coordinates of the second node of the spring
//  l0 - the undeformed length of the spring
//  stiffness - the stiffness constant for this spring
//Output:
//  V - potential energy for the spring
void V_spring_particle_particle(double &V, Eigen ::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d> q1, double l0, double stiffness) {

    double l=(q1-q0).norm();
    V=0.5*stiffness*(l-l0)*(l-l0);

}