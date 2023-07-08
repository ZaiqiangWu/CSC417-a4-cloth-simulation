#include <dV_spring_particle_particle_dq.h>
#include <iostream>
using namespace std;
void dV_spring_particle_particle_dq(Eigen::Ref<Eigen::Vector6d> f, Eigen::Ref<const Eigen::Vector3d> q0,  Eigen::Ref<const Eigen::Vector3d>     q1, double l0, double stiffness) {

    Eigen::Vector3d f0, f1;
    auto q01=q0-q1;
    double l=q01.norm();
    f0 = -stiffness*(q0-q1).normalized()*(l-l0);
    f1=-f0;
    f.block(0,0,3,1) = f0;
    f.block(3,0,3,1) = f1;
    //cout<<q01.hasNaN()<<endl;
    ;

}