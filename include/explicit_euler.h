//
// Created by Zaiqiang Wu on 2023/07/07.
//

#ifndef A4_CLOTH_SIMULATION_EXPLICIT_EULER_H
#define A4_CLOTH_SIMULATION_EXPLICIT_EULER_H
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include<Eigen/SparseCholesky>
#include <EigenTypes.h>
#include <cassert>
//Input:
//  q - generalized coordinates for the FEM system
//  qdot - generalized velocity for the FEM system
//  dt - the time step in seconds
//  mass - the mass matrix
//  force(f, q, qdot) - a function that computes the force acting on the FEM system. This takes q and qdot as parameters, returns the force in f.
//  stiffness(K, q, qdot) - a function that computes the stiffness (negative second derivative of the potential energy). This takes q and qdot as parameters, returns the stiffness matrix in K.
//  tmp_force - scratch space to collect forces
//  tmp_stiffness - scratch space to collect stiffness matrix
//Output:
//  q - set q to the updated generalized coordinate using linearly implicit time integration
//  qdot - set qdot to the updated generalized velocity using linearly implicit time integration
template<typename FORCE>
inline void explicit_euler(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt,
                                    const Eigen::SparseMatrixd &mass,  FORCE &force,
                                    Eigen::VectorXd &tmp_force) {
    force(tmp_force,q,qdot);
    assert(!tmp_force.hasNaN());
    //dt=dt*0.1;

    const Eigen::SparseMatrixd& A=mass;
    //std::cout<<mass.toDense().determinant()<<std::endl;

    Eigen::SimplicialLDLT<Eigen::SparseMatrixd> ldlt;
    //Eigen::ConjugateGradient<Eigen::SparseMatrixd> ldlt;
    ldlt.compute(A);
    auto qdot_new = ldlt.solve(mass*qdot+dt*tmp_force);
    assert(!qdot_new.hasNaN());
    auto q_new = q+dt*qdot_new;
    q=q_new;
    qdot=qdot_new;
    assert(!q.hasNaN());
    assert(!qdot.hasNaN());

}
template<typename FORCE>
inline void energy_minimization(Eigen::VectorXd &q, Eigen::VectorXd &qdot, double dt,
                           const Eigen::SparseMatrixd &mass,  FORCE &force,
                           Eigen::VectorXd &tmp_force) {
    force(tmp_force,q,qdot);
    assert(!tmp_force.hasNaN());
    //dt=dt*0.1;

    qdot = tmp_force*0.001;
    //q+=dt*tmp_force*1.0;

}
#endif //A4_CLOTH_SIMULATION_EXPLICIT_EULER_H
