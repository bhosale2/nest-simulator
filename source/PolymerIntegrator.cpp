/*
 * PolymerIntegrator.cpp
 *
 *  Created on: Jun 20, 2014
 *      Author: mgazzola
 */

#include "PolymerIntegrator.h"

void PolymerIntegrator::_computeForces(const REAL time, const int step) {
    const unsigned int numRods = rodptrs.size();

    // rods
    for (unsigned int j = 0; j < numRods; j++)
        rodptrs[j]->reset();

    // rods
    for (unsigned int j = 0; j < numRods; j++) {
        rodptrs[j]
            ->computeAllInternalResultingForcesAndTorques();  // (in-place)

        efptrs[j]->applyForces(*rodptrs[j], time);
    }

    // rods + piston
    interptrs[0]->applyForces(
        time);  // make sure that interaction forces are in-place

    // Compute Rod-Rod collision force
    ecptrs[0]->RodRodCollision_Simple(time);
    /* ecptrs[0]->RodRodCollision_flex(time, true); */

    if (time > 0.5)
        for (unsigned int j = 0; j < numRods; j++) {
            rodptrs[j]->applyKineticFrictions();
            rodptrs[j]->applyStaticFrictions();
        }

    for (unsigned int j = 0; j < numRods; j++) {
        // Sum up internal and external forces
        // Next line is equivalent (but in-place) to: rodptrs[j]->totalForces =
        // rodptrs[j]->totalInternalForces + rodptrs[j]->externalForces;
        v_a_plus_b_equal_c(rodptrs[j]->totalInternalForces,
                           rodptrs[j]->externalForces, rodptrs[j]->totalForces);

        // Sum up internal and external torques
        // Next line is equivalent (but in-place) to: rodptrs[j]->totalTorques =
        // rodptrs[j]->totalInternalTorques + rodptrs[j]->externalTorques;
        v_a_plus_b_equal_c(rodptrs[j]->totalInternalTorques,
                           rodptrs[j]->externalTorques,
                           rodptrs[j]->totalTorques);
    }
}

void PolymerIntegrator::_computeAccelerations(const REAL time, const int step) {
    const REAL dt = time - time_old;
    assert(dt > Tolerance::tol());
    time_old = time;

    _computeForces(time, step);  // in-place

    // rods
    for (unsigned int j = 0; j < rodptrs.size(); j++) {
        // Compute linear accelerations given the forces (internal+ezternal)
        // acting on the single rod elements
        v_a_divide_b_equal_c(
            rodptrs[j]->totalForces, rodptrs[j]->m,
            rodptrs[j]->a);  // rodptrs[j]->a = rodptrs[j]->totalForces /
                             // rodptrs[j]->m; Note that m = rho*A0*dS = const

        // Compute angullar accelerations given the toques (internal+ezternal)
        // acting on the single rod elements
        rodptrs[j]->deldt = (rodptrs[j]->e - rodptrs[j]->e_old) / dt;
        rodptrs[j]->wDot =
            rodptrs[j]->e * (rodptrs[j]->J0inv * rodptrs[j]->totalTorques) +
            rodptrs[j]->w * rodptrs[j]->deldt / rodptrs[j]->e;

        // v_a_times_b_equal_c(rodptrs[j]->J0inv, rodptrs[j]->totalTorques,
        // rodptrs[j]->wDot); // rodptrs[j]->wDot = rodptrs[j]->I0inv *
        // rodptrs[j]->totalTorques;
    }
}

void PolymerIntegrator::_v_update_Q(const REAL coeffDt, Rod* rod) {
    // Next lines are equivalent (but in-place) to: rodptrs[j]->Q =  vExp(0.5 *
    // dt * rodptrs[j]->w) * rodptrs[j]->Q;
    v_a_times_b_equal_c(rod->w, coeffDt, rod->tempVV3_n);  // in-place
    vExp(rod->tempVV3_n, rod->tempVM3_n);                  // in-place
    v_timesequal(rod->tempVM3_n, rod->Q);                  // in-place
    // cout << "ReachQ"<< endl;
}

void PolymerIntegrator::_v_update_x(const REAL coeffDt, Rod* rod) {
    // Next lines are equivalent (but in-place) to:  rodptrs[j]->x += coeffDt *
    // rodptrs[j]->v;

    v_plusequal_a_times_b(coeffDt, rod->v, rod->x);
}

void PolymerIntegrator::_v_update_v(const REAL coeffDt, Rod* rod) {
    // Next line is equivalent (but in-place) to: rodptrs[j]->v += coeffDt *
    // rodptrs[j]->a;
    v_plusequal_a_times_b(coeffDt, rod->a, rod->v);
}

void PolymerIntegrator::_v_update_w(const REAL coeffDt, Rod* rod) {
    // Next line is equivalent (but in-place) to: rodptrs[j]->w += coeffDt *
    // rodptrs[j]->wDot;
    v_plusequal_a_times_b(coeffDt, rod->wDot, rod->w);
    // cout << "Reach3"<< endl;
}
