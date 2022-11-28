/*
Detect Rod-Rod collision
Detect Rod-Object collision
*/
#ifndef EXTERNALCONTACT_H_
#define EXTERNALCONTACT_H_
#include "GeometryFunctions.h"
#include "Matrix3.h"
#include "Rod.h"
#include "UsualHeaders.h"
#include "Vector3.h"
using namespace std;
class ExternalContact {
   protected:
    typedef std::vector<Rod*> Vrodptr;
    Vrodptr rodptrs;
    const REAL isotropicFCK;
    const REAL isotropicFCS;
    vector<pair<int, int> > attachpoint;
    REAL angle_old = 0.0;
    REAL angle_new = 0.0;
    REAL angleVelocity = 0.0;
    REAL ml_old = 0.0;
    REAL ml_new = 0.0;
    REAL mVelocity = 0.0;
    REAL angleVelocity_old = 0.0;
    Vector3 externalforce = Vector3(0.0, 0.0, 0.0);
    REAL Localtime = -0.1;
    int SwitchFlag = 0;
    int timecountL = 0;
    int timecountT = 0;
    int nor;
    REAL time_old = 0.0;
    const REAL _linear(const REAL _vel, const REAL _velThresold) const {
        const REAL vel = fabs(_vel);
        const REAL velThresold = fabs(_velThresold);
        const REAL width = 0.5 * velThresold;
        const REAL velDiff = vel - velThresold;
        REAL f = 1.0;
        if (vel > (velThresold))
            f = fabs(1.0 - min(1.0, velDiff / width));
        if (vel > (velThresold + width))
            f = 0.0;
        return f;
    }

   public:
    ExternalContact(Vrodptr rodptrs,
                    const REAL _isotropicFCK,
                    const REAL _isotropicFCS,
                    vector<pair<int, int> > _attachpoint)
        : rodptrs(rodptrs),
          isotropicFCK(_isotropicFCK),
          isotropicFCS(_isotropicFCS),
          attachpoint(_attachpoint) {}
    ~ExternalContact(){};

    /* sets default near list for rods as the entire list */
    void setRodNearList() {
        for (auto ri = 0; ri < (rodptrs.size() - 1); ri++)
            for (auto rj = 0; rj < rodptrs.size(); rj++)
                rodptrs[ri]->nearbyRods.push_back(rj);
    }
    /* collision detection for flexible rods */
    /* near list can be turned on by flag*/
    void RodRodCollision_flex(const REAL time, const bool opt) {
        const REAL zetasoft = 1000;
        const REAL zeta_nu_soft = 50;
        const REAL zetahard = 1e5;
        const REAL zeta_nu_hard = 1000;
        const REAL timeThresold = 0.2;
        const REAL vStatic = 1e-4;
        const int near_tol =
            4;  // multiples of sum_radii for near list detection
        const int numRods = rodptrs.size();
        const REAL zeta = (time >= timeThresold) ? zetahard : zetasoft;
        const REAL zeta_nu =
            (time >= timeThresold) ? zeta_nu_hard : zeta_nu_soft;
        // set cord no, contact points to zero before starting to check
        rodptrs[0]->cord = 0.0;
        rodptrs[0]->intercontacts.clear();
        rodptrs[0]->intercontacts_index.clear();
        rodptrs[0]->intercontacts_load.clear();
        /* check for near list update after timer */
        if (opt && rodptrs[0]->timer == 1) {
            for (auto ri = 0; ri < (numRods - 1); ri++) {
                rodptrs[ri]->nearbyRods.clear();
                for (auto rj = ri + 1; rj < numRods; rj++) {
                    const int sizei = rodptrs[ri]->n;
                    const int sizej = rodptrs[rj]->n;
                    for (auto bri = 0; bri < sizei; bri++) {
                        for (auto brj = 0; brj < sizej; brj++) {
                            const Vector3 x1 = rodptrs[ri]->x[bri];
                            const Vector3 x2 = rodptrs[rj]->x[brj];
                            /* const REAL l1 = rodptrs[ri]->l[bri]; */
                            /* const REAL l2 = rodptrs[rj]->l[brj]; */
                            const REAL sum_r1r2 =
                                rodptrs[ri]->r[bri] + rodptrs[rj]->r[brj];
                            /* check distance set as 5 * r1+r2 */
                            if ((x1 - x2).length() <=
                                near_tol * sum_r1r2 + (sum_r1r2)) {
                                rodptrs[ri]->nearbyRods.push_back(rj);
                                bri = sizei;
                                brj = sizej;
                            }
                        }
                    }
                }
            }
            rodptrs[0]->timer = 0;
        }
        /* only check for near list interaction */
        for (unsigned int ri = 0; ri < (numRods - 1); ri++) {
            for (unsigned int rj : rodptrs[ri]->nearbyRods) {
                const int sizei = rodptrs[ri]->n;
                const int sizej = rodptrs[rj]->n;
                for (auto bri = 0; bri < sizei; bri++) {
                    for (auto brj = 0; brj < sizej; brj++) {
                        const Vector3 x1 = rodptrs[ri]->x[bri];
                        const Vector3 x2 = rodptrs[rj]->x[brj];
                        const REAL l1 = rodptrs[ri]->l[bri];
                        const REAL l2 = rodptrs[rj]->l[brj];
                        const REAL sum_r1r2 =
                            rodptrs[ri]->r[bri] + rodptrs[rj]->r[brj];
                        if ((x1 - x2).length() <= (sum_r1r2 + l1 + l2)) {
                            const Vector3 edge1 =
                                rodptrs[ri]->x[bri + 1] - rodptrs[ri]->x[bri];
                            const Vector3 edge2 =
                                rodptrs[rj]->x[brj + 1] - rodptrs[rj]->x[brj];
                            const vector<Vector3> minVectors =
                                findMinDistVectors(x1, edge1, x2, edge2);
                            const Vector3 dVector =
                                minVectors[1];  // pointing to J
                            const Vector3 dVectorDirection = dVector.unitize();
                            const Vector3 pVector =
                                minVectors[0];  // the position on i
                            // gamma tells you whether ther is overlap
                            const REAL gamma = (sum_r1r2 - dVector.length());
                            const bool yesCollision = (gamma > 0.0);
                            if (yesCollision) {
                                // Compute contact force
                                const REAL cForceContact = zeta * gamma;
                                // Compute damping force
                                const Vector3 vInteri =
                                    0.5 * (rodptrs[ri]->v[bri] +
                                           rodptrs[ri]->v[bri + 1]);
                                const Vector3 vInterj =
                                    0.5 * (rodptrs[rj]->v[brj] +
                                           rodptrs[rj]->v[brj + 1]);
                                const Vector3 vInterpenetration =
                                    vInterj - vInteri;
                                const REAL vNorm =
                                    vInterpenetration % dVectorDirection;
                                const REAL cForceDamping = -zeta_nu * vNorm;
                                const Vector3 cForce =
                                    yesCollision *
                                    (cForceContact + cForceDamping) *
                                    dVectorDirection;
                                // Remenber that the first and last points have
                                // half the mass! Hence the strange coefficient
                                // in front of cForce
                                rodptrs[ri]->externalForces[bri] -=
                                    ((bri == 0) ? 0.5 : 1.0) * cForce;
                                rodptrs[ri]->externalForces[bri + 1] -=
                                    (((bri + 1) == sizei) ? 0.5 : 1.0) * cForce;
                                rodptrs[rj]->externalForces[brj] +=
                                    ((brj == 0) ? 0.5 : 1.0) * cForce;
                                rodptrs[rj]->externalForces[brj + 1] +=
                                    (((brj + 1) == sizej) ? 0.5 : 1.0) * cForce;
                                const REAL cForcem = cForce.length();
                                /* push contact data */
                                rodptrs[0]->intercontacts.push_back(
                                    pVector + dVector / 2);
                                rodptrs[0]->intercontacts_load.push_back(
                                    cForcem);
                                rodptrs[0]->cord += 2;
                                rodptrs[0]->intercontacts_index.push_back(
                                    {ri, rj});
                                /* only kinetic sliding friction implemented */
                                const Vector3 slipVel =
                                    vInterpenetration -
                                    (vNorm)*dVectorDirection;
                                const REAL vel = fabs(slipVel.length());
                                const REAL f = _linear(vel, vStatic);
                                const Vector3 kineticFrictionForce =
                                    -(1.0 - f) * isotropicFCK * cForcem *
                                    slipVel.unitize();
                                rodptrs[ri]->kineticFrictionsForce[bri] -=
                                    ((bri == 0) ? 0.5 : 1.0) *
                                    kineticFrictionForce;
                                rodptrs[ri]->kineticFrictionsForce[bri + 1] -=
                                    ((bri + 1 == sizei) ? 0.5 : 1.0) *
                                    kineticFrictionForce;
                                rodptrs[rj]->kineticFrictionsForce[brj] +=
                                    ((brj == 0) ? 0.5 : 1.0) *
                                    kineticFrictionForce;
                                rodptrs[rj]->kineticFrictionsForce[brj + 1] +=
                                    ((brj + 1 == sizej) ? 0.5 : 1.0) *
                                    kineticFrictionForce;
                            }  // collision
                        }      // near elements
                    }          // brj
                }              // bri
            }                  // j
        }                      // i
        rodptrs[0]->cord /= numRods;
    }  // function
    // Nest----Nest----Nest----Nest----Nest----Details:
    // near list maintained and updated every 0.1 sec, beware of hidden
    // timescale!!!! if rod considerably straight, treated as 1 element or
    // at max 2 if bent Caution!! ---- fails for flexible rods
    // assumes unifom radii of rods and same number of elements in each rod
    void RodRodCollision_Simple(const REAL time) {
        const REAL zetasoft = 1500;
        const REAL zeta_nu_soft = 100;
        const REAL zetahard = 10e5;
        const REAL zeta_nu_hard = 1e3;
        /* doesnt work */
        /* const REAL zetasoft = 15000; */
        /* const REAL zeta_nu_soft = 1e3; */
        /* const REAL zetahard = 1e6; */
        /* const REAL zeta_nu_hard = 1e4; */
        const REAL timeThresold = 0.2;
        /* const REAL timeThresold = 0.0; */
        const REAL vStatic = 1e-6;
        const int numRods = rodptrs.size();
        const int size = rodptrs[0]->n;
        const REAL zeta = (time >= timeThresold) ? zetahard : zetasoft;
        const REAL zeta_nu =
            (time >= timeThresold) ? zeta_nu_hard : zeta_nu_soft;
        /* assuming Poisson ratio = 0.3 */
        const REAL zeta_tang =
            (time >= timeThresold) ? 0.8 * zetahard : 0.8 * zetasoft;
        const REAL zeta_nu_tang =
            (time >= timeThresold) ? 1.0 * zeta_nu_hard : 1.0 * zeta_nu_soft;
        // loop over every rod after 0.1 sec
        if (rodptrs[0]->timer == 1) {
            for (unsigned int ri = 0; ri < (numRods - 1); ri++) {
                rodptrs[ri]->nearbyRods.clear();
                const REAL radius_i = rodptrs[ri]->r[0];
                for (unsigned int rj = ri + 1; rj < numRods; rj++) {
                    // part added for bent rods
                    const REAL radius_j = rodptrs[rj]->r[0];
                    const Vector3 rodi =
                        rodptrs[ri]->x[size] - rodptrs[ri]->x[0];
                    const Vector3 rodj =
                        rodptrs[rj]->x[size] - rodptrs[rj]->x[0];
                    const Vector3 straightrodi =
                        (rodptrs[ri]->x[1] - rodptrs[ri]->x[0]).unitize();
                    const Vector3 straightrodj =
                        (rodptrs[rj]->x[1] - rodptrs[rj]->x[0]).unitize();
                    const REAL bendrodi = (rodi * straightrodi).length();
                    const REAL bendrodj = (rodj * straightrodj).length();
                    int tempsize_i = size;
                    int tempsize_j = size;
                    if (bendrodi > 4 * radius_i)
                        tempsize_i = size / 2;
                    if (bendrodj > 4 * radius_j)
                        tempsize_j = size / 2;
                    // new approach  - close cylinder
                    for (unsigned int bri = 0; bri < size; bri += tempsize_i) {
                        for (unsigned int brj = 0; brj < size;
                             brj += tempsize_j) {
                            const Vector3 x1 = rodptrs[ri]->x[bri];
                            const Vector3 x2 = rodptrs[rj]->x[brj];
                            const Vector3 edge1 =
                                rodptrs[ri]->x[bri + tempsize_i] -
                                rodptrs[ri]->x[bri];
                            const Vector3 edge2 =
                                rodptrs[rj]->x[brj + tempsize_j] -
                                rodptrs[rj]->x[brj];
                            const vector<Vector3> minVectors =
                                findMinDistVectors(x1, edge1, x2, edge2);
                            const Vector3 dVector = minVectors[1];
                            if (dVector.length() <= 3 * (radius_i + radius_j)) {
                                rodptrs[ri]->nearbyRods.push_back(rj);
                                bri = size;
                                brj = size;
                            }
                        }  // rj elements
                    }      // ri elements
                }          // for j
            }              // for i
            rodptrs[0]->timer = 0;
        }
        // set cord no, contact points diagnostics to zero before starting to
        // check
        rodptrs[0]->cord = 0.0;
        rodptrs[0]->intercontacts.clear();
        rodptrs[0]->intercontacts_index.clear();
        rodptrs[0]->intercontacts_index_local.clear();
        rodptrs[0]->intercontacts_delta_n.clear();
        rodptrs[0]->intercontacts_delta_t.clear();
        rodptrs[0]->intercontacts_load.clear();
        rodptrs[0]->intercontacts_slip_vel.clear();
        for (unsigned int ri = 0; ri < (numRods - 1); ri++) {
            const REAL radius_i = rodptrs[ri]->r[0];

            /* clear up near contacts and tang_delta */
            rodptrs[ri]->nearcontacts_index_old =
                rodptrs[ri]->nearcontacts_index;
            rodptrs[ri]->tang_delta_old = rodptrs[ri]->tang_delta;
            rodptrs[ri]->nearcontacts_index.clear();
            rodptrs[ri]->tang_delta.clear();

            for (unsigned int j = 0; j < rodptrs[ri]->nearbyRods.size() &&
                                     rodptrs[ri]->nearbyRods.size() > 0;
                 j++) {
                unsigned int rj = rodptrs[ri]->nearbyRods[j];
                const REAL radius_j = rodptrs[rj]->r[0];
                // part added for bent rods
                const Vector3 rodi = rodptrs[ri]->x[size] - rodptrs[ri]->x[0];
                const Vector3 rodj = rodptrs[rj]->x[size] - rodptrs[rj]->x[0];
                const Vector3 straightrodi =
                    (rodptrs[ri]->x[1] - rodptrs[ri]->x[0]).unitize();
                const Vector3 straightrodj =
                    (rodptrs[rj]->x[1] - rodptrs[rj]->x[0]).unitize();
                const REAL bendrodi = (rodi * straightrodi).length();
                const REAL bendrodj = (rodj * straightrodj).length();
                int tempsize_i = size;
                int tempsize_j = size;
                if (bendrodi > 4 * radius_i)
                    tempsize_i = size / 2;
                if (bendrodj > 4 * radius_j)
                    tempsize_j = size / 2;
                for (unsigned int bri = 0; bri < size; bri += tempsize_i) {
                    for (unsigned int brj = 0; brj < size; brj += tempsize_j) {
                        const Vector3 x1 = rodptrs[ri]->x[bri];
                        const Vector3 x2 = rodptrs[rj]->x[brj];
                        const Vector3 edge1 = rodptrs[ri]->x[bri + tempsize_i] -
                                              rodptrs[ri]->x[bri];
                        const Vector3 edge2 = rodptrs[rj]->x[brj + tempsize_j] -
                                              rodptrs[rj]->x[brj];
                        const Vector3 x1_shiftnormal =
                            (rodptrs[ri]->x[bri + tempsize_i / 2] -
                             (rodptrs[ri]->x[bri + tempsize_i] +
                              rodptrs[ri]->x[bri]) *
                                 0.5) *
                            0.5;
                        const Vector3 x2_shiftnormal =
                            (rodptrs[rj]->x[brj + tempsize_j / 2] -
                             (rodptrs[rj]->x[brj + tempsize_j] +
                              rodptrs[rj]->x[brj]) *
                                 0.5) *
                            0.5;
                        const Vector3 x1_shift = x1_shiftnormal + x1;
                        const Vector3 x2_shift = x2_shiftnormal + x2;
                        const REAL sum_r1r2 = radius_i + radius_j;
                        const vector<Vector3> minVectors = findMinDistVectors(
                            x1_shift, edge1, x2_shift, edge2);
                        const Vector3 dVector = minVectors[1];  // pointing to J
                        const Vector3 dVectorDirection = dVector.unitize();
                        if ((dVector.length() - sum_r1r2) < 0) {
                            const bool yesCollision = true;
                            const Vector3 pVector =
                                minVectors[0];  // the position on i
                            // Calculating the contact point previous point:
                            unsigned int ContactIndexI =
                                bri + (pVector - x1_shift).length() /
                                          edge1.length() * tempsize_i;
                            unsigned int ContactIndexJ =
                                brj +
                                ((pVector + dVector) - x2_shift).length() /
                                    edge2.length() * tempsize_j;
                            ContactIndexI = (ContactIndexI == size)
                                                ? size - 1
                                                : ContactIndexI;
                            ContactIndexJ = (ContactIndexJ == size)
                                                ? size - 1
                                                : ContactIndexJ;
                            const REAL gamma = (sum_r1r2 - dVector.length());
                            // Compute contact force
                            /* const REAL cForceContact = zeta * gamma; */
                            const REAL cForceContact = zeta * gamma * pow(fabs(gamma), 0.5);
                            // Compute damping force
                            const Vector3 vInteri =
                                0.5 * (rodptrs[ri]->v[ContactIndexI] +
                                       rodptrs[ri]->v[ContactIndexI + 1]);
                            const Vector3 vInterj =
                                0.5 * (rodptrs[rj]->v[ContactIndexJ] +
                                       rodptrs[rj]->v[ContactIndexJ + 1]);
                            const Vector3 vInterpenetration = vInterj - vInteri;
                            const REAL vNorm =
                                vInterpenetration % dVectorDirection;
                            /* const REAL cForceDamping = -zeta_nu * vNorm; */
                            const REAL cForceDamping = -zeta_nu * vNorm * pow(fabs(gamma), 0.5);
                            const Vector3 cForce =
                                yesCollision * (cForceContact + cForceDamping) *
                                dVectorDirection;
                            // Remenber that the first and last points have half
                            // the mass! Hence the strange coefficient in front
                            // of cForce
                            rodptrs[ri]->externalForces[ContactIndexI] -=
                                ((ContactIndexI == 0) ? 0.5 : 1.0) * cForce;
                            rodptrs[ri]->externalForces[ContactIndexI + 1] -=
                                (((ContactIndexI + 1) == size) ? 0.5 : 1.0) *
                                cForce;
                            rodptrs[rj]->externalForces[ContactIndexJ] +=
                                ((ContactIndexJ == 0) ? 0.5 : 1.0) * cForce;
                            rodptrs[rj]->externalForces[ContactIndexJ + 1] +=
                                (((ContactIndexJ + 1) == size) ? 0.5 : 1.0) *
                                cForce;
                            const REAL cForcem = cForce.length();
                            /* push contact data */
                            rodptrs[0]->intercontacts.push_back(pVector +
                                                                dVector / 2);
                            rodptrs[0]->intercontacts_delta_n.push_back(gamma);
                            rodptrs[0]->intercontacts_load.push_back(cForcem);
                            rodptrs[0]->cord += 2;
                            rodptrs[0]->intercontacts_index.push_back({ri, rj});
                            rodptrs[0]->intercontacts_index_local.push_back(
                                {ContactIndexI, ContactIndexJ});
                            // Compute axial velocity
                            const Vector3 slipVel =
                                vInterpenetration - (vNorm)*dVectorDirection;
                            const REAL vel = slipVel.length();
                            rodptrs[0]->intercontacts_slip_vel.push_back(vel);
                            {
                                const REAL f = _linear(vel, vStatic);
                                const Vector3 kineticFrictionForce =
                                    -(1.0 - f) * isotropicFCK * cForcem *
                                    slipVel.unitize();

                                /* calculating spring friction */
                                bool is_old_contact = false;
                                unsigned int old_contact_serial;
                                Vector3 delta_t;
                                std::array<unsigned, 3> contact_id = {
                                    rj, ContactIndexI, ContactIndexJ};
                                /* find if it is an old contact */
                                for (auto k = 0;
                                     k < rodptrs[ri]
                                             ->nearcontacts_index_old.size() &&
                                     rodptrs[ri]
                                             ->nearcontacts_index_old.size() >
                                         0;
                                     k++) {
                                    std::array<unsigned, 3> old_contact_id =
                                        rodptrs[ri]->nearcontacts_index_old[k];
                                    if (contact_id == old_contact_id) {
                                        is_old_contact = true;
                                        old_contact_serial = k;
                                        k = rodptrs[ri]
                                                ->nearcontacts_index_old.size();
                                    }
                                }
                                const REAL dt = time - time_old;
                                time_old = time;
                                delta_t =
                                    (is_old_contact ==
                                         true &&
                                         rodptrs[ri]->tang_delta_old.size() > 0)
                                        ? (dt * slipVel +
                                           rodptrs[ri]->tang_delta_old
                                               [old_contact_serial])
                                        : Vector3(0.0, 0.0, 0.0);
                                rodptrs[ri]->nearcontacts_index.push_back(
                                    contact_id);

                                /* const Vector3 spring_force = */
                                /*     -zeta_tang * delta_t; */
                                /* const Vector3 damping_force = */
                                /*     -zeta_nu_tang * slipVel; */
                                /* Vector3 kineticFrictionForce; */
                                /* if ((spring_force + damping_force).length() > */
                                /*     (isotropicFCK * cForcem)) { */
                                /*     kineticFrictionForce = -isotropicFCK * */
                                /*                            cForcem * */
                                /*                            slipVel.unitize(); */
                                /*     delta_t = (-isotropicFCK * cForcem * */
                                /*                    slipVel.unitize() - */
                                /*                damping_force) / */
                                /*               zeta_tang; */
                                /* } else */
                                /*     kineticFrictionForce = */
                                /*         -(spring_force + damping_force) */
                                /*              .length() * */
                                /*         slipVel.unitize(); */

                                rodptrs[ri]->tang_delta.push_back(delta_t);
                                rodptrs[0]->intercontacts_delta_t.push_back(
                                    delta_t.length());
                                rodptrs[ri]
                                    ->kineticFrictionsForce[ContactIndexI] -=
                                    ((ContactIndexI == 0) ? 0.5 : 1.0) *
                                    kineticFrictionForce;
                                rodptrs[ri]
                                    ->kineticFrictionsForce[ContactIndexI +
                                                            1] -=
                                    ((ContactIndexI + 1 == size) ? 0.5 : 1.0) *
                                    kineticFrictionForce;
                                rodptrs[rj]
                                    ->kineticFrictionsForce[ContactIndexJ] +=
                                    ((ContactIndexJ == 0) ? 0.5 : 1.0) *
                                    kineticFrictionForce;
                                rodptrs[rj]
                                    ->kineticFrictionsForce[ContactIndexJ +
                                                            1] +=
                                    ((ContactIndexJ + 1 == size) ? 0.5 : 1.0) *
                                    kineticFrictionForce;
                            }  // Friction in axial direction
                        }      // collision
                    }          // brj
                }              // bri
            }                  // j
        }                      // i
        rodptrs[0]->cord /= numRods;
    }  // function
};
#endif
