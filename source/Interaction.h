#ifndef INTERACTION_H_
#define INTERACTION_H_
#include <assert.h>
#include <limits>
#include "MathFunctions.h"
#include "Rod.h"
// Class defining interaction forces between rod and rod, substrate, etc
class Interaction {
   protected:
    typedef std::vector<Vector3> VV3;
    typedef std::vector<Matrix3> VM3;
    typedef std::vector<REAL> VREAL;
    typedef std::vector<bool> VBOOL;
    typedef std::vector<int> VINT;
    typedef std::vector<Rod*> Vrodptr;
    typedef std::vector<Interaction*> Vinterptr;
    Vrodptr rodptrs;

   public:
    Interaction(Vrodptr rodptrs) : rodptrs(rodptrs) {}
    virtual ~Interaction() {}
    virtual void applyForces(const REAL time) = 0;
};
class FrictionPlaneInteraction : public Interaction {
   protected:
    const REAL
        kPlane;  // Elastic response to prevent interpenetration rod-plane
    const REAL etaPlane;  // Dissipative term appliet to the elastice response
                          // rod-plane
    const REAL muKineticPlane;
    const REAL muStaticPlane;
    const REAL vStaticPlane;
    const Vector3 normalPlane;
    const Vector3 originPlane;
    const REAL dt;                  // time step - (s)
    const REAL stroke_start = 1.0;  // stroke beginning
    const REAL cyl_radius;          // radius of the cylinder - mm
    const REAL perturb_start =
        0.5;  // start time of cylindrical radius perturbation
    const REAL perturb_end =
        stroke_start;  // end time of cylindrical radius perturbation
    const REAL perturb_amp = 0.0;   // amplitude of perturbation
    REAL cyl_radius_perturb = 0.0;  // perturbation in the cylindrical radius
    /* const REAL top_start = 218.0;   // position of top wall at the start - mm
     */
    /* const REAL top_start = 136;   // position of top wall at the start - mm
     */
    const REAL top_start = 93;  // position of top wall at the start - mm
    /* const REAL top_start = 89;   // position of top wall at the start - mm */
    /* const REAL stroke_start = 100.0;     // stroke beginning */
    int bdirection = -1;        // direction of stroke: -1:downward, +1:upward
    REAL top_dynH = top_start;  // top dynamic height - mm
    REAL base_dynH = 0;         // base dynamic height - mm
    const REAL cycle_strain_rate =
        5;  // 1.0 rate of movement of top wall - mm/sec
    const REAL init_strain_rate =
        5;                   // faster rate before compression begins- mm/sec
    REAL strain_rate;        // dynamic strain rate
    REAL base_stress = 0.0;  // stress at the base - Pa
    /* const REAL stress_low = 143.66 + 10;    // lower limit of stress - Pa */
    /* const REAL stress_high = 143.66 + 800;  // higher limit of stress - Pa */
    /* const REAL stress_low = 91 + 20;    // lower limit of stress - Pa */
    /* const REAL stress_high = 91 + 1600;  // higher limit of stress - Pa */
    const REAL stress_low = 91 + 30;  // lower limit of stress - Pa
    /* const REAL stress_high = 91 + 2500;  // higher limit of stress - Pa */
    /* const REAL stress_high = 91 + 4000;  // higher limit of stress - Pa */
    /* const REAL stress_high = 91 + 7500;  // higher limit of stress - Pa */
    const REAL stress_high = 91 + 1e4;  // higher limit of stress - Pa
    /* const REAL stress_low = 91 + 50;    // lower limit of stress - Pa */
    /* const REAL stress_high = 91 + 4400;  // higher limit of stress - Pa */
    bool lid_contact = false;  // checks if the lid touches the rods
    const REAL max_cycle_height = 78;

    const REAL _linear(const REAL _vel, const REAL _velThresold) const {
        const REAL vel = fabs(_vel);
        const REAL velThresold = fabs(_velThresold);
        const REAL width = 0.5 * velThresold;
        const REAL velDiff = vel - velThresold;
        REAL f = 1.0;
        if (vel > (velThresold))
            // f = 1.0-velDiff/width;
            f = fabs(1.0 - min(1.0, velDiff / width));
        if (vel > (velThresold + width))
            f = 0.0;
        return f;
    }
    const REAL _sigmoid(const REAL _vel, const REAL _velThresold) const {
        const REAL vel = fabs(_vel);
        const REAL velThresold = fabs(_velThresold);
        const REAL velDiff = vel - 1.5 * velThresold;
        return 0.5 + 0.5 * erf(-(8.0 / velThresold) * velDiff);
    }

   public:
    FrictionPlaneInteraction(Vrodptr& rodptrs,
                             const Vector3 _normalPlane,
                             const Vector3 _originPlane,
                             const REAL _kPlane,
                             const REAL _etaPlane,
                             const REAL _muKineticPlane,
                             const REAL _muStaticPlane,
                             const REAL _cyl_radius,
                             const REAL _dt,
                             const REAL _vStatic = 1e-4)
        : Interaction(rodptrs),
          normalPlane(_normalPlane),
          originPlane(_originPlane),
          kPlane(_kPlane),
          etaPlane(_etaPlane),
          muStaticPlane(_muStaticPlane),
          muKineticPlane(_muKineticPlane),
          vStaticPlane(_vStatic),
          cyl_radius(_cyl_radius),
          dt(_dt) {}
    void applyForces(const REAL time) {
        const int nor = rodptrs.size();
        /* setting top wall height to save time */
        REAL MaxH = 0.0;
        for (auto j = 0; j < nor; j++)
            MaxH = std::max(MaxH, std::max(rodptrs[j]->x.back().z,
                                           rodptrs[j]->x.front().z));
        if (base_stress > stress_high) {
            bdirection = 1;
            strain_rate = cycle_strain_rate;
        }
        /* else if (base_stress < stress_low){ */
        else if (base_stress < stress_low && top_dynH > max_cycle_height) {
            /* else if (base_stress < stress_low && lid_contact == false){ */
            bdirection = -1;
            strain_rate = cycle_strain_rate;
        } else
            strain_rate = cycle_strain_rate;
        /* if (time > stroke_start && top_dynH < MaxH + 2 * rodptrs[0]->r[0]) */
        if (time > stroke_start)
            top_dynH += bdirection * strain_rate * dt;
        /* else */
        /*     top_dynH = MaxH + rodptrs[0]->r[0]; */

        rodptrs[0]->TopHeight = top_dynH;
        base_stress = 0.0;
        REAL rad_stress = 0.0;
        if (time < perturb_start || time > perturb_end)
            cyl_radius_perturb = 0.0;
        else
            cyl_radius_perturb = -(time - perturb_start) /
                                 (perturb_end - perturb_start) * perturb_amp *
                                 cyl_radius;

        lid_contact = false;
        for (auto j = 0; j < nor; j++) {
            Rod* Rod = rodptrs[j];
            for (auto i = 0; i < Rod->x.size();
                 /* i += Rod->x.size() - 1)  // direct jump */
                 i++) {
                // Assumed base doesnt move and is at z = 0
                if (Rod->x[i].z < Rod->r[0] + base_dynH) {
                    const Vector3 elementV = Rod->v[i];
                    const Vector3 normalPlane = Vector3(0.0, 0.0, 1.0);
                    const Vector3 overallRodForces =
                        Rod->totalInternalForces[i] + Rod->externalForces[i];

                    const REAL currentRodForcesInNormaldirectionSign =
                        (overallRodForces % normalPlane);
                    const Vector3 currentRodForcesInNormaldirection =
                        currentRodForcesInNormaldirectionSign * normalPlane;
                    const Vector3 currentRodForcesInplanardirection =
                        overallRodForces - currentRodForcesInNormaldirection;

                    /* const Vector3 elasticPlaneResponse = */
                    /*     -kPlane * (Rod->x[i].z - Rod->r[0]) * normalPlane; */
                    /* const Vector3 dampingForcePlane = */
                    /*     -etaPlane * (elementV % normalPlane) * normalPlane;
                     */
                    const REAL gamma =
                        fabs(Rod->x[i].z - Rod->r[0] - base_dynH);
                    const Vector3 elasticPlaneResponse =
                        -kPlane * (Rod->x[i].z - Rod->r[0] - base_dynH) *
                        pow(gamma, 0.5) * normalPlane;
                    const Vector3 dampingForcePlane =
                        -etaPlane * (elementV % normalPlane) * pow(gamma, 0.5) *
                        normalPlane;
                    Rod->externalForces[i] +=
                        elasticPlaneResponse + dampingForcePlane;
                    const REAL Forcem =
                        (elasticPlaneResponse + dampingForcePlane).length();
                    base_stress += Forcem / M_PI / cyl_radius / cyl_radius;
                    // Compute axial velocity
                    const Vector3 slipVel =
                        elementV - (elementV % normalPlane) * normalPlane;
                    // Friction in axial direction
                    {
                        const REAL vel = fabs(slipVel.length());
                        const REAL f = _linear(vel, vStaticPlane);
                        const Vector3 kineticFrictionForce =
                            -(1.0 - f) * muKineticPlane * Forcem *
                            slipVel.unitize();
                        const Vector3 staticFrictionForce =
                            f * muStaticPlane * Forcem *
                            currentRodForcesInplanardirection.unitize();
                        Rod->kineticFrictionsForce[i] += kineticFrictionForce;
                        Rod->staticFrictionsAxialForceForward[i] +=
                            staticFrictionForce;
                        Rod->staticFrictionsAxialForceBackward[i] +=
                            staticFrictionForce;
                    }  // Friction in axial direction
                }      // if within range

                // For wall=======  unchanged
                const REAL distanceFromAxes =
                    sqrt(Rod->x[i].x * Rod->x[i].x + Rod->x[i].y * Rod->x[i].y);
                if (distanceFromAxes >
                    cyl_radius + cyl_radius_perturb - Rod->r[0]) {
                    const Vector3 elementV = Rod->v[i];
                    const Vector3 normalPlane =
                        -Vector3(Rod->x[i].x, Rod->x[i].y, 0.0).unitize();
                    /* const Vector3 elasticPlaneResponse = */
                    /*     -kPlane * */
                    /*     (cyl_radius + cyl_radius_perturb - Rod->r[0] - */
                    /*      distanceFromAxes) * */
                    /*     normalPlane; */
                    /* const Vector3 dampingForcePlane = */
                    /*     -etaPlane * (elementV % normalPlane) * normalPlane;
                     */
                    const REAL gamma = fabs(cyl_radius + cyl_radius_perturb -
                                            Rod->r[0] - distanceFromAxes);
                    const Vector3 elasticPlaneResponse =
                        -kPlane *
                        (cyl_radius + cyl_radius_perturb - Rod->r[0] -
                         distanceFromAxes) *
                        pow(gamma, 0.5) * normalPlane;
                    const Vector3 dampingForcePlane =
                        -etaPlane * (elementV % normalPlane) * pow(gamma, 0.5) *
                        normalPlane;
                    Rod->externalForces[i] +=
                        elasticPlaneResponse + dampingForcePlane;
                    const REAL Forcem =
                        (elasticPlaneResponse + dampingForcePlane).length();
                    /* rad_stress += Forcem / M_PI / cyl_radius / top_dynH; */
                    const Vector3 overallRodForces =
                        Rod->totalInternalForces[i] + Rod->externalForces[i];

                    const REAL currentRodForcesInNormaldirectionSign =
                        (overallRodForces % normalPlane);
                    const Vector3 currentRodForcesInNormaldirection =
                        currentRodForcesInNormaldirectionSign * normalPlane;
                    const Vector3 currentRodForcesInplanardirection =
                        overallRodForces - currentRodForcesInNormaldirection;
                    // Compute axial velocity
                    const Vector3 slipVel =
                        elementV - (elementV % normalPlane) * normalPlane;
                    // Friction in axial direction
                    {
                        const REAL vel = fabs(slipVel.length());
                        const REAL f = _linear(vel, vStaticPlane);
                        const Vector3 kineticFrictionForce =
                            -(1.0 - f) * muKineticPlane * Forcem *
                            slipVel.unitize();
                        Rod->kineticFrictionsForce[i] += kineticFrictionForce;
                        const Vector3 staticFrictionForce =
                            f * muStaticPlane * Forcem *
                            currentRodForcesInplanardirection.unitize();
                        Rod->staticFrictionsAxialForceForward[i] +=
                            staticFrictionForce;
                        Rod->staticFrictionsAxialForceBackward[i] +=
                            staticFrictionForce;
                    }  // Friction in axial direction
                }
                //    =========For top cover =======
                if (Rod->x[i].z > top_dynH - Rod->r[0]) {
                    lid_contact = true;
                    const Vector3 normalPlane = Vector3(0.0, 0.0, -1.0);
                    const Vector3 elementV =
                        (Rod->v[i] -
                         bdirection * strain_rate * Vector3(0.0, 0.0, 1.0));
                    /* const Vector3 elasticPlaneResponse = */
                    /*     -kPlane * (top_dynH - Rod->r[0] - Rod->x[i].z) * */
                    /*     normalPlane; */
                    /* const Vector3 dampingForcePlane = */
                    /*     -etaPlane * (elementV % normalPlane) * normalPlane;
                     */
                    const REAL gamma = fabs(top_dynH - Rod->r[0] - Rod->x[i].z);
                    const Vector3 elasticPlaneResponse =
                        -kPlane * (top_dynH - Rod->r[0] - Rod->x[i].z) *
                        pow(gamma, 0.5) * normalPlane;
                    const Vector3 dampingForcePlane =
                        -etaPlane * (elementV % normalPlane) * pow(gamma, 0.5) *
                        normalPlane;
                    Rod->externalForces[i] +=
                        elasticPlaneResponse + dampingForcePlane;
                    const REAL Forcem =
                        (elasticPlaneResponse + dampingForcePlane).length();
                    rad_stress += Forcem / M_PI / cyl_radius / cyl_radius;
                    const Vector3 overallRodForces =
                        Rod->totalInternalForces[i] + Rod->externalForces[i];

                    const REAL currentRodForcesInNormaldirectionSign =
                        (overallRodForces % normalPlane);
                    const Vector3 currentRodForcesInNormaldirection =
                        currentRodForcesInNormaldirectionSign * normalPlane;
                    const Vector3 currentRodForcesInplanardirection =
                        overallRodForces - currentRodForcesInNormaldirection;
                    // Compute axial velocity
                    const Vector3 slipVel =
                        elementV - (elementV % normalPlane) * normalPlane;
                    // Friction in axial direction
                    {
                        const REAL vel = fabs(slipVel.length());
                        const REAL f = _linear(vel, vStaticPlane);
                        const Vector3 kineticFrictionForce =
                            -(1.0 - f) * muKineticPlane * Forcem *
                            slipVel.unitize();
                        Rod->kineticFrictionsForce[i] += kineticFrictionForce;
                        const Vector3 staticFrictionForce =
                            f * muStaticPlane * Forcem *
                            currentRodForcesInplanardirection.unitize();
                        Rod->staticFrictionsAxialForceForward[i] +=
                            staticFrictionForce;
                        Rod->staticFrictionsAxialForceBackward[i] +=
                            staticFrictionForce;
                    }  // Friction in axial direction
                }      // if within range
            }          // for loop all elements
        }              // for loop all rods
        rodptrs[0]->base_stress = base_stress;
        rodptrs[0]->rad_stress = rad_stress;
    }  // apply force
};
#endif
