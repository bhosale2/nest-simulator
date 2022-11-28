#include "Mindlin_friction.h"
Mindlin_friction::Mindlin_friction(const int argc, const char** argv) {}
void Mindlin_friction::run() {
    NOR = 2;                      // 1250 : number of sticks
    const REAL real_dens = 8.4e-4;  // g/mm^3
    const REAL mult = 10;           // scaled up density
    const REAL g = 9810 / mult;     // mm/sec^2
    // Rod parameters             mm
    const int n = 8;                        // number of elements in a rod
    const REAL L0 = 130.8;                  // 159.7, 130.8, 76.1, 55.1, 29.9 mm
    const REAL r0 = 1.2;                    // mm
    const REAL density = real_dens * mult;  // 8.233e-4//5e-3//10e-3//20e-3
    const REAL volume = M_PI * r0 * r0 * L0;
    const REAL totalMass = density * volume;  // grams
    const REAL E = 1e09;                      // 12e09 GPa
    const REAL dt = 3e-6;  // 2.53e-7//6.2e-7//8.8e-7//12.5e-7
    const REAL timeSimulation = 2;
    // Dumping frequencies (number of frames/dumps per unit time)
    const REAL diagPerUnitTime = 50;
    const REAL povrayPerUnitTime = 20;
    const REAL dL0 = L0 / (double)n;  // length of cross-section element
    const REAL A0 = M_PI * r0 * r0;
    const REAL poissonRatio = 0.305;  // Incompressible
    const REAL G = E / (poissonRatio + 1.0);
    // Define rod ----Random rods
    vector<Vector3> directionRod;
    vector<Vector3> normalRod;
    vector<Vector3> originRod;
    originRod.push_back(Vector3(-L0/2, 0.0, 10.0 + 2 * r0));
    directionRod.push_back(Vector3(1.0, 0.0, 0.0));
    normalRod.push_back(Vector3(0.0, 0.0, 1.0));
    originRod.push_back(Vector3(0.0, -L0/2, 10.0));
    directionRod.push_back(Vector3(0.0, 1.0, 0.0));
    normalRod.push_back(Vector3(0.0, 0.0, 1.0));
    // Second moment of area for disk cross section
    const REAL I0_1 = A0 * A0 / (4.0 * M_PI);
    const REAL I0_2 = I0_1;
    const REAL I0_3 = 2.0 * I0_1;
    const Matrix3 I0 = Matrix3(I0_1, 0.0, 0.0, 0.0, I0_2, 0.0, 0.0, 0.0, I0_3);
    // Mass inertia matrix for disk cross section
    const Matrix3 J0 = density * dL0 * I0;
    // Bending matrix
    Matrix3 B0 =
        Matrix3(E * I0_1, 0.0, 0.0, 0.0, E * I0_2, 0.0, 0.0, 0.0, G * I0_3);
    // Shear matrix
    Matrix3 S0 = Matrix3((4.0 / 3.0) * G * A0, 0.0, 0.0, 0.0,
                         (4.0 / 3.0) * G * A0, 0.0, 0.0, 0.0, E * A0);
    // Initialize straight rod and pack it into a vector of pointers to rod
    const REAL initialTotalTwist = 0.0;
    const REAL nu = 0.2;
    const REAL relaxationNu = 0.0;
    const bool useSelfContact = false;
    vector<Rod*> rodPtrs;
    for (unsigned int i = 0; i < NOR; i++) {
        Rod* rod = RodInitialConfigurations::multstraightRod(
            n, totalMass, r0, J0, B0, S0, L0, initialTotalTwist, originRod[i],
            directionRod[i], normalRod[i], nu, relaxationNu, useSelfContact);
        rodPtrs.push_back(rod);
        rod->update(0.0);
        rod->computeEnergies();
    }
    // Pack boundary conditions
    FreeBC freeBC = FreeBC();
    RestrictfullBC  resBC  = RestrictfullBC(originRod[1], directionRod[1], L0);
    vector<RodBC*> boundaryConditionsPtrs;
    boundaryConditionsPtrs.push_back(&freeBC);
    boundaryConditionsPtrs.push_back(&resBC);
    // Define plane
    const REAL angle = 0.0;
    const Vector3 originPlane = Vector3(0.0, 0.0, 0.0);
    const Vector3 normalPlane = Vector3(0.0, sin(angle), cos(angle)).unitize();
    // Set up substrate properties
    const REAL kPlane = 8e4;   // stiffness of ground
    const REAL nuPlane = 500;  // viscous damping of ground
    const REAL muKinetic = 0.35 * 0.5;
    const REAL muStatic = 0.35;
    const REAL vStatic = 1e-6;
    const REAL cyl_radius  = 1e6;
    FrictionPlaneInteraction frictionPlane = FrictionPlaneInteraction(
        rodPtrs, normalPlane, originPlane, kPlane, nuPlane, muKinetic, muStatic,
        cyl_radius, dt, vStatic);
    vector<Interaction*> substrateInteractionsPtrs;
    substrateInteractionsPtrs.push_back(&frictionPlane);

    // Pack all forces together
    GravityForce gravity = GravityForce(Vector3(0.0, 0.0, -g));
    GravityForce pull_gravity = GravityForce(0.0 * Vector3(g, 0.0, 0.0));
    MultipleForces multipleForces, multipleForces2;
    multipleForces.add(&gravity);
    multipleForces.add(&pull_gravity);
    MultipleForces* multipleForcesPtr = multipleForces.get();
    MultipleForces* multipleForcesPtr2 = multipleForces2.get();
    vector<ExternalForces*> externalForcesPtrs;
    externalForcesPtrs.push_back(multipleForcesPtr);
    externalForcesPtrs.push_back(multipleForcesPtr2);

    // Set up External Contact
    const REAL muisotropicK = 0.45;
    const REAL muisotropicS = 0.45;
    vector<pair<int, int> > attachpoint;
    vector<ExternalContact*> externalcontactPtrs;
    ExternalContact externalcontact =
        ExternalContact(rodPtrs, muisotropicK, muisotropicS, attachpoint);
    externalcontactPtrs.push_back(&externalcontact);
    externalcontactPtrs[0]->setRodNearList();
    PolymerIntegrator* integrator = new PositionVerlet2nd(
        rodPtrs, externalForcesPtrs, boundaryConditionsPtrs,
        substrateInteractionsPtrs, externalcontactPtrs);
    // Instantiate simulator
    Polymer poly = Polymer(integrator);
    // Run simulation
    string outfileName = string("prova");
    const bool goodRun = poly.simulate(timeSimulation, dt, diagPerUnitTime,
                                       povrayPerUnitTime, outfileName);
    // Throw exception if something went wrong
    if (!goodRun)
        throw "not good run in localized helical buckling, what is going on?";
    exit(0);
}
