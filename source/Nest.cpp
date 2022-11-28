#include "Nest.h"
Nest::Nest(const int argc, const char** argv) {}
void Nest::run() {
    /* const REAL cyl_radius = 0.5 * 241.3;     // Diameter 247 mm */
    /* const REAL cyl_radius = 0.5 * 190.5;     // Diameter 247 mm */
    const REAL cyl_radius = 0.5 * 139.7;  // Diameter 247 mm
    /* const REAL cyl_radius = 0.5 * 114.3;     // Diameter 247 mm */
    /* NOR = 1250;                      // 300, 400, 500 : number of sticks */
    /* NOR = 625;                      // 300, 400, 500 : number of sticks */
    NOR = 458;  // 300, 400, 500 : number of sticks
    /* NOR = 375;                      // 300, 400, 500 : number of sticks */
    const REAL real_dens = 8.94e-4;  // assuming mass of 0.5395 g: g/mm^3
    const REAL mult = 10;            // scaled up density
    const REAL g = 9810 / mult;      // mm/sec^2
    // Rod parameters             mm
    const int n = 8;  // number of elements in a rod
    /* const REAL L0 = 130.83;                  // 159.7,
     * 130.8, 76.1, 55.1, 29.9 mm */
    /* const REAL L0 = 103.29;                  // 159.7,
     * 130.8, 76.1, 55.1, 29.9 mm */
    const REAL L0 = 75.74;  // 159.7, 130.8, 76.1, 55.1, 29.9 mm
    /* const REAL L0 = 61.97;                  // 159.7, 130.8, 76.1, 55.1, 29.9
     * mm */
    const REAL r0 = 0.5 * 2.415;            // mm
    const REAL density = real_dens * mult;  // 8.233e-4//5e-3//10e-3//20e-3
    const REAL volume = M_PI * r0 * r0 * L0;
    const REAL totalMass = density * volume;  // grams
    /* const REAL E = 1e09;                      // 12e09 GPa */
    /* const REAL dt = 3e-6;  // 2.53e-7//6.2e-7//8.8e-7//12.5e-7 */
    const REAL E = 12e09;    // 12e09 GPa
    const REAL dt = 8.8e-7;  // 2.53e-7//6.2e-7//8.8e-7//12.5e-7
    const REAL timeSimulation = 5e4;
    // Dumping frequencies (number of frames/dumps per unit time)
    const REAL diagPerUnitTime = 50;
    const REAL povrayPerUnitTime = 20;
    const REAL dL0 = L0 / (double)n;  // length of cross-section element
    const REAL A0 = M_PI * r0 * r0;
    const REAL poissonRatio = 0.305;  // Incompressible
    const REAL G = E / (poissonRatio + 1.0);
    // Cylinder parameters
    const REAL R =
        1.0 *
        cyl_radius;  // controls radial extent of rod origin initialisation
    /* const REAL H = */
    /*     80.0;  // controls height range of rod origin initialisation */
    const REAL H = 30.0;  // controls height range of rod origin initialisation
    // Define rod ----Random rods
    vector<Vector3> directionRod;
    vector<Vector3> normalRod;
    vector<Vector3> originRod;
    srand(time(NULL));
    for (unsigned int i = 0; i < NOR; i++) {
        int Length = 100 * R;
        REAL a = (rand() % Length) * 0.01;  // distance from axes
        REAL b = (rand() % 360) * M_PI /
                 180;  // rotaroin angle of the starting point
        int r = 1000 * H;
        REAL c = (rand() % r + 1000 * r0) * 0.001;  // distance from ground
        Vector3 RodOrigin =
            Vector3(a * cos(b), a * sin(b), c);  // from 0.003 to c
        originRod.push_back(RodOrigin);
        REAL alpha = (rand() % 360) * M_PI / 180;
        REAL beta = (rand() % 179 + 1) * M_PI / 180;
        Vector3 RelativeEndpoint =
            L0 *
            Vector3(cos(beta) * cos(alpha), cos(beta) * sin(alpha), sin(beta));
        Vector3 Endpoint = RodOrigin + RelativeEndpoint;
        while ((Endpoint.x * Endpoint.x + Endpoint.y * Endpoint.y) >
               cyl_radius * cyl_radius) {
            if (beta < (M_PI) / 2) {
                beta += 0.05;
            } else {
                beta -= 0.05;
            }
            /* alpha = (rand() % 360) * M_PI / 180; */
            /* beta = (rand() % 179 + 1) * M_PI / 180; */
            RelativeEndpoint =
                Vector3(L0 * cos(beta) * cos(alpha),
                        L0 * cos(beta) * sin(alpha), L0 * sin(beta));
            Endpoint = RodOrigin + RelativeEndpoint;
        }
        Vector3 direction = Vector3(L0 * cos(beta) * cos(alpha),
                                    L0 * cos(beta) * sin(alpha), L0 * sin(beta))
                                .unitize();
        directionRod.push_back(direction);
        const Vector3 normal =
            Vector3(1, 1, (-direction.x - direction.y) / direction.z).unitize();
        normalRod.push_back(normal);
    }
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
    vector<RodBC*> boundaryConditionsPtrs;
    for (unsigned int i = 0; i < NOR; i++) {
        boundaryConditionsPtrs.push_back(&freeBC);
    }
    // Define plane
    const REAL angle = 0.0;
    const Vector3 originPlane = Vector3(0.0, 0.0, 0.0);
    const Vector3 normalPlane = Vector3(0.0, sin(angle), cos(angle)).unitize();
    // Set up substrate properties
    const REAL kPlane = 10e5;  // stiffness of ground
    const REAL nuPlane = 1e3;  // viscous damping of ground
    const REAL muStatic = 0.4;
    /* const REAL muStatic = 0.1; */
    const REAL muKinetic = muStatic * 0.5;
    const REAL vStatic = 1e-6;
    FrictionPlaneInteraction frictionPlane = FrictionPlaneInteraction(
        rodPtrs, normalPlane, originPlane, kPlane, nuPlane, muKinetic, muStatic,
        cyl_radius, dt, vStatic);
    vector<Interaction*> substrateInteractionsPtrs;
    substrateInteractionsPtrs.push_back(&frictionPlane);
    // Pack all forces together
    GravityForce gravity = GravityForce(Vector3(0.0, 0.0, -g));
    MultipleForces multipleForces;
    multipleForces.add(&gravity);
    MultipleForces* multipleForcesPtr = multipleForces.get();
    vector<ExternalForces*> externalForcesPtrs;
    for (unsigned int i = 0; i < NOR; i++) {
        externalForcesPtrs.push_back(multipleForcesPtr);
    }
    // Set up External Contact
    /* const REAL muisotropicK = 0.4; */
    /* const REAL muisotropicS = 0.4; */
    const REAL muisotropicK = 0.4;
    const REAL muisotropicS = 0.4;
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
}
