#include "Collision_test.h"
Coll_test::Coll_test(const int argc, const char** argv) {}
void Coll_test::run() {
    NOR = 2;                      // 1250 : number of sticks
    const REAL real_dens = 8.4e-4;  // g/mm^3
    const REAL mult = 10;           // scaled up density
    // Rod parameters             mm
    const int n = 8;                        // number of elements in a rod
    const REAL L0 = 100;                  // 159.7, 130.8, 76.1, 55.1, 29.9 mm
    const REAL r0 = 2;                    // mm
    const REAL density = real_dens * mult;  // 8.233e-4//5e-3//10e-3//20e-3
    const REAL volume = M_PI * r0 * r0 * L0;
    const REAL totalMass = density * volume;  // grams
    const REAL E = 1e09;                      // 12e09 GPa
    // Dumping frequencies (number of frames/dumps per unit time)
    const REAL dL0 = L0 / (double)n;  // length of cross-section element
    const REAL A0 = M_PI * r0 * r0;
    const REAL poissonRatio = 0.305;  // Incompressible
    const REAL G = E / (poissonRatio + 1.0);
    // Cylinder parameters
    // Define rod ----Random rods
    vector<Vector3> directionRod;
    vector<Vector3> normalRod;
    vector<Vector3> originRod;
    /* rod 1 */
    directionRod.push_back(Vector3(1.0, 0.0, 0.0));
    normalRod.push_back(Vector3(0.0, 0.0, 1.0));
    originRod.push_back(Vector3(0.0, 0.0, 0.0));
    /* rod 2 */
    directionRod.push_back(Vector3(0.0, 1.0, 0.0));
    normalRod.push_back(Vector3(0.0, 0.0, 1.0));
    originRod.push_back(Vector3(0.0, 0.0, 2.0));
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
    }
    // Set up External Contact
    const REAL muisotropicK = 0.45 * 0.5;
    const REAL muisotropicS = 0.45;
    vector<pair<int, int> > attachpoint;
    vector<ExternalContact*> externalcontactPtrs;
    ExternalContact externalcontact =
        ExternalContact(rodPtrs, muisotropicK, muisotropicS, attachpoint);
    externalcontactPtrs.push_back(&externalcontact);
    externalcontactPtrs[0]->setRodNearList();
    externalcontactPtrs[0]->RodRodCollision_Simple(0.0);
    cout << rodPtrs[0]->intercontacts;
}
