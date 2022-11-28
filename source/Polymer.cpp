/*
 * Polymer.cpp
 *
 *  Created on: Jun 20, 2014
 *      Author: mgazzola
 */
#include "Polymer.h"
void Polymer::computeEnergies() {
    bendingEnergy = 0.0;
    shearEnergy = 0.0;
    translationalEnergy = 0.0;
    rotationalEnergy = 0.0;
    totalInternalEnergy = 0;
    for (int i = 0; i < numRods; i++) {
        rodptrs[i]->computeEnergies();
        bendingEnergy += rodptrs[i]->bendingEnergy;
        shearEnergy += rodptrs[i]->shearEnergy;
        translationalEnergy += rodptrs[i]->translationalEnergy;
        rotationalEnergy += rodptrs[i]->rotationalEnergy;
        totalInternalEnergy += rodptrs[i]->totalInternalEnergy;
    }
}
bool Polymer::simulate(const REAL simulationTime,
                       const REAL dt0,
                       const unsigned int diagPerUnitTime,
                       const unsigned int povrayFramesPerUnitTime,
                       const string diagnostics,
                       const string integrationType,
                       const REAL CFL) {
#ifdef SNAKE_POV
    long unsigned int povrayStep = 0;
#endif
    /* unsigned int N = 1; */
    long unsigned int step = 0;
    double time = 0.0;
    double diagFotoTimer = 0.0;
    double povrayFotoTimer = 0.0;
    while (time <= simulationTime) {
        // Dump diagnostics

#ifdef SNAKE_POV
        {
            const double fotoDT =
                1.0 / ((diagPerUnitTime > 0) ? diagPerUnitTime : 1e-6);
            if ((diagPerUnitTime != 0) &&
                (diagFotoTimer > fotoDT || step == 0)) {
                rodptrs[0]->timer = 1;  // update list of nearby rods
                diagFotoTimer = 0.0;
            }
        }
#endif

// Dump Povray files
#ifdef SNAKE_POV
        const double fotoDT =
            1.0 /
            ((povrayFramesPerUnitTime > 0) ? povrayFramesPerUnitTime : 1e-6);
        if ((povrayFramesPerUnitTime != 0) &&
            (povrayFotoTimer > fotoDT || step == 0)) {
            povrayFotoTimer = 0.0;
            /* char bufPov[1000]; */
            /* sprintf(bufPov, "snake_%07d.pov", int(povrayStep)); */
            /* char bufData[1000]; */
            /* sprintf(bufData, "data_%07d.inc", int(povrayStep)); */
            // Dump Povray file
            /* fstream f; */
            /* f.open(bufPov, fstream::out); */
            /* f << "#include \"scenepovray.inc\"" */
            /*   << "\n"; */
            /* f << "#include \"" << bufData << "\"" */
            /*   << "\n"; */
            /* f << "\n"; */
            /* f.close(); */
            // Dump data to be rendered in the pov file
            /* for (unsigned int i = 0; i < numRods; i++) { */
            /*     rodptrs[i]->dumpPovray(bufPov, bufData, i, time); */
            /* } */
            /* for (unsigned int i = 0; i < numPistons; i++) { */
            /*     pistonptrs[i]->dumpPovray(bufPov, bufData, i, time); */
            /* } */
            /* dump energies */
            printEnergies(step, povrayStep);
            printRod_X(step, povrayStep);
            print_Nest_Stuff(step, time, povrayStep);
            povrayStep++;
        }
#endif
#ifdef SNAKE_VIZ
        if (step % 5000 == 0) {
            _paint(snake);
            cout << "time = " << time << endl;
        }
#endif
        // Integrate
        REAL dt = 0.0;
        dt = pint->integrate(time, dt0, step);
        // Update cumulative quantities
        time += dt;
        diagFotoTimer += dt;
        povrayFotoTimer += dt;
        step += 1;
        /* #ifndef NDEBUG */
        if (nanCheck()) {
            cout << "Found a NaN. Ending simulation." << endl;
            return false;
        }
        /* #endif */
    }
    // Compute energies
    /* computeEnergies(); */
    return true;
}

bool Polymer::nanCheck() {
    bool foundNan = false;
    for (int i = 0; i < numRods; i++)
        if (rodptrs[i]->nanCheck())
            foundNan = true;
    return foundNan;
}

/* prints different energies with time for all rods in same file */
void Polymer::printEnergies(const int step, const int povrayStep) {
    char buffer[1000];
    sprintf(buffer, "rodEnergies_%05d.csv", povrayStep);
    FILE* outfile = fopen(buffer, (step == 0) ? "w" : "a");
    for (unsigned int i = 0; i < numRods; i++) {
        rodptrs[i]->computeEnergies();
        fprintf(outfile, "%.8f,%.8f,%.8f,%.8f,%.8f\n",
                rodptrs[i]->totalInternalEnergy,
                rodptrs[i]->translationalEnergy, rodptrs[i]->rotationalEnergy,
                rodptrs[i]->bendingEnergy, rodptrs[i]->shearEnergy);
    }
    fclose(outfile);
}

/* prints all rod positions in the same file for a time */
/* need to dump N for each rod later */
void Polymer::printRod_X(const int step, const int povrayStep) {
    char buffer[1000];
    sprintf(buffer, "rodX_%05d.csv", povrayStep);
    FILE* outfile = fopen(buffer, (step == 0) ? "w" : "a");
    for (auto i = 0; i < numRods; i++) {
        for (auto j = 0; j < rodptrs[i]->x.size(); j++)
            fprintf(outfile, "%f,%f,%f\n", rodptrs[i]->x[j].x,
                    rodptrs[i]->x[j].y, rodptrs[i]->x[j].z);
    }
    fclose(outfile);
}

/* prints nest related stuff like nest height, coordination number, etc */
void Polymer::print_Nest_Stuff(const int step,
                               const REAL time,
                               const int povrayStep) {
    REAL MinH = 300;
    REAL MaxH = 0;
    for (auto i = 0; i < numRods; i++) {
        MaxH = std::max(
            MaxH, std::max(rodptrs[i]->x.back().z, rodptrs[i]->x.front().z));
        MinH = std::min(
            MinH, std::min(rodptrs[i]->x.back().z, rodptrs[i]->x.front().z));
    }
    char buffer[1000];
    sprintf(buffer, "nest.csv");
    FILE* outfile = fopen(buffer, (step == 0) ? "w" : "a");
    fprintf(outfile, "%f,%f,%f,%f,%f,%f,%f\n", time, rodptrs[0]->TopHeight,
            MinH, MaxH, rodptrs[0]->base_stress, rodptrs[0]->rad_stress,
            rodptrs[0]->cord);
    fclose(outfile);

    /* dump contacts, load and their indices */
    sprintf(buffer, "contacts_%05d.csv", povrayStep);
    FILE* confile = fopen(buffer, (step == 0) ? "w" : "a");
    for (auto i = 0; i < rodptrs[0]->intercontacts.size(); i++) {
        fprintf(confile, "%f,%f,%f,%f,%u,%u,%u,%u,%f,%f,%f\n",
                rodptrs[0]->intercontacts[i].x, rodptrs[0]->intercontacts[i].y,
                rodptrs[0]->intercontacts[i].z,
                rodptrs[0]->intercontacts_load[i],
                rodptrs[0]->intercontacts_index[i][0],
                rodptrs[0]->intercontacts_index[i][1],
                rodptrs[0]->intercontacts_index_local[i][0],
                rodptrs[0]->intercontacts_index_local[i][1],
                rodptrs[0]->intercontacts_slip_vel[i],
                rodptrs[0]->intercontacts_delta_n[i],
                rodptrs[0]->intercontacts_delta_t[i]);
    }
    fclose(confile);
    /* dump average orientation of rods */
    sprintf(buffer, "orient_%05d.csv", povrayStep);
    FILE* ornfile = fopen(buffer, (step == 0) ? "w" : "a");
    Vector3 orient;
    for (auto i = 0; i < numRods; i++) {
        orient = (rodptrs[i]->x.back() - rodptrs[i]->x.front()).unitize();
        fprintf(ornfile, "%f,%f,%f\n", orient.x, orient.y, orient.z);
    }
    fclose(ornfile);
}

void Polymer::printX(const int step,
                     const REAL time,
                     const string outfilename) {
    for (unsigned int i = 0; i < rodptrs.size(); i++) {
        char buffer[1000];
        sprintf(buffer, "%s_rodX_%05d", outfilename.c_str(), i);
        FILE* outfile = fopen(buffer, (step == 0) ? "w" : "a");
        fprintf(outfile, "%1.10e ", time);
        for (unsigned int j = 0; j < rodptrs[i]->x.size(); j++)
            fprintf(outfile, "%1.10e ", rodptrs[i]->x[j].x);
        fprintf(outfile, "\n");
        fclose(outfile);
    }
    for (unsigned int i = 0; i < rodptrs.size(); i++) {
        const Matrix3 Q0 = rodptrs[i]->Q.front();
        const Vector3 axisOfRotation = Vector3(0, 0, 1);
        char buffer[1000];
        sprintf(buffer, "%s_rodX_angle_%05d", outfilename.c_str(), i);
        FILE* outfile = fopen(buffer, (step == 0) ? "w" : "a");
        fprintf(outfile, "%1.10e ", time);
        for (unsigned int j = 0; j < rodptrs[i]->Q.size(); j++) {
            const Vector3 angvector = (Q0 * rodptrs[i]->Q[j].T()).log();
            const REAL length = angvector.length();
            const int sign = floor((angvector % (length * axisOfRotation)) /
                                   (length * length));
            const REAL angle = -sign * length;
            fprintf(outfile, "%1.10e ", angle);
        }
        fprintf(outfile, "\n");
        fclose(outfile);
    }
    for (unsigned int i = 0; i < rodptrs.size(); i++) {
        const Matrix3 Q0 = rodptrs[i]->Q.front();
        const Vector3 axisOfRotation = Vector3(1, 0, 0);
        char buffer[1000];
        sprintf(buffer, "%s_rodX_anglebend_%05d", outfilename.c_str(), i);
        FILE* outfile = fopen(buffer, (step == 0) ? "w" : "a");
        fprintf(outfile, "%1.10e ", time);
        for (unsigned int j = 0; j < rodptrs[i]->Q.size(); j++) {
            const Vector3 angvector = (Q0 * rodptrs[i]->Q[j].T()).log();
            const REAL length = angvector.length();
            const int sign = floor((angvector % (length * axisOfRotation)) /
                                   (length * length));
            const REAL angle = -sign * length;
            fprintf(outfile, "%1.10e ", angle);
        }
        fprintf(outfile, "\n");
        fclose(outfile);
    }
    for (unsigned int i = 0; i < rodptrs.size(); i++) {
        char buffer[1000];
        sprintf(buffer, "%s_rodX_sigma_%05d", outfilename.c_str(), i);
        FILE* outfile = fopen(buffer, (step == 0) ? "w" : "a");
        fprintf(outfile, "%1.10e ", time);
        for (unsigned int j = 0; j < rodptrs[i]->Q.size(); j++)
            fprintf(outfile, "%1.10e ", rodptrs[i]->shearStrain0[j].length());
        fprintf(outfile, "\n");
        fclose(outfile);
    }
}
void Polymer::printXV(const int step,
                      const REAL time,
                      const string outfilename) {
    for (unsigned int i = 0; i < rodptrs.size(); i++) {
        char buffer[1000];
        sprintf(buffer, "%s_rodXV_%05d", outfilename.c_str(), i);
        FILE* outfile = fopen(buffer, (step == 0) ? "w" : "a");
        fprintf(outfile, "%1.10e ", time);
        for (unsigned int j = 0; j < rodptrs[i]->x.size(); j++)
            fprintf(outfile, "%1.10e ", rodptrs[i]->x[j].x);
        fprintf(outfile, "\n");
        fprintf(outfile, "%1.10e ", time);
        for (unsigned int j = 0; j < rodptrs[i]->x.size(); j++)
            fprintf(outfile, "%1.10e ", rodptrs[i]->x[j].y);
        fprintf(outfile, "\n");
        fprintf(outfile, "%1.10e ", time);
        for (unsigned int j = 0; j < rodptrs[i]->x.size(); j++)
            fprintf(outfile, "%1.10e ", rodptrs[i]->x[j].z);
        fprintf(outfile, "\n");
        fprintf(outfile, "%1.10e ", time);
        for (unsigned int j = 0; j < rodptrs[i]->x.size(); j++)
            fprintf(outfile, "%1.10e ", rodptrs[i]->v[j].x);
        fprintf(outfile, "\n");
        fprintf(outfile, "%1.10e ", time);
        for (unsigned int j = 0; j < rodptrs[i]->x.size(); j++)
            fprintf(outfile, "%1.10e ", rodptrs[i]->v[j].y);
        fprintf(outfile, "\n");
        fprintf(outfile, "%1.10e ", time);
        for (unsigned int j = 0; j < rodptrs[i]->x.size(); j++)
            fprintf(outfile, "%1.10e ", rodptrs[i]->v[j].z);
        fprintf(outfile, "\n");
        fclose(outfile);
    }
}
void Polymer::print_s_internalTorques(const string outfilename) {
    for (unsigned int i = 0; i < rodptrs.size(); i++) {
        char buffer[1000];
        sprintf(buffer, "%s_%05d", outfilename.c_str(), i);
        FILE* outfile = fopen(buffer, "w");
        const vector<REAL> s = vCumSum(rodptrs[i]->l);
        assert(s.size() == rodptrs[i]->bendingInternalTorques0.size() + 1);
        for (unsigned int j = 0; j < rodptrs[i]->bendingInternalTorques0.size();
             j++)
            fprintf(outfile, "%1.10e\t%1.10e\t%1.10e\t%1.10e\n", s[j],
                    rodptrs[i]->bendingInternalTorques0[j].x,
                    rodptrs[i]->bendingInternalTorques0[j].y,
                    rodptrs[i]->bendingInternalTorques0[j].z);
        fprintf(outfile, "\n");
        fclose(outfile);
    }
}
void Polymer::print_s_coordinates(const string outfilename) {
    for (unsigned int i = 0; i < rodptrs.size(); i++) {
        char buffer[1000];
        sprintf(buffer, "%s_%05d", outfilename.c_str(), i);
        FILE* outfile = fopen(buffer, "w");
        const vector<REAL> s = vCumSum(rodptrs[i]->l);
        assert(s.size() + 1 == rodptrs[i]->x.size());
        assert(s.size() > 1);
        const REAL zero = 0.0;
        fprintf(outfile, "%1.10e\t%1.10e\t%1.10e\t%1.10e\n", zero,
                rodptrs[i]->x[0].x, rodptrs[i]->x[0].y, rodptrs[i]->x[0].z);
        for (unsigned int j = 1; j < rodptrs[i]->x.size(); j++) {
            assert((j - 1) < s.size());
            fprintf(outfile, "%1.10e\t%1.10e\t%1.10e\t%1.10e\n", s[j - 1],
                    rodptrs[i]->x[j].x, rodptrs[i]->x[j].y, rodptrs[i]->x[j].z);
        }
        fprintf(outfile, "\n");
        fclose(outfile);
    }
}
void Polymer::print_s_internalShears(const string outfilename) {
    for (unsigned int i = 0; i < rodptrs.size(); i++) {
        char buffer[1000];
        sprintf(buffer, "%s_%05d", outfilename.c_str(), i);
        FILE* outfile = fopen(buffer, "w");
        vector<REAL> s = vector<REAL>(rodptrs[i]->shearInternalForces0.size());
        assert(s.size() == rodptrs[i]->shearInternalForces0.size());
        s[0] = rodptrs[i]->l[0] / 2.0;
        for (unsigned int j = 1; j < s.size(); j++)
            s[j] = s[j - 1] +
                   (rodptrs[i]->l[j - 1] / 2.0 + rodptrs[i]->l[j] / 2.0);
        for (unsigned int j = 0; j < rodptrs[i]->shearInternalForces0.size();
             j++) {
            const Vector3 shearInLabFrame =
                rodptrs[i]->Q[j].T() * rodptrs[i]->shearInternalForces0[j];
            fprintf(outfile, "%1.10e\t%1.10e\t%1.10e\t%1.10e\n", s[j],
                    shearInLabFrame.x, shearInLabFrame.y, shearInLabFrame.z);
        }
        fprintf(outfile, "\n");
        fclose(outfile);
    }
}
void Polymer::print_s_curvatures(const string outfilename) {
    for (unsigned int i = 0; i < rodptrs.size(); i++) {
        char buffer[1000];
        sprintf(buffer, "%s_%05d", outfilename.c_str(), i);
        FILE* outfile = fopen(buffer, "w");
        const vector<REAL> s = vCumSum(rodptrs[i]->l);
        assert(s.size() == rodptrs[i]->k0.size() + 1);
        for (unsigned int j = 0; j < rodptrs[i]->k0.size(); j++) {
            const REAL ed = rodptrs[i]->ed[j];
            fprintf(outfile, "%1.10e\t%1.10e\t%1.10e\t%1.10e\n", s[j],
                    rodptrs[i]->k0[j].x / ed, rodptrs[i]->k0[j].y / ed,
                    rodptrs[i]->k0[j].z / ed);
        }
        fprintf(outfile, "\n");
        fclose(outfile);
    }
}
#ifdef SNAKE_VIZ
void Polymer::_paint(Rod* snake) {
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
    glPushAttrib(GL_ENABLE_BIT);
    snake->paint();
    glPopAttrib();
    glutSwapBuffers();
}
#endif
