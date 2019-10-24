/* -------------------------------------------------------------------------- *
 *                                   OpenMM                                   *
 * -------------------------------------------------------------------------- *
 * This is part of the OpenMM molecular simulation toolkit originating from   *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "ExampleForce.h"
#include "MDIServer.h"
#include "internal/ExampleForceImpl.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/AssertionUtilities.h"
#include "mdi.h"
//////////
#include "openmm/internal/ContextImpl.h"
#include "openmm/NonbondedForce.h"
//////////

using namespace ExamplePlugin;
using namespace OpenMM;
using namespace std;

MDIServer::MDIServer(OpenMM::System& system) : system(system) {
}

void MDIServer::init(string mdi_options) {
    printf("BBBBBBBBBBBBBBBBBBBBBBB\n");
}

void MDIServer::listen(ContextImpl& context, string node) {
    printf("CCCCCCCCCCCCCCCCCCCCCCC\n");
    const OpenMM::System& system = context.getSystem();

    // <COORDS
    vector<Vec3> positions;
    context.getPositions(positions);
    printf("      pos: %f %f %f\n",positions[0][0],positions[0][1],positions[0][2]);

    // <VELOCITIES
    vector<Vec3> velocities;
    context.getVelocities(velocities);
    printf("      vel: %f %f %f\n",velocities[0][0],velocities[0][1],velocities[0][2]);

    // <FORCES
    vector<Vec3> forces;
    context.getForces(forces);
    printf("      for: %f %f %f\n",forces[0][0],forces[0][1],forces[0][2]);

    // <TIME
    double time = context.getTime();
    printf("      time: %f\n",time);

    // <NATOMS
    int natom = system.getNumParticles();
    printf("      natoms: %d\n",natom);

    // <CHARGES
    // identify the NonbondedForce, if any
    int nforce = system.getNumForces();
    printf("      nforces: %d\n",nforce);
    int nbnd_index = -1;
    for (int iforce=0; iforce < nforce; iforce++) {
      const Force& force = system.getForce(iforce);
      if ( dynamic_cast<const NonbondedForce*>( &force ) ){
	printf("      %d TRUE\n",iforce);
	nbnd_index = iforce;
      }
      else {
	printf("      %d FALSE\n",iforce);
      }
    }
    const Force& nbnd_temp = system.getForce(nbnd_index);
    const NonbondedForce* nbnd_force = dynamic_cast<const NonbondedForce*>( &nbnd_temp );
    double charge = 0.0;
    double sigma = 0.0;
    double epsilon = 0.0;
    nbnd_force->getParticleParameters(0, charge, sigma, epsilon);
    printf("      charge: %f\n",charge);
    //nbnd_force->setParticleParameters(0, 0.2, sigma, epsilon);

    // <DIMENSIONS
    bool periodic = nbnd_force->usesPeriodicBoundaryConditions();
    printf("      periodic: %d\n",periodic);

    // <CELL
    Vec3 cell1;
    Vec3 cell2;
    Vec3 cell3;
    context.getPeriodicBoxVectors(cell1, cell2, cell3);

    // >CELL
    context.getOwner().setPeriodicBoxVectors(cell1, cell2, cell3);

    // <ENERGY
    // Compiles, but should be called outside the time step
    //State state = context.getOwner().getState( State::Energy );
    //double ke = state.getKineticEnergy();
    //double pe = state.getPotentialEnergy();

    // <MASSES
    double mass;
    mass = system.getParticleMass(0);
    printf("      mass: %f\n",mass);


}
