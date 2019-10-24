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

#ifdef WIN32
  #define _USE_MATH_DEFINES // Needed to get M_PI
#endif
#include "internal/ExampleForceImpl.h"
#include "ExampleKernels.h"
#include "openmm/OpenMMException.h"
#include "openmm/internal/ContextImpl.h"
#include <cmath>
#include <map>
#include <set>
#include <sstream>
//////////
#include "MDIServer.h"
#include "openmm/NonbondedForce.h"
//////////

using namespace ExamplePlugin;
using namespace OpenMM;
using namespace std;

ExampleForceImpl::ExampleForceImpl(const ExampleForce& owner) : owner(owner) {
}

ExampleForceImpl::~ExampleForceImpl() {
}

void ExampleForceImpl::initialize(ContextImpl& context) {
    kernel = context.getPlatform().createKernel(CalcExampleForceKernel::Name(), context);
    kernel.getAs<CalcExampleForceKernel>().initialize(context.getSystem(), owner);
}

double ExampleForceImpl::calcForcesAndEnergy(ContextImpl& context, bool includeForces, bool includeEnergy, int groups) {
    printf("   @FORCES\n");
    const OpenMM::System& system = context.getSystem();
    MDIServer& server = owner.getServer();

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


    if ((groups&(1<<owner.getForceGroup())) != 0)
        return kernel.getAs<CalcExampleForceKernel>().execute(context, includeForces, includeEnergy);
    return 0.0;
}

void ExampleForceImpl::updateContextState(OpenMM::ContextImpl& context, bool& forcesInvalid) {
    printf("   @UPDATE\n");
}


std::vector<std::string> ExampleForceImpl::getKernelNames() {
    std::vector<std::string> names;
    names.push_back(CalcExampleForceKernel::Name());
    return names;
}

vector<pair<int, int> > ExampleForceImpl::getBondedParticles() const {
    int numBonds = owner.getNumBonds();
    vector<pair<int, int> > bonds(numBonds);
    for (int i = 0; i < numBonds; i++) {
        double length, k;
        owner.getBondParameters(i, bonds[i].first, bonds[i].second, length, k);
    }
    return bonds;
}

void ExampleForceImpl::updateParametersInContext(ContextImpl& context) {
    kernel.getAs<CalcExampleForceKernel>().copyParametersToContext(context, owner);
}
