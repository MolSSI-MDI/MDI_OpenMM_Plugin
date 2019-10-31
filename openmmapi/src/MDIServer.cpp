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
//#include "mdi.h"
//////////
#include "ExampleKernels.h"
#include "openmm/internal/ContextImpl.h"
#include <mpi.h>
//////////

using namespace ExamplePlugin;
using namespace OpenMM;
using namespace std;

MDIServer::MDIServer() {
}

void MDIServer::init(string mdi_options) {
    // Initialize MPI
    //MPI_Init(NULL, NULL);

    // Initialize MDI
    //const char *options = mdi_options.c_str();
    //MDI_Init(options, NULL);

    // Accept the MDI communicator
    //this->mdi_comm = MDI_Accept_Communicator();
    //printf("AAAAAAA: %d\n",this->mdi_comm);

    // register commands
    // this->register("<FORCES",["@GLOBAL","@UPDATE"]);

    // register callbacks
    // this->callbacks("@FORCES",["FORCES"]);
}

void MDIServer::listen(ContextImpl& context, Kernel& kernel, string node, MDI_Comm mdi_comm) {
    const OpenMM::System& system = context.getSystem();
    printf("CCCCCCCC %d\n",mdi_comm);

    // <COORDS
    this->send_coords(context, mdi_comm);

    // <VELOCITIES
    this->send_velocities(context, mdi_comm);

    // <FORCES
    this->send_forces(context, mdi_comm);

    // <TIME
    this->send_time(context, mdi_comm);

    // <NATOMS
    this->send_natoms(context, mdi_comm);

    // <CHARGES
    this->send_charges(context, mdi_comm);

    // <DIMENSIONS
    this->send_dimensions(context, mdi_comm);

    // <CELL
    this->send_cell(context, mdi_comm);

    // >CELL
    this->recv_cell(context, mdi_comm);

    // <ENERGY
    //this->send_energy(context, mdi_comm);

    // <MASSES
    this->send_masses(context, mdi_comm);

    // >FORCES
    this->add_forces(context, kernel, mdi_comm);
}



vector<Vec3> MDIServer::send_coords(ContextImpl& context, MDI_Comm mdi_comm) {
    const OpenMM::System& system = context.getSystem();
    int natoms = system.getNumParticles();
    vector<Vec3> positions;
    context.getPositions(positions);
    printf("      pos: %f %f %f\n",positions[0][0],positions[0][1],positions[0][2]);
    MDI_Send(&positions, 3*natoms, MDI_DOUBLE, mdi_comm);
    return positions;
}

vector<Vec3> MDIServer::send_velocities(ContextImpl& context, MDI_Comm mdi_comm) {
    const OpenMM::System& system = context.getSystem();
    int natoms = system.getNumParticles();
    vector<Vec3> velocities;
    context.getVelocities(velocities);
    printf("      vel: %f %f %f\n",velocities[0][0],velocities[0][1],velocities[0][2]);
    //MDI_Send(&velocities, 3*natoms, MDI_DOUBLE, this->mdi_comm);
    return velocities;
}

vector<Vec3> MDIServer::send_forces(ContextImpl& context, MDI_Comm mdi_comm) {
    const OpenMM::System& system = context.getSystem();
    int natoms = system.getNumParticles();
    vector<Vec3> forces;
    context.getForces(forces);
    printf("      for: %f %f %f\n",forces[0][0],forces[0][1],forces[0][2]);
    //MDI_Send(&forces, 3*natoms, MDI_DOUBLE, this->mdi_comm);
    return forces;
}

double MDIServer::send_time(ContextImpl& context, MDI_Comm mdi_comm) {
    double time = context.getTime();
    printf("      time: %f\n",time);
    //MDI_Send(&time, 1, MDI_DOUBLE, this->mdi_comm);
    return time;
}

int MDIServer::send_natoms(ContextImpl& context, MDI_Comm mdi_comm) {
    const OpenMM::System& system = context.getSystem();
    int natoms = system.getNumParticles();
    printf("      natoms: %d\n",natoms);
    //MDI_Send(&natoms, 1, MDI_INT, this->mdi_comm);
    return natoms;
}

vector<double> MDIServer::send_charges(ContextImpl& context, MDI_Comm mdi_comm) {
    const OpenMM::System& system = context.getSystem();
    int natoms = system.getNumParticles();

    // identify the NonbondedForce, if any
    const NonbondedForce* force = this->get_nonbonded_force(context);

    double charge = 0.0;
    double sigma = 0.0;
    double epsilon = 0.0;
    vector<double> charges;
    for (int iatom=0; iatom<natoms; iatom++) {
      force->getParticleParameters(iatom, charge, sigma, epsilon);
      charges.push_back(charge);
    }
    printf("      charge: %f\n",charge);
    //MDI_Send(&charges, natoms, MDI_DOUBLE, this->mdi_comm);
    return charges;
}

const NonbondedForce* MDIServer::get_nonbonded_force(ContextImpl& context) {
    const OpenMM::System& system = context.getSystem();

    // identify the NonbondedForce, if any
    int nforce = system.getNumForces();
    printf("      nforces: %d\n",nforce);
    int nbnd_index = -1;
    for (int iforce=0; iforce < nforce; iforce++) {
      const Force& force = system.getForce(iforce);
      if ( dynamic_cast<const NonbondedForce*>( &force ) ){
	nbnd_index = iforce;
      }
    }
    /*
    if (nbnd_index == -1) {
      throw "MDI Unable to find a nonbonded force";
    }
    */
    const Force& nbnd_temp = system.getForce(nbnd_index);
    const NonbondedForce* nbnd_force = dynamic_cast<const NonbondedForce*>( &nbnd_temp );
    return nbnd_force;
}

vector<int> MDIServer::send_dimensions(ContextImpl& context, MDI_Comm mdi_comm) {
    const OpenMM::System& system = context.getSystem();

    // identify the NonbondedForce, if any
    const NonbondedForce* force = this->get_nonbonded_force(context);

    bool periodic = force->usesPeriodicBoundaryConditions();
    printf("      periodic: %d\n",periodic);

    vector<int> dimensions;
    for (int idim=0; idim < 3; idim++) {
      if ( periodic ) {
	dimensions.push_back(2);
      }
      else {
	dimensions.push_back(1);
      }
    }

    return dimensions;
}

vector<double> MDIServer::send_cell(ContextImpl& context, MDI_Comm mdi_comm) {
    Vec3 cell1, cell2, cell3;
    context.getPeriodicBoxVectors(cell1, cell2, cell3);
    vector <double> cell { 
        cell1[0], cell1[1], cell1[2],
        cell2[0], cell2[1], cell2[2],
	cell3[0], cell3[1], cell3[2],
	0.0, 0.0, 0.0
    };
    return cell;
}

void MDIServer::recv_cell(ContextImpl& context, MDI_Comm mdi_comm) {
    // get the cell vectors
    // TEMPORARY
    vector<double> cell = this->send_cell(context, mdi_comm);
    Vec3 cell1, cell2, cell3;
    cell1[0] = cell[0];
    cell1[1] = cell[1];
    cell1[2] = cell[2];
    cell2[0] = cell[3];
    cell2[1] = cell[4];
    cell2[2] = cell[5];
    cell3[0] = cell[6];
    cell3[1] = cell[7];
    cell3[2] = cell[8];

    context.setPeriodicBoxVectors(cell1, cell2, cell3);
}

double MDIServer::send_energy(ContextImpl& context, MDI_Comm mdi_comm) {
    // Compiles, but should be called outside the time step
    State state = context.getOwner().getState( State::Energy );
    double ke = state.getKineticEnergy();
    double pe = state.getPotentialEnergy();
    double energy = ke + pe;
    return energy;
}

vector<double> MDIServer::send_masses(ContextImpl& context, MDI_Comm mdi_comm) {
    const OpenMM::System& system = context.getSystem();
    int natoms = system.getNumParticles();
    vector<double> masses;
    for (int iatom=0; iatom < natoms; iatom++) {
      masses.push_back( system.getParticleMass(iatom) );
    }
    printf("      mass: %f\n",masses[0]);
    return masses;
}

void MDIServer::add_forces(ContextImpl& context, Kernel& kernel, MDI_Comm mdi_comm) {
    // Tell the kernel to add a value to the forces
    const OpenMM::System& system = context.getSystem();
    int natoms = system.getNumParticles();

    kernel.getAs<CalcExampleForceKernel>().setAction(1);
    vector<double>* kernel_forces_ptr = kernel.getAs<CalcExampleForceKernel>().getForcesPtr();
    vector<double> &kernel_forces = *kernel_forces_ptr;

    for (int i = 0; i < 3*natoms; i++) {
      kernel_forces[i] = 100.0;
    }
}
