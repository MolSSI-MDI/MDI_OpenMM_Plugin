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

    // <COORDS
    vector<Vec3> positions = this->send_coords(context, mdi_comm);

    // >COORDS
    this->recv_coords(context, mdi_comm, &positions);

    // <VELOCITIES
    vector<Vec3> velocities = this->send_velocities(context, mdi_comm);

    // >VELOCITIES
    this->recv_velocities(context, mdi_comm, &velocities);

    // <FORCES
    this->send_forces(context, mdi_comm);

    // <TIME
    this->send_time(context, mdi_comm);

    // <NATOMS
    double natoms = this->send_natoms(context, mdi_comm);

    // <CHARGES
    this->send_charges(context, mdi_comm);

    // <DIMENSIONS
    this->send_dimensions(context, mdi_comm);

    // <CELL
    this->send_cell(context, mdi_comm);

    // >CELL
    vector<double> cell = this->send_cell(context, mdi_comm);
    this->recv_cell(context, mdi_comm, &cell);

    // <ENERGY
    //this->send_energy(context, mdi_comm);

    // <KE
    //this->send_ke(context, mdi_comm);

    // <KE_NUC
    //this->send_ke_nuc(context, mdi_comm);

    // <PE
    //this->send_pe(context, mdi_comm);

    // <PE_NUC
    //this->send_pe_nuc(context, mdi_comm);

    // <MASSES
    this->send_masses(context, mdi_comm);

    // >FORCES
    double conv = MDI_Conversion_Factor("atomic_unit_of_energy","kilojoule_per_mol");
    vector<double> forces;
    for (int i = 0; i < 3*natoms; i++) {
      forces.push_back( 100.0 / conv );
    }
    this->add_forces(context, kernel, mdi_comm, &forces);
}



vector<Vec3> MDIServer::send_coords(ContextImpl& context, MDI_Comm mdi_comm) {
    const OpenMM::System& system = context.getSystem();
    int natoms = system.getNumParticles();
    vector<Vec3> positions;
    context.getPositions(positions);
    double conv = MDI_Conversion_Factor("nanometer","atomic_unit_of_length");
    for (int iatom = 0; iatom < natoms; iatom++) {
      positions[iatom][0] *= conv;
      positions[iatom][1] *= conv;
      positions[iatom][2] *= conv;
    }
    printf("      pos: %f %f %f\n",positions[0][0],positions[0][1],positions[0][2]);
    MDI_Send(&positions, 3*natoms, MDI_DOUBLE, mdi_comm);
    return positions;
}

void MDIServer::recv_coords(ContextImpl& context, MDI_Comm mdi_comm, vector<Vec3>* coords_in) {
    const OpenMM::System& system = context.getSystem();
    int natoms = system.getNumParticles();

    vector<Vec3> positions;
    positions.resize(natoms);
    if ( coords_in == nullptr ) {
      MDI_Recv(&positions, 3*natoms, MDI_DOUBLE, mdi_comm);
    }
    else {
      for (int i = 0; i < natoms; i++) {
        positions[i] = (*coords_in)[i];
      }
    }

    double conv = MDI_Conversion_Factor("atomic_unit_of_length","nanometer");
    for (int i = 0; i < natoms; i++) {
      positions[i][0] *= conv;
      positions[i][1] *= conv;
      positions[i][2] *= conv;
    }
    context.setPositions(positions);
}

vector<Vec3> MDIServer::send_velocities(ContextImpl& context, MDI_Comm mdi_comm) {
    const OpenMM::System& system = context.getSystem();
    int natoms = system.getNumParticles();
    vector<Vec3> velocities;
    context.getVelocities(velocities);
    double conv = MDI_Conversion_Factor("nanometer","atomic_unit_of_length");
    conv /= MDI_Conversion_Factor("picosecond","atomic_unit_of_time");
    for (int iatom = 0; iatom < natoms; iatom++) {
      velocities[iatom][0] *= conv;
      velocities[iatom][1] *= conv;
      velocities[iatom][2] *= conv;
    }
    printf("      vel: %f %f %f\n",velocities[0][0],velocities[0][1],velocities[0][2]);
    MDI_Send(&velocities, 3*natoms, MDI_DOUBLE, mdi_comm);
    return velocities;
}

void MDIServer::recv_velocities(ContextImpl& context, MDI_Comm mdi_comm, vector<Vec3>* velocities_in) {
    const OpenMM::System& system = context.getSystem();
    int natoms = system.getNumParticles();

    vector<Vec3> velocities;
    velocities.resize(natoms);
    if ( velocities_in == nullptr ) {
      MDI_Recv(&velocities, 3*natoms, MDI_DOUBLE, mdi_comm);
    }
    else {
      for (int i = 0; i < natoms; i++) {
        velocities[i] = (*velocities_in)[i];
      }
    }

    double conv = MDI_Conversion_Factor("atomic_unit_of_length","nanometer");
    conv /= MDI_Conversion_Factor("atomic_unit_of_time","picosecond");
    for (int i = 0; i < natoms; i++) {
      velocities[i][0] *= conv;
      velocities[i][1] *= conv;
      velocities[i][2] *= conv;
    }
    context.setVelocities(velocities);
}

vector<Vec3> MDIServer::send_forces(ContextImpl& context, MDI_Comm mdi_comm) {
    const OpenMM::System& system = context.getSystem();
    int natoms = system.getNumParticles();
    vector<Vec3> forces;
    context.getForces(forces);
    double conv = MDI_Conversion_Factor("kilojoule_per_mol","atomic_unit_of_energy");
    conv *= MDI_Conversion_Factor("nanometer","atomic_unit_of_length");
    for (int iatom = 0; iatom < natoms; iatom++) {
      forces[iatom][0] *= conv;
      forces[iatom][1] *= conv;
      forces[iatom][2] *= conv;
    }
    printf("      for: %f %f %f\n",forces[0][0],forces[0][1],forces[0][2]);
    MDI_Send(&forces, 3*natoms, MDI_DOUBLE, mdi_comm);
    return forces;
}

double MDIServer::send_time(ContextImpl& context, MDI_Comm mdi_comm) {
    double time = context.getTime();
    double conv = MDI_Conversion_Factor("picosecond","atomic_unit_of_time");
    time *= conv;
    printf("      time: %f\n",time);
    MDI_Send(&time, 1, MDI_DOUBLE, mdi_comm);
    return time;
}

int MDIServer::send_natoms(ContextImpl& context, MDI_Comm mdi_comm) {
    const OpenMM::System& system = context.getSystem();
    int natoms = system.getNumParticles();
    printf("      natoms: %d\n",natoms);
    MDI_Send(&natoms, 1, MDI_INT, mdi_comm);
    return natoms;
}

vector<double> MDIServer::send_charges(ContextImpl& context, MDI_Comm mdi_comm) {
    const OpenMM::System& system = context.getSystem();
    int natoms = system.getNumParticles();

    // identify the NonbondedForce, if any
    const NonbondedForce* force = get_nonbonded_force(context);

    double charge = 0.0;
    double sigma = 0.0;
    double epsilon = 0.0;
    vector<double> charges;
    for (int iatom=0; iatom<natoms; iatom++) {
      force->getParticleParameters(iatom, charge, sigma, epsilon);
      charges.push_back(charge);
    }
    printf("      charge: %f\n",charge);
    MDI_Send(&charges, natoms, MDI_DOUBLE, mdi_comm);
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
    if (nbnd_index == -1) {
      throw "MDI Unable to find a nonbonded force";
    }
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
    MDI_Send(&dimensions, 3, MDI_INT, mdi_comm);

    return dimensions;
}

vector<double> MDIServer::send_cell(ContextImpl& context, MDI_Comm mdi_comm) {
    Vec3 cell1, cell2, cell3;
    double conv = MDI_Conversion_Factor("nanometer","atomic_unit_of_length");
    context.getPeriodicBoxVectors(cell1, cell2, cell3);
    vector <double> cell { 
        cell1[0] * conv, cell1[1] * conv, cell1[2] * conv,
        cell2[0] * conv, cell2[1] * conv, cell2[2] * conv,
	cell3[0] * conv, cell3[1] * conv, cell3[2] * conv,
	0.0, 0.0, 0.0
    };
    MDI_Send(&cell, 12, MDI_DOUBLE, mdi_comm);
    return cell;
}

double MDIServer::send_energy(ContextImpl& context, MDI_Comm mdi_comm) {
    // Compiles, but should be called outside the time step
    State state = context.getOwner().getState( State::Energy );
    double ke = state.getKineticEnergy();
    double pe = state.getPotentialEnergy();
    double energy = ke + pe;
    double conv = MDI_Conversion_Factor("kilojoule_per_mol","atomic_unit_of_energy");
    energy *= conv;
    MDI_Send(&energy, 1, MDI_DOUBLE, mdi_comm);
    return energy;
}

double MDIServer::send_ke(ContextImpl& context, MDI_Comm mdi_comm) {
    // Compiles, but should be called outside the time step
    State state = context.getOwner().getState( State::Energy );
    double ke = state.getKineticEnergy();
    double conv = MDI_Conversion_Factor("kilojoule_per_mol","atomic_unit_of_energy");
    ke *= conv;
    MDI_Send(&ke, 1, MDI_DOUBLE, mdi_comm);
    return ke;
}

double MDIServer::send_ke_nuc(ContextImpl& context, MDI_Comm mdi_comm) {
    return this->send_ke(context, mdi_comm);
}

double MDIServer::send_pe(ContextImpl& context, MDI_Comm mdi_comm) {
    // Compiles, but should be called outside the time step
    State state = context.getOwner().getState( State::Energy );
    double pe = state.getPotentialEnergy();
    double conv = MDI_Conversion_Factor("kilojoule_per_mol","atomic_unit_of_energy");
    pe *= conv;
    MDI_Send(&pe, 1, MDI_DOUBLE, mdi_comm);
    return pe;
}

double MDIServer::send_pe_nuc(ContextImpl& context, MDI_Comm mdi_comm) {
    return this->send_pe(context, mdi_comm);
}

vector<double> MDIServer::send_masses(ContextImpl& context, MDI_Comm mdi_comm) {
    const OpenMM::System& system = context.getSystem();
    int natoms = system.getNumParticles();
    vector<double> masses;
    for (int iatom=0; iatom < natoms; iatom++) {
      masses.push_back( system.getParticleMass(iatom) );
    }
    printf("      mass: %f\n",masses[0]);
    MDI_Send(&masses, natoms, MDI_DOUBLE, mdi_comm);
    return masses;
}

void MDIServer::recv_cell(ContextImpl& context, MDI_Comm mdi_comm, vector<double>* cell_in) {
    // get the cell vectors
    vector<double> cell;
    cell.resize(12);
    if ( cell_in == nullptr ) {
      MDI_Recv(&cell, 12, MDI_DOUBLE, mdi_comm);
    }
    else {
      for (int i = 0; i < cell.size(); i++ ) {
	cell[i] = (*cell_in)[i];
      }
    }
    double conv = MDI_Conversion_Factor("atomic_unit_of_length","nanometer");
    Vec3 cell1, cell2, cell3;
    cell1[0] = cell[0] * conv;
    cell1[1] = cell[1] * conv;
    cell1[2] = cell[2] * conv;
    cell2[0] = cell[3] * conv;
    cell2[1] = cell[4] * conv;
    cell2[2] = cell[5] * conv;
    cell3[0] = cell[6] * conv;
    cell3[1] = cell[7] * conv;
    cell3[2] = cell[8] * conv;

    context.setPeriodicBoxVectors(cell1, cell2, cell3);
}

void MDIServer::add_forces(ContextImpl& context, Kernel& kernel, MDI_Comm mdi_comm, vector<double>* forces_in) {
    const OpenMM::System& system = context.getSystem();
    int natoms = system.getNumParticles();

    // Tell the kernel to add a value to the forces
    kernel.getAs<CalcExampleForceKernel>().setAction(1);
    vector<double>* kernel_forces_ptr = kernel.getAs<CalcExampleForceKernel>().getForcesPtr();
    vector<double> &kernel_forces = *kernel_forces_ptr;

    vector<double> forces;
    forces.resize(3*natoms);
    if ( forces_in == nullptr ) {
      MDI_Recv(&kernel_forces, 3*natoms, MDI_DOUBLE, mdi_comm);
    }
    else {
      for (int i = 0; i < 3*natoms; i++) {
	kernel_forces[i] = (*forces_in)[i];
      }
    }

    double conv = MDI_Conversion_Factor("atomic_unit_of_energy","kilojoule_per_mol");
    for (int i = 0; i < 3*natoms; i++) {
      kernel_forces[i] *= conv;
    }
}
