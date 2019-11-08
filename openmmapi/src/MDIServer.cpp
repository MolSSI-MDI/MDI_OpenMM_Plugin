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
#include "ExampleKernels.h"
#include "openmm/internal/ContextImpl.h"
#include <mpi.h>

using namespace ExamplePlugin;
using namespace OpenMM;
using namespace std;

MDIServer::MDIServer(string mdi_options) {
    // Initialize MPI
    MPI_Init(NULL, NULL);
    MPI_Comm world_comm = MPI_COMM_WORLD;

    // Initialize MDI
    int ierr;
    const char *options = mdi_options.c_str();
    printf("   Engine calling mdi_init\n");
    ierr = MDI_Init(options, &world_comm);
    if ( ierr != 0 ) {
      throw OpenMMException("Unable to initialize MDI\n");
    }
    // Register the @GLOBAL node
    MDI_Register_Node("@GLOBAL");
    MDI_Register_Command("@GLOBAL", "<COORDS");
    MDI_Register_Command("@GLOBAL", ">COORDS");
    MDI_Register_Command("@GLOBAL", "EXIT");
    MDI_Register_Command("@GLOBAL", "<VELOCITIES");
    MDI_Register_Command("@GLOBAL", ">VELOCITIES");
    MDI_Register_Command("@GLOBAL", "<FORCES");
    MDI_Register_Command("@GLOBAL", "<TIME");
    MDI_Register_Command("@GLOBAL", "<NATOMS");
    MDI_Register_Command("@GLOBAL", "<CHARGES");
    MDI_Register_Command("@GLOBAL", "<DIMENSIONS");
    MDI_Register_Command("@GLOBAL", "<CELL");
    MDI_Register_Command("@GLOBAL", ">CELL");
    MDI_Register_Command("@GLOBAL", "<MASSES");
    MDI_Register_Command("@GLOBAL", "@");
    MDI_Register_Command("@GLOBAL", "@INIT_MD");

    // Register the @INIT_MD node
    MDI_Register_Node("@INIT_MD");
    MDI_Register_Command("@INIT_MD", "<COORDS");
    MDI_Register_Command("@INIT_MD", ">COORDS");
    MDI_Register_Command("@INIT_MD", "EXIT");
    MDI_Register_Command("@INIT_MD", "<VELOCITIES");
    MDI_Register_Command("@INIT_MD", ">VELOCITIES");
    MDI_Register_Command("@INIT_MD", "<FORCES");
    MDI_Register_Command("@INIT_MD", "<TIME");
    MDI_Register_Command("@INIT_MD", "<NATOMS");
    MDI_Register_Command("@INIT_MD", "<CHARGES");
    MDI_Register_Command("@INIT_MD", "<DIMENSIONS");
    MDI_Register_Command("@INIT_MD", "<CELL");
    MDI_Register_Command("@INIT_MD", ">CELL");
    MDI_Register_Command("@INIT_MD", "<MASSES");
    MDI_Register_Command("@INIT_MD", "@");
    MDI_Register_Command("@INIT_MD", "@ENERGY");
    MDI_Register_Command("@INIT_MD", "@GLOBAL");
    MDI_Register_Command("@INIT_MD", "@FORCES");
    MDI_Register_Command("@INIT_MD", "@UPDATE");

    // Register the @UPDATE node
    MDI_Register_Node("@UPDATE");
    MDI_Register_Command("@UPDATE", "<COORDS");
    MDI_Register_Command("@UPDATE", ">COORDS");
    MDI_Register_Command("@UPDATE", "EXIT");
    MDI_Register_Command("@UPDATE", "<VELOCITIES");
    MDI_Register_Command("@UPDATE", ">VELOCITIES");
    MDI_Register_Command("@UPDATE", "<FORCES");
    MDI_Register_Command("@UPDATE", "<TIME");
    MDI_Register_Command("@UPDATE", "<NATOMS");
    MDI_Register_Command("@UPDATE", "<CHARGES");
    MDI_Register_Command("@UPDATE", "<DIMENSIONS");
    MDI_Register_Command("@UPDATE", "<CELL");
    MDI_Register_Command("@UPDATE", ">CELL");
    MDI_Register_Command("@UPDATE", "<MASSES");
    MDI_Register_Command("@UPDATE", "@");
    MDI_Register_Command("@UPDATE", "@ENERGY");
    MDI_Register_Command("@UPDATE", "@FORCES");
    MDI_Register_Command("@UPDATE", "@GLOBAL");
    MDI_Register_Command("@UPDATE", "@UPDATE");

    // Register the @FORCES node
    MDI_Register_Node("@FORCES");
    MDI_Register_Command("@FORCES", "<COORDS");
    MDI_Register_Command("@FORCES", "EXIT");
    MDI_Register_Command("@FORCES", "<VELOCITIES");
    MDI_Register_Command("@FORCES", "<NATOMS");
    MDI_Register_Command("@FORCES", "<CHARGES");
    MDI_Register_Command("@FORCES", "<DIMENSIONS");
    MDI_Register_Command("@FORCES", "<CELL");
    MDI_Register_Command("@FORCES", "<MASSES");
    MDI_Register_Command("@FORCES", "+FORCES");
    MDI_Register_Command("@FORCES", "@");
    MDI_Register_Command("@FORCES", "@ENERGY");
    MDI_Register_Command("@FORCES", "@FORCES");
    MDI_Register_Command("@FORCES", "@GLOBAL");
    MDI_Register_Command("@FORCES", "@UPDATE");

    // Register the @ENERGY node
    MDI_Register_Node("@ENERGY");
    MDI_Register_Command("@ENERGY", "<ENERGY");
    MDI_Register_Command("@ENERGY", "EXIT");
    MDI_Register_Command("@ENERGY", "<KE");
    MDI_Register_Command("@ENERGY", "<KE_NUC");
    MDI_Register_Command("@ENERGY", "<PE");
    MDI_Register_Command("@ENERGY", "<PE_NUC");
    MDI_Register_Command("@ENERGY", "@");
    MDI_Register_Command("@ENERGY", "@ENERGY");
    MDI_Register_Command("@ENERGY", "@FORCES");
    MDI_Register_Command("@ENERGY", "@GLOBAL");
    MDI_Register_Command("@ENERGY", "@UPDATE");

    // Accept the MDI communicator
    printf("   Engine calling mdi_accept_communicator\n");
    ierr = MDI_Accept_Communicator(&this->mdi_comm);
    if ( ierr != 0 ) {
      throw OpenMMException("Unable to accept communicator from the driver\n");
    }
    printf("   Engine finished init\n");

    // Initialize the target node
    this->target_node = new char[MDI_COMMAND_LENGTH];
    target_node[0] = '\0';
}

MDIServer::~MDIServer() {
    MPI_Finalize();
}

void MDIServer::setActive(bool active) {
  this->is_active = active;
}

std::string MDIServer::getTargetNode() {
  std::string target( this->target_node );
  return target;
}

std::string MDIServer::getPreviousNode() {
  return this->previous_node;
}

std::string MDIServer::listen(string node, ContextImpl& context, Kernel& kernel) {
    printf("   Engine at node: %s %d\n",node.c_str(), this->is_active);

    // Only enter server mode if the server is activated
    if ( not this->is_active ) {
      return "";
    }

    const OpenMM::System& system = context.getSystem();
    int supported;
    int ierr;
    char *command = new char[MDI_COMMAND_LENGTH];

    // If there is a target node, check if this is the target node
    if ( strcmp( this->target_node, "" ) != 0 ) {
      if ( strcmp( this->target_node, "@" ) == 0 ) {
	// The engine was searching for the next node, so this node counts
	target_node[0] = '\0';
      }
      else if ( strcmp( this->target_node, node.c_str() ) == 0 ) {
	// This is the node the engine was searching for
	target_node[0] = '\0';
      }
    }

    // Enter server mode, listening for commands from the driver
    while ( strcmp( this->target_node, "" ) == 0 ) {

      // Receive a new command from the driver
      ierr = MDI_Recv_Command(command, this->mdi_comm);
      if ( ierr != 0 ) {
	throw OpenMMException("Unable to receive command from the driver\n");
      }
      printf("   MDI COMMAND: %s\n",command);

      // Confirm that this command is supported
      MDI_Check_Command_Exists(node.c_str(), command, MDI_NULL_COMM, &supported);
      if ( supported != 1 ) {
	throw OpenMMException("Received unsupported MDI command\n");
      }

      // Respond to the command
      if ( strcmp( command, "<CELL" )  == 0 ) {
	this->send_cell(context);
      }
      else if ( strcmp( command, ">CELL" )  == 0 ) {
	this->recv_cell(context);
      }
      else if ( strcmp( command, "<CHARGES" ) == 0 ) {
	this->send_charges(context);
      }
      else if ( strcmp( command, "<COORDS" )  == 0 ) {
	this->send_coords(context);
      }
      else if ( strcmp( command, ">COORDS" ) == 0 ) {
	this->recv_coords(context);
      }
      else if ( strcmp( command, "<DIMENSIONS" )  == 0 ) {
	this->send_dimensions(context);
      }
      else if ( strcmp( command, "<ENERGY" ) == 0 ) {
	this->send_energy(context);
      }
      else if ( strcmp( command, "EXIT" ) == 0 ) {
	strcpy( this->target_node, command );
      }
      else if ( strcmp( command, "<FORCES" ) == 0 ) {
	this->send_forces(context);
      }
      else if ( strcmp( command, "+FORCES" ) == 0 ) {
	this->add_forces(context, kernel);
      }
      else if ( strcmp( command, "<KE" ) == 0 ) {
	this->send_ke(context);
      }
      else if ( strcmp( command, "<KE_NUC" ) == 0 ) {
	this->send_ke_nuc(context);
      }
      else if ( strcmp( command, "<MASSES" ) == 0 ) {
	this->send_masses(context);
      }
      else if ( strcmp( command, "<NATOMS" ) == 0 ) {
	this->send_natoms(context);
      }
      else if ( strcmp( command, "<PE" ) == 0 ) {
	this->send_pe(context);
      }
      else if ( strcmp( command, "<PE_NUC" ) == 0 ) {
	this->send_pe_nuc(context);
      }
      else if ( strcmp( command, "<TIME" ) == 0 ) {
	this->send_time(context);
      }
      else if ( strcmp( command, "<VELOCITIES" ) == 0 ) {
        this->send_velocities(context);
      }
      else if ( strcmp( command, ">VELOCITIES" ) == 0 ) {
	this->recv_velocities(context);
      }
      else if ( strcmp( command, "@" ) == 0 ) {
	strcpy( this->target_node, command );
      }
      else if ( strcmp( command, "@ENERGY" ) == 0 ) {
	strcpy( this->target_node, command );
      }
      else if ( strcmp( command, "@FORCES" ) == 0 ) {
	strcpy( this->target_node, command );
      }
      else if ( strcmp( command, "@GLOBAL" ) == 0 ) {
	strcpy( this->target_node, command );
      }
      else if ( strcmp( command, "@INIT_MD" ) == 0 ) {
	strcpy( this->target_node, command );
      }
      else if ( strcmp( command, "@UPDATE" ) == 0 ) {
	strcpy( this->target_node, command );
      }
    }

    this->previous_node = node;

    delete [] command;

    return this->target_node;
}



vector<Vec3> MDIServer::send_coords(ContextImpl& context) {
    const OpenMM::System& system = context.getSystem();
    int natoms = system.getNumParticles();
    vector<Vec3> positions;
    context.getPositions(positions);
    double conv;
    MDI_Conversion_Factor("nanometer","atomic_unit_of_length",&conv);
    //vector<double> coords;
    //coords.resize( 3 * natoms );
    for (int iatom = 0; iatom < natoms; iatom++) {
      positions[iatom][0] *= conv;
      positions[iatom][1] *= conv;
      positions[iatom][2] *= conv;
      /*
      coords[3*iatom+0] = conv * positions[iatom][0];
      coords[3*iatom+1] = conv * positions[iatom][1];
      coords[3*iatom+2] = conv * positions[iatom][2];
      */
    }
    printf("      pos: %f %f %f\n",positions[0][0],positions[0][1],positions[0][2]);
    MDI_Send(&positions[0][0], 3*natoms, MDI_DOUBLE, mdi_comm);
    return positions;
}

void MDIServer::recv_coords(ContextImpl& context, vector<Vec3>* coords_in) {
    const OpenMM::System& system = context.getSystem();
    int natoms = system.getNumParticles();

    vector<Vec3> positions;
    positions.resize(natoms);
    if ( coords_in == nullptr ) {
      MDI_Recv(&positions[0][0], 3*natoms, MDI_DOUBLE, mdi_comm);
    }
    else {
      for (int i = 0; i < natoms; i++) {
        positions[i] = (*coords_in)[i];
      }
    }

    double conv;
    MDI_Conversion_Factor("atomic_unit_of_length","nanometer",&conv);
    for (int i = 0; i < natoms; i++) {
      positions[i][0] *= conv;
      positions[i][1] *= conv;
      positions[i][2] *= conv;
    }
    context.setPositions(positions);
}

vector<Vec3> MDIServer::send_velocities(ContextImpl& context) {
    const OpenMM::System& system = context.getSystem();
    int natoms = system.getNumParticles();
    vector<Vec3> velocities;
    context.getVelocities(velocities);
    double conv;
    MDI_Conversion_Factor("nanometer","atomic_unit_of_length",&conv);
    double len_conv;
    MDI_Conversion_Factor("picosecond","atomic_unit_of_time",&len_conv);
    conv /= len_conv;
    for (int iatom = 0; iatom < natoms; iatom++) {
      velocities[iatom][0] *= conv;
      velocities[iatom][1] *= conv;
      velocities[iatom][2] *= conv;
    }
    printf("      vel: %f %f %f\n",velocities[0][0],velocities[0][1],velocities[0][2]);
    MDI_Send(&velocities[0][0], 3*natoms, MDI_DOUBLE, mdi_comm);
    return velocities;
}

void MDIServer::recv_velocities(ContextImpl& context, vector<Vec3>* velocities_in) {
    const OpenMM::System& system = context.getSystem();
    int natoms = system.getNumParticles();

    vector<Vec3> velocities;
    velocities.resize(natoms);
    if ( velocities_in == nullptr ) {
      MDI_Recv(&velocities[0][0], 3*natoms, MDI_DOUBLE, mdi_comm);
    }
    else {
      for (int i = 0; i < natoms; i++) {
        velocities[i] = (*velocities_in)[i];
      }
    }

    double conv;
    MDI_Conversion_Factor("atomic_unit_of_length","nanometer",&conv);
    double time_conv;
    MDI_Conversion_Factor("atomic_unit_of_time","picosecond",&time_conv);
    conv /= time_conv;
    for (int i = 0; i < natoms; i++) {
      velocities[i][0] *= conv;
      velocities[i][1] *= conv;
      velocities[i][2] *= conv;
    }
    context.setVelocities(velocities);
}

vector<Vec3> MDIServer::send_forces(ContextImpl& context) {
    const OpenMM::System& system = context.getSystem();
    int natoms = system.getNumParticles();
    vector<Vec3> forces;
    context.getForces(forces);
    double conv;
    MDI_Conversion_Factor("kilojoule_per_mol","atomic_unit_of_energy",&conv);
    double len_conv;
    MDI_Conversion_Factor("nanometer","atomic_unit_of_length",&len_conv);
    conv *= len_conv;
    for (int iatom = 0; iatom < natoms; iatom++) {
      forces[iatom][0] *= conv;
      forces[iatom][1] *= conv;
      forces[iatom][2] *= conv;
    }
    printf("      for: %f %f %f\n",forces[0][0],forces[0][1],forces[0][2]);
    MDI_Send(&forces[0][0], 3*natoms, MDI_DOUBLE, mdi_comm);
    return forces;
}

double MDIServer::send_time(ContextImpl& context) {
    double time = context.getTime();
    double conv;
    MDI_Conversion_Factor("picosecond","atomic_unit_of_time",&conv);
    time *= conv;
    printf("      time: %f\n",time);
    MDI_Send(&time, 1, MDI_DOUBLE, mdi_comm);
    return time;
}

int MDIServer::send_natoms(ContextImpl& context) {
    const OpenMM::System& system = context.getSystem();
    int natoms = system.getNumParticles();
    printf("      natoms: %d\n",natoms);
    MDI_Send(&natoms, 1, MDI_INT, mdi_comm);
    return natoms;
}

vector<double> MDIServer::send_charges(ContextImpl& context) {
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
    MDI_Send(&charges[0], natoms, MDI_DOUBLE, mdi_comm);
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

vector<int> MDIServer::send_dimensions(ContextImpl& context) {
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
    MDI_Send(&dimensions[0], 3, MDI_INT, mdi_comm);

    return dimensions;
}

vector<double> MDIServer::send_cell(ContextImpl& context) {
    Vec3 cell1, cell2, cell3;
    double conv;
    MDI_Conversion_Factor("nanometer","atomic_unit_of_length",&conv);
    context.getPeriodicBoxVectors(cell1, cell2, cell3);
    vector <double> cell { 
        cell1[0] * conv, cell1[1] * conv, cell1[2] * conv,
        cell2[0] * conv, cell2[1] * conv, cell2[2] * conv,
	cell3[0] * conv, cell3[1] * conv, cell3[2] * conv,
	0.0, 0.0, 0.0
    };
    MDI_Send(&cell[0], 12, MDI_DOUBLE, mdi_comm);
    return cell;
}

double MDIServer::send_energy(ContextImpl& context) {
    // Compiles, but should be called outside the time step
    this->setActive(false);
    State state = context.getOwner().getState( State::Energy );
    double ke = state.getKineticEnergy();
    double pe = state.getPotentialEnergy();
    double energy = ke + pe;
    double conv;
    MDI_Conversion_Factor("kilojoule_per_mol","atomic_unit_of_energy",&conv);
    energy *= conv;
    MDI_Send(&energy, 1, MDI_DOUBLE, mdi_comm);
    this->setActive(true);
    return energy;
}

double MDIServer::send_ke(ContextImpl& context) {
    // Compiles, but should be called outside the time step
    this->setActive(false);
    State state = context.getOwner().getState( State::Energy );
    double ke = state.getKineticEnergy();
    double conv;
    MDI_Conversion_Factor("kilojoule_per_mol","atomic_unit_of_energy",&conv);
    ke *= conv;
    MDI_Send(&ke, 1, MDI_DOUBLE, mdi_comm);
    this->setActive(true);
    return ke;
}

double MDIServer::send_ke_nuc(ContextImpl& context) {
    return this->send_ke(context);
}

double MDIServer::send_pe(ContextImpl& context) {
    // Compiles, but should be called outside the time step
    this->setActive(false);
    State state = context.getOwner().getState( State::Energy );
    double pe = state.getPotentialEnergy();
    double conv;
    MDI_Conversion_Factor("kilojoule_per_mol","atomic_unit_of_energy",&conv);
    pe *= conv;
    MDI_Send(&pe, 1, MDI_DOUBLE, mdi_comm);
    this->setActive(true);
    return pe;
}

double MDIServer::send_pe_nuc(ContextImpl& context) {
    return this->send_pe(context);
}

vector<double> MDIServer::send_masses(ContextImpl& context) {
    const OpenMM::System& system = context.getSystem();
    int natoms = system.getNumParticles();
    vector<double> masses;
    for (int iatom=0; iatom < natoms; iatom++) {
      masses.push_back( system.getParticleMass(iatom) );
    }
    printf("      mass: %f\n",masses[0]);
    MDI_Send(&masses[0], natoms, MDI_DOUBLE, mdi_comm);
    return masses;
}

void MDIServer::recv_cell(ContextImpl& context, vector<double>* cell_in) {
    // get the cell vectors
    vector<double> cell;
    cell.resize(12);
    if ( cell_in == nullptr ) {
      MDI_Recv(&cell[0], 12, MDI_DOUBLE, mdi_comm);
    }
    else {
      for (int i = 0; i < cell.size(); i++ ) {
	cell[i] = (*cell_in)[i];
      }
    }
    double conv;
    MDI_Conversion_Factor("atomic_unit_of_length","nanometer",&conv);
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

void MDIServer::add_forces(ContextImpl& context, Kernel& kernel, vector<double>* forces_in) {
    const OpenMM::System& system = context.getSystem();
    int natoms = system.getNumParticles();

    // Tell the kernel to add a value to the forces
    kernel.getAs<CalcExampleForceKernel>().setAction(1);
    vector<double>* kernel_forces_ptr = kernel.getAs<CalcExampleForceKernel>().getForcesPtr();
    vector<double> &kernel_forces = *kernel_forces_ptr;

    vector<double> forces;
    forces.resize(3*natoms);
    if ( forces_in == nullptr ) {
      MDI_Recv(&kernel_forces[0], 3*natoms, MDI_DOUBLE, mdi_comm);
    }
    else {
      for (int i = 0; i < 3*natoms; i++) {
	kernel_forces[i] = (*forces_in)[i];
      }
    }

    double conv;
    MDI_Conversion_Factor("atomic_unit_of_energy","kilojoule_per_mol",&conv);
    for (int i = 0; i < 3*natoms; i++) {
      kernel_forces[i] *= conv;
    }
}
