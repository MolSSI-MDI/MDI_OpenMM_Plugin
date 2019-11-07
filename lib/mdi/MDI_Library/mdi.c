/*! \file
 *
 * \brief Functions callable by users of the MolSSI Driver Interface
 */



/* ----------------------------------------------------------------------
   MDI - MolSSI Driver Interface
   https://molssi.org/, Molecular Sciences Software Institute
   Taylor Barnes, tbarnes1@vt.edu
-------------------------------------------------------------------------

Contents:
   MDI_Init: Initialize MDI
   MDI_Accept_Communicator: Accepts a new MDI communicator
   MDI_Send: Sends data through the socket
   MDI_Recv: Receives data through the socket
   MDI_Send_Command: Sends a string of length MDI_COMMAND_LENGTH through the
      socket
   MDI_Recv_Command: Receives a string of length MDI_COMMAND_LENGTH through the
      socket
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "mdi.h"
#include "mdi_global.h"
#include "mdi_general.h"
#include "mdi_mpi.h"

/*! \brief MDI version number */
const double MDI_VERSION = 0.5;

/*! \brief length of an MDI command in characters */
const int MDI_COMMAND_LENGTH = 12;

/*! \brief length of an MDI name in characters */
const int MDI_NAME_LENGTH = 12;

/*! \brief value of a null communicator */
const MDI_Comm MDI_NULL_COMM = 0;

// MDI data types
/*! \brief integer data type */
const int MDI_INT          = 0;
/*! \brief double precision float data type */
const int MDI_DOUBLE       = 1;
/*! \brief character data type */
const int MDI_CHAR         = 2;
/*! \brief NumPy integer data type */
const int MDI_INT_NUMPY    = 3;
/*! \brief NumPy double precision float data type */
const int MDI_DOUBLE_NUMPY = 4;

// MDI communication types
/*! \brief TCP/IP communication method */
const int MDI_TCP    = 1;
/*! \brief MPI communication method */
const int MDI_MPI    = 2;
/*! \brief Test communication method */
const int MDI_TEST   = 3;

/*----------------------*/
/* MDI unit conversions */
/*----------------------*/

// length
/*! \brief conversion factor between meters and bohr */
const double MDI_METER_TO_BOHR = 1.88972612546e10;
/*! \brief conversion factor between Angstroms and bohr */
const double MDI_ANGSTROM_TO_BOHR = 1.88972612546;

// time
/*! \brief conversion factor between seconds and atomic units of time */
const double MDI_SECOND_TO_AUT = 4.1341374575751e16;
/*! \brief conversion factor between picoseconds and atomic units of time */
const double MDI_PICOSECOND_TO_AUT = 4.1341374575751e4;

// force
/*! \brief conversion factor between newtons and atomic units of force */
const double MDI_NEWTON_TO_AUF = 1.213780478e7;

// energy
/*! \brief conversion factor between joules and hartrees */
const double MDI_JOULE_TO_HARTREE = 2.29371265835792e17;
/*! \brief conversion factor between kilojoules and hartrees */
const double MDI_KJ_TO_HARTREE = 2.29371265835792e20;
/*! \brief conversion factor between kilojoules/mol and hartrees */
const double MDI_KJPERMOL_TO_HARTREE = 3.80879947807451e-4;
/*! \brief conversion factor between kcal/mol and hartrees */
const double MDI_KCALPERMOL_TO_HARTREE = 1.5941730215480900e-3;
/*! \brief conversion factor between eV and hartrees */
const double MDI_EV_TO_HARTREE = 3.67493266806491e-2;
/*! \brief conversion factor between rydbergs and hartrees */
const double MDI_RYDBERG_TO_HARTREE = 0.5;
/*! \brief conversion factor between Kelvin and hartrees */
const double MDI_KELVIN_TO_HARTREE = 3.16681050847798e-6;


static int is_initialized = 0;


/*! \brief Initialize communication through the MDI library
 *
 * If using the "-method MPI" option, this function must be called by all ranks.
 * The function returns \p 0 on a success.
 *
 * \param [in]       options
 *                   Options describing the communication method used to connect to codes.
 * \param [in, out]  world_comm
 *                   On input, the MPI communicator that spans all of the codes.
 *                   On output, the MPI communicator that spans the single code corresponding to the calling rank.
 *                   Only used if the "-method MPI" option is provided.
 */
int MDI_Init(const char* options, void* world_comm)
{
  if ( is_initialized == 1 ) {
    mdi_error("MDI_Init called after MDI was already initialized");
  }
  general_init(options, world_comm);
  is_initialized = 1;
  return 0;
}


/*! \brief Accept a new MDI communicator
 *
 * The function returns an MDI_Comm that describes a connection between two codes.
 * If no new communicators are available, the function returns \p MDI_NULL_COMM.
 *
 */
MDI_Comm MDI_Accept_Communicator(MDI_Comm* comm)
{
  if ( is_initialized == 0 ) {
    mdi_error("MDI_Accept_Communicator called but MDI has not been initialized");
  }
  *comm = general_accept_communicator();
  return 0;
}


/*! \brief Send data through the MDI connection
 *
 * If running with MPI, this function must be called only by rank \p 0.
 * The function returns \p 0 on a success.
 *
 * \param [in]       buf
 *                   Pointer to the data to be sent.
 * \param [in]       count
 *                   Number of values (integers, double precision floats, characters, etc.) to be sent.
 * \param [in]       datatype
 *                   MDI handle (MDI_INT, MDI_DOUBLE, MDI_CHAR, etc.) corresponding to the type of data to be sent.
 * \param [in]       comm
 *                   MDI communicator associated with the intended recipient code.
 */
int MDI_Send(const void* buf, int count, MDI_Datatype datatype, MDI_Comm comm)
{
  if ( is_initialized == 0 ) {
    mdi_error("MDI_Send called but MDI has not been initialized");
  }
  return general_send(buf, count, datatype, comm);
}


/*! \brief Receive data through the MDI connection
 *
 * If running with MPI, this function must be called only by rank \p 0.
 * The function returns \p 0 on a success.
 *
 * \param [in]       buf
 *                   Pointer to the buffer where the received data will be stored.
 * \param [in]       count
 *                   Number of values (integers, double precision floats, characters, etc.) to be received.
 * \param [in]       datatype
 *                   MDI handle (MDI_INT, MDI_DOUBLE, MDI_CHAR, etc.) corresponding to the type of data to be received.
 * \param [in]       comm
 *                   MDI communicator associated with the connection to the sending code.
 */
int MDI_Recv(void* buf, int count, MDI_Datatype datatype, MDI_Comm comm)
{
  if ( is_initialized == 0 ) {
    mdi_error("MDI_Recv called but MDI has not been initialized");
  }
  return general_recv(buf, count, datatype, comm);
}


/*! \brief Send a command of length \p MDI_COMMAND_LENGTH through the MDI connection
 *
 * If running with MPI, this function must be called only by rank \p 0.
 * The function returns \p 0 on a success.
 *
 * \param [in]       buf
 *                   Pointer to the data to be sent.
 * \param [in]       comm
 *                   MDI communicator associated with the intended recipient code.
 */
int MDI_Send_Command(const char* buf, MDI_Comm comm)
{
  if ( is_initialized == 0 ) {
    mdi_error("MDI_Send_Command called but MDI has not been initialized");
  }
  return general_send_command(buf, comm);
}


/*! \brief Receive a command of length \p MDI_COMMAND_LENGTH through the MDI connection
 *
 * If running with MPI, this function must be called only by rank \p 0.
 * The function returns \p 0 on a success.
 *
 * \param [in]       buf
 *                   Pointer to the buffer where the received data will be stored.
 * \param [in]       comm
 *                   MDI communicator associated with the connection to the sending code.
 */
int MDI_Recv_Command(char* buf, MDI_Comm comm)
{
  if ( is_initialized == 0 ) {
    mdi_error("MDI_Recv_Command called but MDI has not been initialized");
  }
  return general_recv_command(buf, comm);
}


/*! \brief Determine the conversion factor between two units
 *
 * The function determines the conversion factor from \p in_unit to \p out_unit.
 * The function requires that \p in_unit and \p out_unit be members of the same category of unit (\em i.e. charge, energy, force, etc.).
 * For example, calling \p MDI_Conversion_Factor(\p "kilojoule_per_mol",\p "atomic_unit_of_energy") will return the conversion factor from kilojoule/mol to hartrees.
 * The function returns \p 0 on a success.
 *
 * All quantities communicated through MDI must be represented using atomic units.
 * When unit conversions are necessary, this function should be used to obtain the conversion factors, as this will ensure that all drivers and engines use conversion factors that are self-consistent across codes.
 * Use of conversion factors that are not self-consistent can result in numerical instabilities.
 *
 * The following is a list of the unit categories, along with the names of the units associated with each category:
 *
 * - charge
 *    - atomic_unit_of_charge
 *    - coulomb
 * - energy
 *    - atomic_unit_of_energy
 *    - calorie
 *    - electron_volt
 *    - hartree
 *    - inverse_meter_energy
 *    - joule
 *    - kelvin_energy
 *    - kilocalorie
 *    - kilocalorie_per_mol
 *    - kilojoule
 *    - kilojoule_per_mol
 *    - rydberg
 * - force
 *    - atomic_unit_of_force
 *    - newton
 * - length
 *    - angstrom
 *    - atomic_unit_of_length
 *    - bohr
 *    - meter
 *    - nanometer
 *    - picometer
 * - mass
 *    - atomic_mass_unit
 *    - atomic_unit_of_mass
 *    - gram
 *    - kilogram
 * - time
 *    - atomic_unit_of_time
 *    - picosecond
 *    - second
 *
 * All conversion factors were acquired from the <a href="https://physics.nist.gov/cuu/Constants/Table/allascii.txt">NIST CODATA</a>, 
 * except the conversion factors for calorie, kilocalorie, and kilocalorie_per_mol, which were obtained from the <a href="https://www.nist.gov/pml/nist-guide-si-appendix-b9-factors-units-listed-kind-quantity-or-field-science"> NIST Guide to the SI</a>.
 *
 * \param [in]       in_unit
 *                   Name of the unit to convert from.
 * \param [in]       out_unit
 *                   Name of the unit to convert to.
 * \param [out]      conv
 *                   Conversion factor from in_unit to out_unit
 */
int MDI_Conversion_Factor(const char* in_unit, const char* out_unit, double* conv)
{
  // Except where otherwise noted, all values are from:
  //   - https://physics.nist.gov/cuu/Constants/Table/allascii.txt

  // mass
  double atomic_unit_of_mass = 1.0;
  double kilogram = 1.09776910575776e30;
  double gram = 1.09776910575776e27;
  double atomic_mass_unit = 1822.88848621731;

  // charge
  double atomic_unit_of_charge = 1.0;
  double coulomb = 6.24150907446076e18;

  // energy
  double atomic_unit_of_energy = 1.0;
  double hartree = 1.0;
  double joule = 2.29371227839632e17;
  double kilojoule = 2.29371227839632e20;
  double kilojoule_per_mol = 3.80879884713342e-4;
  double calorie = 9.5968921728102e17;
     // from https://www.nist.gov/pml/nist-guide-si-appendix-b9-factors-units-listed-kind-quantity-or-field-science
  double kilocalorie = 9.5968921728102e20;
     // from https://www.nist.gov/pml/nist-guide-si-appendix-b9-factors-units-listed-kind-quantity-or-field-science
  double kilocalorie_per_mol = 1.59360143764062e-3;
     // from https://www.nist.gov/pml/nist-guide-si-appendix-b9-factors-units-listed-kind-quantity-or-field-science
  double electron_volt = 3.67493221756549e-2;
  double rydberg = 0.5;
  double kelvin_energy = 3.16681156345558e-6;
  double inverse_meter_energy = 4.55633525291194e-8;

  // force
  double atomic_unit_of_force = 1.0;
  double newton = 1.21378026608897e7;

  // length
  double atomic_unit_of_length = 1.0;
  double bohr = 1.0;
  double meter = 1.88972612462577e10;
  double nanometer = 1.88972612462577e1;
  double picometer = 1.88972612462577e-2;
  double angstrom = 1.88975437603133;

  // time
  double atomic_unit_of_time = 1.0;
  double second = 4.13413733351821e16;
  double picosecond = 4.13413733351821e4;

  // categories of input and output units (i.e. length, force, time, etc.)
  //    0 - none
  //    1 - mass
  //    2 - charge
  //    3 - energy
  //    4 - force
  //    5 - length
  //    6 - time
  int in_category = 0;
  int out_category = 0;

  // conversion factors for input and output units
  double in_conv = 1.0;
  double out_conv = 1.0;

  // identify the input unit
  if ( strcmp( in_unit, "atomic_unit_of_mass" ) == 0 ) {
    in_category = 1;
    in_conv = atomic_unit_of_mass;
  }
  else if ( strcmp( in_unit, "kilogram" ) == 0 ) {
    in_category = 1;
    in_conv = kilogram;
  }
  else if ( strcmp( in_unit, "gram" ) == 0 ) {
    in_category = 1;
    in_conv = gram;
  }
  else if ( strcmp( in_unit, "atomic_mass_unit" ) == 0 ) {
    in_category = 1;
    in_conv = atomic_mass_unit;
  }
  else if ( strcmp( in_unit, "atomic_unit_of_charge" ) == 0 ) {
    in_category = 2;
    in_conv = atomic_unit_of_charge;
  }
  else if ( strcmp( in_unit, "coulomb" ) == 0 ) {
    in_category = 2;
    in_conv = coulomb;
  }
  else if ( strcmp( in_unit, "atomic_unit_of_energy" ) == 0 ) {
    in_category = 3;
    in_conv = atomic_unit_of_energy;
  }
  else if ( strcmp( in_unit, "hartree" ) == 0 ) {
    in_category = 3;
    in_conv = hartree;
  }
  else if ( strcmp( in_unit, "joule" ) == 0 ) {
    in_category = 3;
    in_conv = joule;
  }
  else if ( strcmp( in_unit, "kilojoule" ) == 0 ) {
    in_category = 3;
    in_conv = kilojoule;
  }
  else if ( strcmp( in_unit, "kilojoule_per_mol" ) == 0 ) {
    in_category = 3;
    in_conv = kilojoule_per_mol;
  }
  else if ( strcmp( in_unit, "calorie" ) == 0 ) {
    in_category = 3;
    in_conv = calorie;
  }
  else if ( strcmp( in_unit, "kilocalorie" ) == 0 ) {
    in_category = 3;
    in_conv = kilocalorie;
  }
  else if ( strcmp( in_unit, "kilocalorie_per_mol" ) == 0 ) {
    in_category = 3;
    in_conv = kilocalorie_per_mol;
  }
  else if ( strcmp( in_unit, "electron_volt" ) == 0 ) {
    in_category = 3;
    in_conv = electron_volt;
  }
  else if ( strcmp( in_unit, "rydberg" ) == 0 ) {
    in_category = 3;
    in_conv = rydberg;
  }
  else if ( strcmp( in_unit, "kelvin_energy" ) == 0 ) {
    in_category = 3;
    in_conv = kelvin_energy;
  }
  else if ( strcmp( in_unit, "inverse_meter_energy" ) == 0 ) {
    in_category = 3;
    in_conv = inverse_meter_energy;
  }
  else if ( strcmp( in_unit, "atomic_unit_of_force" ) == 0 ) {
    in_category = 4;
    in_conv = inverse_meter_energy;
  }
  else if ( strcmp( in_unit, "newton" ) == 0 ) {
    in_category = 4;
    in_conv = newton;
  }
  else if ( strcmp( in_unit, "atomic_unit_of_length" ) == 0 ) {
    in_category = 5;
    in_conv = atomic_unit_of_length;
  }
  else if ( strcmp( in_unit, "bohr" ) == 0 ) {
    in_category = 5;
    in_conv = bohr;
  }
  else if ( strcmp( in_unit, "meter" ) == 0 ) {
    in_category = 5;
    in_conv = meter;
  }
  else if ( strcmp( in_unit, "nanometer" ) == 0 ) {
    in_category = 5;
    in_conv = nanometer;
  }
  else if ( strcmp( in_unit, "picometer" ) == 0 ) {
    in_category = 5;
    in_conv = picometer;
  }
  else if ( strcmp( in_unit, "angstrom" ) == 0 ) {
    in_category = 5;
    in_conv = angstrom;
  }
  else if ( strcmp( in_unit, "atomic_unit_of_time" ) == 0 ) {
    in_category = 6;
    in_conv = atomic_unit_of_time;
  }
  else if ( strcmp( in_unit, "second" ) == 0 ) {
    in_category = 6;
    in_conv = second;
  }
  else if ( strcmp( in_unit, "picosecond" ) == 0 ) {
    in_category = 6;
    in_conv = picosecond;
  }
  else {
    mdi_error("Unit name not recognized");
  }

  // identify the output unit
  if ( strcmp( out_unit, "atomic_unit_of_mass" ) == 0 ) {
    out_category = 1;
    out_conv = atomic_unit_of_mass;
  }
  else if ( strcmp( out_unit, "kilogram" ) == 0 ) {
    out_category = 1;
    out_conv = kilogram;
  }
  else if ( strcmp( out_unit, "gram" ) == 0 ) {
    out_category = 1;
    out_conv = gram;
  }
  else if ( strcmp( out_unit, "atomic_mass_unit" ) == 0 ) {
    out_category = 1;
    out_conv = atomic_mass_unit;
  }
  else if ( strcmp( out_unit, "atomic_unit_of_charge" ) == 0 ) {
    out_category = 2;
    out_conv = atomic_unit_of_charge;
  }
  else if ( strcmp( out_unit, "coulomb" ) == 0 ) {
    out_category = 2;
    out_conv = coulomb;
  }
  else if ( strcmp( out_unit, "atomic_unit_of_energy" ) == 0 ) {
    out_category = 3;
    out_conv = atomic_unit_of_energy;
  }
  else if ( strcmp( out_unit, "hartree" ) == 0 ) {
    out_category = 3;
    out_conv = hartree;
  }
  else if ( strcmp( out_unit, "joule" ) == 0 ) {
    out_category = 3;
    out_conv = joule;
  }
  else if ( strcmp( out_unit, "kilojoule" ) == 0 ) {
    out_category = 3;
    out_conv = kilojoule;
  }
  else if ( strcmp( out_unit, "kilojoule_per_mol" ) == 0 ) {
    out_category = 3;
    out_conv = kilojoule_per_mol;
  }
  else if ( strcmp( out_unit, "calorie" ) == 0 ) {
    out_category = 3;
    out_conv = calorie;
  }
  else if ( strcmp( out_unit, "kilocalorie" ) == 0 ) {
    out_category = 3;
    out_conv = kilocalorie;
  }
  else if ( strcmp( out_unit, "kilocalorie_per_mol" ) == 0 ) {
    out_category = 3;
    out_conv = kilocalorie_per_mol;
  }
  else if ( strcmp( out_unit, "electron_volt" ) == 0 ) {
    out_category = 3;
    out_conv = electron_volt;
  }
  else if ( strcmp( out_unit, "rydberg" ) == 0 ) {
    out_category = 3;
    out_conv = rydberg;
  }
  else if ( strcmp( out_unit, "kelvin_energy" ) == 0 ) {
    out_category = 3;
    out_conv = kelvin_energy;
  }
  else if ( strcmp( out_unit, "inverse_meter_energy" ) == 0 ) {
    out_category = 3;
    out_conv = inverse_meter_energy;
  }
  else if ( strcmp( out_unit, "atomic_unit_of_force" ) == 0 ) {
    out_category = 4;
    out_conv = inverse_meter_energy;
  }
  else if ( strcmp( out_unit, "newton" ) == 0 ) {
    out_category = 4;
    out_conv = newton;
  }
  else if ( strcmp( out_unit, "atomic_unit_of_length" ) == 0 ) {
    out_category = 5;
    out_conv = atomic_unit_of_length;
  }
  else if ( strcmp( out_unit, "bohr" ) == 0 ) {
    out_category = 5;
    out_conv = bohr;
  }
  else if ( strcmp( out_unit, "meter" ) == 0 ) {
    out_category = 5;
    out_conv = meter;
  }
  else if ( strcmp( out_unit, "nanometer" ) == 0 ) {
    out_category = 5;
    out_conv = nanometer;
  }
  else if ( strcmp( out_unit, "picometer" ) == 0 ) {
    out_category = 5;
    out_conv = picometer;
  }
  else if ( strcmp( out_unit, "angstrom" ) == 0 ) {
    out_category = 5;
    out_conv = angstrom;
  }
  else if ( strcmp( out_unit, "atomic_unit_of_time" ) == 0 ) {
    out_category = 6;
    out_conv = atomic_unit_of_time;
  }
  else if ( strcmp( out_unit, "second" ) == 0 ) {
    out_category = 6;
    out_conv = second;
  }
  else if ( strcmp( out_unit, "picosecond" ) == 0 ) {
    out_category = 6;
    out_conv = picosecond;
  }
  else {
    mdi_error("Unit name not recognized");
  }

  // confirm that the two units belong to the same category (i.e. mass, energy, time, etc.)
  if ( in_category != out_category ) {
    mdi_error("The units are of two different types, and conversion is not possible.");
  }

  *conv = in_conv / out_conv;
  return 0;
}


/*! \brief Return order of this code within all codes represented in MPI_COMM_WORLD
 *
 * When using the MPI communication method, all processes across all codes are spawned 
 * as part of the same MPI_COMM_WORLD.
 * This funciton returns the order of the code associated with the calling process 
 * within MPI_COMM_WORLD.
 *
 */
int MDI_Get_MPI_Code_Rank()
{
  return mpi_code_rank;
}

/*! \brief Return the rank of the calling process within its associated code
 *
 * When using the MPI communication method, all processes across all codes are spawned 
 * as part of the same MPI_COMM_WORLD.
 * This funciton returns the rank of the calling process within the subset of processes
 * associated with the same code.
 *
 */
void MDI_Set_MPI_Intra_Rank(int rank)
{
  intra_rank = rank;
}

/*! \brief Set the size of MPI_COMM_WORLD
 *
 * This function is only used if the linked program uses MPI4PY.
 *
 */
void MDI_Set_World_Size(int world_size_in)
{
  set_world_size(world_size_in);
}

/*! \brief Set the rank of this process within MPI_COMM_WORLD
 *
 * This function is only used if the linked program uses MPI4PY.
 *
 */
void MDI_Set_World_Rank(int world_rank_in)
{
  set_world_rank(world_rank_in);
}

/*! \brief Register a node
 *
 * The function returns \p 0 on a success.
 *
 * \param [in]       node_name
 *                   Name of the node.
 */
int MDI_Register_Node(const char* node_name)
{
  if ( is_initialized == 0 ) {
    mdi_error("MDI_Register_Node called but MDI has not been initialized");
  }

  // confirm that the node_name size is not greater than MDI_COMMAND_LENGTH
  if ( strlen(node_name) > COMMAND_LENGTH ) {
    mdi_error("Cannot register name with length greater than MDI_COMMAND_LENGTH");
  }

  // confirm that this node is not already registered
  int node_index = get_node_index(&nodes, node_name);
  if ( node_index != -1 ) {
    mdi_error("This node is already registered");
  }

  node new_node;
  vector* command_vec = malloc(sizeof(vector));
  vector* callback_vec = malloc(sizeof(vector));
  vector_init(command_vec, sizeof(char[COMMAND_LENGTH]));
  vector_init(callback_vec, sizeof(char[COMMAND_LENGTH]));
  new_node.commands = command_vec;
  new_node.callbacks = callback_vec;
  strcpy(new_node.name, node_name);
  vector_push_back(&nodes, &new_node);
  return 0;
}

/*! \brief Check whether a node is supported on a specified engine
 *
 * The function returns \p 0 on a success.
 *
 * \param [in]       node_name
 *                   Name of the node.
 * \param [in]       comm
 *                   MDI communicator of the engine.  If comm is set to 
 *                   MDI_NULL_COMM, the function will check for the calling engine.
 * \param [out]      flag
 *                   On return, 1 if the node is supported and 0 otherwise.
 */
int MDI_Check_Node_Exists(const char* node_name, MDI_Comm comm, int* flag)
{
  if ( is_initialized == 0 ) {
    mdi_error("MDI_Check_Node_Exists called but MDI has not been initialized");
  }

  // confirm that the node_name size is not greater than MDI_COMMAND_LENGTH
  if ( strlen(node_name) > COMMAND_LENGTH ) {
    mdi_error("Node name is greater than MDI_COMMAND_LENGTH");
  }

  // find the node
  int node_index = get_node_index(&nodes, node_name);
  if ( node_index == -1 ) {
    *flag = 0;
  }
  else {
    *flag = 1;
  }
  return 0;
}

/*! \brief Get the number of nodes on a specified engine
 *
 * The function returns \p 0 on a success.
 *
 * \param [in]       comm
 *                   MDI communicator of the engine.  If comm is set to 
 *                   MDI_NULL_COMM, the function will check for the calling engine.
 * \param [out]      nnodes
 *                   On return, the number of nodes supported by the engine.
 */
int MDI_Get_NNodes(MDI_Comm comm, int* nnodes)
{
  if ( is_initialized == 0 ) {
    mdi_error("MDI_Get_NNodes called but MDI has not been initialized");
  }

  *nnodes = nodes.size;
  return 0;
}

/*! \brief Get the name of a node on a specified engine
 *
 * The function returns \p 0 on a success.
 *
 * \param [in]       index
 *                   Index of the node on the specified engine.
 * \param [in]       comm
 *                   MDI communicator of the engine.  If comm is set to 
 *                   MDI_NULL_COMM, the function will check for the calling engine.
 * \param [out]      name
 *                   On return, the name of the node
 */
int MDI_Get_Node(int index, MDI_Comm comm, char* name)
{
  if ( is_initialized == 0 ) {
    mdi_error("MDI_Get_Node called but MDI has not been initialized");
  }

  node* ret_node = vector_get(&nodes, index);
  strcpy(name, &ret_node->name[0]);
  return 0;
}

/*! \brief Register a command on a specified node
 *
 * The function returns \p 0 on a success.
 *
 * \param [in]       node_name
 *                   Name of the node on which the command will be registered.
 * \param [in]       command_name
 *                   Name of the command.
 */
int MDI_Register_Command(const char* node_name, const char* command_name)
{
  if ( is_initialized == 0 ) {
    mdi_error("MDI_Register_Command called but MDI has not been initialized");
  }

  // confirm that the node_name size is not greater than MDI_COMMAND_LENGTH
  if ( strlen(node_name) > COMMAND_LENGTH ) {
    mdi_error("Node name is greater than MDI_COMMAND_LENGTH");
  }

  // confirm that the command_name size is not greater than MDI_COMMAND_LENGTH
  if ( strlen(command_name) > COMMAND_LENGTH ) {
    mdi_error("Cannot register name with length greater than MDI_COMMAND_LENGTH");
  }

  // find the node
  int node_index = get_node_index(&nodes, node_name);
  if ( node_index == -1 ) {
    mdi_error("Attempting to register a command on an unregistered node");
  }
  node* target_node = vector_get(&nodes, node_index);

  // confirm that this command is not already registered
  int command_index = get_command_index(target_node, command_name);
  if ( command_index != -1 ) {
    mdi_error("This command is already registered for this node");
  }

  // register this command
  char new_command[COMMAND_LENGTH];
  strcpy(new_command, command_name);
  vector_push_back( target_node->commands, &new_command );

  return 0;
}

/*! \brief Check whether a command is supported on specified node on a specified engine
 *
 * The function returns \p 0 on a success.
 *
 * \param [in]       node_name
 *                   Name of the command's node.
 * \param [in]       command_name
 *                   Name of the command.
 * \param [in]       comm
 *                   MDI communicator of the engine.  If comm is set to 
 *                   MDI_NULL_COMM, the function will check for the calling engine.
 * \param [out]      flag
 *                   On return, 1 if the command is supported and 0 otherwise.
 */
int MDI_Check_Command_Exists(const char* node_name, const char* command_name, MDI_Comm comm, int* flag)
{
  if ( is_initialized == 0 ) {
    mdi_error("MDI_Check_Command_Exists called but MDI has not been initialized");
  }

  // confirm that the node_name size is not greater than MDI_COMMAND_LENGTH
  if ( strlen(node_name) > COMMAND_LENGTH ) {
    mdi_error("Node name is greater than MDI_COMMAND_LENGTH");
  }

  // confirm that the command_name size is not greater than MDI_COMMAND_LENGTH
  if ( strlen(command_name) > COMMAND_LENGTH ) {
    mdi_error("Cannot register name with length greater than MDI_COMMAND_LENGTH");
  }

  // find the node
  int node_index = get_node_index(&nodes, node_name);
  if ( node_index == -1 ) {
    mdi_error("Could not find the node");
  }
  node* target_node = vector_get(&nodes, node_index);

  // find the command
  int command_index = get_command_index(target_node, command_name);
  if ( command_index == -1 ) {
    *flag = 0;
  }
  else {
    *flag = 1;
  }
  return 0;
}

/*! \brief Get the number of commands supported for a specified node on a specified engine
 *
 * The function returns \p 0 on a success.
 *
 * \param [in]       node_name
 *                   Name of the node.
 * \param [in]       comm
 *                   MDI communicator of the engine.  If comm is set to 
 *                   MDI_NULL_COMM, the function will check for the calling engine.
 * \param [out]      nnodes
 *                   On return, the number of commands supported on the specified engine
 *                   on the specified node.
 */
int MDI_Get_NCommands(const char* node_name, MDI_Comm comm, int* ncommands)
{
  if ( is_initialized == 0 ) {
    mdi_error("MDI_Get_NCommands called but MDI has not been initialized");
  }

  // confirm that the node_name size is not greater than MDI_COMMAND_LENGTH
  if ( strlen(node_name) > COMMAND_LENGTH ) {
    mdi_error("Node name is greater than MDI_COMMAND_LENGTH");
  }

  // find the node
  int node_index = get_node_index(&nodes, node_name);
  if ( node_index == -1 ) {
    mdi_error("Could not find the node");
  }
  node* target_node = vector_get(&nodes, node_index);

  *ncommands = target_node->commands->size;
  return 0;
}

/*! \brief Get the name of a command on a specified node on a specified engine
 *
 * The function returns \p 0 on a success.
 *
 * \param [in]       node_name
 *                   Name of the node on which the command is located.
 * \param [in]       index
 *                   Index of the command on the specified node.
 * \param [in]       comm
 *                   MDI communicator of the engine.  If comm is set to 
 *                   MDI_NULL_COMM, the function will check for the calling engine.
 * \param [out]      name
 *                   On return, the name of the command
 */
int MDI_Get_Command(const char* node_name, int index, MDI_Comm comm, char* name)
{
  if ( is_initialized == 0 ) {
    mdi_error("MDI_Get_Command called but MDI has not been initialized");
  }

  // find the node
  int node_index = get_node_index(&nodes, node_name);
  if ( node_index == -1 ) {
    mdi_error("MDI_Get_Command could not find the requested node");
  }
  node* target_node = vector_get(&nodes, node_index);

  if ( target_node->commands->size <= index ) {
    mdi_error("MDI_Get_Command failed because the command does not exist");
  }

  char* target_command = vector_get( target_node->commands, index );
  strcpy(name, target_command);
  return 0;
}

/*! \brief Register a callback on a specified node
 *
 * The function returns \p 0 on a success.
 *
 * \param [in]       node_name
 *                   Name of the node on which the callback will be registered.
 * \param [in]       callback_name
 *                   Name of the callback.
 */
int MDI_Register_Callback(const char* node_name, const char* callback_name)
{
  if ( is_initialized == 0 ) {
    mdi_error("MDI_Register_Callback called but MDI has not been initialized");
  }

  // confirm that the node_name size is not greater than MDI_COMMAND_LENGTH
  if ( strlen(node_name) > COMMAND_LENGTH ) {
    mdi_error("Node name is greater than MDI_COMMAND_LENGTH");
  }

  // confirm that the callback_name size is not greater than MDI_COMMAND_LENGTH
  if ( strlen(callback_name) > COMMAND_LENGTH ) {
    mdi_error("Cannot register name with length greater than MDI_COMMAND_LENGTH");
  }

  // find the node
  int node_index = get_node_index(&nodes, node_name);
  if ( node_index == -1 ) {
    mdi_error("Attempting to register a callback on an unregistered node");
  }
  node* target_node = vector_get(&nodes, node_index);

  // confirm that this callback is not already registered
  int callback_index = get_callback_index(target_node, callback_name);
  if ( callback_index != -1 ) {
    mdi_error("This callback is already registered for this node");
  }

  // register this callback
  char new_callback[COMMAND_LENGTH];
  strcpy(new_callback, callback_name);
  vector_push_back( target_node->callbacks, &new_callback );

  return 0;
}

/*! \brief Check whether a callback exists on specified node on a specified engine
 *
 * The function returns \p 0 on a success.
 *
 * \param [in]       node_name
 *                   Name of the callbacks's node.
 * \param [in]       command_name
 *                   Name of the callback.
 * \param [in]       comm
 *                   MDI communicator of the engine.  If comm is set to 
 *                   MDI_NULL_COMM, the function will check for the calling engine.
 * \param [out]      flag
 *                   On return, 1 if the callback is supported and 0 otherwise.
 */
int MDI_Check_Callback_Exists(const char* node_name, const char* callback_name, MDI_Comm comm, int* flag)
{
  if ( is_initialized == 0 ) {
    mdi_error("MDI_Check_Callback_Exists called but MDI has not been initialized");
  }

  // confirm that the node_name size is not greater than MDI_COMMAND_LENGTH
  if ( strlen(node_name) > COMMAND_LENGTH ) {
    mdi_error("Node name is greater than MDI_COMMAND_LENGTH");
  }

  // confirm that the callback_name size is not greater than MDI_COMMAND_LENGTH
  if ( strlen(callback_name) > COMMAND_LENGTH ) {
    mdi_error("Cannot register name with length greater than MDI_COMMAND_LENGTH");
  }

  // find the node
  int node_index = get_node_index(&nodes, node_name);
  if ( node_index == -1 ) {
    mdi_error("Could not find the node");
  }
  node* target_node = vector_get(&nodes, node_index);

  // find the callback
  int callback_index = get_callback_index(target_node, callback_name);
  if ( callback_index == -1 ) {
    *flag = 0;
  }
  else {
    *flag = 1;
  }
  return 0;
}

/*! \brief Get the number of callbacks on a specified node on a specified engine
 *
 * The function returns \p 0 on a success.
 *
 * \param [in]       node_name
 *                   Name of the node.
 * \param [in]       comm
 *                   MDI communicator of the engine.  If comm is set to 
 *                   MDI_NULL_COMM, the function will check for the calling engine.
 * \param [out]      ncallbacks
 *                   On return, the number of callbacks on the specified node
 *                   on the specified engine.
 */
int MDI_Get_NCallbacks(const char* node_name, MDI_Comm comm, int* ncallbacks)
{
  if ( is_initialized == 0 ) {
    mdi_error("MDI_Get_NCallbacks called but MDI has not been initialized");
  }

  // confirm that the node_name size is not greater than MDI_COMMAND_LENGTH
  if ( strlen(node_name) > COMMAND_LENGTH ) {
    mdi_error("Node name is greater than MDI_COMMAND_LENGTH");
  }

  // find the node
  int node_index = get_node_index(&nodes, node_name);
  if ( node_index == -1 ) {
    mdi_error("Could not find the node");
  }
  node* target_node = vector_get(&nodes, node_index);

  *ncallbacks = target_node->callbacks->size;
  return 0;
}

/*! \brief Get the name of a callback on a specified node on a specified engine
 *
 * The function returns \p 0 on a success.
 *
 * \param [in]       node_name
 *                   Name of the node on which the callback is located.
 * \param [in]       index
 *                   Index of the callback on the specified node.
 * \param [in]       comm
 *                   MDI communicator of the engine.  If comm is set to 
 *                   MDI_NULL_COMM, the function will check for the calling engine.
 * \param [out]      name
 *                   On return, the name of the callback
 */
int MDI_Get_Callback(const char* node_name, int index, MDI_Comm comm, char* name)
{
  if ( is_initialized == 0 ) {
    mdi_error("MDI_Get_Callback called but MDI has not been initialized");
  }

  // find the node
  int node_index = get_node_index(&nodes, node_name);
  if ( node_index == -1 ) {
    mdi_error("MDI_Get_Command could not find the requested node");
  }
  node* target_node = vector_get(&nodes, node_index);

  if ( target_node->callbacks->size <= index ) {
    mdi_error("MDI_Get_Command failed because the command does not exist");
  }

  char* target_callback = vector_get( target_node->callbacks, index );
  strcpy(name, target_callback);
  return 0;
}
