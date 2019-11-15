/*! \file
 *
 * \brief Generic MDI function calls
 */

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "mdi.h"
#include "mdi_general.h"
#include "mdi_mpi.h"
#include "mdi_tcp.h"
#include "mdi_test.h"

//this is the number of communicator handles that have been returned by MDI_Accept_Connection()
static int returned_comms = 0;

//name of the driver/engine
static char* name;

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
int general_init(const char* options, void* world_comm) {
  returned_comms = 0;
  vector_init(&communicators, sizeof(communicator));

  // initialize the vector of nodes
  vector_init(&nodes, sizeof(node));

  char* strtol_ptr;
  int i;

  int mpi_initialized = 0;

  // values acquired from the input options
  char* role;
  char* method;
  name = malloc(MDI_NAME_LENGTH+1);
  char* hostname;
  int port;
  char* output_file;
  char* language = ((char*)"");
  int has_role = 0;
  int has_method = 0;
  int has_name = 0;
  int has_hostname = 0;
  int has_port = 0;
  int has_output_file = 0;

  // calculate argc
  char* argv_line = strdup(options);
  char* token = strtok(argv_line, " ");
  int argc = 0;
  while (token != NULL) {
    argc++;
    token = strtok(NULL," ");
  }

  // calculate argv
  char** argv = malloc( argc * sizeof(char*) );
  argv_line = strdup(options);
  token = strtok(argv_line, " ");
  for (i=0; i<argc; i++) {
    argv[i] = token;
    token = strtok(NULL," ");
  }

  // read options
  int iarg = 0;
  while (iarg < argc) {

    //-role
    if (strcmp(argv[iarg],"-role") == 0){
      if (iarg+2 > argc) {
	mdi_error("Argument missing from -role option");
      }
      role = argv[iarg+1];
      has_role = 1;
      iarg += 2;
    }
    //-method
    else if (strcmp(argv[iarg],"-method") == 0) {
      if (iarg+2 > argc) {
	mdi_error("Argument missing from -method option");
      }
      method = argv[iarg+1];
      has_method = 1;
      iarg += 2;
    }
    //-name
    else if (strcmp(argv[iarg],"-name") == 0){
      if (iarg+2 > argc) {
	mdi_error("Argument missing from -name option");
      }
      if ( strlen(argv[iarg+1]) > MDI_NAME_LENGTH ) {
	mdi_error("Name argument length exceeds MDI_NAME_LENGTH");
      }
      strcpy(name, argv[iarg+1]);
      has_name = 1;
      iarg += 2;
    }
    //-hostname
    else if (strcmp(argv[iarg],"-hostname") == 0){
      if (iarg+2 > argc) {
	mdi_error("Argument missing from -hostname option");
      }
      hostname = argv[iarg+1];
      has_hostname = 1;
      iarg += 2;
    }
    //-port
    else if (strcmp(argv[iarg],"-port") == 0) {
      if (iarg+2 > argc) {
	mdi_error("Argument missing from -port option");
      }
      port = strtol( argv[iarg+1], &strtol_ptr, 10 );
      has_port = 1;
      iarg += 2;
    }
    //-ipi
    else if (strcmp(argv[iarg],"-ipi") == 0) {
      ipi_compatibility = 1;
      iarg += 1;
    }
    //-out
    else if (strcmp(argv[iarg],"-out") == 0) {
      if (iarg+2 > argc) {
	mdi_error("Argument missing from -out option");
      }
      has_output_file = 1;
      output_file = argv[iarg+1];
      iarg += 2;
    }
    //_language
    else if (strcmp(argv[iarg],"_language") == 0) {
      if (iarg+2 > argc) {
	mdi_error("Argument missing from _language option");
      }
      language = argv[iarg+1];
      iarg += 2;
    }
    else {
      mdi_error("Unrecognized option");
    }
  }

  // get the MPI rank
  MPI_Comm mpi_communicator;
  int mpi_rank = 0;
  if ( world_comm == NULL ) {
    mpi_communicator = 0;
    mpi_rank = 0;
  }
  else {
    if ( world_rank == -1 ) {
      if ( strcmp(language, "Fortran") == 0 ) {
	mpi_communicator = MPI_Comm_f2c( *(MPI_Fint*) world_comm );
      }
      else {
	mpi_communicator = *(MPI_Comm*) world_comm;
      }
      MPI_Comm_rank(mpi_communicator, &mpi_rank);
    }
    else {
      mpi_rank = 0;
    }
  }

  // redirect the standard output
  if ( has_output_file == 1 ) {
    freopen(output_file, "w", stdout);
  }

  // ensure the -role option was provided
  if ( has_role == 0 ) {
    mdi_error("Error in MDI_Init: -role option not provided");
  }

  // ensure the -name option was provided
  if ( has_name == 0 ) {
    mdi_error("Error in MDI_Init: -name option not provided");
  }

  // ensure the -method option was provided
  if ( has_method == 0 ) {
    mdi_error("Error in MDI_Init: -method option not provided");
  }

  // determine whether the intra-code MPI communicator should be split by mpi_init_mdi
  int do_split = 1;
  if ( strcmp(language, "Python") == 0 ) {
    do_split = 0;
  }

  if ( strcmp(role, "DRIVER") == 0 ) {
    // initialize this code as a driver

    if ( strcmp(method, "MPI") == 0 ) {
      mpi_identify_codes("", do_split, mpi_communicator);
      mpi_initialized = 1;
    }
    else if ( strcmp(method, "TCP") == 0 ) {
      if ( has_port == 0 ) {
	mdi_error("Error in MDI_Init: -port option not provided");
      }
      if ( mpi_rank == 0 ) {
	tcp_listen(port);
      }
    }
    else if ( strcmp(method, "TEST") == 0 ) {
      test_initialize();
    }
    else {
      mdi_error("Error in MDI_Init: method not recognized");
    }

  }
  else if ( strcmp(role,"ENGINE") == 0 ) {
    // initialize this code as an engine

    if ( strcmp(method, "MPI") == 0 ) {
      mpi_identify_codes(name, do_split, mpi_communicator);
      mpi_initialized = 1;
    }
    else if ( strcmp(method, "TCP") == 0 ) {
      if ( has_hostname == 0 ) {
	mdi_error("Error in MDI_Init: -hostname option not provided");
      }
      if ( has_port == 0 ) {
	mdi_error("Error in MDI_Init: -port option not provided");
      }
      if ( mpi_rank == 0 ) {
	tcp_request_connection(port, hostname);
      }
    }
    else if ( strcmp(method, "TEST") == 0 ) {
      test_initialize();
    }
    else {
      mdi_error("Error in MDI_Init: method not recognized");
    }

    
  }
  else {
    mdi_error("Error in MDI_Init: role not recognized");
  }

  // set the MPI communicator correctly
  if ( mpi_initialized == 1 ) {
    if ( do_split == 1 ) {
      if ( strcmp(language, "Fortran") == 0 ) {
	mpi_communicator = MPI_Comm_f2c( *(MPI_Fint*) world_comm );
	mpi_update_world_comm( (void*) &mpi_communicator);
        MPI_Fint f_comm = MPI_Comm_c2f( mpi_communicator );
	world_comm = (void*) &f_comm;
      }
      else {
	mpi_update_world_comm(world_comm);
      }
    }
  }

  free( argv_line );
  free( argv );

  return 0;
}


/*! \brief Accept a new MDI communicator
 *
 * The function returns an MDI_Comm that describes a connection between two codes.
 * If no new communicators are available, the function returns \p MDI_NULL_COMM.
 *
 */
int general_accept_communicator() {
  // if MDI hasn't returned some connections, do that now
  if ( returned_comms < communicators.size ) {
    returned_comms++;
    return returned_comms;
  }

  // check for any production codes connecting via TCP
  if ( tcp_socket > 0 ) {

    //accept a connection via TCP
    tcp_accept_connection();

    // if MDI hasn't returned some connections, do that now
    if ( returned_comms < communicators.size ) {
      returned_comms++;
      return (MDI_Comm)returned_comms;
    }

  }

  // unable to accept any connections
  return MDI_NULL_COMM;
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
int general_send(const void* buf, int count, MDI_Datatype datatype, MDI_Comm comm) {
  if ( intra_rank != 0 ) {
    mdi_error("Called MDI_Send with incorrect rank");
  }

  //communicator_send(buf, count, datatype, comm);

  communicator* this = vector_get(&communicators, comm-1);

  if ( this->method == MDI_MPI ) {
    mpi_send(buf, count, datatype, comm);
  }
  else if ( this->method == MDI_TCP ) {
    tcp_send(buf, count, datatype, comm);
  }
  else if ( this->method == MDI_TEST ) {
    test_send(buf, count, datatype, comm);
  }
  else {
    mdi_error("MDI method not recognized in communicator_send");
  }


  return 0;
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
int general_recv(void* buf, int count, MDI_Datatype datatype, MDI_Comm comm) {
  if ( intra_rank != 0 ) {
    mdi_error("Called MDI_Recv with incorrect rank");
  }

  //communicator_recv(buf, count, datatype, comm);

  communicator* this = vector_get(&communicators, comm-1);

  if ( this->method == MDI_MPI ) {
    mpi_recv(buf, count, datatype, comm);
  }
  else if ( this->method == MDI_TCP ) {
    tcp_recv(buf, count, datatype, comm);
  }
  else if ( this->method == MDI_TEST ) {
    test_recv(buf, count, datatype, comm);
  }
  else {
    mdi_error("MDI method not recognized in communicator_send");
  }


  return 0;
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
int general_send_command(const char* buf, MDI_Comm comm) {
  if ( intra_rank != 0 ) {
    mdi_error("Called MDI_Send_Command with incorrect rank");
  }
  int count = MDI_COMMAND_LENGTH;
  char* command = malloc( MDI_COMMAND_LENGTH * sizeof(char) );

  strncpy(command, buf, MDI_COMMAND_LENGTH);
  int ret = general_send( command, count, MDI_CHAR, comm );
  free( command );
  return ret;
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
int general_recv_command(char* buf, MDI_Comm comm) {
  if ( intra_rank != 0 ) {
    mdi_error("Called MDI_Recv_Command with incorrect rank");
  }
  int count = MDI_COMMAND_LENGTH;
  int datatype = MDI_CHAR;

  int ret = general_recv( buf, count, datatype, comm );

  // check if this command corresponds to one of MDI's standard built-in commands
  if ( strcmp( buf, "<NAME" ) == 0 ) {
    MDI_Send(name, MDI_NAME_LENGTH, MDI_CHAR, comm);
    return general_recv_command(buf, comm);
  }
  else if ( strcmp( buf, "<VERSION" ) == 0 ) {
    MDI_Send(&MDI_VERSION, 1, MDI_DOUBLE, comm);
    return general_recv_command(buf, comm);
  }
  else if ( strcmp( buf, "<COMMANDS" ) == 0 ) {
    send_command_list(comm);
    return general_recv_command(buf, comm);
  }
  else if ( strcmp( buf, "<CALLBACKS" ) == 0 ) {
    send_callback_list(comm);
    return general_recv_command(buf, comm);
  }
  else if ( strcmp( buf, "<NODES" ) == 0 ) {
    send_node_list(comm);
    return general_recv_command(buf, comm);
  }
  else if ( strcmp( buf, "<NCOMMANDS" ) == 0 ) {
    send_ncommands(comm);
    return general_recv_command(buf, comm);
  }
  else if ( strcmp( buf, "<NCALLBACKS" ) == 0 ) {
    send_ncallbacks(comm);
    return general_recv_command(buf, comm);
  }
  else if ( strcmp( buf, "<NNODES" ) == 0 ) {
    send_nnodes(comm);
    return general_recv_command(buf, comm);
  }

  return ret;
}


/*! \brief Register a node
 *
 * The function returns \p 0 on a success.
 *
 * \param [in]       node_vec
 *                   Vector of nodes, into which the new node will be added.
 * \param [in]       node_name
 *                   Name of the node.
 */
int register_node(vector* node_vec, const char* node_name)
{
  // confirm that the node_name size is not greater than MDI_COMMAND_LENGTH
  if ( strlen(node_name) > COMMAND_LENGTH ) {
    mdi_error("Cannot register name with length greater than MDI_COMMAND_LENGTH");
  }

  // confirm that this node is not already registered
  int node_index = get_node_index(node_vec, node_name);
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
  vector_push_back(node_vec, &new_node);
  return 0;
}


/*! \brief Register a command on a specified node
 *
 * The function returns \p 0 on a success.
 *
 * \param [in]       node_vec
 *                   Vector of nodes, into which the new node will be added.
 * \param [in]       node_name
 *                   Name of the node on which the command will be registered.
 * \param [in]       command_name
 *                   Name of the command.
 */
int register_command(vector* node_vec, const char* node_name, const char* command_name)
{
  // confirm that the node_name size is not greater than MDI_COMMAND_LENGTH
  if ( strlen(node_name) > COMMAND_LENGTH ) {
    mdi_error("Node name is greater than MDI_COMMAND_LENGTH");
  }

  // confirm that the command_name size is not greater than MDI_COMMAND_LENGTH
  if ( strlen(command_name) > COMMAND_LENGTH ) {
    mdi_error("Cannot register name with length greater than MDI_COMMAND_LENGTH");
  }

  // find the node
  int node_index = get_node_index(node_vec, node_name);
  if ( node_index == -1 ) {
    mdi_error("Attempting to register a command on an unregistered node");
  }
  node* target_node = vector_get(node_vec, node_index);

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


/*! \brief Register a callback on a specified node
 *
 * The function returns \p 0 on a success.
 *
 * \param [in]       node_vec
 *                   Vector of nodes, into which the new node will be added.
 * \param [in]       node_name
 *                   Name of the node on which the callback will be registered.
 * \param [in]       callback_name
 *                   Name of the callback.
 */
int register_callback(vector* node_vec, const char* node_name, const char* callback_name)
{
  // confirm that the node_name size is not greater than MDI_COMMAND_LENGTH
  if ( strlen(node_name) > COMMAND_LENGTH ) {
    mdi_error("Node name is greater than MDI_COMMAND_LENGTH");
  }

  // confirm that the callback_name size is not greater than MDI_COMMAND_LENGTH
  if ( strlen(callback_name) > COMMAND_LENGTH ) {
    mdi_error("Cannot register name with length greater than MDI_COMMAND_LENGTH");
  }

  // find the node
  int node_index = get_node_index(node_vec, node_name);
  if ( node_index == -1 ) {
    mdi_error("Attempting to register a callback on an unregistered node");
  }
  node* target_node = vector_get(node_vec, node_index);

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


/*! \brief Send the list of supported commands
 *
 * If running with MPI, this function must be called only by rank \p 0.
 * The function returns \p 0 on a success.
 *
 * \param [in]       comm
 *                   MDI communicator associated with the intended recipient code.
 */
int send_command_list(MDI_Comm comm) {
  if ( intra_rank != 0 ) {
    mdi_error("Attempting to send command information from the incorrect rank");
  }
  int ncommands = 0;
  int nnodes = nodes.size;
  int inode, icommand;
  int stride = MDI_COMMAND_LENGTH + 1;

  // determine the number of commands
  for (inode = 0; inode < nnodes; inode++) {
    node* this_node = vector_get(&nodes, inode);
    ncommands += this_node->commands->size;
  }

  // allocate memory for the commands list
  int count = ( ncommands + nnodes ) * stride;
  char* commands = malloc( count * sizeof(char) );

  // form the list of commands
  int islot = 0;
  for (inode = 0; inode < nnodes; inode++) {
    // add the name of this node to the list
    node* this_node = vector_get(&nodes, inode);
    int length = strlen(this_node->name);
    strcpy(&commands[ islot * stride ], this_node->name);
    int ichar;
    for (ichar = length; ichar < stride-1; ichar++) {
      strcpy( &commands[ islot * stride + ichar ], " " );
    }
    if ( this_node->commands->size > 0 ) {
      strcpy( &commands[ islot * stride + stride - 1 ], "," );
    }
    else {
      strcpy( &commands[ islot * stride + stride - 1 ], ";" );
    }
    islot++;

    // add the commands for this node
    for (icommand = 0; icommand < this_node->commands->size; icommand++) {
      char* command = vector_get(this_node->commands, icommand);
      length = strlen(command);
      strcpy(&commands[ islot * stride ], command);
      for (ichar = length; ichar < stride-1; ichar++) {
	strcpy( &commands[ islot * stride + ichar ], " " );
      }
      if ( icommand == this_node->commands->size - 1 ) {
	strcpy( &commands[ islot * stride + stride - 1], ";");
      }
      else {
	strcpy( &commands[ islot * stride + stride - 1], ",");
      }
      islot++;
    }
  }
  //printf("---COMMANDS: %d %s\n",count, commands);

  int ret = general_send( commands, count, MDI_CHAR, comm );
  free( commands );
  return ret;
}


/*! \brief Send the list of callbacks
 *
 * If running with MPI, this function must be called only by rank \p 0.
 * The function returns \p 0 on a success.
 *
 * \param [in]       comm
 *                   MDI communicator associated with the intended recipient code.
 */
int send_callback_list(MDI_Comm comm) {
  if ( intra_rank != 0 ) {
    mdi_error("Attempting to send callback information from the incorrect rank");
  }
  int ncallbacks = 0;
  int nnodes = nodes.size;
  int inode, icallback;
  int stride = MDI_COMMAND_LENGTH + 1;

  // determine the number of callbakcs
  for (inode = 0; inode < nnodes; inode++) {
    node* this_node = vector_get(&nodes, inode);
    ncallbacks += this_node->callbacks->size;
  }

  // allocate memory for the callbacks list
  int count = ( ncallbacks + nnodes ) * stride;
  char* callbacks = malloc( count * sizeof(char) );

  // form the list of callbacks
  int islot = 0;
  for (inode = 0; inode < nnodes; inode++) {
    // add the name of this node to the list
    node* this_node = vector_get(&nodes, inode);
    int length = strlen(this_node->name);
    strcpy(&callbacks[ islot * stride ], this_node->name);
    int ichar;
    for (ichar = length; ichar < stride-1; ichar++) {
      strcpy( &callbacks[ islot * stride + ichar ], " " );
    }
    if ( this_node->callbacks->size > 0 ) {
      strcpy( &callbacks[ islot * stride + stride - 1 ], "," );
    }
    else {
      strcpy( &callbacks[ islot * stride + stride - 1 ], ";" );
    }
    islot++;

    // add the callbacks for this node
    for (icallback = 0; icallback < this_node->callbacks->size; icallback++) {
      char* callback = vector_get(this_node->callbacks, icallback);
      length = strlen(callback);
      strcpy(&callbacks[ islot * stride ], callback);
      for (ichar = length; ichar < stride-1; ichar++) {
	strcpy( &callbacks[ islot * stride + ichar ], " " );
      }
      if ( icallback == this_node->callbacks->size - 1 ) {
	strcpy( &callbacks[ islot * stride + stride - 1], ";");
      }
      else {
	strcpy( &callbacks[ islot * stride + stride - 1], ",");
      }
      islot++;
    }
  }

  int ret = general_send( callbacks, count, MDI_CHAR, comm );
  free( callbacks );
  return ret;
}


/*! \brief Send the list of nodes
 *
 * If running with MPI, this function must be called only by rank \p 0.
 * The function returns \p 0 on a success.
 *
 * \param [in]       comm
 *                   MDI communicator associated with the intended recipient code.
 */
int send_node_list(MDI_Comm comm) {
  if ( intra_rank != 0 ) {
    mdi_error("Attempting to send node information from the incorrect rank");
  }
  int nnodes = nodes.size;
  int inode;
  int stride = MDI_COMMAND_LENGTH + 1;

  // allocate memory for the node list
  int count = nnodes * stride;
  char* node_list = malloc( count * sizeof(char) );

  // form the list of nodes
  for (inode = 0; inode < nnodes; inode++) {
    // add the name of this node to the list
    node* this_node = vector_get(&nodes, inode);
    int length = strlen(this_node->name);
    strcpy(&node_list[ inode * stride ], this_node->name);
    int ichar;
    for (ichar = length; ichar < stride-1; ichar++) {
      strcpy( &node_list[ inode * stride + ichar ], " " );
    }
    strcpy( &node_list[ inode * stride + stride - 1 ], "," );
  }

  int ret = general_send( node_list, count, MDI_CHAR, comm );
  free( node_list );
  return ret;
}


/*! \brief Send the number of supported commands
 *
 * If running with MPI, this function must be called only by rank \p 0.
 * The function returns \p 0 on a success.
 *
 * \param [in]       comm
 *                   MDI communicator associated with the intended recipient code.
 */
int send_ncommands(MDI_Comm comm) {
  if ( intra_rank != 0 ) {
    mdi_error("Attempting to send command information from the incorrect rank");
  }
  int ncommands = 0;
  int nnodes = nodes.size;
  int inode;
  int stride = MDI_COMMAND_LENGTH + 1;

  // determine the number of commands
  for (inode = 0; inode < nnodes; inode++) {
    node* this_node = vector_get(&nodes, inode);
    ncommands += this_node->commands->size;
  }

  int ret = general_send( &ncommands, 1, MDI_INT, comm );
  return ret;
}


/*! \brief Send the number of supported callbacks
 *
 * If running with MPI, this function must be called only by rank \p 0.
 * The function returns \p 0 on a success.
 *
 * \param [in]       comm
 *                   MDI communicator associated with the intended recipient code.
 */
int send_ncallbacks(MDI_Comm comm) {
  if ( intra_rank != 0 ) {
    mdi_error("Attempting to send callback information from the incorrect rank");
  }
  int ncallbacks = 0;
  int nnodes = nodes.size;
  int inode;
  int stride = MDI_COMMAND_LENGTH + 1;

  // determine the number of callbacks
  for (inode = 0; inode < nnodes; inode++) {
    node* this_node = vector_get(&nodes, inode);
    ncallbacks += this_node->callbacks->size;
  }

  int ret = general_send( &ncallbacks, 1, MDI_INT, comm );
  return ret;
}


/*! \brief Send the number of supported nodes
 *
 * If running with MPI, this function must be called only by rank \p 0.
 * The function returns \p 0 on a success.
 *
 * \param [in]       comm
 *                   MDI communicator associated with the intended recipient code.
 */
int send_nnodes(MDI_Comm comm) {
  if ( intra_rank != 0 ) {
    mdi_error("Attempting to send callback information from the incorrect rank");
  }
  int nnodes = nodes.size;
  int ret = general_send( &nnodes, 1, MDI_INT, comm );
  return ret;
}


/*! \brief Get information about the nodes of a particular code
 *
 * If running with MPI, this function must be called only by rank \p 0.
 * The function returns \p 0 on a success.
 *
 * \param [in]       comm
 *                   MDI communicator associated with the connection to the sending code.
 */
int get_node_info(MDI_Comm comm) {
  printf("~~~Getting node information\n");
  size_t stride = MDI_COMMAND_LENGTH + 1;
  communicator* this = vector_get(&communicators, comm-1);
  char* current_node = malloc( MDI_COMMAND_LENGTH * sizeof(char) );

  // get the number of nodes
  size_t nnodes;
  MDI_Send_Command("<NNODES",comm);
  MDI_Recv(&nnodes, 1, MDI_INT, comm);

  // get the nodes
  size_t list_size = nnodes * stride * sizeof(char);
  char* node_list = malloc( list_size );
  MDI_Send_Command("<NODES",comm);
  MDI_Recv(node_list, nnodes * stride, MDI_CHAR, comm);
  printf("DRIVER NODE_LIST: %s\n",node_list);

  // register the nodes
  int inode;
  int ichar;
  for (inode = 0; inode < nnodes; inode++) {
    // find the end of the node name
    char* name_start = &node_list[ inode * stride ];
    char* name_end = strchr( name_start, ' ' );
    int name_length;
    if (name_end == NULL) {
      name_length = MDI_COMMAND_LENGTH;
    }
    else {
      name_length = name_end - name_start;
    }
    if ( name_length > MDI_COMMAND_LENGTH ) {
      mdi_error("Error obtaining node information: could not parse node name");
    }

    // construct the name of the node
    char* node_name = malloc( MDI_COMMAND_LENGTH * sizeof(char) );
    strncpy( node_name, name_start, name_length );
    for (ichar = name_length; ichar < MDI_COMMAND_LENGTH; ichar++) {
      node_name[ichar] = '\0';
    }
    printf("DRIVER NODE: %d %d %s\n",inode,name_length,node_name);

    // register this node
    register_node(this->nodes, node_name);

    // free the memory for the node name
    free( node_name );
  }

  // get the number of commands
  int ncommands;
  MDI_Send_Command("<NCOMMANDS",comm);
  MDI_Recv(&ncommands, 1, MDI_INT, comm);

  // get the list of commands
  int count = ( ncommands + nnodes ) * stride;
  char* commands = malloc( count * sizeof(char) );
  MDI_Send_Command("<COMMANDS",comm);
  MDI_Recv(commands, count, MDI_CHAR, comm);
  printf("~~~COMMANDS: %d %s\n",ncommands,commands);

  // register the commands
  int node_flag = 1;
  for (inode = 0; inode < nnodes + ncommands; inode++) {
    // find the end of the name
    char* name_start = &commands[ inode * stride ];
    char* name_end = strchr( name_start, ' ' );
    int name_length;
    if (name_end == NULL) {
      name_length = MDI_COMMAND_LENGTH;
    }
    else {
      name_length = name_end - name_start;
    }
    if ( name_length > MDI_COMMAND_LENGTH ) {
      mdi_error("Error obtaining node information: could not parse command name");
    }

    // construct the name
    char* command_name = malloc( MDI_COMMAND_LENGTH * sizeof(char) );
    strncpy( command_name, name_start, name_length );
    for (ichar = name_length; ichar < MDI_COMMAND_LENGTH; ichar++) {
      command_name[ichar] = '\0';
    }
    printf("DRIVER COMMAND: %d %d %s\n",inode,name_length,command_name);

    if ( node_flag == 1 ) { // node
      // store the name of the current node
      strcpy( current_node, command_name );
    }
    else { // command
      // register this command
      register_command(this->nodes, current_node, command_name);
    }

    // determine whether the next name is for a node or a command
    if ( name_start[stride - 1] == ';' ) {
      printf("NODE\n");
    }
    else if ( name_start[stride - 1] == ',' ) {
      printf("COMMAND\n");
    }
    else {
      mdi_error("Error obtaining node information: could not parse delimiter");
    }

    // free the memory for the node name
    free( command_name );
  }


  // get the number of callbacks
  int ncallbacks;
  MDI_Send_Command("<NCALLBACKS",comm);
  MDI_Recv(&ncallbacks, 1, MDI_INT, comm);

  // get the list of callbacks
  count = ( ncallbacks + nnodes ) * stride;
  char* callbacks = malloc( count * sizeof(char) );
  MDI_Send_Command("<CALLBACKS",comm);
  MDI_Recv(callbacks, count, MDI_CHAR, comm);
  printf("~~~CALLBACKS: %d %s\n",ncallbacks,callbacks);

  // register the callbacks
  node_flag = 1;
  for (inode = 0; inode < nnodes + ncallbacks; inode++) {
    // find the end of the name
    char* name_start = &callbacks[ inode * stride ];
    char* name_end = strchr( name_start, ' ' );
    int name_length;
    if (name_end == NULL) {
      name_length = MDI_COMMAND_LENGTH;
    }
    else {
      name_length = name_end - name_start;
    }
    if ( name_length > MDI_COMMAND_LENGTH ) {
      mdi_error("Error obtaining node information: could not parse callback name");
    }

    // construct the name
    char* callback_name = malloc( MDI_COMMAND_LENGTH * sizeof(char) );
    strncpy( callback_name, name_start, name_length );
    for (ichar = name_length; ichar < MDI_COMMAND_LENGTH; ichar++) {
      callback_name[ichar] = '\0';
    }
    printf("DRIVER CALLBACK: %d %d %s\n",inode,name_length,callback_name);

    if ( node_flag == 1 ) { // node
      // store the name of the current node
      strcpy( current_node, callback_name );
    }
    else { // callback
      // register this callback
      register_callback(this->nodes, current_node, callback_name);
    }

    // determine whether the next name is for a node or a callback
    if ( name_start[stride - 1] == ';' ) {
      printf("NODE\n");
    }
    else if ( name_start[stride - 1] == ',' ) {
      printf("CALLBACK\n");
    }
    else {
      mdi_error("Error obtaining node information: could not parse delimiter");
    }

    // free the memory for the node name
    free( callback_name );
  }
  
  // free the memory
  free( node_list );
  free( commands );
  free( callbacks );
  free( current_node );

  return 0;
}


/*! \brief Get the node vector associated with a particular communicator
 *
 * The function returns the node vector for the communicator.
 *
 * \param [in]       comm
 *                   MDI communicator of the engine.  If comm is set to 
 *                   MDI_NULL_COMM, the function will return the node vector for the calling engine.
 */
vector* get_node_vector(MDI_Comm comm) {
  // get the vector of nodes associated with the communicator
  vector* node_vec;
  if ( comm == MDI_NULL_COMM ) {
    node_vec = &nodes;
  }
  else {
    communicator* this = vector_get(&communicators, comm-1);
    if ( this->nodes->size == 0 ) {
      // acquire node information for this communicator
      get_node_info(comm);
    }
    node_vec = this->nodes;
  }
  return node_vec;
}
