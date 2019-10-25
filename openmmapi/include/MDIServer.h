#ifndef OPENMM_MDISERVER_H_
#define OPENMM_MDISERVER_H_

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
#include "openmm/Context.h"
#include "openmm/System.h"
#include "openmm/Force.h"
#include "openmm/Kernel.h"
#include <vector>
#include "internal/windowsExportExample.h"
#include "openmm/NonbondedForce.h"

namespace ExamplePlugin {

/**
 * This class implements an MDI server that listens for commands from an external driver.
 */

class OPENMM_EXPORT_EXAMPLE MDIServer {
public:
    /**
     * Create an MDIServer.
     */
    MDIServer();
    /**
     * Initialize the MDIServer.
     */
    void init(std::string mdi_options);
    /**
     * Listen for commands from the external driver.
     */
    void listen(OpenMM::ContextImpl& context, OpenMM::Kernel& kernel, std::string node);
    /**
     * Get the NonbondedForce.
     */
    const OpenMM::NonbondedForce* get_nonbonded_force(OpenMM::ContextImpl& context);
    /**
     * Respond to <COORDS.
     */
    std::vector<OpenMM::Vec3> send_coords(OpenMM::ContextImpl& context);
    /**
     * Respond to <VELOCITIES.
     */
    std::vector<OpenMM::Vec3> send_velocities(OpenMM::ContextImpl& context);
    /**
     * Respond to <FORCES.
     */
    std::vector<OpenMM::Vec3> send_forces(OpenMM::ContextImpl& context);
    /**
     * Respond to <TIME.
     */
    double send_time(OpenMM::ContextImpl& context);
    /**
     * Respond to <NATOMS.
     */
    int send_natoms(OpenMM::ContextImpl& context);
    /**
     * Respond to <CHARGES.
     */
    std::vector<double> send_charges(OpenMM::ContextImpl& context);
    /**
     * Respond to <DIMENSIONS.
     */
    std::vector<int> send_dimensions(OpenMM::ContextImpl& context);
    /**
     * Respond to <CELL.
     */
    std::vector<double> send_cell(OpenMM::ContextImpl& context);
    /**
     * Respond to >CELL.
     */
    void recv_cell(OpenMM::ContextImpl& context);
    /**
     * Respond to <ENERGY.
     */
    double send_energy(OpenMM::ContextImpl& context);
    /**
     * Respond to <MASSES.
     */
    std::vector<double> send_masses(OpenMM::ContextImpl& context);
    /**
     * Respond to >FORCES.
     */
    void recv_forces(OpenMM::ContextImpl& context, OpenMM::Kernel& kernel);
    /**
     * Additional responses needed:
     * @
     * <@
     * <COMMANDS
     * @COORDS
     * >COORDS
     * <KE
     * <PE
     * <KE_NUC
     * <PE_NUC
     * EXIT
     * @FORCES
     * +FORCES
     * >FORCES
     * <GLOBAL
     * @INIT_MD
     * @INIT_OPTG
     * <NCOMMANDS
     * @PRE-FORCES
     * >VELOCITIES
     */
    //private:
    //OpenMM::Kernel *kernel;
};


} // namespace ExamplePlugin

#endif /*OPENMM_MDISERVER_H_*/
