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

#include "OpenCLExampleKernels.h"
#include "OpenCLExampleKernelSources.h"
#include "openmm/internal/ContextImpl.h"
#include "openmm/opencl/OpenCLBondedUtilities.h"
#include "openmm/opencl/OpenCLForceInfo.h"

using namespace ExamplePlugin;
using namespace OpenMM;
using namespace std;

class OpenCLExampleForceInfo : public OpenCLForceInfo {
public:
    OpenCLExampleForceInfo(const ExampleForce& force) : OpenCLForceInfo(0), force(force) {
    }
    int getNumParticleGroups() {
        return force.getNumBonds();
    }
    void getParticlesInGroup(int index, vector<int>& particles) {
        int particle1, particle2;
        double length, k;
        force.getBondParameters(index, particle1, particle2, length, k);
        particles.resize(2);
        particles[0] = particle1;
        particles[1] = particle2;
    }
    bool areGroupsIdentical(int group1, int group2) {
        int particle1, particle2;
        double length1, length2, k1, k2;
        force.getBondParameters(group1, particle1, particle2, length1, k1);
        force.getBondParameters(group2, particle1, particle2, length2, k2);
        return (length1 == length2 && k1 == k2);
    }
private:
    const ExampleForce& force;
};

OpenCLCalcExampleForceKernel::~OpenCLCalcExampleForceKernel() {
    if (params != NULL)
        delete params;
    if (MDIForces != NULL)
        delete MDIForces;
}

void OpenCLCalcExampleForceKernel::initialize(const System& system, const ExampleForce& force) {
    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cl.getContextIndex()+1)*force.getNumBonds()/numContexts;
    numBonds = endIndex-startIndex;
    map<string, string> replacements;

    // OBSOLETE: OLD EXAMPLE CODE
    vector<vector<int> > atoms(numBonds, vector<int>(2));
    if (numBonds > 0) {
      params = OpenCLArray::create<mm_float2>(cl, numBonds, "bondParams");
      vector<mm_float2> paramVector(numBonds);
      for (int i = 0; i < numBonds; i++) {
        double length, k;
        force.getBondParameters(startIndex+i, atoms[i][0], atoms[i][1], length, k);
        paramVector[i] = mm_float2((cl_float) length, (cl_float) k);
      }
      params->upload(paramVector);
      replacements["PARAMS"] = cl.getBondedUtilities().addArgument(params->getDeviceBuffer(), "float2");
    }

    // Create addForces kernel
    map<string, string> defines;
    defines["NUM_ATOMS"] = cl.intToString(cl.getNumAtoms());
    defines["PADDED_NUM_ATOMS"] = cl.intToString(cl.getPaddedNumAtoms());
    cl::Program program = cl.createProgram(OpenCLExampleKernelSources::mdiForce, defines);
    addForcesKernel = cl::Kernel(program, "addForces");

    // Upload MDIForces
    MDIForces = OpenCLArray::create<cl_float>(cl, 3*system.getNumParticles(), "MDIForces");
    vector<cl_float> MDIForcesVector(3*system.getNumParticles());
    for (int i = 0; i < 3*system.getNumParticles(); i++) {
      MDIForcesVector[i] = (cl_float) 0.0;
    }
    MDIForces->upload(MDIForcesVector);
    replacements["MDIADD"] = cl.getBondedUtilities().addArgument(MDIForces->getDeviceBuffer(), "float");

    // OBSOLETE: OLD EXAMPLE CODE
    if (numBonds > 0) {
      cl.getBondedUtilities().addInteraction(atoms, cl.replaceStrings(OpenCLExampleKernelSources::exampleForce, replacements), force.getForceGroup());
    }
    cl.addForce(new OpenCLExampleForceInfo(force));
}

double OpenCLCalcExampleForceKernel::execute(ContextImpl& context, bool includeForces, bool includeEnergy) {
    // Upload the forces to the device.

  /*
    int numParticles = context.getSystem().getNumParticles();
    if (cl.getUseDoublePrecision()) {
      double* buffer = (double*) cl.getPinnedBuffer();
      for (int i = 0; i < numParticles; ++i) {
	//const Vec3& p = forces[i];
	Vec3 p;
	p[0] = 0.0;
	p[1] = 0.0;
	p[2] = 0.0;
	buffer[3*i] = p[0];
	buffer[3*i+1] = p[1];
	buffer[3*i+2] = p[2];
      }
    }
    else {
      float* buffer = (float*) cl.getPinnedBuffer();
      for (int i = 0; i < numParticles; ++i) {
	//const Vec3& p = forces[i];
	Vec3 p;
	p[0] = 0.0;
	p[1] = 0.0;
	p[2] = 0.0;
	buffer[3*i] = (float) p[0];
	buffer[3*i+1] = (float) p[1];
	buffer[3*i+2] = (float) p[2];
      }
    }
  */

    addForcesKernel.setArg<cl::Buffer>(0, MDIForces->getDeviceBuffer());
    addForcesKernel.setArg<cl::Buffer>(1, cl.getForceBuffers().getDeviceBuffer());
    addForcesKernel.setArg<cl::Buffer>(2, cl.getAtomIndexArray().getDeviceBuffer());
    cl.executeKernel(addForcesKernel, cl.getNumAtoms());

    return 0.0;
}

void OpenCLCalcExampleForceKernel::copyParametersToContext(ContextImpl& context, const ExampleForce& force) {
    int numContexts = cl.getPlatformData().contexts.size();
    int startIndex = cl.getContextIndex()*force.getNumBonds()/numContexts;
    int endIndex = (cl.getContextIndex()+1)*force.getNumBonds()/numContexts;
    if (numBonds != endIndex-startIndex)
        throw OpenMMException("updateParametersInContext: The number of bonds has changed");
    if (numBonds == 0)
        return;
    
    // Record the per-bond parameters.
    
    vector<mm_float2> paramVector(numBonds);
    for (int i = 0; i < numBonds; i++) {
        int atom1, atom2;
        double length, k;
        force.getBondParameters(startIndex+i, atom1, atom2, length, k);
        paramVector[i] = mm_float2((cl_float) length, (cl_float) k);
    }
    params->upload(paramVector);
    
    // Mark that the current reordering may be invalid.
    
    cl.invalidateMolecules();
}

