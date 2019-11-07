import simtk.openmm as mm
import simtk.openmm.app as mmapp
import simtk.unit as mmunit
from .openmmmdi import ExampleForce, MDIServer

class MDISimulation(mmapp.Simulation):
    def __init__(self, mdiOptions, topology, system, integrator, platform=None, platformProperties=None, state=None):
        ## Create an MDI server object
        self.server = MDIServer()
        self.server.init(mdiOptions)

        ## Add the MDI force
        force = ExampleForce(mdiOptions, self.server)
        #for i in range(1000):
        #    force.addBond(i, i, 1.0, 10.0)
        system.addForce(force)

        ## NOTE: AT THIS POINT, SHOULD SET A FLAG THAT PREVENTS THE FORCES FROM ACTUALLY LISTENING FOR COMMANDS
        ## MDIFORCE SHOULDN'T LISTEN FOR COMMANDS UNTIL RUNMDI IS CALLED

        mmapp.Simulation.__init__(self, topology=topology, system=system, integrator=integrator, 
                                  platform=platform, platformProperties=platformProperties, state=state)

    def runMDI(self):
        self.server.run()
        for istep in range(10):
            #new_state = self.context.getState(getEnergy = True)
            #print("      ------------------- " + str(new_state.getKineticEnergy()))
            self.step(1)
