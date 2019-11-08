import simtk.openmm as mm
import simtk.openmm.app as mmapp
import simtk.unit as mmunit
from .openmmmdi import ExampleForce, MDIServer

class MDISimulation(mmapp.Simulation):
    def __init__(self, mdiOptions, topology, system, integrator, platform=None, platformProperties=None, state=None):
        ## Create an MDI server object
        #self.server = MDIServer()
        #self.server.init(mdiOptions)

        ## Add the MDI force
        self.mdi_force = ExampleForce(mdiOptions)
        #for i in range(1000):
        #    force.addBond(i, i, 1.0, 10.0)
        system.addForce(self.mdi_force)

        ## NOTE: AT THIS POINT, SHOULD SET A FLAG THAT PREVENTS THE FORCES FROM ACTUALLY LISTENING FOR COMMANDS
        ## MDIFORCE SHOULDN'T LISTEN FOR COMMANDS UNTIL RUNMDI IS CALLED

        mmapp.Simulation.__init__(self, topology=topology, system=system, integrator=integrator, 
                                  platform=platform, platformProperties=platformProperties, state=state)

    def runMDI(self):
        print("runMDI")
        #self.server.run()
        #for istep in range(10):
        #    #new_state = self.context.getState(getEnergy = True)
        #    #print("      ------------------- " + str(new_state.getKineticEnergy()))
        #    self.step(1)

        #self.step(1)
        #new_state = self.context.getState(getEnergy = True)
        #new_state = self.context.getState(getEnergy = True)
        #new_state = self.context.getState(getEnergy = True)
        #new_state = self.context.getState(getEnergy = True)

        self.mdi_force.setActive(True, self.context)

        command = "@GLOBAL"
        current_simulation = ""
        while command != "EXIT":

            if command == "@GLOBAL":
                current_simulation = ""
                command = self.mdi_force.mdiListen("@GLOBAL", self.context)

            elif command == "@INIT_MD":
                current_simulation = "MD"
                command = self.mdi_force.mdiListen("@INIT_MD", self.context)

            elif command == "@":
                if current_simulation == "MD":
                    self.step(1)
                else:
                    raise Exception("Cannot proceed to the next node unless an MD simulation has been started")

            elif command == "@ENERGY":
                if current_simulation == "MD":
                    command = self.mdi_force.mdiListen("@ENERGY", self.context)
                else:
                    raise Exception("Cannot proceed to the @ENERGY node unless an MD simulation has been started")

            elif command == "@UPDATE":
                if current_simulation == "MD":
                    self.step(1)
                else:
                    raise Exception("Cannot proceed to the @UPDATE node unless an MD simulation has been started")

            elif command == "@FORCES":
                if current_simulation == "MD":
                    self.step(1)
                else:
                    raise Exception("Cannot proceed to the @FORCES node unless an MD simulation has been started")

            elif command == "EXIT":
                pass

            else:
                raise Exception("Unsupported MDI command: " + str(command))

            print( "   Engine outer command: " + str(command) )
