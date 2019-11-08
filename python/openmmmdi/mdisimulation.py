import simtk.openmm as mm
import simtk.openmm.app as mmapp
import simtk.unit as mmunit
from .openmmmdi import ExampleForce, MDIServer

class MDISimulation(mmapp.Simulation):
    def __init__(self, mdiOptions, topology, system, integrator, platform=None, platformProperties=None, state=None):
        ## Add the MDI force
        self.mdi_force = ExampleForce(mdiOptions)
        system.addForce(self.mdi_force)

        ## Create the simulation
        mmapp.Simulation.__init__(self, topology=topology, system=system, integrator=integrator, 
                                  platform=platform, platformProperties=platformProperties, state=state)

    def runMDI(self):
        self.mdi_force.setActive(True, self.context)

        command = "@GLOBAL"
        current_simulation = ""
        while command != "EXIT":

            if command == "@GLOBAL":
                current_simulation = ""
                self.mdi_force.mdiListen("@GLOBAL", self.context)

            elif command == "@INIT_MD":
                current_simulation = "MD"
                self.mdi_force.mdiListen("@INIT_MD", self.context)

            elif command == "@":
                if current_simulation == "MD":
                    # Place an @ENERGY node between steps
                    if previous_node != "@ENERGY":
                        self.mdi_force.mdiListen("@ENERGY", self.context)
                    else:
                        self.step(1)
                else:
                    raise Exception("Cannot proceed to the next node unless an MD simulation has been started")

            elif command == "@ENERGY":
                if current_simulation == "MD":
                    if previous_node == "@ENERGY":
                        # The engine was already at @ENERGY, so do another step first
                        self.step(1)
                    self.mdi_force.mdiListen("@ENERGY", self.context)
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

            previous_node = self.mdi_force.getPreviousNode(self.context)
            command = self.mdi_force.getTargetNode(self.context)
