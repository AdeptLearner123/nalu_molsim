from openmm import *
from openmm.app import *
from openmm.unit import *

forcefield = ForceField('amber99sb.xml', 'amber99_obc.xml')
system = forcefield.createSystem()

# one particle, no forces
#system = System()
system.addParticle(0.0)  # massless
system.addForce(CustomExternalForce("0"))  # dummy force so Context is valid

#integrator = VerletIntegrator(0.002)
integrator = LangevinIntegrator(300*kelvin, 1/picosecond, 0.002*picoseconds)

simulation = Simulation(Topology(), system, integrator)
# add a position
simulation.context.setPositions([[0,0,0]*nanometer])

for step in range(10):
    pos = simulation.context.getState(getPositions=True).getPositions()[0]
    print(step, pos)
    simulation.step(1)