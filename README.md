# Pythoncodes


A) Microcanonical Ensemble Simulation using Lennard-Jones Potential:-

Introduction: -
This Python code implementation is for simulating a microcanonical ensemble (constant energy, volume, and number of particles) using the Lennard-Jones potential to model interactions between particles in a cubic simulation box with periodic boundary conditions.
Key Points of the Code
1. Constants & Initialization:
Physical constants such as Lennard-Jones parameters (σ, ε), particle mass (m), and cutoff radius (rc) are defined. The simulation initializes a fixed number of particles (N = 512) in a cubic box with dimensions Lx, Ly, and Lz. Initial particle positions are arranged in a lattice structure, and velocities are initialized randomly while ensuring a net-zero center of mass velocity.
2. Force Calculation:
The Lennard-Jones potential is used to compute forces between particles within the cutoff distance (rc). Periodic boundary conditions are applied to account for particles leaving the box, ensuring continuity.
3. Velocity-Verlet Integration:
The Velocity-Verlet algorithm is used to iteratively update particle positions and velocities, providing numerical stability. At each step, the kinetic energy (KE) and potential energy (PE) of the system are calculated.
4. Simulation:
The simulation runs for a specified number of iterations (ni = 1500). It records energy values and updates particle positions and velocities while ensuring energy conservation throughout the process.
5. Visualization:
Several visual outputs are generated, including:
- 3D plots of the initial and final positions of particles.
- A graph showing the evolution of kinetic, potential, and total energy over the iterations.
- A histogram of particle speeds to illustrate the velocity distribution.


B) Canonical Ensemble Simulation:-

This code simulates the dynamics of a system of 512 particles in a 3D box using molecular dynamics (MD) with the Lennard-Jones potential. The simulation tracks particle positions, velocities, and forces, while applying periodic boundary conditions and a Berendsen thermostat to control temperature. Key steps include:

1. Initialization:
Particle positions are set on a 3D grid.
Velocities are initialized randomly and adjusted to have zero net momentum.
2. Force Calculation:
The Lennard-Jones potential is used to compute forces and potential energy between particles within a cutoff radius.
3. Integration:
Positions and velocities are updated iteratively using the Verlet integration scheme, incorporating the thermostat.
4. Energy Tracking:
Potential, kinetic, and total energies are computed and tracked.
5. Visualization:
Initial and final particle positions are plotted in 3D.
Energy evolution and speed distribution are visualized.
This provides insights into the system's energy conservation and thermalization dynamics.




