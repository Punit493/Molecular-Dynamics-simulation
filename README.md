# Pythoncodes


Microcanonical Ensemble Simulation using Lennard-Jones Potential

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

