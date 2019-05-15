# Quantum Scripts

Scripts of/for simple quantum system simulations using MATLAB.

### Simple Quantum Mechanics

- gaussianwave: Create a normalized gaussian wavefunction;
- get_fs_mat: Create finite stencil differential matrix for kinetic energy operator;
- get_laplacian: Create multidimensional Laplacian operator (only the kernel) for kinetic energy;
- get_potential: Quickly create some common potentials;
- runhm: Propogating wavefunction with given Hamiltonian by eigenvector decomposition;
- runho: Propogating wavefunction with given Hamiltonian by ODE (for extermely large matrix);
- wavepacket: Example of wavepacket evolution;
- wavepacket2d: Example of wavepacket evolution in 2D;
- wavepacketms: Example of wavepacket evolution with multiple potential surfaces;
- tunnel: Example of quantum tunneling using Fourier transform;
- scatter: Example of scattering state simulation;

### Relaxation

- get_env: Generate a list of environment energy with fixed interval;
- get_gamma: Calculate relaxation rate (Gamma) and energy shift from a given environment;
- relaxation: Example of relaxation of a quantum state;

### Boson Dynamics

- get_env_boson: Generate a boson environment, with equal-interval energy and corresponding population;
- runam: Propogating boson system with given Hamiltonian by eigenvector decomposition;
- runao: Propogating boson system with given Hamiltonian by ODE;
- relaxation_boson: Example of relaxation of boson system;
