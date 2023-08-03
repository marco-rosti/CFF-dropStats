Ianto Cannon, 25 July 2023

Code adapted from Giovanni Soligo. Uses a flood-fill algorithm to count the number of droplets in a 3D binary file such as: dataVOF000001700000. The code uses mpi to run in parallel. Each process searches a box-shaped region in the full domain. An mpi window is used to share large fields like the velocity or surfactant concentration. Calculates statistics such as centre of mass, moments of inertia, velocity, area, surfactant concentration, and the Euler characteristic of each drop.
Run with command 'sbatch launchDropCount.sh'.

Okinawa Instiute of Science and Technology, Complex Fluids and Flows unit - Rosti unit
