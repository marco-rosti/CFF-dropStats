2021 Aug 11
Ianto Cannon

Code from Giovanni Soligo which counts the number of droplets in a 3D binary file such as: dataVOF000001700000
Edited by Ianto Cannon to calculate droplet morphology
Each process is assigned a separate timestep, so ntasks in launchDropCount.sh should be less than or equal to the number of timesteps.

Then run with command 'sbatch launchDropCount.sh'
