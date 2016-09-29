# particletrack

Diffusing particle simulation. Particles move randomly inside a 3D volume and boundaries can be periodic or absorbing. Flux boundaries can also be used in order to create concentration gradients across the volume.

Cells can be placed within the simulation volume. The number of particles within a cell can be counted in order to assign a cell a polarization vector.

## Simulation Procedure

First, the system is allowed to reach an equilibrium. Particles diffuse and are produced for `ntItl` time-steps.

Next, further time evolution of the system can be done. Particles diffuse and are produced for `ntTotal` time-steps. During each time-step cell polarization is calculated based on the diffusing molecule population within each cell. The total cluster polarization is stored for that time-step.

## Output

The time-averaged total cluster polarization is output from the program by using the `wrtPlrTime` subroutine. During each time-step the total cluster polarization is stored in the array `timePolar(1:3,nt)`. At the end of a run `wrtPlrTime` is called and the mean of each component of the polarization vector is calculated. The mean is written to the file `mean###.200`.
