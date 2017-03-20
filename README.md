# particletrack

[![DOI](https://zenodo.org/badge/64701444.svg)](https://zenodo.org/badge/latestdoi/64701444)

Diffusing particle simulation. Particles move randomly inside a 3D volume and boundaries can be periodic or absorbing. Flux boundaries can also be used in order to create concentration gradients across the volume.

The most up-to-date version of the code can be found at the code's [GitHub repo](https://github.com/varennes/particletrack).

Cells can be placed within the simulation volume. The particles within a cell are tracked in order to assign a cell a polarization vector. Two principal methods of polarization are used. One referred to as *many wrongs* (also known as individual-based chemotaxis) in which cells individually make a weighted average of the number particles within their bodies. Particles at the front of the cell are given positive weights, whereas particles at the back are weighted negatively. The other is referred to as *emergent chemotaxis* in which cells count the all encolsed particles with equal weight and polarize in the direction away from all their neighbors.

Currently, the master branch can be used for creating linear particles concentration profiles. The [**`exp`**](https://github.com/varennes/particletrack/tree/exp) branch can be used for creating exponential concentration profiles.

## Simulation Procedure

First, the system is allowed to reach an equilibrium. Particles diffuse and are produced for `ntItl` time-steps.

Next, further time evolution of the system can be done. Particles diffuse and are produced for `ntTotal` time-steps. During each time-step cell polarization is calculated based on the diffusing particle population within each cell. The total collective polarization is stored for that time-step. This procedure can be repeated for many simulation instances for `runTotal` number of times.

## Output

At each time-step `nt` the emergent chemotaxis and many wrongs total polarization are stored in `timePolarEC(:,nt)` and `timePolarMW(:,nt)`, respectively. The averaged over all `ntTotal` time-steps is output from the program by the `wrtPlrECMW` subroutine. At the end of each instance the subroutine `wrtPlrECMW` is called, and the total mean polarization vector for emergent chemotaxis and many wrongs are written to the files `ec###.dat` and `mw###.dat`.
