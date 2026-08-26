# ClusterConstrict

### lptt-log-conform.h vs lptt-log-conform-embed.h

The "embed" variant is the same linear-PTT log-conformation solver, adapted
for embedded boundaries (embed.h):

* In the lptt-log-conform-embed.h variant, every stress update loop is protected by a loop that skips all solid/cut cells part of the embedded boundary.

# LPTTValidation : 2-d Axisymmetric Pipe flow Solver

A Basilisk C solver for simulating fully developed axisymmetric pipe flow of a viscoelastic matrix described by the L-PTT model with a Newtonian solvent contribution. 

This solver uses the log-conformation formulation to maintain numerical stability at high weissenberg numbers. The implementation is derived frimm the log-conform.h header present in the Basilisk source code.

In particular validation against analytical solutions for fully developed velocity profile described by Cruz et al., (2005)[1], is done using this solver.



The solver is fully non-dimensionalized.





## Compilation & Execution
 **Compile and run the solver:**
   ```bash
   qcc -O2 -Wall -fopenmp code.c -o run -lm
   
   export OMP_NUM_THREADS=12 
   
   ./run

  Number of threads subject to change based on the machine used to operate the solver.
 ```


 ## DropInChannel : 2-D Axisymmetric Channel With Viscoelastic Droplet


This Basilisk C script simulates the axisymmetric hydrodynamic migration of a droplet through a cylindrical channel. Configured in a purely viscous limit ($Wi = 0$, $\beta = 1.0$), it serves as a baseline validation setup to reproduce Newtonian droplet deformation and kinematic results from Nath et al. (2017)[2].


The current non-dimensional parameters replicate the results found in figure 11.


## Compilation and execution

```bash
qcc -g -Wall -O2 -fopenmp -disable-dimensions DropChannel.c     -o run     -L$HOME/basilisk/src/gl     -lglutils -lfb_tiny -lm

export OMP_NUM_THREADS=12 

./run

  Number of threads subject to change based on the machine used to operate the solver.
  
 ```


## Ongoing_Force_Model

This folder contains two .mp4 files.
* Force_test.mp4 : Shows a test case where 2 droplets are initialised at a distance of 1 unit away from each other, and given an initial velocity of 1 unit. The video demonstrates the current results of the force model, which is currently being validated. 

* Channel_geometry_2Drop.mp4 : Shows the proposed geometry for investigation. Two drops are initialised and travel through the constriction number, at present the force model is not integrated. The $Ca$ number is 0.05.

# References

[1] D. O. A. Cruz, F. T. Pinho, and P. J. Oliveira, "Analytical solutions for fully developed laminar flow of some viscoelastic liquids with a Newtonian solvent contribution," *Journal of Non-Newtonian Fluid Mechanics*, 132, 28–35, 2005.

[2] B. Nath, G. Biswas, A. Dalal, and K. C. Sahu, "Migration of a droplet in a cylindrical tube in the creeping flow regime," *Physical Review E*, 95, 033110, 2017. doi:10.1103/PhysRevE.95.033110.
