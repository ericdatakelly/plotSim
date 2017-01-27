# plotSim
Calibration and visualization for Crystallize3D


This figure shows simulated crystals resulting from Crystallize3D modeling, displayed with plotSim.

<img src="https://github.com/ericdavidkelly/plotSim/blob/master/SimulationResults/test_160907_1709_Plot3.png" width="600"\>



Crystallize3D is a finite-difference element model that simulates solid-state diffusion representing conditions of crystal growth in the Earth's crust (Ketcham and Carlson, 2012).  The chemical reactions that produce the crystals are literally impossible to reproduce in physical lab experiments due to time and temperature limitations (too fast and hot in the lab), so we use numerical techniques to overcome that problem.  plotSim enhances the sensitivity of Crystallize3D and provides analysis of simulation results.  


A sample is scanned using computed tomography (CT scan), and slices of the sample are stacked into a 3D rendering of the volume, which allows quantification of the spatial features of the sample, notably the locations and sizes of all crystals.


<img src="https://github.com/ericdavidkelly/plotSim/blob/master/Extras/sample.jpg" width="400"\>


<img src="https://github.com/ericdavidkelly/plotSim/blob/master/Extras/ct_slice.jpg" width="250"\>


<img src="https://github.com/ericdavidkelly/plotSim/blob/master/Extras/ct_rendering.jpg" width="600"\>




This figure shows 2D images of a model through four time steps.

<img src="https://github.com/ericdavidkelly/plotSim/blob/master/Extras/2D_model_images.png"/>


In the past, only the final state of the system was explicitly fit by a simulation, primarily using statistical analysis of the resulting crystal locations and sizes.  One of the measures is the degree of ordering of crystal centers (vs. pseudorandom locations), which can be determined using correlation functions that evaluate the nearest-neighbor distances of all crystals in the sample and the simulation.  

To visualize this, imagine trees in a forest and their root systems competing for water.  A seedling is unlikely to get enough water and sunlight to grow if it is next to a big tree because the big tree's roots grab up too much of the water and its branches block sunlight.  However, if the seed begins growth away from the big tree, it may have more water and light, allowing it to grow to full size.  Diffusion gradients that exist around crystals have an analogous effect, so the degree of crystal ordering (neighbors are consistently separated) relative to randomness (neighbors are close and far away) can constrain diffusion gradients critical for understanding crystal growth.  The sizes of the crystals can also be used to statistically test the model fit.  This figure shows an example.


<img src="https://github.com/ericdavidkelly/plotSim/blob/master/Extras/cfs.png"/>




Although these measures of fit are robust for constraining the final state of the system, several paths to the final state are possible, so constraints on the start and intermediate steps are highly desirable for producing much more sensitive models.  

The additional constraints provided by plotSim result in an elegant solution.  Let’s use another tree analogy, this time from the rings of a tree.  Starting at the core of the tree, each consecutive ring represents additional growth stages.  These rings contain abundant data on the timing of growth, the rate of growth, and the environment in which the tree grew, including climate conditions.

<img src="https://github.com/ericdavidkelly/plotSim/blob/master/Extras/Tree_rings.jpg"/>
(image source unknown)

In crystals, chemical changes are recorded in zones within each crystal, much like the rings of a tree.  These reveal the time at which each crystal began to grow and how fast it grew at any point in its history, allowing many more properties of the system to be accessed (e.g., flux of material supplied to each crystal, sample-wide changes in growth rates, etc.).  By simply fitting chemical profiles to a group of representative crystals using plotSim, we can simulate crystallization in the highest fidelity to date for these types of samples.  This allows extremely well-constrained determinations of diffusion and other kinetic parameters of crystallization that are fundamental to understanding Earth processes, such as ore formation (ore provides stuff to make cell phones and other things).


The images below represent output from plotSim, including the simple chemical zoning calibration (shown in "Measured garnet zoning profile", spline fit) and a key measure of model fit (shown in "Relative nucleation time").

<img src="https://github.com/ericdavidkelly/plotSim/blob/master/SimulationResults/test_160907_1709_Plot4.png" width="700"\>

<img src="https://github.com/ericdavidkelly/plotSim/blob/master/SimulationResults/PM4_120619_1554_rnd1_Plot1.png" width="700"\>

