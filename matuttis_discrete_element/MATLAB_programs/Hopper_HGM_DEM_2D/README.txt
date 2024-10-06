This directory contains the MATLAB m-files for a discrete-element 
code with polygonal particles as companion program for 

Understanding the Discrete Element Method: 
Simulation of Non-Spherical Particles for 
Granular and Multi-body Systems

Hans-Georg Matuttis, Jian Chen

Wiley 2014

!!!!!  THIS PROGRAM AND ITS FUNCTIONS COME WITHOUT WARRANTY  !!!!
% REVISION HISTORY: 26-May-2014 H.-G. Matuttis

The purpose of this program is not to supply a perfectly running code 
(though it should be sufficiently bug-free to use it for simulations),
but to illustrate the concepts from the book, including
various programming approaches, documentation examples etc.
 
PROGRAM AUTHORS: Shi Han Ng, Hans-Georg Matuttis

The main program is dr_particledry_hg.m
The program flow corresponds to the typical program-flow outlined in 
the book. The page- and section numbers in the comments refer to our book.

The simulation geometry (a hopper) is contained in the file
walldata.dat.
It should be relatively easy to modify that file for other geometries
with non-moving walls: Don't forget to change the number of the
walls at the beginning of the files when you add or remove walls.

README.txt           this file

dr_HGM_DEM_2D.m      the main file

demInitial.m         initializes the positions, masses etc. for particles and masses, 
                     calls mass_momentinertia.m, read the geometry-data from  walldata.dat

draw_particle.m      draws the particles

gearPredictor.m      Predictor-Corrector Pair after GEAR (Backward-difference formular, BDF)
gearCorrector.m

updatevertices.m     updates the particle corners/vertices after a predictor step, calls 
                     rotation.m

neighborhood_algorithm.m Neighborhood-algorithm, calls  init_neigh_coord_poly.m at
                     the first timestep and after that update_boundbox.m, update_a_list.m
                     start_search_x2D.m, start_search_y2D.m, inc_sort_x2D.m, inc_sort_y2D.m
                  	 in every timestep

deminteraction.m 	 computes the interaction, calls first  overlap.m (which calls centerOfMass.m) 
                     and then findForce.m  

                
               
              
        
               
                  
               
                
             
               

 