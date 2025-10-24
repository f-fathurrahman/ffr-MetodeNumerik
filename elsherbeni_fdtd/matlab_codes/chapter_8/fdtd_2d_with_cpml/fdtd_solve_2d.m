       % initialize the matlab workspace
 clear all; close all; clc;

% define the problem 
define_problem_space_parameters_2d;
define_geometry_2d;
define_sources_2d;
define_output_parameters_2d;

% initialize the problem space and parameters 
initialize_fdtd_material_grid_2d;
initialize_fdtd_parameters_and_arrays_2d;
initialize_sources_2d;
initialize_updating_coefficients_2d;
initialize_boundary_conditions_2d;
initialize_output_parameters_2d;
initialize_display_parameters_2d;

% draw the objects in the problem space 
draw_objects_2d;

% FDTD time marching loop
run_fdtd_time_marching_loop_2d;

% display simulation results
post_process_and_display_results_2d;
