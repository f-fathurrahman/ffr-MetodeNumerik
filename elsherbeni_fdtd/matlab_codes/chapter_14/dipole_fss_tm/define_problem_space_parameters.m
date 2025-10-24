disp('defining the problem space parameters'); 

% maximum number of time steps to run FDTD simulation 
number_of_time_steps = 5000;

% A factor that determines duration of a time step wrt CFL limit
courant_factor = 0.9;

% Dimensions of a unit cell in x, y, and z directions (meters)
dx = 0.5e-3;    
dy = 0.5e-3;   
dz = 0.5e-3;

% ==<periodic boundary simulation parameters>========
disp('Defining periodic boundary simulation'); 

periodic_boundary.kx = 20;
periodic_boundary.ky = 7.8;
periodic_boundary.mode = 'TM';
periodic_boundary.H_phi = 1;
periodic_boundary.E_x = 1;
periodic_boundary.E_y = 0.5;
periodic_boundary.source_z = 18e-3;
periodic_boundary.reflection_z = 16e-3;
periodic_boundary.transmission_z = -9e-3;

% ==<boundary conditions>========
% Here we define the boundary conditions parameters 
% 'pec' : perfect electric conductor
% 'cpml' : convolutional PML 
% 'pbc' : Periodic Boundary Condition 
% if cpml_number_of_cells is less than zero
% CPML extends inside of the domain rather than outwards

boundary.type_xn = 'pbc';
boundary.air_buffer_number_of_cells_xn = 0;
boundary.cpml_number_of_cells_xn = 8;

boundary.type_xp = 'pbc';
boundary.air_buffer_number_of_cells_xp = 0;
boundary.cpml_number_of_cells_xp = 8;

boundary.type_yn = 'pbc';
boundary.air_buffer_number_of_cells_yn = 0;
boundary.cpml_number_of_cells_yn = 10;

boundary.type_yp = 'pbc';
boundary.air_buffer_number_of_cells_yp = 0;
boundary.cpml_number_of_cells_yp = 8;

boundary.type_zn = 'cpml';
boundary.air_buffer_number_of_cells_zn = 20;
boundary.cpml_number_of_cells_zn = 10;

boundary.type_zp = 'cpml';
boundary.air_buffer_number_of_cells_zp = 40;
boundary.cpml_number_of_cells_zp = 10;

boundary.cpml_order = 3;
boundary.cpml_sigma_factor = 1.3;
boundary.cpml_kappa_max = 7;
boundary.cpml_alpha_min = 0;
boundary.cpml_alpha_max = 0.05;

% ===<material types>============
% Here we define and initialize the arrays of material types
% eps_r   : relative permittivity
% mu_r    : relative permeability
% sigma_e : electric conductivity
% sigma_m : magnetic conductivity

% air
material_types(1).eps_r   = 1;
material_types(1).mu_r    = 1;
material_types(1).sigma_e = 0;
material_types(1).sigma_m = 0;
material_types(1).color   = [1 1 1];

% PEC : perfect electric conductor;
material_types(2).eps_r   = 1;
material_types(2).mu_r    = 1;
material_types(2).sigma_e = 1e10;
material_types(2).sigma_m = 0;
material_types(2).color   = [1 0 0];

% a material type
material_types(3).eps_r   = 1;
material_types(3).mu_r    = 1;
material_types(3).sigma_e = 0;
material_types(3).sigma_m = 1e10;
material_types(3).color   = [0 1 0];

% substrate 
material_types(4).eps_r   = 2.2;
material_types(4).mu_r    = 1;
material_types(4).sigma_e = 0;
material_types(4).sigma_m = 0; 
material_types(4).color   = [1 1 0];

