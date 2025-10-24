disp('defining the problem space parameters'); 

% maximum number of time steps to run FDTD simulation 
number_of_time_steps = 200000;

% A factor that determines duration of a time step wrt CFL limit
courant_factor = 0.9;

% Dimensions of a unit cell in x, y, and z directions (meters)
dx = 0.1e-3;
dy = 0.1e-3;
dz = 0.127e-3;

% ==<boundary conditions>========
% Here we define the boundary conditions parameters 
% 'pec' : perfect electric conductor
% 'cpml' : convolutional PML 
% if cpml_number_of_cells is less than zero
% CPML extends inside of the domain rather than outwards

boundary.type_xn = 'pec';
boundary.air_buffer_number_of_cells_xn = 0;
boundary.cpml_number_of_cells_xn = 8;

boundary.type_xp = 'pec';
boundary.air_buffer_number_of_cells_xp = 0;
boundary.cpml_number_of_cells_xp = 8;

boundary.type_yn = 'pec';
boundary.air_buffer_number_of_cells_yn = 0;
boundary.cpml_number_of_cells_yn = 10;

boundary.type_yp = 'pec';
boundary.air_buffer_number_of_cells_yp = 0;
boundary.cpml_number_of_cells_yp = 8;

boundary.type_zn = 'pec';
boundary.air_buffer_number_of_cells_zn = 0;
boundary.cpml_number_of_cells_zn = 8;

boundary.type_zp = 'pec';
boundary.air_buffer_number_of_cells_zp = 10;
boundary.cpml_number_of_cells_zp = 8;

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
material_types(3).eps_r   = 10.8;
material_types(3).mu_r    = 1;
material_types(3).sigma_e = 0;
material_types(3).sigma_m = 0;
material_types(3).color   = [0.94118 0.41176 0.88235];

