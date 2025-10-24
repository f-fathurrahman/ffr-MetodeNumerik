disp('defining the problem space parameters'); 

% maximum number of time steps to run FDTD simulation 
number_of_time_steps = 330;

% A factor that determines duration of a time step wrt CFL limit
courant_factor = 0.9;

% Dimensions of a unit cell in x, y, and z directions (meters)
dx = 0.001;
dy = 0.001;
dz = 0.001;

% ==<boundary conditions>========
% Here we define the boundary conditions parameters 
% 'pec' : perfect electric conductor
% 'cpml' : convolutional PML 
% if cpml_number_of_cells is less than zero
% CPML extends inside of the domain rather than outwards

boundary.type_xn = 'cpml';
boundary.air_buffer_number_of_cells_xn = 10;
boundary.cpml_number_of_cells_xn = 8;

boundary.type_xp = 'cpml';
boundary.air_buffer_number_of_cells_xp = 10;
boundary.cpml_number_of_cells_xp = 8;

boundary.type_yn = 'cpml';
boundary.air_buffer_number_of_cells_yn = 10;
boundary.cpml_number_of_cells_yn = 10;

boundary.type_yp = 'cpml';
boundary.air_buffer_number_of_cells_yp = 10;
boundary.cpml_number_of_cells_yp = 8;

boundary.type_zn = 'cpml';
boundary.air_buffer_number_of_cells_zn = 10;
boundary.cpml_number_of_cells_zn = 8;

boundary.type_zp = 'cpml';
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
material_types(3).eps_r   = 3;
material_types(3).mu_r    = 1;
material_types(3).sigma_e = 0;
material_types(3).sigma_m = 1e10;
material_types(3).color   = [0.94118 0.41176 0.88235];

% a material type
material_types(4).eps_r   = 4;
material_types(4).mu_r    = 1;
material_types(4).sigma_e = 0;
material_types(4).sigma_m = 1e10;
material_types(4).color   = [0.2549 0.88235 0.5098];

