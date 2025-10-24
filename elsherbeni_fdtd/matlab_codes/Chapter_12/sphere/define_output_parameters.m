disp('defining output parameters');

sampled_electric_fields = [];
sampled_magnetic_fields = [];
sampled_voltages = [];
sampled_currents = [];
ports = [];
farfield.frequencies = [];

% figure refresh rate
plotting_step = 10;

% mode of operation
run_simulation = true;
show_material_mesh = true;
show_problem_space = true;

% far field calculation parameters
farfield.frequencies(1) = 1.0e9;
farfield.number_of_cells_from_outer_boundary = 13;

% frequency domain parameters
frequency_domain.start = 20e6;
frequency_domain.end   = 4e9;
frequency_domain.step  = 20e6;

% define sampled electric fields
% component: vector component ’x ’,’ y ’,’ z ’, or magnitude ’m’
% display plot = true, in order to plot field during simulation
sampled_electric_fields(1).x = 0;
sampled_electric_fields(1).y = 0;
sampled_electric_fields(1).z = 0;
sampled_electric_fields(1).component = 'x';
sampled_electric_fields(1).display_plot = false;
 
% define animation
% field_type shall be 'e' or 'h'
% plane cut shall be 'xy', yz, or zx
% component shall be 'x', 'y', 'z', or 'm;
animation(1).field_type = 'e';
animation(1).component = 'm';
animation(1).plane_cut(1).type = 'zx';
animation(1).plane_cut(1).position  = 0;
animation(1).display_grid = false;
animation(1).display_objects = true;
animation(1).save_movie = true;
animation(1).view_angles = [15 15];
animation(1).zoom = 0.8;
animation(1).enable = true;
