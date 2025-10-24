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

% frequency domain parameters
frequency_domain.start = 2e6;
frequency_domain.end   = 2e9;
frequency_domain.step  = 1e6;

% define sampled electric fields
% component: vector component ’x ’,’ y ’,’ z ’, or magnitude ’m’
% display plot = true, in order to plot field during simulation
sampled_electric_fields(1).x = 0;
sampled_electric_fields(1).y = 0;
sampled_electric_fields(1).z = -0.025;
sampled_electric_fields(1).component = 'x';
sampled_electric_fields(1).display_plot = false;

sampled_electric_fields(2).x = 0;
sampled_electric_fields(2).y = 0;
sampled_electric_fields(2).z = 0.225;
sampled_electric_fields(2).component = 'x';
sampled_electric_fields(2).display_plot = false;

animation(1).field_type = 'e';
animation(1).component = 'x';
animation(1).plane_cut(1).type = 'zx';
animation(1).plane_cut(1).position  = 0;
animation(1).enable = false;
animation(1).display_grid = false;
animation(1).display_objects = true;

% display problem space parameters
problem_space_display.labels = true;
problem_space_display.axis_at_origin = false;
problem_space_display.axis_outside_domain = true;
problem_space_display.grid_xn = false;
problem_space_display.grid_xp = false;
problem_space_display.grid_yn = false;
problem_space_display.grid_yp = false;
problem_space_display.grid_zn = false;
problem_space_display.grid_zp = false;
problem_space_display.outer_boundaries = true;
problem_space_display.cpml_boundaries = true;
