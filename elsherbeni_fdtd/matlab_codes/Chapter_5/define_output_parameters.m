disp('defining output parameters');

sampled_electric_fields = [];
sampled_magnetic_fields = [];
sampled_voltages = [];
sampled_currents = [];

% figure refresh rate
plotting_step = 10;

% mode of operation
run_simulation = true;
show_material_mesh = true;
show_problem_space = true;

% frequency domain parameters
frequency_domain.start = 1e7;
frequency_domain.end   = 1e9;
frequency_domain.step  = 1e7;

% define sampled electric fields
% component: vector component 'x','y','z', or magnitude 'm'
% display_plot = true, in order to plot field during simulation 
sampled_electric_fields(1).x = 30*dx;
sampled_electric_fields(1).y = 30*dy;
sampled_electric_fields(1).z = 10*dz;
sampled_electric_fields(1).component = 'x';
sampled_electric_fields(1).display_plot = true;

% define sampled magnetic fields
% component: vector component 'x','y','z', or magnitude 'm'
% display_plot = true, in order to plot field during simulation 
sampled_magnetic_fields(1).x = 30*dx;
sampled_magnetic_fields(1).y = 30*dy;
sampled_magnetic_fields(1).z = 10*dz;
sampled_magnetic_fields(1).component = 'm';
sampled_magnetic_fields(1).display_plot = true;

% define sampled voltages
sampled_voltages(1).min_x = 5.0e-3;
sampled_voltages(1).min_y = 0;
sampled_voltages(1).min_z = 0;
sampled_voltages(1).max_x = 5.0e-3;
sampled_voltages(1).max_y = 2.0e-3;
sampled_voltages(1).max_z = 4.0e-3;
sampled_voltages(1).direction = 'zp';
sampled_voltages(1).display_plot = true;

% define sampled currents
sampled_currents(1).min_x = 5.0e-3;
sampled_currents(1).min_y = 0;
sampled_currents(1).min_z = 4.0e-3;
sampled_currents(1).max_x = 5.0e-3;
sampled_currents(1).max_y = 2.0e-3;
sampled_currents(1).max_z = 4.0e-3;
sampled_currents(1).direction = 'xp';
sampled_currents(1).display_plot = true;
