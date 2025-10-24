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

% define sampled voltages
sampled_voltages(1).min_x = 2.0e-3;
sampled_voltages(1).min_y = 1.0e-3;
sampled_voltages(1).min_z = 0;
sampled_voltages(1).max_x = 2.0e-3;
sampled_voltages(1).max_y = 1.0e-3;
sampled_voltages(1).max_z = 1.0e-3;
sampled_voltages(1).direction = 'zp';
sampled_voltages(1).display_plot = true;