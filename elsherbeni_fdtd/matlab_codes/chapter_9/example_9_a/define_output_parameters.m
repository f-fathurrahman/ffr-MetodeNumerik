disp('defining output parameters');

sampled_electric_fields = [];
sampled_magnetic_fields = [];
sampled_voltages = [];
sampled_currents = [];
ports = [];
farfield.frequencies = [];

% figure refresh rate
plotting_step = 100;

% mode of operation
run_simulation = true;
show_material_mesh = true;
show_problem_space = true;

% far field calculation parameters
farfield.frequencies(1) = 2.4e9;
farfield.frequencies(2) = 5.8e9;
farfield.number_of_cells_from_outer_boundary = 13;

% frequency domain parameters
frequency_domain.start = 20e6;
frequency_domain.end   = 10e9;
frequency_domain.step  = 20e6;

% define sampled voltages
sampled_voltages(1).min_x = -0.787e-3;
sampled_voltages(1).min_y = 0;
sampled_voltages(1).min_z = 24.4e-3;
sampled_voltages(1).max_x = 0;
sampled_voltages(1).max_y = 0;
sampled_voltages(1).max_z = 26.4e-3;
sampled_voltages(1).direction = 'xp';
sampled_voltages(1).display_plot = false;

% define sampled currents
sampled_currents(1).min_x = -0.39e-3;
sampled_currents(1).min_y = 0;
sampled_currents(1).min_z = 24e-3;
sampled_currents(1).max_x = -0.39e-3;
sampled_currents(1).max_y = 0;
sampled_currents(1).max_z = 26.4e-3;
sampled_currents(1).direction = 'xp';
sampled_currents(1).display_plot = false;

% define ports
ports(1).sampled_voltage_index = 1;
ports(1).sampled_current_index = 1;
ports(1).impedance = 50;
ports(1).is_source_port = true;

