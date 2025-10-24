disp('defining output parameters');

sampled_electric_fields = [];
sampled_magnetic_fields = [];
sampled_transient_E_planes = [];
sampled_frequency_E_planes = [];

% figure refresh rate
plotting_step = 10;

% frequency domain parameters
frequency_domain.start = 20e6;
frequency_domain.end   = 20e9;
frequency_domain.step  = 20e6;

% define sampled magnetic fields
sampled_magnetic_fields(1).x = 8e-3;
sampled_magnetic_fields(1).y = 8e-3;
sampled_magnetic_fields(1).component = 'z';
sampled_magnetic_fields(1).display_plot = false;

draw_domain.enable = true;
draw_domain.draw_axis = true;
draw_domain.draw_grid = true;
