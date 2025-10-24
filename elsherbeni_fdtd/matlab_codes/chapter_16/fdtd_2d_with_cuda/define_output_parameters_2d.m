disp('defining output parameters');

sampled_electric_fields = [];
sampled_magnetic_fields = [];
sampled_transient_E_planes = [];
sampled_frequency_E_planes = [];
sampled_transient_H_planes = [];

% figure refresh rate
plotting_step = 1;

% frequency domain parameters
frequency_domain.start = 20e6;
frequency_domain.end   = 20e9;
frequency_domain.step  = 20e6;

% % define sampled electric fields
sampled_electric_fields(1).x = 0.7;
sampled_electric_fields(1).y = 0.5;
sampled_electric_fields(1).component = 'z';
sampled_electric_fields(1).display_plot = false;

sampled_transient_E_planes(1).component = 'z';

draw_domain.enable = false;
draw_domain.draw_axis = false;
draw_domain.draw_grid = false;

