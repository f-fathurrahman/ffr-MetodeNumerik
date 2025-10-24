disp('defining sources and lumped element components');

voltage_sources = [];
current_sources = [];
diodes = [];
resistors = [];
inductors = [];
capacitors = [];

% define source waveform types and parameters
waveforms.gaussian(1).number_of_cells_per_wavelength = 0; 
waveforms.gaussian(2).number_of_cells_per_wavelength = 15; 

% voltage sources 
% direction: 'xp', 'xn', 'yp', 'yn', 'zp', or 'zn'
% resistance : ohms, magitude   : volts
voltage_sources(1).min_x = 24*dx;
voltage_sources(1).min_y = 0;
voltage_sources(1).min_z = 0;
voltage_sources(1).max_x = 36*dx;
voltage_sources(1).max_y = 0;
voltage_sources(1).max_z = 6*dz;
voltage_sources(1).direction = 'zp';
voltage_sources(1).resistance = 50;
voltage_sources(1).magnitude = 1;
voltage_sources(1).waveform_type = 'gaussian';
voltage_sources(1).waveform_index = 1;
