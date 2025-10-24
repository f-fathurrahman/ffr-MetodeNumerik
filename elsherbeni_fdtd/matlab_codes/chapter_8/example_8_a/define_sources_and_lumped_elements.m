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
voltage_sources(1).min_x = 14*dx;
voltage_sources(1).min_y = 0;
voltage_sources(1).min_z = 0;
voltage_sources(1).max_x = 20*dx;
voltage_sources(1).max_y = 0;
voltage_sources(1).max_z = 3*dz;
voltage_sources(1).direction = 'zp';
voltage_sources(1).resistance = 50;
voltage_sources(1).magnitude = 1;
voltage_sources(1).waveform_type = 'gaussian';
voltage_sources(1).waveform_index = 1;

% resistors
% direction: 'x', 'y', or 'z'
% resistance : ohms
resistors(1).min_x = 30*dx;
resistors(1).min_y = 46*dy;
resistors(1).min_z = 0;
resistors(1).max_x = 36*dx;
resistors(1).max_y = 46*dy;
resistors(1).max_z = 3*dz;
resistors(1).direction = 'z';
resistors(1).resistance = 50;
