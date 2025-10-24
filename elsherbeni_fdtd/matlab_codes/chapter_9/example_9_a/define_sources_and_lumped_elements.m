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
voltage_sources(1).min_x = -0.787e-3;
voltage_sources(1).min_y = 0;
voltage_sources(1).min_z =  24e-3;
voltage_sources(1).max_x = 0;
voltage_sources(1).max_y = 0;
voltage_sources(1).max_z =  26.4e-3;
voltage_sources(1).direction = 'xp';
voltage_sources(1).resistance = 50;
voltage_sources(1).magnitude = 1;
voltage_sources(1).waveform_type = 'gaussian';
voltage_sources(1).waveform_index = 1;
