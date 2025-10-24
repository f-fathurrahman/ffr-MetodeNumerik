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
voltage_sources(1).min_x = -15*dx;
voltage_sources(1).min_y = -25*dy;
voltage_sources(1).min_z = 0;
voltage_sources(1).max_x = -9*dx;
voltage_sources(1).max_y = -25*dy;
voltage_sources(1).max_z = 3*dz;
voltage_sources(1).direction = 'zp';
voltage_sources(1).resistance = 50;
voltage_sources(1).magnitude = 1;
voltage_sources(1).waveform_type = 'gaussian';
voltage_sources(1).waveform_index = 1;

% resistors
% direction: 'x', 'y', or 'z'
% resistance : ohms
resistors(1).min_x = 9*dx;
resistors(1).min_y = -25*dy;
resistors(1).min_z = 0;
resistors(1).max_x = 15*dx;
resistors(1).max_y = -25*dy;
resistors(1).max_z = 3*dz;
resistors(1).direction = 'z';
resistors(1).resistance = 50;

% resistors
% direction: 'x', 'y', or 'z'
% resistance : ohms
resistors(2).min_x = 9*dx;
resistors(2).min_y = 25*dy;
resistors(2).min_z = 0;
resistors(2).max_x = 15*dx;
resistors(2).max_y = 25*dy;
resistors(2).max_z = 3*dz;
resistors(2).direction = 'z';
resistors(2).resistance = 50;

% resistors
% direction: 'x', 'y', or 'z'
% resistance : ohms
resistors(3).min_x = -15*dx;
resistors(3).min_y = 25*dy;
resistors(3).min_z = 0;
resistors(3).max_x = -9*dx;
resistors(3).max_y = 25*dy;
resistors(3).max_z = 3*dz;
resistors(3).direction = 'z';
resistors(3).resistance = 50;

