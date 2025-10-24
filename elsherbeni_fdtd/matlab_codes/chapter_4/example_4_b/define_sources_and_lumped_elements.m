disp('defining sources and lumped element components');

voltage_sources = [];
current_sources = [];
diodes = [];
resistors = [];
inductors = [];
capacitors = [];

% define source waveform types and parameters
waveforms.sinusoidal(1).frequency = 1e9; 
waveforms.sinusoidal(2).frequency = 5e8; 
waveforms.unit_step(1).start_time_step = 50; 

% voltage sources 
% direction: 'xp', 'xn', 'yp', 'yn', 'zp', or 'zn'
% resistance : ohms, magitude   : volts
voltage_sources(1).min_x = 0;
voltage_sources(1).min_y = 0;
voltage_sources(1).min_z = 0;
voltage_sources(1).max_x = 0.0e-3;
voltage_sources(1).max_y = 2.0e-3;
voltage_sources(1).max_z = 1.0e-3;
voltage_sources(1).direction = 'zp';
voltage_sources(1).resistance = 50;
voltage_sources(1).magnitude = 10;
voltage_sources(1).waveform_type = 'sinusoidal';
voltage_sources(1).waveform_index = 2;
 
% diodes
% direction: 'xp', 'xn', 'yp', 'yn', 'zp', or 'zn'
diodes(1).min_x = 2.0e-3;
diodes(1).min_y = 1.0e-3;
diodes(1).min_z = 0.0e-3;
diodes(1).max_x = 2.0e-3;
diodes(1).max_y = 1.0e-3;
diodes(1).max_z = 1.0e-3;
diodes(1).direction = 'zn';
