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
waveforms.gaussian(1).number_of_cells_per_wavelength = 0; 
waveforms.gaussian(2).number_of_cells_per_wavelength = 15; 
waveforms.derivative_gaussian(1).number_of_cells_per_wavelength = 20; 
waveforms.cosine_modulated_gaussian(1).bandwidth = 1e9; 
waveforms.cosine_modulated_gaussian(1).modulation_frequency = 2e9; 

% voltage sources 
% direction: 'xp', 'xn', 'yp', 'yn', 'zp', or 'zn'
% resistance : ohms, magitude   : volts
voltage_sources(1).min_x = 0;
voltage_sources(1).min_y = 0;
voltage_sources(1).min_z = 0;
voltage_sources(1).max_x = 1.0e-3;
voltage_sources(1).max_y = 2.0e-3;
voltage_sources(1).max_z = 4.0e-3;
voltage_sources(1).direction = 'zp';
voltage_sources(1).resistance = 50;
voltage_sources(1).magnitude = 1;
voltage_sources(1).waveform_type = 'gaussian';
voltage_sources(1).waveform_index = 2;
