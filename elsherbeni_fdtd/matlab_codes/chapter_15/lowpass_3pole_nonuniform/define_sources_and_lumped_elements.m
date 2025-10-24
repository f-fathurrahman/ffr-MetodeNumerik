disp('defining sources and lumped element components'); 

% voltage sources
% direction: 'xp', 'xn', 'yp', 'yn', 'zp', or 'zn' 
% resistance : ohms, magitude   : volts 

% a voltage source
voltage_sources(1).min_x = 0;
voltage_sources(1).min_y = 7.8e-3;
voltage_sources(1).min_z = 0;
voltage_sources(1).max_x = 0.4e-3;
voltage_sources(1).max_y = 9e-3;
voltage_sources(1).max_z = 0;
voltage_sources(1).direction = 'xp';
voltage_sources(1).resistance = 50;
voltage_sources(1).magnitude = 1;
voltage_sources(1).waveform_type = 'gaussian';
voltage_sources(1).number_of_cells_per_wavelength = 20;

% a resistor
resistors(1).min_x = 29.6e-3;
resistors(1).min_y = 7.8e-3;
resistors(1).min_z = 0e-3;
resistors(1).max_x = 30e-3;
resistors(1).max_y = 9e-3;
resistors(1).max_z = 0;
resistors(1).direction = 'x';
resistors(1).resistance = 50;





