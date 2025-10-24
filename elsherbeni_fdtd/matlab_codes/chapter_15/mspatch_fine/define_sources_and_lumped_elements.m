disp('defining sources and lumped element components'); 

% voltage sources
% direction: 'xp', 'xn', 'yp', 'yn', 'zp', or 'zn' 
% resistance : ohms, magitude   : volts 

% a voltage source
voltage_sources(1).min_x = 30e-3;
voltage_sources(1).min_y = 10e-3;
voltage_sources(1).min_z = 0;
voltage_sources(1).max_x = 30e-3;
voltage_sources(1).max_y = 10e-3;
voltage_sources(1).max_z = 1.9e-3;
voltage_sources(1).direction = 'zp';
voltage_sources(1).resistance = 50;
voltage_sources(1).magnitude = 1;
voltage_sources(1).waveform_type = 'gaussian';
voltage_sources(1).number_of_cells_per_wavelength = 20;





