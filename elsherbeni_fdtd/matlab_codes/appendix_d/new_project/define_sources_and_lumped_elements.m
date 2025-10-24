disp('defining sources and lumped element components'); 

% voltage sources
% direction: 'xp', 'xn', 'yp', 'yn', 'zp', or 'zn' 
% resistance : ohms, magitude   : volts 

% a voltage source
voltage_sources(1).min_x = 0.1;
voltage_sources(1).min_y = 0.1;
voltage_sources(1).min_z = 0;
voltage_sources(1).max_x = 0.2;
voltage_sources(1).max_y = 0.2;
voltage_sources(1).max_z = 0.1;
voltage_sources(1).direction = 'zp';
voltage_sources(1).resistance = 50;
voltage_sources(1).magnitude = 1;
voltage_sources(1).waveform_type = 'gaussian';
voltage_sources(1).number_of_cells_per_wavelength = 20;

% a voltage source
voltage_sources(2).min_x = 0.2;
voltage_sources(2).min_y = 0.2;
voltage_sources(2).min_z = 0;
voltage_sources(2).max_x = 0.3;
voltage_sources(2).max_y = 0.3;
voltage_sources(2).max_z = 0.1;
voltage_sources(2).direction = 'zp';
voltage_sources(2).resistance = 50;
voltage_sources(2).magnitude = 1;
voltage_sources(2).waveform_type = 'gaussian';
voltage_sources(2).number_of_cells_per_wavelength = 20;



% Define incident plane wave, angles are in degrees
incident_plane_wave.E_theta = 1;
incident_plane_wave.E_phi = 0;
incident_plane_wave.theta_incident = 0;
incident_plane_wave.phi_incident = 0;
incident_plane_wave.waveform_type = 'gaussian';
incident_plane_wave.number_of_cells_per_wavelength = 20;
incident_plane_wave.enable = false;

% resistors
% direction: 'x', 'y', or 'z' 
% resistance : ohms

% a resistor
resistors(1).min_x = 0.1;
resistors(1).min_y = 0.1;
resistors(1).min_z = 0;
resistors(1).max_x = 0.2;
resistors(1).max_y = 0.2;
resistors(1).max_z = 0.1;
resistors(1).direction = 'z';
resistors(1).resistance = 50;

% a resistor
resistors(2).min_x = 0.2;
resistors(2).min_y = 0.2;
resistors(2).min_z = 0;
resistors(2).max_x = 0.3;
resistors(2).max_y = 0.3;
resistors(2).max_z = 0.1;
resistors(2).direction = 'z';
resistors(2).resistance = 50;





