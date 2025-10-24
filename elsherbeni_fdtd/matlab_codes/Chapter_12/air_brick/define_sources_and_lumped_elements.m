disp('defining sources and lumped element components'); 

% Define incident plane wave, angles are in degrees
incident_plane_wave.E_theta = 1;
incident_plane_wave.E_phi = 0;
incident_plane_wave.theta_incident = 30;
incident_plane_wave.phi_incident = 90;
incident_plane_wave.waveform_type = 'gaussian';
incident_plane_wave.number_of_cells_per_wavelength = 20;
incident_plane_wave.enable = true;






