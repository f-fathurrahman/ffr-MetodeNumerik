% Calculate total radiated power
if incident_plane_wave.enabled == false
    return;
end
x = incident_plane_wave.waveform;    
time_shift = 0;
[X] = time_to_frequency_domain(x,dt,farfield.frequencies,time_shift);
incident_plane_wave.frequency_domain_value = X;
incident_plane_wave.frequencies = frequency_array;
E_t = incident_plane_wave.E_theta;
E_p = incident_plane_wave.E_phi;
eta_0 = sqrt(mu_0/eps_0);
incident_plane_wave.incident_power = ...
    (0.5/eta_0) * (E_t^2 + E_p^2) * abs(X)^2; 
