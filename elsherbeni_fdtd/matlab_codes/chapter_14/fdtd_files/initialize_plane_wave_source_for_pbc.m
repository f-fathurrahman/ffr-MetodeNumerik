
% Here the incident plane wave is assumed to propagate in the -z direction
kx = periodic_boundary.kx;
ky = periodic_boundary.ky;
periodic_boundary.kxy = sqrt(kx^2+ky^2);
phi_incident = atan2(ky,kx);

if strcmp(periodic_boundary.mode,'TE')
    phi_incident = atan2(ky,kx);
    E_phi = periodic_boundary.E_phi;
    periodic_boundary.Exi0 = -E_phi * sin(phi_incident);
    periodic_boundary.Eyi0 = E_phi * cos(phi_incident);
    Exi = zeros(nx, nyp1); 
    Eyi = zeros(nxp1, ny); 
end
if strcmp(periodic_boundary.mode,'TM')
    phi_incident = atan2(ky,kx);
    H_phi = periodic_boundary.H_phi;
    periodic_boundary.Hxi0 = -H_phi * sin(phi_incident);
    periodic_boundary.Hyi0 =  H_phi * cos(phi_incident);
    Hxi = zeros(nxp1, ny); 
    Hyi = zeros(nx, nyp1); 
end
if strcmp(periodic_boundary.mode,'TEM')
    phi_incident = 0;
    periodic_boundary.Exi0 = periodic_boundary.E_x;
    periodic_boundary.Eyi0 = periodic_boundary.E_y;
    Exi = zeros(nx, nyp1); 
    Eyi = zeros(nxp1, ny); 
end

periodic_boundary.phi_incident = phi_incident;
f_min = c*periodic_boundary.kxy/(2*pi);
f_max = c/(40*max([dx,dy,dz]));%40 cells per wavelength
bandwidth = f_max-f_min;
frequency = (f_max+f_min)/2;
tau = (2*4.29/pi)./bandwidth;  %for edge of Gaussian 40 dB below max
t_0 = 4.5 * tau;
periodic_boundary.frequency_start = f_min;
periodic_boundary.frequency_end = f_max;
periodic_boundary.modulation_frequency = frequency;
periodic_boundary.bandwidth = bandwidth;
periodic_boundary.tau = tau;
periodic_boundary.t_0 = t_0;
periodic_boundary.waveform = ...
    cos(2*pi*frequency*(time - t_0)).*exp(-((time - t_0)/tau).^2);
