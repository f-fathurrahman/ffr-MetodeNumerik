% set simulation mode as PBC if needed
% all four sides (xn, xp, yn, yp) should be PBC  
if (strcmp(boundary.type_xn, 'pbc') && strcmp(boundary.type_xp, 'pbc')...
    && strcmp(boundary.type_xn, 'pbc') && strcmp(boundary.type_xp, 'pbc'))
    periodic_boundary_simulation = true;
else
    periodic_boundary_simulation = false;
    return;
end
    
periodic_boundary.source_ks = ...
    1+round((periodic_boundary.source_z - ...
    fdtd_domain.min_z)/fdtd_domain.dz);

if isfield(periodic_boundary,'reflection_z')
    periodic_boundary.calculate_reflection = true;
    periodic_boundary.reflection_ks = ...
    1+round((periodic_boundary.reflection_z - ...
    fdtd_domain.min_z)/fdtd_domain.dz);
    periodic_boundary.reflection_ex = zeros(1,number_of_time_steps);
    periodic_boundary.reflection_ey = zeros(1,number_of_time_steps);
    periodic_boundary.reflection_ez = zeros(1,number_of_time_steps);
    periodic_boundary.reflection_hx = zeros(1,number_of_time_steps);
    periodic_boundary.reflection_hy = zeros(1,number_of_time_steps);
    periodic_boundary.reflection_hz = zeros(1,number_of_time_steps);
else
        periodic_boundary.calculate_reflection = false;
end

if isfield(periodic_boundary,'transmission_z')
    periodic_boundary.calculate_transmission = true;
    periodic_boundary.transmission_ks = ...
    1+round((periodic_boundary.transmission_z - ...
    fdtd_domain.min_z)/fdtd_domain.dz);
    periodic_boundary.transmission_ex = zeros(1,number_of_time_steps);
    periodic_boundary.transmission_ey = zeros(1,number_of_time_steps);
    periodic_boundary.transmission_ez = zeros(1,number_of_time_steps);
    periodic_boundary.transmission_hx = zeros(1,number_of_time_steps);
    periodic_boundary.transmission_hy = zeros(1,number_of_time_steps);
    periodic_boundary.transmission_hz = zeros(1,number_of_time_steps);
    periodic_boundary.reflection_transmission_distance = ...
        dz*(periodic_boundary.reflection_ks - ...
        periodic_boundary.transmission_ks);
else
    periodic_boundary.calculate_transmission = false;
end

if strcmp(periodic_boundary.mode,'TEM')
    periodic_boundary.kx = 0;
    periodic_boundary.ky = 0;
end

kx = periodic_boundary.kx;
ky = periodic_boundary.ky;

pex = zeros(nx,ny);
pey = zeros(nx,ny);
pez = zeros(nx,ny);
phz = zeros(nx,ny);

for mi=1:nx
    for mj=1:ny
        pex(mi,mj) = exp(j*kx*(mi-0.5)*dx)*exp(j*ky*(mj-1)*dy);
        pey(mi,mj) = exp(j*kx*(mi-1)*dx)*exp(j*ky*(mj-0.5)*dy);
        pez(mi,mj) = exp(j*kx*(mi-1)*dx)*exp(j*ky*(mj-1)*dy);
        phz(mi,mj) = exp(j*kx*(mi-0.5)*dx)*exp(j*ky*(mj-0.5)*dy);
    end
end

periodic_boundary.phase_correction_ex = pex;
periodic_boundary.phase_correction_ey = pey;
periodic_boundary.phase_correction_ez = pez;

periodic_boundary.phase_correction_hx = pey;
periodic_boundary.phase_correction_hy = pex;
periodic_boundary.phase_correction_hz = phz;
    
periodic_boundary.Px = nx*dx;
periodic_boundary.Py = ny*dy;

initialize_plane_wave_source_for_pbc;





        
