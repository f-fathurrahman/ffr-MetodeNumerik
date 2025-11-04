
nxy = nx*ny;

if periodic_boundary.calculate_reflection

    ks = periodic_boundary.reflection_ks;   
    
    ex = Ex(1:nx,1:ny,ks).*periodic_boundary.phase_correction_ex;
    ey = Ey(1:nx,1:ny,ks).*periodic_boundary.phase_correction_ey;
    ez = 0.5*(Ez(1:nx,1:ny,ks)+Ez(1:nx,1:ny,ks-1)) ...
            .*periodic_boundary.phase_correction_ez;

    hx = 0.5*(Hx(1:nx,1:ny,ks)+Hx(1:nx,1:ny,ks-1)) ...
        .*periodic_boundary.phase_correction_hx;
    hy = 0.5*(Hy(1:nx,1:ny,ks)+Hy(1:nx,1:ny,ks-1)) ...
        .*periodic_boundary.phase_correction_hy;
    hz = Hz(1:nx,1:ny,ks).*periodic_boundary.phase_correction_hz;

    periodic_boundary.reflection_ex(time_step) = sum(sum(ex))/nxy;
    periodic_boundary.reflection_ey(time_step) = sum(sum(ey))/nxy;
    periodic_boundary.reflection_ez(time_step) = sum(sum(ez))/nxy;

    periodic_boundary.reflection_hx(time_step) = sum(sum(hx))/nxy;
    periodic_boundary.reflection_hy(time_step) = sum(sum(hy))/nxy;
    periodic_boundary.reflection_hz(time_step) = sum(sum(hz))/nxy;

end

if periodic_boundary.calculate_transmission
    ks = periodic_boundary.transmission_ks;   
   
    ex = Ex(1:nx,1:ny,ks).*periodic_boundary.phase_correction_ex;
    ey = Ey(1:nx,1:ny,ks).*periodic_boundary.phase_correction_ey;
    ez = 0.5*(Ez(1:nx,1:ny,ks)+Ez(1:nx,1:ny,ks-1)) ...
            .*periodic_boundary.phase_correction_ez;

    hx = 0.5*(Hx(1:nx,1:ny,ks)+Hx(1:nx,1:ny,ks-1)) ...
        .*periodic_boundary.phase_correction_hx;
    hy = 0.5*(Hy(1:nx,1:ny,ks)+Hy(1:nx,1:ny,ks-1)) ...
        .*periodic_boundary.phase_correction_hy;
    hz = Hz(1:nx,1:ny,ks).*periodic_boundary.phase_correction_hz;

    periodic_boundary.transmission_ex(time_step) = sum(sum(ex))/nxy;
    periodic_boundary.transmission_ey(time_step) = sum(sum(ey))/nxy;
    periodic_boundary.transmission_ez(time_step) = sum(sum(ez))/nxy;

    periodic_boundary.transmission_hx(time_step) = sum(sum(hx))/nxy;
    periodic_boundary.transmission_hy(time_step) = sum(sum(hy))/nxy;
    periodic_boundary.transmission_hz(time_step) = sum(sum(hz))/nxy;   
end

