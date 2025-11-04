% update electric field source for TE and TEM mode PBC

if strcmp(periodic_boundary.mode,'TE')
    ks = periodic_boundary.source_ks;

    Exi(:,:) = periodic_boundary.Exi0 ...
        * periodic_boundary.waveform(time_step);
    Eyi(:,:) = periodic_boundary.Eyi0 ...
        * periodic_boundary.waveform(time_step);

    for m = 1:nx
        Eyi(m,:) = Eyi(m,:) * exp(-j*kx*(m-1)*dx);
        Exi(m,:) = Exi(m,:) * exp(-j*kx*(m-0.5)*dx);
    end
    Eyi(nxp1,:) = Eyi(nxp1,:) * exp(-j*kx*nx*dx);

    for m = 1:ny
        Eyi(:,m) = Eyi(:,m) * exp(-j*ky*(m-0.5)*dy);
        Exi(:,m) = Exi(:,m) * exp(-j*ky*(m-1)*dy);
    end
    Exi(:,nyp1) = Exi(:,nyp1) * exp(-j*ky*ny*dy);

    Ex(:,:,ks) = Ex(:,:,ks) +  Exi;
    Ey(:,:,ks) = Ey(:,:,ks) +  Eyi;

end   
if strcmp(periodic_boundary.mode,'TEM')
    ks = periodic_boundary.source_ks;

    Exi(:,:) = periodic_boundary.Exi0 ...
        * periodic_boundary.waveform(time_step);
    Eyi(:,:) = periodic_boundary.Eyi0 ...
        * periodic_boundary.waveform(time_step);

    Ex(:,:,ks) = Ex(:,:,ks) +  Exi;
    Ey(:,:,ks) = Ey(:,:,ks) +  Eyi;

end
