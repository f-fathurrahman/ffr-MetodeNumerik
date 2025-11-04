% update magnetic field source for TM mode PBC
if strcmp(periodic_boundary.mode,'TM')

    ks = periodic_boundary.source_ks;

    Hxi(:,:) = periodic_boundary.Hxi0 ...
        * periodic_boundary.waveform(time_step);
    Hyi(:,:) = periodic_boundary.Hyi0 ...
        * periodic_boundary.waveform(time_step);

    for m = 1:nx
        Hyi(m,:) = Hyi(m,:) * exp(-j*kx*(m-0.5)*dx);
        Hxi(m,:) = Hxi(m,:) * exp(-j*kx*(m-1)*dx);
    end
    Hxi(nxp1,:) = Hxi(nxp1,:)*exp(-j*kx*nx*dx);

    for m = 1:ny
        Hyi(:,m) = Hyi(:,m) * exp(-j*ky*(m-1)*dy);
        Hxi(:,m) = Hxi(:,m) * exp(-j*ky*(m-0.5)*dy);
    end
    Hyi(:,nyp1) = Hyi(:,nyp1) * exp(-j*ky*ny*dy);

    Hx(:,:,ks) = Hx(:,:,ks) +  Hxi;
    Hy(:,:,ks) = Hy(:,:,ks) +  Hyi;

end    