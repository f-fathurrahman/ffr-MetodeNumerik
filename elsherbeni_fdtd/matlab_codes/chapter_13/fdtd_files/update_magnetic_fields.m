% update magnetic fields

current_time  = current_time + dt/2;
    
Hx = Chxh.*Hx+Chxey.*(Ey(1:nxp1,1:ny,2:nzp1)-Ey(1:nxp1,1:ny,1:nz)) ...
    + Chxez.*(Ez(1:nxp1,2:nyp1,1:nz)-Ez(1:nxp1,1:ny,1:nz)); 
                  
Hy = Chyh.*Hy+Chyez.*(Ez(2:nxp1,1:nyp1,1:nz)-Ez(1:nx,1:nyp1,1:nz)) ...
    + Chyex.*(Ex(1:nx,1:nyp1,2:nzp1)-Ex(1:nx,1:nyp1,1:nz)); 
               
Hz = Chzh.*Hz+Chzex.*(Ex(1:nx,2:nyp1,1:nzp1)-Ex(1:nx,1:ny,1:nzp1))  ...
    + Chzey.*(Ey(2:nxp1,1:ny,1:nzp1)-Ey(1:nx,1:ny,1:nzp1)); 

if incident_plane_wave.enable
   xui = incident_plane_wave.ui;
   xuj = incident_plane_wave.uj;
   xuk = incident_plane_wave.uk;
   xli = incident_plane_wave.li;
   xlj = incident_plane_wave.lj;
   xlk = incident_plane_wave.lk;
    % Hx_yn
        Hx(xli:xui,xlj-1,xlk:xuk-1) = ...
        Hx(xli:xui,xlj-1,xlk:xuk-1) + dt/(mu_0*dy) * Ezi_yn(:,1,:); 
    % Hx_yp
        Hx(xli:xui,xuj,xlk:xuk-1) = ...
        Hx(xli:xui,xuj,xlk:xuk-1) - dt/(mu_0*dy) * Ezi_yp(:,1,:); 
    % Hx_zn
        Hx(xli:xui,xlj:xuj-1,xlk-1) = ...
        Hx(xli:xui,xlj:xuj-1,xlk-1) - dt/(mu_0*dz) * Eyi_zn(:,:,1); 
    % Hx_zp
        Hx(xli:xui,xlj:xuj-1,xuk) = ...
        Hx(xli:xui,xlj:xuj-1,xuk) + dt/(mu_0*dz) * Eyi_zp(:,:,1); 

    % Hy_xn
        Hy(xli-1,xlj:xuj,xlk:xuk-1) = ...
        Hy(xli-1,xlj:xuj,xlk:xuk-1) - dt/(mu_0*dx) * Ezi_xn(1,:,:);
    % Hy_xp
        Hy(xui,xlj:xuj,xlk:xuk-1) = ...
        Hy(xui,xlj:xuj,xlk:xuk-1) + dt/(mu_0*dx) * Ezi_xp(1,:,:);
    % Hy_zn
        Hy(xli:xui-1,xlj:xuj,xlk-1) = ...
        Hy(xli:xui-1,xlj:xuj,xlk-1) + dt/(mu_0*dz) * Exi_zn(:,:,1); 
    % Hy_zp
        Hy(xli:xui-1,xlj:xuj,xuk) = ...
        Hy(xli:xui-1,xlj:xuj,xuk) - dt/(mu_0*dz) * Exi_zp(:,:,1); 

    % Hz_yn
        Hz(xli:xui-1, xlj-1,xlk:xuk) = ...
        Hz(xli:xui-1, xlj-1,xlk:xuk) - dt/(mu_0*dy) * Exi_yn(:,1,:);
    % Hz_yp
        Hz(xli:xui-1, xuj,xlk:xuk) = ...
        Hz(xli:xui-1, xuj,xlk:xuk) + dt/(mu_0*dy) * Exi_yp(:,1,:);
    % Hz_xn
        Hz(xli-1,xlj:xuj-1,xlk:xuk) = ...
        Hz(xli-1,xlj:xuj-1,xlk:xuk) + dt/(mu_0*dx) * Eyi_xn(1,:,:); 
    % Hz_xp
        Hz(xui,xlj:xuj-1,xlk:xuk) = ...
        Hz(xui,xlj:xuj-1,xlk:xuk) - dt/(mu_0*dx) * Eyi_xp(1,:,:); 
end
