% initializing the incident plane wave

if ~incident_plane_wave.enable
    return;
end

% leaving three cells between the outer boundary and the TF/SF boundary
li = 3;     lj = 3;     lk = 3; 
ui = nx-3;  uj = ny-3;  uk = nz-3;
if strcmp(boundary.type_xn, 'cpml')
    li = boundary.cpml_number_of_cells_xn + 3; 
end
if strcmp(boundary.type_yn, 'cpml')
    lj = boundary.cpml_number_of_cells_yn + 3; 
end
if strcmp(boundary.type_zn, 'cpml')
    lk = boundary.cpml_number_of_cells_zn + 3; 
end
if strcmp(boundary.type_xp, 'cpml')
    ui = nx-boundary.cpml_number_of_cells_xp - 3;
end
if strcmp(boundary.type_yp, 'cpml')
    uj = ny-boundary.cpml_number_of_cells_yp - 3; 
end
if strcmp(boundary.type_zp, 'cpml')
    uk = nz-boundary.cpml_number_of_cells_zp - 3;
end

Exi_yn = zeros(ui-li,1,uk-lk+1); 
Exi_yp = zeros(ui-li,1,uk-lk+1); 
Exi_zn = zeros(ui-li,uj-lj+1,1); 
Exi_zp = zeros(ui-li,uj-lj+1,1); 

Eyi_xn = zeros(1,uj-lj,uk-lk+1); 
Eyi_xp = zeros(1,uj-lj,uk-lk+1); 
Eyi_zn = zeros(ui-li+1,uj-lj,1); 
Eyi_zp = zeros(ui-li+1,uj-lj,1); 

Ezi_xn = zeros(1,uj-lj+1,uk-lk); 
Ezi_xp = zeros(1,uj-lj+1,uk-lk); 
Ezi_yn = zeros(ui-li+1,1,uk-lk); 
Ezi_yp = zeros(ui-li+1,1,uk-lk); 

Hzi_yn = zeros(ui-li,1,uk-lk+1); 
Hzi_yp = zeros(ui-li,1,uk-lk+1); 
Hyi_zn = zeros(ui-li,uj-lj+1,1); 
Hyi_zp = zeros(ui-li,uj-lj+1,1); 

Hzi_xn = zeros(1,uj-lj,uk-lk+1); 
Hzi_xp = zeros(1,uj-lj,uk-lk+1); 
Hxi_zn = zeros(ui-li+1,uj-lj,1); 
Hxi_zp = zeros(ui-li+1,uj-lj,1); 

Hyi_xn = zeros(1,uj-lj+1,uk-lk); 
Hyi_xp = zeros(1,uj-lj+1,uk-lk); 
Hxi_yn = zeros(ui-li+1,1,uk-lk); 
Hxi_yp = zeros(ui-li+1,1,uk-lk); 

k_dot_r_ex_yn = zeros(ui-li,1,uk-lk+1); 
k_dot_r_ex_yp = zeros(ui-li,1,uk-lk+1); 
k_dot_r_ex_zn = zeros(ui-li,uj-lj+1,1); 
k_dot_r_ex_zp = zeros(ui-li,uj-lj+1,1); 

k_dot_r_ey_xn = zeros(1,uj-lj,uk-lk+1); 
k_dot_r_ey_xp = zeros(1,uj-lj,uk-lk+1); 
k_dot_r_ey_zn = zeros(ui-li+1,uj-lj,1); 
k_dot_r_ey_zp = zeros(ui-li+1,uj-lj,1); 

k_dot_r_ez_xn = zeros(1,uj-lj+1,uk-lk); 
k_dot_r_ez_xp = zeros(1,uj-lj+1,uk-lk); 
k_dot_r_ez_yn = zeros(ui-li+1,1,uk-lk); 
k_dot_r_ez_yp = zeros(ui-li+1,1,uk-lk); 

k_dot_r_hz_yn = zeros(ui-li,1,uk-lk+1); 
k_dot_r_hz_yp = zeros(ui-li,1,uk-lk+1); 
k_dot_r_hy_zn = zeros(ui-li,uj-lj+1,1); 
k_dot_r_hy_zp = zeros(ui-li,uj-lj+1,1); 

k_dot_r_hz_xn = zeros(1,uj-lj,uk-lk+1); 
k_dot_r_hz_xp = zeros(1,uj-lj,uk-lk+1); 
k_dot_r_hx_zn = zeros(ui-li+1,uj-lj,1); 
k_dot_r_hx_zp = zeros(ui-li+1,uj-lj,1); 

k_dot_r_hy_xn = zeros(1,uj-lj+1,uk-lk); 
k_dot_r_hy_xp = zeros(1,uj-lj+1,uk-lk); 
k_dot_r_hx_yn = zeros(ui-li+1,1,uk-lk); 
k_dot_r_hx_yp = zeros(ui-li+1,1,uk-lk); 
  
% calculate the amplitude factors for field components
theta_incident = incident_plane_wave.theta_incident*pi/180;
phi_incident = incident_plane_wave.phi_incident*pi/180;
E_theta = incident_plane_wave.E_theta;
E_phi = incident_plane_wave.E_phi;
eta_0 = sqrt(mu_0/eps_0);
Exi0 = E_theta * cos(theta_incident) * cos(phi_incident) ...
    - E_phi * sin(phi_incident);
Eyi0 = E_theta * cos(theta_incident) * sin(phi_incident) ...
    + E_phi * cos(phi_incident);
Ezi0 = -E_theta * sin(theta_incident);
Hxi0 = (-1/eta_0)*(E_phi * cos(theta_incident) ...
    * cos(phi_incident) + E_theta * sin(phi_incident));
Hyi0 = (-1/eta_0)*(E_phi * cos(theta_incident) ...
    * sin(phi_incident) - E_theta * cos(phi_incident));
Hzi0 = (1/eta_0)*(E_phi * sin(theta_incident));

% calculate spatial shift, l_0, required for incident plane wave
r0 =[fdtd_domain.min_x fdtd_domain.min_y fdtd_domain.min_z;
fdtd_domain.min_x fdtd_domain.min_y fdtd_domain.max_z;
fdtd_domain.min_x fdtd_domain.max_y fdtd_domain.min_z;
fdtd_domain.min_x fdtd_domain.max_y fdtd_domain.max_z;
fdtd_domain.max_x fdtd_domain.min_y fdtd_domain.min_z;
fdtd_domain.max_x fdtd_domain.min_y fdtd_domain.max_z;
fdtd_domain.max_x fdtd_domain.max_y fdtd_domain.min_z;
fdtd_domain.max_x fdtd_domain.max_y fdtd_domain.max_z;];

k_vec_x =  sin(theta_incident)*cos(phi_incident);    
k_vec_y =  sin(theta_incident)*sin(phi_incident);    
k_vec_z =  cos(theta_incident);

k_dot_r0 = k_vec_x * r0(:,1) ...
    + k_vec_y * r0(:,2) ...
    + k_vec_z * r0(:,3);

l_0 = min(k_dot_r0)/c;

min_x = fdtd_domain.min_x;
min_y = fdtd_domain.min_y;
min_z = fdtd_domain.min_z;

% calculate k.r for every field component on the TF/SF boundary
for mk=lk:uk
        for mi=li:ui-1
        mii = mi - li + 1;
        mjj = 1;
        mkk = mk - lk + 1;

        mj = lj;
        x = min_x + dx * (mi-0.5);
        y = min_y + dy * (mj-1);
        z = min_z + dz * (mk-1);
        
        k_dot_r_ex_yn(mii,mjj,mkk) = ...
            (x*k_vec_x + y*k_vec_y + z*k_vec_z)/c;

        y = min_y + dy * (mj-1.5);
        k_dot_r_hz_yn(mii,mjj,mkk) = ...
            (x*k_vec_x + y*k_vec_y + z*k_vec_z)/c;

        mj = uj;
        y = min_y + dy * (mj-1);
        k_dot_r_ex_yp(mii,mjj,mkk) = ...
            (x*k_vec_x + y*k_vec_y + z*k_vec_z)/c;

        y = min_y + dy * (mj-0.5);
        k_dot_r_hz_yp(mii,mjj,mkk) = ...
            (x*k_vec_x + y*k_vec_y + z*k_vec_z)/c;

    end 
end 

for mj=lj:uj
    for mi=li:ui-1
        mii = mi - li + 1;
        mjj = mj - lj + 1;
        mkk = 1;

        mk = lk;                        
        x = min_x + dx * (mi-0.5);
        y = min_y + dy * (mj-1);
        z = min_z + dz * (mk-1);
        
        k_dot_r_ex_zn(mii,mjj,mkk) = ...
            (x*k_vec_x + y*k_vec_y + z*k_vec_z)/c;

        z = min_z + dz * (mk-1.5);
        k_dot_r_hy_zn(mii,mjj,mkk) = ...
            (x*k_vec_x + y*k_vec_y + z*k_vec_z)/c;

        mk = uk;                        
        z = min_z + dz * (mk-1);
        k_dot_r_ex_zp(mii,mjj,mkk) = ...
            (x*k_vec_x + y*k_vec_y + z*k_vec_z)/c;

        z = min_z + dz * (mk-0.5);
        k_dot_r_hy_zp(mii,mjj,mkk) = ...
            (x*k_vec_x + y*k_vec_y + z*k_vec_z)/c;
    end 
end 


for mj=lj:uj-1
    for mk=lk:uk
        mii = 1;
        mjj = mj - lj + 1;
        mkk = mk - lk + 1;

        mi = li;                        
        x = min_x + dx * (mi-1);
        y = min_y + dy * (mj-0.5);
        z = min_z + dz * (mk-1);
        
        k_dot_r_ey_xn(mii,mjj,mkk) = ...
            (x*k_vec_x + y*k_vec_y + z*k_vec_z)/c;

        x = min_x + dx * (mi-1.5);
        k_dot_r_hz_xn(mii,mjj,mkk) = ...
            (x*k_vec_x + y*k_vec_y + z*k_vec_z)/c;

        mi = ui;                        
        x = min_x + dx * (mi-1);
        k_dot_r_ey_xp(mii,mjj,mkk) = ...
            (x*k_vec_x + y*k_vec_y + z*k_vec_z)/c;

        x = min_x + dx * (mi-0.5);
        k_dot_r_hz_xp(mii,mjj,mkk) = ...
            (x*k_vec_x + y*k_vec_y + z*k_vec_z)/c;
    end 
end 

for mj=lj:uj-1
    for mi=li:ui
        mii = mi - li + 1;
        mjj = mj - lj + 1;
        mkk = 1;

        mk = lk;                        
        x = min_x + dx * (mi-1);
        y = min_y + dy * (mj-0.5);
        z = min_z + dz * (mk-1);
        k_dot_r_ey_zn(mii,mjj,mkk) = ...
            (x*k_vec_x + y*k_vec_y + z*k_vec_z)/c;

        z = min_z + dz * (mk-1.5);
        k_dot_r_hx_zn(mii,mjj,mkk) = ...
            (x*k_vec_x + y*k_vec_y + z*k_vec_z)/c;

        mk = uk;                        
        z = min_z + dz * (mk-1);
        k_dot_r_ey_zp(mii,mjj,mkk) = ...
            (x*k_vec_x + y*k_vec_y + z*k_vec_z)/c;

        z = min_z + dz * (mk-0.5);
        k_dot_r_hx_zp(mii,mjj,mkk) = ...
            (x*k_vec_x + y*k_vec_y + z*k_vec_z)/c;
    end 
end

for mj=lj:uj
    for mk=lk:uk-1
        mii = 1;
        mjj = mj - lj + 1;
        mkk = mk - lk + 1;

        mi = li;                        
        x = min_x + dx * (mi-1);
        y = min_y + dy * (mj-1);
        z = min_z + dz * (mk-0.5);
        k_dot_r_ez_xn(mii,mjj,mkk) = ...
            (x*k_vec_x + y*k_vec_y + z*k_vec_z)/c;

        x = min_x + dx * (mi-1.5);
        k_dot_r_hy_xn(mii,mjj,mkk) = ...
            (x*k_vec_x + y*k_vec_y + z*k_vec_z)/c;

        mi = ui;                        
        x = min_x + dx * (mi-1);
        k_dot_r_ez_xp(mii,mjj,mkk) = ...
            (x*k_vec_x + y*k_vec_y + z*k_vec_z)/c;

        x = min_x + dx * (mi-0.5);
        k_dot_r_hy_xp(mii,mjj,mkk) = ...
            (x*k_vec_x + y*k_vec_y + z*k_vec_z)/c;

    end 
end 


for mi=li:ui
    for mk=lk:uk-1
        mjj = 1;
        mii = mi - li + 1;
        mkk = mk - lk + 1;

        mj = lj;                        
        x = min_x + dx * (mi-1);
        y = min_y + dy * (mj-1);
        z = min_z + dz * (mk-0.5);
        k_dot_r_ez_yn(mii,mjj,mkk) = ...
            (x*k_vec_x + y*k_vec_y + z*k_vec_z)/c;

        y = min_y + dy * (mj-1.5);
        k_dot_r_hx_yn(mii,mjj,mkk) = ...
            (x*k_vec_x + y*k_vec_y + z*k_vec_z)/c;

        mj = uj;                        
        y = min_y + dy * (mj-1);
        k_dot_r_ez_yp(mii,mjj,mkk) = ...
            (x*k_vec_x + y*k_vec_y + z*k_vec_z)/c;

        y = min_y + dy * (mj-0.5);
        k_dot_r_hx_yp(mii,mjj,mkk) = ...
            (x*k_vec_x + y*k_vec_y + z*k_vec_z)/c;
    end
end

% embed spatial shift in k.r
k_dot_r_ex_yn = k_dot_r_ex_yn - l_0;   
k_dot_r_ex_yp = k_dot_r_ex_yp - l_0;   
k_dot_r_ex_zn = k_dot_r_ex_zn - l_0;   
k_dot_r_ex_zp = k_dot_r_ex_zp - l_0;   

k_dot_r_ey_zn = k_dot_r_ey_zn - l_0;   
k_dot_r_ey_zp = k_dot_r_ey_zp - l_0;   
k_dot_r_ey_xn = k_dot_r_ey_xn - l_0;   
k_dot_r_ey_xp = k_dot_r_ey_xp - l_0;   

k_dot_r_ez_xn = k_dot_r_ez_xn - l_0;   
k_dot_r_ez_xp = k_dot_r_ez_xp - l_0;   
k_dot_r_ez_yn = k_dot_r_ez_yn - l_0;   
k_dot_r_ez_yp = k_dot_r_ez_yp - l_0;   

k_dot_r_hx_yn = k_dot_r_hx_yn - l_0;   
k_dot_r_hx_yp = k_dot_r_hx_yp - l_0;   
k_dot_r_hx_zn = k_dot_r_hx_zn - l_0;   
k_dot_r_hx_zp = k_dot_r_hx_zp - l_0;   

k_dot_r_hy_zn = k_dot_r_hy_zn - l_0;   
k_dot_r_hy_zp = k_dot_r_hy_zp - l_0;   
k_dot_r_hy_xn = k_dot_r_hy_xn - l_0;   
k_dot_r_hy_xp = k_dot_r_hy_xp - l_0;   

k_dot_r_hz_xn = k_dot_r_hz_xn - l_0;   
k_dot_r_hz_xp = k_dot_r_hz_xp - l_0;   
k_dot_r_hz_yn = k_dot_r_hz_yn - l_0;   
k_dot_r_hz_yp = k_dot_r_hz_yp - l_0;

current_object = incident_plane_wave;
calculate_waveform;
incident_plane_wave.waveform = current_object.waveform;
incident_plane_wave.tau = current_object.tau;
incident_plane_wave.t_0 = current_object.t_0;
   
incident_plane_wave.ui = ui;
incident_plane_wave.uj = uj;
incident_plane_wave.uk = uk;
incident_plane_wave.li = li;
incident_plane_wave.lj = lj;
incident_plane_wave.lk = lk;
