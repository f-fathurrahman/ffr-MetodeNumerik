% initializing the incident plane wave

if incident_plane_wave.enable
    % create incident field arrays for current time step
    Hxic = zeros(nxp1,ny,nz);
    Hyic = zeros(nx,nyp1,nz);
    Hzic = zeros(nx,ny,nzp1);
    Exic = zeros(nx,nyp1,nzp1);
    Eyic = zeros(nxp1,ny,nzp1);
    Ezic = zeros(nxp1,nyp1,nz);
    % create incident field arrays for previous time step
    Hxip = zeros(nxp1,ny,nz);
    Hyip = zeros(nx,nyp1,nz);
    Hzip = zeros(nx,ny,nzp1);
    Exip = zeros(nx,nyp1,nzp1);
    Eyip = zeros(nxp1,ny,nzp1);
    Ezip = zeros(nxp1,nyp1,nz);
    
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
    
    % calculate k.r for every field component

    % Create position arrays indicating the coordinates of the nodes
    [x_pos, y_pos, z_pos] = ndgrid(node_coordinates_xh(2:nxp1), ...
        node_coordinates_ye, node_coordinates_ze); 

    k_dot_r_ex = (x_pos * k_vec_x ...
                + y_pos * k_vec_y ...
                + z_pos * k_vec_z)/c;
    
    [x_pos, y_pos, z_pos] = ndgrid(node_coordinates_xe, ...
        node_coordinates_yh(2:nyp1), node_coordinates_ze); 

    k_dot_r_ey = (x_pos * k_vec_x ...
                + y_pos * k_vec_y ...
                + z_pos * k_vec_z)/c;
    
    [x_pos, y_pos, z_pos] = ndgrid(node_coordinates_xe, ...
        node_coordinates_ye, node_coordinates_zh(2:nzp1)); 
    
    k_dot_r_ez = (x_pos * k_vec_x ...
                + y_pos * k_vec_y ...
                + z_pos * k_vec_z)/c;
    
    [x_pos, y_pos, z_pos] = ndgrid(node_coordinates_xe, ...
        node_coordinates_yh(2:nyp1), node_coordinates_zh(2:nzp1)); 

    k_dot_r_hx = (x_pos * k_vec_x ...
                + y_pos * k_vec_y ...
                + z_pos * k_vec_z)/c;

    [x_pos, y_pos, z_pos] = ndgrid(node_coordinates_xh(2:nxp1), ...
        node_coordinates_ye, node_coordinates_zh(2:nzp1)); 

    k_dot_r_hy = (x_pos * k_vec_x ...
                + y_pos * k_vec_y ...
                + z_pos * k_vec_z)/c;
    
    [x_pos, y_pos, z_pos] = ndgrid(node_coordinates_xh(2:nxp1), ...
        node_coordinates_yh(2:nyp1), node_coordinates_ze); 

    k_dot_r_hz = (x_pos * k_vec_x ...
                + y_pos * k_vec_y ...
                + z_pos * k_vec_z)/c;
   
    % embed spatial shift in k.r
    k_dot_r_ex = k_dot_r_ex - l_0;   
    k_dot_r_ey = k_dot_r_ey - l_0;   
    k_dot_r_ez = k_dot_r_ez - l_0;   
    k_dot_r_hx = k_dot_r_hx - l_0;   
    k_dot_r_hy = k_dot_r_hy - l_0;   
    k_dot_r_hz = k_dot_r_hz - l_0;   

    current_object = incident_plane_wave;
    calculate_waveform;
    incident_plane_wave.waveform = current_object.waveform;
    incident_plane_wave.tau = current_object.tau;
    incident_plane_wave.t_0 = current_object.t_0;
    
    clear x_pos y_pos z_pos; 
end
