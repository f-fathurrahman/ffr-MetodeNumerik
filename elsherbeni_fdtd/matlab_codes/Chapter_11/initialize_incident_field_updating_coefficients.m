% initialize incident field updating coefficients

if incident_plane_wave.enabled == false
    return;
end

% Coeffiecients updating Ex
Cexeic= (2*(1-eps_r_x)*eps_0-dt*sigma_e_x) ...
    ./(2*eps_r_x*eps_0+dt*sigma_e_x);
Cexeip=-(2*(1-eps_r_x)*eps_0+dt*sigma_e_x) ...
    ./(2*eps_r_x*eps_0+dt*sigma_e_x);

% Coeffiecients updating Ey
Ceyeic= (2*(1-eps_r_y)*eps_0-dt*sigma_e_y) ...
    ./(2*eps_r_y*eps_0+dt*sigma_e_y);
Ceyeip=-(2*(1-eps_r_y)*eps_0+dt*sigma_e_y) ...
    ./(2*eps_r_y*eps_0+dt*sigma_e_y);

% Coeffiecients updating Ez
Cezeic= (2*(1-eps_r_z)*eps_0-dt*sigma_e_z) ...
    ./(2*eps_r_z*eps_0+dt*sigma_e_z);
Cezeip=-(2*(1-eps_r_z)*eps_0+dt*sigma_e_z) ...
    ./(2*eps_r_z*eps_0+dt*sigma_e_z);

% Coeffiecients updating Hx
Chxhic= (2*(1-mu_r_x)*mu_0-dt*sigma_m_x)./(2*mu_r_x*mu_0+dt*sigma_m_x);
Chxhip=-(2*(1-mu_r_x)*mu_0+dt*sigma_m_x)./(2*mu_r_x*mu_0+dt*sigma_m_x);

% Coeffiecients updating Hy
Chyhic= (2*(1-mu_r_y)*mu_0-dt*sigma_m_y)./(2*mu_r_y*mu_0+dt*sigma_m_y);
Chyhip=-(2*(1-mu_r_y)*mu_0+dt*sigma_m_y)./(2*mu_r_y*mu_0+dt*sigma_m_y);

% Coeffiecients updating Hz
Chzhic=(2*(1-mu_r_z)*mu_0-dt*sigma_m_z)./(2*mu_r_z*mu_0+dt*sigma_m_z);
Chzhip=-(2*(1-mu_r_z)*mu_0+dt*sigma_m_z)./(2*mu_r_z*mu_0+dt*sigma_m_z);
