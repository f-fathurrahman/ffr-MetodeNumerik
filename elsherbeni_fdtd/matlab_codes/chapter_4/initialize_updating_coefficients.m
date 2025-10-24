disp('initializing general updating coefficients');

% General electric field updating coefficients
% Coeffiecients updating Ex
Cexe  =  (2*eps_r_x*eps_0 - dt*sigma_e_x) ...
    ./(2*eps_r_x*eps_0 + dt*sigma_e_x);
Cexhz =  (2*dt/dy)./(2*eps_r_x*eps_0 + dt*sigma_e_x);
Cexhy = -(2*dt/dz)./(2*eps_r_x*eps_0 + dt*sigma_e_x);

% Coeffiecients updating Ey
Ceye  =  (2*eps_r_y*eps_0 - dt*sigma_e_y) ...
    ./(2*eps_r_y*eps_0 + dt*sigma_e_y);
Ceyhx =  (2*dt/dz)./(2*eps_r_y*eps_0 + dt*sigma_e_y);
Ceyhz = -(2*dt/dx)./(2*eps_r_y*eps_0 + dt*sigma_e_y);

% Coeffiecients updating Ez
Ceze  =  (2*eps_r_z*eps_0 - dt*sigma_e_z) ...
    ./(2*eps_r_z*eps_0 + dt*sigma_e_z);
Cezhy =  (2*dt/dx)./(2*eps_r_z*eps_0 + dt*sigma_e_z);
Cezhx = -(2*dt/dy)./(2*eps_r_z*eps_0 + dt*sigma_e_z);

% General magnetic field updating coefficients
% Coeffiecients updating Hx
Chxh  =  (2*mu_r_x*mu_0 - dt*sigma_m_x) ...
    ./(2*mu_r_x*mu_0 + dt*sigma_m_x);
Chxez = -(2*dt/dy)./(2*mu_r_x*mu_0 + dt*sigma_m_x);
Chxey =  (2*dt/dz)./(2*mu_r_x*mu_0 + dt*sigma_m_x);

% Coeffiecients updating Hy
Chyh  =  (2*mu_r_y*mu_0 - dt*sigma_m_y) ...
    ./(2*mu_r_y*mu_0 + dt*sigma_m_y);
Chyex = -(2*dt/dz)./(2*mu_r_y*mu_0 + dt*sigma_m_y);
Chyez =  (2*dt/dx)./(2*mu_r_y*mu_0 + dt*sigma_m_y);

% Coeffiecients updating Hz
Chzh  =  (2*mu_r_z*mu_0 - dt*sigma_m_z) ...
    ./(2*mu_r_z*mu_0 + dt*sigma_m_z);
Chzey = -(2*dt/dx)./(2*mu_r_z*mu_0 + dt*sigma_m_z);
Chzex =  (2*dt/dy)./(2*mu_r_z*mu_0 + dt*sigma_m_z);

% Initialize coeffiecients for lumped element components
initialize_voltage_source_updating_coefficients;
initialize_current_source_updating_coefficients;
initialize_resistor_updating_coefficients;
initialize_capacitor_updating_coefficients;
initialize_inductor_updating_coefficients;
initialize_diode_updating_coefficients;

