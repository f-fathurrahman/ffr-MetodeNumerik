% update incident fields for the current time step
if incident_plane_wave.enable == false
    return;
end

tm = current_time + dt/2;
te = current_time + dt;
 
% if waveform is Gaussian waveforms
tau = incident_plane_wave.tau;
t_0 = incident_plane_wave.t_0;

Exi_yn = Exi0 * exp(-((tm - t_0 - k_dot_r_ex_yn )/tau).^2);
Exi_yp = Exi0 * exp(-((tm - t_0 - k_dot_r_ex_yp )/tau).^2);
Exi_zn = Exi0 * exp(-((tm - t_0 - k_dot_r_ex_zn )/tau).^2);
Exi_zp = Exi0 * exp(-((tm - t_0 - k_dot_r_ex_zp )/tau).^2);

Eyi_zn = Eyi0 * exp(-((tm - t_0 - k_dot_r_ey_zn )/tau).^2);
Eyi_zp = Eyi0 * exp(-((tm - t_0 - k_dot_r_ey_zp )/tau).^2);
Eyi_xn = Eyi0 * exp(-((tm - t_0 - k_dot_r_ey_xn )/tau).^2);
Eyi_xp = Eyi0 * exp(-((tm - t_0 - k_dot_r_ey_xp )/tau).^2);

Ezi_xn = Ezi0 * exp(-((tm - t_0 - k_dot_r_ez_xn )/tau).^2);
Ezi_xp = Ezi0 * exp(-((tm - t_0 - k_dot_r_ez_xp )/tau).^2);
Ezi_yn = Ezi0 * exp(-((tm - t_0 - k_dot_r_ez_yn )/tau).^2);
Ezi_yp = Ezi0 * exp(-((tm - t_0 - k_dot_r_ez_yp )/tau).^2);

Hxi_yn = Hxi0 * exp(-((te - t_0 - k_dot_r_hx_yn )/tau).^2);
Hxi_yp = Hxi0 * exp(-((te - t_0 - k_dot_r_hx_yp )/tau).^2);
Hxi_zn = Hxi0 * exp(-((te - t_0 - k_dot_r_hx_zn )/tau).^2);
Hxi_zp = Hxi0 * exp(-((te - t_0 - k_dot_r_hx_zp )/tau).^2);

Hyi_zn = Hyi0 * exp(-((te - t_0 - k_dot_r_hy_zn )/tau).^2);
Hyi_zp = Hyi0 * exp(-((te - t_0 - k_dot_r_hy_zp )/tau).^2);
Hyi_xn = Hyi0 * exp(-((te - t_0 - k_dot_r_hy_xn )/tau).^2);
Hyi_xp = Hyi0 * exp(-((te - t_0 - k_dot_r_hy_xp )/tau).^2);

Hzi_xn = Hzi0 * exp(-((te - t_0 - k_dot_r_hz_xn )/tau).^2);
Hzi_xp = Hzi0 * exp(-((te - t_0 - k_dot_r_hz_xp )/tau).^2);
Hzi_yn = Hzi0 * exp(-((te - t_0 - k_dot_r_hz_yn )/tau).^2); 
Hzi_yp = Hzi0 * exp(-((te - t_0 - k_dot_r_hz_yp )/tau).^2);
