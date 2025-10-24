% update incident fields for the current time step
if incident_plane_wave.enabled == false
    return;
end

tm = current_time + dt/2;
te = current_time + dt;

% update incident fields for previous time step
Hxip = Hxic; Hyip = Hyic; Hzip = Hzic;
Exip = Exic; Eyip = Eyic; Ezip = Ezic;

wt_str = incident_plane_wave.waveform_type;    
wi = incident_plane_wave.waveform_index; 
    
% if waveform is Gaussian waveforms
if strcmp(incident_plane_wave.waveform_type,'gaussian')
    tau = waveforms.gaussian(wi).tau;
    t_0 = waveforms.gaussian(wi).t_0;
    Exic = Exi0 * exp(-((te - t_0 - k_dot_r_ex )/tau).^2);
    Eyic = Eyi0 * exp(-((te - t_0 - k_dot_r_ey )/tau).^2);
    Ezic = Ezi0 * exp(-((te - t_0 - k_dot_r_ez )/tau).^2);
    Hxic = Hxi0 * exp(-((tm - t_0 - k_dot_r_hx )/tau).^2);
    Hyic = Hyi0 * exp(-((tm - t_0 - k_dot_r_hy )/tau).^2);
    Hzic = Hzi0 * exp(-((tm - t_0 - k_dot_r_hz )/tau).^2);
end

% if waveform is derivative of Gaussian 
if strcmp(incident_plane_wave.waveform_type,'derivative_gaussian')
    tau = waveforms.derivative_gaussian(wi).tau;
    t_0 = waveforms.derivative_gaussian(wi).t_0;
    Exic = Exi0 * (-sqrt(2*exp(1))/tau)*(te - t_0 - k_dot_r_ex) ...
           .*exp(-((te - t_0 - k_dot_r_ex)/tau).^2);
    Eyic = Eyi0 * (-sqrt(2*exp(1))/tau)*(te - t_0 - k_dot_r_ey) ...
           .*exp(-((te - t_0 - k_dot_r_ey)/tau).^2);
    Ezic = Ezi0 * (-sqrt(2*exp(1))/tau)*(te - t_0 - k_dot_r_ez) ...
           .*exp(-((te - t_0 - k_dot_r_ez)/tau).^2);
    Hxic = Hxi0 * (-sqrt(2*exp(1))/tau)*(tm - t_0 - k_dot_r_hx) ...
           .*exp(-((tm - t_0 - k_dot_r_hx)/tau).^2);
    Hyic = Hyi0 * (-sqrt(2*exp(1))/tau)*(tm - t_0 - k_dot_r_hy) ...
           .*exp(-((tm - t_0 - k_dot_r_hy)/tau).^2);
    Hzic = Hzi0 * (-sqrt(2*exp(1))/tau)*(tm - t_0 - k_dot_r_hz) ...
           .*exp(-((tm - t_0 - k_dot_r_hz)/tau).^2);
end

% if waveform is cosine modulated Gaussian 
if strcmp(incident_plane_wave.waveform_type, ...
    'cosine_modulated_gaussian')
    f = waveforms.cosine_modulated_gaussian(wi).modulation_frequency;
    tau = waveforms.cosine_modulated_gaussian(wi).tau;
    t_0 = waveforms.cosine_modulated_gaussian(wi).t_0;
    Exic = Exi0 * cos(2*pi*f*(te - t_0 - k_dot_r_ex)) ...
           .*exp(-((te - t_0 - k_dot_r_ex)/tau).^2);
    Eyic = Eyi0 * cos(2*pi*f*(te - t_0 - k_dot_r_ey)) ...
           .*exp(-((te - t_0 - k_dot_r_ey)/tau).^2);
    Ezic = Ezi0 * cos(2*pi*f*(te - t_0 - k_dot_r_ez)) ...
           .*exp(-((te - t_0 - k_dot_r_ez)/tau).^2);
    Hxic = Hxi0 * cos(2*pi*f*(tm - t_0 - k_dot_r_hx)) ...
           .*exp(-((tm - t_0 - k_dot_r_hx)/tau).^2);
    Hyic = Hyi0 * cos(2*pi*f*(tm - t_0 - k_dot_r_hy)) ...
           .*exp(-((tm - t_0 - k_dot_r_hy)/tau).^2);
    Hzic = Hyi0 * cos(2*pi*f*(tm - t_0 - k_dot_r_hz)) ...
           .*exp(-((tm - t_0 - k_dot_r_hz)/tau).^2);
end
