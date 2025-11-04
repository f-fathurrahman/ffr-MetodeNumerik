if ~periodic_boundary.calculate_reflection
    return;
end
frequencies = linspace(periodic_boundary.frequency_start, ...
    periodic_boundary.frequency_end, 200);

time_shift = 0;
reflection_ex_dft = time_to_frequency_domain( ...
    periodic_boundary.reflection_ex, dt, frequencies, time_shift);
reflection_ey_dft = time_to_frequency_domain( ...
    periodic_boundary.reflection_ey, dt, frequencies, time_shift);
reflection_ez_dft = time_to_frequency_domain( ...
    periodic_boundary.reflection_ez, dt, frequencies, time_shift);
time_shift = -dt/2;
reflection_hx_dft = time_to_frequency_domain( ...
    periodic_boundary.reflection_hx, dt, frequencies, time_shift);
reflection_hy_dft = time_to_frequency_domain( ...
    periodic_boundary.reflection_hy, dt, frequencies, time_shift);
reflection_hz_dft = time_to_frequency_domain( ...
    periodic_boundary.reflection_hz, dt, frequencies, time_shift);    
    
if periodic_boundary.calculate_transmission
    time_shift = 0;
    transmission_ex_dft = time_to_frequency_domain( ...
        periodic_boundary.transmission_ex, dt, frequencies, time_shift);
    transmission_ey_dft = time_to_frequency_domain( ...
        periodic_boundary.transmission_ey, dt, frequencies, time_shift);
    transmission_ez_dft = time_to_frequency_domain( ...
        periodic_boundary.transmission_ez, dt, frequencies, time_shift);
    time_shift = -dt/2;
    transmission_hx_dft = time_to_frequency_domain( ...
        periodic_boundary.transmission_hx, dt, frequencies, time_shift);
    transmission_hy_dft = time_to_frequency_domain( ...
        periodic_boundary.transmission_hy, dt, frequencies, time_shift);
    transmission_hz_dft = time_to_frequency_domain( ...
        periodic_boundary.transmission_hz, dt, frequencies, time_shift);    
end

% Free space wavenumber
k_0 = frequencies*2*pi/c;
% kz wavenumber
kx = periodic_boundary.kx;
ky = periodic_boundary.ky;
kxy = sqrt(kx^2+ky^2);
kz = conj(sqrt(k_0.^2-kxy^2));
phi_incident = atan2(ky,kx);

eta_0 = sqrt(mu_0/eps_0);

if strcmp(periodic_boundary.mode,'TE')
    Eco = (reflection_ey_dft*kx/kxy-reflection_ex_dft*ky/kxy);
    Hco = (reflection_hx_dft*kx/kxy+reflection_hy_dft*ky/kxy);
    Ecr = (reflection_ex_dft*kx/kxy+reflection_ey_dft*ky/kxy);
    
    eta = Eco./Hco; 
    Zi_TE = eta;
    Z0_TE = eta_0 * k_0./kz;

    Gamma_co  = (Zi_TE-Z0_TE)./(Zi_TE+Z0_TE);
    Eco_inc   = Eco./(Gamma_co+1);
    aGamma_co = abs(Gamma_co);                    % Magnitude
    pGamma_co = angle(Gamma_co)/pi*180;           % Phase in degree
    
    
    Gamma_cr  = Ecr./Eco_inc;
    aGamma_cr = abs(Gamma_cr);                    % Magnitude
    pGamma_cr = angle(Gamma_cr)/pi*180;           % Phase in degree
    
    if periodic_boundary.calculate_transmission
        Etco = (transmission_ey_dft*kx/kxy-transmission_ex_dft*ky/kxy);

        % Transmission coefficient
        T  = Etco./Eco_inc;
        aT = abs(T);                    % Magnitude
        pT = angle(T)/pi*180;           % Phase in degree
    end
end    
    
if strcmp(periodic_boundary.mode,'TM')
    Eco = (-reflection_ex_dft*kx/k_t-reflection_ey_dft*ky/k_t);
    Hco = (reflection_hy_dft*kx/k_t-reflection_hx_dft*ky/k_t);
    Ecr  = (reflection_ex_dft*ky/k_t-reflection_ey_dft*kx/k_t);

    eta = Eco./Hco;
    Zi_TM = eta;
    Z0_TM = eta_0 * kz./k0;
    
    Gamma_co  = -(Zi_TM-Z0_TM)./(Zi_TM+Z0_TM);
    aGamma_co = abs(Gamma_co);                    % Magnitude
    pGamma_co = angle(Gamma_co)/pi*180;           % Phase in degree
  
    Eco_inc = Eco./(Gamma_co+1);
    Hco_inc = Hco./(Gamma_co+1);
    Gamma_cr = Ecr./Eco_inc;
    aGamma_cr = abs(Rx);                    % Magnitude
    pGamma_cr = angle(Rx)/pi*180;           % Phase in degree
    
    if periodic_boundary.calculate_transmission
        Htco = (transmission_hy_dft*kx/k_t-transmission_hx_dft*ky/k_t);
        Etco = (-transmission_ex_dft*kx/k_t-transmission_ey_dft*ky/k_t);

        % Transmission coefficient
        T  = Htco./Hco_inc;
        aT = abs(T);                    % Magnitude
        pT = angle(T)/pi*180;           % Phase in degree
    end
    
end

if strcmp(periodic_boundary.mode,'TEM')
    
    Exi0 = periodic_boundary.Exi0;
    Eyi0 = periodic_boundary.Eyi0;
    Ei0 = sqrt(Exi0^2+Eyi0^2);
    phi = atan2(Eyi0,Exi0);
    
    Eco = reflection_ex_dft*cos(phi) + reflection_ey_dft*sin(phi);
    Hco = reflection_hx_dft*sin(phi) - reflection_hy_dft*cos(phi);
    Ecr = reflection_ex_dft*sin(phi) - reflection_ey_dft*cos(phi);
    
    eta = Eco./Hco; 
    Zi_TE = eta;
    Z0_TE = eta_0 * k_0./kz;

    Gamma_co  = (Zi_TE-Z0_TE)./(Zi_TE+Z0_TE);
    Eco_inc   = Eco./(Gamma_co+1);
    aGamma_co = abs(Gamma_co);                    % Magnitude
    pGamma_co = angle(Gamma_co)/pi*180;           % Phase in degree
    
    
    Gamma_cr  = Ecr./Eco_inc;
    aGamma_cr = abs(Gamma_cr);                    % Magnitude
    pGamma_cr = angle(Gamma_cr)/pi*180;           % Phase in degree
    
    if periodic_boundary.calculate_transmission
        Etco = transmission_ex_dft*cos(phi) + transmission_ey_dft*sin(phi);

        % Transmission coefficient
        T  = Etco./Eco_inc;
        aT = abs(T);                    % Magnitude
        pT = angle(T)/pi*180;           % Phase in degree
    end
end    

fGHz = frequencies*1e-9;

figure;
plot(fGHz,aGamma_co,'b-',fGHz,aGamma_cr,'r--','linewidth',2);
legend('\Gamma co-pol','\Gamma x-pol');
if periodic_boundary.calculate_transmission
    hold on;
    plot(fGHz,aT,'g-.','linewidth',2);
    legend('\Gamma co-pol','\Gamma x-pol','T');
end
grid on;
xlabel('Frequency [GHz]');
ylabel('Magnitude');

figure;
plot(fGHz,pGamma_co,'b-',fGHz,pGamma_cr,'r--','linewidth',2);
legend('\Gamma co-pol','\Gamma x-pol');
if periodic_boundary.calculate_transmission
    hold on;
    plot(fGHz,pT,'g-.','linewidth',2);
    legend('\Gamma co-pol','\Gamma x-pol','T');
end
grid on;
xlabel('Frequency [GHz]');
ylabel('Phase [degrees]');

