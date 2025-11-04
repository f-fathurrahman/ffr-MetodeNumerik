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

% wavenumber components
k = frequencies*2*pi/c;
kx = periodic_boundary.kx;
ky = periodic_boundary.ky;
kh = sqrt(kx^2+ky^2);
kz = sqrt(k.^2-kh^2);

phi_incident = atan2(ky,kx);

eta_0 = sqrt(mu_0/eps_0);

if strcmp(periodic_boundary.mode,'TE')
    Erp_co = (reflection_ey_dft*kx/kh-reflection_ex_dft*ky/kh);
    Hrp_co = (reflection_hx_dft*kx/kh+reflection_hy_dft*ky/kh);
    Eref_cr = (reflection_ex_dft*kx/kh+reflection_ey_dft*ky/kh);
    
    Einc = (Erp_co+eta_0*Hrp_co.*k./kz)/2;
    Eref_co = Erp_co-Einc;

    Gamma_co = Eref_co./Einc;
    aGamma_co = abs(Gamma_co);                    % Magnitude
    pGamma_co = angle(Gamma_co)/pi*180;           % Phase in degree
    
    Gamma_cr  = Eref_cr./Einc;
    aGamma_cr = abs(Gamma_cr);                    % Magnitude
    pGamma_cr = angle(Gamma_cr)/pi*180;           % Phase in degree
        
    if periodic_boundary.calculate_transmission
        Etra_co = (transmission_ey_dft*kx/kh-transmission_ex_dft*ky/kh);
        Etra_cr = (transmission_ex_dft*kx/kh+transmission_ey_dft*ky/kh);
        
        Einc = Einc .* ...
            exp(j*kz*periodic_boundary.reflection_transmission_distance);
        
        T_co = Etra_co./Einc;
        aT_co = abs(T_co);                    % Magnitude
        pT_co = angle(T_co)/pi*180;           % Phase in degree

        T_cr = Etra_cr./Einc;
        aT_cr = abs(T_cr);                    % Magnitude
        pT_cr = angle(T_cr)/pi*180;           % Phase in degree
    end
end

if strcmp(periodic_boundary.mode,'TM')
    Erp_co  = (-reflection_ex_dft*kx/kh-reflection_ey_dft*ky/kh);
    Hrp_co  = (reflection_hy_dft*kx/kh-reflection_hx_dft*ky/kh);
    Href_cr = (reflection_hx_dft*kx/kh+reflection_hy_dft*ky/kh);
        
    Hinc = (Hrp_co+Erp_co.*k./(eta_0*kz))/2;
    Href_co = Hrp_co-Hinc;

    Gamma_co = Href_co./Hinc;
    aGamma_co = abs(Gamma_co);                    % Magnitude
    pGamma_co = angle(Gamma_co)/pi*180;           % Phase in degree
  
    Gamma_cr = Href_cr./Hinc;
    aGamma_cr = abs(Gamma_cr);                    % Magnitude
    pGamma_cr = angle(Gamma_cr)/pi*180;           % Phase in degree
           
    if periodic_boundary.calculate_transmission
        Htra_co = (transmission_hy_dft*kx/kh-transmission_hx_dft*ky/kh);
        Htra_cr = (transmission_hx_dft*kx/kh+transmission_hy_dft*ky/kh);

        Hinc = Hinc .* ...
            exp(j*kz*periodic_boundary.reflection_transmission_distance);

        T_co = Htra_co./Hinc;
        aT_co = abs(T_co);                    % Magnitude
        pT_co = angle(T_co)/pi*180;           % Phase in degree

        T_cr = Htra_cr./Hinc;
        aT_cr = abs(T_cr);                    % Magnitude
        pT_cr = angle(T_cr)/pi*180;           % Phase in degree
    end
end

if strcmp(periodic_boundary.mode,'TEM')
    
    Exi0 = periodic_boundary.Exi0;
    Eyi0 = periodic_boundary.Eyi0;
    Ei0 = sqrt(Exi0^2+Eyi0^2);
    phi = atan2(Eyi0,Exi0);
    
    Erp_co  = reflection_ex_dft*cos(phi) + reflection_ey_dft*sin(phi);
    Hrp_co  = reflection_hx_dft*sin(phi) - reflection_hy_dft*cos(phi);
    Eref_cr = reflection_ex_dft*sin(phi) - reflection_ey_dft*cos(phi);
    
    Einc = (Erp_co+eta_0*Hrp_co)/2;
    Eref_co = Erp_co-Einc;

    Gamma_co = Eref_co./Einc;
    aGamma_co = abs(Gamma_co);                    % Magnitude
    pGamma_co = angle(Gamma_co)/pi*180;           % Phase in degree
    
    Gamma_cr  = Eref_cr./Einc;
    aGamma_cr = abs(Gamma_cr);                    % Magnitude
    pGamma_cr = angle(Gamma_cr)/pi*180;           % Phase in degree

    if periodic_boundary.calculate_transmission
        Etra_co = transmission_ex_dft*cos(phi) + transmission_ey_dft*sin(phi);
        Etra_cr = transmission_ex_dft*sin(phi) - transmission_ey_dft*cos(phi);

        Einc = Einc .* ...
            exp(j*kz*periodic_boundary.reflection_transmission_distance);
        
        T_co = Etra_co./Einc;
        aT_co = abs(T_co);                    % Magnitude
        pT_co = angle(T_co)/pi*180;           % Phase in degree

        T_cr = Etra_cr./Einc;
        aT_cr = abs(T_cr);                    % Magnitude
        pT_cr = angle(T_cr)/pi*180;           % Phase in degree
    end
end    

fGHz = frequencies*1e-9;

figure;
plot(fGHz,aGamma_co,'b-',fGHz,aGamma_cr,'r--','linewidth',1.5);
    legend('|\Gamma_{co}|','|\Gamma_{cr}|');
if periodic_boundary.calculate_transmission
    hold on;
    plot(fGHz,aT_co,'g-.',fGHz,aT_cr,'m:','linewidth',1.5);
    legend('|\Gamma_{co}|','|\Gamma_{cr}|','|T_{co}|', '|T_{cr}|');
end
grid on;
xlabel('frequency [GHz]','fontsize',12);
ylabel('magnitude','fontsize',12);

figure;
plot(fGHz,pGamma_co,'b-',fGHz,pGamma_cr,'r--','linewidth',1.5);
legend('|\Gamma_{co}|','|\Gamma_{cr}|');
if periodic_boundary.calculate_transmission
    hold on;
    plot(fGHz,pT_co,'g-.',fGHz,pT_cr,'m:','linewidth',1.5);
    legend('|\Gamma_{co}|','|\Gamma_{cr}|','|T_{co}|', '|T_{cr}|');
end
grid on;
xlabel('frequency [GHz]','fontsize',12);
ylabel('phase [degrees]','fontsize',12);
set(gca,'fontsize',12);

figure;
timens = time*1e9;
plot(timens, abs(periodic_boundary.reflection_ex),'b-', ...
     timens, abs(periodic_boundary.reflection_ey),'r--', ...
     timens, abs(periodic_boundary.reflection_ez),'k-.', ...
    'linewidth',1.5);
xlabel('time [ns]','fontsize',12);
ylabel('magnitude','fontsize',12);
grid on;
legend('reflection E_{x}','reflection E_{y}','reflection E_{z}');
set(gca,'fontsize',12);

if periodic_boundary.calculate_transmission
    figure;
    plot(timens, abs(periodic_boundary.transmission_ex),'b-', ...
         timens, abs(periodic_boundary.transmission_ey),'r--', ...
         timens, abs(periodic_boundary.transmission_ez),'k-.', ...
        'linewidth',1.5);
    xlabel('time [ns]','fontsize',12);
    ylabel('magnitude','fontsize',12);
    grid on;
    legend('transmission E_{x}','transmission E_{y}','transmission E_{z}');
    set(gca,'fontsize',12);
end
