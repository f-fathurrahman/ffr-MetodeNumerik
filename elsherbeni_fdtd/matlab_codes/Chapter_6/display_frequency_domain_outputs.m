disp('plotting the frequency domain parameters');

% figures for sampled electric fields
for ind=1:number_of_sampled_electric_fields 
    frequencies = sampled_electric_fields(ind).frequencies*1e-9;
    fd_value = sampled_electric_fields(ind).frequency_domain_value;
    figure;
    subplot(2,1,1);
    plot(frequencies, abs(fd_value),'b-','linewidth',1.5);
    title(['sampled electric field [' num2str(ind) ']'],'fontsize',12);
    xlabel('frequency (GHz)','fontsize',12);
    ylabel('magnitude','fontsize',12);
    grid on;
    subplot(2,1,2);
    plot(frequencies, angle(fd_value)*180/pi,'r-','linewidth',1.5);
    xlabel('frequency (GHz)','fontsize',12);
    ylabel('phase (degrees)','fontsize',12);
    grid on;
    drawnow;
end

% figures for sampled magnetic fields
for ind=1:number_of_sampled_magnetic_fields 
    frequencies = sampled_magnetic_fields(ind).frequencies*1e-9;
    fd_value = sampled_magnetic_fields(ind).frequency_domain_value;
    figure;
    subplot(2,1,1);
    plot(frequencies, abs(fd_value),'b-','linewidth',1.5);
    title(['sampled magnetic field [' num2str(ind) ']'],'fontsize',12);
    xlabel('frequency (GHz)','fontsize',12);
    ylabel('magnitude','fontsize',12);
    grid on;
    subplot(2,1,2);
    plot(frequencies, angle(fd_value)*180/pi,'r-','linewidth',1.5);
    xlabel('frequency (GHz)','fontsize',12);
    ylabel('phase (degrees)','fontsize',12);
    grid on;
    drawnow;
end

% figures for sampled voltages
for ind=1:number_of_sampled_voltages 
    frequencies = sampled_voltages(ind).frequencies*1e-9;
    fd_value = sampled_voltages(ind).frequency_domain_value;
    figure;
    subplot(2,1,1);
    plot(frequencies, abs(fd_value),'b-','linewidth',1.5);
    title(['sampled voltage [' num2str(ind) ']'],'fontsize',12);
    xlabel('frequency (GHz)','fontsize',12);
    ylabel('magnitude','fontsize',12);
    grid on;
    subplot(2,1,2);
    plot(frequencies, angle(fd_value)*180/pi,'r-','linewidth',1.5);
    xlabel('frequency (GHz)','fontsize',12);
    ylabel('phase (degrees)','fontsize',12);
    grid on;
    drawnow;
end

% figures for sampled currents
for ind=1:number_of_sampled_currents 
    frequencies = sampled_currents(ind).frequencies*1e-9;
    fd_value = sampled_currents(ind).frequency_domain_value;
    figure;
    subplot(2,1,1);
    plot(frequencies, abs(fd_value),'b-','linewidth',1.5);
    title(['sampled current [' num2str(ind) ']'],'fontsize',12);
    xlabel('frequency (GHz)','fontsize',12);
    ylabel('magnitude','fontsize',12);
    grid on;
    subplot(2,1,2);
    plot(frequencies, angle(fd_value)*180/pi,'r-','linewidth',1.5);
    xlabel('frequency (GHz)','fontsize',12);
    ylabel('phase (degrees)','fontsize',12);
    grid on;
    drawnow;
end

% figures for voltage sources
for ind=1:number_of_voltage_sources 
    frequencies = voltage_sources(ind).frequencies*1e-9;
    fd_value = voltage_sources(ind).frequency_domain_value;
    figure;
    subplot(2,1,1);
    plot(frequencies, abs(fd_value),'b-','linewidth',1.5);
    title(['voltage source [' num2str(ind) ']'],'fontsize',12);
    xlabel('frequency (GHz)','fontsize',12);
    ylabel('magnitude','fontsize',12);
    grid on;
    subplot(2,1,2);
    plot(frequencies, angle(fd_value)*180/pi,'r-','linewidth',1.5);
    xlabel('frequency (GHz)','fontsize',12);
    ylabel('phase (degrees)','fontsize',12);
    grid on;
    drawnow;
end

% figures for current sources
for ind=1:number_of_current_sources 
    frequencies = current_sources(ind).frequencies*1e-9;
    fd_value = current_sources(ind).frequency_domain_value;
    figure;
    subplot(2,1,1);
    plot(frequencies, abs(fd_value),'b-','linewidth',1.5);
    title(['current source [' num2str(ind) ']'],'fontsize',12);
    xlabel('frequency (GHz)','fontsize',12);
    ylabel('magnitude','fontsize',12);
    grid on;
    subplot(2,1,2);
    plot(frequencies, angle(fd_value)*180/pi,'r-','linewidth',1.5);
    xlabel('frequency (GHz)','fontsize',12);
    ylabel('phase (degrees)','fontsize',12);
    grid on;
    drawnow;
end

% figures for S-parameters
for ind=1:number_of_ports 
    if ports(ind).is_source_port == true
        frequencies = ports(ind).frequencies*1e-9;
        for oind=1:number_of_ports
            S = ports(ind).S(oind).values;
            Sdb = 20 * log10(abs(S));
            Sphase = angle(S)*180/pi;
            figure;
            subplot(2,1,1);
            plot(frequencies, Sdb,'b-','linewidth',1.5);
            title(['S' num2str(oind) num2str(ind)],'fontsize',12);
            xlabel('frequency (GHz)','fontsize',12);
            ylabel('magnitude (dB)','fontsize',12);
            grid on;
            subplot(2,1,2);
            plot(frequencies, Sphase,'r-','linewidth',1.5);
            xlabel('frequency (GHz)','fontsize',12);
            ylabel('phase (degrees)','fontsize',12);
            grid on;
            drawnow;
        end
    end
end
