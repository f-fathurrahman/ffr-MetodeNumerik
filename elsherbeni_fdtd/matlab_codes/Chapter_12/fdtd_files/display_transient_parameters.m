disp('plotting the transient parameters');

% figures for sampled electric fields
for ind=1:number_of_sampled_electric_fields 
    if sampled_electric_fields(ind).display_plot == false
        sampled_electric_fields(ind).figure_number = figure;
        xlabel('time (ns)','fontsize',12);
        ylabel('(volt/meter)','fontsize',12);
        title(['sampled electric field [' num2str(ind) ']'],'fontsize',12);
        grid on; hold on;
    else
        figure(sampled_electric_fields(ind).figure_number);
        delete(sampled_electric_fields(ind).plot_handle);
    end
    sampled_time = sampled_electric_fields(ind).time(1:time_step)*1e9;
    sampled_value = sampled_electric_fields(ind).sampled_value(1:time_step);
    plot(sampled_time, sampled_value(1:time_step),'b-','linewidth',1.5);
    set(gca,'fontsize',12);
    drawnow;
end

% figures for sampled magnetic fields
for ind=1:number_of_sampled_magnetic_fields
    if sampled_magnetic_fields(ind).display_plot == false
        sampled_magnetic_fields(ind).figure_number = figure;
        xlabel('time (ns)','fontsize',12);
        ylabel('(ampere/meter)','fontsize',12);
        title(['sampled magnetic field [' num2str(ind) ']'],'fontsize',12);
        grid on; hold on;
    else
        figure(sampled_magnetic_fields(ind).figure_number);
        delete(sampled_magnetic_fields(ind).plot_handle);
    end
    sampled_time = sampled_magnetic_fields(ind).time(1:time_step)*1e9;
    sampled_value = sampled_magnetic_fields(ind).sampled_value(1:time_step);
    plot(sampled_time, sampled_value(1:time_step),'b-','linewidth',1.5);
    drawnow;
end

% figures for sampled voltages
for ind=1:number_of_sampled_voltages  
    if sampled_voltages(ind).display_plot == false
        sampled_voltages(ind).figure_number = figure;
        xlabel('time (ns)','fontsize',12);
        ylabel('(volt)','fontsize',12);
        title(['sampled voltage [' num2str(ind) ']'],'fontsize',12);
        grid on; hold on;
    else
        figure(sampled_voltages(ind).figure_number);
        delete(sampled_voltages(ind).plot_handle);
    end
    sampled_time = sampled_voltages(ind).time(1:time_step)*1e9;
    sampled_value = sampled_voltages(ind).sampled_value(1:time_step);
    plot(sampled_time, sampled_value(1:time_step),'b-','linewidth',1.5);
    drawnow;
end

% figures for sampled currents
for ind=1:number_of_sampled_currents  
    if sampled_currents(ind).display_plot == false
        sampled_currents(ind).figure_number = figure;
        xlabel('time (ns)','fontsize',12);
        ylabel('(ampere)','fontsize',12);
        title(['sampled current [' num2str(ind) ']'],'fontsize',12);
        grid on; hold on;
    else
        figure(sampled_currents(ind).figure_number);    
        delete(sampled_currents(ind).plot_handle);
    end
    sampled_time = sampled_currents(ind).time(1:time_step)*1e9;
    sampled_value = sampled_currents(ind).sampled_value(1:time_step);
    plot(sampled_time, sampled_value(1:time_step),'b-','linewidth',1.5);
    drawnow;
end

% figures for voltage sources
for ind=1:number_of_voltage_sources 
    voltage_sources(ind).figure_number = figure;
    sampled_time = time(1:time_step)*1e9;
    sampled_value = voltage_sources(ind).waveform(1:time_step);
    plot(sampled_time, sampled_value(1:time_step),'r-','linewidth',1.5);
    xlabel('time (ns)','fontsize',12);
    ylabel('(volt)','fontsize',12);
    title(['Voltage Source [' num2str(ind) ']'],'fontsize',12);
    grid on; 
    drawnow;
end

% figures for current sources
for ind=1:number_of_current_sources 
    current_sources(ind).figure_number = figure;
    sampled_time = time(1:time_step)*1e9;
    sampled_value = current_sources(ind).waveform(1:time_step);
    plot(sampled_time, sampled_value(1:time_step),'r-','linewidth',1.5);
    xlabel('time (ns)','fontsize',12);
    ylabel('(ampere)','fontsize',12);
    title(['Current Source [' num2str(ind) ']'],'fontsize',12);
    grid on; 
    drawnow;
end

% figure for incident plane wave 
if incident_plane_wave.enable == true
    figure;
    xlabel('time (ns)','fontsize',12);
    ylabel('(volt/meter)','fontsize',12);
    title('Incident electric field','fontsize',12);
    grid on; hold on;
    sampled_value_theta = incident_plane_wave.E_theta * ...
        incident_plane_wave.waveform;
    sampled_value_phi = incident_plane_wave.E_phi * ...
        incident_plane_wave.waveform;
    sampled_time = time(1:time_step)*1e9;
    plot(sampled_time, sampled_value_theta,'b-',...
        sampled_time, sampled_value_phi,'r:','linewidth',1.5);
    legend('E_{\theta,inc}','E_{\phi,inc}'); drawnow;
end
