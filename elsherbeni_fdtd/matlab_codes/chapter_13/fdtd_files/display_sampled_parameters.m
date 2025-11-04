% displaying sampled parameters

if mod(time_step,plotting_step) ~= 0   
    return;
end

remaining_time = (number_of_time_steps-time_step) ...
    *(cputime-start_time)/(60*time_step);
disp([num2str(time_step) ' of ' ...
    num2str(number_of_time_steps) ' is completed, ' ...
    num2str(remaining_time) ' minutes remaining']); 

% display sampled electric fields
for ind = 1:number_of_sampled_electric_fields
    if sampled_electric_fields(ind).display_plot == true
    sampled_time = sampled_electric_fields(ind).time(1:time_step)*1e9;
    sampled_value = sampled_electric_fields(ind).sampled_value(1:time_step);
    figure(sampled_electric_fields(ind).figure_number);
    delete(sampled_electric_fields(ind).plot_handle);
    sampled_electric_fields(ind).plot_handle = ...
    plot(sampled_time, sampled_value(1:time_step),'b-','linewidth',1.5);
    drawnow;
    end
end

% display sampled magnetic fields
for ind = 1:number_of_sampled_magnetic_fields
    if sampled_magnetic_fields(ind).display_plot == true
    sampled_time = sampled_magnetic_fields(ind).time(1:time_step)*1e9;
    sampled_value = sampled_magnetic_fields(ind).sampled_value(1:time_step);
    figure(sampled_magnetic_fields(ind).figure_number);
    delete(sampled_magnetic_fields(ind).plot_handle);
    sampled_magnetic_fields(ind).plot_handle = ...
    plot(sampled_time, sampled_value(1:time_step),'b-','linewidth',1.5);
    drawnow;
    end
end

% display sampled voltages
for ind = 1:number_of_sampled_voltages
    if sampled_voltages(ind).display_plot == true
    sampled_time = sampled_voltages(ind).time(1:time_step)*1e9;
    sampled_value = sampled_voltages(ind).sampled_value(1:time_step);
    figure(sampled_voltages(ind).figure_number);
    delete(sampled_voltages(ind).plot_handle);
    sampled_voltages(ind).plot_handle = ...
    plot(sampled_time, sampled_value(1:time_step),'b-','linewidth',1.5);
    drawnow;
    end
end

% display sampled currents
for ind = 1:number_of_sampled_currents
    if sampled_currents(ind).display_plot == true
    sampled_time = sampled_currents(ind).time(1:time_step)*1e9;
    sampled_value = sampled_currents(ind).sampled_value(1:time_step);
    figure(sampled_currents(ind).figure_number);
    delete(sampled_currents(ind).plot_handle);
    sampled_currents(ind).plot_handle = ...
    plot(sampled_time, sampled_value(1:time_step),'b-','linewidth',1.5);
    drawnow;
    end
end

% display animated fields
display_animation;
