disp (['Starting the time marching loop']);
disp(['Total number of time steps : ' ...
    num2str(number_of_time_steps)]); 

current_time = 0;
start_time = tic; 
if computation_platform ~=1
    launch_executable;
else
    for time_step = 1:number_of_time_steps  
        update_magnetic_fields_2d;
        update_impressed_M;
        update_magnetic_fields_for_CPML_2d;
        capture_sampled_magnetic_fields_2d;
        update_electric_fields_2d;
        update_impressed_J;
        update_electric_fields_for_CPML_2d;
        capture_sampled_electric_fields_2d;
        display_sampled_parameters_2d;
    end 
end
elapsed_time = toc(start_time);
total_time_in_minutes = elapsed_time/60;
disp(['Total simulation time is ' ...
    num2str(total_time_in_minutes) ' minutes.']);

