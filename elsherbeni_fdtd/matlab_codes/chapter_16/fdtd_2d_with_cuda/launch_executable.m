disp('Launching executable to run time marching on exe!');

% save the data to a file to process on gpu
exe_data_file_name = 'exe_data_file.dat';
save_project_data_to_file;

cmd = ['fdtd2d.exe ' exe_data_file_name];
status = system(cmd);
if status
    disp('Calculations on with executable are not successful!');
    return;
end

time_step = number_of_time_steps;

result_data_file_name = ['result_' exe_data_file_name];
fid = fopen(result_data_file_name, 'r');

for ind = 1:number_of_sampled_electric_fields
    sampled_electric_fields(ind).sampled_value = ...
        fread(fid, number_of_time_steps, 'float32').';
end

for ind = 1:number_of_sampled_magnetic_fields
    sampled_magnetic_fields(ind).sampled_value = ...
        fread(fid, number_of_time_steps, 'float32').';
end

fclose(fid);