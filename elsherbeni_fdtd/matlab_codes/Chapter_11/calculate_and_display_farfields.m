% This file calls the routines necessary for calculating 
% farfield patterns in xy, xz, and yz plane cuts, and displays them. 
% The display can be modified as desired.  
% You will find the instructions the formats for
% radiation pattern plots can be set by the user.

if number_of_farfield_frequencies == 0
    return;
end

calculate_radiated_power;
calculate_incident_plane_wave_power;

j = sqrt(-1);
number_of_angles = 360;

% parameters used by polar plotting functions
step_size = 10;       % increment between the rings in the polar grid
Nrings = 4;           % number of rings in the polar grid
line_style1 = 'b-';   % line style for theta component
line_style2 = 'r--';  % line style for phi component
scale_type = 'dB';    % linear or dB

if incident_plane_wave.enabled == false
    plot_type = 'D';
else
    plot_type = 'RCS';
end
% xy plane
% ===============================================
farfield_theta = zeros(number_of_angles, 1);
farfield_phi   = zeros(number_of_angles, 1);
farfield_theta = farfield_theta + pi/2;
farfield_phi = (pi/180)*[-180:1:179].';
const_theta = 90; % used for plot

% calculate farfields
calculate_farfields_per_plane;

% plotting the farfield data
for mi=1:number_of_farfield_frequencies 
 f = figure;
 pat1 = farfield_dataTheta(mi,:).';
 pat2 = farfield_dataPhi(mi,:).';

 % if scale_type is db use these, otherwise comment these two lines
 pat1 = 10*log10(pat1); 
 pat2 = 10*log10(pat2);

 max_val = max(max([pat1 pat2]));
 max_val = step_size * ceil(max_val/step_size);

 legend_str1 = ...
 [plot_type '_{\theta}, f=' num2str(farfield.frequencies(mi)*1e-9) ' GHz'];
 legend_str2 = ...
 [plot_type '_{\phi}, f=' num2str(farfield.frequencies(mi)*1e-9) ' GHz'];

 polar_plot_constant_theta(farfield_phi,pat1,pat2,max_val, ...
        step_size, Nrings,line_style1,line_style2,const_theta, ...
        legend_str1,legend_str2,scale_type);
end

% xz plane
% ===============================================
farfield_theta = zeros(number_of_angles, 1);
farfield_phi   = zeros(number_of_angles, 1);
farfield_theta = (pi/180)*[-180:1:179].';
const_phi = 0; % used for plot

% calculate farfields
calculate_farfields_per_plane;

% plotting the farfield data
for mi=1:number_of_farfield_frequencies 
 f = figure;
 pat1 = farfield_dataTheta(mi,:).';
 pat2 = farfield_dataPhi(mi,:).';

% if scale_type is db use these, otherwise comment these two lines
 pat1 = 10*log10(pat1); 
 pat2 = 10*log10(pat2);

 max_val = max(max([pat1 pat2]));
 max_val = step_size * ceil(max_val/step_size);
 
 legend_str1 = ...
 [plot_type '_{\theta}, f=' num2str(farfield.frequencies(mi)*1e-9) ' GHz'];
 legend_str2 = ...
 [plot_type '_{\phi}, f=' num2str(farfield.frequencies(mi)*1e-9) ' GHz'];

 polar_plot_constant_phi(farfield_theta,pat1,pat2,max_val, ...
        step_size, Nrings,line_style1,line_style2,const_phi, ...
        legend_str1,legend_str2,scale_type);
end

% yz plane
% ===============================================
farfield_theta = zeros(number_of_angles, 1);
farfield_phi   = zeros(number_of_angles, 1);
farfield_phi = farfield_phi + pi/2;
farfield_theta = (pi/180)*[-180:1:179].';
const_phi = 90; % used for plot

% calculate farfields
calculate_farfields_per_plane;

% plotting the farfield data
for mi=1:number_of_farfield_frequencies 
 f = figure;
 pat1 = farfield_dataTheta(mi,:).';
 pat2 = farfield_dataPhi(mi,:).';

% if scale_type is db use these, otherwise comment these two lines
 pat1 = 10*log10(pat1); 
 pat2 = 10*log10(pat2);

 max_val = max(max([pat1 pat2]));
 max_val = step_size * ceil(max_val/step_size);
 
 legend_str1 = ...
 [plot_type '_{\theta}, f=' num2str(farfield.frequencies(mi)*1e-9) ' GHz'];
 legend_str2 = ...
 [plot_type '_{\phi}, f=' num2str(farfield.frequencies(mi)*1e-9) ' GHz'];

 polar_plot_constant_phi(farfield_theta,pat1,pat2,max_val, ...
        step_size, Nrings,line_style1,line_style2,const_phi, ...
        legend_str1,legend_str2,scale_type);
end
