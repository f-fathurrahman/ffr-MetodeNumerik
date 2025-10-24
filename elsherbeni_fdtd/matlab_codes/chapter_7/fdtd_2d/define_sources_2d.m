disp('defining sources');

impressed_J = [];
impressed_M = [];

% define source waveform types and parameters
waveforms.gaussian(1).number_of_cells_per_wavelength = 0; 
waveforms.gaussian(2).number_of_cells_per_wavelength = 15; 

% electric current sources
% direction: 'xp', 'xn', 'yp', 'yn', 'zp', or 'zn'
impressed_J(1).min_x = -0.1e-3;
impressed_J(1).min_y = -0.1e-3;
impressed_J(1).max_x = 0.1e-3;
impressed_J(1).max_y = 0.1e-3;
impressed_J(1).direction = 'zp';
impressed_J(1).magnitude = 1;
impressed_J(1).waveform_type = 'gaussian';
impressed_J(1).waveform_index = 1;

% % magnetic current sources
% % direction: 'xp', 'xn', 'yp', 'yn', 'zp', or 'zn'
% impressed_M(1).min_x = -0.1e-3;
% impressed_M(1).min_y = -0.1e-3;
% impressed_M(1).max_x = 0.1e-3;
% impressed_M(1).max_y = 0.1e-3;
% impressed_M(1).direction = 'zp';
% impressed_M(1).magnitude = 1;
% impressed_M(1).waveform_type = 'gaussian';
% impressed_M(1).waveform_index = 1;
