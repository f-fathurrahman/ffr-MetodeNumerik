% calculate the waveform for a source 

% initialize sinusoidal waveform with a gaussian start
if strcmp(current_object.waveform_type,'sinusoidal')
    % gaussian transition
    nc = 20;
    maximum_frequency = c/(nc*max([dx,dy,dz]));
    tau = (nc*max([dx,dy,dz]))/(2*c);
    t_0 = 4.5 * tau;
    gaussian_waveform = exp(-((time - t_0)/tau).^2);
    I = find(time>t_0);
    gaussian_waveform(I) = 1;
    % sinusiodal waveform   
    if isfield(current_object, 'phase')
        w_phase = current_object.phase*pi/180;
    else
        w_phase = 0;
    end
    current_object.waveform = ...
        sin(2 * pi * current_object.frequency * time + w_phase) ...
        .* gaussian_waveform;        
end

% initialize unit step waveform with a gaussian start
if strcmp(current_object.waveform_type,'unit_step')
    % gaussian transition
    nc = 20;
    maximum_frequency = c/(nc*max([dx,dy,dz]));
    tau = (nc*max([dx,dy,dz]))/(2*c);
    t_0 = 4.5 * tau;
    gaussian_waveform = exp(-((time - t_0)/tau).^2);
    I = find(time>t_0);
    gaussian_waveform(I) = 1;
    current_object.waveform = gaussian_waveform;        
end

% initialize Gaussian waveform
if strcmp(current_object.waveform_type,'gaussian')
    if isfield(current_object, 'number_of_cells_per_wavelength')
        nc = current_object.number_of_cells_per_wavelength;
    else
        nc = 20;
    end
    current_object.maximum_frequency = c/(nc*max([dx,dy,dz]));
    tau = (nc*max([dx,dy,dz]))/(2*c);
    t_0 = 4.5 * tau;
    current_object.waveform = exp(-((time - t_0)/tau).^2);
    current_object.tau = tau;        
    current_object.t_0 = t_0;        
end

% initialize derivative of Gaussian waveform
if strcmp(current_object.waveform_type,'derivative_gaussian')
    if isfield(current_object, 'number_of_cells_per_wavelength')
        nc = current_object.number_of_cells_per_wavelength;
    else
        nc = 20;
    end
    current_object.maximum_frequency = c/(nc*max([dx,dy,dz]));
    tau = (nc*max([dx,dy,dz]))/(2*c);
    t_0 = 4.5 * tau;
    current_object.waveform = ...
        -(sqrt(2*exp(1))/tau)*(time - t_0).*exp(-((time - t_0)/tau).^2);
    current_object.tau = tau;        
    current_object.t_0 = t_0;        
end

% initialize cosine modulated Gaussian waveform
if strcmp(current_object.waveform_type,'cosine_modulated_gaussian')
    frequency = ...
        current_object.modulation_frequency;
    tau = 0.966/current_object.bandwidth;
    t_0 = 4.5 * tau;
    current_object.waveform = ...
        cos(2*pi*frequency*(time - t_0)).*exp(-((time - t_0)/tau).^2);
    current_object.tau = tau;        
    current_object.t_0 = t_0;        
end
