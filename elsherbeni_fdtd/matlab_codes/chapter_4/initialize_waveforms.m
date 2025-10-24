disp('initializing source waveforms');

% initialize sinusoidal waveforms
for ind=1:size(waveforms.sinusoidal,2)
    waveforms.sinusoidal(ind).waveform = ...
        sin(2 * pi * waveforms.sinusoidal(ind).frequency * time);
end

% initialize unit step waveforms
for ind=1:size(waveforms.unit_step,2)
    start_index = waveforms.unit_step(ind).start_time_step;
    waveforms.unit_step(ind).waveform(1:number_of_time_steps) = 1;
    waveforms.unit_step(ind).waveform(1:start_index-1) = 0;
end
