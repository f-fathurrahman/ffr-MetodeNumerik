% matlab script sho.m
% finite difference scheme for simple harmonic oscillator

%%%%%% begin global parameters

SR = 44100;                  % sample rate (Hz)
f0 = 1000;                   % fundamental frequency (Hz)
TF = 1.0;                    % duration of simulation (s)
u0 = 0.3;                    % initial displacement
v0 = 0.0;                    % initial velocity

%%%%%% end global parameters

% check that stability condition is satisfied

if(SR<=pi*f0)
    error('Stability condition violated');
end

% derived parameters

k = 1/SR;                    % time step
coef = 2-k^2*(2*pi*f0)^2;    % scheme update coefficient
NF = floor(TF*SR);           % duration of simulation (samples)

% initialize state of scheme

u1 = u0+k*v0;                % last value of time series
u2 = u0;                     % one before last value of time series

% initialize readout

out = zeros(NF,1); out(1) = u2; out(2) = u1;

%%%%%% start main loop

for n=3:NF
    u=coef*u1-u2;            % difference scheme calculation
    out(n) = u;              % read value to output vector
    u2 = u1; u1 = u;         % update state of difference scheme
end

%%%%%% end main loop

% play sound

soundsc(out,SR);  

% plot output waveform

plot([0:NF-1]*k, out, 'k'); 
xlabel('t'); ylabel('u'); title('SHO: Scheme Output'); axis tight


                  
    
