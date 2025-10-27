% matlab script waveeq1dfd.m
% finite difference scheme for the 1D wave equation 
% fixed boundary conditions
% raised cosine initial conditions

%%%%%% begin global parameters

SR = 44100;                       % sample rate (Hz)
f0 = 441;                         % fundamental frequency (Hz)
TF = 1;                           % duration of simulation (s)
ctr = 0.7; wid = 0.1;             % center location/width of excitation
u0 = 1; v0 = 0;                   % maximum initial displacement/velocity
rp = 0.3;                         % position of readout (0-1)
lambda = 1;                       % Courant number

%%%%%% end global parameters

% begin derived parameters

gamma = 2*f0;                     % wave equation free parameter
k = 1/SR;                         % time step
NF = floor(SR*TF);                % duration of simulation (samples)

% stability condition/scheme parameters

h = gamma*k/lambda; N = floor(1/h); h = 1/N; lambda = gamma*k/h;         
s0 = 2*(1-lambda^2); s1 = lambda^2;       

% readout interpolation parameters

rp_int = 1+floor(N*rp);           % rounded grid index for readout
rp_frac = 1+rp/h-rp_int;          % fractional part of readout location

% create raised cosine

xax = [0:N]'*h;
ind = sign(max(-(xax-ctr-wid/2).*(xax-ctr+wid/2),0));
rc = 0.5*ind.*(1+cos(2*pi*(xax-ctr)/wid));

% initialize grid functions and output

u2 = u0*rc; u1 = (u0+k*v0)*rc; u = zeros(N+1,1); out = zeros(NF,1);

%%%%%% start main loop

for n=3:NF
    u(2:N) = -u2(2:N)+s0*u1(2:N)+s1*(u1(1:N-1)+u1(3:N+1));  % scheme calculation
    out(n) = (1-rp_frac)*u(rp_int)+rp_frac*u(rp_int+1);     % readout
    u2 = u1; u1 = u;                                        % update of grid variables
end

%%%%%% end main loop

% plot output waveform

plot([0:NF-1]*k, out, 'k'); 
xlabel('t'); ylabel('u'); title('1D Wave Equation: FD Output'); axis tight

% play sound

soundsc(out,SR);                    
    
