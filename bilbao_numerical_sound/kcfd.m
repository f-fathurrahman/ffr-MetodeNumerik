% matlab script kcfd.m
% finite difference scheme for the Kirchhoff-Carrier equation 
% fixed boundary conditions
% triangular initial conditions

%%%%%% begin global parameters

SR = 44100;                              % sample rate (Hz)
f0 = 200;                                % fundamental frequency (Hz)
alpha = 10;                              % nonlinear string parameter
TF = 0.03;                               % duration of simulation (s)
ctr = 0.5;                               % center location of excitation (0-1)
u0 = 0.05;                               % maximum initial displacement
rp = 0.5;                                % position of readout (0-1)
lambda = 0.7;                            % Courant number

%%%%%% end global parameters

% begin derived parameters

gamma = 2*f0;                            % wave equation free parameter
k = 1/SR;                                % time step
NF = floor(SR*TF);                       % duration of simulation (samples)

% stability condition

h = gamma*k/lambda; N = floor(1/h); h = 1/N; lambda = gamma*k/h;     

% readout interpolation parameters

rp_int = 1+floor(N*rp);                  % rounded grid index for readout
rp_frac = 1+rp/h-rp_int;                 % fractional part of readout location

% create triangular function

xax = [0:N]'*h;
tri = min(xax/ctr-1,0)+1+min((1-xax)/(1-ctr)-1,0);

% initialize grid functions and output

u2 = u0*tri; u1 = u2; u = zeros(N+1,1); out = zeros(NF,1); 

%%%%%% start main loop

for n=1:NF
    % calculate nonlinearity g
    u1x = (u1(2:N+1)-u1(1:N))/h; u1xx = (u1x(2:N)-u1x(1:N-1))/h;
    g = (1+0.5*alpha^2*h*sum(u1x.*u1x))/...
        (1+0.25*alpha^2*k^2*gamma^2*h*sum(u1xx.*u1xx)); 
    % scheme update
    u(2:N) = 2*u1(2:N)-u2(2:N)+g*gamma^2*k^2*u1xx(1:N-1);       % calculation
    out(n) = (1-rp_frac)*u(rp_int)+rp_frac*u(rp_int+1);         % readout
    u2 = u1; u1 = u;                                            % update
end

%%%%% end main loop

% plot output waveform

plot([0:NF-1]*k, out, 'k'); 
xlabel('t'); ylabel('u'); title('Kirchhoff-Carrier Equation: FD Output');
axis tight
                
    
