% matlab script plateloss.m
% finite difference scheme for the thin plate equation with loss 
% simply-supported boundary conditions
% raised cosine initial conditions
% vector/matrix update form
% zeroth order interpolation

%%%%%% begin global parameters

SR = 44100;                         % sample rate(Hz)
K = 20;                             % plate stiffness parameter (1/s)
T60 = 8;                            % loss [freq.(Hz), T60(s), freq.(Hz), T60(s)]
epsilon = 1.2;                      % domain aspect ratio
TF = 2;                             % duration of simulation(s)
ctr = [0.8 0.9]; wid = 0.3;         % center location/width of excitation
u0 = 0; v0 = 1;                     % maximum initial displacement/velocity
rp = [0.05 0.7];                    % position of readout([0-1,0-1])
mu = 0.25;                          % scheme free parameter

%%%%%% end global parameters

% begin derived parameters

k = 1/SR;                           % time step
NF = floor(SR*TF);                  % duration of simulation (samples)
sig0 = 6*log(10)/T60;               % loss parameter

% stability condition/scheme parameters

h = sqrt(K*k/mu);                   % find grid spacing
Nx = floor(sqrt(epsilon)/h);        % number of x-subdivisions of spatial domain
Ny = floor(1/(sqrt(epsilon)*h));    % number of y-subdivisions of spatial domain
h = sqrt(epsilon)/Nx;                         
ss = (Nx-1)*(Ny-1);                 % total grid size

% generate difference matrix/scheme matrices

Dxx = sparse(toeplitz([-2/h^2;1/h^2;zeros(Nx-3,1)]));
Dyy = sparse(toeplitz([-2/h^2;1/h^2;zeros(Ny-3,1)]));
D = kron(eye(Nx-1), Dyy)+kron(Dxx, eye(Ny-1)); DD = D*D; 
B = sparse((2*eye(ss)-K^2*k^2*DD)/(1+sig0*k));
C = ((1-sig0*k)/(1+sig0*k))*sparse(eye(ss));  

% readout interpolation parameters

rp_index = (Ny-1)*floor(rp(1)*Nx)+floor(rp(2)*Ny);

% create 2D raised cosine

[X, Y] = meshgrid([1:Nx-1]*h, [1:Ny-1]*h); dist = sqrt((X-ctr(1)*sqrt(epsilon)).^2+(Y-ctr(2)/sqrt(epsilon)).^2);
ind = sign(max(-dist+wid/2,0)); rc = 0.5*ind.*(1+cos(2*pi*dist/wid));
rc = reshape(rc, ss,1);

% set initial conditions/initialize output

u2 = u0*rc; u1 = (u0+k*v0)*rc; u = zeros(ss,1);out = zeros(NF,1); 

%%%%%% start main loop

for n=3:NF
    u = B*u1-C*u2;
    u2 = u1; u1 = u;
    out(n) = u(rp_index);   
end

%%%%%% end main loop

% plot output waveform

plot([0:NF-1]*k, out, 'k'); xlabel('t'); ylabel('u'); 
title('Thin Plate Equation with Loss: FD Output'); axis tight

% play sound

soundsc(out,SR);                    
    
