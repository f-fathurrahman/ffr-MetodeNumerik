% matlab script waveeq2dloss.m
% finite difference scheme for the 2D wave equation with loss 
% fixed boundary conditions
% raised cosine initial conditions
% bilinear interpolation

%%%%%% begin global parameters

SR = 16000;                         % sample rate(Hz)
gamma = 200;                        % wave speed (1/s)
T60 = 8;                            % loss [freq.(Hz), T60(s), freq.(Hz), T60(s)]
epsilon = 1.3;                      % domain aspect ratio
TF = 2;                             % duration of simulation(s)
ctr = [0.3 0.5]; wid = 0.15;        % center location/width of excitation
u0 = 0; v0 = 1;                     % maximum initial displacement/velocity
rp = [0.5 0.6];                     % position of readout([0-1,0-1])
lambda = 1/sqrt(2);                 % Courant number

%%%%%% end global parameters

% begin derived parameters

k = 1/SR;                           % time step
NF = floor(SR*TF);                  % duration of simulation (samples)
sig0 = 6*log(10)/T60;               % loss parameter

% stability condition/scheme parameters

h = gamma*k/lambda;                 % find grid spacing
Nx = floor(sqrt(epsilon)/h);        % number of x-subdivisions of spatial domain
Ny = floor(1/(sqrt(epsilon)*h));    % number of y-subdivisions of spatial domain
h = sqrt(epsilon)/Nx; lambda = gamma*k/h;                        % reset Courant number
s0 = (2-4*lambda^2)/(1+sig0*k); s1 = lambda^2/(1+sig0*k); t0 = -(1-sig0*k)/(1+sig0*k);

% readout interpolation parameters

rp_int = 1+floor([Nx Ny].*rp); rp_frac = 1+rp/h-rp_int;                   

% create 2D raised cosine

[X, Y] = meshgrid([0:Nx]*h, [0:Ny]*h); dist = sqrt((X-ctr(1)).^2+(Y-ctr(2)).^2);
ind = sign(max(-dist+wid/2,0)); rc = 0.5*ind'.*(1+cos(2*pi*dist'/wid));

% set initial conditions

u2 = u0*rc; u1 = (u0+k*v0)*rc; u = zeros(Nx+1,Ny+1); out = zeros(NF,2); 

%%%%%% start main loop

for n=3:NF
    u(2:Nx,2:Ny) = s1*(u1(3:Nx+1,2:Ny)+u1(1:Nx-1,2:Ny)+u1(2:Nx,3:Ny+1)+u1(2:Nx,1:Ny-1))...
        +s0*u1(2:Nx,2:Ny)+t0*u2(2:Nx,2:Ny);
    out(n,:) = (1-rp_frac(1))*(1-rp_frac(2))*u(rp_int(1),rp_int(2))+...
        (1-rp_frac(1))*rp_frac(2)*u(rp_int(1),rp_int(2)+1)+...
        rp_frac(1)*(1-rp_frac(2))*u(rp_int(1)+1,rp_int(2))+...
        rp_frac(1)*rp_frac(2)*u(rp_int(1)+1,rp_int(2)+1);
    u2 = u1; u1 = u;                                         
end

%%%%% end main loop

% plot output waveform

plot([0:NF-1]*k, out, 'k'); xlabel('t'); ylabel('u'); 
title('2D Wave Equation with Loss: FD Output'); axis tight

% play sound

soundsc(out,SR);                    
    
