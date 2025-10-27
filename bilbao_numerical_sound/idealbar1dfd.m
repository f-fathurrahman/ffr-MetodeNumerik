% matlab script idealbarfd.m
% finite difference scheme for the ideal bar equation 
% *fixed boundary conditions
% *raised cosine initial conditions

%%%%%% begin global parameters

SR = 44100;                                     % sample rate (Hz)
K = 1;                                          % stiffness parameter
TF = 1;                                         % duration of simulation (s)
ctr = 0.5;                                      % center location of excitation (0-1)
wid = 0.1;                                      % width of excitation 
u0 = 1;                                         % maximum initial displacement
v0 = 0;                                         % maximum initial velocity
mu = 0.5;                                       % scheme free parameter
rp = 0.5;                                      % position of readout (0-1)

%%%%%% end global parameters

%%%%%% begin derived parameters

k = 1/SR;                                       % time step
NF = floor(SR*TF);                              % duration of simulation (samples)
h = sqrt(K*k/mu);                               % grid spacing
N = floor(1/h);                                 % number of subdivisions of spatial domain
h = 1/N;                                        % reset h 
mu = K*k/h^2;                                   % reset Courant number
s0 = 2*(1-3*mu^2); s1 = 4*mu^2; s2 = -mu^2;     % scheme parameters  
rp_int = 1+floor(N*rp);                         % rounded grid index for readout
rp_frac = 1+rp/h-rp_int;                        % fractional part of readout location

%%%%%% initialize grid functions and output

u = zeros(N+1,1); u1 = zeros(N+1,1); u2 = zeros(N+1,1); 
out = zeros(NF,1); 

%%%%%% create raised cosine

rc = zeros(N+1,1);
for qq=1:N+1
    pos = (qq-1)*h; dist = ctr-pos; 
    if(abs(dist)<=wid/2)
        rc(qq) = 0.5*(1+cos(2*pi*dist/wid));
    end
end

%%%%%% set initial conditions

u2 = u0*rc; u1 = (u0+k*v0)*rc; 

%%%%%% start main loop

for n=3:NF
    u(3:N-1) = -u2(3:N-1)+s0*u1(3:N-1)+s1*(u1(2:N-2)+u1(4:N))...
        +s2*(u1(1:N-3)+u1(5:N+1));                          % scheme calculation
    out(n) = (1-rp_frac)*u(rp_int)+rp_frac*u(rp_int+1);     % readout
    u2 = u1; u1 = u;                                        % update of grid variables
end

%%%%% end main loop

% plot output waveform

plot([0:NF-1]*k, out, 'k'); 
xlabel('t'); ylabel('u'); title('Ideal Bar Equation: Difference Scheme Output');
axis tight

% play sound

soundsc(out,SR);                    
    
