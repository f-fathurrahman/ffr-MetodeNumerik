% matlab script ssfd.m
% finite difference scheme for the stiff string 
% clamped boundary conditions
% raised cosine initial conditions
% stereo output
% implicit scheme: matrix update form
% two-parameter frequency-dependent loss

%%%%%% begin global parameters

SR = 44100;                                % sample rate(Hz)
B = 0.001;                                 % inharmonicity parameter (>0)
f0 = 100;                                  % fundamental(Hz)
TF = 2;                                    % duration of simulation(s)
ctr = 0.1; wid = 0.05;                     % center location/width of excitation
u0 = 1; v0 = 0;                            % maximum initial displacement/velocity
rp = [0.3 0.7];                            % positions of readout(0-1)
loss = [100, 10; 1000, 8];                 % loss [freq.(Hz), T60(s), freq.(Hz), T60(s)]
theta = 1.0;                               % implicit scheme free parameter (>0.5)

%%%%%% end global parameters

% begin derived parameters

k = 1/SR;                                  % time step
NF = floor(SR*TF);                         % duration of simulation (samples)
gamma = 2*f0; K = sqrt(B)*(gamma/pi);      % set parameters

% stability conditions

h = sqrt((gamma^2*k^2+sqrt(gamma^4*k^4+16*K^2*k^2*(2*theta-1)))/(2*(2*theta-1)));
N = floor(1/h); h = 1/N; mu = K*k/h^2; lambda = gamma*k/h;          

% readout interpolation parameters

rp_int = 1+floor(N*rp);           % rounded grid index for readout
rp_frac = 1+rp/h-rp_int;          % fractional part of readout location

% set scheme loss parameters

zeta1 = (-gamma^2+sqrt(gamma^4+4*K^2*(2*pi*loss(1,1))^2))/(2*K^2);
zeta2 = (-gamma^2+sqrt(gamma^4+4*K^2*(2*pi*loss(2,1))^2))/(2*K^2);
sig0 = 6*log(10)*(-zeta2/loss(1,2)+zeta1/loss(2,2))/(zeta1-zeta2);
sig1 = 6*log(10)*(1/loss(1,2)-1/loss(2,2))/(zeta1-zeta2);

% create update matrices

M = sparse(toeplitz([theta (1-theta)/2 zeros(1,N-3)]));
A = M+sparse(toeplitz([sig1*k/(h^2)+sig0*k/2 -sig1*k/(2*h^2) zeros(1,N-3)]));
C = M+sparse(toeplitz([-sig1*k/(h^2)-sig0*k/2 sig1*k/(2*h^2) zeros(1,N-3)]));
B = 2*M+sparse(toeplitz([-2*lambda^2-6*mu^2 lambda^2+4*mu^2 -mu^2 zeros(1,N-4)]));

% create raised cosine

xax = [1:N-1]'*h;
ind = sign(max(-(xax-ctr-wid/2).*(xax-ctr+wid/2),0));
rc = 0.5*ind.*(1+cos(2*pi*(xax-ctr)/wid));

% set initial conditions

u2 = u0*rc; u1 = (u0+k*v0)*rc; u = zeros(N+1,1); out = zeros(NF,2);

%%%%%% start main loop

for n=3:NF
    u = A\(B*u1-C*u2);                   
    out(n,:) = (1-rp_frac).*u(rp_int)'+rp_frac.*u(rp_int+1)'; % readout
    u2 = u1; u1 = u;                                          % update 
end

%%%%%% end main loop

% plot output waveform

subplot(2,1,1); plot([0:NF-1]*k, out(:,1), 'k');
xlabel('t'); ylabel('u'); title('Stiff String Equation: FD Output (left)');
subplot(2,1,2); plot([0:NF-1]*k, out(:,2), 'k');
xlabel('t'); ylabel('u'); title('Stiff String Equation: FD Output (right)');
axis tight

% play sound

soundsc(out,SR);                    
    
