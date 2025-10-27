% matlab script waveeq1ddw.m
% digital waveguide method for the 1D wave equation 
% fixed boundary conditions
% raised cosine initial conditions

%%%%%% begin global parameters

SR = 44100;                       % sample rate (Hz)
f0 = 441;                         % fundamental frequency (Hz)
TF = 1;                           % duration of simulation (s)
ctr = 0.7; wid = 0.1;             % center location/width of excitation
u0 = 1; v0 = 0;                   % maximum initial displacement/velocity
rp = 0.3;                         % position of readout (0-1)

%%%%%% end global parameters

% begin derived parameters

k = 1/SR;                         % time step
NF = floor(SR*TF);                % duration of simulation (samples)
N = floor(0.5*SR/f0);             % length of delay lines
rp_int = 1+floor(N*rp);           % rounded grid index for readout
rp_frac = 1+rp*N-rp_int;          % fractional part of readout location

% initialize delay lines and output

wleft = zeros(N,1); wright = zeros(N,1);  
out = zeros(NF,1); 

% create raised cosine and integral 

xax = ([1:N]'-1/2)/N;
ind = sign(max(-(xax-ctr-wid/2).*(xax-ctr+wid/2),0));
rc = 0.5*ind.*(1+cos(2*pi*(xax-ctr)/wid));
rcint = zeros(N,1);
for qq=2:N
    rcint(qq) = rcint(qq-1)+rc(qq)/N;
end

% set initial conditions

wleft = 0.5*(u0*rc+v0*rcint/(2*f0)); 
wright = 0.5*(u0*rc-v0*rcint/(2*f0));

%%%%%% start main loop

for n=3:NF
    temp1 = wright(N); temp2 = wleft(1);
    wright(2:N) = wright(1:N-1); wleft(1:N-1) = wleft(2:N);
    wright(1) = -temp2; wleft(N) = -temp1;
    % readout
    out(n) = (1-rp_frac)*(wleft(rp_int)+wright(rp_int))...
        +rp_frac*(wleft(rp_int+1)+wright(rp_int+1));     
end

%%%%%% end main loop

% plot output waveform

plot([0:NF-1]*k, out, 'k'); 
xlabel('t'); ylabel('u'); title('1D Wave Equation: Digital Waveguide Synthesis Output');
axis tight

% play sound

soundsc(out,SR);                    
    
