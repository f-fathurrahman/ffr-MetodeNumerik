% matlab script waveeq1dmod.m
% modal synthesis method for the 1D wave equation 
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

% begin derived parameters/temporary storage

k = 1/SR;                         % time step
NF = floor(SR*TF);                % duration of simulation (samples)
N = floor(0.5*SR/f0);             % number of modes
temp = 2*pi*[1:N]'*f0/SR; coeff = 2*cos(temp); outexp = sin([1:N]*pi*rp);

% initialize grid functions and output

U = zeros(N,1); U1 = zeros(N,1); U2 = zeros(N,1); 
out2 = zeros(NF,1); 

% create raised cosine and find Fourier coefficients

xax = [0:N-1]'/N;
ind = sign(max(-(xax-ctr-wid/2).*(xax-ctr+wid/2),0));
rc = 0.5*ind.*(1+cos(2*pi*(xax-ctr)/wid));
rcfs = -imag(fft([rc; zeros(N,1)])); rcfs = 2*rcfs(2:N+1)/N;

% set initial conditions

U2(1:N) = u0*rcfs; 
U1(1:N) = (u0*cos(temp)+v0*sin(temp)./(2*pi*[1:N]'*f0)).*rcfs; 

%%%%%% start main loop

for n=3:NF
    U = -U2+coeff.*U1;            % scheme calculation
    out(n) = outexp*U;            % readout
    U2 = U1; U1 = U;              % update of modal weights
end

%%%%%% end main loop

% plot output waveform

plot([0:NF-1]*k, out, 'k'); 
xlabel('t'); ylabel('u'); title('1D Wave Equation: Modal Synthesis Output');
axis tight

% play sound

soundsc(out,SR);                    
    
