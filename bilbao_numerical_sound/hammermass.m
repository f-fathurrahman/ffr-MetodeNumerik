% matlab script hammermass.m
% hammer collision with a mass-spring system

%%%%%% begin global parameters

SR = 44100;                     % sample rate (Hz)
xH0 = -0.001; vH0 = 2;          % initial conditions of hammer
TF = 0.05;                      % duration of simulation (s)
w0 = 2000;                      % angular frequency of mass-spring system
MR = 10;                        % hammer/target mass ratio
wH = 1000;                      % stiffness parameter for hammer
alpha = 2;                      % hammer stiffness nonlinearity exponent

%%%%%% end global parameters

% derived parameters

k = 1/SR;
NF = floor(TF*SR); 

% initialization

uH2 = xH0; uH1 = xH0+k*vH0;         % hammer
u2 = 0; u1 = 0;                     % mass-spring system
out = zeros(NF,1); f = zeros(NF,1);
out(1) = u2; out(2) = u1;

%%%%%% start main loop

for n=3:NF
    if(uH1>u1)
        f(n-1) = wH^(1+alpha)*(uH1-u1)^alpha;
    else f(n-1) = 0;
    end
    uH = 2*uH1-uH2-k^2*f(n-1);
    u = 2*u1-u2-w0^2*k^2*u1+MR*k^2*f(n-1);
    out(n) = u;
    u2 = u1; u1 = u;
    uH2 = uH1; uH1 = uH;
end

%%%%%% end main loop

% plots of displacement of target mass and force

subplot(2,1,1)
plot([0:NF-1]*k, out, 'k'); title('Position of Target Mass'); xlabel('t');
axis tight
subplot(2,1,2)
plot([0:NF-1]*k, f, 'k'); title('Hammer Force/Mass'); xlabel('t');
axis tight    

