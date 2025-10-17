function input = initAC()
% Time factor = exp(-iwt);                  % Source wavelet s rho_0 Q = 1
% Spatial units is in m
global nDIM;  nDIM = 2;                     % set dimension of space

input.rho_0   = 2500;                       % mass density in embedding
input.rho_sct = 1500;                       % mass density in scatterer
input.c_0     = 1500;                       % wave speed in embedding             
% wave speed in scatterer 
input.c_sct   = input.c_0 * sqrt(input.rho_0/input.rho_sct); 
 
f             = 50;                         % temporal frequency
wavelength    = input.c_0 / f;              % wavelength
s             = 1e-16 - 1i*2*pi*f;          % LaPlace parameter
input.gamma_0 = s/input.c_0;                % propagation coefficient                                         
disp(['wavelength = ' num2str(wavelength)]);
   
% add input data to structure array 'input'
  input = initSourceReceiver(input);     % add location of source/receiver
  
  input = initGrid(input);               % add grid in either 1D, 2D or 3D
  
  input = initFFTGreenfun(input);        % compute FFT of Green function
   
  input = initAcousticContrast(input);   % add contrast distribution
  
  input.Errcri = 1e-5;
  
  input.Noise = 0;
 
  if input.Noise == 1  % generate random numbers that are repeatable
    sigma  = 0.2;       % sigma = standard deviation and mu =0 
    randn('state',1);   ReRand = randn(input.NR,input.NS) * sigma;
    randn('state',2);   ImRand = randn(input.NR,input.NS) * sigma;
    input.Rand = ReRand + 1i*ImRand;
  end