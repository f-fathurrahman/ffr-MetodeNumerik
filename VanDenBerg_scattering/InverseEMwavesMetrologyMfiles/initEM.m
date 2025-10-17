function input = initEM()
% Time factor = exp(-iwt) 
% Spatial units is in m
% Source wavelet M Z_0 / gamma_0  = 1   (Z_0 M = gamma_0)

global nDIM;  nDIM = 2;               % set dimension of space
if nDIM ==1;  disp('nDIM should be either 2 or 3');  return; end;
input.c_0     = 3e8;                  % wave speed in embedding
input.eps_sct = 1.55;                 % relative permittivity of scatterer
input.mu_sct  = 1;                    % relative permeability of scatterer

f             = 10e6;                 % temporal frequency
wavelength    = input.c_0 / f;        % wavelength
s             = 1e-16 - 1i*2*pi*f;    % LaPlace parameter
input.gamma_0 = s/input.c_0;          % propagation coefficient
disp(['wavelength = ' num2str(wavelength)]);
   
% add input data to structure array 'input'
  input = initSourceReceiver(input);  % add location of source/receiver 
  
  input = initEMgrid(input);          % add grid in either 2D or 3D
   
  input = initFFTGreen(input);        % compute FFT of Green function
  
  input = initEMcontrast(input);      % add contrast distribution
  
  input.Errcri = 1e-3;
  
  input.Noise = 0;
 
  if input.Noise == 1;  % generate random numbers that are repeatable
    sigma  = 0.2;       % sigma = standard deviation and mu =0 
    randn('state',1);   ReRand = randn(input.NR,input.NS) * sigma;
    randn('state',2);   ImRand = randn(input.NR,input.NS) * sigma;
    input.Rand = ReRand + 1i*ImRand;
  end;
  
  