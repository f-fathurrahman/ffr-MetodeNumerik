function input = init()
% Time factor = exp(-iwt) 
% Spatial units is in m
% Source wavelet  Q = 1
 
global nDIM;  nDIM = 2;                     % set dimension of space
         
input.c_0     = 1500;                       % wave speed in embedding
input.c_sct   = 3000;                       % wave speed in scatterer

f             = 50;                         % temporal frequency
wavelength    = input.c_0 / f;              % wavelength
s             = 1e-16 - 1i*2*pi*f;          % LaPlace parameter
input.gamma_0 = s/input.c_0;                % propagation coefficient
                                           
disp(['wavelength = ' num2str(wavelength)]);
   
% add input data to structure array 'input'
  input = initSourceReceiver(input);   % add location of source/receiver
  
  input = initGrid(input);             % add grid in either 1D, 2D or 3D
  
  input = initFFTGreen(input);         % compute FFT of Green function
  
  input = initContrast(input);         % add contrast distribution

  input.Errcri = 1e-3;