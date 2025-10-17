function input = init()
% Time factor = exp(-iwt); Spatial units is in m; Source wavelet  Q = 1
  global nDIM;  nDIM = 2;                     % dimension of space
         
  input.c_0     = 1500;         % wave speed in embedding
 
  input.c_sct   = 2000;         % wave speed in scatterer (CHI =   0.4375)
% input.c_sct   = 1200;         % wave speed in scatterer (CHI = - 0.5625)

  f             = 50;                         % temporal frequency
  wavelength    = input.c_0 / f;              % wavelength
  s             = 1e-16 - 1i*2*pi*f;          % LaPlace parameter
  input.gamma_0 = s/input.c_0;                % propagation coefficient
                                           
  disp(['wavelength = ' num2str(wavelength)]);
   
% add input data to structure array 'input'
  input = initSourceReceiver(input);   % add location of source/receiver
  
  input = initGrid(input);             % add grid in 2D
  
  input = initFFTGreen(input);         % compute FFT of Green function
  
  input = initContrast(input);         % add contrast distribution

  input.Errcri = 1e-5;
  
  input.Noise = 0;
  if input.Noise == 1  % generate random numbers that are repeatable
    sigma  = 0.2;       % sigma = standard deviation and mu =0 
    randn('state',1);   ReRand = randn(input.NR,input.NS) * sigma;
    randn('state',2);   ImRand = randn(input.NR,input.NS) * sigma;
    input.Rand = ReRand + 1i*ImRand;
  end