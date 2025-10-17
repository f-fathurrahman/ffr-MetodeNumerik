clear all; clc; close all; clear workspace;
input = initEM(); 

input = WaveletE(input); 
Wavelfrq = input.Wavelfrq;   
Wavelmax = max(abs(Wavelfrq(:)));

% Redefine frequency-independent parameters
  input.N1 = 600;    N1 = input.N1;             % number of samples in x_1             
  input.N2 = 600;    N2 = input.N2;             % number of samples in x_2 
  input.dx = 1;                              % grid size
    x1 = -(input.N1+1)*input.dx/2 + (1:input.N1)*input.dx;   
    x2 = -(input.N2+1)*input.dx/2 + (1:input.N2)*input.dx;
  [input.X1,input.X2] = ndgrid(x1,x2);
   input.xS  = [0 ,-170];                     % source position      
  [input] = initEMcontrast(input);            % Input contrast
  [input] = initEMcontrastIntf(input);        % Input interface contrast
  
Errcrr_0 = input.Errcri;         Ef = cell(2,2);
                                 Ef{2,1} = zeros(N1,N2,input.Nfft);
                                 Ef{1,2} = zeros(N1,N2,input.Nfft);
for f = 2 : input.fsamples  
%  (0) Make error criterion frequency dependent
       factor = abs(Wavelfrq(f))/Wavelmax;
       Errcri = min([Errcrr_0 / factor, 0.999]);     
       input.Errcri = Errcri;
%  (1) Redefine  frequency parameters
       freq = (f-1) * input.df;        disp(['freq sample: ',num2str(f)]);
       s = 1e-16 - 1i*2*pi*freq;           % LaPlace parameter
       input.gamma_0 = s/input.c_0;        % propagation coefficient
       input = initFFTGreen(input);        % compute FFT of Green function 
%  (2) Compute incident field ---------------------------------------------      
       [Einc_tau] = IncEtaufield(input);        
%  (3) Solve integral equation for contrast source with FFT ---------------
       [dw] = ITERBiCGSTABEdw(Einc_tau,input);
%  (4) Compute total  wave field on grid and add to frequency components
       E_sct = KopEdw(dw,input);    
       % take average value at midpoints between interfaces
       E_sct{2,1}(2:N1,:) = (E_sct{2,1}(1:N1-1,:) + E_sct{2,1}(2:N1,:))/2;
       E_sct{1,2}(:,2:N2) = (E_sct{1,2}(:,1:N2-1) + E_sct{1,2}(:,2:N2))/2;
       Ef{2,1}(:,:,f) = Wavelfrq(f).*(Einc_tau{2,1}(:,:)+E_sct{2,1}(:,:));
       Ef{1,2}(:,:,f) = Wavelfrq(f).*(Einc_tau{1,2}(:,:)+E_sct{1,2}(:,:));            
end; % frequency loop

Et{2,1} = 2 * input.df * real(fft(Ef{2,1},[],3)); 
Et{1,2} = 2 * input.df * real(fft(Ef{1,2},[],3));            clear Ef;
Etot = sqrt(abs(Et{2,1}).^2+ abs(Et{1,2}).^2 );              clear Et;
Etime(:,:,1:input.fsamples) = Etot(:,:,1:input.fsamples);    clear Etot;

SnapshotE;                        % Make snapshots for a few time instants