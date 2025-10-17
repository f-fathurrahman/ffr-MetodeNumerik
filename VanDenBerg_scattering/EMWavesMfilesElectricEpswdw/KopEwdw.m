function [Kw,Kdw] = KopEwdw(w,dw,input)   
global nDIM;                

% Correction factor for interface integrals
factor = 2 / (input.gamma_0^2 * input.dx); 

if nDIM == 2;      
   Kw{1} = Kop(w{1},input.FFTG);           %---------------------%   K00_1
  Kdw{1} = extrapolate(Kw{1},1,0,input);                         %   K10
   Kw{2} = Kop(w{2},input.FFTG);           %---------------------%   K00_2
  Kdw{2} = extrapolate(Kw{2},2,0,input);                         %   K20   
           KOPdw  = factor * Kop(dw{1},input.FFTG); %----------------------
   Kw{1} = Kw{1}  + gradj(extrapolate(KOPdw,0,1,input),1,input); % + d1K01
  Kdw{1} = Kdw{1} + gradj(            KOPdw,           1,input); % + d1K11
   Kw{2} = Kw{2}  + gradj(extrapolate(KOPdw,0,1,input),2,input); % + d2K01
  Kdw{2} = Kdw{2} + gradj(interpolate(KOPdw,2,1,input),2,input); % + d2K21 
           KOPdw  = factor * Kop(dw{2},input.FFTG); %----------------------
   Kw{1} = Kw{1}  + gradj(extrapolate(KOPdw,0,2,input),1,input); % + d1K02
  Kdw{1} = Kdw{1} + gradj(interpolate(KOPdw,1,2,input),1,input); % + d1K12
   Kw{2} = Kw{2}  + gradj(extrapolate(KOPdw,0,2,input),2,input); % + d2K02
  Kdw{2} = Kdw{2} + gradj(            KOPdw,           2,input); % + d2K22
                   
elseif nDIM == 3;   
   Kw{1} = Kop(w{1},input.FFTG);           %---------------------%   K00_1
  Kdw{1} = extrapolate(Kw{1},1,0,input);                         %   K10
   Kw{2} = Kop(w{2},input.FFTG);           %---------------------%   K00_2
  Kdw{2} = extrapolate(Kw{2},2,0,input);                         %   K20  
   Kw{3} = Kop(w{3},input.FFTG);           %---------------------%   K00_3
  Kdw{3} = extrapolate(Kw{3},3,0,input);                         %   K30  
           KOPdw  = factor * Kop(dw{1},input.FFTG); %----------------------
   Kw{1} = Kw{1}  + gradj(extrapolate(KOPdw,0,1,input),1,input); % + d1K01
  Kdw{1} = Kdw{1} + gradj(            KOPdw,           1,input); % + d1K11
   Kw{2} = Kw{2}  + gradj(extrapolate(KOPdw,0,1,input),2,input); % + d2K01
  Kdw{2} = Kdw{2} + gradj(interpolate(KOPdw,2,1,input),2,input); % + d2K21
   Kw{3} = Kw{3}  + gradj(extrapolate(KOPdw,0,1,input),3,input); % + d3K01
  Kdw{3} = Kdw{3} + gradj(interpolate(KOPdw,3,1,input),3,input); % + d3K31
           KOPdw  = factor * Kop(dw{2},input.FFTG); %----------------------
   Kw{1} = Kw{1}  + gradj(extrapolate(KOPdw,0,2,input),1,input); % + d1K02
  Kdw{1} = Kdw{1} + gradj(interpolate(KOPdw,1,2,input),1,input); % + d1K12
   Kw{2} = Kw{2}  + gradj(extrapolate(KOPdw,0,2,input),2,input); % + d2K02
  Kdw{2} = Kdw{2} + gradj(            KOPdw,           2,input); % + d2K22
   Kw{3} = Kw{3}  + gradj(extrapolate(KOPdw,0,2,input),3,input); % + d3K02
  Kdw{3} = Kdw{3} + gradj(extrapolate(KOPdw,3,2,input),3,input); % + d3K32
           KOPdw  = factor * Kop(dw{3},input.FFTG); %----------------------
   Kw{1} = Kw{1}  + gradj(extrapolate(KOPdw,0,3,input),1,input); % + d1K03
  Kdw{1} = Kdw{1} + gradj(interpolate(KOPdw,1,3,input),1,input); % + d1K13
   Kw{2} = Kw{2}  + gradj(extrapolate(KOPdw,0,3,input),2,input); % + d2K03 
  Kdw{2} = Kdw{2} + gradj(extrapolate(KOPdw,2,3,input),2,input); % + d2K23 
   Kw{3} = Kw{3}  + gradj(extrapolate(KOPdw,0,3,input),3,input); % + d3K03 
  Kdw{3} = Kdw{3} + gradj(            KOPdw,           3,input); % + d3K33   
end;   