function [Kw,Kdw] = KOPwdw(w,dw,input)   
global nDIM;               

FFTG = input.FFTG;  FFTdG = cell(1,nDIM);
for n=1: nDIM;   FFTdG{n} = input.FFTdG{n} * 2 / input.dx;  end

if nDIM == 1   
    
   Kdiag = Kop(w,FFTG);  % ----------------------------------------  % K00 
      Kw = Kdiag;  Kdw{1} = extrapolate(Kdiag,1,0,input);            % K10                                                                                       
   Kdiag = Kop(dw{1},FFTdG{1}); % ----------------------------------------                                         
                   Kdw{1} = Kdw{1} + Kdiag;                          % K11
      Kw = Kw + extrapolate(Kdiag,0,1,input);                        % K01

elseif nDIM == 2  
    
   Kdiag = Kop(w,FFTG);  % ----------------------------------------  % K00
      Kw = Kdiag;  Kdw{1} = extrapolate(Kdiag,1,0,input);            % K10
                   Kdw{2} = extrapolate(Kdiag,2,0,input);            % K20                                     
   Kdiag = Kop(dw{1},FFTdG{1}); % ----------------------------------------                                      
                   Kdw{1} = Kdw{1} + Kdiag;                          % K11
      Kw = Kw + extrapolate(Kdiag,0,1,input);                        % K01
                   Kdw{2} = Kdw{2} + interpolate(Kdiag,2,1,input);   % K21
   Kdiag = Kop(dw{2},FFTdG{2}); % ----------------------------------------                                     
                   Kdw{2} = Kdw{2} + Kdiag;                          % K22
      Kw = Kw + extrapolate(Kdiag,0,2,input);                        % K02
                   Kdw{1} = Kdw{1} + interpolate(Kdiag,1,2,input);   % K12
           
elseif nDIM == 3  
    
   Kdiag = Kop(w,FFTG); % -------------------------------------------- K00
      Kw = Kdiag;  Kdw{1} = extrapolate(Kdiag,1,0,input);            % K10
                   Kdw{2} = extrapolate(Kdiag,2,0,input);            % K20
                   Kdw{3} = extrapolate(Kdiag,3,0,input);            % K30                                                    
   Kdiag = Kop(dw{1},FFTdG{1}); % ---------------------------------------- 
                   Kdw{1} = Kdw{1} + Kdiag;                          % K11
      Kw = Kw + extrapolate(Kdiag,0,1,input);                        % K01
                   Kdw{2} = Kdw{2} + interpolate(Kdiag,2,1,input);   % K21
                   Kdw{3} = Kdw{3} + interpolate(Kdiag,3,1,input);   % K31  
   Kdiag = Kop(dw{2},FFTdG{2}); % ---------------------------------------- 
                   Kdw{2} = Kdw{2} + Kdiag;                          % K22          
      Kw = Kw + extrapolate(Kdiag,0,2,input);                        % K02
                   Kdw{1} = Kdw{1} + interpolate(Kdiag,1,2,input);   % K12
                   Kdw{3} = Kdw{3} + interpolate(Kdiag,3,2,input);   % K32
   Kdiag = Kop(dw{3},FFTdG{3}); % ---------------------------------------- 
                   Kdw{3} = Kdw{3} + Kdiag;                          % K33
      Kw = Kw + extrapolate(Kdiag,0,3,input);                        % K03
                   Kdw{1} = Kdw{1} + interpolate(Kdiag,1,3,input);   % K13
                   Kdw{2} = Kdw{2} + interpolate(Kdiag,2,3,input);   % K23

end   