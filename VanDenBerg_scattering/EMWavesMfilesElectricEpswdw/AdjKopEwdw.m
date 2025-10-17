 function [Kf,Kdf] = AdjKopEwdw(f,df,input)
global nDIM;

%correction factor for  interface integrals
 factor = 2 / (conj(input.gamma_0)^2 * input.dx); 

if nDIM == 2; 
    
   arg = f{1} + extrapolate(df{1},0,1,input);                 % x^1-to-x^0
 Kf{1} = AdjKop(arg,input.FFTG);    % ----------------- convolution in x^0
   arg = f{2} + extrapolate(df{2},0,2,input);                 % x^2-to-x^0
 Kf{2} = AdjKop(arg,input.FFTG);    % ----------------- convolution in x^0
 
   arg = extrapolate(f{1},1,0,input) + df{1};
 Kdf{1}=  - factor * gradj(AdjKop(arg,input.FFTG),1,input);         % convolution in x^1        
   arg = extrapolate(f{2},1,0,input) + interpolate(df{2},1,2,input);     
 Kdf{1}= Kdf{1} - factor * gradj(AdjKop(arg,input.FFTG),2,input);   % convolution in x^1

   arg = extrapolate( f{1},2,0,input) + interpolate(df{1},2,1,input);  
 Kdf{2}= - factor * gradj( AdjKop(arg,input.FFTG),1,input);         % convolution in x^2
   arg = extrapolate(f{2},2,0,input) + df{2};                     
 Kdf{2}= Kdf{2} - factor * gradj(AdjKop(arg,input.FFTG),2,input);   % convolution in x^2  
   

elseif nDIM == 3; 
    
   arg = f{1} + extrapolate(df{1},0,1,input);                 % x^1-to-x^0
 Kf{1} = AdjKop(arg,input.FFTG);    % ----------------- convolution in x^0
   arg = f{2} + extrapolate(df{2},0,2,input);                 % x^2-to-x^0
 Kf{2} = AdjKop(arg,input.FFTG);    % ----------------- convolution in x^0
   arg = f{3} + extrapolate(df{3},0,3,input);                 % x^3-to-x^0
 Kf{3} = AdjKop(arg,input.FFTG);    % ----------------- convolution in x^0
 
   arg = - extrapolate(gradj( f{1},1,input),1,0,input) ...    % x^0-to-x^1
         - extrapolate(gradj( f{2},2,input),1,0,input) ...    % x^0-to-x^1
         - gradj(df{1},1,input)                        ...    % x^1-to-x^1
         - interpolate(gradj(df{2},2,input),1,2,input);       % x^2-to-x^1
         - interpolate(gradj(df{3},3,input),1,3,input);       % x^3-to-x^1   
 Kdf{1}= factor * AdjKop(arg,input.FFTG);   % --------- convolution in x^1
   arg = - gradj(df{2},2,input)                        ...    % x^2-to-x^2
         - extrapolate(gradj( f{1},1,input),2,0,input) ...    % x^0-to-x^2
         - extrapolate(gradj( f{2},2,input),2,0,input) ...    % x^0-to-x^2
         - interpolate(gradj(df{1},1,input),2,1,input);       % x^1-to-x^2
         - interpolate(gradj(df{3},3,input),2,3,input);       % x^3-to-x^2   
 Kdf{2}= factor * AdjKop(arg,input.FFTG);   % --------- convolution in x^2
 
   arg = - gradj(df{3},3,input)                        ...    % x^3-to-x^3
         - extrapolate(gradj( f{1},1,input),3,0,input) ...    % x^0-to-x^3
         - extrapolate(gradj( f{3},3,input),3,0,input) ...    % x^0-to-x^3
         - interpolate(gradj(df{1},1,input),3,1,input);       % x^1-to-x^3
         - interpolate(gradj(df{2},2,input),3,2,input);       % x^2-to-x^3   
 Kdf{3}= factor * AdjKop(arg,input.FFTG);   % --------- convolution in x^3
end;  



     