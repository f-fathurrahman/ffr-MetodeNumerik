function [Kf,Kdf] = AdjKOPwdw(f,df,input)
global nDIM;

 FFTG = input.FFTG;  FFTdG = cell(1,nDIM);
 for n=1: nDIM;   FFTdG{n} = input.FFTdG{n} * 2 / input.dx;  end

if nDIM == 1
    
    arg = f + extrapolate(df{1},0,1,input);     % from x^(1) back to x^(0)                         
     Kf = AdjKop(arg,FFTG);     % ------------- convolution in x^(0) space
      
    arg = df{1} + extrapolate(f,1,0,input);     % from x^(0) towards x^(1)             
 Kdf{1}= AdjKop(arg,FFTdG{1});  % ------------- convolution in x^(1) space  
  
elseif nDIM == 2 
    
   arg = f   + extrapolate(df{1},0,1,input);    % from x^(1) back to x^(0) 
   arg = arg + extrapolate(df{2},0,2,input);    % from x^(2) back to x^(0)                  
    Kf = AdjKop(arg,FFTG);       % ------------ convolution in x^(0) space
     
   arg = df{1} + extrapolate(f,1,0,input);     % from x^(0) towards x^(1) 
   arg = arg   + interpolate(df{2},1,2,input); % from x^(2) towards x^(1)          
 Kdf{1}= AdjKop(arg,FFTdG{1});   % ------------ convolution in x^(1) space
 
   arg = df{2} + extrapolate(f,2,0,input);     % from x^(0) towards x^(2)
   arg = arg   + interpolate(df{1},2,1,input); % from x^(1) towards x^(2)
 Kdf{2}= AdjKop(arg,FFTdG{2});   % ------------ convolution in x^(2) space 
           
elseif nDIM == 3
    
   arg = f   + extrapolate(df{1},0,1,input);    % from x^(1) back to x^(0)
   arg = arg + extrapolate(df{2},0,2,input);    % from x^(2) back to x^(0)
   arg = arg + extrapolate(df{3},0,3,input);    % from x^(3) back to x^(0)
    Kf = AdjKop(arg,FFTG);       % ------------ convolution in x^(0) space
    
   arg = df{1} + extrapolate(f,1,0,input);      % from x^(0) towards x^(1) 
   arg = arg   + interpolate(df{2},1,2,input);  % from x^(2) towards x^(1) 
   arg = arg   + interpolate(df{3},1,3,input);  % from x^(3) towards x^(1)
 Kdf{1}= AdjKop(arg,FFTdG{1});   % ------------ convolution in x^(1) space
 
   arg = df{2} + extrapolate(f,2,0,input);      % from x^(0) towards x^(2) 
   arg = arg   + interpolate(df{1},2,1,input);  % from x^(1) towards x^(2) 
   arg = arg   + interpolate(df{3},2,3,input);  % from x^(3) towards x^(2)
 Kdf{2}= AdjKop(arg,FFTdG{2});   % ------------ convolution in x^(2) space
 
   arg = df{3} + extrapolate(f,3,0,input);      % from x^(0) towards x^(3) 
   arg = arg   + interpolate(df{1},3,1,input);  % from x^(1) towards x^(3) 
   arg = arg   + interpolate(df{2},3,2,input);  % from x^(2) towards x^(3)
 Kdf{3}= AdjKop(arg,FFTdG{3});   % ------------ convolution in x^(3) space 

end      
     