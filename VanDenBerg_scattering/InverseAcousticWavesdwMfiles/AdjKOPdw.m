function [Kdf] = AdjKOPdw(df,input)
global nDIM;
FFTdG = cell(1,nDIM);
for n = 1:nDIM;   FFTdG{n} = input.FFTdG{n} * 2 / input.dx;  end

  arg  = df{1} + interpolate(df{2},1,2,input);          
Kdf{1} = AdjKop(arg,FFTdG{1});  
 
  arg  = df{2} + interpolate(df{1},2,1,input);
Kdf{2} = AdjKop(arg,FFTdG{2});                