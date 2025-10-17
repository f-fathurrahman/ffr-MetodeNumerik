function CheckAdjointE(r_E,FFTG,input) 
global nDIM;  dummy = cell(1,nDIM);
CHI_eps = input.CHI_eps;

% check adjoint operator [I-K*CHI] via inner product ----------------------
for n = 1:nDIM;  
    dummy{n} = conj(CHI_eps).*r_E{n};  
end;
[KE] = AdjKopE(dummy,FFTG,input); 
Result1  = 0;
for n = 1:nDIM;  
   dummy{n} = r_E{n} - KE{n};  
   Result1  = Result1 + sum( r_E{n}(:) .* conj(dummy{n}(:)) );  
end;  
[KE] = KopE(r_E,FFTG,input); 
Result2  = 0;
for n = 1:nDIM;
   dummy{n} = r_E{n} - CHI_eps .* KE{n};                                
   Result2  = Result2 + sum( dummy{n}(:) .* conj(r_E{n}(:)) ); 
end;  

fprintf('Check adjoint: %e\n',abs(Result1-Result2));