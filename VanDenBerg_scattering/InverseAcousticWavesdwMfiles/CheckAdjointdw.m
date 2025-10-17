function CheckAdjointdw(dr,input)
global nDIM;
   Rfl = input.Rfl;
dummy = cell(1,nDIM);

% Adjoint operator [I - K*conj(CHI)] via inner product --------------------
for n = 1:nDIM  
    dummy{n} = conj(Rfl{n}).*dr{n};  
end
[dKf] = AdjKOPdw(dummy,input); 
Result1  = 0;
for n = 1:nDIM  
   dummy{n} = dr{n} - dKf{n};  
   Result1  = Result1 + sum( dr{n}(:) .* conj(dummy{n}(:)) );  
end  
[Kdf] = KOPdw(dr,input); 

% Operator via [I - CHI K] via inner product ------------------------------
Result2  = 0;
for n = 1:nDIM
   dummy{n} = dr{n} - Rfl{n} .* Kdf{n};                                
   Result2  = Result2 + sum( dummy{n}(:) .* conj(dr{n}(:)) ); 
end  

% Print the difference ----------------------------------------------------
fprintf('Check adjoint: %e\n',abs(Result1-Result2));