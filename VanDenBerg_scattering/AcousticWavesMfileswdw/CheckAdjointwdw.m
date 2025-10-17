function CheckAdjointwdw(r,dr,input)
global nDIM;
  CHI = input.CHI;  Rfl = input.Rfl;
dummy = cell(1,nDIM);

% Adjoint operator [I - K*conj(CHI)] via inner product --------------------
for n = 1:nDIM  
    dummy{n} = conj(Rfl{n}).*dr{n};  
end
[Kf,dKf] = AdjKOPwdw(conj(CHI).*r,dummy,input); 
Result1  = sum( r(:) .* conj(r(:)-Kf(:)) );
for n = 1:nDIM  
   dummy{n} = dr{n} - dKf{n};  
   Result1  = Result1 + sum( dr{n}(:) .* conj(dummy{n}(:)) );  
end  
[Kf,Kdf] = KOPwdw(r,dr,input); 

% Operator via [I - CHI K] via inner product ------------------------------
Result2  = sum( (r(:)-CHI(:).*Kf(:)) .* conj(r(:)));
for n = 1:nDIM
   dummy{n} = dr{n} - Rfl{n} .* Kdf{n};                                
   Result2  = Result2 + sum( dummy{n}(:) .* conj(dr{n}(:)) ); 
end  

% Print the difference ----------------------------------------------------
fprintf('Check adjoint: %e\n',abs(Result1-Result2));