function CheckAdjointpZv(r_p,r_Zv,input) 
global nDIM;  
CHI_kap = input.CHI_kap;   CHI_rho = input.CHI_rho; 
  dummy = cell(1,nDIM);

% Adjoint operator [I - K*conj(CHI)] via inner product --------------------
for n = 1:nDIM
    dummy{n} = conj(CHI_rho).*r_Zv{n};  
end
[Kp,KZv] = AdjKOPpZv(conj(CHI_kap).*r_p,dummy,input); 
Result1  = sum( r_p(:) .* conj(r_p(:)-Kp(:)) );
for n = 1:nDIM 
   dummy{n} = r_Zv{n} - KZv{n};  
   Result1  = Result1 + sum( r_Zv{n}(:) .* conj(dummy{n}(:)) );  
end

% Operator [I - CHI K] via inner product ----------------------------------
[Kp,KZv] = KOPpZv(r_p,r_Zv,input);
Result2  = sum( (r_p(:)-CHI_kap(:).*Kp(:)) .* conj(r_p(:)));
for n = 1:nDIM
   dummy{n} = r_Zv{n} - CHI_rho .* KZv{n};                                
   Result2  = Result2 + sum( dummy{n}(:) .* conj(r_Zv{n}(:)) ); 
end

% Print the difference ----------------------------------------------------
fprintf('Check adjoint: %e\n',abs(Result1-Result2));