function CheckAdjointwcdrho(r_c,r_drho,input) 
global nDIM; 
CHI_c = input.CHI_c;  CHI_drho = input.CHI_drho;   dummy = cell(1,nDIM);

% Adjoint operator [I - K*conj(CHI)] via inner product --------------------
for n = 1:nDIM 
    dummy{n} = conj(CHI_drho{n}).*r_drho{n};  
end
[Kf,dKf] = AdjKOPwcdrho(conj(CHI_c).*r_c,dummy,input); 
Result1  = sum( r_c(:) .* conj(r_c(:)-Kf(:)) );
for n = 1:nDIM  
   dummy{n} = r_drho{n} - dKf{n};  
   Result1  = Result1 + sum( r_drho{n}(:) .* conj(dummy{n}(:)) );  
end

% Adjoint operator [I - CHI K] via inner product -------------------------- 
[Kwcdrho] = KOPwcdrho(r_c,r_drho,input); 
Result2  = sum( (r_c(:)-CHI_c(:).*Kwcdrho(:)) .* conj(r_c(:)));
for n = 1:nDIM
   dummy{n} = r_drho{n} - CHI_drho{n} .* Kwcdrho;                                
   Result2  = Result2 + sum( dummy{n}(:) .* conj(r_drho{n}(:)) ); 
end

% Print the differences
fprintf('Check adjoint: %e\n',abs(Result1-Result2));