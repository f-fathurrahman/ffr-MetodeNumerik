function CheckAdjointEwdw(r,dr,input) 
global nDIM;  
 CHI  = input.CHI_eps;   Rfl = input.Rfl;  
dummy = cell(1,nDIM);  ddummy = cell(1,nDIM);


% Adjoint operator  via inner product -------------------------------------
for n = 1:nDIM;  
    dummy{n} = conj(CHI) .*  r{n}; 
   ddummy{n} = conj(Rfl{n}) .* dr{n};  
end;

[Kf,Kdf] = AdjKopEwdw(dummy,ddummy,input); 
Result1  = 0; 
for n = 1:nDIM; 
    dummy{n} =  r{n} -  Kf{n};
    ddummy{n} = dr{n} - Kdf{n};
   Result1   = Result1 + sum( r{n}(:).*conj( dummy{n}(:))) ...
                       + sum(dr{n}(:).*conj(ddummy{n}(:)));
end;

% Operator via inner product ------------------------------------------


[Kf,Kdf] =  KopEwdw(r,dr,input); 
Result2  = 0; 
for n = 1:nDIM; 
    dummy{n} =  r{n} -  CHI  .*  Kf{n};
   ddummy{n} = dr{n} -  Rfl{n} .* Kdf{n};
   Result2   = Result2 + sum( dummy{n}(:).*conj(r{n}(:) )) ...
                       + sum(ddummy{n}(:).*conj(dr{n}(:)));
end;
 
% Print the difference
fprintf('Check adjoint: %e\n',abs(Result1-Result2));