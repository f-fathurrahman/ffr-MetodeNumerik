% updating electric field components 
% associated with the diodes

for ind = 1:number_of_diodes
    fi = diodes(ind).field_indices;  
    B  = diodes(ind).B;
    switch (diodes(ind).direction(1))
    case 'x' 
        E = diodes(ind).Exn;    
        C = -Ex(fi) + diodes(ind).Cexd;
        A = -diodes(ind).Cexd * exp(B * E);
        E = solve_diode_equation(A, B, C, E);
        Ex(fi) = E;
        diodes(ind).Exn = E;
    case 'y' 
        E = diodes(ind).Eyn;    
        C = -Ey(fi) + diodes(ind).Ceyd;
        A = -diodes(ind).Ceyd * exp(B * E);
        E = solve_diode_equation(A, B, C, E);
        Ey(fi) = E;
        diodes(ind).Eyn = E;
    case 'z' 
        E = diodes(ind).Ezn;    
        C = -Ez(fi) + diodes(ind).Cezd;
        A = -diodes(ind).Cezd * exp(B * E);
        E = solve_diode_equation(A, B, C, E);
        Ez(fi) = E;
        diodes(ind).Ezn = E;
    end
end
