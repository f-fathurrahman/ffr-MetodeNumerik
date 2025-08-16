function d = dcombined(P, comb)
    %   Compute signed distance function for a composite object
    %   Original code: DISTMESH 2004-2012 Per-Olof Persson

    Ptemp(:, 1) = P(:, 1)-comb.X(1);
    Ptemp(:, 2) = P(:, 2)-comb.Y(1); 
    if strcmp(comb.type(1), 'c')            
        d = dcircle(Ptemp, comb.R(1));
    else
        d = drectangle(Ptemp, comb.L(1), comb.W(1));
    end    
    
    for m = 2:length(comb.type)
        Ptemp(:, 1) = P(:, 1)-comb.X(m);
        Ptemp(:, 2) = P(:, 2)-comb.Y(m);        
        if strcmp(comb.type(m), 'c')
            dtemp = dcircle(Ptemp, comb.R(m));
        else
            dtemp = drectangle(Ptemp, comb.L(m), comb.W(m));  
        end
        if strcmp(comb.boolean(m), 'u')
            d = min(d, dtemp);
        end
        if strcmp(comb.boolean(m), 'd')            
            d = max(d, -dtemp);          
        end
        if strcmp(comb.boolean(m), 'i')
            d = max(d, dtemp);
        end        
    end    
end