switch animation(mi).field_type 
    case 'e'
        minvalx = min(min(min(Ex)));    
        minvaly = min(min(min(Ey)));    
        minvalz = min(min(min(Ez)));
        maxvalx = max(max(max(Ex)));    
        maxvaly = max(max(max(Ey)));    
        maxvalz = max(max(max(Ez)));
    case 'h'
        minvalx = min(min(min(Hx)));    
        minvaly = min(min(min(Hy)));    
        minvalz = min(min(min(Hz)));
        maxvalx = max(max(max(Hx)));    
        maxvaly = max(max(max(Hy)));    
        maxvalz = max(max(max(Hz)));
    case 'j'
        minvalx = min(min(min(Jx)));    
        minvaly = min(min(min(Jy)));    
        minvalz = min(min(min(Jz)));
        maxvalx = max(max(max(Jx)));    
        maxvaly = max(max(max(Jy)));    
        maxvalz = max(max(max(Jz)));
    case 'm'
        minvalx = min(min(min(Mx)));    
        minvaly = min(min(min(My)));    
        minvalz = min(min(min(Mz)));
        maxvalx = max(max(max(Mx)));    
        maxvaly = max(max(max(My)));    
        maxvalz = max(max(max(Mz)));
end
    minval  = min([minvalx minvaly minvalz]);
    maxval  = max([maxvalx maxvaly maxvalz]);
    maxval = 1.1*max([abs(minval) abs(maxval)]);
    if animation(current_animation_index).component == 'm'
        caxis([0 maxval]);
    else
        caxis([-maxval maxval]);
    end
