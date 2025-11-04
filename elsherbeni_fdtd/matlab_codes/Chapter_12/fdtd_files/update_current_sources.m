% updating electric field components 
% associated with the current sources

for ind = 1:number_of_current_sources
    fi = current_sources(ind).field_indices;
    switch (current_sources(ind).direction(1))
    case 'x' 
        Ex(fi) = Ex(fi) + current_sources(ind).Cexs ...
            * current_sources(ind).current_per_e_field(time_step);
    case 'y'
        Ey(fi) = Ey(fi) + current_sources(ind).Ceys ...
            * current_sources(ind).current_per_e_field(time_step);
    case 'z'
        Ez(fi) = Ez(fi) + current_sources(ind).Cezs ...
            * current_sources(ind).current_per_e_field(time_step);
    end
end
