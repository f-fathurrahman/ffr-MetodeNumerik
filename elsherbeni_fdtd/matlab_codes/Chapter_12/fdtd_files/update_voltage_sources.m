% updating electric field components 
% associated with the voltage sources

for ind = 1:number_of_voltage_sources
    fi = voltage_sources(ind).field_indices;
    switch (voltage_sources(ind).direction(1))
    case 'x' 
        Ex(fi) = Ex(fi) + voltage_sources(ind).Cexs ...
            * voltage_sources(ind).voltage_per_e_field(time_step);
    case 'y'
        Ey(fi) = Ey(fi) + voltage_sources(ind).Ceys ...
            * voltage_sources(ind).voltage_per_e_field(time_step);
    case 'z'
        Ez(fi) = Ez(fi) + voltage_sources(ind).Cezs ...
            * voltage_sources(ind).voltage_per_e_field(time_step);
    end
end
