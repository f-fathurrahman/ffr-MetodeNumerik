% Capturing sampled voltages

for ind=1:number_of_sampled_voltages
    fi   = sampled_voltages(ind).field_indices;
    Csvf = sampled_voltages(ind).Csvf;
    switch (sampled_voltages(ind).direction(1))
    case 'x'
        sampled_value = Csvf * sum(Ex(fi)); 
    case 'y'
        sampled_value = Csvf * sum(Ey(fi)); 
    case 'z'
        sampled_value = Csvf * sum(Ez(fi)); 
    end
    sampled_voltages(ind).sampled_value(time_step) = sampled_value;
end

