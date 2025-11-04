% Capturing sampled voltages

for ind=1:number_of_sampled_voltages
    fi   = sampled_voltages(ind).field_indices;
    Csvf = sampled_voltages(ind).Csvf;
    switch (sampled_voltages(ind).direction(1))
    case 'x'
        sampled_value = sum(Csvf .* Ex(fi)); 
    case 'y'
        sampled_value = sum(Csvf .* Ey(fi)); 
    case 'z'
        sampled_value = sum(Csvf .* Ez(fi)); 
    end
    sampled_voltages(ind).sampled_value(time_step) = sampled_value;
end

