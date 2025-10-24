% Capturing electric fields

for ind=1:number_of_sampled_electric_fields
    is = sampled_electric_fields(ind).is;
    js = sampled_electric_fields(ind).js;
    ks = sampled_electric_fields(ind).ks;

    switch (sampled_electric_fields(ind).component)
        case 'x'
            sampled_value = 0.5 * sum(Ex(is-1:is,js,ks)); 
        case 'y'
            sampled_value = 0.5 * sum(Ey(is,js-1:js,ks)); 
        case 'z'
            sampled_value = 0.5 * sum(Ez(is,js,ks-1:ks)); 
        case 'm'
            svx = 0.5 * sum(Ex(is-1:is,js,ks)); 
            svy = 0.5 * sum(Ey(is,js-1:js,ks)); 
            svz = 0.5 * sum(Ez(is,js,ks-1:ks)); 
            sampled_value = sqrt(svx^2 + svy^2 + svz^2);
    end
    sampled_electric_fields(ind).sampled_value(time_step) = sampled_value;
end

