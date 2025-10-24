% Capturing magnetic fields

for ind=1:number_of_sampled_magnetic_fields
    is = sampled_magnetic_fields(ind).is;
    js = sampled_magnetic_fields(ind).js;
    ks = sampled_magnetic_fields(ind).ks;

    switch (sampled_magnetic_fields(ind).component)
        case 'x'
            sampled_value = 0.25 * sum(sum(Hx(is,js-1:js,ks-1:ks))); 
        case 'y'
            sampled_value = 0.25 * sum(sum(Hy(is-1:is,js,ks-1:ks))); 
        case 'z'
            sampled_value = 0.25 * sum(sum(Hz(is-1:is,js-1:js,ks))); 
        case 'm'
            svx = 0.25 * sum(sum(Hx(is,js-1:js,ks-1:ks))); 
            svy = 0.25 * sum(sum(Hy(is-1:is,js,ks-1:ks))); 
            svz = 0.25 * sum(sum(Hz(is-1:is,js-1:js,ks))); 
            sampled_value = sqrt(svx^2 + svy^2 + svz^2);
    end
    sampled_magnetic_fields(ind).sampled_value(time_step) = ...
        sampled_value;
end

