% Capturing sampled currents

for ind=1:number_of_sampled_currents
    ni = sampled_currents(ind).node_indices;
    is = ni.is; js = ni.js; ks = ni.ks; 
    ie = ni.ie; je = ni.je; ke = ni.ke; 

    switch (sampled_currents(ind).direction(1))
    case 'x'
        fihy1 = sampled_currents(ind).fihy1;
        fihy2 = sampled_currents(ind).fihy2;
        fihz1 = sampled_currents(ind).fihz1;
        fihz2 = sampled_currents(ind).fihz2;
        sampled_value = ... 
            sum(sampled_currents(ind).Cschy .* (Hy(fihy1) - Hy(fihy2))) ...
          + sum(sampled_currents(ind).Cschz .* (Hz(fihz1) - Hz(fihz2)));        
    case 'y'
        fihz1 = sampled_currents(ind).fihz1;
        fihz2 = sampled_currents(ind).fihz2;
        fihx1 = sampled_currents(ind).fihx1;
        fihx2 = sampled_currents(ind).fihx2;
        sampled_value = ... 
            sum(sampled_currents(ind).Cschz .* (Hz(fihz1) - Hz(fihz2))) ...
          + sum(sampled_currents(ind).Cschx .* (Hx(fihx1) - Hx(fihx2)));
    case 'z'
        fihx1 = sampled_currents(ind).fihx1;
        fihx2 = sampled_currents(ind).fihx2;
        fihy1 = sampled_currents(ind).fihy1;
        fihy2 = sampled_currents(ind).fihy2;
        sampled_value = ... 
            sum(sampled_currents(ind).Cschx .* (Hx(fihx1) - Hx(fihx2))) ...
          + sum(sampled_currents(ind).Cschy .* (Hy(fihy1) - Hy(fihy2)));
    end
    sampled_currents(ind).sampled_value(time_step) = sampled_value;
end

