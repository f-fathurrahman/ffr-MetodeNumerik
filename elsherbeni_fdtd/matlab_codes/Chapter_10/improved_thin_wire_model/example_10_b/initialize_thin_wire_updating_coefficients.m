disp('initializing thin wire updating coefficients');

dtm = dt/mu_0;

for ind = 1:number_of_thin_wires
    is = round((thin_wires(ind).min_x - fdtd_domain.min_x)/dx)+1;
    js = round((thin_wires(ind).min_y - fdtd_domain.min_y)/dy)+1;
    ks = round((thin_wires(ind).min_z - fdtd_domain.min_z)/dz)+1;
    ie = round((thin_wires(ind).max_x - fdtd_domain.min_x)/dx)+1;
    je = round((thin_wires(ind).max_y - fdtd_domain.min_y)/dy)+1;
    ke = round((thin_wires(ind).max_z - fdtd_domain.min_z)/dz)+1;
    r_o = thin_wires(ind).radius;
    
    switch (thin_wires(ind).direction(1))
        case 'x'
            khy = (dz/dy)*atan(dy/dz);
            khz = (dy/dz)*atan(dz/dy);
            Cexe (is:ie-1,js,ks) = 0;
            Cexhy(is:ie-1,js,ks) = 0;
            Cexhz(is:ie-1,js,ks) = 0;
            Chyh (is:ie-1,js,ks-1:ks) = 1;
            Chyex(is:ie-1,js,ks-1:ks) = -2 * dtm * khy...
                ./ (mu_r_y(is:ie-1,js,ks-1:ks) * dz * log(dz/r_o));
            Chzh( is:ie-1,js-1:js,ks) = 1;
            Chzex(is:ie-1,js-1:js,ks) = 2 * dtm * khz...
                ./ (mu_r_z(is:ie-1,js-1:js,ks) * dy * log(dy/r_o));
        case 'y'
            khz = (dx/dz)*atan(dz/dx);
            khx = (dz/dx)*atan(dx/dz);
            Ceye (is,js:je-1,ks) = 0;
            Ceyhx(is,js:je-1,ks) = 0;
            Ceyhz(is,js:je-1,ks) = 0;
            Chzh (is-1:is,js:je-1,ks) = 1;
            Chzey(is-1:is,js:je-1,ks) = -2 * dtm * khz...
                ./ (mu_r_z(is-1:is,js:je-1,ks) * dx * log(dx/r_o));
            Chxh (is,js:je-1,ks-1:ks) = 1;
            Chxey(is,js:je-1,ks-1:ks) = 2 * dtm *khx...
                ./ (mu_r_x(is,js:je-1,ks-1:ks) * dz * log(dz/r_o));
        case 'z'
            khx = (dy/dx)*atan(dx/dy);
            khy = (dx/dy)*atan(dy/dx);
            Ceze (is,js,ks:ke-1) = 0;
            Cezhx(is,js,ks:ke-1) = 0;
            Cezhy(is,js,ks:ke-1) = 0;
            Chxh (is,js-1:js,ks:ke-1) = 1;
            Chxez(is,js-1:js,ks:ke-1) = -2 * dtm * khx...
                ./ (mu_r_x(is,js-1:js,ks:ke-1) * dy * log(dy/r_o));
            Chyh (is-1:is,js,ks:ke-1) = 1;
            Chyez(is-1:is,js,ks:ke-1) = 2 * dtm * khy...
                ./ (mu_r_y(is-1:is,js,ks:ke-1) * dx * log(dx/r_o));
    end
end
