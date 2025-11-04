disp('initializing thin wire updating coefficients');

dtm = dt/mu_0;

for ind = 1:number_of_thin_wires
    ni = thin_wires(ind).node_indices;
    is = ni.is; js = ni.js; ks = ni.ks; 
    ie = ni.ie; je = ni.je; ke = ni.ke; 

    r_o = thin_wires(ind).radius;
    
    switch (thin_wires(ind).direction(1))
        case 'x'
            
            csx = cell_sizes_xe(is:ie-1);
            csy = cell_sizes_ye(js);
            csz = cell_sizes_ze(ks);
            [DX,DY,DZ] = ndgrid(csx, csy, csz);            
            DX = reshape(DX,[],1);
            DY = reshape(DY,[],1);
            DZ = reshape(DZ,[],1);

            Cexe (is:ie-1,js,ks) = 0;
            Cexhy(is:ie-1,js,ks) = 0;
            Cexhz(is:ie-1,js,ks) = 0;
            Chyh (is:ie-1,js,ks-1:ks) = 1;
            Chyez(is:ie-1,js,ks-1:ks) = dtm ...
                ./ (mu_r_y(is:ie-1,js,ks-1:ks) .* DX);
            Chyex(is:ie-1,js,ks-1:ks) = -2 * dtm ...
                ./ (mu_r_y(is:ie-1,js,ks-1:ks) .* DZ .* log(DZ/r_o));
            Chzh( is:ie-1,js-1:js,ks) = 1;
            Chzex(is:ie-1,js-1:js,ks) = 2 * dtm ...
                ./ (mu_r_z(is:ie-1,js-1:js,ks) .* DY .* log(DY/r_o));
            Chzey(is:ie-1,js-1:js,ks) = -dtm ...
                ./ (mu_r_z(is:ie-1,js-1:js,ks) .* DX);
        case 'y'
            csx = cell_sizes_xe(is);
            csy = cell_sizes_ye(js:je-1);
            csz = cell_sizes_ze(ks);
            [DX,DY,DZ] = ndgrid(csx, csy, csz);            
            DX = reshape(DX,[],1);
            DY = reshape(DY,[],1);
            DZ = reshape(DZ,[],1);
            
            Ceye (is,js:je-1,ks) = 0;
            Ceyhx(is,js:je-1,ks) = 0;
            Ceyhz(is,js:je-1,ks) = 0;
            Chzh (is-1:is,js:je-1,ks) = 1;
            Chzex(is-1:is,js:je-1,ks) = dtm ...
                ./ (mu_r_z(is-1:is,js:je-1,ks) .* DY);
            Chzey(is-1:is,js:je-1,ks) = -2 * dtm ...
                ./ (mu_r_z(is-1:is,js:je-1,ks) .* DX .* log(DX/r_o));
            Chxh (is,js:je-1,ks-1:ks) = 1;
            Chxey(is,js:je-1,ks-1:ks) = 2 * dtm ...
                ./ (mu_r_x(is,js:je-1,ks-1:ks) .* DZ .* log(DZ/r_o));
            Chxez(is,js:je-1,ks-1:ks) = -dtm ...
                ./ (mu_r_x(is,js:je-1,ks-1:ks) .* DY);
        case 'z'
            
            csx = cell_sizes_xe(is);
            csy = cell_sizes_ye(js);
            csz = cell_sizes_ze(ks:ke-1);
            [DX,DY,DZ] = ndgrid(csx, csy, csz);            
            DX = reshape(DX,[],1);
            DY = reshape(DY,[],1);
            DZ = reshape(DZ,[],1);

            Ceze (is,js,ks:ke-1) = 0;
            Cezhx(is,js,ks:ke-1) = 0;
            Cezhy(is,js,ks:ke-1) = 0;
            Chxh (is,js-1:js,ks:ke-1) = 1;
            Chxey(is,js-1:js,ks:ke-1) = dtm ...
                ./ (mu_r_x(is,js-1:js,ks:ke-1) .* DZ);
            Chxez(is,js-1:js,ks:ke-1) = -2 * dtm ...
                ./ (mu_r_x(is,js-1:js,ks:ke-1) .* DY .* log(DY/r_o));
            Chyh (is-1:is,js,ks:ke-1) = 1;
            Chyez(is-1:is,js,ks:ke-1) = 2 * dtm ...
                ./ (mu_r_y(is-1:is,js,ks:ke-1) .* DX .* log(DX/r_o));
            Chyex(is-1:is,js,ks:ke-1) = -dtm ...
                ./ (mu_r_y(is-1:is,js,ks:ke-1) .* DZ);
    end
end
