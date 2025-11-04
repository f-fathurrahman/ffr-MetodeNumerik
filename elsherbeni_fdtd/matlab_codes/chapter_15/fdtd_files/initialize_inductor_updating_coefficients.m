disp('initializing inductor updating coefficients');

for ind = 1:number_of_inductors
    ni = inductors(ind).node_indices;
    is = ni.is; js = ni.js; ks = ni.ks; 
    ie = ni.ie; je = ni.je; ke = ni.ke; 

    L = inductors(ind).inductance_per_component;
        
    switch (inductors(ind).direction(1))
    case 'x'
        fi = create_linear_index_list(eps_r_x,is:ie-1,js:je,ks:ke);
        
        csx = cell_sizes_xe(is:ie-1);
        csy = cell_sizes_yh(js:je);
        csz = cell_sizes_zh(ks:ke);
        [DX,DY,DZ] = ndgrid(csx, csy, csz);            
        DX = reshape(DX,[],1);
        DY = reshape(DY,[],1);
        DZ = reshape(DZ,[],1);

        inductors(ind).Cexj = -(2*dt) ...
            ./ (2*eps_r_x(fi)*eps_0+dt*sigma_e_x(fi));      
        inductors(ind).Jix = zeros(size(fi));
        inductors(ind).Cjex = (dt*DX)./(L*DY.*DZ);      
    case 'y'
        fi = create_linear_index_list(eps_r_y,is:ie,js:je-1,ks:ke);
     
        csx = cell_sizes_xh(is:ie);
        csy = cell_sizes_ye(js:je-1);
        csz = cell_sizes_zh(ks:ke);
        [DX,DY,DZ] = ndgrid(csx, csy, csz);            
        DX = reshape(DX,[],1);
        DY = reshape(DY,[],1);
        DZ = reshape(DZ,[],1);
        
        inductors(ind).Ceyj = -(2*dt) ...
            ./ (2*eps_r_y(fi)*eps_0+dt*sigma_e_y(fi));      
        inductors(ind).Jiy = zeros(size(fi));
        inductors(ind).Cjey = (dt*DY)./(L*DZ.*DX);      
    case 'z'
        fi = create_linear_index_list(eps_r_z,is:ie,js:je,ks:ke-1);
        
        csx = cell_sizes_xh(is:ie);
        csy = cell_sizes_yh(js:je);
        csz = cell_sizes_ze(ks:ke-1);
        [DX,DY,DZ] = ndgrid(csx, csy, csz);            
        DX = reshape(DX,[],1);
        DY = reshape(DY,[],1);
        DZ = reshape(DZ,[],1);
        
        inductors(ind).Cezj = -(2*dt) ...
            ./ (2*eps_r_z(fi)*eps_0+dt*sigma_e_z(fi));      
        inductors(ind).Jiz = zeros(size(fi));
        inductors(ind).Cjez = (dt*DZ)./(L*DX.*DY);      
    end
    inductors(ind).field_indices = fi;
end
