disp('initializing the output parameters');

number_of_sampled_electric_fields  = 0;
number_of_sampled_magnetic_fields  = 0;
number_of_sampled_voltages  = 0;
number_of_sampled_currents  = 0;
number_of_ports             = 0;

if exist('sampled_electric_fields','var')
    number_of_sampled_electric_fields  = size(sampled_electric_fields,2);
end
if exist('sampled_magnetic_fields','var')
    number_of_sampled_magnetic_fields  = size(sampled_magnetic_fields,2);
end
if exist('sampled_voltages','var')
    number_of_sampled_voltages  = size(sampled_voltages,2);
end
if exist('sampled_currents','var')
    number_of_sampled_currents  = size(sampled_currents,2);
end
if exist('ports','var')
    number_of_ports = size(ports,2);
end

% intialize frequency domain parameters
frequency_domain.frequencies = [frequency_domain.start: ...
    frequency_domain.step:frequency_domain.end];
frequency_domain.number_of_frequencies = ...
    size(frequency_domain.frequencies,2);

% initialize sampled electric field terms
for ind=1:number_of_sampled_electric_fields
    node_indices = get_node_indices(sampled_electric_fields(ind), fdtd_domain);
    sampled_electric_fields(ind).is = node_indices.is;
    sampled_electric_fields(ind).js = node_indices.js;
    sampled_electric_fields(ind).ks = node_indices.ks;
    sampled_electric_fields(ind).sampled_value = ...
        zeros(1, number_of_time_steps);
    sampled_electric_fields(ind).time = ([1:number_of_time_steps])*dt;
end

% initialize sampled magnetic field terms
for ind=1:number_of_sampled_magnetic_fields  
    node_indices = get_node_indices(sampled_magnetic_fields(ind), fdtd_domain);
    sampled_magnetic_fields(ind).is = node_indices.is;
    sampled_magnetic_fields(ind).js = node_indices.js;
    sampled_magnetic_fields(ind).ks = node_indices.ks;
    sampled_magnetic_fields(ind).sampled_value = ...
        zeros(1, number_of_time_steps);
    sampled_magnetic_fields(ind).time = ([1:number_of_time_steps]-0.5)*dt;
end

% initialize sampled voltage terms
for ind=1:number_of_sampled_voltages 
    node_indices = get_node_indices(sampled_voltages(ind), fdtd_domain);
    is = node_indices.is; js = node_indices.js; ks = node_indices.ks; 
    ie = node_indices.ie; je = node_indices.je; ke = node_indices.ke; 
    sampled_voltages(ind).node_indices = node_indices;
    
    sampled_voltages(ind).sampled_value = ...
                        zeros(1, number_of_time_steps);

    switch (sampled_voltages(ind).direction(1))
    case 'x'
        fi = create_linear_index_list(Ex,is:ie-1,js:je,ks:ke);
        csx = cell_sizes_xe(is:ie-1);
        csy = cell_sizes_yh(js:je);
        csz = cell_sizes_zh(ks:ke);
        [DX,DY,DZ] = ndgrid(csx, csy, csz);            
        DX = reshape(DX,[],1);
        sampled_voltages(ind).Csvf = -DX/((je-js+1)*(ke-ks+1));
    case 'y'
        fi = create_linear_index_list(Ey,is:ie,js:je-1,ks:ke);

        csx = cell_sizes_xh(is:ie);
        csy = cell_sizes_ye(js:je-1);
        csz = cell_sizes_zh(ks:ke);
        [DX,DY,DZ] = ndgrid(csx, csy, csz);            
        DY = reshape(DY,[],1);
        
        sampled_voltages(ind).Csvf = -DY/((ke-ks+1)*(ie-is+1));
    case 'z'
        fi = create_linear_index_list(Ez,is:ie,js:je,ks:ke-1);
        
        csx = cell_sizes_xh(is:ie);
        csy = cell_sizes_yh(js:je);
        csz = cell_sizes_ze(ks:ke-1);
        [DX,DY,DZ] = ndgrid(csx, csy, csz);            
        DZ = reshape(DZ,[],1);

        sampled_voltages(ind).Csvf = -DZ/((ie-is+1)*(je-js+1));
    end    
    if strcmp(sampled_voltages(ind).direction(2),'n')
        sampled_voltages(ind).Csvf =  -1 * sampled_voltages(ind).Csvf;
    end
    sampled_voltages(ind).field_indices = fi;
    sampled_voltages(ind).time = ([1:number_of_time_steps])*dt;
end

% initialize sampled current terms
for ind=1:number_of_sampled_currents  
    node_indices = get_node_indices(sampled_currents(ind), fdtd_domain);
    is = node_indices.is; js = node_indices.js; ks = node_indices.ks; 
    ie = node_indices.ie; je = node_indices.je; ke = node_indices.ke; 

    sampled_currents(ind).sampled_value = ...
                        zeros(1, number_of_time_steps);
    sampled_currents(ind).time =([1:number_of_time_steps]-0.5)*dt;

    sgn = 1;
    if strcmp(sampled_currents(ind).direction(2),'n')
        sgn = -1;
    end    
    
    switch (sampled_currents(ind).direction(1))
    case 'x'
        if (is==ie) 
            is = is-1;  node_indices.is = is;
            ie = ie+1;  node_indices.ie = ie;
        end

        sampled_currents(ind).fihy1 = create_linear_index_list(Hy,is:ie-1,js:je,ks-1);
        sampled_currents(ind).fihy2 = create_linear_index_list(Hy,is:ie-1,js:je,ke);
        [DX, DY, DZ] = ndgrid(cell_sizes_xe, cell_sizes_yh, cell_sizes_ze); 
        Cschy = DY(is:ie-1,js:je,ke)/(ie-is);
        sampled_currents(ind).Cschy = sgn * reshape(Cschy,[],1);
        
        sampled_currents(ind).fihz1 = create_linear_index_list(Hz,is:ie-1,je,ks:ke);
        sampled_currents(ind).fihz2 = create_linear_index_list(Hz,is:ie-1,js-1,ks:ke);

        [DX, DY, DZ] = ndgrid(cell_sizes_xe, cell_sizes_ye, cell_sizes_zh); 
        Cschz = DZ(is:ie-1,je,ks:ke)/(ie-is); 
        sampled_currents(ind).Cschz = sgn * reshape(Cschz,[],1);
        
    case 'y'
        if (js==je)
            js = js-1;  node_indices.js = js;
            je = je+1;  node_indices.je = je;
        end

        sampled_currents(ind).fihz1 = create_linear_index_list(Hz,is-1,js:je-1,ks:ke);
        sampled_currents(ind).fihz2 = create_linear_index_list(Hz,ie,js:je-1,ks:ke);
        [DX, DY, DZ] = ndgrid(cell_sizes_xe, cell_sizes_ye, cell_sizes_zh); 
        Cschz = DZ(is,js:je-1,ks:ke)/(je-js);        
        sampled_currents(ind).Cschz = sgn * reshape(Cschz,[],1);
        
        sampled_currents(ind).fihx1 = create_linear_index_list(Hx,is:ie,js:je-1,ke);
        sampled_currents(ind).fihx2 = create_linear_index_list(Hx,is:ie,js:je-1,ks-1);
        [DX, DY, DZ] = ndgrid(cell_sizes_xh, cell_sizes_ye, cell_sizes_ze); 
        Cschx = DX(is:ie,js:je-1,ks)/(je-js);        
        sampled_currents(ind).Cschx = sgn * reshape(Cschx,[],1);
        
    case 'z'
        if (ks==ke) 
            ks = ks-1;  node_indices.ks = ks;
            ke = ke+1;  node_indices.ke = ke;
        end
        
        sampled_currents(ind).fihx1 = create_linear_index_list(Hx,is:ie,js-1,ks:ke-1);
        sampled_currents(ind).fihx2 = create_linear_index_list(Hx,is:ie,je,ks:ke-1);
        [DX, DY, DZ] = ndgrid(cell_sizes_xh, cell_sizes_ye, cell_sizes_ze); 
        Cschx = DX(is:ie,je,ks:ke-1)/(ke-ks);        
        sampled_currents(ind).Cschx = sgn * reshape(Cschx,[],1);
        
        sampled_currents(ind).fihy1 = create_linear_index_list(Hy,ie,js:je,ks:ke-1);
        sampled_currents(ind).fihy2 = create_linear_index_list(Hy,is-1,js:je,ks:ke-1);
        [DX, DY, DZ] = ndgrid(cell_sizes_xe, cell_sizes_yh, cell_sizes_ze); 
        Cschy = DY(ie,js:je,ks:ke-1)/(ke-ks);        
        sampled_currents(ind).Cschy = sgn * reshape(Cschy,[],1);
    end
    sampled_currents(ind).node_indices = node_indices;    
end
