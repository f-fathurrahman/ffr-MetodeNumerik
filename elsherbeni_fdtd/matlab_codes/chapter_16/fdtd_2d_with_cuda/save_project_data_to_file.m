% save project data to a file
% an executable reads the data from thid file and
% runs the time marching loop on GPU

% extend the domain and adjust number of cells for gpu
nxx = (floor(nx/16)+1)*16; 
nyy = (floor(ny/16)+1)*16; 

fid = fopen(exe_data_file_name, 'w');

fwrite(fid, computation_platform, 'int32');

fwrite(fid, number_of_time_steps, 'int32');

fwrite(fid, plotting_step, 'int32');

fwrite(fid,dt, 'float32');
fwrite(fid,dx, 'float32');
fwrite(fid,dy, 'float32');
fwrite(fid,nx, 'int32');
fwrite(fid,ny, 'int32');
fwrite(fid,nxx, 'int32');
fwrite(fid,nyy, 'int32');

fwrite(fid,fdtd_domain.min_x, 'float32');
fwrite(fid,fdtd_domain.min_y, 'float32');
fwrite(fid,fdtd_domain.max_x, 'float32');
fwrite(fid,fdtd_domain.max_y, 'float32');

fwrite(fid, number_of_circles, 'int32');
for i=1:number_of_circles
    fwrite(fid,circles(i).center_x, 'float32');
    fwrite(fid,circles(i).center_y, 'float32');
    fwrite(fid,circles(i).radius, 'float32');    
end
fwrite(fid, number_of_rectangles, 'int32');
for i=1:number_of_rectangles
    fwrite(fid,rectangles(i).min_x, 'float32');
    fwrite(fid,rectangles(i).min_y, 'float32');
    fwrite(fid,rectangles(i).max_x, 'float32');
    fwrite(fid,rectangles(i).max_y, 'float32');
end

fwrite(fid, is_TEz, 'bit8');
fwrite(fid, is_TMz, 'bit8');

Ceze(1,:) = 0;
Ceze(:,1) = 0;
Ceze(nxp1,:) = 0;
Ceze(:,nyp1) = 0;

Cezhx(1,:) = 0;
Cezhx(:,1) = 0;
Cezhx(nxp1,:) = 0;
Cezhx(:,nyp1) = 0;

Cezhy(1,:) = 0;
Cezhy(:,1) = 0;
Cezhy(nxp1,:) = 0;
Cezhy(:,nyp1) = 0;

Cexe(:,1) = 0;
Cexe(:,nyp1) = 0;
Cexhz(:,1) = 0;
Cexhz(:,nyp1) = 0;

Ceye(1,:) = 0;
Ceye(nxp1,:) = 0;
Ceyhz(1,:) = 0;
Ceyhz(nxp1,:) = 0;

fwrite(fid,append_pads(Ex, nxx, nyy+2, 0, 1),  'float32');
fwrite(fid,append_pads(Ey, nxx, nyy+2, 0, 1),  'float32');
fwrite(fid,append_pads(Ez, nxx, nyy+2, 0, 1),  'float32');
fwrite(fid,append_pads(Hx, nxx, nyy+2, 0, 1),  'float32');
fwrite(fid,append_pads(Hy, nxx, nyy+2, 0, 1),  'float32');
fwrite(fid,append_pads(Hz, nxx, nyy+2, 0, 1),  'float32');

fwrite(fid,append_pads(Cexe , nxx, nyy+2, 0, 1),  'float32');
fwrite(fid,append_pads(Cexhz, nxx, nyy+2, 0, 1),  'float32');
fwrite(fid,append_pads(Ceye , nxx, nyy+2, 0, 1),  'float32');
fwrite(fid,append_pads(Ceyhz, nxx, nyy+2, 0, 1),  'float32');
fwrite(fid,append_pads(Ceze , nxx, nyy+2, 0, 1),  'float32');
fwrite(fid,append_pads(Cezhy, nxx, nyy+2, 0, 1),  'float32');
fwrite(fid,append_pads(Cezhx, nxx, nyy+2, 0, 1),  'float32');

fwrite(fid,append_pads(Chxh , nxx, nyy+2, 0, 1),  'float32');
fwrite(fid,append_pads(Chxez, nxx, nyy+2, 0, 1),  'float32');
fwrite(fid,append_pads(Chyh , nxx, nyy+2, 0, 1),  'float32');
fwrite(fid,append_pads(Chyez, nxx, nyy+2, 0, 1),  'float32');
fwrite(fid,append_pads(Chzh , nxx, nyy+2, 0, 1),  'float32');
fwrite(fid,append_pads(Chzex, nxx, nyy+2, 0, 1),  'float32');
fwrite(fid,append_pads(Chzey, nxx, nyy+2, 0, 1),  'float32');

fwrite(fid,is_any_side_cpml, 'bit8');
fwrite(fid,is_cpml_xn, 'bit8');
fwrite(fid,is_cpml_xp, 'bit8');
fwrite(fid,is_cpml_yn, 'bit8');
fwrite(fid,is_cpml_yp, 'bit8');

if is_any_side_cpml
    fwrite(fid,n_cpml_xn, 'int32');
    fwrite(fid,n_cpml_xp, 'int32');
    fwrite(fid,n_cpml_yn, 'int32');
    fwrite(fid,n_cpml_yp, 'int32');
end
if is_cpml_xn
    fwrite(fid,append_pads(cpml_a_ex_xn, 1, 16, 0, 1), 'float32');
    fwrite(fid,append_pads(cpml_b_ex_xn, 1, 16, 0, 1), 'float32');
    fwrite(fid,append_pads(cpml_a_mx_xn, 1, 16, 0, 0), 'float32');
    fwrite(fid,append_pads(cpml_b_mx_xn, 1, 16, 0, 0), 'float32');
    if is_TEz
        fwrite(fid,append_pads( Psi_hzx_xn, 16, nyy+2, 0, 1), 'float32');
        fwrite(fid,append_pads(CPsi_hzx_xn, 16, nyy+2, 0, 1), 'float32');
        fwrite(fid,append_pads( Psi_eyx_xn, 16, nyy+2, 1, 1), 'float32');
        fwrite(fid,append_pads(CPsi_eyx_xn, 16, nyy+2, 1, 1), 'float32');
    end
    if is_TMz
        fwrite(fid,append_pads( Psi_ezx_xn, 16, nyy+2, 1, 1), 'float32');
        fwrite(fid,append_pads(CPsi_ezx_xn, 16, nyy+2, 1, 1), 'float32');
        fwrite(fid,append_pads( Psi_hyx_xn, 16, nyy+2, 0, 1), 'float32');
        fwrite(fid,append_pads(CPsi_hyx_xn, 16, nyy+2, 0, 1), 'float32');
    end
end

if is_cpml_xp
    npad = nxx-nx;
    nshift = 32-n_cpml_xp-npad;   
    fwrite(fid,append_pads(cpml_a_ex_xp, 1, 32, 0, nshift), 'float32');
    fwrite(fid,append_pads(cpml_b_ex_xp, 1, 32, 0, nshift), 'float32');
    fwrite(fid,append_pads(cpml_a_mx_xp, 1, 32, 0, nshift), 'float32');
    fwrite(fid,append_pads(cpml_b_mx_xp, 1, 32, 0, nshift), 'float32');
    if is_TEz
        fwrite(fid,append_pads( Psi_hzx_xp, 32, nyy+2, nshift, 1), 'float32');
        fwrite(fid,append_pads(CPsi_hzx_xp, 32, nyy+2, nshift, 1), 'float32');
        fwrite(fid,append_pads( Psi_eyx_xp, 32, nyy+2, nshift, 1), 'float32');
        fwrite(fid,append_pads(CPsi_eyx_xp, 32, nyy+2, nshift, 1), 'float32');
    end
    if is_TMz
        fwrite(fid,append_pads( Psi_ezx_xp, 32, nyy+2, nshift, 1), 'float32');
        fwrite(fid,append_pads(CPsi_ezx_xp, 32, nyy+2, nshift, 1), 'float32');
        fwrite(fid,append_pads( Psi_hyx_xp, 32, nyy+2, nshift, 1), 'float32');
        fwrite(fid,append_pads(CPsi_hyx_xp, 32, nyy+2, nshift, 1), 'float32');
    end    
end
if is_cpml_yn
    fwrite(fid,append_pads(cpml_a_ey_yn, 1, 16, 0, 1), 'float32');
    fwrite(fid,append_pads(cpml_b_ey_yn, 1, 16, 0, 1), 'float32');
    fwrite(fid,append_pads(cpml_a_my_yn, 1, 16, 0, 0), 'float32');
    fwrite(fid,append_pads(cpml_b_my_yn, 1, 16, 0, 0), 'float32');
    if is_TEz
        fwrite(fid,append_pads( Psi_hzy_yn, nxx, 16, 0, 0), 'float32');
        fwrite(fid,append_pads(CPsi_hzy_yn, nxx, 16, 0, 0), 'float32');
        fwrite(fid,append_pads( Psi_exy_yn, nxx, 16, 0, 1), 'float32');
        fwrite(fid,append_pads(CPsi_exy_yn, nxx, 16, 0, 1), 'float32');
    end
    if is_TMz
        fwrite(fid,append_pads( Psi_hxy_yn, nxx, 16, 0, 0), 'float32');
        fwrite(fid,append_pads(CPsi_hxy_yn, nxx, 16, 0, 0), 'float32');
        fwrite(fid,append_pads( Psi_ezy_yn, nxx, 16, 0, 1), 'float32');
        fwrite(fid,append_pads(CPsi_ezy_yn, nxx, 16, 0, 1), 'float32');
    end
end
if is_cpml_yp
    npad = nyy-ny;
    nshift = 32-n_cpml_yp-npad;   
    fwrite(fid,append_pads(cpml_a_ey_yp, 1, 32, 0, nshift), 'float32');
    fwrite(fid,append_pads(cpml_b_ey_yp, 1, 32, 0, nshift), 'float32');
    fwrite(fid,append_pads(cpml_a_my_yp, 1, 32, 0, nshift), 'float32');
    fwrite(fid,append_pads(cpml_b_my_yp, 1, 32, 0, nshift), 'float32');
    if is_TEz
        fwrite(fid,append_pads( Psi_hzy_yp, nxx, 32, 0, nshift), 'float32');
        fwrite(fid,append_pads(CPsi_hzy_yp, nxx, 32, 0, nshift), 'float32');
        fwrite(fid,append_pads( Psi_exy_yp, nxx, 32, 0, nshift), 'float32');
        fwrite(fid,append_pads(CPsi_exy_yp, nxx, 32, 0, nshift), 'float32');
    end
    if is_TMz
        fwrite(fid,append_pads( Psi_hxy_yp, nxx, 32, 0, nshift), 'float32');
        fwrite(fid,append_pads(CPsi_hxy_yp, nxx, 32, 0, nshift), 'float32');
        fwrite(fid,append_pads( Psi_ezy_yp, nxx, 32, 0, nshift), 'float32');
        fwrite(fid,append_pads(CPsi_ezy_yp, nxx, 32, 0, nshift), 'float32');
    end
end

fwrite(fid,number_of_sampled_electric_fields, 'int32');
fwrite(fid,number_of_sampled_magnetic_fields, 'int32');
fwrite(fid,number_of_impressed_J, 'int32');
fwrite(fid,number_of_impressed_M, 'int32');

for ind = 1:number_of_sampled_electric_fields
    fwrite(fid,sampled_electric_fields(ind).component, 'char');
    fwrite(fid,sampled_electric_fields(ind).is, 'int32');
    fwrite(fid,sampled_electric_fields(ind).js, 'int32');
    fwrite(fid,sampled_electric_fields(ind).sampled_value, 'float32');
end

for ind = 1:number_of_sampled_magnetic_fields
    fwrite(fid,sampled_magnetic_fields(ind).component, 'char');
    fwrite(fid,sampled_magnetic_fields(ind).is, 'int32');
    fwrite(fid,sampled_magnetic_fields(ind).js, 'int32');
    fwrite(fid,sampled_magnetic_fields(ind).sampled_value, 'float32');
end

number_of_cej_components = [];

for ind = 1:number_of_impressed_J
    switch (impressed_J(ind).direction(1))
    case 'x'
        cej = impressed_J(ind).Cexj;
    case 'y'
        cej = impressed_J(ind).Ceyj;
    case 'z'
        cej = impressed_J(ind).Cezj;
    end
    sz = size(cej);
    number_of_cej_components(ind) = sz(1)*sz(2);
end

fwrite(fid,number_of_cej_components, 'int32');
    
for ind = 1:number_of_impressed_J
    fwrite(fid,impressed_J(ind).direction(1), 'char');
    fwrite(fid,impressed_J(ind).is, 'int32');
    fwrite(fid,impressed_J(ind).ie, 'int32');
    fwrite(fid,impressed_J(ind).js, 'int32');
    fwrite(fid,impressed_J(ind).je, 'int32');
    fwrite(fid,impressed_J(ind).min_x, 'float32');
    fwrite(fid,impressed_J(ind).min_y, 'float32');
    fwrite(fid,impressed_J(ind).max_x, 'float32');
    fwrite(fid,impressed_J(ind).max_y, 'float32');
    
    fwrite(fid,impressed_J(ind).waveform, 'float32');
    switch (impressed_J(ind).direction(1))
    case 'x'
        fwrite(fid,impressed_J(ind).Cexj, 'float32');
    case 'y'
        fwrite(fid,impressed_J(ind).Ceyj, 'float32');
    case 'z'
        fwrite(fid,impressed_J(ind).Cezj, 'float32');
    end
end

number_of_cem_components = [];

for ind = 1:number_of_impressed_M
    switch (impressed_M(ind).direction(1))
    case 'x'
        cem = impressed_M(ind).Chxm;
    case 'y'
        cem = impressed_M(ind).Chym;
    case 'z'
        cem = impressed_M(ind).Chzm;
    end
    sz = size(cem);
    number_of_cem_components(ind) = sz(1)*sz(2);
end

fwrite(fid,number_of_cem_components, 'int32');
    
for ind = 1:number_of_impressed_M
    fwrite(fid,impressed_M(ind).direction(1), 'char');
    fwrite(fid,impressed_M(ind).is, 'int32');
    fwrite(fid,impressed_M(ind).ie, 'int32');
    fwrite(fid,impressed_M(ind).js, 'int32');
    fwrite(fid,impressed_M(ind).je, 'int32');
    fwrite(fid,impressed_M(ind).min_x, 'float32');
    fwrite(fid,impressed_M(ind).min_y, 'float32');
    fwrite(fid,impressed_M(ind).max_x, 'float32');
    fwrite(fid,impressed_M(ind).max_y, 'float32');

    fwrite(fid,impressed_M(ind).waveform, 'float32');
    switch (impressed_M(ind).direction(1))
    case 'x'
        fwrite(fid,impressed_M(ind).Chxm, 'float32');
    case 'y'
        fwrite(fid,impressed_M(ind).Chym, 'float32');
    case 'z'
        fwrite(fid,impressed_M(ind).Chzm, 'float32');
    end
end

show_Hz = false;
show_Ez = false;
if is_TEz
    for mi=1:number_of_sampled_transient_H_planes
        if sampled_transient_H_planes(mi).component == 'z'
            show_Hz = true;
        end
    end
end
if is_TMz
    for mi=1:number_of_sampled_transient_E_planes
        if sampled_transient_E_planes(mi).component == 'z'
            show_Ez = true;
        end
    end
end

fwrite(fid, show_Hz, 'bit8');
fwrite(fid, show_Ez, 'bit8');
    
fclose(fid);

