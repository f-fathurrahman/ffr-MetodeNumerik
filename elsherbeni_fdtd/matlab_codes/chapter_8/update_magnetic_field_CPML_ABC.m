% apply CPML to magnetic field components
if is_cpml_xn
    for i = 1: n_cpml_xn
        Psi_hyx_xn(i,:,:) = cpml_b_mx_xn(i) * Psi_hyx_xn(i,:,:) ...
            + cpml_a_mx_xn(i)*(Ez(i+1,:,:)-Ez(i,:,:)); 
        Psi_hzx_xn(i,:,:) = cpml_b_mx_xn(i) * Psi_hzx_xn(i,:,:) ...
            + cpml_a_mx_xn(i)*(Ey(i+1,:,:)-Ey(i,:,:)); 
    end
        Hy(1:n_cpml_xn,:,:) = Hy(1:n_cpml_xn,:,:) ...
            + CPsi_hyx_xn(:,:,:) .* Psi_hyx_xn(:,:,:);
        Hz(1:n_cpml_xn,:,:) = Hz(1:n_cpml_xn,:,:) ...
            + CPsi_hzx_xn(:,:,:) .* Psi_hzx_xn(:,:,:);    
end

if is_cpml_xp
    n_st = nx - n_cpml_xp;
    for i = 1:n_cpml_xp
        Psi_hyx_xp(i,:,:) = cpml_b_mx_xp(i) * Psi_hyx_xp(i,:,:) ...
            + cpml_a_mx_xp(i)*(Ez(i+n_st+1,:,:)-Ez(i+n_st,:,:)); 
        Psi_hzx_xp(i,:,:) = cpml_b_mx_xp(i) * Psi_hzx_xp(i,:,:) ...
            + cpml_a_mx_xp(i)*(Ey(i+n_st+1,:,:)-Ey(i+n_st,:,:)); 
    end

        Hy(n_st+1:nx,:,:) = Hy(n_st+1:nx,:,:) ...
            + CPsi_hyx_xp(:,:,:) .* Psi_hyx_xp(:,:,:);
        Hz(n_st+1:nx,:,:) = Hz(n_st+1:nx,:,:) ...
            + CPsi_hzx_xp(:,:,:) .* Psi_hzx_xp(:,:,:);    
end

if is_cpml_yn
    for i = 1:n_cpml_yn
        Psi_hzy_yn(:,i,:) = cpml_b_my_yn(i) * Psi_hzy_yn(:,i,:) ...
            + cpml_a_my_yn(i)*(Ex(:,i+1,:)-Ex(:,i,:)); 
        Psi_hxy_yn(:,i,:) = cpml_b_my_yn(i) * Psi_hxy_yn(:,i,:) ...
            + cpml_a_my_yn(i)*(Ez(:,i+1,:)-Ez(:,i,:)); 
    end
        Hz(:,1:n_cpml_yn,:) = Hz(:,1:n_cpml_yn,:) ...
            + CPsi_hzy_yn(:,:,:) .* Psi_hzy_yn(:,:,:);
        Hx(:,1:n_cpml_yn,:) = Hx(:,1:n_cpml_yn,:) ...
            + CPsi_hxy_yn(:,:,:) .* Psi_hxy_yn(:,:,:);    
end

if is_cpml_yp
    n_st = ny - n_cpml_yp;
    for i = 1:n_cpml_yp
        Psi_hzy_yp(:,i,:) = cpml_b_my_yp(i) * Psi_hzy_yp(:,i,:) ...
            + cpml_a_my_yp(i)*(Ex(:,i+n_st+1,:)-Ex(:,i+n_st,:)); 
        Psi_hxy_yp(:,i,:) = cpml_b_my_yp(i) * Psi_hxy_yp(:,i,:) ...
            + cpml_a_my_yp(i)*(Ez(:,i+n_st+1,:)-Ez(:,i+n_st,:)); 
    end

        Hz(:,n_st+1:ny,:) = Hz(:,n_st+1:ny,:) ...
            + CPsi_hzy_yp(:,:,:) .* Psi_hzy_yp(:,:,:);
        Hx(:,n_st+1:ny,:) = Hx(:,n_st+1:ny,:) ...
            + CPsi_hxy_yp(:,:,:) .* Psi_hxy_yp(:,:,:);    
end

if is_cpml_zn
    for i = 1:n_cpml_zn
        Psi_hxz_zn(:,:,i) = cpml_b_mz_zn(i) * Psi_hxz_zn(:,:,i) ...
            + cpml_a_mz_zn(i)*(Ey(:,:,i+1)-Ey(:,:,i)); 
        Psi_hyz_zn(:,:,i) = cpml_b_mz_zn(i) * Psi_hyz_zn(:,:,i) ...
            + cpml_a_mz_zn(i)*(Ex(:,:,i+1)-Ex(:,:,i)); 
    end
        Hx(:,:,1:n_cpml_zn) = Hx(:,:,1:n_cpml_zn) ...
            + CPsi_hxz_zn(:,:,:) .* Psi_hxz_zn(:,:,:);
        Hy(:,:,1:n_cpml_zn) = Hy(:,:,1:n_cpml_zn) ...
            + CPsi_hyz_zn(:,:,:) .* Psi_hyz_zn(:,:,:);    
end

if is_cpml_zp
    n_st = nz - n_cpml_zp;
    for i = 1:n_cpml_zp
        Psi_hxz_zp(:,:,i) = cpml_b_mz_zp(i) * Psi_hxz_zp(:,:,i) ...
            + cpml_a_mz_zp(i)*(Ey(:,:,i+n_st+1)-Ey(:,:,i+n_st)); 
        Psi_hyz_zp(:,:,i) = cpml_b_mz_zp(i) * Psi_hyz_zp(:,:,i) ...
            + cpml_a_mz_zp(i)*(Ex(:,:,i+n_st+1)-Ex(:,:,i+n_st)); 
    end
        Hx(:,:,n_st+1:nz) = Hx(:,:,n_st+1:nz) ...
            + CPsi_hxz_zp(:,:,:) .* Psi_hxz_zp(:,:,:);
        Hy(:,:,n_st+1:nz) = Hy(:,:,n_st+1:nz) ...
            + CPsi_hyz_zp(:,:,:) .* Psi_hyz_zp(:,:,:);    
end
