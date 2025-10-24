% apply CPML to electric field components
if is_cpml_xn
    for i = 1:n_cpml_xn
        Psi_eyx_xn(i,:,:) = cpml_b_ex_xn(i) * Psi_eyx_xn(i,:,:) ...
            + cpml_a_ex_xn(i)*(Hz(i+1,:,:)-Hz(i,:,:)); 
        Psi_ezx_xn(i,:,:) = cpml_b_ex_xn(i) * Psi_ezx_xn(i,:,:) ...
            + cpml_a_ex_xn(i)*(Hy(i+1,:,:)-Hy(i,:,:)); 
    end
    Ey(2:n_cpml_xn+1,:,:) = Ey(2:n_cpml_xn+1,:,:) ...
            + CPsi_eyx_xn .* Psi_eyx_xn;
    Ez(2:n_cpml_xn+1,:,:) = Ez(2:n_cpml_xn+1,:,:) ...
            + CPsi_ezx_xn .* Psi_ezx_xn;    
end

if is_cpml_xp 
    n_st = nx - n_cpml_xp;
    for i = 1:n_cpml_xp
        Psi_eyx_xp(i,:,:) = cpml_b_ex_xp(i) * Psi_eyx_xp(i,:,:) ...
            + cpml_a_ex_xp(i)*(Hz(i+n_st,:,:)-Hz(i+n_st-1,:,:)); 
        Psi_ezx_xp(i,:,:) = cpml_b_ex_xp(i) * Psi_ezx_xp(i,:,:) ...
            + cpml_a_ex_xp(i)*(Hy(i+n_st,:,:)-Hy(i+n_st-1,:,:)); 
    end
    Ey(n_st+1:nx,:,:) = Ey(n_st+1:nx,:,:) ...
        + CPsi_eyx_xp .* Psi_eyx_xp;
    Ez(n_st+1:nx,:,:) = Ez(n_st+1:nx,:,:) ...
        + CPsi_ezx_xp .* Psi_ezx_xp;    
end

if is_cpml_yn
    for i = 1:n_cpml_yn
        Psi_ezy_yn(:,i,:) = cpml_b_ey_yn(i) * Psi_ezy_yn(:,i,:) ...
            + cpml_a_ey_yn(i)*(Hx(:,i+1,:)-Hx(:,i,:)); 
        Psi_exy_yn(:,i,:) = cpml_b_ey_yn(i) * Psi_exy_yn(:,i,:) ...
            + cpml_a_ey_yn(i)*(Hz(:,i+1,:)-Hz(:,i,:)); 
    end
    Ez(:,2:n_cpml_yn+1,:) = Ez(:,2:n_cpml_yn+1,:) ...
        + CPsi_ezy_yn .* Psi_ezy_yn;
    Ex(:,2:n_cpml_yn+1,:) = Ex(:,2:n_cpml_yn+1,:) ...
        + CPsi_exy_yn .* Psi_exy_yn;    
end

if is_cpml_yp
    n_st = ny - n_cpml_yp;
    for i = 1:n_cpml_yp
        Psi_ezy_yp(:,i,:) = cpml_b_ey_yp(i) * Psi_ezy_yp(:,i,:) ...
            + cpml_a_ey_yp(i)*(Hx(:,i+n_st,:)-Hx(:,i+n_st-1,:)); 
        Psi_exy_yp(:,i,:) = cpml_b_ey_yp(i) * Psi_exy_yp(:,i,:) ...
            + cpml_a_ey_yp(i)*(Hz(:,i+n_st,:)-Hz(:,i+n_st-1,:)); 
    end
    Ez(:,n_st+1:ny,:) = Ez(:,n_st+1:ny,:) ...
        + CPsi_ezy_yp .* Psi_ezy_yp;
    Ex(:,n_st+1:ny,:) = Ex(:,n_st+1:ny,:) ...
        + CPsi_exy_yp .* Psi_exy_yp;    
end

if is_cpml_zn
    for i = 1:n_cpml_zn
        Psi_exz_zn(:,:,i) = cpml_b_ez_zn(i) * Psi_exz_zn(:,:,i) ...
            + cpml_a_ez_zn(i)*(Hy(:,:,i+1)-Hy(:,:,i)); 
        Psi_eyz_zn(:,:,i) = cpml_b_ez_zn(i) * Psi_eyz_zn(:,:,i) ...
            + cpml_a_ez_zn(i)*(Hx(:,:,i+1)-Hx(:,:,i)); 
    end
    Ex(:,:,2:n_cpml_zn+1) = Ex(:,:,2:n_cpml_zn+1) ...
        + CPsi_exz_zn .* Psi_exz_zn;
    Ey(:,:,2:n_cpml_zn+1) = Ey(:,:,2:n_cpml_zn+1) ...
        + CPsi_eyz_zn .* Psi_eyz_zn;    
end

if is_cpml_zp 
    n_st = nz - n_cpml_zp;
    for i = 1:n_cpml_zp
        Psi_exz_zp(:,:,i) = cpml_b_ez_zp(i) * Psi_exz_zp(:,:,i) ...
            + cpml_a_ez_zp(i)*(Hy(:,:,i+n_st)-Hy(:,:,i+n_st-1)); 
        Psi_eyz_zp(:,:,i) = cpml_b_ez_zp(i) * Psi_eyz_zp(:,:,i) ...
            + cpml_a_ez_zp(i)*(Hx(:,:,i+n_st)-Hx(:,:,i+n_st-1)); 
    end
    Ex(:,:,n_st+1:nz) = Ex(:,:,n_st+1:nz) ...
        + CPsi_exz_zp .* Psi_exz_zp;
    Ey(:,:,n_st+1:nz) = Ey(:,:,n_st+1:nz) ...
        + CPsi_eyz_zp .* Psi_eyz_zp;    
end
