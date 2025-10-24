% TMz
% apply CPML to magnetic field components
if is_cpml_xn
    for i = 1: n_cpml_xn
        Psi_hyx_xn(i,:,:) = cpml_b_mx_xn(i) * Psi_hyx_xn(i,:,:) ...
            + cpml_a_mx_xn(i)*(Ez(i+1,:,:)-Ez(i,:,:)); 
    end
        Hy(1:n_cpml_xn,:,:) = Hy(1:n_cpml_xn,:,:) ...
            + CPsi_hyx_xn(:,:,:) .* Psi_hyx_xn(:,:,:);
end

if is_cpml_xp
    n_st = nx - n_cpml_xp;
    for i = 1:n_cpml_xp
        Psi_hyx_xp(i,:,:) = cpml_b_mx_xp(i) * Psi_hyx_xp(i,:,:) ...
            + cpml_a_mx_xp(i)*(Ez(i+n_st+1,:,:)-Ez(i+n_st,:,:)); 
    end

        Hy(n_st+1:nx,:,:) = Hy(n_st+1:nx,:,:) ...
            + CPsi_hyx_xp(:,:,:) .* Psi_hyx_xp(:,:,:);
end

if is_cpml_yn
    for i = 1:n_cpml_yn
        Psi_hxy_yn(:,i,:) = cpml_b_my_yn(i) * Psi_hxy_yn(:,i,:) ...
            + cpml_a_my_yn(i)*(Ez(:,i+1,:)-Ez(:,i,:)); 
    end
        Hx(:,1:n_cpml_yn,:) = Hx(:,1:n_cpml_yn,:) ...
            + CPsi_hxy_yn(:,:,:) .* Psi_hxy_yn(:,:,:);    
end

if is_cpml_yp
    n_st = ny - n_cpml_yp;
    for i = 1:n_cpml_yp
        Psi_hxy_yp(:,i,:) = cpml_b_my_yp(i) * Psi_hxy_yp(:,i,:) ...
            + cpml_a_my_yp(i)*(Ez(:,i+n_st+1,:)-Ez(:,i+n_st,:)); 
    end
        Hx(:,n_st+1:ny,:) = Hx(:,n_st+1:ny,:) ...
            + CPsi_hxy_yp(:,:,:) .* Psi_hxy_yp(:,:,:);    
end

