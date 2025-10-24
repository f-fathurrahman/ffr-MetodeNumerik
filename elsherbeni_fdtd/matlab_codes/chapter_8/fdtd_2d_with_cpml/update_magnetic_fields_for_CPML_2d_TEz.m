% TEz: apply CPML to magnetic field components
if is_cpml_xn
    for i = 1: n_cpml_xn
        Psi_hzx_xn(i,:,:) = cpml_b_mx_xn(i) * Psi_hzx_xn(i,:,:) ...
            + cpml_a_mx_xn(i)*(Ey(i+1,:,:)-Ey(i,:,:)); 
    end
        Hz(1:n_cpml_xn,:,:) = Hz(1:n_cpml_xn,:,:) ...
            + CPsi_hzx_xn(:,:,:) .* Psi_hzx_xn(:,:,:);    
end

if is_cpml_xp
    n_st = nx - n_cpml_xp;
    for i = 1:n_cpml_xp
        Psi_hzx_xp(i,:,:) = cpml_b_mx_xp(i) * Psi_hzx_xp(i,:,:) ...
            + cpml_a_mx_xp(i)*(Ey(i+n_st+1,:,:)-Ey(i+n_st,:,:)); 
    end

        Hz(n_st+1:nx,:,:) = Hz(n_st+1:nx,:,:) ...
            + CPsi_hzx_xp(:,:,:) .* Psi_hzx_xp(:,:,:);    
end

if is_cpml_yn
    for i = 1:n_cpml_yn
        Psi_hzy_yn(:,i,:) = cpml_b_my_yn(i) * Psi_hzy_yn(:,i,:) ...
            + cpml_a_my_yn(i)*(Ex(:,i+1,:)-Ex(:,i,:)); 
    end
        Hz(:,1:n_cpml_yn,:) = Hz(:,1:n_cpml_yn,:) ...
            + CPsi_hzy_yn(:,:,:) .* Psi_hzy_yn(:,:,:);
end

if is_cpml_yp
    n_st = ny - n_cpml_yp;
    for i = 1:n_cpml_yp
        Psi_hzy_yp(:,i,:) = cpml_b_my_yp(i) * Psi_hzy_yp(:,i,:) ...
            + cpml_a_my_yp(i)*(Ex(:,i+n_st+1,:)-Ex(:,i+n_st,:)); 
    end

        Hz(:,n_st+1:ny,:) = Hz(:,n_st+1:ny,:) ...
            + CPsi_hzy_yp(:,:,:) .* Psi_hzy_yp(:,:,:);
end

