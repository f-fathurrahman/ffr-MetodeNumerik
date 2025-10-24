% update electric fields at the PML regions
% TMz
% apply CPML to electric field components
if is_cpml_xn
    for i = 1:n_cpml_xn
        Psi_ezx_xn(i,:,:) = cpml_b_ex_xn(i) * Psi_ezx_xn(i,:,:) ...
            + cpml_a_ex_xn(i)*(Hy(i+1,:,:)-Hy(i,:,:)); 
    end
    Ez(2:n_cpml_xn+1,:,:) = Ez(2:n_cpml_xn+1,:,:) ...
            + CPsi_ezx_xn .* Psi_ezx_xn;    
end

if is_cpml_xp 
    n_st = nx - n_cpml_xp;
    for i = 1:n_cpml_xp
        Psi_ezx_xp(i,:,:) = cpml_b_ex_xp(i) * Psi_ezx_xp(i,:,:) ...
            + cpml_a_ex_xp(i)*(Hy(i+n_st,:,:)-Hy(i+n_st-1,:,:)); 
    end
    Ez(n_st+1:nx,:,:) = Ez(n_st+1:nx,:,:) ...
        + CPsi_ezx_xp .* Psi_ezx_xp;    
end

if is_cpml_yn
    for i = 1:n_cpml_yn
        Psi_ezy_yn(:,i,:) = cpml_b_ey_yn(i) * Psi_ezy_yn(:,i,:) ...
            + cpml_a_ey_yn(i)*(Hx(:,i+1,:)-Hx(:,i,:)); 
    end
    Ez(:,2:n_cpml_yn+1,:) = Ez(:,2:n_cpml_yn+1,:) ...
        + CPsi_ezy_yn .* Psi_ezy_yn;
end

if is_cpml_yp
    n_st = ny - n_cpml_yp;
    for i = 1:n_cpml_yp
        Psi_ezy_yp(:,i,:) = cpml_b_ey_yp(i) * Psi_ezy_yp(:,i,:) ...
            + cpml_a_ey_yp(i)*(Hx(:,i+n_st,:)-Hx(:,i+n_st-1,:)); 
    end
    Ez(:,n_st+1:ny,:) = Ez(:,n_st+1:ny,:) ...
        + CPsi_ezy_yp .* Psi_ezy_yp;
end
