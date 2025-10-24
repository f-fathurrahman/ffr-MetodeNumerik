% update electric fields at the PML regions
% TEz
% apply CPML to electric field components
if is_cpml_xn
    for i = 1:n_cpml_xn
        Psi_eyx_xn(i,:,:) = cpml_b_ex_xn(i) * Psi_eyx_xn(i,:,:) ...
            + cpml_a_ex_xn(i)*(Hz(i+1,:,:)-Hz(i,:,:)); 
    end
    Ey(2:n_cpml_xn+1,:,:) = Ey(2:n_cpml_xn+1,:,:) ...
            + CPsi_eyx_xn .* Psi_eyx_xn;
end

if is_cpml_xp 
    n_st = nx - n_cpml_xp;
    for i = 1:n_cpml_xp
        Psi_eyx_xp(i,:,:) = cpml_b_ex_xp(i) * Psi_eyx_xp(i,:,:) ...
            + cpml_a_ex_xp(i)*(Hz(i+n_st,:,:)-Hz(i+n_st-1,:,:)); 
    end
    Ey(n_st+1:nx,:,:) = Ey(n_st+1:nx,:,:) ...
        + CPsi_eyx_xp .* Psi_eyx_xp;
end

if is_cpml_yn
    for i = 1:n_cpml_yn
        Psi_exy_yn(:,i,:) = cpml_b_ey_yn(i) * Psi_exy_yn(:,i,:) ...
            + cpml_a_ey_yn(i)*(Hz(:,i+1,:)-Hz(:,i,:)); 
    end
    Ex(:,2:n_cpml_yn+1,:) = Ex(:,2:n_cpml_yn+1,:) ...
        + CPsi_exy_yn .* Psi_exy_yn;    
end

if is_cpml_yp
    n_st = ny - n_cpml_yp;
    for i = 1:n_cpml_yp
        Psi_exy_yp(:,i,:) = cpml_b_ey_yp(i) * Psi_exy_yp(:,i,:) ...
            + cpml_a_ey_yp(i)*(Hz(:,i+n_st,:)-Hz(:,i+n_st-1,:)); 
    end
    Ex(:,n_st+1:ny,:) = Ex(:,n_st+1:ny,:) ...
        + CPsi_exy_yp .* Psi_exy_yp;    
end
