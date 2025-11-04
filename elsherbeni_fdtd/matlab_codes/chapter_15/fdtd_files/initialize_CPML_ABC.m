% Initialize CPML boundary condition    

p_order = boundary.cpml_order; % order of the polynomial distribution
sigma_ratio = boundary.cpml_sigma_factor;
kappa_max = boundary.cpml_kappa_max;
alpha_min = boundary.cpml_alpha_min;
alpha_max = boundary.cpml_alpha_max;

% Initialize cpml for xn region
if is_cpml_xn

    % define one-dimensional temporary cpml parameter arrays 
    sigma_max = sigma_ratio  * (p_order+1)/(150*pi*dx);
    ncells = n_cpml_xn;
    rho_e = ([ncells:-1:1]-0.75)/ncells;
    rho_m = ([ncells:-1:1]-0.25)/ncells;
    sigma_pex_xn = sigma_max * rho_e.^p_order;
    sigma_pmx_xn = sigma_max * rho_m.^p_order;
    sigma_pmx_xn = (mu_0/eps_0) * sigma_pmx_xn;
    kappa_ex_xn = 1 + (kappa_max - 1) * rho_e.^p_order;
    kappa_mx_xn = 1 + (kappa_max - 1) * rho_m.^p_order;
    alpha_ex_xn = alpha_min + (alpha_max - alpha_min) * (1-rho_e);
    alpha_mx_xn = alpha_min + (alpha_max - alpha_min) * (1-rho_m);
    alpha_mx_xn = (mu_0/eps_0) * alpha_mx_xn;
    
    % define one-dimensional cpml parameter arrays 
    cpml_b_ex_xn = exp((-dt/eps_0) ...
        *((sigma_pex_xn./kappa_ex_xn)+ alpha_ex_xn)); 
    cpml_a_ex_xn = (1/dx)*(cpml_b_ex_xn-1.0).* sigma_pex_xn ...
        ./(kappa_ex_xn.*(sigma_pex_xn+kappa_ex_xn.*alpha_ex_xn));
    cpml_b_mx_xn = exp((-dt/mu_0) ...
        *((sigma_pmx_xn./kappa_mx_xn)+ alpha_mx_xn)); 
    cpml_a_mx_xn = (1/dx)*(cpml_b_mx_xn-1.0) .* sigma_pmx_xn ...
        ./(kappa_mx_xn.*(sigma_pmx_xn+kappa_mx_xn.*alpha_mx_xn));

    % Create and initialize two-dimensional cpml convolution parameters 
    Psi_eyx_xn = zeros(ncells,ny,nzp1); 
    Psi_ezx_xn = zeros(ncells,nyp1,nz); 
    Psi_hyx_xn = zeros(ncells,nyp1,nz); 
    Psi_hzx_xn = zeros(ncells,ny,nzp1); 

    % Create and initialize two-dimensional cpml convolution coefficients 
    % Notice that Ey(1,:,:) and Ez(1,:,:) are not updated by cmpl 
    CPsi_eyx_xn = Ceyhz(2:ncells+1,:,:)*dx;
    CPsi_ezx_xn = Cezhy(2:ncells+1,:,:)*dx;
    CPsi_hyx_xn = Chyez(1:ncells,:,:)*dx;
    CPsi_hzx_xn = Chzey(1:ncells,:,:)*dx;
    
    % Adjust FDTD coefficients in the CPML region 
    % Notice that Ey(1,:,:) and Ez(1,:,:) are not updated by cmpl 
    for i = 1: ncells                                           
        Ceyhz(i+1,:,:) = Ceyhz(i+1,:,:)/kappa_ex_xn(i);
        Cezhy(i+1,:,:) = Cezhy(i+1,:,:)/kappa_ex_xn(i);
        Chyez(i,:,:) = Chyez(i,:,:)/kappa_mx_xn(i);
        Chzey(i,:,:) = Chzey(i,:,:)/kappa_mx_xn(i);
    end
    
    % Delete temporary arrays. These arrays will not be used any more.
    clear sigma_pex_xn sigma_pmx_xn;
    clear kappa_ex_xn kappa_mx_xn;
    clear alpha_ex_xn alpha_mx_xn;
end

% Initialize cpml for xp region
if is_cpml_xp

    % define one-dimensional temporary cpml parameter arrays 
    sigma_max = sigma_ratio  * (p_order+1)/(150*pi*dx);
    ncells = n_cpml_xp;
    rho_e = ([1:ncells]-0.75)/ncells;
    rho_m = ([1:ncells]-0.25)/ncells;
    sigma_pex_xp = sigma_max * rho_e.^p_order;
    sigma_pmx_xp = sigma_max * rho_m.^p_order;
    sigma_pmx_xp = (mu_0/eps_0) * sigma_pmx_xp;
    kappa_ex_xp = 1 + (kappa_max - 1) * rho_e.^p_order;
    kappa_mx_xp = 1 + (kappa_max - 1) * rho_m.^p_order;
    alpha_ex_xp = alpha_min + (alpha_max - alpha_min) * (1-rho_e);
    alpha_mx_xp = alpha_min + (alpha_max - alpha_min) * (1-rho_m);
    alpha_mx_xp = (mu_0/eps_0) * alpha_mx_xp;
    
    % define one-dimensional cpml parameter arrays 
    cpml_b_ex_xp = exp((-dt/eps_0) ...
        *((sigma_pex_xp./kappa_ex_xp)+ alpha_ex_xp)); 
    cpml_a_ex_xp = (1/dx)*(cpml_b_ex_xp-1.0).* sigma_pex_xp ...
        ./(kappa_ex_xp.*(sigma_pex_xp+kappa_ex_xp.*alpha_ex_xp));
    cpml_b_mx_xp = exp((-dt/mu_0) ...
        *((sigma_pmx_xp./kappa_mx_xp)+ alpha_mx_xp)); 
    cpml_a_mx_xp = (1/dx)*(cpml_b_mx_xp-1.0) .* sigma_pmx_xp ...
        ./(kappa_mx_xp.*(sigma_pmx_xp+kappa_mx_xp.*alpha_mx_xp));

    % Create and initialize two-dimensional cpml convolution parameters 
    Psi_eyx_xp = zeros(ncells,ny,nzp1); 
    Psi_ezx_xp = zeros(ncells,nyp1,nz); 
    Psi_hyx_xp = zeros(ncells,nyp1,nz); 
    Psi_hzx_xp = zeros(ncells,ny,nzp1); 

    % Create and initialize two-dimensional cpml convolution coefficients 
    % Notice that Ey(nxp1,:,:) and Ez(nxp1,:,:) are not updated by cmpl 
    CPsi_eyx_xp = Ceyhz(nxp1-ncells:nx,:,:)*dx;
    CPsi_ezx_xp = Cezhy(nxp1-ncells:nx,:,:)*dx;
    CPsi_hyx_xp = Chyez(nxp1-ncells:nx,:,:)*dx;
    CPsi_hzx_xp = Chzey(nxp1-ncells:nx,:,:)*dx;
    
    % Adjust FDTD coefficients in the CPML region 
    % Notice that Ey(nxp1,:,:) and Ez(nxp1,:,:) are not updated by cmpl 
    for i = 1: ncells                                           
        Ceyhz(nx-ncells+i,:,:) = Ceyhz(nx-ncells+i,:,:)/kappa_ex_xp(i);
        Cezhy(nx-ncells+i,:,:) = Cezhy(nx-ncells+i,:,:)/kappa_ex_xp(i);
        Chyez(nx-ncells+i,:,:) = Chyez(nx-ncells+i,:,:)/kappa_mx_xp(i);
        Chzey(nx-ncells+i,:,:) = Chzey(nx-ncells+i,:,:)/kappa_mx_xp(i);
    end
    
    % Delete temporary arrays. These arrays will not be used any more.
    clear sigma_pex_xp sigma_pmx_xp;
    clear kappa_ex_xp kappa_mx_xp;
    clear alpha_ex_xp alpha_mx_xp;
end

% Initialize cpml for yn region
if is_cpml_yn

    % define one-dimensional temporary cpml parameter arrays 
    sigma_max = sigma_ratio  * (p_order+1)/(150*pi*dy);
    ncells = n_cpml_yn;
    rho_e = ([ncells:-1:1]-0.75)/ncells;
    rho_m = ([ncells:-1:1]-0.25)/ncells;
    sigma_pey_yn = sigma_max * rho_e.^p_order;
    sigma_pmy_yn = sigma_max * rho_m.^p_order;
    sigma_pmy_yn = (mu_0/eps_0) * sigma_pmy_yn;
    kappa_ey_yn = 1 + (kappa_max - 1) * rho_e.^p_order;
    kappa_my_yn = 1 + (kappa_max - 1) * rho_m.^p_order;
    alpha_ey_yn = alpha_min + (alpha_max - alpha_min) * (1-rho_e);
    alpha_my_yn = alpha_min + (alpha_max - alpha_min) * (1-rho_m);
    alpha_my_yn = (mu_0/eps_0) * alpha_my_yn;
    
    % define one-dimensional cpml parameter arrays 
    cpml_b_ey_yn = exp((-dt/eps_0) ...
        *((sigma_pey_yn./kappa_ey_yn)+ alpha_ey_yn)); 
    cpml_a_ey_yn = (1/dy)*(cpml_b_ey_yn-1.0).* sigma_pey_yn ...
        ./(kappa_ey_yn.*(sigma_pey_yn+kappa_ey_yn.*alpha_ey_yn));
    cpml_b_my_yn = exp((-dt/mu_0) ...
        *((sigma_pmy_yn./kappa_my_yn)+ alpha_my_yn)); 
    cpml_a_my_yn = (1/dy)*(cpml_b_my_yn-1.0) .* sigma_pmy_yn ...
        ./(kappa_my_yn.*(sigma_pmy_yn+kappa_my_yn.*alpha_my_yn));

    % Create and initialize two-dimensional cpml convolution parameters 
    Psi_ezy_yn = zeros(nxp1,ncells,nz); 
    Psi_exy_yn = zeros(nx,ncells,nzp1); 
    Psi_hzy_yn = zeros(nx,ncells,nzp1); 
    Psi_hxy_yn = zeros(nxp1,ncells,nz); 

    % Create and initialize two-dimensional cpml convolution coefficients 
    % Notice that Ez(:,1,:) and Ex(:,1,:) are not updated by cmpl 
    CPsi_ezy_yn = Cezhx(:,2:ncells+1,:)*dy;
    CPsi_exy_yn = Cexhz(:,2:ncells+1,:)*dy;
    CPsi_hzy_yn = Chzex(:,1:ncells,:)*dy;
    CPsi_hxy_yn = Chxez(:,1:ncells,:)*dy;
    
    % Adjust FDTD coefficients in the CPML region 
    % Notice that Ez(:,1,:) and Ex(:,1,:) are not updated by cmpl 
    for i = 1: ncells                                           
        Cezhx(:,i+1,:) = Cezhx(:,i+1,:)/kappa_ey_yn(i);
        Cexhz(:,i+1,:) = Cexhz(:,i+1,:)/kappa_ey_yn(i);
        Chzex(:,i,:) = Chzex(:,i,:)/kappa_my_yn(i);
        Chxez(:,i,:) = Chxez(:,i,:)/kappa_my_yn(i);
    end
    
    % Delete temporary arrays. These arrays will not be used any more.
    clear sigma_pey_yn sigma_pmy_yn;
    clear kappa_ey_yn kappa_my_yn;
    clear alpha_ey_yn alpha_my_yn;
end

% Initialize cpml for yp region
if is_cpml_yp

    % define one-dimensional temporary cpml parameter arrays 
    sigma_max = sigma_ratio  * (p_order+1)/(150*pi*dy);
    ncells = n_cpml_yp;
    rho_e = ([1:ncells]-0.75)/ncells;
    rho_m = ([1:ncells]-0.25)/ncells;
    sigma_pey_yp = sigma_max * rho_e.^p_order;
    sigma_pmy_yp = sigma_max * rho_m.^p_order;
    sigma_pmy_yp = (mu_0/eps_0) * sigma_pmy_yp;
    kappa_ey_yp = 1 + (kappa_max - 1) * rho_e.^p_order;
    kappa_my_yp = 1 + (kappa_max - 1) * rho_m.^p_order;
    alpha_ey_yp = alpha_min + (alpha_max - alpha_min) * (1-rho_e);
    alpha_my_yp = alpha_min + (alpha_max - alpha_min) * (1-rho_m);
    alpha_my_yp = (mu_0/eps_0) * alpha_my_yp;
    
    % define one-dimensional cpml parameter arrays 
    cpml_b_ey_yp = exp((-dt/eps_0) ...
        *((sigma_pey_yp./kappa_ey_yp)+ alpha_ey_yp)); 
    cpml_a_ey_yp = (1/dy)*(cpml_b_ey_yp-1.0).* sigma_pey_yp ...
        ./(kappa_ey_yp.*(sigma_pey_yp+kappa_ey_yp.*alpha_ey_yp));
    cpml_b_my_yp = exp((-dt/mu_0) ...
        *((sigma_pmy_yp./kappa_my_yp)+ alpha_my_yp)); 
    cpml_a_my_yp = (1/dy)*(cpml_b_my_yp-1.0) .* sigma_pmy_yp ...
        ./(kappa_my_yp.*(sigma_pmy_yp+kappa_my_yp.*alpha_my_yp));

    % Create and initialize two-dimensional cpml convolution parameters 
    Psi_ezy_yp = zeros(nxp1,ncells,nz); 
    Psi_exy_yp = zeros(nx,ncells,nzp1); 
    Psi_hzy_yp = zeros(nx,ncells,nzp1); 
    Psi_hxy_yp = zeros(nxp1,ncells,nz); 

    % Create and initialize two-dimensional cpml convolution coefficients 
    % Notice that Ez(:,nyp1,:) and Ex(:,nyp1,:) are not updated by cmpl 
    CPsi_ezy_yp = Cezhx(:,nyp1-ncells:ny,:)*dy;
    CPsi_exy_yp = Cexhz(:,nyp1-ncells:ny,:)*dy;
    CPsi_hzy_yp = Chzex(:,nyp1-ncells:ny,:)*dy;
    CPsi_hxy_yp = Chxez(:,nyp1-ncells:ny,:)*dy;
    
    % Adjust FDTD coefficients in the CPML region 
    % Notice that Ez(:,nyp1,:) and Ex(:,nyp1,:) are not updated by cmpl 
    for i = 1: ncells                                           
        Cezhx(:,ny-ncells+i,:) = Cezhx(:,ny-ncells+i,:)/kappa_ey_yp(i);
        Cexhz(:,ny-ncells+i,:) = Cexhz(:,ny-ncells+i,:)/kappa_ey_yp(i);
        Chzex(:,ny-ncells+i,:) = Chzex(:,ny-ncells+i,:)/kappa_my_yp(i);
        Chxez(:,ny-ncells+i,:) = Chxez(:,ny-ncells+i,:)/kappa_my_yp(i);
    end
    
    % Delete temporary arrays. These arrays will not be used any more.
    clear sigma_pey_yp sigma_pmy_yp;
    clear kappa_ey_yp kappa_my_yp;
    clear alpha_ey_yp alpha_my_yp;
end

% Initialize cpml for zn region
if is_cpml_zn

    % define one-dimensional temporary cpml parameter arrays 
    sigma_max = sigma_ratio  * (p_order+1)/(150*pi*dz);
    ncells = n_cpml_zn;
    rho_e = ([ncells:-1:1]-0.75)/ncells;
    rho_m = ([ncells:-1:1]-0.25)/ncells;
    sigma_pez_zn = sigma_max * rho_e.^p_order;
    sigma_pmz_zn = sigma_max * rho_m.^p_order;
    sigma_pmz_zn = (mu_0/eps_0) * sigma_pmz_zn;
    kappa_ez_zn = 1 + (kappa_max - 1) * rho_e.^p_order;
    kappa_mz_zn = 1 + (kappa_max - 1) * rho_m.^p_order;
    alpha_ez_zn = alpha_min + (alpha_max - alpha_min) * (1-rho_e);
    alpha_mz_zn = alpha_min + (alpha_max - alpha_min) * (1-rho_m);
    alpha_mz_zn = (mu_0/eps_0) * alpha_mz_zn;
    
    % define one-dimensional cpml parameter arrays 
    cpml_b_ez_zn = exp((-dt/eps_0) ...
        *((sigma_pez_zn./kappa_ez_zn)+ alpha_ez_zn)); 
    cpml_a_ez_zn = (1/dz)*(cpml_b_ez_zn-1.0).* sigma_pez_zn ...
        ./(kappa_ez_zn.*(sigma_pez_zn+kappa_ez_zn.*alpha_ez_zn));
    cpml_b_mz_zn = exp((-dt/mu_0) ...
        *((sigma_pmz_zn./kappa_mz_zn)+ alpha_mz_zn)); 
    cpml_a_mz_zn = (1/dz)*(cpml_b_mz_zn-1.0) .* sigma_pmz_zn ...
        ./(kappa_mz_zn.*(sigma_pmz_zn+kappa_mz_zn.*alpha_mz_zn));

    % Create and initialize two-dimensional cpml convolution parameters 
    Psi_exz_zn = zeros(nx,nyp1,ncells); 
    Psi_eyz_zn = zeros(nxp1,ny,ncells); 
    Psi_hxz_zn = zeros(nxp1,ny,ncells); 
    Psi_hyz_zn = zeros(nx,nyp1,ncells); 

    % Create and initialize two-dimensional cpml convolution coefficients 
    % Notice that Ex(:,:,1) and Ey(:,:,1) are not updated by cmpl 
    CPsi_exz_zn = Cexhy(:,:,2:ncells+1)*dz;
    CPsi_eyz_zn = Ceyhx(:,:,2:ncells+1)*dz;
    CPsi_hxz_zn = Chxey(:,:,1:ncells)*dz;
    CPsi_hyz_zn = Chyex(:,:,1:ncells)*dz;
    
    % Adjust FDTD coefficients in the CPML region 
    % Notice that Ex(:,:,1) and Ey(:,:,1) are not updated by cmpl 
    for i = 1: ncells                                           
        Cexhy(:,:,i+1) = Cexhy(:,:,i+1)/kappa_ez_zn(i);
        Ceyhx(:,:,i+1) = Ceyhx(:,:,i+1)/kappa_ez_zn(i);
        Chxey(:,:,i) = Chxey(:,:,i)/kappa_mz_zn(i);
        Chyex(:,:,i) = Chyex(:,:,i)/kappa_mz_zn(i);
    end
    
    % Delete temporary arrays. These arrays will not be used any more.
    clear sigma_pez_zn sigma_pmz_zn;
    clear kappa_ez_zn kappa_mz_zn;
    clear alpha_ez_zn alpha_mz_zn;
end

% Initialize cpml for zp region
if is_cpml_zp

    % define one-dimensional temporary cpml parameter arrays 
    sigma_max = sigma_ratio  * (p_order+1)/(150*pi*dz);
    ncells = n_cpml_zp;
    rho_e = ([1:ncells]-0.75)/ncells;
    rho_m = ([1:ncells]-0.25)/ncells;
    sigma_pez_zp = sigma_max * rho_e.^p_order;
    sigma_pmz_zp = sigma_max * rho_m.^p_order;
    sigma_pmz_zp = (mu_0/eps_0) * sigma_pmz_zp;
    kappa_ez_zp = 1 + (kappa_max - 1) * rho_e.^p_order;
    kappa_mz_zp = 1 + (kappa_max - 1) * rho_m.^p_order;
    alpha_ez_zp = alpha_min + (alpha_max - alpha_min) * (1-rho_e);
    alpha_mz_zp = alpha_min + (alpha_max - alpha_min) * (1-rho_m);
    alpha_mz_zp = (mu_0/eps_0) * alpha_mz_zp;
    
    % define one-dimensional cpml parameter arrays 
    cpml_b_ez_zp = exp((-dt/eps_0) ...
        *((sigma_pez_zp./kappa_ez_zp)+ alpha_ez_zp)); 
    cpml_a_ez_zp = (1/dz)*(cpml_b_ez_zp-1.0).* sigma_pez_zp ...
        ./(kappa_ez_zp.*(sigma_pez_zp+kappa_ez_zp.*alpha_ez_zp));
    cpml_b_mz_zp = exp((-dt/mu_0) ...
        *((sigma_pmz_zp./kappa_mz_zp)+ alpha_mz_zp)); 
    cpml_a_mz_zp = (1/dz)*(cpml_b_mz_zp-1.0) .* sigma_pmz_zp ...
        ./(kappa_mz_zp.*(sigma_pmz_zp+kappa_mz_zp.*alpha_mz_zp));

    % Create and initialize two-dimensional cpml convolution parameters 
    Psi_exz_zp = zeros(nx,nyp1,ncells); 
    Psi_eyz_zp = zeros(nxp1,ny,ncells); 
    Psi_hxz_zp = zeros(nxp1,ny,ncells); 
    Psi_hyz_zp = zeros(nx,nyp1,ncells); 

    % Create and initialize two-dimensional cpml convolution coefficients 
    % Notice that Ex(:,:,nzp1) and Ey(:,:,nzp1) are not updated by cmpl 
    CPsi_exz_zp = Cexhy(:,:,nzp1-ncells:nz)*dz;
    CPsi_eyz_zp = Ceyhx(:,:,nzp1-ncells:nz)*dz;
    CPsi_hxz_zp = Chxey(:,:,nzp1-ncells:nz)*dz;
    CPsi_hyz_zp = Chyex(:,:,nzp1-ncells:nz)*dz;
    
    % Adjust FDTD coefficients in the CPML region 
    % Notice that Ex(:,:,nzp1) and Ey(:,:,nzp1) are not updated by cmpl 
    for i = 1: ncells                                           
        Cexhy(:,:,nz-ncells+i) = Cexhy(:,:,nz-ncells+i)/kappa_ez_zp(i);
        Ceyhx(:,:,nz-ncells+i) = Ceyhx(:,:,nz-ncells+i)/kappa_ez_zp(i);
        Chxey(:,:,nz-ncells+i) = Chxey(:,:,nz-ncells+i)/kappa_mz_zp(i);
        Chyex(:,:,nz-ncells+i) = Chyex(:,:,nz-ncells+i)/kappa_mz_zp(i);
    end
    
    % Delete temporary arrays. These arrays will not be used any more.
    clear sigma_pez_zp sigma_pmz_zp;
    clear kappa_ez_zp kappa_mz_zp;
    clear alpha_ez_zp alpha_mz_zp;
end