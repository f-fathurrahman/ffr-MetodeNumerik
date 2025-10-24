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

    % Create and initialize cpml convolution parameters 
    Psi_eyx_xn = zeros(ncells,ny); 
    Psi_hzx_xn = zeros(ncells,ny); 

    % Create and initialize cpml convolution coefficients 
    CPsi_eyx_xn = Ceyhz(2:ncells+1,:)*dx;
    CPsi_hzx_xn = Chzey(1:ncells,:)*dx;
    
    % Adjust FDTD coefficients in the CPML region 
    for i = 1: ncells                                           
        Ceyhz(i+1,:) = Ceyhz(i+1,:)/kappa_ex_xn(i);
        Chzey(i,:) = Chzey(i,:)/kappa_mx_xn(i);
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

    % Create and initialize cpml convolution parameters 
    Psi_eyx_xp = zeros(ncells,ny); 
    Psi_hzx_xp = zeros(ncells,ny); 

    % Create and initialize cpml convolution coefficients 
    CPsi_eyx_xp = Ceyhz(nxp1-ncells:nx,:)*dx;
    CPsi_hzx_xp = Chzey(nxp1-ncells:nx,:)*dx;
    
    % Adjust FDTD coefficients in the CPML region 
    for i = 1: ncells                                           
        Ceyhz(nx-ncells+i,:) = Ceyhz(nx-ncells+i,:)/kappa_ex_xp(i);
        Chzey(nx-ncells+i,:) = Chzey(nx-ncells+i,:)/kappa_mx_xp(i);
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

    % Create and initialize cpml convolution parameters 
    Psi_exy_yn = zeros(nx,ncells); 
    Psi_hzy_yn = zeros(nx,ncells); 

    % Create and initialize cpml convolution coefficients 
    CPsi_exy_yn = Cexhz(:,2:ncells+1)*dy;
    CPsi_hzy_yn = Chzex(:,1:ncells)*dy;
    
    % Adjust FDTD coefficients in the CPML region 
    for i = 1: ncells                                           
        Cexhz(:,i+1) = Cexhz(:,i+1)/kappa_ey_yn(i);
        Chzex(:,i) = Chzex(:,i)/kappa_my_yn(i);
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

    % Create and initialize cpml convolution parameters 
    Psi_exy_yp = zeros(nx,ncells); 
    Psi_hzy_yp = zeros(nx,ncells); 

    % Create and initialize cpml convolution coefficients 
    CPsi_exy_yp = Cexhz(:,nyp1-ncells:ny)*dy;
    CPsi_hzy_yp = Chzex(:,nyp1-ncells:ny)*dy;
    
    % Adjust FDTD coefficients in the CPML region 
    for i = 1: ncells                                           
        Cexhz(:,ny-ncells+i) = Cexhz(:,ny-ncells+i)/kappa_ey_yp(i);
        Chzex(:,ny-ncells+i) = Chzex(:,ny-ncells+i)/kappa_my_yp(i);
    end
    
    % Delete temporary arrays. These arrays will not be used any more.
    clear sigma_pey_yp sigma_pmy_yp;
    clear kappa_ey_yp kappa_my_yp;
    clear alpha_ey_yp alpha_my_yp;
end

