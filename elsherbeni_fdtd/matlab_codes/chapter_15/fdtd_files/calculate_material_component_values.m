disp('filling material components arrays');

% creating temporary 1D arrays for storing 
% parameter values of material types
for ind = 1:size(material_types,2)
    t_eps_r(ind)   = material_types(ind).eps_r;
    t_mu_r(ind)    = material_types(ind).mu_r;
    t_sigma_e(ind) = material_types(ind).sigma_e;
    t_sigma_m(ind) = material_types(ind).sigma_m;
end

% assign negligibly small values to t_mu_r and t_sigma_m where they are
% zero in order to prevent division by zero error 
t_mu_r(find(t_mu_r==0)) = 1e-20;
t_sigma_m(find(t_sigma_m==0)) = 1e-20;

% create arrays to store cell crosssection areas
[X, Y, Z] = ndgrid(cell_sizes_xe, cell_sizes_ye, cell_sizes_ze); 
ca_xy = X.*Y;
ca_yz = Y.*Z;
ca_zx = Z.*X;

disp('Calculating eps_r_x');
% eps_r_x(i,j,k) is cell area weighted average of four cells 

area_1 = ca_yz(1:nx,2:ny,2:nz);
area_2 = ca_yz(1:nx,1:ny-1,2:nz);
area_3 = ca_yz(1:nx,2:ny,1:nz-1);
area_4 = ca_yz(1:nx,1:ny-1,1:nz-1);
area_total = area_1 + area_2 + area_3 + area_4;

eps_1 = t_eps_r(material_3d_space(1:nx,2:ny,2:nz));
eps_2 = t_eps_r(material_3d_space(1:nx,1:ny-1,2:nz));
eps_3 = t_eps_r(material_3d_space(1:nx,2:ny,1:nz-1));
eps_4 = t_eps_r(material_3d_space(1:nx,1:ny-1,1:nz-1));

eps_r_x(1:nx,2:ny,2:nz) = (1./area_total) ...
    .* (area_1.* eps_1 + area_2.* eps_2 + area_3.* eps_3 + area_4.* eps_4);

disp('Calculating eps_r_y');
% eps_r_y(i,j,k) is cell area weighted average of four cells  

area_1 = ca_zx(2:nx,1:ny,2:nz);
area_2 = ca_zx(1:nx-1,1:ny,2:nz);
area_3 = ca_zx(2:nx,1:ny,1:nz-1);
area_4 = ca_zx(1:nx-1,1:ny,1:nz-1);
area_total = area_1 + area_2 + area_3 + area_4;

eps_1 = t_eps_r(material_3d_space(2:nx,1:ny,2:nz));
eps_2 = t_eps_r(material_3d_space(1:nx-1,1:ny,2:nz));
eps_3 = t_eps_r(material_3d_space(2:nx,1:ny,1:nz-1));
eps_4 = t_eps_r(material_3d_space(1:nx-1,1:ny,1:nz-1));

eps_r_y(2:nx,1:ny,2:nz) = (1./area_total) ...
    .* (area_1.* eps_1 + area_2.* eps_2 + area_3.* eps_3 + area_4.* eps_4);
                    
disp('Calculating eps_r_z');
% eps_r_z(i,j,k) is cell area weighted average of four cells 

area_1 = ca_xy(2:nx,2:ny,1:nz);
area_2 = ca_xy(1:nx-1,2:ny,1:nz);
area_3 = ca_xy(2:nx,1:ny-1,1:nz);
area_4 = ca_xy(1:nx-1,1:ny-1,1:nz);
area_total = area_1 + area_2 + area_3 + area_4;

eps_1 = t_eps_r(material_3d_space(2:nx,2:ny,1:nz));
eps_2 = t_eps_r(material_3d_space(1:nx-1,2:ny,1:nz));
eps_3 = t_eps_r(material_3d_space(2:nx,1:ny-1,1:nz));
eps_4 = t_eps_r(material_3d_space(1:nx-1,1:ny-1,1:nz));

eps_r_z(2:nx,2:ny,1:nz) = (1./area_total) ...
    .* (area_1.* eps_1 + area_2.* eps_2 + area_3.* eps_3 + area_4.* eps_4);

disp('Calculating sigma_e_x');
% sigma_e_x(i,j,k) is cell area weighted average of four cells 

area_1 = ca_yz(1:nx,2:ny,2:nz);
area_2 = ca_yz(1:nx,1:ny-1,2:nz);
area_3 = ca_yz(1:nx,2:ny,1:nz-1);
area_4 = ca_yz(1:nx,1:ny-1,1:nz-1);
area_total = area_1 + area_2 + area_3 + area_4;

sigma_1 = t_sigma_e(material_3d_space(1:nx,2:ny,2:nz));
sigma_2 = t_sigma_e(material_3d_space(1:nx,1:ny-1,2:nz));
sigma_3 = t_sigma_e(material_3d_space(1:nx,2:ny,1:nz-1));
sigma_4 = t_sigma_e(material_3d_space(1:nx,1:ny-1,1:nz-1));

sigma_e_x(1:nx,2:ny,2:nz) = (1./area_total) ...
    .* (area_1.* sigma_1 + area_2.* sigma_2 + area_3.* sigma_3 + area_4.* sigma_4);

disp('Calculating sigma_e_y');
% sigma_e_y(i,j,k) is cell area weighted average of four cells  

area_1 = ca_zx(2:nx,1:ny,2:nz);
area_2 = ca_zx(1:nx-1,1:ny,2:nz);
area_3 = ca_zx(2:nx,1:ny,1:nz-1);
area_4 = ca_zx(1:nx-1,1:ny,1:nz-1);
area_total = area_1 + area_2 + area_3 + area_4;

sigma_1 = t_sigma_e(material_3d_space(2:nx,1:ny,2:nz));
sigma_2 = t_sigma_e(material_3d_space(1:nx-1,1:ny,2:nz));
sigma_3 = t_sigma_e(material_3d_space(2:nx,1:ny,1:nz-1));
sigma_4 = t_sigma_e(material_3d_space(1:nx-1,1:ny,1:nz-1));

sigma_e_y(2:nx,1:ny,2:nz) = (1./area_total) ...
    .* (area_1.* sigma_1 + area_2.* sigma_2 + area_3.* sigma_3 + area_4.* sigma_4);
                    
disp('Calculating sigma_e_z');
% sigma_e_z(i,j,k) is cell area weighted average of four cells 

area_1 = ca_xy(2:nx,2:ny,1:nz);
area_2 = ca_xy(1:nx-1,2:ny,1:nz);
area_3 = ca_xy(2:nx,1:ny-1,1:nz);
area_4 = ca_xy(1:nx-1,1:ny-1,1:nz);
area_total = area_1 + area_2 + area_3 + area_4;

sigma_1 = t_sigma_e(material_3d_space(2:nx,2:ny,1:nz));
sigma_2 = t_sigma_e(material_3d_space(1:nx-1,2:ny,1:nz));
sigma_3 = t_sigma_e(material_3d_space(2:nx,1:ny-1,1:nz));
sigma_4 = t_sigma_e(material_3d_space(1:nx-1,1:ny-1,1:nz));

sigma_e_z(2:nx,2:ny,1:nz) = (1./area_total) ...
    .* (area_1.* sigma_1 + area_2.* sigma_2 + area_3.* sigma_3 + area_4.* sigma_4);

% clear temporary arrays
clear X Y Z ca_xy ca_yz ca_zx;
clear area_1 area_2 area_3 area_4;
clear sigma_1 sigma_2 sigma_3 sigma_4;
clear eps_1 eps_2 eps_3 eps_4;

% create arrays to store cell lengths
[cl_x, cl_y, cl_z] = ndgrid(cell_sizes_xe, cell_sizes_ye, cell_sizes_ze); 

disp('Calculating mu_r_x');
% mu_r_x(i,j,k) is average of two cells (i,j,k),(i-1,j,k)
length_1 = cl_x(2:nx,1:ny,1:nz);
length_2 = cl_x(1:nx-1,1:ny,1:nz);
mu_1 = t_mu_r(material_3d_space(2:nx,1:ny,1:nz));
mu_2 = t_mu_r(material_3d_space(1:nx-1,1:ny,1:nz));

mu_r_x(2:nx,1:ny,1:nz) = ...
    ((length_1 + length_2) .* mu_1 .* mu_2) ...
    ./(length_2.*mu_1 + length_1.*mu_2);
                    
disp('Calculating mu_r_y');
% mu_r_y(i,j,k) is average of two cells (i,j,k),(i,j-1,k)
length_1 = cl_y(1:nx,2:ny,1:nz);
length_2 = cl_y(1:nx,1:ny-1,1:nz);
mu_1 = t_mu_r(material_3d_space(1:nx,2:ny,1:nz));
mu_2 = t_mu_r(material_3d_space(1:nx,1:ny-1,1:nz));

mu_r_y(1:nx,2:ny,1:nz) = ...
    ((length_1 + length_2) .* mu_1 .* mu_2) ...
    ./(length_2.*mu_1 + length_1.*mu_2);
                    
disp('Calculating mu_r_z');
% mu_r_z(i,j,k) is average of two cells (i,j,k),(i,j,k-1)

length_1 = cl_z(1:nx,1:ny,2:nz);
length_2 = cl_z(1:nx,1:ny,1:nz-1);
mu_1 = t_mu_r(material_3d_space(1:nx,1:ny,2:nz));
mu_2 = t_mu_r(material_3d_space(1:nx,1:ny,1:nz-1));

mu_r_z(1:nx,1:ny,2:nz) = ...
    ((length_1 + length_2) .* mu_1 .* mu_2) ...
    ./(length_2.*mu_1 + length_1.*mu_2);
                    
disp('Calculating sigma_m_x');
% sigma_m_x(i,j,k) is average of two cells (i,j,k),(i-1,j,k)
length_1 = cl_x(2:nx,1:ny,1:nz);
length_2 = cl_x(1:nx-1,1:ny,1:nz);
sigma_1 = t_sigma_m(material_3d_space(2:nx,1:ny,1:nz));
sigma_2 = t_sigma_m(material_3d_space(1:nx-1,1:ny,1:nz));

sigma_m_x(2:nx,1:ny,1:nz) = ...
    ((length_1 + length_2) .* sigma_1 .* sigma_2) ...
    ./(length_2.*sigma_1 + length_1.*sigma_2);
                    
disp('Calculating sigma_m_y');
% sigma_m_y(i,j,k) is average of two cells (i,j,k),(i,j-1,k)
length_1 = cl_y(1:nx,2:ny,1:nz);
length_2 = cl_y(1:nx,1:ny-1,1:nz);
sigma_1 = t_sigma_m(material_3d_space(1:nx,2:ny,1:nz));
sigma_2 = t_sigma_m(material_3d_space(1:nx,1:ny-1,1:nz));

sigma_m_y(1:nx,2:ny,1:nz) = ...
    ((length_1 + length_2) .* sigma_1 .* sigma_2) ...
    ./(length_2.*sigma_1 + length_1.*sigma_2);
                    
disp('Calculating sigma_m_z');
% sigma_m_z(i,j,k) is average of two cells (i,j,k),(i,j,k-1)

length_1 = cl_z(1:nx,1:ny,2:nz);
length_2 = cl_z(1:nx,1:ny,1:nz-1);
sigma_1 = t_sigma_m(material_3d_space(1:nx,1:ny,2:nz));
sigma_2 = t_sigma_m(material_3d_space(1:nx,1:ny,1:nz-1));

sigma_m_z(1:nx,1:ny,2:nz) = ...
    ((length_1 + length_2) .* sigma_1 .* sigma_2) ...
    ./(length_2.*sigma_1 + length_1.*sigma_2);

% clear temporary arrays
clear length_1 length_2 mu_1 mu_2 sigma_1 sigma_2
clear t_eps_r t_mu_r t_sigma_e t_sigma_m;
