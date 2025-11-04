disp('filling material components arrays');

% creating temporary 1D arrays for storing 
% parameter values of material types
for ind = 1:size(material_types,2)
    t_eps_r(ind)   = material_types(ind).eps_r;
    t_mu_r(ind)    = material_types(ind).mu_r;
    t_sigma_e(ind) = material_types(ind).sigma_e;
    t_sigma_m(ind) = material_types(ind).sigma_m;
end

disp('Calculating eps_r_x');
% eps_r_x(i,j,k) is average of four cells 
% (i,j,k),(i,j-1,k), (i,j,k-1), (i,j-1,k-1) 
eps_r_x(1:nx,2:ny,2:nz) = 0.25 * (t_eps_r(material_3d_space(1:nx,2:ny,2:nz)) ...
                        + t_eps_r(material_3d_space(1:nx,1:ny-1,2:nz)) ...
                        + t_eps_r(material_3d_space(1:nx,2:ny,1:nz-1)) ...
                        + t_eps_r(material_3d_space(1:nx,1:ny-1,1:nz-1)));

disp('Calculating eps_r_y');
% eps_r_y(i,j,k) is average of four cells 
% (i,j,k),(i-1,j,k), (i,j,k-1), (i-1,j,k-1) 
eps_r_y(2:nx,1:ny,2:nz) = 0.25 * (t_eps_r(material_3d_space(2:nx,1:ny,2:nz)) ...
                        + t_eps_r(material_3d_space(1:nx-1,1:ny,2:nz)) ...
                        + t_eps_r(material_3d_space(2:nx,1:ny,1:nz-1)) ...
                        + t_eps_r(material_3d_space(1:nx-1,1:ny,1:nz-1)));
                    
disp('Calculating eps_r_z');
% eps_r_z(i,j,k) is average of four cells 
% (i,j,k),(i-1,j,k), (i,j-1,k), (i-1,j-1,k) 
eps_r_z(2:nx,2:ny,1:nz) = 0.25 * (t_eps_r(material_3d_space(2:nx,2:ny,1:nz)) ...
                        + t_eps_r(material_3d_space(1:nx-1,2:ny,1:nz)) ...
                        + t_eps_r(material_3d_space(2:nx,1:ny-1,1:nz)) ...
                        + t_eps_r(material_3d_space(1:nx-1,1:ny-1,1:nz)));

disp('Calculating sigma_e_x');
% sigma_e_x(i,j,k) is average of four cells 
% (i,j,k),(i,j-1,k), (i,j,k-1), (i,j-1,k-1) 
sigma_e_x(1:nx,2:ny,2:nz) = 0.25 * (t_sigma_e(material_3d_space(1:nx,2:ny,2:nz)) ...
                        + t_sigma_e(material_3d_space(1:nx,1:ny-1,2:nz)) ...
                        + t_sigma_e(material_3d_space(1:nx,2:ny,1:nz-1)) ...
                        + t_sigma_e(material_3d_space(1:nx,1:ny-1,1:nz-1)));

disp('Calculating sigma_e_y');
% sigma_e_y(i,j,k) is average of four cells
% (i,j,k),(i-1,j,k), (i,j,k-1), (i-1,j,k-1) 
sigma_e_y(2:nx,1:ny,2:nz) = 0.25 * (t_sigma_e(material_3d_space(2:nx,1:ny,2:nz)) ...
                        + t_sigma_e(material_3d_space(1:nx-1,1:ny,2:nz)) ...
                        + t_sigma_e(material_3d_space(2:nx,1:ny,1:nz-1)) ...
                        + t_sigma_e(material_3d_space(1:nx-1,1:ny,1:nz-1)));
                    
disp('Calculating sigma_e_z');
% sigma_e_z(i,j,k) is average of four cells
% (i,j,k),(i-1,j,k), (i,j-1,k), (i-1,j-1,k) 
sigma_e_z(2:nx,2:ny,1:nz) = 0.25 * (t_sigma_e(material_3d_space(2:nx,2:ny,1:nz)) ...
                        + t_sigma_e(material_3d_space(1:nx-1,2:ny,1:nz)) ...
                        + t_sigma_e(material_3d_space(2:nx,1:ny-1,1:nz)) ...
                        + t_sigma_e(material_3d_space(1:nx-1,1:ny-1,1:nz)));

                    
disp('Calculating mu_r_x');
% mu_r_x(i,j,k) is average of two cells (i,j,k),(i-1,j,k)
mu_r_x(2:nx,1:ny,1:nz) = 0.5 * (t_mu_r(material_3d_space(2:nx,1:ny,1:nz)) ...
                        + t_mu_r(material_3d_space(1:nx-1,1:ny,1:nz)));
                    
disp('Calculating mu_r_y');
% mu_r_y(i,j,k) is average of two cells (i,j,k),(i,j-1,k)
mu_r_y(1:nx,2:ny,1:nz) = 0.5 * (t_mu_r(material_3d_space(1:nx,2:ny,1:nz)) ...
                        + t_mu_r(material_3d_space(1:nx,1:ny-1,1:nz)));
                    
disp('Calculating mu_r_z');
% mu_r_z(i,j,k) is average of two cells (i,j,k),(i,j,k-1)
mu_r_z(1:nx,1:ny,2:nz) = 0.5 * (t_mu_r(material_3d_space(1:nx,1:ny,2:nz)) ...
                        + t_mu_r(material_3d_space(1:nx,1:ny,1:nz-1)));
                    
disp('Calculating sigma_m_x');
% sigma_m_x(i,j,k) is average of two cells (i,j,k),(i-1,j,k)
sigma_m_x(2:nx,1:ny,1:nz) = 0.5 * (t_sigma_m(material_3d_space(2:nx,1:ny,1:nz)) ...
                        + t_sigma_m(material_3d_space(1:nx-1,1:ny,1:nz)));
                    
disp('Calculating sigma_m_y');
% sigma_m_y(i,j,k) is average of two cells (i,j,k),(i,j-1,k)
sigma_m_y(1:nx,2:ny,1:nz) = 0.5 * (t_sigma_m(material_3d_space(1:nx,2:ny,1:nz)) ...
                        + t_sigma_m(material_3d_space(1:nx,1:ny-1,1:nz)));
                    
disp('Calculating sigma_m_z');
% sigma_m_z(i,j,k) is average of two cells (i,j,k),(i,j,k-1)
sigma_m_z(1:nx,1:ny,2:nz) = 0.5 * (t_sigma_m(material_3d_space(1:nx,1:ny,2:nz)) ...
                        + t_sigma_m(material_3d_space(1:nx,1:ny,1:nz-1)));

clear t_eps_r t_mu_r t_sigma_e t_sigma_m;
