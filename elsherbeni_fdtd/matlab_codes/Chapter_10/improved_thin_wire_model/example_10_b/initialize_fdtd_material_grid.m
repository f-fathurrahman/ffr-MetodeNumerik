disp('initializing FDTD material grid');

% calculate problem space size based on the object 
% locations and boundary conditions
calculate_domain_size;

% Array to store material type indices for every cell
% in the problem space. By default the space is filled 
% with air by initializing the array by ones
material_3d_space = ones(nx, ny, nz); 

% Create the 3D objects in the problem space by
% assigning indices of material types in the cells 
% to material_3d_space  

% create spheres
create_spheres;

% create bricks
create_bricks;

% Material component arrays for a problem space
% composed of (nx, ny, nz) cells
eps_r_x     = ones (nx  , nyp1 , nzp1);
eps_r_y     = ones (nxp1, ny   , nzp1);
eps_r_z     = ones (nxp1, nyp1 , nz);
mu_r_x      = ones (nxp1, ny   , nz);
mu_r_y      = ones (nx  , nyp1 , nz);
mu_r_z      = ones (nx  , ny   , nzp1);
sigma_e_x   = zeros(nx  , nyp1 , nzp1);
sigma_e_y   = zeros(nxp1, ny   , nzp1);
sigma_e_z   = zeros(nxp1, nyp1 , nz);
sigma_m_x   = zeros(nxp1, ny   , nz);
sigma_m_y   = zeros(nx  , nyp1 , nz);
sigma_m_z   = zeros(nx  , ny   , nzp1);

% calculate material component values by averaging
calculate_material_component_values;

% create zero thickness PEC plates
create_PEC_plates;

