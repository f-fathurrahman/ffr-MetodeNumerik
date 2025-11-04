% adjust cell sizes and number of cells in the computation space

number_of_subregions = 0;
if exist('subregions','var')
    number_of_subregions = size(subregions,2);
end

% number of cells in x, y, and z directions
nx = round(fdtd_domain.size_x/dx);  
ny = round(fdtd_domain.size_y/dy);
nz = round(fdtd_domain.size_z/dz);
node_coordinates_xe = (0:nx)*dx + fdtd_domain.min_x;
node_coordinates_ye = (0:ny)*dy + fdtd_domain.min_y;
node_coordinates_ze = (0:nz)*dz + fdtd_domain.min_z;

node_coordinates_xh = (0:nx+1)*dx + fdtd_domain.min_x - dx/2;
node_coordinates_yh = (0:ny+1)*dy + fdtd_domain.min_y - dy/2;
node_coordinates_zh = (0:nz+1)*dz + fdtd_domain.min_z - dz/2;
    
for ind = 1:number_of_subregions
    switch (subregions(ind).direction)
        case 'x'
        [node_coordinates_xe, node_coordinates_xh] ...
            = insert_subregion_into_domain ...
            (node_coordinates_xe, node_coordinates_xh, subregions(ind), dx);
        case 'y'
        [node_coordinates_ye, node_coordinates_yh] ...
            = insert_subregion_into_domain ...
            (node_coordinates_ye, node_coordinates_yh, subregions(ind), dy);
        case 'z'
        [node_coordinates_ze, node_coordinates_zh] ...
            = insert_subregion_into_domain ...
            (node_coordinates_ze, node_coordinates_zh, subregions(ind), dz);
    end    
end

if isfield(boundary,'is_nonuniform_cpml')
    if (boundary.is_nonuniform_cpml)    
        set_nonuniform_cpml_grid;   
    end
end

% reset number of cells
nx = size(node_coordinates_xe, 2)-1;
ny = size(node_coordinates_ye, 2)-1;
nz = size(node_coordinates_ze, 2)-1;

% create cell size arrays
cell_sizes_xe = node_coordinates_xe(2:nx+1) - node_coordinates_xe(1:nx);
cell_sizes_ye = node_coordinates_ye(2:ny+1) - node_coordinates_ye(1:ny);
cell_sizes_ze = node_coordinates_ze(2:nz+1) - node_coordinates_ze(1:nz);

% create cell size arrays conforming the magnetic field grid
cell_sizes_xh = node_coordinates_xh(2:nx+2) - node_coordinates_xh(1:nx+1);
cell_sizes_yh = node_coordinates_yh(2:ny+2) - node_coordinates_yh(1:ny+1);
cell_sizes_zh = node_coordinates_zh(2:nz+2) - node_coordinates_zh(1:nz+1);

fdtd_domain.node_coordinates_xe = node_coordinates_xe;
fdtd_domain.node_coordinates_ye = node_coordinates_ye;
fdtd_domain.node_coordinates_ze = node_coordinates_ze;

fdtd_domain.node_coordinates_xh = node_coordinates_xh;
fdtd_domain.node_coordinates_yh = node_coordinates_yh;
fdtd_domain.node_coordinates_zh = node_coordinates_zh;

% some frequently used auxiliary parameters 
nxp1 = nx+1;    nyp1 = ny+1;    nzp1 = nz+1;
nxm1 = nx-1;    nxm2 = nx-2;	nym1 = ny-1;
nym2 = ny-2;	nzm1 = nz-1;	nzm2 = nz-2;

% create arrays to store cell center coordinates
cell_center_coordinates_x = 0.5 * (node_coordinates_xe(2:nx+1) ...
    + node_coordinates_xe(1:nx));
cell_center_coordinates_y = 0.5 * (node_coordinates_ye(2:ny+1) ...
    + node_coordinates_ye(1:ny));
cell_center_coordinates_z = 0.5 * (node_coordinates_ze(2:nz+1) ...
    + node_coordinates_ze(1:nz));

fdtd_domain.min_x = fdtd_domain.node_coordinates_xe(1);
fdtd_domain.min_y = fdtd_domain.node_coordinates_ye(1);
fdtd_domain.min_z = fdtd_domain.node_coordinates_ze(1);

fdtd_domain.max_x = fdtd_domain.node_coordinates_xe(nxp1);
fdtd_domain.max_y = fdtd_domain.node_coordinates_ye(nyp1);
fdtd_domain.max_z = fdtd_domain.node_coordinates_ze(nzp1);

% Determining the problem space size
fdtd_domain.size_x = fdtd_domain.max_x - fdtd_domain.min_x;
fdtd_domain.size_y = fdtd_domain.max_y - fdtd_domain.min_y;
fdtd_domain.size_z = fdtd_domain.max_z - fdtd_domain.min_z;


