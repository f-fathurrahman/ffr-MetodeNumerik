% Create dispersive objects

number_of_material_types = size(material_types,2);
for ind = 1:number_of_material_types
    if ~isfield(material_types(ind),'lorentz_dispersive')
        material_types(ind).lorentz_dispersive = false;
    end            
end

% Three dimensional arrays to store field coordinates
% Used to identify field components that are inside a sphere
Ex_x_coordinates = zeros(nx, nyp1, nzp1);
Ex_y_coordinates = zeros(nx, nyp1, nzp1);
Ex_z_coordinates = zeros(nx, nyp1, nzp1);

Ey_x_coordinates = zeros(nxp1, ny, nzp1);
Ey_y_coordinates = zeros(nxp1, ny, nzp1);
Ey_z_coordinates = zeros(nxp1, ny, nzp1);

Ez_x_coordinates = zeros(nxp1, nyp1, nz);
Ez_y_coordinates = zeros(nxp1, nyp1, nz);
Ez_z_coordinates = zeros(nxp1, nyp1, nz);

% Arrays to store material indices distribution of dispersive objects
dispersive_Ex = zeros(nx, nyp1, nzp1);
dispersive_Ey = zeros(nxp1, ny, nzp1);
dispersive_Ez = zeros(nxp1, nyp1, nz);

% calculate field coordinates
for ind = 1:nx
    Ex_x_coordinates(ind,:,:) = (ind - 0.5) * dx + fdtd_domain.min_x;
end
for ind = 1:nyp1
    Ex_y_coordinates(:,ind,:) = (ind - 1) * dy + fdtd_domain.min_y;
end
for ind = 1:nzp1
    Ex_z_coordinates(:,:,ind) = (ind - 1) * dz + fdtd_domain.min_z;
end

for ind = 1:nxp1
    Ey_x_coordinates(ind,:,:) = (ind - 1) * dx + fdtd_domain.min_x;
end
for ind = 1:ny
    Ey_y_coordinates(:,ind,:) = (ind - 0.5) * dy + fdtd_domain.min_y;
end
for ind = 1:nzp1
    Ey_z_coordinates(:,:,ind) = (ind - 1) * dz + fdtd_domain.min_z;
end

for ind = 1:nxp1
    Ez_x_coordinates(ind,:,:) = (ind - 1) * dx + fdtd_domain.min_x;
end
for ind = 1:nyp1
    Ez_y_coordinates(:,ind,:) = (ind - 1) * dy + fdtd_domain.min_y;
end
for ind = 1:nz
    Ez_z_coordinates(:,:,ind) = (ind - 0.5) * dz + fdtd_domain.min_z;
end

disp('creating dispersive spheres');
for ind=1:number_of_spheres
    mat_ind = spheres(ind).material_type;
    material = material_types(mat_ind);
    if material.lorentz_dispersive 
        % distance of the Ex components from the center of the sphere
        distance = sqrt((spheres(ind).center_x - Ex_x_coordinates).^2 ...
                + (spheres(ind).center_y - Ex_y_coordinates).^2 ...
                + (spheres(ind).center_z - Ex_z_coordinates).^2);
        I = find(distance<=spheres(ind).radius);
        dispersive_Ex(I) = mat_ind;

        % distance of the Ey components from the center of the sphere
        distance = sqrt((spheres(ind).center_x - Ey_x_coordinates).^2 ...
                + (spheres(ind).center_y - Ey_y_coordinates).^2 ...
                + (spheres(ind).center_z - Ey_z_coordinates).^2);
        I = find(distance<=spheres(ind).radius);
        dispersive_Ey(I) = mat_ind;

        % distance of the Ez components from the center of the sphere
        distance = sqrt((spheres(ind).center_x - Ez_x_coordinates).^2 ...
                + (spheres(ind).center_y - Ez_y_coordinates).^2 ...
                + (spheres(ind).center_z - Ez_z_coordinates).^2);
        I = find(distance<=spheres(ind).radius);
        dispersive_Ez(I) = mat_ind;
    end
end

disp('creating dispersive bricks');
for ind = 1:number_of_bricks
    mat_ind = bricks(ind).material_type;
    material = material_types(mat_ind);
    if material.lorentz_dispersive 
        % convert brick end coordinates to node indices 
        ni = get_node_indices(bricks(ind), fdtd_domain);
        is = ni.is; js = ni.js; ks = ni.ks; 
        ie = ni.ie; je = ni.je; ke = ni.ke; 

        dispersive_Ex (is:ie-1, js:je, ks:ke) = mat_ind; 
        dispersive_Ey (is:ie, js:je-1, ks:ke) = mat_ind;  
        dispersive_Ez (is:ie, js:je, ks:ke-1) = mat_ind;  
    end
end

% create a new structure array to store indices of field components 
% that need to be updated using Lorentz formulation
number_of_lorentz = 0;
for ind = 1:number_of_material_types
    material = material_types(ind);
    if material.lorentz_dispersive
        number_of_lorentz = number_of_lorentz + 1;

        lorentz(number_of_lorentz).material = ind;
        
        I = find(dispersive_Ex == ind);
        lorentz(number_of_lorentz).Ex_indices = I;
        eps_r_x(I) = material.eps_r;
        sigma_e_x(I) = material.sigma_e;   
        
        I = find(dispersive_Ey == ind);
        lorentz(number_of_lorentz).Ey_indices = I;
        eps_r_y(I) = material.eps_r;
        sigma_e_y(I) = material.sigma_e;   
    
        I = find(dispersive_Ez == ind);
        lorentz(number_of_lorentz).Ez_indices = I;
        eps_r_z(I) = material.eps_r;
        sigma_e_z(I) = material.sigma_e;
    end
end

clear Ex_x_coordinates Ex_y_coordinates Ex_z_coordinates 
clear Ey_x_coordinates Ey_y_coordinates Ey_z_coordinates 
clear Ez_x_coordinates Ez_y_coordinates Ez_z_coordinates 
clear dispersive_Ex dispersive_Ey dispersive_Ez 
