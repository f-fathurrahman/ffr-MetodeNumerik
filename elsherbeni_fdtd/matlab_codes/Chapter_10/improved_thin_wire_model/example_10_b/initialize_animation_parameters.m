if ~exist('animation','var')
    return;
end

% xy plane 
%==========================================================================
% x coordinate of vertices
x = meshgrid(0:nx,1:ny+1)*dx;
x = reshape(x,1,[]).';
vertices_xy_x = reshape(x.',1,[]).';
vertices_xy_x = vertices_xy_x + fdtd_domain.min_x;

% y coordinate of vertices
y = (meshgrid(0:ny,1:nx+1)*dy).';
y = reshape(y,1,[]).';
vertices_xy_y = reshape(y.',1,[]).';
vertices_xy_y = vertices_xy_y + fdtd_domain.min_y;

% xy faces list
number_of_vertices = (nx+1)*(ny+1);
v_indices = 1:number_of_vertices;
v_indices = reshape(v_indices,ny+1,[]).';
faces_xy = [reshape(v_indices(1:nx,1:ny),1,[]).' reshape(v_indices(2:nx+1,1:ny),1,[]).' ...
    reshape(v_indices(2:nx+1,2:ny+1),1,[]).' reshape(v_indices(1:nx,2:ny+1),1,[]).'];

% yz plane 
%==========================================================================
% y coordinate of vertices
y = meshgrid(0:ny,1:nz+1)*dy;
y = reshape(y,1,[]).';
vertices_yz_y = reshape(y.',1,[]).';
vertices_yz_y = vertices_yz_y + fdtd_domain.min_y;

% z coordinate of vertices
z = (meshgrid(0:nz,1:ny+1)*dz).';
z = reshape(z,1,[]).';
vertices_yz_z = reshape(z.',1,[]).';
vertices_yz_z = vertices_yz_z + fdtd_domain.min_z;

% yz faces list
number_of_vertices = (ny+1)*(nz+1);
v_indices = 1:number_of_vertices;
v_indices = reshape(v_indices,nz+1,[]).';
faces_yz = [reshape(v_indices(1:ny,1:nz),1,[]).' reshape(v_indices(2:ny+1,1:nz),1,[]).' ...
    reshape(v_indices(2:ny+1,2:nz+1),1,[]).' reshape(v_indices(1:ny,2:nz+1),1,[]).'];

% zx plane 
%==========================================================================
% z coordinate of vertices
z = meshgrid(0:nz,1:nx+1)*dz;
z = reshape(z,1,[]).';
vertices_zx_z = reshape(z.',1,[]).';
vertices_zx_z = vertices_zx_z + fdtd_domain.min_z;

% x coordinate of vertices
x = (meshgrid(0:nx,1:nz+1)*dx).';
x = reshape(x,1,[]).';
vertices_zx_x = reshape(x.',1,[]).';
vertices_zx_x = vertices_zx_x + fdtd_domain.min_x;

% zx faces list
number_of_vertices = (nz+1)*(nx+1);
v_indices = 1:number_of_vertices;
v_indices = reshape(v_indices,nx+1,[]).';
faces_zx = [reshape(v_indices(1:nz,1:nx),1,[]).' reshape(v_indices(2:nz+1,1:nx),1,[]).' ...
    reshape(v_indices(2:nz+1,2:nx+1),1,[]).' reshape(v_indices(1:nz,2:nx+1),1,[]).'];

%===================================
for i=size(animation,2):-1:1
    if animation(i).enable == false
        disp(['animation ' num2str(i) ' is disabled and will be deleted.']);
        animation(i) = [];
        continue;
    end
    for j=size(animation(i).plane_cut,2):-1:1
        if animation(i).plane_cut(j).type == 'xy';
            animation(i).plane_cut(j).position_index = round((animation(i).plane_cut(j).position-fdtd_domain.min_z)/dz)+1;
            x = vertices_xy_x;
            y = vertices_xy_y;
            z = vertices_xy_y*0+(animation(i).plane_cut(j).position_index-1) * dz + fdtd_domain.min_z;
            animation(i).plane_cut(j).vertices = [x y z];
            animation(i).plane_cut(j).faces = faces_xy;
            nc = nz;
        end
        if animation(i).plane_cut(j).type == 'yz';
            animation(i).plane_cut(j).position_index = round((animation(i).plane_cut(j).position-fdtd_domain.min_x)/dx)+1;
            y = vertices_yz_y;
            z = vertices_yz_z;
            x = vertices_yz_z*0+(animation(i).plane_cut(j).position_index-1) * dx + fdtd_domain.min_x;
            animation(i).plane_cut(j).vertices = [x y z];
            animation(i).plane_cut(j).faces = faces_yz;
            nc = nx;
        end
        if animation(i).plane_cut(j).type == 'zx';
            animation(i).plane_cut(j).position_index = round((animation(i).plane_cut(j).position-fdtd_domain.min_y)/dy)+1;
            z = vertices_zx_z;
            x = vertices_zx_x;
            y = vertices_zx_x*0+(animation(i).plane_cut(j).position_index-1) * dy + fdtd_domain.min_y;
            animation(i).plane_cut(j).vertices = [x y z];
            animation(i).plane_cut(j).faces = faces_zx;
            nc = ny;
        end
        if (animation(i).plane_cut(j).position_index > nc) || (animation(i).plane_cut(j).position_index < 2)
            disp(['animation ' num2str(i) '- plane cut ' num2str(j) ' is invalid and will be deleted.']);
            animation(i).plane_cut(j) = [];
        end
    end
end
%=============================
% rearrange vertices and faces
for i=1:size(animation,2)
    vertices = [];
    faces = [];
    for j=1:size(animation(i).plane_cut,2)
        faces = [faces; size(vertices,1)+animation(i).plane_cut(j).faces];
        vertices = [vertices; animation(i).plane_cut(j).vertices];
        animation(i).plane_cut(j).vertices = [];
        animation(i).plane_cut(j).faces = [];
    end
    animation(i).vertices = vertices;
    animation(i).faces = faces;
end

%=================================================
for i=1:size(animation,2)
    if animation(i).enable
         animation(i).figure_number = figure;
         cameratoolbar(animation(i).figure_number);
         axis vis3d;
         axis equal;
         view(20,30);
         axis([fdtd_domain.min_x fdtd_domain.max_x ...
             fdtd_domain.min_y fdtd_domain.max_y ...
             fdtd_domain.min_z fdtd_domain.max_z]);
         axis off;
         if animation(i).field_type == 'e';
             str = 'Electric field';
         else
             str = 'Magnetic field';
         end
        switch animation(i).component
            case 'x'
                str = [str ', x-component'];
            case 'y'
                str = [str ', y-component'];
            case 'z'
                str = [str ', z-component'];
            case 'm'
                str = [str ', magnitude'];
        end
          set(animation(i).figure_number,'name',str);
    end
    if animation(i).display_grid 
        animation(i).edgecolor = 'k';
    else
        animation(i).edgecolor = 'none';
    end    
end
