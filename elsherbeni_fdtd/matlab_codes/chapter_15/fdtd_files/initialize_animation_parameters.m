if ~exist('animation','var')
    return;
end

% xy plane 
%==========================================================================
% x and y coordinates of vertices
x = meshgrid(0:nx,1:ny+1)*dx;
[x, y] = meshgrid(node_coordinates_xe, node_coordinates_ye);
x = reshape(x,1,[]).';
y = reshape(y,1,[]).';
vertices_xy_x = reshape(x.',1,[]).';
vertices_xy_y = reshape(y.',1,[]).';

% xy faces list
number_of_vertices = (nx+1)*(ny+1);
v_indices = 1:number_of_vertices;
v_indices = reshape(v_indices,ny+1,[]).';
faces_xy = [reshape(v_indices(1:nx,1:ny),1,[]).' reshape(v_indices(2:nx+1,1:ny),1,[]).' ...
    reshape(v_indices(2:nx+1,2:ny+1),1,[]).' reshape(v_indices(1:nx,2:ny+1),1,[]).'];

% yz plane 
%==========================================================================
% y and z coordinates of vertices
[y, z] = meshgrid(node_coordinates_ye, node_coordinates_ze);
y = reshape(y,1,[]).';
z = reshape(z,1,[]).';
vertices_yz_y = reshape(y.',1,[]).';
vertices_yz_z = reshape(z.',1,[]).';

% yz faces list
number_of_vertices = (ny+1)*(nz+1);
v_indices = 1:number_of_vertices;
v_indices = reshape(v_indices,nz+1,[]).';
faces_yz = [reshape(v_indices(1:ny,1:nz),1,[]).' reshape(v_indices(2:ny+1,1:nz),1,[]).' ...
    reshape(v_indices(2:ny+1,2:nz+1),1,[]).' reshape(v_indices(1:ny,2:nz+1),1,[]).'];

% zx plane 
%==========================================================================
% z and x coordinates of vertices
[z, x] = meshgrid(node_coordinates_ze, node_coordinates_xe);
z = reshape(z,1,[]).';
x = reshape(x,1,[]).';
vertices_zx_z = reshape(z.',1,[]).';
vertices_zx_x = reshape(x.',1,[]).';

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
            [V, I] = min(abs(animation(i).plane_cut(j).position - node_coordinates_ze));
            animation(i).plane_cut(j).position_index = I(1);    
            x = vertices_xy_x;
            y = vertices_xy_y;
            z = vertices_xy_y*0+node_coordinates_ze(animation(i).plane_cut(j).position_index);
            animation(i).plane_cut(j).vertices = [x y z];
            animation(i).plane_cut(j).faces = faces_xy;
            nc = nz;
        end
        if animation(i).plane_cut(j).type == 'yz';            
            [V, I] = min(abs(animation(i).plane_cut(j).position - node_coordinates_xe));
            animation(i).plane_cut(j).position_index = I(1);    
            y = vertices_yz_y;
            z = vertices_yz_z;
            x = vertices_yz_z*0+node_coordinates_xe(animation(i).plane_cut(j).position_index);
            animation(i).plane_cut(j).vertices = [x y z];
            animation(i).plane_cut(j).faces = faces_yz;
            nc = nx;
        end
        if animation(i).plane_cut(j).type == 'zx';
            [V, I] = min(abs(animation(i).plane_cut(j).position - node_coordinates_ye));
            animation(i).plane_cut(j).position_index = I(1);    
            z = vertices_zx_z;
            x = vertices_zx_x;
            y = vertices_zx_x*0+node_coordinates_ye(animation(i).plane_cut(j).position_index);
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
         if isfield(animation,'view_angles')
             view(animation(i).view_angles(1), animation(i).view_angles(2));
         else
            view(20,30);
         end
         if isfield(animation,'zoom')
             zoom(animation(i).zoom);
         end
%          axis([fdtd_domain.min_x fdtd_domain.max_x ...
%              fdtd_domain.min_y fdtd_domain.max_y ...
%              fdtd_domain.min_z fdtd_domain.max_z]);
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
        if animation(i).save_movie
            set(gca,'nextplot','replacechildren');          
            animation(i).frame_number = 0;
        end
    end
    if animation(i).display_grid 
        animation(i).edgecolor = 'k';
    else
        animation(i).edgecolor = 'none';
    end    
    if animation(i).display_objects
        display_objects_mesh_in_animation;
    end
end
