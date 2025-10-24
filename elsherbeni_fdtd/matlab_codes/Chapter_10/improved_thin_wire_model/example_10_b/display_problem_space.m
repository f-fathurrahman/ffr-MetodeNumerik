% This subroutine displays the objects in the FDTD problem space

if ~exist('show_problem_space','var') 
    show_problem_space = false;
end
if ~show_problem_space
    return;
end

disp('drawing objects');

% default display settings
if ~exist('problem_space_display','var')
    problem_space_display.labels = true;
    problem_space_display.axis_at_origin = false;
    problem_space_display.axis_outside_domain = true;
    problem_space_display.grid_xn = false;
    problem_space_display.grid_xp = false;
    problem_space_display.grid_yn = false;
    problem_space_display.grid_yp = false;
    problem_space_display.grid_zn = false;
    problem_space_display.grid_zp = false;
    problem_space_display.outer_boundaries = true;
    problem_space_display.cpml_boundaries = true;
end

if ~isfield(problem_space_display,'labels')
    problem_space_display.labels = true;
end
if ~isfield(problem_space_display,'axis_at_origin')
    problem_space_display.axis_at_origin = false;
end
if ~isfield(problem_space_display,'axis_outside_domain')
    problem_space_display.axis_outside_domain = true;
end    
if ~isfield(problem_space_display,'grid_xn')
    problem_space_display.grid_xn = false;
end
if ~isfield(problem_space_display,'grid_xp')
    problem_space_display.grid_xp = false;
end
if ~isfield(problem_space_display,'grid_yn')
    problem_space_display.grid_yn = false;
end
if ~isfield(problem_space_display,'grid_yp')
    problem_space_display.grid_yp = false;
end
if ~isfield(problem_space_display,'grid_zn')
    problem_space_display.grid_zn = false;
end
if ~isfield(problem_space_display,'grid_zp')
    problem_space_display.grid_zp = false;
end    
if ~isfield(problem_space_display,'outer_boundaries')
    problem_space_display.outer_boundaries = true;
end
if ~isfield(problem_space_display,'cpml_boundaries')
    problem_space_display.cpml_boundaries = true;
end



%==========================================================================
% Display 3D view of problem geometry
figure_number = figure;

% display spheres
if exist('number_of_spheres','var')
    for i=1:number_of_spheres
        [X, Y, Z] = sphere(36);
        X = X * spheres(i).radius;
        Y = Y * spheres(i).radius;
        Z = Z * spheres(i).radius;
        X = X + spheres(i).center_x;
        Y = Y + spheres(i).center_y;
        Z = Z + spheres(i).center_z;
        RGB = material_types(spheres(i).material_type).color;
        patch(surf2patch(X,Y,Z,Z), 'FaceColor', RGB, ...
            'FaceAlpha', 0.5, 'facelighting', 'flat', ...
            'EdgeColor', 'none', 'facealpha', 0.5);
    end
end

% display bricks
if exist('number_of_bricks','var')
    for i=1:number_of_bricks
        vertices=[0 0 0; 1 0 0; 1 1 0; 0 1 0; 0 0 1; 1 0 1; 1 1 1; 0 1 1];
        faces=[1 2 6 5;2 3 7 6;3 4 8 7;4 1 5 8;1 2 3 4;5 6 7 8 ];
        vertices(:,1)=(vertices(:,1) * (bricks(i).max_x ...
            - bricks(i).min_x)) + bricks(i).min_x;
        vertices(:,2)=(vertices(:,2) * (bricks(i).max_y ...
            - bricks(i).min_y)) + bricks(i).min_y;
        vertices(:,3)=(vertices(:,3) * (bricks(i).max_z ...
            - bricks(i).min_z)) + bricks(i).min_z;
        RGB = material_types(bricks(i).material_type).color;
        patch('Vertices', vertices, 'Faces', faces,...
            'FaceColor', RGB, 'facealpha', 0.5);
    end
end

display_3D_geometry;

% enable the camera toolbar
cameratoolbar(figure_number);
axis equal;
axis off;
light;
camlight;
view(20,30);

%==========================================================================
% Display cell view of problem geometry
figure_number = figure;

for ind =1:size(material_types,2)
    matcol(ind,:) = material_types(ind).color;
end

m3d  = ones(nx+2, ny+2, nz+2);
tm3d = zeros(nx+2, ny+2, nz+2);
m3d(2:nx+1, 2:ny+1, 2:nz+1) = material_3d_space;
tm3d(2:nx+1, 2:ny+1, 2:nz+1) = m3d(2:nx+1, 2:ny+1, 2:nz+1) ...
    - m3d(1:nx, 2:ny+1, 2:nz+1);
I1 = find(tm3d>0);
tm3d(2:nx+1, 2:ny+1, 2:nz+1) = m3d(2:nx+1, 2:ny+1, 2:nz+1) ...
    - m3d(3:nx+2, 2:ny+1, 2:nz+1);
I2 = find(tm3d>2);
tm3d(2:nx+1, 2:ny+1, 2:nz+1) = m3d(2:nx+1, 2:ny+1, 2:nz+1) ...
    - m3d(2:nx+1, 1:ny, 2:nz+1);
I3 = find(tm3d>2);
tm3d(2:nx+1, 2:ny+1, 2:nz+1) = m3d(2:nx+1, 2:ny+1, 2:nz+1) ...
    - m3d(2:nx+1, 3:ny+2, 2:nz+1);
I4 = find(tm3d>2);
tm3d(2:nx+1, 2:ny+1, 2:nz+1) = m3d(2:nx+1, 2:ny+1, 2:nz+1) ...
    - m3d(2:nx+1, 2:ny+1, 1:nz);
I5 = find(tm3d>2);
tm3d(2:nx+1, 2:ny+1, 2:nz+1) = m3d(2:nx+1, 2:ny+1, 2:nz+1) ...
    - m3d(2:nx+1, 2:ny+1, 3:nz+2);
I6 = find(tm3d>2);
clear tm3d;

xx = zeros(nx+2,ny+2,nz+2);
yy = zeros(nx+2,ny+2,nz+2);
zz = zeros(nx+2,ny+2,nz+2);

for ind = -1:nx
    xx(ind+2,:,:) = ind * dx;
end
xx = xx + fdtd_domain.min_x;
for ind = -1:ny
    yy(:,ind+2,:) = ind * dy;
end
yy = yy + fdtd_domain.min_y;
for ind = -1:nz
    zz(:,:,ind+2) = ind * dz;
end
zz = zz + fdtd_domain.min_z;

len = size(I1,1);
vx = [xx(I1) xx(I1) xx(I1) xx(I1)];
vy = [yy(I1) yy(I1)+dy  yy(I1)+dy yy(I1)];
vz = [zz(I1) zz(I1) zz(I1)+dz zz(I1)+dz];
vx  = reshape(vx.',4*len,1);    
vy  = reshape(vy.',4*len,1);    
vz  = reshape(vz.',4*len,1);    
v1 = [vx vy vz];
f1 = [(1:4:4*(len-1)+1).' (2:4:4*(len-1)+2).' (3:4:4*(len-1)+3).' (4:4:4*(len-1)+4).'];
RGB1 = matcol(m3d(I1),:);

len = size(I2,1);
vx = [xx(I2) xx(I2) xx(I2) xx(I2)]+dx;
vy = [yy(I2) yy(I2)+dy  yy(I2)+dy yy(I2)];
vz = [zz(I2) zz(I2) zz(I2)+dz zz(I2)+dz];
vx  = reshape(vx.',4*len,1);    
vy  = reshape(vy.',4*len,1);    
vz  = reshape(vz.',4*len,1);    
v2 = [vx vy vz];
f2 = [(1:4:4*(len-1)+1).' (2:4:4*(len-1)+2).' (3:4:4*(len-1)+3).' (4:4:4*(len-1)+4).'];
RGB2 = matcol(m3d(I2),:);

len = size(I3,1);
vx = [xx(I3) xx(I3)+dx xx(I3)+dx xx(I3)];
vy = [yy(I3) yy(I3) yy(I3) yy(I3)];
vz = [zz(I3) zz(I3) zz(I3)+dz zz(I3)+dz];
vx  = reshape(vx.',4*len,1);    
vy  = reshape(vy.',4*len,1);    
vz  = reshape(vz.',4*len,1);    
v3 = [vx vy vz];
f3 = [(1:4:4*(len-1)+1).' (2:4:4*(len-1)+2).' (3:4:4*(len-1)+3).' (4:4:4*(len-1)+4).'];
RGB3 = matcol(m3d(I3),:);

len = size(I4,1);
vx = [xx(I4) xx(I4)+dx xx(I4)+dx xx(I4)];
vy = [yy(I4) yy(I4) yy(I4) yy(I4)] + dy;
vz = [zz(I4) zz(I4) zz(I4)+dz zz(I4)+dz];
vx  = reshape(vx.',4*len,1);    
vy  = reshape(vy.',4*len,1);    
vz  = reshape(vz.',4*len,1);    
v4 = [vx vy vz];
f4 = [(1:4:4*(len-1)+1).' (2:4:4*(len-1)+2).' (3:4:4*(len-1)+3).' (4:4:4*(len-1)+4).'];
RGB4 = matcol(m3d(I4),:);

len = size(I5,1);
vx = [xx(I5) xx(I5) xx(I5)+dx xx(I5)+dx];
vy = [yy(I5) yy(I5)+dy  yy(I5)+dy yy(I5)];
vz = [zz(I5) zz(I5) zz(I5) zz(I5)];
vx  = reshape(vx.',4*len,1);    
vy  = reshape(vy.',4*len,1);    
vz  = reshape(vz.',4*len,1);    
v5 = [vx vy vz];
f5 = [(1:4:4*(len-1)+1).' (2:4:4*(len-1)+2).' (3:4:4*(len-1)+3).' (4:4:4*(len-1)+4).'];
RGB5 = matcol(m3d(I5),:);

len = size(I6,1);
vx = [xx(I6) xx(I6) xx(I6)+dx xx(I6)+dx];
vy = [yy(I6) yy(I6)+dy  yy(I6)+dy yy(I6)];
vz = [zz(I6) zz(I6) zz(I6) zz(I6)] + dz;
vx  = reshape(vx.',4*len,1);    
vy  = reshape(vy.',4*len,1);    
vz  = reshape(vz.',4*len,1);    
v6 = [vx vy vz];
f6 = [(1:4:4*(len-1)+1).' (2:4:4*(len-1)+2).' (3:4:4*(len-1)+3).' (4:4:4*(len-1)+4).'];
RGB6 = matcol(m3d(I6),:);

vertices = [v1;v2;v3;v4;v5;v6];

faces = f1; 
faces = [faces; f2+ max(max(faces))];
faces = [faces; f3+ max(max(faces))];
faces = [faces; f4+ max(max(faces))];
faces = [faces; f5+ max(max(faces))];
faces = [faces; f6+ max(max(faces))];
RGB  = [RGB1;RGB2;RGB3;RGB4;RGB5;RGB6];

if ~isempty(vertices)
    patch('Vertices',vertices,'Faces',faces,'facevertexcdata',RGB,'FaceColor','flat');
end

% display bricks with zero thickness
if exist('number_of_bricks','var')
    for ind = 1:number_of_bricks

        mtype = bricks(ind).material_type;
        sigma_pec = material_types(mtype).sigma_e;

        % convert coordinates to node indices on the FDTD grid
        blx = round((bricks(ind).min_x - fdtd_domain.min_x)/dx)+1; 
        bly = round((bricks(ind).min_y - fdtd_domain.min_y)/dy)+1; 
        blz = round((bricks(ind).min_z - fdtd_domain.min_z)/dz)+1; 

        bux = round((bricks(ind).max_x - fdtd_domain.min_x)/dx)+1; 
        buy = round((bricks(ind).max_y - fdtd_domain.min_y)/dy)+1; 
        buz = round((bricks(ind).max_z - fdtd_domain.min_z)/dz)+1; 

        % find the zero thickness bricks
        if (blx == bux)
            len = (buy-bly) * (buz-blz);
            vx = xx(blx,bly:buy-1,blz:buz-1)+dx; 
            vy = yy(blx,bly:buy-1,blz:buz-1)+dy; 
            vz = zz(blx,bly:buy-1,blz:buz-1)+dz; 
            vx = reshape(vx,len,1);
            vy = reshape(vy,len,1);
            vz = reshape(vz,len,1);
            vx = [vx vx vx vx];
            vy = [vy vy+dy vy+dy vy];
            vz = [vz vz vz+dz vz+dz];
            vx  = reshape(vx.',4*len,1);    
            vy  = reshape(vy.',4*len,1);    
            vz  = reshape(vz.',4*len,1);    
            v = [vx vy vz];
            f = [(1:4:4*(len-1)+1).' (2:4:4*(len-1)+2).' (3:4:4*(len-1)+3).' (4:4:4*(len-1)+4).'];
            RGB = matcol(mtype,:);
            patch('Vertices',v,'Faces',f,'facevertexcdata',RGB,'FaceColor','flat');
        end
        if (bly == buy)
            len = (bux-blx) * (buz-blz);
            vx = xx(blx:bux-1,bly,blz:buz-1)+dx; 
            vy = yy(blx:bux-1,bly,blz:buz-1)+dy; 
            vz = zz(blx:bux-1,bly,blz:buz-1)+dz; 
            vx = reshape(vx,len,1);
            vy = reshape(vy,len,1);
            vz = reshape(vz,len,1);
            vx = [vx vx+dx vx+dx vx];
            vy = [vy vy vy vy];
            vz = [vz vz vz+dz vz+dz];
            vx  = reshape(vx.',4*len,1);    
            vy  = reshape(vy.',4*len,1);    
            vz  = reshape(vz.',4*len,1);    
            v = [vx vy vz];
            f = [(1:4:4*(len-1)+1).' (2:4:4*(len-1)+2).' (3:4:4*(len-1)+3).' (4:4:4*(len-1)+4).'];
            RGB = matcol(mtype,:);
            patch('Vertices',v,'Faces',f,'facevertexcdata',RGB,'FaceColor','flat');
        end
        if (blz == buz)
            len = (bux-blx) * (buy-bly);
            vx = xx(blx:bux-1,bly:buy-1,blz)+dx; 
            vy = yy(blx:bux-1,bly:buy-1,blz)+dy; 
            vz = zz(blx:bux-1,bly:buy-1,blz)+dz; 
            vx = reshape(vx,len,1);
            vy = reshape(vy,len,1);
            vz = reshape(vz,len,1);
            vx = [vx vx+dx vx+dx vx];
            vy = [vy vy vy+dy vy+dy];
            vz = [vz vz vz vz];
            vx  = reshape(vx.',4*len,1);    
            vy  = reshape(vy.',4*len,1);    
            vz  = reshape(vz.',4*len,1);    
            v = [vx vy vz];
            f = [(1:4:4*(len-1)+1).' (2:4:4*(len-1)+2).' (3:4:4*(len-1)+3).' (4:4:4*(len-1)+4).'];
            RGB = matcol(mtype,:);
            patch('Vertices',v,'Faces',f,'facevertexcdata',RGB,'FaceColor','flat');
        end
    end
end
 
display_3D_geometry;

% enable the camera toolbar
cameratoolbar(figure_number);
axis equal;
axis off;
view(20,30);

clear m3d xx yy zz vx vy vz v I1 I2 I3 I4 I5 I6;
clear v1 v2 v3 v4 v5 v6 f1 f2 f3 f4 f5 f6 RGB1 RGB2 RGB3 RGB4 RGB5 RGB6 RGB;
 