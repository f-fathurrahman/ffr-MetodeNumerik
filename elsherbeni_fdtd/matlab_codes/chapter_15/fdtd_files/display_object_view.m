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