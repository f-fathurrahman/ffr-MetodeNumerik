disp('creating spheres');

cx = fdtd_domain.cell_center_coordinates_x;
cy = fdtd_domain.cell_center_coordinates_y;

for ind=1:number_of_circles
% distance of the centers of the cells from the center of the circle
    distance = sqrt((circles(ind).center_x - cx).^2 ...
            + (circles(ind).center_y - cy).^2);
    I = find(distance<=circles(ind).radius);
    material_2d_space(I) = circles(ind).material_type;
end
clear cx cy;
