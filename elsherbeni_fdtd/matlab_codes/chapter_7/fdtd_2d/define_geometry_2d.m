disp('defining the problem geometry');

rectangles  = [];
circles = [];
% 
% define a rectangle
rectangles(1).min_x = 0;
rectangles(1).min_y = 1e-3;
rectangles(1).max_x = 1e-3;
rectangles(1).max_y = 2e-3;
rectangles(1).material_type = 4; 

% define a rectangle
rectangles(2).min_x = -1e-3;
rectangles(2).min_y = -2e-3;
rectangles(2).max_x = 0;
rectangles(2).max_y = -1e-3;
rectangles(2).material_type = 2; 

% define a circle
circles(1).center_x = -2e-3;
circles(1).center_y = 1e-3;
circles(1).radius = 1e-3;
circles(1).material_type = 4; 

% define a circle
circles(2).center_x = 2e-3;
circles(2).center_y = -1e-3;
circles(2).radius = 1e-3;
circles(2).material_type = 5; 
