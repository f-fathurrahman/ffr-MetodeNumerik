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

display_object_view;

display_cell_view;

 