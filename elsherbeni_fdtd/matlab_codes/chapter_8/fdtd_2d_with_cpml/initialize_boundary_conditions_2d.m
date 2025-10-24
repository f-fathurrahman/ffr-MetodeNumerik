% initialize boundary parameters

disp('initializing boundary conditions');

% define logical parameters for the conditions that  will be used often 
is_cpml_xn = false; is_cpml_xp = false; 
is_cpml_yn = false; is_cpml_yp = false; 
is_any_side_cpml = false;

if strcmp(boundary.type_xn, 'cpml')
    is_cpml_xn = true;
    n_cpml_xn = abs(boundary.cpml_number_of_cells_xn);
end
if strcmp(boundary.type_xp, 'cpml')
    is_cpml_xp = true;
    n_cpml_xp = abs(boundary.cpml_number_of_cells_xp);
end
if strcmp(boundary.type_yn, 'cpml')
    is_cpml_yn = true;
    n_cpml_yn = abs(boundary.cpml_number_of_cells_yn);
end
if strcmp(boundary.type_yp, 'cpml')
    is_cpml_yp = true;
    n_cpml_yp = abs(boundary.cpml_number_of_cells_yp);
end

if (is_cpml_xn || is_cpml_xp || is_cpml_yn || is_cpml_yp)
    is_any_side_cpml = true;
end

% Call CPML initialization routine if any side is CPML
if is_any_side_cpml
    if is_TEz
        initialize_cpml_boundary_conditions_2d_TEz;
    end
    if is_TMz
        initialize_cpml_boundary_conditions_2d_TMz;
    end
end
