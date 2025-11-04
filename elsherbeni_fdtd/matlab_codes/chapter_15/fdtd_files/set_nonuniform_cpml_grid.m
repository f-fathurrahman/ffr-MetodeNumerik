        %xn and xp
        % Combine e grid with h grid
        node_coordinates_e = node_coordinates_xe;
        node_coordinates_h = node_coordinates_xh;
        
        n_nodes_e = size(node_coordinates_e,2); 
        n_nodes_h = size(node_coordinates_h,2); 
        n_nodes   = n_nodes_e + n_nodes_h;
        node_coordinates_eh =  zeros(1,n_nodes);

        node_coordinates_eh(1:2:n_nodes) = node_coordinates_h;
        node_coordinates_eh(2:2:n_nodes) = node_coordinates_e;
        if strcmp(boundary.type_xn, 'cpml') && ...
                (boundary.cpml_number_of_cells_xn>0)
                ncpml = boundary.cpml_number_of_cells_xn;
                ds = node_coordinates_eh (2*ncpml+1) - node_coordinates_eh (2*ncpml);                
                for mi=2*ncpml:-1:1
                    node_coordinates_eh(mi) = node_coordinates_eh(mi+1) - ds;
                    ds = ds * boundary.nonuniform_rate;
                end            
        end
        if strcmp(boundary.type_xp, 'cpml') && ...
                (boundary.cpml_number_of_cells_xp>0)
                ncpml = boundary.cpml_number_of_cells_xp;
                ds = node_coordinates_eh (n_nodes-2*ncpml+1) - node_coordinates_eh (n_nodes-2*ncpml);                
                for mi=n_nodes-2*ncpml+1:n_nodes
                    node_coordinates_eh(mi) = node_coordinates_eh(mi-1) + ds;
                    ds = ds * boundary.nonuniform_rate;
                end            
        end

        n_nodes = size(node_coordinates_eh, 2);
        node_coordinates_h = node_coordinates_eh(1:2:n_nodes);
        node_coordinates_e = node_coordinates_eh(2:2:n_nodes);

        node_coordinates_xe = node_coordinates_e;
        node_coordinates_xh = node_coordinates_h;

        %yn and yp
        % Combine e grid with h grid
        node_coordinates_e = node_coordinates_ye;
        node_coordinates_h = node_coordinates_yh;
        
        n_nodes_e = size(node_coordinates_e,2); 
        n_nodes_h = size(node_coordinates_h,2); 
        n_nodes   = n_nodes_e + n_nodes_h;
        node_coordinates_eh =  zeros(1,n_nodes);

        node_coordinates_eh(1:2:n_nodes) = node_coordinates_h;
        node_coordinates_eh(2:2:n_nodes) = node_coordinates_e;
        if strcmp(boundary.type_yn, 'cpml') && ...
                (boundary.cpml_number_of_cells_yn>0)
                ncpml = boundary.cpml_number_of_cells_yn;
                ds = node_coordinates_eh (2*ncpml+1) - node_coordinates_eh (2*ncpml);                
                for mi=2*ncpml:-1:1
                    node_coordinates_eh(mi) = node_coordinates_eh(mi+1) - ds;
                    ds = ds * boundary.nonuniform_rate;
                end            
        end
        if strcmp(boundary.type_yp, 'cpml') && ...
                (boundary.cpml_number_of_cells_yp>0)
                ncpml = boundary.cpml_number_of_cells_yp;
                ds = node_coordinates_eh (n_nodes-2*ncpml+1) - node_coordinates_eh (n_nodes-2*ncpml);                
                for mi=n_nodes-2*ncpml+1:n_nodes
                    node_coordinates_eh(mi) = node_coordinates_eh(mi-1) + ds;
                    ds = ds * boundary.nonuniform_rate;
                end            
        end
        n_nodes = size(node_coordinates_eh, 2);
        node_coordinates_h = node_coordinates_eh(1:2:n_nodes);
        node_coordinates_e = node_coordinates_eh(2:2:n_nodes);

        node_coordinates_ye = node_coordinates_e;
        node_coordinates_h = node_coordinates_h;
        
        %zn and zp
        % Combine e grid with h grid
        node_coordinates_e = node_coordinates_ze;
        node_coordinates_h = node_coordinates_zh;
        
        n_nodes_e = size(node_coordinates_e,2); 
        n_nodes_h = size(node_coordinates_h,2); 
        n_nodes   = n_nodes_e + n_nodes_h;
        node_coordinates_eh =  zeros(1,n_nodes);

        node_coordinates_eh(1:2:n_nodes) = node_coordinates_h;
        node_coordinates_eh(2:2:n_nodes) = node_coordinates_e;
        if strcmp(boundary.type_zn, 'cpml') && ...
                (boundary.cpml_number_of_cells_zn>0)
                ncpml = boundary.cpml_number_of_cells_zn;
                ds = node_coordinates_eh (2*ncpml+1) - node_coordinates_eh (2*ncpml);                
                for mi=2*ncpml:-1:1
                    node_coordinates_eh(mi) = node_coordinates_eh(mi+1) - ds;
                    ds = ds * boundary.nonuniform_rate;
                end            
        end
        if strcmp(boundary.type_zp, 'cpml') && ...
                (boundary.cpml_number_of_cells_zp>0)
                ncpml = boundary.cpml_number_of_cells_zp;
                ds = node_coordinates_eh (n_nodes-2*ncpml+1) - node_coordinates_eh (n_nodes-2*ncpml);                
                for mi=n_nodes-2*ncpml+1:n_nodes
                    node_coordinates_eh(mi) = node_coordinates_eh(mi-1) + ds;
                    ds = ds * boundary.nonuniform_rate;
                end            
        end
        n_nodes = size(node_coordinates_eh, 2);
        node_coordinates_h = node_coordinates_eh(1:2:n_nodes);
        node_coordinates_e = node_coordinates_eh(2:2:n_nodes);
        
        node_coordinates_ze = node_coordinates_e;
        node_coordinates_zh = node_coordinates_h;
        
        
