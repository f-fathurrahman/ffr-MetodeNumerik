if ~exist('animation','var')
    return;
end

% display animated fields
for mi=1:size(animation,2)
    current_animation_index = mi;
    f = figure(animation(mi).figure_number);
    axis equal;
    po = findobj(gca,'type','patch');
    set_color_axis_scaling;
    color_data = [];
    for mj = 1:size(animation(mi).plane_cut,2)
        position_index = animation(mi).plane_cut(mj).position_index;
        switch animation(mi).field_type  
            case 'e'
                switch animation(mi).plane_cut(mj).type 
                    case 'xy'
                        plot_e_xy;
                    case 'yz'
                        plot_e_yz;
                    case 'zx'
                        plot_e_zx;
                end
            case 'h'
                switch animation(mi).plane_cut(mj).type 
                    case 'xy'
                        plot_h_xy;
                    case 'yz'
                        plot_h_yz;
                    case 'zx'
                        plot_h_zx;
                end
        end
    end
    patch('faces',animation(mi).faces,'vertices', ...
    animation(mi).vertices,'facecolor','interp',...
    'FaceLighting','flat','FaceVertexCData',color_data,...
    'edgecolor',animation(mi).edgecolor);

    delete(po);
    if animation(mi).display_objects
        display_objects_mesh_in_animation;
    end
    colorbar;
    drawnow;
end
