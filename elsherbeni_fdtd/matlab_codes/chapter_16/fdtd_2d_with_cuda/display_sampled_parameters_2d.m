% displaying sampled parameters

if mod(time_step,plotting_step) ~= 0   
    return;
end

disp([num2str(time_step) ' of ' ...
    num2str(number_of_time_steps) ' is completed.']) ; 

% display sampled electric fields
for ind = 1:number_of_sampled_electric_fields
    if sampled_electric_fields(ind).display_plot == true
    sampled_time = sampled_electric_fields(ind).time(1:time_step)*1e9;
    sampled_value = sampled_electric_fields(ind).sampled_value(1:time_step);
    figure(sampled_electric_fields(ind).figure_number);
    delete(sampled_electric_fields(ind).plot_handle);
    sampled_electric_fields(ind).plot_handle = ...
    plot(sampled_time, sampled_value(1:time_step),'b-','linewidth',1.5);
    drawnow;
    end
end

% display sampled magnetic fields
for ind = 1:number_of_sampled_magnetic_fields
    if sampled_magnetic_fields(ind).display_plot == true
    sampled_time = sampled_magnetic_fields(ind).time(1:time_step)*1e9;
    sampled_value = sampled_magnetic_fields(ind).sampled_value(1:time_step);
    figure(sampled_magnetic_fields(ind).figure_number);
    delete(sampled_magnetic_fields(ind).plot_handle);
    sampled_magnetic_fields(ind).plot_handle = ...
    plot(sampled_time, sampled_value(1:time_step),'b-','linewidth',1.5);
    drawnow;
    end
end

% display sampled electric fields on a plane
for ind=1:number_of_sampled_transient_E_planes  
    figure(sampled_transient_E_planes(ind).figure);
    Es = zeros(nxp1, nyp1);    
    component = sampled_transient_E_planes(ind).component;
    switch (component)
        case 'x'
            Es(2:nx,:) = 0.5 * (Ex(1:nx-1,:) + Ex(2:nx,:)); 
        case 'y'
            Es(:,2:ny) = 0.5 * (Ey(:,1:ny-1) + Ey(:,2:ny)); 
        case 'z'
            Es = Ez;
        case 'm'
            Exs(2:nx,:) = 0.5 * (Ex(1:nx-1,:) + Ex(2:nx,:)); 
            Eys(:,2:ny) = 0.5 * (Ey(:,1:ny-1) + Ey(:,2:ny)); 
            Ezs = Ez;
            Es = sqrt(Exs.^2 + Eys.^2 + Ezs.^2);
    end
    imagesc(xcoor,ycoor,Es.');
    axis equal; axis xy; colorbar; 
    title(['Electric field <' component '>[' num2str(ind) ']']);
    drawnow;
end

% display sampled magnetic fields on a plane
for ind=1:number_of_sampled_transient_H_planes  
    figure(sampled_transient_H_planes(ind).figure);
    Hs = zeros(nxp1, nyp1);    
    component = sampled_transient_H_planes(ind).component;
    switch (component)
        case 'x'
            Hs(:,2:ny) = 0.5 * (Hx(:,1:ny-1) + Hx(:,2:ny)); 
        case 'y'
            Hs(2:nx,:) = 0.5 * (Hy(1:nx-1,:) + Hy(2:nx,:)); 
        case 'z'
            Hs(2:nx,2:ny) = 0.25 * (Hz(2:nx,2:ny) + Hz(1:nx-1,2:ny) + Hz(2:nx,1:ny-1) + Hz(1:nx-1,1:ny-1)); 
        case 'm'
            Hsx(:,2:ny) = 0.5 * (Hx(:,1:ny-1) + Hx(:,2:ny)); 
            Hys(2:nx,:) = 0.5 * (Hy(1:nx-1,:) + Hy(2:nx,:)); 
            Hzs(2:nx,2:y) = 0.25 * (Hz(2:nx,2:ny) + Hz(1:nx-1,2:ny) + Hz(2:nx,1:ny-1) + Hz(1:nx-1,1:ny-1)); 
            Hs = sqrt(Hxs.^2 + Hys.^2 + Hzs.^2);
    end
    imagesc(xcoor,ycoor,Hs.');
    axis equal; axis xy; colorbar; 
    title(['Magnetic field <' component '>[' num2str(ind) ']']);
    drawnow;
end
