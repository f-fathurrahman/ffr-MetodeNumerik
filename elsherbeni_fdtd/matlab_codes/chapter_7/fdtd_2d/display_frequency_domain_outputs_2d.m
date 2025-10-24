disp('plotting the frequency domain parameters');

% figures for sampled electric fields
for ind=1:number_of_sampled_electric_fields 
    frequencies = sampled_electric_fields(ind).frequencies;
    fd_value = sampled_electric_fields(ind).frequency_domain_value;
    figure;
    subplot(2,1,1);
    plot(frequencies*1e-9, abs(fd_value),'b-','linewidth',1.5);
    title(['sampled electric field [' num2str(ind) ']'],'fontsize',12);
    xlabel('frequency (GHz)','fontsize',12);
    ylabel('magnitude','fontsize',12);
    grid on;
    subplot(2,1,2);
    plot(frequencies*1e-9, angle(fd_value)*180/pi,'r-','linewidth',1.5);
    xlabel('frequency (GHz)','fontsize',12);
    ylabel('phase (degrees)','fontsize',12);
    grid on;
    drawnow;
end

% figures for sampled magnetic fields
for ind=1:number_of_sampled_magnetic_fields 
    frequencies = sampled_magnetic_fields(ind).frequencies;
    fd_value = sampled_magnetic_fields(ind).frequency_domain_value;
    figure;
    subplot(2,1,1);
    plot(frequencies*1e-9, abs(fd_value),'b-','linewidth',1.5);
    title(['sampled magnetic field [' num2str(ind) ']'],'fontsize',12);
    xlabel('frequency (GHz)','fontsize',12);
    ylabel('magnitude','fontsize',12);
    grid on;
    subplot(2,1,2);
    plot(frequencies*1e-9, angle(fd_value)*180/pi,'r-','linewidth',1.5);
    xlabel('frequency (GHz)','fontsize',12);
    ylabel('phase (degrees)','fontsize',12);
    grid on;
    drawnow;
end

% display sampled time harmonic electric fields on a plane
for ind=1:number_of_sampled_frequency_E_planes  
    figure;
    f = sampled_frequency_E_planes(ind).frequency;
    component = sampled_frequency_E_planes(ind).component;
    Es = abs(sampled_frequency_E_planes(ind).sampled_field);
    Es = Es/max(max(Es));
    imagesc(xcoor,ycoor,Es.');
    axis equal; axis xy; colorbar; 
    title(['Electric field at f = ' ...
        num2str(f*1e-9) ' GHz, <' component '>[' num2str(ind) ']']);
    drawnow;
end
