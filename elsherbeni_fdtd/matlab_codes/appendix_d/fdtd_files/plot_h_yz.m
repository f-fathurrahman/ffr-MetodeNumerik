
if animation(current_animation_index).component == 'y'
    % plot y components;
    pdata = squeeze(0.5*(Hy(position_index,:,:)+Hy(position_index-1,:,:)));
    pdata(:,nz+1) = pdata(:,nz);
    pdata(:,2:nz) = 0.5 * (pdata(:,1:nz-1) + pdata(:,2:nz));
    pdata = reshape(pdata.',1,[]).';
end

if animation(current_animation_index).component == 'z'
    % plot z components;
    pdata = squeeze(0.5*(Hz(position_index,:,:)+Hz(position_index-1,:,:)));
    pdata(ny+1,:) = pdata(ny,:);
    pdata(2:ny,:) = 0.5 * (pdata(1:ny-1,:) + pdata(2:ny,:));
    pdata = reshape(pdata.',1,[]).';
end

if animation(current_animation_index).component == 'x'
    % plot x components;
    pdata = squeeze(Hx(position_index,:,:));
    pdata(:,nz+1) = pdata(:,nz);
    pdata(ny+1,:) = pdata(ny,:);
    pdata(2:ny,2:nz) = 0.25 * (pdata(1:ny-1,2:nz) + pdata(2:ny,2:nz) ...       
                            + pdata(1:ny-1,1:nz-1) + pdata(2:ny,1:nz-1));
    pdata = reshape(pdata.',1,[]).';
end

if animation(current_animation_index).component == 'm'
    % plot magnitude;

    pdatay = squeeze(0.5*(Hy(position_index,:,:)+Hy(position_index-1,:,:)));
    pdatay(:,nz+1) = pdatay(:,nz);
    pdatay(:,2:nz) = 0.5 * (pdatay(:,1:nz-1) + pdatay(:,2:nz));
    pdatay = reshape(pdatay.',1,[]).';
    
    pdataz = squeeze(0.5*(Hz(position_index,:,:)+Hz(position_index-1,:,:)));
    pdataz(ny+1,:) = pdataz(ny,:);
    pdataz(2:ny,:) = 0.5 * (pdataz(1:ny-1,:) + pdataz(2:ny,:));
    pdataz = reshape(pdataz.',1,[]).';
    
    pdatax = squeeze(Hx(position_index,:,:));
    pdatax(:,nz+1) = pdatax(:,nz);
    pdatax(ny+1,:) = pdatax(ny,:);
    pdatax(2:ny,2:nz) = 0.25 * (pdatax(1:ny-1,2:nz) + pdatax(2:ny,2:nz) ...       
                            + pdatax(1:ny-1,1:nz-1) + pdatax(2:ny,1:nz-1));
    pdatax = reshape(pdatax.',1,[]).';
    pdata = (pdatay.^2+pdataz.^2+pdatax.^2).^0.5;
end

color_data = [color_data;pdata];
