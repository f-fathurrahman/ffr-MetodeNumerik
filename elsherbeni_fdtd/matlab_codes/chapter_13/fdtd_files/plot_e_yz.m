
if animation(current_animation_index).component == 'y'
    % plot y components;
    pdata = squeeze(Ey(position_index,:,:));
    pdata(ny+1,:) = pdata(ny,:);
    pdata(2:ny,:) = 0.5 * (pdata(1:ny-1,:) + pdata(2:ny,:));
    pdata = reshape(pdata.',1,[]).';
end

if animation(current_animation_index).component == 'z'
    % plot z components;
    pdata = squeeze(Ez(position_index,:,:));
    pdata(:,nz+1) = pdata(:,nz);
    pdata(:,2:nz) = 0.5 * (pdata(:,1:nz-1) + pdata(:,2:nz));
    pdata = reshape(pdata.',1,[]).';
end

if animation(current_animation_index).component == 'x'
    % plot x components;
    pdata = 0.5*squeeze((Ex(position_index,:,:)+Ex(position_index-1,:,:)));
    pdata = reshape(pdata.',1,[]).';
end

if animation(current_animation_index).component == 'm'
    % plot magnitude;
    pdatay = squeeze(Ey(position_index,:,:));
    pdatay(ny+1,:) = pdatay(ny,:);
    pdatay(2:ny,:) = 0.5 * (pdatay(1:ny-1,:) + pdatay(2:ny,:));
    pdatay = reshape(pdatay.',1,[]).';

    pdataz = squeeze(Ez(position_index,:,:));
    pdataz(:,nz+1) = pdataz(:,nz);
    pdataz(:,2:nz) = 0.5 * (pdataz(:,1:nz-1) + pdataz(:,2:nz));
    pdataz = reshape(pdataz.',1,[]).';

    pdatax = 0.5*squeeze((Ex(position_index,:,:)+Ex(position_index-1,:,:)));
    pdatax = reshape(pdatax.',1,[]).';
   
    pdata = (pdatay.^2+pdataz.^2+pdatax.^2).^0.5;

end

color_data = [color_data;pdata];
