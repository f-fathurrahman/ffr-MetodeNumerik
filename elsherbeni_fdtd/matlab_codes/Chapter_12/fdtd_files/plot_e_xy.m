
if animation(current_animation_index).component == 'x'
    % plot x components;
    pdata = Ex(:,:,position_index);
    pdata(nx+1,:) = pdata(nx,:);
    pdata(2:nx,:) = 0.5 * (pdata(1:nx-1,:) + pdata(2:nx,:));
    pdata = reshape(pdata.',1,[]).';
end

if animation(current_animation_index).component == 'y'
    % plot y components;
    pdata = Ey(:,:,position_index);
    pdata(:,ny+1) = pdata(:,ny);
    pdata(:,2:ny) = 0.5 * (pdata(:,1:ny-1) + pdata(:,2:ny));
    pdata = reshape(pdata.',1,[]).';
end

if animation(current_animation_index).component == 'z'
    % plot z components;
    pdata = 0.5*(Ez(:,:,position_index)+Ez(:,:,position_index-1));
    pdata = reshape(pdata.',1,[]).';
end

if animation(current_animation_index).component == 'm'
    % plot magnitude;
    pdatax = Ex(:,:,position_index);
    pdatax(nx+1,:) = pdatax(nx,:);
    pdatax(2:nx,:) = 0.5 * (pdatax(1:nx-1,:) + pdatax(2:nx,:));
    pdatax = reshape(pdatax.',1,[]).';

    pdatay = Ey(:,:,position_index);
    pdatay(:,ny+1) = pdatay(:,ny);
    pdatay(:,2:ny) = 0.5 * (pdatay(:,1:ny-1) + pdatay(:,2:ny));
    pdatay = reshape(pdatay.',1,[]).';

    pdataz = 0.5*(Ez(:,:,position_index)+Ez(:,:,position_index-1));
    pdataz = reshape(pdataz.',1,[]).';

    pdata = (pdatax.^2+pdatay.^2+pdataz.^2).^0.5;

end
color_data = [color_data;pdata];
