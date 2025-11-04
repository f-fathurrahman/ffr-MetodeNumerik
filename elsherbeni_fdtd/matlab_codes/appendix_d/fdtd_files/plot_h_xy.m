
if animation(current_animation_index).component == 'x'
    % plot x components;
    pdata = 0.5*(Hx(:,:,position_index)+Hx(:,:,position_index-1));
    pdata(:,ny+1) = pdata(:,ny);
    pdata(:,2:ny) = 0.5 * (pdata(:,1:ny-1) + pdata(:,2:ny));
    pdata = reshape(pdata.',1,[]).';
end

if animation(current_animation_index).component == 'y'
    % plot y components;
    pdata = 0.5*(Hy(:,:,position_index)+Hy(:,:,position_index-1));
    pdata(nx+1,:) = pdata(nx,:);
    pdata(2:nx,:) = 0.5 * (pdata(1:nx-1,:) + pdata(2:nx,:));
    pdata = reshape(pdata.',1,[]).';
end

if animation(current_animation_index).component == 'z'
    % plot z components;
    pdata = Hz(:,:,position_index);
    pdata(:,ny+1) = pdata(:,ny);
    pdata(nx+1,:) = pdata(nx,:);
    pdata(2:nx,2:ny) = 0.25 * (pdata(1:nx-1,2:ny) + pdata(2:nx,2:ny) ...       
                            + pdata(1:nx-1,1:ny-1) + pdata(2:nx,1:ny-1));
    pdata = reshape(pdata.',1,[]).';
end

if animation(current_animation_index).component == 'm'
    % plot magnitude;
    pdatax = 0.5*(Hx(:,:,position_index)+Hx(:,:,position_index-1));
    pdatax(:,ny+1) = pdatax(:,ny);
    pdatax(:,2:ny) = 0.5 * (pdatax(:,1:ny-1) + pdatax(:,2:ny));
    pdatax = reshape(pdatax.',1,[]).';

    pdatay = 0.5*(Hy(:,:,position_index)+Hy(:,:,position_index-1));
    pdatay(nx+1,:) = pdatay(nx,:);
    pdatay(2:nx,:) = 0.5 * (pdatay(1:nx-1,:) + pdatay(2:nx,:));
    pdatay = reshape(pdatay.',1,[]).';

    pdataz = Hz(:,:,position_index);
    pdataz(:,ny+1) = pdataz(:,ny);
    pdataz(nx+1,:) = pdataz(nx,:);
    pdataz(2:nx,2:ny) = 0.25 * (pdataz(1:nx-1,2:ny) + pdataz(2:nx,2:ny) ...       
                            + pdataz(1:nx-1,1:ny-1) + pdataz(2:nx,1:ny-1));
    pdataz = reshape(pdataz.',1,[]).';

    pdata = (pdatax.^2+pdatay.^2+pdataz.^2).^0.5;

end
color_data = [color_data;pdata];
