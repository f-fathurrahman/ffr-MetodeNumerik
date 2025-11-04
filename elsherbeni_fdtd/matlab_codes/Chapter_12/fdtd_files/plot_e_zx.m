
if animation(current_animation_index).component == 'z'
    % plot z components;
    pdata = squeeze(Ez(:,position_index,:));
    pdata(:,nz+1) = pdata(:,nz);
    pdata(:,2:nz) = 0.5 * (pdata(:,1:nz-1) + pdata(:,2:nz));
    pdata = reshape(pdata,1,[]).';
end

if animation(current_animation_index).component == 'x'
    % plot x components;
    pdata = squeeze(Ex(:,position_index,:));
    pdata(nx+1,:) = pdata(nx,:);
    pdata(2:nx,:) = 0.5 * (pdata(1:nx-1,:) + pdata(2:nx,:));
    pdata = reshape(pdata,1,[]).';
end

if animation(current_animation_index).component == 'y'
    % plot y components;
    pdata = 0.5*squeeze((Ey(:,position_index,:)+Ey(:,position_index-1,:)));
    pdata = reshape(pdata,1,[]).';
end

if animation(current_animation_index).component == 'm'
    % plot magnitude;
    pdataz = squeeze(Ez(:,position_index,:));
    pdataz(:,nz+1) = pdataz(:,nz);
    pdataz(:,2:nz) = 0.5 * (pdataz(:,1:nz-1) + pdataz(:,2:nz));
    pdataz = reshape(pdataz,1,[]).';
   
    pdatax = squeeze(Ex(:,position_index,:));
    pdatax(nx+1,:) = pdatax(nx,:);
    pdatax(2:nx,:) = 0.5 * (pdatax(1:nx-1,:) + pdatax(2:nx,:));
    pdatax = reshape(pdatax,1,[]).';

    pdatay = 0.5*squeeze((Ey(:,position_index,:)+Ey(:,position_index-1,:)));
    pdatay = reshape(pdatay,1,[]).';

    pdata = (pdataz.^2+pdatax.^2+pdatay.^2).^0.5;

end

color_data = [color_data;pdata];
