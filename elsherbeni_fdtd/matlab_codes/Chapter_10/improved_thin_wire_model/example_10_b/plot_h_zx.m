
if animation(current_animation_index).component == 'z'
    % plot z components;
    pdata = squeeze((Hz(:,position_index,:)+Hz(:,position_index-1,:)));
    pdata(nx+1,:) = pdata(nx,:);
    pdata(2:nx,:) = 0.5 * (pdata(1:nx-1,:) + pdata(2:nx,:));
    pdata = reshape(pdata,1,[]).';
end

if animation(current_animation_index).component == 'x'
    % plot x components;
    pdata = squeeze((Hx(:,position_index,:)+Hx(:,position_index-1,:)));
    pdata(:,nz+1) = pdata(:,nz);
    pdata(:,2:nz) = 0.5 * (pdata(:,1:nz-1) + pdata(:,2:nz));
    pdata = reshape(pdata,1,[]).';
end

if animation(current_animation_index).component == 'y'
    % plot y components;
    pdata = 0.5*squeeze((Hy(:,position_index,:)+Hy(:,position_index-1,:)));
    pdata(:,nz+1) = pdata(:,nz);
    pdata(nx+1,:) = pdata(nx,:);
    pdata(2:nx,2:nz) = 0.25 * (pdata(1:nx-1,2:nz) + pdata(2:nx,2:nz) ...       
                            + pdata(1:nx-1,1:nz-1) + pdata(2:nx,1:nz-1));
    pdata = reshape(pdata,1,[]).';
end

if animation(current_animation_index).component == 'm'
    % plot z components;
    pdataz = squeeze((Hz(:,position_index,:)+Hz(:,position_index-1,:)));
    pdataz(nx+1,:) = pdataz(nx,:);
    pdataz(2:nx,:) = 0.5 * (pdataz(1:nx-1,:) + pdataz(2:nx,:));
    pdataz = reshape(pdataz,1,[]).';

    pdatax = squeeze((Hx(:,position_index,:)+Hx(:,position_index-1,:)));
    pdatax(:,nz+1) = pdatax(:,nz);
    pdatax(:,2:nz) = 0.5 * (pdatax(:,1:nz-1) + pdatax(:,2:nz));
    pdatax = reshape(pdatax,1,[]).';

    pdatay = 0.5*squeeze((Hy(:,position_index,:)+Hy(:,position_index-1,:)));
    pdatay(:,nz+1) = pdatay(:,nz);
    pdatay(nx+1,:) = pdatay(nx,:);
    pdatay(2:nx,2:nz) = 0.25 * (pdatay(1:nx-1,2:nz) + pdatay(2:nx,2:nz) ...       
                            + pdatay(1:nx-1,1:nz-1) + pdatay(2:nx,1:nz-1));
    pdatay = reshape(pdatay,1,[]).';
    
    pdata = (pdataz.^2+pdatax.^2+pdatay.^2).^0.5;
end

color_data = [color_data;pdata];
