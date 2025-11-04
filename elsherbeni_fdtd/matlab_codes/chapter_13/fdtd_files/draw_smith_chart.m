function draw_smith_chart(varargin)

resolution = 0.5; 
numbers_on = true;
number_bcolor = 'none';
show_admittance = false;

for k = 1:2:size(varargin,2)
    arg = varargin{k}; 
    val = varargin{k+1}; 
    if strcmp(arg, 'numbers')
        if strcmp(val,'off')
            numbers_on = false;
        end
    end
    if strcmp(arg, 'show_admittance')
        if strcmp(val,'on')
            show_admittance = true;
        end
    end
    if strcmp(arg, 'number_bcolor')
        number_bcolor = val;
    end    
    if strcmp(arg, 'resolution')
        resolution = val;
    end    
end

n = round(2/resolution);
hold on;

if show_admittance
    adcolor = [0 1 0];
    yr = logspace(-3,4,100);
    for mi=-n:n
        yi = resolution * mi;
        y = (yr+j*yi);
        g = (1./y-1)./(1./y+1);
        plot(real(g),imag(g),'color',adcolor);
    end
    yi = logspace(-3,4,100);
    for mi=0:n
        yr = resolution * mi;
        y = (yr+j*yi);
        g = (1./y-1)./(1./y+1);
        plot(real(g),imag(g),'color',adcolor);
    end
    yi = -yi;
    for mi=0:n
        yr = resolution * mi;
        y = (yr+j*yi);
        g = (1./y-1)./(1./y+1);
        plot(real(g),imag(g),'color',adcolor);
    end
end

zr = logspace(-3,4,100);
for mi=-n:n
    zi = resolution * mi;
    z = (zr+j*zi);
    g = (z-1)./(z+1);
    plot(real(g),imag(g),'k');
end
zi = logspace(-3,4,100);
for mi=0:n
    zr = resolution * mi;
    z = (zr+j*zi);
    g = (z-1)./(z+1);
    plot(real(g),imag(g),'k');
end
zi = -zi;
for mi=0:n
    zr = resolution * mi;
    z = (zr+j*zi);
    g = (z-1)./(z+1);
    plot(real(g),imag(g),'k');
end

if numbers_on 
    for mi=0:n
        z = resolution * mi;
        g = (z-1)./(z+1);
        text(real(g),imag(g),num2str(z),'BackgroundColor',number_bcolor,'HorizontalAlignment','center','color','b','fontweight','bold');
    end
    for mi=-n:n
        z = j*resolution * mi;
        g = (z-1)./(z+1);
        text(real(g),imag(g),num2str(imag(z)),'BackgroundColor',number_bcolor,'HorizontalAlignment','center','color','b','fontweight','bold');
    end
end



hold off;
axis equal;
axis off;
