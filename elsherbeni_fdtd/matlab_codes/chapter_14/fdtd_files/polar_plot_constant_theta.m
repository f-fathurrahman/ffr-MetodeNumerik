 function polar_plot_constant_theta(phi,pattern_1,pattern_2, ...
        max_val, step_size, number_of_rings,...
        line_style_1, line_style_2,constant_theta, ...
        legend_1,legend_2,scale_type)

% this function plots two polar plots in the same figure
plot_range = step_size * number_of_rings;
min_val = max_val - plot_range; 

hold on;
th = 0:(pi/50):2*pi; circle_x = cos(th); circle_y = sin(th);
for mi = 1:number_of_rings
    r = (1/number_of_rings) * mi; 
    plot(r*circle_x,r*circle_y,':','color','k','linewidth',1);
    text(0.04,r,[num2str(min_val+step_size*mi)],...
        'verticalalignment','bottom','color','k',...
        'fontweight','demi','fontsize',10);
end

r=[0:0.1:1];
for mi = 0:11
    th=mi*pi/6;
    plot(r*cos(th),r*sin(th),':','color','k','linewidth',1);
text(1.1*cos(th),1.1*sin(th),[num2str(30*mi)],...
        'horizontalalignment','center','color','k',...
        'fontweight','demi','fontsize',10);
end

pattern_1(find(pattern_1 < min_val)) = min_val;  
pattern_1 = (pattern_1 - min_val)/plot_range;
pattern_2(find(pattern_2 < min_val)) = min_val;  
pattern_2 = (pattern_2 - min_val)/plot_range;

% transform data to Cartesian coordinates
x1 = pattern_1.*cos(phi);
y1 = pattern_1.*sin(phi);

x2 = pattern_2.*cos(phi);
y2 = pattern_2.*sin(phi);

% plot data on top of grid
p = plot(x1,y1,line_style_1,x2,y2,line_style_2,'linewidth',2);
text(1.2*cos(pi/4),1.2*sin(pi/4),...
    ['\theta = ' num2str(constant_theta) '^o'],...
    'color','b','fontweight','demi');
legend(p,legend_1,legend_2,'location','southeast');
text(-1, -1.1, scale_type,'fontsize',12);
text(1.02 * 1.1,  0.13 * 1.1,'\uparrow',...
    'color','b','fontweight','demi');
text(1.08 * 1.1,  0.13 * 1.1,'\phi',...
    'fontname','arial','color','b','fontweight','demi','fontsize',12);

if constant_theta == 90
    text(1.2,  0.06,'x','fontname','arial',...
        'color','b','fontweight','demi');
    text(1.2, 0,'\rightarrow','color','b','fontweight','demi');
    text(0.06,1.23,'y','fontname','arial',...
        'color','b','fontweight','demi');
    text(0,1.23,'\uparrow','color','b','fontweight','demi');
    text(1.2*cos(pi/4),1.18*sin(pi/4)-0.12,...
        'xy plane','color','b','fontweight','demi');
end

axis([-1.2 1.2 -1.2 1.2]);
axis('equal');axis('off');
hold off;
set(gcf,'PaperPositionMode','auto');
set(gca,'fontsize',12);
       

    