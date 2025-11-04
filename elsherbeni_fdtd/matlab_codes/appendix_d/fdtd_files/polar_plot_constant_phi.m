 function polar_plot_constant_phi(theta,pattern_1,pattern_2, ...
        max_val, step_size, number_of_rings,...
        line_style_1, line_style_2,constant_phi, ...
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

r=[-1:0.1:1];
for mi = 3:8
    th=mi*pi/6;
    plot(r*cos(th),r*sin(th),':','color','k','linewidth',1);
    text(1.1*cos(th),1.1*sin(th),[num2str(30*(mi-3))],...
        'horizontalalignment','center','color','k',...
        'fontweight','demi','fontsize',10);
    text(-1.1*cos(th),1.1*sin(th),[num2str(30*(mi-3))],...
        'horizontalalignment','center','color','k',...
        'fontweight','demi','fontsize',10);    
end
    text(0,-1.1,'180',...
        'horizontalalignment','center','color','k',...
        'fontweight','demi','fontsize',10);    
    
pattern_1(find(pattern_1 < min_val)) = min_val;  
pattern_1 = (pattern_1 - min_val)/plot_range;
pattern_2(find(pattern_2 < min_val)) = min_val;  
pattern_2 = (pattern_2 - min_val)/plot_range;

% transform data to Cartesian coordinates
x1 = -pattern_1.*cos(theta+pi/2);
y1 = pattern_1.*sin(theta+pi/2);

x2 = -pattern_2.*cos(theta+pi/2);
y2 = pattern_2.*sin(theta+pi/2);

% plot data on top of grid
p = plot(x1,y1,line_style_1,x2,y2,line_style_2,'linewidth',2);
text(1.2*cos(pi/4),1.2*sin(pi/4),...
    ['\phi = ' num2str(constant_phi) '^o'],...
    'color','b','fontweight','demi');
legend(p,legend_1,legend_2,'location','southeast');
text(-1, -1.1, scale_type,'fontsize',12);
text(0.2,1.02,'\rightarrow','color','b','fontweight','demi');
text(0.2, 1.08,'\theta','fontname','arial','color','b',...
    'fontweight','demi','fontsize',12);
text(-0.21,1.02,'\leftarrow','color','b','fontweight','demi');
text(-0.2, 1.08,'\theta','fontname','arial','color','b',...
    'fontweight','demi','fontsize',12);

if constant_phi == 0
    text(1.2,0.06,'x','fontname','arial','color','b','fontweight','demi');
    text(1.2,0,'\rightarrow','color','b','fontweight','demi');
    text(0.06,1.23,'z','fontname','arial','color','b','fontweight','demi');
    text(0,1.23,'\uparrow','color','b','fontweight','demi');
    text(1.2*cos(pi/4),1.18*sin(pi/4)-0.12,'xz plane',...
        'color','b','fontweight','demi');
end
if constant_phi == 90
    text(1.2,0.06,'y','fontname','arial','color','b','fontweight','demi');
    text(1.2, 0,'\rightarrow','color','b','fontweight','demi');
    text(0.06,1.23,'z','fontname','arial','color','b','fontweight','demi');
    text(0,1.23,'\uparrow','color','b','fontweight','demi');
    text(1.2*cos(pi/4),1.18*sin(pi/4)-0.12,'yz plane',...
        'color','b','fontweight','demi');
end

axis([-1.2 1.2 -1.2 1.2]);
axis('equal');axis('off');
hold off;
set(gcf,'PaperPositionMode','auto');
set(gca,'fontsize',12);
