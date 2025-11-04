function show_s_on_smith(g,varargin)

numbers_on = true;
number_bcolor = 'y';

for k = 1:2:size(varargin,2)
    arg = varargin{k}; 
    val = varargin{k+1}; 
    if strcmp(arg, 'numbers')
        if strcmp(val,'off')
            numbers_on = false;
        end
    end
    if strcmp(arg, 'number_bcolor')
        number_bcolor = val;
    end    
end

hold on;

for mi=1:size(g,2)
    plot(real(g),imag(g),'--r','LineWidth',2,...
                    'MarkerEdgeColor','k',...
                    'MarkerFaceColor','g',...
                    'MarkerSize',8);
    if numbers_on 
        text(real(g),imag(g),num2str(g(mi)),'BackgroundColor',number_bcolor,'verticalalignment','bottom');
    end
end

hold off;

