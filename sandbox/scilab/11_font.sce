x = linspace(0,3,100)
y = sin(x)
xlfont("reset")
xlfont("Times New Roman",10)
plot(x,y)
xstring(0.5,0.5,"A Text from ffr")
figure_entity = gcf();
axes_entity = figure_entity.children
title_entity = axes_entity.children
title_entity.font_style = 10


