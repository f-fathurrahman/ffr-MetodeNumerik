function draw_particle(s, x, y)
%  PURPOSE: Plot the polygonal particles and walls
%  REVISION HISTORY: S. H. Ng, HG Matuttis 09-May-2014
%  TODO: Plot particles and walls in different color

n=length(x);

hold off
hold on       
for i = 1:n
  ncorner = s(i);
  patch(x(1:ncorner,i),y(1:ncorner,i),'g');
% remove the comment to annotate the corner number 
% for each particle
%  text(x(1:ncorner,i),y(1:ncorner,i),num2str([1:ncorner]'));
end

return