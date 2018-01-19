// Calling Sequence
//
// grayplot(x,y,z,[strf,rect,nax])
// grayplot(x,y,z,<opt_args>)
//
// Arguments
// 
//   x,y
// real row vectors of size n1 and n2.
//
//   z
// real matrix of size (n1,n2). z(i,j) is the value of the surface at the point (x(i),y(j)).
// 
//   <opt_args>
// This represents a sequence of statements key1=value1, key2=value2 ,... where key1, key2,...
// can be one of the following: rect, nax, strf, logflag or axesflag and frameflag (see plot2d).
// strf,rect,nax
// see plot2d.
//
// Description
// grayplot makes a 2D plot of the surface given by z on a grid defined by x and y.
// Each rectangle on the grid is filled with a gray or color level depending on the average
// value of z on the corners of the rectangle. If z contains %nan values, the suroundin
// rectangles are not displayed.

clf()
x = -10:10;
y = -10:10;
m = rand(21,21);
grayplot( x, y, m, rect=[-20,-20,20,20] )

xs2pdf( gcf(), "images/06_grayplot_v1.pdf" )

clf()
grayplot( x, y, m )
xs2pdf( gcf(), "images/06_grayplot_v2.pdf" )

if getscilabmode() ~= "STD"
  quit()
end

