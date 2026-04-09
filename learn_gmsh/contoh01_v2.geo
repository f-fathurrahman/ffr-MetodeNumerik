lc = 0.1;

Point(1) = {0.0, 0.0, 0, lc};
Point(2) = {0.2, 0.0, 0, lc};
Point(3) = {0.3, 0.2, 0, lc};
Point(4) = {0.1, 0.4, 0, lc};
Point(5) = {-0.1, 0.2, 0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 1};

Curve Loop(1) = {1, 2, 3, 4, 5};
Plane Surface(1) = {1};

//
Physical Curve(5) = {1, 2, 4};
Physical Surface("My surface") = {1};

