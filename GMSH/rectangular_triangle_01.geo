// define new variable
lc = 0.1;

Point(1) = {0.0, 0.0, 0.0, lc};
Point(2) = {3.0, 0.0, 0.0, lc};
Point(3) = {3.0, 3.0, 0.0, lc};
Point(4) = {0.0, 3.0, 0.0, lc};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {4, 1, 2, 3};

Point(5) = {1.0, 1.0, 0.0, lc};
Point(6) = {1.0, 2.0, 0.0, lc};
Point(7) = {2.0, 2.0, lc};

Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 5};

Curve Loop(2) = {5, 6, 7};

Plane Surface(1) = {1, 2};

//Physical Curve(5) = {1, 2, 4};
//Physical Surface("My surface") = {1};

// Mesh 2;
// Save "TEMP_t1.msh";
// Save "TEMP_t1.m"; // Matlab script

