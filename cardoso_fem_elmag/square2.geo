// characteristic length
lc = 0.05;

// These are defined in 2d, so z coord is 0.0
Point(1) = {0.0, 0.0, 0.0, lc};
Point(2) = {1.0, 0.0, 0.0, lc};
Point(3) = {1.0, 1.0, 0.0, lc};
Point(4) = {0.0, 1.0, 0.0, lc};

Line(1) = {1, 2}; // bottom
Line(2) = {2, 3}; // right
Line(3) = {3, 4}; // top
Line(4) = {4, 1}; // left

Line Loop(10) = {1, 2, 3, 4};
Plane Surface(20) = {10};

// Physical groups
Physical Surface("Air", 1) = {20}; // 20 -> Plane Surface(20)

// Physical lines for boundary conditions
Physical Line("Bottom", 11) = {1};
Physical Line("Right", 12) = {2};
Physical Line("Top", 13) = {3};
Physical Line("Left", 14) = {4};
