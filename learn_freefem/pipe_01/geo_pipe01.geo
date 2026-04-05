// Parameters
R = 1.0;   // Bend radius
L = 2.0;   // Vertical length
h = 0.1;   // Mesh size

// Points (clockwise)
Point(1) = {0, 0, 0, h};                  // Inlet bottom
Point(2) = {0, R, 0, h};                  // Inlet top
Point(3) = {2*R, R, 0, h};                // Mid bend right
Point(4) = {2*R, R+L, 0, h};              // Top vertical connector
Point(5) = {0, R+L, 0, h};                // Top left
Point(6) = {0, 2*R + L, 0, h};            // Outlet top
Point(7) = {2*R, 2*R + L, 0, h};          // Outlet right

// Arcs and segments
Circle(1) = {2, {R, R, 0}, 3};            // Bottom bend
Line(2) = {3, 4};                          // Vertical up
Circle(3) = {4, {R, R+L, 0}, 5};          // Top bend
Line(4) = {5, 6};                          // Vertical down (left)
Line(5) = {6, 7};                          // Top outlet
Line(6) = {7, 4};                          // Down right

// Line loop and surface
Line Loop(10) = {1, 2, 3, 4, 5, 6};
Plane Surface(11) = {10};

// Physical groups
Physical Surface("Fluid") = {11};
Physical Line("Inlet") = {4};
Physical Line("Outlet") = {5};
Physical Line("Walls") = {1, 2, 3, 6};
