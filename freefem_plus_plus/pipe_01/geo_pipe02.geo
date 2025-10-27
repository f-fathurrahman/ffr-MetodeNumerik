// Parameters
R = 1.0;    // Radius of semicircles
L = 2.0;    // Vertical connector length
h = 0.1;    // Mesh size

// Points for bottom semicircle
Point(1) = {0, 0, 0, h};         // Inlet bottom
Point(2) = {0, R, 0, h};         // Inlet top
Point(3) = {2*R, R, 0, h};       // Bottom arc end

// Vertical segment up
Point(4) = {2*R, R+L, 0, h};     // Connector top

// Top semicircle (right to left)
Point(5) = {0, R+L, 0, h};       // Top arc start
Point(6) = {0, 2*R + L, 0, h};   // Outlet top
Point(7) = {2*R, 2*R + L, 0, h}; // Outlet end

// Arc centers
Point(10) = {R, R, 0, h};        // Center of bottom arc
Point(11) = {R, R+L, 0, h};      // Center of top arc

// Geometry: arcs and lines
Circle(1) = {2, 10, 3};         // Bottom arc
Line(2) = {3, 4};               // Vertical up
Circle(3) = {4, 11, 5};         // Top arc
Line(4) = {5, 6};               // Vertical down (left side)
Line(5) = {6, 7};               // Top outlet
Line(6) = {7, 4};               // Down on right

// Close loop (reversed lines to make it clockwise)
Line Loop(100) = {1, 2, 3, 4, 5, 6};
Plane Surface(200) = {100};

// Physical groups
Physical Surface("Fluid") = {200};
Physical Line("Inlet") = {4};     // Line 4: vertical left
Physical Line("Outlet") = {5};    // Line 5: top
Physical Line("Walls") = {1,2,3,6};
