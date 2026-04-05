// Characteristic mesh sizes
lc_inner = 0.1; // Finer mesh for the capacitor
lc_outer = 0.5; // Coarser mesh for the external domain

// Points for the inner rectangle (capacitor)
Point(1) = {-2, -0.5, 0, lc_inner};
Point(2) = { 2, -0.5, 0, lc_inner};
Point(3) = { 2, 0.5, 0, lc_inner};
Point(4) = {-2, 0.5, 0, lc_inner};

// Points for the outer rectangle (background domain)
Point(5) = {-6, -3, 0, lc_outer};
Point(6) = { 6, -3, 0, lc_outer};
Point(7) = { 6, 3, 0, lc_outer};
Point(8) = {-6, 3, 0, lc_outer};

// Lines for the inner rectangle
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5}; // Inner surface (capacitor)

// Lines for the outer rectangle
Line(7) = {5, 6};
Line(8) = {6, 7};
Line(9) = {7, 8};
Line(10) = {8, 5};
Line Loop(11) = {7, 8, 9, 10};
Plane Surface(12) = {11, 5}; // External domain minus capacitor (fragmentation)

// Physical surfaces
Physical Surface("Insulator", 1) = {6}; // Capacitor dielectric
Physical Surface("External", 2) = {12}; // Surrounding space

// Physical lines for boundary conditions
Physical Line("V0", 13) = {3}; // Top plate of capacitor
Physical Line("GND", 14) = {1}; // Bottom plate of capacitor
