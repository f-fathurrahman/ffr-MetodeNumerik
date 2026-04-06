SetFactory("OpenCASCADE");

// Mesh sizes
lc_core = 0.5;
lc_coil = 1.0;

// Iron core
Point(1) = {-17, -24.8, 0, lc_core};
Point(2) = {23, -24.8, 0, lc_core};
Point(3) = {23, -0.3, 0, lc_core};
Point(4) = {13, -0.3, 0, lc_core};
Point(5) = {13, -15.3, 0, lc_core};
Point(6) = {-7, -15.3, 0, lc_core};
Point(7) = {-7, 14.7, 0, lc_core};
Point(8) = {13, 14.7, 0, lc_core};
Point(9) = {13, 1.7, 0, lc_core};
Point(10) = {23, 1.7, 0, lc_core};
Point(11) = {23, 24.7, 0, lc_core};
Point(12) = {-17, 24.7, 0, lc_core};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 1};

Curve Loop(1) = {1:12};
Plane Surface(1) = {1};

// -------------------------------
// Positive coil (Coil+)
// -------------------------------
Point(13) = {-7, -14.8, 0, lc_coil};


Point(14) = {-2, -14.8, 0, lc_coil};
Point(15) = {-2, 14.7, 0, lc_coil};
Point(16) = {-7, 14.7, 0, lc_coil};
Line(13) = {13, 14};
Line(14) = {14, 15};
Line(15) = {15, 16};
Line(16) = {16, 13};
Curve Loop(2) = {13,14,15,16};
Plane Surface(2) = {2};

// Negative coil (Coil-)
Point(17) = {-22, -15.3, 0, lc_coil};
Point(18) = {-17, -15.3, 0, lc_coil};
Point(19) = {-17, 14.7, 0, lc_coil};
Point(20) = {-22, 14.7, 0, lc_coil};
Line(17) = {17, 18};
Line(18) = {18, 19};
Line(19) = {19, 20};
Line(20) = {20, 17};
Curve Loop(3) = {17,18,19,20};
Plane Surface(3) = {3};

// Outer air domain
Disk(100) = {0, 0, 0, 100, 100}; // Surface ID 100

// Fragmentation
out[] = BooleanFragments{
    Surface{100};
    Delete;
}{
    Surface{1,2,3};
    Delete;
};

// Physical groups — adjust as needed after inspection
Physical Surface("Air") = { 4 };
Physical Surface("Iron") = { 1 };
Physical Surface("Coil+") = { 3 };
Physical Surface("Coil-") = { 2 };

// Outer boundary: only curves belonging to the "Air" region
// Manually defined here based on visual inspection of the geometry
// Adjust the IDs according to your final geometry
Physical Curve("A0") = { 19 };