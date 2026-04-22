SetFactory("OpenCASCADE");

// Lattice vectors
a1[] = {10.0, 0.0, 0.0};
a2[] = {5.0, 10.0, 0.0};
a3[] = {2.0, 3.0, 10.0};

// Base point
p0 = {0, 0, 0};

// Other vertices
p1 = {a1[0], a1[1], a1[2]};
p2 = {a2[0], a2[1], a2[2]};
p3 = {a3[0], a3[1], a3[2]};
p4 = {a1[0]+a2[0], a1[1]+a2[1], a1[2]+a2[2]};
p5 = {a1[0]+a3[0], a1[1]+a3[1], a1[2]+a3[2]};
p6 = {a2[0]+a3[0], a2[1]+a3[1], a2[2]+a3[2]};
p7 = {a1[0]+a2[0]+a3[0], a1[1]+a2[1]+a3[1], a1[2]+a2[2]+a3[2]};

// Create points
Point(1) = {p0[0], p0[1], p0[2], 1.0};
Point(2) = {p1[0], p1[1], p1[2], 1.0};
Point(3) = {p2[0], p2[1], p2[2], 1.0};
Point(4) = {p4[0], p4[1], p4[2], 1.0};

Point(5) = {p3[0], p3[1], p3[2], 1.0};
Point(6) = {p5[0], p5[1], p5[2], 1.0};
Point(7) = {p6[0], p6[1], p6[2], 1.0};
Point(8) = {p7[0], p7[1], p7[2], 1.0};

// Bottom face
Line(1) = {1,2};
Line(2) = {2,4};
Line(3) = {4,3};
Line(4) = {3,1};

Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

// Top face
Line(5) = {5,6};
Line(6) = {6,8};
Line(7) = {8,7};
Line(8) = {7,5};

Curve Loop(2) = {5,6,7,8};
Plane Surface(2) = {2};

// Vertical edges
Line(9)  = {1,5};
Line(10) = {2,6};
Line(11) = {4,8};
Line(12) = {3,7};

// Side surfaces
Surface(3) = {1,10,-5,-9};
Surface(4) = {2,11,-6,-10};
Surface(5) = {3,12,-7,-11};
Surface(6) = {4,9,-8,-12};

// Volume
Surface Loop(1) = {1,2,3,4,5,6};
Volume(1) = {1};

Periodic Surface {2} = {1} Translate {a1[0], a1[1], a1[2]};
Periodic Surface {4} = {3} Translate {a2[0], a2[1], a2[2]};
Periodic Surface {6} = {5} Translate {a3[0], a3[1], a3[2]};

Mesh 3;

Transfinite Line "*";
Transfinite Surface "*";
Transfinite Volume "*";
Recombine Volume "*";


