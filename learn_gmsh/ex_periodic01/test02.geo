SetFactory("OpenCASCADE");

// Lattice vectors
a1[] = {10.0, 0.0, 0.0};
a2[] = {5.0, 10.0, 0.0};
a3[] = {2.0, 3.0, 10.0};

// Base rectangle (spanned by a1, a2)
Point(1) = {0,0,0,1.0};
Point(2) = {a1[0], a1[1], a1[2],1.0};
Point(3) = {a1[0]+a2[0], a1[1]+a2[1], a1[2]+a2[2],1.0};
Point(4) = {a2[0], a2[1], a2[2],1.0};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Curve Loop(1) = {1,2,3,4};
Plane Surface(1) = {1};

out[] = Extrude {a3[0], a3[1], a3[2]} {
  Surface{1};
};

// a3 direction (top ↔ bottom)
Periodic Surface {out[0]} = {1} Translate {a3[0], a3[1], a3[2]};


Periodic Surface {side2} = {side1} Translate {a1[]};
Periodic Surface {side4} = {side3} Translate {a2[]};


