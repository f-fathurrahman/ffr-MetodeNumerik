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

cellVol = out[1];

r1 = 2.5; // corner sphere
r2 = 2;  // center sphere

// Corner sphere
Sphere(100) = {0, 0, 0, r1};

// Center position
xc = 0.5*(a1[0] + a2[0] + a3[0]);
yc = 0.5*(a1[1] + a2[1] + a3[1]);
zc = 0.5*(a1[2] + a2[2] + a3[2]);

Sphere(200) = {xc, yc, zc, r2};

cornerSpheres[] = {100};

For i In {0:1}
  For j In {0:1}
    For k In {0:1}
      If(i != 0 || j != 0 || k != 0)
        dx = i*a1[0] + j*a2[0] + k*a3[0];
        dy = i*a1[1] + j*a2[1] + k*a3[1];
        dz = i*a1[2] + j*a2[2] + k*a3[2];

        newv = newv;
        Sphere(newv) = {dx, dy, dz, r1};
        cornerSpheres[] += {newv};
      EndIf
    EndFor
  EndFor
EndFor

allSpheres[] = cornerSpheres[];
allSpheres[] += {200};

out[] = BooleanFragments{
  Volume{cellVol};
}{
  Volume{allSpheres[]};
};

Physical Volume("Matrix") = {out[0]}; // cell minus spheres
Physical Volume("Inclusions") = {out[1:#out[]-1]};


