SetFactory("OpenCASCADE");

// ==============================
// 1. LATTICE VECTORS
// ==============================
a1[] = {1.0, 0.0, 0.0};
a2[] = {0.5, 1.0, 0.0};
a3[] = {0.2, 0.3, 1.0};

// ==============================
// 2. BASE PARALLELOGRAM (a1, a2)
// ==============================
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

// ==============================
// 3. EXTRUDE → UNIT CELL
// ==============================
out[] = Extrude {a3[0], a3[1], a3[2]} {
  Surface{1};
};

cellVol = out[1];

// ==============================
// 4. SPHERES
// ==============================

// Radii
r_corner = 0.25;
r_center = 0.2;

// ---- Corner sphere (origin)
Sphere(100) = {0, 0, 0, r_corner};

// ---- Center sphere
xc = 0.5*(a1[0] + a2[0] + a3[0]);
yc = 0.5*(a1[1] + a2[1] + a3[1]);
zc = 0.5*(a1[2] + a2[2] + a3[2]);

Sphere(200) = {xc, yc, zc, r_center};

// ==============================
// 5. REPLICATE CORNER SPHERE
// ==============================
cornerSpheres[] = {100};

For i In {-1:1}
  For j In {-1:1}
    For k In {-1:1}
      If(i != 0 || j != 0 || k != 0)
        dx = i*a1[0] + j*a2[0] + k*a3[0];
        dy = i*a1[1] + j*a2[1] + k*a3[1];
        dz = i*a1[2] + j*a2[2] + k*a3[2];

        s = newv;
        Sphere(s) = {dx, dy, dz, r_corner};
        cornerSpheres[] += {s};
      EndIf
    EndFor
  EndFor
EndFor

allSpheres[] = cornerSpheres[];
allSpheres[] += {200};

// ==============================
// 6. CUT SPHERES WITH CELL
// ==============================

// Keep only parts INSIDE the cell
inclusions[] = BooleanIntersection{
  Volume{allSpheres[]};
}{
  Volume{cellVol};
};

// Remove original full spheres
Recursive Delete {
  Volume{allSpheres[]};
}

// ==============================
// 7. MATRIX (CELL MINUS SPHERES)
// ==============================
matrix[] = BooleanDifference{
  Volume{cellVol};
}{
  Volume{inclusions[]};
};

// ==============================
// 8. CONFORMAL FRAGMENT
// ==============================
final[] = BooleanFragments{
  Volume{matrix[], inclusions[]};
}{};

// ==============================
// 9. PERIODIC BOUNDARIES
// ==============================

// NOTE: surface IDs depend on extrusion
// You should verify in GUI (Visibility tool)

// a3 direction (top ↔ bottom)
Periodic Surface {out[0]} = {1} Translate {a3[0], a3[1], a3[2]};

// For a1 and a2:
// Identify surfaces visually and replace IDs below

// Example placeholders:
// Periodic Surface {Sx2} = {Sx1} Translate {a1[0], a1[1], a1[2]};
// Periodic Surface {Sy2} = {Sy1} Translate {a2[0], a2[1], a2[2]};

// ==============================
// 10. PHYSICAL GROUPS
// ==============================
Physical Volume("Matrix") = {matrix[]};
Physical Volume("Inclusions") = {inclusions[]};

// ==============================
// 11. MESH
// ==============================
Mesh 3;