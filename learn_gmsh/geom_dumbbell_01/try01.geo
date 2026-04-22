SetFactory("OpenCASCADE");

// Parameters
r_sphere = 1.0;
r_cyl    = 0.4;
L        = 3.0;   // distance between sphere centers

// Sphere centers
x1 = -L/2;
x2 =  L/2;

// --------------------------------------------------
// Create primitives
// --------------------------------------------------

// Two spheres
Sphere(1) = {x1, 0, 0, r_sphere};
Sphere(2) = {x2, 0, 0, r_sphere};

// Cylinder connecting them
Cylinder(3) = {x1, 0, 0, L, 0, 0, r_cyl};

// --------------------------------------------------
// Remove overlaps using BooleanFragments
// (this automatically splits and removes intersections)
// --------------------------------------------------

BooleanFragments{ Volume{1,2,3}; Delete; }{}

// --------------------------------------------------
// Optional: fuse into a single volume (clean dumbbell)
// --------------------------------------------------

v[] = Volume{:};
BooleanUnion{ Volume{v[]}; Delete; }{}

// --------------------------------------------------
// Mesh settings
// --------------------------------------------------

Mesh.CharacteristicLengthMin = 0.3;
Mesh.CharacteristicLengthMax = 0.4;

// Generate 3D mesh
Mesh 3;

Save "dumbbell.vtk";
