// ------------------------------------------------------------
// Dumbbell: two spheres + connecting cylinder (OpenCASCADE)
// ------------------------------------------------------------
SetFactory("OpenCASCADE");

// Parameters
r_sphere = 0.5;      // sphere radius
r_cyl    = 0.3;      // cylinder radius
L_cyl    = 1.5;      // distance between sphere centers

// Sphere centers (cylinder will connect them)
x_left  = -L_cyl/2;
x_right =  L_cyl/2;

// Create spheres
Sphere(1) = {x_left,  0, 0, r_sphere};
Sphere(2) = {x_right, 0, 0, r_sphere};

// Create cylinder: from left sphere center to right sphere center
// (so the cylinder fully penetrates both spheres, ensuring overlap)
Cylinder(3) = {x_left, 0, 0, L_cyl, 0, 0, r_cyl};

// Boolean union – combine all three into a single volume
BooleanUnion(4) = { Volume{1}; Volume{2}; Volume{3}; };

// Optional: delete original volumes to keep the workspace clean
Delete { Volume{1}; Volume{2}; Volume{3}; }

// ------------------------------------------------------------
// Mesh size control
// ------------------------------------------------------------
Mesh.CharacteristicLengthMax = 0.1;

// Refine near the joints (optional – improves quality)
// Field[1] = Distance;
// Field[1].NodesList = {x_left,0,0, x_right,0,0};
// Field[2] = Threshold;
// Field[2].InField = 1;
// Field[2].SizeMin = 0.05;
// Field[2].SizeMax = 0.1;
// Field[2].DistMin = 0.2;
// Field[2].DistMax = 0.5;
// Background Field = 2;

// ------------------------------------------------------------
// Generate 3D tetrahedral mesh
// ------------------------------------------------------------
Mesh 3;

// Save mesh in GMSH format (version 2.2 for compatibility)
// The file will be written in the current directory.
Save "dumbbell.msh" Version 2;

// Optional: launch the GUI to view the mesh
// View "Geometry" = 1;
// View "Mesh" = 1;
// View[0].Show = 1;