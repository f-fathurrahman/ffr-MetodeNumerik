SetFactory("OpenCASCADE");

r_sphere = 1.0;
r_cyl    = 0.4;
L        = 3.0;

x1 = -L/2;
x2 =  L/2;

// Geometry
Sphere(1) = {x1, 0, 0, r_sphere};
Sphere(2) = {x2, 0, 0, r_sphere};
Cylinder(3) = {x1, 0, 0, L, 0, 0, r_cyl};

// Fragment and clean
BooleanFragments{ Volume{1,2,3}; Delete; }{}

// --------------------------------------------------
// Identify surfaces via bounding boxes
// --------------------------------------------------

// Left sphere surfaces
s_left[] = Surface In BoundingBox{
    x1 - r_sphere - 1e-3, -r_sphere, -r_sphere,
    x1 + r_sphere + 1e-3,  r_sphere,  r_sphere
};

// Right sphere surfaces
s_right[] = Surface In BoundingBox{
    x2 - r_sphere - 1e-3, -r_sphere, -r_sphere,
    x2 + r_sphere + 1e-3,  r_sphere,  r_sphere
};

// Cylinder surfaces (middle region)
s_cyl[] = Surface In BoundingBox{
    x1 + 0.2, -r_cyl-1e-3, -r_cyl-1e-3,
    x2 - 0.2,  r_cyl+1e-3,  r_cyl+1e-3
};

// --------------------------------------------------
// Define physical groups (CRUCIAL for VTK)
// --------------------------------------------------

Physical Surface("left_sphere")  = {s_left[]};
Physical Surface("right_sphere") = {s_right[]};
Physical Surface("cylinder")     = {s_cyl[]};

Physical Volume("domain") = Volume{:};

// Mesh
Mesh.CharacteristicLengthMin = 0.1;
Mesh.CharacteristicLengthMax = 0.2;

Mesh 3;
Save "dumbbell.vtk";
