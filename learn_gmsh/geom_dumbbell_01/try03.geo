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

// --------------------------------------------------
// KEY FIX: true union (removes internal overlaps)
// --------------------------------------------------
BooleanUnion{ Volume{1,2,3}; Delete; }{}
Coherence;
RemoveAllDuplicates;


// --------------------------------------------------
// Identify boundary surfaces
// --------------------------------------------------

s_left[] = Surface In BoundingBox{
    x1 - r_sphere - 1e-3, -r_sphere, -r_sphere,
    x1 + r_sphere + 1e-3,  r_sphere,  r_sphere
};

s_right[] = Surface In BoundingBox{
    x2 - r_sphere - 1e-3, -r_sphere, -r_sphere,
    x2 + r_sphere + 1e-3,  r_sphere,  r_sphere
};

s_cyl[] = Surface In BoundingBox{
    x1 + 0.2, -r_cyl-1e-3, -r_cyl-1e-3,
    x2 - 0.2,  r_cyl+1e-3,  r_cyl+1e-3
};

// Physical groups
Physical Surface("left_sphere")  = {s_left[]};
Physical Surface("right_sphere") = {s_right[]};
Physical Surface("cylinder")     = {s_cyl[]};

Physical Volume("domain") = Volume{:};

// Mesh
Mesh.CharacteristicLengthMin = 0.1;
Mesh.CharacteristicLengthMax = 0.2;

Mesh 3;
Save "dumbbell.vtk";

