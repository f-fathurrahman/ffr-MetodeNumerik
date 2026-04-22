SetFactory("OpenCASCADE");

r_sphere = 1.0;
r_cyl    = 0.4;
L        = 3.0;

x1 = -L/2;
x2 =  L/2;

// Primitives
Sphere(1) = {x1, 0, 0, r_sphere};
Sphere(2) = {x2, 0, 0, r_sphere};
Cylinder(3) = {x1, 0, 0, L, 0, 0, r_cyl};

// --------------------------------------------------
// STEP 1: Fragment everything
// --------------------------------------------------
v() = BooleanFragments{ Volume{1,2,3}; Delete; }{};

// --------------------------------------------------
// STEP 2: Keep only the OUTER volume
// --------------------------------------------------
v_keep[] = Volume In BoundingBox{
    -L, -r_sphere, -r_sphere,
     L,  r_sphere,  r_sphere
};

// Fuse to remove internal partitions
BooleanUnion{ Volume{v_keep[]}; Delete; }{}

// Clean topology
Coherence;
RemoveAllDuplicates;

Mesh.RemoveDuplicateElements = 1;
Mesh.RemoveDuplicateNodes = 1;


Mesh 3;
Save "dumbbell.vtk";

