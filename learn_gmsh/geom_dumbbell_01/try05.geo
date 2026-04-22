SetFactory("OpenCASCADE");

r_sphere = 1.0;
r_cyl    = 0.4;
L        = 3.0;

x1 = -L/2;
x2 =  L/2;

// --------------------------------------------------
// Step 1: Create primitives
// --------------------------------------------------

Sphere(1) = {x1, 0, 0, r_sphere};
Sphere(2) = {x2, 0, 0, r_sphere};
Cylinder(3) = {x1, 0, 0, L, 0, 0, r_cyl};

// --------------------------------------------------
// Step 2: CUT cylinder by spheres
// (this removes cylinder parts inside spheres)
// --------------------------------------------------

cyl_cut[] = BooleanDifference{
    Volume{3}; Delete;
}{
    Volume{1,2};
};

// --------------------------------------------------
// Step 3: Fuse everything into one domain
// --------------------------------------------------

BooleanUnion{
    Volume{1,2,cyl_cut[]}; Delete;
}{}

// --------------------------------------------------
// Clean geometry
// --------------------------------------------------

Coherence;
RemoveAllDuplicates;

// Mesh
Mesh.CharacteristicLengthMin = 0.1;
Mesh.CharacteristicLengthMax = 0.2;

Mesh 3;
Save "dumbbell.vtk";