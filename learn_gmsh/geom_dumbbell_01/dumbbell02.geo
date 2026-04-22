SetFactory("OpenCASCADE");

r_sphere = 0.5;
r_cyl    = 0.3;
L_cyl    = 1.5;

x_left  = -L_cyl/2;
x_right =  L_cyl/2;

// Create spheres
Sphere(1) = {x_left,  0, 0, r_sphere};
Sphere(2) = {x_right, 0, 0, r_sphere};

// Create a longer cylinder that fully passes through both spheres
// (extend beyond sphere centers to ensure cutting works)
Cylinder(3) = {x_left - 0.2, 0, 0, L_cyl + 0.4, 0, 0, r_cyl};

// Cut the cylinder with the spheres: keep only the part that is OUTSIDE the spheres
// First, get the part of cylinder that is inside sphere 1
BooleanIntersection(4) = { Volume{3}; Volume{1}; };
// Then, get the part of cylinder that is inside sphere 2
BooleanIntersection(5) = { Volume{3}; Volume{2}; };
// Subtract both inside parts from the original cylinder
BooleanDifference(6) = { Volume{3}; Delete; Volume{4,5}; };
// Now Volume 6 is the cylinder with both ends trimmed to the sphere surfaces

// Union the trimmed cylinder with the two spheres
BooleanUnion(7) = { Volume{1}; Volume{2}; Volume{6}; };

Delete { Volume{1,2,3,4,5,6}; }

Mesh.CharacteristicLengthMax = 0.1;
Mesh 3;
Save "dumbbell_trimmed.msh" Version 2;

// After Mesh 3 and Save
// Launch the GUI with the mesh loaded
// (Only works when running from command line without -nopopup)
View "Mesh" = 1;
View[0].Show = 1;      // Show the mesh view
View[0].Visible = 1;
View[0].ShowElements = 3;  // Show 3D elements (tetrahedra)

