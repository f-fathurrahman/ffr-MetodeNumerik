SetFactory("OpenCASCADE");
r_sphere = 0.5; r_cyl = 0.3; L_cyl = 1.5;
x_left = -L_cyl/2; x_right = L_cyl/2;

Sphere(1) = {x_left,0,0,r_sphere};
Sphere(2) = {x_right,0,0,r_sphere};
Cylinder(3) = {x_left,0,0,L_cyl,0,0,r_cyl};

// Fragment everything – splits into regions: left sphere only, right sphere only,
// cylinder only (outside spheres), and overlapping parts.
//BooleanFragments(4) = { Volume{1,2,3}; Delete; };

// Now you have multiple volumes. Delete the ones that are inside the spheres.
// You can identify them by their centroids.
// For example, volume 5 might be the left sphere part, volume 6 the right sphere part,
// volume 7 the cylinder outside, and volumes 8 and 9 the cylinder inside spheres.
// Delete those internal fragments.
//Delete { Volume{8,9}; } // adjust indices based on your run

// Then union the remaining volumes into one.
//BooleanUnion(10) = { Volume{5,6,7}; };

Mesh 3;

