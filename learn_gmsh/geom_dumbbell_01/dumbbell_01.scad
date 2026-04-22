// Dumbbell parameters
r_sphere = 5;      // Radius of the spheres
r_bar    = 3;      // Radius of the connecting bar
len_bar  = 15;     // Total length of the bar (from sphere to sphere)

// Calculate sphere positions based on the bar length
x_left  = -len_bar / 2;
x_right =  len_bar / 2;

// Create the dumbbell using a union
union() {
    // Left sphere
    translate([x_left, 0, 0]) sphere(r = r_sphere);
    
    // Right sphere
    translate([x_right, 0, 0]) sphere(r = r_sphere);
    
    // Connecting cylinder (bar)
    // Cylinder along X‑axis: rotate the default Z‑axis cylinder by 90° around Y
    rotate([0, 90, 0]) cylinder(h = len_bar, r = r_bar, center = true);
}

