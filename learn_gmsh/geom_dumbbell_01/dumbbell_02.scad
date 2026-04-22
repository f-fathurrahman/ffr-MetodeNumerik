// Parameters
r_sphere = 5;      // sphere radius
r_bar    = 3;      // bar radius
len_bar  = 15;     // distance between sphere centers

x_left  = -len_bar / 2;
x_right =  len_bar / 2;

// Final dumbbell: spheres + (cylinder trimmed by spheres)
union() {
    // Left and right spheres
    translate([x_left, 0, 0]) sphere(r = r_sphere);
    translate([x_right, 0, 0]) sphere(r = r_sphere);
    
    // Cylinder with the parts inside the spheres cut away
    difference() {
        // Cylinder along X axis
        rotate([0, 90, 0]) 
            cylinder(h = len_bar, r = r_bar, center = true);
        
        // Remove the volume that lies inside the left sphere
        translate([x_left, 0, 0]) sphere(r = r_sphere);
        
        // Remove the volume that lies inside the right sphere
        translate([x_right, 0, 0]) sphere(r = r_sphere);
    }
}