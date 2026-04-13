// Define points of a non-orthogonal parallelepiped
p1 = 1; Point(p1) = {0, 0, 0, l};
p2 = 2; Point(p2) = {1, 0.2, 0, l}; // Non-orthogonal skew
p3 = 3; Point(p3) = {1.2, 1.2, 0.1, l};
p4 = 4; Point(p4) = {0.2, 1, 0.1, l};

// Add top points (p5-p8) similarly...
// Create lines and surface loops to form a Volume
