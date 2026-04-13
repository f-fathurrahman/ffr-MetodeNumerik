SetFactory("OpenCASCADE");
// Create a cube of size 1x1x1
Box(1) = {0, 0, 0, 1, 1, 1};

// Define Periodic Boundaries
// Surface 1 (bottom), Surface 3 (top) - Transate along Z
Periodic Surface {3} = {1} Translate {0, 0, 1}; 

// Surface 2 (front), Surface 4 (back) - Translate along Y
Periodic Surface {4} = {2} Translate {0, 1, 0}; 

// Surface 5 (left), Surface 6 (right) - Translate along X
Periodic Surface {6} = {5} Translate {1, 0, 0}; 

// Define Mesh Size and Generate
MeshSize {:} = 0.1;
Mesh 3;

