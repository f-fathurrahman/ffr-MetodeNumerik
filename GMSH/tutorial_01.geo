// define new variable
lc = 1e-2;

// Use the variable in other construction
// Define a Point
// Point is uniquely identified by a tag (must be positive number)
// and defined by a list of four numbers: three coords, x, y, z and target
// mesh size (lc) close to the point
Point(1) = {0.0, 0.0, 0.0, lc};

// If no target mesh size of provided, a default uniform coarse size will be
// used for the model, based on the overall model size.

// We can then define some additional points. All points should have different
// tags:
Point(2) = {0.1, 0.0, 0.0, lc};
Point(3) = {0.5, 0.3, 0.0, lc};
Point(4) = {0.0, 0.3, 0.0, lc};

// Curves are Gmsh’s second type of elementary entities, and, amongst curves,
// straight lines are the simplest. A straight line is identified by a tag and
// is defined by a list of two point tags. In the commands below, for example,
// the line 1 starts at point 1 and ends at point 2.
//
// Note that curve tags are separate from point tags - hence we can reuse tag
// ‘1’ for our first curve. And as a general rule, elementary entity tags in
// Gmsh have to be unique per geometrical dimension

Line(1) = {1, 2};
Line(2) = {3, 2}; // why not {2,3} ?
Line(3) = {3, 4};
Line(4) = {4, 1};


// The third elementary entity is the surface. In order to define a simple
// rectangular surface from the four curves defined above, a curve loop has
// first to be defined. A curve loop is also identified by a tag (unique amongst
// curve loops) and defined by an ordered list of connected curves, a sign being
// associated with each curve (depending on the orientation of the curve to form
// a loop):
Curve Loop(1) = {4, 1, -2, 3};
// why using -2 instead of 2? See definition of Line(2)

// We can then define the surface as a list of curve loops (only one here,
// representing the external contour, since there are no holes--see ‘t4.geo’ for
// an example of a surface with a hole):
Plane Surface(1) = {1};

// At this level, Gmsh knows everything to display the rectangular surface 1 and
// to mesh it. An optional step is needed if we want to group elementary
// geometrical entities into more meaningful groups, e.g. to define some
// mathematical ("domain", "boundary"), functional ("left wing", "fuselage") or
// material ("steel", "carbon") properties.

// Such groups are called "Physical Groups" in Gmsh. By default, if physical
// groups are defined, Gmsh will export in output files only mesh elements that
// belong to at least one physical group. (To force Gmsh to save all elements,
// whether they belong to physical groups or not, set ‘Mesh.SaveAll=1;’, or
// specify ‘-save_all’ on the command line.) Physical groups are also identified
// by tags, i.e. strictly positive integers, that should be unique per dimension
// (0D, 1D, 2D or 3D). Physical groups can also be given names.
// Here we define a physical curve that groups the left, bottom and right curves
// in a single group (with prescribed tag 5); and a physical surface with name
// "My surface" (with an automatic tag) containing the geometrical surface 1:

// Here we define a physical curve that groups the left, bottom and right curves
// in a single group (with prescribed tag 5); and a physical surface with name
// "My surface" (with an automatic tag) containing the geometrical surface 1:
Physical Curve(5) = {1, 2, 4};
Physical Surface("My surface") = {1};

// Mesh 2;
// Save "TEMP_t1.msh";
// Save "TEMP_t1.m"; // Matlab script

