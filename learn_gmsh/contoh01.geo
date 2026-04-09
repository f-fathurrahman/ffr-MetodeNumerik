// -----------------------------------------------------------------------------
//
//  Gmsh GEO tutorial 1
//
// Geometry basics, elementary entities, physical groups
// - Geometri dasar
// - entitas elementer (elementary entities),
// - physical groups
//
// -----------------------------------------------------------------------------

// Istilah yang digunakan adalah `affectation`
// Mungkin ini sama dengan pendefinisian variabel baru.
lc = 1e-2;

// Variabel ini kemudian digunakan untuk mendefinisikan entitas elementer
// paling sederhana yang ada pada Gmsh, yaitu titik atau `Point`.
// Definisi `Point` memerlukan: tag (integer positif), tiga angka (X, Y, Z)
// untuk koordinat dan target ukuran mesh di dekat titik ini. 
Point(1) = {0, 0, 0, lc};
// Distribusi ukuran mesh kemudian dapat diperoleh dengan interpolasi ukuran
// mesh tersebut diseluruh geometri.
// Ada beberapa cara lain yang dapat digunakan untuk spesifikasi ukuran mesh.
// Jika ukura mesh tidak diberikan maka ukuran default yang akan digunakan.

// Definisikan beberapa titik lain yang harus memiliki tag yang berbeda.
Point(2) = {0.1, 0.0, 0, lc};
Point(3) = {0.1, 0.3, 0, lc};
Point(4) = {0.0, 0.3, 0, lc};
// Apakah tag harus berurutan?


// Entitas geometri Gmsh yang lain adalah `Curve`. Ada beberapa jenis `Curve`
// yang paling sederhana dalah garis lurus `Line`. Suatu line didefinisikan
// dengan suatu tag, serta
// dua `Point` yang dirujuk dengan menggunakan tag dari `Point` tersebut.
Line(1) = {1, 2};
Line(2) = {3, 2};
Line(3) = {3, 4};
Line(4) = {4, 1};
// Tag untuk `Line` terpisah dari tag untuk `Point` maupun entitas geometri
// yang lain.
// Note that curve tags are separate from point tags - hence we can reuse tag
// `1' for our first curve. And as a general rule, elementary entity tags in
// Gmsh have to be unique per geometrical dimension.
// ffr: what is geometrical dimension ?


// The third elementary entity is the surface. In order to define a simple
// rectangular surface from the four curves defined above, a curve loop has
// first to be defined. A curve loop is also identified by a tag (unique amongst
// curve loops) and defined by an ordered list of connected curves, a sign being
// associated with each curve (depending on the orientation of the curve to form
// a loop):
// Entitas geometri Gmsh yang lain adalah `surface`. Sebagai contoh, kita
// akan mendefinisikan permukanaan persegi yang didefinisikan dengan suatu
// kurva tertutup atau `Loop`
Curve Loop(1) = {4, 1, -2, 3};

// We can then define the surface as a list of curve loops (only one here,
// representing the external contour, since there are no holes--see `t4.geo' for
// an example of a surface with a hole):
Plane Surface(1) = {1};

// At this level, Gmsh knows everything to display the rectangular surface 1 and
// to mesh it. An optional step is needed if we want to group elementary
// geometrical entities into more meaningful groups, e.g. to define some
// mathematical ("domain", "boundary"), functional ("left wing", "fuselage") or
// material ("steel", "carbon") properties.
//
// Such groups are called "Physical Groups" in Gmsh. By default, if physical
// groups are defined, Gmsh will export in output files only mesh elements that
// belong to at least one physical group. (To force Gmsh to save all elements,
// whether they belong to physical groups or not, set `Mesh.SaveAll=1;', or
// specify `-save_all' on the command line.) Physical groups are also identified
// by tags, i.e. strictly positive integers, that should be unique per dimension
// (0D, 1D, 2D or 3D). Physical groups can also be given names.
//
// Here we define a physical curve that groups the left, bottom and right curves
// in a single group (with prescribed tag 5); and a physical surface with name
// "My surface" (with an automatic tag) containing the geometrical surface 1:

Physical Curve(5) = {1, 2, 4};
Physical Surface("My surface") = {1};

// Now that the geometry is complete, you can
// - either open this file with Gmsh and select `2D' in the `Mesh' module to
//   create a mesh; then select `Save' to save it to disk in the default format
//   (or use `File->Export' to export in other formats);
// - or run `gmsh t1.geo -2` to mesh in batch mode on the command line.

// You could also uncomment the following lines in this script:
//
//   Mesh 2;
//   Save "t1.msh";
//
// which would lead Gmsh to mesh and save the mesh every time the file is
// parsed. (To simply parse the file from the command line, you can use `gmsh
// t1.geo -')

// By default, Gmsh saves meshes in the latest version of the Gmsh mesh file
// format (the `MSH' format). You can save meshes in other mesh formats by
// specifying a filename with a different extension in the GUI, on the command
// line or in scripts. For example
//
//   Save "t1.unv";
//
// will save the mesh in the UNV format. You can also save the mesh in older
// versions of the MSH format:
//
// - In the GUI: open `File->Export', enter your `filename.msh' and then pick
//   the version in the dropdown menu.
// - On the command line: use the `-format' option (e.g. `gmsh file.geo -format
//   msh2 -2').
// - In a `.geo' script: add `Mesh.MshFileVersion = x.y;' for any version
//   number `x.y'.
// - As an alternative method, you can also not specify the format explicitly,
//   and just choose a filename with the `.msh2' or `.msh4' extension.

// Note that starting with Gmsh 3.0, models can be built using other geometry
// kernels than the default built-in kernel. By specifying
//
//   SetFactory("OpenCASCADE");
//
// any subsequent command in the `.geo' file would be handled by the OpenCASCADE
// geometry kernel instead of the built-in kernel. Different geometry kernels
// have different features. With OpenCASCADE, instead of defining the surface by
// successively defining 4 points, 4 curves and 1 curve loop, one can define the
// rectangular surface directly with
//
//   Rectangle(2) = {.2, 0, 0, .1, .3};
//
// The underlying curves and points could be accessed with the `Boundary' or
// `CombinedBoundary' operators.
//
// See e.g. `t16.geo', `t18.geo', `t19.geo' or `t20.geo' for complete examples
// based on OpenCASCADE, and `examples/boolean' for more.