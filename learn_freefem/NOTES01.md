Visualisasi objek setelah membuat objek `mesh` pada FreeFEM:
```ff++
Mesh = buildmesh(b1(meshSize) + b2(meshSize) + b3(meshSize) + b4(meshSize));
plot( Mesh, wait=true );
```

Alternatif eksport ke format VTK, yang dapat divisualisasi dengan ParaView,
atau PyVista dan juga banyak tools visualisasi lainnya.
```
load "iovtk"
savevtk("IMG_mesh1.vtk", Mesh);
```