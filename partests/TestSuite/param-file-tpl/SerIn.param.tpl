#=============================================================
#    DROPS parameter file for    TestMGSerPar
#    Testroutines for serialization
#=============================================================

Brick {
  BasicRefX   = 4                   # number of basic refines in x direction
  BasicRefY   = 4                   # number of basic refines in y direction
  BasicRefZ   = 4                   # number of basic refines in z direction
  dim         = 1 1 1               # dx, dy, dz
  orig        = 0 0 0               # origin of the brick
}

Refining {
  MarkAll     = __REFINE__
  MarkDrop    = 1
  MarkCorner  = 1
  MarkingProc = -1
}

MGSerialization {
  Serialization = 1
  Overwrite     = 1
  SerDir        = ./geometry/
}

VTK{
    VTKOut      = 0
    VTKDir      = vtk
    VTKName     = Serialization
    Binary      = 0
}

Misc{
  Mode        = 1                   # mode=0: write serialized multigrid into a file, mode=1: read serialized multigrid
  Unknowns    = 1                   # test with or without unknowns
}
