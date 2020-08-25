// Copyright (c) Lawrence Livermore National Security, LLC and other VisIt
// Project developers.  See the top-level LICENSE file for dates and other
// details.  No copyright assignment is required to contribute to VisIt.

#ifndef SIMV2_CURVILINEARMESH_H
#define SIMV2_CURVILINEARMESH_H
#include "export.h"
#include "VisItDataTypes.h"

// C-callable implementation of front end functions
#ifdef __cplusplus
extern "C"
{
#endif

V_LIBSIMXPORT int simv2_CurvilinearMesh_alloc(visit_handle*);
V_LIBSIMXPORT int simv2_CurvilinearMesh_free(visit_handle);

V_LIBSIMXPORT int simv2_CurvilinearMesh_setCoordsXY(visit_handle obj, int dims[2], 
                                                visit_handle x, visit_handle y);
V_LIBSIMXPORT int simv2_CurvilinearMesh_setCoordsXYZ(visit_handle obj, int dims[3], 
                                                 visit_handle x, visit_handle y, 
                                                 visit_handle z);
V_LIBSIMXPORT int simv2_CurvilinearMesh_setCoords2(visit_handle obj, int dims[2],
                                               visit_handle c);
V_LIBSIMXPORT int simv2_CurvilinearMesh_setCoords3(visit_handle obj, int dims[3],
                                               visit_handle c);
V_LIBSIMXPORT int simv2_CurvilinearMesh_setRealIndices(visit_handle obj, int min[3], int max[3]);
V_LIBSIMXPORT int simv2_CurvilinearMesh_setBaseIndex(visit_handle obj, int base_index[3]);
V_LIBSIMXPORT int simv2_CurvilinearMesh_setGhostCells(visit_handle h, visit_handle gz);
V_LIBSIMXPORT int simv2_CurvilinearMesh_setGhostNodes(visit_handle h, visit_handle gn);

V_LIBSIMXPORT int simv2_CurvilinearMesh_getCoords(visit_handle h, int *ndims, int dims[3],
                                              int *coordMode, 
                                              visit_handle *x,
                                              visit_handle *y,
                                              visit_handle *z,
                                              visit_handle *c);
V_LIBSIMXPORT int simv2_CurvilinearMesh_getRealIndices(visit_handle h, int min[3], int max[3]);
V_LIBSIMXPORT int simv2_CurvilinearMesh_getBaseIndex(visit_handle h, int base_index[3]);
V_LIBSIMXPORT int simv2_CurvilinearMesh_getGhostCells(visit_handle h, visit_handle *gz);
V_LIBSIMXPORT int simv2_CurvilinearMesh_getGhostNodes(visit_handle h, visit_handle *gn);

#ifdef __cplusplus
}
#endif


// Callable from within the runtime and SimV2
V_LIBSIMXPORT int simv2_CurvilinearMesh_check(visit_handle h);

#endif
