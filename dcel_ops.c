// dcel_ops.c
// DCEL Operations

// Handles the following:
// - Creating a DCEL and its components
//       (vertices, faces, edges, half-edges i.e. hedges)
// - Constructing a polygon
// - Splitting a face with a line segment between two given points on the boundary
// - Checking if a point lies in a "clockwise" half-plane of two points
// - Checking if a point lies inside a face
// - Finding the perpendicular bisector of two points
// - Finding the diameter of a face
// - Incremental Voronoi algorithm
// - Finding the type of intersection of two line segments and point of intersection (if exists)
// - Printing bisector equations, splitting information and DCEL
// - Freeing a DCEL and its components

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <math.h>
#include "dcel_ops.h"

// Need to know about watchTowers to store them in the face
// and store face in them
#include "wt_ops.h"

#define V_START_SIZE 4
#define F_START_SIZE 4
#define E_START_SIZE 4
#define H_START_SIZE 4

// Labelling the outside face simplifies things significantly
// And it is also counted in the formula V + F = E + 2 for planar graphs.
#define EXTERIOR_FACE -1

#define NOT_SPLIT -1
#define EPS 1e-9
#define DELETED_VERTEX -1
#define DELETED_FACE -2
#define DELETED_EDGE -1

//==============================================================================
// The following 5 functions create a DCEL and its elements
//==============================================================================

// Returns pointer to a newly created DCEL
dcel_t *CreateDcel() {
    dcel_t *dcel;
    dcel = (dcel_t*)malloc(sizeof(dcel_t));
    assert(dcel);

    dcel->num_vertex = 0;
    dcel->size_vertex_list = V_START_SIZE;
    dcel->vertex_list = (vertex_t**)malloc(sizeof(vertex_t*)*V_START_SIZE);
    assert(dcel->vertex_list);

    dcel->num_face = 0;
    dcel->size_face_list = F_START_SIZE;
    dcel->face_list = (face_t**)malloc(sizeof(face_t*)*F_START_SIZE);
    assert(dcel->face_list);

    dcel->num_edge = 0;
    dcel->size_edge_list = E_START_SIZE;
    dcel->edge_list = (edge_t**)malloc(sizeof(edge_t*)*E_START_SIZE);
    assert(dcel->edge_list);

    return dcel;
}

// Adds a vertex to an existing DCEL
void AddVertex(dcel_t *dcel, double x, double y) {
    int *n = &dcel->num_vertex;
    int *size = &dcel->size_vertex_list;
    // reallocate space if not enough
    if (*n == *size) {
        *size*=2;
        dcel->vertex_list = (vertex_t**)
        realloc(dcel->vertex_list, sizeof(vertex_t*)*(*size));
        assert(dcel->vertex_list);
    }

    vertex_t *v;
    v = (vertex_t*)malloc(sizeof(vertex_t));
    assert(v);

    v->index = *n;
    v->x = x;
    v->y = y;
    dcel->vertex_list[*n] = v;
    *n+=1;
}

// Adds a face to an existing DCEL
void AddFace(dcel_t *dcel) {
    int *n = &dcel->num_face;
    int *size = &dcel->size_face_list;
    // reallocate space if not enough
    if (*n == *size) {
        *size*=2;
        dcel->face_list = (face_t**)
        realloc(dcel->face_list, sizeof(face_t*)*(*size));
        assert(dcel->face_list);
    }

    face_t *f;
    f = (face_t*)malloc(sizeof(face_t));
    assert(f);

    f->index = *n;
    f->hedge = NULL;
    dcel->face_list[*n] = f;
    *n+=1;
}

// Adds an edge to an existing DCEL
// Also adds 2 half-edges, v_s-->v_e having face f1 and v_e-->v_s having face f2
// This new edge is automatically linked to the first half edge.
// The two half-edges are automatically linked as twins
void AddEdge(dcel_t *dcel, int v_s, int v_e, int f1, int f2) {
    int *n = &dcel->num_edge;
    int *size = &dcel->size_edge_list;
    // reallocate space if not enough
    if (*n == *size) {
        *size*=2;
        dcel->edge_list = (edge_t**)
        realloc(dcel->edge_list, sizeof(edge_t*)*(*size));
        assert(dcel->edge_list);
    }

    edge_t *e;
    e = (edge_t*)malloc(sizeof(edge_t));
    assert(e);

    e->index = *n;
    e->hedge = AddHedge(v_s, v_e, f1, *n);
    e->hedge->twin = AddHedge(v_e, v_s, f2, *n);
    e->hedge->twin->twin = e->hedge;
    dcel->edge_list[*n] = e;
    *n+=1;
}

// Returns pointer to a newly created half edge
hedge_t *AddHedge(int v_s, int v_e, int f, int e) {
    hedge_t *hedge;
    hedge = (hedge_t*)malloc(sizeof(hedge_t));
    assert(hedge);

    hedge->v_start = v_s;
    hedge->v_end = v_e;
    hedge->face = f;
    hedge->edge = e;
    hedge->next = NULL;
    hedge->prev = NULL;
    hedge->twin = NULL;
    return hedge;
}

//==============================================================================
// Reads a file containing the vertices of the polygon
// and constructs a DCEL
//==============================================================================
void FirstPolygon(dcel_t *dcel, FILE *file) {
    // Add all vertices
    double X, Y;
    while (fscanf(file, "%lf %lf", &X, &Y) > 0) {
        AddVertex(dcel, X, Y);
    }

    // Add 1 face and all edges
    AddFace(dcel);
    int i;
    int n = dcel->num_vertex;
    for (i=0;i<n;i++) {
        AddEdge(dcel, i, (i+1)%n, 0, EXTERIOR_FACE);
    }

    // Link half edges, assuming that the polygon's vertices are oriented clockwise
    // %n as a neat trick to "wrap around" between indices 0 and n-1
    dcel->face_list[0]->hedge = dcel->edge_list[0]->hedge;
    for (i=0;i<n;i++) {
        dcel->edge_list[i]->hedge->next = dcel->edge_list[(i+1)%n]->hedge;
        dcel->edge_list[i]->hedge->prev = dcel->edge_list[(n+i-1)%n]->hedge;
        dcel->edge_list[i]->hedge->twin->next = dcel->edge_list[(n+i-1)%n]->hedge->twin;
        dcel->edge_list[i]->hedge->twin->prev = dcel->edge_list[(i+1)%n]->hedge->twin;
    }
}

//==============================================================================
// A special case for splitting, where both edges are bisected.
// Returns this splitting information
//==============================================================================
split_info_t TwoMidpoints(dcel_t *dcel, int e1, int e2) {
    split_info_t split;
    split.face = -1;
    hedge_t *h1 = dcel->edge_list[e1]->hedge;
    hedge_t *h2 = dcel->edge_list[e2]->hedge;

    // Temporary variables to help compute
    // the midpoints of AB and CD
    // (AB and CD could be oriented incorrectly, but that does not matter)
    int e1A = h1->v_start;
    int e1B = h1->v_end;
    double Ax = dcel->vertex_list[e1A]->x;
    double Ay = dcel->vertex_list[e1A]->y;
    double Bx = dcel->vertex_list[e1B]->x;
    double By = dcel->vertex_list[e1B]->y;
    int e2A = h2->v_start;
    int e2B = h2->v_end;
    double Cx = dcel->vertex_list[e2A]->x;
    double Cy = dcel->vertex_list[e2A]->y;
    double Dx = dcel->vertex_list[e2B]->x;
    double Dy = dcel->vertex_list[e2B]->y;

    split.edge_start = e1;
    split.x_start = (Ax + Bx) / 2;
    split.y_start = (Ay + By) / 2;
    split.edge_end = e2;
    split.x_end = (Cx + Dx) / 2;
    split.y_end = (Cy + Dy) / 2;
    split.v_start = NOT_SPLIT;
    split.v_end = NOT_SPLIT;

    return split;
}

// Reverse a split's direction
void ReverseSplit(split_info_t *split) {
    int temp_int;
    double temp_double;
    temp_int = split->v_start;
    split->v_start = split->v_end;
    split->v_end = temp_int;
    temp_int = split->edge_start;
    split->edge_start = split->edge_end;
    split->edge_end = temp_int;
    temp_double = split->x_start;
    split->x_start = split->x_end;
    split->x_end = temp_double;
    temp_double = split->y_start;
    split->y_start = split->y_end;
    split->y_end = temp_double;
}

//==============================================================================
// Splits a face with the provided splitting information.
// This adds a face and the edge which splits the face.
// For each vertex that doesn't already exist (0, 1 or 2), 
// a new vertex and edge will also be created.
// Also makes appropriate changes to the half edges.
//==============================================================================
void SplitFace(dcel_t *dcel, split_info_t *split) {
    int m_i, n_i, B_index, C_index;
    hedge_t *AM, *MB, *CN, *ND, *MN;
    hedge_t *old_next_hedge, *old_prev_hedge, *old_next_hedge_twin, *old_prev_hedge_twin;
    int start_exists = (split->v_start != NOT_SPLIT);
    int end_exists = (split->v_end != NOT_SPLIT);

    int e1 = split->edge_start;
    int e2 = split->edge_end;
    hedge_t *h1 = dcel->edge_list[e1]->hedge;
    hedge_t *h2 = dcel->edge_list[e2]->hedge;
    // First identify the face where the split is occuring.
    // Then make both edges point to the half-edge that is part of the face, 
    // rather than point to the half-edge that is not part of the face.
    // This ensures that the new edges are created in the correct position and order.
    if (split->face >= 0) {
        if (h1->face != split->face) {
            dcel->edge_list[e1]->hedge = dcel->edge_list[e1]->hedge->twin;
        }
        if (h2->face != split->face) {
            dcel->edge_list[e2]->hedge = dcel->edge_list[e2]->hedge->twin;
        }
    } else {
        // create a temporary test point P for testing face
        double Px = (split->x_start + split->x_end) / 2;
        double Py = (split->y_start + split->y_end) / 2;
        if (HalfPlane(dcel, h1->v_start, h1->v_end, Px, Py) == -1) {
            dcel->edge_list[e1]->hedge = dcel->edge_list[e1]->hedge->twin;
        }
        if (HalfPlane(dcel, h2->v_start, h2->v_end, Px, Py) == -1) {
            dcel->edge_list[e2]->hedge = dcel->edge_list[e2]->hedge->twin;
        }
    }
    // Now both edges point to half-edges which are in the face we are splitting

    // Now let first edge = A-->B, second edge = C-->D (clockwise order)
    // Add 2 vertices
    // Let this be M, midpoint of A-->B
    if (start_exists) {
        m_i = split->v_start;
    } else {
        m_i = dcel->num_vertex;
        AddVertex(dcel, split->x_start, split->y_start);
    }
    // Let this be N, midpoint of C-->D
    if (end_exists) {
        n_i = split->v_end;
    } else {
        n_i = dcel->num_vertex;
        AddVertex(dcel, split->x_end, split->y_end);
    }

    // Give names to edges
    AM = dcel->edge_list[e1]->hedge;
    ND = dcel->edge_list[e2]->hedge;
    if (start_exists) {
        if (AM->v_start == m_i) {
            AM = AM->prev;
        }
    }
    if (end_exists) {
        if (ND->v_end == n_i) {
            ND = ND->next;
        }
    }

    // face indices
    int face_old = AM->face;
    int face_out1 = AM->twin->face;
    int face_out2 = ND->twin->face;
    int face_new = dcel->num_face;
    AddFace(dcel);

    // Add 3 edges

    // 1st edge, add M-->N
    AddEdge(dcel, m_i, n_i, face_old, face_new);
    MN = dcel->edge_list[dcel->num_edge-1]->hedge;
    dcel->face_list[face_old]->hedge = MN;
    dcel->face_list[face_new]->hedge = MN->twin;

    // Pre-processing for 2nd edge, change A-->B into A-->M
    if (start_exists) {
        MB = AM->next;
    } else {
        B_index = AM->v_end;
        old_next_hedge = AM->next;
        old_prev_hedge_twin = AM->twin->prev;
        AM->v_end = m_i;
        AM->twin->v_start = m_i;
    }
    AM->next = MN;
    MN->prev = AM;

    // 2nd edge, add M-->B
    if (!start_exists) {
        AddEdge(dcel, m_i, B_index, face_new, face_out1);
        MB = dcel->edge_list[dcel->num_edge-1]->hedge;
        MB->next = old_next_hedge;
        old_next_hedge->prev = MB;
        MB->twin->prev = old_prev_hedge_twin;
        old_prev_hedge_twin->next = MB->twin;
        MB->twin->next = AM->twin;
        AM->twin->prev = MB->twin;
    }
    MB->prev = MN->twin;
    MN->twin->next = MB;

    // Pre-processing for 3rd edge, change C-->D into N-->D
    if (end_exists) {
        CN = ND->prev;
    } else {
        C_index = ND->v_start;
        old_prev_hedge = dcel->edge_list[e2]->hedge->prev;
        old_next_hedge_twin = dcel->edge_list[e2]->hedge->twin->next;
        ND->v_start = n_i;
        ND->twin->v_end = n_i;
    }
    ND->prev = MN;
    MN->next = ND;

    // 3rd edge, add C-->N
    if (!end_exists) {
        AddEdge(dcel, C_index, n_i, face_new, face_out2);
        CN = dcel->edge_list[dcel->num_edge-1]->hedge;
        CN->prev = old_prev_hedge;
        old_prev_hedge->next = CN;
        CN->twin->next = old_next_hedge_twin;
        old_next_hedge_twin->prev = CN->twin;
        CN->twin->prev = ND->twin;
        ND->twin->next = CN->twin;
    }
    CN->next = MN->twin;
    MN->twin->prev = CN;

    // update faces
    hedge_t *h = MB;
    while (h != MN->twin) {
        h->face = face_new;
        h = h->next;
    }
}

// Merges the second face into the first face, given that they share exactly
// one consecutive chain of edges, cleaning up unused geometry in the process.
void MergeFaces(dcel_t *dcel, int f1, int f2) {
    // We want to locate these 4 important half edges
    // they describe the half edges before/after the chain (going clockwise)
    hedge_t *f1_before, *f1_after, *f2_before, *f2_after;

    hedge_t *h = dcel->face_list[f1]->hedge;
    // Stops when current edge is not common and next edge is common
    while ((h->twin->face == f2) || (h->next->twin->face != f2)) {
        h = h->next;
    }
    f1_before = h;
    f2_after = f1_before->next->twin->next;
    // Stops when current edge is not common and previous edge is common
    while ((h->twin->face == f2) || (h->prev->twin->face != f2)) {
        h = h->next;
    }
    f1_after = h;
    f2_before = f1_after->prev->twin->prev;

    h = f1_before->next;
    while (h != f1_after) {
        // This vertex is redundant because it is in the interior of the chain
        if (h != f1_before->next) {
            dcel->vertex_list[h->v_start]->index = DELETED_VERTEX;
        }

        dcel->edge_list[h->edge]->index = DELETED_EDGE;
        h = h->next;
    }

    f1_before->next = f2_after;
    f2_after->prev = f1_before;
    f1_after->prev = f2_before;
    f2_before->next = f1_after;
    h = f2_after;
    while (h != f1_after) {
        h->face = f1;
        h = h->next;
    }
    dcel->face_list[f2]->index = DELETED_FACE;
}

// Returns 1 if P is in the clockwise half plane as the half edge v1-->v2,
// returns -1 if P is in the anticlockwise half plane instead,
// returns 0 if P lies on v1-->v2.
int HalfPlane(dcel_t *dcel, int v1, int v2, double Px, double Py) {
    vertex_t *A = dcel->vertex_list[v1];
    vertex_t *B = dcel->vertex_list[v2];
    return AreaSign(A->x, A->y, Px, Py, B->x, B->y);
}

// Returns -1, 0 or 1, based on the area enclosed by the three points. 0 corresponds
// to no area enclosed.
int AreaSign(double sx, double sy, double ex, double ey, double x, double y){
    double areaSq;
    /* |AB x AC|^2, squared area */
    /* See https://mathworld.wolfram.com/CrossProduct.html */
    areaSq = (ex - sx) * (y  - sy) -
             (x  - sx) * (ey - sy);
    
    if(areaSq > 0.0){
        return 1;
    } else if(areaSq == 0.0){
        return 0;
    } else {
        return -1;
    }
}

// Returns 1 if the given x,y point is inside or on the given face, 0 otherwise
int InFace(dcel_t *dcel, int face, double x, double y) {
    if (dcel->face_list[face]->index == DELETED_FACE) {
        return 0;
    }
    int passed = 0;
    hedge_t *hedge = dcel->face_list[face]->hedge;
    while (hedge != dcel->face_list[face]->hedge || !passed) {
        passed = 1;
        if ((HalfPlane(dcel, hedge->v_start, hedge->v_end, x, y)) == -1) {
            return 0;
        }
        hedge = hedge->next;
    }
    return 1;
}

// Obtains information which describes the perpendicular bisector of two points
void GetBisector(bisector_t *b, double x1, double y1, double x2, double y2) {
    b->M_x = (x1 + x2) / 2;
    b->M_y = (y1 + y2) / 2;
    if (y1 == y2) {
        b->notVertical = 0;
        return;
    }
    b->notVertical = 1;
    b->m = (x1 - x2)/(y2 - y1);
    if ((b->m) == 0) {
        b->m = (double)0;
    }
}

// Gets a point at least abs(distance) away from the midpoint of the bisector given.
void GetBisectorPoint(double distance, bisector_t *b, double *x, double *y) {
    if (!b->notVertical) {
        *x = b->M_x;
        *y = b->M_y + distance;
        return;
    }
    *x = b->M_x + distance;
    *y = b->m * (*x - b->M_x) + b->M_y;
}

// Obtain information which describes a split in a specified face
// (find_start, find_end) combinations:
// (1,1): First/Second intersection found becomes start/end of split respectively
// (1,0): First intersection found becomes start of split
// (0,1): First intersection found becomes end of split
void GetSplitInfo(dcel_t *dcel, int face, double distance, bisector_t *b, 
                    split_info_t *split, int find_start, int find_end) {
    double Px, Py;
    int points_found = (!find_start);
    hedge_t *hedge = dcel->face_list[face]->hedge;
    while (points_found < (2 - !find_end) ) {
        if ((intersects(dcel, b, hedge, distance, &Px, &Py) == INTERSECT)
         || (intersects(dcel, b, hedge, distance, &Px, &Py) == ENDS_OVERLAP)) {
            points_found++;
            if (points_found == 1) {
                split->v_start = NOT_SPLIT;
                split->edge_start = hedge->edge;
                split->x_start = Px;
                split->y_start = Py;
            } else {
                split->v_end = NOT_SPLIT;
                split->edge_end = hedge->edge;
                split->x_end = Px;
                split->y_end = Py;
            }
        }
        hedge = hedge->next;
    }
    split->face = face;
}

// Calculate the diameter of the given face
double GetDiameter(dcel_t *dcel, int face) {
    double diameter = 0;
    double dist, x1, y1, x2, y2;
    hedge_t *start = dcel->face_list[face]->hedge;
    hedge_t *end = start->prev;
    hedge_t *h1 = start, *h2 = start->next;

    while (h1 != end) {
        if (h2 == start) {
            h1 = h1->next;
            h2 = h1->next;
        }
        while (h2 != start) {
            x1 = dcel->vertex_list[h1->v_start]->x;
            y1 = dcel->vertex_list[h1->v_start]->y;
            x2 = dcel->vertex_list[h2->v_start]->x;
            y2 = dcel->vertex_list[h2->v_start]->y;

            dist = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
            if (dist > diameter) {
                diameter = dist;
            }
            h2 = h2->next;
        }
    }

    return diameter;
}

// Returns 1 if a is almost equal to b, 0 otherwise
int AlmostEqual(double a, double b) {
    return ((-EPS < a-b) && (a-b < EPS));
}

// Adds the watchtower to the Voronoi diagram represented by the given DCEL, 
// applying required splits and setting the watchtower as required.
void IncrementalVoronoi(dcel_t *dcel, wt_info_t *new_wt, double max_distance) {
    face_t *startFace, *currentFace;
    hedge_t *startHedge, *hedge;
    bisector_t bisector;
    split_info_t split;
    vertex_t *startVertex;
    int currentV_index;
    int newFace1, newFace2;
    int face = 0, looped = 0;
    double x, y;

    // First locate a face which the new watchtower belongs to
    // Assuming new watchtower is not in the exterior of the DCEL
    while (!InFace(dcel, face, new_wt->x, new_wt->y)) {
        face++;
        if (face == dcel->num_face) {
            printf("Error: New watchtower is not in the interior\n");
            exit(EXIT_FAILURE);
        }
    }
    startFace = dcel->face_list[face];

    // Find first split
    GetBisector(&bisector, startFace->wt->x, startFace->wt->y, new_wt->x, new_wt->y);
    GetSplitInfo(dcel, startFace->index, max_distance, &bisector, &split, 1, 1);

    // Make the new watchtower be in the new face
    if (AreaSign(split.x_start, split.y_start, new_wt->x, new_wt->y, 
        split.x_end, split.y_end) == 1) {
        ReverseSplit(&split);
    }
    
    SplitFace(dcel, &split);
    dcel->face_list[dcel->num_face-1]->wt = new_wt;
    new_wt->face = dcel->num_face-1;

    // Algorithm:
    // 0. Do first split
    // 1. Traverse the newest split, clockwise (anticlockwise if exterior face)
    // 2. When it meets an edge at point V1, go to its other side
    // 3.   2 cases:
    // 3a. Interior face: Find next split which goes from V1 to V2
    //     We know this is always true because V1 must be the intersection
    //     of 3 perpendicular bisectors, which are Voronoi edges.
    // 3b. Exterior face: Find the next place where we enter the interior
    //     Instead of being a split, it'll just be traversing the boundary
    // Merge faces and clean up the unused geometry as we go
    // 4. Repeat 1-3 until we get back to the original split
    startHedge = dcel->face_list[startFace->index]->hedge;
    startVertex = dcel->vertex_list[startHedge->v_start];
    hedge = startHedge->next->twin->next;
    newFace1 = startHedge->twin->face;
    newFace2 = newFace1;

    while (!looped) {        
        currentV_index = hedge->v_start;

        if (hedge->face != EXTERIOR_FACE) { // Interior face

            currentFace = dcel->face_list[hedge->face];
            // Want to begin intersection search past currentV_index
            dcel->face_list[hedge->face]->hedge = hedge->next;

            // Find bisector between new_wt and current face's wt, and split
            GetBisector(&bisector, currentFace->wt->x, currentFace->wt->y, new_wt->x, new_wt->y);
            GetSplitInfo(dcel, currentFace->index, max_distance, &bisector, &split, 0, 1);
            split.v_start = currentV_index;
            split.edge_start = hedge->edge;

            // If original starting point, done
            if (AlmostEqual(split.x_end, startVertex->x) &&
                AlmostEqual(split.y_end, startVertex->y)) {
                looped = 1;
                split.v_end = startVertex->index;
            }
            
            SplitFace(dcel, &split);
            hedge = dcel->face_list[currentFace->index]->hedge->next->twin->next;
            MergeFaces(dcel, newFace2, dcel->num_face-1);
        } else { // Exterior face

            if (currentV_index == startVertex->index) { // original starting point, done
                looped = 1;
            } else if (hedge->twin->face < new_wt->face) { // check if we re-enter the interior here
                // Check if the bisector (between new_wt and interior face's wt)
                // intersects current edge. Otherwise, move on
                currentFace = dcel->face_list[hedge->twin->face];
                GetBisector(&bisector, currentFace->wt->x, currentFace->wt->y, new_wt->x, new_wt->y);

                if ((intersects(dcel, &bisector, hedge, max_distance, &x, &y) == INTERSECT)
                || intersects(dcel, &bisector, hedge, max_distance, &x, &y) == ENDS_OVERLAP) {
                    
                    // found interior re-entry, same algorithm as main case with some changes
                    GetSplitInfo(dcel, currentFace->index, max_distance, &bisector, &split, 1, 1);
                    if (AreaSign(split.x_start, split.y_start, new_wt->x, new_wt->y, 
                        split.x_end, split.y_end) == 1) {
                        ReverseSplit(&split);
                    }

                    // If original starting point, done
                    if (AlmostEqual(split.x_end, startVertex->x) &&
                        AlmostEqual(split.y_end, startVertex->y)) {
                        looped = 1;
                        split.v_end = startVertex->index;
                        split.edge_end = startHedge->prev->edge;
                    }
                    
                    SplitFace(dcel, &split);
                    hedge = dcel->face_list[currentFace->index]->hedge->next->twin->next;
                    newFace2 = dcel->num_face-1;
                } else {
                    hedge = hedge->next;
                }

            } else { // this face is already part of the new face, so move on
                hedge = hedge->next;
            }

        }
    }

    if (newFace1 != newFace2) {
        MergeFaces(dcel, newFace1, newFace2);
    }
}

/* 
This intersection is based on code by Joseph O'Rourke and is provided for use in 
COMP20003 Assignment 2.

The approach for intersections is:
- Use the bisector to construct a finite segment and test it against the half-edge.
- Use O'Rourke's segseg intersection (https://hydra.smith.edu/~jorourke/books/ftp.html)
    to check if the values overlap.
*/

/* Returns 1 if the point (x, y) is in the line from s(x, y) to e(x, y), 0 otherwise. */
int collinear(double sx, double sy, double ex, double ey, double x, double y){
    /* Work out area of parallelogram - if it's 0, points are in the same line. */
    if (AreaSign(sx, sy, ex, ey, x, y) == 0){
        return 1;
    } else {
        return 0;
    }
}

/* Returns 1 if point (x, y) is between (sx, sy) and (se, se) */
int between(double sx, double sy, double ex, double ey, double x, double y){
    if(sx != ex){
        /* If not vertical, check whether between x. */
        if((sx <= x && x <= ex) || (sx >= x && x >= ex)){
            return 1;
        } else {
            return 0;
        }
    } else {
        /* Vertical, so can't check _between_ x-values. Check y-axis. */
        if((sy <= y && y <= ey) || (sy >= y && y >= ey)){
            return 1;
        } else {
            return 0;
        }
    }
}

/* 
    Generates a segment with each end at least minLength away in each direction 
    from the bisector midpoint. Returns 1 if b intersects the given half-edge
    on this segment, 0 otherwise. Sets the intersection point to the given x, y
    positions.
*/
enum intersectType intersects(dcel_t *dcel, bisector_t *b, hedge_t *he,
    double minLength, double *x, double *y){
    /* Half-edge x, y pair */
    double heSx = dcel->vertex_list[he->v_start]->x;
    double heSy = dcel->vertex_list[he->v_start]->y;
    double heEx = dcel->vertex_list[he->v_end]->x;
    double heEy = dcel->vertex_list[he->v_end]->y;
    
    /* Bisector x, y pair */
    double bSx;
    double bSy;
    double bEx;
    double bEy;
    GetBisectorPoint(-minLength, b, &bSx, &bSy);
    GetBisectorPoint(minLength, b, &bEx, &bEy);
    
    /* Parametric equation parameters */
    double t1;
    double t2;
    /* Numerators for X and Y coordinate of intersection. */
    double numeratorX;
    double numeratorY;
    /* Denominators of intersection coordinates. */
    double denominator;
    
    /*
    See http://www.cs.jhu.edu/~misha/Spring20/15.pdf
    for explanation and intuition of the algorithm here.
    x_1 = heSx, y_1 = heSy    |    p_1 = heS
    x_2 = heEx, y_2 = heEy    |    q_1 = heE
    x_3 = bSx , y_3 = bSy     |    p_2 =  bS
    x_4 = bEx , y_4 = bEy     |    q_2 =  bE
    ----------------------------------------
    So the parameters t1 and t2 are given by:
    | t1 |   | heEx - heSx  bSx - bEx | -1  | bSx - heSx |
    |    | = |                        |     |            |
    | t2 |   | heEy - heSy  bSy - bEy |     | bSy - heSy |
    
    Hence:
    | t1 |       1     | bSy - bEy        bEx - bSx |  | bSx - heSx |
    |    | = --------- |                            |  |            |
    | t2 |    ad - bc  | heSy - heEy    heEx - heSx |  | bSy - heSy |
    
        where 
        a = heEx - heSx
        b = bSx  -  bEx
        c = heEy - heSy
        d = bSy  -  bEy
    */
    
    /* Here we calculate ad - bc */
    denominator = heSx * (bEy  -  bSy) +
                  heEx * (bSy  -  bEy) +
                  bEx  * (heEy - heSy) +
                  bSx  * (heSy - heEy);
    
    if(denominator == 0){
        /* In this case the two are parallel */
        return parallelIntersects(heSx, heSy, heEx, heEy, bSx, bSy, bEx, bEy, x, y);
    }
    
    /*
    Here we calculate the top row.
    | bSy - bEy        bEx - bSx |  | bSx - heSx |
    |                            |  |            |
    |                            |  | bSy - heSy |
    */
    numeratorX = heSx * (bEy  -  bSy) +
                 bSx  * (heSy -  bEy) +
                 bEx  * (bSy  - heSy);
    
    /*
    Here we calculate the bottom row.
    |                            |  | bSx - heSx |
    |                            |  |            |
    | heSy - heEy    heEx - heSx |  | bSy - heSy |
    */
    numeratorY = -(heSx * (bSy  -  heEy) +
                   heEx * (heSy -  bSy) +
                   bSx  * (heEy  - heSy));
    
    /* Use parameters to convert to the intersection point */
    t1 = numeratorX/denominator;
    t2 = numeratorY/denominator;
    *x = heSx + t1 * (heEx - heSx);
    *y = heSy + t1 * (heEy - heSy);
    
    /* Make final decision - if point is on segments, parameter values will be
    between 0, the start of the line segment, and 1, the end of the line segment.
    */
    if (0.0 < t1 && t1 < 1.0 && 0.0 < t2 && t2 < 1.0){
        return INTERSECT;
    } else if(t1 < 0.0 || 1.0 < t1 || t2 < 0.0 || 1.0 < t2){
        /* s or t outside of line segment. */
        return DOESNT_INTERSECT;
    } else {
        /* 
        ((numeratorX == 0) || (numeratorY == 0) || 
         (numeratorX == denominator) || (numeratorY == denominator))
        */
        return ENDS_OVERLAP;
    }
}

enum intersectType parallelIntersects(double heSx, double heSy, double heEx, double heEy,
    double bSx, double bSy, double bEx, double bEy, double *x, double *y){
    if(!collinear(heSx, heSy, heEx, heEy, bSx, bSy)){
        /* Parallel, no intersection so don't set (x, y) */
        return DOESNT_INTERSECT;
    }
    /* bS between heS and heE */
    if(between(heSx, heSy, heEx, heEy, bSx, bSy)){
        *x = bSx; 
        *y = bSy;
        return SAME_LINE_OVERLAP;
    }
    /* bE between heS and heE */
    if(between(heSx, heSy, heEx, heEy, bEx, bEy)){
        *x = bEx;
        *y = bEy;
        return SAME_LINE_OVERLAP;
    }
    /* heS between bS and bE */
    if(between(bSx, bSy, bEx, bEy, heSx, heSy)){
        *x = heSx;
        *y = heSy;
        return SAME_LINE_OVERLAP;
    }
    /* heE between bS and bE */
    if(between(bSx, bSy, bEx, bEy, heEx, heEy)){
        *x = heEx; 
        *y = heEy;
        return SAME_LINE_OVERLAP;
    }
    
    return DOESNT_INTERSECT;
}

//==============================================================================
// Printing and Freeing
//==============================================================================

// Print the bisector equation into a file
void PrintBisector(FILE *file, bisector_t *b) {
    if (!b->notVertical) {
        fprintf(file, "x = %f\n", b->M_x);
        return;
    }
    fprintf(file, "y = %f * (x - %f) + %f\n", b->m, b->M_x, b->M_y);
}

// Print a information about a split into a file
void PrintSplitInfo(FILE *file, split_info_t *s) {
    fprintf(file, "From Edge %d (%f, %f) to Edge %d (%f, %f)\n", s->edge_start, 
            s->x_start, s->y_start, s->edge_end, s->x_end, s->y_end);
}

// Prints out DCEL for debugging
void PrintDcel(FILE *file, dcel_t *dcel) {
    int i, count, passed;
    vertex_t *v;
    face_t *f;
    hedge_t *h;

    count = 0;
    for (i=0; i<dcel->num_vertex; i++) {
        if (dcel->vertex_list[i] && dcel->vertex_list[i]->index != DELETED_VERTEX) {
            count++;
        }
    }
    fprintf(file, "Number of vertices: %d + %d deleted\n", count, dcel->num_vertex-count);

    count = 0;
    for (i=0; i<dcel->num_face; i++) {
        if (dcel->face_list[i] && dcel->face_list[i]->index != DELETED_FACE) {
            count++;
        }
    }
    fprintf(file, "Number of faces: %d + %d deleted\n", count, dcel->num_face-count);

    count = 0;
    for (i=0; i<dcel->num_edge; i++) {
        if (dcel->edge_list[i] && dcel->edge_list[i]->index != DELETED_EDGE) {
            count++;
        }
    }
    fprintf(file, "Number of edges: %d + %d deleted\n", count, dcel->num_edge-count);
    for (i=0; i<dcel->num_edge; i++) {
        PrintEdge(file, dcel, i);
    }

    fprintf(file, "Edge cycles:\n");
    for (i=0;i<dcel->num_face;i++) {
        f = dcel->face_list[i];
        if (!f || f->index == DELETED_FACE) {
            continue;
        }
        passed = 0;
        h = f->hedge;
        fprintf(file, "Face %d: %d", f->index, h->edge);
        while (h != f->hedge || !passed) {
            h = h->next;
            passed=1;
            fprintf(file, "->%d", h->edge);
        }
        fprintf(file, "\n");
    }

    for (i=0;i<dcel->num_vertex;i++) {
        v = dcel->vertex_list[i];
        if (!v || v->index == DELETED_VERTEX) {
            continue;
        }
        fprintf(file, "P_{%d}=(%f,%f)\n", v->index, v->x, v->y);
    }
    fprintf(file, "Vertex cycles:\n");
    for (i=0;i<dcel->num_face;i++) {
        f = dcel->face_list[i];
        if (!f || f->index == DELETED_FACE) {
            continue;
        }
        h = f->hedge;
        fprintf(file, "W_{%d}=(%f,%f)\n", f->index, f->wt->x, f->wt->y);
        fprintf(file, "\\operatorname{polygon}(P_{%d}", h->v_start);
        while (h != f->hedge->prev) {
            h = h->next;
            fprintf(file, ",P_{%d}", h->v_start);
        }
        fprintf(file, ")\n");
    }
}

void PrintEdge(FILE *file, dcel_t *dcel, int e) {
    fprintf(file, "Edge %d:\n", e);
    if (dcel->edge_list[e] && dcel->edge_list[e]->index == DELETED_EDGE) {
        fprintf(file, "\tDeleted\n");
        return;
    }
    hedge_t *h1 = dcel->edge_list[e]->hedge;
    hedge_t *h2 = h1->twin;
    fprintf(file, "\tVertices: %d->%d, Face: %d, Edge: %d\n", h1->v_start, h1->v_end, h1->face, h1->edge);
    fprintf(file, "\tVertices: %d->%d, Face: %d, Edge: %d\n", h2->v_start, h2->v_end, h2->face, h2->edge);
}

// Free the DCEL and its components
void FreeDcel(dcel_t *dcel) {
    int i;

    for (i=0;i<dcel->num_vertex;i++) {
        FreeVertex(dcel->vertex_list[i]);
    }
    free(dcel->vertex_list);

    for (i=0;i<dcel->num_face;i++) {
        FreeFace(dcel->face_list[i]);
    }
    free(dcel->face_list);

    for (i=0;i<dcel->num_edge;i++) {
        FreeEdge(dcel->edge_list[i]);
    }
    free(dcel->edge_list);

    free(dcel);
}

void FreeVertex(vertex_t *vertex) {
    if (vertex) free(vertex);
}

void FreeFace(face_t *face) {
    if (face) free(face);
}

// Assuming that the half edges only get freed when the edge is freed
void FreeEdge(edge_t *edge) {
    if (edge) {
        free(edge->hedge->twin);
        free(edge->hedge);
        free(edge);
    }
}
