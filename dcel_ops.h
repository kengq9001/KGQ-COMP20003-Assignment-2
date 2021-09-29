#ifndef _DCEL_OPS_H_
#define _DCEL_OPS_H_
/* Need to know about watchTowers to store them in the face. */
#include "wt_ops.h"

// hedge meaning half-edge
typedef struct hedge hedge_t;

// The following 5 structs are for the DCEL components and the DCEL
struct hedge {
    int v_start;
    int v_end;
    int face;
    int edge;
    hedge_t *next;
    hedge_t *prev;
    hedge_t *twin;
};

typedef struct {
    int index;
    double x;
    double y;
} vertex_t;

typedef struct {
    int index;
    hedge_t *hedge;
    wt_info_t *wt;
} face_t;

typedef struct {
    int index;
    hedge_t *hedge;
} edge_t;

typedef struct {
    int num_vertex;
    int size_vertex_list;
    vertex_t **vertex_list;
    int num_face;
    int size_face_list;
    face_t **face_list;
    int num_edge;
    int size_edge_list;
    edge_t **edge_list;
} dcel_t;

// Splitting information struct
// If a vertex already exists, it will split using that vertex.
// Otherwise, use the new information to create a new vertex.
typedef struct {
    int face; // if we are confident which face the split occurs in, otherwise set to -1
    int v_start;
    int edge_start;
    double x_start;
    double y_start;
    int v_end;
    int edge_end;
    double x_end;
    double y_end;
} split_info_t;

// (Perpendicular) Bisector information struct
// The midpoint has coordinates (M_x, M_y)
// If the line is not vertical, then the equation is y = m(x - M_x) + M_y
// Otherwise, the equation is x = M_x.
typedef struct {
    int notVertical;
    double M_x;
    double M_y;
    double m;    
} bisector_t;

typedef struct {
    double x;
    double y;
} intersection_t;

enum intersectType;

enum intersectType {
    DOESNT_INTERSECT  = 0, // Doesn't intersect
    INTERSECT         = 1, // Intersects
    SAME_LINE_OVERLAP = 2, // Lines are the same
    ENDS_OVERLAP      = 3  // Intersects at exactly one point (endpoint)
};

// The following 5 functions create a DCEL and its elements
dcel_t *CreateDcel();

void AddVertex(dcel_t *dcel, double x, double y);

void AddFace(dcel_t *dcel);

void AddEdge(dcel_t *dcel, int v_s, int v_e, int f1, int f2);

hedge_t *AddHedge(int v_s, int v_e, int f, int e);

// Constructs a DCEL that describes one polygon
void FirstPolygon(dcel_t *dcel, FILE *file);

// Returns splitting information which describes two edges being bisected
split_info_t TwoMidpoints(dcel_t *dcel, int e1, int e2);

// Reverse a split's direction
void ReverseSplit(split_info_t *split);

// Splits a face with the provided splitting information
void SplitFace(dcel_t *dcel, split_info_t *split);

// Merges the second face into the first face, given that they share exactly
// one consecutive chain of edges, cleaning up unused geometry in the process.
void MergeFaces(dcel_t *dcel, int f1, int f2);

// Returns 1, 0, -1 depending on where P is in relation to v1-->v2
int HalfPlane(dcel_t *dcel, int v1, int v2, double P_x, double P_y);

// Returns -1, 0 or 1, based on the area enclosed by the three points. 0 corresponds
// to no area enclosed.
int AreaSign(double sx, double sy, double ex, double ey, double x, double y);

// Returns 1 if the given x,y point is inside or on the given face, 0 otherwise
int InFace(dcel_t *dcel, int face, double x, double y);

// Obtains information which describes the perpendicular bisector of two points
void GetBisector(bisector_t *b, double x1, double y1, double x2, double y2);

// Gets a point at least abs(distance) away from the midpoint of the bisector given.
void GetBisectorPoint(double distance, bisector_t *b, double *x, double *y);

// Obtain information which describes a split in a specified face
void GetSplitInfo(dcel_t *dcel, int face, double distance, bisector_t *b, 
                    split_info_t *split, int find_start, int find_end);

// Calculate the diameter of the given face
double GetDiameter(dcel_t *dcel, int face);

// Returns 1 if a is almost equal to b, 0 otherwise
int AlmostEqual(double a, double b);

// Adds the watchtower to the Voronoi diagram represented by the given DCEL, 
// applying required splits and setting the watchtower as required.
void IncrementalVoronoi(dcel_t *dcel, wt_info_t *wt, double max_distance);

/* 
This intersection is based on code by Joseph O'Rourke and is provided for use in 
COMP20003 Assignment 2.

The approach for intersections is:
- Use the bisector to construct a finite segment and test it against the half-edge.
- Use O'Rourke's segseg intersection (https://hydra.smith.edu/~jorourke/books/ftp.html)
    to check if the values overlap.
*/

/* Returns 1 if the point (x, y) is in the line from s(x, y) to e(x, y), 0 otherwise. */
int collinear(double sx, double sy, double ex, double ey, double x, double y);

/* Returns 1 if point (x, y) is between (sx, sy) and (se, se) */
int between(double sx, double sy, double ex, double ey, double x, double y);

/* 
    Generates a segment with each end at least minLength away in each direction 
    from the bisector midpoint. Returns 1 if b intersects the given half-edge
    on this segment, 0 otherwise. Sets the intersection point to the given x, y
    positions.
*/
enum intersectType intersects(dcel_t *dcel, bisector_t *b, hedge_t *he, 
    double minLength, double *x, double *y);

enum intersectType parallelIntersects(double heSx, double heSy, double heEx, double heEy,
    double bSx, double bSy, double bEx, double bEy, double *x, double *y);

// Print the bisector equation into a file
void PrintBisector(FILE *file, bisector_t *b);

// Print a information about a split into a file
void PrintSplitInfo(FILE *file, split_info_t *split);

// Prints out DCEL for debugging
void PrintDcel(FILE *file, dcel_t *dcel);
void PrintEdge(FILE *file, dcel_t *dcel, int e);

// Free the DCEL and its components
void FreeDcel(dcel_t *dcel);
void FreeVertex(vertex_t *vertex);
void FreeFace(face_t *face);
void FreeEdge(edge_t *edge);
#endif