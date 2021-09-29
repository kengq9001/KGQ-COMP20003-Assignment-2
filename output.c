// Functions used for processing output

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "wt_ops.h"
#include "dcel_ops.h"

// Compute and store diameters
void ComputeDiameters(dcel_t *dcel, wt_info_t **wts, int nWt) {
    for (int i = 0; i < nWt; i++) {
        if (wts[i]->face != NO_FACE) {
            wts[i]->diameter = GetDiameter(dcel, wts[i]->face);
        }
    }
}

// Prints all watchtower information into a file
void PrintWtsAndDiameters(FILE *file, dcel_t *dcel, wt_info_t **wts, int nWt) {
    for (int i = 0; i < nWt; i++) {
        PrintWtInfo(file, wts[i]);
    }
}

// Sort watchtowers by diameter of cells using insertion sort
void InsertionSortByDiameters(wt_info_t **wts, int nWt) {
    if (nWt <= 1) {
        return;
    }
    int i, j;
    wt_info_t *wt;
    for (i=1;i<nWt;i++) {
        wt = wts[i];
        j = i - 1;
        while ((j>=0) && (wts[j]->diameter > wt->diameter)) {
            wts[j+1] = wts[j];
            j--;
        }
        wts[j+1] = wt;
    }
}