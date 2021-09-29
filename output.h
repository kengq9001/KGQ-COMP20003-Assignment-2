// Compute and store diameters
void ComputeDiameters(dcel_t *dcel, wt_info_t **wts, int nWt);

// Prints all watchtower information into a file
void PrintWtsAndDiameters(FILE *file, dcel_t *dcel, wt_info_t **wts, int nWt);

// Sort watchtowers by diameter of cells using insertion sort
void InsertionSortByDiameters(wt_info_t **wts, int nWt);