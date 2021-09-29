// Ken Gene Quah
// COMP20003 Assignment 2

// main.c
// Passes the files to the necessary functions/operations
// and creates the desired output.

// Using parts of base code from
// https://edstem.org/au/courses/6413/discussion/576843
// https://edstem.org/au/courses/6413/workspaces/p8zK87PpJ43kbuGCtxLm51c1G6iBg5Tc

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "wt_ops.h"
#include "dcel_ops.h"
#include "output.h"

#define STAGE1 1 
#define STAGE2 2 
#define STAGE3 3 
#define STAGE4 4

#define MAX_DISTANCE 200

void ProcessArg(int argc, char **argv, int *stage, char **towerFName,
    char **polygonFName, char **outFName, char **pointsFName); 
void DoStage1(char *pointsFName, FILE *outFile);
void DoStage2(char *pointsFName, char *polygonFName, FILE *outFile);
void DoStage34(int stage, char *towerFName, char *polygonFName, FILE *outFile);

int main(int argc, char **argv) {
    char *towerFName, 
         *polygonFName,
         *outFName,
	     *pointsFName;
	int stage;
	ProcessArg(argc, argv, &stage, &towerFName, &polygonFName, &outFName, &pointsFName);
   
	FILE *outFile= fopen(outFName, "w");

	if (stage==STAGE1) {
		DoStage1(pointsFName, outFile);
	} else if (stage==STAGE2) {
		DoStage2(pointsFName, polygonFName, outFile);
	} else {
		DoStage34(stage, towerFName, polygonFName, outFile);
	}

	fclose(outFile);
    return 0;
}

//==============================================================================
// Anh Vo's Assignment 2 Base Code
// The following functions handle the stages
//==============================================================================

// Stage 1:
// Outputs the bisectors for the given file pointsFName	
void DoStage1(char *pointsFName, FILE *outFile) {
	FILE *inFile = fopen(pointsFName, "r");
    if (!inFile) {
        printf("Error: File not found");
        exit(EXIT_FAILURE);
    }

	// For each line, get and print the bisector equation
	double x1, y1, x2, y2;
	bisector_t bisector;
	while (fscanf(inFile, "%lf %lf %lf %lf", &x1, &y1, &x2, &y2) == 4) {
		GetBisector(&bisector, x1, y1, x2, y2);
		PrintBisector(outFile, &bisector);
	}
	fclose(inFile); 
}

// Stage 2:
// 1. Makes the first face from the file polygonFName
// 2. Outputs the intersections of the bisectors of points in pointsFName
//    with the polygon
void DoStage2(char *pointsFName, char *polygonFName, FILE *outFile) {
    FILE *inFile = fopen(polygonFName, "r");
    if (!inFile) {
        printf("Error: File not found");
        exit(EXIT_FAILURE);
    }
    dcel_t *DCEL;
    DCEL = CreateDcel();
    FirstPolygon(DCEL, inFile);
	fclose(inFile);

	inFile = fopen(pointsFName, "r");
    if (!inFile) {
        printf("Error: File not found");
        exit(EXIT_FAILURE);
    }

    double x1, y1, x2, y2;
	bisector_t bisector;
    split_info_t split;
	// Start from edge 0 because we want the output to be in order (already is)
	// DCEL->face_list[0]->hedge = DCEL->edge_list[0]->hedge;

	// For each line, get the bisector, then get the split information
	// between the initial face and the bisector.
	while (fscanf(inFile, "%lf %lf %lf %lf", &x1, &y1, &x2, &y2) == 4) {
        GetBisector(&bisector, x1, y1, x2, y2);
        GetSplitInfo(DCEL, 0, MAX_DISTANCE, &bisector, &split, 1, 1);
        PrintSplitInfo(outFile, &split);
	}
	fclose(inFile);
	FreeDcel(DCEL);
}

// Stage 3/4:
// 1. Makes the first face from the file polygonFName
// 2. Genarates voronopi cell for each tower from file twoerFName
// 3. (stage 4 only) sorts the tower by diameter of their faces
// 4. Outputs the towers
void DoStage34(int stage, char *towerFName, char *polygonFName, FILE *outFile) {
	FILE *inFile = fopen(polygonFName, "r");
    if (!inFile) {
        printf("Error: File not found");
        exit(EXIT_FAILURE);
    }
    dcel_t *DCEL;
    DCEL = CreateDcel();
    FirstPolygon(DCEL, inFile);
	fclose(inFile);

	inFile = fopen(towerFName, "r");
    if (!inFile) {
        printf("Error: File not found");
        exit(EXIT_FAILURE);
    }
	wt_info_t **watchtowers; int watchtowers_num=0;
	watchtowers = ReadWtInfo(inFile, &watchtowers_num);
	fclose(inFile);

	// Build the Voronoi diagram tower by tower
	int i;
	if (watchtowers_num) {
		DCEL->face_list[0]->wt = watchtowers[0];
		watchtowers[0]->face = 0;
	}
	for (i=1; i<watchtowers_num; i++) {
		IncrementalVoronoi(DCEL, watchtowers[i], MAX_DISTANCE);
	}
	
	ComputeDiameters(DCEL, watchtowers, watchtowers_num);
	if (stage==STAGE4) {
		InsertionSortByDiameters(watchtowers, watchtowers_num);
	}

	// Print output
	PrintWtsAndDiameters(outFile, DCEL, watchtowers, watchtowers_num);
	FILE *graphFile = fopen("graph.txt", "w");
	PrintDcel(graphFile, DCEL);
	fclose(graphFile);
	FreeWts(watchtowers, &watchtowers_num);
	FreeDcel(DCEL);
}

// Processing Arguments
void ProcessArg(int argc, char **argv, int *stage, char **towerFName,
	char **polygonFName, char **outFName, char **pointsFName) {
	int valid= argc >1 && atoi(argv[1])>=1 && atoi(argv[1])<=4;
	if (valid) {
		*stage= atoi(argv[1]);
		valid = (*stage==STAGE1 && argc == 4) || 
		        (*stage>=STAGE2 && *stage<=STAGE4 && argc==5);
	}
	if (!valid){
		fprintf(stderr, "Usage:\n");
		fprintf(stderr, "\t%s %d point_pairs_file output_file\n"
	                "\t%s %d point_pairs_file polygon_file output_file\n"
	                "\t%s %d csv_data_file polygon_file output_file\n"
	                "\t%s %d csv_data_file polygon_file output_file\n",
	                argv[0], STAGE1, argv[0], STAGE2, argv[0], STAGE3, argv[0],STAGE4);
		exit(EXIT_FAILURE);
	}
	*towerFName= *pointsFName= *polygonFName= NULL;
	if (*stage==STAGE1) {
		*pointsFName= argv[2];
		*outFName= argv[3];
	} else {
		*polygonFName= argv[3];
		*outFName= argv[4];
		if (*stage==2) 
			*pointsFName= argv[2];
		else
			*towerFName= argv[2];
	} 
}