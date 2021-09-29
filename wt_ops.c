// wt_ops.c
// Reads, prints, frees watchtower information

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include "wt_ops.h"

#define WT_START_SIZE 10

// Reads watchtower information from the csv file and returns an array
// of wt_info_t pointers.
// Also updates the number of watchtowers as the file is read.
// realloc usage idea from https://people.eng.unimelb.edu.au/ammoffat/ppsaa/c/realloc.c
wt_info_t **ReadWtInfo(FILE *file, int *n) {
    wt_info_t **watchtowers;
    int currentSize = WT_START_SIZE;
    watchtowers = (wt_info_t**)malloc(sizeof(wt_info_t*)*currentSize);
    assert(watchtowers);

    char *line = NULL;
    size_t lineBufferLength = 0;

    // first row which is just headings
    getline(&line, &lineBufferLength, file);

    // remaining rows
    while (getline(&line, &lineBufferLength, file) > 0) {
        // reallocate space if not enough
        if (*n == currentSize) {
            currentSize*=2;
            watchtowers = (wt_info_t**)
            realloc(watchtowers, sizeof(wt_info_t*)*currentSize);
            assert(watchtowers);
        }
        watchtowers[*n] = SplitLine(line);
        *n+=1;
    }
    
    free(line);
    return watchtowers;
}

// Split each line into 6 components and store each in the struct pointer
// Assuming the headings and their order is always the same
// strtok usage idea from https://www.tutorialspoint.com/c_standard_library/c_function_strtok.htm
wt_info_t *SplitLine(char *line) {
    wt_info_t *info;
    info = (wt_info_t*)malloc(sizeof(wt_info_t));
    assert(info);
    char *token;

    token = strtok(line, ",");
    info->ID = OneString(token);

    token = strtok(NULL, ",");
    info->postcode = OneString(token);

    token = strtok(NULL, ",");
    info->population = atoi(token);

    token = strtok(NULL, ",");
    info->contact_name = OneString(token);

    token = strtok(NULL, ",");
    info->x = atof(token);

    token = strtok(NULL, "\n");
    info->y = atof(token);

    info->face = NO_FACE;
    info->diameter = NO_DIAMETER;
    
    return info;
}

// Further string processing for string entries
char *OneString(char *line) {
    int len;
    len = strlen(line);
    
    char *temp;
    temp = (char*)malloc(sizeof(char)*(len+1));
    assert(temp);

    strcpy(temp, line);
    temp[len] = '\0';
    return temp;
}

// Outputs one watchtower information into a file
void PrintWtInfo(FILE *file, wt_info_t *entry) {
    fprintf(file, "Watchtower ID: %s, Postcode: %s, Population Served: %d, "
           "Watchtower Point of Contact Name: %s, x: %f, y: %f",
    entry->ID,entry->postcode,entry->population,entry->contact_name,entry->x,entry->y);
    if (entry->diameter != NO_DIAMETER) {
        fprintf(file, ", Diameter of Cell: %f", entry->diameter);
    }
    fprintf(file, "\n");
}

// Frees the string entries, watchtower entries and entire array
void FreeWts(wt_info_t **watchtowers, int *n) {
    int i;
    for (i=0;i<*n;i++) {
        wt_info_t *entry = watchtowers[i];
        free(entry->ID);
        free(entry->postcode);
        free(entry->contact_name);
        free(entry);
    }
    *n = 0;
    free(watchtowers);
}