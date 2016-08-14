
#ifndef SIFTDETECTION_H
#define SIFTDETECTION_H



#include <stdlib.h>
#include <stdio.h> 
#include <string.h> 
#include <math.h>
#include "stash.h"
#include "stdarg.h"

typedef unsigned char BYTE;

const int octaves = 4;

const int intervals = 2;

static unsigned int * byte2int (BYTE *image_b, int width, int height);


typedef struct KeypointSt {
  int x, y;                                  /* location of keypoint. */
  float feature_orentation, scale;           /* Scale and orientation (range [-PI,PI]) */
  float *feature_descriptor;                 /* Vector of descriptor values */
  struct KeypointSt *next;                   /* Pointer to next keypoint in list. */
} *Keypoint;


float *normalize(unsigned int *image, int width, int height);

float *filtering(float *image, int width, int height, float sigma, int length = 5 );

float **gaussian_coeff(float sigma, int length );

float bspline_func(const float x);

float *bicubic_zooming(float *image, int width, int height, int new_width, int new_height);

float ***sift_extraction(BYTE *image_b, int width, int height, 
					 float contrast_threshold = 0.03, float curvature_threshold = 10.0 );

void fix_edges(float *image, int edge_width, int width, int height);

float *bilinear_zooming(float *image, int width, int height, int new_width, int new_height);

int min_bi(int a, int b);


Stash ***keypoint_location(float ***differenceofgaussian, int width, int height,  
							float contrast_threshold = 0.03, float curvature_threshold = 10.0);

float findmax_inmatrix(float matrix[3][3][3]);

float findmin_inmatrix(float matrix[3][3][3]);


Stash *orientation_assignment(float ***gaussain_pyrmiad, Stash ***pos, float ***abs_sigma, int width, int height);

int finding_peak(float *peaks, int num_bins);

float parabolic_interp(float a, float fa, float b, float fb, float c, float fc);


void local_descriptor_calculation(float ***gaussain_pyrmiad, Stash *orientation_pos, int width, int height);

static float interp2(float *image, int width, int height, float fx, float fy);

Keypoint ReadKeys(FILE *fp);

int MatchCalculation(float *keypoint, Keypoint klist);

void FatalError(char *fmt, ...);

void SiftDrawLine(BYTE *image, int width, int height, int x1, int y1, int x2, int y2);


void SiftFindMatches();


///////////////////////////////////////////////////////////
//   Test section
///////////////////////////////////////////////////////////

float orientation_calc_test(BYTE *image, int width, int height);

void local_descriptor_calc(BYTE *image, int width, int height);

#endif