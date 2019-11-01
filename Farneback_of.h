#pragma once
#include "define.h"
#include <math.h>
#define POLY_EXP_SAMPLE_SIZE 15
#define DE_SAMPLE_SIZE 7
#define MAXSIZE (WIDTH)*(HEIGHT)
#define SCALING_FACTOR 5

/** @brief Resize a mat into another one according to the scale factor.
@param in  input image.
@param out output image.
@param width parameter. specifying the input image's width
@param height parameter. specifying the input image's height
@param scale parameter. specifying the image scale (>1) to produce a new image which both width
and height are scale times smaller than the input image 
*/
void Resize(pix_t in[MAXSIZE], pix_t out[MAXSIZE], int width, int height, int scale = 1);

/** @brief Smooth the input image with a converlutional kernel
@param in  input image.
@param out output image.
@param width parameter. specifying the input image's width
@param height parameter. specifying the input image's height
*/
void Smooth(pix_t in[MAXSIZE], pix_t out[MAXSIZE], int width, int height);

/** @brief Calculate the polynomial expression  
@param in  input image.
@param out output polynomial expression which contains five dimensions
@param width parameter. specifying the input image's width
@param height parameter. specifying the input image's height
*/
void Poly_Exp(pix_t in[MAXSIZE], data_t out[MAXSIZE][5], int width, int height);

/** @brief Naive displacement estimation. Search the neighbors of the aim point and select the point
best fits the target point
@param src_poly input source polynomial expression
@param dst_poly input destination polynomial expression
@param flow_in input the priori optical flow
@param flow_out output the result flow of such function
@param width parameter. specifying the input image's width
@param height parameter. specifying the input image's height
@param scale parameter. show the scale factor between flow_in and the source images
*/
void Displacement_Est(data_t src_poly[MAXSIZE][5], data_t dst_poly[MAXSIZE][5], data_t flow_in[MAXSIZE][2], data_t flow_out[MAXSIZE][2], int width, int height, int scale = 1);
