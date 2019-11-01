#pragma once
#include "define.h"
#include <math.h>
#define POLY_EXP_SAMPLE_SIZE 15
#define DE_SAMPLE_SIZE 7
#define MAXSIZE (WIDTH)*(HEIGHT)
#define SCALING_FACTOR 5




void Resize(pix_t in[MAXSIZE], pix_t out[MAXSIZE], int width, int height, int scale = 1);
void Smooth(pix_t in[MAXSIZE], pix_t out[MAXSIZE], int width, int height);
void Poly_Exp(pix_t in[MAXSIZE], data_t out[MAXSIZE][5], int width, int height);
void Displacement_Est(data_t src_poly[MAXSIZE][5], data_t dst_poly[MAXSIZE][5], data_t flow_in[MAXSIZE][2], data_t flow_out[MAXSIZE][2], int width, int height, int scale = 1);
