#pragma once
#define WIDTH 640
#define HEIGHT 480

#define POLY_EXP_SAMPLE_SIZE 7
#define MAXPOLY ((POLY_EXP_SAMPLE_SIZE) * (WIDTH) * (HEIGHT))
#define DE_SAMPLE_SIZE 11
#define MAXSIZE ((WIDTH)*(HEIGHT))
#define RAM_LENGTH 6348//(((HEIGHT-1)%POLY_EXP_SAMPLE_SIZE+1) * ((WIDTH-1)%POLY_EXP_SAMPLE_SIZE+1))

typedef unsigned char pix_t;
typedef float data_t;
#define SMALL_NUM 1e-3

