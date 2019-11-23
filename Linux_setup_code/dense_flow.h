#pragma once

#include <opencv2/opencv.hpp>
#include "opencv2/highgui.hpp"

using namespace cv;

void calcOpticalFlowFarneback(const Mat& prev0, const Mat& next0,
	Mat& flow0, double pyr_scale, int levels, int winsize,
	int iterations, int poly_n, double poly_sigma, int flags);

void FarnebackPolyExp(const Mat& src, Mat& dst, int n, double sigma);

void FarnebackUpdateMatrices(const Mat& _R0, const Mat& _R1, const Mat& _flow, Mat& _M, int _y0, int _y1);

void FarnebackUpdateFlow_Blur(const Mat& _R0, const Mat& _R1, Mat& _flow, Mat& _M, int block_size, bool update_matrices);
