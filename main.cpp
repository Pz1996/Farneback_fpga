#include <opencv2/opencv.hpp>
#include "opencv2/highgui.hpp"
#include "ref.h"
#include "opencv2/optflow.hpp"
#include <string>
#include <iostream>

#include "define.h"
#include "Farneback_of.h"

using namespace std;
using namespace cv;

#define BOUND 6400

#define L 4
#define R 6


void Mat_to_Array(Mat& m, pix_t* dst, int width, int height) {
	int count = 0;
	for (int i = 0; i < HEIGHT; i++) {
		for (int j = 0; j < WIDTH; j++) {
			dst[count++] = m.at<unsigned char>(i, j);
		}
	}

}

void DrawOptFlowMap(data_t **flow, Mat& cflowmap, int step, const Scalar& color) {
	for (int y = 0; y < cflowmap.rows; y += step)
		for (int x = 0; x < cflowmap.cols; x += step)
		{
			data_t fx = flow[y*cflowmap.cols + x][1];
			data_t fy = flow[y*cflowmap.cols + x][0];
			if (fx * fx + fy * fy > 1) {
				line(cflowmap, Point(x, y), Point(cvRound(x + fx), cvRound(y + fy)), color);
				circle(cflowmap, Point(cvRound(x + fx), cvRound(y + fy)), 1, color, -1);
			}
		}
}

void DrawOptFlowMap(const Mat& flow, Mat& cflowmap, int step, const Scalar& color) {
	for (int y = 0; y < cflowmap.rows; y += step)
		for (int x = 0; x < cflowmap.cols; x += step)
		{
			const Point2f& fxy = flow.at< Point2f>(y, x);
			if (fxy.x * fxy.x + fxy.y + fxy.y > 1) {
				line(cflowmap, Point(x, y), Point(cvRound(x + fxy.x), cvRound(y + fxy.y)), color);
				circle(cflowmap, Point(cvRound(x + fxy.x), cvRound(y + fxy.y)), 1, color, -1);
			}
		}
}

int main() {
	Mat pre_rgb = imread("pre.jpg", IMREAD_COLOR);
	Mat aft_rgb = imread("aft.jpg", IMREAD_COLOR);
	Mat pre, aft, pre_rgb2;
	pre_rgb.copyTo(pre_rgb2);
	cvtColor(pre_rgb, pre, CV_BGR2GRAY);
	cvtColor(aft_rgb, aft, CV_BGR2GRAY);;
	pix_t pre_img_1[MAXSIZE], aft_img_1[MAXSIZE], pre_img_2[MAXSIZE], aft_img_2[MAXSIZE];
	//data_t pre_poly[MAXSIZE][5], aft_poly[MAXSIZE][5], flow_in[MAXSIZE][2], flow_out[MAXSIZE][2];
	data_t ** pre_poly, ** aft_poly, **flow_in, **flow_out;
	pre_poly = new data_t *[MAXSIZE];
	aft_poly = new data_t *[MAXSIZE];
	flow_in = new data_t *[MAXSIZE];
	flow_out = new data_t *[MAXSIZE];
	for (int i = 0; i < MAXSIZE; i++) {
		pre_poly[i] = new data_t[5];
		aft_poly[i] = new data_t[5];
		flow_in[i] = new data_t[2];
		flow_out[i] = new data_t[2];
	}



	Mat_to_Array(pre, pre_img_1, WIDTH, HEIGHT);
	Mat_to_Array(aft, aft_img_1, WIDTH, HEIGHT);
	Resize(pre_img_1, pre_img_2, WIDTH, HEIGHT, 2);
	Resize(aft_img_1, aft_img_2, WIDTH, HEIGHT, 2);
	Smooth(pre_img_2, pre_img_1, WIDTH / 2, HEIGHT / 2);
	Smooth(aft_img_2, aft_img_1, WIDTH / 2, HEIGHT / 2);
	Poly_Exp(pre_img_1, pre_poly, WIDTH / 2, HEIGHT / 2);
	Poly_Exp(aft_img_1, aft_poly, WIDTH / 2, HEIGHT / 2);
	Displacement_Est(pre_poly, aft_poly, flow_out, flow_in, WIDTH / 2, HEIGHT / 2, 0);
	
	Mat_to_Array(pre, pre_img_2, WIDTH, HEIGHT);
	Mat_to_Array(aft, aft_img_2, WIDTH, HEIGHT);
	Smooth(pre_img_2, pre_img_1, WIDTH, HEIGHT);
	Smooth(aft_img_2, aft_img_1, WIDTH, HEIGHT);
	Poly_Exp(pre_img_1, pre_poly, WIDTH, HEIGHT);
	Poly_Exp(aft_img_1, aft_poly, WIDTH, HEIGHT);
	Displacement_Est(pre_poly, aft_poly, flow_in, flow_out, WIDTH, HEIGHT, 2);
	
	Mat flow;
	calcOpticalFlowFarneback(pre, aft, flow, 0.5, 2, 11, 3, 7, 1.5, 0);


	Mat k = aft_rgb - pre_rgb;


	DrawOptFlowMap(flow_out, pre_rgb, 5, CV_RGB(0, 255, 0));
	DrawOptFlowMap(flow, pre_rgb2, 5, CV_RGB(0, 255, 0));

	imshow("flow", k);
	imshow("fpga", pre_rgb);
	imshow("software", pre_rgb2);
	//imshow("aft", aft_rgb);

	waitKey(0);
	return 0;
}

/*
int main() {
	VideoCapture cap(0);
	if (!cap.isOpened())
		return -1;
	Mat pre, pre2, aft, pre_rgb, aft_rgb;
	while (1) {
		cap >> pre_rgb;
		waitKey(100);
		cap >> aft_rgb;
		pre_rgb.copyTo(pre2);
		cvtColor(pre_rgb, pre, CV_BGR2GRAY);
		cvtColor(aft_rgb, aft, CV_BGR2GRAY);;

		imwrite("pre.jpg", pre_rgb);
		imwrite("aft.jpg", aft_rgb);


		pix_t pre_img_1[MAXSIZE], aft_img_1[MAXSIZE], pre_img_2[MAXSIZE], aft_img_2[MAXSIZE];
		data_t pre_poly[MAXSIZE][5], aft_poly[MAXSIZE][5], flow_in[MAXSIZE][2], flow_out[MAXSIZE][2];
		
		Mat_to_Array(pre, pre_img_1, WIDTH, HEIGHT);
		Mat_to_Array(aft, aft_img_1, WIDTH, HEIGHT);
		Resize(pre_img_1, pre_img_2, WIDTH, HEIGHT, 2);
		Resize(aft_img_1, aft_img_2, WIDTH, HEIGHT, 2);
		Smooth(pre_img_2, pre_img_1, WIDTH / 2, HEIGHT / 2);
		Smooth(aft_img_2, aft_img_1, WIDTH / 2, HEIGHT / 2);
		Poly_Exp(pre_img_1, pre_poly, WIDTH / 2, HEIGHT / 2);
		Poly_Exp(aft_img_1, aft_poly, WIDTH / 2, HEIGHT / 2);
		Displacement_Est(pre_poly,aft_poly,flow_out,flow_in, WIDTH / 2, HEIGHT / 2,0);

		Mat_to_Array(pre, pre_img_2, WIDTH, HEIGHT);
		Mat_to_Array(aft, aft_img_2, WIDTH, HEIGHT);
		Smooth(pre_img_2, pre_img_1, WIDTH, HEIGHT);
		Smooth(aft_img_2, aft_img_1, WIDTH, HEIGHT);
		Poly_Exp(pre_img_1, pre_poly, WIDTH, HEIGHT);
		Poly_Exp(aft_img_1, aft_poly, WIDTH, HEIGHT);
		Displacement_Est(pre_poly, aft_poly, flow_in, flow_out, WIDTH, HEIGHT, 2);
		
		Mat flow;
		calcOpticalFlowFarneback(pre, aft, flow, 0.707, 5, 10, 3, 7, 1.5, cv::OPTFLOW_FARNEBACK_GAUSSIAN);

		DrawOptFlowMap(flow_out, pre_rgb, 5, CV_RGB(0, 255, 0));
		DrawOptFlowMap(flow, pre2, 5, CV_RGB(0, 255, 0));
		
		imshow("fpga", pre_rgb);
		imshow("software", pre2);
		imshow("aft", aft_rgb);

		waitKey(0);
	}
	return 0;
}
*/