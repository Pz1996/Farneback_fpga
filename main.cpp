#include <opencv2/opencv.hpp>
#include "opencv2/highgui.hpp"
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

void DrawOptFlowMap(data_t flow[MAXSIZE][2], Mat& cflowmap, int step, const Scalar& color) {
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
	
	VideoCapture cap(0);
	if (!cap.isOpened())
		return -1;
	Mat pre, aft, pre_rgb;
	while (1) {
		cap >> pre_rgb;
		waitKey(100);
		cap >> aft;
		cvtColor(pre_rgb, pre, CV_BGR2GRAY);
		cvtColor(aft, aft, CV_BGR2GRAY);
		cvtColor(pre_rgb, pre_rgb, CV_BGR2GRAY);
		/*
		pre.copyTo(aft);
		for (int i = 0; i < HEIGHT - L; i++) {
			for (int j = 0; j < WIDTH - R; j++) {
				aft.at<unsigned char>(i, j) = pre.at<unsigned char>(i + L, j + R);
			}
		}
		*/

		for (int i = 100; i < 120; i++) {
			for (int j = 100; j < 120; j++)
				cout << (int)(pre.at<unsigned char>(i, j)) << " ";
			cout << endl;
		}
		cout << endl;
		for (int i = 100; i < 120; i++) {
			for (int j = 100; j < 120; j++)
				cout << (int)(aft.at<unsigned char>(i, j)) << " ";
			cout << endl;
		}
		cout << endl;

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
		calcOpticalFlowFarneback(pre, aft, flow, 0.5, 2, 1, 1, 7, 1.5, cv::OPTFLOW_FARNEBACK_GAUSSIAN);


		DrawOptFlowMap(flow_out, pre_rgb, 5, CV_RGB(0, 255, 0));
		DrawOptFlowMap(flow, pre, 5, CV_RGB(0, 255, 0));
		
		imshow("fpga", pre_rgb);
		imshow("software", pre);
		imshow("aft", aft);

		waitKey(0);
	}
	
	return 0;
}


/*
Mat pre = imread("2hh.bmp", IMREAD_COLOR), aft;
cvtColor(pre, pre, CV_BGR2GRAY);
pre.copyTo(aft);
for (int i = 0; i < HEIGHT - L; i++) {
for (int j = 0; j < WIDTH - R; j++) {
aft.at<unsigned char>(i, j) = pre.at<unsigned char>(i + L, j + R);
}
}
//pix_t pre_in[MAXSIZE], aft_in[MAXSIZE];
pix_t pre_1[MAXSIZE], aft_1[MAXSIZE];
data_t pre_poly_1[MAXSIZE][5], aft_poly_1[MAXSIZE][5], flow_1[MAXSIZE][2];
pix_t pre_2[MAXSIZE], aft_2[MAXSIZE];
data_t pre_poly_2[MAXSIZE][5], aft_poly_2[MAXSIZE][5], flow_2[MAXSIZE][2];

for (int i = 0; i < HEIGHT * WIDTH; i++) {
pre_1[i] = pre.data[i];
aft_1[i] = aft.data[i];
}
for (int i = 0; i < HEIGHT / SCALING_FACTOR; i++) {
for (int j = 0; j < WIDTH / SCALING_FACTOR; j++) {
pre_2[i*WIDTH / SCALING_FACTOR + j] = (pre_1[i * SCALING_FACTOR * WIDTH + j * SCALING_FACTOR] + pre_1[i * SCALING_FACTOR * WIDTH + j * SCALING_FACTOR + WIDTH] + pre_1[i * SCALING_FACTOR * WIDTH + j * SCALING_FACTOR + 1] + pre_1[i * SCALING_FACTOR * WIDTH + j * SCALING_FACTOR + WIDTH + 1]) / SCALING_FACTOR /  SCALING_FACTOR;
aft_2[i*WIDTH / SCALING_FACTOR + j] = (aft_1[i * SCALING_FACTOR * WIDTH + j * SCALING_FACTOR] + aft_1[i * SCALING_FACTOR * WIDTH + j * SCALING_FACTOR + WIDTH] + aft_1[i * SCALING_FACTOR * WIDTH + j * SCALING_FACTOR + 1] + aft_1[i * SCALING_FACTOR * WIDTH + j * SCALING_FACTOR + WIDTH + 1]) / SCALING_FACTOR /  SCALING_FACTOR;
}
}
Poly_Exp(pre_2, pre_poly_2, WIDTH / SCALING_FACTOR, HEIGHT / SCALING_FACTOR);
Poly_Exp(aft_2, aft_poly_2, WIDTH / SCALING_FACTOR, HEIGHT / SCALING_FACTOR);
Displacement_Est(pre_poly_2, aft_poly_2, NULL, flow_1, WIDTH / SCALING_FACTOR, HEIGHT / SCALING_FACTOR, 0);
Poly_Exp(pre_1, pre_poly_1, WIDTH, HEIGHT);
Poly_Exp(aft_1, aft_poly_1, WIDTH, HEIGHT);
Displacement_Est(pre_poly_1, aft_poly_1, flow_1, flow_2, WIDTH, HEIGHT, SCALING_FACTOR);
for (int i = 10; i < HEIGHT - 10; i++) {
for (int j = 10; j < WIDTH - 10; j++) {
cout << flow_2[i*WIDTH + j][0] <<" "<< flow_2[i*WIDTH + j][1] << endl;
}
}
*/