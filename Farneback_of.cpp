#include "Farneback_of.h"
#include<iostream>
using namespace std;

void Resize(pix_t in[MAXSIZE], pix_t out[MAXSIZE], int width, int height, int scale)
{
	for (int i = 0; i < (height - 1) / scale + 1; i++) {
		for (int j = 0; j < (width - 1) / scale + 1; j++) {  // divide the original image into pieces size scale*scale
			int count = 0, sum = 0;
			for (int ii = i*scale; ii < (i + 1)*scale && ii < height; ii++) {
				for (int jj = j*scale; jj < (j + 1)*scale && jj < width; jj++) {
					count++;
					sum += in[ii * width + jj];
				}
			}
			sum = sum / count;  // calculate all the point and get the average value to stand for the block
			out[i * ((width - 1) / scale + 1) + j] = sum;
		}
	}
}

void Smooth(pix_t in[MAXSIZE], pix_t out[MAXSIZE], int width, int height)
{
	int coeff[25] = {
		1,	4,	6,	4,	1,
		4,	16,	24,	16,	4,
		6,	24,	36,	24,	6,
		4,	16,	24,	16,	4,
		1,	4,	6,	4,	1
	};
	// use a kernel to converlute the image
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			int cnt = 0;
			int sum = 0;
			for (int m = -2; m < 3; m++) {
				for (int n = -2; n < 3; n++) {
					// get the valid point position
					int px, py;
					px = i + m;
					py = j + n;
					if (px < 0) px = 0;
					if (px >= height) px = height - 1;
					if (py < 0) py = 0;
					if (py >= width) py = width - 1;
					sum += in[px * width + py] * coeff[cnt++];
				}
			}
			out[i*width + j] = (pix_t)(sum / 256);
		}
	}
}

void Poly_Exp(pix_t in[MAXSIZE], data_t out[MAXSIZE][5], int width, int height)
{
	//initialize the matrix
	int n = (POLY_EXP_SAMPLE_SIZE - 1) / 2;
	data_t coeff[5][POLY_EXP_SAMPLE_SIZE][POLY_EXP_SAMPLE_SIZE];
	data_t sigma = 1.5;
	data_t m[POLY_EXP_SAMPLE_SIZE], *g, s = 0;
	g = m + n;
	
	for (int x = -n; x <= n; x++) {
		g[x] = exp(-x*x / (2 * sigma*sigma));
		s += g[x];
	}
	// calculate a Gaussian distribution and normalize it, store in g[]
	for (int x = -n; x <= n; x++)
		g[x] = g[x] / s;
	
	data_t a, b, c, d, X;
	data_t ig00, ig11, ig03, ig33, ig55;
	// delta(d) = [B_ W B]^-1 B_ W f
	// B_ W B = G

	a = b = c = d = 0;
	// invG:
	// [ x        e  e    ]
	// [    y             ]
	// [       y          ]
	// [ e        z       ]
	// [ e           z    ]
	// [                u ]
	
	//calculate mat G
	for (int y = -n; y <= n; y++)
		for (int x = -n; x <= n; x++)
		{
			a += g[y] * g[x];
			b += g[y] * g[x] * x*x;
			c += g[y] * g[x] * x*x*x*x;
			d += g[y] * g[x] * x*x*y*y;
		}
	// calculate mat G_inv, according to the special structure
	X = a*c*c + b*b*d * 2 - a*d*d - b*b*c * 2;
	ig11 = 1 / b;
	ig33 = (a*c - b*b) / X;
	ig55 = 1 / d;
	ig03 = (b*d - c*b) / X;

	//calculate the coeff, which makes the poly_exp like 5 converlutional operate
	for (int i = -n; i <= n; i++) {
		for (int j = -n; j <= n; j++) {
			data_t b1, b2, b3, b4, b5, b6;
			b1 = g[i] * g[j];
			b2 = i * g[i] * g[j];
			b3 = j * g[i] * g[j];
			b4 = i * i * g[i] * g[j];
			b5 = j * j * g[i] * g[j];
			b6 = i * j * g[i] * g[j];
			coeff[0][i + n][j + n] = ig11 * b2;
			coeff[1][i + n][j + n] = ig11 * b3;
			coeff[2][i + n][j + n] = ig33 * b4 + ig03 * b1;
			coeff[3][i + n][j + n] = ig33 * b5 + ig03 * b1;
			coeff[4][i + n][j + n] = ig55 * b6;
		}
	}

	// dot_mul to get the poly_exp result
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			data_t sum[5];
			for (int k = 0; k < 5; k++)
				sum[k] = 0;
			for (int ii = -n; ii <= n; ii++) {
				for (int jj = -n; jj <= n; jj++) {
					int px, py;
					px = i + ii;
					px = (px > 0) ? px:0;
					px = (px < height) ? px : height - 1;
					py = j + jj;
					py = (py > 0) ? py : 0;
					py = (py < width) ? py:width - 1;
					for (int k = 0; k < 5; k++)
						sum[k] += coeff[k][ii+n][jj+n] * in[px * width + py];
				}
			}
			for (int k = 0; k < 5; k++)
				out[i*width + j][k] = sum[k];
		}
	}
}

void Displacement_Est(data_t src_poly[MAXSIZE][5], data_t dst_poly[MAXSIZE][5], data_t flow_in[MAXSIZE][2], data_t flow_out[MAXSIZE][2], int width, int height, int scale)
{
	data_t M[MAXSIZE][5], flow_t[MAXSIZE][2];
	UpdateMat(src_poly, dst_poly, flow_in, M, width, height, scale);
	UpdateFlow(M, flow_t, width, height);
	UpdateMat(src_poly, dst_poly, flow_t, M, width, height, 1);
	UpdateFlow(M, flow_out, width, height);


	/*
	SmoothFlow(flow_t, flow_in, width, height);
	UpdateMat(src_poly, dst_poly, flow_in, M, width, height, 1);
	UpdateFlow(M, flow_t, width, height);
	SmoothFlow(flow_t, flow_in, width, height);
	UpdateMat(src_poly, dst_poly, flow_in, M, width, height, 1);
	UpdateFlow(M, flow_t, width, height);
	
	SmoothFlow(flow_t, flow_out, width, height);
	*/
}

void UpdateMat(data_t src_poly[MAXSIZE][5], data_t dst_poly[MAXSIZE][5], data_t flow_in[MAXSIZE][2], data_t M[MAXSIZE][5], int width, int height, int scale)
{
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			int fx, fy, dx = 0, dy = 0;
			if (scale != 0) {// get the priori flow
				data_t dx_f, dy_f;
				dx_f = flow_in[(i / scale)*(width / scale) + (j / scale)][0] * scale;
				dy_f = flow_in[(i / scale)*(width / scale) + (j / scale)][1] * scale;
				dx = (dx_f >= 0) ? (int)(dx_f + 0.5) : (int)(dx_f - 0.5);
				dy = (dy_f >= 0) ? (int)(dy_f + 0.5) : (int)(dy_f - 0.5);
			}
			// get the valid dstination point
			fx = i + dx;
			fy = j + dy;
			fx = (fx < 0) ? 0 : (fx >= height) ? height - 1 : fx;
			fy = (fy < 0) ? 0 : (fy >= width) ? width - 1 : fy;
			data_t r[5], _r[5], a00, a01, a11, b0, b1;
			for (int k = 0; k < 5; k++) {
				r[k] = src_poly[i*width + j][k];
				_r[k] = dst_poly[fx * width + fy][k];
			}
			a00 = (r[2] + _r[2]) / 2; //r4
			a01 = (r[4] + _r[4]) / 4; //r6
			a11 = (r[3] + _r[3]) / 2; //r5
			b0 = (r[0] - _r[0]) / 2; //r2
			b1 = (r[1] - _r[1]) / 2; //r3

			b0 += a00 * dx + a01 * dy;
			b1 += a01 * dx + a11 * dy;

			//if ((dx != 0 || dy != 0)) {
			//	cout << dx << " " << dy << endl;
			//}


			//Calculate G and h
			M[i * width + j][0] = a00*a00 + a01*a01; // G(0, 0)
			M[i * width + j][1] = a01*(a00 + a11); // G(0, 1)
			M[i * width + j][2] = a11*a11 + a01*a01; // G(1, 1)
			M[i * width + j][3] = a00*b0 + a01*b1; // H(0)
			M[i * width + j][4] = a01*b0 + a11*b1; // H(1)
			/*
			data_t G00, G11;
			G00 = a00*a00 + a01*a01;
			G11 = a11*a11 + a01*a01;
			
			if (G00 > 100000 * G11 || G11 > 100000 * G00) {
				cout << "(" << i << "," << j << ")" << endl;
				for (int k = 0; k < 5; k++)
					cout << r[k] << " ";
				cout << endl;
				for (int k = 0; k < 5; k++)
					cout << _r[k] << " ";
				cout << endl;
			}
			*/

		}
	}
}

void UpdateFlow(data_t M[MAXSIZE][5], data_t flow_out[MAXSIZE][2], int width, int height)
{
	int n = (DE_SAMPLE_SIZE - 1) / 2;
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			int count = 0;
			data_t r[5];
			for (int k = 0; k < 5; k++)
				r[k] = 0;
			for (int ii = -n; ii <= n; ii++) {
				for (int jj = -n; jj <= n; jj++) {
					// find a neighbor and use the average parameter as G and h
					int fx, fy;
					fx = i + ii;
					fy = j + jj;
					if (fx >= 0 && fx < height && fy >= 0 && fy < width) {
						for (int k = 0; k < 5; k++)
							r[k] += M[fx*width + fy][k];
						count++;
					}
				}
			}
			for (int k = 0; k < 5; k++)
				r[k] = r[k] / count;
			data_t g00, g01, g11, h0, h1;
			g00 = r[0]; 
			g01 = r[1];
			g11 = r[2];
			h0 = r[3];
			h1 = r[4];

			data_t t_x, t_y, idet;
			idet = 1.0f / (g00*g11 - g01*g01 + SMALL_NUM);
			t_x = (g11*h0 - g01*h1) * idet;
			t_y = (g00*h1 - g01*h0) * idet;

			flow_out[i*width + j][0] = t_x;
			flow_out[i*width + j][1] = t_y;
		}
	}
}

void SmoothFlow(data_t flow_in[MAXSIZE][2], data_t flow_out[MAXSIZE][2], int width, int height)
{
	int n = 1;
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			data_t sum[2], max[2], score;
			int count = 0;
			sum[0] = sum[1] = max[0] = max[1] = score = 0;
			for (int ii = -n; ii <= n; ii++) {
				for (int jj = -n; jj <= n; jj++) {
					int fx = i + ii;
					int fy = j + jj;
					if (fx >= 0 && fx < height&&fy >= 0 && fy < width) {
						data_t flow_x, flow_y;

						if (fx < POLY_EXP_SAMPLE_SIZE / 2 || fx + POLY_EXP_SAMPLE_SIZE / 2 > height ||
							fy < POLY_EXP_SAMPLE_SIZE / 2 || fy + POLY_EXP_SAMPLE_SIZE / 2 > width){
							flow_x = 0;
							flow_y = 0;
						}
						else {
							flow_x = flow_in[fx*width + fy][0];
							flow_y = flow_in[fx*width + fy][1];
						}
							
						if (flow_x * flow_x + flow_y *flow_y > score) {
							max[0] = flow_x;
							max[1] = flow_y;
							score = flow_x * flow_x + flow_y *flow_y;
						}
						sum[0] += flow_x;
						sum[1] += flow_y;
						count++;
					}
				}
			}
			count--;
			flow_out[i*width + j][0] = (sum[0] - max[0]) / count;
			flow_out[i*width + j][1] = (sum[1] - max[1]) / count;
		}
	}

}




