#include "Farneback_of.h"
#include<iostream>
using namespace std;



void Resize(pix_t in[MAXSIZE], pix_t out[MAXSIZE], int width, int height, int scale)
{
	for (int i = 0; i < (height - 1) / scale + 1; i++) {
		for (int j = 0; j < (width - 1) / scale + 1; j++) {
			int count = 0, sum = 0;
			for (int ii = i*scale; ii < (i + 1)*scale && ii < height; ii++) {
				for (int jj = j*scale; jj < (j + 1)*scale && jj < width; jj++) {
					count++;
					sum += in[ii * width + jj];
				}
			}
			sum = sum / count;
			out[i * ((width - 1) / scale + 1) + j] = sum;
		}
	}
}

void Smooth(pix_t in[MAXSIZE], pix_t out[MAXSIZE], int width, int height)
{
	int coeff[25] = {
		1,	2,	4,	2,	1,
		2,	4,	8,	4,	2,
		4,	8,	16,	8,	4,
		2,	4,	8,	4,	2,
		1,	2,	4,	2,	1
	};

	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			int cnt = 0;
			int sum = 0;
			for (int m = -2; m < 3; m++) {
				for (int n = -2; n < 3; n++) {
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
			out[i*width + j] = (pix_t)(sum / 100);
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
	for (int x = -n; x <= n; x++)
		g[x] = g[x] / s;
	
	data_t a, b, c, d, X;
	data_t ig00, ig11, ig03, ig33, ig55;
	a = b = c = d = 0;
	// invG:
	// [ x        e  e    ]
	// [    y             ]
	// [       y          ]
	// [ e        z       ]
	// [ e           z    ]
	// [                u ]
	for (int y = -n; y <= n; y++)
		for (int x = -n; x <= n; x++)
		{
			a += g[y] * g[x];
			b += g[y] * g[x] * x*x;
			c += g[y] * g[x] * x*x*x*x;
			d += g[y] * g[x] * x*x*y*y;
		}
	X = a*c*c + b*b*d * 2 - a*d*d - b*b*c * 2;
	//ig00 = (c*c - d*d) / X;
	ig11 = 1 / b;
	ig33 = (a*c - b*b) / X;
	ig55 = 1 / d;
	ig03 = (b*d - c*b) / X;
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
	for (int i = 0; i < height; i++) {
		for (int j = 0; j < width; j++) {
			int fx, fy, dx = 0, dy = 0;
			int flow_x = 0, flow_y = 0;
			data_t score = 1e10;
			int n = (DE_SAMPLE_SIZE - 1) / 2;
			data_t sigma = 1.5;
			data_t m[DE_SAMPLE_SIZE], *g, s = 0;
			g = m + n;
			for (int x = -n; x <= n; x++) {
				g[x] = exp(-x*x / (2 * sigma*sigma));
				s += g[x];
			}
			for (int x = -n; x <= n; x++)
				g[x] = g[x] / s;


			if (scale != 0) {
				data_t dx_f, dy_f;
				dx_f = flow_in[(i / scale)*(width / scale) + (j / scale)][0] * scale;
				dy_f = flow_in[(i / scale)*(width / scale) + (j / scale)][1] * scale;
				dx = (dx_f >= 0) ? (int)(dx_f + 0.5) : (int)(dx_f - 0.5);
				dy = (dy_f >= 0) ? (int)(dy_f + 0.5) : (int)(dy_f - 0.5);
			}
			fx = i + dx;
			fy = j + dy;
			if (fx < n) fx = n;
			if (fx >= height - n) fx = height - n - 1;
			if (fy < n) fy = n;
			if (fy >= width - n) fy = width - n - 1;

			data_t sum_x = 0, sum_y = 0;
			int count = 0;

			for (int ii = -n; ii <= n; ii++) {
				for (int jj = -n; jj <= n; jj++) {
					data_t r[5], _r[5], a00, a11, a01, b0, b1, _a00, _a01, _a11;
					data_t t_x, t_y;
					for (int k = 0; k < 5; k++) {
						r[k] = src_poly[i * width + j][k];
						_r[k] = dst_poly[(fx + ii) * width + fy + jj][k];
					}
					a00 = (r[2] + _r[2]) / 2;
					a01 = (r[4] + _r[4]) / 4;
					a11 = (r[3] + _r[3]) / 2;
					b0 = (r[0] - _r[0]) / 2;
					b1 = (r[1] - _r[1]) / 2;
					//如果周围一片都类似，那么可能随机选点
					if (b0<1e-9 && b0 >  -1e-9)
						b0 = 1e-9;
					if (b1<1e-9 && b1 >  -1e-9)
						b1 = 1e-9;

					_a00 = a11 / (a00*a11 - a01*a01);
					_a01 = -a01 / (a00*a11 - a01*a01);
					_a11 = a00 / (a00*a11 - a01*a01);
					t_x = _a00 * b0 - _a01 * b1;
					t_y = -_a01 * b0 + _a11 * b1;

					if (t_x*t_x + t_y*t_y + (ii*ii+jj*jj)/1e3 < score) {
						//cout << t_x*t_x + t_y*t_y << " "<< (ii*ii + jj*jj) / 1e2 << endl;
						score = t_x*t_x + t_y*t_y + (ii*ii + jj*jj) / 1e3;
						flow_x = fx + ii - i + t_x;
						flow_y = fy + jj - j + t_y;
					}
				}

			}
		
			flow_out[i * width + j][0] = flow_x;
			flow_out[i * width + j][1] = flow_y;
		}
	}
}


