#include "dense_flow.h"
#include <iostream>
using namespace std;

void calcOpticalFlowFarneback(const Mat & prev0, const Mat & next0, Mat & flow0, double pyr_scale, int levels, int winsize, int iterations, int poly_n, double poly_sigma, int flags)
{
	const int min_size = 32;
	const Mat* img[2] = { &prev0, &next0 };
	Mat fimg;

	int i, k;
	double scale;
	Mat prevFlow, flow;

	CV_Assert(prev0.size() == next0.size() && prev0.channels() == next0.channels() &&
		prev0.channels() == 1);
	flow0.create(prev0.size(), CV_32FC2);

	for (k = 0, scale = 1; k < levels; k++)
	{
		scale *= pyr_scale;
		if (prev0.cols*scale < min_size || prev0.rows*scale < min_size)
			break;
	}

	levels = k;

	for (k = levels; k >= 0; k--)
	{
		for (i = 0, scale = 1; i < k; i++)
			scale *= pyr_scale;

		double sigma = (1. / scale - 1)*0.5;
		int smooth_sz = cvRound(sigma * 5) | 1;
		smooth_sz = std::max(smooth_sz, 3);

		int width = cvRound(prev0.cols*scale);
		int height = cvRound(prev0.rows*scale);

		if (k > 0)
			flow.create(height, width, CV_32FC2);
		else
			flow = flow0;

		if (!prevFlow.data)
		{
			if (flags & OPTFLOW_USE_INITIAL_FLOW)
			{
				resize(flow0, flow, Size(width, height), 0, 0, INTER_AREA);
				flow *= scale;
			}
			else
				flow = Mat::zeros(height, width, CV_32FC2);
		}
		else
		{
			resize(prevFlow, flow, Size(width, height), 0, 0, INTER_LINEAR);
			flow *= 1. / pyr_scale;
		}

		Mat R[2], I, M;
		for (i = 0; i < 2; i++)
		{
			img[i]->convertTo(fimg, CV_32F);
			GaussianBlur(fimg, fimg, Size(smooth_sz, smooth_sz), sigma, sigma);
			resize(fimg, I, Size(width, height), CV_INTER_LINEAR);
			FarnebackPolyExp(I, R[i], poly_n, poly_sigma);
		}

		FarnebackUpdateMatrices(R[0], R[1], flow, M, 0, flow.rows);

		for (i = 0; i < iterations; i++)
			FarnebackUpdateFlow_Blur(R[0], R[1], flow, M, winsize, i < iterations - 1);

		prevFlow = flow;
	}
}

//#define DEBUG_EXP
#ifdef DEBUG_EXP
void FarnebackPolyExp(const Mat & src, Mat & dst, int n, double sigma) {
	pix_t in[MAXSIZE];
	data_t out[MAXSIZE][5];
	int width = src.cols;
	int height = src.rows;
	float* srca, *dsta;
	for (int i = 0; i < height; i++) {
		srca = (float*)(src.data + i*src.step);
		for (int j = 0; j < width; j++)
			in[i*width + j] = (pix_t)srca[j];
	}
	Poly_Exp(in, out, width, height);
	dst.create(height, width, CV_32FC(5));
	for (int i = 0; i < height; i++) {
		dsta = (float*)(dst.data + i*dst.step);
		for (int j = 0; j < width; j++) {
			dsta[j * 5] = out[i*width + j][0];
			dsta[j * 5 + 1] = out[i*width + j][1];
			dsta[j * 5 + 2] = out[i*width + j][2];
			dsta[j * 5 + 3] = out[i*width + j][3];
			dsta[j * 5 + 4] = out[i*width + j][4];
		}
	}
}
#else
void FarnebackPolyExp(const Mat & src, Mat & dst, int n, double sigma)
{
	int k, x, y;

	assert(src.type() == CV_32FC1);
	int width = src.cols;
	int height = src.rows;
	AutoBuffer<float> kbuf(n * 6 + 3), _row((width + n * 2) * 3);
	float* g = kbuf + n;
	float* xg = g + n * 2 + 1;
	float* xxg = xg + n * 2 + 1;
	float *row = (float*)_row + n * 3;

	if (sigma < FLT_EPSILON)
		sigma = n*0.3;

	double s = 0.;
	for (x = -n; x <= n; x++)
	{
		g[x] = (float)std::exp(-x*x / (2 * sigma*sigma));
		s += g[x];
	}

	s = 1. / s;
	for (x = -n; x <= n; x++)
	{
		g[x] = (float)(g[x] * s);
		xg[x] = (float)(x*g[x]);
		xxg[x] = (float)(x*x*g[x]);
	}

	Mat_<double> G = Mat_<double>::zeros(6, 6);

	for (y = -n; y <= n; y++)
		for (x = -n; x <= n; x++)
		{
			G(0, 0) += g[y] * g[x];
			G(1, 1) += g[y] * g[x] * x*x;
			G(3, 3) += g[y] * g[x] * x*x*x*x;
			G(5, 5) += g[y] * g[x] * x*x*y*y;
		}

	//G[0][0] = 1.;
	G(2, 2) = G(0, 3) = G(0, 4) = G(3, 0) = G(4, 0) = G(1, 1);
	G(4, 4) = G(3, 3);
	G(3, 4) = G(4, 3) = G(5, 5);

	// invG:
	// [ x        e  e    ]
	// [    y             ]
	// [       y          ]
	// [ e        z       ]
	// [ e           z    ]
	// [                u ]
	Mat_<double> invG = G.inv(DECOMP_CHOLESKY);
	double ig11 = invG(1, 1), ig03 = invG(0, 3), ig33 = invG(3, 3), ig55 = invG(5, 5);

	dst.create(height, width, CV_32FC(5));

	for (y = 0; y < height; y++)
	{
		float g0 = g[0], g1, g2;
		float *srow0 = (float*)(src.data + src.step*y), *srow1 = 0;
		float *drow = (float*)(dst.data + dst.step*y);

		// vertical part of convolution
		for (x = 0; x < width; x++)
		{
			row[x * 3] = srow0[x] * g0;
			row[x * 3 + 1] = row[x * 3 + 2] = 0.f;
		}

		for (k = 1; k <= n; k++)
		{
			g0 = g[k]; g1 = xg[k]; g2 = xxg[k];
			srow0 = (float*)(src.data + src.step*std::max(y - k, 0));
			srow1 = (float*)(src.data + src.step*std::min(y + k, height - 1));

			for (x = 0; x < width; x++)
			{
				float p = srow0[x] + srow1[x];
				float t0 = row[x * 3] + g0*p;
				float t1 = row[x * 3 + 1] + g1*(srow1[x] - srow0[x]);
				float t2 = row[x * 3 + 2] + g2*p;

				row[x * 3] = t0;
				row[x * 3 + 1] = t1;
				row[x * 3 + 2] = t2;
			}
		}

		// horizontal part of convolution
		for (x = 0; x < n * 3; x++)
		{
			row[-1 - x] = row[2 - x];
			row[width * 3 + x] = row[width * 3 + x - 3];
		}

		for (x = 0; x < width; x++)
		{
			g0 = g[0];
			// r1 ~ 1, r2 ~ x, r3 ~ y, r4 ~ x^2, r5 ~ y^2, r6 ~ xy
			double b1 = row[x * 3] * g0, b2 = 0, b3 = row[x * 3 + 1] * g0,
				b4 = 0, b5 = row[x * 3 + 2] * g0, b6 = 0;

			for (k = 1; k <= n; k++)
			{
				double tg = row[(x + k) * 3] + row[(x - k) * 3];
				g0 = g[k];
				b1 += tg*g0;
				b4 += tg*xxg[k];
				b2 += (row[(x + k) * 3] - row[(x - k) * 3])*xg[k];
				b3 += (row[(x + k) * 3 + 1] + row[(x - k) * 3 + 1])*g0;
				b6 += (row[(x + k) * 3 + 1] - row[(x - k) * 3 + 1])*xg[k];
				b5 += (row[(x + k) * 3 + 2] + row[(x - k) * 3 + 2])*g0;
			}

			// do not store r1
			drow[x * 5 + 1] = (float)(b2*ig11);
			drow[x * 5] = (float)(b3*ig11);
			drow[x * 5 + 3] = (float)(b1*ig03 + b4*ig33);
			drow[x * 5 + 2] = (float)(b1*ig03 + b5*ig33);
			drow[x * 5 + 4] = (float)(b6*ig55);
		}
	}

	row -= n * 3;
}
#endif

//#define DEBUG_MAT
#ifdef DEBUG_MAT
void FarnebackUpdateMatrices(const Mat & _R0, const Mat & _R1, const Mat & _flow, Mat & _M, int _y0, int _y1) {
	data_t src_poly[MAXSIZE][5];
	data_t dst_poly[MAXSIZE][5];
	data_t flow_in[MAXSIZE][2];
	data_t M[MAXSIZE][5];
	int width;
	int height;

	height = _R0.rows;
	width = _R0.cols;
	float *R0, *R1, *flow, *Ma;

	for (int i = 0; i < height; i++) {
		R0 = (float*)(_R0.data + i*_R0.step);
		R1 = (float*)(_R1.data + i*_R1.step);
		flow = (float*)(_flow.data + i*_flow.step);
		for (int j = 0; j < width; j++) {
			src_poly[i*width + j][0] = R0[j * 5 + 0];
			src_poly[i*width + j][1] = R0[j * 5 + 1];
			src_poly[i*width + j][2] = R0[j * 5 + 2];
			src_poly[i*width + j][3] = R0[j * 5 + 3];
			src_poly[i*width + j][4] = R0[j * 5 + 4];
			dst_poly[i*width + j][0] = R1[j * 5 + 0];
			dst_poly[i*width + j][1] = R1[j * 5 + 1];
			dst_poly[i*width + j][2] = R1[j * 5 + 2];
			dst_poly[i*width + j][3] = R1[j * 5 + 3];
			dst_poly[i*width + j][4] = R1[j * 5 + 4];
			flow_in[i*width + j][0] = flow[j * 2 + 1];
			flow_in[i*width + j][1] = flow[j * 2 ];
		}
	}

	UpdateMat(src_poly, dst_poly, flow_in, M, width, height, 1);
	_M.create(height, width, CV_32FC(5));
	for (int i = 0; i < height; i++) {
		Ma = (float*)(_M.data + i*_M.step);
		for (int j = 0; j < width; j++) {
			Ma[j * 5 + 0] = M[i*width + j][0];
			Ma[j * 5 + 1] = M[i*width + j][1];
			Ma[j * 5 + 2] = M[i*width + j][2];
			Ma[j * 5 + 3] = M[i*width + j][3];
			Ma[j * 5 + 4] = M[i*width + j][4];

		}
	}

}

#else
void FarnebackUpdateMatrices(const Mat & _R0, const Mat & _R1, const Mat & _flow, Mat & _M, int _y0, int _y1)
{
	const int BORDER = 5;
	static const float border[BORDER] = { 0.14f, 0.14f, 0.4472f, 0.4472f, 0.4472f };

	int x, y, width = _flow.cols, height = _flow.rows;
	const float* R1 = (float*)_R1.data;
	size_t step1 = _R1.step / sizeof(R1[0]);

	_M.create(height, width, CV_32FC(5));

	for (y = _y0; y < _y1; y++)
	{
		const float* flow = (float*)(_flow.data + y*_flow.step);
		const float* R0 = (float*)(_R0.data + y*_R0.step);
		float* M = (float*)(_M.data + y*_M.step);

		for (x = 0; x < width; x++)
		{
			float dx = flow[x * 2], dy = flow[x * 2 + 1];
			float fx = x + dx, fy = y + dy;

			int x1 = cvRound(fx), y1 = cvRound(fy);
			const float* ptr = R1 + y1*step1 + x1 * 5;
			float r2, r3, r4, r5, r6;

			if ((unsigned)x1 < (unsigned)width &&
				(unsigned)y1 < (unsigned)height)
			{
				r2 = ptr[0];
				r3 = ptr[1];
				r4 = (R0[x * 5 + 2] + ptr[2])*0.5f;
				r5 = (R0[x * 5 + 3] + ptr[3])*0.5f;
				r6 = (R0[x * 5 + 4] + ptr[4])*0.25f;
			}
			else
			{
				r2 = r3 = 0.f;
				r4 = R0[x * 5 + 2];
				r5 = R0[x * 5 + 3];
				r6 = R0[x * 5 + 4] * 0.5f;
			}

			r2 = (R0[x * 5] - r2)*0.5f;
			r3 = (R0[x * 5 + 1] - r3)*0.5f;

			r2 += r4*dy + r6*dx;
			r3 += r6*dy + r5*dx;

			if ((unsigned)(x - BORDER) >= (unsigned)(width - BORDER * 2) ||
				(unsigned)(y - BORDER) >= (unsigned)(height - BORDER * 2))
			{
				float scale = (x < BORDER ? border[x] : 1.f)*
					(x >= width - BORDER ? border[width - x - 1] : 1.f)*
					(y < BORDER ? border[y] : 1.f)*
					(y >= height - BORDER ? border[height - y - 1] : 1.f);

				r2 *= scale; r3 *= scale; r4 *= scale;
				r5 *= scale; r6 *= scale;
			}

			M[x * 5] = r4*r4 + r6*r6; // G(1,1)
			M[x * 5 + 1] = (r4 + r5)*r6;  // G(1,2)=G(2,1)
			M[x * 5 + 2] = r5*r5 + r6*r6; // G(2,2)
			M[x * 5 + 3] = r4*r2 + r6*r3; // h(1)
			M[x * 5 + 4] = r6*r2 + r5*r3; // h(2)
		}
	}
}
#endif

//#define DEBUG_FLOW
#ifdef DEBUG_FLOW
void FarnebackUpdateFlow_Blur(const Mat & _R0, const Mat & _R1, Mat & _flow, Mat & _M, int block_size, bool update_matrices) {
	data_t M[MAXSIZE][5];
	data_t flow_out[MAXSIZE][2];
	int width = _R0.cols;
	int height = _R0.rows;
	float *Ma, *flow;
	for (int i = 0; i < height; i++) {
		Ma = (float*)(_M.data + i*_M.step);
		for (int j = 0; j < width; j++) {
			for (int k = 0; k < 5; k++)
				M[i*width + j][k] = Ma[j * 5 + k];
		}
	}
	UpdateFlow(M, flow_out,width, height);
	for (int i = 0; i < height; i++) {
		flow = (float*)(_flow.data + i*_flow.step);
		for (int j = 0; j < width; j++) {
			flow[j * 2] = flow_out[i*width + j][1];
			flow[j * 2 + 1] = flow_out[i*width + j][0];
		}
	}
}
#else
void FarnebackUpdateFlow_Blur(const Mat & _R0, const Mat & _R1, Mat & _flow, Mat & _M, int block_size, bool update_matrices)
{
	int x, y, width = _flow.cols, height = _flow.rows;
	int m = block_size / 2;
	int y0 = 0, y1;
	int min_update_stripe = std::max((1 << 10) / width, block_size);
	double scale = 1. / (block_size*block_size);

	AutoBuffer<double> _vsum((width + m * 2 + 2) * 5);
	double* vsum = _vsum + (m + 1) * 5;

	// init vsum
	const float* srow0 = (const float*)_M.data;
	for (x = 0; x < width * 5; x++)
		vsum[x] = srow0[x] * (m + 2);

	for (y = 1; y < m; y++)
	{
		srow0 = (float*)(_M.data + _M.step*std::min(y, height - 1));
		for (x = 0; x < width * 5; x++)
			vsum[x] += srow0[x];
	}

	// compute blur(G)*flow=blur(h)
	for (y = 0; y < height; y++)
	{
		double g11, g12, g22, h1, h2;
		float* flow = (float*)(_flow.data + _flow.step*y);

		srow0 = (const float*)(_M.data + _M.step*std::max(y - m - 1, 0));
		const float* srow1 = (const float*)(_M.data + _M.step*std::min(y + m, height - 1));

		// vertical blur
		for (x = 0; x < width * 5; x++)
			vsum[x] += srow1[x] - srow0[x];

		// update borders
		for (x = 0; x < (m + 1) * 5; x++)
		{
			vsum[-1 - x] = vsum[4 - x];
			vsum[width * 5 + x] = vsum[width * 5 + x - 5];
		}

		// init g** and h*
		g11 = vsum[0] * (m + 2);
		g12 = vsum[1] * (m + 2);
		g22 = vsum[2] * (m + 2);
		h1 = vsum[3] * (m + 2);
		h2 = vsum[4] * (m + 2);

		for (x = 1; x < m; x++)
		{
			g11 += vsum[x * 5];
			g12 += vsum[x * 5 + 1];
			g22 += vsum[x * 5 + 2];
			h1 += vsum[x * 5 + 3];
			h2 += vsum[x * 5 + 4];
		}

		// horizontal blur
		for (x = 0; x < width; x++)
		{
			g11 += vsum[(x + m) * 5] - vsum[(x - m) * 5 - 5];
			g12 += vsum[(x + m) * 5 + 1] - vsum[(x - m) * 5 - 4];
			g22 += vsum[(x + m) * 5 + 2] - vsum[(x - m) * 5 - 3];
			h1 += vsum[(x + m) * 5 + 3] - vsum[(x - m) * 5 - 2];
			h2 += vsum[(x + m) * 5 + 4] - vsum[(x - m) * 5 - 1];

			double g11_ = g11*scale;
			double g12_ = g12*scale;
			double g22_ = g22*scale;
			double h1_ = h1*scale;
			double h2_ = h2*scale;

			double idet = 1. / (g11_*g22_ - g12_*g12_ + 1e-3);

			flow[x * 2] = (float)((g11_*h2_ - g12_*h1_)*idet);
			flow[x * 2 + 1] = (float)((g22_*h1_ - g12_*h2_)*idet);
		}

		y1 = y == height - 1 ? height : y - block_size;
		if (update_matrices && (y1 == height || y1 >= y0 + min_update_stripe))
		{
			FarnebackUpdateMatrices(_R0, _R1, _flow, _M, y0, y1);
			y0 = y1;
		}
	}
}
#endif
