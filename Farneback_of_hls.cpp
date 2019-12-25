#include "Farneback_of_hls.h"

void Smooth_hls(hls::stream<pix_t>& in, hls::stream<pix_t>&out, int width, int height){
#pragma HLS DATAFLOW
	const int K = 5;
	const unsigned char coeff[5] = {1, 4, 6, 4, 1};
	unsigned int hwin[5];
	hls::stream<unsigned int> hconv("hconv");
	unsigned int linebuf[4][WIDTH];
#pragma HLS ARRAY_PARTITION variable=linebuf dim=1 complete
	hls::stream<unsigned int> vconv("vconv");
	const int vconv_xlim = width + 1 - K;

	// Horizontal convolution
	HConvH:for(int i=0;i<height;i++){
#pragma HLS LOOP_TRIPCOUNT min=240 max=480
		HConvW:for(int j=0;j<width;j++){
#pragma HLS LOOP_TRIPCOUNT min=320 max=640
#pragma HLS PIPELINE II=5
			unsigned int in_val = in.read();
			unsigned int out_val = 0;
			Hconv:for(int k=0;k<K;k++){
				hwin[k] = k < K-1 ? hwin[k+1] : in_val;
				out_val = out_val + hwin[k] * coeff[k];
			}
			if(j >= K-1)
				hconv.write(out_val);
		}
	}

	// Vertical convolution
	VConvH:for(int i=0;i<height;i++){
#pragma HLS LOOP_TRIPCOUNT min=240 max=480
		VconvW:for(int j=0;j<vconv_xlim;j++){
#pragma HLS LOOP_TRIPCOUNT min=320 max=640
#pragma HLS DEPENDENCE variable=linebuf inter false
#pragma HLS PIPELINE II=5
			unsigned int in_val = hconv.read();
			unsigned int out_val = 0;
			Vconv:for(int k=0; k<K; k++){
				unsigned int vwin_val = k < K-1? linebuf[k][j]: in_val;
				out_val = out_val + vwin_val * coeff[k];
				if(k > 0)
					linebuf[k-1][j] = vwin_val;
			}
			if(i >= K-1)
				vconv.write(out_val);
		}
	}

	const int border_width = int(K/2);
	unsigned int borderbuf[WIDTH + 1 -K];

	BorderH:for(int i=0; i<height; i++){
#pragma HLS LOOP_TRIPCOUNT min=240 max=480
		BorderW:for(int j=0; j<width; j++){
#pragma HLS PIPELINE II=5
#pragma HLS LOOP_TRIPCOUNT min=320 max=640
			unsigned int d_in, l_edge, r_edge, d_out;
			if (i == 0 || (i > border_width && i < height - border_width)) {
				if (j < width - (K - 1)) {
					d_in = vconv.read();
			        borderbuf[j] = d_in;
				}
				if (j == 0) {
					l_edge = d_in;
				}
				if (j == width - K) {
					r_edge = d_in;
				}
			}
			if (j <= border_width) {
				d_out = l_edge;
			} else if (j >= width - border_width - 1) {
				d_out = r_edge;
			} else {
				d_out = borderbuf[j - border_width];
			}
			out.write(d_out / 256);

		}
	}

}

void Poly_Exp_hls_strm(hls::stream<pix_t>& in, hls::stream<Data_5> &out, int width, int height){
#pragma HLS DATAFLOW
	//initialize the matrix
	int n = (POLY_EXP_SAMPLE_SIZE - 1) / 2;
	/*
	data_t sigma = 1.5;
	data_t m[POLY_EXP_SAMPLE_SIZE], *g, s = 0;
#pragma HLS ARRAY_PARTITION variable=m complete dim=1
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
	a = b = c = d = 0;
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
	 */

	data_t m[POLY_EXP_SAMPLE_SIZE] = {0.0366328, 0.111281, 0.216745, 0.270682, 0.216745, 0.111281, 0.0366328 };
#pragma HLS ARRAY_PARTITION variable=m complete dim=1
	data_t *g = m + n;
	data_t ig11 = 0.504254;
	data_t ig33 = 0.166772;
	data_t ig55 = 0.254272;
	data_t ig03 = -0.330731;
	pix_t hwin[POLY_EXP_SAMPLE_SIZE];
	hls::stream<Data_3> hconv("hconv");

	Data_3 linebuf[POLY_EXP_SAMPLE_SIZE-1][WIDTH];
#pragma HLS ARRAY_PARTITION variable=linebuf dim=1 complete
	hls::stream<Data_5> vconv("vconv");
	Data_5 borderbuf[WIDTH + 1 -POLY_EXP_SAMPLE_SIZE];

	HConvH:for(int i=0;i<height;i++){
#pragma HLS LOOP_TRIPCOUNT min=240 max=480
		HConvW:for(int j=0;j<width;j++){
#pragma HLS PIPELINE II=5
#pragma HLS LOOP_TRIPCOUNT min=320 max=640
			pix_t in_val = in.read();
			Data_3 out_val = {0,0,0};
			HConv:for(int k=0;k<POLY_EXP_SAMPLE_SIZE;k++){
				hwin[k] = k < (POLY_EXP_SAMPLE_SIZE - 1) ? hwin[k+1]:in_val;
				int p = k-n;
				out_val.r += hwin[k] * g[p];
				out_val.xr += hwin[k] * p * g[p];
				out_val.xxr += hwin[k] * p * p * g[p];
			}
			if(j >= POLY_EXP_SAMPLE_SIZE - 1)
				hconv << out_val;
		}
	}

	VConvH:for(int i=0;i<height;i++){
#pragma HLS LOOP_TRIPCOUNT min=240 max=480
		VConvW:for(int j=0;j<width + 1 - POLY_EXP_SAMPLE_SIZE;j++){
#pragma HLS DEPENDENCE variable=linebuf inter false
#pragma HLS PIPELINE II=5
#pragma HLS LOOP_TRIPCOUNT min=320 max=640
			Data_3 in_val = hconv.read();
			data_t b1, b2, b3, b4, b5, b6;
			Data_5 out_val;
			b1 = b2 = b3 = b4 = b5 = b6 = 0;
			VConv:for(int k=0;k<POLY_EXP_SAMPLE_SIZE;k++){
				Data_3 vwin_val = k < (POLY_EXP_SAMPLE_SIZE - 1) ? linebuf[k][j]:in_val;
				int p = k-n;
				b1 += vwin_val.r * g[p];
				b2 += vwin_val.r * p * g[p];
			    b3 += vwin_val.xr * g[p];
				b4 += vwin_val.r * p * p * g[p];
				b5 += vwin_val.xxr * g[p];
				b6 += vwin_val.xr * p * g[p];
				if(k > 0)
					linebuf[k-1][j] = vwin_val;

			}
			if(i >= POLY_EXP_SAMPLE_SIZE - 1){
				out_val.r0 = ig11 * b2;
				out_val.r1 = ig11 * b3;
				out_val.r2 = ig33 * b4 + ig03 * b1;
				out_val.r3 = ig33 * b5 + ig03 * b1;
				out_val.r4 = ig55 * b6;
				vconv.write(out_val);
			}
		}
	}

	const int border_width = int(POLY_EXP_SAMPLE_SIZE/2);
	BorderH:for(int i=0; i<height; i++){
#pragma HLS LOOP_TRIPCOUNT min=240 max=480
		BorderW:for(int j=0; j<width; j++){
#pragma HLS PIPELINE II=5
#pragma HLS LOOP_TRIPCOUNT min=320 max=640
			Data_5 d_in, l_edge, r_edge, d_out;
			if (i == 0 || (i > border_width && i < height - border_width)) {
				if (j < width - (POLY_EXP_SAMPLE_SIZE - 1)) {
					d_in = vconv.read();
			        borderbuf[j] = d_in;
				}
				if (j == 0) {
					l_edge = d_in;
				}
				if (j == width - POLY_EXP_SAMPLE_SIZE) {
					r_edge = d_in;
				}
			}
			if (j <= border_width) {
				d_out = l_edge;
			} else if (j >= width - border_width - 1) {
				d_out = r_edge;
			} else {
				d_out = borderbuf[j - border_width];
			}
			out.write(d_out);
		}
	}
}

void Poly_Exp_hls_bram(pix_t img_in[POLY_EXP_SAMPLE_SIZE*POLY_EXP_SAMPLE_SIZE][RAM_LENGTH],
		hls::stream<Data_2u>& flow_in, hls::stream<Data_5> &out, int width, int height){
#pragma HLS ARRAY_PARTITION variable=img_in complete dim=1
#pragma HLS DATAFLOW
	int n = (POLY_EXP_SAMPLE_SIZE - 1) / 2;
	int row_blocks = (width - 1) / POLY_EXP_SAMPLE_SIZE + 1;
	/*
	data_t sigma = 1.5;
	data_t m[POLY_EXP_SAMPLE_SIZE], xm[POLY_EXP_SAMPLE_SIZE], xxm[POLY_EXP_SAMPLE_SIZE];
	data_t *g, *xg, *xxg, s = 0;
#pragma HLS ARRAY_PARTITION variable=m complete dim=1
#pragma HLS ARRAY_PARTITION variable=xm complete dim=1
#pragma HLS ARRAY_PARTITION variable=xxm complete dim=1
	g = m + n;
	xg = xm + n;
	xxg = xxm + n;

	for (int x = -n; x <= n; x++) {
		g[x] = exp(-x*x / (2 * sigma*sigma));
		s += g[x];
	}
	// calculate a Gaussian distribution and normalize it, store in g[]
	for (int x = -n; x <= n; x++){
		g[x] = g[x] / s;
		xg[x] = x * g[x];
		xxg[x] = x * x * g[x];
	}
	data_t a, b, c, d, X;
	data_t ig00, ig11, ig03, ig33, ig55;
	a = b = c = d = 0;
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

	cout << "data_t m[POLY_EXP_SAMPLE_SIZE] = {";
	for(int i=0;i<POLY_EXP_SAMPLE_SIZE;i++)
		cout << m[i]<<", ";
	cout << "};"<<endl;

	cout << "data_t xm[POLY_EXP_SAMPLE_SIZE] = {";
	for(int i=0;i<POLY_EXP_SAMPLE_SIZE;i++)
		cout << xm[i]<<", ";
	cout << "};"<<endl;

	cout << "data_t xxm[POLY_EXP_SAMPLE_SIZE] = {";
	for(int i=0;i<POLY_EXP_SAMPLE_SIZE;i++)
		cout << xxm[i]<<", ";
	cout << "};"<<endl;

	cout << "data_t ig11 = "<< ig11 << ";" << endl;
	cout << "data_t ig33 = "<< ig33 << ";" << endl;
	cout << "data_t ig55 = "<< ig55 << ";" << endl;
	cout << "data_t ig03 = "<< ig03 << ";" << endl;
	*/
	data_t m[POLY_EXP_SAMPLE_SIZE] = {0.0366328, 0.111281, 0.216745, 0.270682, 0.216745, 0.111281, 0.0366328 };
	data_t xm[POLY_EXP_SAMPLE_SIZE] = {-0.109899, -0.222562, -0.216745, 0, 0.216745, 0.222562, 0.109899 };
	data_t xxm[POLY_EXP_SAMPLE_SIZE] = {0.329696, 0.445123, 0.216745, 0, 0.216745, 0.445123, 0.329696 };
#pragma HLS ARRAY_PARTITION variable=m complete dim=1
#pragma HLS ARRAY_PARTITION variable=xm complete dim=1
#pragma HLS ARRAY_PARTITION variable=xxm complete dim=1
	data_t ig11 = 0.504254;
	data_t ig33 = 0.166772;
	data_t ig55 = 0.254272;
	data_t ig03 = -0.330731;

	for(int i=0;i<height;i++){
#pragma HLS LOOP_TRIPCOUNT min=240 max=480
		for(int j=0;j<width;j++){
#pragma HLS LOOP_TRIPCOUNT min=320 max=640
#pragma HLS PIPELINE II=5
			Data_2u in_val = flow_in.read();
			data_t b1,b2,b3,b4,b5,b6;
			b1 = b2 = b3 = b4 = b5 = b6 = 0;
			int fx, fy;
			fx = i + in_val.r0 - n;
			fy = j + in_val.r1 - n;
			fx = (fx < 0) ? 0 : (fx > height -POLY_EXP_SAMPLE_SIZE) ? height - POLY_EXP_SAMPLE_SIZE: fx;
			fy = (fy < 0) ? 0 : (fy > width -POLY_EXP_SAMPLE_SIZE) ? width - POLY_EXP_SAMPLE_SIZE: fy;
			for(int k=0;k<POLY_EXP_SAMPLE_SIZE*POLY_EXP_SAMPLE_SIZE;k++){
				int dx, dy;
				dx = (MAXPOLY + k/POLY_EXP_SAMPLE_SIZE - fx) % POLY_EXP_SAMPLE_SIZE;
				dy = (MAXPOLY + k%POLY_EXP_SAMPLE_SIZE - fy) % POLY_EXP_SAMPLE_SIZE;
				int addr = (fx + dx) / POLY_EXP_SAMPLE_SIZE * row_blocks + (fy + dy)/POLY_EXP_SAMPLE_SIZE;
				pix_t in_val = img_in[k][addr];
				b1 += in_val * m[dx] * m[dy];
				b2 += in_val * xm[dx] * m[dy];
				b3 += in_val * m[dx] * xm[dy];
				b4 += in_val * xxm[dx] * m[dy];
				b5 += in_val * m[dx] * xxm[dy];
				b6 += in_val * xm[dx] * xm[dy];
			}

			Data_5 out_val;
			out_val.r0 = ig11 * b2;
			out_val.r1 = ig11 * b3;
			out_val.r2 = ig33 * b4 + ig03 * b1;
			out_val.r3 = ig33 * b5 + ig03 * b1;
			out_val.r4 = ig55 * b6;
			out.write(out_val);
		}
	}

}

void UpdateMat_0_hls(hls::stream<Data_5> &src_poly, hls::stream<Data_5> &dst_poly, hls::stream<Data_5>& M,
		int width, int height){
	UpdataMat_0H:for(int i=0;i<height;i++){
#pragma HLS LOOP_TRIPCOUNT min=240 max=480
		UpdataMat_0W:for(int j=0;j<width;j++){
#pragma HLS LOOP_TRIPCOUNT min=320 max=640
#pragma HLS PIPELINE II=5
			Data_5 src, dst, m;
			data_t a00, a01, a11, b0, b1;
			src = src_poly.read();
			dst = dst_poly.read();
			a00 = (src.r2 + dst.r2) / 2; //r4
			a01 = (src.r4 + dst.r4) / 4; //r6
			a11 = (src.r3 + dst.r3) / 2; //r5
			b0 = (src.r0 - dst.r0) / 2; //r2
			b1 = (src.r1 - dst.r1) / 2; //r3

			m.r0 = a00*a00 + a01*a01; // G(0, 0)
			m.r1 = a01*(a00 + a11); // G(0, 1)
			m.r2 = a11*a11 + a01*a01; // G(1, 1)
			m.r3 = a00*b0 + a01*b1; // H(0)
			m.r4 = a01*b0 + a11*b1; // H(1)
			M.write(m);
		}
	}
}

void UpdateMat_hls(hls::stream<Data_5> &src_poly, hls::stream<Data_5> &dst_poly,hls::stream<Data_2u> &flow,
		hls::stream<Data_5>& M, int width, int height){
	UpdataMat_H:for(int i=0;i<height;i++){
#pragma HLS LOOP_TRIPCOUNT min=240 max=480
		UpdataMat_W:for(int j=0;j<width;j++){
#pragma HLS LOOP_TRIPCOUNT min=320 max=640
#pragma HLS PIPELINE II=5
			Data_5 src, dst, m;
			data_t a00, a01, a11, b0, b1;
			Data_2u flow_val;
			flow_val = flow.read();
			int dx, dy;
			dx = flow_val.r0;
			dy = flow_val.r1;
			src = src_poly.read();
			dst = dst_poly.read();
			a00 = (src.r2 + dst.r2) / 2; //r4
			a01 = (src.r4 + dst.r4) / 4; //r6
			a11 = (src.r3 + dst.r3) / 2; //r5
			b0 = (src.r0 - dst.r0) / 2; //r2
			b1 = (src.r1 - dst.r1) / 2; //r3

			b0 += a00 * dx + a01 * dy;
			b1 += a01 * dx + a11 * dy;

			m.r0 = a00*a00 + a01*a01; // G(0, 0)
			m.r1 = a01*(a00 + a11); // G(0, 1)
			m.r2 = a11*a11 + a01*a01; // G(1, 1)
			m.r3 = a00*b0 + a01*b1; // H(0)
			m.r4 = a01*b0 + a11*b1; // H(1)
			M.write(m);
		}
	}
}

void UpdateFlow_hls(hls::stream<Data_5>&M, hls::stream<Data_2>&flow_out, int width, int height){
#pragma HLS DATAFLOW
	const int K = DE_SAMPLE_SIZE;
	// Horizontal pixel window
	Data_5 hwin[DE_SAMPLE_SIZE];
	hls::stream<Data_5> hconv("hconv");
	// Vertical pixel window
	static Data_5 linebuf[DE_SAMPLE_SIZE-1][WIDTH];
#pragma HLS ARRAY_PARTITION variable=linebuf dim=1 complete
	hls::stream<Data_5> vconv("vconv");
	const int vconv_xlim = width - (K - 1);

	// Horizontal convolution
	HConvH:for(int i=0;i<height;i++){
#pragma HLS LOOP_TRIPCOUNT min=240 max=480
		HconvW:for(int j=0;j<width;j++){
#pragma HLS LOOP_TRIPCOUNT min=320 max=640
#pragma HLS PIPELINE II=5
			Data_5 in_val = M.read();
			Data_5 out_val = {0,0,0,0,0};
			Hconv:for(int k=0; k<K; k++){
				hwin[k] = k < K-1 ? hwin[k+1] : in_val;
				out_val = out_val + hwin[k];
			}
			if(j >= K-1)
				hconv << out_val;
		}
	}

	// Vertical convolution
	VConvH:for(int i=0;i<height;i++){
#pragma HLS LOOP_TRIPCOUNT min=240 max=480
		VconvW:for(int j=0;j<vconv_xlim;j++){
#pragma HLS LOOP_TRIPCOUNT min=320 max=640
#pragma HLS DEPENDENCE variable=linebuf inter false
#pragma HLS PIPELINE II=5
			Data_5 in_val = hconv.read();
			Data_5 out_val = {0,0,0,0,0};
			Vconv:for(int k=0; k<K; k++){
				Data_5 vwin_val = k < K-1? linebuf[k][j]: in_val;
				out_val = out_val + vwin_val;
				if(k > 0)
					linebuf[k-1][j] = vwin_val;
			}
			if(i >= K-1)
				vconv << out_val;
		}
	}

	for(int i=0;i<height;i++){
#pragma HLS LOOP_TRIPCOUNT min=240 max=480
		for(int j=0;j<width;j++){
#pragma HLS LOOP_TRIPCOUNT min=320 max=640
#pragma HLS PIPELINE II=5
			Data_2 out = {0,0};
			if((unsigned)(i - K/2) < height-K+1 && (unsigned)(j - K/2) < width-K+1){
				Data_5 m_in = vconv.read()/(K*K);
				data_t g00, g01, g11, h0, h1;
				g00 = m_in.r0;
				g01 = m_in.r1;
				g11 = m_in.r2;
				h0 = m_in.r3;
				h1 = m_in.r4;

				data_t idet;
				idet = 1.0f / (g00*g11 - g01*g01 + SMALL_NUM);
				out.r0 = (g11*h0 - g01*h1) * idet;
				out.r1 = (g00*h1 - g01*h0) * idet;
			}
			flow_out.write(out);
		}
	}
}

void Farneback_top(volatile pix_t* mig_in, volatile data_t* mig_out){
#pragma HLS INTERFACE s_axilite port=return
#pragma HLS INTERFACE ap_ctrl_hs port=return
#pragma HLS INTERFACE m_axi depth=640000 port=mig_in
#pragma HLS INTERFACE m_axi depth=3200000 port=mig_out
#pragma HLS DATAFLOW
	hls::stream<pix_t> src_img_orig_strm("src_img_orig_strm");
	hls::stream<pix_t> dst_img_orig_strm("dst_img_orig_strm");
	hls::stream<pix_t> src_img_strm("src_img_strm");
	hls::stream<pix_t> dst_img_strm("dst_img_strm");
	hls::stream<Data_5> src_poly("src_poly");
#pragma HLS DATA_PACK variable=src_poly
	hls::stream<Data_5> dst_poly("dst_poly");
#pragma HLS DATA_PACK variable=dst_poly
	hls::stream<Data_5> M("M");
#pragma HLS DATA_PACK variable=M
	hls::stream<Data_2> flow("flow");
#pragma HLS DATA_PACK variable=flow

	for(int i=0;i<MAXSIZE;i++){
#pragma HLS PIPELINE II=5
		pix_t tmp = mig_in[i*2];
		pix_t tmp2 = mig_in[i*2+1];
		src_img_orig_strm << tmp;
		dst_img_orig_strm << tmp2;
	}

	//Smooth_hls(src_img_orig_strm, src_img_strm, WIDTH, HEIGHT);
	//Smooth_hls(dst_img_orig_strm, dst_img_strm, WIDTH, HEIGHT);
	Poly_Exp_hls_strm(src_img_orig_strm,src_poly,WIDTH,HEIGHT);
	Poly_Exp_hls_strm(dst_img_orig_strm,dst_poly,WIDTH,HEIGHT);
	UpdateMat_0_hls(src_poly, dst_poly, M, WIDTH,HEIGHT);
	UpdateFlow_hls(M, flow, WIDTH, HEIGHT);

	for(int i=0; i<MAXSIZE; i++){
#pragma HLS PIPELINE II=5
		Data_2 tmp = flow.read();
		mig_out[i*2] = tmp.r0;
		mig_out[i*2+1] = tmp.r1;
	}
}

void Farneback_core(hls::stream<pix_t> &src_img, pix_t dst_img[POLY_EXP_SAMPLE_SIZE*POLY_EXP_SAMPLE_SIZE][RAM_LENGTH],
		hls::stream<Data_2u>& flow_in, hls::stream<Data_2>& flow_out, const pix_t con){
#pragma HLS ARRAY_PARTITION variable=dst_img complete dim=1
#pragma HLS DATAFLOW
	//hls::stream<pix_t> src_img_strm("src_img_strm");
	hls::stream<Data_5> src_poly("src_poly");
	hls::stream<Data_5> dst_poly("dst_poly");
	hls::stream<Data_2u> flow_in1("flow_in1");
	hls::stream<Data_2u> flow_in2("flow_in2");
#pragma HLS STREAM variable=flow_in2 depth=192 dim=1
	hls::stream<Data_5> M("M");

	int scale = 1 << ((con&0x3f)-1);
	int height, width;
	height = HEIGHT / scale;
	width = WIDTH / scale;

	for(int i=0;i<width*height;i++){
#pragma HLS LOOP_TRIPCOUNT min=76800 max=307200
		Data_2u in_val = flow_in.read();
		flow_in1.write(in_val);
		flow_in2.write(in_val);
	}

	//Smooth_hls(src_img, src_img_strm, width, height);
	Poly_Exp_hls_strm(src_img, src_poly, width, height);
	Poly_Exp_hls_bram(dst_img, flow_in1, dst_poly, width, height);
	UpdateMat_hls(src_poly, dst_poly, flow_in2, M, width,height);
	UpdateFlow_hls(M, flow_out, width, height);

}

void Farneback_core_wrapper(pix_t *src_in, pix_t dst_in[POLY_EXP_SAMPLE_SIZE*POLY_EXP_SAMPLE_SIZE][RAM_LENGTH],
		Data_2u flow_in[MAXSIZE], Data_2u* flow_out, const pix_t con){
#pragma HLS ARRAY_PARTITION variable=dst_in complete dim=1
#pragma HLS DATAFLOW
	//control bit sample 1bit|useflow 1bit|scale 6bit
	int scale = 1 << ((con&0x3f)-1);
	int height, width;
	bool sampleFlow = (con&0x80) != 0;
	height = HEIGHT / scale;
	width = WIDTH / scale;
	int src_base_addr = (scale * scale - 1) * 8 / 3 * MAXSIZE / scale / scale;
	int br = (width - 1) / POLY_EXP_SAMPLE_SIZE + 1;

	hls::stream<pix_t> src_img("src_img");
	hls::stream<Data_2u> core_flow_in("core_flow_in");
	hls::stream<Data_2> core_flow_out("core_flow_out");

	IMG_src_rd_loop:for(int i=0; i<height; i++){
		for(int j=0; j<width; j++){
#pragma HLS PIPELINE
			pix_t img_val = src_in[src_base_addr++];
			src_img.write(img_val);
		}
	}
	Flow_src_rd_loop:for(int i=0; i<height; i++){
		for(int j=0; j<width; j++){
#pragma HLS PIPELINE
			Data_2u flow_val;
			int flow_addr;
			if(sampleFlow)
				flow_addr = (i&0xffe) * scale * WIDTH + (j&0xffe) * scale;
			else
				flow_addr = i * scale * WIDTH + j * scale;
			flow_val = flow_in[flow_addr];
			flow_val.r0 = flow_val.r0/scale;
			flow_val.r1 = flow_val.r1/scale;
			core_flow_in.write(flow_val);
		}
	}

	Farneback_core(src_img, dst_in, core_flow_in, core_flow_out, con);
	Out_loop:for(int i=0;i<height;i++){
		for(int j=0;j<width;j++){
#pragma HLS PIPELINE II=2
			Data_2 out_val;
			Data_2u res;
			out_val = core_flow_out.read();
			res.r0 = (char)(out_val.r0 * scale * 2)/2 + (char)(out_val.r0 * scale * 2)%2;
			res.r1 = (char)(out_val.r1 * scale * 2)/2 + (char)(out_val.r1 * scale * 2)%2;
			flow_out[i*scale*WIDTH + j*scale] = res;
		}
	}
}

void Farneback_top_v2_1(pix_t* mig_in, Data_2u* mig_out, const pix_t con){
#pragma HLS DATA_PACK variable=mig_out
#pragma HLS INTERFACE ap_vld port=con
#pragma HLS INTERFACE s_axilite port=con bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite port=return bundle=CTRL_BUS
#pragma HLS INTERFACE ap_ctrl_hs port=return
//#pragma HLS INTERFACE m_axi depth=576 port=mig_in
//#pragma HLS INTERFACE m_axi depth=192 port=mig_out
#pragma HLS INTERFACE m_axi depth=960000 port=mig_in
#pragma HLS INTERFACE m_axi depth=640000 port=mig_out
	Data_2u flow_in[MAXSIZE];
	pix_t dst_img[POLY_EXP_SAMPLE_SIZE*POLY_EXP_SAMPLE_SIZE][RAM_LENGTH];
#pragma HLS ARRAY_PARTITION variable=dst_img complete dim=1
	/* ddr3 img memory
	 * 640 * 480 src_img
	 * 640 * 480 dst_img
	 * 320 * 240 src_img_2
	 * 320 * 240 dst_img_2
	 * */
	int scale = 1 << ((con&0x3f)-1);
	int height, width;
	bool useFlow = (con&0x40) != 0;
	bool sampleFlow = (con&0x80) != 0;
	height = HEIGHT / scale;
	width = WIDTH / scale;
	int dst_base_addr = ((scale * scale - 1) * 8 / 3 + 1)* MAXSIZE / scale / scale;
	int br = (width - 1) / POLY_EXP_SAMPLE_SIZE + 1;
	IMG_dst_rd_loop:for(int i=0; i<height; i++){
		for(int j=0; j<width; j++){
#pragma HLS PIPELINE
			int sel = i % POLY_EXP_SAMPLE_SIZE * POLY_EXP_SAMPLE_SIZE + j % POLY_EXP_SAMPLE_SIZE;
			int addr = i / POLY_EXP_SAMPLE_SIZE * br + j / POLY_EXP_SAMPLE_SIZE;
			dst_img[sel][addr] = mig_in[dst_base_addr++];
			int flow_addr;
			if(sampleFlow)
				flow_addr = (i&0xffe) * scale * WIDTH + (j&0xffe) * scale;
			else
				flow_addr = i * scale * WIDTH + j * scale;
			if(useFlow){
				Data_2u val = mig_out[flow_addr];
				flow_in[flow_addr] = val;
			}
			else{
				Data_2u val = {0,0};
				flow_in[flow_addr] = val;
			}
		}
	}

	Farneback_core_wrapper(mig_in, dst_img, flow_in, mig_out, con);
//	void Farneback_core_wrapper(pix_t *src_in, pix_t dst_in[POLY_EXP_SAMPLE_SIZE*POLY_EXP_SAMPLE_SIZE][RAM_LENGTH],
//			Data_2u flow_in[MAXSIZE], Data_2u* flow_out, pix_t con)
}

void Farneback_top_v2(volatile pix_t* mig_in_src, volatile pix_t* mig_in_dst,volatile data_t* mig_out, unsigned int con){
#pragma HLS INTERFACE ap_vld port=con
#pragma HLS INTERFACE s_axilite port=con bundle=CTRL_BUS
#pragma HLS INTERFACE s_axilite port=return bundle=CTRL_BUS
#pragma HLS INTERFACE ap_ctrl_hs port=return
#pragma HLS INTERFACE m_axi depth=640000 port=mig_in_src
#pragma HLS INTERFACE m_axi depth=640000 port=mig_in_dst
#pragma HLS INTERFACE m_axi depth=3200000 port=mig_out
	static Data_2u flow[HEIGHT * WIDTH];
	pix_t dst_img[POLY_EXP_SAMPLE_SIZE*POLY_EXP_SAMPLE_SIZE][RAM_LENGTH];
#pragma HLS ARRAY_PARTITION variable=dst_img complete dim=1
	/* ddr3 img memory
	 * 640 * 480 src_img
	 * 640 * 480 dst_img
	 * 320 * 240 src_img_2
	 * 320 * 240 dst_img_2
	 * */
	int scale = 1 << ((con & 7) - 1); //scale: 2 power integer
	int height, width;
	bool useFlow = (con & 8) != 0;
	bool smallFlow = (con & 16) != 0;
	bool outFlow = (con & 32) != 0;
	height = HEIGHT / scale;
	width = WIDTH / scale;
	int src_base_addr = (scale * scale - 1) * 8 / 3 * MAXSIZE / scale / scale;
	int dst_base_addr = src_base_addr+ MAXSIZE/scale/scale;
	int br = (width - 1) / POLY_EXP_SAMPLE_SIZE + 1;

	IMG_dst_rd_loop:for(int i=0; i<height; i++){
		for(int j=0; j<width; j++){
#pragma HLS PIPELINE
			int sel = i % POLY_EXP_SAMPLE_SIZE * POLY_EXP_SAMPLE_SIZE + j % POLY_EXP_SAMPLE_SIZE;
			int addr = i / POLY_EXP_SAMPLE_SIZE * br + j / POLY_EXP_SAMPLE_SIZE;
			dst_img[sel][addr] = mig_in_dst[dst_base_addr++];
		}
	}

	Core:{
#pragma HLS DATAFLOW
	hls::stream<pix_t> src_img("src_img");
	hls::stream<Data_2u> flow_in("flow_in");
	hls::stream<Data_2> flow_out("flow_out");
#pragma HLS DEPENDENCE variable=flow inter false
	IMG_src_rd_loop:for(int i=0; i<height; i++){
		for(int j=0; j<width; j++){
#pragma HLS PIPELINE
			pix_t img_val = mig_in_src[src_base_addr++];
			src_img.write(img_val);
			Data_2u flow_val;
			if(useFlow){
				int flow_addr;
				if(smallFlow)
					flow_addr = (i&0xffe) * scale * WIDTH + (j&0xffe) * scale;
				else
					flow_addr = i * scale * WIDTH + j * scale;
				flow_val = flow[flow_addr];
				flow_val.r0 = flow_val.r0/scale;
				flow_val.r1 = flow_val.r1/scale;
			}
			else{
				flow_val.r0 = 0;
				flow_val.r1 = 0;
			}
			flow_in.write(flow_val);
		}
	}


	Farneback_core(src_img, dst_img, flow_in, flow_out, con);
	int flow_out_base = 0;
	Out_loop:for(int i=0;i<height;i++){
		for(int j=0;j<width;j++){
#pragma HLS PIPELINE II=5
			Data_2 out_val;
			out_val = flow_out.read();
			if(outFlow){
				mig_out[flow_out_base] = out_val.r0;
				mig_out[flow_out_base+1] = out_val.r1;
				flow_out_base += 2;
			}
			else{
				Data_2u out_buf;
				out_buf.r0 = (pix_t)(out_val.r0 * scale * 2)/2 + (pix_t)(out_val.r0 * scale * 2)%2;
				out_buf.r1 = (pix_t)(out_val.r1 * scale * 2)/2 + (pix_t)(out_val.r1 * scale * 2)%2;
				flow[i*scale*WIDTH + j*scale] = out_buf;
			}
		}
	}
	}//Core
}
/*

void UpdataMat_2_1_hls(hls::stream<Data_5> &src_poly, Data_5 dst_poly[MAXSIZE], hls::stream<Data_2>& flow_in,
		hls::stream<Data_5>& M, int width, int height){
	Data_2 flow_buf[WIDTH], data2_tmp;
	UpdataMat_2_1_hls_label2:for(int i=0;i<height;i++){
#pragma HLS LOOP_TRIPCOUNT min=240 max=480
		UpdataMat_2_1_hls_label1:for(int j=0;j<width;j++){
#pragma HLS LOOP_TRIPCOUNT min=320 max=640
#pragma HLS PIPELINE
			int fx, fy, dx, dy;
			data_t dx_f, dy_f;
			if(i%2 == 0 && j%2 ==0){
				data2_tmp = flow_in.read();
				flow_buf[j/2] = data2_tmp;
			}
			else
				data2_tmp = flow_buf[j/2];
			dx_f = data2_tmp.r0 * 2;
			dy_f = data2_tmp.r1 * 2;
			dx = (dx_f >= 0) ? (int)(dx_f + 0.5) : (int)(dx_f - 0.5);
			dy = (dy_f >= 0) ? (int)(dy_f + 0.5) : (int)(dy_f - 0.5);
			fx = i + dx;
			fy = j + dy;
			fx = (fx < 0) ? 0 : (fx >= height) ? height - 1 : fx;
			fy = (fy < 0) ? 0 : (fy >= width) ? width - 1 : fy;
			Data_5 src, dst, m;
			data_t a00, a01, a11, b0, b1;
			src = src_poly.read();
			dst = dst_poly[fx*width + fy];
			a00 = (src.r2 + dst.r2) / 2; //r4
			a01 = (src.r4 + dst.r4) / 4; //r6
			a11 = (src.r3 + dst.r3) / 2; //r5
			b0 = (src.r0 - dst.r0) / 2; //r2
			b1 = (src.r1 - dst.r1) / 2; //r3

			b0 += a00 * dx + a01 * dy;
			b1 += a01 * dx + a11 * dy;

			m.r0 = a00*a00 + a01*a01; // G(0, 0)
			m.r1 = a01*(a00 + a11); // G(0, 1)
			m.r2 = a11*a11 + a01*a01; // G(1, 1)
			m.r3 = a00*b0 + a01*b1; // H(0)
			m.r4 = a01*b0 + a11*b1; // H(1)

			M.write(m);
		}
	}
}
*/
