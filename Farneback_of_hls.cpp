#include "Farneback_of_hls.h"

void Smooth_hls(hls::stream<pix_t>& in, hls::stream<pix_t>&out, int width, int height){
#pragma HLS DATAFLOW
	const int K = 5;
	const unsigned char coeff[5] = {1, 4, 6, 4, 1};
	unsigned short hwin[5];
	hls::stream<unsigned short> hconv("hconv");
	static unsigned short linebuf[4][WIDTH];
#pragma HLS ARRAY_PARTITION variable=linebuf dim=1 complete
	hls::stream<unsigned short> vconv("vconv");
	const int vconv_xlim = width + 1 - K;

	// Horizontal convolution
	HConvH:for(int i=0;i<height;i++){
#pragma HLS LOOP_TRIPCOUNT min=240 max=480
		HConvW:for(int j=0;j<width;j++){
#pragma HLS LOOP_TRIPCOUNT min=320 max=640
#pragma HLS PIPELINE II=2
			unsigned short in_val = in.read();
			unsigned short out_val = 0;
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
#pragma HLS PIPELINE II=2
			unsigned short in_val = hconv.read();
			unsigned short out_val = 0;
			Vconv:for(int k=0; k<K; k++){
				unsigned short vwin_val = k < K-1? linebuf[k][j]: in_val;
				out_val = out_val + vwin_val * coeff[k];
				if(k > 0)
					linebuf[k-1][j] = vwin_val;
			}
			if(i >= K-1)
				vconv.write(out_val);
		}
	}

	const int border_width = int(K/2);
	unsigned short borderbuf[WIDTH + 1 -K];

	BorderH:for(int i=0; i<height; i++){
#pragma HLS LOOP_TRIPCOUNT min=240 max=480
		BorderW:for(int j=0; j<width; j++){
#pragma HLS PIPELINE II=2
#pragma HLS LOOP_TRIPCOUNT min=320 max=640
			unsigned short d_in, l_edge, r_edge, d_out;
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

void Poly_Exp_hls_strm(hls::stream<pix_t>& in, hls::stream<Data_5> &out1, hls::stream<Data_5> &out2, int width, int height){
#pragma HLS DATAFLOW
	//initialize the matrix
	int n = (POLY_EXP_SAMPLE_SIZE - 1) / 2;
	data_t coeff[5][POLY_EXP_SAMPLE_SIZE][POLY_EXP_SAMPLE_SIZE];
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

	pix_t hwin[POLY_EXP_SAMPLE_SIZE];
	hls::stream<Data_3> hconv("hconv");

	Data_3 linebuf[POLY_EXP_SAMPLE_SIZE-1][WIDTH];
#pragma HLS ARRAY_PARTITION variable=linebuf dim=1 complete
	hls::stream<Data_5> vconv("vconv");
	Data_5 borderbuf[WIDTH + 1 -POLY_EXP_SAMPLE_SIZE];

	HConvH:for(int i=0;i<height;i++){
#pragma HLS LOOP_TRIPCOUNT min=240 max=480
		HConvW:for(int j=0;j<width;j++){
#pragma HLS PIPELINE II=2
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
#pragma HLS PIPELINE II=2
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
#pragma HLS PIPELINE II=2
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
			out1.write(d_out);
			out2.write(d_out);
		}
	}
}

void UpdataMat_2_1_hls(hls::stream<Data_5> &src_poly, Data_5 dst_poly[MAXSIZE], hls::stream<Data_2>& flow_in,
		hls::stream<Data_5>& M, short width, short height){
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

void UpdateMat_0_hls(hls::stream<Data_5> &src_poly, hls::stream<Data_5> &dst_poly, hls::stream<Data_5>& M,
		short width, short height){
	UpdataMat_0H:for(int i=0;i<height;i++){
#pragma HLS LOOP_TRIPCOUNT min=240 max=480
		UpdataMat_0W:for(int j=0;j<width;j++){
#pragma HLS LOOP_TRIPCOUNT min=320 max=640
#pragma HLS PIPELINE II=2
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
#pragma HLS PIPELINE II=2
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
#pragma HLS PIPELINE II=2
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
#pragma HLS PIPELINE II=2
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

void Farneback_core(hls::stream<pix_t>& img_in, hls::stream<Data_5>& d5_in, hls::stream<Data_5>& d5_out, hls::stream<Data_2>& d2_out, short width, short height){
#pragma HLS DATAFLOW
	hls::stream<Data_5> poly("poly");
	hls::stream<Data_5> M("M");
	Poly_Exp_hls_strm(img_in, poly, d5_out, width, height);
	UpdateMat_0_hls(poly, d5_in, M, width, height);
	UpdateFlow_hls(M, d2_out, width, height);
}


void Farneback_top(volatile pix_t* mig_in, volatile data_t* mig_out){
#pragma HLS INTERFACE s_axilite port=return
#pragma HLS INTERFACE ap_ctrl_hs port=return
#pragma HLS INTERFACE m_axi depth=640000 port=mig_in
#pragma HLS INTERFACE m_axi depth=3200000 port=mig_out
#pragma HLS DATAFLOW
	hls::stream<pix_t> src_img_strm("src_img_strm");
#pragma HLS DATA_PACK variable=src_img_strm
	hls::stream<pix_t> dst_img_strm("dst_img_strm");
#pragma HLS DATA_PACK variable=dst_img_strm
	hls::stream<Data_5> src_poly("src_poly");
#pragma HLS DATA_PACK variable=src_poly
	hls::stream<Data_5> dst_poly("dst_poly");
#pragma HLS DATA_PACK variable=dst_poly
	hls::stream<Data_5> M("M");
#pragma HLS DATA_PACK variable=M
	hls::stream<Data_2> flow("flow");
#pragma HLS DATA_PACK variable=flow

	hls::stream<Data_5> src_poly2("src_poly2");
	hls::stream<Data_5> dst_poly2("dst_poly2");

	for(int i=0;i<MAXSIZE;i++){
#pragma HLS PIPELINE II=2
		pix_t tmp = mig_in[i*2];
		pix_t tmp2 = mig_in[i*2+1];
		src_img_strm << tmp;
		dst_img_strm << tmp2;
	}

	Poly_Exp_hls_strm(src_img_strm,src_poly,src_poly2,WIDTH,HEIGHT);
	Poly_Exp_hls_strm(dst_img_strm,dst_poly,dst_poly2,WIDTH,HEIGHT);
	UpdateMat_0_hls(src_poly, dst_poly, M, WIDTH,HEIGHT);
	UpdateFlow_hls(M, flow, WIDTH, HEIGHT);

	for(int i=0; i<MAXSIZE; i++){
#pragma HLS PIPELINE II=2
		src_poly2.read();
		dst_poly2.read();
		Data_2 tmp = flow.read();
		mig_out[i*2] = tmp.r0;
		mig_out[i*2+1] = tmp.r1;
	}
}


