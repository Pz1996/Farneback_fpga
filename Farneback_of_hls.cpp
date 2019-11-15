#include "Farneback_of_hls.h"

void Poly_Exp_hls_strm(hls::stream<pix_t>& in, hls::stream<Data_5> &out, int width, int height){

}

/*
Data_5 Cal_poly_kernel(data_t in[POLY_EXP_SAMPLE_SIZE * POLY_EXP_SAMPLE_SIZE * 5]){
#pragma HLS ARRAY_PARTITION variable=in complete dim=1
#pragma HLS ARRAY_RESHAPE variable=in cyclic factor=5 dim=1
#pragma HLS PIPELINE
	Data_5 out_val = {0,0,0,0,0};
	int n = (POLY_EXP_SAMPLE_SIZE - 1) / 2;
	data_t coeff[POLY_EXP_SAMPLE_SIZE * POLY_EXP_SAMPLE_SIZE * 5];
#pragma HLS ARRAY_PARTITION variable=coeff complete dim=1
#pragma HLS ARRAY_RESHAPE variable=coeff cyclic factor=5 dim=1
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
				int base_addr = ((i+n) * POLY_EXP_SAMPLE_SIZE +(j+n))* 5;
				coeff[base_addr] = ig11 * b2;
				coeff[base_addr + 1] = ig11 * b3;
				coeff[base_addr + 2] = ig33 * b4 + ig03 * b1;
				coeff[base_addr + 3] = ig33 * b5 + ig03 * b1;
				coeff[base_addr + 4] = ig55 * b6;
			}
		}

	for(int i=0;i<POLY_EXP_SAMPLE_SIZE * POLY_EXP_SAMPLE_SIZE; i++){
		out_val.r0 += coeff[i * 5 + 0] * in[i * 5 + 0];
		out_val.r1 += coeff[i * 5 + 1] * in[i * 5 + 1];
		out_val.r2 += coeff[i * 5 + 2] * in[i * 5 + 2];
		out_val.r3 += coeff[i * 5 + 3] * in[i * 5 + 3];
		out_val.r4 += coeff[i * 5 + 4] * in[i * 5 + 4];
	}
	return out_val;


}




void Poly_Exp_hls_strm(hls::stream<pix_t>& in, hls::stream<Data_5> &out, int width, int height){
	//initialize the matrix
	int n = (POLY_EXP_SAMPLE_SIZE - 1) / 2;
	data_t coeff[5][POLY_EXP_SAMPLE_SIZE][POLY_EXP_SAMPLE_SIZE];
#pragma HLS ARRAY_PARTITION variable=coeff complete dim=0
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

	pix_t window[POLY_EXP_SAMPLE_SIZE][POLY_EXP_SAMPLE_SIZE];
#pragma HLS ARRAY_PARTITION variable=window complete dim=0
	pix_t buf[POLY_EXP_SAMPLE_SIZE-1][WIDTH];
	for(int i=0;i<height;i++){
		for(int j=0;i<width;j++){
#pragma HLS PIPELINE
			pix_t in_val = in.read();
			if(i == 0){
				for(int k=0; k<POLY_EXP_SAMPLE_SIZE-1; k++){
					buf[k][0] = in_val;
				}
			}
			if(j == 0){
				for(int ii=0;ii<POLY_EXP_SAMPLE_SIZE - 1;ii++){
					for(int jj=0;jj<POLY_EXP_SAMPLE_SIZE;jj++){
						window[ii][jj] = buf[ii][0];
					}
				}
				for(int k=0; k<POLY_EXP_SAMPLE_SIZE; k++){
					window[POLY_EXP_SAMPLE_SIZE - 1][k] = in_val;
				}
			}
			else{
				for(int ii=0;ii<POLY_EXP_SAMPLE_SIZE;ii++){
					for(int jj=0;jj<POLY_EXP_SAMPLE_SIZE - 1;jj++){
						window[ii][jj] = window[ii][jj+1];
					}
				}
				for(int k=0; k<POLY_EXP_SAMPLE_SIZE-1; k++)
					window[k][POLY_EXP_SAMPLE_SIZE - 1] = buf[k][j];
				window[POLY_EXP_SAMPLE_SIZE - 1][POLY_EXP_SAMPLE_SIZE - 1] = in_val;
			}
			for(int k=0; k<POLY_EXP_SAMPLE_SIZE-1; k++)
				buf[k][j] = buf[k+1][j];
			buf[POLY_EXP_SAMPLE_SIZE-1][j] = in_val;
			data_t r[5];
#pragma HLS ARRAY_PARTITION variable=r complete dim=1
			for(int k=0;k<5;k++)
				r[k] = 0;
			for(int ii=0;ii<POLY_EXP_SAMPLE_SIZE;ii++){
				for(int jj=0;jj<POLY_EXP_SAMPLE_SIZE;jj++){
					for(int k=0;k<5;k++){
						r[k]+=window[ii][jj] * coeff[k][ii][jj];
					}
				}
			}
			Data_5 out_val = {r[0],r[1],r[2],r[3],r[4]};
			out.write(out_val);
		}
	}
}
*/

void UpdataMat_2_1_hls(hls::stream<Data_5> &src_poly, Data_5 dst_poly[MAXSIZE], hls::stream<Data_2>& flow_in,
		hls::stream<Data_5>& M, short width, short height){
	Data_2 flow_buf[WIDTH], data2_tmp;
	UpdataMat_2_1_hls_label2:for(int i=0;i<height;i++){
		UpdataMat_2_1_hls_label1:for(int j=0;j<width;j++){
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

void UpdateFlow_hls(hls::stream<Data_5>&M, hls::stream<Data_2>&flow_out, int width, int height){
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
#pragma HLS PIPELINE
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
#pragma HLS PIPELINE
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
#pragma HLS PIPELINE
			Data_2 out = {0, 0};
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
