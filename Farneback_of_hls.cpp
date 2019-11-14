#include "Farneback_of_hls.h"

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
			if((unsigned)(i - K/2) < height-K && (unsigned)(j - K/2) < width-K){
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
