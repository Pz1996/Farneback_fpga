#pragma once
#include "define.h"
#include "hls_stream.h"
#include <math.h>

struct Data_3{
	data_t r, xr, xxr;
};

struct Data_5{
	data_t r0, r1, r2, r3, r4;
	Data_5 operator+(const Data_5& in){
		Data_5 out;
		out.r0 = in.r0 + r0;
		out.r1 = in.r1 + r1;
		out.r2 = in.r2 + r2;
		out.r3 = in.r3 + r3;
		out.r4 = in.r4 + r4;
		return out;
	}

	Data_5 operator/(pix_t coeff){
		Data_5 out;
		out.r0 = r0 / coeff;
		out.r1 = r1 / coeff;
		out.r2 = r2 / coeff;
		out.r3 = r3 / coeff;
		out.r4 = r4 / coeff;
		return out;
	}
};

struct Data_2{
	data_t r0, r1;
	Data_2 operator+(const Data_2& in){
		Data_2 out;
		out.r0 = in.r0 + r0;
		out.r1 = in.r1 + r1;
		return out;
	}

	Data_2 operator/(pix_t coeff){
		Data_2 out;
		out.r0 = r0 / coeff;
		out.r1 = r1 / coeff;
		return out;
	}
};

void Poly_Exp_hls_strm(hls::stream<pix_t>& in, hls::stream<Data_5> &out, int width, int height);

void UpdataMat_2_1_hls(hls::stream<Data_5> &src_poly, Data_5 dst_poly[MAXSIZE], hls::stream<Data_2>& flow_in,
		hls::stream<Data_5>& M, short width, short height);

void UpdateMat_0_hls(hls::stream<Data_5> &src_poly, hls::stream<Data_5> &dst_poly, hls::stream<Data_5>& M,
		short width, short height);

void UpdateFlow_hls(hls::stream<Data_5>&M, hls::stream<Data_2>&flow_out, int width, int height);

void Farneback_top(volatile pix_t* mig_in, volatile data_t* mig_out);
