 #include "Farneback_of.h"
#include "Farneback_of_hls.h"
#include <stdio.h>
#include <iostream>
using namespace std;

int test_mat_2_1(){
	data_t **src_poly;//[MAXSIZE][5];
	data_t **dst_poly;//[MAXSIZE][5];
	data_t **flow_in;//[MAXSIZE][2];
	data_t **M;//[MAXSIZE][5];

	int scale = 2;
	hls::stream<Data_5> hls_src_poly("hls_src_poly");
	Data_5* hls_dst_poly;//[MAXSIZE];
	hls::stream<Data_2> hls_flow_in("hls_flow_in");
	hls::stream<Data_5>hls_M("hls_M");

	src_poly = new data_t*[MAXSIZE];
	dst_poly = new data_t*[MAXSIZE];
	flow_in = new data_t*[MAXSIZE];
	M = new data_t*[MAXSIZE];
	hls_dst_poly = new Data_5[MAXSIZE];

	for(int i=0;i<MAXSIZE;i++){
		src_poly[i] = new data_t[5];
		dst_poly[i] = new data_t[5];
		flow_in[i] = new data_t[2];
		M[i] = new data_t[5];
	}


	//Generate input data
	for(int i=0;i<WIDTH*HEIGHT;i++){
		Data_5 d5;
		src_poly[i][0] = d5.r0 = i%255;
		src_poly[i][1] = d5.r1 = i%233;
		src_poly[i][2] = d5.r2 = i%244;
		src_poly[i][3] = d5.r3 = i%222;
		src_poly[i][4] = d5.r4 = i%211;
		hls_src_poly.write(d5);
		dst_poly[i][0] = hls_dst_poly[i].r0 = i%155;
		dst_poly[i][1] = hls_dst_poly[i].r1 = i%133;
		dst_poly[i][2] = hls_dst_poly[i].r2 = i%144;
		dst_poly[i][3] = hls_dst_poly[i].r3 = i%122;
		dst_poly[i][4] = hls_dst_poly[i].r4 = i%111;
	}
	for(int i=0;i<WIDTH*HEIGHT/4; i++){
		Data_2 d2;
		flow_in[i][0] = d2.r0 = i%2;
		flow_in[i][1] = d2.r1 = i%3;
		hls_flow_in.write(d2);
	}

	//ref model
	UpdateMat(src_poly, dst_poly, flow_in, M, WIDTH, HEIGHT, 2);
	//DUT
	UpdataMat_2_1_hls(hls_src_poly, hls_dst_poly, hls_flow_in, hls_M, WIDTH, HEIGHT);

	//check
	int err_cnt = 0;
	for(int i=0;i<WIDTH*HEIGHT;i++){
		Data_5 d5;
		d5 = hls_M.read();
		if(d5.r0 != M[i][0]) err_cnt++;
		if(d5.r1 != M[i][1]) err_cnt++;
		if(d5.r2 != M[i][2]) err_cnt++;
		if(d5.r3 != M[i][3]) err_cnt++;
		if(d5.r4 != M[i][4]) err_cnt++;
	}

	if(err_cnt == 0)
		cout << "*** TEST PASSED ***" << endl;
	else
		cout << err_cnt << " errors are detected!\n"<< "*** TEST FAILED ***" << endl;
	return (err_cnt == 0)? 0 : -1;
}

int test_flow(){
	data_t **M;
	data_t **flow_out;
	hls::stream<Data_5> hls_M("hls_M");
	hls::stream<Data_2> hls_flow_out("hls_flow_out");
	M = new data_t*[MAXSIZE];
	flow_out = new data_t*[MAXSIZE];
	for(int i=0; i<MAXSIZE; i++){
		M[i] = new data_t[5];
		flow_out[i] = new data_t[2];
	}

	//Generate input data
	for(int i=0;i<WIDTH*HEIGHT;i++){
		Data_5 d5;
		M[i][0] = d5.r0 = i%255;
		M[i][1] = d5.r1 = i%233;
		M[i][2] = d5.r2 = i%244;
		M[i][3] = d5.r3 = i%222;
		M[i][4] = d5.r4 = i%211;
		hls_M.write(d5);

	}

	//ref model
	UpdateFlow(M, flow_out, WIDTH, HEIGHT);
	//dut model
	UpdateFlow_hls(hls_M, hls_flow_out, WIDTH, HEIGHT);

	//check
	int err_cnt = 0;
	for(int i=0;i<WIDTH*HEIGHT;i++){
		Data_2 d2;
		d2 = hls_flow_out.read();
		if((d2.r0 - flow_out[i][0])*(d2.r1 - flow_out[i][1]) > 1 || (d2.r0 - flow_out[i][0])*(d2.r1 - flow_out[i][1]) < -1){
			err_cnt++;
			printf("%d:(%f,%f) (%f,%f)\n",i, d2.r0, d2.r1, flow_out[i][0],flow_out[i][1]);
		}
	}

	if(err_cnt == 0)
		cout << "*** TEST PASSED ***" << endl;
	else
		cout << err_cnt << " errors are detected!\n"<< "*** TEST FAILED ***" << endl;
	return (err_cnt == 0)? 0 : -1;
}

int test_poly(){
	pix_t *in;
	data_t **out;
	in = new pix_t[MAXSIZE];
	out = new data_t*[MAXSIZE];
	for(int i=0;i<MAXSIZE;i++)
		out[i] = new data_t[5];
	hls::stream<pix_t> in_hls("in_hls");
	hls::stream<Data_5> out_hls("out_hls");

	for(int i=0;i<MAXSIZE;i++){
		in[i] = i % 234;
		in_hls.write(i % 234);
	}

	Poly_Exp(in, out, WIDTH, HEIGHT);
	Poly_Exp_hls_strm(in_hls, out_hls, WIDTH, HEIGHT);

	//check
	int err_cnt = 0;
	for(int i=0;i<MAXSIZE;i++){
		Data_5 d5;
		d5 = out_hls.read();
		data_t sq_sum = 0;
		sq_sum += (out[i][0] - d5.r0)*(out[i][0] - d5.r0);
		sq_sum += (out[i][1] - d5.r1)*(out[i][1] - d5.r1);
		sq_sum += (out[i][2] - d5.r2)*(out[i][2] - d5.r2);
		sq_sum += (out[i][3] - d5.r3)*(out[i][3] - d5.r3);
		sq_sum += (out[i][4] - d5.r4)*(out[i][4] - d5.r4);

		if(sq_sum >= 1){
			cout << out[i][0] << "\t"<< out[i][1] << "\t"<< out[i][2] << "\t"<< out[i][3] << "\t"<< out[i][4] << endl;
			cout << d5.r0 << "\t"<< d5.r1 << "\t"<<d5.r2 << "\t"<<d5.r3<< "\t"<< d5.r4 << endl << endl;
			err_cnt++;
		}
	}
	if(err_cnt == 0)
		cout << "*** TEST PASSED! ***" << endl;
	else
		cout << err_cnt << " errors are detected!\n"<< "*** TEST FAILED! ***" << endl;
	return (err_cnt == 0)? 0 : -1;


	for(int i=0;i<MAXSIZE;i++)
		delete[]out[i];
	delete[] in;
	delete[] out;
}

int main(){
	return test_poly();
}
