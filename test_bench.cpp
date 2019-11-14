#include "Farneback_of.h"
#include "Farneback_of_hls.h"
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


}

int main(){
	return test_mat_2_1();
}
