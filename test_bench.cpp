#include "Farneback_of.h"
#include "Farneback_of_hls.h"
#include <stdio.h>
#include <iostream>
using namespace std;

int test_smooth(){
	int height=48, width = 64;
	pix_t *in, *out;
	hls::stream<pix_t> in_hls("in_hls");
	hls::stream<pix_t> out_hls("out_hls");
	in = new pix_t[height * width];
	out = new pix_t[height * width];

	//Generate input data
	for(int i=0;i<height*width;i++){
		in[i] = i%123 + i%23;
		in_hls.write(i%123 + i%23);
	}

	Smooth_hls(in_hls, out_hls, width, height);
	Smooth(in, out, width, height);

	//check
	int err_cnt = 0;
	for(int i=0;i<width*height;i++){
		pix_t in_val = out_hls.read();
		if(in_val != out[i]) err_cnt++;

	}

	if(err_cnt == 0)
		cout << "*** SMOOTH TEST PASSED ***" << endl;
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

int test_top_v2_case0(){
	//void Farneback_top_v2_1(pix_t* mig_in, Data_2u* mig_out, pix_t con)
	pix_t* hls_in;
	pix_t* in;
	Data_2u* hls_out;
	in = new pix_t[MAXSIZE * 2];
	hls_in = new pix_t[MAXSIZE * 2];
	hls_out = new Data_2u[MAXSIZE];
	data_t **src_poly, **dst_poly, **M, **flow_out, **flow_in;
	src_poly = new data_t*[MAXSIZE];
	dst_poly = new data_t*[MAXSIZE];
	M = new data_t*[MAXSIZE];
	flow_in = new data_t*[MAXSIZE];
	flow_out = new data_t*[MAXSIZE];
	for(int i=0;i<MAXSIZE;i++){
		src_poly[i] = new data_t[5];
		dst_poly[i] = new data_t[5];
		M[i] = new data_t[5];
		flow_out[i] = new data_t[2];
		flow_in[i] = new data_t[2];
	}
	for(int i=0;i<MAXSIZE;i++){
		flow_in[i][0] = hls_out[i].r0 = (i*i)%5 - 3;
		flow_in[i][1] = hls_out[i].r1 = (i*i)%7 - 4;
	}
	//source code
	for(int i=0;i<MAXSIZE;i++){
		hls_in[i] = in[i] = (i*i) % 179;
		hls_in[i + MAXSIZE] = in[i + MAXSIZE] = (i*i) %123;
	}
	pix_t con = 0x41;
	//dut
	Farneback_top_v2_1(hls_in, hls_out, con);
	//ref model
	Poly_Exp(in, src_poly, WIDTH, HEIGHT);
	Poly_Exp(in+MAXSIZE, dst_poly, WIDTH, HEIGHT);
	UpdateMat(src_poly, dst_poly, flow_in, M, WIDTH, HEIGHT, 1);
	UpdateFlow(M, flow_out, WIDTH, HEIGHT);
	//check
	int err_cnt = 0;
	for(int i=0;i<MAXSIZE;i++){
		Data_2u out_val = hls_out[i];
		data_t score = (flow_out[i][0] - out_val.r0) * (flow_out[i][0] - out_val.r0) +
				(flow_out[i][1] - out_val.r1) * (flow_out[i][1] - out_val.r1);
		if(score > 1){
			err_cnt++;
			printf("[%d,%d]: (%f,%f)(%d,%d)\n",i/WIDTH, i%WIDTH, flow_out[i][0], flow_out[i][1], (int)out_val.r0, (int)out_val.r1);
		}
	}
	if(err_cnt == 0)
		cout << "*** V2 TOP TEST PASSED! ***" << endl;
	else
		cout << err_cnt << " errors are detected!\n"<< "*** TEST FAILED! ***" << endl;

	for(int i=0;i<MAXSIZE;i++){
		delete[]src_poly[i];
		delete[]dst_poly[i];
		delete[]M[i];
		delete[]flow_out[i];
	}

	delete[] hls_in;
	delete[] in;
	delete[] hls_out;
	delete[] src_poly;
	delete[] dst_poly;
	delete[] M;
	delete[] flow_out;
	return (err_cnt == 0)? 0 : -1;


}

int test_top(){
	pix_t* hls_in;
	pix_t* in;
	data_t* hls_out;
	in = new pix_t[MAXSIZE * 2];
	hls_in = new pix_t[MAXSIZE * 2];
	hls_out = new data_t[MAXSIZE * 2];
	data_t **src_poly, **dst_poly, **M, **flow;
	src_poly = new data_t*[MAXSIZE];
	dst_poly = new data_t*[MAXSIZE];
	M = new data_t*[MAXSIZE];
	flow = new data_t*[MAXSIZE];
	for(int i=0;i<MAXSIZE;i++){
		src_poly[i] = new data_t[5];
		dst_poly[i] = new data_t[5];
		M[i] = new data_t[5];
		flow[i] = new data_t[2];
	}

	//source code
	for(int i=0;i<MAXSIZE;i++){
		hls_in[i*2] = in[i] = (i*i) % 179;
		hls_in[i*2+1] = in[i + MAXSIZE] = (i*i) %123;
	}

	//dut
	Farneback_top(hls_in, hls_out);

	//ref model
	Poly_Exp(in, src_poly, WIDTH, HEIGHT);
	Poly_Exp(in+MAXSIZE, dst_poly, WIDTH, HEIGHT);
	UpdateMat(src_poly, dst_poly, NULL, M, WIDTH, HEIGHT, 0);
	UpdateFlow(M, flow, WIDTH, HEIGHT);

	//check
	int err_cnt = 0;
	for(int i=0;i<MAXSIZE;i++){
		data_t score = (flow[i][0] - hls_out[i*2]) * (flow[i][0] - hls_out[i*2]) +
				(flow[i][1] - hls_out[i*2+1]) * (flow[i][1] - hls_out[i*2+1]);
		if(score > 1){
			err_cnt++;
			printf("[%d,%d]: (%f,%f)(%f,%f)\n",i/WIDTH, i%WIDTH, flow[i][0], flow[i][1], hls_out[i*2], hls_out[i*2+1]);
		}
	}
	if(err_cnt == 0)
		cout << "*** TOP TEST PASSED! ***" << endl;
	else
		cout << err_cnt << " errors are detected!\n"<< "*** TEST FAILED! ***" << endl;

	for(int i=0;i<MAXSIZE;i++){
		delete[]src_poly[i];
		delete[]dst_poly[i];
		delete[]M[i];
		delete[]flow[i];
	}

	delete[] hls_in;
	delete[] in;
	delete[] hls_out;
	delete[] src_poly;
	delete[] dst_poly;
	delete[] M;
	delete[] flow;
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
		cout << "*** POLY STRM TEST PASSED! ***" << endl;
	else
		cout << err_cnt << " errors are detected!\n"<< "*** TEST FAILED! ***" << endl;
	return (err_cnt == 0)? 0 : -1;


	for(int i=0;i<MAXSIZE;i++)
		delete[]out[i];
	delete[] in;
	delete[] out;
}

int test_poly_bram(){
	pix_t *in;
	data_t **out;
	data_t **flow_in;
	int rb = (WIDTH - 1) / POLY_EXP_SAMPLE_SIZE + 1;
	in = new pix_t[MAXSIZE];
	out = new data_t*[MAXSIZE];
	flow_in = new data_t*[MAXSIZE];
	for(int i=0;i<MAXSIZE;i++){
		out[i] = new data_t[5];
		flow_in[i] = new data_t[2];
	}
	pix_t in_hls[POLY_EXP_SAMPLE_SIZE*POLY_EXP_SAMPLE_SIZE][RAM_LENGTH];
	hls::stream<Data_2u> flow_in_hls("flow_in_hls");
	hls::stream<Data_5> out_hls("out_hls");

	for(int i=0;i<HEIGHT;i++){
		for(int j=0;j<WIDTH;j++){
			pix_t in_val = rand() % 234;
			pix_t fx, fy;
			fx = rand()%5 - 2;
			fy = rand()%5 - 2;
			int sel = i % POLY_EXP_SAMPLE_SIZE * POLY_EXP_SAMPLE_SIZE + j % POLY_EXP_SAMPLE_SIZE;
			int addr = i / POLY_EXP_SAMPLE_SIZE * rb + j / POLY_EXP_SAMPLE_SIZE;

			in[i*WIDTH+j] = in_val;
			in_hls[sel][addr] = in_val;
			Data_2u d2;
			flow_in[i*WIDTH+j][0] = d2.r0 = fx;
			flow_in[i*WIDTH+j][1] = d2.r1 = fy;
			flow_in_hls.write(d2);
		}
	}

	Poly_Exp(in, out, WIDTH, HEIGHT);
	Poly_Exp_hls_bram(in_hls, flow_in_hls, out_hls, WIDTH, HEIGHT);

	//check
	int err_cnt = 0;
	for(int i=0;i<HEIGHT;i++){
		for(int j=0;j<WIDTH;j++){
			int n = (POLY_EXP_SAMPLE_SIZE - 1) / 2;
			int fx, fy;
			int buf = i;
			fx = (int)(i + flow_in[i*WIDTH+j][0]);
			fy = (int)(j + flow_in[i*WIDTH+j][1]);
			fx = (fx < n) ? n : (fx >= HEIGHT - n) ? HEIGHT - n - 1: fx;
			fy = (fy < n) ? n : (fy >= WIDTH - n) ? WIDTH - n - 1: fy;
			Data_5 d5;
			d5 = out_hls.read();
			i = fx*WIDTH + fy;
			data_t sq_sum = 0;
			sq_sum += (out[i][0] - d5.r0)*(out[i][0] - d5.r0);
			sq_sum += (out[i][1] - d5.r1)*(out[i][1] - d5.r1);
			sq_sum += (out[i][2] - d5.r2)*(out[i][2] - d5.r2);
			sq_sum += (out[i][3] - d5.r3)*(out[i][3] - d5.r3);
			sq_sum += (out[i][4] - d5.r4)*(out[i][4] - d5.r4);

			if(sq_sum >= 1){
				cout <<"("<<buf<<","<<j<<")"<<endl;
				cout << out[i][0] << "\t"<< out[i][1] << "\t"<< out[i][2] << "\t"<< out[i][3] << "\t"<< out[i][4] << endl;
				cout << d5.r0 << "\t"<< d5.r1 << "\t"<<d5.r2 << "\t"<<d5.r3<< "\t"<< d5.r4 << endl << endl;
				err_cnt++;
			}
			i = buf;
		}
	}

	if(err_cnt == 0)
		cout << "*** POLY BRAM TEST PASSED! ***" << endl;
	else
		cout << err_cnt << " errors are detected!\n"<< "*** TEST FAILED! ***" << endl;

	for(int i=0;i<MAXSIZE;i++){
		delete[]out[i];
		delete[]flow_in[i];
	}
	delete[] flow_in;
	delete[] in;
	delete[] out;

	return (err_cnt == 0)? 0 : -1;
}

int main(){
	pix_t* mig_in;
	Data_2u* mig_out;
	pix_t con = 0x1;
	mig_in = new pix_t[MAXSIZE * 3];
	mig_out = new Data_2u[MAXSIZE];
	Farneback_top_v2_1(mig_in, mig_out, con);
	return 0;
	//void Farneback_top_v2_1(pix_t* mig_in, Data_2u* mig_out, pix_t con)
	//return test_top_v2_case0();
	/*
	hls::stream<pix_t> src_img;
	pix_t dst_img[POLY_EXP_SAMPLE_SIZE*POLY_EXP_SAMPLE_SIZE][RAM_LENGTH];
	hls::stream<Data_2u> flow_in;
	hls::stream<Data_2> flow_out;
	int width, height;
	width = 16;
	height = 12;
	const int p = POLY_EXP_SAMPLE_SIZE;
	const int br = (width-1)/p+1;
	for(int i=0;i<height;i++){
		for(int j=0;j<width;j++){
			int sel = i%p*p+j%p;
			int addr = i/p*br * j/p;
			dst_img[sel][addr] = (i*j)%123;
		}
	}


	for(int i=0;i<height;i++){
		for(int j=0;j<width;j++){
			src_img.write((i+j)%123);
			Data_2u val = {i%3-2, j%5-3};
			flow_in.write(val);
		}
	}
	Farneback_core(src_img, dst_img, flow_in, flow_out, width, height);

	for(int i=0;i<height;i++){
		for(int j=0;j<width;j++){
			flow_out.read();
		}
	}
	 */
}

/*
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
		cout << "*** MAT 2 1 TEST PASSED ***" << endl;
	else
		cout << err_cnt << " errors are detected!\n"<< "*** TEST FAILED ***" << endl;
	return (err_cnt == 0)? 0 : -1;
}

int test_mat(){
	data_t **src_poly;//[MAXSIZE][5];
	data_t **dst_poly;//[MAXSIZE][5];
	data_t **flow_in;//[MAXSIZE][2];
	data_t **M;//[MAXSIZE][5];

	int scale = 2;
	hls::stream<Data_5> hls_src_poly("hls_src_poly");
	hls::stream<Data_5> hls_dst_poly("hls_dst_poly");//[MAXSIZE];
	hls::stream<Data_2> hls_flow_in("hls_flow_in");
	hls::stream<Data_5>hls_M("hls_M");

	src_poly = new data_t*[MAXSIZE];
	dst_poly = new data_t*[MAXSIZE];
	flow_in = new data_t*[MAXSIZE];
	M = new data_t*[MAXSIZE];

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
		dst_poly[i][0] = i%155;
		dst_poly[i][1] = i%133;
		dst_poly[i][2] = i%144;
		dst_poly[i][3] = i%122;
		dst_poly[i][4] = i%111;
	}
	for(int i=0;i<WIDTH*HEIGHT; i++){
		Data_2 d2;
		data_t dx_f, dy_f;
		dx_f = flow_in[i][0] = d2.r0 = (i%4)/1.6-2;
		dy_f = flow_in[i][1] = d2.r1 = (i%3)/0.8-1;
		hls_flow_in.write(d2);
		int fx, fy, dx, dy;
		dx = (dx_f >= 0) ? (int)(dx_f + 0.5) : (int)(dx_f - 0.5);
		dy = (dy_f >= 0) ? (int)(dy_f + 0.5) : (int)(dy_f - 0.5);
		fx = i/WIDTH + dx;
		fy = i%WIDTH + dy;
		fx = (fx < 0) ? 0 : (fx >= HEIGHT) ? HEIGHT - 1 : fx;
		fy = (fy < 0) ? 0 : (fy >= WIDTH) ? WIDTH - 1 : fy;
		Data_5 d5;
		d5.r0 = dst_poly[fx*WIDTH+fy][0];
		d5.r1 = dst_poly[fx*WIDTH+fy][1];
		d5.r2 = dst_poly[fx*WIDTH+fy][2];
		d5.r3 = dst_poly[fx*WIDTH+fy][3];
		d5.r4 = dst_poly[fx*WIDTH+fy][4];
		hls_dst_poly.write(d5);
	}

	//ref model
	UpdateMat(src_poly, dst_poly, flow_in, M, WIDTH, HEIGHT, 1);
	//DUT
	UpdateMat_hls(hls_src_poly, hls_dst_poly, hls_flow_in, hls_M, WIDTH, HEIGHT);

	//check
	int err_cnt = 0;
	for(int i=0;i<WIDTH*HEIGHT;i++){
		Data_5 d5;
		d5 = hls_M.read();
		if(d5.r0 - M[i][0] > 0.1 ||d5.r0 - M[i][0] < -0.1 || d5.r1 - M[i][1] > 0.1 ||d5.r1 - M[i][1] < -0.1 ||
				d5.r2 - M[i][2] > 0.1 ||d5.r2 - M[i][2] < -0.1 ||d5.r3 - M[i][3] > 0.1 ||d5.r3 - M[i][3] < -0.1 ||
				d5.r4 - M[i][4] > 0.1 ||d5.r4 - M[i][4] < -0.1)
		{
			err_cnt++;
			printf("[%d:%d] REF(%f,%f,%f,%f,%f)DUT(%f,%f,%f,%f,%f)\n", i/WIDTH, i%WIDTH, M[i][0], M[i][1], M[i][2], M[i][3], M[i][4], d5.r0, d5.r1, d5.r2, d5.r3, d5.r4);
		}

	}

	if(err_cnt == 0)
		cout << "*** MAT TEST PASSED ***" << endl;
	else
		cout << err_cnt << " errors are detected!\n"<< "*** TEST FAILED ***" << endl;
	return (err_cnt == 0)? 0 : -1;
}
*/
