#include <assert.h>
#include <fcntl.h>
#include <getopt.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>

#include <sys/mman.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

#include <sys/time.h>
#include <math.h>
#include <iostream>
#include <opencv2/opencv.hpp>
#include "dense_flow.h"

using namespace std;
using namespace cv;

void DrawOptFlowMap(const Mat& flow, Mat& cflowmap, int step, const Scalar& color){
    for(int y=0;y<cflowmap.rows;y+=step){
        for(int x=0;x<cflowmap.cols;x+=step){
            const Point2f& fxy = flow.at<Point2f>(y,x);
            if(fxy.x*fxy.x + fxy.y*fxy.y > 4){            
                line(cflowmap, Point(x,y), Point(cvRound(x + fxy.x), cvRound(y + fxy.y)), color);
                circle(cflowmap, Point(cvRound(x + fxy.x), cvRound(y + fxy.y)), 1, color, -1);
            }
        }
    }
}

void DrawOptFlowMap(float* flow, Mat& cflowmap, int step, const Scalar& color){
    for(int y=0; y<cflowmap.rows; y+=step){
        for(int x=0; x<cflowmap.cols; x+=step){
            float fx = flow[(y*cflowmap.cols + x)*2];
            float fy = flow[(y*cflowmap.cols + x)*2 + 1];
            if(fx*fx + fy*fy > 4){
                line(cflowmap, Point(x,y), Point(cvRound(x+fx), cvRound(y+fy)),color);
                circle(cflowmap, Point(cvRound(x+fx), cvRound(y+fy)), 1, color, -1);
            }    
        }    
    }
}

int main(int argc, char *argv[]){    
    VideoCapture capture(0);
    if(!capture.isOpened()){
        printf("No device!\n");
        return -1;
    }
    cout << "fps=" << capture.get(CV_CAP_PROP_FPS) << endl
         << "rows=" << capture.get(CV_CAP_PROP_FRAME_WIDTH) << endl
         << "cols=" << capture.get(CV_CAP_PROP_FRAME_HEIGHT) << endl << endl;
    
    Mat pre, aft, pre_rgb, aft_rgb, flow;
    Mat fpga_img(480, 640, CV_8UC1);
    clock_t s_time, e_time;    
    uchar *in_buf, *out_buf;
    posix_memalign((void **)&in_buf, 4096/*alignment*/, 640*480*2 + 4096);
    posix_memalign((void **)&out_buf, 4096/*alignment*/, 640*480*2*sizeof(float) + 4096);
    float* flow_buf = new float[640*480*2];
    int fpga_w = open("/dev/xdma0_h2c_1", O_RDWR);
    int fpga_r = open("/dev/xdma0_c2h_0", O_RDWR|O_NONBLOCK);
    int rc;

    string in;
    cout <<"Please input the mode (f - FPGA/s - software):";
    cin >> in;
    bool isFPGA = (in == string("f"));

    while(true){
        cout << "start!" << endl;    
        capture >> pre_rgb;
        capture >> aft_rgb;
        cvtColor(pre_rgb, pre, CV_BGR2GRAY);
        cvtColor(aft_rgb, aft, CV_BGR2GRAY);

        if(isFPGA){
            cout << "FPGA mode:" << endl;
            s_time = clock();            
            for(int i=0;i<640*480;i++){
                in_buf[i*2] = pre.data[i];
                in_buf[i*2+1] = aft.data[i];       
            }
            unsigned int ret = lseek(fpga_w, 0x80000000, SEEK_SET);
            rc = write(fpga_w, in_buf, 640*480*2);
            //cout << "DMA write successfully!" << endl;
            assert(rc == 640*480*2);
            lseek(fpga_w, 0xC0000000, SEEK_SET);
            in_buf[0] = 1;
            rc = write(fpga_w, in_buf, 8);
            //cout << "write control signal successfully!" << endl;
            while(true){
                lseek(fpga_r, 0xC0000000, SEEK_SET);
                rc = read(fpga_r, out_buf, 8);
                if(*out_buf == 0)
                    break;
            }
            lseek(fpga_r, 0x90000000,SEEK_SET);       
            rc = read(fpga_r, out_buf, 640 * 480 * 2 * sizeof(float));
            //cout << "read successfully!" << endl;
            memcpy(flow_buf, out_buf, 640*480*2*sizeof(float));
            e_time = clock();  
            DrawOptFlowMap(flow_buf,pre_rgb,5, CV_RGB(0,255,0));
        }
        else{
            cout << "software mode:" << endl;
            s_time = clock();
            calcOpticalFlowFarneback(pre, aft, flow, 0.5, 0, 11, 1, 7, 1.5, 0);
            e_time = clock(); 
            DrawOptFlowMap(flow,pre_rgb, 5, CV_RGB(0,255,0));
        }
        cout <<"Total time:"<<(double)(e_time-s_time)/CLOCKS_PER_SEC<<"s"<<endl<<endl;      
        imshow("pic", pre_rgb);
        waitKey(1);
    }
    delete []in_buf;
    delete []out_buf;
    delete []flow_buf;
    return 0;    
}
