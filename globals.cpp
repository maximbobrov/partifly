
#include "globals.h"
#include<xmmintrin.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <math.h>


int num_thread=1;
uint32_t xr[32],yr[32],zr[32],wr[32];
int NX;
int NY;
int NZ;
void rand_init()
{
    for (int i=0;i<num_thread;i++)
    {
        xr[i] = 123456789+i;
        yr[i] = 362436069+2*i;
        zr[i] = 521288629+3*i;
        wr[i] = 88675123+5*i;
    }
}

float my_rand(int i) {
    uint32_t t;
    t = xr[i] ^ (xr[i] << 11);
    xr[i] = yr[i]; yr[i] = zr[i]; zr[i] = wr[i];
    wr[i] = wr[i] ^ (wr[i] >> 19) ^ (t ^ (t >> 8));
    return wr[i]*0.5/0x7FFFFFFF;
}

 bool move_particles=false;


int itn=0;
int clear_w=1.0;
//int t=0;

double dt=5.0;
float rx=0;
float ry=0;
int mx0,my0;
int rotate=0;
float rx0=0;
float ry0=0;
double d_temp;
double mouse_x,mouse_y;

double r_2g=0.0;//cg res
