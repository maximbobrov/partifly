

//just a test comment
#ifndef GLOBALS_H
#define GLOBALS_H

 //for multigrid
 //#define N_X 1025
 //#define N_Y 513
//#define N_Y_DIEL 255

#define M_PI 3.1415926535

 #define W_WIDTH 1400
 #define W_HEIGHT 900

 #define RES2_MIN 0.000001

#define VIEW_U 0
#define VIEW_V 1
#define VIEW_W 2

void rand_init();

extern int NX;
extern int NY;
extern int NZ;

float my_rand(int i);


extern bool move_particles;

extern  int itn;
extern  int clear_w;


extern  float rx;
extern  float ry;
extern  int mx0,my0;
extern  int rotate;
extern  float rx0;
extern  float ry0;
extern double d_temp;
extern double mouse_x,mouse_y;

extern double r_2g;//cg res
extern double dt;
#endif
