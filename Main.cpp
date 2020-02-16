
#include <stdio.h>
#include <stdlib.h>

#include  <GL/gl.h>
#include  <GL/glu.h>
#include  <GL/glut.h>/* glut.h includes gl.h and glu.h*/

//#include <my_include/gl.h>
//#include <my_include/glu.h>
//#include <my_include/glut.h>
#include  <math.h>
#include <time.h>
#include "globals.h"
#include <iostream>
#include <vector>

#include <dirent.h>
#include <iostream>
#include <fstream>
#include <vtkXMLStructuredGridWriter.h>
#include <vtkStructuredGridWriter.h>
#include <vtkStructuredGridReader.h>
#include <vtkStructuredGrid.h>
#include <vtkDataSet.h>
#include <vtkPointData.h>
#include <vtkCellData.h>
#include <vtkCellArray.h>
#include <vtkDoubleArray.h>

template<typename T>
class vec3 {
public:
    T x;
    T y;
    T z;
};

template<typename T>
class vec4 {
public:
    T x;
    T y;
    T z;
    T w;
};

vec3<float> *bodyAccel;
vec4<float> *bodyPos;
vec3<float> *bodyVel;

double x_save[4096][1000];
double y_save[4096][1000];
double z_save[4096][1000];

double u_save[4096][1000];
double v_save[4096][1000];
double w_save[4096][1000];

int save_num = 0;

int kCur= 0;
int* dims;
double* x;
double* y;
double* z;
double*** u;
double*** v;
double*** w;
double*** T;

double*** Wx;
double*** Wy;
double*** Wz;


double*** um;
double*** vm;
double*** wm;

double*** tilt;
double*** str;
double*** adv;

bool created=false;


vtkStructuredGridReader *reader;

vtkStructuredGrid *structured_grid;


void save_fields();
void load_fields();




void display(void);
void sweep_init();
void init();

void fmm_step(double dt);
void sweep();

double sc=1.0;
int view=VIEW_U;

int redr=0;

double ck=2.0;

double cv=0.001;

bool clearc=true;

const int maxParticles=4096;
int numParticles=0;//4096;

int i_tick=0;

void display(void)
{

    //printf ("%e %e %e  %e %e %e \n",x[0],x[NX-1],y[0],y[NY-1],z[0],z[NZ-1]);
    if (redr==1)
    {
        for (int i=0;i<1;i++)
            sweep();
    }

    if (clearc)
        glClear(GL_COLOR_BUFFER_BIT);


    glLoadIdentity();

    glTranslatef(x[NX/2],z[NZ/2],y[NY/2]);
    glScalef(sc,sc,sc);
    glRotatef(-90,1,0,0);
    glRotatef(ry,1.0,0,0);
    glRotatef(rx,0.0,1.0,0);
    //glRotatef(90,1,0,0);
    glTranslatef(-x[NX/2],-z[NZ/2],-y[NY/2]);


    glColor3f(1,1,1);




    double l_2;

    for (int i=0;i<NX-1;i++)
    {
        /*glBegin(GL_TRIANGLE_STRIP);


        for (int j=0;j<NY;j++)
        {
            if (view==VIEW_U)
                l_2=ck*(u[i][j][kCur]);
            if (view==VIEW_V)
                l_2=ck*(v[i][j][kCur]);
            if (view==VIEW_W)
                l_2=ck*(w[i][j][kCur]);
            glColor3f(l_2,l_2,-l_2);
            glVertex3f(x[i],y[j],z[kCur]);

            if (view==VIEW_U)
                l_2=ck*(u[i+1][j][kCur]);
            if (view==VIEW_V)
                l_2=ck*(v[i+1][j][kCur]);
            if (view==VIEW_W)
                l_2=ck*(w[i+1][j][kCur]);

            glColor3f(l_2,l_2,-l_2);
            glVertex3f(x[i+1],y[j],z[kCur]);
        }


        glEnd();*/

    }

    glEnable(GL_BLEND);

    glEnable(GL_POINT_SMOOTH);
    glPointSize(1.5);



    double vel_scale=0.0;
    for( int i=0; i<numParticles; i++ ) {
        vel_scale+=bodyVel[i].x*bodyVel[i].x+bodyVel[i].y*bodyVel[i].y+bodyVel[i].z*bodyVel[i].z;

    }

    vel_scale=1*sqrt(vel_scale/numParticles+0.0001);
    double leng_sacle=0.01*(x[NX-1]-x[0]);
    /*glBegin(GL_LINES);

    for( int i=0; i<numParticles; i++ )
    {
        // float pot=get_nearwall_potential(bodyPos[i].x,bodyPos[i].y);
        glColor3f(1.0*fabs(bodyVel[i].x)/vel_scale,4.0*fabs(bodyVel[i].y)/vel_scale,3.0*fabs(bodyVel[i].z)/vel_scale);
        // glColor3f(ck*pot,ck*pot,-ck*pot);
        glVertex3f(bodyPos[i].x,bodyPos[i].y,bodyPos[i].z);
        glColor3f(0.0,0.0,0.0);
        glVertex3f(bodyPos[i].x+leng_sacle*bodyVel[i].x/vel_scale,bodyPos[i].y+leng_sacle*bodyVel[i].y/vel_scale,bodyPos[i].z+leng_sacle*bodyVel[i].z/vel_scale);
    }
    glEnd();*/


    printf("ii=%d \n",save_num);
    for( int i=0; i<numParticles; i++ )
    {
        glBegin(GL_LINE_STRIP);


        for( int j=1; j<save_num-1; j++ )
        {
            // float pot=get_nearwall_potential(bodyPos[i].x,bodyPos[i].y);
            glColor3f(1.0*fabs(u_save[i][j])/vel_scale,4.0*fabs(v_save[i][j])/vel_scale,3.0*fabs(w_save[i][j])/vel_scale);
            // glColor3f(ck*pot,ck*pot,-ck*pot);
            glVertex3f(x_save[i][j],z_save[i][j],y_save[i][j]);
            //glColor3f(0.0,0.0,0.0);
            //glVertex3f(bodyPos[i].x+leng_sacle*bodyVel[i].x/vel_scale,bodyPos[i].y+leng_sacle*bodyVel[i].y/vel_scale,bodyPos[i].z+leng_sacle*bodyVel[i].z/vel_scale);
        }
        glEnd();
    }

    glLineWidth(2);


    glColor3f(0.5,0.5,0.5);

    glBegin(GL_LINE_LOOP);


    glVertex3f(x[0],z[0],y[0]);
    glVertex3f(x[NX-1],z[0],y[0]);
    glVertex3f(x[NX-1],z[0],y[NY-1]);
    glVertex3f(x[0],z[0],y[NY-1]);
    glEnd();

    glColor3f(1,1,1);

    glBegin(GL_LINE_LOOP);

    glVertex3f(x[0],z[NZ-1],y[0]);
    glVertex3f(x[NX-1],z[NZ-1],y[0]);
    glVertex3f(x[NX-1],z[NZ-1],y[NY-1]);
    glVertex3f(x[0],z[NZ-1],y[NY-1]);
    glEnd();
    glPointSize(3.0);
    glLineWidth(1.0);
    glutSwapBuffers();
    if (redr==1) glutPostRedisplay();

}

void m_m(int x,int y) //mouse move
{
    if (rotate==1)
    {
        rx=rx0+0.05*(x-mx0);
        ry=ry0+0.05*(y-my0);
    }
    glutPostRedisplay();
}



void m_d(int button, int state,int x, int y)  //mouse down
{
    if (state==GLUT_UP)
    {
        rotate=0;
        rx0=rx;
        ry0=ry;
    }
    if (state==GLUT_DOWN)
    {
        rotate=1;
        mx0=x;
        my0=y;
    }

    mouse_x=(1.0*x)/W_WIDTH;
    mouse_y=(W_HEIGHT-(1.0*y))/W_HEIGHT;
    glutPostRedisplay();
}




vec4<float> init_pos()
{
    vec4<float> res;

    res.x = rand()*1.0/RAND_MAX*x[NX-1];
    res.y = rand()*1.0/RAND_MAX*y[NY-1];
    res.z = 100;
    res.w = 0.5+my_rand(0)*0.5;

    return res;

}

vec3<float> init_vel()
{
    vec3<float> vel;
    vel.x = my_rand(0)*0.1-0.05;
    vel.y = my_rand(0)*0.1-0.05;
    vel.z = my_rand(0)*0.1-0.05;

    return vel;
}



void create_random_particles(int threadIdx, vec4<float> *bodyPos_, vec3<float> *bodyVel_,vec3<float> *bodyAccel_)
{
    int curNum = numParticles;
    int numToAdd = std::min(int(my_rand(threadIdx) * 300.0), maxParticles - numParticles-1);
    numParticles += numToAdd;
    for (int i = curNum; i < curNum + numToAdd; ++i)
    {
        vec4<float> pos=init_pos();
        vec3<float> vel=init_vel();
        bodyPos_[i].x=pos.x;
        bodyPos_[i].y=pos.y;
        bodyPos_[i].z=pos.z;
        bodyPos_[i].w=pos.w;
        bodyVel_[i]=vel;
        bodyAccel_[i].x=0.0;
        bodyAccel_[i].y=0.0;
        bodyAccel_[i].z=0.0;
    }

}


/*void delete_particle(int threadIdx, int particlesIdx, vec3<float> *bodyAccel_, vec4<float> *bodyPos_, vec3<float> *bodyVel_)
{
    bodyPos_[particlesIdx] = bodyPos_[numParticles-1];
    bodyVel_[particlesIdx] = bodyVel_[numParticles-1];
    bodyAccel_[particlesIdx] = bodyAccel_[numParticles-1];
    numParticles -= 1;
    //float r = my_rand(threadIdx);
    //randMax  = randMax > r ? randMax : r;
    // randMin  = randMin < r ? randMin : r;
}

*/

vec3<float> getV(vec4<float> bodyP)
{
    int i_,j_,k_;
    double a,b,c;

    i_=(int)((NX-1) * (bodyP.x - x[0]) / (x[NX-1] - x[0]));
    j_=(int)((NY-1) * (bodyP.y - y[0]) / (y[NY-1] - y[0]));

    k_=0;
    while ((z[k_] < bodyP.z) && (k_ < NZ-2))
        k_++;

    a=(bodyP.x - x[i_]) / (x[i_+1] - x[i_]);
    b=(bodyP.y - y[j_]) / (y[j_+1] - y[j_]);
    c=(bodyP.z - z[k_]) / (z[k_+1] - z[k_]);

    vec3<float> res;

    i_=fmax(1,i_);
    i_=fmin(NX-2,i_);
    j_=fmax(1,j_);
    j_=fmin(NY-2,j_);
    k_=fmax(1,k_);
    k_=fmin(NZ-2,k_);



    a=fmax(0.0,a);
    a=fmin(1.0,a);
    b=fmax(0.0,b);
    b=fmin(1.0,b);
    c=fmax(0.0,c);
    c=fmin(1.0,c);


    printf("i=%d j=%d k=%d  %f %f %f\n", i_,j_,k_,a,b,c);

    res.x=  ((u[i_][j_][k_]*(1.0-a)  +u[i_+1][j_][k_]  *(a))*(1.0-b) + (u[i_][j_+1][k_]*(1.0-a)  +u[i_+1][j_+1][k_]*a)*b)*(1.0-c)+
            ((u[i_][j_][k_+1]*(1.0-a)+u[i_+1][j_][k_+1]*(a))*(1.0-b) + (u[i_][j_+1][k_+1]*(1.0-a)+u[i_+1][j_+1][k_+1]*a)*b)*(c);

    res.y=((v[i_][j_][k_]*(1.0-a)+v[i_+1][j_][k_]*(a))*(1.0-b) + (v[i_][j_+1][k_]*(1.0-a)+v[i_+1][j_+1][k_]*(a))*b)*(1.0-c)+
            ((v[i_][j_][k_+1]*(1.0-a)+v[i_+1][j_][k_+1]*(a))*(1.0-b) + (v[i_][j_+1][k_+1]*(1.0-a)+v[i_+1][j_+1][k_+1]*(a))*b)*(c);

    res.z=((w[i_][j_][k_]*(1.0-a)+w[i_+1][j_][k_]*(a))*(1.0-b) + (w[i_][j_+1][k_]*(1.0-a)+w[i_+1][j_+1][k_]*(a))*b)*(1.0-c)+
            ((w[i_][j_][k_+1]*(1.0-a)+w[i_+1][j_][k_+1]*(a))*(1.0-b) + (w[i_][j_+1][k_+1]*(1.0-a)+w[i_+1][j_+1][k_+1]*(a))*b)*(c);

    return res;
}

int counter=0;
void fmm_step(double dt)
{
    int i;


    counter++;

    if (counter>10)
    {
    for( i=0; i<numParticles; i++ )
    {
        x_save[i][save_num]=bodyPos[i].x;
        y_save[i][save_num]=bodyPos[i].y;
        z_save[i][save_num]=bodyPos[i].z;
        u_save[i][save_num]=bodyVel[i].x;
        v_save[i][save_num]=bodyVel[i].y;
        w_save[i][save_num]=bodyVel[i].z;
    }
    save_num++;
    counter=0;
    }
    for( i=0; i<numParticles; i++ )
    {
        float magn=0.001;//qe/Me;//1e-1;
        vec3<float> ev = getV(bodyPos[i]);



        bodyVel[i].x -= magn*(bodyVel[i].x - ev.x)*dt;//-= magn*(ev.x )*dt;
        bodyVel[i].y -= magn*(bodyVel[i].y - ev.y)*dt;//-= magn*(ev.y )*dt;
        bodyVel[i].z -= magn*(bodyVel[i].z - ev.z)*dt;//-= magn*(ev.z )*dt;

        printf("vx=%e vy=%e vz=%e  \n",ev.x,ev.y,ev.z);
    }



    for( i=0; i<numParticles; i++ )
    {
        bodyPos[i].x += dt*bodyVel[i].x;
        bodyPos[i].y += dt*bodyVel[i].y;
        bodyPos[i].z += dt*bodyVel[i].z;

    }



}

void fmm_init()
{
    int i;
    rand_init();
    bodyAccel = new vec3<float>[maxParticles];
    bodyVel = new vec3<float>[maxParticles];
    bodyPos = new vec4<float>[maxParticles];

    create_random_particles(0, bodyPos, bodyVel,bodyAccel);



    /*   for( i=0; i<maxParticles; i++ ) {
        bodyPos[i].x = rand()/(float) RAND_MAX*(w_x1-w_x0)*0.95+w_x0;
        bodyPos[i].y = rand()/(float) RAND_MAX*(w_y1-0)*0.95+0;
        bodyPos[i].z = rand()/(float) RAND_MAX*(w_z1-w_z0)+w_z0;
        bodyPos[i].w = 1.0;

        bodyVel[i].x = 300.0*(0.4+rand()/(float) RAND_MAX*0.1);
        bodyVel[i].y = 300.0*(-0.005+rand()/(float) RAND_MAX*0.01-0.005);
        bodyVel[i].z = 300.0*(rand()/(float) RAND_MAX*0.01-0.005);
    }

    for (i=0; i<N_X;i++)
    {
        BoundaryLayer[i]=0.000001;
    }
*/

}




void kb(unsigned char key, int x, int y)
{
    int i,j,k,nn,n;
    double m,sum;
    double max_err=0.0;
    if (key==']')
    {
        //ck*=1.1;
        if (kCur<NZ)
            kCur++;


    }


    if (key=='[')
    {
        if (kCur>0)
            kCur--;
        //ck/=1.1;
    }


    if (key=='q')
    {
        sc*=1.1;
    }

    if (key=='e')
    {
        sc/=1.1;
    }


    if (key=='.')
    {
        ck*=1.1;


    }


    if (key==',')
    {
        ck/=1.1;
    }


    if (key==']')
    {
        dt*=1.1;
        printf("dt=%e \n",dt);
    }


    if (key=='[')
    {


        dt/=1.1;
        printf("dt=%e \n",dt);

    }


    if (key=='1')
    {

        view=VIEW_U;
        printf("viewing PHI \n");

    }


    if (key=='2')
    {

        view=VIEW_V;
        printf("viewing E \n");

    }


    if (key=='m')
    {

        move_particles=!move_particles;

    }

    if (key=='3')
    {

        view=VIEW_W;
        printf("viewing P \n");

    }


    if (key=='n')
    {
        clearc=!clearc;
    }

    if (key=='s')
    {
        //   for(int i=0;i<10;i++)
        sweep();
        //   fmm_step(0.0000001);
    }

    if (key==' ')
    {
        redr=!redr;
        // filter_conv(iglob,1,NULL,false);
        //iglob++;
    }


    glutPostRedisplay();
}

void sweep_init()
{
    fmm_init();
}

void sweep()
{
    //if (move_particles)
    {
        fmm_step(dt);
        printf("numParticles=%d \n",numParticles);
    }


}


void init()
{
    glClearColor (0.0, 0.0, 0.0, 0.0);
    glColor3f(1.0, 1.0, 1.0);
    glMatrixMode (GL_PROJECTION);
    glLoadIdentity ();
    glOrtho(x[0]-1000, x[NX-1]+1000,  z[0]-1000, z[NZ-1]+1000,y[0]-100000, y[NY-1]+100000);
    glMatrixMode (GL_MODELVIEW);
    rand_init();
    sweep_init();

    // fmm_init();


}

/*void fmm_step(double dt)
{
    numParticles ++;
    bodyPos_[numParticles - 1].x = x[NX/2];
    bodyPos_[numParticles - 1].y = y[NY/2];
    bodyPos_[numParticles - 1].z = 30;
    bodyPos_[numParticles - 1].w = 1;
    bodyVel_[numParticles - 1] = vel;
    bodyAccel_[numParticles - 1].x = 0.0;
    bodyAccel_[numParticles - 1].y = 0.0;
    bodyAccel_[numParticles - 1].z = 0.0;

}*/

void loadVTK(char *name)
{
    if(!created)
    {
        reader = vtkStructuredGridReader::New();
        //structured_grid = vtkStructuredGrid::New();
    }

    vtkDataArray *Temp;
    vtkDataArray *vel;

    reader->ReadAllScalarsOn();
    reader->ReadAllVectorsOn();
    reader->SetFileName(name);

    reader->Update();
    Temp=reader->GetOutput()->GetPointData()->GetArray("Temperature");
    vel=reader->GetOutput()->GetPointData()->GetArray("Velocity");
    structured_grid =reader->GetOutput();
    dims = structured_grid->GetDimensions();

    if(!created)
    {

        NX = dims[0];
        NY = dims[1];
        NZ =  dims[2];
        x=new double [dims[0]];
        y=new double [dims[1]];
        z=new double [dims[2]];

        u=new double **[dims[0]];
        v=new double **[dims[0]];
        w=new double **[dims[0]];


        um=new double **[dims[0]];
        vm=new double **[dims[0]];
        wm=new double **[dims[0]];

        T=new double **[dims[0]];

        Wx=new double **[dims[0]];
        Wy=new double **[dims[0]];
        Wz=new double **[dims[0]];

        tilt=new double **[dims[0]];
        str=new double **[dims[0]];
        adv=new double **[dims[0]];

        for (int i=0;i<dims[0];i++)
        {
            u[i]=new double *[dims[1]];
            v[i]=new double *[dims[1]];
            w[i]=new double *[dims[1]];

            um[i]=new double *[dims[1]];
            vm[i]=new double *[dims[1]];
            wm[i]=new double *[dims[1]];
            T[i]=new double *[dims[1]];

            Wx[i]=new double *[dims[1]];
            Wy[i]=new double *[dims[1]];
            Wz[i]=new double *[dims[1]];


            tilt[i]=new double *[dims[1]];
            str[i]=new double *[dims[1]];
            adv[i]=new double *[dims[1]];

            for (int j=0;j<dims[1];j++)
            {
                u[i][j]=new double [dims[2]];
                v[i][j]=new double [dims[2]];
                w[i][j]=new double [dims[2]];

                um[i][j]=new double [dims[2]];
                vm[i][j]=new double [dims[2]];
                wm[i][j]=new double [dims[2]];

                T[i][j]=new double [dims[2]];

                Wx[i][j]=new double [dims[2]];
                Wy[i][j]=new double [dims[2]];
                Wz[i][j]=new double [dims[2]];

                tilt[i][j]=new double [dims[2]];
                str[i][j]=new double [dims[2]];
                adv[i][j]=new double [dims[2]];

                for (int k=0;k<dims[2];k++)
                {
                    u[i][j][k]=0.0;
                    v[i][j][k]=0.0;
                    w[i][j][k]=0.0;
                    um[i][j][k]=0.0;
                    vm[i][j][k]=0.0;
                    wm[i][j][k]=0.0;

                    T[i][j][k]=0.0;

                    Wx[i][j][k]=0.0;
                    Wy[i][j][k]=0.0;
                    Wz[i][j][k]=0.0;

                    tilt[i][j][k]=0.0;
                    str[i][j][k]=0.0;
                    adv[i][j][k]=0.0;
                }
            }
        }
        created=true;
        printf("done alloc %d %d %d \n",dims[0],dims[1],dims[2]);
    }


    reader->CloseVTKFile();
    double p[3];

    for (int k = 0; k < dims[2]; k++)
    {
        for (int j = 0; j < dims[1]; j++)
        {
            for (int i = 0; i < dims[0]; i++)
            {
                structured_grid->GetPoint(i, j, k, p);
                x[i]=p[0];
            }
            y[j]=p[1];
        }
        z[k]=p[2];
    }


    for (int k = 0; k < dims[2]; k++)
    {
        for (int j = 0; j < dims[1]; j++)
        {
            for (int i = 0; i < dims[0]; i++)
            {
                int ijk=i+j*dims[0]+k*dims[1]*dims[1];
                double tmp;
                u[i][j][k]=vel->GetComponent(ijk, 0);
                v[i][j][k]=vel->GetComponent(ijk, 1);
                w[i][j][k]=vel->GetComponent(ijk, 2);
                Temp->GetTuple(ijk, &tmp);
                T[i][j][k]+=tmp;

                um[i][j][k]+=u[i][j][k];//(w[i][j+1][k]-w[i][j-1][k])/(y[j+1]-y[j-1])-(v[i][j][k+1]-v[i][j][k-1])/(z[k+1]-z[k-1]);
                vm[i][j][k]+=v[i][j][k];//(u[i][j][k+1]-u[i][j][k-1])/(z[k+1]-z[k-1])-(w[i+1][j][k]-v[i-1][j][k])/(x[i+1]-x[i-1]);
                wm[i][j][k]+=w[i][j][k];//(v[i+1][j][k]-v[i-1][j][k])/(x[i+1]-x[i-1])-(u[i][j+1][k]-u[i][j-1][k])/(y[j+1]-y[j-1]);
            }
        }
    }

    // structured_grid->Delete();
    // reader->Delete();
}

void calcNewFields()
{
    for (int k = 1; k < dims[2]-1; k++)
    {
        for (int j = 1; j < dims[1]-1; j++)
        {
            for (int i = 1; i < dims[0]-1; i++)
            {
                um[i][j][k]+=u[i][j][k];//(w[i][j+1][k]-w[i][j-1][k])/(y[j+1]-y[j-1])-(v[i][j][k+1]-v[i][j][k-1])/(z[k+1]-z[k-1]);
                vm[i][j][k]+=v[i][j][k];//(u[i][j][k+1]-u[i][j][k-1])/(z[k+1]-z[k-1])-(w[i+1][j][k]-v[i-1][j][k])/(x[i+1]-x[i-1]);
                wm[i][j][k]+=w[i][j][k];//(v[i+1][j][k]-v[i-1][j][k])/(x[i+1]-x[i-1])-(u[i][j+1][k]-u[i][j-1][k])/(y[j+1]-y[j-1]);
            }
        }
    }

}


void calcOmega_gen()
{
    for (int k = 1; k < dims[2]-1; k++)
    {
        for (int j = 1; j < dims[1]-1; j++)
        {
            for (int i = 1; i < dims[0]-1; i++)
            {
                Wx[i][j][k]=(w[i][j+1][k]-w[i][j-1][k])/(y[j+1]-y[j-1])-(v[i][j][k+1]-v[i][j][k-1])/(z[k+1]-z[k-1]);
                Wy[i][j][k]=(u[i][j][k+1]-u[i][j][k-1])/(z[k+1]-z[k-1])-(w[i+1][j][k]-v[i-1][j][k])/(x[i+1]-x[i-1]);
                Wz[i][j][k]=(v[i+1][j][k]-v[i-1][j][k])/(x[i+1]-x[i-1])-(u[i][j+1][k]-u[i][j-1][k])/(y[j+1]-y[j-1]);
            }
        }
    }

    for (int k = 1; k < dims[2]-1; k++)
    {
        for (int j = 1; j < dims[1]-1; j++)
        {
            for (int i = 1; i < dims[0]-1; i++)
            {
                double wx_duzdx=Wx[i][j][k]*(w[i+1][j][k]-w[i-1][j][k])/(x[i+1]-x[i-1]);
                double wy_duzdy=Wy[i][j][k]*(w[i][j+1][k]-w[i][j+1][k])/(y[j+1]-y[j-1]);
                double wz_duzdz=Wz[i][j][k]*(w[i][j][k+1]-w[i][j][k-1])/(z[k+1]-z[k-1]);

                double ux_dwzdx=u[i][j][k]*(Wz[i+1][j][k]-Wz[i-1][j][k])/(x[i+1]-x[i-1]);
                double uy_dwzdy=v[i][j][k]*(Wz[i][j+1][k]-Wz[i][j+1][k])/(y[j+1]-y[j-1]);
                double uz_dwzdz=w[i][j][k]*(Wz[i][j][k+1]-Wz[i][j][k-1])/(z[k+1]-z[k-1]);



                tilt[i][j][k]+=Wz[i][j][k]*(wx_duzdx+wy_duzdy);
                str[i][j][k]+=Wz[i][j][k]*(wz_duzdz);
                adv[i][j][k]+=-Wz[i][j][k]*(ux_dwzdx+uy_dwzdy+uz_dwzdz);
            }
        }
    }
}


void saveVTK(char *name,int num)
{
    vtkStructuredGridWriter *writer = vtkStructuredGridWriter::New();
    vtkStructuredGrid *structured_grid_new = vtkStructuredGrid::New();
    structured_grid_new->SetDimensions(dims[0],dims[1],dims[2]);
    vtkPoints *nodes = vtkPoints::New();
    nodes->Allocate(dims[0]*dims[1]*dims[2]);
    vtkDoubleArray* velArray = vtkDoubleArray::New();
    velArray->SetNumberOfComponents(3);
    velArray->SetName("Velocity");
    vtkDoubleArray* vortArray = vtkDoubleArray::New();
    vortArray->SetNumberOfComponents(3);
    vortArray->SetName("Vel_mean");
    vtkDoubleArray* TArray = vtkDoubleArray::New();
    TArray->SetNumberOfComponents(1);
    TArray->SetName("T_mean");

    vtkDoubleArray* tilt_a = vtkDoubleArray::New();
    tilt_a->SetNumberOfComponents(1);
    tilt_a->SetName("wz_tilt");
    vtkDoubleArray* str_a = vtkDoubleArray::New();
    str_a->SetNumberOfComponents(1);
    str_a->SetName("wz_str");
    vtkDoubleArray* adv_a = vtkDoubleArray::New();
    adv_a->SetNumberOfComponents(1);
    adv_a->SetName("wz_adv");


    for (int k = 0; k < dims[2]; k++)
    {
        for (int j = 0; j < dims[1]; j++)
        {
            for (int i = 0; i < dims[0]; i++)
            {
                nodes->InsertNextPoint(x[i], y[j], z[k]);
                velArray->InsertNextTuple3(u[i][j][k], v[i][j][k], w[i][j][k]);
                TArray->InsertNextTuple1(T[i][j][k]/num);

                tilt_a->InsertNextTuple1(tilt[i][j][k]/num);
                str_a->InsertNextTuple1(str[i][j][k]/num);
                adv_a->InsertNextTuple1(adv[i][j][k]/num);


                vortArray->InsertNextTuple3(um[i][j][k]/num, vm[i][j][k]/num, wm[i][j][k]/num);
            }
        }
    }

    structured_grid_new->SetPoints(nodes);
    structured_grid_new->GetPointData()->AddArray(velArray);
    structured_grid_new->GetPointData()->AddArray(TArray);
    structured_grid_new->GetPointData()->AddArray(tilt_a);
    structured_grid_new->GetPointData()->AddArray(str_a);
    structured_grid_new->GetPointData()->AddArray(adv_a);

    structured_grid_new->GetPointData()->AddArray(vortArray);


#if (VTK_MAJOR_VERSION >=6)
    writer->SetInputData(structured_grid_new);
#else
    writer->SetInput(structured_grid_new);
#endif

    writer->SetFileName(name);
    writer->SetFileTypeToBinary();
    writer->Write();

    structured_grid_new->Delete();
    writer->Delete();
    nodes->Delete();
    velArray->Delete();
    TArray->Delete();
    tilt_a->Delete();
    str_a->Delete();
    adv_a->Delete();
    vortArray->Delete();

}

int main(int argc, char** argv)
{
    //srand(time(NULL));

    std::string  nameFile,nameForOpen,nameForSave;
    std::string nameDir("/home/user/RFFI_Petya/ferroPlasma/vtk/");//("../");
    DIR *mydir = opendir(nameDir.data());
    if(mydir == NULL) {
        perror("opendir");
        return -1;
    }
    printf("start\n");
    struct dirent *entry;
    int fnum=0;
    int i=0;
    while ((entry = readdir(mydir))&&(i<3))
    {
        i++;
        //entry = readdir(mydir);
        int len = strlen (entry->d_name);
        if (len >= 4) {
            if (strcmp (".vtk", entry->d_name + len - 4) == 0) {
                fnum++;
                nameFile=entry->d_name;
                nameForOpen = nameDir+ nameFile;
                printf("%s \n",nameForOpen.data());
                loadVTK( const_cast<char*>(nameForOpen.data()));
                //   calcOmega_gen();
                //                calcNewFields();

                //                nameForSave = "vtk/"+ nameFile.substr(0,nameFile.size()-4)+"_new.vtk";
                //                saveVTK( const_cast<char*>(nameForSave.data()));
            }
        }
        //break;
    }
    closedir(mydir);
    printf("end\n");
    glutInit(&argc,argv);
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB);
    glutInitWindowSize(W_HEIGHT,W_HEIGHT);
    glutInitWindowPosition(0,0);
    glutCreateWindow("simple");
    glutDisplayFunc(display);
    glutMotionFunc(m_m);
    glutMouseFunc(m_d);
    glutKeyboardFunc(kb);
    init();
    glutMainLoop();
    delete x;
    delete y;
    delete z;
    delete u;
    delete v;
    delete w;
    delete T;
}
