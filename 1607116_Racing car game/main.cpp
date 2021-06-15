#include<windows.h>
#include <GL/glut.h>
#include<bits/stdc++.h>
#include <stdlib.h>
#define rad (3.1416/180)
#define EN_SIZE 20
#include "BmpLoder.h"

const double PI = 3.14159265389;

int anglex= 0, angley = 0, anglez = 0;          //rotation angles
int window;
int wired=0;
int shcpt=1;
int animat = 0;
const int L=5;
const int L1=14;
const int dgre=3;
int ncpt=L+1;
int clikd=0;
const int nt = 40;				//number of slices along x-direction
const int ntheta = 30;

GLfloat ctrlpoints[L+1][3] =
{
     { 0.0, 0.0, 0.0},{0.1,1.5,0},{0.5,2.0,0},{0.9,3.0,0},{1.4,4,0},{1.9,5,0}

    /*{ 0.0, 0.0, 0.0}, { -0.3, 0.5, 0.0},
    { 0.1, 1.7, 0.0},{ 0.5, 1.5, 0.0},
    {1.0, 1.5, 0.0}, {1.4, 1.4, 0.0},
    {1.8, 0.4, 0.0},{2.2, 0.4, 0.0},
    {2.6, 1.5, 0.0}, {3.0, 1.4, 0.0},
    {3.4, 1.4, 0.0},{3.8, 1.4, 0.0},
    {4.2, 1.0, 0.0},{4.6, 1.0, 0.0},
    {5.0, 1.0, 0.0},{5.4, 1.0, 0.0},
    {5.8, 0.5, 0.0},{6.2, 0.5, 0.0},
    {6.6, 0.5, 0.0},{7.2, 0.2, 0.0},
    {6.8, 0.52, 0.0}*/
};

GLfloat ctrlpoints_tunnel[L1+1][3]=
{
    {73.2249, 23.684, 0},
    {73.2249, 23.684, 0},
    {-29.9166, 23.602, 0},
    {-29.9166, 23.602, 0},
    {51.3486, 23.4378, 0},
    {51.3486, 23.4378, 0},
    {58.2295, 23.2681, 0},
    {58.2295, 23.2681, 0},
    {-21.0119, 23.0602, 0},
    {-21.0119, 23.0602, 0},
    {-17.9402, 22.9124, 0},
    {-17.9402, 22.9124, 0},
    {60.5347, 22.7647, 0},
    {60.5347, 22.7647, 0},
    {48.1226, 22.6935, 0}
};

float wcsClkDn[3],wcsClkUp[3];
///////////////////////////////
class point1
{
public:
    point1()
    {
        x=0;
        y=0;
    }
    int x;
    int y;
} clkpt[2];
int flag=0;
GLint viewport[4]; //var to hold the viewport info
GLdouble modelview[16]; //var to hold the modelview info
GLdouble projection[16]; //var to hold the projection matrix info


void scsToWcs(float sx,float sy, float wcsv[3] );
void processMouse(int button, int state, int x, int y);
void matColor(float kdr, float kdg, float kdb,  float shiny, int frnt_Back=0, float ambFactor=1.0, float specFactor=1.0);
///////////////////////////


void scsToWcs(float sx,float sy, float wcsv[3] )
{

    GLfloat winX, winY, winZ; //variables to hold screen x,y,z coordinates
    GLdouble worldX, worldY, worldZ; //variables to hold world x,y,z coordinates

    //glGetDoublev( GL_MODELVIEW_MATRIX, modelview ); //get the modelview info
    glGetDoublev( GL_PROJECTION_MATRIX, projection ); //get the projection matrix info
    glGetIntegerv( GL_VIEWPORT, viewport ); //get the viewport info

    winX = sx;
    winY = (float)viewport[3] - (float)sy;
    winZ = 0;

    //get the world coordinates from the screen coordinates
    gluUnProject( winX, winY, winZ, modelview, projection, viewport, &worldX, &worldY, &worldZ);
    wcsv[0]=worldX;
    wcsv[1]=worldY;
    wcsv[2]=worldZ;


}
void processMouse(int button, int state, int x, int y)
{
    if(button==GLUT_LEFT_BUTTON && state==GLUT_DOWN)
    {
        if(flag!=1)
        {
            flag=1;
            clkpt[0].x=x;
            clkpt[0].y=y;
        }


        scsToWcs(clkpt[0].x,clkpt[0].y,wcsClkDn);
        //cout<<"\nD: "<<x<<" "<<y<<" wcs: "<<wcsClkDn[0]<<" "<<wcsClkDn[1];
    }
    else if(button==GLUT_LEFT_BUTTON && state==GLUT_UP)
    {
        if (flag==1)
        {
            clkpt[1].x=x;
            clkpt[1].y=y;
            flag=0;
        }
        float wcs[3];
        scsToWcs(clkpt[1].x,clkpt[1].y,wcsClkUp);
        //cout<<"\nU: "<<x<<" "<<y<<" wcs: "<<wcsClkUp[0]<<" "<<wcsClkUp[1];

        clikd=!clikd;
    }
}

//control points
long long nCr(int n, int r)
{
    if(r > n / 2) r = n - r; // because C(n, r) == C(n, n - r)
    long long ans = 1;
    int i;

    for(i = 1; i <= r; i++)
    {
        ans *= n - r + i;
        ans /= i;
    }

    return ans;
}

//polynomial interpretation for N points
void BezierCurve ( double t,  float xy[2])
{
    double y=0;
    double x=0;
    t=t>1.0?1.0:t;
    for(int i=0; i<=L; i++)
    {
        int ncr=nCr(L,i);
        double oneMinusTpow=pow(1-t,double(L-i));
        double tPow=pow(t,double(i));
        double coef=oneMinusTpow*tPow*ncr;
        x+=coef*ctrlpoints[i][0];
        y+=coef*ctrlpoints[i][1];

    }
    xy[0] = float(x);
    xy[1] = float(y);

    //return y;
}

void BezierCurvet ( double t,  float xy[2])
{
    double y=0;
    double x=0;
    t=t>1.0?1.0:t;
    for(int i=0; i<=L1; i++)
    {
        int ncr=nCr(L1,i);
        double oneMinusTpow=pow(1-t,double(L1-i));
        double tPow=pow(t,double(i));
        double coef=oneMinusTpow*tPow*ncr;
        x+=coef*ctrlpoints_tunnel[i][0];
        y+=coef*ctrlpoints_tunnel[i][1];

    }
    xy[0] = float(x);
    xy[1] = float(y);

    //return y;
}

///////////////////////
void setNormal(GLfloat x1, GLfloat y1,GLfloat z1, GLfloat x2, GLfloat y2,GLfloat z2, GLfloat x3, GLfloat y3,GLfloat z3)
{
    GLfloat Ux, Uy, Uz, Vx, Vy, Vz, Nx, Ny, Nz;

    Ux = x2-x1;
    Uy = y2-y1;
    Uz = z2-z1;

    Vx = x3-x1;
    Vy = y3-y1;
    Vz = z3-z1;

    Nx = Uy*Vz - Uz*Vy;
    Ny = Uz*Vx - Ux*Vz;
    Nz = Ux*Vy - Uy*Vx;

    glNormal3f(-Nx,-Ny,-Nz);
}

void tunnelBezier()
{
    int i, j;
    float x, y, z, r;				//current coordinates
    float x1, y1, z1, r1;			//next coordinates
    float theta;

    const float startx = 0, endx = ctrlpoints_tunnel[L1][0];
    //number of angular slices
    const float dx = (endx - startx) / nt;	//x step size
    const float dtheta = 2*PI / ntheta;		//angular step size

    float t=0;
    float dt=1.0/nt;
    float xy[2];
    BezierCurvet( t,  xy);
    x = xy[0];
    r = xy[1];
    //rotate about z-axis
    float p1x,p1y,p1z,p2x,p2y,p2z;
    for ( i = 0; i < nt; ++i )  			//step through x
    {
        theta = 0;
        t+=dt;
        BezierCurvet( t,  xy);
        x1 = xy[0];
        r1 = xy[1];

        //draw the surface composed of quadrilaterals by sweeping theta
        glBegin( GL_QUAD_STRIP );
        for ( j = 0; j <= ntheta/2; ++j )
        {
            theta += dtheta;
            double cosa = cos( theta );
            double sina = sin ( theta );
            y = r * cosa;
            y1 = r1 * cosa;	//current and next y
            z = r * sina;
            z1 = r1 * sina;	//current and next z

            //edge from point at x to point at next x
            glVertex3f (x, y, z);

            if(j>0)
            {
                setNormal(p1x,p1y,p1z,p2x,p2y,p2z,x, y, z);
            }
            else
            {
                p1x=x;
                p1y=y;
                p1z=z;
                p2x=x1;
                p2y=y1;
                p2z=z1;

            }
            glVertex3f (x1, y1, z1);

            //forms quad with next pair of points with incremented theta value
        }
        glEnd();
        x = x1;
        r = r1;
    } //for i

}

void station_top()
{
    int i, j;
    float x, y, z, r;				//current coordinates
    float x1, y1, z1, r1;			//next coordinates
    float theta;

    const float startx = 0, endx = ctrlpoints[L][0];
    //number of angular slices
    const float dx = (endx - startx) / nt;	//x step size
    const float dtheta = 2*PI / ntheta;		//angular step size

    float t=0;
    float dt=1.0/nt;
    float xy[2];
    BezierCurve( t,  xy);
    x = xy[0];
    r = xy[1];
    //rotate about z-axis
    float p1x,p1y,p1z,p2x,p2y,p2z;
    for ( i = 0; i < 30; ++i )  			//step through x
    {
        theta = 0;
        t+=dt;
        BezierCurve( t,  xy);
        x1 = xy[0];
        r1 = xy[1];

        //draw the surface composed of quadrilaterals by sweeping theta
        glBegin( GL_QUAD_STRIP );
        //glBegin( GL_QUADS );

        for ( j = 0; j <= ntheta; ++j )
        {
            theta += dtheta;
            double cosa = cos( theta );
            double sina = sin ( theta );
            y = r * cosa;
            y1 = r1 * cosa;	//current and next y
            z = r * sina;
            z1 = r1 * sina;	//current and next z

            //edge from point at x to point at next x
            glVertex3f (x, y, z);

            if(j>0)
            {
                setNormal(p1x,p1y,p1z,p2x,p2y,p2z,x, y, z);
            }
            else
            {
                p1x=x;
                p1y=y;
                p1z=z;
                p2x=x1;
                p2y=y1;
                p2z=z1;

            }
            glVertex3f (x1, y1, z1);

            //forms quad with next pair of points with incremented theta value
        }
        glEnd();
        x = x1;
        r = r1;
    } //for i

}
void showControlPoints()
{
    glPointSize(5.0);
    glColor3f(1.0, 0.0, 1.0);
    glBegin(GL_POINTS);
    for (int i = 0; i <=L; i++)
        glVertex3fv(&ctrlpoints[i][0]);
    glEnd();
}


using namespace std;
unsigned int ID;

int score=0;
float spt_cutoff = 20;
bool sp_flag= false;
float zoom=4;
int tola[5000][5000];
float tX=0,tY=0,tZ=-8,rX=0,rY=0,rZ=4;
float tZ1=-20,tZ2=-40,tZ3=-60,tZ4=-80,tZ5=-100,tZ6=-120;
float rotX=0,rotY=0,rotZ=0;
float cosX=0,cosY=1,cosZ=0;
float angle=0;
//gluLookAt(	0.0, 14.5, 30.0,0, 4, 0,0, 1.0f, 0.0f);
float eye_x=0,eye_y=14.5,eye_z=30,c_x=0,c_y=4,c_z=0,up_x=0,up_y=1.0f,up_z=0.0f;
float xEye=0.0f,yEye=5.0f,zEye=30.0f;
float cenX=0,cenY=0,cenZ=0,roll=0;
float radius=0;
float theta=0,slope=0;
float speed = 0.0;
float angleBackFrac = 0.2;
bool saheedMinarVisible = false;
//float r[] = {0.1,0.4,0.0,0.9,0.2,0.5,0.0,0.7,0.5,0.0};
//float g[] = {0.2,0.0,0.4,0.5,0.2,0.0,0.3,0.9,0.0,0.2};
//float b[] = {0.4,0.5,0.0,0.7,0.9,0.0,0.1,0.2,0.5,0.0};
int TIME=0;
bool START = false,START1=false;
float torusPosX[7] = {1,-2,3,-4,-2,0,2};
float torusPosY[7] = {2,3,10,6,7,4,1};
bool car_flage= true, car_flage1=false, car_flage2=false,car_flage3=false;

bool rot = false;

///
//light violet
GLfloat mat_ambient1[] = { 1.0,0.8,1.0, 1.0 };
GLfloat mat_diffuse1[] = { 1.0,0.8,1.0, 1.0 };
GLfloat mat_specular1[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat mat_shininess1[] = {60};
// coffee color
GLfloat mat_ambient2[] = { 0.4,0.2,0.2, 1.0 };
GLfloat mat_diffuse2[] = { 0.4,0.2,0.2, 1.0 };
GLfloat mat_specular2[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat mat_shininess2[] = {60};
//white
GLfloat mat_ambient3[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat mat_diffuse3[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat mat_specular3[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat mat_shininess3[] = {60};
// black
GLfloat mat_ambient4[] = { 0.0, 0.0, 0.0, 1.0 };
GLfloat mat_diffuse4[] = { 0.0, 0.0, 0.0, 1.0 };
GLfloat mat_specular4[] = { 0.0, 0.0, 0.0, 1.0 };
GLfloat mat_shininess4[] = {60};
//
GLfloat mat_ambient5[] = { 0.50, 0.0, 0.0, 1.0 };
GLfloat mat_diffuse5[] = { 0.8, 0.0, 0.0, 1.0 };
GLfloat mat_specular5[] = { 1.0, 0.0, 0.0, 1.0 };
GLfloat mat_shininess5[] = {60};

//upper body color
GLfloat mat_ambient6[] = { 0.6, 0.0, 0.2, 1.0 };
GLfloat mat_diffuse6[] = { 0.6, 0.0, 0.2, 1.0 };
GLfloat mat_specular6[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat mat_shininess6[] = {60};

//color for window body
GLfloat mat_ambient7[] = { 1.0, 1.0, 0.8, 1.0 };
GLfloat mat_diffuse7[] = { 1.0, 1.0, 0.8, 1.0 };
GLfloat mat_specular7[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat mat_shininess7[] = {60};

//color for carring body
GLfloat mat_ambient8[] = { 0.8, 1.0, 1.0, 1.0 };
GLfloat mat_diffuse8[] = { 0.8, 1.0, 1.0, 1.0 };
GLfloat mat_specular8[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat mat_shininess8[] = {60};

//color for carring object
GLfloat mat_ambient9[] = { 0.5, 1.0, 0.0, 1.0 };
GLfloat mat_diffuse9[] = { 0.5, 1.0, 0.0, 1.0 };
GLfloat mat_specular9[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat mat_shininess9[] = {60};


GLfloat mat_ambient10[] = { 0.0, 0.4, 0.0, 1.0 };
GLfloat mat_diffuse10[] = { 0.0, 0.4, 0.0, 1.0 };
GLfloat mat_specular10[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat mat_shininess10[] = {60};
//deep violet
GLfloat mat_ambient11[] = { 0.4, 0.0, 0.4, 1.0 };
GLfloat mat_diffuse11[] = { 0.4, 0.0, 0.4, 1.0 };
GLfloat mat_specular11[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat mat_shininess11[] = {60};
// yellow
GLfloat mat_ambient12[] = { 1.0, 0.6, 0.0, 1.0 };
GLfloat mat_diffuse12[] = { 1.0, 0.6, 0.0, 1.0 };
GLfloat mat_specular12[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat mat_shininess12[] = {60};
// off blue
GLfloat mat_ambient13[] = { 0.6, 0.8, 1.0, 1.0 };
GLfloat mat_diffuse13[] = { 0.6, 0.8, 1.0, 1.0 };
GLfloat mat_specular13[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat mat_shininess13[] = {60};
// off fox
GLfloat mat_ambient14[] = { 1.0, 0.6, 0.6, 1.0 };
GLfloat mat_diffuse14[] = { 1.0, 0.6, 0.6, 1.0 };
GLfloat mat_specular14[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat mat_shininess14[] = {60};
// building color fox
GLfloat mat_ambient15[] = { 0.74, 0.71, 0.4, 1.0 };
GLfloat mat_diffuse15[] = { 0.74, 0.71, 0.4, 1.0 };
GLfloat mat_specular15[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat mat_shininess15[] = {60};
// resturen building color
GLfloat mat_ambient16[] = { 0.0, 0.98, 0.60, 1.0 };
GLfloat mat_diffuse16[] = { 0.0, 0.98, 0.60, 1.0 };
GLfloat mat_specular16[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat mat_shininess16[] = {60};



///
//glColor3d(0.4,0.2,0.2);

static void resize(int width, int height)
{
    const float ar = (float) width / (float) height;

    glViewport(0, 0, width, height);
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glFrustum(-ar, ar, -1.0, 1.0, 2.0, 1000.0);

    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
}
static GLfloat v_cube[8][3] =
{
    {0,0,0},
    {0,0,1},
    {0,1,0},
    {0,1,1},

    {1,0,0},
    {1,0,1},
    {1,1,0},
    {1,1,1}
};

static GLubyte c_ind[6][4] =
{
    {0,2,6,4},
    {1,5,7,3},
    {0,4,5,1},
    {2,3,7,6},
    {0,1,3,2},
    {4,6,7,5}
};

static void getNormal3p(GLfloat x1, GLfloat y1, GLfloat z1, GLfloat x2, GLfloat y2, GLfloat z2, GLfloat x3, GLfloat y3, GLfloat z3)
{
    GLfloat Ux, Uy, Uz, Vx, Vy, Vz, Nx, Ny, Nz;

    Ux = x2-x1;
    Uy = y2-y1;
    Uz = z2-z1;

    Vx = x3-x1;
    Vy = y3-y1;
    Vz = z3-z1;

    Nx = Uy*Vz - Uz*Vy;
    Ny = Uz*Vx - Ux*Vz;
    Nz = Ux*Vy - Uy*Vx;

    glNormal3f(Nx,Ny,Nz);
}

//void cube(float R=0.5, float G=0.5, float B=0.5, bool e=false, float alpha=1)
void cube()
{

    glBegin(GL_QUADS);
    for (GLint i = 0; i <6; i++)
    {
        getNormal3p(v_cube[c_ind[i][0]][0], v_cube[c_ind[i][0]][1], v_cube[c_ind[i][0]][2],
                    v_cube[c_ind[i][1]][0], v_cube[c_ind[i][1]][1], v_cube[c_ind[i][1]][2],
                    v_cube[c_ind[i][2]][0], v_cube[c_ind[i][2]][1], v_cube[c_ind[i][2]][2]);

        for (GLint j=0; j<4; j++)
        {
            glVertex3fv(&v_cube[c_ind[i][j]][0]);
        }
    }
    glEnd();
}
void Cube()
{

    glBegin(GL_QUADS);
    for (GLint i = 0; i <6; i++)
    {
        getNormal3p(v_cube[c_ind[i][0]][0], v_cube[c_ind[i][0]][1], v_cube[c_ind[i][0]][2],
                    v_cube[c_ind[i][1]][0], v_cube[c_ind[i][1]][1], v_cube[c_ind[i][1]][2],
                    v_cube[c_ind[i][2]][0], v_cube[c_ind[i][2]][1], v_cube[c_ind[i][2]][2]);

        for (GLint j=0; j<4; j++)
        {
            glVertex3fv(&v_cube[c_ind[i][0]][0]);glTexCoord2f(1,1);
            glVertex3fv(&v_cube[c_ind[i][1]][0]);glTexCoord2f(1,0);
            glVertex3fv(&v_cube[c_ind[i][2]][0]);glTexCoord2f(0,0);
            glVertex3fv(&v_cube[c_ind[i][3]][0]);glTexCoord2f(0,1);
        }
    }
    glEnd();
}


void car_wheel()
{
    glBegin(GL_POLYGON);

    for(int i=0; i<360; i=i+5)
    {
        float x1,y1;
        x1 = 1.5 * cos(((-3.1416)*i)/180);
        y1 = 1.5 * sin(((-3.1416)*i)/180);
        glVertex3f(x1, y1, -.1);
    }
    glEnd();
    glBegin(GL_POLYGON);

    for(int i=0; i<360; i=i+5)
    {
        float x2,y2;
        x2 = 1.5 * cos(((3.1416)*i)/180);
        y2 = 1.5 * sin(((3.1416)*i)/180);
        glVertex3f(x2, y2, .1);
    }

    glEnd();

    glBegin(GL_QUAD_STRIP);

    for(int i=0; i<=360; i=i+5)
    {
        float x3,y3;
        x3= 1.5 * cos(((-3.1416)*i)/180);
        y3 = 1.5 * sin(((-3.1416)*i)/180);
        glVertex3f(x3,y3, -.1);
        glVertex3f(x3,y3, .1);
        glColor3d(0,0,0);
    }
    glEnd();
}

float mx1=2,my1=0.2,mz1=1,ux1=0.5,uy1=0.5,uz1=1;
void car7()
{
    const double t = glutGet(GLUT_ELAPSED_TIME) / 1000.0;
    const double a = t*180.0;
    glPushMatrix();
        glTranslated(0,-1.55,0);
        glScaled(1.5,1.5,1.5);

        /// Main body
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient16);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse16);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular16);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess16);
        glPushMatrix();
            glTranslated(0,0,0);
            glScaled(mx1,my1,mz1);
            cube();
        glPopMatrix();
        ///top
        glPushMatrix();
            glTranslated(0.6,0.2,0);
            glScaled(0.8,0.3,1);
            cube();
        glPopMatrix();
        /// front part
        glPushMatrix();
            /// front shape right
            glPushMatrix();
                glTranslated(0.25,0.2,0);
                glRotated(-45,0,0,1);
                glScaled(0.4,0.45,0.01);
                cube();
            glPopMatrix();
            /// front shape left
            glPushMatrix();
                glTranslated(0.25,0.2,0.99);
                glRotated(-45,0,0,1);
                glScaled(0.4,0.45,0.01);
                cube();
            glPopMatrix();

            /// front shape top cube
            glPushMatrix();
                glTranslated(0.3,0.2,0.0);
                glRotated(45,0,0,1);
                glScaled(0.4,0.01,1);
                cube();
            glPopMatrix();

        glPopMatrix();// front  end

        /// back part
        glPushMatrix();
            glTranslated(0.8,0,0);
            /// front shape right
            glPushMatrix();
                glTranslated(0.25,0.2,0);
                glRotated(-45,0,0,1);
                glScaled(0.4,0.45,0.01);
                cube();
            glPopMatrix();
            /// front shape left
            glPushMatrix();
                glTranslated(0.25,0.2,0.99);
                glRotated(-45,0,0,1);
                glScaled(0.4,0.45,0.01);
                cube();
            glPopMatrix();
        glPopMatrix();
        /// front shape top cube
        glPushMatrix();
            glTranslated(1.4,0.5,0.0);
            glRotated(-45,0,0,1);
            glScaled(0.4,0.01,1);
            cube();
        glPopMatrix();

        /// window
        /// front window 1
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
        ///front 2 side window
        glPushMatrix();
            glTranslated(0.6,0.2,0);
        ///front right
            glPushMatrix();
                glTranslated(0,0,-0.001);
                glScaled(0.25,0.25,0.01);
                cube();
            glPopMatrix();
          /// front left
            glPushMatrix();
                glTranslated(0,0,0.995);
                glScaled(0.25,0.25,0.01);
                cube();
            glPopMatrix();
        glPopMatrix();

        ///back 2 side window
        glPushMatrix();
        glTranslated(1.15,0.2,0);
        ///back right
            glPushMatrix();
                glTranslated(0,0,-0.001);
                glScaled(0.25,0.25,0.01);
                cube();
            glPopMatrix();
          /// back left
            glPushMatrix();
                glTranslated(0,0,0.995);
                glScaled(0.25,0.25,0.01);
                cube();
            glPopMatrix();
        glPopMatrix();

        ///  front light
        glPushMatrix();
        /// light 1
            glPushMatrix();
                glTranslated(0.0,0.08,0.2);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
            /// light 2
            glPushMatrix();
                glTranslated(0.0,0.08,0.8);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
        glPopMatrix();

        ///  front light
        glPushMatrix();
            glTranslated(2,0,0);
        /// light 1
            glPushMatrix();
                glTranslated(0.0,0.08,0.2);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
            /// light 2
            glPushMatrix();
                glTranslated(0.0,0.08,0.8);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
        glPopMatrix();

        /// car wheel

        /// front right
        glPushMatrix();
            //glColor3d(0,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(mx1-0.5,0,0);
            glScalef(0.1,0.1,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();
        ///front left
        glPushMatrix();
            //glColor3d(0,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(mx1-0.5,0,mz1);
            glScalef(0.1,0.1,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();

        ///back right
        glPushMatrix();
            //glColor3d(0,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(mx1-1.5,0,0);
            glScalef(0.1,0.1,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();
        ///front left
        glPushMatrix();
            //glColor3d(0,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(mx1-1.5,0,mz1);
            glScalef(0.1,0.1,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();


    glPopMatrix();

}
void car6()
{
    const double t = glutGet(GLUT_ELAPSED_TIME) / 1000.0;
    const double a = t*180.0;
    glPushMatrix();
        glTranslated(0,-1.55,0);
        glScaled(1.5,1.5,1.5);
        //glRotated(90,0,1,0);

        /// Main body
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient13);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse13);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular13);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess13);
        glPushMatrix();
            glTranslated(0,0,0);
            glScaled(mx1,my1,mz1);
            //glutSolidSphere(1,30,30);
            cube();
        glPopMatrix();
        ///top
        glPushMatrix();
            glTranslated(0.6,0.2,0);
            glScaled(0.8,0.3,1);
            //glutSolidSphere(1,30,30);
            cube();
        glPopMatrix();
        /// front part
        glPushMatrix();
            /// front shape right
            glPushMatrix();
                glTranslated(0.25,0.2,0);
                glRotated(-45,0,0,1);
                glScaled(0.4,0.45,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
            /// front shape left
            glPushMatrix();
                glTranslated(0.25,0.2,0.99);
                glRotated(-45,0,0,1);
                glScaled(0.4,0.45,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();

            /// front shape top cube
            glPushMatrix();
                glTranslated(0.3,0.2,0.0);
                glRotated(45,0,0,1);
                glScaled(0.4,0.01,1);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();

        glPopMatrix();// front  end

        /// back part
        glPushMatrix();
            glTranslated(0.8,0,0);
            /// front shape right
            glPushMatrix();
                glTranslated(0.25,0.2,0);
                glRotated(-45,0,0,1);
                glScaled(0.4,0.45,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
            /// front shape left
            glPushMatrix();
                glTranslated(0.25,0.2,0.99);
                glRotated(-45,0,0,1);
                glScaled(0.4,0.45,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
        glPopMatrix();
        /// front shape top cube
        glPushMatrix();
            glTranslated(1.4,0.5,0.0);
            glRotated(-45,0,0,1);
            glScaled(0.4,0.01,1);
            //glutSolidSphere(1,30,30);
            cube();
        glPopMatrix();

        /// window
        /// front window 1
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
        ///front 2 side window
        glPushMatrix();
            glTranslated(0.6,0.2,0);
        ///front right
            glPushMatrix();
                glTranslated(0,0,-0.001);
                //glRotated(-45,0,0,1);
                glScaled(0.25,0.25,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
          /// front left
            glPushMatrix();
                glTranslated(0,0,0.995);
                //glRotated(-45,0,0,1);
                glScaled(0.25,0.25,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
        glPopMatrix();

        ///back 2 side window
        glPushMatrix();
        glTranslated(1.15,0.2,0);
        ///back right
            glPushMatrix();
                glTranslated(0,0,-0.001);
                //glRotated(-45,0,0,1);
                glScaled(0.25,0.25,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
          /// back left
            glPushMatrix();
                glTranslated(0,0,0.995);
                //glRotated(-45,0,0,1);
                glScaled(0.25,0.25,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
        glPopMatrix();

        ///  front light
        glPushMatrix();
        /// light 1
            glPushMatrix();
                glTranslated(0.0,0.08,0.2);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
            /// light 2
            glPushMatrix();
                glTranslated(0.0,0.08,0.8);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
        glPopMatrix();

        ///  front light
        glPushMatrix();
            glTranslated(2,0,0);
        /// light 1
            glPushMatrix();
                glTranslated(0.0,0.08,0.2);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
            /// light 2
            glPushMatrix();
                glTranslated(0.0,0.08,0.8);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
        glPopMatrix();

        /// car wheel

        /// front right
        glPushMatrix();
            //glColor3d(0,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(mx1-0.5,0,0);
            glScalef(0.1,0.1,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();
        ///front left
        glPushMatrix();
            //glColor3d(0,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(mx1-0.5,0,mz1);
            glScalef(0.1,0.1,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();

        ///back right
        glPushMatrix();
            //glColor3d(0,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(mx1-1.5,0,0);
            glScalef(0.1,0.1,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();
        ///front left
        glPushMatrix();
            //glColor3d(0,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(mx1-1.5,0,mz1);
            glScalef(0.1,0.1,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();


    glPopMatrix();

}

void car5()
{
    const double t = glutGet(GLUT_ELAPSED_TIME) / 1000.0;
    const double a = t*180.0;
    glPushMatrix();
        glTranslated(0,-1.55,0);
        glScaled(1.5,1.5,1.5);
        //glRotated(90,0,1,0);

        /// Main body
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient12);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse12);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular12);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess12);
        glPushMatrix();
            glTranslated(0,0,0);
            glScaled(mx1,my1,mz1);
            //glutSolidSphere(1,30,30);
            cube();
        glPopMatrix();
        ///top
        glPushMatrix();
            glTranslated(0.6,0.2,0);
            glScaled(0.8,0.3,1);
            //glutSolidSphere(1,30,30);
            cube();
        glPopMatrix();
        /// front part
        glPushMatrix();
            /// front shape right
            glPushMatrix();
                glTranslated(0.25,0.2,0);
                glRotated(-45,0,0,1);
                glScaled(0.4,0.45,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
            /// front shape left
            glPushMatrix();
                glTranslated(0.25,0.2,0.99);
                glRotated(-45,0,0,1);
                glScaled(0.4,0.45,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();

            /// front shape top cube
            glPushMatrix();
                glTranslated(0.3,0.2,0.0);
                glRotated(45,0,0,1);
                glScaled(0.4,0.01,1);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();

        glPopMatrix();// front  end

        /// back part
        glPushMatrix();
            glTranslated(0.8,0,0);
            /// front shape right
            glPushMatrix();
                glTranslated(0.25,0.2,0);
                glRotated(-45,0,0,1);
                glScaled(0.4,0.45,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
            /// front shape left
            glPushMatrix();
                glTranslated(0.25,0.2,0.99);
                glRotated(-45,0,0,1);
                glScaled(0.4,0.45,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
        glPopMatrix();
        /// front shape top cube
        glPushMatrix();
            glTranslated(1.4,0.5,0.0);
            glRotated(-45,0,0,1);
            glScaled(0.4,0.01,1);
            //glutSolidSphere(1,30,30);
            cube();
        glPopMatrix();

        /// window
        /// front window 1
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
        ///front 2 side window
        glPushMatrix();
            glTranslated(0.6,0.2,0);
        ///front right
            glPushMatrix();
                glTranslated(0,0,-0.001);
                //glRotated(-45,0,0,1);
                glScaled(0.25,0.25,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
          /// front left
            glPushMatrix();
                glTranslated(0,0,0.995);
                //glRotated(-45,0,0,1);
                glScaled(0.25,0.25,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
        glPopMatrix();

        ///back 2 side window
        glPushMatrix();
        glTranslated(1.15,0.2,0);
        ///back right
            glPushMatrix();
                glTranslated(0,0,-0.001);
                //glRotated(-45,0,0,1);
                glScaled(0.25,0.25,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
          /// back left
            glPushMatrix();
                glTranslated(0,0,0.995);
                //glRotated(-45,0,0,1);
                glScaled(0.25,0.25,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
        glPopMatrix();

        ///  front light
        glPushMatrix();
        /// light 1
            glPushMatrix();
                glTranslated(0.0,0.08,0.2);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
            /// light 2
            glPushMatrix();
                glTranslated(0.0,0.08,0.8);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
        glPopMatrix();

        ///  front light
        glPushMatrix();
            glTranslated(2,0,0);
        /// light 1
            glPushMatrix();
                glTranslated(0.0,0.08,0.2);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
            /// light 2
            glPushMatrix();
                glTranslated(0.0,0.08,0.8);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
        glPopMatrix();

        /// car wheel

        /// front right
        glPushMatrix();
            //glColor3d(0,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(mx1-0.5,0,0);
            glScalef(0.1,0.1,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();
        ///front left
        glPushMatrix();
            //glColor3d(0,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(mx1-0.5,0,mz1);
            glScalef(0.1,0.1,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();

        ///back right
        glPushMatrix();
            //glColor3d(0,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(mx1-1.5,0,0);
            glScalef(0.1,0.1,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();
        ///front left
        glPushMatrix();
            //glColor3d(0,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(mx1-1.5,0,mz1);
            glScalef(0.1,0.1,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();


    glPopMatrix();

}
void car4()
{
    const double t = glutGet(GLUT_ELAPSED_TIME) / 1000.0;
    const double a = t*180.0;
    glPushMatrix();
        glTranslated(0,-1.55,0);
        glScaled(1.5,1.5,1.5);
        //glRotated(90,0,1,0);

        /// Main body
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient11);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse11);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular11);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess11);
        glPushMatrix();
            glTranslated(0,0,0);
            glScaled(mx1,my1,mz1);
            //glutSolidSphere(1,30,30);
            cube();
        glPopMatrix();
        ///top
        glPushMatrix();
            glTranslated(0.6,0.2,0);
            glScaled(0.8,0.3,1);
            //glutSolidSphere(1,30,30);
            cube();
        glPopMatrix();
        /// front part
        glPushMatrix();
            /// front shape right
            glPushMatrix();
                glTranslated(0.25,0.2,0);
                glRotated(-45,0,0,1);
                glScaled(0.4,0.45,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
            /// front shape left
            glPushMatrix();
                glTranslated(0.25,0.2,0.99);
                glRotated(-45,0,0,1);
                glScaled(0.4,0.45,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();

            /// front shape top cube
            glPushMatrix();
                glTranslated(0.3,0.2,0.0);
                glRotated(45,0,0,1);
                glScaled(0.4,0.01,1);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();

        glPopMatrix();// front  end

        /// back part
        glPushMatrix();
            glTranslated(0.8,0,0);
            /// front shape right
            glPushMatrix();
                glTranslated(0.25,0.2,0);
                glRotated(-45,0,0,1);
                glScaled(0.4,0.45,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
            /// front shape left
            glPushMatrix();
                glTranslated(0.25,0.2,0.99);
                glRotated(-45,0,0,1);
                glScaled(0.4,0.45,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
        glPopMatrix();
        /// front shape top cube
        glPushMatrix();
            glTranslated(1.4,0.5,0.0);
            glRotated(-45,0,0,1);
            glScaled(0.4,0.01,1);
            //glutSolidSphere(1,30,30);
            cube();
        glPopMatrix();

        /// window
        /// front window 1
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
        ///front 2 side window
        glPushMatrix();
            glTranslated(0.6,0.2,0);
        ///front right
            glPushMatrix();
                glTranslated(0,0,-0.001);
                //glRotated(-45,0,0,1);
                glScaled(0.25,0.25,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
          /// front left
            glPushMatrix();
                glTranslated(0,0,0.995);
                //glRotated(-45,0,0,1);
                glScaled(0.25,0.25,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
        glPopMatrix();

        ///back 2 side window
        glPushMatrix();
        glTranslated(1.15,0.2,0);
        ///back right
            glPushMatrix();
                glTranslated(0,0,-0.001);
                //glRotated(-45,0,0,1);
                glScaled(0.25,0.25,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
          /// back left
            glPushMatrix();
                glTranslated(0,0,0.995);
                //glRotated(-45,0,0,1);
                glScaled(0.25,0.25,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
        glPopMatrix();

        ///  front light
        glPushMatrix();
        /// light 1
            glPushMatrix();
                glTranslated(0.0,0.08,0.2);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
            /// light 2
            glPushMatrix();
                glTranslated(0.0,0.08,0.8);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
        glPopMatrix();

        ///  front light
        glPushMatrix();
            glTranslated(2,0,0);
        /// light 1
            glPushMatrix();
                glTranslated(0.0,0.08,0.2);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
            /// light 2
            glPushMatrix();
                glTranslated(0.0,0.08,0.8);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
        glPopMatrix();

        /// car wheel

        /// front right
        glPushMatrix();
            //glColor3d(0,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(mx1-0.5,0,0);
            glScalef(0.1,0.1,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();
        ///front left
        glPushMatrix();
            //glColor3d(0,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(mx1-0.5,0,mz1);
            glScalef(0.1,0.1,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();

        ///back right
        glPushMatrix();
            //glColor3d(0,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(mx1-1.5,0,0);
            glScalef(0.1,0.1,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();
        ///front left
        glPushMatrix();
            //glColor3d(0,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(mx1-1.5,0,mz1);
            glScalef(0.1,0.1,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();


    glPopMatrix();

}
void car3()
{
    const double t = glutGet(GLUT_ELAPSED_TIME) / 1000.0;
    const double a = t*180.0;
    glPushMatrix();
        glTranslated(0,-1.55,0);
        glScaled(1.5,1.5,1.5);
        //glRotated(90,0,1,0);

        /// Main body
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient5);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse5);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular5);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess5);
        glPushMatrix();
            glTranslated(0,0,0);
            glScaled(mx1,my1,mz1);
            //glutSolidSphere(1,30,30);
            cube();
        glPopMatrix();
        ///top
        glPushMatrix();
            glTranslated(0.6,0.2,0);
            glScaled(0.8,0.3,1);
            //glutSolidSphere(1,30,30);
            cube();
        glPopMatrix();
        /// front part
        glPushMatrix();
            /// front shape right
            glPushMatrix();
                glTranslated(0.25,0.2,0);
                glRotated(-45,0,0,1);
                glScaled(0.4,0.45,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
            /// front shape left
            glPushMatrix();
                glTranslated(0.25,0.2,0.99);
                glRotated(-45,0,0,1);
                glScaled(0.4,0.45,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();

            /// front shape top cube
            glPushMatrix();
                glTranslated(0.3,0.2,0.0);
                glRotated(45,0,0,1);
                glScaled(0.4,0.01,1);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();

        glPopMatrix();// front  end

        /// back part
        glPushMatrix();
            glTranslated(0.8,0,0);
            /// front shape right
            glPushMatrix();
                glTranslated(0.25,0.2,0);
                glRotated(-45,0,0,1);
                glScaled(0.4,0.45,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
            /// front shape left
            glPushMatrix();
                glTranslated(0.25,0.2,0.99);
                glRotated(-45,0,0,1);
                glScaled(0.4,0.45,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
        glPopMatrix();
        /// front shape top cube
        glPushMatrix();
            glTranslated(1.4,0.5,0.0);
            glRotated(-45,0,0,1);
            glScaled(0.4,0.01,1);
            //glutSolidSphere(1,30,30);
            cube();
        glPopMatrix();

        /// window
        /// front window 1
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
        ///front 2 side window
        glPushMatrix();
            glTranslated(0.6,0.2,0);
        ///front right
            glPushMatrix();
                glTranslated(0,0,-0.001);
                //glRotated(-45,0,0,1);
                glScaled(0.25,0.25,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
          /// front left
            glPushMatrix();
                glTranslated(0,0,0.995);
                //glRotated(-45,0,0,1);
                glScaled(0.25,0.25,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
        glPopMatrix();

        ///back 2 side window
        glPushMatrix();
        glTranslated(1.15,0.2,0);
        ///back right
            glPushMatrix();
                glTranslated(0,0,-0.001);
                //glRotated(-45,0,0,1);
                glScaled(0.25,0.25,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
          /// back left
            glPushMatrix();
                glTranslated(0,0,0.995);
                //glRotated(-45,0,0,1);
                glScaled(0.25,0.25,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
        glPopMatrix();

        ///  front light
        glPushMatrix();
        /// light 1
            glPushMatrix();
                glTranslated(0.0,0.08,0.2);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
            /// light 2
            glPushMatrix();
                glTranslated(0.0,0.08,0.8);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
        glPopMatrix();

        ///  front light
        glPushMatrix();
            glTranslated(2,0,0);
        /// light 1
            glPushMatrix();
                glTranslated(0.0,0.08,0.2);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
            /// light 2
            glPushMatrix();
                glTranslated(0.0,0.08,0.8);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
        glPopMatrix();

        /// car wheel

        /// front right
        glPushMatrix();
            //glColor3d(0,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(mx1-0.5,0,0);
            glScalef(0.1,0.1,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();
        ///front left
        glPushMatrix();
            //glColor3d(0,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(mx1-0.5,0,mz1);
            glScalef(0.1,0.1,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();

        ///back right
        glPushMatrix();
            //glColor3d(0,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(mx1-1.5,0,0);
            glScalef(0.1,0.1,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();
        ///front left
        glPushMatrix();
            //glColor3d(0,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(mx1-1.5,0,mz1);
            glScalef(0.1,0.1,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();


    glPopMatrix();

}
void car2()
{
    const double t = glutGet(GLUT_ELAPSED_TIME) / 1000.0;
    const double a = t*180.0;
    glPushMatrix();
        glTranslated(0,-1.55,0);
        glScaled(1.5,1.5,1.5);
        //glRotated(90,0,1,0);

        /// Main body
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient2);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse2);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular2);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess2);
        glPushMatrix();
            glTranslated(0,0,0);
            glScaled(mx1,my1,mz1);
            //glutSolidSphere(1,30,30);
            cube();
        glPopMatrix();
        ///top
        glPushMatrix();
            glTranslated(0.6,0.2,0);
            glScaled(0.8,0.3,1);
            //glutSolidSphere(1,30,30);
            cube();
        glPopMatrix();
        /// front part
        glPushMatrix();
            /// front shape right
            glPushMatrix();
                glTranslated(0.25,0.2,0);
                glRotated(-45,0,0,1);
                glScaled(0.4,0.45,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
            /// front shape left
            glPushMatrix();
                glTranslated(0.25,0.2,0.99);
                glRotated(-45,0,0,1);
                glScaled(0.4,0.45,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();

            /// front shape top cube
            glPushMatrix();
                glTranslated(0.3,0.2,0.0);
                glRotated(45,0,0,1);
                glScaled(0.4,0.01,1);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();

        glPopMatrix();// front  end

        /// back part
        glPushMatrix();
            glTranslated(0.8,0,0);
            /// front shape right
            glPushMatrix();
                glTranslated(0.25,0.2,0);
                glRotated(-45,0,0,1);
                glScaled(0.4,0.45,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
            /// front shape left
            glPushMatrix();
                glTranslated(0.25,0.2,0.99);
                glRotated(-45,0,0,1);
                glScaled(0.4,0.45,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
        glPopMatrix();
        /// front shape top cube
        glPushMatrix();
            glTranslated(1.4,0.5,0.0);
            glRotated(-45,0,0,1);
            glScaled(0.4,0.01,1);
            //glutSolidSphere(1,30,30);
            cube();
        glPopMatrix();

        /// window
        /// front window 1
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
        ///front 2 side window
        glPushMatrix();
            glTranslated(0.6,0.2,0);
        ///front right
            glPushMatrix();
                glTranslated(0,0,-0.001);
                //glRotated(-45,0,0,1);
                glScaled(0.25,0.25,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
          /// front left
            glPushMatrix();
                glTranslated(0,0,0.995);
                //glRotated(-45,0,0,1);
                glScaled(0.25,0.25,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
        glPopMatrix();

        ///back 2 side window
        glPushMatrix();
        glTranslated(1.15,0.2,0);
        ///back right
            glPushMatrix();
                glTranslated(0,0,-0.001);
                //glRotated(-45,0,0,1);
                glScaled(0.25,0.25,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
          /// back left
            glPushMatrix();
                glTranslated(0,0,0.995);
                //glRotated(-45,0,0,1);
                glScaled(0.25,0.25,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
        glPopMatrix();

        ///  front light
        glPushMatrix();
        /// light 1
            glPushMatrix();
                glTranslated(0.0,0.08,0.2);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
            /// light 2
            glPushMatrix();
                glTranslated(0.0,0.08,0.8);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
        glPopMatrix();

        ///  front light
        glPushMatrix();
            glTranslated(2,0,0);
        /// light 1
            glPushMatrix();
                glTranslated(0.0,0.08,0.2);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
            /// light 2
            glPushMatrix();
                glTranslated(0.0,0.08,0.8);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
        glPopMatrix();

        /// car wheel

        /// front right
        glPushMatrix();
            //glColor3d(0,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(mx1-0.5,0,0);
            glScalef(0.1,0.1,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();
        ///front left
        glPushMatrix();
            //glColor3d(0,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(mx1-0.5,0,mz1);
            glScalef(0.1,0.1,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();

        ///back right
        glPushMatrix();
            //glColor3d(0,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(mx1-1.5,0,0);
            glScalef(0.1,0.1,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();
        ///front left
        glPushMatrix();
            //glColor3d(0,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(mx1-1.5,0,mz1);
            glScalef(0.1,0.1,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();


    glPopMatrix();

}
void car1()
{
    const double t = glutGet(GLUT_ELAPSED_TIME) / 1000.0;
    const double a = t*180.0;
    glPushMatrix();
        glTranslated(0,-1.55,0);
        glScaled(1.5,1.5,1.5);
        //glRotated(90,0,1,0);

        /// Main body
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient10);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse10);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular10);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess10);
        glPushMatrix();
            glTranslated(0,0,0);
            glScaled(mx1,my1,mz1);
            //glutSolidSphere(1,30,30);
            cube();
        glPopMatrix();
        ///top
        glPushMatrix();
            glTranslated(0.6,0.2,0);
            glScaled(0.8,0.3,1);
            //glutSolidSphere(1,30,30);
            cube();
        glPopMatrix();
        /// front part
        glPushMatrix();
            /// front shape right
            glPushMatrix();
                glTranslated(0.25,0.2,0);
                glRotated(-45,0,0,1);
                glScaled(0.4,0.45,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
            /// front shape left
            glPushMatrix();
                glTranslated(0.25,0.2,0.99);
                glRotated(-45,0,0,1);
                glScaled(0.4,0.45,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();

            /// front shape top cube
            glPushMatrix();
                glTranslated(0.3,0.2,0.0);
                glRotated(45,0,0,1);
                glScaled(0.4,0.01,1);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();

        glPopMatrix();// front  end

        /// back part
        glPushMatrix();
            glTranslated(0.8,0,0);
            /// front shape right
            glPushMatrix();
                glTranslated(0.25,0.2,0);
                glRotated(-45,0,0,1);
                glScaled(0.4,0.45,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
            /// front shape left
            glPushMatrix();
                glTranslated(0.25,0.2,0.99);
                glRotated(-45,0,0,1);
                glScaled(0.4,0.45,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
        glPopMatrix();
        /// front shape top cube
        glPushMatrix();
            glTranslated(1.4,0.5,0.0);
            glRotated(-45,0,0,1);
            glScaled(0.4,0.01,1);
            //glutSolidSphere(1,30,30);
            cube();
        glPopMatrix();

        /// window
        /// front window 1
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
        ///front 2 side window
        glPushMatrix();
            glTranslated(0.6,0.2,0);
        ///front right
            glPushMatrix();
                glTranslated(0,0,-0.001);
                //glRotated(-45,0,0,1);
                glScaled(0.25,0.25,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
          /// front left
            glPushMatrix();
                glTranslated(0,0,0.995);
                //glRotated(-45,0,0,1);
                glScaled(0.25,0.25,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
        glPopMatrix();

        ///back 2 side window
        glPushMatrix();
        glTranslated(1.15,0.2,0);
        ///back right
            glPushMatrix();
                glTranslated(0,0,-0.001);
                //glRotated(-45,0,0,1);
                glScaled(0.25,0.25,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
          /// back left
            glPushMatrix();
                glTranslated(0,0,0.995);
                //glRotated(-45,0,0,1);
                glScaled(0.25,0.25,0.01);
                //glutSolidSphere(1,30,30);
                cube();
            glPopMatrix();
        glPopMatrix();

        ///  front light
        glPushMatrix();
        /// light 1
            glPushMatrix();
                glTranslated(0.0,0.08,0.2);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
            /// light 2
            glPushMatrix();
                glTranslated(0.0,0.08,0.8);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
        glPopMatrix();

        ///  front light
        glPushMatrix();
            glTranslated(2,0,0);
        /// light 1
            glPushMatrix();
                glTranslated(0.0,0.08,0.2);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
            /// light 2
            glPushMatrix();
                glTranslated(0.0,0.08,0.8);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
        glPopMatrix();

        /// car wheel

        /// front right
        glPushMatrix();
            //glColor3d(0,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(mx1-0.5,0,0);
            glScalef(0.1,0.1,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();
        ///front left
        glPushMatrix();
            //glColor3d(0,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(mx1-0.5,0,mz1);
            glScalef(0.1,0.1,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();

        ///back right
        glPushMatrix();
            //glColor3d(0,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(mx1-1.5,0,0);
            glScalef(0.1,0.1,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();
        ///front left
        glPushMatrix();
            //glColor3d(0,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(mx1-1.5,0,mz1);
            glScalef(0.1,0.1,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();

    glPopMatrix();

}

float mx=3,my=0.3,mz=1,ux=0.5,uy=0.5,uz=1;
void car()
{

    const double t = glutGet(GLUT_ELAPSED_TIME) / 1000.0;
    const double a = t*180.0;
    glPushMatrix();
        glTranslated(0,-1.55,0);
        glScaled(1.5,1.5,1.5);
        /// Main body
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient6);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse6);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular6);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess6);
        glPushMatrix();
            glTranslated(0,0,0);
            glScaled(mx,my,mz);
            cube();
        glPopMatrix();

        glColor3d(0,0,0);

        ///upper body///
        glPushMatrix();
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient6);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse6);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular6);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess6);
            glTranslatef(mx-ux-0.2,my,0);
            glScalef(ux,uy,uz);
            cube();
        glPopMatrix();

          ///front window///
        glPushMatrix();
            //glColor3d(1,1,0.8);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient7);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse7);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular7);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess7);

            glTranslatef(mx-0.2,my+0.05,0.1);
            glScalef(0.0097,uy*0.8,uz*0.8);
            cube();
        glPopMatrix();

        /// back window///
        glPushMatrix();
            //glColor3d(1,1,0.8);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient7);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse7);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular7);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess7);
            glTranslatef(mx-ux-0.25,my+uy-0.25,0.4);
            glScalef(0.0097,0.2,0.2);
            cube();
        glPopMatrix();

        /// right window///
        glPushMatrix();
            //glColor3d(1,1,0.8);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient7);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse7);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular7);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess7);
            glTranslatef(mx-ux-0.1,my+0.1,-0.01);
            glScalef(ux*0.6,uy*0.6,0.0097);
            cube();
        glPopMatrix();

        ///left window///
        glPushMatrix();
            //glColor3d(1,1,0.8);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient7);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse7);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular7);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess7);
            glTranslatef(mx-ux-0.1,my+0.1,1.01);
            glScalef(ux*0.6,uy*0.6,0.0097);
            cube();
        glPopMatrix();

        ///  back light
        glPushMatrix();
        /// light 1
            glPushMatrix();
                glTranslated(0.0,0.13,0.2);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
            /// light 2
            glPushMatrix();
                glTranslated(0.0,0.13,0.8);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
        glPopMatrix();

        ///  front light
        glPushMatrix();
            glTranslated(3,0,0);
        /// light 1
            glPushMatrix();
                glTranslated(0.0,0.08,0.2);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
            /// light 2
            glPushMatrix();
                glTranslated(0.0,0.08,0.8);
                glutSolidSphere(0.05,20,20);
            glPopMatrix();
        glPopMatrix();



       ///back right wheel 1///
        glPushMatrix();
            //glColor3d(0,0,1);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(.5,0,0);
            //glScalef(1.8,0.8,0.097);
            glScalef(0.15,0.15,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();

        ///back left wheel 2///
        glPushMatrix();
            //glColor3d(0,0,1);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(.5,0,mz);
            //glScalef(1.8,0.8,0.097);
            glScalef(0.15,0.15,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();

        ///front right wheel 1///
        glPushMatrix();
            //glColor3d(0,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(mx-0.5,0,0);
            glScalef(0.15,0.15,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();

        ///front left wheel 2///
        glPushMatrix();
            //glColor3d(0,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslatef(mx-0.5,0,mz);
            glScalef(0.15,0.15,1);
            glRotated(-a,0,0,1);
            car_wheel();
        glPopMatrix();

        ///front carring part///
        glPushMatrix();
            //glColor3d(0.8,1,1
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient8);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse8);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular8);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess8);
            glTranslatef(mx-ux-0.3,my,0);
            glScalef(0.0097,uy*0.5,uz);
            cube();
        glPopMatrix();

        ///back carring box///
        //glEnable(GL_TEXTURE_2D);
        glPushMatrix();
            //glColor3d(0.8,1,1);
            glTranslatef(0,my,0);
            glScalef(0.0097,uy*0.5,uz);
            cube();
        glPopMatrix();
        //glDisable(GL_TEXTURE_2D);

        ///right carring part///
        glPushMatrix();
            //glColor3d(0.8,1,1);
            glTranslatef(0,my,0);
            glScalef(mx-ux-0.3,uy*0.5,0.0097);
            cube();
        glPopMatrix();

        ///left carring part///
        glPushMatrix();
            //glColor3d(0.8,1,1);
            glTranslatef(0,my,mz);
            glScalef(mx-ux-0.3,uy*0.5,0.0097);
            cube();
        glPopMatrix();

        ///carring object 1
        glPushMatrix();
            //glColor3d(0.5,1,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient9);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse9);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular9);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess9);
            glTranslated(1.1,my,0.3);
            glScaled(1,0.1,0.2);
            glutSolidSphere(1,30,30);
        glPopMatrix();

        ///carring object 2
        glPushMatrix();
            //glColor3d(0.5,1,0);
            glTranslated(1.1,my,0.8);
            glScaled(1,0.1,0.2);
            glutSolidSphere(1,30,30);
        glPopMatrix();



    glPopMatrix();
}

void car_station()
{
        glPushMatrix();
            glTranslated(-2,0,0);
            glPushMatrix();
                glTranslated(0,1.55,0);
                glScaled(0.05,1.45,0.05);
                Cube();
            glPopMatrix();
            glPushMatrix();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient11);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse11);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular11);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess11);
                glTranslated(0,3.88,0);
                glScaled(0.75,0.750,0.5);
                glRotated(-90,0,0,1);
                station_top();
            glPopMatrix();
        glPopMatrix();

        glPushMatrix();
            ///1
            glPushMatrix();
                glTranslated(1,2.1 ,0);
                glRotated(-90,0,1,0);
                glScaled(0.3,0.3,0.3);
                car();
            glPopMatrix();

            ///2
            glPushMatrix();
                glTranslated(0,2.1 ,0);
                glRotated(-90,0,1,0);
                glScaled(0.3,0.3,0.3);
                car();
            glPopMatrix();

            ///3
            glPushMatrix();
                glTranslated(-1,2.1 ,0);
                glRotated(-90,0,1,0);
                glScaled(0.3,0.3,0.3);
                car7();
            glPopMatrix();
            ///4
            glPushMatrix();
                glTranslated(-2,2.1 ,0);
                glRotated(-90,0,1,0);
                glScaled(0.3,0.3,0.3);
                car7();
            glPopMatrix();
            ///5
            glPushMatrix();
                glTranslated(-3,2.1 ,0);
                glRotated(-90,0,1,0);
                glScaled(0.3,0.3,0.3);
                car6();
            glPopMatrix();
            ///6
            glPushMatrix();
                glTranslated(-4,2.1 ,0);
                glRotated(-90,0,1,0);
                glScaled(0.3,0.3,0.3);
                car6();
            glPopMatrix();

            ///7
            glPushMatrix();
                glTranslated(-5,2.1 ,0);
                glRotated(-90,0,1,0);
                glScaled(0.3,0.3,0.3);
                car5();
            glPopMatrix();

        glPopMatrix();



        ///
        glPushMatrix();
            glTranslated(0,0,-2);
            ///1
            glPushMatrix();
                glTranslated(1,2.1 ,0);
                glRotated(-90,0,1,0);
                glScaled(0.3,0.3,0.3);
                car4();
            glPopMatrix();

            ///2
            glPushMatrix();
                glTranslated(0,2.1 ,0);
                glRotated(-90,0,1,0);
                glScaled(0.3,0.3,0.3);
                car4();
            glPopMatrix();

            ///3
            glPushMatrix();
                glTranslated(-1,2.1 ,0);
                glRotated(-90,0,1,0);
                glScaled(0.3,0.3,0.3);
                car3();
            glPopMatrix();
            ///4
            glPushMatrix();
                glTranslated(-2,2.1 ,0);
                glRotated(-90,0,1,0);
                glScaled(0.3,0.3,0.3);
                car3();
            glPopMatrix();
            ///5
            glPushMatrix();
                glTranslated(-3,2.1 ,0);
                glRotated(-90,0,1,0);
                glScaled(0.3,0.3,0.3);
                car2();
            glPopMatrix();
            ///6
            glPushMatrix();
                glTranslated(-4,2.1 ,0);
                glRotated(-90,0,1,0);
                glScaled(0.3,0.3,0.3);
                car2();
            glPopMatrix();

            ///7
            glPushMatrix();
                glTranslated(-5,2.1 ,0);
                glRotated(-90,0,1,0);
                glScaled(0.3,0.3,0.3);
                car1();
            glPopMatrix();

        glPopMatrix();
}

void Man()
{
    glPushMatrix();
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient2);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse2);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular2);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess2);
        ///body
        glPushMatrix();
            glTranslated(0,1.75,0);
            glScaled(0.05,0.2,0.14);
            cube();
        glPopMatrix();

        ///leg left
        glPushMatrix();
            glTranslated(0,1.55,0);
            glScaled(0.025,0.2,0.025);
            cube();
        glPopMatrix();

        ///leg right
        glPushMatrix();
            glTranslated(0,1.55,0.115);
            glScaled(0.025,0.2,0.025);
            cube();
        glPopMatrix();
        /// head connector
        glPushMatrix();
            glTranslated(0,1.95,0.065);
            glScaled(0.025,0.05,0.025);
            cube();
        glPopMatrix();

        ///head
        glPushMatrix();
            glTranslated(0,2.025,0.065);
            glScaled(0.5,0.5,0.5);
            //cube();
            glutSolidSphere(0.1,20,20);
        glPopMatrix();

        ///Eye left
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
        glPushMatrix();
            glTranslated(0.02,2.05,0.090);
            glScaled(0.15,0.15,0.15);
            //cube();
            glutSolidSphere(0.1,20,20);
        glPopMatrix();
        ///eye right
        glPushMatrix();
            glTranslated(0.02,2.05,0.040);
            glScaled(0.15,0.15,0.15);
            //cube();
            glutSolidSphere(0.1,20,20);
        glPopMatrix();
        /// mouth
        glPushMatrix();
            glTranslated(0.035,2.02,0.040);
            glScaled(0.015,0.015,0.05);
            cube();
        glPopMatrix();
        /// hand left
        glPushMatrix();
            glTranslated(0.02,1.95,-0.015);
            //glRotated(-45,0,1,0);
            glRotated(-135,0,0,1);
            glScaled(0.015,0.15,0.015);
            cube();
        glPopMatrix();

        /// hand right
        glPushMatrix();
            glTranslated(0.02,1.95,0.14);
            //glRotated(-45,0,1,0);
            glRotated(-135,0,0,1);
            glScaled(0.015,0.15,0.015);
            cube();
        glPopMatrix();

    glPopMatrix();

}

void drawShohidMinar(){
    ///base 1
    //glColor3d(0.4,0.2,0.2);
	glPushMatrix();
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient2);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse2);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular2);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess2);
        glTranslated(0,1.55,0);
        glScaled(2,0.05,1.5);
        glutSolidCube(1);
        //cube();
    glPopMatrix();
  /// top of base 1
    //base 2
    //glColor3d(0.4,0.2,0.2);
	glPushMatrix();
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient2);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse2);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular2);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess2);
        glTranslated(0,1.6,0);
        glScaled(1.9,0.05,1.4);
        glutSolidCube(1);
    glPopMatrix();
    /// top of base 2
    // base 3
    //glColor3d(0.4,0.2,0.2);
	glPushMatrix();
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient2);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse2);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular2);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess2);
        glTranslated(0,1.65,0);
        glScaled(1.8,0.05,1.3);
        glutSolidCube(1);
    glPopMatrix();

     /// pataton

    //glColor3d(1,1,1);
	glPushMatrix();
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
        glTranslated(0,1.68,-0.4);
        glScaled(0.5,0.02,0.08);
        glutSolidCube(1);
    glPopMatrix();

    /// Piller

	glPushMatrix();
        //glTranslated(0,1.99,-0.4);
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
        glTranslated(0,2.04,-0.4);
        glScaled(0.06,0.7,0.04);
        glutSolidCube(1);
    glPopMatrix();

   /// ROD
   glPushMatrix();// back side 3 rod

        //glColor3d(0,0,0);

        glPushMatrix();
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslated(0.07,2.04,-0.4);
            glScaled(0.003,0.7,0.003);
            glutSolidCube(1);
        glPopMatrix();

        glPushMatrix();
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslated(0.11,2.04,-0.4);
            glScaled(0.003,0.7,0.003);
            glutSolidCube(1);
        glPopMatrix();

        glPushMatrix();
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslated(0.15,2.04,-0.4);
            glScaled(0.003,0.7,0.003);
            glutSolidCube(1);
        glPopMatrix();

    glPopMatrix();

    /// ROD
    glPushMatrix();// front side rod
    glTranslated(-0.22,0,0);
        glPushMatrix();
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslated(0.07,2.04,-0.4);
            glScaled(0.003,0.7,0.003);
            glutSolidCube(1);
        glPopMatrix();

        glPushMatrix();
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslated(0.11,2.04,-0.4);
            glScaled(0.003,0.7,0.003);
            glutSolidCube(1);
        glPopMatrix();

        glPushMatrix();
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glTranslated(0.15,2.04,-0.4);
            glScaled(0.003,0.7,0.003);
            glutSolidCube(1);
        glPopMatrix();
    glPopMatrix();

   /* ///Horizontal rod
        glPushMatrix();
            glTranslated(2.2,0,-0.1);
            glScaled(4.2,1,1);
                glColor3d(0,0,0);
                glPushMatrix();
                    glTranslated(-0.528,1.85,-0.3);
                    glScaled(0.1,0.003,0.003);
                    glutSolidCube(1);
                glPopMatrix();

                glPushMatrix();
                    glTranslated(-0.528,2.02,-0.3);
                    glScaled(0.1,0.003,0.003);
                    glutSolidCube(1);
                glPopMatrix();

                glPushMatrix();
                    glTranslated(-0.528,2.18,-0.3);
                    glScaled(0.1,0.003,0.003);
                    glutSolidCube(1);
                glPopMatrix();
                glColor3d(1,1,1);
            glPopMatrix();
*/
   /// ROD END
    /// front side piller of middele
    glColor3d(1,1,1);
    glPushMatrix();
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
        glTranslated(-0.22,2.04,-0.4);
        glScaled(0.06,0.7,0.04);
        glutSolidCube(1);
    glPopMatrix();
    /// back side piller of middele
    glPushMatrix();
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
        glTranslated(0.22,2.04,-0.4);
        glScaled(0.06,0.7,0.04);
        glutSolidCube(1);
    glPopMatrix();

    /// Uporer piller

    glPushMatrix();
        glTranslated(0,0.743,-1.424);
        glRotated(45,1,0,0);

        //glColor3d(1,0,1);
        //middle
        glPushMatrix();
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
            glTranslated(0,2.04,-0.4);
            glScaled(0.06,0.3,0.04);
            glutSolidCube(1);
        glPopMatrix();
        //glColor3d(1,1,1);

        glPushMatrix();

            glTranslated(-0.22,2.04,-0.4);
            glScaled(0.06,0.3,0.04);
            glutSolidCube(1);
        glPopMatrix();

        glPushMatrix();
            glTranslated(0.22,2.04,-0.4);
            glScaled(0.06,0.3,0.04);
            glutSolidCube(1);
        glPopMatrix();

        glPushMatrix();
            glTranslated(0,2.20,-0.4);
            glScaled(0.5,0.04,0.04);
            glutSolidCube(1);
        glPopMatrix();

        /// ROD back left
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);

        glPushMatrix();
            glPushMatrix();
                glTranslated(0.07,2.04,-0.4);
                glScaled(0.003,0.277,0.003);
                glutSolidCube(1);
            glPopMatrix();

            glPushMatrix();
                glTranslated(0.11,2.04,-0.4);
                glScaled(0.003,0.277,0.003);
                glutSolidCube(1);
            glPopMatrix();

            glPushMatrix();
                glTranslated(0.15,2.04,-0.4);
                glScaled(0.003,0.277,0.003);
                glutSolidCube(1);
            glPopMatrix();
        glPopMatrix();

        ///

        ///rod front left
        glPushMatrix();
            glTranslated(-0.22,0,0);
                glPushMatrix();
                    glTranslated(0.07,2.04,-0.4);
                    glScaled(0.003,0.277,0.003);
                    glutSolidCube(1);
                glPopMatrix();

                glPushMatrix();
                    glTranslated(0.11,2.04,-0.4);
                    glScaled(0.003,0.277,0.003);
                    glutSolidCube(1);
                glPopMatrix();

                glPushMatrix();
                    glTranslated(0.15,2.04,-0.4);
                    glScaled(0.003,0.277,0.003);
                    glutSolidCube(1);
                glPopMatrix();
        glPopMatrix();
        /// ROD END

       /* ///Horizontal rod
        glPushMatrix();
            glTranslated(2.2,0,-0.1);
            glScaled(4.2,1,1);
                glColor3d(0,0,0);
                glPushMatrix();
                    glTranslated(-0.528,1.85,-0.3);
                    glScaled(0.1,0.003,0.003);
                    glutSolidCube(1);
                glPopMatrix();

                glPushMatrix();
                    glTranslated(-0.528,2,-0.3);
                    glScaled(0.1,0.003,0.003);
                    glutSolidCube(1);
                glPopMatrix();

                glPushMatrix();
                    glTranslated(-0.528,2.15,-0.3);
                    glScaled(0.1,0.003,0.003);
                    glutSolidCube(1);
                glPopMatrix();
                glColor3d(1,1,1);
        glPopMatrix();*/

    glPopMatrix();


    /// pasher piller left 1

    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
    glPushMatrix();
        glTranslated(0.1,0,-0.4);
        glRotated(45,0,1,0);

        glPushMatrix();
            glTranslated(-0.605,2.02,-0.3);
            glScaled(0.045,0.65,0.03);
            glutSolidCube(1);
        glPopMatrix();

        glPushMatrix();
            glTranslated(-0.45,2.02,-0.3);
            glScaled(0.045,0.65,0.03);
            glutSolidCube(1);
        glPopMatrix();

        glPushMatrix();
            glTranslated(-0.528,2.35,-0.3);
            glScaled(0.199,0.04,0.03);
            glutSolidCube(1);
        glPopMatrix();

        glPushMatrix();
            glTranslated(-0.528,1.68,-0.3);
            glScaled(0.199,0.02,0.06);
            glutSolidCube(1);
        glPopMatrix();

        /// ROD
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
        glPushMatrix();
            glTranslated(-0.64,-0.05,0.1);
            glScaled(1,1.02,1);
            glPushMatrix();
                glTranslated(0.078,2.02,-0.4);
                glScaled(0.003,0.56,0.003);
                glutSolidCube(1);
            glPopMatrix();

            glPushMatrix();
                glTranslated(0.11,2.02,-0.4);
                glScaled(0.003,0.56,0.003);
                glutSolidCube(1);
            glPopMatrix();

            glPushMatrix();
                glTranslated(0.145,2.02,-0.4);
                glScaled(0.003,0.56,0.003);
                glutSolidCube(1);
            glPopMatrix();
        glPopMatrix();

        ///
        //glColor3d(1,1,1);

        /*///Horizontal rod
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
            glPushMatrix();
                glTranslated(-0.528,1.85,-0.3);
                glScaled(0.1,0.003,0.003);
                glutSolidCube(1);
            glPopMatrix();

            glPushMatrix();
                glTranslated(-0.528,2,-0.3);
                glScaled(0.1,0.003,0.003);
                glutSolidCube(1);
            glPopMatrix();

            glPushMatrix();
                glTranslated(-0.528,2.15,-0.3);
                glScaled(0.1,0.003,0.003);
                glutSolidCube(1);
            glPopMatrix();
            glColor3d(1,1,1);*/

    glPopMatrix();

    /// pasher piller left 2
    glPushMatrix();
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);

        glTranslated(0.65,0,0.3);
        glRotated(-45,0,1,0);

        glPushMatrix();
            glTranslated(-0.605,2.02,-0.3);
            glScaled(0.045,0.65,0.03);
            glutSolidCube(1);
        glPopMatrix();

        glPushMatrix();
            glTranslated(-0.45,2.02,-0.3);
            glScaled(0.045,0.65,0.03);
            glutSolidCube(1);
        glPopMatrix();

        glPushMatrix();
            glTranslated(-0.528,2.35,-0.3);
            glScaled(0.199,0.04,0.03);
            glutSolidCube(1);
        glPopMatrix();

        glPushMatrix();
            glTranslated(-0.528,1.68,-0.3);
            glScaled(0.199,0.02,0.06);
            glutSolidCube(1);
        glPopMatrix();

        ///ROD
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);

        glPushMatrix();
            glTranslated(-0.64,-0.05,0.1);
            glScaled(1,1.02,1);
            glPushMatrix();
                glTranslated(0.078,2.02,-0.4);
                glScaled(0.003,0.56,0.003);
                glutSolidCube(1);
            glPopMatrix();

            glPushMatrix();
                glTranslated(0.11,2.02,-0.4);
                glScaled(0.003,0.56,0.003);
                glutSolidCube(1);
            glPopMatrix();

            glPushMatrix();
                glTranslated(0.145,2.02,-0.4);
                glScaled(0.003,0.56,0.003);
                glutSolidCube(1);
            glPopMatrix();
        glPopMatrix();
        glColor3d(1,1,1);
        ///ROD

      /*  ///Horizontal rod
            glColor3d(0,0,0);
            glPushMatrix();
                glTranslated(-0.528,1.85,-0.3);
                glScaled(0.1,0.003,0.003);
                glutSolidCube(1);
            glPopMatrix();

            glPushMatrix();
                glTranslated(-0.528,2,-0.3);
                glScaled(0.1,0.003,0.003);
                glutSolidCube(1);
            glPopMatrix();

            glPushMatrix();
                glTranslated(-0.528,2.15,-0.3);
                glScaled(0.1,0.003,0.003);
                glutSolidCube(1);
            glPopMatrix();
            glColor3d(1,1,1);*/

    glPopMatrix();


    ///Choto pillers

    glPushMatrix();
        /// pasher piller left 1
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);

        glTranslated(0.06,0,0.14);
        glPushMatrix();

            glTranslated(-0.2,0,-0.31);
            glRotated(45,0,1,0);

            glPushMatrix();
                glTranslated(-0.605,1.88,-0.3);
                glScaled(0.045,0.4,0.03);
                glutSolidCube(1);
            glPopMatrix();

            glPushMatrix();
                glTranslated(-0.45,1.88,-0.3);
                glScaled(0.045,0.4,0.03);
                glutSolidCube(1);
            glPopMatrix();

            glPushMatrix();
                glTranslated(-0.528,2.08,-0.3);
                glScaled(0.2,0.04,0.03);
                glutSolidCube(1);
            glPopMatrix();

            glPushMatrix();
                glTranslated(-0.528,1.68,-0.3);
                glScaled(0.199,0.02,0.06);
                glutSolidCube(1);
            glPopMatrix();

        ///ROD
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);

        glPushMatrix();
        glTranslated(-0.641,0.43,0.1);
        glScaled(1,0.73,1);
            glPushMatrix();
                glTranslated(0.078,1.99,-0.4);
                glScaled(0.003,0.56,0.003);
                glutSolidCube(1);
            glPopMatrix();

            glPushMatrix();
                glTranslated(0.11,1.99,-0.4);
                glScaled(0.003,0.56,0.003);
                glutSolidCube(1);
            glPopMatrix();

            glPushMatrix();
                glTranslated(0.145,1.99,-0.4);
                glScaled(0.003,0.56,0.003);
                glutSolidCube(1);
            glPopMatrix();
        glPopMatrix();

     /*   ///Horizontal rod
            glColor3d(0,0,0);
            glPushMatrix();
                glTranslated(-0.528,1.8,-0.3);
                glScaled(0.1,0.003,0.003);
                glutSolidCube(1);
            glPopMatrix();

            glPushMatrix();
                glTranslated(-0.528,1.96,-0.3);
                glScaled(0.1,0.003,0.003);
                glutSolidCube(1);
            glPopMatrix();
            glColor3d(1,1,1);
        ///ROD
        */

        glPopMatrix();

        /// pasher piller left 2
        glPushMatrix();
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);

            glTranslated(0.83,0,0.39);
            glRotated(-45,0,1,0);

            glPushMatrix();
                glTranslated(-0.605,1.88,-0.3);
                glScaled(0.045,0.4,0.03);
                glutSolidCube(1);
            glPopMatrix();

            glPushMatrix();
                glTranslated(-0.45,1.88,-0.3);
                glScaled(0.045,0.4,0.03);
                glutSolidCube(1);
            glPopMatrix();

            glPushMatrix();
                glTranslated(-0.528,2.1,-0.3);
                glScaled(0.199,0.04,0.03);
                glutSolidCube(1);
            glPopMatrix();

            glPushMatrix();
                glTranslated(-0.528,1.68,-0.3);
                glScaled(0.199,0.02,0.06);
                glutSolidCube(1);
            glPopMatrix();




        ///ROD
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);

        glPushMatrix();
        glTranslated(-0.641,0.43,0.1);
        glScaled(1,0.73,1);
            glPushMatrix();
                glTranslated(0.078,1.99,-0.4);
                glScaled(0.003,0.56,0.003);
                glutSolidCube(1);
            glPopMatrix();

            glPushMatrix();
                glTranslated(0.11,1.99,-0.4);
                glScaled(0.003,0.56,0.003);
                glutSolidCube(1);
            glPopMatrix();

            glPushMatrix();
                glTranslated(0.145,1.99,-0.4);
                glScaled(0.003,0.56,0.003);
                glutSolidCube(1);
            glPopMatrix();
        glPopMatrix();
        glColor3d(1,1,1);
        ///ROD

        glPopMatrix();


    glPopMatrix();

    /// Circle

    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient5);
    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse5);
    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular5);
    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess5);

    glPushMatrix();
        glTranslated(0,2.1,-0.44);
        glScaled(0.35,0.35,0.01);
        glutSolidSphere(1,50,50);
    glPopMatrix();

    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
    glPushMatrix();
        glTranslated(-0.18,1.9,-0.45);
        glScaled(0.01,0.5,0.01);
        glutSolidCube(1);
    glPopMatrix();

    glPushMatrix();
        glTranslated(0.18,1.9,-0.45);
        glScaled(0.01,0.5,0.01);
        glutSolidCube(1);
    glPopMatrix();


}


void soheedMinarEnv(){
    /// Ground

    //glColor3d(0,0.5,0.1);
    glPushMatrix();
        //glColor3d(1,0.8,1);
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient1);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse1);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular1);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess1);
        glTranslated(0,0,0);
        glScaled(EN_SIZE*2,0.3,EN_SIZE*2);
        glutSolidCube(1);
        //Cube();
    glPopMatrix();

    //tunnel
    glPushMatrix();
    glScaled(1,1,0.6);
        glPushMatrix();
            //glColor3d(1,0.8,1);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient14);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse14);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular14);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess14);
            glTranslated(0,-2.0,3);
            glRotated(110,0,0,1);
            glRotated(90,0,1,0);
            glScaled(.3,.3,.3);
            tunnelBezier();
        glPopMatrix();
    glPopMatrix();


    glPushMatrix();
        glTranslated(-9,-2.7,-5);
        glRotated(65,0,1,0);
        //glRotated(15,0,1,0);
        glScaled(2,2,2);
        drawShohidMinar();
    glPopMatrix();

    glPushMatrix();
        glTranslated(9,-2.7,-5);
        glRotated(-65,0,1,0);
        //glRotated(15,0,1,0);
        glScaled(2,2,2);
        drawShohidMinar();
    glPopMatrix();
}

void football_ground(){
    /// people
    //1
    glPushMatrix();
        glTranslated(0,0,-0.5);
        //glScaled(2,1,2);
        Man();
    glPopMatrix();
    //2
    glPushMatrix();
        glTranslated(0,0,-1);
        //glScaled(2,1,2);
        Man();
    glPopMatrix();
        //3
    glPushMatrix();
        glTranslated(0,0,-1.5);
        //glScaled(2,1,2);
        Man();
    glPopMatrix();
        //4
    glPushMatrix();
        glTranslated(0,0,-2);
        //glScaled(2,1,2);
        Man();
    glPopMatrix();
    //5
    glPushMatrix();
        glTranslated(-1,0,-0.5);
        //glScaled(2,1,2);
        Man();
    glPopMatrix();
    //6
    glPushMatrix();
        glTranslated(-1,0,-1);
        //glScaled(2,1,2);
        Man();
    glPopMatrix();
        //7
    glPushMatrix();
        glTranslated(-1,0,-1.5);
        //glScaled(2,1,2);
        Man();
    glPopMatrix();
        //8
    glPushMatrix();
        glTranslated(-1,0,-2);
        //glScaled(2,1,2);
        Man();
    glPopMatrix();
    ///goal kepper
    glPushMatrix();
        glTranslated(-1.5,0,-1.5);
        //glScaled(2,1,2);
        Man();
    glPopMatrix();


	/// ground
	glPushMatrix();
        glTranslated(-1.5,0,-2.8);
        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D,1);
        glPushMatrix();
            //glColor3d(1,1.0,1.0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
            //glTranslated(-1.5,1.55,-2.5);
            glTranslated(0,1.55,0);
            glScaled(4,0.05,3);
            //glutSolidCube(1);
            Cube();
        glPopMatrix();
        glDisable(GL_TEXTURE_2D);

        ///goalbar front
        ///f1
        glPushMatrix();
            glTranslated(0,1.55,1.2);
            glScaled(0.02,0.5,0.02);
            cube();
        glPopMatrix();

        ///f2
        glPushMatrix();
            glTranslated(0,1.55,1.8);
            glScaled(0.02,0.5,0.02);
            cube();
        glPopMatrix();

        ///f3
        glPushMatrix();
            glTranslated(0,2.05,1.2);
            glScaled(0.02,0.02,0.6);
            //glutSolidCube(1);
            cube();
        glPopMatrix();

        ///goalbar back
        ///b1
        glPushMatrix();
            //glColor3d(1,1.0,1.0);
            //glTranslated(-1.5,1.55,-2.5);
            glTranslated(4,1.55,1.2);
            glScaled(0.02,0.5,0.02);
            //glutSolidCube(1);
            cube();
        glPopMatrix();

        ///b2
        glPushMatrix();
            //glColor3d(1,1.0,1.0);
            //glTranslated(-1.5,1.55,-2.5);
            glTranslated(4,1.55,1.8);
            glScaled(0.02,0.5,0.02);
            //glutSolidCube(1);
            cube();
        glPopMatrix();

        ///b3
        glPushMatrix();
            //glColor3d(1,1.0,1.0);
            //glTranslated(-1.5,1.55,-2.5);
            glTranslated(4,2.05,1.2);
            glScaled(0.02,0.02,0.6);
            //glutSolidCube(1);
            cube();
        glPopMatrix();

    glPopMatrix();//finished ground

    /// stand
    glPushMatrix();
	    //glColor3d(0,0.0,0.0);
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient4);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse4);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular4);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess4);
        //glTranslated(-1.5,1.55,-2.5);
        glTranslated(0,1.55,0.5);
        glScaled(0.05,1,0.05);
        //glutSolidCube(1);
        cube();
    glPopMatrix();

    ///board
    glEnable(GL_TEXTURE_2D);
	glBindTexture(GL_TEXTURE_2D,2);
    glPushMatrix();
	    //glColor3d(1.0,1.0,1.0);
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
        //glTranslated(-1.5,1.55,-2.5);
        glTranslated(-0.5,2.55,0.5);
        glScaled(1,1,0.05);
        //glRotated(80,0,1,0);
        //glutSolidCube(1);
        Cube();
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);


}



void road(){
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,4);
    //glColor3d(1,1,1);
    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
    Cube();
    glDisable(GL_TEXTURE_2D);

}

void banch1(){
        //1
    glPushMatrix();
        glTranslated(0.0,0,1.0);
        glRotated(-90,0,1,0);
        Man();
    glPopMatrix();

    //2
    glPushMatrix();
        glTranslated(-1.5,0,1.0);
        glRotated(-90,0,1,0);
        Man();
    glPopMatrix();

    //3
    glPushMatrix();
        glTranslated(-4.0,0,1.0);
        glRotated(-90,0,1,0);
        Man();
    glPopMatrix();

    /// moshjid
	glPushMatrix();

        glPushMatrix();// start banch1
            glTranslated(-2,0,0.29);
            ///banch start 1
            /// base
            glPushMatrix();
                //glColor3d(1,1.0,1.0);
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(0,1.55,0);
                glScaled(3,0.25,0.5);
                cube();
           glPopMatrix();
           /// back support
           glPushMatrix();
                glTranslated(0,1.8,0);
                glScaled(3,0.25,0.1);
                cube();
           glPopMatrix();

       glPopMatrix();// finished banch 1

       glPushMatrix();// start banch 2
            glTranslated(-7,0,0.29);
            ///banch start 1
            /// base
            glPushMatrix();
                //glColor3d(1,1.0,1.0);
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(0,1.55,0);
                glScaled(3,0.25,0.5);
                cube();
           glPopMatrix();
           /// back support
           glPushMatrix();
                glTranslated(0,1.8,0);
                glScaled(3,0.25,0.1);
                //glutSolidCube(1);
                cube();
           glPopMatrix();

       glPopMatrix();// finished banch 2

        glTranslated(-7,0,0);
        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D,5);
        glPushMatrix();
            glTranslated(0,1.55,0);
            glScaled(10,2,0.1);
            Cube();
       glPopMatrix();
       glDisable(GL_TEXTURE_2D);

    glPopMatrix();


}

void house_building3()
{
    glPushMatrix();
        glPushMatrix();
            glTranslated(0,1.55,0);
            /// base building
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                cube();
            glPopMatrix();
        /// window front
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(-0.01,0.3,0.3);
                glScaled(0.01,0.4,0.4);
                cube();
            glPopMatrix();

            /// window back
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(1.01,0.3,0.3);
                glScaled(0.01,0.4,0.4);
                cube();
            glPopMatrix();

            /// door front
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(0.3,0.15,1.01);
                glScaled(0.4,0.7,0.01);
                cube();
            glPopMatrix();

            /// top floor of building;
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(-0.1,1,-0.1);
                glScaled(1.2,0.02,1.2);
                cube();
            glPopMatrix();

        glPopMatrix();

        ///2nd floor
        glPushMatrix();
            glTranslated(0,2.57,0);
            /// base building
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                cube();
            glPopMatrix();
        /// window front
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(-0.01,0.3,0.3);
                glScaled(0.01,0.4,0.4);
                cube();
            glPopMatrix();

            /// window back
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(1.01,0.3,0.3);
                glScaled(0.01,0.4,0.4);
                cube();
            glPopMatrix();

            /// door front
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(0.3,0.15,1.01);
                glScaled(0.4,0.7,0.01);
                cube();
            glPopMatrix();

            /// top floor of building;
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(-0.1,1,-0.1);
                glScaled(1.2,0.02,1.2);
                cube();
            glPopMatrix();

        glPopMatrix();

        ///3rd floor
        glPushMatrix();
            glTranslated(0,3.59,0);
            /// base building
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                cube();
            glPopMatrix();
        /// window front
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(-0.01,0.3,0.3);
                glScaled(0.01,0.4,0.4);
                cube();
            glPopMatrix();

            /// window back
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(1.01,0.3,0.3);
                glScaled(0.01,0.4,0.4);
                cube();
            glPopMatrix();

            /// door front
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(0.3,0.15,1.01);
                glScaled(0.4,0.7,0.01);
                cube();
            glPopMatrix();

            /// top floor of building;
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(-0.1,1,-0.1);
                glScaled(1.2,0.02,1.2);
                cube();
            glPopMatrix();

        glPopMatrix();

    glPopMatrix();




}
void house_building2()
{
    glPushMatrix();
        glPushMatrix();
            glTranslated(0,1.55,0);
            /// base building
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient14);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse14);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular14);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess14);
                cube();
            glPopMatrix();
        /// window front
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(-0.01,0.3,0.3);
                glScaled(0.01,0.4,0.4);
                cube();
            glPopMatrix();

            /// window back
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(1.01,0.3,0.3);
                glScaled(0.01,0.4,0.4);
                cube();
            glPopMatrix();

            /// door front
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(0.3,0.15,1.01);
                glScaled(0.4,0.7,0.01);
                cube();
            glPopMatrix();

            /// top floor of building;
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(-0.1,1,-0.1);
                glScaled(1.2,0.02,1.2);
                cube();
            glPopMatrix();

        glPopMatrix();

        ///2nd floor
        glPushMatrix();
            glTranslated(0,2.57,0);
            /// base building
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient14);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse14);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular14);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess14);
                cube();
            glPopMatrix();
        /// window front
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(-0.01,0.3,0.3);
                glScaled(0.01,0.4,0.4);
                cube();
            glPopMatrix();

            /// window back
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(1.01,0.3,0.3);
                glScaled(0.01,0.4,0.4);
                cube();
            glPopMatrix();

            /// door front
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(0.3,0.15,1.01);
                glScaled(0.4,0.7,0.01);
                cube();
            glPopMatrix();

            /// top floor of building;
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(-0.1,1,-0.1);
                glScaled(1.2,0.02,1.2);
                cube();
            glPopMatrix();

        glPopMatrix();
    glPopMatrix();
}

void house_building1()
{
    glPushMatrix();
        glPushMatrix();
            glTranslated(0,1.55,0);
            /// base building
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient13);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse13);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular13);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess13);
                cube();
            glPopMatrix();
        /// window front
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(-0.01,0.3,0.3);
                glScaled(0.01,0.4,0.4);
                cube();
            glPopMatrix();

            /// window back
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(1.01,0.3,0.3);
                glScaled(0.01,0.4,0.4);
                cube();
            glPopMatrix();

            /// door front
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(0.3,0.15,1.01);
                glScaled(0.4,0.7,0.01);
                cube();
            glPopMatrix();

            /// top floor of building;
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(-0.1,1,-0.1);
                glScaled(1.2,0.02,1.2);
                cube();
            glPopMatrix();

        glPopMatrix();

        ///2nd floor
        glPushMatrix();
            glTranslated(0,2.57,0);
            /// base building
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient13);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse13);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular13);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess13);
                cube();
            glPopMatrix();
        /// window front
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(-0.01,0.3,0.3);
                glScaled(0.01,0.4,0.4);
                cube();
            glPopMatrix();

            /// window back
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(1.01,0.3,0.3);
                glScaled(0.01,0.4,0.4);
                cube();
            glPopMatrix();

            /// door front
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(0.3,0.15,1.01);
                glScaled(0.4,0.7,0.01);
                cube();
            glPopMatrix();

            /// top floor of building;
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(-0.1,1,-0.1);
                glScaled(1.2,0.02,1.2);
                cube();
            glPopMatrix();

        glPopMatrix();

        ///3rd floor
        glPushMatrix();
            glTranslated(0,3.59,0);
            /// base building
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient13);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse13);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular13);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess13);
                cube();
            glPopMatrix();
        /// window front
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(-0.01,0.3,0.3);
                glScaled(0.01,0.4,0.4);
                cube();
            glPopMatrix();

            /// window back
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(1.01,0.3,0.3);
                glScaled(0.01,0.4,0.4);
                cube();
            glPopMatrix();

            /// door front
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(0.3,0.15,1.01);
                glScaled(0.4,0.7,0.01);
                cube();
            glPopMatrix();

            /// top floor of building;
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(-0.1,1,-0.1);
                glScaled(1.2,0.02,1.2);
                cube();
            glPopMatrix();

        glPopMatrix();

    glPopMatrix();




}

void house_building()
{
    /// people
    //1
    glPushMatrix();
        glTranslated(0,0,1);
        //glScaled(2,1,2);
        Man();
    glPopMatrix();
    //2
    glPushMatrix();
        glTranslated(-1,0,1);
        //glScaled(2,1,2);
        Man();
    glPopMatrix();

    //3
    glPushMatrix();
        glTranslated(-2,0,1);
        //glRotated(180,0,1,0);
        //glScaled(2,1,2);
        Man();
    glPopMatrix();
    ///

    glPushMatrix();
        glTranslated(0,0,-0.25);
        house_building1();
    glPopMatrix();

    glPushMatrix();
        glTranslated(-1.5,0,-0.25);
        house_building2();
    glPopMatrix();

    glPushMatrix();
        glTranslated(-3,0,-0.25);
        house_building3();
    glPopMatrix();

    glPushMatrix();
        glTranslated(-4.5,0,-0.25);
        house_building1();
    glPopMatrix();


    glPushMatrix();
        glTranslated(0,0,-1.75);
        house_building3();
    glPopMatrix();

    glPushMatrix();
        glTranslated(-1.5,0,-1.75);
        house_building1();
    glPopMatrix();

    glPushMatrix();
        glTranslated(-3,0,-1.75);
        house_building2();
    glPopMatrix();

    glPushMatrix();
        glTranslated(-4.5,0,-1.75);
        house_building3();
    glPopMatrix();

}

void football(){
    /// Ground
    glPushMatrix();
        glPushMatrix();
            //glColor3d(1,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
            glTranslated(-(EN_SIZE*0.33),0,-EN_SIZE);
            glScaled(EN_SIZE*0.66,0.3,EN_SIZE*2);
            road();
            //cube();
        glPopMatrix();

        glPushMatrix();
        glTranslated(0,2.1,0);
            glRotated(180,0,1,0);
            car7();
        glPopMatrix();

        glPushMatrix();
            glTranslated(-2,2.1,-5);
            glRotated(90,0,1,0);
            car5();
        glPopMatrix();

        glPushMatrix();
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient1);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse1);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular1);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess1);
            glTranslated(0,0,0);
            glScaled(EN_SIZE*2,0.3,EN_SIZE*2);
            glutSolidCube(1);
        glPopMatrix();
    glPopMatrix();


    glPushMatrix();
        glTranslated(-8,-2.7,-5);
        glRotated(90,0,1,0);
        glTranslated(-4,0,0);
        glScaled(3,2,2);
        football_ground();
    glPopMatrix();

    glPushMatrix();
        glTranslated(8,-2.7,-5);
        glRotated(-90,0,1,0);
        glScaled(2,2,2);
        glTranslated(6,0,0);
        house_building();
    glPopMatrix();
}

void Duplex()
{
    glPushMatrix();
            glTranslated(-1.5,0,-1.75);
            glPushMatrix();
                glTranslated(0,1.55,0);
                /// base building
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient14);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse14);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular14);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess14);
                    glScaled(2,1,2);
                    cube();
                glPopMatrix();
            /// window front 1
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                    glTranslated(-0.01,0.3,0.2);
                    glScaled(0.01,0.4,0.4);
                    cube();
                glPopMatrix();

                ///door between 2 window
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                    glTranslated(-0.01,0.05,0.8);
                    glScaled(0.01,0.7,0.4);
                    cube();
                glPopMatrix();

                /// window front 2
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                    glTranslated(-0.01,0.3,1.4);
                    glScaled(0.01,0.4,0.4);
                    cube();
                glPopMatrix();

                /// window back 1
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                    glTranslated(2.01,0.3,0.2);
                    glScaled(0.01,0.4,0.4);
                    cube();
                glPopMatrix();

                /// window back 2
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                    glTranslated(2.01,0.3,1.4);
                    glScaled(0.01,0.4,0.4);
                    cube();
                glPopMatrix();

                ///  rode side window front
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                    glTranslated(0.2,0.3,2.01);
                    glScaled(0.4,0.4,0.01);
                    cube();
                glPopMatrix();

                ///  rode side door front 1
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                    glTranslated(0.8,0.05,2.01);
                    glScaled(0.4,0.7,0.01);
                    cube();
                glPopMatrix();

                ///  rode side window front 2
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                    glTranslated(1.4,0.3,2.01);
                    glScaled(0.4,0.4,0.01);
                    cube();
                glPopMatrix();




            /// ground floor of building;
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient10);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse10);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular10);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess10);
                glTranslated(-0.5,0,0.0);
                glScaled(2.5,0.02,2.5);
                cube();
            glPopMatrix();

            glPushMatrix();///ralling ground floor

                /// piler of building 1;
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(-0.5,0,0.0);
                    glScaled(0.1,1,0.1);
                    cube();
                glPopMatrix();

                /// piler of building 2;
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(-0.5,0,2.4);
                    glScaled(0.1,1,0.1);
                    cube();
                glPopMatrix();

                /// piler of building 3;
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(1.9,0,2.4);
                    glScaled(0.1,1,0.1);
                    cube();
                glPopMatrix();

                ///grill front left part small
                    glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(-0.5,0.35,0);
                    glScaled(0.05,0.05,0.8);
                    cube();
                glPopMatrix();



                ///grill front left part large
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(-0.5,0.35,1.2);
                    glScaled(0.05,0.05,1.3);
                    cube();
                glPopMatrix();

                ///rod 1
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(-0.5,0.0,0.25);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();

                ///rod 2
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(-0.5,0.0,0.5);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();

                ///rod 3
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(-0.5,0.0,0.75);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();
                /// large part rod
                ///rod 4
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(-0.5,0.0,1.2);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();
                ///rod 5
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(-0.5,0.0,1.45);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();
                ///rod 6
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(-0.5,0.0,1.7);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();

                ///rod 7
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(-0.5,0.0,1.95);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();

                ///rod 8
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(-0.5,0.0,2.2);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();

                /// rode side grill large
                    glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(-0.5,0.35,2.45);
                    glScaled(1.3,0.05,0.05);
                    cube();
                glPopMatrix();

                /// rode side grill small
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(1.2,0.35,2.45);
                    glScaled(0.8,0.05,0.05);
                    cube();
                glPopMatrix();

                ///rod 9
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(-0.25,0.0,2.45);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();

                ///rod 10
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(0.0,0.0,2.45);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();

                ///rod 11
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(0.25,0.0,2.45);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();

                ///rod 12
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(0.5,0.0,2.45);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();

                ///rod 13
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(0.75,0.0,2.45);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();
                ///small
                ///rod 14
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(1.2,0.0,2.45);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();

                ///rod 15
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(1.45,0.0,2.45);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();

                ///rod 16
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(1.75,0.0,2.45);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();

            glPopMatrix();/// finished ground floor railing

            glPushMatrix();///start 2nd floor ralling
                glTranslated(0,1.01,0);
                /// piler of building 1;
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(-0.5,0,0.0);
                    glScaled(0.1,1,0.1);
                    cube();
                glPopMatrix();

                /// piler of building 2;
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(-0.5,0,2.4);
                    glScaled(0.1,1,0.1);
                    cube();
                glPopMatrix();

                /// piler of building 3;
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(1.90,0,2.4);
                    glScaled(0.1,1,0.1);
                    cube();
                glPopMatrix();

                ///grill front left part small
                    glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(-0.5,0.35,0);
                    glScaled(0.05,0.05,0.8);
                    cube();
                glPopMatrix();
                ///grill front left part large
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(-0.5,0.35,1.2);
                    glScaled(0.05,0.05,1.3);
                    cube();
                glPopMatrix();

                ///rod 1
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(-0.5,0.0,0.25);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();

                ///rod 2
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(-0.5,0.0,0.5);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();

                ///rod 3
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(-0.5,0.0,0.75);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();
                /// large part rod
                ///rod 4
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(-0.5,0.0,1.2);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();
                ///rod 5
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(-0.5,0.0,1.45);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();
                ///rod 6
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(-0.5,0.0,1.7);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();

                ///rod 7
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(-0.5,0.0,1.95);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();

                ///rod 8
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(-0.5,0.0,2.2);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();

                /// rode side grill large
                    glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(-0.5,0.35,2.45);
                    glScaled(1.3,0.05,0.05);
                    cube();
                glPopMatrix();

                /// rode side grill small
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(1.2,0.35,2.45);
                    glScaled(0.8,0.05,0.05);
                    cube();
                glPopMatrix();

                ///rod 9
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(-0.25,0.0,2.45);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();

                ///rod 10
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(0.0,0.0,2.45);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();

                ///rod 11
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(0.25,0.0,2.45);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();

                ///rod 12
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(0.5,0.0,2.45);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();

                ///rod 13
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(0.75,0.0,2.45);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();
                ///small
                ///rod 14
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(1.2,0.0,2.45);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();

                ///rod 15
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(1.45,0.0,2.45);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();

                ///rod 16
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                    glTranslated(1.75,0.0,2.45);
                    glScaled(0.05,0.35,0.05);
                    cube();
                glPopMatrix();

            glPopMatrix();///finished top floor railing

   glPopMatrix();

///1st floor
            glPushMatrix();
                glTranslated(0,2.56,0);
                /// base building
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient14);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse14);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular14);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess14);
                    glScaled(2,1,2);
                    cube();
                glPopMatrix();
            /// window front 1
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                    glTranslated(-0.01,0.3,0.2);
                    glScaled(0.01,0.4,0.4);
                    cube();
                glPopMatrix();

                ///door between 2 window
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                    glTranslated(-0.01,0.05,0.8);
                    glScaled(0.01,0.7,0.4);
                    cube();
                glPopMatrix();

                /// window front 2
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                    glTranslated(-0.01,0.3,1.4);
                    glScaled(0.01,0.4,0.4);
                    cube();
                glPopMatrix();

                /// window back 1
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                    glTranslated(2.01,0.3,0.2);
                    glScaled(0.01,0.4,0.4);
                    cube();
                glPopMatrix();

                /// window back 2
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                    glTranslated(2.01,0.3,1.4);
                    glScaled(0.01,0.4,0.4);
                    cube();
                glPopMatrix();

                ///  rode side window front
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                    glTranslated(0.2,0.3,2.01);
                    glScaled(0.4,0.4,0.01);
                    cube();
                glPopMatrix();

                ///  rode side door front 1
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                    glTranslated(0.8,0.05,2.01);
                    glScaled(0.4,0.7,0.01);
                    cube();
                glPopMatrix();

                ///  rode side window front 2
                glPushMatrix();
                    //house_building1();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                    glTranslated(1.4,0.3,2.01);
                    glScaled(0.4,0.4,0.01);
                    cube();
                glPopMatrix();

            /// 1st floor of building;
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient10);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse10);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular10);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess10);
                glTranslated(-0.5,0,0.0);
                glScaled(2.5,0.02,2.5);
                cube();
            glPopMatrix();

            /// top floor of building;
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient10);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse10);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular10);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess10);
                glTranslated(-0.5,1.01,0.0);
                glScaled(2.5,0.02,2.5);
                cube();
            glPopMatrix();
            glPopMatrix();

        glPopMatrix();
}

void building(){

    //people
    //1
    glPushMatrix();
        glTranslated(0.0,0,1.0);
        Man();
    glPopMatrix();
    //2
    glPushMatrix();
        glTranslated(-1.5,0,1.0);
        Man();
    glPopMatrix();
    //3
    glPushMatrix();
        glTranslated(-4.0,0,1.0);
        Man();
    glPopMatrix();

    glPushMatrix();
        Duplex();
    glPopMatrix();
    /// moshjid
	glPushMatrix();
        glTranslated(-4.5,0,-2.0);
        glEnable(GL_TEXTURE_2D);
        glBindTexture(GL_TEXTURE_2D,3);
        glPushMatrix();
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
            glTranslated(-1,1.55,0);
            glScaled(1.5,1,1.5);
            Cube();
       glPopMatrix();
       glDisable(GL_TEXTURE_2D);





    glPopMatrix();


}

void banch(){
    /// Ground
    glPushMatrix();
        glPushMatrix();
            //glColor3d(1,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
            glTranslated(-(EN_SIZE*0.33),0,-EN_SIZE);
            glScaled(EN_SIZE*0.66,0.3,EN_SIZE*2);
            road();
            //cube();
        glPopMatrix();

        glPushMatrix();
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient12);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse12);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular12);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess12);
            glTranslated(0,1.3,0);
            //glRotated(180,0,1,0);
            glutSolidSphere(0.7,30,30);
        glPopMatrix();

        glPushMatrix();
            glTranslated(0,0,0);
            glScaled(EN_SIZE*2,0.3,EN_SIZE*2);
            //glutSolidCube(1);
        glPopMatrix();
    glPopMatrix();


    glPushMatrix();
        glTranslated(-8,-2.7,-5);
        glRotated(90,0,1,0);
        glScaled(2,2,2);
         building();
    glPopMatrix();

    glPushMatrix();
        glTranslated(8,-2.7,-5);
        glRotated(-90,0,1,0);
        glTranslated(10,0,0);
        glScaled(2,2,2);
        banch1();
    glPopMatrix();
}

void tree2(){
	/// ground
	glPushMatrix();
    ///base body;
	glPushMatrix();
	    //glColor3d(0,0.40,0.0);
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient10);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse10);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular10);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess10);
        //glTranslated(-2,0,3.5);
        glTranslated(0,1.55,0);
        glScaled(4,0.05,3);
        glRotated(-90,1,0,0);;
        GLUquadricObj *quadratic;
        quadratic = gluNewQuadric();
        gluCylinder(quadratic, 0.0125, 0.0125, 20, 20, 8);
    glPopMatrix();

    ///top leaf stand cylinder;
    glPushMatrix();// start leaf with stand;
    glTranslated(0,2.55,0);
    //stand
        glPushMatrix();
            glScaled(4,0.05,3);
            glRotated(-90,1,0,0);
            quadratic = gluNewQuadric();
            gluCylinder(quadratic, 0.005, 0.005,3 , 20, 8);
        glPopMatrix();
        ///top leaf sphare
            glPushMatrix();
            glTranslated(0,.30,0);
            glScaled(0.18,0.4,0.18);
            glutSolidSphere(1,30,30);
        glPopMatrix();

    glPopMatrix();// finished leaf with stand;

    glPushMatrix();// down layer of leaf start

        ///road side leaf stand cylinder;
        glPushMatrix();// start leaf with stand;
             glTranslated(0,2.0,0);
             glRotated(45,1,0,0);
            /// stand
            glPushMatrix();
                glScaled(4,0.05,3);
                glRotated(-90,1,0,0);
                quadratic = gluNewQuadric();
                gluCylinder(quadratic, 0.005, 0.005,3 , 20, 8);
            glPopMatrix();
            ///top leaf sphare
                glPushMatrix();
                glTranslated(0,.30,0);
                glScaled(0.18,0.4,0.18);
                glutSolidSphere(1,30,30);
            glPopMatrix();

        glPopMatrix();// finished leaf with stand;

        ///road back side leaf stand cylinder;
        glPushMatrix();// start leaf with stand;
             glTranslated(0,2.0,0);
             glRotated(-45,1,0,0);
            /// stand
            glPushMatrix();
                glScaled(4,0.05,3);
                glRotated(-90,1,0,0);
                quadratic = gluNewQuadric();
                gluCylinder(quadratic, 0.005, 0.005,3 , 20, 8);
            glPopMatrix();
            ///top leaf sphare
                glPushMatrix();
                glTranslated(0,.30,0);
                glScaled(0.18,0.4,0.18);
                glutSolidSphere(1,30,30);
            glPopMatrix();

        glPopMatrix();// finished leaf with stand;

        ///road 1st side leaf stand cylinder;
        glPushMatrix();// start leaf with stand;
             glTranslated(0,2.0,0);
             glRotated(45,0,0,1);
            /// stand
            glPushMatrix();
                glScaled(4,0.05,3);
                glRotated(-90,1,0,0);
                quadratic = gluNewQuadric();
                gluCylinder(quadratic, 0.005, 0.005,3 , 20, 8);
            glPopMatrix();
            ///top leaf sphare
                glPushMatrix();
                glTranslated(0,.30,0);
                glScaled(0.18,0.4,0.18);
                glutSolidSphere(1,30,30);
            glPopMatrix();

        glPopMatrix();// finished leaf with stand;

        ///road 2st side leaf stand cylinder;
        glPushMatrix();// start leaf with stand;
             glTranslated(0,2.0,0);
             glRotated(-45,0,0,1);
            /// stand
            glPushMatrix();
                glScaled(4,0.05,3);
                glRotated(-90,1,0,0);

                quadratic = gluNewQuadric();
                gluCylinder(quadratic, 0.005, 0.005,3 , 20, 8);
            glPopMatrix();
            ///top leaf sphare
                glPushMatrix();

                glTranslated(0,.30,0);
                glScaled(0.18,0.4,0.18);
                glutSolidSphere(1,30,30);
            glPopMatrix();

        glPopMatrix();// finished leaf with stand;

    glPopMatrix(); //down layer finished

    ///

    glPushMatrix();// up layer of leaf start
    glTranslated(0,0.35,0);

        ///road side leaf stand cylinder;
        glPushMatrix();// start leaf with stand;
             glTranslated(0,2.0,0);
             glRotated(45,1,0,0);
            /// stand
            glPushMatrix();
                glScaled(4,0.05,3);
                glRotated(-90,1,0,0);
                quadratic = gluNewQuadric();
                gluCylinder(quadratic, 0.005, 0.005,3 , 20, 8);
            glPopMatrix();
            ///top leaf sphare
                glPushMatrix();
                glTranslated(0,.30,0);
                glScaled(0.18,0.18,0.18);
                glutSolidSphere(1,30,30);
            glPopMatrix();

        glPopMatrix();// finished leaf with stand;

        ///road back side leaf stand cylinder;
        glPushMatrix();// start leaf with stand;
             glTranslated(0,2.0,0);
             glRotated(-45,1,0,0);
            /// stand
            glPushMatrix();
                glScaled(4,0.05,3);
                glRotated(-90,1,0,0);
                quadratic = gluNewQuadric();
                gluCylinder(quadratic, 0.005, 0.005,3 , 20, 8);
            glPopMatrix();
            ///top leaf sphare
            glPushMatrix();
                glTranslated(0,.30,0);
                glScaled(0.18,0.18,0.18);
                glutSolidSphere(1,30,30);
            glPopMatrix();

        glPopMatrix();// finished leaf with stand;

        ///road 1st side leaf stand cylinder;
        glPushMatrix();// start leaf with stand;
             glTranslated(0,2.0,0);
             glRotated(45,0,0,1);
            /// stand
            glPushMatrix();
                glScaled(4,0.05,3);
                glRotated(-90,1,0,0);
                quadratic = gluNewQuadric();
                gluCylinder(quadratic, 0.005, 0.005,3 , 20, 8);
            glPopMatrix();
            ///top leaf sphare
                glPushMatrix();
                glTranslated(0,.30,0);
                glScaled(0.18,0.18,0.18);
                glutSolidSphere(1,30,30);
            glPopMatrix();

        glPopMatrix();// finished leaf with stand;

        ///road 2st side leaf stand cylinder;
        glPushMatrix();// start leaf with stand;
             glTranslated(0,2.0,0);
             glRotated(-45,0,0,1);
            /// stand
            glPushMatrix();
                glScaled(4,0.05,3);
                glRotated(-90,1,0,0);
                quadratic = gluNewQuadric();
                gluCylinder(quadratic, 0.005, 0.005,3 , 20, 8);
            glPopMatrix();
            ///top leaf sphare
                glPushMatrix();
                glTranslated(0,.30,0);
                glScaled(0.18,0.18,0.18);
                glutSolidSphere(1,30,30);
            glPopMatrix();

        glPopMatrix();// finished leaf with stand;




    glPopMatrix(); //up layer finished


    glPopMatrix();

}

void tree1(){


    glPushMatrix();
    ///
        glPushMatrix();//3rd
            glTranslated(1,0,0.5);
            tree2();
        glPopMatrix();
     ///
        glPushMatrix();//4th
            glTranslated(0,0,0.5);
            tree2();
        glPopMatrix();
     ///
        glPushMatrix();//5th
            glTranslated(-1,0,0.5);
            tree2();
        glPopMatrix();

     ///
        glPushMatrix();//6th
            glTranslated(-2,0,0.5);
            tree2();
        glPopMatrix();

         ///
        glPushMatrix();//7th
            glTranslated(-3,0,0.5);
            tree2();
        glPopMatrix();
      ///
        glPushMatrix();//8th
            glTranslated(-4,0,0.5);
            tree2();
        glPopMatrix();

         ///
        glPushMatrix();//9th
            glTranslated(-5,0,0.5);
            //tree2();
        glPopMatrix();
    glPopMatrix();


    /// 2nd layer

    glPushMatrix();
        glTranslated(0,0,-1.5);
    ///
        glPushMatrix();//3rd
            glTranslated(1,0,0.5);
            tree2();
        glPopMatrix();
     ///
        glPushMatrix();//4th
            glTranslated(0,0,0.5);
            tree2();
        glPopMatrix();
     ///
        glPushMatrix();//5th
            glTranslated(-1,0,0.5);
            tree2();
        glPopMatrix();

     ///
        glPushMatrix();//6th
            glTranslated(-2,0,0.5);
            tree2();
        glPopMatrix();

         ///
        glPushMatrix();//7th
            glTranslated(-3,0,0.5);
            tree2();
        glPopMatrix();
      ///
        glPushMatrix();//8th
            glTranslated(-4,0,0.5);
            tree2();
        glPopMatrix();

         ///
        glPushMatrix();//9th
            glTranslated(-5,0,0.5);
            tree2();
        glPopMatrix();
    glPopMatrix();

    /// 3rd layer

    glPushMatrix();
        glTranslated(0,0,-3);
    ///
        glPushMatrix();//3rd
            glTranslated(1,0,0.5);
            tree2();
        glPopMatrix();
     ///
        glPushMatrix();//4th
            glTranslated(0,0,0.5);
            tree2();
        glPopMatrix();
     ///
        glPushMatrix();//5th
            glTranslated(-1,0,0.5);
            tree2();
        glPopMatrix();

     ///
        glPushMatrix();//6th
            glTranslated(-2,0,0.5);
            tree2();
        glPopMatrix();

         ///
        glPushMatrix();//7th
            glTranslated(-3,0,0.5);
            tree2();
        glPopMatrix();
      ///
        glPushMatrix();//8th
            glTranslated(-4,0,0.5);
            tree2();
        glPopMatrix();

         ///
        glPushMatrix();//9th
            glTranslated(-5,0,0.5);
            tree2();
        glPopMatrix();
    glPopMatrix();




}

void tree(){
    /// Ground

    //glColor3d(0,0.5,0.1);

    glPushMatrix();
        glPushMatrix();
            glColor3d(1,0,0);
            glTranslated(-(EN_SIZE*0.33),0,-EN_SIZE);
            glScaled(EN_SIZE*0.66,0.3,EN_SIZE*2);
            road();
            //cube();
        glPopMatrix();

        glPushMatrix();
            //glColor3d(0,0.5,0.1);
            glColor3d(1,0.8,1);
            glTranslated(0,0,0);
            glScaled(EN_SIZE*2,0.3,EN_SIZE*2);
            glutSolidCube(1);
        glPopMatrix();
        //Cube();

        glPushMatrix();
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient12);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse12);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular12);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess12);
            glTranslated(-4,1.55,0);
            glScaled(0.5,0.5,0.5);
            glutSolidTorus(.5,2,20,20);
        glPopMatrix();

        glPushMatrix();
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient5);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse5);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular5);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess5);
            glTranslated(3,1.55,0);
            glScaled(0.5,0.5,0.5);
            glutSolidTorus(.5,2,20,20);
        glPopMatrix();

        glPushMatrix();
            glTranslated(3,2.1,-3);
            glRotated(90,0,1,0);
            car3();
        glPopMatrix();

    glPopMatrix();


    glPushMatrix();
        glTranslated(-8,-2.7,-5);
        glRotated(90,0,1,0);
        glTranslated(-4,0,0);
        glScaled(3,2,2);
         tree1();
    glPopMatrix();

    glPushMatrix();
        glTranslated(8,-2.7,-5);
        glRotated(-90,0,1,0);
        glScaled(2,2,2);
        glTranslated(5,0,0);
        tree1();
    glPopMatrix();
}

void truck(){
    /// Ground

    //glColor3d(0,0.5,0.1);

    glPushMatrix();
        glPushMatrix();
            glColor3d(1,0,0);
            glTranslated(-(EN_SIZE*0.33),0,-EN_SIZE);
            glScaled(EN_SIZE*0.66,0.3,EN_SIZE*2);
            road();
            //cube();
        glPopMatrix();

        glPushMatrix();
            //glColor3d(0,0.5,0.1);
            glColor3d(1,0.8,1);
            glTranslated(0,0,0);
            glScaled(EN_SIZE*2,0.3,EN_SIZE*2);
            glutSolidCube(1);
        glPopMatrix();
        //Cube();

        glPushMatrix();
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient5);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse5);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular5);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess5);
            glTranslated(-4,1.55,0);
            glScaled(0.5,0.5,0.5);
            glutSolidTorus(.5,2,20,20);
        glPopMatrix();

        glPushMatrix();
            glTranslated(-3,2.1,-3);
            glRotated(90,0,1,0);
            car4();
        glPopMatrix();

        glPushMatrix();
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient12);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse12);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular12);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess12);
            glTranslated(3,1.55,0);
            glScaled(0.5,0.5,0.5);
            glutSolidTorus(.5,2,20,20);
        glPopMatrix();

    glPopMatrix();


    glPushMatrix();
        glTranslated(-8,-2.7,-5);
        glRotated(90,0,1,0);
        glTranslated(-4,0,0);
        glScaled(3,2,2);
         car_station();
    glPopMatrix();

    glPushMatrix();
        glTranslated(8,-2.7,-5);
        glRotated(-90,0,1,0);
        glScaled(2,2,2);
        glTranslated(8,0,0);
        car_station();
    glPopMatrix();
}

void home3(){
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,7);
    //glColor3d(1,1.0,1.0);
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
	glPushMatrix();
        glTranslated(0,1.55,0);
        glScaled(1,2,0.8);
        //glutSolidCube(1);
        Cube();
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);
}

void home2(){
    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D,6);
    //glColor3d(1,1.0,1.0);
        glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
        glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
        glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
        glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
	glPushMatrix();
        glTranslated(0,1.55,0);
        glScaled(1,2,0.8);
        //glutSolidCube(1);
        Cube();
    glPopMatrix();
    glDisable(GL_TEXTURE_2D);
}

void home1(){
    ///1
    glPushMatrix();
        home2();
    glPopMatrix();

    ///1 back
    glPushMatrix();
        glTranslated(0,0,-1.3);
        home2();
    glPopMatrix();

    ///2
    glPushMatrix();
        glTranslated(-2,0,0);
        home3();
    glPopMatrix();

    ///2 back
    glPushMatrix();
        glTranslated(-2,0,-1.3);
        home3();
    glPopMatrix();

    ///3
    glPushMatrix();
        glTranslated(-4,0,0);
        home2();
    glPopMatrix();
    ///3 back
    glPushMatrix();
        glTranslated(-4,0,-1.3);
        home2();
    glPopMatrix();

    ///4
    glPushMatrix();
        glTranslated(-6,0,0);
        home3();
    glPopMatrix();

    ///4 back
    glPushMatrix();
        glTranslated(-6,0,-1.3);
        home3();
    glPopMatrix();

    ///5
    glPushMatrix();
        glTranslated(-8,0,0);
        home3();
    glPopMatrix();
    ///5 back
    glPushMatrix();
        glTranslated(-8,0,-1.3);
        home3();
    glPopMatrix();

}

void station_unit()
{
    glPushMatrix();
            glTranslated(0,1.55,0);
            /// base building
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient16);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse16);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular16);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess16);
                cube();
            glPopMatrix();
        /// window front
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(-0.01,0.3,0.3);
                glScaled(0.01,0.4,0.4);
                cube();
            glPopMatrix();

            /// window back
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(1.01,0.3,0.3);
                glScaled(0.01,0.4,0.4);
                cube();
            glPopMatrix();

            /// door front
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(0.3,0.15,1.01);
                glScaled(0.4,0.7,0.01);
                cube();
            glPopMatrix();

            /// top floor of building;
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(-0.1,1,-0.1);
                glScaled(1.2,0.02,1.2);
                cube();
            glPopMatrix();

            /// front left border
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(0.0,0,-0.01);
                glScaled(0.01,1,0.02);
                cube();
            glPopMatrix();

            /// front right border
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(0.0,0,1.01);
                glScaled(0.01,1,0.02);
                cube();
            glPopMatrix();

            /// back left border
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(1.01,0,-0.01);
                glScaled(0.01,1,0.02);
                cube();
            glPopMatrix();

            /// back right border
            glPushMatrix();
                //house_building1();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
                glTranslated(1.01,0,1.01);
                glScaled(0.01,1,0.02);
                cube();
            glPopMatrix();


        glPopMatrix();
}

station_rest()
{
        if(wired)
    {
        glPolygonMode( GL_FRONT, GL_LINE ) ;
        glPolygonMode( GL_BACK, GL_LINE ) ;

    }
    else
    {
        glPolygonMode( GL_FRONT,GL_FILL ) ;
        glPolygonMode( GL_BACK, GL_FILL ) ;
    }

    glPushMatrix();
        glTranslated(0,0,-0.6);
        glPushMatrix();
                glTranslated(0,1.55,0.0);
                ///banch start 1
                glPushMatrix();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient15);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse15);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular15);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess15);
                glRotated(90,1,0,0);
                glScaled(0.75,0.75,0.75);
                glutSolidTorus(0.3,1.5,20,20);

               glPopMatrix();
        glPopMatrix();

        glPushMatrix();
            glPushMatrix();
                glTranslated(0,1.55,0);
                glScaled(0.05,1.25,0.05);
                Cube();
            glPopMatrix();
            glPushMatrix();
                    glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient11);
                    glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse11);
                    glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular11);
                    glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess11);
                glTranslated(0,3.4,0);
                glScaled(0.25,0.50,0.25);
                glRotated(-90,0,0,1);
                station_top();
            glPopMatrix();
        glPopMatrix();

        glPushMatrix();
            glTranslated(-0.25,0,-0.25);
            ///table base
            glPushMatrix();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient11);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse11);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular11);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess11);
                glTranslated(-0.05,2.0,-0.05);
                glScaled(0.6,0.02,0.6);
                cube();
            glPopMatrix();
             ///front left leg
            glPushMatrix();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient11);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse11);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular11);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess11);
                glTranslated(0,1.55,0);
                glScaled(0.02,0.45,0.02);
                cube();
            glPopMatrix();

            ///front right leg
            glPushMatrix();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient11);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse11);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular11);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess11);
                glTranslated(0,1.55,0.5);
                glScaled(0.02,0.45,0.02);
                cube();
            glPopMatrix();

            ///back left leg
            glPushMatrix();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient11);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse11);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular11);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess11);
                glTranslated(0.5,1.55,0);
                glScaled(0.02,0.45,0.02);
                cube();
            glPopMatrix();

            ///right right leg
            glPushMatrix();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient11);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse11);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular11);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess11);
                glTranslated(0.5,1.55,0.5);
                glScaled(0.02,0.45,0.02);
                cube();
            glPopMatrix();


        glPopMatrix();
    glPopMatrix();


    if(wired)
    {
        glPolygonMode( GL_FRONT, GL_LINE ) ;
        glPolygonMode( GL_BACK, GL_LINE ) ;

    }
    else
    {
        glPolygonMode( GL_FRONT,GL_FILL ) ;
        glPolygonMode( GL_BACK, GL_FILL ) ;
    }

}

void station()
{
    ///1st floor
    glPushMatrix();
        ///widthwise
        ///1st
       glPushMatrix();
            glTranslated(0,0,-.25);
            //glRotated(-90,0,1,0);
            station_unit();
       glPopMatrix();

        ///2nd
       glPushMatrix();
            glTranslated(1,0,-1.25);
            glRotated(-90,0,1,0);
            station_unit();
       glPopMatrix();

       ///3rd
       glPushMatrix();
            glTranslated(1,0,-2.25);
            glRotated(-90,0,1,0);
            station_unit();
       glPopMatrix();

        ///4th
       glPushMatrix();
            glTranslated(1,0,-3.25);
            glRotated(-90,0,1,0);
            station_unit();
        glPopMatrix();

        ///length wise

        ///5th
       glPushMatrix();
            glTranslated(-1,0,-4.25);
            //glRotated(-90,0,1,0);
            station_unit();
       glPopMatrix();

        ///6th
       glPushMatrix();
            glTranslated(-2,0,-4.25);
            //glRotated(-90,0,1,0);
            station_unit();
       glPopMatrix();

        ///7th
       glPushMatrix();
            glTranslated(-3,0,-4.25);
            //glRotated(-90,0,1,0);
            station_unit();
       glPopMatrix();

        ///8th
       glPushMatrix();
            glTranslated(-4,0,-4.25);
            //glRotated(-90,0,1,0);
            station_unit();
       glPopMatrix();

        ///9th
       glPushMatrix();
            glTranslated(-5,0,-4.25);
            //glRotated(-90,0,1,0);
            station_unit();
       glPopMatrix();

       ///

       ///wide
       ///10th
          glPushMatrix();
            glTranslated(-6,0,-2.25);
            glRotated(90,0,1,0);
            station_unit();
       glPopMatrix();

          ///11th
          glPushMatrix();
            glTranslated(-6,0,-1.25);
            glRotated(90,0,1,0);
            station_unit();
       glPopMatrix();

          ///12th
          glPushMatrix();
            glTranslated(-6,0,-0.25);
            glRotated(90,0,1,0);
            station_unit();
       glPopMatrix();

          ///13th
          glPushMatrix();
            glTranslated(-6,0,0.75);
            glRotated(90,0,1,0);
            station_unit();
       glPopMatrix();

    glPopMatrix();

    ///2nd floor
    glPushMatrix();
        glTranslated(0,1.02,0);
        ///widthwise
        ///1st
       glPushMatrix();
            glTranslated(0,0,-.25);
            //glRotated(-90,0,1,0);
            station_unit();
       glPopMatrix();

        ///2nd
       glPushMatrix();
            glTranslated(1,0,-1.25);
            glRotated(-90,0,1,0);
            station_unit();
       glPopMatrix();

       ///3rd
       glPushMatrix();
            glTranslated(1,0,-2.25);
            glRotated(-90,0,1,0);
            station_unit();
       glPopMatrix();

        ///4th
       glPushMatrix();
            glTranslated(1,0,-3.25);
            glRotated(-90,0,1,0);
            station_unit();
        glPopMatrix();

        ///length wise

        ///5th
       glPushMatrix();
            glTranslated(-1,0,-4.25);
            //glRotated(-90,0,1,0);
            station_unit();
       glPopMatrix();

        ///6th
       glPushMatrix();
            glTranslated(-2,0,-4.25);
            //glRotated(-90,0,1,0);
            station_unit();
       glPopMatrix();

        ///7th
       glPushMatrix();
            glTranslated(-3,0,-4.25);
            //glRotated(-90,0,1,0);
            station_unit();
       glPopMatrix();

        ///8th
       glPushMatrix();
            glTranslated(-4,0,-4.25);
            //glRotated(-90,0,1,0);
            station_unit();
       glPopMatrix();

        ///9th
       glPushMatrix();
            glTranslated(-5,0,-4.25);
            //glRotated(-90,0,1,0);
            station_unit();
       glPopMatrix();

       ///

       ///wide
       ///10th
          glPushMatrix();
            glTranslated(-6,0,-2.25);
            glRotated(90,0,1,0);
            station_unit();
       glPopMatrix();

          ///11th
          glPushMatrix();
            glTranslated(-6,0,-1.25);
            glRotated(90,0,1,0);
            station_unit();
       glPopMatrix();

          ///12th
          glPushMatrix();
            glTranslated(-6,0,-0.25);
            glRotated(90,0,1,0);
            station_unit();
       glPopMatrix();

          ///13th
          glPushMatrix();
            glTranslated(-6,0,0.75);
            glRotated(90,0,1,0);
            station_unit();
       glPopMatrix();

    glPopMatrix();



    glPushMatrix();
        glTranslated(-2,0,0);
        station_rest();
   glPopMatrix();

    //people
    //people
    glPushMatrix();
        glTranslated(-2.75,-1.6,-0.5);
        //glRotated(180,0,1,0);
        glScaled(1.5,2,1.5);
        //glScaled(2,1,2);
        Man();
    glPopMatrix();
    glPushMatrix();
        glTranslated(-2,-1.6,-1.25);
        glRotated(-90,0,1,0);
        glScaled(1.5,2,1.5);
        //glScaled(2,1,2);
        Man();
    glPopMatrix();

        //people
    glPushMatrix();
        glTranslated(-1.25,-1.6,-0.5);
        glRotated(180,0,1,0);
        glScaled(1.5,2,1.5);
        //glScaled(2,1,2);
        Man();
    glPopMatrix();
}

void home(){
    /// Ground

    //glColor3d(0,0.5,0.1);

    glPushMatrix();
        glPushMatrix();
            //glColor3d(1,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
        //glTranslated(-2,0,3.5);
            glTranslated(-(EN_SIZE*0.33),0,-EN_SIZE);
            glScaled(EN_SIZE*0.66,0.3,EN_SIZE*2);
            road();
            //cube();
        glPopMatrix();

        glPushMatrix();
            glTranslated(0,2.1,0);
            //glRotated(180,0,1,0);
            car6();
        glPopMatrix();

        glPushMatrix();
            glTranslated(0,2.1,-5);
            glRotated(90,0,1,0);
            car4();
        glPopMatrix();

        glPushMatrix();
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient1);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse1);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular1);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess1);
            glTranslated(0,0,0);
            glScaled(EN_SIZE*2,0.3,EN_SIZE*2);
            glutSolidCube(1);
        glPopMatrix();
        //Cube();
    glPopMatrix();


    glPushMatrix();
        glTranslated(-8,-2.7,-5);
        glRotated(90,0,1,0);
        glTranslated(-4,0,0);
        glScaled(2,2,2);
         station();
    glPopMatrix();

    glPushMatrix();
        glTranslated(8,-2.7,-5);
        glRotated(-90,0,1,0);
        glScaled(2,2,2);
        glTranslated(8,0,0);
         station();
    glPopMatrix();
}

void banch2(){
    /// Ground
    glPushMatrix();
        glPushMatrix();
            //lColor3d(1,0,0);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient3);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse3);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular3);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess3);
            glTranslated(-(EN_SIZE*0.33),0,-EN_SIZE);
            glScaled(EN_SIZE*0.66,0.3,EN_SIZE*2);
            road();
            //cube();
        glPopMatrix();
        glPushMatrix();
            glTranslated(0,0,0);
            glScaled(EN_SIZE*2,0.3,EN_SIZE*2);
            glutSolidCube(1);
        glPopMatrix();
        //Cube();
        glPushMatrix();
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient10);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse10);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular10);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess10);
            glTranslated(0,1.0,0);
            glScaled(0.5,0.5,0.5);
            glutSolidTorus(0.5,1,20,20);
        glPopMatrix();


    glPopMatrix();


    glPushMatrix();
        glTranslated(-8,-2.7,-5);
        glRotated(90,0,1,0);
        glScaled(2,2,2);
         house_building();
    glPopMatrix();

    glPushMatrix();
        glTranslated(8,-2.7,-5);
        glRotated(-90,0,1,0);
        glTranslated(10,0,0);
        glScaled(2,2,2);
        home1();
    glPopMatrix();
}

void draw(){
    double t = glutGet(GLUT_ELAPSED_TIME) / 1000.0;
    double a = t*90.0;

    TIME = t;

    ///Plane
    if(rotX>11)rotX=11;
    if(rotX<-11)rotX=-11;
    if(rotZ>10)rotZ=10;
    if(rotZ<-15)rotZ=-15;

    glPushMatrix();
        //glTranslated(0,1,0);
        glTranslated(0,1,2);
        glRotated(90,0,1,0);
        glRotated(5,0,0,1);
        glRotated(rotX,1,0,0);
        glRotated(rotY,0,1,0);
        glRotated(rotZ,0,0,1);

        glScaled(0.4,0.4,0.4);
        //plane();
        if(car_flage==true)
            {
               car();
            }
        else if(car_flage1==true)
        {
            car1();
        }
        else if(car_flage2==true)
        {
            car2();
        }
        else if(car_flage3==true)
        {
            car3();
        }

    glPopMatrix();

    ///Environment
    if(tX>=6){START=false;START1= true;tX=0;}//tX=4.1;
    if(tX<=-5.1) {START=false;START1= true;tX=0;}

    glPushMatrix();
        glTranslated(tX,tY,tZ);
        football();// field not here
    glPopMatrix();

    glPushMatrix();
        glTranslated(tX,tY,tZ1);
        soheedMinarEnv();
    glPopMatrix();

    glPushMatrix();
        glTranslated(tX,tY,tZ2);
        home();
    glPopMatrix();

    glPushMatrix();
        glTranslated(tX,tY,tZ3);
        banch();
    glPopMatrix();

    glPushMatrix();
        glTranslated(tX,tY,tZ4);
        banch2();
    glPopMatrix();

    glPushMatrix();
        glTranslated(tX,tY,tZ5);
        tree();
    glPopMatrix();

    glPushMatrix();
        glTranslated(tX,tY,tZ6);
        truck();
    glPopMatrix();

    /// point collect portion
    //for 1

     if( tZ<0.5 && tZ>-1)
    {
        if(tX<3.5 && tX> 0.1)
        {
            START= false;
            START1= true;
            tX=0;
        }

    }
    ///
        if( tZ2<0.5 && tZ2>-1)
    {
        if(tX<0.5 && tX> -3.5)
        {
            START= false;
            START1= true;
            tX=3;
        }

    }
    ///

        if( tZ3<0.5 && tZ3>-2)
    {
        if(tX<1 && tX> -0.1)
        {
            START= false;
            START1= true;
            tX=3;
        }

    }
    ///
    if( tZ4<0.3 && tZ4>0)
    {
        if(tX<1 && tX> -0.1)
        {
            score=score+5;
        }

    }
    ///

    if( tZ5<0.3 && tZ5>0)
    {
        if((tX<6 && tX> 4))
        {
            score=score+5;
        }
        else if((tX<-1.5 && tX> -5.1))
        {
            START= false;
            START1= true;
            tX=0;

        }

    }

        if( tZ6<0.3 && tZ6>0)
    {
        if((tX<6 && tX> 4))
        {
            START= false;
            START1= true;
            tX=0;
        }
        else if((tX<-1.5 && tX> -5.1))
        {
            score=score+5;
        }

    }


    tZ+=speed;
    tZ1+=speed;
    tZ2+=speed;
    tZ3+=speed;
    tZ4+=speed;
    tZ5+=speed;
    tZ6+=speed;
///front gare
    if(tZ>=20)tZ=-110;
    if(tZ1>=20)tZ1=-110;
    if(tZ2>=20)tZ2=-110;
    if(tZ3>=20)tZ3=-110;
    if(tZ4>=20)tZ4=-110;
    if(tZ5>=20)tZ5=-110;
    if(tZ6>=20)tZ6=-110;
 /// back gare


    if(rotX>0)rotX-=angleBackFrac;
    if(rotX<0)rotX+=angleBackFrac;
    if(rotY>0)rotY-=angleBackFrac;
    if(rotY<0)rotY+=angleBackFrac;
    if(rotZ>0)rotZ-=angleBackFrac;
    if(rotZ<0)rotZ+=angleBackFrac;

}



void drawBitmapText(char *str,float x,float y,float z)
{
	char *c;
	glRasterPos3f(x,y+8,z);

	for (c=str; *c != '\0'; c++)
	{
		glutBitmapCharacter(GLUT_BITMAP_TIMES_ROMAN_10, *c);
	}
}

void drawStrokeText(char* str,int x,int y,int z)
{
	  char *c;
	  glPushMatrix();
	  glTranslatef(x, y+8,z);
	  glScalef(0.003f,0.003f,z);

	  for (c=str; *c != '\0'; c++)
	  {
    		glutStrokeCharacter(GLUT_STROKE_ROMAN , *c);
	  }
	  glPopMatrix();
}

void drawStrokeText2(char* str,int x,int y,int z)
{
	  char *c;
	  glPushMatrix();
	  glTranslatef(x, y+8,z);
	  glScalef(0.005f,0.005f,z);

	  for (c=str; *c != '\0'; c++)
	  {
    		glutStrokeCharacter(GLUT_STROKE_ROMAN , *c);
	  }
	  glPopMatrix();
}
void drawStrokeChar(char c,float x,float y,float z)
{
	  glPushMatrix();
          glTranslatef(x, y+8,z);
          glScalef(0.003f,0.003f,z);
          glutStrokeCharacter(GLUT_STROKE_ROMAN , c);
	  glPopMatrix();
}

///light 02
GLfloat light_ambients[] = {0.3, 0.3, 0.3, 1};
GLfloat light_diffuses[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat light_speculars[] = { 1.0, 1.0, 1.0, 1.0 };
GLfloat light_positions[] = { 0,0.5,4.5 ,1.0 };

static void display(void)
{
    const double t = glutGet(GLUT_ELAPSED_TIME) / 1000.0;
    double a = t*90.0;
    double aa=a;

    if(!rot){
        a=0;
    }

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();

    gluLookAt(eye_x,eye_y,eye_z,c_x,c_y,c_z,up_x,up_y,up_z);
    //gluLookAt(	0.0, 14.5, 30.0,0, 4, 0,0, 1.0f, 0.0f);



    if(START==true){
    //if(1){

        glPushMatrix();
            glTranslated(0,0,0);
            glScaled(zoom,zoom,zoom);
            glRotated(a,0,1,0);
            draw();
        glPopMatrix();
                glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient16);
                glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse16);
                glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular16);
                glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess16);
        drawStrokeText("  Speed UP: 1, Speed DOWN: 2, Break: 3, LEFT: A, RIGHT: D, MAIN MENU: G",-14,9,2);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient5);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse5);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular5);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess5);

        drawStrokeText("Score : ",3,9,2);
        //drawStrokeChar(score+48,4,9,2);
        int pp=score, mod,number=0;

        float tmp=0;
        //drawStrokeChar(10,4+tmp,9,2);
        //while(score){
        for( int i=0;i<3;i++){
            //mod=Time%10;
            mod=score%10;
            drawStrokeChar(mod+48,4.6-tmp,9,2);
            //drawStrokeChar(score+48,4+tmp,9,2);
            score/=10;

            tmp+=0.2;
        }
        score=pp;
    }
    else if(START1== true)
    {
        glPushMatrix();
            glTranslated(0,3,0);
            glRotated(aa,0,1,0);
            glScaled(1.5,1.5,1.5);
            //plane();
            car();
        glPopMatrix();

        drawStrokeText2("Press G to Restart the game again",-1,-1,0);
        //drawStrokeText2("Plane Game",-2,0,0);
    }
    else{

        glPushMatrix();
            glTranslated(0,3,0);
            glRotated(aa,0,1,0);
            glScaled(1.5,1.5,1.5);
            //plane();
            car();
        glPopMatrix();

        drawStrokeText("Press G to Start",0,-1,2);
        drawStrokeText2("Racing Car Game",0,0,2);
            glMaterialfv( GL_FRONT, GL_AMBIENT, mat_ambient16);
            glMaterialfv( GL_FRONT, GL_DIFFUSE, mat_diffuse16);
            glMaterialfv( GL_FRONT, GL_SPECULAR, mat_specular16);
            glMaterialfv( GL_FRONT, GL_SHININESS, mat_shininess16);
        drawStrokeText2("Submitted By :",-10,3.5,3);
        drawStrokeText2("________",-10,3,3);
        drawStrokeText2("Md. Mushfiqur Rahman Fakir",-10,2,3);
        drawStrokeText2("Roll: 1607116",-10,1,3);
        drawStrokeText2("Dept. of CSE",-10,0,3);

        drawStrokeText2("Submitted To :",8,3.5,3);
        drawStrokeText2("________",8,3,3);
        drawStrokeText2("Dr. Sk. Md. Masudul Ahsan",8,2,3);
        drawStrokeText2("Professor",8,1,3);
        drawStrokeText2("Dept. of CSE",8,0,3);
        drawStrokeText2("Kazi Saeed Alam",8,-1,3);
        drawStrokeText2("Lecturer",8,-2,3);
        drawStrokeText2("Dept. of CSE",8,-3,3);
    }

    glLightfv(GL_LIGHT1, GL_AMBIENT,  light_ambients);
    glLightfv(GL_LIGHT1, GL_DIFFUSE,  light_diffuses);
    glLightfv(GL_LIGHT1, GL_SPECULAR, light_speculars);
    glLightfv(GL_LIGHT1, GL_POSITION, light_positions);
    //glLightfv(GL_LIGHT1, GL_POSITION, light_positions);
    if(sp_flag==true)
    {
        glEnable(GL_LIGHT1);
        GLfloat l_spt[] = {0.0,0,-1,0.0};
        GLfloat spt_ct[] = {spt_cutoff};
        glLightfv(GL_LIGHT1, GL_SPOT_DIRECTION, l_spt);
        glLightfv(GL_LIGHT1, GL_SPOT_CUTOFF, spt_ct);
    }

    glutSwapBuffers();
}
///light 01
const GLfloat no_light[]  = { 0.0f, 0.0f, 0.0f, 1.0f };
const GLfloat light_ambient[]  = { 0.30f, 0.30f, 0.30f, 1.0f };
const GLfloat light_diffuse[]  = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat light_specular[] = { 1.0f, 1.0f, 1.0f, 1.0f };
const GLfloat light_position[] = { 2.0f, 5.0f, 2.0f, 0.0f };

static void idle(void)
{
    glutPostRedisplay();
}

void LoadTexture(const char*filename)
{
    glGenTextures(1, &ID);
    glBindTexture(GL_TEXTURE_2D, ID);
    glPixelStorei(GL_UNPACK_ALIGNMENT, ID);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    BmpLoader bl(filename);
    gluBuild2DMipmaps(GL_TEXTURE_2D, GL_RGB, bl.iWidth, bl.iHeight, GL_RGB, GL_UNSIGNED_BYTE, bl.textureData );
}
static void key(unsigned char key, int x, int y)
{
    float frac = 0.3;
    float rotFrac = 1;
    switch (key)
    {
        case 27 :
        case 'q':
            exit(0);
            break;
        case 'r':
            rot=true;
            break;
        case 't':
            rot=false;
            break;
        case 'z':
            zoom+=0.05;
            break;
        case 'Z':
            zoom-=0.05;
        case 'w':
            tY-=frac;
            rotZ+=rotFrac;
            break;
        case 's':
            tY+=frac;
            rotZ-=rotFrac;
            break;
        case 'a':
            tX+=frac;
            break;
        case 'd':
            tX-=frac;
            break;

        case 'g':
            START=true;
            START1=false;
            break;

        case 'G':
            START=false;
            START1=false;
            break;
        case 'm':
            START=false;
            break;
        case '7':
            eye_x+=0.1;
            break;
       case '4':
            eye_x-=0.1;
            break;

        case 'P':
            eye_x=0,eye_y=14.5,eye_z=30,c_x=0,c_y=4,c_z=0,up_x=0,up_y=1.0f,up_z=0.0f;
            break;
        case 'p':
            eye_x=0,eye_y=3,eye_z=1,c_x=0,c_y=4,c_z=-10,up_x=0,up_y=1.0f,up_z=0.0f;
            break;
        case '6':
            eye_y-=0.1;
            break;
        case '9':
            eye_y+=0.1;
            break;
        case '8':
            eye_z-=.1;
            break;
        case '5':
            eye_z+=.1;
            break;
        case '1':
            speed +=0.02;
            if(speed > 1.5) speed =1.5;
            break;
        case '2':
            speed -=.02;
            if(speed <-1) speed =0;
            break;
        case '3':
            speed =0;
            break;

        case 'v':
            glDisable(GL_LIGHT0);
            break;
        case 'V':
            glEnable(GL_LIGHT0);
            break;
        case 'b':
            glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ambient);
            glLightfv(GL_LIGHT0, GL_DIFFUSE,  no_light);
            glLightfv(GL_LIGHT0, GL_SPECULAR, no_light);
            break;

        case 'n':
            glLightfv(GL_LIGHT0, GL_AMBIENT,  no_light);
            glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diffuse);
            glLightfv(GL_LIGHT0, GL_SPECULAR, no_light);
            break;
        case 'M':
            glLightfv(GL_LIGHT0, GL_AMBIENT,  no_light);
            glLightfv(GL_LIGHT0, GL_DIFFUSE,  no_light);
            glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);;
            break;
        case 'l':
            sp_flag=true;
            break;

        case 'L':
            glDisable(GL_LIGHT1);
            sp_flag=false;
            break;
    case 'h':
    case 'H':
        wired=!wired;
        break;
    case 'c':
        car_flage= true;
        car_flage1=false;
        car_flage2=false;
        car_flage3=false;
        break;
    case 'C':
        car_flage= false;
        car_flage1=true;
        car_flage2=false;
        car_flage3=false;
        break;
    case 'f':
        car_flage= false;
        car_flage1=false;
        car_flage2=true;
        car_flage3=false;
        break;
    case 'F':
        car_flage= false;
        car_flage1=false;
        car_flage2=false;
        car_flage3=true;
        break;

    }

    glutPostRedisplay();
}

int main(int argc, char *argv[])
{
    glutInit(&argc, argv);
    glutInitWindowPosition(0,0);
	glutInitWindowSize(1366,720);
    glutInitDisplayMode(GLUT_RGB | GLUT_DOUBLE | GLUT_DEPTH | GLUT_RGBA);

    glutCreateWindow("GLUT Shapes");
    LoadTexture("C:\\Users\\User\\Desktop\\Online class 4_2\\Graphics lab\\Project\\Project with Texture\\football.bmp");
    LoadTexture("C:\\Users\\User\\Desktop\\Online class 4_2\\Graphics lab\\Project\\Project with Texture\\player.bmp");
    LoadTexture("C:\\Users\\User\\Desktop\\Online class 4_2\\Graphics lab\\Project\\Project with Texture\\building.bmp");
    LoadTexture("C:\\Users\\User\\Desktop\\Online class 4_2\\Graphics lab\\Project\\Project with Texture\\road.bmp");
    LoadTexture("C:\\Users\\User\\Desktop\\Online class 4_2\\Graphics lab\\Project\\Project with Texture\\brick.bmp");
    LoadTexture("C:\\Users\\User\\Desktop\\Online class 4_2\\Graphics lab\\Project\\Project with Texture\\home2.bmp");
    //LoadTexture("C:\\Users\\User\\Desktop\\Online class 4_2\\Graphics lab\\Project\\Project with Texture\\home3.bmp");
    LoadTexture("C:\\Users\\User\\Desktop\\Online class 4_2\\Graphics lab\\Project\\Project with Texture\\f2.bmp");


    glutReshapeFunc(resize);
    glutDisplayFunc(display);
    glutKeyboardFunc(key);
    glutMouseFunc(processMouse);
    glutIdleFunc(idle);

    //PlaySound("starwars.wav", NULL, SND_ASYNC|SND_FILENAME|SND_LOOP);

    //glClearColor(1,1,1,1);
    glEnable(GL_CULL_FACE);
    glCullFace(GL_BACK);

    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LESS);

    glEnable(GL_LIGHT0);
    glEnable(GL_NORMALIZE);
    //glEnable(GL_COLOR_MATERIAL);
    glEnable(GL_LIGHTING);

    glLightfv(GL_LIGHT0, GL_AMBIENT,  light_ambient);
    glLightfv(GL_LIGHT0, GL_DIFFUSE,  light_diffuse);
    glLightfv(GL_LIGHT0, GL_SPECULAR, light_specular);
    glLightfv(GL_LIGHT0, GL_POSITION, light_position);


    glutMainLoop();

    return EXIT_SUCCESS;
}
