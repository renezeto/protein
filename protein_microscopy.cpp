#include <iostream>
using namespace std;
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "protein.h"
#include "MersenneTwister.h"
#include <cassert>

const double difD = 2.5; // (um)^2 s^- 1
const double difE = 2.5; // (um)^2 s^-1
const double rate_ADP_ATP = 1; // s^-1
const double rate_D = .025; // um s^-1
const double rate_dD = .0015; // (um)^3 s^-1
const double rate_de = .7; // s^-1
const double rate_E = .093; // (um)^3 s^-1

const double nATP_starting_density = 1000.0/(M_PI*0.5*0.5); //proteins per micrometer^3 (values from paper)
const double nE_starting_density = 350.0/(M_PI*0.5*0.5); // proteins per micrometer
double density_factor;

int area_rating_flag = 0;
int slice_flag = 0;
int dump_flag = 0;

char* hires_flag_str = new char[1024];
char* slice_flag_str = new char[1024];

double dx; //grid spacing
double tot_time; //total simulation time
double time_step; //simulation time step
int iter; //total # of iterations
int printout_iterations; //iterations between file printout

int box_divider_left;
int box_divider_right;

double x, y, z;
int Nx, Ny, Nz; //number of gridpoints in each direction

string mem_f_shape; //cell shape argument
double A, B, C, D; //specific shape parameters, set by command line args

//N denotes protein number, n denotes protein number density
double *nATP; //min D bound to an ATP
double *nADP; //min D bound to an ADP
double *nE; //loose min E in cytoplasm
double *ND; //min D bound to ATP on the wall
double *NDE; //min D bound to ATP and min E on the wall
double *f_mem;

const int starting_num_guassians=20;
double guass[3*starting_num_guassians]; //stores y,z, and sigma for each guassian when creating random cell wall
const int random_num_guassians=5;
int rand_seed=0; //=14; at this point I have this passed in from the command line as the D argument
double Norm = 15.0; //This is the height of the guassians that make the cell wall

//begin randst shape functions
double rand_dis(double d0,double d_fac,int i) {
  int fac = 1;
  srand(i+rand_seed);
  double x = (rand()%1000);
  x = x/1000.0;
  if ( (rand()%1000) < 500) fac = -1;
  return d0*(1 + fac*(-d_fac*log(1-x)));
}

void randomize_cell_wall(double guass[]){
  //double X = Nx*dx; unused
  double Y = Ny*dx;
  double Z = Nz*dx;
  guass[0]=Y/2.0;
  guass[1]=Z/2.0;
  guass[2] = rand_dis(0.2,.13,2);
  for (int i = 1; random_num_guassians; i++){
    srand(i+rand_seed);
    double sigma = rand_dis(0.2,.13,i);
    guass[3*i+2]=sigma;
    double d = rand_dis(0.4,.13,i);
    double theta = fmod(rand(),2*M_PI);
    double y_change = d*sin(theta);
    double z_change = d*cos(theta);
    guass[i*3] = guass[(i-1)*3]+y_change;
    guass[i*3+1] = guass[(i-1)*3+1]+z_change;
  }
}

double f_2D_TIE_fighter(double y, double z){
  double f = 0;double f1 = 0;double f2 = 0;double f3 = 0;
  double Y = Ny*dx;
  double Z = Nz*dx;
  if (y<2.6){
    f1 = (z-Z/2.0)*(z-Z/2.0)/4 + (y-1.5)*(y-1.5)/0.6 - 1.0;
  } else {f1 = 0.1;}
  if (Y-y < 2.6){
    f2 = (z-Z/2.0)*(z-Z/2.0)/4 + (y-(Y-1.5))*(y-(Y-1.5))/0.6 - 1.0;
  } else {f2 = 0.1;}
  if (abs(z-Z/2.0) < 0.8 && abs(y-Y/2.0) < 2.0){
    f3 = abs(z-Z/2.0)-.5;
  } else {f3 = 0.1;}
  f = f1+f2+f3;
  return f;
}

bool only_print_once = true;
double f_2D_triangle(double y, double z){
  double Y = Ny*dx; double Z = Nz*dx; // total height and width of grid
  double z1 = A/2.0+2*dx; double y1 = A/2.0+2*dx; // lower left corner of triangle
  if (y < y1) {
    return 0.1; // it's too low to be in the triangle
  }

  double z2 = Z-A/2.0-2*dx; double y2 = y1; // lower right corner of triangle
  //Using law of cosines from lengths of sides we get top corner:
  double cos_theta = (B*B+D*D-C*C)/(2*B*D);
  double z3 = A/2.0+2*dx+D*cos_theta; double y3 = Y-A/2.0-2*dx; // top corner of triangle
  if (only_print_once==true){
    printf("z1 = %g y1 = %g\nz2 = %g y2 = %g\nz3 = %g y3 = %g\n",z1,y1,z2,y2,z3,y3);
    printf("cos_theta = %g\nZ = %g Y = %g A = %g",cos_theta,Z,Y,A/2.0);
    only_print_once = false;
  }

  if (z >= z3) {
    double fac = (z-z3)/(z2-z3); //how far along the line the z coordinate is
    double y_line = fac*(y2-y3)+y3; //the y coordinate of the line at the z point
    if (y > y_line) {
      return 0.1; // it's outside the triangle on the right side
    }
  }
  if (z < z3) {
    double fac = (z-z1)/(z3-z1); //how far along the line the z coordinate is
    double y_line = fac*(y3-y1)+y1; //the y coordinate of the line at the z point
    if (y > y_line) {
      return 0.1; // its outside the triangle on the left side
    }
  }
  //double rad = 1.75*(z2-z1)*sqrt(3.0)/6.0;
  //double y_circle = Y/2.0; double z_circle = z1 + sqrt(3)*(y2-y1)/6.0;
  //if ((z < zl1) && (z < zl3) && (z > z1) && ((z-z_circle)*(z-z_circle) + (y-y_circle)*(y-y_circle)) < rad*rad){
  //  return -0.1;
  //} else {
  return -0.1;
}

double f_2D_randst(double y, double z){
  int num_guassians = 0;
  for(int i=0;i<starting_num_guassians;i++){
    if (guass[i*3]!=0.0 || guass[i*3+1]!=0.0 || guass[i*3+2]!=0.0){
      num_guassians++;
    }
  }
  double f=0;
  for (int i = 0; i<num_guassians; i++){
    double arg = ((y-guass[3*i])*(y-guass[3*i]) + (z-guass[3*i+1])*(z-guass[3*i+1]))/(guass[3*i+2]*guass[3*i+2]);
    double fi = -((Norm*exp(-arg))-exp(-1.0));
    f+=fi;
  }
  return f;
}
//end randst shape functions

//mem_f produces a scalar field on the grid. points where mem_f = 0 are cell wall.
double mem_f(double x, double y, double z) {
  if(mem_f_shape=="randst" || mem_f_shape=="TIE_fighter" || mem_f_shape=="triangle"){
    double f = 0;
    double f0 = 0;
    double X = Nx*dx;
    double x1 = (X-A)/2.0;
    double x2 = (X+A)/2.0;
    if (x< x2 && x > x1){
      if(mem_f_shape=="randst") f = f_2D_randst(y,z);
      else if(mem_f_shape=="TIE_fighter") f = f_2D_TIE_fighter(y,z);
      else if(mem_f_shape=="triangle") f = f_2D_triangle(y,z);
      else {
        printf("somethings wrong with the shape argument!!!");
        exit(1);
      }
      if (f<=0) {
        f = abs(2*(x-(X/2))/A) - 1;
        return f;
      }
      fflush(stdout);
      if (f>0) {
        double closest_y0 = -100.0;
        double closest_z0 = -100.0;
        //bool there_is_closest_point=0; unused
        for (double y0 = y-(A/2.0+2.0*dx); y0<y+(A/2.0+2.0*dx); y0+=dx) {
          for (double z0 = z-(A/2.0+2.0*dx); z0<z+(A/2.0+2.0*dx); z0+=dx) {
            if ( (y-y0)*(y-y0)+(z-z0)*(z-z0) < (y-closest_y0)*(y-closest_y0)+(z-closest_z0)*(z-closest_z0) ) {
              if(mem_f_shape=="randst") f0 = f_2D_randst(y0,z0);
              if(mem_f_shape=="TIE_fighter") f0 = f_2D_TIE_fighter(y0,z0);
              if(mem_f_shape=="triangle") f0 = f_2D_triangle(y0,z0);
              if (f0 <= 0) {              //there_is_closest_point = 1; unused
                closest_y0 = y0;
                closest_z0 = z0;
              }
            }
          }
        }
        double dis = sqrt((y-closest_y0)*(y-closest_y0) + (z-closest_z0)*(z-closest_z0) + (x-X/2.0)*(x-X/2.0));
        if (dis <= A) {
          f = 2*(dis/A-.5);
          return f;
        } else {
          return 1;
        }
      }
    } else if (x >= x2) {
      f = 2*(x-x2)/A;
      return f;
    } else {
      f = 2*(x1-x)/A;
      return f;
    }
  }

  if (mem_f_shape=="p"){
    //A = length, B = radius of endcap and cylinder
    double f;
    double X = Nx*dx;
    double Y = Ny*dx;
    double Z = Nz*dx;
    double z1 = (Z-A)/2;
    double z2 = (A+(Z-A)/2);
    double x1 = X/2;
    double y1 = Y/2;
    f = sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1))-B;
    if (z < z1) {
      f = sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1)+(z-z1)*(z-z1))-B;
    }
    if (z > z2) {
      f = sqrt((x-x1)*(x-x1)+(y-y1)*(y-y1)+(z-z2)*(z-z2))-B;
    }
    return f;
  }

  if (mem_f_shape=="b"){
    //A (x),B (z),C (y) lengths
    double f;
    double X = Nx*dx;
    double Y = Ny*dx;
    double Z = Nz*dx;
    if ( (x >= (X-A)/2) && (x <= A+(X-A)/2) && (z >= (Z-B)/2) && (z <= B+(Z-B)/2) && (y >= (Y-C)/2) && (y <= C+(Y-C)/2)) {
      f = -1;
    }
    else {
      f = 1;
    }
    return f;
  }

  if (mem_f_shape=="c"){
    //A = length of cone, B = radius of base
    double f;
    double X = Nx*dx;
    double Y = Ny*dx;
    double Z = Nz*dx;
    double x1=X/2;
    double y1=Y/2;
    double z1=(Z-A)/2;
    double z2=(A+(Z-A)/2);
    if ((z > z1) && (z < z2)){
      f = sqrt((x-x1)*(x-x1) + (y-y1)*(y-y1))/(z-z1) - B/(z2-z1);
    }
    else { f = 1; }
    return f;
  }

  if (mem_f_shape=="st"){
    //A = length, B = x axis radius radius, C = y axis radius radius, D = z axis radius radius
    double f;
    double X = Nx*dx;
    double Y = Ny*dx;
    double Z = Nz*dx;
    double z1 = (Z-A)/2;
    double z2 = (A+(Z-A)/2);
    double x1 = X/2;
    double y1 = Y/2;
    f = sqrt((x-x1)*(x-x1)/(B*B)+(y-y1)*(y-y1)/(C*C))-1;
    if (z < z1) {
      f = sqrt((x-x1)*(x-x1)/(B*B)+(y-y1)*(y-y1)/(C*C)+(z-z1)*(z-z1)/(D*D))-1;
    }
    if (z > z2) {
      f = sqrt((x-x1)*(x-x1)/(B*B)+(y-y1)*(y-y1)/(C*C)+(z-z2)*(z-z2)/(D*D))-1;
    }
    return f;
  }

  if (mem_f_shape=="sp"){
    // A = radius
    double X = Nx*dx;
    double Y = Ny*dx;
    double Z = Nz*dx;
    double x1 = X/2;
    double y1 = Y/2;
    double z1 = Z/2;
    double f;
    f = sqrt((x-x1)*(x-x1) + (y-y1)*(y-y1) + (z-z1)*(z-z1)) - A;
    return f;
  }

  if (mem_f_shape=="e"){
    //B = x axis radius radius, C = y axis radius radius, A = z axis radius radius
    double X = Nx*dx;
    double Y = Ny*dx;
    double Z = Nz*dx;
    double x1 = X/2;
    double y1 = Y/2;
    double z1 = Z/2;
    double f = sqrt( (x-x1)*(x-x1)/(B*B) + (y-y1)*(y-y1)/(C*C)+ (z-z1)*(z-z1)/(A*A) ) - 1;
    return f;
  }

  else {
    double f = 1;
  	return f;
  }
}


//struct for storing plot information
struct protein {
  char* name;

  //time_map
  double* sum;

  //box_plot
  double* numLeft;
  double* numMid;
  double* numRight;

  double* numRightUp;
  double* numRightDown;
  double* numLeftUp;
  double* numLeftDown;

  //arrow_plot
  double* maxloc;
};

char* print_filename(const char* plotname, const char* proteinname) {
  char* filename = new char[1024];
  sprintf(filename,"data/shape-%s/%s%s%s-%s-%s-%1.2f-%1.2f-%1.2f-%1.2f-%1.2f.dat",mem_f_shape.c_str(),hires_flag_str,slice_flag_str,plotname,proteinname,mem_f_shape.c_str(),A,B,C,D,density_factor);
  return filename;
}

bool only_once = true;
string triangle_section (double y, double z) {
  //needs editing!!!!!
  double Y = Ny*dx; double Z = Nz*dx; // total height and width of grid
  double z1 = A/2.0+2*dx; double y1 = A/2.0+2*dx; // lower left corner of triangle
  double z2 = Z-A/2.0-2*dx; double y2 = y1; // lower right corner of triangle
  //Using law of cosines from lengths of sides we get top corner:
  double cos_theta = (B*B+D*D-C*C)/(2*B*D);
  double z3 = A/2.0+2*dx+D*cos_theta; double y3 = Y-A/2.0-2*dx; // top corner of triangle
  if (only_once == true) {
    printf("z1 = %g y1 = %g\nz2 = %g y2 = %g\nz3 %g y3 = %g\n",z1,y1,z2,y2,z3,y3);
  }

  //get bisecting points and lines:
  double y_21 = (y1 + y2)/2.0; double z_21 = (z2 + z1)/2.0;
  double y_32 = (y3 + y2)/2.0; double z_32 = (z3 + z2)/2.0;
  double y_13 = (y1 + y3)/2.0; double z_13 = (z1 + z3)/2.0;
  double slope1 = (y_32-y1)/(z_32-z1); // from left corner to right line
  double slope2 = (y2-y_13)/(z2-z_13); //from right corner to left line
  double slope3 = (y_21-y3)/(z_21-z3); //from top corner to bottom
  //running into nan issues when z3 is same as z_21, so I brute force a large negative slope:
  if (z_21-z3 < .000001){
    slope3 = -1000000;
  }
  if (only_once==true){
    printf("slope1 = %g slope2 = %g slope3 = %g\n",slope1,slope2,slope3);
  }
  //find centroid, which is where all three lines above intersect:
  double z_cen = (y3 - y1 + slope1*z1 - slope3*z3)/(slope1 -slope3);
  double y_cen = slope1*(z_cen-z1) + y1;
  if (only_once ==true){
    printf("z_cen = %g y_cen = %g\n",z_cen,y_cen);
    only_once = false;
  }
  if (z > z_cen) {
    if (y > slope1*(z-z_cen)+y_cen) {
      return "Mid";
    }
  } else {
    if (y > slope2*(z-z_13)+y_13) {
      return "Mid";
    }
  }
  if (y > slope3*(z-z_cen)+y_cen) {
    return "Right";
  }
  return "Left";
}

int main (int argc, char *argv[]) {
  //command line parameters
  mem_f_shape = argv[1];
  A = atof(argv[2]);
  B = atof(argv[3]);
  C = atof(argv[4]);
  D = atof(argv[5]);
  density_factor = atof(argv[6]);
  dx=.05;

  //flag checking
  for (int i=0; i<argc; i++) {
    if (strcmp(argv[i],"-area")==0) {
      area_rating_flag = 1;
      printf("Area rating printout.\n");
    }
    if (strcmp(argv[i],"-hires")==0) {
      dx=.025;
      printf("Using high resolution.\n");
      sprintf(hires_flag_str,"hires-");
    }
    if (strcmp(argv[i],"-slice")==0) {
      slice_flag = 1;
      printf("Printing middle slice data.\n");
      sprintf(slice_flag_str,"slice-");
    }
    if (strcmp(argv[i],"-dump")==0) {
      dump_flag = 1;
      printf("Printing all 501 data files.\n");
    }
  }

  //fixed simulation parameters
  tot_time = 2500; //sec
  time_step = .1*dx*dx/difD;//sec
  iter = int(20*1000);
  //iter = int(tot_time/time_step);
  printout_iterations = int(5.0/time_step);
  printf("%d\n",printout_iterations);//approximately 5 seconds between each printout
  double dV = dx*dx*dx;

  //compute grid size based on cell parameters
  if (mem_f_shape=="p") {
    Nx = ceil(2*B/dx) + 4;
    Ny = ceil(2*B/dx) + 4;
    Nz = ceil((A + 2*B)/dx) + 4;
    box_divider_left = int(Nz/3);
    box_divider_right = int(2*Nz/3);
  }
  if (mem_f_shape=="b") {
    Nx = ceil(A/dx) + 4;
    Ny = ceil(C/dx) + 4;
    Nz = ceil(B/dx) + 4;
  }
  if (mem_f_shape=="c") {
    Nx = ceil(2*B/dx) + 4;
    Ny = ceil(2*B/dx) + 4;
    Nz = ceil(A/dx) + 4;
  }
  if (mem_f_shape=="randst") {
    Nx = ceil(A/dx) + 4;
    Ny = ceil(B/dx) + 4;
    Nz = ceil(C/dx) + 4;
    rand_seed = int(D);
    if (D != round(D)) {
      printf("WARNING!!! When using randst the last argument, the rand_seed, should be an integer!  For now I've truncated it!!!\n");
    }
  }
  if (mem_f_shape=="TIE_fighter") {
    Nx = ceil(A/dx) + 4;
    Ny = ceil(B/dx) + 4;
    Nz = ceil(C/dx) + 4;
  }
  if (mem_f_shape=="triangle") {
    Nx = ceil(A/dx) + 4;
    Nz = ceil((A+B)/dx) + 4;
    //Using law of cosines we get height of triangle:
    double theta = acos((B*B+D*D-C*C)/(2*B*D));
    Ny = ceil((A+D*sin(theta))/dx) + 4;
  }
  if (mem_f_shape=="st") {
    Nx = ceil(2*B/dx) + 4;
    Ny = ceil(2*C/dx) + 4;
    Nz = ceil((A+2*D)/dx) + 4;
  }
  if (mem_f_shape=="sp") {
    Nx = ceil(2*A/dx) + 4;
    Ny = ceil(2*A/dx) + 4;
    Nz = ceil(2*A/dx) + 4;
  }
  if (mem_f_shape=="e") {
    Nx = ceil(1/dx) + 4;
    Ny = ceil(2*A/dx) + 4;
    Nz = ceil(2*B/dx) + 4;
  }

  //open out file to begin recording info about simulation
  char * out_file_name = new char[1024];
  sprintf(out_file_name,"data/shape-%s/out_files/%s-%4.02f-%4.02f-%4.02f-%4.02f-%4.02f.out",mem_f_shape.c_str(),mem_f_shape.c_str(),A,B,C,D,density_factor);
  FILE * out_file = fopen((const char *)out_file_name,"w");

  //random stuff that needs to be contained
  fprintf(out_file,"Nx=%d\nNy=%d\nNz=%d\nX=%f\nY=%f\nZ=%f\n",Nx,Ny,Nz,(Nx*dx),(Ny*dx),(Nz*dx));
  for (int i=0;i<3*starting_num_guassians;i++){
    guass[i]=0;
  }

  //In the following, for every set of three numbers, the 1st is y and he 2nd is z and the 3rd is quassian width
  double guass99[] = {2.0,2.2,.50,3,3,.50,4.0,3.6,.50,3,4.2,.50,2.0,5,.50};
  double guass98[] = {2.0,2.0,.3,3,3,.6,4.2,3.4,.3,4.6,4.6,.6,3.4,5.6,.6};
  double guass97[] = {1.4,3,.4,1.8,3,.4,2.2,3,.4,2.6,3,.4,3,3,.4,3.4,3,.4,
                      3.8,3,.4,4.2,3,.4,4.6,3,.4,5,3,.4,5.4,3,.4,3.4,2.4,.6};
  double guass96[] = {1.3,1.3,.7,2.1,2,.7,3,2,.7,3.9,2,.7,4.7,1.3,.7,4,2.1,.7,4,3,.7,4,3.9,.7,
                      4.7,4.7,.7,3.9,4,.7,3,4,.7,2.3,4,.7,1.3,4.7,.7,2.1,3.9,.7,3,3.9,.7,2.1,3.9,.7};
  if (rand_seed == 99){
    for (int i=0;i<3*5;i++){
      guass[i]=guass99[i];
      fprintf(out_file,"rand_seed is 99!");
      fflush(stdout);
    }
  } else if (rand_seed == 98){
    for (int i=0;i<3*5;i++){
      guass[i]=guass98[i];
    }
  } else if (rand_seed == 97){
    for (int i=0;i<3*12;i++){
      guass[i]=guass97[i];
    }
  } else if (rand_seed == 96){
    for (int i=0;i<3*16;i++){
      guass[i]=guass96[i];
    }
  } else {
    if (mem_f_shape == "randst"){
      fprintf(out_file,"rand_seed is not 99!");
      fflush(stdout);
      randomize_cell_wall(guass);
    }
  }
  //end random stuff

  //global arrays for storing simulation data
  nATP = new double[Nx*Ny*Nz];
  nADP = new double[Nx*Ny*Nz];
  nE = new double[Nx*Ny*Nz];
  ND = new double[Nx*Ny*Nz];
  NDE = new double[Nx*Ny*Nz];
  f_mem = new double[Nx*Ny*Nz];
  fprintf(out_file,"For this simulation,\ndx = %f\ntot_time = %f\ntimestep = %f\ntotal iterations = %d\niter at five sec = %d\n",
          dx, tot_time, time_step, iter, printout_iterations);
  double *JxATP = new double[Nx*Ny*Nz];
  double *JyATP = new double[Nx*Ny*Nz];
  double *JzATP = new double[Nx*Ny*Nz];
  double *JxADP = new double[Nx*Ny*Nz];
  double *JyADP = new double[Nx*Ny*Nz];
  double *JzADP = new double[Nx*Ny*Nz];
  double *JxE = new double[Nx*Ny*Nz];
  double *JyE = new double[Nx*Ny*Nz];
  double *JzE = new double[Nx*Ny*Nz];
  double *mem_A = new double[Nx*Ny*Nz];
  double *curvature = new double[Nx*Ny*Nz];
  bool *insideArr = new bool[Nx*Ny*Nz];
  for (int i=0;i<Nx*Ny*Nz;i++){mem_A[i] = 0;}

  const int numProteins = 5;

  protein* nATP_plot = new protein;
  protein* nE_plot = new protein;
  protein* nADP_plot = new protein;
  protein* NDE_plot = new protein;
  protein* ND_plot = new protein;

  protein* proteinList[numProteins] = { nATP_plot, nADP_plot, nE_plot, ND_plot, NDE_plot };
  double* accessGlobals[numProteins] = { nATP, nADP, nE, ND, NDE };

  int print_denominator = 1000;

  //initialize things
  for (int pNum=0; pNum<numProteins; pNum++) {
    proteinList[pNum]->sum = new double[Ny*Nz];
    proteinList[pNum]->name = new char[1024];

    proteinList[pNum]->numLeft = new double[iter];
    proteinList[pNum]->numMid = new double[iter];
    proteinList[pNum]->numRight = new double[iter];

    proteinList[pNum]->numRightUp = new double[iter];
    proteinList[pNum]->numRightDown = new double[iter];
    proteinList[pNum]->numLeftUp = new double[iter];
    proteinList[pNum]->numLeftDown = new double[iter];
  }

  sprintf(proteinList[0]->name,"D_nATP");
  sprintf(proteinList[1]->name,"D_nADP");
  sprintf(proteinList[2]->name,"E_nE");
  sprintf(proteinList[3]->name,"D_ND");
  sprintf(proteinList[4]->name,"D_E_NDE");

    set_membrane(out_file, mem_f, mem_A);
  set_curvature(mem_A,curvature);

  //begin area rating
  // printf("Starting with the area rating stuff!!!!\n");
  // fflush(stdout);
  // char *curvature_out = new char[1024];
  // if (dx==.05) {
  //   sprintf(curvature_out, "data/shape-%s/hires-curvature-%4.02f-%4.02f-%4.02f-%4.02f-%4.02f.dat",mem_f_shape.c_str(),A,B,C,D,density_factor);
  // }
  // else {
  //   sprintf(curvature_out, "data/shape-%s/curvature-%4.02f-%4.02f-%4.02f-%4.02f-%4.02f.dat",mem_f_shape.c_str(),A,B,C,D,density_factor);
  // }
  // FILE *curvature_file = fopen((const char *)curvature_out,"w");
  // if (curvature_file == NULL){
  //   printf("Error: curvature_file == null \n");
  // }

  // char *area_rating_out = new char[1024];
  // if (dx==.05) {
  //   sprintf(area_rating_out, "data/shape-%s/hires-area_rating-%4.02f-%4.02f-%4.02f-%4.02f-%4.02f.dat",mem_f_shape.c_str(),A,B,C,D,density_factor);
  // }
  // else {
  //   sprintf(area_rating_out, "data/shape-%s/area_rating-%4.02f-%4.02f-%4.02f-%4.02f-%4.02f.dat",mem_f_shape.c_str(),A,B,C,D,density_factor);
  // }
  // FILE *area_rating_file = fopen((const char *)area_rating_out,"w");
  // if (area_rating_file == NULL){
  //   printf("Error: area_rating_file == null \n");
  // }

  // char * area_rating_out_two = new char[1024];
  // if (dx==.05) {
  //   sprintf(area_rating_out_two, "data/shape-%s/hires-area_rating_two-%4.02f-%4.02f-%4.02f-%4.02f-%4.02f.dat",mem_f_shape.c_str(),A,B,C,D,density_factor);
  // }
  // else {
  //   sprintf(area_rating_out_two, "data/shape-%s/area_rating_two-%4.02f-%4.02f-%4.02f-%4.02f-%4.02f.dat",mem_f_shape.c_str(),A,B,C,D,density_factor);
  // }
  // FILE *area_rating_file_two = fopen((const char *)area_rating_out_two,"w");
  // if (area_rating_file_two == NULL){
  //   printf("Error: area_rating_file_two == null \n");
  // }

  // char * area_rating_out_three = new char[1024];
  // if (dx==.05) {
  //   sprintf(area_rating_out_three, "data/shape-%s/hires-area_rating_three-%4.02f-%4.02f-%4.02f-%4.02f-%4.02f.dat",mem_f_shape.c_str(),A,B,C,D,density_factor);
  // }
  // else {
  //   sprintf(area_rating_out_three, "data/shape-%s/area_rating_three-%4.02f-%4.02f-%4.02f-%4.02f-%4.02f.dat",mem_f_shape.c_str(),A,B,C,D,density_factor);
  // }
  // FILE *area_rating_file_three = fopen((const char *)area_rating_out_three,"w");
  // if (area_rating_file_three == NULL){
  //   printf("Error: area_rating_file_three == null \n");
  // }

  // fprintf(out_file,"Finished opening area_rating file.\n");
  set_insideArr(insideArr);
  // fprintf(out_file,"Finished with insideArr function.\n");

  double total_cell_volume = 0;
  double total_cell_area = 0;

  for (int i=0;i<Nx*Ny*Nz;i++){
    total_cell_area += mem_A[i];
    if (insideArr[i]==true) {
      total_cell_volume += dx*dx*dx;
    }
  }
  printf("total_cell_volume = %g\n",total_cell_volume);
  fprintf(out_file, "total_cell_volume = %g\n",total_cell_volume);

  //add cell params file

  // int i=int(Nx/2);
  // printf("Hello!!!Nx=%d, x=%g, Nx/2=%d, x/2=%g\n",Nx,Nx*dx,i,Nx*dx/2.0);
  // fflush(stdout);
  // for (int j=0;j<Ny;j++){
  //   for (int k=0;k<Nz;k++){

  //     double area_rating = 0;
  //     double area_rating_two = 0;
  //     double area_rating_three = 0;
  //     fprintf(curvature_file,"%g\t%g\t%g\t%g\t%g\n",i*dx,j*dx,k*dx,curvature[i*Ny*Nz+j*Nz+k],mem_A[i*Ny*Nz+j*Nz+k]);
  //     if (insideArr[i*Ny*Nz+j*Nz+k]==true){
  //       for (int i2=0;i2<Nx;i2++){
  //         for (int j2=0;j2<Ny;j2++){
  //           for (int k2=0;k2<Nz;k2++){
  //             if(i2!=i && j2!=j && k2!=k){
  //               double dis = dx*sqrt((i-i2)*(i-i2)+(j-j2)*(j-j2)+(k-k2)*(k-k2));
  //               area_rating += curvature[i2*Ny*Nz+j2*Nz+k2]/(dis*dis);
  //               if (dis<1.5*A){
  //                 area_rating_two += curvature[i2*Ny*Nz+j2*Nz+k2]/(dis*dis);
  //                 area_rating_three += curvature[i2*Ny*Nz+j2*Nz+k2];
  //               }
  //             }
  //           }
  //         }
  //       }
  //       fprintf(area_rating_file,"%g\t%g\t%g\t%g\t%g\n",i*dx,j*dx,k*dx,area_rating,mem_A[i*Ny*Nz+j*Nz+k]);
  //       fprintf(area_rating_file_two,"%g\t%g\t%g\t%g\t%g\n",i*dx,j*dx,k*dx,area_rating_two,mem_A[i*Ny*Nz+j*Nz+k]);
  //       fprintf(area_rating_file_three,"%g\t%g\t%g\t%g\n",i*dx,j*dx,k*dx,area_rating_three);
  //     } else {
  //       fprintf(area_rating_file,"%g\t%g\t%g\t%g\n",i*dx,j*dx,k*dx,0.0);
  //       fprintf(area_rating_file_two,"%g\t%g\t%g\t%g\n",i*dx,j*dx,k*dx,0.0);
  //       fprintf(area_rating_file_three,"%g\t%g\t%g\t%g\n",i*dx,j*dx,k*dx,0.0);
  //     }
  //   }
  // }
  // printf("Done with the area rating stuff!!!!\n");
  // fflush(stdout);

  // fclose(curvature_file);
  // fclose(area_rating_file_three);
  // fclose(area_rating_file_two);
  // fclose(area_rating_file);

  // if(area_rating_flag==1) {
  //   printf("the area rating flag was on so we're quitting the simulation now!!!\n");
  //   exit(0);
  // }
  // fprintf(out_file,"Area_rating_two file is using %g as the radius of the sphere.",A);
  // fprintf(out_file,"Finished writing to the area_rating file.\n");
  // fprintf(out_file,"Total cell volume = %g.\nTotal cell area = %g.\n",total_cell_volume,total_cell_area);
  // fflush(out_file);
  //end area rating

  fflush(out_file);

  //begin membrane printing - need to change this to mem_f instead of 1's and 0's
  printf("HELLLLOOOOOOO %s\n",mem_f_shape.c_str());
  char* outfilename = new char[1024];
  sprintf(outfilename,"data/shape-%s/membrane_files/membrane-%4.02f-%4.02f-%4.02f-%4.02f-%4.02f.dat",mem_f_shape.c_str(),A,B,C,D,density_factor);
  FILE *out = fopen((const char *)outfilename,"w");
  if (out==0){
    printf ("couldn't print outfile\n");
    //exit();
  }
  double marker;
  //  double inmarker; unused
  //  double zt = A/2; double yt = B/2; double xt = C/2; unused
  //  double ft = mem_f(zt,yt,xt); unused
  for (int j=0;j<Ny;j++){
    for (int i=0;i<Nz;i++){
      if (mem_A[(int(Nx/2))*Ny*Nz+j*Nz+i]!=0) {
        marker = 1;
      }
      else {
        marker = 0;
      }
      fprintf(out, "%g  ", marker);
    }
    fprintf(out, "\n");
  }
  fflush(stdout);
  fclose(out);
  printf("\nMembrane file printed.\n");

  char* outfilename_sections = new char[1024];
  sprintf(outfilename_sections, "data/shape-%s/membrane_files/sections-%4.02f-%4.02f-%4.02f-%4.02f-%4.02f.dat",mem_f_shape.c_str(),A,B,C,D,density_factor);
  FILE *outfile_sections = fopen((const char*)outfilename_sections,"w");
  for (int j=0;j<Ny;j++){
    for (int i=0;i<Nz;i++) {
      if (triangle_section(j*dx,i*dx)=="Right"){
        marker = 1;
      }
      if (triangle_section(j*dx,i*dx)=="Mid"){
        marker = 2;
      }
      if (triangle_section(j*dx,i*dx)=="Left"){
        marker = 3;
      }
      fprintf(outfile_sections,"%g ",marker);
    }
    fprintf(outfile_sections,"\n");
  }
  fflush(stdout);
  fclose(outfile_sections);
  printf("\nMembrane sections file printed.\n");

//end membrane printing

  //eventually uncommend and replace membrane.dat prints with mem_f prints for extrema.py
  //begin mem_f printing for randst, tie fighter, triangle - possibly do this for other shapes if needed -- NEEDS UPDATE.
  // if (mem_f_shape == "randst"||mem_f_shape == "TIE_fighter"||mem_f_shape == "triangle") {
  //   char *f_file_name = new char[1024];
  //   if(f_file_name==NULL) {
  //     fprintf(out_file,"Error: f_file_name is null.");
  //     exit(1);
  //   }
  //   sprintf(f_file_name,"data/shape-%s/membrane_files/f_membrane-%4.02f-%4.02f-%4.02f-%4.02f-%4.02f.dat", mem_f_shape.c_str(),A,B,C,D,density_factor);
  //   FILE *f_file = fopen((const char *)f_file_name,"w");
  //   double x = Nx/2.0*dx;
  //   for (int i=0;i<Ny;i++) {
  //     for (int j=0;j<Nz;j++) {
  //       fprintf(f_file,"%g\t%g\t%g\n",i*dx,j*dx,mem_f(x,i*dx,j*dx));
  //     }
  //   }
  //   fclose(f_file);
  //   fprintf(out_file,"Finished printing out the mem_f_shape function.\n");
  //   fflush(stdout);
  // }
  // fprintf (out_file,"Membrane set with density in it.\n");
  // fflush(out_file);
  //end mem_f printing

  set_density(nATP, nE, mem_A);
  fflush(out_file);
  printf("density set.\n");

  //   Reference for creating the randst dividers:
  //   the arrays are structured in the order of y,z,width, etc.  of the guassians that define the shape
  //   double guass99[] = {2.0,2.2,.50,3,3,.50,4.0,3.6,.50,3,4.2,.50,2.0,5,.50};
  //   double guass98[] = {2.0,2.0,.3,3,3,.6,4.2,3.4,.3,4.6,4.6,.6,3.4,5.6,.6};
  //   double guass97[] = {1.4,3,.4,1.8,3,.4,2.2,3,.4,2.6,3,.4,3,3,.4,3.4,3,.4,
  //                       3.8,3,.4,4.2,3,.4,4.6,3,.4,5,3,.4,5.4,3,.4,3.4,2.4,.6};
  //   double guass96[] = {1.3,1.3,.7,2.1,2,.7,3,2,.7,3.9,2,.7,4.7,1.3,.7,4,2.1,.7,4,3,.7,4,3.9,.7,
  //                       4.7,4.7,.7,3.9,4,.7,3,4,.7,2.3,4,.7,1.3,4.7,.7,2.1,3.9,.7,3,3.9,.7,2.1,3.9,.7};


  double vert_div;
  double vert_div_two;
  double hor_div;
  double hor_div_two;


  if (mem_f_shape=="randst") {
    if (rand_seed == 99) {
      vert_div = 2.6/dx;
      vert_div_two = 4.6/dx;
    }
    if (rand_seed == 98) {
      vert_div = 3.0/dx;
      vert_div_two = 4.8/dx;
    }
    if (rand_seed == 97) {
      vert_div = 2.6/dx;
      hor_div = 1.8/dx;
      hor_div_two = 3.0/dx;
    }
    if (rand_seed == 96) {
      vert_div = 1.7/dx;
      hor_div = 1.7/dx;
    }
  }


  //begin simulation
  for (int i=0;i<iter;i++){
    get_J(difD, nATP, nADP, nE, JxATP, JyATP,
          JzATP, JxADP, JyADP, JzADP, JxE, JyE, JzE);
    get_next_density(mem_A, insideArr, nATP, nADP, nE, ND, NDE, JxATP, JyATP, JzATP,
                     JxADP, JyADP, JzADP, JxE, JyE, JzE);
    if (i%print_denominator==0) {
      printf("Finished sim loop # i=%d, We're %1.2f percent done\n",i,double(100*i/iter));
    }

    //capture plot information at each time step -- need to work after first 10%
    for (int pNum=0; pNum<numProteins; pNum++) {

      //time map
      for (int a=0; a<Ny; a++) {
        for (int b=0; b<Nz; b++) {
          if (slice_flag==0) {
            for (int c=0; c<Nx; c++) {
              proteinList[pNum]->sum[a*Nz+b] += accessGlobals[pNum][c*Ny*Nz+a*Nz+b];
            }
          }
          else {
            proteinList[pNum]->sum[a*Nz+b] += accessGlobals[pNum][int(Nx/2)*Ny*Nz+a*Nz+b];
          }
        }
      }

      //box plot ...
      if (i%print_denominator==0) {
        if ((strcmp(proteinList[pNum]->name,"D_nATP")==0) || (strcmp(proteinList[pNum]->name,"E_nE")==0) || (strcmp(proteinList[pNum]->name,"D_nADP")==0)) {
          dV = dx*dx*dx;
        }
        else {
          dV = 1;
        }
        int i_dat = i/print_denominator;
        for (int a=0; a<Ny; a++) {
          for (int b=0; b<Nz; b++) {
            for (int c=0; c<Nx; c++) {
              if (mem_f_shape == "p") {
                if (b < box_divider_left) {
                  proteinList[pNum]->numLeft[i_dat] += accessGlobals[pNum][c*Ny*Nz+a*Nz+b]*dV;
                }
                if (b > box_divider_right) {
                  proteinList[pNum]->numRight[i_dat] += accessGlobals[pNum][c*Ny*Nz+a*Nz+b]*dV;
                }
                if ((b <= box_divider_right) && (b >= box_divider_left)) {
                  proteinList[pNum]->numMid[i_dat] += accessGlobals[pNum][c*Ny*Nz+a*Nz+b]*dV;
                }
              }

              if (mem_f_shape == "randst") {
                if (rand_seed == 99 || rand_seed == 98) {
                  if (b < vert_div) {
                    proteinList[pNum]->numLeft[i_dat] += accessGlobals[pNum][c*Ny*Nz+a*Nz+b]*dV;
                  }
                  else if (b > vert_div_two) {
                    proteinList[pNum]->numRight[i_dat] += accessGlobals[pNum][c*Ny*Nz+a*Nz+b]*dV;
                  }
                  else {
                    proteinList[pNum]->numMid[i_dat] += accessGlobals[pNum][c*Ny*Nz+a*Nz+b]*dV;
                  }
                }
                if (rand_seed == 97) {
                  if (b < vert_div) {
                    proteinList[pNum]->numLeft[i_dat] += accessGlobals[pNum][c*Ny*Nz+a*Nz+b]*dV;
                  }
                  else if (b > vert_div && a > hor_div_two) {
                    proteinList[pNum]->numRightUp[i_dat] += accessGlobals[pNum][c*Ny*Nz+a*Nz+b]*dV;
                  }
                  else if (b > vert_div && a < hor_div) {
                    proteinList[pNum]->numRightDown[i_dat] += accessGlobals[pNum][c*Ny*Nz+a*Nz+b]*dV;
                  }
                  else {
                    proteinList[pNum]->numMid[i_dat] += accessGlobals[pNum][c*Ny*Nz+a*Nz+b]*dV;
                  }
                }
                if (rand_seed == 96) {
                  if (b < vert_div && a < hor_div) {
                    proteinList[pNum]->numLeftDown[i_dat] += accessGlobals[pNum][c*Ny*Nz+a*Nz+b]*dV;
                  }
                  else if (b < vert_div_two && a > hor_div) {
                    proteinList[pNum]->numLeftUp[i_dat] += accessGlobals[pNum][c*Ny*Nz+a*Nz+b]*dV;
                  }
                  else if (b > vert_div_two && a < hor_div) {
                    proteinList[pNum]->numRightDown[i_dat] += accessGlobals[pNum][c*Ny*Nz+a*Nz+b]*dV;
                  }
                  else {
                    proteinList[pNum]->numRightUp[i_dat] += accessGlobals[pNum][c*Ny*Nz+a*Nz+b]*dV;
                  }
                }
              }
              if (mem_f_shape == "triangle") {
                if (triangle_section(a*dx,b*dx) == "Left") {
                  proteinList[pNum]->numLeft[i_dat] += accessGlobals[pNum][c*Ny*Nz+a*Nz+b]*dV;
                }
                else if (triangle_section(a*dx,b*dx) == "Right") {
                  proteinList[pNum]->numRight[i_dat] += accessGlobals[pNum][c*Ny*Nz+a*Nz+b]*dV;
                }
                else {
                    proteinList[pNum]->numMid[i_dat] += accessGlobals[pNum][c*Ny*Nz+a*Nz+b]*dV;
                }
              }
            }
          }
        }
      }
    }

    //begin file printing
    if ((dump_flag == 1) && (i%printout_iterations == 0)) {
      dV = dx*dx*dx;
      fprintf(out_file,"Printing at iteration number = %d\n",i);
      int k = i/printout_iterations;

      //begin nATP printing.
      char *outfilenameATP = new char[1024];
      sprintf(outfilenameATP, "data/shape-%s/%s%snATP-%s-%03.2f-%03.2f-%03.2f-%03.2f-%03.2f-%03d.dat", argv[1],hires_flag_str,slice_flag_str,argv[1],A,B,C,D,density_factor,k);
      FILE *nATPfile = fopen((const char *)outfilenameATP,"w");
      delete[] outfilenameATP;

      if (slice_flag==1) {
        for (int a=0;a<Ny;a++){
          for (int b=0;b<Nz;b++){
            fprintf(nATPfile, "%1.2f ", nATP[(int(Nx/2))*Ny*Nz+a*Nz+b]);
          }
          fprintf(nATPfile, "\n");
        }
        fclose(nATPfile);
      }

      else {
        for (int a=0;a<Ny;a++){
          for (int b=0;b<Nz;b++){
            double nATPsum = 0;
            for (int c=0;c<Nx;c++){
              nATPsum += nATP[c*Ny*Nz+a*Nz+b];
            }
            fprintf(nATPfile, "%1.2f ", nATPsum);
          }
          fprintf(nATPfile, "\n");
        }
        fclose(nATPfile);
      }
      //end nATP printing

      //nE printing
      char *outfilenameE = new char[1000];
      sprintf(outfilenameE, "data/shape-%s/%s%snE-%s-%03.2f-%03.2f-%03.2f-%03.2f-%03.2f-%03d.dat", argv[1],hires_flag_str,slice_flag_str,argv[1],A,B,C,D,density_factor,k);
      FILE *nEfile = fopen((const char *)outfilenameE,"w");
      delete[] outfilenameE;

      if (slice_flag==1) {
        for (int a=0;a<Ny;a++){
          for (int b=0;b<Nz;b++){
            fprintf(nEfile, "%1.2f ", nE[(int(Nx/2))*Ny*Nz+a*Nz+b]);
          }
          fprintf(nEfile, "\n");
        }
        fclose(nEfile);
      }

      else {
        for (int a=0;a<Ny;a++){
          for (int b=0;b<Nz;b++){
            double nEsum = 0;
            for (int c=0;c<Nx;c++){
              nEsum += nE[c*Ny*Nz+a*Nz+b];
            }
            fprintf(nEfile, "%1.2f ", nEsum);
          }
          fprintf(nEfile, "\n");
        }
        fclose(nEfile);
      }
      //end nE printing


      //nADP printing
      char *outfilenameADP = new char[1000];
      sprintf(outfilenameADP, "data/shape-%s/%s%snADP-%s-%03.2f-%03.2f-%03.2f-%03.2f-%03.2f-%03d.dat", argv[1],hires_flag_str,slice_flag_str,argv[1],A,B,C,D,density_factor,k);
      FILE *nADPfile = fopen((const char *)outfilenameADP,"w");
      delete[] outfilenameADP;

      if (slice_flag==1) {
        for (int a=0;a<Ny;a++){
          for (int b=0;b<Nz;b++){
            fprintf(nADPfile, "%1.2f ", nADP[(int(Nx/2))*Ny*Nz+a*Nz+b]);
          }
          fprintf(nADPfile, "\n");
        }
        fclose(nADPfile);
      }

      else {
        for (int a=0;a<Ny;a++){
          for (int b=0;b<Nz;b++){
            double nADPsum = 0;
            for (int c=0;c<Nx;c++){
              nADPsum += nADP[c*Ny*Nz+a*Nz+b];
            }
            fprintf(nADPfile, "%1.2f ", nADPsum);
          }
          fprintf(nADPfile, "\n");
        }
        fclose(nADPfile);
      }
      //end nADP printing

      //begin ND printing
      char *outfilenameD = new char[1000];
      sprintf(outfilenameD, "data/shape-%s/%s%sND-%s-%03.2f-%03.2f-%03.2f-%03.2f-%03.2f-%03d.dat", argv[1],hires_flag_str,slice_flag_str,argv[1],A,B,C,D,density_factor,k);
      FILE *NDfile = fopen((const char *)outfilenameD,"w");
      delete[] outfilenameD;

      if (slice_flag==1) {
        for (int a=0;a<Ny;a++){
          for (int b=0;b<Nz;b++){
            fprintf(NDfile, "%1.2f ", ND[(int(Nx/2))*Ny*Nz+a*Nz+b]);
          }
          fprintf(NDfile, "\n");
        }
        fclose(NDfile);
      }

      else {
        for (int a=0;a<Ny;a++){
          for (int b=0;b<Nz;b++){
            double NDsum = 0;
            for (int c=0;c<Nx;c++){
              NDsum += ND[c*Ny*Nz+a*Nz+b];
            }
            fprintf(NDfile, "%1.2f ", NDsum);
          }
          fprintf(NDfile, "\n");
        }
        fclose(NDfile);
      }
      //end ND printing

      //begin NDE printing
      char *outfilenameDE = new char[1000];
      sprintf(outfilenameDE, "data/shape-%s/%s%sNDE-%s-%03.2f-%03.2f-%03.2f-%03.2f-%03.2f-%03d.dat", argv[1],hires_flag_str,slice_flag_str,argv[1],A,B,C,D,density_factor,k);
      FILE *NDEfile = fopen((const char *)outfilenameDE,"w");
      delete[] outfilenameDE;

      if (slice_flag==1) {
        for (int a=0;a<Ny;a++){
          for (int b=0;b<Nz;b++){
            fprintf(NDEfile, "%1.2f ", NDE[(int(Nx/2))*Ny*Nz+a*Nz+b]);
          }
          fprintf(NDEfile, "\n");
        }
        fclose(NDEfile);
      }

      else {
        for (int a=0;a<Ny;a++){
          for (int b=0;b<Nz;b++){
            double NDEsum = 0;
            for (int c=0;c<Nx;c++){
              NDEsum += NDE[c*Ny*Nz+a*Nz+b];
            }
            fprintf(NDEfile, "%1.2f ", NDEsum);
          }
          fprintf(NDEfile, "\n");
        }
        fclose(NDEfile);
      }
      //end NDE printing

      //begin NflE printing
      char *outfilenameflE = new char[1000];
      sprintf(outfilenameflE, "data/shape-%s/%s%sNflE-%s-%03.2f-%03.2f-%03.2f-%03.2f-%03.2f-%03d.dat", argv[1],hires_flag_str,slice_flag_str,argv[1],A,B,C,D,density_factor,k);
      FILE *NflEfile = fopen((const char *)outfilenameflE,"w");
      delete[] outfilenameflE;

      if (slice_flag==1) {
        for (int a=0;a<Ny;a++){
          for (int b=0;b<Nz;b++){
            fprintf(NflEfile, "%1.2f ", nE[(int(Nx/2))*Ny*Nz+a*Nz+b]*dV + NDE[(int(Nx/2))*Ny*Nz+a*Nz+b]);
          }
          fprintf(NflEfile, "\n");
        }
        fclose(NflEfile);
      }

      else {
        for (int a=0;a<Ny;a++){
          for (int b=0;b<Nz;b++){
            double NflEsum = 0;
            for (int c=0;c<Nx;c++){
              NflEsum += nE[c*Ny*Nz+a*Nz+b]*dV + NDE[c*Ny*Nz+a*Nz+b];
            }
            fprintf(NflEfile, "%1.2f ", NflEsum);
          }
          fprintf(NflEfile, "\n");
        }
        fclose(NflEfile);
      }
      //end NflE printing

      //begin NflD printing
      char *outfilenameflD = new char[1000];
      sprintf(outfilenameflD, "data/shape-%s/%s%sNflD-%s-%03.2f-%03.2f-%03.2f-%03.2f-%03.2f-%03d.dat", argv[1],hires_flag_str,slice_flag_str,argv[1],A,B,C,D,density_factor,k);
      FILE *NflDfile = fopen((const char *)outfilenameflD,"w");
      delete[] outfilenameflD;

      if (slice_flag==1) {
        for (int a=0;a<Ny;a++){
          for (int b=0;b<Nz;b++){
            fprintf(NflDfile, "%1.2f ", NDE[(int(Nx/2))*Ny*Nz+a*Nz+b] + nADP[(int(Nx/2))*Ny*Nz+a*Nz+b]*dV + nATP[(int(Nx/2))*Ny*Nz+a*Nz+b]*dV + ND[(int(Nx/2))*Ny*Nz+a*Nz+b]);
          }
          fprintf(NflDfile, "\n");
        }
        fclose(NflDfile);
      }

      else {
        for (int a=0;a<Ny;a++){
          for (int b=0;b<Nz;b++){
            double NflDsum = 0;
            for (int c=0;c<Nx;c++){
              NflDsum += NDE[c*Ny*Nz+a*Nz+b] + nADP[c*Ny*Nz+a*Nz+b]*dV + nATP[c*Ny*Nz+a*Nz+b]*dV + ND[c*Ny*Nz+a*Nz+b];
            }
            fprintf(NflDfile, "%1.2f ", NflDsum);
          }
          fprintf(NflDfile, "\n");
        }
        fclose(NflDfile);
      }
      //end NflD printing
      k++;
      fflush(out_file);
    }
  }
  //end file printing
  //end simulation

  fclose(out_file);

  //printing plot information
  for (int pNum=0; pNum<numProteins; pNum++) {

    //time map
    char *timename = new char[1024];
    sprintf(timename,"%s",print_filename("time-map",proteinList[pNum]->name));
    FILE* time_map = fopen(timename,"w");

    for (int a=0; a<Ny; a++) {
      for (int b=0; b<Nz; b++) {
        fprintf(time_map,"%1.2f\t",(proteinList[pNum]->sum[a*Nz+b])/((double)iter));
      }
      fprintf(time_map,"\n");
    }

    fclose(time_map);
    delete[] timename;
  }
  char *boxname = new char[1024];
  sprintf(boxname,"%s",print_filename("box-plot",""));
  FILE* box_plot = fopen(boxname,"w");

  for (int pNum=0; pNum<numProteins; pNum++) {

    if (mem_f_shape == "p" || rand_seed == 99 || rand_seed == 98 || mem_f_shape == "triangle") {
      fprintf(box_plot,"%s\tleft\t",proteinList[pNum]->name);
      for (int i_dat=0; i_dat<iter/print_denominator; i_dat++) {
        fprintf(box_plot,"%1.2f\t",(proteinList[pNum]->numLeft[i_dat]));
      }
      fprintf(box_plot,"\n");
      fprintf(box_plot,"%s\tmid\t",proteinList[pNum]->name);
      for (int i_dat=0; i_dat<iter/print_denominator; i_dat++) {
        fprintf(box_plot,"%1.2f\t",(proteinList[pNum]->numMid[i_dat]));
      }
      fprintf(box_plot,"\n");
      fprintf(box_plot,"%s\tright\t",proteinList[pNum]->name);
      for (int i_dat=0; i_dat<iter/print_denominator; i_dat++) {
        fprintf(box_plot,"%1.2f\t",(proteinList[pNum]->numRight[i_dat]));
      }
      fprintf(box_plot,"\n");
      fprintf(box_plot,"\n");
    }
    if (rand_seed == 97) {
      fprintf(box_plot,"%s\tleft\t",proteinList[pNum]->name);
      for (int i_dat=0; i_dat<iter/print_denominator; i_dat++) {
        fprintf(box_plot,"%1.2f\t",(proteinList[pNum]->numLeft[i_dat]));
      }
      fprintf(box_plot,"\n");
      fprintf(box_plot,"%s\trightup\t",proteinList[pNum]->name);
      for (int i_dat=0; i_dat<iter/print_denominator; i_dat++) {
        fprintf(box_plot,"%1.2f\t",(proteinList[pNum]->numRightUp[i_dat]));
      }
      fprintf(box_plot,"%s\tmid\t",proteinList[pNum]->name);
      fprintf(box_plot,"\n");
      for (int i_dat=0; i_dat<iter/print_denominator; i_dat++) {
        fprintf(box_plot,"%1.2f\t",(proteinList[pNum]->numMid[i_dat]));
      }
      fprintf(box_plot,"\n");
      fprintf(box_plot,"%s\trightdown\t",proteinList[pNum]->name);
      for (int i_dat=0; i_dat<iter/print_denominator; i_dat++) {
        fprintf(box_plot,"%1.2f\t",(proteinList[pNum]->numRightDown[i_dat]));
      }
      fprintf(box_plot,"\n");
      fprintf(box_plot,"\n");
    }
    if (rand_seed == 96) {
      fprintf(box_plot,"%s\trightup\t",proteinList[pNum]->name);
      for (int i_dat=0; i_dat<iter/print_denominator; i_dat++) {
        fprintf(box_plot,"%1.2f\t",(proteinList[pNum]->numRightUp[i_dat]));
      }
      fprintf(box_plot,"\n");
      fprintf(box_plot,"%s\tleftup\t",proteinList[pNum]->name);
      for (int i_dat=0; i_dat<iter/print_denominator; i_dat++) {
        fprintf(box_plot,"%1.2f\t",(proteinList[pNum]->numLeftUp[i_dat]));
      }
      fprintf(box_plot,"\n");
      fprintf(box_plot,"%s\tleftdown\t",proteinList[pNum]->name);
      for (int i_dat=0; i_dat<iter/print_denominator; i_dat++) {
        fprintf(box_plot,"%1.2f\t",(proteinList[pNum]->numLeftDown[i_dat]));
      }
      fprintf(box_plot,"\n");
      fprintf(box_plot,"%s\trightdown\t",proteinList[pNum]->name);
      for (int i_dat=0; i_dat<iter/print_denominator; i_dat++) {
        fprintf(box_plot,"%1.2f\t",(proteinList[pNum]->numRightDown[i_dat]));
      }
      fprintf(box_plot,"\n");
      fprintf(box_plot,"\n");
    }
  }

  fclose(box_plot);
  delete[] boxname;

  char *avename = new char[1024];
  sprintf(avename,"%s",print_filename("ave_plot",""));
  FILE* ave_plot = fopen(avename,"w");

  double left_area_total = 0;
  double middle_area_total = 0;
  double right_area_total = 0;
  for (int a=0;a<Ny; a++) {
    for (int b=0; b<Nz; b++) {
      for (int c=0; c<Nx; c++) {
        if (b < box_divider_left) {
          left_area_total += mem_A[c*Ny*Nz+a*Nz+b];
        }
        else if (b > box_divider_right) {
          right_area_total += mem_A[c*Ny*Nz+a*Nz+b];
        }
        else {
          middle_area_total += mem_A[c*Ny*Nz+a*Nz+b];
        }
      }
    }
  }

  if (mem_f_shape == "p"){
    for (int pNum=3; pNum<numProteins; pNum++) {

      fprintf(ave_plot,"%s\tleft\t",proteinList[pNum]->name);
      for (int i_dat=0; i_dat<iter/print_denominator; i_dat++) {
        fprintf(ave_plot,"%1.2f\t",(proteinList[pNum]->numLeft[i_dat]/left_area_total));
      }
      fprintf(ave_plot,"\n");
      fprintf(ave_plot,"%s\tmid\t",proteinList[pNum]->name);
      for (int i_dat=0; i_dat<iter/print_denominator; i_dat++) {
        fprintf(ave_plot,"%1.2f\t",(proteinList[pNum]->numMid[i_dat]/middle_area_total));
      }
      fprintf(ave_plot,"\n");
      fprintf(ave_plot,"%s\tright\t",proteinList[pNum]->name);
      for (int i_dat=0; i_dat<iter/print_denominator; i_dat++) {
        fprintf(ave_plot,"%1.2f\t",(proteinList[pNum]->numRight[i_dat]/right_area_total));
      }
      fprintf(ave_plot,"\n");
      fprintf(ave_plot,"\n");
    }
  }
  fclose(ave_plot);
  delete[] avename;

  //printing to the project directory so we have a shortlist of what we've done.
  char *fname = new char[1024];
  sprintf(fname,"catalog.txt");
  FILE * catalog;
  int catalog_exists;
  catalog = fopen(fname,"r");
  if (catalog==NULL) {
    catalog_exists=0;
  }
  else {
    catalog_exists=1;
    fclose(catalog);
  }
  if (catalog_exists==1) {
    catalog=fopen(fname,"a+b");
  }
  else {
    catalog=fopen(fname,"w+b");
  }
  if (catalog!=NULL) {
    fprintf(catalog,"%s %1.2f %1.2f %1.2f %1.2f %1.2f", mem_f_shape.c_str(),A,B,C,D,density_factor);
    if (dx==.05) {
      fprintf(catalog," -hires\n");
    }
    else {
      fprintf(catalog,"\n");
    }
    fclose(catalog);
  }
  delete[] fname;
  //end catalog

  return 0;
}


void set_membrane(FILE * out_file, double (*mem_f)(double x, double y, double z),
                  double mem_A[]) {
  clock_t old_time = clock();
  for(int xi=0;xi<Nx;xi++){
    clock_t time = clock();
    fprintf(out_file, "x row %d in set_membrane took %4.02f seconds",xi, (time-old_time)/double(CLOCKS_PER_SEC));
    fflush(stdout);
    old_time = time;
    for(int yi=0;yi<Ny;yi++){
      for(int zi=0;zi<Nz;zi++){
        double fXYZ = mem_f((xi+0.5)*dx, (yi+0.5)*dx, (zi+0.5)*dx);
        double fXYz = mem_f((xi+0.5)*dx, (yi+0.5)*dx, (zi-0.5)*dx);
        double fXyZ = mem_f((xi+0.5)*dx, (yi-0.5)*dx, (zi+0.5)*dx);
        double fXyz = mem_f((xi+0.5)*dx, (yi-0.5)*dx, (zi-0.5)*dx);
        double fxYZ = mem_f((xi-0.5)*dx, (yi+0.5)*dx, (zi+0.5)*dx);
        double fxYz = mem_f((xi-0.5)*dx, (yi+0.5)*dx, (zi-0.5)*dx);
        double fxyZ = mem_f((xi-0.5)*dx, (yi-0.5)*dx, (zi+0.5)*dx);
        double fxyz = mem_f((xi-0.5)*dx, (yi-0.5)*dx, (zi-0.5)*dx);
        double f = mem_f(xi*dx, yi*dx, zi*dx);
        mem_A[xi*Ny*Nz+yi*Nz+zi] = find_intersection(fXYZ, fXYz, fXyZ, fXyz, fxYZ, fxYz, fxyZ, fxyz, f);
      }
    }
  }
}

void set_curvature(double mem_A[], double curvature[]){
  printf("doing set curvature!!!\n");
  fflush(stdout);
  double X = Nx*dx;
  double x1 = (X-A)/2.0;
  double x2 = (X+A)/2.0;
  for(int xi=0;xi<Nx;xi++){
    for(int yi=0;yi<Ny;yi++){
      for(int zi=0;zi<Nz;zi++){
        if (mem_A[xi*Ny*Nz+yi*Nz+zi]==0 || (xi*dx+0.05)>x2 || (xi*dx-0.05)<x1 ){
          curvature[xi*Ny*Nz+yi*Nz+zi]=0;
        } else {
          double fX = mem_f((xi+0.5)*dx, yi*dx, zi*dx);
          double fx = mem_f((xi-0.5)*dx, yi*dx, zi*dx);
          double fY = mem_f(xi*dx, (yi+0.5)*dx, zi*dx);
          double fy = mem_f(xi*dx, (yi-0.5)*dx, zi*dx);
          double fZ = mem_f(xi*dx, yi*dx, (zi+0.5)*dx);
          double fz = mem_f(xi*dx, yi*dx, (zi-0.5)*dx);
          double f = mem_f(xi*dx, yi*dx, zi*dx);

          double df_dx = (fX-fx)/dx;
          double df_dy = (fY-fy)/dx;
          double df_dz = (fZ-fz)/dx;
          double constant = sqrt((df_dx)*(df_dx) + (df_dy)*(df_dy) + (df_dz)*(df_dz));

          curvature[xi*Ny*Nz+yi*Nz+zi] = 4*(fX + fx + fY + fy + fZ + fz - 6*f)/dx/dx/constant;
        }
      }
    }
  }
  fflush(stdout);
}

double find_intersection(const double fXYZ, const double fXYz, const double fXyZ, const double fXyz,
                         const double fxYZ, const double fxYz, const double fxyZ, const double fxyz,
                         const double f) {
  double dA = 0;
  double *ptsx = new double[8];
  double *ptsy = new double[8];
  double *ptsz = new double[8];
  for (int i=0;i<8;i++) ptsx[i] = 0;
  for (int i=0;i<8;i++) ptsy[i] = 0;
  for (int i=0;i<8;i++) ptsz[i] = 0;
  int np = 0;
  double df_dx = (fXYZ + fXYz + fXyZ + fXyz - fxYZ - fxyZ - fxYz - fxyz)/(4*dx);
  double df_dy = (fXYZ + fXYz + fxYZ + fxYz - fXyZ - fxyZ - fXyz - fxyz)/(4*dx);
  double df_dz = (fXYZ + fXyZ + fxYZ + fxyZ - fXYz - fXyz - fxYz - fxyz)/(4*dx);
  for (double j=-0.5; j<1.0; j++){
    for (double k=-0.5; k<1.0; k++){
      if ((-f - j*dx*df_dy - k*dx*df_dz)/(df_dx*dx) < 0.5 && (-f - j*dx*df_dy - k*dx*df_dz)/(df_dx*dx) > -0.5){
        ptsx[np] = (-f - j*dx*df_dy - k*dx*df_dz)/df_dx;
        ptsy[np] = j*dx;
        ptsz[np] = k*dx;
        np++;
      }
    }
  }
  for (double i=-0.5; i<1.0; i++){
    for (double k=-0.5; k<1.0; k++){
      if ((-f - i*dx*df_dx - k*dx*df_dz)/(df_dy*dx) < 0.5 && (-f - i*dx*df_dx - k*dx*df_dz)/(df_dy*dx) > -0.5){
        ptsy[np] = (-f - i*dx*df_dx - k*dx*df_dz)/df_dy;
        ptsx[np] = i*dx;
        ptsz[np] = k*dx;
        np++;
      }
    }
  }
  for (double j=-0.5; j<1.0; j++){
    for (double i=-0.5; i<1.0; i++){
      if ((-f - j*dx*df_dy - i*dx*df_dx)/(df_dz*dx) < 0.5 && (-f - j*dx*df_dy - i*dx*df_dx)/(df_dz*dx) > -0.5){
        ptsz[np] = (-f - j*dx*df_dy - i*dx*df_dx)/df_dz;
        ptsy[np] = j*dx;
        ptsx[np] = i*dx;
        np++;
      }
    }
  }
  int nz0 = 0;
  int ny0 = 0;
  int nx0 = 0;
  double *ptx = new double[8];
  double *pty = new double[8];
  double *ptz = new double[8];
  for (int i=0;i<8;i++) ptx[i] = 0;
  for (int i=0;i<8;i++) pty[i] = 0;
  for (int i=0;i<8;i++) ptz[i] = 0;
  double *line = new double[8];
  for (int i=0;i<8;i++) line[i] = 0;
  double *cos = new double[8];
  double cos_max = -2;
  int as=0;
  int bs=0;
  int cs=0;
  int ds=0;
  int s = 1;
  double eline = 0;
  double p = 0;
  for (int i=0;i<np;i++){
    if (ptsz[i] == -0.5*dx) {ptz[nz0]=ptsz[i]; ptx[nz0]=ptsx[i]; pty[nz0]=ptsy[i]; nz0++;}
  }
  if (nz0 == 2){
    line[0] = sqrt((ptx[1]-ptx[0])*(ptx[1]-ptx[0]) + (pty[1]-pty[0])*(pty[1]-pty[0])
                   + (ptz[1]-ptz[0])*(ptz[1]-ptz[0]));
    for (int i=0;i<np;i++){
      if (ptsz[i] != -0.5*dx){
        line[s] = sqrt((ptsx[i]-ptx[0])*(ptsx[i]-ptx[0]) + (ptsy[i]-pty[0])*(ptsy[i]-pty[0])
                       + (ptsz[i]-ptz[0])*(ptsz[i]-ptz[0]));
        cos[s] = ((ptsx[i]-ptx[0])*(ptx[1]-ptx[0]) + (ptsy[i]-pty[0])*(pty[1]-pty[0])
                  + (ptsz[i]-ptz[0])*(ptz[1]-ptz[0])) / (line[s]*line[0]);
        ptz[s+1] = ptsz[i]; ptx[s+1] = ptsx[i]; pty[s+1] = ptsy[i];
        s++;
      }
    }
  } else {
    for (int i=0;i<np;i++){
      if (ptsx[i] == -0.5*dx) {ptz[nx0]=ptsz[i]; ptx[nx0]=ptsx[i]; pty[nx0]=ptsy[i]; nx0++;}
    }
    if (nx0 == 2){
      line[0] = sqrt((ptx[1]-ptx[0])*(ptx[1]-ptx[0]) + (pty[1]-pty[0])*(pty[1]-pty[0])
                     + (ptz[1]-ptz[0])*(ptz[1]-ptz[0]));
      for (int i=0;i<np;i++){
        if (ptsx[i] != -0.5*dx){
          line[s] = sqrt((ptsx[i]-ptx[0])*(ptsx[i]-ptx[0]) + (ptsy[i]-pty[0])*(ptsy[i]-pty[0])
                         + (ptsz[i]-ptz[0])*(ptsz[i]-ptz[0]));
          cos[s] = ((ptsx[i]-ptx[0])*(ptx[1]-ptx[0]) + (ptsy[i]-pty[0])*(pty[1]-pty[0])
                    + (ptsz[i]-ptz[0])*(ptz[1]-ptz[0])) / (line[s]*line[0]);
          ptz[s+1] = ptsz[i]; ptx[s+1] = ptsx[i]; pty[s+1] = ptsy[i];
          s++;
        }
      }
    } else {
      for (int i=0;i<np;i++){
        if (ptsy[i] == -0.5*dx) {ptz[ny0]=ptsz[i]; ptx[ny0]=ptsx[i]; pty[ny0]=ptsy[i]; ny0++;}
      }
      if (ny0 == 2){
        line[0] = sqrt((ptx[1]-ptx[0])*(ptx[1]-ptx[0]) + (pty[1]-pty[0])*(pty[1]-pty[0])
                       + (ptz[1]-ptz[0])*(ptz[1]-ptz[0]));
        for (int i=0;i<np;i++){
          if (ptsy[i] != -0.5*dx){
            line[s] = sqrt((ptsx[i]-ptx[0])*(ptsx[i]-ptx[0]) + (ptsy[i]-pty[0])*(ptsy[i]-pty[0])
                           + (ptsz[i]-ptz[0])*(ptsz[i]-ptz[0]));
            cos[s] = ((ptsx[i]-ptx[0])*(ptx[1]-ptx[0]) + (ptsy[i]-pty[0])*(pty[1]-pty[0])
                      + (ptsz[i]-ptz[0])*(ptz[1]-ptz[0])) / (line[s]*line[0]);
            ptz[s+1] = ptsz[i]; ptx[s+1] = ptsx[i]; pty[s+1] = ptsy[i];
            s++;
          }
        }
      }
    }
  }
  for (int i=1;i<s;i++){
    if (cos[i] > cos_max){cos_max = cos[i]; as=i;}
  }
  eline = sqrt((ptx[as+1]-ptx[1])*(ptx[as+1]-ptx[1]) + (pty[as+1]-pty[1])*(pty[as+1]-pty[1])
               + (ptz[as+1]-ptz[1])*(ptz[as+1]-ptz[1]));
  p = (line[as]+line[0]+eline)/2;
  dA = sqrt(p*(p-line[as])*(p-line[0])*(p-eline));
  cos[as] = -2; cos_max = -2;

  if (np==4 || np==5 || np==6){
    for (int i=1;i<s;i++){
      if (cos[i] > cos_max){cos_max = cos[i]; bs=i;}
    }
    eline = sqrt((ptx[bs+1]-ptx[as+1])*(ptx[bs+1]-ptx[as+1]) + (pty[bs+1]-pty[as+1])*(pty[bs+1]-pty[as+1])
                 + (ptz[bs+1]-ptz[as+1])*(ptz[bs+1]-ptz[as+1]));
    p = (line[bs]+line[as]+eline)/2;
    dA += sqrt(p*(p-line[bs])*(p-line[as])*(p-eline));
    cos[bs] = -2; cos_max = -2;
  }

  if (np==5 || np==6){
    for (int i=1;i<s;i++){
      if (cos[i] > cos_max){cos_max = cos[i]; cs=i;}
    }
    eline = sqrt((ptx[cs+1]-ptx[bs+1])*(ptx[cs+1]-ptx[bs+1]) + (pty[cs+1]-pty[bs+1])*(pty[cs+1]-pty[bs+1])
                 + (ptz[cs+1]-ptz[bs+1])*(ptz[cs+1]-ptz[bs+1]));
    p = (line[cs]+line[bs]+eline)/2;
    dA += sqrt(p*(p-line[cs])*(p-line[bs])*(p-eline));
    cos[cs] = -2; cos_max = -2;
  }

  if (np == 6) {
    for (int i=1;i<s;i++){
      if (cos[i] > cos_max){cos_max = cos[i]; ds=i;}
    }
    eline = sqrt((ptx[ds+1]-ptx[cs+1])*(ptx[ds+1]-ptx[cs+1]) + (pty[ds+1]-pty[cs+1])*(pty[ds+1]-pty[cs+1])
                 + (ptz[ds+1]-ptz[cs+1])*(ptz[ds+1]-ptz[cs+1]));
    p = (line[ds]+line[cs]+eline)/2;
    dA += sqrt(p*(p-line[ds])*(p-line[cs])*(p-eline));
  }
  delete[] ptsx;
  delete[] ptsy;
  delete[] ptsz;
  delete[] ptx;
  delete[] pty;
  delete[] ptz;
  delete[] line;
  delete[] cos;
  return dA;
}

void set_insideArr(bool *insideArr){
  for(int xi=0;xi<Nx;xi++){
    for(int yi=0;yi<Ny;yi++){
      for(int zi=0;zi<Nz;zi++){
        insideArr[xi*Ny*Nz+yi*Nz+zi] = inside(xi,yi,zi);
      }
    }
  }
}

bool inside (int xi, int yi, int zi){
  if (mem_f((xi-0.5)*dx,(yi-0.5)*dx,(zi-0.5)*dx) <= 0) {return true;}
  if (mem_f((xi+0.5)*dx,(yi-0.5)*dx,(zi-0.5)*dx) <= 0) {return true;}
  if (mem_f((xi-0.5)*dx,(yi+0.5)*dx,(zi-0.5)*dx) <= 0) {return true;}
  if (mem_f((xi-0.5)*dx,(yi-0.5)*dx,(zi+0.5)*dx) <= 0) {return true;}
  if (mem_f((xi-0.5)*dx,(yi+0.5)*dx,(zi+0.5)*dx) <= 0) {return true;}
  if (mem_f((xi+0.5)*dx,(yi-0.5)*dx,(zi+0.5)*dx) <= 0) {return true;}
  if (mem_f((xi+0.5)*dx,(yi+0.5)*dx,(zi-0.5)*dx) <= 0) {return true;}
  if (mem_f((xi+0.5)*dx,(yi+0.5)*dx,(zi+0.5)*dx) <= 0) {return true;}
  return false;
}

int get_J(double difD, double *nATP, double *nADP, double *nE,
          double *JxATP, double *JyATP, double *JzATP,
          double *JxADP, double *JyADP, double *JzADP,
          double *JxE, double *JyE, double *JzE){
  for(int xi=0;xi<Nx-1;xi++){
    for(int yi=0;yi<Ny-1;yi++){
      for(int zi=0;zi<Nz-1;zi++){
        JzATP[xi*Ny*Nz+yi*Nz+zi] = -difD*(nATP[xi*Ny*Nz+yi*Nz+zi+1]-nATP[xi*Ny*Nz+yi*Nz+zi])/dx;
        JyATP[xi*Ny*Nz+yi*Nz+zi] = -difD*(nATP[xi*Ny*Nz+(yi+1)*Nz+zi]-nATP[xi*Ny*Nz+yi*Nz+zi])/dx;
        JxATP[xi*Ny*Nz+yi*Nz+zi] = -difD*(nATP[(xi+1)*Ny*Nz+yi*Nz+zi]-nATP[xi*Ny*Nz+yi*Nz+zi])/dx;
        JzADP[xi*Ny*Nz+yi*Nz+zi] = -difD*(nADP[xi*Ny*Nz+yi*Nz+zi+1]-nADP[xi*Ny*Nz+yi*Nz+zi])/dx;
        JyADP[xi*Ny*Nz+yi*Nz+zi] = -difD*(nADP[xi*Ny*Nz+(yi+1)*Nz+zi]-nADP[xi*Ny*Nz+yi*Nz+zi])/dx;
        JxADP[xi*Ny*Nz+yi*Nz+zi] = -difD*(nADP[(xi+1)*Ny*Nz+yi*Nz+zi]-nADP[xi*Ny*Nz+yi*Nz+zi])/dx;
        JzE[xi*Ny*Nz+yi*Nz+zi] = -difD*(nE[xi*Ny*Nz+yi*Nz+zi+1]-nE[xi*Ny*Nz+yi*Nz+zi])/dx;
        JyE[xi*Ny*Nz+yi*Nz+zi] = -difD*(nE[xi*Ny*Nz+(yi+1)*Nz+zi]-nE[xi*Ny*Nz+yi*Nz+zi])/dx;
        JxE[xi*Ny*Nz+yi*Nz+zi] = -difD*(nE[(xi+1)*Ny*Nz+yi*Nz+zi]-nE[xi*Ny*Nz+yi*Nz+zi])/dx;
      }
    }
  }
  return 0;
}

int get_next_density(double *mem_A, bool *insideArr, double *nATP, double *nADP,
                     double *nE, double *Nd, double *Nde,
                     double *JxATP, double *JyATP, double *JzATP,
                     double *JxADP, double *JyADP, double *JzADP,
                     double *JxE, double *JyE, double *JzE){
  for(int xi=0;xi<Nx-1;xi++){
    for(int yi=0;yi<Ny-1;yi++){
      for(int zi=0;zi<Nz-1;zi++){
        //for the diffusion terms, we use dn/dt = dA*J/dV thinking in terms of tot #, but dA/dV=1/dx
        if (insideArr[(xi+1)*Ny*Nz+yi*Nz+zi] && insideArr[xi*Ny*Nz+yi*Nz+zi]){
          nADP[(xi+1)*Ny*Nz+yi*Nz+zi] += JxADP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
          nADP[xi*Ny*Nz+yi*Nz+zi] -= JxADP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
          nATP[(xi+1)*Ny*Nz+yi*Nz+zi] += JxATP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
          nATP[xi*Ny*Nz+yi*Nz+zi] -= JxATP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
          nE[(xi+1)*Ny*Nz+yi*Nz+zi] += JxE[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
          nE[xi*Ny*Nz+yi*Nz+zi] -= JxE[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
        }
        if (insideArr[xi*Ny*Nz+(yi+1)*Nz+zi] && insideArr[xi*Ny*Nz+yi*Nz+zi]){
          nADP[xi*Ny*Nz+(yi+1)*Nz+zi] += JyADP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
          nADP[xi*Ny*Nz+yi*Nz+zi] -= JyADP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
          nATP[xi*Ny*Nz+(yi+1)*Nz+zi] += JyATP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
          nATP[xi*Ny*Nz+yi*Nz+zi] -= JyATP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
          nE[xi*Ny*Nz+(yi+1)*Nz+zi] += JyE[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
          nE[xi*Ny*Nz+yi*Nz+zi] -= JyE[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
        }
        if (insideArr[xi*Ny*Nz+yi*Nz+(zi+1)] && insideArr[xi*Ny*Nz+yi*Nz+zi]){
          nADP[xi*Ny*Nz+yi*Nz+(zi+1)] += JzADP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
          nADP[xi*Ny*Nz+yi*Nz+zi] -= JzADP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
          nATP[xi*Ny*Nz+yi*Nz+(zi+1)] += JzATP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
          nATP[xi*Ny*Nz+yi*Nz+zi] -= JzATP[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
          nE[xi*Ny*Nz+yi*Nz+(zi+1)] += JzE[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
          nE[xi*Ny*Nz+yi*Nz+zi] -= JzE[xi*Ny*Nz+yi*Nz+zi]/dx*time_step;
        }
      }
    }
  }
  double ADP_to_ATP;
  double de_to_ADP_E;
  double ATP_to_d;
  double E_d_to_de;
  for(int xi=0;xi<Nx;xi++){
    for(int yi=0;yi<Ny;yi++){
      for(int zi=0;zi<Nz;zi++){
        ADP_to_ATP = rate_ADP_ATP*nADP[xi*Ny*Nz+yi*Nz+zi]*time_step;
        de_to_ADP_E = rate_de*NDE[xi*Ny*Nz+yi*Nz+zi]/mem_A[xi*Ny*Nz+yi*Nz+zi]*time_step;
        ATP_to_d = (rate_D + rate_dD*(ND[xi*Ny*Nz+yi*Nz+zi] + NDE[xi*Ny*Nz+yi*Nz+zi])/mem_A[xi*Ny*Nz+yi*Nz+zi])
          *nATP[xi*Ny*Nz+yi*Nz+zi]*time_step;
        E_d_to_de = rate_E*ND[xi*Ny*Nz+yi*Nz+zi]/mem_A[xi*Ny*Nz+yi*Nz+zi]*nE[xi*Ny*Nz+yi*Nz+zi]*time_step;
        //Jeff!  remember that when you gain cyto density and lose the same amount of wall density,
        //the numbers of proteins gained/lost will be different, and it's the numbers that you want to be the same!!
        //also, keep thinking about the issue below where all additions are divided and then mult by the same mem_A fun
        nADP[xi*Ny*Nz+yi*Nz+zi] -= ADP_to_ATP;
        nATP[xi*Ny*Nz+yi*Nz+zi] += ADP_to_ATP;
        if (mem_A[xi*Ny*Nz+yi*Nz+zi] != 0){
          NDE[xi*Ny*Nz+yi*Nz+zi] += -de_to_ADP_E*mem_A[xi*Ny*Nz+yi*Nz+zi];
          nADP[xi*Ny*Nz+yi*Nz+zi] += de_to_ADP_E*mem_A[xi*Ny*Nz+yi*Nz+zi]/(dx*dx*dx);
          nE[xi*Ny*Nz+yi*Nz+zi] += de_to_ADP_E*mem_A[xi*Ny*Nz+yi*Nz+zi]/(dx*dx*dx);

          nATP[xi*Ny*Nz+yi*Nz+zi] -= ATP_to_d*mem_A[xi*Ny*Nz+yi*Nz+zi]/(dx*dx*dx);
          ND[xi*Ny*Nz+yi*Nz+zi] += ATP_to_d*mem_A[xi*Ny*Nz+yi*Nz+zi];

          nE[xi*Ny*Nz+yi*Nz+zi] -= E_d_to_de*mem_A[xi*Ny*Nz+yi*Nz+zi]/(dx*dx*dx);
          ND[xi*Ny*Nz+yi*Nz+zi] -= E_d_to_de*mem_A[xi*Ny*Nz+yi*Nz+zi];
          NDE[xi*Ny*Nz+yi*Nz+zi] += E_d_to_de*mem_A[xi*Ny*Nz+yi*Nz+zi];
        }
      }
    }
  }
  return 0;
}

int set_density(double *nATP, double *nE, double *mem_A){
  int right_most_point_z=0; //left and right most points for z
  int left_most_point_z=Nz;
  int right_most_point_y=0; //"left" and "right" most points for y, in terms of magnitude (right = larger y value)
  int left_most_point_y=Ny;
  for (int i=0;i<Nx;i++){
    for (int j=0;j<Ny;j++){
      for (int k=0;k<Nz;k++){
        if (inside(i,j,k)){
          if (k>right_most_point_z){
            right_most_point_z = k;
          }
          if(k<left_most_point_z){
            left_most_point_z = k;
          }
          if (j>right_most_point_y){
            right_most_point_y = j;
          }
          if (j<left_most_point_y){
            left_most_point_y = j;
          }
        }
      }
    }
  }

  int density_divider_z = int(right_most_point_z - (right_most_point_z - left_most_point_z)/3);

  //get total gridpoints, gridpoints left of divide, gridpoints right of divide for protein count
  int gridpoints_left = 0;
  int gridpoints_right = 0;
  int gridpoints_total = 0;
  for (int i=0;i<Nx;i++){
    for (int j=0;j<Ny;j++){
      for (int k=0;k<Nz;k++){
        if (inside(i,j,k)){
          gridpoints_total++;
          if (k<=density_divider_z) {
            gridpoints_left++;
          }
          else {
            gridpoints_right++;
          }
        }
      }
    }
  }

  //compute density scale factors left and right of divide (to ensure correct protein #)
  double density_factor_left = gridpoints_total/(gridpoints_left + density_factor*gridpoints_right);
  double density_factor_right = density_factor*gridpoints_total/(gridpoints_left + density_factor*gridpoints_right);

  printf("Density factors: left: %f, right: %f, ratio: %f\n", density_factor_left, density_factor_right, density_factor_right/density_factor_left);

  printf("Gridpoints left of the divider: %d\n",gridpoints_left);
  printf("Gridpoints right of the divider:%d\n",gridpoints_right);
  printf("Gridpoints total: %d\n",gridpoints_total);

  //begin setting density at each gridpoint:
  for (int i=0;i<Nx;i++){
    for (int j=0;j<Ny;j++){
      for (int k=0;k<Nz;k++){
        if (inside(i,j,k)){
          if(k>density_divider_z){
            nATP[i*Ny*Nz+j*Nz+k] = nATP_starting_density*density_factor_right;
            nE[i*Ny*Nz+j*Nz+k] = nE_starting_density*density_factor_right;
            nADP[i*Ny*Nz+j*Nz+k] =0;
            ND[i*Ny*Nz+j*Nz+k] =0;
            NDE[i*Ny*Nz+j*Nz+k] = 0;
          }
          else {
            nATP[i*Ny*Nz+j*Nz+k] = nATP_starting_density*density_factor_left;
            nE[i*Ny*Nz+j*Nz+k] = nE_starting_density*density_factor_left;
            nADP[i*Ny*Nz+j*Nz+k] =0;
            ND[i*Ny*Nz+j*Nz+k] =0;
            NDE[i*Ny*Nz+j*Nz+k] = 0;
          }
        }
        else {
          nATP[i*Ny*Nz+j*Nz+k] = 0;
          nE[i*Ny*Nz+j*Nz+k] = 0;
          nADP[i*Ny*Nz+j*Nz+k] =0;
          ND[i*Ny*Nz+j*Nz+k] =0;
          NDE[i*Ny*Nz+j*Nz+k] = 0;
        }
      }
    }
  }
  return 0;
}
