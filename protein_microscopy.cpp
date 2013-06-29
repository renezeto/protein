#include <iostream>
using namespace std;
#include <stdio.h>
#include <time.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "test.h"
#include "protein.h"
#include "MersenneTwister.h"
#include <cassert>
//  #include <ctime>

//notes:
//variables commented out and flagged as unused were causing compiler errors
//since they were unused.

const double difD = 2.5; // (um)^2 s^- 1
const double difE = 2.5; // (um)^2 s^-1
const double rate_ADP_ATP = 1; // s^-1
const double rate_D = .025; // um s^-1
const double rate_dD = .0015; // (um)^3 s^-1 
const double rate_de = .7; // s^-1
const double rate_E = .093; // (um)^3 s^-1

const double nATP_starting_density = 1000.0; //proteins per micrometer
const double nE_starting_density = 350.0; //proteins per micrometer
double density_factor;

const int n = 706; //what is this?
int area_rating_flag;

const double dx=0.15;
const double tot_time = 186;
const double time_step = .1*dx*dx/difD;
const int iter = int(tot_time/time_step)+3;
const int iter_at_half_sec = int(0.5/time_step)+1;

double x, y, z;
int Nx, Ny, Nz;

string mem_f_shape;
double A;
double B;
double C;
double D;

double *nATP;
double *nADP;
double *nE;
double *Nd;
double *Nde;
double *f_mem;

const int starting_num_guassians=20;
double guass[3*starting_num_guassians]; //stores y,z, and sigma for each guassian when creating random cell wall
const int random_num_guassians=5;
int rand_seed=0; //=14; at this point I have this passed in from the command line as the D argument
double Norm = 15.0; //This is the height of the guassians that make the cell wall

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

double f_2D_triangle(double y, double z){
  double f = 0; double Y = Ny*dx; double Z = Nz*dx;
  double y1 = A + 2*dx; double z1 = A + 2*dx;
  double y2 = Y-A-2*dx; double z2 = z1;
  double y3 = Y/2.0; double z3 = Z-A-2*dx;
  double b1 = (z3-z1)/(y3-y1); double a1 = z1 - b1*y1;
  double b3 = (z2-z3)/(y2-y3); double a3 = z3 - b3*y3;
  double zl1 = b1*y +a1;
  double zl3 = b3*y +a3;
  double rad = 1.75*(y2-y1)*sqrt(3.0)/6.0;
  double y_circle = Y/2.0; double z_circle = z1 + sqrt(3)*(y2-y1)/6.0;
  if ((z < zl1) && (z < zl3) && (z > z1) && ((z-z_circle)*(z-z_circle) + (y-y_circle)*(y-y_circle)) < rad*rad){
    f = -0.1;
  } else {
    f = 0.1;
  }
  return f;
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
      if (f>0) {
        double closest_y0 = -100.0; 
        double closest_z0 = -100.0; 
        //bool there_is_closest_point=0; unused
        for (double y0 = y-A; y0<y+A; y0+=dx) {
          for (double z0 = z-A; z0<z+A; z0+=dx) {
            if(mem_f_shape=="randst") f0 = f_2D_randst(y0,z0);
            if(mem_f_shape=="TIE_fighter") f0 = f_2D_TIE_fighter(y0,z0);
            if(mem_f_shape=="triangle") f0 = f_2D_triangle(y0,z0);
            if (f0 <= 0) {
              //printf("f0 = %g\n",f0);
              fflush(stdout);
              //there_is_closest_point = 1; unused
              if ( (y-y0)*(y-y0)+(z-z0)*(z-z0) < (y-closest_y0)*(y-closest_y0)+(z-closest_z0)*(z-closest_z0) ) {
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

int main (int argc, char *argv[]) {
  mem_f_shape = argv[1];
  A = atof(argv[2]);
  B = atof(argv[3]);
  C = atof(argv[4]);
  D = atof(argv[5]);
  density_factor = atof(argv[6]);

  // input extra command line argument equal to 1 for area_rating check only
  if(argc > 7) { 
    area_rating_flag = atoi(argv[7]);
  }

  //compute grid size based on cell parameters
  if (mem_f_shape=="p") {
    Nx = ceil(2*B/dx) + 4;
    Ny = ceil(2*B/dx) + 4;
    Nz = ceil((A + 2*B)/dx) + 4;
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
    Ny = ceil(B/dx) + 4;
    Nz = ceil(C/dx) + 4;
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
  char * out_file_name = new char[1024];
  sprintf(out_file_name,"data/shape-%s/out_files/%s-%4.02f-%4.02f-%4.02f-%4.02f-%4.02f.out",mem_f_shape.c_str(),mem_f_shape.c_str(),A,B,C,D,density_factor);
  FILE * out_file = fopen((const char *)out_file_name,"w");
  
  time_t t = time(0);   // get time now
  struct tm * now = localtime( & t );
  char * timevar = new char[1024];
  sprintf(timevar, "%d/%d/%d at %d hours and %d minutes", now->tm_mon +1,now->tm_mday,now->tm_year +1900,now->tm_hour,now->tm_min);
  fprintf(out_file,"This simulation was run on %s\n",timevar);


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
  nATP = new double[Nx*Ny*Nz];
  nADP = new double[Nx*Ny*Nz];
  nE = new double[Nx*Ny*Nz];
  Nd = new double[Nx*Ny*Nz];
  Nde = new double[Nx*Ny*Nz];
  f_mem = new double[Nx*Ny*Nz];
  fprintf(out_file,"For this simulation,\ndx = %f\ntot_time = %f\ntimestep = %f\ntotal iterations = %d\niter at five sec = %d\n",
         dx, tot_time, time_step, iter, iter_at_half_sec);
  double *JxATP = new double[Nx*Ny*Nz];
  double *JyATP = new double[Nx*Ny*Nz];
  double *JzATP = new double[Nx*Ny*Nz];
  double *JxADP = new double[Nx*Ny*Nz];
  double *JyADP = new double[Nx*Ny*Nz];
  double *JzADP = new double[Nx*Ny*Nz];
  double *JxE = new double[Nx*Ny*Nz];
  double *JyE = new double[Nx*Ny*Nz];
  double *JzE = new double[Nx*Ny*Nz];
  double *mem_A = new double[Nx*Ny*Nz]; //area of membrane in each cube
  bool *insideArr = new bool[Nx*Ny*Nz]; //whether each cube is inside at all
  for (int i=0;i<Nx*Ny*Nz;i++){mem_A[i] = 0;}
  bool force_to_generate_new_memA = true;
  if (mem_f_shape=="randst") {
    char* memA_name = new char[1024];
    sprintf(memA_name,"data/shape-randst/membrane_files/memA-%4.02f-%4.02f-%4.02f-%4.02f-%4.02f.dat",A,B,C,D,density_factor);
    FILE *memAin = fopen(memA_name,"r");
    if (!memAin || force_to_generate_new_memA) {
      if (memAin && force_to_generate_new_memA) fclose(memAin);
      fprintf(out_file,"There is evidently no file called %s,\n so we're going to create are own and fill it with memA information for future use\n",memA_name);
      set_membrane(out_file, mem_f, mem_A);
      fprintf (out_file,"\nFinished with set_membrane function now we have a mem_A\n");
      char* memA_out = new char[1024];
      sprintf(memA_out,"data/shape-randst/membrane_files/memA-%4.02f-%4.02f-%4.02f-%4.02f-%4.0f.dat",A,B,C,D,density_factor);
      FILE *memAout = fopen((const char *)memA_out,"w");
      for (int i=0;i<Nx*Ny*Nz;i++) {
        fprintf(memAout, "%g\t",mem_A[i]);
      }
      fclose(memAout);
      delete[] memA_out;
      fprintf(out_file,"\nfinished printing the memA file, now we're moving on with simulation\n");
    } else {
      fprintf(out_file,"We're taking the memA info from a file that already exists\n");
      for (int i=0;i<Nx*Ny*Nz;i++) {
        if (fscanf(memAin, "%lg\t",&mem_A[i])!=1) {
          fprintf(out_file,"There was a problem in trying to read into the mem_A array! RUN!!!!\n");
            exit(1);
        }
      }
    }
  } else {
    set_membrane(out_file, mem_f, mem_A);
    fprintf (out_file,"\nFinished with set_membrane function now we have a mem_A and its not randst and Nx is = %d\n",Nx);
  }

  //begin area rating
  char *area_rating_out = new char[1024];
  sprintf(area_rating_out, "data/shape-%s/area_rating-%4.02f-%4.02f-%4.02f-%4.02f-%4.02f.dat",mem_f_shape.c_str(),A,B,C,D,density_factor);
  FILE *area_rating_file = fopen((const char *)area_rating_out,"w");
  if (area_rating_file == NULL){
    printf("WAAAAAAAAAAAAAAAAAAAAAA - area_rating_file == null \n");
  }

  char * area_rating_out_two = new char[1024];
  sprintf(area_rating_out_two, "data/shape-%s/area_rating_two-%4.02f-%4.02f-%4.02f-%4.02f-%4.02f.dat",mem_f_shape.c_str(),A,B,C,D,density_factor);
  FILE *area_rating_file_two = fopen((const char *)area_rating_out_two,"w");
  if (area_rating_file_two == NULL){
    printf("WAAAAAAAAAAAAAAAAAAAAAA - area_rating_file_two == null \n");
  }

  fprintf(out_file,"finished opening area_rating file\n");
  set_insideArr(insideArr);
  fprintf(out_file,"Finished with inside Arr function\n");

  double total_cell_volume = 0;
  double total_cell_area = 0;
  for (int i=0;i<Nx*Ny*Nz;i++){
    if (insideArr[i]==true) {
      total_cell_volume += dx*dx*dx;
    }
  }
  for (int i=0;i<Nx;i++){
    for (int j=0;j<Ny;j++){
      for (int k=0;k<Nz;k++){
        total_cell_area += mem_A[i*Ny*Nz+j*Nz+k];
        if(i==int(Nx/2)){  //FLAG: fix this for microscopy file
          if (insideArr[i*Ny*Nz+j*Nz+k]==true){
            double area_rating = 0;
            double area_rating_two = 0;
            for (int i2=0;i2<Nx;i2++){
              for (int j2=0;j2<Ny;j2++){
                for (int k2=0;k2<Nz;k2++){
                  if(i2!=i && j2!=j && k2!=k){
                    double dis = dx*sqrt((i-i2)*(i-i2)+(j-j2)*(j-j2)+(k-k2)*(k-k2));
                    area_rating += mem_A[i2*Ny*Nz+j2*Nz+k2]/(dis*dis);
                    if (dis<1.5*A){
                      area_rating_two += mem_A[i2*Ny*Nz+j2*Nz+k2];
                    }
                  }
                }
              }
            }
            fprintf(area_rating_file,"%g\t%g\t%g\n",j*dx,k*dx,area_rating);
            fprintf(area_rating_file_two,"%g\t%g\t%g\n",j*dx,k*dx,area_rating_two);
          } else {
            fprintf(area_rating_file,"%g\t%g\t%g\n",j*dx,k*dx,0.0);
            fprintf(area_rating_file_two,"%g\t%g\t%g\n",j*dx,k*dx,0.0);
          }
        }
      }
    }
  }
  fclose(area_rating_file_two);
  fclose(area_rating_file);

  if(area_rating_flag==1) {
    exit(0);
  }
  fprintf(out_file,"area_rating_two file is using %g as the radius of the sphere its looking at areas in",A);
  fprintf(out_file,"Finished writing to the area_rating file!\n");
  fprintf(out_file,"Total cell volume = %g\nTotal cell area = %g\n",total_cell_volume,total_cell_area);
  fflush(out_file);
  //end area rating
  
  //begin membrane printing
  char* outfilename = new char[1024];
  sprintf(outfilename,"data/shape-%s/membrane_files/membrane-%4.02f-%4.02f-%4.02f-%4.02f-%4.02f.dat",mem_f_shape.c_str(),A,B,C,D,density_factor);
  FILE *out = fopen((const char *)outfilename,"w");
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
  fprintf(out_file,"\nMEMBRANE FILE PRINTED\n");
  //end membrane printing

  //begin mem_f printing for randst, tie fighter, triangle - possibly do this for other shapes if needed -- NEEDS UPDATE.
  if (mem_f_shape == "randst"||mem_f_shape == "TIE_fighter"||mem_f_shape == "triangle") {
    char *f_file_name = new char[1024];
    if(f_file_name==NULL){fprintf(out_file,"OOOOOOOOOOOOOOOH no.");exit(1);}
    sprintf(f_file_name,"data/shape-%s/membrane_files/f_membrane-%4.02f-%4.02f-%4.02f-%4.02f-%4.02f.dat", mem_f_shape.c_str(),A,B,C,D,density_factor);
    FILE *f_file = fopen((const char *)f_file_name,"w");
    double x = Nx/2.0*dx;
    for (int i=0;i<Ny;i++) {
      for (int j=0;j<Nz;j++) {
        fprintf(f_file,"%g\t%g\t%g\n",i*dx,j*dx,mem_f(x,i*dx,j*dx));
      }
    }
    fclose(f_file);
    fprintf(out_file,"Finished printing out the mem_f_shape function!\n");
    fflush(stdout);
  }
  fprintf (out_file,"membrane set with density in it!\n");
  fflush(out_file);
  //end mem_f printing

  set_density(nATP, nE, mem_A);
  printf("nATP_starting_density = %g \n nE_starting_density = %g \n Density_factor = %g",
         nATP_starting_density, nE_starting_density, density_factor);
  fprintf (out_file,"nATP_starting density = %g and nE_starting_density = %g and Density_factor = %g",
           nATP_starting_density, nE_starting_density, density_factor);
  fflush(out_file);

  double bef_total_NATP=0;
  double bef_total_NADP=0;
  double bef_total_NE=0;
  double bef_total_Nde=0;
  double bef_total_Nd=0;
  double bef_total_N=0;
  double total_NATP=0;
  double total_NADP=0;
  double total_NE=0;
  double total_Nde=0;
  double total_Nd=0;
  double total_N=0;
  for (int i=0;i<Nx*Ny*Nz;i++){
    bef_total_NATP += nATP[i]*dx*dx*dx;
    bef_total_NADP += nADP[i]*dx*dx*dx;
    bef_total_NE += nE[i]*dx*dx*dx;
    bef_total_Nde += Nde[i];
    bef_total_Nd += Nd[i];
  }
  bef_total_N = bef_total_NATP*2 + bef_total_NADP*2 + bef_total_NE + bef_total_Nde*3 + bef_total_Nd*2;

  int percent = int(iter/100);
  double time_for_percent;
  bool check = true;
  //  clock_t newtime; unused
  clock_t start = clock();
  int k=0;
  for (int i=0;i<iter;i++){
    get_J(difD, nATP, nADP, nE, JxATP, JyATP,
          JzATP, JxADP, JyADP, JzADP, JxE, JyE, JzE);
    get_next_density(mem_A, insideArr, nATP, nADP, nE, Nd, Nde, JxATP, JyATP, JzATP,
                     JxADP, JyADP, JzADP, JxE, JyE, JzE);
    
    //begin status update
    if (i%percent == 0){
      if (i!=0){
        clock_t newtime = clock();
        if (check){
          time_for_percent = double(newtime - start)/CLOCKS_PER_SEC;
          check = false;
        }
        int percents_to_go = int(iter/percent - i/percent);
        if(percents_to_go%10==0 || percents_to_go == 99){
          fprintf(out_file,"We are %d percent complete and have %f seconds to go!\n",
                 i/percent, percents_to_go*time_for_percent);
        }
      }
    }
    //end status update

    //begin file printing
    if (i%iter_at_half_sec == 0) {
      fprintf(out_file,"******this is printing at iteration number = %d\n\n",i);
      //if(i>30){exit(1);}
      //int k = i/iter_at_half_sec;
      char *outfilenameATP = new char[1000];
      sprintf(outfilenameATP, "data/shape-%s/m_natp-%s-%03.2f-%03.2f-%03.2f-%03.2f-%03.2f-%03d.dat", argv[1],argv[1],A,B,C,D,density_factor,k);
      FILE *nATPfile = fopen((const char *)outfilenameATP,"w");
      delete[] outfilenameATP;
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
      fprintf(out_file,"printed out new file = natp\n");
      char *outfilenameE = new char[1000];
      sprintf(outfilenameE, "data/shape-%s/m_ne-%s-%03.2f-%03.2f-%03.2f-%03.2f-%03.2f-%03d.dat", argv[1],argv[1],A,B,C,D,density_factor,k);
      FILE *nEfile = fopen((const char *)outfilenameE,"w");
      delete[] outfilenameE;
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
      fprintf(out_file,"printed out new file = nadp\n");
      char *outfilenameADP = new char[1000];
      sprintf(outfilenameADP, "data/shape-%s/m_nadp-%s-%03.2f-%03.2f-%03.2f-%03.2f-%03.2f-%03d.dat", argv[1],argv[1],A,B,C,D,density_factor,k);
      FILE *nADPfile = fopen((const char *)outfilenameADP,"w");
      delete[] outfilenameADP;
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
      fprintf(out_file,"printed out new file = nadp\n");
      char *outfilenameD = new char[1000];
      sprintf(outfilenameD, "data/shape-%s/m_nd-%s-%03.2f-%03.2f-%03.2f-%03.2f-%03.2f-%03d.dat", argv[1],argv[1],A,B,C,D,density_factor,k);
      FILE *nDfile = fopen((const char *)outfilenameD,"w");
      delete[] outfilenameD;
      for (int a=0;a<Ny;a++){
        for (int b=0;b<Nz;b++){
          double nDsum = 0;
          for (int c=0;c<Nx;c++){
            nDsum += Nd[c*Ny*Nz+a*Nz+b];
          }
          fprintf(nDfile, "%1.2f ", nDsum);
        }
        fprintf(nDfile, "\n");
      }
      fclose(nDfile);
      fprintf(out_file,"printed out new file = nd\n");
      k++;
      fflush(out_file);
    }
  }
  //end file printing

  for (int i=0;i<Nx*Ny*Nz;i++){
    total_NATP += nATP[i]*dx*dx*dx;
    total_NADP += nADP[i]*dx*dx*dx;
    total_NE += nE[i]*dx*dx*dx;
    total_Nde += Nde[i];
    total_Nd += Nd[i];
  }
  total_N = total_NATP*2 + total_NADP*2 + total_NE + total_Nde*3 + total_Nd*2;
  fprintf(out_file,"total before NATP is = %f\n",bef_total_NATP);
  fprintf(out_file,"total before NADP is = %f\n",bef_total_NADP);
  fprintf(out_file,"total before NE is = %f\n",bef_total_NE);
  fprintf(out_file,"total before Nd is = %f\n",bef_total_Nd);
  fprintf(out_file,"total before Nde is = %f\n",bef_total_Nde);
  fprintf(out_file,"total before N is = %f\n",bef_total_N);
  fprintf(out_file,"total after NATP is = %f\n",total_NATP);
  fprintf(out_file,"total after NADP is = %f\n",total_NADP);
  fprintf(out_file,"total after NE is = %f\n",total_NE);
  fprintf(out_file,"total after Nd is = %f\n",total_Nd);
  fprintf(out_file,"total after Nde is = %f\n",total_Nde);
  fprintf(out_file,"total after N is = %f\n",total_N);
  fflush(out_file);
  cout << "Program has Run!!\n";
  fclose(out_file);

  //printing to a catalog for each directory so we know what we've run.
  char *fname = new char[1024];
  sprintf(fname,"data/shape-%s/catalog.txt",mem_f_shape.c_str());
  FILE * catalog;
  int catalog_exists;
  time_t rawtime;
  //  struct tm * timeinfo;
  time (&rawtime);
  //  timeinfo = localtime (&rawtime);
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
    fprintf(catalog," %s %1.2f %1.2f %1.2f %1.2f %1.2f\n", mem_f_shape.c_str(),A,B,C,D,density_factor);
    fclose(catalog);
  }
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
        if (xi == int(Nx/2)){
          //printf("x = %g y = %g, z = %g, f = %g\n",xi*dx,yi*dx,zi*dx,f);
        }
        mem_A[xi*Ny*Nz+yi*Nz+zi] = find_intersection(fXYZ, fXYz, fXyZ, fXyz, fxYZ, fxYz, fxyZ, fxyz, f);
        //printf(" x =%g y = %g z = %g f = %g\n",xi*dx,yi*dx,zi*dx,f);
        //fflush(stdout);
      }
    }
  }
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
  if (np > 6){cout << "There are more than six points in the cube!!!\n";}
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
  //  double F; unused
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
        de_to_ADP_E = rate_de*Nde[xi*Ny*Nz+yi*Nz+zi]/mem_A[xi*Ny*Nz+yi*Nz+zi]*time_step;
        ATP_to_d = (rate_D + rate_dD*(Nd[xi*Ny*Nz+yi*Nz+zi] + Nde[xi*Ny*Nz+yi*Nz+zi])/mem_A[xi*Ny*Nz+yi*Nz+zi])
        *nATP[xi*Ny*Nz+yi*Nz+zi]*time_step;
        E_d_to_de = rate_E*Nd[xi*Ny*Nz+yi*Nz+zi]/mem_A[xi*Ny*Nz+yi*Nz+zi]*nE[xi*Ny*Nz+yi*Nz+zi]*time_step;
        //Jeff!  remember that when you gain cyto density and lose the same amount of wall density,
        //the numbers of proteins gained/lost will be different, and it's the numbers that you want to be the same!!
        //also, keep thinking about the issue below where all additions are divided and then mult by the same mem_A fun
        nADP[xi*Ny*Nz+yi*Nz+zi] -= ADP_to_ATP;
        nATP[xi*Ny*Nz+yi*Nz+zi] += ADP_to_ATP;
        if (mem_A[xi*Ny*Nz+yi*Nz+zi] != 0){
          Nde[xi*Ny*Nz+yi*Nz+zi] += -de_to_ADP_E*mem_A[xi*Ny*Nz+yi*Nz+zi];
          nADP[xi*Ny*Nz+yi*Nz+zi] += de_to_ADP_E*mem_A[xi*Ny*Nz+yi*Nz+zi]/(dx*dx*dx);
          nE[xi*Ny*Nz+yi*Nz+zi] += de_to_ADP_E*mem_A[xi*Ny*Nz+yi*Nz+zi]/(dx*dx*dx);

          nATP[xi*Ny*Nz+yi*Nz+zi] -= ATP_to_d*mem_A[xi*Ny*Nz+yi*Nz+zi]/(dx*dx*dx);
          Nd[xi*Ny*Nz+yi*Nz+zi] += ATP_to_d*mem_A[xi*Ny*Nz+yi*Nz+zi];

          nE[xi*Ny*Nz+yi*Nz+zi] -= E_d_to_de*mem_A[xi*Ny*Nz+yi*Nz+zi]/(dx*dx*dx);
          Nd[xi*Ny*Nz+yi*Nz+zi] -= E_d_to_de*mem_A[xi*Ny*Nz+yi*Nz+zi];
          Nde[xi*Ny*Nz+yi*Nz+zi] += E_d_to_de*mem_A[xi*Ny*Nz+yi*Nz+zi];
        }
      }
    }
  }
  return 0;
}

double ran(){
	const long unsigned int x=0;
	static MTRand my_mtrand(x); // always use the same random number generator (for debugging)!
	return my_mtrand.randExc(); // which is the range of [0,1)
}

int set_density(double *nATP, double *nE, double *mem_A){
  int count_inside = 0;
  for (int i=0;i<Nx;i++){
    for (int j=0;j<Ny;j++){
      for (int k=0;k<Nz;k++){
        if (inside(i,j,k)){
          count_inside++; 
        }
      }
    }
  }

  //  double NE_per_cell = 1000*dx*dx*dx; //possible change to 350
  //  double NATP_per_cell = 350*dx*dx*dx; //possible change to 1000

  //  double NE_stdev = sqrt(NE_per_cell); unused
  //  double NATP_stdev = sqrt(NATP_per_cell); unused

  //  printf("Total inside = %d\nTotal nE should be = %f\n NE per cell = %f\n", count_inside,count_inside*NE_per_cell, NE_per_cell);
  //  printf("Total inside = %d\nTotal nATP should be = %f\n NE per cell = %f\n", count_inside,count_inside*NE_per_cell, NE_per_cell);

  //double r2,U,V; unused
  int right_most_point_z=0;
  int left_most_point_z=Nz;
  int right_most_point_y=0;
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
  int density_divider_z = int(right_most_point_z - (right_most_point_z - left_most_point_z)/2.0);
  //  int density_divider_y = int(right_most_point_y - (right_most_point_y - left_most_point_y)/2.0); unused
  for (int i=0;i<Nx;i++){
    for (int j=0;j<Ny;j++){
      for (int k=0;k<Nz;k++){
        if (inside(i,j,k)){
          if(k>density_divider_z){
            nATP[i*Ny*Nz+j*Nz+k] = nATP_starting_density*density_factor;
          } else {
            nATP[i*Ny*Nz+j*Nz+k] = nATP_starting_density;
          }
        }
      }
    }
  }
  for (int i=0;i<Nx;i++){
    for (int j=0;j<Ny;j++){
      for (int k=0;k<Nz;k++){
        if (inside(i,j,k)){
          if(k>density_divider_z){
            nE[i*Ny*Nz+j*Nz+k] = nE_starting_density*density_factor;
          }
          else {
            nE[i*Ny*Nz+j*Nz+k] = nE_starting_density;
          }
        }
      }
    }
  }
  for (int i=0;i<Nx*Ny*Nz;i++) {
    nADP[i] = 0;
    Nd[i] = 0;
    Nde[i] = 0;
  }
  return 0;
}
