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


double difD = 2.5; double difE = 2.5;
double rate_ADP_ATP = 1;
double rate_D = .025;
double rate_dD = .0015;
double rate_de = .7;
double rate_E = .093;

const int n = 706;

const double dx=0.05;
const double tot_time = 86;
const double time_step = .1*dx*dx/difD;
const int iter = int(tot_time/time_step)+3;
const int iter_at_five_sec = int(5.0/time_step)+1;

double x, y, z;

int Nx;
int Ny;
int Nz;

double *nATP;
double *nADP;
double *nE;
double *Nd;
double *Nde;
double *f_mem;

string mem_f_shape;
double A;
double B;
double C;
double D;

double mem_f(double x, double y, double z) {
  if (mem_f_shape=="p"){ //pill
    //A = length, B = radius of endcap and cylinder, C = ???
    double f;
    double X = Nx*dx;
    double Y = Ny*dx;
    double Z = Nz*dx;
    double z1 = (Z-A)/2;
    //printf("z1 = %f and Z = %f and A = %f and B = %f\n",z1,Z,A,B);
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
  if (mem_f_shape=="c"){ //cone
    //A = length of cone, B = radius of base, C = ???
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
  if (mem_f_shape=="sp"){ //sphere
    // A = radius, B = ???, C = ???
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
	if (mem_f_shape=="e"){ //ellipsoid B = x axis radius radius, C = y axis radius radius, A = z axis radius radius
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
  if (argc != 2){
    printf("usage: %s mem_f_shape\n", argv[0]);
  }
  mem_f_shape = argv[1];
  A = atof(argv[2]);
  B = atof(argv[3]);
  C = atof(argv[4]);
  D = atof(argv[5]);
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
  printf("Nx=%d\nNy=%d\nNz=%d\nX=%f\nY=%f\nZ=%f\n",Nx,Ny,Nz,(Nx*dx),(Ny*dx),(Nz*dx));
  nATP = new double[Nx*Ny*Nz];
  nADP = new double[Nx*Ny*Nz];
  nE = new double[Nx*Ny*Nz];
  Nd = new double[Nx*Ny*Nz];
  Nde = new double[Nx*Ny*Nz];
  f_mem = new double[Nx*Ny*Nz];
  printf("For this simulation,\ndx = %f\ntot_time = %f\ntimestep = %f\ntotal iterations = %d\niter at five sec = %d\n",
         dx, tot_time, time_step, iter, iter_at_five_sec);
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
  printf("here\n");
  set_membrane(mem_f, mem_A);
  printf ("jheretwo\n");
  set_insideArr(insideArr);
  printf("herethree\n");
  set_density(nATP,nE, mem_A);
  for (int a=0;a<Ny;a++){
    for (int b=0;b<Nz;b++){
      //printf("nATP starts as %g\n", nATP[(int(Nx/2))*Ny*Nz+a*Nz+b]);
    }
  }
  printf ("membrane set with density in it!\n");
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
  printf("herefive\n");
  for (int i=0;i<Nx*Ny*Nz;i++){
    bef_total_NATP += nATP[i]*dx*dx*dx;
    bef_total_NADP += nADP[i]*dx*dx*dx;
    bef_total_NE += nE[i]*dx*dx*dx;
    bef_total_Nde += Nde[i];
    bef_total_Nd += Nd[i];
  }
  printf("heresix\n");
  bef_total_N = bef_total_NATP*2 + bef_total_NADP*2 + bef_total_NE + bef_total_Nde*3 + bef_total_Nd*2;
  //moved membrane
  printf("herefour\n");
  const char* outfilename = "membrane.dat";
	FILE *out = fopen((const char *)outfilename,"w");
  double marker;
  double inmarker;
  double zt = A/2; double yt = B/2; double xt = C/2;
  double ft = mem_f(zt,yt,xt);
  printf("in the box has mem_f function = %f\n",ft);
  printf("what is Nx/2 = %f\n",Nz/2.0);
  for (int j=0;j<Ny;j++){ // possible source of membrane.dat problem?
    for (int i=0;i<Nz;i++){
      if (insideArr[(int(Nx/2))*Ny*Nz+j*Nz+i]==true) {inmarker = 1;}
      else {inmarker = 0;}
      if (mem_A[(int(Nx/2))*Ny*Nz+j*Nz+i]!=0) {marker = 1;}
      else {marker = 0;}
      fprintf(out, "%g  ", inmarker);
    }
    fprintf(out, "\n");
  }
  fflush(stdout);
  fclose(out);
  //end of membrane
  printf("\nMEMBRANE FILE PRINTED\n");
  char *outfilenameStart = new char[1000];
  sprintf(outfilenameStart, "starting_natp.dat");
  FILE *nATPStartfile = fopen((const char *)outfilenameStart,"w");
  delete[] outfilenameStart;
  for (int a=0;a<Ny;a++){
    for (int b=0;b<Nz;b++){
      //printf("nATP is %g\n", nATP[(int(Nx/2))*Ny*Nz+a*Nz+b]);
      fprintf(nATPStartfile, "%1.2f ", nATP[(int(Nx/2))*Ny*Nz+a*Nz+b]);
    }
    fprintf(nATPStartfile, "\n");
  }
  fclose(nATPStartfile);
  int percent = int(iter/100);
  double time_for_percent;
  bool check = true;
  clock_t newtime;
  clock_t start = clock();
  int k=0;
  for (int i=0;i<iter;i++){
    get_J(difD, nATP, nADP, nE, JxATP, JyATP,
          JzATP, JxADP, JyADP, JzADP, JxE, JyE, JzE);
    get_next_density(mem_A, insideArr, nATP, nADP, nE, Nd, Nde, JxATP, JyATP, JzATP,
                     JxADP, JyADP, JzADP, JxE, JyE, JzE);
    if (i%percent == 0){
      if (i!=0){
        clock_t newtime = clock();
        if (check){
          time_for_percent = double(newtime - start)/CLOCKS_PER_SEC;
          check = false;
        }
        int percents_to_go = int(iter/percent - i/percent);
        if(percents_to_go%10==0 || percents_to_go == 99){
          printf("We are %d percent complete and have %f seconds to go!\n",
                 i/percent, percents_to_go*time_for_percent);
        }
      }
    }
    if (i%iter_at_five_sec == 0){printf("did this work????????????????? = %d\n",i);}
    //printf("iter_at_five_sec = %d\n\n",iter_at_five_sec);
    if (i%iter_at_five_sec == 0 || ((i*125)%iter_at_five_sec == 0 && i<iter_at_five_sec)) {
      printf("******this is printing at iteration number = %d\n\n",i);
      fflush(stdout);
      //if(i>30){exit(1);}
      //int k = i/iter_at_five_sec;
      char *outfilenameATP = new char[1000];
      sprintf(outfilenameATP, "shape-%s/natp-%s-%03.1f-%03.1f-%03.1f-%03.1f-%03d.dat", argv[1],argv[1],A,B,C,D,k);
      FILE *nATPfile = fopen((const char *)outfilenameATP,"w");
      delete[] outfilenameATP;
      for (int a=0;a<Ny;a++){
        for (int b=0;b<Nz;b++){
          fprintf(nATPfile, "%1.2f ", nATP[(int(Nx/2))*Ny*Nz+a*Nz+b]);
        }
        fprintf(nATPfile, "\n");
      }
      fclose(nATPfile);
      printf("printed out new file = natp\n");
      char *outfilenameE = new char[1000];
      sprintf(outfilenameE, "shape-%s/ne-%s-%03.1f-%03.1f-%03.1f-%03.1f-%03d.dat", argv[1],argv[1],A,B,C,D,k);
      FILE *nEfile = fopen((const char *)outfilenameE,"w");
      delete[] outfilenameE;
      for (int a=0;a<Ny;a++){
        for (int b=0;b<Nz;b++){
          fprintf(nEfile, "%1.2f ", nE[(int(Nx/2))*Ny*Nz+a*Nz+b]);
        }
        fprintf(nEfile, "\n");
      }
      fclose(nEfile);
      printf("printed out new file = nadp\n");
      char *outfilenameADP = new char[1000];
      sprintf(outfilenameADP, "shape-%s/nadp-%s-%03.1f-%03.1f-%03.1f-%03.1f-%03d.dat", argv[1],argv[1],A,B,C,D,k);
      FILE *nADPfile = fopen((const char *)outfilenameADP,"w");
      delete[] outfilenameADP;
      for (int a=0;a<Ny;a++){
        for (int b=0;b<Nz;b++){
          fprintf(nADPfile, "%1.2f ", nADP[(int(Nx/2))*Ny*Nz+a*Nz+b]);
        }
        fprintf(nADPfile, "\n");
      }
      fclose(nADPfile);
      printf("printed out new file = nadp\n");
      char *outfilenameD = new char[1000];
      sprintf(outfilenameD, "shape-%s/nd-%s-%03.1f-%03.1f-%03.1f-%03.1f-%03d.dat", argv[1],argv[1],A,B,C,D,k);
      FILE *nDfile = fopen((const char *)outfilenameD,"w");
      delete[] outfilenameD;
      for (int a=0;a<Ny;a++){
        for (int b=0;b<Nz;b++){
          fprintf(nDfile, "%1.2f ", Nd[(int(Nx/2))*Ny*Nz+a*Nz+b]);
        }
        fprintf(nDfile, "\n");
      }
      fclose(nDfile);
      printf("printed out new file = nd\n");
      k++;
    }
  }
  for (int i=0;i<Nx*Ny*Nz;i++){
    total_NATP += nATP[i]*dx*dx*dx;
    total_NADP += nADP[i]*dx*dx*dx;
    total_NE += nE[i]*dx*dx*dx;
    total_Nde += Nde[i];
    total_Nd += Nd[i];
  }
  total_N = total_NATP*2 + total_NADP*2 + total_NE + total_Nde*3 + total_Nd*2;
  printf("total before NATP is = %f\n",bef_total_NATP);
  printf("total before NADP is = %f\n",bef_total_NADP);
  printf("total before NE is = %f\n",bef_total_NE);
  printf("total before Nd is = %f\n",bef_total_Nd);
  printf("total before Nde is = %f\n",bef_total_Nde);
  printf("total before N is = %f\n",bef_total_N);
  printf("total after NATP is = %f\n",total_NATP);
  printf("total after NADP is = %f\n",total_NADP);
  printf("total after NE is = %f\n",total_NE);
  printf("total after Nd is = %f\n",total_Nd);
  printf("total after Nde is = %f\n",total_Nde);
  printf("total after N is = %f\n",total_N);

  cout << "Program has Run!!\n";
  return 0;
}

//checks corners of each gridpoint, creates mem_A
void set_membrane(double (*mem_f)(double x, double y, double z),
		 double mem_A[]) {
  for(int xi=0;xi<Nx;xi++){
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

//unsure
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

//unsure
void set_insideArr(bool *insideArr){
  for(int xi=0;xi<Nx;xi++){
    for(int yi=0;yi<Ny;yi++){
      for(int zi=0;zi<Nz;zi++){
        insideArr[xi*Ny*Nz+yi*Nz+zi] = inside(xi,yi,zi);
      }
    }
  }
}

//checks midddle of each gridpoint to see if it is inside the cell or not
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

//unsure
int get_J(double difD, double *nATP, double *nADP, double *nE,
	  double *JxATP, double *JyATP, double *JzATP,
	  double *JxADP, double *JyADP, double *JzADP,
	  double *JxE, double *JyE, double *JzE){
  double F;
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
}


//unsure
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


//random # generator
double ran(){
	const long unsigned int x=0;
	static MTRand my_mtrand(x); // always use the same random number generator (for debugging)!
	return my_mtrand.randExc(); // which is the range of [0,1)
}


//density intializer
int set_density(double *nATP, double *nE, double *mem_A){
  int count_inside = 0;
  for (int i=0;i<Nx;i++){
    for (int j=0;j<Ny;j++){
      for (int k=0;k<Nz;k++){
        if (inside(i,j,k)){
          count_inside++; // counts # of inside gridpoints
        }
      }
    }
  }
  double NE_per_cell = 1000*dx*dx*dx;
  double NATP_per_cell = 350*dx*dx*dx;
  double NE_stdev = sqrt(NE_per_cell);
  double NATP_stdev = sqrt(NATP_per_cell);
  printf("NATP_per_cell = %g and NATP_stdev = %g\n", NATP_per_cell, NATP_stdev);
  printf("total inside = %d\nTotal nE should be = %f\nE_per_cell = %f\n", count_inside,
  count_inside*NE_per_cell, NE_per_cell);
  double r2,U,V;
  for (int i=0;i<Nx;i++){
    for (int j=0;j<Ny;j++){
      for (int k=0;k<Nz;k++){
        if (inside(i,j,k)){
          if (NATP_per_cell <= 1){
            if (ran() <= NATP_per_cell){
              nATP[i*Ny*Nz+j*Nz+k] = 1/(dx*dx*dx);
            }
            else {
              nATP[i*Ny*Nz+j*Nz+k] = 0;
            }
            if(k>2*Nz/3){
              if (ran() <= 2.5*NATP_per_cell){
                nATP[i*Ny*Nz+j*Nz+k] = 1/(dx*dx*dx);
              }
              else {
                nATP[i*Ny*Nz+j*Nz+k] = 0;
              }
            }
          } else {exit(1);}
        }
      }
    }
  }
  for (int i=0;i<Nx;i++){
    for (int j=0;j<Ny;j++){
      for (int k=0;k<Nz;k++){
        if (inside(i,j,k)){
         if (NE_per_cell <= 1){
            if (ran() <= NE_per_cell){
              nE[i*Ny*Nz+j*Nz+k] = 1/(dx*dx*dx);
            }
            else {
              nE[i*Ny*Nz+j*Nz+k] = 0;
            }
            if(k>2*Nz/3){
              if (ran() <= 2.5*NE_per_cell){
                nE[i*Ny*Nz+j*Nz+k] = 1/(dx*dx*dx);
              }
              else {
                nE[i*Ny*Nz+j*Nz+k] = 0;
              }
            }
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
