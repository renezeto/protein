extern double *nATP;
extern double *nADP;
extern double *nE;
extern double *Nd;
extern double *Nde;
extern double *f_mem;
#include <stdio.h>


int set_density(double *nATP, double *nE, double *mem_A);
int get_next_density(double *mem_A, bool *insideArr, double *nATP, double *nADP,
                     double *nE, double *Nd, double *Nde,
                     double *JxATP, double *JyATP, double *JzATP,
                     double *JxADP, double *JyADP, double *JzADP,
                     double *JxE, double *JyE, double *JzE);
int get_J(double difD, double *nATP, double *nADP, double *nE,
          double *JxATP, double *JyATP, double *JzATP,
          double *JxADP, double *JyADP, double *JzADP,
          double *JxE, double *JyE, double *JzE);
void set_membrane(FILE * out_file, double (*mem_f)(double x, double y, double z), double mem_A[]);
void set_curvature(double mem_A[], double curvature[]);
void set_insideArr(bool *insideArr);
bool inside(int xi, int yi, int zi);
double find_intersection(const double fXYZ, const double fXYz, const double fXyZ, const double fxYZ,
		       const double fxyZ, const double fxYz, const double fXyz, const double fxyz,
		       const double f_minus_C);
