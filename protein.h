

/* struct densities { */
/* public: */
/*   densities(int xnum, int ynum, int znum, double mydx); // allocates the arrays */
/*   //densities(double size, double mydx); // allocates the arrays */
/*   ~densities(); // free all the arrays */
/*   int Nx, Ny, Nz; */
/*   double dx; */
/*   double *Nd, *Nde, *nATP; */
/* private: */
/*   density(const density &); */
/*   void operator=(const density &); */
/* }; */
/*
extern const int Nx;
extern const int Ny;
extern const int Nz;
*/
extern double *nATP;
extern double *nADP;
extern double *nE;
extern double *Nd;
extern double *Nde;
extern double *f_mem;


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
void set_membrane(double (*mem_f)(double x, double y, double z), double mem_A[]);
void set_insideArr(bool *insideArr);
bool inside(double xi, double yi, double zi);
double find_intersection(const double fXYZ, const double fXYz, const double fXyZ, const double fxYZ,
		       const double fxyZ, const double fxYz, const double fXyz, const double fxyz,
		       const double f_minus_C);
