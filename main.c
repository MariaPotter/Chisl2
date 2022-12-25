#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

_Pragma ("GCC diagnostic push")
_Pragma ("GCC diagnostic ignored \"-Wunused-parameter\"")
static double x_d  (double t, double x, double y, double p_x, double p_y, double mult)        
{
    return y;
}

static double y_d (double t, double x, double y, double p_x, double p_y, double mult)
{
    return p_y/2.-sqrt(2.)*x*exp(-mult*t);
}

static double px_d(double t, double x, double y, double p_x, double p_y, double mult)
{
    return 2.*x+sqrt(2.)*p_y*exp(-mult*t);
}

static double py_d (double t, double x, double y, double p_x, double p_y, double mult)
{
    return -p_x;
}

_Pragma ("GCC diagnostic pop")
    
double rd (FILE* inpf)        // reading floating point number from a/b format
{
    int chisl = 0;
    int znam  = 1;
    char tmp;
    _Pragma ("GCC diagnostic push");
    _Pragma ("GCC diagnostic ignored \"-Wunused-result\"");
    fscanf (inpf, "%d", &chisl);
    if (fscanf (inpf, "%c", &tmp) && tmp == '/') fscanf (inpf, "%d", &znam);
    _Pragma ("GCC diagnostic pop");
    return ((double)chisl) / znam;
}

// Runge–Kutta in general form (one step)
void RK (double t,
         double* x0,
         double* y0,
		 double* px0,
         double* py0,
         double h,
         unsigned int s,
         double** k,
         double** cab,
         double f (double, double, double, double, double, double),
         double g (double, double, double, double, double, double),
		 double u (double, double, double, double, double, double),
         double v (double, double, double, double, double, double),
         double mult,
         bool check)
{
    double temp_x;
    double temp_y;
	double temp_px;
    double temp_py;

    for (unsigned int i = 0; i < s; i++)
    {
        temp_x = 0;
        temp_y = 0;
        temp_px = 0;
        temp_py = 0;
        for (unsigned int j = 0; j < i; j++)
        {
            temp_x += cab[i + 1][j] * k[0][j];
            temp_y += cab[i + 1][j] * k[1][j];
			temp_px += cab[i + 1][j] * k[2][j];
			temp_py += cab[i + 1][j] * k[3][j];
        }
        // printf("%f\n",f(1,0));
        k[0][i] = g (t + h * cab[0][i],
					*x0 + h * temp_x,
					*y0 + h * temp_y,
					*px0 + h * temp_px,
					*py0 + h * temp_py,
                    mult);
					
        k[1][i] = f (t + h * cab[0][i],
					*x0 + h * temp_x,
					*y0 + h * temp_y,
					*px0 + h * temp_px,
					*py0 + h * temp_py,
                    mult);
					
		k[2][i] = u (t + h * cab[0][i],
					*x0 + h * temp_x,
					*y0 + h * temp_y,
					*px0 + h * temp_px,
					*py0 + h * temp_py,
                    mult);
					
		k[3][i] = v (t + h * cab[0][i],
					*x0 + h * temp_x,
					*y0 + h * temp_y,
					*px0 + h * temp_px,
					*py0 + h * temp_py,
                    mult);
    }

    temp_x = 0;
    temp_y = 0;
	temp_px = 0;
    temp_py = 0;
	
    if (check)
    {
        for (unsigned int i = 0; i < s; i++)
        {
            temp_x += cab[s + 2][i] * k[0][i];
            temp_y += cab[s + 2][i] * k[1][i];
			temp_px += cab[s + 2][i] * k[2][i];
			temp_py += cab[s + 2][i] * k[3][i];
        }
    }
    else
    {
        for (unsigned int i = 0; i < s; i++)
        {
            temp_x += cab[s + 1][i] * k[0][i];
            temp_y += cab[s + 1][i] * k[1][i];
			temp_px += cab[s + 1][i] * k[2][i];
			temp_py += cab[s + 1][i] * k[3][i];
        }
    }

    *x0 += h * temp_x;
    *y0 += h * temp_y;
	*px0 += h * temp_px;
    *py0 += h * temp_py;
}

// Adaptive Runge-Kutta
double astep (double T,
              double* x,
              double* y,
			  double* px,
			  double* py,
			  
              double* x_,
              double* y_,
			  double* px_,
			  double* py_,
			  
              long long unsigned* i,
              long long unsigned* j,
			  
              unsigned int p,
              unsigned int s,
			  
              double** k,
              double** cab,
			  
              double tol,
			  
              double f (double, double, double, double, double, double),
			  double g (double, double, double, double, double, double),
			  double u (double, double, double, double, double, double),
			  double v (double, double, double, double, double, double),
              double mult)
{
    double temp_x, temp_y, temp_px, temp_py, dist, h, temp, fac, err, norm;
    *x_ = temp_x = *x;
    *y_ = temp_y = *y;
    
	*px_ = temp_px = *px;
    *py_ = temp_py = *py;
    
    dist        = 0;
    h           = 0.01;
    fac         = 1.7;
    err         = 0;
    norm        = 0;
	
    printf("%.2e    %.6e  |  %.6e    %.2e\n", *x, *y, *px, *py);
	
    for (*i = *j = 0; T - dist > pow(10,-17);)
    {   
        
        if (dist + h > T) h = T - dist;

        RK (dist, &temp_x, &temp_y, &temp_px, &temp_py, h, s, k, cab, f, g, u, v, mult, false);
        RK (dist, x_, y_, px_, py_, h, s, k, cab, f, g, u, v, mult, true);
        
        
        norm = fmax( fmax( fabs (temp_x - *x_), fabs (temp_y - *y_)), fmax(fabs (temp_px - *px_), fabs (temp_py - *py_)));
        temp = h;
        h *= fmin (fac,
                   fmax (0.7,
                         pow (0.98 * tol / norm,
                              1. / (p + 1))));
        if (h < pow (10, -18))
        {
            printf ("\nSomething goes wrong...\n");
            return 0;
        }
        if (norm > tol)
        {
            *x_ = temp_x = *x;
            *y_ = temp_y = *y;

            *px_ = temp_px = *px;
            *py_ = temp_py = *py;
			
            fac         = 1;
            *j += 1;
            continue;
        }

        err += norm;
        dist += temp;
        *x  = temp_x;
        *y  = temp_y;
		*px = temp_px;
		*py = temp_py;
        fac = 1.7;
        *i += 1;
    }
    //printf("%.2e    %.2e  |  %.2e    %.2e\n", *x, *y, *px, *py);
    return err;
}

// discrepancy counting
double error(double T,
	        double x,
	        double y,
	        double px,
	        double py,
	        double* err,
	        unsigned int p,
	        unsigned int s,

	        double** k,
	        double** cab,

	        double tol,

	        double f (double, double, double, double, double, double),
	        double g (double, double, double, double, double, double),
	        double u (double, double, double, double, double, double),
	        double v (double, double, double, double, double, double),
            double mult)
{
	double x_ = x;
	double y_ = y;
	double px_ = px;
	double py_ = py;
	long long unsigned i = 0;
	long long unsigned j = 0;
    double delt;

	delt = astep(T, &x, &y, &px, &py, &x_, &y_, &px_, &py_, &i, &j, p, s, k, cab, tol, f, g, u, v, mult);

	err[0] = x;
	err[1] = py;
    
    return delt;
}

/* 
void inverse_matrix(double** A) {
    double temp = 0.0;
    double MATR[2][2];
    for(int i=0;i<2;i++)
        for(int j=0;j<2;j++)
            MATR[i][j]=A[i][j];
    temp = 1/(MATR[0][0]*MATR[1][1] - MATR[0][1]*MATR[1][0]);
    A[0][0] = temp * MATR[1][1];
    A[0][1] = -temp * MATR[0][1];
    A[1][0] = -temp * MATR[1][0];
    A[1][1] = temp * MATR[0][0];
}

double fedorenko_norm(double** A, double* B) {
    return sqrt(pow(B[0],2)/(pow(A[0][0],2) + pow(A[0][1],2)) + pow(B[1],2)/(pow(A[1][0],2) + pow(A[1][1],2)));
}
*/   

// the shooting method
void shooting_method(unsigned int max_iterations,
                     double T,
	                 double x,
	                 double py,
	                 double* err,
	                 unsigned int p,
	                 unsigned int s,

	                 double** k,
	                 double** cab,

	                 double tol,

	                 double f (double, double, double, double, double, double),
	                 double g (double, double, double, double, double, double),
	                 double u (double, double, double, double, double, double),
	                 double v (double, double, double, double, double, double),
                     double mult)
{
    double alpha_1 = 1; //x
    double alpha_2 = 1; //py
    double** A = (double**)malloc (2*sizeof(double*));
    A[0] = (double*)malloc (2*sizeof(double));
    A[1] = (double*)malloc (2*sizeof(double));
    double* B = (double*)malloc (2*sizeof(double));
    double* C = (double*)malloc (2*sizeof(double));
    double delta = pow(10, -6);
    double temp_alpha[2];

    for(unsigned int iteration = 0; iteration < max_iterations; iteration++){
        error(T, x, alpha_1+delta, alpha_2, py, err, p, s, k, cab, tol, f, g, u, v, mult);
        temp_alpha[0] = err[0];  // y
        temp_alpha[1] = err[1];  // px
        error(T, x, alpha_1-delta, alpha_2, py, err, p, s, k, cab, tol, f, g, u, v, mult);
        A[0][0] = (temp_alpha[0]-err[0])/(2*delta);
        A[1][0] = (temp_alpha[1]-err[1])/(2*delta);

        error(T, x, alpha_1, alpha_2+delta, py, err, p, s, k, cab, tol, f, g, u, v, mult);
        temp_alpha[0] = err[0];  // y
        temp_alpha[1] = err[1];  // px
        error(T, x, alpha_1, alpha_2-delta, py, err, p, s, k, cab, tol, f, g, u, v, mult);
        A[0][1] = (temp_alpha[0]-err[0])/(2*delta);
        A[1][1] = (temp_alpha[1]-err[1])/(2*delta);

        error(T, x, alpha_1, alpha_2, py, err, p, s, k, cab, tol, f, g, u, v, mult);
        B[0] = err[0];  //y
        B[1] = err[1];  //px
        
        printf("%.2e    %.2e  |  %.2e\n%.2e    %.2e  |  %.2e\n", A[0][0], A[0][1], B[0], A[1][0], A[1][1], B[1]);
        
        if(sqrt(pow(B[0],2)+pow(B[1],2)) < tol*pow(10,2)){
            //printf and break
            printf("Shooting iteration: %ud,\nAlpha: %.3f,\nX(0)=%.3e,\nY(0)=%.3e,\nDiscrepancy: %.3e\n", iteration, mult,alpha_1,alpha_2,sqrt(pow(B[0],2)+pow(B[1],2)));
            break;
        }
        
        printf("Shooting iteration: %ud,\nAlpha: %.3f,\nX(0)=%.3e,\nY(0)=%.3e,\nDiscrepancy: %.3e\n\n", iteration, mult,alpha_1,alpha_2,sqrt(pow(B[0],2)+pow(B[1],2)));
        
        C[0]=B[0]-A[0][0]*alpha_1-A[0][1]*alpha_2;
        C[1]=B[1]-A[1][0]*alpha_1-A[1][1]*alpha_2;

        alpha_1 = (A[1][1]*C[0]-A[0][1]*C[1])/(A[1][0]*A[0][1]-A[0][0]*A[1][1]);
        alpha_2 = (A[1][1]*C[1]-A[0][1]*C[0])/(A[1][0]*A[0][1]-A[0][0]*A[1][1]);
                
        //changing alpha_1, alpha_2
    }
    free (A[0]);
    free (A[1]);
    free (A);
    free (B);
    free (C);
}

int main()
{
    double T             = 1;
    double tol           = pow (10, -7);
    long long unsigned i = 0;
    long long unsigned j = 0;

    double x=1;
	//double y;
	//double px;
	double py=0;
	
    double*  mult; // alpha из условия
    double* alpha; // под невязки

    unsigned int s;
    unsigned int p;
    double** k;
    double** cab;

    time_t start, end;

    /////////////////////////////////////////////////////////////////////////////////////////

    // file check

    time (&start);

    FILE* inpf = fopen ("koef (8).txt", "r");
    if (inpf == NULL)
    {
        printf ("File doen`t exist\n");
        return -1;
    }
    if (!fscanf (inpf, "%ud", &p) || !fscanf (inpf, "%ud", &s))
    {
        printf ("File isn`t correct\n");
        fclose (inpf);
        return -2;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    // CAB matrix initialization

    cab = (double**)malloc ((s + 3) * sizeof (double*));
    for (i = 0; i < s + 3; i++)
    {
        cab[i] = (double*)malloc (s * sizeof (double));
        for (j = 0; j < s; j++)
        {
            cab[i][j] = 0;
        }
    }
    
    mult = (double*)malloc (4 * sizeof (double));
    mult[0] = 0.0;
    mult[1] = 0.01;
    mult[2] = 0.1;
    mult[3] = 1;

	alpha = (double*)malloc (2 * sizeof (double));
    alpha[0] = 1;
    alpha[1] = 1;

    k    = (double**)malloc (4 * sizeof (double*));
    k[0] = (double*)malloc (s * sizeof (double));
    k[1] = (double*)malloc (s * sizeof (double));
	k[2] = (double*)malloc (s * sizeof (double));
    k[3] = (double*)malloc (s * sizeof (double));
	
    for (i = 0; i < s; i++)
    {
        k[0][i] = 0;
        k[1][i] = 0;
		k[2][i] = 0;
        k[3][i] = 0;
    }

    /////////////////////////////////////////////////////////////////////////////////////////

    // CAB matrix reading

    for (i = 0; i < s; i++)
    {
        cab[0][i] = rd (inpf);
    }

    for (i = 2; i < s + 3; i++)
        for (j = 0; j + 1 < i && j < s && !feof (inpf); j++)
            cab[i][j] = rd (inpf);

    // CAB matrix printing
    /*
            printf ("\n");
            for (unsigned int i = 0; i < s + 3; i++)
            {
                for (unsigned int j = 0; j < s; j++)
                    printf ("%6.3f ", cab[i][j]);
                printf ("\n");
            }
            printf ("\n");
        */
    ///////////////////////////////////////////////////////////////////////////////////////// 
    
    // the shooting method
	
    for(unsigned int l = 0; l < 1; l++)
        shooting_method(200, T, x, py, alpha, p, s, k, cab, tol, x_d, y_d, px_d, py_d, mult[l]);

    /////////////////////////////////////////////////////////////////////////////////////////

    // cleaning :)

    fclose (inpf);
    for (i = 0; i < s + 3; i++)
    {
        free (cab[i]);
    }
    free (cab);
	free (mult);
	free (alpha);
    free (k[0]);
    free (k[1]);
	free (k[2]);
    free (k[3]);
    free (k);

    time (&end);
    printf ("%.f seconds\n", difftime (end, start));
    return 0;
}
