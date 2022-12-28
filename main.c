#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define pi 3.14159265358979323846264338328

_Pragma ("GCC diagnostic push")
_Pragma ("GCC diagnostic ignored \"-Wunused-parameter\"")
static double x_d  (double t, double x, double y, double p_x, double p_y, double mult)        
{
    return y;
}

static double y_d (double t, double x, double y, double p_x, double p_y, double mult)
{
    return p_y;
}

static double px_d(double t, double x, double y, double p_x, double p_y, double mult)
{
    return -24./(1.+mult*pow(y,2.));
}

static double py_d (double t, double x, double y, double p_x, double p_y, double mult)
{
    return -p_x +48.*x*y*mult/pow((1.+mult*pow(y,2.)),2.);
}

static double B_d  (double t, double x, double y, double p_x, double p_y, double mult)        
{
    return pow(p_y,2.) -48.*x/(1.+mult*pow(y,2.));
}

static double test_x_d  (double t, double x, double y, double p_x, double p_y, double mult)        
{
    return -y;
}

static double test_y_d (double t, double x, double y, double p_x, double p_y, double mult)
{
    return x;
}

static double test_px_d(double t, double x, double y, double p_x, double p_y, double mult)
{
    return -p_y;
}

static double test_py_d (double t, double x, double y, double p_x, double p_y, double mult)
{
    return p_x;
}

static double test_B_d  (double t, double x, double y, double p_x, double p_y, double mult)        
{
    return sqrt(pow(x-p_x,2) + pow(y-p_y,2));
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
         double* b0,
         double h,
         unsigned int s,
         double** k,
         double** cab,
         double f (double, double, double, double, double, double),
         double g (double, double, double, double, double, double),
		 double u (double, double, double, double, double, double),
         double v (double, double, double, double, double, double),
         double l (double, double, double, double, double, double),
         double mult,
         bool check)
{
    double temp_x;
    double temp_y;
	double temp_px;
    double temp_py;
    double temp_b;
    

    for (unsigned int i = 0; i < s; i++)
    {
        temp_x = 0;
        temp_y = 0;
        temp_px = 0;
        temp_py = 0;
        temp_b = 0;
        for (unsigned int j = 0; j < i; j++)
        {
            temp_x += cab[i + 1][j] * k[0][j];
            temp_y += cab[i + 1][j] * k[1][j];
			temp_px += cab[i + 1][j] * k[2][j];
			temp_py += cab[i + 1][j] * k[3][j];
			temp_b += cab[i + 1][j] * k[4][j];
        }
        // printf("%f\n",f(1,0));
        k[0][i] = f (t + h * cab[0][i],
					*x0 + h * temp_x,
					*y0 + h * temp_y,
					*px0 + h * temp_px,
					*py0 + h * temp_py,
                    mult);
					
        k[1][i] = g (t + h * cab[0][i],
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
        
        k[4][i] = l (t + h * cab[0][i],
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
    temp_b = 0;
    
	
    if (check)
    {
        for (unsigned int i = 0; i < s; i++)
        {
            temp_x += cab[s + 2][i] * k[0][i];
            temp_y += cab[s + 2][i] * k[1][i];
			temp_px += cab[s + 2][i] * k[2][i];
			temp_py += cab[s + 2][i] * k[3][i];
			temp_b += cab[s + 2][i] * k[4][i];
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
			temp_b += cab[s + 1][i] * k[4][i];           
        }
    }

    *x0 += h * temp_x;
    *y0 += h * temp_y;
	*px0 += h * temp_px;
    *py0 += h * temp_py;      
    *b0 += h * temp_b;       
}

// Adaptive Runge-Kutta
double astep (double T,
              double* x,
              double* y,
			  double* px,
			  double* py,
			  double* b,
			  
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
			  double l (double, double, double, double, double, double),
              double mult)
{
    double temp_x, temp_y, temp_px, temp_py, temp_b, x_, y_, px_, py_, b_, dist, h, temp, fac, err, norm;
    x_ = temp_x = *x;
    y_ = temp_y = *y;
    
	px_ = temp_px = *px;
    py_ = temp_py = *py;
    b_ = temp_b = *b;
    
    
    dist        = 0;
    h           = 0.01;
    fac         = 1.7;
    err         = 0;
    norm        = 0;
	
    // printf("%.2e    %.6e  |  %.6e    %.2e\n", *x, *y, *px, *py);
	
    for (*i = *j = 0; T - dist > tol*pow(10,-2);)
    {   
        
        if (dist + h > T) h = T - dist;

        RK (dist, &temp_x, &temp_y, &temp_px, &temp_py, &temp_b, h, s, k, cab, f, g, u, v, l, mult, false);
        RK (dist, &x_, &y_, &px_, &py_, &b_, h, s, k, cab, f, g, u, v, l, mult, true);
        
        
        // norm = fmax( fmax( fabs (temp_x - x_), fabs (temp_y - y_)), fmax(fabs (temp_px - px_), fabs (temp_py - py_)));
        norm = sqrt(pow(temp_x - x_,2)+pow(temp_y - y_,2)+pow(temp_px - px_,2)+pow(temp_py - py_,2));
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
            x_ = temp_x = *x;
            y_ = temp_y = *y;

            px_ = temp_px = *px;
            py_ = temp_py = *py;
			b_ = temp_b = *b;
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
		*b = temp_b;
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
	        double b,
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
	        double l (double, double, double, double, double, double),
            double mult)
{
	long long unsigned i = 0;
	long long unsigned j = 0;

	astep(T, &x, &y, &px, &py, &b, &i, &j, p, s, k, cab, tol, f, g, u, v, l, mult);

	err[0] = x-0;//
	err[1] = py-0;
    
    return b;
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
	                 double y,
	                 double px,
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
	                 double l (double, double, double, double, double, double),
                     double mult)
{
    double alpha_1 = -1.5; //x
    double alpha_2 = -0.5; //py
    double** A = (double**)malloc (2*sizeof(double*));
    A[0] = (double*)malloc (2*sizeof(double));
    A[1] = (double*)malloc (2*sizeof(double));
    double* B = (double*)malloc (2*sizeof(double));
    double delta = pow(10, -6);
    double integral = 0;

    for(unsigned int iteration = 0; iteration < max_iterations; iteration++){
        error(T, alpha_1+delta, y, px, alpha_2, integral, err, p, s, k, cab, tol, f, g, u, v, l, mult);
        B[0] = err[0];  // x
        B[1] = err[1];  // py
        error(T, alpha_1-delta, y, px, alpha_2, integral, err, p, s, k, cab, tol, f, g, u, v, l, mult);
        // printf("\n%.7e    %.7e    %.7e\n",B[0], err[0], (B[0]-err[0])/delta);
        // printf("%.7e    %.7e    %.7e\n\n",B[1], err[1], (B[1]-err[1])/delta);
        A[0][0] = (B[0]-err[0])/(2*delta);
        A[1][0] = (B[1]-err[1])/(2*delta);

        error(T, alpha_1, y, px, alpha_2+delta, integral, err, p, s, k, cab, tol, f, g, u, v, l, mult);
        B[0] = err[0];  // x
        B[1] = err[1];  // py
        error(T, alpha_1, y, px, alpha_2-delta, integral, err, p, s, k, cab, tol, f, g, u, v, l, mult);
        A[0][1] = (B[0]-err[0])/(2*delta);
        A[1][1] = (B[1]-err[1])/(2*delta);

        integral = error(T, alpha_1, y, px, alpha_2, integral, err, p, s, k, cab, tol, f, g, u, v, l, mult);
        B[0] = err[0];  //x
        B[1] = err[1];  //py
        
        // printf("%10.2e    %10.2e  |  %10.2e\n%10.2e    %10.2e  |  %10.2e\n", A[0][0], A[0][1], B[0], A[1][0], A[1][1], B[1]);
        
        if(sqrt(pow(B[0],2)+pow(B[1],2)) < tol*pow(10,2)){
            //printf and break
            printf("Shooting iteration: %u,\nAlpha: %.3f,\nX(0)=%.14f,\nPy(0)=%.14f,\nB=%.14f,\nDiscrepancy: %.3e\n\n", iteration, mult, alpha_1, alpha_2, integral, sqrt(pow(B[0], 2) + pow(B[1], 2)));

            free (A[0]);
            free (A[1]);
            free (A);
            free (B);
            
            return;
        }
        
        // printf("Shooting iteration: %ud,\nAlpha: %.3f,\nX(0)=%.3e,\nY(0)=%.3e,\nDiscrepancy: %.3e\n\n", iteration, mult,alpha_1,alpha_2,sqrt(pow(B[0],2)+pow(B[1],2)));

        alpha_1 -= (A[0][0]*B[0]-A[1][0]*B[1])/(A[1][0]*A[0][1]-A[0][0]*A[1][1]);
        alpha_2 += (A[0][0]*B[1]-A[1][0]*B[0])/(A[1][0]*A[0][1]-A[0][0]*A[1][1]);
        integral = 0;
    }
    free (A[0]);
    free (A[1]);
    free (A);
    free (B);
}

int main()
{
    double T             = 1;
    double tol           = pow (10, -17);
    long long unsigned i = 0;
    long long unsigned j = 0;

    double x=1;
	double y=0;
	double px=1;
	double py=0;
    double b=0;
	
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

    FILE* inpf = fopen ("koef (8).txt", "r"); //207
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
    alpha[0] = 0;
    alpha[1] = 0;

    k    = (double**)malloc (5 * sizeof (double*));
    k[0] = (double*)malloc (s * sizeof (double));
    k[1] = (double*)malloc (s * sizeof (double));
	k[2] = (double*)malloc (s * sizeof (double));
    k[3] = (double*)malloc (s * sizeof (double));
    k[4] = (double*)malloc (s * sizeof (double));
	
    for (i = 0; i < s; i++)
    {
        k[0][i] = 0;
        k[1][i] = 0;
		k[2][i] = 0;
        k[3][i] = 0;
        k[4][i] = 0;
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
	astep(pi*pow(10,3), &x, &y, &px, &py, &b, &i, &j, p, s, k, cab, tol, test_x_d, test_y_d, test_px_d, test_py_d, test_B_d, 0);
    printf("The Runge-Kutta test:\n%.2e    %.2e  |  %.2e    %.2e  |  %.7e  |  10^3*Pi\n\n", x - 1, y, px - 1, py, b);
    
	// x=1;
	y=0;
	px=0;
	// py=0;
    
    for(unsigned int l = 0; l < 4; l++)
        shooting_method(1500, T, y, px, alpha, p, s, k, cab, tol, x_d, y_d, px_d, py_d, B_d, mult[l]);

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
