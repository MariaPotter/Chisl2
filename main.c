#include <math.h>
#include <stdbool.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#define pi 3.14159265358979323846264338328
#define EPS pow(10,-17)
_Pragma ("GCC diagnostic push")
_Pragma ("GCC diagnostic ignored \"-Wunused-parameter\"")
static double x_d  (double t, double x, double y, double p_x, double p_y, double mult, double contol(double,double,double,double,double,double))        
{
    return y;
}

static double y_d (double t, double x, double y, double p_x, double p_y, double mult, double contol(double,double,double,double,double,double))
{
    return p_y/2-sqrt(2.)*x*exp(-mult*t);
}

static double px_d(double t, double x, double y, double p_x, double p_y, double mult, double contol(double,double,double,double,double,double))
{
    return 2*x+sqrt(2.)*p_y*exp(-mult*t);
}

static double py_d (double t, double x, double y, double p_x, double p_y, double mult, double contol(double,double,double,double,double,double))
{
    return -p_x;
}

static double B_d  (double t, double x, double y, double p_x, double p_y, double mult, double contol(double,double,double,double,double,double))        
{
    return pow(x,2.) + pow(p_y,2.)/4.;
}

static double test_x_d  (double t, double x, double y, double p_x, double p_y, double mult, double contol(double,double,double,double,double,double))        
{
    return -y;
}

static double test_y_d (double t, double x, double y, double p_x, double p_y, double mult, double contol(double,double,double,double,double,double))
{
    return x;
}

static double test_px_d(double t, double x, double y, double p_x, double p_y, double mult, double contol(double,double,double,double,double,double))
{
    return -p_y;
}

static double test_py_d (double t, double x, double y, double p_x, double p_y, double mult, double contol(double,double,double,double,double,double))
{
    return p_x;
}

static double test_B_d  (double t, double x, double y, double p_x, double p_y, double mult, double contol(double,double,double,double,double,double))        
{
    return sqrt(pow(x-p_x,2) + pow(y-p_y,2));
}

static double contol_n (double t, double x, double y, double p_x, double p_y, double mult)
{
    if(p_y > 0)
        return (p_y/cos(mult*x)<24)?p_y/cos(mult*x):24;
    else
        return (p_y/cos(mult*x)>-24)?p_y/cos(mult*x):-24;
}
/*
static double contol_p (double t, double x, double y, double p_x, double p_y, double mult)
{
    return 24;
}

static double contol_m (double t, double x, double y, double p_x, double p_y, double mult)
{
    return -24;
}
*/
_Pragma ("GCC diagnostic pop")
    
    /*
void diagonal() //finding crossing point
{    
    double lx = x, ly = y, lpx = px, lpy = py, lb = b, ldist = *dist;
    double _x = x, _y = y, _dist = *dist, _h = temp;
                
    for (; fabs (_x - *startx) > EPSM;)        // finding point of crossing
    {
        if ((_x - *startx) * (tempx - *startx) < 0)
        {
            lx = _x;
            ly = _y;
            ldist = _dist;
            _h = fabs (cos(mult*_x)) / fabs (cos(mult*tempx) - cos(mult*_x)) * (*dist + temp - _dist);        // step from _x to startx

            RK (_dist, &_x, &_y, _h, s, k, cab, f, g, false);
            _dist += _h;
        }
        else
        {
            tempx = _x;
            tempy = _y;
            temp = _dist - *dist;

            _h = fabs (cos(mult*lx)) / fabs (cos(mult*lx) - cos(mult*_x)) * (_dist - ldist);        // step from lx to startx

            _x    = lx;
            _y    = ly;
            _dist = ldist;

            RK (_dist, &_x, &_y, _h, s, k, cab, f, g, false);

            _dist += _h;
        }
    }

    tempx = x;
    tempy = y;
    temp  = _dist - *dist;

    RK (*dist, &tempx, &tempy, temp, s, k, cab, f, g, false);        // counting starty point
    RK (dist, &temp_x, &temp_y, &temp_px, &temp_py, &temp_b, h, s, k, cab, f, g, u, v, l, mult, *c, false);
}
    */
    
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
         double f (double, double, double, double, double, double, double (*)(double,  double,  double,  double,  double,  double)),
         double g (double, double, double, double, double, double, double (*)(double,  double,  double,  double,  double,  double)),
		 double u (double, double, double, double, double, double, double (*)(double,  double,  double,  double,  double,  double)),
         double v (double, double, double, double, double, double, double (*)(double,  double,  double,  double,  double,  double)),
         double l (double, double, double, double, double, double, double (*)(double,  double,  double,  double,  double,  double)),
         double mult,
         double c (double, double, double, double, double, double),
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
                    mult, c);
					
        k[1][i] = g (t + h * cab[0][i],
					*x0 + h * temp_x,
					*y0 + h * temp_y,
					*px0 + h * temp_px,
					*py0 + h * temp_py,
                    mult, c);
					
		k[2][i] = u (t + h * cab[0][i],
					*x0 + h * temp_x,
					*y0 + h * temp_y,
					*px0 + h * temp_px,
					*py0 + h * temp_py,
                    mult, c);
					
		k[3][i] = v (t + h * cab[0][i],
					*x0 + h * temp_x,
					*y0 + h * temp_y,
					*px0 + h * temp_px,
					*py0 + h * temp_py,
                    mult, c);
        
        k[4][i] = l (t + h * cab[0][i],
					*x0 + h * temp_x,
					*y0 + h * temp_y,
					*px0 + h * temp_px,
					*py0 + h * temp_py,
                    mult, c);
        
        
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
			  
              double f (double, double, double, double, double, double, double (*)(double,  double,  double,  double,  double,  double)),
			  double g (double, double, double, double, double, double, double (*)(double,  double,  double,  double,  double,  double)),
			  double u (double, double, double, double, double, double, double (*)(double,  double,  double,  double,  double,  double)),
			  double v (double, double, double, double, double, double, double (*)(double,  double,  double,  double,  double,  double)),
			  double l (double, double, double, double, double, double, double (*)(double,  double,  double,  double,  double,  double)),
              double mult,
              bool test)
{
    double temp_x, temp_y, temp_px, temp_py, temp_b, x_, y_, px_, py_, b_, dist, h, temp, fac, err, norm;
    double (*c) (double, double, double, double, double, double);
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
	c = (test)?0:contol_n;
    // printf("%.2e    %.6e  |  %.6e    %.2e\n", *x, *y, *px, *py);
	
    for (*i = *j = 0; T - dist > EPS;)
    {   
        
        if (dist + h > T) h = T - dist;

        RK (dist, &temp_x, &temp_y, &temp_px, &temp_py, &temp_b, h, s, k, cab, f, g, u, v, l, mult, *c, false);
        RK (dist, &x_, &y_, &px_, &py_, &b_, h, s, k, cab, f, g, u, v, l, mult, *c, true);
        
        
        // norm = fmax( fmax( fabs (temp_x - x_), fabs (temp_y - y_)), fmax(fabs (temp_px - px_), fabs (temp_py - py_)));
        norm = sqrt(pow(temp_x - x_,2)+pow(temp_y - y_,2)+pow(temp_px - px_,2)+pow(temp_py - py_,2));
        temp = h;
        h *= fmin (fac,
                   fmax (0.7,
                         pow (0.98 * tol / norm,
                              1. / (p + 1))));
        if (h < pow (10, -18))
        {
            // printf ("\nSomething goes wrong...\n");
            return -1000000;
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

	        double f (double, double, double, double, double, double, double (*)(double,  double,  double,  double,  double,  double)),
	        double g (double, double, double, double, double, double, double (*)(double,  double,  double,  double,  double,  double)),
	        double u (double, double, double, double, double, double, double (*)(double,  double,  double,  double,  double,  double)),
	        double v (double, double, double, double, double, double, double (*)(double,  double,  double,  double,  double,  double)),
	        double l (double, double, double, double, double, double, double (*)(double,  double,  double,  double,  double,  double)),
            double mult)
{
	long long unsigned i = 0;
	long long unsigned j = 0;

	if(astep(T, &x, &y, &px, &py, &b, &i, &j, p, s, k, cab, tol, f, g, u, v, l, mult, 0) + 1000000 < EPS)
        return -1000000;

	err[0] = px-0;//
	err[1] = py-0;
    
    return b;
}

// the shooting method
int shooting_method(unsigned int max_iterations,
                     double T,
	                 double x,
	                 double py,
	                 double* alpha,
	                 unsigned int p,
	                 unsigned int s,

	                 double** k,
	                 double** cab,

	                 double tol,

	                 double f (double, double, double, double, double, double, double (*)(double,  double,  double,  double,  double,  double)),
	                 double g (double, double, double, double, double, double, double (*)(double,  double,  double,  double,  double,  double)),
	                 double u (double, double, double, double, double, double, double (*)(double,  double,  double,  double,  double,  double)),
	                 double v (double, double, double, double, double, double, double (*)(double,  double,  double,  double,  double,  double)),
	                 double l (double, double, double, double, double, double, double (*)(double,  double,  double,  double,  double,  double)),
                     double mult)
{
    double alpha_1; //x
    double alpha_2; //py
    double betta;   //dempher
    unsigned int iteration;
    double** A = (double**)malloc (2*sizeof(double*));
    A[0] = (double*)malloc (2*sizeof(double));
    A[0][0] = 0;
    A[0][1] = 0;
    A[1] = (double*)malloc (2*sizeof(double));
    A[1][0] = 0;
    A[1][1] = 0;
    double* B = (double*)malloc (6*sizeof(double));
    for(int i = 0; i < 6; i++) B[i] = 0;
    double* newton = (double*)malloc (2*sizeof(double));
    newton[0] = 0;
    newton[1] = 0;
	double* err = (double*)malloc (2*sizeof(double));
    err[0] = 0;
    err[1] = 0;
    double delta = pow(2, -11);
    double step_m = pow(10, -1);
    double step_T = pow(10, -1);
    double integral;
    for(;1;){
        integral = 0;
        alpha_1 = alpha[2]; //y
        alpha_2 = alpha[3]; //px
         // printf("mult_0=%13.7e  |  mult_k=%13.7e  |  mult=%13.7e  |  T_0=%13.7e  |  T_k=%13.7e  |  T=%13.7e  |  step_m=%13.7e  |  step_T=%13.7e\n",alpha[0], alpha[4], mult, alpha[1], alpha[5], T, step_m, step_T);

        for(iteration = 0; ; iteration++){
            // printf("T=%14.7e  |  Y(0)=%14.7e  |  Px(0)=%14.7e  | Mult=%14.7e\n",alpha[5], alpha_1, alpha_2, alpha[4]);
            if(iteration == max_iterations)
            {
                if(fabs(alpha[4] - (alpha[0] + step_m)) < EPS)
                    step_m *= 0.1;
                if(fabs(alpha[5] - (alpha[1] + step_T)) < EPS)
                    step_T *= 0.1;
                break;
            }
            if(iteration % 5 == 0)
            {
                //dy
                
                if(error(alpha[5], x, alpha_1+delta, alpha_2, py, integral, err, p, s, k, cab, tol, f, g, u, v, l, alpha[4]) + 1000000 < EPS){
                    if(fabs(alpha[4] - (alpha[0] + step_m)) < EPS)
                        step_m *= 0.1;
                    if(fabs(alpha[5] - (alpha[1] + step_T)) < EPS)
                        step_T *= 0.1;
                    break;
                }

                B[0] = err[0];  // x
                B[3] = err[1];  // py

                if(error(alpha[5], x, alpha_1-delta, alpha_2, py, integral, err, p, s, k, cab, tol, f, g, u, v, l, alpha[4]) + 1000000 < EPS){
                    if(fabs(alpha[4] - (alpha[0] + step_m)) < EPS)
                        step_m *= 0.1;
                    if(fabs(alpha[5] - (alpha[1] + step_T)) < EPS)
                        step_T *= 0.1;
                    break;
                }
                
                B[1] = err[0];  // x
                B[4] = err[1];  // py

                if(error(alpha[5], x, alpha_1+2*delta, alpha_2, py, integral, err, p, s, k, cab, tol, f, g, u, v, l, alpha[4]) + 1000000 < EPS){
                    if(fabs(alpha[4] - (alpha[0] + step_m)) < EPS)
                        step_m *= 0.1;
                    if(fabs(alpha[5] - (alpha[1] + step_T)) < EPS)
                        step_T *= 0.1;
                    break;
                }

                B[2] = err[0];  // x
                B[5] = err[1];  // py

                if(error(alpha[5], x, alpha_1-2*delta, alpha_2, py, integral, err, p, s, k, cab, tol, f, g, u, v, l, alpha[4]) + 1000000 < EPS){
                    if(fabs(alpha[4] - (alpha[0] + step_m)) < EPS)
                        step_m *= 0.1;
                    if(fabs(alpha[5] - (alpha[1] + step_T)) < EPS)
                        step_T *= 0.1;
                    break;
                }
                
                A[0][0] = (2.*B[0]-2.*B[1]-B[2]/4.+err[0]/4.)/(3*delta);
                A[1][0] = (2.*B[3]-2.*B[4]-B[5]/4.+err[1]/4.)/(3*delta);

                // dpx
                
                if(error(alpha[5], x, alpha_1, alpha_2+delta, py, integral, err, p, s, k, cab, tol, f, g, u, v, l, alpha[4]) + 1000000 < EPS){
                    if(fabs(alpha[4] - (alpha[0] + step_m)) < EPS)
                        step_m *= 0.1;
                    if(fabs(alpha[5] - (alpha[1] + step_T)) < EPS)
                        step_T *= 0.1;
                    break;
                }

                B[0] = err[0];  // x
                B[3] = err[1];  // py

                if(error(alpha[5], x, alpha_1, alpha_2-delta, py, integral, err, p, s, k, cab, tol, f, g, u, v, l, alpha[4]) + 1000000 < EPS){
                    if(fabs(alpha[4] - (alpha[0] + step_m)) < EPS)
                        step_m *= 0.1;
                    if(fabs(alpha[5] - (alpha[1] + step_T)) < EPS)
                        step_T *= 0.1;
                    break;
                }

                B[1] = err[0];  // x
                B[4] = err[1];  // py

                if(error(alpha[5], x, alpha_1, alpha_2+2.*delta, py, integral, err, p, s, k, cab, tol, f, g, u, v, l, alpha[4]) + 1000000 < EPS){
                    if(fabs(alpha[4] - (alpha[0] + step_m)) < EPS)
                        step_m *= 0.1;
                    if(fabs(alpha[5] - (alpha[1] + step_T)) < EPS)
                        step_T *= 0.1;
                    break;
                }

                B[2] = err[0];  // x
                B[5] = err[1];  // py

                if(error(alpha[5], x, alpha_1, alpha_2-2.*delta, py, integral, err, p, s, k, cab, tol, f, g, u, v, l, alpha[4]) + 1000000 < EPS){
                    if(fabs(alpha[4] - (alpha[0] + step_m)) < EPS)
                        step_m *= 0.1;
                    if(fabs(alpha[5] - (alpha[1] + step_T)) < EPS)
                        step_T *= 0.1;
                    break;
                }

                A[0][1] = (2.*B[0]-2.*B[1]-B[2]/4.+err[0]/4.)/(3*delta);
                A[1][1] = (2.*B[3]-2.*B[4]-B[5]/4.+err[1]/4.)/(3*delta);
            }
            
            integral = error(alpha[5], x, alpha_1, alpha_2, py, integral, err, p, s, k, cab, tol, f, g, u, v, l, alpha[4]);
            if(integral + 1000000 < EPS){
                if(fabs(alpha[4] - (alpha[0] + step_m)) < EPS)
                    step_m *= 0.1;
                if(fabs(alpha[5] - (alpha[1] + step_T)) < EPS)
                    step_T *= 0.1;
                break;
            }
            
            B[0] = err[0];  //x
            B[1] = err[1];  //py

            // printf("%10.2e    %10.2e  |  %10.2e\n%10.2e    %10.2e  |  %10.2e\n", A[0][0], A[0][1], B[0], A[1][0], A[1][1], B[1]);
            
            B[5] = sqrt(pow(B[0],2)+pow(B[1],2));
            if(B[5] < tol*pow(10, 3)){
                // printf and break
                // printf("Shooting iteration: %u,\nAlpha: %.4f,\nDistance: %.4f,\nY(0) = %.14f,\nPx(0) = %.14f,\nB = %.14f,\nDiscrepancy: %.3e\n\n", iteration, alpha[4], alpha[5], alpha_1, alpha_2, integral, sqrt(pow(B[0], 2) + pow(B[1], 2)));
                
                alpha[0] = alpha[4];
                alpha[1] = alpha[5];
                alpha[2] = alpha_1;
                alpha[3] = alpha_2;

                break;
            }

            // printf("Shooting iteration: %ud,\nAlpha: %.3f,\nY(0) = %.3e,\nPx(0) = %.3e,\nDiscrepancy: %.3e\n\n", iteration, mult,alpha_1,alpha_2,sqrt(pow(B[0],2)+pow(B[1],2)));
            
            newton[0] = (A[0][0]*B[0]-A[1][0]*B[1])/(A[1][0]*A[0][1]-A[0][0]*A[1][1]);
            newton[1] = (A[0][0]*B[1]-A[1][0]*B[0])/(A[1][0]*A[0][1]-A[0][0]*A[1][1]);
            
            //demph
            
            for(betta = 1;;){
                if(error(alpha[5], x, alpha_1 - betta * newton[0], alpha_2 + betta * newton[1], py, integral, err, p, s, k, cab, tol, f, g, u, v, l, alpha[4]) + 1000000 < EPS){
                    betta *= 1.01;
                    break;
                }
                if(sqrt(pow(err[0],2)+pow(err[1],2))-B[5]<EPS && betta > pow(2, -15))
                {
                    B[5] = sqrt(pow(err[0],2)+pow(err[1],2));
                    betta /= 1.01;
                }else{
                    betta *= 1.01;
                    break;
                }
            }
            
            
            alpha_1 -= betta * newton[0];
            alpha_2 += betta * newton[1];
            integral = 0;
        }
        
        if(fabs(alpha[0] - alpha[4]) < EPS && fabs(alpha[1] - alpha[5]) < EPS )
        {
            step_m *= 1.1;
            step_T *= 1.1;
            if(mult - alpha[4] < step_m)
                step_m = mult - alpha[4];
            
            if(T - alpha[5] < step_T)
                step_T = T - alpha[5];
        }
        
        if(fabs(mult - alpha[0]) > EPS)
            alpha[4] = alpha[0] + step_m;
        if(fabs(T - alpha[1]) > EPS)
            alpha[5] = alpha[1] + step_T;
        
        
        if(fabs(alpha[0] - mult) < EPS && fabs(alpha[1] - T) < EPS) {
            printf("Shooting iteration: %u,\nAlpha: %.4f,\nDistance: %.4f,\nY(0) = %.14f,\nPx(0) = %.14f,\nB = %.14f,\nDiscrepancy: %.3e\n\n", iteration, alpha[4], alpha[5], alpha_1, alpha_2, integral, sqrt(pow(B[0], 2) + pow(B[1], 2)));
            break;
        }
        
        if((fabs(alpha[5] - alpha[1]) < pow(10, -4) && fabs(alpha[5] - alpha[1]) > EPS) || (fabs(alpha[4] - alpha[0]) < pow(10, -4) && fabs(alpha[4] - alpha[0]) > EPS)) {
                printf("Newton's method does not converge!\n");
                integral = 0;
                integral = error(alpha[1], x, alpha[2], alpha[3], py, integral, err, p, s, k, cab, tol, f, g, u, v, l, alpha[0]);
                printf("Shooting iteration: %u,\nAlpha: %16.4f,\nRequested alpha: %.4f,\nDistance: %16.4f,\nRequested distance: %.4f,\nY(0) = %.14f,\nPx(0) = %.14f,\nB = %.14f,\nDiscrepancy: %.3e\nSteps: %.4e   %.4e\n\n", iteration, alpha[0], mult, alpha[1], T, alpha[2], alpha[3], integral, sqrt(pow(err[0], 2) + pow(err[1], 2)), step_m, step_T);
    
                free (A[0]);
                free (A[1]);
                free (A);
                free (B);
                free (err);
                
                return -1;
        }
        iteration = 0;
    }
    
    free (A[0]);
    free (A[1]);
    free (A);
    free (B);
    free (err);
    
    return 0;
}

int main()
{
    double tol           = pow (10, -17);
    long long unsigned i = 0;
    long long unsigned j = 0;

    double x=1;
	double y=0;
	double px=1;
	double py=0;
    double b=0;
	
    double*     T; // T из условия
    unsigned int T_c = 2;
    double*  mult; // alpha из условия
    unsigned int mult_c = 8;
    double* alpha; // под невязки

    unsigned int s;
    unsigned int p;
    double** k;
    double** cab;

    time_t start, end;

    /////////////////////////////////////////////////////////////////////////////////////////

    // file check


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
    
    T = (double*)malloc (T_c * sizeof (double));
    for(i=0; i < T_c; i++){
        switch (i) {
            case 0:
                T[i] = 1;
                break;
            case 1:
                T[i] = 3;
                break;
            case 2:
                T[i] = 3;
                break;
            case 3:
                T[i] = 3.5;
                break;
            case 4:
                T[i] = 10;
                break;
            case 5:
                T[i] = 0.1;
                break;
            default:
                T[i] = 1;
        }
    }
    
    mult = (double*)malloc (mult_c * sizeof (double));
	/*
    for(i=0; i < mult_c; i++){
        switch (i) {
            case 0:
                mult[i] = 0.0;
                break;
            default:
                mult[i] = (i<=4)?pow(10,(double)i - 4.):5.*((double)i - 4.);
        }
    }
	*/
    mult[0] = 0.0;
    mult[1] = 0.01;
    mult[2] = 0.1;
    mult[3] = 1;
    mult[4] = 2;
    mult[5] = 5;
    mult[6] = 10;
    mult[7] = 50;

	alpha = (double*)malloc (6 * sizeof (double));
    alpha[0] = mult[0];
    alpha[1] = T[0];
    alpha[2] = 5.;
    alpha[3] = -12.;

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
	astep(pi*pow(10,3), &x, &y, &px, &py, &b, &i, &j, p, s, k, cab, tol, test_x_d, test_y_d, test_px_d, test_py_d, test_B_d, 0, 1);
    printf("The Runge-Kutta test:\n%.2e    %.2e  |  %.2e    %.2e  |  %.7e  |  10^3*Pi\n\n", x - 1, y, px - 1, py, b);
    
    time (&start);
    
	x=1;
	y=0;
	px=0;
	py=0;
    
    for(unsigned int m = 0; m < T_c; m++)
        for(unsigned int l = 0; l < mult_c; l++){
            alpha[4] = mult[l];
            alpha[5] = T[m];
        
            if(shooting_method(500, T[m], x, py, alpha, p, s, k, cab, tol, x_d, y_d, px_d, py_d, B_d, mult[l]) == -1){  
                alpha[4] = mult[l+1];
                alpha[5] = T[m];
                if(shooting_method(500, T[m], x, py, alpha, p, s, k, cab, tol, x_d, y_d, px_d, py_d, B_d, mult[l+1]) == -1)
                    break;
            }
        }

    /////////////////////////////////////////////////////////////////////////////////////////

    // cleaning :)

    fclose (inpf);
    for (i = 0; i < s + 3; i++)
    {
        free (cab[i]);
    }
    free (cab);
	free (T);
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
