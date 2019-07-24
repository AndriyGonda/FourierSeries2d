#include <stdio.h>
#include <math.h>
#include <time.h>
double integral(double ,double ,unsigned int , double (*) (double));
double an(double ,double ,unsigned int , unsigned int ,double (*) (double));
double bn(double ,double ,unsigned int , unsigned int ,double (*) (double));
double f1(double);
double fourier_series(double , double, double, unsigned int , unsigned int,double(*) (double));
int main()
 {
   /*----------------------------------*/
   unsigned int H = 100,integral_it=200; // numbers of harmonics and number for integral iterations
   double a = -3.14,b = 3.14; //interval min and max values
   double num_of_points = 100; 
   /*----------------------------------*/
   FILE *data = fopen("data.txt","w");
   double x = a,dx=(b-a)/num_of_points;
   char s[100];
    printf("The number of harmonics: %d\nThe number of iterations for the integral: %d\n",H,integral_it);
    printf("\t     x:\t\ty=f(x):\t\t  y=f_F(x):\n");
    clock_t time = clock();
    while(x<=b)
     {
        sprintf(s,"\t%12.6f\t%12.6f\t%12.6f\n",x,f1(x),fourier_series(a,b,x,H,integral_it,f1));
        puts(s);
        fputs(s,data);
       x+=dx;
     }
     time = clock() - time;
    fclose(data);
    printf("Calculations time : %f sec \n", (double)time/CLOCKS_PER_SEC);
    FILE *gr  = popen("gnuplot -persist", "w");
    fprintf(gr,"plot 'data.txt' using 1:2 w p lc 3 title 'Start function ', 'data.txt' using 1:3 title 'Fourier function'  w l lc 2\n");
    fclose(gr);
   return 0;
 }
 
 //Sympson integral method
 double integral(double a,double b ,unsigned int iteration , double (*f) (double x))
  {
    double h = (b-a)/(double )iteration;
    double sum = 0;
    double x0 = a;
    double x1 = a+h;
    for (register unsigned int i = 0; i<iteration; i++)
     {
       sum += f(x0) + 4*f(x0+h/2.) + f(x1);
       x0+=h;
       x1+=h;
     }
    return sum*h/6.;
  }

// Start function
 double f1(double x)
  {
    return x;
  }

// An 
 double an(double a,double b,unsigned int n , unsigned int iteration,double (*f) (double x))
  {
     double h = (b-a)/(double )iteration;
    double sum = 0;
    double x0 = a;
    double x1 = a+h;
    for (unsigned int i = 0; i<iteration; i++)
     {
       sum += f(x0)*cos((2.*3.14*n*x0)/(b-a)) + 4*f(x0+h/2.)*cos(2.*3.14*n*(x0+h/2.)/(b-a)) + f(x1)*cos(2.*3.14*n*x1/(b-a));
       x0+=h;
       x1+=h;
     }
    return (2./(b-a)) * sum*h/6.;
  }

// Bn 
 double bn(double a,double b,unsigned int n , unsigned int iteration,double (*f) (double x))
  {
     double h = (b-a)/(double)iteration;
    double sum = 0;
    double x0 = a;
    double x1 = a+h;
    for (unsigned int i = 0; i<iteration; i++)
     {
       sum += f(x0)*sin(2.*3.14*n*x0/(b-a)) + 4*f(x0+h/2.)*sin(2.*3.14*n*(x0+h/2.)/(b-a)) + f(x1)*sin(2.*3.14*n*x1/(b-a));
       x0+=h;
       x1+=h;
     }
    return (2./(b-a)) * sum*h/6.;
  }

 double  fourier_series(double a, double b, double x,unsigned int series_size, unsigned int intgr_iteration,double(*f)(double))
  {
    double interval = b - a;
     double a0 = 2./interval * integral (a , b, intgr_iteration,f);
        double sum = 0;
        for(unsigned int n =1 ; n<=series_size; n++)
         {
           sum += an(a,b,n,intgr_iteration,f)*cos((2*3.14*n*x)/interval) + bn(a,b,n,intgr_iteration,f)*sin((2*3.14*n*x)/interval) ;
         }
    return a0/2.+sum;
  }
 
 // Created by Andriy Gonda 
