#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define IA 16807
#define IM 2147483647
#define AM (1.0/IM)
#define IQ 127773
#define IR 2836
#define NTAB 32
#define NDIV (1+(IM-1)/NTAB)
#define EPS 1.2e-7
#define RNMX (1.0-EPS)

float ran1(long *idum);

#define h 0.01
#define N 40000 //iteracoes
#define transiente  38000
#define Neqs 3
//parametros do modelo
#define a 40.0
#define c 33.0
//~ #define b 2.75

FILE *input;
FILE *input1;
FILE *input2;

void integrador(double y[]);
double equacoes(double y[],int );
double t, y[Neqs];
int j;
double b;

#define NR_END 1
#define FREE_ARG char*
void nrerror(char error_text[]);
double *dvector(long nl,long nh);
void free_dvector(double *v, long nl, long nh);

int main(){
	double *newx, *newy, *newz;
    
    newx=dvector(1,N+10);
    newy=dvector(1,N+10);
    newz=dvector(1,N+10);
	
	input=fopen("bif_indo_x.dat","wt"); 
	input1=fopen("bif_indo_y.dat","wt"); 
	input2=fopen("bif_indo_z.dat","wt"); 

	y[0]=1.0;
	y[1]=1.0; 
	y[2]=0.0; 
	for(b=2.5;b<=4.5;b=b+0.0001){
	t=0.0;
	for(j=0;j<=N;j++){
		t=t+h;
		integrador(y);
		
		newx[j]=y[0];
		newy[j]=y[1];
		newz[j]=y[2];
		if(j>transiente){
			//~ fprintf(input,"%.12f %.12f %.12f %.12f %.12f\n", b, t, y[0], y[1],y[2]);
			if(newx[j-1]>newx[j] && newx[j-1]>newx[j-2]){
			fprintf(input,"%.6f %.6f\n", b, newx[j-1]);
			}
			if(newy[j-1]>newy[j] && newy[j-1]>newy[j-2]){
			fprintf(input1,"%.6f %.6f\n", b, newy[j-1]);
			}
			if(newz[j-1]>newz[j] && newz[j-1]>newz[j-2]){
			fprintf(input2,"%.6f %.6f\n", b, newz[j-1]);
			}
			
			}
		}
	}
	fclose(input);
	fclose(input1);
	fclose(input2);
	
	free_dvector(newx,1,N+10);
	free_dvector(newy,1,N+10);
	free_dvector(newz,1,N+10);
	
	return 0;
}

double equacoes(double y[],int l)
 {
  double resp;
  
  if(l==0)//x=y[0]
    resp=a*(-y[0] + y[1]);
  if(l==1)//y=y[1]
    resp=-y[2]*y[0] + (c-a)*y[0] + c*y[1];
  if(l==2)//z=y[2]
    resp=y[0]*y[1] - b*y[2];
  
  return(resp);
 }

void integrador(double y[])
 {
  double t1[Neqs],t2[Neqs],t3[Neqs],k1[Neqs],k2[Neqs],k3[Neqs],k4[Neqs];
  int n1;

  for(n1=0;n1<Neqs;n1++){
	  k1[n1]=h*equacoes(y,n1);
	  t1[n1]=y[n1]+0.5*k1[n1];
	  }
  
  for(n1=0;n1<Neqs;n1++){
	  k2[n1]=h*equacoes(t1,n1);
	  t2[n1]=y[n1]+0.5*k2[n1];
	  }
	  
  for(n1=0;n1<Neqs;n1++){
	  k3[n1]=h*equacoes(t2,n1);
	  t3[n1]=y[n1]+k3[n1];
	  }
	  
  for(n1=0;n1<Neqs;n1++){
	  k4[n1]=h*equacoes(t3,n1);
	  }
	  
  for(n1=0;n1<Neqs;n1++){
	   y[n1]=y[n1]+(k1[n1]+2.0*k2[n1]+2.0*k3[n1]+k4[n1])/6.0; 
	  }  
 }
 
double *dvector(long nl,long nh)
{
   double *v;
   
   v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
   if (!v) nrerror("allocation failure in dvector()");
   return v-nl+NR_END;
}


void free_dvector(double *v, long nl, long nh)
{
   free((FREE_ARG) (v+nl-NR_END));
}

void nrerror(char error_text[])
{
  fprintf(stderr,"Numerical Recipes run-time error...\n");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"...now exiting to system...\n");
  exit(1);
}

/////////////////////////////////////////
float ran1(long *idum)
{
 int j; long k;
 static long iy=0;
 static long iv[NTAB];
 float temp;
 if(*idum<=0 || !iy)
  {
   if(-(*idum)<1) *idum=1;
   else *idum = -(*idum);
   for(j=NTAB+7;j>=0;j--)
      {
       k=(*idum)/IQ;
       *idum=IA*(*idum-k*IQ)-IR*k;
       if(*idum<0) *idum +=IM;
       if(j<NTAB) iv[j]=*idum;
      }
      iy=iv[0];
   }
 k=(*idum)/IQ;
 *idum=IA*(*idum-k*IQ)-IR*k;
 if(*idum<0) *idum += IM;
 j=iy/NDIV;
 iy=iv[j];
 iv[j]=*idum;
 if((temp=AM*iy)>RNMX) return RNMX;
 else return temp;
}
