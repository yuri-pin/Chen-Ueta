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
#define transiente  39500
#define transiente2  39000
#define Neqs 3
//parametros do modelo
#define a 40.0
#define c 33.0
//~ #define b 2.75
#define ncis 40

FILE *input;
FILE *ler_cis;
FILE *input_estx;

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
    double delta[Neqs],media[Neqs],M2[Neqs],variancia[Neqs],desvio[Neqs];
    double xmax, xmin;
	double ymax, ymin;
	double zmax, zmin;
	int i,iteracoes,ci;
	
    newx=dvector(1,N+10);
    newy=dvector(1,N+10);
    newz=dvector(1,N+10);
    
    /********Leitura CIs*************/ 
	ler_cis=fopen("cis_yuri.dat","r");
	int linha;
	double id[ncis],x0[ncis],y0[ncis],z0[ncis];
	for(linha=1;linha<=ncis;linha++){ //linhas
		  fscanf(ler_cis,
           "%lf %lf %lf %lf",
           &id[linha],
           &x0[linha],
           &y0[linha],
           &z0[linha]);
		}
	 fclose(ler_cis);
	 /********Fim Leitura CIs*************/
	
	/***abertura dos arquivos***/
	int variable;
	char output_filename[200];      
	char output_filenamex2[200];      
	scanf("%d", &variable);
	sprintf(output_filename,"bif_x_%d.dat",variable);
	input=fopen(output_filename,"wt");
	sprintf(output_filenamex2,"est_x_%d.dat",variable);
	input_estx=fopen(output_filenamex2,"wt");
	
	for(b=3.0;b<=4.0;b=b+0.0001){
	y[0]=x0[variable];
	y[1]=y0[variable]; 
	y[2]=z0[variable];
	
	t=0.0;
	
	iteracoes=0;
	xmax = -1.0e30;
	xmin =  1.0e30;

	ymax = -1.0e30;
	ymin =  1.0e30;

	zmax = -1.0e30;
	zmin =  1.0e30;
	for(i=0;i<Neqs;i++){
		delta[i]=0.0;
		media[i]=0.0;
		M2[i]=0.0;
		}
		
	for(j=0;j<=N;j++){
		t=t+h;
		integrador(y);
		
		newx[j]=y[0];
		newy[j]=y[1];
		newz[j]=y[2];
		//calculo das estatisticas
		if(j>transiente2){ //transiente para calcular as estatisticas
				if(newx[j-1]>newx[j] && newx[j-1]>newx[j-2]){
					fprintf(input,"%.6f %.6f %.6f %.6f %.6f %.6f\n", b, newx[j-1], 
					id[variable],x0[variable],y0[variable],z0[variable]);
				}
				
				if(y[0] > xmax) xmax = y[0];
				if(y[0] < xmin) xmin = y[0];

				if(y[1] > ymax) ymax = y[1];
				if(y[1] < ymin) ymin = y[1];

				if(y[2] > zmax) zmax = y[2];
				if(y[2] < zmin) zmin = y[2];
            
				iteracoes++;
				for(i=0;i<Neqs;i++){
					delta[i]=y[i] - media[i];
					media[i]=media[i] + delta[i]/iteracoes;
					M2[i]=M2[i]+delta[i]*(y[i] - media[i]);
				}				
				}//fim transiente2	
		}//fim j
		for(i=0;i<Neqs;i++){
				variancia[i]=M2[i]/iteracoes;
				desvio[i]=sqrt(variancia[i]);
			}
			fprintf(input_estx,"%d %.6f %.6f %.6f %.6f %.6f %.6f %.6f %.6f\n", variable, b, x0[ci],y0[ci],desvio[0],media[0],xmax,xmin,xmax-xmin);
	}//fim b
	fclose(input);
	fclose(input_estx);
	
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
