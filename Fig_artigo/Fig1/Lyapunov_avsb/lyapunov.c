#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define h 0.01    // step
#define N 3
#define NN 12
//~ #define a 40.0
#define c 33.0
//~ #define b 2.666666667
#define transiente 40000
#define tempofinal 120000

FILE *arquivo;

void rk41(double y[]);
double campo1(double y[],int,double b,double a);
void rk42(double y[]);
double campo2(double y[],int,double b,double a);
double y[NN],znorm[N],gsc[N],cum[N],lyap1,lyap2,lyap3,T,t;
int i,j,k,kk;
double a,b;

int main()
 { 
	 double lyap[N],temp;
	int variable;
	char output_filename[200];    
	scanf("%d", &variable);
	sprintf(output_filename,"x_%d.dat",variable);
	arquivo=fopen(output_filename,"wt");
	
	b=2.0+ 0.007*variable;
	
  for(a=38.0;a<=60.0;a=a+0.044)
   {
    y[0]=1.0;
    y[1]=1.0;
    y[2]=0.0;
    t=0.0;  

    while(t<=transiente)
      {
       rk41(y);
       t+=h;       	 
      }

    for(i=N+1;i<=NN;i++) 
       y[i-1]=0.0;
    
    for(i=1;i<=N;i++)
      {
       y[(N+1)*i-1]=1.0; 
       cum[i-1]=0.0;
      }

    lyap1=0.0;
    lyap2=0.0;
    lyap3=0.0;
  
    kk=0;
    T=50.0*h;
    
    while(t<=tempofinal)
      {
       kk++;
       for(j=1;j<=50;j++)
	 {
	  rk42(y);
	  t+=h;
	 }
 	
       znorm[0]=0.0;
     
       for(j=1;j<=N;j++) 
           znorm[0]=znorm[0]+pow(y[N*j],2);
       znorm[0]=sqrt(znorm[0]);

       for(j=1;j<=N;j++) 
          y[N*j]=y[N*j]/znorm[0];
	
       for(j=2;j<=N;j++)
	 {
	  for(k=1;k<=(j-1);k++) 
	    {
	     gsc[k-1]=0.0;
	     for(i=1;i<=N;i++) 
                 gsc[k-1]=gsc[k-1]+y[N*i+j-1]*y[N*i+k-1];
	    } 
	  for(k=1;k<=N;k++)
	     for(i=1;i<=(j-1);i++) 
                 y[N*k+j-1]=y[N*k+j-1]-gsc[i-1]*y[N*k+i-1];
	  znorm[j-1]=0.0;
	    
	  for(k=1;k<=N;k++)  
             znorm[j-1]=znorm[j-1]+pow(y[N*k+j-1],2);
	  znorm[j-1]=sqrt(znorm[j-1]);
	    
	  for(k=1;k<=N;k++) 
             y[N*k+j-1]=y[N*k+j-1]/znorm[j-1];
	 }
				  
       for(k=1;k<=N;k++) 
          cum[k-1]+=log(fabs(znorm[k-1]));
      }
    
    for(i=0;i<N;i++)
		lyap[i]=cum[i]/(kk*T);
	
	 for(i=0;i<N-1;i++)
{
    for(j=i+1;j<N;j++)
     {
         if(lyap[j] > lyap[i])
         {
            temp=lyap[i];
             lyap[i]=lyap[j];
             lyap[j]=temp;
         }
     }
 }

    fprintf(arquivo,"%20.6f %20.6f %20.6f %20.6f %20.6f\n",b,a,lyap[0],lyap[1],lyap[2]);
   }
    
  fclose(arquivo);
  return 0;
 }

void rk41(double y[])
 {
  double t1[N],t2[N],t3[N],k1[N],k2[N],k3[N],k4[N];
  int n1;

  for(n1=0;n1<N;n1++) t1[n1]=y[n1]+0.5*(k1[n1]=h*campo1(y,n1,b,a));
  for(n1=0;n1<N;n1++) t2[n1]=y[n1]+0.5*(k2[n1]=h*campo1(t1,n1,b,a));
  for(n1=0;n1<N;n1++) t3[n1]=y[n1]+(k3[n1]=h*campo1(t2,n1,b,a));
  for(n1=0;n1<N;n1++) k4[n1]=h*campo1(t3,n1,b,a);
  for(n1=0;n1<N;n1++) y[n1]+=(k1[n1]+2.0*k2[n1]+2.0*k3[n1]+k4[n1])/6.0;  
 }

double campo1(double y[],int l,double b, double a)
 {
  double resp;
   
  if(l==0)
    resp=a*(-y[0] + y[1]);
  if(l==1)
    resp=-y[2]*y[0] + (c-a)*y[0] + c*y[1];
  if(l==2)
    resp=y[0]*y[1] - b*y[2];
  
  return(resp);
 }

void rk42(double y[])
 {
  double t1[NN],t2[NN],t3[NN],k1[NN],k2[NN],k3[NN],k4[NN];
  int n1;

  for(n1=0;n1<NN;n1++) t1[n1]=y[n1]+0.5*(k1[n1]=h*campo2(y,n1,b,a));
  for(n1=0;n1<NN;n1++) t2[n1]=y[n1]+0.5*(k2[n1]=h*campo2(t1,n1,b,a));
  for(n1=0;n1<NN;n1++) t3[n1]=y[n1]+(k3[n1]=h*campo2(t2,n1,b,a));
  for(n1=0;n1<NN;n1++) k4[n1]=h*campo2(t3,n1,b,a);
  for(n1=0;n1<NN;n1++) y[n1]+=(k1[n1]+2.0*k2[n1]+2.0*k3[n1]+k4[n1])/6.0;
 }

double campo2(double y[],int l,double b, double a)
 {
  double resp,j[N+5][N+5],matrix[N+5][N+5],soma[N+5][N+5],R[NN+5];
  int i,m,n,cont;

  /////// Jacobiana //////////////

  j[1][1]=-a;
  j[1][2]=a;
  j[1][3]=0.0;   
  j[2][1]=c-a-y[2];
  j[2][2]=c;
  j[2][3]=-y[0];
  j[3][1]=y[1];
  j[3][2]=y[0];
  j[3][3]=-b;

  /////////// Modelo /////////////////////////
  
  if(l==0)
    resp=a*(-y[0] + y[1]);
  if(l==1)
    resp=-y[2]*y[0] + (c-a)*y[0] + c*y[1];
  if(l==2)
    resp=y[0]*y[1] - b*y[2];

  ////////////////////////////////////////////

  cont=N-1;
  for(m=1;m<=N;m++)
    for(n=1;n<=N;n++)
      {
       cont=cont+1;
       matrix[m][n]=y[cont];
      }

  for(m=1;m<=N;m++)
    for(n=1;n<=N;n++)   
       soma[m][n]=0.0; 

  for(m=1;m<=N;m++)
    for(n=1;n<=N;n++)
      for(i=1;i<=N;i++)
         soma[m][n]=soma[m][n]+j[m][i]*matrix[i][n];

  cont=N-1;
  for(m=1;m<=N;m++)
    for(n=1;n<=N;n++) 
      {
       cont=cont+1;	
       R[cont]=soma[m][n];
      }

  for(i=N;i<NN;i++)
    if(l==i)
      resp=R[l];

  return(resp);
 }
