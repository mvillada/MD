/******************************************************************************/
/* Dinámica molecular velocity-Verlet 1 especie                               */
/* Estructura FCC                                                             */
/* Potencial de Morse                                                         */
/* Versión 1                                                                  */
/*                                                                            */
/* Calcula:                                                                   */
/* g(r)  Función de distribución radial                                       */
/*                                                                            */
/* Parámetros de simulación                                                   */
/* Pi           Valor del número pi                                           */
/* Ncx          Número de celdas en x                                         */
/* Ncy          Número de celdas en y                                         */
/* Ncz          Número de celdas en z                                         */
/* frand()      Generador de números aleatorios                               */
/* maxint       Número máximo para el arreglo de la g(r)                      */
/* Temp         Temperatura reducida                                          */
/* dt           Paso de integración                                           */
/* dr           Intervalo de muestreo                                         */
/* intstep      Número de pasos de integración                                */
/* ncg          Frecuencia de cálculo de la g(r)                              */
/* rprint       Frecuencia de impresión                                       */
/* rprintf      Frecuencia de impresión en el archivo                         */
/* a            Parametro de red                                              */
/* b            Parametro de red                                              */
/* c            Parametro de red                                              */
/* U0           Fondo del pozo de potencial                                   */
/* alpha        Parámetro del potencial de Morse                              */
/* Rij_0        Distancia del mínimo                                          */
/* epsilon      Profundidad del pozo                                          */
/* Npt          Número de partículas                                          */
/* intstepinit  Número para iniciar la estadística                            */
/******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#define Pi		3.14159265358979
#define Ncx		4
#define Ncy		4
#define Ncz		4
#define frand()		((double) rand() / (RAND_MAX + 1.0)) 
#define maxint		1000000                       
#define Temp		1.0
#define dt		0.005
#define dr		0.005
#define intstep		1000
#define intstepinit	100
#define ncg		10
#define rprint		5000
#define rprintf		10
#define a		cbrt(4.0/0.8)
#define b		cbrt(4.0/0.8)
#define c		cbrt(4.0/0.8)
#define U0		0.1
#define alpha		20.0
#define Rij_0		1.0
#define sigma		1.0                         
#define Npt		(4*Ncx*Ncy*Ncz)
#define dk		0.05
#define mass		1.0

//Variables de simulación
int Npi,i,j,k,k1,k2,iRij,isample,icel,istep,iNBin,NBin_ad;
double Rx[Npt+1],Ry[Npt+1],Rz[Npt+1],Xi,Yi,Zi,Rijx,Rijy,Rijz,Rij,U,IU;
double Lboxx,Lboxy,Lboxz,Vol,pack,Vx[Npt+1],Vy[Npt+1],Vz[Npt+1];
double KEn,TIKEn,DU,Fijx,Fijy,Fijz,denssd,Fx[Npt+1],Fy[Npt+1],Fz[Npt+1];
double Virx,Viry,Virz,TMass,TMassI,Momx,Momy,Momz,IKEn;
double KENAXS,KENAYS,KENAZS,KENAX,KENAY,KENAZ,Lmin,c_v,PITEn2,ITEn2,TITEn2;
double IpressP,LboxxI,LboxyI,LboxzI,ITemp,Tratio;
double IpressK,IPressT,ITEn,Ttemp,TIU,TIPressT,TITEn,Ptemp,PIKEn,PIU,PIPressT,PITEn;
double NBin,Dens,Dvol,R1,R2,DifC_R;
double RxRe[Npt+1],RyRe[Npt+1],RzRe[Npt+1];
double RxAf[Npt+1],RyAf[Npt+1],RzAf[Npt+1];
double Mass[Npt+1],IMass[Npt+1];
double FxBe[Npt+1],FyBe[Npt+1],FzBe[Npt+1];
double Hist[maxint],Rdf1[maxint],R[maxint];
double intstat,Vx2,Vy2,Vz2;
const double dt2 = dt*dt;
double DRij;

//Funciones
void Init(void);
float grandom(void);
void Force(void);
void Instant(void);
void Prom(void);
void Zero(void);
void Sample(void);
void Rdf(void);
void Move1(void);
void Move2(void);

/*********************Inicia programa principal*******************************/
int main(){

srand ( time(NULL) );

puts("\n\t\t********************");
puts("\t\t INICIANDO PROGRAMA...");
puts("\t\t********************");

FILE *Stat;
Stat=fopen("StatD0844T071.dat","w");

Zero();

Init();

//* Película
FILE *Movie;
Movie=fopen("filmD0844T071.xyz","w");
fprintf(Movie,"%d\n",Npt);
fprintf(Movie,"FRAME 1\n");
for(i=1;i<=Npt;i++){
fprintf(Movie,"C %f %f %f \n",Rx[i],Ry[i],Rz[i]);
}
//Película */

Force();

for(istep=1;istep<=intstep;istep++){
Move1();

//* Película
fprintf(Movie,"%d\n",Npt);
fprintf(Movie,"FRAME %d\n",istep+1);
for(i=1;i<=Npt;i++){
fprintf(Movie,"C %f %f %f \n",Rx[i],Ry[i],Rz[i]);
}
//Película */

Force();
Move2();
Instant();

if (istep%rprintf==0){
fprintf(Stat," %f \t %f \t %f \t %f \t %f \t %f \n",istep*dt,IKEn,IU,ITEn,ITemp,IPressT);
}

if (istep>intstepinit && istep%ncg==0) Sample();

}

Rdf();
Prom();

  puts("\n\t\t*********************");
  puts("\t \t PROGRAMA FINALIZADO");
  puts("\t\t*********************");

fclose(Stat);
fclose(Movie);

return 0;
}
/*********************Termina programa principal******************************/


//Configuración inicial
void Init(void){

FILE *data;
data=fopen("dataD0844T071.xyz","w");

TMass  = (double)Npt;
TMassI = 1.0/TMass;
Momx   = 0.0;
Momy   = 0.0;
Momz   = 0.0;
Npi    = 1;

//Longitudes de la caja, volumen y densidad
Lboxx  = a*(double)Ncx;
Lboxy  = b*(double)Ncy;
Lboxz  = c*(double)Ncz;
LboxxI = 1.0/Lboxx;
LboxyI = 1.0/Lboxy;
LboxzI = 1.0/Lboxz;
Vol    = Lboxx*Lboxy*Lboxz;
Dens   = (double)Npt/Vol;
pack   = Pi*Dens/6.0;

//Posiciones de la celda unitaria
Rx[1] = 0.0-0.5*Lboxx;
Ry[1] = 0.0-0.5*Lboxy;
Rz[1] = 0.0-0.5*Lboxz;

Rx[2] = a/2.0-0.5*Lboxx;
Ry[2] = b/2.0-0.5*Lboxy;
Rz[2] = 0.0-0.5*Lboxz;

Rx[3] = 0.0-0.5*Lboxx;
Ry[3] = b/2.0-0.5*Lboxy;
Rz[3] = c/2.0-0.5*Lboxz;

Rx[4] = a/2.0-0.5*Lboxx;
Ry[4] = 0.0-0.5*Lboxy;
Rz[4] = c/2.0-0.5*Lboxz;

//Crear posiciones de las partículas
for(i=1;i<=Ncx;i++){
 for(j=1;j<=Ncy;j++){
  for(k=1;k<=Ncz;k++){
   for(icel=0;icel<=3;icel++){
   Rx[Npi+icel]=Rx[icel+1]+a*(double)(i-1);
   Ry[Npi+icel]=Ry[icel+1]+b*(double)(j-1);
   Rz[Npi+icel]=Rz[icel+1]+c*(double)(k-1);
   Mass[Npi+icel]=mass;
   IMass[Npi+icel]=1.0/mass;
}
Npi=Npi+4;
}}}

//Condiciones de frontera periodica
for (i=1;i<=Npt;i++){
Rx[i] = Rx[i]-rint(Rx[i]*LboxxI)*Lboxx;
Ry[i] = Ry[i]-rint(Ry[i]*LboxyI)*Lboxy;
Rz[i] = Rz[i]-rint(Rz[i]*LboxzI)*Lboxz;
}

//Archivo de posiciones iniciales
fprintf(data,"%d \n",Npt);
fprintf(data,"FRAME 1\n");

for(i=1;i<=Npt;i++){
fprintf(data,"C %f %f %f\n",Rx[i],Ry[i],Rz[i]);
}

fclose(data);

//Calcula longitud mínima para calcular la g(r)
Lmin     = Lboxx;
NBin    = Lmin/(2.0*dr);
iNBin   = rint(NBin);
NBin_ad = iNBin;

if(iNBin%2==1) NBin_ad=NBin_ad+1;

  puts("\n\t\t********************");
  puts("\t\t INICIANDO PROGRAMA...");
  puts("\t\t********************");

printf("\ta            %f \n",a);
printf("\tb            %f \n",b);
printf("\tc            %f \n",c);
printf("\tTemp         %f \n",Temp);
printf("\tmass         %f \n",mass);
printf("\tTMass        %f \n",TMass);
printf("\tTMassI       %f \n",TMassI);
printf("\tLboxx        %f \n",Lboxx);
printf("\tLboxy        %f \n",Lboxy);
printf("\tLboxz        %f \n",Lboxz);
printf("\tVol          %f \n",Vol);
printf("\tDens         %f \n",Dens);
printf("\tpack         %f \n",pack);
printf("\tdt           %f \n",dt);
printf("\tdr           %f \n",dr);
printf("\tdk           %f \n",dk);
printf("\tintstep      %d \n",intstep);
printf("\tintstepinit  %d \n",intstepinit);
printf("\tLboxxI       %f \n",LboxxI);
printf("\tLboxyI       %f \n",LboxyI);
printf("\tLboxzI       %f \n",LboxzI);
printf("\tLmin         %f \n",Lmin);
printf("\tNBin        %f \n",NBin);
printf("\tiNBin       %d \n",iNBin);
  
  puts("\n\t\t********************");
  puts("\t\t INICIALIZACION...");
  puts("\t\t********************");
  printf("\tNUMERO DE PARTICULAS        %d \n", Npt);
  puts("\n\t\t POSICIONES INICIALIZADAS...\n");


for(i=1;i<=Npt;i++){
Vx[i] = sqrt(Temp*IMass[i])*grandom();
Vy[i] = sqrt(Temp*IMass[i])*grandom();
Vz[i] = sqrt(Temp*IMass[i])*grandom();
}

//Quito velocidad del centro de masa
for(i=1;i<=Npt;i++){
Momx = Momx+Vx[i];
Momy = Momy+Vy[i];
Momz = Momz+Vz[i];
}

Momx = Momx*TMassI;
Momy = Momy*TMassI;
Momz = Momz*TMassI;

for(i=1;i<=Npt;i++){
Vx[i] = Vx[i]-Momx;
Vy[i] = Vy[i]-Momy;
Vz[i] = Vz[i]-Momz;
}

puts("\t\t VELOCIDADES INICIALIZADAS...\n");

}


//Generador de números aleatorios
float grandom( void ){

float x,y,gx,gy,g;

x=frand();y=frand();

gx = sqrt( -2.0*log(1.0-(x)) )*cos(2.0*Pi*(y));
gy = sqrt( -2.0*log(1.0-(x)) )*sin(2.0*Pi*(y));

g = (frand() < 0.5)? gx : gy ;
g = (frand() < 0.5)? -g : g ;

return g;
}


//Calcular la fuerza
void Force(void){

for(i=1;i<=Npt;i++){
Fx[i] = 0.0;
Fy[i] = 0.0;
Fz[i] = 0.0;
}
Virx  = 0.0;
Viry  = 0.0;
Virz  = 0.0;

IU    = 0.0;

for(i=1;i<=Npt-1;i++){
Xi = Rx[i];
Yi = Ry[i];
Zi = Rz[i];
for(j=i+1;j<=Npt;j++){
Rijx = Xi-Rx[j];
Rijy = Yi-Ry[j];
Rijz = Zi-Rz[j];

Rijx = Rijx-rint(Rijx*LboxxI)*Lboxx;
Rijy = Rijy-rint(Rijy*LboxyI)*Lboxy;
Rijz = Rijz-rint(Rijz*LboxzI)*Lboxz;

Rij  = sqrt(Rijx*Rijx+Rijy*Rijy+Rijz*Rijz);
DRij = Rij-Rij_0;

U  = exp(-2.0*alpha*DRij)-2.0*exp(-alpha*DRij);
U  = U*U0;
IU = IU+U;

DU = exp(-2.0*alpha*DRij)-exp(-alpha*DRij);
DU = -2.0*alpha*U0*DU;

Fijx = -DU*Rijx/Rij;
Fijy = -DU*Rijy/Rij;
Fijz = -DU*Rijz/Rij;

Fx[i] = Fx[i] + Fijx;
Fy[i] = Fy[i] + Fijy;
Fz[i] = Fz[i] + Fijz;
Fx[j] = Fx[j] - Fijx;
Fy[j] = Fy[j] - Fijy;
Fz[j] = Fz[j] - Fijz;

Virx = Virx + Fijx*Rijx;
Viry = Viry + Fijy*Rijy;
Virz = Virz + Fijz*Rijz;

}}

}


//Calcular valores instantáneos
void Instant(void){

IpressP = (Virx+Viry+Virz)/(3.0*Vol);
IpressK = (KENAX+KENAY+KENAZ)/(3.0*Vol);
IPressT = IpressP+IpressK;

IKEn  = IKEn/(double)Npt;
IU    = IU/(double)Npt;
ITEn  = IKEn+IU;
ITEn2 = ITEn*ITEn;

if(istep==1){     
printf("\tU INICIAL        %f \n",IU);
printf("\tK INICIAL        %f \n",IKEn);   
printf("\tE INICIAL        %f \n",IKEn+IU);
printf("\n istep*dt \t\t IKEn \t\t IU \t\t ITEn \t\t ITemp \t\t IPressT\n");
}

//Imprimir valores a la pantalla
if (istep%rprint==0) printf(" %f \t\t %f \t %f \t %f \t %f \t %f \n",istep*dt,IKEn,IU,ITEn,ITemp,IPressT);

if (istep>intstepinit){
Ttemp    = Ttemp+ITemp;
TIKEn    = TIKEn+IKEn;
TIU      = TIU+IU;
TIPressT = TIPressT+IPressT;
TITEn    = TITEn+ITEn;
TITEn2   = TITEn2+ITEn2;
}

}


//Calcular valores promedio
void Prom(void){

intstat = (double)(intstep-intstepinit);

Ptemp    = Ttemp/intstat;
PIKEn    = TIKEn/intstat;
PIU      = TIU/intstat;
PIPressT = TIPressT/intstat;
PITEn    = TITEn/intstat;
PITEn2   = TITEn2/intstat;
c_v      = (PITEn2-PITEn*PITEn)/Temp/Temp;

printf("\n\t\t< T >    %f \n",Ptemp);
printf("\t\t< K >    %f \n",PIKEn);
printf("\t\t< U >    %f \n",PIU);
printf("\t\t< E >    %f \n",PITEn);
printf("\t\t< E*E >  %f \n",PITEn2);
printf("\t\t< P >    %f \n",PIPressT);
printf("\t\t  D      %f \n",DifC_R);
printf("\t\t  c_v    %f \n",c_v);
printf("\t\t  C_v    %f \n",c_v*(double)Npt);

FILE *promv;
promv = fopen("PromvD0844T071.dat","w");

fprintf(promv,"\ta            %f \n",a);
fprintf(promv,"\tb            %f \n",b);
fprintf(promv,"\tc            %f \n",c);
fprintf(promv,"\tTemp         %f \n",Temp);
fprintf(promv,"\tmass         %f \n",mass);
fprintf(promv,"\tTMass        %f \n",TMass);
fprintf(promv,"\tTMassI       %f \n",TMassI);
fprintf(promv,"\tLboxx        %f \n",Lboxx);
fprintf(promv,"\tLboxy        %f \n",Lboxy);
fprintf(promv,"\tLboxz        %f \n",Lboxz);
fprintf(promv,"\tVol          %f \n",Vol);
fprintf(promv,"\tDens         %f \n",Dens);
fprintf(promv,"\tpack         %f \n",pack);
fprintf(promv,"\tdt           %f \n",dt);
fprintf(promv,"\tdr           %f \n",dr);
fprintf(promv,"\tdk           %f \n",dk);
fprintf(promv,"\tintstep      %d \n",intstep);
fprintf(promv,"\tintstepinit  %d \n",intstepinit);
fprintf(promv,"\tLboxxI       %f \n",LboxxI);
fprintf(promv,"\tLboxyI       %f \n",LboxyI);
fprintf(promv,"\tLboxzI       %f \n",LboxzI);
fprintf(promv,"\tLmin         %f \n",Lmin);
fprintf(promv,"\tNBin        %f \n",NBin);
fprintf(promv,"\tiNBin       %d \n",iNBin);
fprintf(promv,"< T >          %f \n",Ptemp);
fprintf(promv,"< K >          %f \n",PIKEn);
fprintf(promv,"< U >          %f \n",PIU);
fprintf(promv,"< E >          %f \n",PITEn);
fprintf(promv,"< E*E >        %f \n",PITEn2);
fprintf(promv,"< P >          %f \n",PIPressT);
fprintf(promv,"  D            %f \n",DifC_R);
fprintf(promv,"  c_v          %f \n",c_v);
fprintf(promv,"  C_v          %f \n",c_v*(double)Npt);

fclose(promv);
}


//Vectores y variables inicializan con 0
void Zero(void){

for(i=0;i<=Npt;i++){
Rx[i]     = Ry[i]     = Rz[i]     = 0.0;
Vx[i]     = Vy[i]     = Vz[i]     = 0.0;
Fx[i]     = Fy[i]     = Fz[i]     = 0.0;
RxRe[i]   = RyRe[i]   = RzRe[i]   = 0.0;
RxAf[i]   = RyAf[i]   = RzAf[i]   = 0.0;
FxBe[i]   = FyBe[i]   = FzBe[i]   = 0.0;

Mass[i]   = 0.0;
IMass[i]  = 0.0;
}

for(i=0;i<=maxint-1;i++){
Hist[i] = 0.0;
Rdf1[i] = 0.0;
R[i]    = 0.0;
}

Ttemp    = 0.0;
TIKEn    = 0.0;
TIU      = 0.0;
TIPressT = 0.0;
TITEn    = 0.0;
TITEn2   = 0.0;
isample  = 0;
KENAX    = 0.0;
KENAY    = 0.0;
KENAZ    = 0.0;
k1       = 0;
k2       = 0;
}


//Calcular r(t+dt)
void Move1(void){

if(istep>intstepinit) k1=k1+1;

for(i=1;i<=Npt;i++){
RxAf[i] = Rx[i]+dt*Vx[i]+0.5*dt2*Fx[i]*IMass[i];
RyAf[i] = Ry[i]+dt*Vy[i]+0.5*dt2*Fy[i]*IMass[i];
RzAf[i] = Rz[i]+dt*Vz[i]+0.5*dt2*Fz[i]*IMass[i];
}

for (i=1;i<=Npt;i++){
Rx[i] = RxAf[i]-rint(RxAf[i]*LboxxI)*Lboxx;
Ry[i] = RyAf[i]-rint(RyAf[i]*LboxyI)*Lboxy;
Rz[i] = RzAf[i]-rint(RzAf[i]*LboxzI)*Lboxz;

FxBe[i] = Fx[i];
FyBe[i] = Fy[i];
FzBe[i] = Fz[i];
}

//Archivo de posiciones finales
if(istep==intstep){
	FILE *dataf;
	dataf = fopen("datafD0844T071.xyz","w");
	fprintf(dataf,"%d \n",Npt);
	fprintf(dataf,"FRAME %d\n",istep+1);
	for(i=1;i<=Npt;i++){
		fprintf(dataf,"C %f %f %f\n",Rx[i],Ry[i],Rz[i]);
	}
	fclose(dataf);
}

}


//Calcular v(t+dt)
void Move2(void){

IKEn   = 0.0;
KENAXS = 0.0;
KENAYS = 0.0;
KENAZS = 0.0;

if(istep>intstepinit) k2=k2+1;

for(i=1;i<=Npt;i++){
Vx[i] = Vx[i]+0.5*dt*(Fx[i]+FxBe[i])*IMass[i];
Vy[i] = Vy[i]+0.5*dt*(Fy[i]+FyBe[i])*IMass[i];
Vz[i] = Vz[i]+0.5*dt*(Fz[i]+FzBe[i])*IMass[i];
}

//Quito velocidad del centro de masa
Momx = 0.0;
Momy = 0.0;
Momz = 0.0; 

for(i=1;i<=Npt;i++){
Momx = Momx+Vx[i];
Momy = Momy+Vy[i];
Momz = Momz+Vz[i];
}

Momx = Momx*TMassI;
Momy = Momy*TMassI;
Momz = Momz*TMassI;

for(i=1;i<=Npt;i++){
Vx[i] = Vx[i]-Momx;
Vy[i] = Vy[i]-Momy;
Vz[i] = Vz[i]-Momz;
}

for(i=1;i<=Npt;i++){
Vx2 = Vx[i]*Vx[i];
Vy2 = Vy[i]*Vy[i];
Vz2 = Vz[i]*Vz[i];

IKEn = IKEn+0.5*Mass[i]*(Vx2+Vy2+Vz2);

KENAXS = KENAXS+(Mass[i]*Vx2);
KENAYS = KENAYS+(Mass[i]*Vy2);
KENAZS = KENAZS+(Mass[i]*Vz2);
}

KENAX = KENAXS;
KENAY = KENAYS;
KENAZ = KENAZS;

ITemp = 2.0*IKEn/(3.0*(double)Npt-3.0);

Tratio = sqrt(Temp/ITemp);

for(i=1;i<=Npt;i++){
Vx[i] = Vx[i]*Tratio;
Vy[i] = Vy[i]*Tratio;
Vz[i] = Vz[i]*Tratio;
}

}


//Calcular histograma
void Sample(void){

isample = isample+1;

for(i=1;i<=Npt-1;i++){
Xi = Rx[i];
Yi = Ry[i];
Zi = Rz[i];
for(j=i+1;j<=Npt;j++){
Rijx = Xi-Rx[j];
Rijy = Yi-Ry[j];
Rijz = Zi-Rz[j];

Rijx = Rijx-rint(Rijx*LboxxI)*Lboxx;
Rijy = Rijy-rint(Rijy*LboxyI)*Lboxy;
Rijz = Rijz-rint(Rijz*LboxzI)*Lboxz;

Rij = sqrt(Rijx*Rijx+Rijy*Rijy+Rijz*Rijz);

iRij = (int)(Rij/dr)+1;

if(iRij<=NBin) Hist[iRij]=Hist[iRij]+2;
}}

}


//Calcular g(r)
void Rdf(void){

FILE *rad;
rad = fopen("RdfD0844T071.dat","w");

for(i=1;i<=NBin;i++){
R2   = (double)i*dr;
R1   = (double)(i-1)*dr;
Dvol = 4.0*Pi*(R2*R2*R2-R1*R1*R1)/3.0;

Rdf1[i] = (double)Hist[i]/(Dens*(double)Npt*(double)isample*Dvol);

R[i]    = ((double)i-0.5)*dr;
fprintf(rad,"%f %f\n",R[i],Rdf1[i]);
}

fclose(rad);
}
