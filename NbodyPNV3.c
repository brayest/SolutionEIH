#include <stdlib.h> /* Librerias Adicionales */
#include <stdio.h> /* Sistema */
#include <math.h>  /* Librerias de Matematica */


void RK4MATRIX ( double **newy, double **y0, double paso, double time, void(*derivada)(double, double **, double **, double *, int, int, int), int nvar, double *Mases, int CADA, int IT)
  {        
    int i,j,k,w,PRINT;
    int coord = 6;
    double **k0,**k1,**k2,**k3,**z;
    k0= malloc(nvar *sizeof(double));
    k1= malloc(nvar *sizeof(double));
    k2= malloc(nvar *sizeof(double));
    k3= malloc(nvar *sizeof(double));
    z= malloc(nvar *sizeof(double)); 
    
    for ( i = 0; i < nvar; i++) {
     k0[i] = (double *)malloc( coord *sizeof(double));
     k1[i] = (double *)malloc( coord *sizeof(double));
     k2[i] = (double *)malloc( coord *sizeof(double));
     k3[i] = (double *)malloc( coord *sizeof(double));
     z[i] = (double *)malloc( coord *sizeof(double)); 
    }
    
    PRINT = 0;
     if ( IT != -1 && IT % CADA == 0 ) {
      PRINT = 1;
    }
    

    (*derivada)(time,y0,k0,Mases,nvar,PRINT,0);
        
    for (i = 0; i < nvar; i++)  {
     for ( j = 0; j < coord; j++) {
        z[i][j] = y0[i][j] + k0[i][j]*paso/2.0;
     }
    } 

    (*derivada)(time + paso/2.0,z,k1,Mases,nvar,PRINT,0);
    
    for (i = 0; i < nvar; i++) {
     for (j = 0; j < coord; j++) {
      z[i][j] = y0[i][j] + k1[i][j]*paso/2.0;
     }
    } 
      
    (*derivada)(time+paso/2.0,z,k2,Mases,nvar,PRINT,0);
      
    for (i = 0; i < nvar; i++) {
     for ( j = 0; j < coord; j++) {
      z[i][j] = y0[i][j] + k2[i][j]*paso;
     }
    } 
      
    (*derivada)(time+paso,z,k3,Mases,nvar,PRINT,1);
    
    for ( i = 0; i < nvar; i++) {
     for ( j = 0; j < coord; j++) {
      newy[i][j] = y0[i][j] + (paso/6.0)*(k0[i][j] + 2.0*k1[i][j] + 2.0*k2[i][j] + k3[i][j]);
     }
    }
    if ( PRINT == 1 ) {    
      FILE *out1;
      out1 = fopen("SolucionNA.dat","a");
      fprintf(out1, "%.10f         ",time);
      for (i = 0; i < nvar; i++) {
        fprintf(out1, "%.30f     %.30f	%.30f	%.30f ",k3[i][3],k3[i][4],k3[i][5],sqrt( pow(k3[i][3],2) + pow(k3[i][4],2) + pow(k3[i][5],2) ) );
      }
      fprintf(out1,"\n");
      fclose(out1);
    }                                  
    
    free (k0);
    free (k1);
    free (k2);
    free (k3);
    free (z);
  }



void PostNewtonNcuerposV3 ( double x, double **y, double **dvdt, double *M, int N, int P, int SN) {

  double G,c;
  int z,k,w,v,r, coordenadas, relaciones;
  
  coordenadas = 3;
  G = 6.6726e-11;
  c = 299792458;  
  
  double ***Xab,**Rab;
  
  Xab = malloc( N * sizeof(double **));
  Rab = malloc( N * sizeof(double *));
  for ( z = 0; z < N ; z++) {
    Xab[z] = malloc( N * sizeof(double *));
    Rab[z] = malloc(N * sizeof(double));
    for ( v = 0; v < N; v++) {
      Xab[z][v] = malloc( coordenadas * sizeof(double));

    }
  }
    
  
  for ( z = 0; z < N; z++)
    for ( v = 0; v < N; v++) {
      for ( w = 0; w < coordenadas; w++) {
        Xab[z][v][w] = y[z][w] - y[v][w];
      }
      Rab[z][v] = sqrt( pow(Xab[z][v][0],2)+ pow(Xab[z][v][1],2) + pow(Xab[z][v][2],2) );
    }

  int N1,l;
  double *SUM, *CONST;
  N1 = 10;
  SUM = malloc( N1 * sizeof(double));
  CONST = malloc( N1 * sizeof(double));

  FILE *out2;
  if ( P == 1 && SN ==1 ) {
    out2 = fopen("SolucionNC.dat","a");
    fprintf(out2, "%.5f         ",x);    
  }
  
  for ( z = 0; z < N; z++) {
    for ( v = 0; v < 3; v++)
      dvdt[z][v] = y[z][v+3];
    for ( v = 3; v < 6; v++) {
      for( k = 0; k < N; k++) {
        if ( k != z ) {
          SUM[0] = SUM[0] + (G*M[k]*Xab[z][k][v-3])/pow(Rab[z][k],3);
          for ( l = 3; l < 6; l++) {
             CONST[0] = CONST[0] + pow(y[z][l],2); // Va^2
             CONST[1] = CONST[1] + pow(y[k][l],2); // Vb^2
             CONST[2] = CONST[2] + y[z][l]*y[k][l]; // Va*Vb
             CONST[3] = CONST[3] + ((Xab[z][k][l-3]/Rab[z][k])*y[k][l]); //nab*vb
             CONST[4] = CONST[4] + (Xab[z][k][l-3]/Rab[z][k])*(4*y[z][l]-3*y[k][l]); //nab*(4va-3vb)
             }
          if ( N > 2 ) {
            for ( w = 0; w < N; w++ ) {
              if (w != z && w !=k) {
                for ( l=3; l<6; l++) {
                  CONST[5] = CONST[5] + (Xab[z][k][l-3]/Rab[z][k])*(Xab[k][w][l-3]/Rab[k][w]);
                }
                SUM[3] = SUM[3] + (G*M[w])/Rab[z][w]; // G*Mc/Zac
                SUM[4] = SUM[4] + (G*M[w])/Rab[k][w]; // G*Mc/Zbc
                SUM[5] = SUM[5] + ((G*M[w]*Rab[z][k])/pow(Rab[k][w],2))*CONST[5];
                SUM[6] = SUM[6] + (pow(G,2)*M[k]*M[w]*(Xab[k][w][v-3]/Rab[k][w]))/(Rab[z][k]*pow(Rab[k][w],2));
                CONST[5] = 0;
              }
            }
          }
          SUM[1] = SUM[1] + ((G*M[k]*Xab[z][k][v-3])/pow(Rab[z][k],3))*( CONST[0] + 2*CONST[1] - 4*CONST[2] - (3.0/2.0)*pow(CONST[3],2) 
          -((5*G*M[z])/Rab[z][k]) -((4*G*M[k])/Rab[z][k]) - 4*SUM[3] - SUM[4]  + 0.5*SUM[5]);
          SUM[2] = SUM[2] + ((G*M[k])/pow(Rab[z][k],2))*CONST[4]*(y[z][v]-y[k][v]);        
          }
        for ( l=0; l<N1; l++) { CONST[l] = 0; }
        SUM[3] = 0;
        SUM[4] = 0;
        SUM[5] = 0;
      }
      SUM[9] = (1/pow(c,2))*(-SUM[1] + SUM[2] - (7.0/2.0)*SUM[6] );
      if ( P == 1 && SN ==1 ) {
//        fprintf(out2, "%.30f	%.30f	%.30f	%.30f	%.30f  ",SUM[0],SUM[1],SUM[2],(7.0/2.0)*SUM[6],SUM[9]);
        fprintf(out2, "%.30f	",SUM[9]);
      }
      dvdt[z][v] = -SUM[0] + SUM[9];
      for ( l = 0; l < N1; l++) { SUM[l] = 0; }
    }
//    printf("\n");
  }
  if ( P == 1 && SN ==1 ) {
    fprintf(out2,"\n");
    fclose(out2);
  }
  
  free(Xab);
  free(Rab);
  free(SUM);
  free(CONST);

}



void NewtonNcuerpos ( double x, double **y, double **dvdt, double *M, int N, int P, int SN) {
  double G,c;
  int z,k,w,v,r, coordenadas, relaciones;
  
  coordenadas = 3;
  G = 6.6726e-11;
  c = 299792458;  
  
  double ***Xab,**Rab;
  
  Xab = malloc( N * sizeof(double **));
  Rab = malloc( N * sizeof(double *));
  for ( z = 0; z < N ; z++) {
    Xab[z] = malloc( N * sizeof(double *));
    Rab[z] = malloc(N * sizeof(double));
    for ( v = 0; v < N; v++) {
      Xab[z][v] = malloc( coordenadas * sizeof(double));

    }
  }
    
  
  for ( z = 0; z < N; z++)
    for ( v = 0; v < N; v++) {
      for ( w = 0; w < coordenadas; w++) {
        Xab[z][v][w] = y[z][w] - y[v][w];
      }
      Rab[z][v] = sqrt( pow(Xab[z][v][0],2)+ pow(Xab[z][v][1],2) + pow(Xab[z][v][2],2) );
    }

  int N1,l;
  double *SUM;
  N1 = 10;
  SUM = malloc( N1 * sizeof(double));

    
  
  for ( z = 0; z < N; z++) {
    for ( v = 0; v < 3; v++)
      dvdt[z][v] = y[z][v+3];
    for ( v = 3; v < 6; v++) {
      for( k = 0; k < N; k++) {
        if ( k != z ) {
          SUM[0] = SUM[0] + (G*M[k]*Xab[z][k][v-3])/pow(Rab[z][k],3);
        }
      }
      dvdt[z][v] = -SUM[0];
      for ( l = 0; l < N1; l++) {
        SUM[l] = 0;  
      }
    }
  }
  
  free(Xab);
  free(Rab);
  free(SUM);

}




double EnergiaN (int N, double **Move, double *mases) {

    int v,z,j,k,w,coordenadas;
    coordenadas = 3;
    double ***Xab,**Rab,*AP,*SUMAS,Energia,G,c;
    
    
    G = 6.6726e-11;
    c = 299792458;  
    
    SUMAS = malloc( 10 * sizeof(double));
    AP = malloc( 10 * sizeof(double));
  
  
    Xab = malloc( N * sizeof(double **));
    Rab = malloc( N * sizeof(double *));
    for ( z = 0; z < N ; z++) {
      Xab[z] = malloc( N * sizeof(double *));
      Rab[z] = malloc(N * sizeof(double *));
      for ( v = 0; v < N; v++) {
        Xab[z][v] = malloc( coordenadas * sizeof(double));
      }
    }
    
    
    for ( j = 0; j < N; j++) {
      for ( k = 0; k < N; k++) {
        for ( w = 0; w < coordenadas; w++) {
          Xab[j][k][w] = Move[j][w] - Move[k][w];
        }
        Rab[j][k] = sqrt( pow(Xab[j][k][0],2)+ pow(Xab[j][k][1],2) + pow(Xab[j][k][2],2) );
      }
    }
    
    
    for ( j = 0; j < N; j++) {
      for ( k = 0; k < coordenadas; k++) {
        AP[0] = AP[0] + pow(Move[j][k+3],2); //Va^2
      }
      SUMAS[0] = SUMAS[0] + mases[j]*AP[0];      //maVa^2
      for ( k = 0; k < N; k++) {
        if ( k!=j ) {
          SUMAS[1] = SUMAS[1] + (0.5*G*mases[j]*mases[k])/Rab[j][k]; //G*ma*mb/Rab
        }
      }
      AP[0] = 0;
    }
    
    Energia = 0.5*SUMAS[0] - SUMAS[1];
    
    
    for ( w = 0; w < 10; w++) {
      SUMAS[w] = 0;
      AP[w] = 0;
    }
    

    free(Xab);
    free(Rab);
    free(SUMAS);
    free(AP);
    
    return Energia;

}

double EnergiaPN2 (int N, double **Move, double *mases) {

    int v,z,j,k,w,coordenadas,l;
    coordenadas = 3;
    double ***Xab,**Rab,*AP,*SUMAS,Energia,G,c;
    
    
    G = 6.6726e-11;
    c = 299792458;  
    
    SUMAS = malloc( 10 * sizeof(double));
    AP = malloc( 10 * sizeof(double));
  
  
    Xab = malloc( N * sizeof(double **));
    Rab = malloc( N * sizeof(double *));
    for ( z = 0; z < N ; z++) {
      Xab[z] = malloc( N * sizeof(double *));
      Rab[z] = malloc(N * sizeof(double *));
      for ( v = 0; v < N; v++) {
        Xab[z][v] = malloc( coordenadas * sizeof(double));
      }
    }
    
    
    for ( j = 0; j < N; j++) {
      for ( k = 0; k < N; k++) {
        for ( w = 0; w < coordenadas; w++) {
          Xab[j][k][w] = Move[j][w] - Move[k][w];
        }
        Rab[j][k] = sqrt( pow(Xab[j][k][0],2)+ pow(Xab[j][k][1],2) + pow(Xab[j][k][2],2) );
      }
    }
    
    
    for ( j = 0; j < N; j++) {
      for ( k = 0; k < coordenadas; k++) {
        AP[0] = AP[0] + pow(Move[j][k+3],2); //Va^2
      }
      SUMAS[0] = SUMAS[0] + mases[j]*AP[0];      //maVa^2
      SUMAS[2] = SUMAS[2] + mases[j]*pow(AP[0],2); //maVa^4
      for ( k = 0; k < N; k++) {
        if ( k!=j ) {
          SUMAS[3] = SUMAS[3] + G*mases[k]/Rab[j][k]; //G*mb/Rab
          SUMAS[1] = SUMAS[1] + (0.5*G*mases[j]*mases[k])/Rab[j][k]; //G*ma*mb/Rab
          for ( w = 0; w < coordenadas; w++ ) {
            AP[1] = AP[1] + Move[j][w+3]*(Xab[j][k][w]/Rab[j][k]); //va*nab
            AP[2] = AP[2] + Move[k][w+3]*(Xab[j][k][w]/Rab[j][k]); //vb*nab
            AP[3] = AP[3] + Move[j][w+3]*Move[k][w+3]; //va*vb
          }
          SUMAS[5] = SUMAS[5] + ((G*mases[k])/Rab[j][k])*(7.0*AP[3]+AP[1]*AP[2]);
          for ( l = 0; l < N; l++ ) {
            if ( l!=j && l!=k  ) {
              SUMAS[7] = SUMAS[7] + (pow(G,2)*mases[k]*mases[l])/(Rab[j][k]*Rab[j][l]); //G^2*mb*mc/rab*rac
            }
          }
          SUMAS[8] = SUMAS[8] + SUMAS[7];
        }
      }
      SUMAS[4] = SUMAS[4] + mases[j]*AP[0]*SUMAS[3];
      SUMAS[6] = SUMAS[6] + mases[j]*SUMAS[5];
      SUMAS[9] = SUMAS[9] + mases[j]*SUMAS[8];
      AP[0] = 0;
      AP[1] = 0;
      AP[2] = 0;
      AP[3] = 0;
    }
    
    Energia = 0.5*SUMAS[0] - SUMAS[1] + (1/pow(c,2))*( (3.0/8.0)*SUMAS[2] +(3.0/2.0)*SUMAS[4] - 0.25*SUMAS[6] + 0.5*SUMAS[9] );
    
    
    for ( w = 0; w < 10; w++) {
      SUMAS[w] = 0;
      AP[w] = 0;
    } 
    
    free(Xab);
    free(Rab);
    free(SUMAS);
    free(AP);
    
    return Energia;

}



void CondIn (int BDS, double *M, double **DAT, double **IN) {
  
  int i, es = 0;  
  double RL = 0;
  double G = 6.6726e-11;
  
  for (i=0; i < BDS; i++) {
    if ( DAT[i][2] != 0 ) {
      es = (int) DAT[i][2];
      if ((int)IN[es][4] == 0) {
        IN[es][0] = (DAT[es][0]*(1-pow(DAT[es][1],2)))/(1+DAT[es][1]); 
        IN[es][1] = 0; 	// Y1
        IN[es][2] = 0;	// Z1
        IN[es][3] = 0;	// Vx
        IN[es][4] = sqrt(G*M[(int)DAT[es][2]]*((2/IN[es][0]) - 1/DAT[es][0]));	// Vy
        IN[es][5] = 0;	// Vz
      }
      IN[i][0] = IN[es][0] + (DAT[i][0]*(1-pow(DAT[i][1],2)))/(1+DAT[i][1]);  // X1 Distancia Sol-Tierra + Distancia Luna-Tierra
      IN[i][1] = 0; 	// Y1
      IN[i][2] = 0;	// Z1
      IN[i][3] = 0;	// Vx
      RL = (DAT[i][0]*(1-pow(DAT[i][1],2)))/(1+DAT[i][1]);
      IN[i][4] = IN[es][4] + sqrt(G*M[es]*((2/RL) - 1/DAT[i][0]));	// Vy Velocidad Tierra+Luna
      IN[i][5] = 0;	// Vz
      
    }
  }
  
  for ( i=1; i < BDS; i++ ) {  
    if ( IN[i][4] == 0  ) {
      IN[i][0] = (DAT[i][0]*(1-pow(DAT[i][1],2)))/(1+DAT[i][1]);  // X1
      IN[i][1] = 0; 	// Y1
      IN[i][2] = 0;	// Z1
      IN[i][3] = 0;	// Vx
      IN[i][4] = sqrt(G*M[(int)DAT[i][2]]*((2/IN[i][0]) - 1/DAT[i][0]));	// Vy
      IN[i][5] = 0;	// Vz
    }	
  }
  
  double momentox, momentoy, momentoz = 0;
  for (i = 1; i < BDS; i++ ) {
    momentoy = momentoy + M[i]*IN[i][4];
  }
  
  // El cuerpo de mayor influencia - Sol
  IN[0][0] = 0;  // X1
  IN[0][1] = 0; 	// Y1
  IN[0][2] = 0;	// Z1
  IN[0][3] = 0;	// Vx
  IN[0][4] = -momentoy/M[0];	// Vy
  IN[0][5] = 0;	// Vz
  
}


void PosIN (double **IN, double **PR, double *M, int M1, int M2, int N, int CADA, int IT ) {

  double *MA,**IN1,**MV1,OR,TM,E,ex,an,T,h,TMID;
  int var,z,i,j,k;

  MA = malloc( 2 *sizeof(double)); 
  IN1 = malloc( 2 *sizeof(double));
  MV1 = malloc( 2 *sizeof(double));
  
  var = 6;

  for ( i = 0; i < 2; i++) {
    IN1[i] = (double *)malloc( var *sizeof(double));
    MV1[i] = (double *)malloc( var *sizeof(double));
  }

  ex = PR[M2][1];
  an = (M_PI*PR[M2][4])/180.00;
  MA[0] = M[M1];
  MA[1] = M[M2];
  for ( j = 0; j < 6; j++ ) {
    IN1[0][j] = 0;
    IN1[1][j] = IN[M2][j] - IN[M1][j];
    MV1[0][j] = 0;
    MV1[1][j] = IN[M2][j] - IN[M1][j];
   }
   
  
  OR = PR[M2][3];
  TMID = 0;
  if (an > M_PI) {
    TMID = OR/2.00;
    an = an - M_PI;
  }
  
  E = 2.0*atan( sqrt((1-ex)/(1+ex))*tan(an/2.0)  ); // Anomalia Excentrica
  TM = fabs((OR/(2*M_PI))*(E-ex*sin(E))) + TMID;
  
  printf("Angulo:	%.5f,	Tiempo:		%.5f	\n",E,TM);
  h = ((TM)/10000.00);
  T = 0;
  while ( T < TM ) {
    RK4MATRIX(MV1,IN1,h,T,NewtonNcuerpos,2,MA,CADA,IT);
    
 //   for ( j = 0; j < 2; j++)
      for ( k = 0; k < var; k++)
        IN1[1][k] = MV1[1][k]; 
      i++;
      T += h;
    }
  T = 0;
  h =0;
  
  for ( j = 0; j < var; j++ ) {
//    IN[M1][j] = MV1[0][j];
    IN[M2][j] = MV1[1][j];
  }
  
  free(MA);
  free(IN1);
  free(MV1);

}

main (int argc, char *argv[]) {
  // Definición numero de cuerpos
  int cuerpos = 3;
  int variables = 6;
  int coordenadas = 3;
  int i,j,k,z,w;
  
  // constantes
  double G,c;  
  G = 6.6726e-11;
  c = 299792458;
  
  // Definicion de variables  respecto a cantidad de cuerpos
  double *masas, **iniciales, **movimiento, MS,AU, *SM,EX,PC,**param;
  
  // Vectores 
  masas = malloc(cuerpos  *sizeof(double)); 
  iniciales = malloc(cuerpos *sizeof(double));
  movimiento = malloc(cuerpos  *sizeof(double));
  param = malloc(cuerpos *sizeof(double)); // parametros de cada cuerpos
   
  
  // Matrices : Con la siguiente forma
  //
  // Masa	x	y	z	vx	vy	vz
  //  1		*	*	*	*	*	*
  //  2		*	*	*	*	*	*
  //  3		*	*	*	*	*	*
  
  for ( i = 0; i < cuerpos; i++) {
    iniciales[i] = (double *)malloc( variables *sizeof(double));
    movimiento[i] = (double *)malloc( variables *sizeof(double));
  }
  
  // Parametros
  //		0		1		2			3		4			5			6
  // Masa	SemiejeMayor	Excentricidad	CuerpoMayorEfectoG	Orbita		Angulo Inicial		Angulo Inclinacion	Argumento Peri
  //	1	*		*		*			*		*			*			*
  //	2	*		*		*			*		*			*			*
  //	3	*		*		*			*		*			*			*

  for ( i = 0; i < cuerpos; i++) {
    param[i] = (double *)malloc( 7 *sizeof(double));
  }
  
  
  // Variables de Control
  double dt, tiempo, Nit,orbita,E0,time;
  int print;
  char *Type = argv[4];
  
  time = strtod(argv[1], NULL); // En periodos orbitales 
  print = atoi(argv[2]);
  Nit = atoi(argv[3]);
  

// Masa Solar 
MS = 1.9891e30;
// Unidad Astronómica
AU = 149597870700.00;
// Parsec
PC = 206265*AU;


//Definicion de parametros 

//Masas  
masas[0] = 1e6*MS; // el cuerpo de mayor influencia es siempre el "0"
masas[1] = 1*MS;
masas[2] = 1.898e27;
//masas[3] = 50000*MS;
//masas[4] = 1*MS;
//masas[5] = 5.684e26;
//masas[6] = 1.345e23;
//masas[7] = 1.25e22;
// Semi eje Mayor respecto al cuerpo de mayor influencia
param[0][0] = 0;	
param[1][0] = 54.4*AU;
param[2][0] = 0.005*AU; 
//param[3][0] = 20*AU;
//param[4][0] = 0.007155 *AU;
//param[5][0] = 9.5549*AU;
//param[6][0] = 0.008168*AU;
//param[7][0] = 39.264*AU; 
// Excentricidad
param[0][1] = 0;
param[1][1] = 0.7;
param[2][1] = 0.001;
//param[3][1] = 0.7;
//param[4][1] = 0.0013;
//param[5][1] = 0.05555;
//param[6][1] = 0.008168;
//param[6][1] = 0,244;
// Mayor influencia de 
param[0][2] = 0;  // Para el cuerpos mas influyente es el mismo y no se utiliza
param[1][2] = 0;  	
param[2][2] = 1;
//param[3][2] = 0;
//param[4][2] = 3;
//param[5][2] = 0;
//param[6][2] = 5;
//param[7][2] = 0;
// Angulo inicial del plano en orbita
param[0][4] = 0;
param[1][4] = 0;
param[2][4] = 0;
//param[3][4] = 0;
//param[4][4] = 0;
//param[5][4] = 111;
//param[6][4] = 0;
//param[7][4] = 0;
// Angulo de Inclinacion
param[0][5] = 0;
param[1][5] = 91;
param[2][5] = 0;
//param[3][5] = 0;
//param[4][5] = 0.20; 
//param[5][5] = 2.485240;
//param[6][5] = 0.34854;
//param[7][5] = 17.2;
//Argumento de Periastro
param[0][6] = 0; 
param[1][6] = 90;
param[2][6] = 90;
//param[3][6] = 90;
//param[4][6] = 0;
//param[5][6] = 339.392;
//param[6][6] = 0;
//param[7][6] = 113.76349;

// N cuerpos: Las relaciones entre masas y distancias deben definirse jerarquicamente. 
// Por ejemplo: Para el sistema Sol-Tierra-Luna la distancia del semiejemayor de la luna es relativo
// al orgien del sistema cartesiano de coordenadas no de la tierra.  Por lo tanto es neceario
// definir los valores inciales jeraquicamente. 
   
  CondIn(cuerpos, masas, param, iniciales);  // Condiciones Iniciales
  
  double BIG1;
  for ( i = 0 ; i< cuerpos; i++ ) {
    param[i][3] = sqrt(4*pow(M_PI,2)*pow(param[i][0],3)/((masas[i]+masas[(int)param[i][2]])*G));    // Periodo Orbital
  }
  
    
  double A,A1,R,R1,X,Z,X1,Z1,Y,Y1,V,V1,Vx,Vy,Vx1,Vy1;
  int base;
  
  // Angulo de Inclinacion de la orbtia
  
  
  for (i=1; i<cuerpos; i++) {
    base = (int)param[i][2];
    if ( base != 0 && iniciales[base][2] == 0 ) {
      R1 = iniciales[i][0] - iniciales[base][0];
      A1 = (M_PI*param[i][5])/180.00;
      X1 = R1*cos(A1);
      Z1 = R1*sin(A1);
      R = iniciales[base][0];
      A = (M_PI*param[base][5])/180.00;
      X = R*cos(A);
      Z = R*sin(A);
      iniciales[base][0] = X;
      iniciales[base][2] = Z;
      iniciales[i][0] = X + X1;
      iniciales[i][2] = Z + Z1;
      A1 = (M_PI*param[i][6])/180.00 - M_PI/2.00;
      R1 = iniciales[i][0] - iniciales[base][0];
      X1 = R1*cos(A1);
      Y1 = R1*sin(A1);
      V1 = iniciales[i][4] - iniciales[base][4];
      Vx1 = -V1*sin(A1);
      Vy1 = V1*cos(A1);
      A = (M_PI*param[base][6])/180.00 - M_PI/2.00;
      R = iniciales[base][0];
      X = R*cos(A);
      Y = R*sin(A);
      V = iniciales[base][4];
      Vx = -V*sin(A);
      Vy = V*cos(A);
      iniciales[base][0] = X;
      iniciales[base][1] = Y;
      iniciales[i][0] = X + X1;
      iniciales[i][1] = Y + Y1;
      iniciales[base][3] = Vx;
      iniciales[base][4] = Vy;
      iniciales[i][3] = Vx + Vx1;
      iniciales[i][4] = Vy + Vy1;  
    }
  }
  
  
  for (i=1; i<cuerpos; i++) {
    base = (int)param[i][2];
    if ( base == 0 && iniciales[i][2] == 0 ) {
      R = iniciales[i][0];
      A = (M_PI*param[i][5])/180.00;
      X = R*cos(A);
      Z = R*sin(A);
      iniciales[i][0] = X;
      iniciales[i][2] = Z;
      A = (M_PI*param[i][6])/180.00 - M_PI/2.0;
      R = iniciales[i][0]; 
      X = R*cos(A);
      Y = R*sin(A);
      iniciales[i][0] = X;
      iniciales[i][1] = Y;
      iniciales[i][3] = -iniciales[i][4]*sin(A);
      iniciales[i][4] = iniciales[i][4]*cos(A);  
    }
  }
  

   // Determinar condiciones iniciales basadas en un angulo de referencia incial de orbita
  for (i = 0; i < cuerpos; i++) {
    for ( j = 0; j < variables; j++) {
      movimiento[i][j] = iniciales[i][j];
    }
  }

  for (i = 1; i < cuerpos; i++ ) {
    if (param[i][4] != 0 ) {
     if ( (int)param[i][2] == 0 ) {
       PosIN(iniciales,param,masas,(int)param[i][2],i,Nit,print,-1); // Angulo inicial 
     } else {
         PosIN(movimiento,param,masas,(int)param[i][2],i,Nit,print,-1);
         for ( j = 0; j < variables; j++) {
           iniciales[i][j] = movimiento[i][j] + iniciales[(int)param[i][2]][j];
         }
      }
     }
  }
  
  
    // Conservacion de Momento
  
  double momentox, momentoy, momentoz = 0;
  for (i = 1; i < cuerpos; i++ ) {
    momentox = momentox + masas[i]*iniciales[i][3];
    momentoy = momentoy + masas[i]*iniciales[i][4];
    momentoz = momentoz + masas[i]*iniciales[i][5];
  }
  
  // El cuerpo de mayor influencia - Sol
  iniciales[0][0] = 0;  // X1
  iniciales[0][1] = 0; 	// Y1
  iniciales[0][2] = 0;	// Z1
  iniciales[0][3] = -momentox/masas[0];	// Vx
  iniciales[0][4] = -momentoy/masas[0];	// Vy
  iniciales[0][5] = -momentoz/masas[0];	// Vz

  // Determinar orbita mayor.
  BIG1 = 0;
  for ( i = 0 ; i< cuerpos; i++ ) {
    if ( param[i][3] > BIG1 ) {
      BIG1 = param[i][3];
    }    
  }
  
  
  // Override
  //masas[0] = 1*MS;
  //masas[1] = 1*MS;
  //masas[2] = 1*MS;
//  masas[3] = 100*MS;
  
  
  /* 
  iniciales[0][0] = 10*AU;
  iniciales[0][1] = 0; 
  iniciales[0][2] = 0; 
  iniciales[0][3] = 0;
  iniciales[0][4] = 10000;
  iniciales[0][5] = 0;
  
  iniciales[1][0] = 1*AU;
  iniciales[1][1] = 0; 
  iniciales[1][2] = 0; 
  iniciales[1][3] = 0;
  iniciales[1][4] = 10000;
  iniciales[1][5] = 0;
  
  iniciales[2][0] = 4*AU;
  iniciales[2][1] = 0; 
  iniciales[2][2] = 0; 
  iniciales[2][3] = 0;
  iniciales[2][4] = -10000;
  iniciales[2][5] = 0;
  
  iniciales[3][0] = -100*AU;
  iniciales[3][1] = 0; 
  iniciales[3][2] = 100*AU; 
  iniciales[3][3] = 0;
  iniciales[3][4] = 0;
  iniciales[3][5] = -3000000;
  */
  

   
//  BIG1 = sqrt((4*pow(M_PI,2)*pow(R2,3))/(G*125*MS));
  
  // el parametro de paso se define acorde al period orbital mas largo
  dt = ((time*BIG1)/Nit);
    
  E0 = EnergiaN(cuerpos,iniciales,masas);
  FILE *out1;
  out1 = fopen("SolucionNA.dat","w");
  fclose(out1);
  FILE *out2;
  out2 = fopen("SolucionNC.dat","w");
  fclose(out2);
  FILE *out3;
  out3 = fopen("SolucionNL.dat","w");
  fclose(out3);
      
  



  
  for (i = 0; i < cuerpos; i++) {
    for ( j = 0; j < variables; j++) {
      movimiento[i][j] = iniciales[i][j];
    }
  }
    
  
  char CR1[6];
  char CR2[6];
  
  CR1[0] = 'X';
  CR1[1] = 'Y';
  CR1[2] = 'Z';
  CR1[3] = 'V';
  CR1[4] = 'V';
  CR1[5] = 'V';
  
  CR2[0] = 'x';
  CR2[1] = 'y';
  CR2[2] = 'z';
  CR2[3] = 'x';
  CR2[4] = 'y';
  CR2[5] = 'z';
  
  FILE *out;
  
  if ( *Type == 'N' ) 
    out = fopen("SolucionNN.dat","w");
  if ( *Type == 'P' )
    out = fopen("SolucionNPN.dat","w");
  
      
  for ( i = 0; i < cuerpos; i++) {
    printf("[%d]	Masa:	%.6f	Masas Solares\n",i,masas[i]/MS);
    printf("	Satélite de: [%.0f]     \n",param[i][2]);
    printf("	Semi eje Mayor(m): %.5f	\n",param[i][0]);
    printf("	Semi eje Mayor(AU): %.5f	\n",param[i][0]/AU);
    printf("	Excentricidad	: %.5f	\n",param[i][1]);
    printf("	Período Orbital:	%.5f\n",param[i][3]);
    printf("	Inclinación Orbital:    %.5f\n",param[i][5]);
    printf("	Argumento de Periastro: %.5f\n",param[i][6]);
    printf("	Anomalía Verdadera:     %.5f\n",param[i][4]);
    for ( j = 0; j < 3; j++) {
      printf("	%c%c:	%.5f	AU\n",CR1[j],CR2[j],iniciales[i][j]/AU);
    }
    for ( j = 3; j < 6; j++) {
      printf("	%c%c:	%.5f	m/s\n",CR1[j],CR2[j],iniciales[i][j]);
    }
    printf("\n");
  }
  printf("\n"); 
  printf("Parametro de Paso(s): %.50f\n",dt);
  printf("Energía N:	%.5f\n",EnergiaN(cuerpos,iniciales,masas));
  printf("Energía PN:	%.5f\n",EnergiaPN2(cuerpos,iniciales,masas));
  printf("Numero de Iteraciones:	%f\n",Nit);
  printf("Mayor periodo Orbital: %.5f\n",BIG1);
  

  double Rx,Ry,Rz,Vz,L,Lx,Ly,Lz,LL,EE,E1,E2,R1x,R1y,R1z,R1T;
  i=0;  
  E1=-337121935767753980208388427515078639616.0000000000;
  if ( *Type == 'N' ) {
    out3 = fopen("SolucionNLN.dat","w");
    printf("Tipo: Newton: %s	\n",Type);
    while ( i < Nit ) {
      
      if ( i % print == 0 ) {
        fprintf(out,"%.10f	%.10f		",tiempo/BIG1,EnergiaN(cuerpos,iniciales,masas));
        for ( j = 0; j < cuerpos; j++) {
          if ( j !=0  && j != 1 ) { 
            base = (int)param[j][2];
            Rx = iniciales[j][0]-iniciales[base][0]; 
            Ry = iniciales[j][1]-iniciales[base][1];
            Rz = iniciales[j][2]-iniciales[base][2];
            Vx = iniciales[j][3]-iniciales[base][3];
            Vy = iniciales[j][4]-iniciales[base][4];
            Vz = iniciales[j][5]-iniciales[base][5];
            Lx = masas[j]*(Ry*Vz-Vy*Rz);
            Ly = masas[j]*(Vx*Rz-Rx*Vz);
            Lz = masas[j]*(Rx*Vy-Vx*Ry);
            LL = sqrt(pow(Lx,2)+pow(Ly,2)+pow(Lz,2)); 
            R = sqrt(pow(Rx,2)+pow(Ry,2)+pow(Rz,2));
            V = sqrt(pow(Vx,2)+pow(Vy,2)+pow(Vz,2));           
            R1 = (pow(V,2)/(G*masas[base])) - (1/R);
            Z = ((Rx*Vx+Ry*Vy+Rz*Vz)/(G*masas[base]));
            EE = 0.5*masas[j]*(pow(Vx,2)+pow(Vz,2)+pow(Vz,2)) - (G*masas[base]*masas[j])/R;
            //EX = sqrt(1+(2*EE*pow(LL,2))/(pow(G*masas[j]*masas[base],2)*masas[j]));
            EX = sqrt(pow(R1*Rx-Z*Vx,2)+pow(R1*Ry-Z*Vy,2)+pow(R1*Rz-Z*Vz,2));
            R1x = iniciales[base][0]-iniciales[0][0];
            R1y = iniciales[base][1]-iniciales[0][1];
            R1z = iniciales[base][2]-iniciales[0][2];
            R1T = sqrt(pow(R1x,2)+pow(R1y,2)+pow(R1z,2));
//            A = acos((fabs(Lx*R1x)+fabs(Ly*R1y)+fabs(Lz*R1z))/(R1T*LL));
            A1 = sqrt(0.5);
            A = fabs(acos((Lx*0.5+Ly*0.5+Lz*sqrt(0.5))/LL) - M_PI/4.0);
            L = sqrt(1-pow(EX,2))*cos(A);
            if ( (EE-E1)/E1 > -0.000050001 ) {
              fprintf(out3,"%.10f         ",tiempo/BIG1);
              fprintf(out3,"%.10f		%.10f	%.10f	",EX,(A*180)/M_PI,L);
              fprintf(out3,"\n");
            }
          }
          for ( k = 0; k < variables; k++) {
            fprintf(out, "%.10f	",movimiento[j][k]);
          }
        }
        fprintf(out,"\n");
      }
      RK4MATRIX(movimiento,iniciales,dt,tiempo,NewtonNcuerpos,cuerpos,masas,print,i);
    
      for ( j = 0; j < cuerpos; j++)
        for ( k = 0; k < variables; k++)
          iniciales[j][k] = movimiento[j][k]; 
              
      i++;
      tiempo += dt;
    }
    fclose(out);
    fclose(out3);
  }
  
  if ( *Type == 'P' ) {
    out3 = fopen("SolucionNLP.dat","w");
    printf("Tipo: PostNewton: %s	\n",Type);
    while ( i < Nit ) {
      if ( i % print == 0 ) {
        fprintf(out,"%.10f	%.10f		",tiempo/BIG1,EnergiaPN2(cuerpos,iniciales,masas));
        for ( j = 0; j < cuerpos; j++) {
          if ( j !=0  && j != 1 ) { 
            base = (int)param[j][2];
            Rx = iniciales[j][0]-iniciales[base][0]; 
            Ry = iniciales[j][1]-iniciales[base][1];
            Rz = iniciales[j][2]-iniciales[base][2];
            Vx = iniciales[j][3]-iniciales[base][3];
            Vy = iniciales[j][4]-iniciales[base][4];
            Vz = iniciales[j][5]-iniciales[base][5];
            Lx = masas[j]*(Ry*Vz-Vy*Rz);
            Ly = masas[j]*(Vx*Rz-Rx*Vz);
            Lz = masas[j]*(Rx*Vy-Vx*Ry);
            LL = sqrt(pow(Lx,2)+pow(Ly,2)+pow(Lz,2)); 
            R = sqrt(pow(Rx,2)+pow(Ry,2)+pow(Rz,2));
            V = sqrt(pow(Vx,2)+pow(Vy,2)+pow(Vz,2));           
            R1 = (pow(V,2)/(G*masas[base])) - (1/R);
            Z = ((Rx*Vx+Ry*Vy+Rz*Vz)/(G*masas[base]));
            EE = 0.5*masas[j]*(pow(Vx,2)+pow(Vz,2)+pow(Vz,2)) - (G*masas[base]*masas[j])/R;
            //EX = sqrt(1+(2*EE*pow(LL,2))/(pow(G*masas[j]*masas[base],2)*masas[j]));
            EX = sqrt(pow(R1*Rx-Z*Vx,2)+pow(R1*Ry-Z*Vy,2)+pow(R1*Rz-Z*Vz,2));
            R1x = iniciales[base][0]-iniciales[0][0];
            R1y = iniciales[base][1]-iniciales[0][1];
            R1z = iniciales[base][2]-iniciales[0][2];
            R1T = sqrt(pow(R1x,2)+pow(R1y,2)+pow(R1z,2));
//            A = acos((fabs(Lx*R1x)+fabs(Ly*R1y)+fabs(Lz*R1z))/(R1T*LL));
            //A1 = sqrt(0.5);
            A = fabs(acos((Lx*0.5+Ly*0.5+Lz*sqrt(0.5))/LL) - M_PI/4.0);
            L = sqrt(1-pow(EX,2))*cos(A);
            if ( (EE-E1)/E1 > -0.0010882 ) {
              fprintf(out3,"%.10f         ",tiempo/BIG1);
              fprintf(out3,"%.10f		%.10f	%.10f	",EX,(A*180)/M_PI,L);
              fprintf(out3,"\n");
            }
          }
          for ( k = 0; k < variables; k++) {
            fprintf(out, "%.10f	",movimiento[j][k]);
          }
        }
        fprintf(out,"\n");
      }
      RK4MATRIX(movimiento,iniciales,dt,tiempo,PostNewtonNcuerposV3,cuerpos,masas,print,i);
    
      for ( j = 0; j < cuerpos; j++)
        for ( k = 0; k < variables; k++)
          iniciales[j][k] = movimiento[j][k]; 
    
      i++;
      tiempo += dt;
    }
    fclose(out);
  }
  
  if ( *Type == 'C' ) {
    printf("Tipo: Convergencia: %s	\n", Type);
    char str[1];
    char command[50];
    int n,Nit2,Cada;
    double **iniciales2;
    iniciales2 = malloc(cuerpos *sizeof(double));
    for ( i = 0; i < cuerpos; i++) {
      iniciales2[i] = (double *)malloc( variables *sizeof(double));
    }
    for ( j = 0; j < cuerpos; j++)
      for ( k = 0; k < variables; k++)
        iniciales2[j][k] = iniciales[j][k]; 
    i=0;    
    for ( n = 0; n < 3; n++) {
      dt = (time*BIG1)*(pow(0.5,n)/Nit);
      sprintf(str, "SolucionCPN-%d.dat", n);
      out = fopen(str,"w");
      if ( n != 0 ) {
        Nit2 = Nit*(2*n);
        Cada = print*(2*n);
      } else { Nit2 = Nit; Cada = print;}
      printf("Convergencia: %d		Iteraciones:	%d	Parametro de Paso: %.5f 	\n",n,Nit2,dt);
      while ( i < Nit2 ) {
        if ( i % Cada == 0 ) {
          fprintf(out,"%.5f	",tiempo);
          for ( j = 0; j < cuerpos; j++) {
            for ( k = 0; k < 1; k++) {
              fprintf(out, "%.20f	",movimiento[j][k]);
            }
          }
          fprintf(out,"\n");
        }
//        RK4MATRIX(movimiento,iniciales,dt,tiempo,PostNewtonNcuerposV3,cuerpos,masas);
        RK4MATRIX(movimiento,iniciales,dt,tiempo,PostNewtonNcuerposV3,cuerpos,masas,print,i);
        
        for ( j = 0; j < cuerpos; j++)
          for ( k = 0; k < variables; k++)
            iniciales[j][k] = movimiento[j][k]; 
    
       i++;
       tiempo += dt;
      }
    i = 0;
    tiempo = 0;
    for ( j = 0; j < cuerpos; j++)
      for ( k = 0; k < variables; k++) {
        movimiento[j][k] = iniciales2[j][k]; 
        iniciales[j][k] = iniciales2[j][k]; 
      }
        
    fclose(out);
    }
    system("cat SolucionCPN-0.dat | awk '{print $1}' > Time.txt");
    system("cat SolucionCPN-0.dat | awk '{print $3}' > Sol1.txt");
    system("cat SolucionCPN-1.dat | awk '{print $3}' > Sol2.txt");
    system("cat SolucionCPN-2.dat | awk '{print $3}' > Sol3.txt");
    system("paste Time.txt Sol1.txt Sol2.txt Sol3.txt > Convergencia.dat");
    
  }
  
  free(masas);
  free(iniciales);
  free(movimiento);
  
}

  

