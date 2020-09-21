//****************************************************
//	ktlattice:	implements the polarization trick to calculate the superfluid density
//			of a system of vortices mapped in a coulomb 2D lattice gas.
//
//	R. Diaz-Mendez 
//
//	Stockholm, November 2016.
//
//****************************************************

//****************************************************
//	ktnn:	clon to study the non-neutral version of the problem, where the polarization trick is not longer good.
//		Instead the net vorticity is recorded.
//
//	Stockholm, January 2017.
//**************************************************** 

//****************************************************
//	kt:	Clon to include all models into a single (versatile) code. 
//
//	Stockholm, January 2017.
//**************************************************** 

//*********************
//  	A reweighting technique is used when Tmin is negative and nTs=1
//  	A number of rwt temperatures is extrapolated uniformly in
//	a vicinity of fabs(Tmin) around Tmax. 
//	 


//********   Last Update: 2017/07/11



#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>

#include "aleat.c"



#define pi acos(-1)
#define Z 4  
#define N (l*l)  
#define kappa (2*pi)

//#define rwt 10	// number of reweighted temperatures


// *** Functions

int 	myrdtsc();
void	memory(int);
void	set_nearest_neighbours();
double	set_energy();
void 	initV(double);
double 	mcstep();
void	initOccFilled(int);
double	set_de(int,int,int);
void	update_filled(int,int);
void 	storer(char,double);
int	set_m2();
void	set_pxy();



// *** Variables

int	model;

int	l;
double	beta;
double	mu;	// chemical potential

int	*occ, *neighbours, *filled;
double	*V;
int 	*mod;
double 	energy;
double 	avem2, avee;

char	*dir;

char 	argL[200]; int argN;

double 	polx=0, poly=0, avep2=0; // if model is not 0 this will never change

//double 	Trw[rwt], averm2num[rwt], averp2num[rwt], averden[rwt];
int rwt;
double 	*Trw, *averm2num, *averp2num, *averden;




// *** Main
int main(int argc, char *argv[]){

	int i,k,j;
	int tr;
	int nm,ns;

	int nTs;
	double T,Tmax,Tmin;

	double cl; // constant for screenin:  lambda=cl*L


	time_t beg, end;
	FILE *opf;

	// starting time
	beg = time(NULL);	


	// input
	argN=13;
	sprintf(argL, "usage: %s <model>\t<l>\t<cl>\t<mu>\t<Tmin>\t<Tmax>\t<nTs>\t<tr>\t<nmesures>\t<nsamples>\t<dir>\t<seed> \n", argv[0]);
        if(argc!=argN) {
                printf("%s",argL);
                exit(1);
        }

	model=(int)atoi(argv[1]);       // 0: neutral case. 1: non-neutral case.        
	l=(int)atoi(argv[2]);           // size    
	cl=(double)atof(argv[3]);	// factor for the screening \lambda=cl*L. If 0: infinite.         
	mu=(double)atof(argv[4]);	// chemical potential           
	Tmin=(double)atof(argv[5]);     // If negative activates the reweighting in a vicinity fabs(Tmin) around Tmax     
	Tmax=(double)atof(argv[6]);          
	nTs=(int)atoi(argv[7]);        	// number of T's (+/- : cooling/heating). If nTs=1 runs only Tmax. 
	tr=(int)atoi(argv[8]);          // relaxation time for each temperature
	nm=(int)atoi(argv[9]);       	// number of measurments at each temperature       
	ns=(int)atoi(argv[10]);         // number of samples (each sample runs all temperatures)
	dir=argv[11];          		// new folder for the outputs    	
	myrand=(unsigned)atoi(argv[12]);// seed for the random numbers generator, 0 means random seed

	// some input-related redefinitions
	if (myrand==0) myrand=myrdtsc(); // random seed from the processor cycles
	if (Tmin<0) {rwt=nTs; nTs=1;} // reweighting mechanism
	strcat(argL, "\t"); for (i=0;i<argN-1;i++) {strcat(argL, argv[i]); strcat(argL, "\t");}  // store "call-1" in argL

	// 

	storer('0',0);	// makes the directory and saves the call

	Init_Random();	// initializing the random number generator
	
	memory(1);	// allocates memory

	set_nearest_neighbours();	// names the neighbors

	cl = (cl>0 ? 1.0/(cl*l) : 0); // trick of the inverse
	cl = (model==2 ? -cl : cl); // trick of model 2 (screening only at k=0) 
	initV(cl);	// builds the potential interaction with screening cl*L (or zero)

	if (Tmin<0) for (i=0;i<rwt;i++) Trw[i]=(Tmax-fabs(Tmin))+i*(2*fabs(Tmin)/rwt); // temperatures for reweighting



////////////////// BYPASS
//		initOccFilled(0);
//		set_pxy();
//		energy=set_energy();
//		beta=(1/0.2115);
//		for (i=0;i<10000;i++) mcstep();
//
//exit(0);
////////////////////

	// loop for the samples
	for (k=0;k<abs(ns);k++){


		initOccFilled(0);
		if (model==0) set_pxy(); // neutral case
		energy=set_energy();


		// loop for the temperatures
		for (T=(nTs>0?Tmin+((Tmax-Tmin)/nTs):Tmax);(T<=Tmax)&&(T>Tmin);T+=(Tmax-Tmin)/nTs){ // cooling or heating depending on the sign of nTs
	
			beta=(T>0?1/T:999999999999.0);
	
			for (i=0;i<tr;i++) mcstep(); // relaxation tr monte carlo steps

			avem2=0; avee=0;	
			avep2=0; // meaningful for neutral case model 0

			if (Tmin<0) for (j=0;j<rwt;j++) {averm2num[j]=0; averp2num[j]=0; averden[j]=0;} // for reweighting 

			for (i=0;i<nm;i++){
		
				//for (j=0;j<tr/nm;j++) mcstep();
				mcstep();

				avem2+=set_m2(); 
				avee+=energy;
				avep2+=(polx*polx+poly*poly); // neutral case 

				if (Tmin<0) for (j=0;j<rwt;j++){	// reweighting
					averm2num[j]+=set_m2()*exp((beta-(1/Trw[j]))*energy);
					averp2num[j]+=(polx*polx+poly*poly)*exp((beta-(1/Trw[j]))*energy);
					averden[j]+=exp((beta-(1/Trw[j]))*energy);
				}
			}

			avem2/=nm; 
			avee/=nm; 

			avep2/=nm; // neutral case
			avep2=1-(kappa*avep2/(2*N*T)); // neutral case: superfluid density

			storer('a',T);

			if (Tmin<0) storer('r',T); // write reweighted sums

			if (k==0) storer('c',T); // snapshot of the last confs of sample 0

		}
		
		storer('A',T);

	}
	



	memory(0);	// free memory


	// writing total time of the run
	end = time(NULL);
	sprintf(argL,"%s/readme.dat",dir);
  	opf = fopen(argL, "a");
    	fprintf(opf,"\n\n  segundos -> %f\n",difftime(end,beg));
    	fprintf(opf,"     horas -> %f\n",difftime(end,beg)/3600);
    	fprintf(opf,"\n");
	fclose(opf);



	return 0;

}


// storing info
void storer(char op, double t){
	FILE *fr;

	char ord[200];

	int i;

	switch (op){

		case '0':
			sprintf(ord, "touch %s; rm -r %s; mkdir %s",dir,dir,dir);
			system(ord);
			sprintf(ord, "mkdir %s/confs",dir);
			system(ord);
			sprintf(ord, "%s/readme.dat",dir);
			fr = fopen (ord, "a"); 
			fprintf(fr,"%s", argL); 
			fprintf(fr,"%i\n", myrand); // before InitRandom !! 	
			fclose(fr);
			
		break;
		case 'a': // temperature - energy - magnetization^2 - superfluid density
			sprintf(ord, "%s/enm2rho.dat",dir); 
			fr = fopen (ord, "a"); 
			fprintf(fr,"%f\t%f\t%f\t%f\n", t, avee/N, avem2, avep2); 
			fclose(fr);
		break;
		case 'A': // double line in file 'a'
			sprintf(ord, "%s/enm2rho.dat",dir);
			fr = fopen (ord, "a"); 
			fprintf(fr,"\n\n"); 
			fclose(fr);
		break;
		case 'r': // rw-temperature - m^2 num - p2num - denominator
			sprintf(ord, "%s/m2p2den.dat",dir); 
			fr = fopen (ord, "a");
 			for (i=0;i<rwt;i++)
				fprintf(fr,"%f\t%f\t%f\t%f\n", Trw[i], averm2num[i], averp2num[i], averden[i]); 
			fprintf(fr,"\n\n"); 
			fclose(fr);
		break;
		case 'c':	
			sprintf(ord, "%s/confs/confT%2.3f.dat",dir,t);
			fr = fopen (ord, "a"); 
			for (i=0;i<N;i++){
                        	if (i%l==0) fprintf(fr,"\n");
                        	fprintf(fr,"%i\t",occ[i]);
                	}
                	fprintf(fr,"\n\n");
			fclose(fr);
		break;
	}


}




// initial configuration
void initOccFilled(int op){

	int i;


	switch (op) {

		case 0 : 	// no particles
 			
			for (i=0;i<N;i++) {occ[i]=0; filled[i]=-1;} 

 			break;

		case 1 :	// two particles

			for (i=0;i<N;i++) {occ[i]=0; filled[i]=-1;} 
 			
			occ[12]=1; occ[22]=-1;
			update_filled(12,22);

 			break;

	}

}




// update filled when occupations are already updated by adding +1 in np and -1 in nm 
// upgrade: to include only one particle 
void update_filled(int np, int nm){
	int i;

	if (np!=-1){
		if (occ[np]==1){
			for (i=0; filled[i]!=-1; i++);
			filled[i]=np;
		}
	
		if (occ[np]==0){
			for (i=0; filled[i]!=np; i++);
			for (i; filled[i]!=-1&&i<(N-1); i++)
				filled[i]=filled[i+1];
			if (i==N-1) filled[i]=-1;
		}
	}	


	if (nm!=-1){
		if (occ[nm]==-1){
			for (i=0; filled[i]!=-1; i++);
			filled[i]=nm;
		}
	
		if (occ[nm]==0){
			for (i=0; filled[i]!=nm; i++);
			for (i; filled[i]!=-1&&i<(N-1); i++)
				filled[i]=filled[i+1];
			if (i==N-1) filled[i]=-1;
		}
	}	

}
	






// *** Montecarlo step (global T, *occ and updated energy)
double mcstep(){
	int ni, nj, ne, i;
	double de, nrate=0.0;
	double dice;

	double dpx, dpy;

	for (i=0;i<N;i++){

		ni=floor(FRANDOM*N);	// positive site
		ne=floor(FRANDOM*Z); 	// neighbor number
		nj=neighbours[Z*ni+ne];	// negative site (neighbor ne of ni)
		
		if (model!=0){ // for the non-neutral case
			dice=FRANDOM;			// just a dice
			if (dice<.25) ni=-1;		// remove the positive charge	
			else if (dice>.75) nj=-1;	// remove the negative charge	
		}	

		de=set_de(ni,nj,ne);

		if ((de<0) || (FRANDOM<exp(-de*beta))){
			if (ni!=-1) occ[ni]+=1; // update occupations
			if (nj!=-1) occ[nj]-=1;	// update occupations
			update_filled(ni,nj);	// update the filled array	
			
			//printf("en %f %f %f %f %i %i %i %i\n",energy, de, energy+de, set_energy(), ni, nj, occ[ni], occ[nj]);

			energy+=de;	// update energy
	
			//printf("en  %f %f %f %i %i %i %i\n",energy, de, set_energy(), ni, nj, occ[ni], occ[nj]);
		

                        if (model==0){ // for the neutral case
				dpx=(ne<2?(ne==1?1:-1):0);
                        	dpy=(ne>1?(ne==2?1:-1):0);
                        	polx+=dpx;      // update polarization (not testable!)
                        	poly+=dpy;
			}



	
			nrate+=1;		


		}

	}
	
	
	return (nrate/N);
}


// change in energy due to adding +1 in np and -1 in nm (which is the neighbor number ne of site np)
double set_de(int np, int nm, int nei){
	int i;
	double qp, qm;
	double 	de_coul=0, de_chem=0;
	int fi,index;

	double dpx, dpy, de_pol=0;

	// coulomb interaction
	qp=(np==-1?0:1);
	qm=(nm==-1?0:-1);
	
	if(qp==1){ 
		for (i=0; fi=filled[i],fi!=-1&&i<N; i++)
			if ((fi!=np)&&(fi!=nm)){
			
				index=abs(np-fi)+((fi<np)?(mod[fi]<=mod[np]?0:l):(mod[np]<=mod[fi]?0:l));
				//index=abs(np-fi)+((fi<np)?((fi%l)<=(np%l)?0:l):((np%l)<=(fi%l)?0:l));

				de_coul+=occ[fi]*qp*V[index];
			}

		de_coul+=((occ[np]*qp)+(qp*qp/2))*V[0];
	}
	if(qm==-1){
		for (i=0; fi=filled[i],fi!=-1&&i<N; i++)
			if ((fi!=np)&&(fi!=nm)){
				
				index=abs(nm-fi)+((fi<nm)?(mod[fi]<=mod[nm]?0:l):(mod[nm]<=mod[fi]?0:l));
				//index=abs(nm-fi)+((fi<nm)?((fi%l)<=(nm%l)?0:l):((nm%l)<=(fi%l)?0:l));

				de_coul+=occ[fi]*qm*V[index];

			}
		de_coul+=((occ[nm]*qm)+(qm*qm/2))*V[0];
	}

	index=abs(nm-np)+((np<nm)?(mod[np]<=mod[nm]?0:l):(mod[nm]<=mod[np]?0:l));
	//index=abs(nm-np)+((np<nm)?((np%l)<=(nm%l)?0:l):((nm%l)<=(np%l)?0:l));
	
	if (qp==1&&qm==-1) de_coul+=(occ[np]*qm+occ[nm]*qp+qp*qm)*V[index];


			


	// chemical potential
	if (qp==1)  de_chem+=fabs(occ[np]+1)-fabs(occ[np]);
	if (qm==-1) de_chem+=fabs(occ[nm]-1)-fabs(occ[nm]);
	de_chem=-mu*de_chem;


        // polarization term
	if (model==0){ // for the neutral case
        	dpx=(nei<2?(nei==1?1:-1):0);
        	dpy=(nei>1?(nei==2?1:-1):0);
        	de_pol=(kappa/(2*N))*((2*polx*dpx)+(dpx*dpx)+(2*poly*dpy)+(dpy*dpy));
	}







	//printf("decoul=%f\n",de_coul);


	return (de_coul+de_chem+de_pol);
}



// *** calculates the energy (global *occ and *cfield)
double set_energy(){

	double 	e_coul=0, e_chem=0;
	int i,j;
	int fi,fj,index;


	// coulomb interaction
	for (i=0; fi=filled[i], fi!=-1&&i<N; i++)
		for (j=0; fj=filled[j],fj!=-1&&j<N; j++){

			//index=abs(fj-fi)+((fi<fj)?((fi%l)<=(fj%l)?0:l):((fj%l)<=(fi%l)?0:l));
			index=abs(fj-fi)+((fi<fj)?(mod[fi]<=mod[fj]?0:l):(mod[fj]<=mod[fi]?0:l));

			e_coul+=.5*occ[fi]*occ[fj]*V[index];

		}



	//printf("%f %f %f\n",polx,poly,e_pol);

	// chemical potential
	for (i=0;i<N;i++) e_chem+=fabs(occ[i]);
	e_chem=-mu*e_chem;


	//printf("en=%f\n",e_coul+e_chem);

	return (e_coul+e_chem);
	//return (e_pol);
}





// *** calculates the squared polarization
int set_m2(){
	int i, cnt_m=0;

	for (i=0; filled[i]!=-1&&i<N; i++)
		cnt_m+=occ[filled[i]];
	
	return (cnt_m*cnt_m);
}



// *** calculates the squared polarization
void set_pxy(){
        int i;
        double  px=0, py=0;
        int     rx, ry;

        for (i=0;i<N;i++){
                rx=i%l;     // can use mod[]
                ry=floor(i/l);
                px+=occ[i]*rx;
                py+=occ[i]*ry;
        }

        polx=px;
        poly=py;

}








// *** generates the neighbours array (global l)
void set_nearest_neighbours()
{
	int i;
	for (i = 0; i < N; i++){	// not the optimal way but happens only once 
		// rigth               
		if (div(i+1, l).rem == 0) neighbours[(Z * i)] = i - (l -1);
		else neighbours[(Z * i)] = i + 1;
		// left
		if (div(i, l).rem == 0) neighbours[(Z * i) + 1] = i + (l -1);
		else neighbours[(Z * i) + 1] = i -1;
		// up
		if (i < l)  neighbours[(Z * i) + 2] = i + ((l -1) * l);
		else neighbours[(Z * i) + 2] = i - l;
		// down
		if ((i + l) >= N) neighbours[(Z * i) + 3] = i -((l - 1) * l);
		else neighbours[(Z * i) + 3] = i + l;
    	}
}



// *** memory (global N)
void memory(int order){

	int i;

	if (order==1){
		occ = (int *) malloc(N*sizeof(int));		
		neighbours = (int *) malloc(N*Z*sizeof(int));		
		filled = (int *) malloc(N*sizeof(int));		

		V=(double *) malloc(N*sizeof(double));
		mod=(int *) malloc(N*sizeof(int));
		for (i=0;i<N;i++) mod[i]=i%l; // definition

		Trw=(double *) malloc(rwt*sizeof(double));
		averp2num=(double *) malloc(rwt*sizeof(double));
		averm2num=(double *) malloc(rwt*sizeof(double));
		averden=(double *) malloc(rwt*sizeof(double));
	}
	else{
		free(occ);
		free(neighbours);
		free(filled);

		free(V);
		free(mod);

		free(Trw);
		free(averm2num);
		free(averp2num);
		free(averden);
	}

}




// returns the last five digits of the total pseudo-cycles since the processor was powered on
int myrdtsc(){
	unsigned long long res;
	unsigned int lo,hi;

	__asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
	res = ((unsigned long long)hi << 32) | lo;

	return (abs(res)%100000);
}








// builts the 2D potential interaction, input: the inverse of the screening lenght 
void initV(double invl){

        //double kappa=2*pi;

        int i, j;
        int dx, dy;

        int nx, ny;
        double kx, ky;


        double **vv;
        int fl2;

        fl2=floor(l/2);

        vv=(double **) malloc(l*sizeof(double *));
                for(i=0;i<l;i++)
                        vv[i]=(double *) malloc(l*sizeof(double));

        for (i=0;i<l;i++)
                for (j=i;j<l;j++){

                        vv[i][j]=0;

                        dx=i-fl2;
                        dy=j-fl2;

                        for (nx=0;nx<l;nx++)
                                for (ny=0;ny<l;ny++){

                                        kx=(2*pi*nx/l);
                                        ky=(2*pi*ny/l);

                                        if (kx!=0 || ky!=0 || invl>0) 
                                                vv[i][j] += cos((kx*dx)+(ky*dy))/(4-(2*cos(kx))-(2*cos(ky))+(invl*invl));
					else if (invl<0) // model 2, with screening only at k=0
                                                vv[i][j] += cos((kx*dx)+(ky*dy))/(4-(2*cos(kx))-(2*cos(ky))+(invl*invl));
                                }
                        vv[i][j]*=kappa/N;

                }


                
	for (j=0;j<N;j++){

      		dx=j%l;
    		dy=floor(j/l);

          	dx=(abs(dx)<l/2 ? dx : -(l-abs(dx)));
          	dy=(abs(dy)<l/2 ? dy : -(l-abs(dy)));

       	        //------
           	dx=dx+fl2;
           	dy=dy+fl2;

        	if (dx<=dy) V[j]=vv[dx][dy];
             	else V[j]=vv[dy][dx];
  	}



        for(i=0;i<l;i++) free(vv[i]);
        free(vv);

}







