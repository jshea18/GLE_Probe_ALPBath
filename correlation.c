/******************************************************************************************************************************************
 * Program:
 * - calculate correlation functions
 * - calculate memory kernels
 ******************************************************************************************************************************************/

//include
#include <fstream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <cstring>
#include <iostream>
#include <stdio.h>
#include <unistd.h>

//read in variables
char *inputFile;
char *meanFile;
int colloids;
int columns;
double dt;
int Ncor;
int Neq;
int N;
int Npart;
int Nav;
int fCG;
int dim;
double M;
double kharm;
std::vector< std::vector<int> > cor_code;
std::vector<char*> cor_filename;
int forward_flag = 0;
std::vector<int> forward_code;
char *forward_filename;
char *forward_filename_2nd;
int backward_flag = 0;
std::vector<int> backward_code;
char *backward_filename;
char *backward_filename_2nd;
int volterra_flag = 0;
int volterraharm_flag = 0;
char *volterra_filename;
char *joint_dist;
char *data_proj_file;
char *dist_proj_file;
char *dist_dot_file;
char *corr_dot_file;
char *vvFile;
char *vfFile;
char *ffFile;
char *sfFile;
char *jd1File;
char *jd2File;
char *sfxFile;
char *sfyFile;
char *sfzFile;
char *msdFile;
int nbin = 5000;
//char *msd_filename;
std::vector<double> data_memory_volterra;
int non_eq, stochastic_force, force_index, linear_response_project, msd, joint, noise;

//global data
std::vector<double> data;			//vector with position and velocity data

//functions
int print_digit(int input, int digit);
void read();
void correlation(int beginA, int beginB, int dimension, const char* outFile);
void memory_volterra();
void memory_volterra_harmonic();
void memory_forward(int beginA, int beginB, int dimensionO, int beginP, int beginF, int dimensionR );
void memory_forwardsecond(int beginA, int beginB, int dimensionO, int beginP, int beginF, int dimensionR );
void memory_backward(int beginA, int beginB, int dimensionO, int beginP, int beginF, int dimensionR );
void memory_backwardsecond(int beginA, int beginB, int dimensionO, int beginP, int beginF, int dimensionR );
void calc_stochastic_force();
void calc_linear_response();
void calc_msd();
void calc_joint_dist();
void calc_proj_noise();

// how to use program
static void show_usage(char* argv[]){
  std::cout << "Usage: " << argv[0] 
    << " inputFile"
    << " columns"
    << " N"
    << " dt"
    << " Neq"
    << " Ncor"
    << " fCG"
    << std::endl;
}

// read in parameters from command line
void read_command(int argc, char* argv[])
{
  if(argc < 8){ 
    show_usage(argv);
    exit(1);
  }
  inputFile = argv[1];
  columns = atoi(argv[2]);
  N = atoi(argv[3]);
  dt = atof(argv[4]);
  Neq = atoi(argv[5]);
  Ncor = atoi(argv[6]);
  fCG = atoi(argv[7]);
  dim = atoi(argv[8]);
  M = atof(argv[9]);
  Npart = atof(argv[10]);
  meanFile = argv[11];
  std::cout << "Parameters:"  
    << "\n inputFile " << inputFile
    << "\n columns " << columns
    << "\n N " << N
    << "\n dt " << dt
    << "\n Neq " << Neq
    << "\n Ncor " << Ncor
    << "\n fCG " << fCG
    << "\n dim " << dim
    << "\n mass " << M
    << std::endl;	

  int opt;
  
  non_eq = 0;
  noise = 0;
  stochastic_force = 0;
  msd = 0;
  joint = 0;

  while ((opt = getopt(argc, argv, "c:f:b:v:ns:l:m:j:p:h")) != -1) {
    std::vector<int> cor_code_local;
    switch (opt) {
      case 'c': 
        for(int i=0; i<3; i++) cor_code_local.push_back(atoi(argv[optind-1+i]));
        cor_code.push_back(cor_code_local);
        cor_filename.push_back(argv[optind+2]);
        break;
      case 'f': 
        forward_flag = 1;
        for(int i=0; i<6; i++) forward_code.push_back(atoi(argv[optind-1+i]));
        forward_filename = argv[optind+5];
        forward_filename_2nd = argv[optind+6];
        break;
      case 'b': 
        backward_flag = 1;
        for(int i=0; i<6; i++) backward_code.push_back(atoi(argv[optind-1+i]));
        backward_filename = argv[optind+5];
        backward_filename_2nd = argv[optind+6];
        break;
      case 'v': 
        volterra_flag = 1;
        // filenames for the correlation functions of the relevant variable and first and second derivative
        vvFile = optarg;
printf("%s \n", vvFile);
        vfFile = argv[optind];
        ffFile = argv[optind+1];
        volterra_filename = argv[optind+2];
        break;
      case 'h': 
        volterraharm_flag = 1;
        // filenames for the correlation functions of the relevant variable and first and second derivative
        vvFile = argv[optind];
printf("%s \n", vvFile);
        vfFile = argv[optind+1];
printf("%s \n", vfFile);
        ffFile = argv[optind+2];
printf("%s \n", ffFile);
        volterra_filename = argv[optind+3];
printf("%s \n", volterra_filename);
        kharm = atof(argv[optind+4]);
printf("%f \n", kharm);
        break;
      case 'n': 
        non_eq = 1;
        break;
      case 's': 
        stochastic_force = 1;
        sfFile = optarg;
        force_index = atoi(argv[optind]);
        printf(" %d\n", force_index);
        break;
      case 'l': 
        linear_response_project = 1;
        break;
      case 'm': 
        msd = 1;
        msdFile = optarg;
        break;
      case 'j': 
        joint = 1;
        jd1File = argv[optind];
        jd2File = argv[optind+1];
        joint_dist = argv[optind+2];
        break;
      case 'p': 
        noise = 1;
        sfxFile = argv[optind];
        sfyFile = argv[optind+1];
        sfzFile = argv[optind+2];
        data_proj_file = argv[optind+3];
        dist_proj_file = argv[optind+4];
        dist_dot_file = argv[optind+5];
        corr_dot_file = argv[optind+6];
        break;
    }
  }

  std::cout << "Correlation functions:";
  for (int c=0; c < cor_code.size(); c++) {
    std::cout << "\n "; 
    for(int i=0; i<3; i++)
      std::cout << cor_code[c][i]; 
      std::cout << " " << cor_filename[c];
    }
    std::cout << std::endl;

  std::cout << "Flags:"  
    << "\n forward " << forward_flag
    << "\n backward " << backward_flag
    << "\n volterra " << volterra_flag
    << "\n non-eq " << non_eq
    << "\n stochastic_force " << stochastic_force
    << "\n linear_response_project " << linear_response_project
    << "\n msd " << msd
    << "\n joint " << joint
    << "\n average noise squared " << noise
    << std::endl;

}



/*****************************************************
 * Main-Method
 * ***************************************************/
int main(int argc, char* argv[]){
  read_command(argc, argv);
  
  N /= fCG;
  Neq /= fCG;
  Ncor /= fCG;
  dt *= fCG;

  printf("New version \n");
  printf("Reading ... \n");
  read();
  printf("Reading ... DONE! \n");

  for (int c=0; c < cor_code.size(); c++) {
    printf("Calculating correlation function %d/%d ... \n",c+1,cor_code.size());
    correlation(cor_code[c][0],cor_code[c][1],cor_code[c][2],cor_filename[c]);
    printf("Calculating correlation function ... DONE! \n");
  }

  if(forward_flag) {
    printf("Calculating memory using forward orthogonal dynamics ...\n");
    memory_forward(forward_code[0],forward_code[1],forward_code[2],forward_code[3],forward_code[4],forward_code[5]);
    printf("Calculating memory using forward orthogonal dynamics ... DONE\n");

    printf("Calculating memory using forward orthogonal dynamics (second order) ...\n");
    memory_forwardsecond(forward_code[0],forward_code[1],forward_code[2],forward_code[3],forward_code[4],forward_code[5]);
    printf("Calculating memory using forward orthogonal dynamics (second order) ... DONE\n");
  }

  if(backward_flag) {
    printf("Calculating memory using backward orthogonal dynamics ...\n");
    memory_backward(backward_code[0],backward_code[1],backward_code[2],backward_code[3],backward_code[4],backward_code[5]);
    printf("Calculating memory using backward orthogonal dynamics ... DONE\n");
    
    printf("Calculating memory using backward orthogonal dynamics (second order) ...\n");
    memory_backwardsecond(backward_code[0],backward_code[1],backward_code[2],backward_code[3],backward_code[4],backward_code[5]);
    printf("Calculating memory using backward orthogonal dynamics (second order) ... DONE\n");
  }

  if(volterra_flag) {
    printf("Calculating memory using inverse volterra ...\n");
    memory_volterra();
    printf("Calculating memory using inverse volterra ... DONE\n");
  }
  
  if(volterraharm_flag) {
    printf("Calculating memory using inverse volterra with harmonic potential ...\n");
    memory_volterra_harmonic();
    printf("Calculating memory using inverse volterra with harmonic potential ... DONE\n");
  }
  
  if(stochastic_force) {
    printf("Calculating stochastic force ...\n");
    calc_stochastic_force();
    printf("Calculating stochastic force ... DONE\n");
  }
  
  if(linear_response_project) {
    printf("Calculating linear response project ...\n");
    calc_linear_response();
    printf("Calculating linear response project ... DONE\n");
  }

  if(msd) {
    printf("Calculating msd ...\n");
    calc_msd();
    printf("Calculating msd ... DONE\n");
  }

  if(joint) {
    printf("Calculating joint distributioin function ...\n");
    calc_joint_dist();
    printf("Calculating joint distribution function ... DONE\n");
  }

  if(noise) {
    printf("Calculating average noise squared ...\n");
    calc_proj_noise();
    printf("Calculating average noise squared ... DONE\n");
  }

}

/********************************************************
 * reads in position,momenta and force from md-simulation
 * ******************************************************/
void read(){
  
  double x,y,z,vx,vy,vz,fx,fy,fz,x2,y2,z2,vx2,vy2,vz2,fx2,fy2,fz2;
  data.resize(N*columns);
  
  FILE * input = fopen(inputFile, "r");
  
  char buffer[100];
  fgets(buffer, 100, input);
  
  double meanV=0.,meanF=0.,meanV2=0.,meanVx=0.,meanVy=0.,meanVz=0.,meanVx2=0.,meanVy2=0.,meanVz2=0.;
  double stdV,stdVx,stdVy,stdVz;
  double meanV_2=0.,meanF_2=0.,meanV2_2=0.,meanVx_2=0.,meanVy_2=0.,meanVz_2=0.,meanVx2_2=0.,meanVy2_2=0.,meanVz2_2=0.;
  double stdV_2,stdVx_2,stdVy_2,stdVz_2;
  
  for (int i=0; i<N*fCG; i++) {
    if (Npart==1) {
      fscanf(input,"%lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&x,&y,&z,&vx,&vy,&vz,&fx,&fy,&fz);
      if (i%fCG==0) {
        data[columns*i/fCG]=x;
        data[columns*i/fCG+1]=y;
        data[columns*i/fCG+2]=z;
        data[columns*i/fCG+3]=vx;
        //if (i> Neq*fCG) meanV+=vx;
        if (i> Neq*fCG) meanVx+=vx;
        if (i> Neq*fCG) meanVx2+=vx*vx;
        data[columns*i/fCG+4]=vy;
        if (i> Neq*fCG) meanVy+=vy;
        if (i> Neq*fCG) meanVy2+=vy*vy;
        data[columns*i/fCG+5]=vz;
        if ((i> Neq*fCG) && dim==3){
          meanVz+=vz;
          meanVz2+=vz*vz;
          meanV+=sqrt(vx*vx+vy*vy+vz*vz);
          meanV2+=vx*vx+vy*vy+vz*vz;
        }
        if ((i> Neq*fCG) && dim==2){ 
          meanV+=sqrt(vx*vx+vy*vy);
          meanV2+=vx*vx+vy*vy;
        }
        data[columns*i/fCG+6]=fx;
        if (i> Neq*fCG) meanF+=fx;
        data[columns*i/fCG+7]=fy;
        data[columns*i/fCG+8]=fz;
      }
    }
        
    if (Npart==2) {
      fscanf(input,"%lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf\n",&x,&y,&z,&vx,&vy,&vz,&fx,&fy,&fz,&x2,&y2,&z2,&vx2,&vy2,&vz2,&fx2,&fy2,&fz2);
      if (i%fCG==0) {
        data[columns*i/fCG]=x;
        data[columns*i/fCG+9]=x2;
        data[columns*i/fCG+1]=y;
        data[columns*i/fCG+10]=y2;
        data[columns*i/fCG+2]=z;
        data[columns*i/fCG+11]=z2;
        data[columns*i/fCG+3]=vx;
        data[columns*i/fCG+12]=vx2;
        //if (i> Neq*fCG) meanV+=vx;
        meanVx+=vx;
        meanVx_2+=vx2;
        meanVx2+=vx*vx;
        meanVx2_2+=vx2*vx2;
        data[columns*i/fCG+4]=vy;
        data[columns*i/fCG+13]=vy2;
        meanVy+=vy;
        meanVy_2+=vy2;
        meanVy2+=vy*vy;
        meanVy2_2+=vy2*vy2;
        data[columns*i/fCG+5]=vz;
        data[columns*i/fCG+14]=vz2;
        if ((i> Neq*fCG) && dim==3){
          meanVz+=vz;
          meanVz_2+=vz2;
          meanVz2+=vz*vz;
          meanVz2_2+=vz2*vz2;
          meanV+=sqrt(vx*vx+vy*vy+vz*vz);
          meanV2+=vx*vx+vy*vy+vz*vz;
          meanV_2+=sqrt(vx2*vx2+vy2*vy2+vz2*vz2);
          meanV2_2+=vx2*vx2+vy2*vy2+vz2*vz2;

        }
        if ((i> Neq*fCG) && dim==2){ 
          meanV+=sqrt(vx*vx+vy*vy);
          meanV2+=vx*vx+vy*vy;
          meanV_2+=sqrt(vx2*vx2+vy2*vy2);
          meanV2_2+=vx2*vx2+vy2*vy2;
        }
        data[columns*i/fCG+6]=fx;
        data[columns*i/fCG+15]=fx2;        
        if (i> Neq*fCG) meanF+=fx;        
        if (i> Neq*fCG) meanF_2+=fx2;
        data[columns*i/fCG+7]=fy;
        data[columns*i/fCG+16]=fy2;
        data[columns*i/fCG+8]=fz;
        data[columns*i/fCG+17]=fz2;
      }
    }
      //printf("%d %f\n",i,y);
  }
  
  // calculate and print mean for non_eq simulations

  stdV=sqrt((meanV2-meanV*meanV/(1.0*(N-Neq)))/(1.0*(N-Neq)));
  meanV/=(1.0*(N-Neq));
  meanF/=(1.0*(N-Neq));
  meanV2/=(1.0*(N-Neq));
  meanV_2/=(1.0*(N-Neq));
  meanF_2/=(1.0*(N-Neq));
  meanV2_2/=(1.0*(N-Neq));
  stdV_2=sqrt(meanV2_2-meanV_2*meanV_2);

  stdVx=sqrt((meanVx2-meanVx*meanVx/(1.0*(N-Neq)))/(1.0*(N-Neq)));    
  meanVx/=(1.0*(N-Neq));
  meanVx2/=(1.0*(N-Neq));
  meanVx_2/=(1.0*(N-Neq));
  meanVx2_2/=(1.0*(N-Neq));
  stdVx_2=sqrt(meanVx2_2-meanVx_2*meanVx_2);

  stdVy=sqrt((meanVy2-meanVy*meanVy/(1.0*(N-Neq)))/(1.0*(N-Neq)));
  meanVy/=(1.0*(N-Neq));
  meanVy2/=(1.0*(N-Neq));
  meanVy_2/=(1.0*(N-Neq));
  meanVy2_2/=(1.0*(N-Neq));
  stdVy_2=sqrt(meanVy2_2-meanVy_2*meanVy_2);

  if (dim==3) {
    stdVz=sqrt((meanVz2-meanVz*meanVz/(1.0*(N-Neq)))/(1.0*(N-Neq)));  
    meanVz/=(1.0*(N-Neq));
    meanVz2/=(1.0*(N-Neq));
    meanVz_2/=(1.0*(N-Neq));
    meanVz2_2/=(1.0*(N-Neq));
    stdVz_2=sqrt(meanVz2_2-meanVz_2*meanVz_2);
  }

  if (dim==2) stdVz=0.,stdVz_2=0.;
        
  FILE * in;
  in = fopen(meanFile,"w");
  if (Npart==1) fprintf(in,"%lf \n %lf %lf \n %lf %lf \n %lf %lf \n",meanV2,meanVx,stdVx,meanVy,stdVy,meanVz,stdVz);
  if (Npart==2) fprintf(in,"%lf %lf \n %lf %lf \n %lf %lf \n %lf %lf \n \n \n",meanV,stdV,meanVx,stdVx,meanVy,stdVy,meanVz,stdVz,meanV_2,stdV_2,meanVx_2,stdVx_2,meanVy_2,stdVy_2,meanVz_2,stdVz_2);
  fclose(in);
    
  if (non_eq) {    
    for (int i=0; i<N; i++) {
      data[columns*i+3]-=meanVx;
      //data[columns*i+6]-=meanF;
    }
  }
  
  // substract mean from correlation functions
  

  Nav = N-Ncor;

}

/***********************************************************
 * calculate an autocorrelation function <B(t)A>
 * int beginA: index of first column to include for A
 * int beginB: index of first column to include for B
 * int dimension: dimension of the vectors A and B
 * ********************************************************/
void correlation(int beginA, int beginB, int dimension, const char* outFile){

  int n,m,d;
  double result = 0;
  double error = 0;

  std::vector<double> Azero(dimension*N);
  std::vector<double> cor(Ncor);
  std::vector<double> cor_error(Ncor);
	
  n=0;
  // special case n=0 

  for (m=Neq; m<N; m++) {
    for (d=0; d<dimension; d++){
      Azero[dimension*m+d]=data[columns*m+d+beginA];
      //printf("partial: %f\n",Azero[dimension*m+d]);
    }
    double partial_result = 0;
    for (d=0; d<dimension; d++){
      partial_result+=Azero[dimension*m+d]*data[columns*m+d+beginB];
    }
    result += partial_result;
    error += partial_result*partial_result;
  }
  cor[0] = result/(N-Neq);
  //include correlation (about 2*dt steps) into error calculation
  cor_error[0] = sqrt(error/(N-Neq)-cor[0]*cor[0])/sqrt((N-Neq)/2*dt);

  for (n=1; n<Ncor; n++) {
    result = 0;
    error = 0;
    for (m=Neq;m<N-n;m++) {
      double partial_result = 0;
      for (d=0; d<dimension; d++){
	partial_result+=Azero[dimension*m+d]*data[columns*m+columns*n+d+beginB];
      }
      result += partial_result;
      error += partial_result*partial_result;
    }

    cor[n]=result/(N-Neq-n);
    cor_error[n] = sqrt(error/(N-Neq)-cor[n]*cor[n])/sqrt((N-Neq)/2*dt);

  }
	
  FILE * out;
  out = fopen(outFile, "w");
  for(n=0; n<Ncor; n++)
    fprintf(out, "%f %f %f\n",n*dt, cor[n]/dimension, cor_error[n]/dimension);
  fclose(out); 
}

/***********************************************************
 * calculate memory function using inverse volterra method
 * ref: Shin et. al, Chem.Phys. 375, 316 (2010) 
 * ********************************************************/
void memory_volterra(){

  double vv0;
  std::vector<double> vf;	
  std::vector<double> ff;
  data_memory_volterra.resize(Ncor);
  double dummy_t,dummy_v,dummy_e;

  //read in
  FILE * in;
  in = fopen(vvFile,"r");
  fscanf(in,"%lf %lf %lf\n",&dummy_t,&vv0,&dummy_e);
  fclose(in);
  vv0*=M*M;

  in = fopen(vfFile,"r");
  for(int i=0;i<Ncor;i++){
    fscanf(in,"%lf %lf %lf\n",&dummy_t,&dummy_v,&dummy_e);
    vf.push_back(dummy_v*M);
  }
  fclose(in);

  in = fopen(ffFile,"r");
  for(int i=0;i<Ncor;i++){
    fscanf(in,"%lf %lf %lf\n",&dummy_t,&dummy_v,&dummy_e);
    ff.push_back(dummy_v);
  }
  fclose(in);

  //calc memory
  data_memory_volterra[0] = ff[0]/vv0;
  for(int i=1;i<Ncor;i++){
    double inv = 1/(vv0+dt*0.5*vf[0]);
    double sum=0;
    for(int j=0;j<i;j++){
      double wj = (j==0) ? 0.5 : 1.0;
      sum += wj*vf[i-j]*data_memory_volterra[j];
    }
    data_memory_volterra[i]=inv*(ff[i]-dt*sum);
  }

  //print memory
  int n;
  FILE * out;
  out = fopen(volterra_filename, "w");
  for(n=0; n<Ncor; n++) {
    data_memory_volterra[n]=data_memory_volterra[n]*M;
    fprintf(out, "%f %f\n",n*dt,data_memory_volterra[n]);
  }
  fclose(out); 

}

/***********************************************************
 * calculate memory function using inverse volterra method with harmonic potential
 * ref: Shin et. al, Chem.Phys. 375, 316 (2010) 
 * ********************************************************/
void memory_volterra_harmonic(){

  double vv0;
  std::vector<double> vv;
  std::vector<double> vf;	
  std::vector<double> ff;
  data_memory_volterra.resize(Ncor);
  double dummy_t,dummy_v,dummy_e;

  //read in
  FILE * in;
  in = fopen(vvFile,"r");
  for(int i=0;i<Ncor;i++){
    fscanf(in,"%lf %lf %lf\n",&dummy_t,&dummy_v,&dummy_e);
    vv.push_back(dummy_v*M*M);
  }
  fclose(in);
  vv0=vv[0];

  in = fopen(vfFile,"r");
  for(int i=0;i<Ncor;i++){
    fscanf(in,"%lf %lf %lf\n",&dummy_t,&dummy_v,&dummy_e);
    vf.push_back(dummy_v*M);
  }
  fclose(in);

  in = fopen(ffFile,"r");
  for(int i=0;i<Ncor;i++){
    fscanf(in,"%lf %lf %lf\n",&dummy_t,&dummy_v,&dummy_e);
    ff.push_back(dummy_v);
  }
  fclose(in);

  //calc memory
  data_memory_volterra[0] = ff[0]/vv0-kharm/M;
  for(int i=1;i<Ncor;i++){
    double inv = 1/(vv0+dt*0.5*vf[0]);
    double sum=0;
    for(int j=0;j<i;j++){
      double wj = (j==0) ? 0.5 : 1.0;
      sum += wj*vf[i-j]*data_memory_volterra[j];
    }
    data_memory_volterra[i]=inv*(ff[i]-kharm/M*vv[i]-dt*sum);
  }

  //print memory
  int n;
  FILE * out;
  out = fopen(volterra_filename, "w");
  for(n=0; n<Ncor; n++) {
    data_memory_volterra[n]=data_memory_volterra[n]*M;
    fprintf(out, "%f %f\n",n*dt,data_memory_volterra[n]);
  }
  fclose(out); 

}

/*************************************************************************************
 * calculate memory function using forward orthogonal dynamics
 * ref: Rotenberg et. al, J. Chem.Phys. 140, 124103 (2014) 
 * int beginA: index of first column to include for A
 * int beginB: index of first column to include for B
 * int dimensionO: dimension of the vectors A and B
 * int beginP: index of first column to include for relevant variable
 * int beginF: index of first column to include for derivative of relevant variable
 * int dimensionR: dimension of the relevant variable
 * ***********************************************************************************/
void memory_forward(int beginA, int beginB, int dimensionO, int beginP, int beginF, int dimensionR ){

  int n,m,dO,dR;
  double result;

  std::vector<double> Bplus(dimensionO*N,0);
  std::vector<double> memory(Ncor);
	
  //printf("beginA=%d beginB=%d dimension0=%d beginP=%d beginF=%d dimensionR=%d\n",beginA, beginB, dimensionO, beginP, beginF, dimensionR);

  n=0;
  // special case n=0
  result=0;
  for (m=Neq;m<N-Ncor;m++) {
    double partial_result = 0;
    for (dO=0; dO<dimensionO; dO++){
      //determine new Bplus
      Bplus[dimensionO*m+dO]=data[columns*(N-m)+dO+beginB];
      partial_result += data[columns*(N-m)+dO+beginA]*Bplus[dimensionO*m+dO];
    }
    //calculate correlation function
    result+=partial_result;
  }
  memory[0]=result/(N-Neq-Ncor);

  
  for(n=1;n<Ncor;n++){ //main-loop for t=n*dt
    // calculate beta matrix to evaluate orthogonal dynamics
    std::vector<double> beta(dimensionR*dimensionO,0);
    std::vector<double> norm(dimensionR,0);
    for (m=Neq;m<N-Ncor;m++) {
      for (dR=0; dR<dimensionR; dR++){
	norm[dR] += data[columns*(N-m-n)+dR+beginP]*data[columns*(N-m-n)+dR+beginP];
      }
    }
    double Mf = 3.0/(norm[0]+norm[1]+norm[2])*(N-Ncor-Neq);

    for (m=Neq;m<N-Ncor;m++) {
      for(dO=0; dO<dimensionO; dO++) {
	for (dR=0; dR<dimensionR; dR++){
	  beta[dO*dimensionO+dR] += data[columns*(N-m-n)+dR+beginF]*Bplus[dimensionO*m+dO];
	}
      }
    }

    for(dO=0; dO<dimensionO; dO++) {
      for (dR=0; dR<dimensionR; dR++){
	beta[dO*dimensionO+dR] /= norm[dR]*Mf*Mf;
      }
    }

    result = 0;

    for (m=Neq;m<N-Ncor;m++) {
      double partial_result = 0;
      for(dO=0; dO<dimensionO; dO++) {
	for (dR=0; dR<dimensionR; dR++){
	  Bplus[dimensionO*m+dO]=Bplus[dimensionO*m+dO]+beta[dO*dimensionO+dR]*data[columns*(N-m-n)+dR+beginP]*dt*Mf;
	}
	partial_result += Bplus[dimensionO*m+dO]*data[columns*(N-m-n)+dO+beginA];
      }
    result += partial_result;
    }
  
    memory[n]=result/(N-Ncor-Neq);
  }


  FILE * out;
  out = fopen(forward_filename, "w");
	
  for(int n=0; n<Ncor; n++)
    fprintf(out, "%f %f\n", n*dt, memory[n]);
  fclose(out); 

}

/*************************************************************************************
 * calculate memory function using forward orthogonal dynamics (second order)
 * ref: Rotenberg et. al, (submitted)
 * int beginA: index of first column to include for A
 * int beginB: index of first column to include for B
 * int dimensionO: dimension of the vectors A and B
 * int beginP: index of first column to include for relevant variable
 * int beginF: index of first column to include for derivative of relevant variable
 * int dimensionR: dimension of the relevant variable
 * ***********************************************************************************/
void memory_forwardsecond(int beginA, int beginB, int dimensionO, int beginP, int beginF, int dimensionR ){

  int n,m,dO,dR;
  double result;

  std::vector<double> Bplus(dimensionO*N,0);
  std::vector<double> memory(Ncor);

  n=0;
  // special case n=0
  result=0;
  for (m=Neq;m<N-Ncor;m++) {
    double partial_result = 0;
    for (dO=0; dO<dimensionO; dO++){
      //determine new Bplus
      Bplus[dimensionO*m+dO]=data[columns*(N-m)+dO+beginB];
      partial_result += data[columns*(N-m)+dO+beginA]*Bplus[dimensionO*m+dO];
    }
    //calculate correlation function
    result+=partial_result;
  }
  memory[0]=result/(N-Neq-Ncor);

  for(n=1;n<Ncor;n++){ //main-loop for t=n*dt

    // calculate beta, gamma matrix and kappa and delta to evaluate orthogonal dynamics
    std::vector<double> beta(dimensionR*dimensionO,0);
    std::vector<double> gamma(dimensionR*dimensionO,0);
    std::vector<double> norm(dimensionR,0);
    double kappa = 0;
    double delta = 0;
    for (m=Neq;m<N-Ncor;m++) {
      for (dR=0; dR<dimensionR; dR++){
	norm[dR] += data[columns*(N-m-n)+dR+beginP]*data[columns*(N-m-n)+dR+beginP];
	kappa += data[columns*(N-m-n)+dR+beginF]*data[columns*(N-m-n)+dR+beginP];
	delta += data[columns*(N-m-n)+dR+beginF]*data[columns*(N-m-1-n)+dR+beginP];
      }
    }

    double Mf2 = 3.0/(norm[0]+norm[1]+norm[2])*(N-Ncor-Neq);
    
    for (m=Neq;m<N-Ncor;m++) {
      for(dO=0; dO<dimensionO; dO++) {
	for (dR=0; dR<dimensionR; dR++){
	  beta[dO*dimensionO+dR] += data[columns*(N-m-n)+dR+beginF]*Bplus[dimensionO*m+dO];
	  gamma[dO*dimensionO+dR] += data[columns*(N-m-n)+dR+beginF]*Bplus[dimensionO*(m+1)+dO];
	}
      }
    }

    for(dO=0; dO<dimensionO; dO++) {
      for (dR=0; dR<dimensionR; dR++){
	beta[dO*dimensionO+dR] /= norm[dR]*Mf2*Mf2;
	gamma[dO*dimensionO+dR] /= norm[dR]*Mf2*Mf2;
      }
    }
    kappa /= (norm[0]+norm[1]+norm[2])*Mf2*Mf2;
    delta /= (norm[0]+norm[1]+norm[2])*Mf2*Mf2;

    result = 0;

    for (m=Neq;m<N-Ncor;m++) {
      double partial_result = 0;
      for(dO=0; dO<dimensionO; dO++) {
	for (dR=0; dR<dimensionR; dR++){
	  Bplus[dimensionO*m+dO]=Bplus[dimensionO*m+dO]+beta[dO*dimensionO+dR]*data[columns*(N-m-n)+dR+beginP]*Mf2*dt/2.0+dt/2.0/(1-dt/2.0*kappa*Mf2)*data[columns*(N-m-n+1)+dR+beginP]*Mf2*(gamma[dO*dimensionO+dR]+dt/2.0*delta*Mf2*beta[dO*dimensionO+dR]);
	}
	partial_result += Bplus[dimensionO*m+dO]*data[columns*(N-m-n)+dO+beginA];
      }
      result += partial_result;
    }
    memory[n]=result/(N-Ncor-Neq);
  }

  FILE * out;
  out = fopen(forward_filename_2nd, "w");

  for(int n=0; n<Ncor; n++)
    fprintf(out, "%f %f\n", n*dt, memory[n]);
  fclose(out); 

}

/**************************************************************
 * calculate memory function using backward orthogonal dynamics
 * ref: Rotenberg et. al, J. Chem.Phys. 140, 124103 (2014) 
 * int beginA: index of first column to include for A
 * int beginB: index of first column to include for B
 * int dimensionO: dimension of the vectors A and B
 * int beginP: index of first column to include for relevant variable
 * int beginF: index of first column to include for derivative of relevant variable
 * int dimensionR: dimension of the relevant variable
 * ***********************************************************/
void memory_backward(int beginA, int beginB, int dimensionO, int beginP, int beginF, int dimensionR ){

  int n,m,dO,dR;
  double result;
	
  std::vector<double> Atilde(dimensionO*Nav,0);
  std::vector<double> memory(Ncor);

  n=0;
  // special case n=0
  result=0;
  for (m=Neq;m<Nav;m++) {
    double partial_result = 0;
    for (dO=0; dO<dimensionO; dO++){
      //determine new Bplus
      Atilde[dimensionO*m+dO]=data[columns*m+dO+beginA];
      partial_result += data[columns*m+dO+beginB]*Atilde[dimensionO*m+dO];
    }
    //calculate correlation function
    result+=partial_result;
  }
  memory[0]=result/(Nav-Neq);

  for(n=1;n<Ncor;n++){ //main-loop for t=n*dt
    // calculate beta matrix to evaluate orthogonal dynamics
    std::vector<double> alpha(dimensionR*dimensionO,0);
    std::vector<double> norm(dimensionR,0);
    for (m=Neq;m<Nav;m++) {
      for (dR=0; dR<dimensionR; dR++){
	norm[dR] += data[columns*(m+n-1)+dR+beginP]*data[columns*(m+n-1)+dR+beginP];
      }
    }

    double Mb = 3.0/(norm[0]+norm[1]+norm[2])*(N-Ncor-Neq);

    for (m=Neq;m<Nav;m++) {
      for(dO=0; dO<dimensionO; dO++) {
	for (dR=0; dR<dimensionR; dR++){
	  alpha[dO*dimensionO+dR] += data[columns*(m+n-1)+dR+beginP]*Atilde[dimensionO*m+dO]*Mb;
	}
      }
    }

    for(dO=0; dO<dimensionO; dO++) {
      for (dR=0; dR<dimensionR; dR++){
	alpha[dO*dimensionO+dR] /= norm[dR]*Mb*Mb;
      }
    }

    result = 0;
    for (m=Neq;m<Nav;m++) {
      double partial_result = 0;
      for(dO=0; dO<dimensionO; dO++) {
	for (dR=0; dR<dimensionR; dR++){
	  Atilde[dimensionO*m+dO]=Atilde[dimensionO*m+dO]+alpha[dO*dimensionO+dR]*data[columns*(m+n-1)+dR+beginF]*dt;
	}
	partial_result += Atilde[dimensionO*m+dO]*data[columns*(m+n)+dO+beginB];
      }
      result += partial_result;
    }
    memory[n]=result/(Nav-Neq);
  }

  FILE * out;
  out = fopen(backward_filename, "w");

  for(int n=0; n<Ncor; n++)
    fprintf(out, "%f %f\n", n*dt, memory[n]);
  fclose(out); 

}

/**************************************************************
 * calculate memory function using backward orthogonal dynamics (second order)
 * ref: Rotenberg et. al, J. Chem.Phys. 140, 124103 (2014) / Jung et al., JCTC
 * int beginA: index of first column to include for A
 * int beginB: index of first column to include for B
 * int dimensionO: dimension of the vectors A and B
 * int beginP: index of first column to include for relevant variable
 * int beginF: index of first column to include for derivative of relevant variable
 * int dimensionR: dimension of the relevant variable
 * ***********************************************************/

void memory_backwardsecond(int beginA, int beginB, int dimensionO, int beginP, int beginF, int dimensionR ){
  
  int n,m,dO,dR;
  double result;
	
  std::vector<double> Atilde(dimensionO*Nav,0);
  std::vector<double> memory(Ncor);

  n=0;
  // special case n=0
  result=0;
  for (m=Neq;m<Nav;m++) {
    double partial_result = 0;
    for (dO=0; dO<dimensionO; dO++){
      //determine new Bplus
      Atilde[dimensionO*m+dO]=data[columns*m+dO+beginA];
      partial_result += data[columns*m+dO+beginB]*Atilde[dimensionO*m+dO];
    }
    //calculate correlation function
    result+=partial_result;
  }
  memory[0]=result/(Nav-Neq);

  for(n=1;n<Ncor;n++){ //main-loop for t=n*dt
    // calculate beta matrix to evaluate orthogonal dynamics
    std::vector<double> alpha(dimensionR*dimensionO,0);
    std::vector<double> zeta(dimensionR*dimensionO,0);
    std::vector<double> norm(dimensionR,0);
    double kappa = 0;
    double epsilon = 0;
    for (m=Neq+1;m<Nav;m++) {
      for (dR=0; dR<dimensionR; dR++){
	norm[dR] += data[columns*(m+n-2)+dR+beginP]*data[columns*(m+n-2)+dR+beginP];
	kappa += data[columns*(m+n-2)+dR+beginF]*data[columns*(m+n-2)+dR+beginP];
	epsilon += data[columns*(m+n-2)+dR+beginF]*data[columns*(m+n-1)+dR+beginP];
      }
    }

    double Mb2 = 3.0/(norm[0]+norm[1]+norm[2])*(N-Ncor-Neq);

    for (m=Neq+1;m<Nav;m++) {
      for(dO=0; dO<dimensionO; dO++) {
	for (dR=0; dR<dimensionR; dR++){
	  alpha[dO*dimensionO+dR] += data[columns*(m+n-2)+dR+beginP]*Atilde[dimensionO*(m-1)+dO]*Mb2;
	  zeta[dO*dimensionO+dR] += data[columns*(m+n-1)+dR+beginP]*Atilde[dimensionO*(m-1)+dO]*Mb2;
	}
      }
    }

    for(dO=0; dO<dimensionO; dO++) {
      for (dR=0; dR<dimensionR; dR++){
	alpha[dO*dimensionO+dR] /= norm[dR]*Mb2*Mb2;
	zeta[dO*dimensionO+dR] /= norm[dR]*Mb2*Mb2;
      }
    }
    kappa /= (norm[0]+norm[1]+norm[2])*Mb2*Mb2;
    epsilon /= (norm[0]+norm[1]+norm[2])*Mb2*Mb2;

    result = 0;
    for (m=Neq+1;m<Nav;m++) {
      double partial_result = 0;
      for(dO=0; dO<dimensionO; dO++) {std::vector<char*> cor_filename;
	for (dR=0; dR<dimensionR; dR++){
	  Atilde[dimensionO*(m-1)+dO]=Atilde[dimensionO*(m-1)+dO]+alpha[dO*dimensionO+dR]*data[columns*(m+n-2)+dR+beginF]*dt/2.0+dt/2.0/(1-dt/2.0*kappa*Mb2)*data[columns*(m+n-1)+dR+beginF]*(zeta[dO*dimensionO+dR]+dt/2.0*epsilon*Mb2*alpha[dO*dimensionO+dR]);
	}
	partial_result += Atilde[dimensionO*(m-1)+dO]*data[columns*(m+n-1)+dO+beginB];
      }
      result += partial_result;
    }
    memory[n]=result/(Nav-Neq);
  }
  
  
  FILE * out;
  out = fopen(backward_filename_2nd, "w");

  for(int n=0; n<Ncor; n++)
    fprintf(out, "%f %f\n", n*dt, memory[n]);
  fclose(out); 

}


/******** special functions **********/


/***********************************************************
 * calculate the stoachstic force using the memory kernel and the trajectory
 * ********************************************************/
void calc_stochastic_force() {
  int n,m;
  std::vector<double> data_stochastic_force(N);
  std::vector<double> cor(Ncor);
  std::vector<double> cor_error(Ncor);
  
  for (m=Neq; m<N; m++) {
    data_stochastic_force[m-Neq] = data[columns*m+6+force_index];
    data_stochastic_force[m-Neq] += 0.5*dt*data_memory_volterra[0]*data[columns*m+3+force_index];
    for (n=1; n<Ncor; n++) {
      //printf(" %f\n", n*dt);
      data_stochastic_force[m-Neq] += dt*data_memory_volterra[n]*data[columns*(m-n)+3+force_index];
    }
  }

  // print out data
  FILE * out;
  char array1[] ="data_";
  char * newArray = new char[std::strlen(array1)+std::strlen(sfFile)+1];
  std::strcpy(newArray,array1);
  std::strcat(newArray,sfFile);
  out = fopen(newArray, "w");

  for (m=0; m<N-Neq; m++) 
    fprintf(out, "%f %f\n", m*dt, data_stochastic_force[m]);
  fclose(out); 

//  const char* outFile = "check.dat";
//  out = fopen(outFile, "w");
//  for (m=Neq; m<N; m++) 
//    fprintf(out, "%f %f %f\n", m*dt, data[columns*m+3+force_index],data_stochastic_force[m-Neq]);
//  fclose(out); 

  
  // correlation
    n=0;
  // special case n=0 
  double partial_result = 0;
  double error = 0;
  for (m=0; m<N-Neq; m++) {
    partial_result+=data_stochastic_force[m]*data_stochastic_force[m];
    error += partial_result*partial_result;
  }
  cor[0] = partial_result/(N-Neq);
  //include correlation (about 2*dt steps) into error calculation
  cor_error[0] = sqrt(error/(N-Neq)-cor[0]*cor[0])/sqrt((N-Neq)/2*dt);

  for (n=1; n<Ncor; n++) {
    partial_result = 0;
    error = 0;
    for (m=0;m<N-Neq-n;m++) {
      partial_result+=data_stochastic_force[m]*data_stochastic_force[m+n];
      error += partial_result*partial_result;
    }
    cor[n]=partial_result/(N-Neq-n);
    cor_error[n] = sqrt(error/(N-Neq-n)-cor[n]*cor[n])/sqrt((N-Neq-n)/2*dt);
  }
	
  char array2[] ="correlation_";
  newArray = new char[std::strlen(array2)+std::strlen(sfFile)+1];
  std::strcpy(newArray,array2);
  std::strcat(newArray,sfFile);
  out = fopen(newArray, "w");
  for(n=0; n<Ncor; n++)
    fprintf(out, "%f %f %f\n",n*dt, cor[n], cor_error[n]);
  fclose(out); 
  
  
  // distribution
  double fmax = 1000;
  double vmax = 10.;
  double vavg = 0.;
  double vstd = 0.;
  double abs_v;
  std::vector<int> dist(nbin);
  std::vector<int> dist_v(nbin);
  std::vector<int> dist_absv(nbin);
  
  for (n=0; n<nbin; n++) {
    dist[n] = 0;
    dist_v[n] = 0;
  }
  int ntot = 0;
  int ntot_v = 0;
  int ntot_absv = 0;
  for (m=0; m<N-Neq; m++) {
    double data_f = data_stochastic_force[m] + fmax;
    double data_v = data[columns*(m+Neq)+3+force_index] + vmax;
    if (dim==3) {
      abs_v = sqrt(data[columns*(m+Neq)+3]*data[columns*(m+Neq)+3]+data[columns*(m+Neq)+4]*data[columns*(m+Neq)+4]+data[columns*(m+Neq)+5]*data[columns*(m+Neq)+5]);
    }
    if (dim==2) {
      abs_v = sqrt(data[columns*(m+Neq)+3]*data[columns*(m+Neq)+3]+data[columns*(m+Neq)+4]*data[columns*(m+Neq)+4]);
    }
    vavg += data[columns*(m+Neq)+3+force_index];
    vstd += data[columns*(m+Neq)+3+force_index]*data[columns*(m+Neq)+3+force_index];

    int bin = data_f*((double) nbin)/(2.0*fmax);
    int bin_v = data_v*((double) nbin)/(2.0*vmax);
    int bin_absv = abs_v*((double) nbin)/(vmax);
    if (bin >= 0 && bin <nbin) {
      ntot++;
      dist[bin]++;
    }
    if (bin_v >= 0 && bin_v <nbin) {
      ntot_v++;
      dist_v[bin_v]++;
    }
    if (bin_absv >= 0 && bin_absv <nbin) {
      ntot_absv++;
      dist_absv[bin_absv]++;
    }
  }
  
  char array3[] ="distribution_";
  newArray = new char[std::strlen(array3)+std::strlen(sfFile)+1];
  std::strcpy(newArray,array3);
  std::strcat(newArray,sfFile);
  out = fopen(newArray, "w");

  for(n=0; n<nbin; n++)
    fprintf(out, "%f %f\n",(((double) n)+0.5)*2.0*fmax/((double) nbin)-fmax,  ((double) 2.0*fmax*dist[n])/((double) ntot) );
  fclose(out); 

  char array4[] ="distribution_velocity_";
  newArray = new char[std::strlen(array4)+std::strlen(sfFile)+1];
  std::strcpy(newArray,array4);
  std::strcat(newArray,sfFile);
  out = fopen(newArray, "w");
  for(n=0; n<nbin; n++)
    fprintf(out, "%f %f\n",(((double) n)+0.5)*2.0*vmax/((double) nbin)-vmax, ((double) 2.0*vmax*dist_v[n])/((double) ntot_v) );
  fclose(out); 

  char array5[] ="avg_std_vel_";
  newArray = new char[std::strlen(array5)+std::strlen(sfFile)+1];
  std::strcpy(newArray,array5);
  std::strcat(newArray,sfFile);
  out = fopen(newArray, "w");
  fprintf(out, "%f %f\n",((double) vavg)/((double) (N-Neq)), sqrt( (vstd-vavg*vavg/((double) (N-Neq)))/((double) (N-Neq))) );
  fclose(out); 

  char array6[] ="distribution_absvel_";
  newArray = new char[std::strlen(array6)+std::strlen(sfFile)+1];
  std::strcpy(newArray,array6);
  std::strcat(newArray,sfFile);
  out = fopen(newArray, "w");
  for(n=0; n<nbin; n++)
    fprintf(out, "%f %f\n",(((double) n)+0.5)*vmax/((double) nbin), ((double) vmax*dist_absv[n])/((double) ntot_absv) );
  fclose(out); 

}

/***********************************************************
 * calculate the important quantities for the linear response project
 * ********************************************************/
void calc_linear_response() {
  
  int dimension=1;
  int n,m,d;
  double result = 0;
  double result3 = 0;
  double result4 = 0;
  double error = 0;

  std::vector<double> Azero(dimension*N);
  std::vector<double> cor(Ncor);
  std::vector<double> cor_error(Ncor);
  
  std::vector<double> cor2(Ncor);
  std::vector<double> cor_error2(Ncor);
  
  std::vector<double> cor3(Ncor);
  std::vector<int> countg3(Ncor);
  std::vector<double> cor_error3(Ncor);
  
  std::vector<double> cor4(Ncor);
  std::vector<int> countg4(Ncor);
  std::vector<double> cor_error4(Ncor);
	
  n=0;
  // special case n=0 
  int count3=0;
  int count4=0;
  for (m=Neq; m<N; m++) {
    for (d=0; d<dimension; d++){
      Azero[dimension*m+d]=data[columns*m+d+3];
      //printf("partial: %f\n",Azero[dimension*m+d]);
    }
    double partial_result = 0;
    double partial_result3 = 0;
    double partial_result4 = 0;
    for (d=0; d<dimension; d++){
      partial_result+=Azero[dimension*m+d]*sqrt(data[columns*m+d+4]*data[columns*m+d+4]+data[columns*m+d+5]*data[columns*m+d+5]);
      if (data[columns*m+d+4] > 0.0) {
        partial_result3+=Azero[dimension*m+d]*data[columns*m+d+4];
        count3++;
      } else {
        partial_result4+=Azero[dimension*m+d]*data[columns*m+d+4];
        count4++;
      }
    }
    result += partial_result;
    result3 += partial_result3;
    result4 += partial_result4;
    error += partial_result*partial_result;
    
  }
  cor[0] = result/(N-Neq);
  cor2[0] = result/(N-Neq);
  cor3[0] = result3/count3;
  cor4[0] = result4/count4;
  countg3[0] = count3;
 countg4[0] = count4;
  //include correlation (about 2*dt steps) into error calculation
  cor_error[0] = sqrt(error/(N-Neq)-cor[0]*cor[0])/sqrt((N-Neq)/2*dt);
  cor_error2[0] = sqrt(error/(N-Neq)-cor[0]*cor[0])/sqrt((N-Neq)/2*dt);

  for (n=1; n<Ncor; n++) {
    count3=0;
    count4=0;
    result = 0;
    error = 0;
    double result2 = 0;
    double error2 = 0;
    result3= result4=0.0;
    for (m=Neq;m<N-n;m++) {
      double partial_result = 0;
      double partial_result2 = 0;
      double partial_result3 = 0;
      double partial_result4 = 0;
      for (d=0; d<dimension; d++){
	partial_result+=Azero[dimension*m+d]*sqrt(data[columns*m+columns*n+d+4]*data[columns*m+columns*n+d+4]+data[columns*m+columns*n+d+5]*data[columns*m+columns*n+d+5]);
  partial_result2+=Azero[dimension*m+dimension*n+d]*sqrt(data[columns*m+d+4]*data[columns*m+d+4]+data[columns*m+d+5]*data[columns*m+d+5]);
        if (data[columns*m+columns*n+d+4] > 0.0) {
        partial_result3+=Azero[dimension*m+d]*data[columns*m+columns*n+d+4];
        count3++;
      } else {
        partial_result4+=Azero[dimension*m+d]*data[columns*m+columns*n+d+4];
        count4++;
      }
      }
      result += partial_result;
      error += partial_result*partial_result;
      result2 += partial_result2;
      error2 += partial_result2*partial_result2;
      result3 += partial_result3;
    result4 += partial_result4;
    }

    cor[n]=result/(N-Neq-n);
    cor2[n]=result2/(N-Neq-n);
    cor3[n]=result3/(count3);
    cor4[n]=result4/(count4);
      countg3[n] = count3;
 countg4[n] = count4;
    cor_error[n] = sqrt(error/(N-Neq)-cor[n]*cor[n])/sqrt((N-Neq)/2*dt);
    cor_error2[n] = sqrt(error2/(N-Neq)-cor2[n]*cor2[n])/sqrt((N-Neq)/2*dt);

  }
	
	const char* outFile = "vv_perp_test.dat";
  FILE * out;
  out = fopen(outFile, "w");
  for(n=0; n<Ncor; n++)
    fprintf(out, "%f %f %f %f %f %f %d %f %d\n",n*dt, cor[n]/dimension, cor_error[n]/dimension,cor2[n]/dimension, cor_error2[n]/dimension,cor3[n],countg3[n],cor4[n],countg4[n]);
  fclose(out); 
  
}


/***********************************************************
 * calculate mean squared displacement
 * ********************************************************/
void calc_msd() {
  
  FILE * out;
  out = fopen(msdFile, "w");

  std::vector<double> x(N);
  std::vector<double> y(N);
  std::vector<double> z(N);
  std::vector<double> msd(Ncor);

  data.resize(N*columns);
  
  for (int i=Neq; i<N; i++) {
    if (i%fCG==0) {
      x[i] = data[columns*i/fCG];
      y[i] = data[columns*i/fCG+1];
      z[i] = data[columns*i/fCG+2];
    }
    //printf("%d %f\n",i,y[i]);
  }

  for (int i=1; i<Ncor; i++) {
    double msd_x_sum=0;
    double msd_y_sum=0;
    double msd_z_sum=0;    
    
    for (int j=Neq; j<N-i; j++) {
//printf("%f %f\n",x[j+i],x[j]);
      msd_x_sum += (x[j+i]-x[j])*(x[j+i]-x[j]);
      msd_y_sum += (y[j+i]-y[j])*(y[j+i]-y[j]);
      msd_z_sum += (z[j+i]-z[j])*(z[j+i]-z[j]);
    }

    double msd_x = msd_x_sum/(N-Neq-i);
    double msd_y = msd_y_sum/(N-Neq-i);
    double msd_z = msd_z_sum/(N-Neq-i);
    
    fprintf(out, "%f %f\n", i*dt, msd_x+msd_y+msd_z);

  }

  fclose(out); 

}

/***********************************************************
 * calculate joint probability of the stoachstic force
 * ********************************************************/
void calc_joint_dist() {

  std::vector<double> jd1;	
  std::vector<double> jd2;
  std::vector<double> f1;	
  std::vector<double> f2;
  //data_memory_volterra.resize(Ncor);
  double dummy_f,dummy_p;

  //read in
  FILE * in;

  in = fopen(jd1File,"r");

  for(int i=0;i<nbin;i++){
    fscanf(in,"%lf %lf\n",&dummy_f,&dummy_p);
    jd1.push_back(dummy_p);
    f1.push_back(dummy_f);
    
  }
  fclose(in);

  in = fopen(jd2File,"r");
  for(int i=0;i<nbin;i++){
    fscanf(in,"%lf %lf\n",&dummy_f,&dummy_p);
    jd2.push_back(dummy_p);
    f2.push_back(dummy_f);
    //printf("%lf %lf\n",jd2[i],f2[i]);
  }
  fclose(in);

  FILE * out;
  out = fopen(joint_dist, "w");
  for(int n=0; n<nbin; n++) {
    for(int m=0; m<nbin; m++) {
      if(jd1[n]*jd2[m] != 0) {
        fprintf(out, "%f %f %f\n",f1[n],f2[m],(double (jd1[n]*jd2[m])));
      }
    }
  }

  fclose(out); 

/*
double tot=0;
    for(int n=0; n<nbin; n++) {
      for(int m=0; m<nbin; m++) {
        if(jd1[n]*jd2[m] != 0) {
          tot += jd1[n]*jd2[m];
        }
      }
    }
*/

}

/***********************************************************
 * calculate proj of noise on velocity
 * ********************************************************/
void calc_proj_noise() {

  std::vector<double> sfx;	
  std::vector<double> sfy;
  std::vector<double> sfz;
  std::vector<double> proj(N-Neq);
  std::vector<double> dot(N-Neq);
  double sfxavg = 0.;
  double sfxstd = 0.;		
  double sfyavg = 0.;
  double sfystd = 0.;		
  double sfzavg = 0.;
  double sfzstd = 0.;		
  double dummy_t, dummy_f, mag_v;
  //read in
  FILE * in;
  FILE * out;

  in = fopen(sfxFile,"r");

  for(int i=0;i<N-Neq;i++){
    fscanf(in,"%lf %lf\n",&dummy_t,&dummy_f);
    sfx.push_back(dummy_f);  
  }
  fclose(in);

  in = fopen(sfyFile,"r");

  for(int i=0;i<N-Neq;i++){
    fscanf(in,"%lf %lf\n",&dummy_t,&dummy_f);
    sfy.push_back(dummy_f);  
  }
  fclose(in);

  if (dim==3){

    in = fopen(sfzFile,"r");

    for(int i=0;i<N-Neq;i++){
      fscanf(in,"%lf %lf\n",&dummy_t,&dummy_f);
      sfz.push_back(dummy_f);  
    }
    fclose(in);
  }

  out = fopen(data_proj_file, "w");
  for(int n=0; n<N-Neq; n++) {
    sfxavg += sfx[n];
    sfxstd += sfx[n]*sfx[n];
    sfyavg += sfy[n];
    sfystd += sfy[n]*sfy[n];
    if (dim==3){
      sfzavg += sfz[n];
      sfzstd += sfz[n]*sfz[n];
    }

    //printf("Velocity: %f %f %f\n",data[columns*(n+Neq)/fCG+3],data[columns*(n+Neq)/fCG+4],data[columns*(n+Neq)/fCG+5]);
    //printf("SF: %f %f %f\n",sfx[n],sfy[n],sfz[n]);

    mag_v=data[columns*(n+Neq)/fCG+3]*data[columns*(n+Neq)/fCG+3]+data[columns*(n+Neq)/fCG+4]*data[columns*(n+Neq)/fCG+4];
    if (dim==3) mag_v+=data[columns*(n+Neq)/fCG+5]*data[columns*(n+Neq)/fCG+5];

    dot[n]=data[columns*(n+Neq)/fCG+3]*sfx[n]+data[columns*(n+Neq)/fCG+4]*sfy[n];
    if (dim==3) dot[n]+=data[columns*(n+Neq)/fCG+5]*sfz[n];

    //printf("Mag: %f, Dot: %f \n",mag_v,dot);

    proj[n] = dot[n]/sqrt(mag_v);

    //printf("Proj: %f \n",proj[n]);
    
    fprintf(out, "%f %f\n",dot[n],proj[n]);
  }

  fclose(out); 

  // distribution
  double proj_max = 1000;
  std::vector<int> proj_dist(nbin);
  double dot_max = 50;
  std::vector<int> dot_dist(nbin);
  
  for (int n=0; n<nbin; n++) {
    proj_dist[n] = 0;
    dot_dist[n] = 0;
  }
  int proj_ntot = 0;
  int dot_ntot = 0;
  for (int m=0; m<N-Neq; m++) {
    double data_proj = proj[m] + proj_max;
    double data_dot = dot[m] + dot_max;

    int bin_proj = data_proj*((double) nbin)/(2.0*proj_max);
    int bin_dot = data_dot*((double) nbin)/(2.0*dot_max);
    if (bin_proj >= 0 && bin_proj <nbin) {
      proj_ntot++;
      proj_dist[bin_proj]++;
    }
    if (bin_dot >= 0 && bin_dot <nbin) {
      dot_ntot++;
      dot_dist[bin_dot]++;
    }
  }

  out = fopen(dist_proj_file, "w");
  for(int n=0; n<nbin; n++) {
    fprintf(out, "%f %f\n",(((double) n)+0.5)*2.0*proj_max/((double) nbin)-proj_max, ((double) proj_dist[n])/((double) proj_ntot) );
  }

  fclose(out); 

  out = fopen(dist_dot_file, "w");
  for(int n=0; n<nbin; n++) {
    fprintf(out, "%f %f\n",(((double) n)+0.5)*2.0*dot_max/((double) nbin)-dot_max, ((double) dot_dist[n])/((double) dot_ntot) );
  }

  fclose(out);

  char array5[] ="avg_std_sf_";
  char * newArray = new char[std::strlen(array5)+std::strlen(sfFile)+1];
  std::strcpy(newArray,array5);
  std::strcat(newArray,sfFile);
  out = fopen(newArray, "w");
  fprintf(out, "%f %f %f %f %f %f\n",((double) sfxavg)/((double) (N-Neq)), sqrt( (sfxstd-sfxavg*sfxavg/((double) (N-Neq)))/((double) (N-Neq))),((double) sfyavg)/((double) (N-Neq)), sqrt( (sfystd-sfyavg*sfyavg/((double) (N-Neq)))/((double) (N-Neq))),((double) sfzavg)/((double) (N-Neq)), sqrt( (sfzstd-sfzavg*sfzavg/((double) (N-Neq)))/((double) (N-Neq))) );
  fclose(out);  

}
