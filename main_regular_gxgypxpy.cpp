//#include <omp.h>
#include "fci.hpp"
#include <cstring>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <sys/types.h>
#include <time.h>
#include <cmath>
extern "C"{
#include <stdio.h>
#include <stdlib.h>
}

using namespace std;
int myrank=0;
int psize=1;
int nrep_ed_max=2000;

extern "C"{
  double ran_();
  void initran_(int*);  
}
int model=1;
int main(int argc,char** argv){
  int lx,ly,i,j,kx,ky,thread,sec,ns,sec_gs,kx_gs,ky_gs,ntau,mlanc,px,py,qr,ppxy,zz,perix,periy,sgm1,sgm2,sec_target,px0,py0;
  unsigned long nbasis,nrep;
  double tau,tmp,jcoup,gcoup,hfldz=0,enr_gs,e0;
  time_t t1;
  bool ready;
  
  //t1=omp_get_wtime();
  t1=139;
  srand48((long)t1);
  initran_((int*)&t1);

  for(i=1;i<argc;i++){
    if(strcmp(argv[i],"-lx")==0){
      lx=atoi(argv[i+1]);
      if(myrank==0)cout<<"lx="<<lx<<endl;
    }
    else if(strcmp(argv[i],"-ly")==0){
      ly=atoi(argv[i+1]);
      if(myrank==0)cout<<"ly="<<ly<<endl;
    }
    else if(strcmp(argv[i],"-jcoup")==0){
      jcoup=atof(argv[i+1]);
      if(myrank==0)cout<<"jcoup="<<jcoup<<endl;
    }
    else if(strcmp(argv[i],"-hfldz")==0){
      hfldz=atof(argv[i+1]);
      cout<<"hfldz="<<hfldz<<endl;
    }
    else if(strcmp(argv[i],"-psize")==0){
      thread=atoi(argv[i+1]);
      cout<<"thread="<<thread<<endl;
      psize=thread;
    }
    else if(strcmp(argv[i],"-sec")==0){
      sec_target=atoi(argv[i+1]);
      cout<<"sec_target="<<sec_target<<endl;
    }
    else if(strcmp(argv[i],"-kx")==0){
      kx=atoi(argv[i+1]);
      cout<<"kx="<<kx<<endl;
    }
    else if(strcmp(argv[i],"-ky")==0){
      ky=atoi(argv[i+1]);
      cout<<"ky="<<ky<<endl;
    }
    else if(strcmp(argv[i],"-px")==0){
      px=atoi(argv[i+1]);
      cout<<"px="<<px<<endl;
    }
    else if(strcmp(argv[i],"-py")==0){
      py=atoi(argv[i+1]);
      cout<<"py="<<py<<endl;
    }
    else if(strcmp(argv[i],"-sgm1")==0){
      sgm1=atoi(argv[i+1]);
      cout<<"sgm1="<<sgm1<<endl;
    }
    else if(strcmp(argv[i],"-sgm2")==0){
      sgm2=atoi(argv[i+1]);
      cout<<"sgm2="<<sgm2<<endl;
    }
    else if(strcmp(argv[i],"-zz")==0){
      zz=atoi(argv[i+1]);
      cout<<"zz="<<zz<<endl;
    }
    else if(strcmp(argv[i],"-model")==0){
      model=atoi(argv[i+1]);
      cout<<"model="<<model<<endl;
    }
  }
  ns=lx*ly;
  perix=lx;
  periy=ly;
  nbasis=0;
  fci momentum;
  enr_gs=100;
  //omp_set_num_threads(thread);
  sec=sec_target;
  for(sec=sec_target;sec<=sec_target;sec++){
    if(fabs(hfldz)<1.e-6&&(sec%2)==1)continue;
    if((ky!=0&&ky!=(periy/2))&&(kx!=0&&kx!=(perix/2))&&(!(px==0&&py==0)))continue;
    if((kx==0||kx==(perix/2))&&(ky!=0&&ky!=(periy/2))&&(!(py==0&&abs(px)==1)))continue;
    if((ky==0||ky==(periy/2))&&(kx!=0&&kx!=(perix/2))&&(!(px==0&&abs(py)==1)))continue;
    if((ky==0||ky==periy/2)&&(kx==0||kx==perix/2)&&kx!=ky&&(!(abs(py)==1&&abs(px)==1)))continue;
    if((ky==0||ky==periy/2)&&(kx==0||kx==perix/2)&&kx==ky&&(!(abs(py)==1&&abs(px)==1)))continue;

    if((ky!=0&&ky!=periy/2)&&(kx!=0&&kx!=perix/2)&&(kx+ky)%perix==0&&(!(sgm1==0&&abs(sgm2)==1)))continue;
    if((ky!=0&&ky!=periy/2)&&(kx!=0&&kx!=perix/2)&&kx==ky&&(!(abs(sgm1)==1&&sgm2==0)))continue;	
    if(((kx+ky)%perix!=0)&&kx!=ky&&(!(sgm1==0&&sgm2==0)))continue;
    if(kx==ky&&(kx==0||kx==perix/2)&&px!=py&&(!(sgm1==0&&sgm2==0)))continue;
    if(kx==ky&&(kx==0||kx==perix/2)&&px==py&&(!(abs(sgm1)==1&&abs(sgm2)==1)))continue;
    
    if(sec==0&&zz==0)continue;
    else if(sec!=0&&zz!=0)continue;
    qr=0;
    cout<<"kx="<<kx<<"\tky="<<ky<<"\tsec="<<sec<<"\tpx="<<px<<"\tpy="<<py<<"\tz="<<zz<<"\tsgm1="<<sgm1<<"\tsgm2="<<sgm2<<"\tqr="<<qr<<endl;
    if(sec%2==1)
      ready=momentum.setup(lx,ly,jcoup,hfldz,kx,ky,ns/2-(sec+1)/2,px,py,zz,sgm1,sgm2,qr);
    else if(sec%2==0)
      ready=momentum.setup(lx,ly,jcoup,hfldz,kx,ky,ns/2+sec/2,px,py,zz,sgm1,sgm2,qr);
    nrep=momentum.get_nrep();
	    /*
	      if(fabs(hfldz)<1.e-6&&sec>0)
	      nbasis+=2*momentum.get_nrep();
	      else
	      nbasis+=momentum.get_nrep();
	      momentum.clean();
	      continue;
	    */
    if(ready&&nrep<nrep_ed_max){
      if((kx%(perix/2)==0)&&(ky%(periy/2)==0)&&qr%2==0){ //use real hamiltonian matrix
	momentum.hamiltonianr();
	momentum.checkhermitianr();
      }
      else{
	momentum.hamiltonian();
	cout<<"done buid complex hamiltonian "<<endl;
	momentum.checkhermitian();
      }
      momentum.diagonalize(1,0);
      momentum.measure();
      momentum.save_eigval();
    }
    else if(ready&&nrep>nrep_ed_max){
      cout<<"initialize lanczos"<<endl;
      momentum.initialize_lanczos();
      if((kx%(perix/2)==0)&&(ky%(periy/2)==0)&&qr%2==0){ //use real hamiltonian matrix
	mlanc=momentum.lanczos1real();
	//momentum.measure();
	momentum.save_eigval();
      }
      else{
	mlanc=momentum.lanczos1();
	//momentum.measure();
	momentum.save_eigval();
      }
    }
    if(fabs(hfldz)<1.e-6&&sec>0)
      nbasis+=2*momentum.get_nrep();
    else
      nbasis+=momentum.get_nrep();
    momentum.clean();
    cout<<"total number of basis="<<nbasis<<endl;
  }
}
