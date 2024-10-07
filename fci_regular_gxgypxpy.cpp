#include "fci.hpp"
//#include <omp.h>
#include <complex>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstring>
#include <stdlib.h>
#include <math.h>
#include <time.h>

using namespace std;
extern "C"{
  void dsyevx_(char*,char*,char*,int*,double*,int*,double*,double*,int*,int*,double*,int*,double*,double*,int*,double*,int*,int*,int*,int*);
  void dsyev_(char*,char*,int*,double*,int*,double*,double*,int*,int*);
  void zheev_(char*,char*,int*,dcmplex*,int*,double*,dcmplex*,int*,double*,int*);
  void zheevx_(char*,char*,char*,int*,dcmplex*,int*,double*,double*,int*,int*,double*,int*,double*,dcmplex*,int*,dcmplex*,int*,double*,int*,int*,int*);
  void dstev_(char*,int*,double*,double*,double*,int*,double*,int*);
  double ran_();
}
extern int nrep_ed_max,psize,model;
void idsort(unsigned long,unsigned long*,double*);
dcmplex inner_prod(unsigned long, dcmplex*,dcmplex*);
double  inner_prod(unsigned long, double*,double*);
void obtain_symmetric_matrix_eigenvector1(double*,double*,int,int);
void obtain_hermitian_matrix_eigenvector(dcmplex*,double*,int,int);
void obtain_hermitian_matrix_eigenvector(dcmplex*,double*,int);
void obtain_hermitian_matrix_eigenvector1(dcmplex*,double*,int,int);

//--------------------------------------------------------------------------------------
fci::fci(){
//--------------------------------------------------------------------------------------
  nrep=0;
  repr=NULL;
  normfac=NULL;
  eigval=NULL;
  totspn=NULL;
  f4func=NULL;
  mat=NULL;
  bst=NULL;
  plq=NULL;
  txmap=NULL;
  tymap=NULL;
  pxmap=NULL;
  pymap=NULL;
  sigma1map=NULL;
  sigma2map=NULL;
  rotmap=NULL;
  mapptr=NULL;
  transtep=NULL;
  ntransite=NULL;
  transite=NULL;
  tranbit=NULL;
  //lanczos space
  eigvec=NULL;
  lanczos_vec=NULL;
  aal=NULL;
  nnl=NULL;
  ff=NULL;
  lookuptable=NULL;
  matr=NULL;
  eigvecr=NULL;
  ffr=NULL;
  latticemaps=NULL;
}

//--------------------------------------------------------------------------------------
fci::fci(int x,int y,double j,double z,int kkx, int kky,int nup,int ppx,int ppy,int zz,int sigma1,int sigma2,int rr):lx(x),ly(y),jcoup(j),hfldz(z),kx(kkx),ky(kky),px(ppx),py(ppy),qz(zz),sgm1(sigma1),sgm2(sigma2),qr(rr){
//--------------------------------------------------------------------------------------
  setup(x,y,j,z,kkx,kky,nup,ppx,ppy,zz,sigma1,sigma2,rr);
}

//--------------------------------------------------------------------------------------
bool fci::setup(int x,int y,double j,double z,int kkx, int kky, int nup, int ppx, int ppy,int zz, int sigma1,int sigma2,int rr){
//--------------------------------------------------------------------------------------
  unsigned long nmat,s0;
  int i,k;
  bool check;
  double tm1,tm2;
  eps=1.e-12;
  lx=x;
  ly=y;
  if(ly!=lx){
    cout<<"couldn't process such a cluster lx="<<lx<<"\tly="<<ly<<endl;
    exit(0);
  }
  perix=lx;
  periy=ly;
  ns=lx*ly;
  jcoup=j;
  //if(lx==4&&ly==4&&model==2)jcoup*=2.;
  hfldz=z;
  kx=kkx;
  ky=kky;
  px=ppx;
  py=ppy;
  sgm1=sigma1;
  sgm2=sigma2;
  qz=zz;
  qr=rr;
  sec=(nup-(ns-nup))/2;
  if(model==1)
    build_hamiltonian_j1j2();
  else if(model==2)
    build_hamiltonian_j1j3();
  check=makefgfunctions();
  if(check==false)return check;
  build_latticemaps();
  //check_latticemaps();
  //tm1=omp_get_wtime();
  build_basis(nup);
  //tm2=omp_get_wtime();
  cout<<"build basis time="<<tm2-tm1<<endl;
  //return true;
  if(nrep>0){
    repr=new unsigned long[nrep];
    normfac=new double[nrep];
    //tm1=omp_get_wtime();
    build_basis1(nup);
    //tm2=omp_get_wtime();
    //for(i=0;i<nrep;i++)
    //cout<<"i="<<i<<"\t"<<repr[i]<<"\t"<<normfac[i]<<endl;
    cout<<"build basis1 time="<<tm2-tm1<<endl;
    if(nrep<nrep_ed_max){
      if(kx%(perix/2)==0&&ky%(periy/2)==0&&qr%2==0) //use real hamiltonian matrix
	matr=new double[nrep*nrep];	
      else   //use complex hamiltonian matrix
	mat=new dcmplex[nrep*nrep];
      eigval=new double[nrep];
      totspn=new double[nrep];
    }
    else{
      mlanc=200;
      neig=5;
      aal=new double[mlanc];
      nnl=new double[mlanc];
      eigval=new double[mlanc];
      totspn=new double[mlanc];
      lanczos_vec=new double[mlanc*mlanc];
      if(kx%(perix/2)==0&&ky%(periy/2)==0&&qr%2==0){ //use real hamiltonian matrix
	eigvecr=new double*[neig];
	eigvecr[0]=new double[nrep*(unsigned long)neig];
	for(i=1;i<neig;i++)
	  eigvecr[i]=&(eigvecr[0][(unsigned long)i*nrep]);
	ffr=new double*[mlanc];
	nmat=(unsigned long)mlanc*nrep;
	ffr[0]=new double[nmat];
	for(i=1;i<mlanc;i++){
	  nmat=(unsigned long)i*nrep;
	  ffr[i]=&(ffr[0][nmat]);
	}
      }
      else{
	eigvec=new dcmplex*[neig];
	eigvec[0]=new dcmplex[nrep*(unsigned long)neig];
	for(i=1;i<neig;i++)
	  eigvec[i]=&(eigvec[0][(unsigned long)i*nrep]);
	ff=new dcmplex*[mlanc];
	nmat=(unsigned long)mlanc*nrep;
	ff[0]=new dcmplex[nmat];
	for(i=1;i<mlanc;i++){
	  nmat=(unsigned long)i*nrep;
	  ff[i]=&(ff[0][nmat]);
	}
      }
    }
    return true;
  }
  else return false;
}

//--------------------------------------------------------------------------------------
fci::~fci(){
//--------------------------------------------------------------------------------------
  clean();
}

//--------------------------------------------------------------------------------------
void fci::clean(){
//--------------------------------------------------------------------------------------
  int i,j,k,l,m,tx,ty,tpx,tpy,tsx,tsy;
  if(repr!=NULL)
    delete []repr;
  if(normfac!=NULL)
    delete []normfac;
  if(eigval!=NULL)
    delete []eigval;
  if(totspn!=NULL)
    delete []totspn;
  if(f4func!=NULL){
    for(j=0;j<ly*2+1;j++){
      for(i=0;i<lx*2+1;i++)
	delete []f4func[j][i];
      delete []f4func[j];      
    }
    delete []f4func;
  }
  if(bst!=NULL){
    delete []bst[0];
    delete []bst;
  }
  if(mat!=NULL)
    delete []mat;
  if(plq!=NULL){
    delete []plq[0];
    delete []plq;
  }
  if(eigvec!=NULL){
    delete []eigvec[0];
    delete []eigvec;
  }
  if(lanczos_vec!=NULL)
    delete []lanczos_vec;
  if(aal!=NULL)
    delete []aal;
  if(nnl!=NULL)
    delete []nnl;
  if(ff!=NULL){
    delete []ff[0];
    delete []ff;
  }
  if(txmap!=NULL)
    delete []txmap;
  if(tymap!=NULL)
    delete []tymap;
  if(sigma1map!=NULL)
    delete []sigma1map;
  if(sigma2map!=NULL)
    delete []sigma2map;
  if(rotmap!=NULL)
    delete []rotmap;
  if(lookuptable!=NULL)
    delete []lookuptable;
  if(matr!=NULL)
    delete []matr;
  if(ffr!=NULL){
    delete []ffr[0];
    delete []ffr;
  }
  if(eigvecr!=NULL){
    delete []eigvecr[0];
    delete []eigvecr;
  }
  if(latticemaps!=NULL){
    for(tx=0;tx<perix;tx+=1){
      for(ty=0;ty<periy;ty+=1){
	for(tpx=0;tpx<2;tpx++){
	  for(tpy=0;tpy<2;tpy++){
	    for(tsx=0;tsx<2;tsx++){
	      for(tsy=0;tsy<2;tsy++)
		delete []latticemaps[tx][ty][tpx][tpy][tsx][tsy];
	      delete []latticemaps[tx][ty][tpx][tpy][tsx];
	    }
	    delete []latticemaps[tx][ty][tpx][tpy];
	  }
	  delete []latticemaps[tx][ty][tpx];
	}
	delete []latticemaps[tx][ty];
      }
      delete []latticemaps[tx];
    }
    delete []latticemaps;
  }
  if(mapptr!=NULL)delete []mapptr;
  if(transtep!=NULL){
    delete []transtep[0];
    delete []transtep;
  }
  if(tranbit!=NULL)delete []tranbit;
  if(ntransite!=NULL)delete []ntransite;
  if(transite!=NULL){
    for(i=0;i<ntrans;i++)
      delete []transite[i];
    delete []transite;
  }
  ntransite=NULL;
  transite=NULL;
  tranbit=NULL;
  mapptr=NULL;
  transtep=NULL;
  latticemaps=NULL;
  nrep=0;
  repr=NULL;
  normfac=NULL;
  eigval=NULL;
  totspn=NULL;
  f4func=NULL;
  mat=NULL;
  bst=NULL;
  plq=NULL;
  txmap=NULL;
  tymap=NULL;
  pxmap=NULL;
  pymap=NULL;
  sigma1map=NULL;
  sigma2map=NULL;
  rotmap=NULL;
  //lanczos space
  eigvec=NULL;
  lanczos_vec=NULL;
  aal=NULL;
  nnl=NULL;
  ff=NULL;
  lookuptable=NULL;
  matr=NULL;
  eigvecr=NULL;
  ffr=NULL;
}

//--------------------------------------------------------------------------------------
void fci::build_hamiltonian(){
//--------------------------------------------------------------------------------------
  int i,j,i1,i2,j1,k,l,s[4];
  ns=lx*ly;
  nb=4*ns;
  bst=new int*[ns*4];
  bst[0]=new int[ns*4*2];
  for(i=1;i<ns*4;i++)
    bst[i]=&(bst[0][i*2]);
  for(i=0;i<lx;i++)
    for(j=0;j<ly;j++){
      i1=(i+1)%lx;
      j1=(j+1)%ly;
      i2=(i+lx-1)%lx;
      l=i+j*lx;
      s[0]=i1+j*lx;
      s[1]=i+j1*lx;
      s[2]=i1+j1*lx;
      s[3]=i2+j1*lx;
      for(k=0;k<4;k++)
	if(s[k]>l){
	  bst[l+k*ns][0]=l;
	  bst[l+k*ns][1]=s[k];
	}
	else{
	  bst[l+k*ns][0]=s[k];
	  bst[l+k*ns][1]=l;
	}
    }
}
 
//--------------------------------------------------------------------------------------
void fci::build_hamiltonian_shastry_sutherland_32(){
//--------------------------------------------------------------------------------------
  int i,j,k,l,m,i0,i1,i2,i3,i4,cnt;
  ofstream fout;
  ns=lx*ly;
  nb=2*ns+ns/2;
  bst=new int*[nb];
  bst[0]=new int[nb*2];
  for(i=1;i<nb;i++)
    bst[i]=&(bst[0][i*2]);
  plq=new int*[ns];
  plq[0]=new int[ns*4];
  for(i=1;i<ns;i++)
    plq[i]=&(plq[0][i*4]);
  cnt=0;
  for(j=0;j<ly;j++)
    for(i=0;i<lx;i++){
      k=i+j*lx;
      if(j%2==0){
	plq[k][0]=i+j*lx;
	plq[k][1]=i+((j+1)%ly)*lx;
	plq[k][2]=(i+1)%lx+((j+1)%ly)*lx;
	plq[k][3]=i+((j+2)%ly)*lx;
	l=i+((j+1)%ly)*lx;
	m=(i+1)%lx+((j+1)%ly)*lx;
      }
      else if(j%2==1){
	plq[k][0]=i+j*lx;
	plq[k][1]=(i+lx-1)%lx+((j+1)%ly)*lx;
	plq[k][2]=i+((j+1)%ly)*lx;
	plq[k][3]=i+((j+2)%ly)*lx;
	l=(i+lx-1)%lx+((j+1)%ly)*lx;
	m=i+((j+1)%ly)*lx;
      }
      bst[k][0]=k;
      bst[k][1]=l;
      bst[k+ns][0]=k;
      bst[k+ns][1]=m;
      if((j/2)%2==0&&(j%2)==0&&(i%2)==1){
	l=i+((j+2)%ly)*lx;
	bst[cnt+ns*2][0]=k;
	bst[cnt+ns*2][1]=l;
	cnt++;
      }
      else if((j/2)%2==0&&(j%2)==1&&(i%2)==0){
	l=(i+1)%lx+j*lx;
	bst[cnt+ns*2][0]=k;
	bst[cnt+ns*2][1]=l;
	cnt++;
      }
      else if((j/2)%2==1&&(j%2)==0&&(i%2)==0){
	l=i+((j+2)%ly)*lx;
	bst[cnt+ns*2][0]=k;
	bst[cnt+ns*2][1]=l;
	cnt++;
      }
      else if((j/2)%2==1&&(j%2)==1&&(i%2)==1){
	l=(i+1)%lx+j*lx;
	bst[cnt+ns*2][0]=k;
	bst[cnt+ns*2][1]=l;
	cnt++;
      }
    }	
  txmap=new int[ns*2*lx];
  tymap=new int[ns*2*ly];
  pxmap=&(tymap[ns*ly]);
  pymap=&(txmap[ns*lx]);
  sigma1map=new int[ns*2];
  sigma2map=new int[ns*2];
  rotmap=new int[ns*4];
  for(i=0;i<lx;i++)
    for(j=0;j<ly;j++){
      k=i+j*lx;
      if(j%2==0)
	i1=((lx-1-i)%lx)+((j+2)%ly)*lx;
      else if(j%2==1)
	i1=((lx-i)%lx)+((j+2)%ly)*lx;
      i2=((i+1)%lx)+((ly-j)%ly)*lx;
      i3=((i+1)%lx)+j*lx;
      i4=i+((j+2)%ly)*lx;
      pxmap[ns+k]=i1;
      pymap[ns+k]=i2;
      txmap[ns+k]=i3;
      tymap[ns+k]=i4;
      pxmap[k]=k;
      pymap[k]=k;
      txmap[k]=k;
      tymap[k]=k;
      //reflection map
      i1=i+((ly-(j-1)+1)%ly)*lx;
      sigma2map[k]=k;
      sigma2map[k+ns]=i1;
      if(j%2==1)
	i1=((lx-(i-1))%lx)+j*lx;
      else if(j%2==0)
	i1=((lx-i)%lx)+j*lx;
      sigma1map[k]=k;
      sigma1map[k+ns]=i1;
    }
  for(i=2;i<lx;i++)
    for(j=0;j<ns;j++){
      txmap[i*ns+j]=txmap[ns+txmap[(i-1)*ns+j]];
      pymap[i*ns+j]=pymap[ns+pymap[(i-1)*ns+j]];
    }
  for(i=2;i<ly;i++)
    for(j=0;j<ns;j++){
      tymap[i*ns+j]=tymap[ns+tymap[(i-1)*ns+j]];
      pxmap[i*ns+j]=pxmap[ns+pxmap[(i-1)*ns+j]];
    }
}

//--------------------------------------------------------------------------------------
void fci::build_hamiltonian_shastry_sutherland(){
//--------------------------------------------------------------------------------------
  int i,j,i0,i1,i2,i3,i4;
  ofstream fout;
  ns=lx*ly;
  nb=2*ns+ns/2;
  bst=new int*[nb];
  bst[0]=new int[nb*2];
  for(i=1;i<nb;i++)
    bst[i]=&(bst[0][i*2]);
  for(i=0;i<lx;i++){
    for(j=0;j<ly;j++){
      i0=i+j*lx;
      i1=(i+1)%lx+j*lx;
      i2=i+((j+1)%ly)*lx;
      bst[i0][0]=i0;
      bst[i0][1]=i1;
      bst[i0+ns][0]=i0;
      bst[i0+ns][1]=i2;
      if(j%2==0&&i%2==0){
	i1=(i+1)%lx+((j+1)%ly)*lx;
	i3=i+(j/2)*lx;
	bst[i3+ns*2][0]=i0;
	bst[i3+ns*2][1]=i1;
      }
      else if(j%2==0&&i%2==1){
	i1=(i+1)%lx+((j+ly-1)%ly)*lx;
	i3=i+(j/2)*lx;
	bst[i3+ns*2][0]=i0;
	bst[i3+ns*2][1]=i1;
      }
    }
  }
  plq=new int*[ns];
  plq[0]=new int[ns*4];
  for(i=1;i<ns;i++)
    plq[i]=&(plq[0][i*4]);
  for(j=0;j<ly;j++)
    for(i=0;i<lx;i++){
      i0=i+j*lx;
      i1=(i+1)%lx+j*lx;
      i2=(i+1)%lx+((j+1)%ly)*lx;
      i3=i+((j+1)%ly)*lx;
      plq[i0][0]=i0;
      plq[i0][1]=i1;
      plq[i0][2]=i2;
      plq[i0][3]=i3;
    }
  txmap=new int[ns*2*lx];
  tymap=new int[ns*2*ly];
  pxmap=&(tymap[ns*ly]);//gliding reflection x map
  pymap=&(txmap[ns*lx]);//gliding reflection y map
  sigma1map=new int[ns*2];
  sigma2map=new int[ns*2];
  rotmap=new int[ns*4];
  for(i=0;i<lx;i++)
    for(j=0;j<ly;j++){
      i0=i+j*lx;
      i1=((lx-i)%lx)+((j+1)%ly)*lx;
      i2=((i+1)%lx)+((ly-j)%ly)*lx;
      i3=((i+1)%lx)+j*lx;
      i4=i+((j+1)%ly)*lx;
      pxmap[ns+i0]=i1;
      pymap[ns+i0]=i2;
      txmap[ns+i0]=i3;
      tymap[ns+i0]=i4;
      sigma1map[ns+i0]=j+i*lx;
      sigma2map[ns+i0]=((lx-1-j)%lx)+((ly-1-i)%ly)*lx;      
      pxmap[i0]=i0;
      pymap[i0]=i0;
      txmap[i0]=i0;
      tymap[i0]=i0;
      sigma1map[i0]=i0;
      sigma2map[i0]=i0;
    }
  for(i=2;i<lx;i++)
    for(j=0;j<ns;j++){
      txmap[i*ns+j]=txmap[ns+txmap[(i-1)*ns+j]];
      pymap[i*ns+j]=pymap[ns+pymap[(i-1)*ns+j]];
    }
  for(i=2;i<ly;i++)
    for(j=0;j<ns;j++){
      tymap[i*ns+j]=tymap[ns+tymap[(i-1)*ns+j]];
      pxmap[i*ns+j]=pxmap[ns+pxmap[(i-1)*ns+j]];
    }
  for(i=0;i<lx;i++)
    for(j=0;j<ly;j++){
      //i0=(lx-1-j)+i*lx;
      i0=(lx-j)%lx+((i+1)%ly)*lx;
      i1=i+j*lx;
      rotmap[i1]=i1;
      rotmap[ns+i1]=i0;
      //cout<<"i1="<<i1<<"\trotmap="<<i0<<endl;
    }
  for(i=2;i<4;i++)
    for(j=0;j<ns;j++)
      rotmap[i*ns+j]=rotmap[ns+rotmap[(i-1)*ns+j]];
}

//--------------------------------------------------------------------------------------
void fci::build_hamiltonian_j1j2(){
//--------------------------------------------------------------------------------------
  int i,j,i0,i1,i2,i3,i4;
  ofstream fout;
  ns=lx*ly;
  nb=4*ns;
  bst=new int*[nb];
  bst[0]=new int[nb*2];
  for(i=1;i<nb;i++)
    bst[i]=&(bst[0][i*2]);
  for(i=0;i<lx;i++){
    for(j=0;j<ly;j++){
      i0=i+j*lx;
      i1=(i+1)%lx+j*lx;
      i2=i+((j+1)%ly)*lx;
      bst[i0][0]=i0;
      bst[i0][1]=i1;
      bst[i0+ns][0]=i0;
      bst[i0+ns][1]=i2;
      i1=(i+1)%lx+((j+1)%ly)*lx;
      i2=(i+1)%lx+((j+ly-1)%ly)*lx;
      bst[i0+ns*2][0]=i0;
      bst[i0+ns*2][1]=i1;
      bst[i0+ns*3][0]=i0;
      bst[i0+ns*3][1]=i2;
    }
  }
  plq=new int*[ns];
  plq[0]=new int[ns*4];
  for(i=1;i<ns;i++)
    plq[i]=&(plq[0][i*4]);
  for(j=0;j<ly;j++)
    for(i=0;i<lx;i++){
      i0=i+j*lx;
      i1=(i+1)%lx+j*lx;
      i2=(i+1)%lx+((j+1)%ly)*lx;
      i3=i+((j+1)%ly)*lx;
      plq[i0][0]=i0;
      plq[i0][1]=i1;
      plq[i0][2]=i2;
      plq[i0][3]=i3;
    }
  txmap=new int[ns*(lx+2)];
  tymap=new int[ns*(ly+2)];
  pxmap=&(txmap[ns*lx]);//gliding reflection x map
  pymap=&(tymap[ns*ly]);//gliding reflection y map
  sigma1map=new int[ns*2];
  sigma2map=new int[ns*2];
  rotmap=new int[ns*4];
  for(i=0;i<lx;i++)
    for(j=0;j<ly;j++){
      i0=i+j*lx;
      i1=((lx-1-i)%lx)+j*lx;
      i2=i+((ly-1-j)%ly)*lx;
      i3=((i+1)%lx)+j*lx;
      i4=i+((j+1)%ly)*lx;
      pxmap[ns+i0]=i1;
      pymap[ns+i0]=i2;
      txmap[ns+i0]=i3;
      tymap[ns+i0]=i4;
      i1=j+i*lx;
      i2=((lx-1-j)%lx)+((ly-1-i)%ly)*lx;
      sigma1map[ns+i0]=i1;
      sigma2map[ns+i0]=i2;
      //cout<<"i0="<<i0<<"\t"<<i1<<"\t"<<i2<<endl;
      pxmap[i0]=i0;
      pymap[i0]=i0;
      txmap[i0]=i0;
      tymap[i0]=i0;
      sigma1map[i0]=i0;
      sigma2map[i0]=i0;
    }
  for(i=2;i<lx;i++)
    for(j=0;j<ns;j++)
      txmap[i*ns+j]=txmap[ns+txmap[(i-1)*ns+j]];
  for(i=2;i<ly;i++)
    for(j=0;j<ns;j++)
      tymap[i*ns+j]=tymap[ns+tymap[(i-1)*ns+j]];
  for(i=0;i<lx;i++)
    for(j=0;j<ly;j++){
      i0=(lx-1-j)%lx+i*lx;
      //i0=(lx-j)%lx+((i+1)%ly)*lx;//this is for ss model, rotate and TxTy
      i1=i+j*lx;
      rotmap[i1]=i1;
      rotmap[ns+i1]=i0;
      //cout<<"i1="<<i1<<"\trotmap="<<i0<<endl;
    }
  for(i=2;i<4;i++)
    for(j=0;j<ns;j++)
      rotmap[i*ns+j]=rotmap[ns+rotmap[(i-1)*ns+j]];
}

//--------------------------------------------------------------------------------------
void fci::build_hamiltonian_j1j3(){
//--------------------------------------------------------------------------------------
  int i,j,i0,i1,i2,i3,i4;
  ofstream fout;
  ns=lx*ly;
  nb=4*ns;
  bst=new int*[nb];
  bst[0]=new int[nb*2];
  for(i=1;i<nb;i++)
    bst[i]=&(bst[0][i*2]);
  for(i=0;i<lx;i++){
    for(j=0;j<ly;j++){
      i0=i+j*lx;
      i1=(i+1)%lx+j*lx;
      i2=i+((j+1)%ly)*lx;
      bst[i0][0]=i0;
      bst[i0][1]=i1;
      bst[i0+ns][0]=i0;
      bst[i0+ns][1]=i2;
      i1=(i+2)%lx+((j+0)%ly)*lx;
      i2=(i+0)%lx+((j+2)%ly)*lx;
      bst[i0+ns*2][0]=i0;
      bst[i0+ns*2][1]=i1;
      bst[i0+ns*3][0]=i0;
      bst[i0+ns*3][1]=i2;
    }
  }
  plq=new int*[ns];
  plq[0]=new int[ns*4];
  for(i=1;i<ns;i++)
    plq[i]=&(plq[0][i*4]);
  for(j=0;j<ly;j++)
    for(i=0;i<lx;i++){
      i0=i+j*lx;
      i1=(i+1)%lx+j*lx;
      i2=(i+1)%lx+((j+1)%ly)*lx;
      i3=i+((j+1)%ly)*lx;
      plq[i0][0]=i0;
      plq[i0][1]=i1;
      plq[i0][2]=i2;
      plq[i0][3]=i3;
    }
  txmap=new int[ns*(lx+2)];
  tymap=new int[ns*(ly+2)];
  pxmap=&(txmap[ns*lx]);//gliding reflection x map
  pymap=&(tymap[ns*ly]);//gliding reflection y map
  sigma1map=new int[ns*2];
  sigma2map=new int[ns*2];
  rotmap=new int[ns*4];
  for(i=0;i<lx;i++)
    for(j=0;j<ly;j++){
      i0=i+j*lx;
      i1=((lx-1-i)%lx)+j*lx;
      i2=i+((ly-1-j)%ly)*lx;
      i3=((i+1)%lx)+j*lx;
      i4=i+((j+1)%ly)*lx;
      pxmap[ns+i0]=i1;
      pymap[ns+i0]=i2;
      txmap[ns+i0]=i3;
      tymap[ns+i0]=i4;
      i1=j+i*lx;
      i2=((lx-1-j)%lx)+((ly-1-i)%ly)*lx;
      sigma1map[ns+i0]=i1;
      sigma2map[ns+i0]=i2;
      //cout<<"i0="<<i0<<"\t"<<i1<<"\t"<<i2<<endl;
      pxmap[i0]=i0;
      pymap[i0]=i0;
      txmap[i0]=i0;
      tymap[i0]=i0;
      sigma1map[i0]=i0;
      sigma2map[i0]=i0;
    }
  for(i=2;i<lx;i++)
    for(j=0;j<ns;j++)
      txmap[i*ns+j]=txmap[ns+txmap[(i-1)*ns+j]];
  for(i=2;i<ly;i++)
    for(j=0;j<ns;j++)
      tymap[i*ns+j]=tymap[ns+tymap[(i-1)*ns+j]];
  for(i=0;i<lx;i++)
    for(j=0;j<ly;j++){
      i0=(lx-1-j)%lx+i*lx;
      //i0=(lx-j)%lx+((i+1)%ly)*lx;//this is for ss model, rotate and TxTy
      i1=i+j*lx;
      rotmap[i1]=i1;
      rotmap[ns+i1]=i0;
      //cout<<"i1="<<i1<<"\trotmap="<<i0<<endl;
    }
  for(i=2;i<4;i++)
    for(j=0;j<ns;j++)
      rotmap[i*ns+j]=rotmap[ns+rotmap[(i-1)*ns+j]];
}

//--------------------------------------------------------------------------------------
void fci::hamiltonian(){
//--------------------------------------------------------------------------------------
  unsigned long i,j,t0,t1,sa,sb,sbb,b0,b1;
  int k,s0,s1,sx,sy,qx,qy,rz,qs1,qs2,tr;
  double ez;
  ofstream fout;
  for(i=0;i<nrep*nrep;i++)mat[i]=0;
  for(i=0;i<nrep;i++){
    sa=repr[i];
    ez=0;
    for(k=0;k<nb;k++){
      s0=bst[k][0];
      s1=bst[k][1];
      t0=(unsigned long)1<<s0;
      t1=(unsigned long)1<<s1;
      b0=sa&t0;
      b1=sa&t1;
      b0>>=s0;
      b1>>=s1;
      if(b0==b1)
	if(k>=2*ns)ez+=0.25*jcoup;
	else ez+=0.25;
      else{
	if(k>=2*ns)ez-=0.25*jcoup;
	else ez-=0.25;
	sb=(sa^t0)^t1;
	findrepresentative(sb,sbb,sx,sy,qx,qy,rz,qs1,qs2,tr);
	findstate(sbb,j);
	if(j!=nrep){
	  if(k>=2*ns)mat[j*nrep+i]+=0.5*jcoup*matrixelement(i,j,sx,sy,qx,qy,rz,qs1,qs2,tr);
	  else mat[j*nrep+i]+=0.5*matrixelement(i,j,sx,sy,qx,qy,rz,qs1,qs2,tr);
	}
      }
    }
    mat[i*nrep+i]+=ez;
    mat[i*nrep+i]-=hfldz*sec;
  }
}

//--------------------------------------------------------------------------------------
void fci::hamiltonianr(){
//--------------------------------------------------------------------------------------
  unsigned long i,j,t0,t1,sa,sb,sbb,b0,b1;
  int k,s0,s1,sx,sy,qx,qy,rz,qs1,qs2,tr;
  double ez;
  for(i=0;i<nrep*nrep;i++)matr[i]=0;

  for(i=0;i<nrep;i++){
    sa=repr[i];
    ez=0;
    for(k=0;k<nb;k++){
      s0=bst[k][0];
      s1=bst[k][1];
      t0=(unsigned long)1<<s0;
      t1=(unsigned long)1<<s1;
      b0=sa&t0;
      b1=sa&t1;
      b0>>=s0;
      b1>>=s1;
      if(b0==b1)
	if(k>=2*ns)ez+=0.25*jcoup;
	else ez+=0.25;
      else{
	if(k>=2*ns)ez-=0.25*jcoup;
	else ez-=0.25;
	sb=(sa^t0)^t1;
	findrepresentative(sb,sbb,sx,sy,qx,qy,rz,qs1,qs2,tr);
	findstate(sbb,j);
	if(j!=nrep){
	  if(k>=2*ns)matr[j*nrep+i]+=0.5*jcoup*matrixelementr(i,j,sx,sy,qx,qy,rz,qs1,qs2,tr);
	  else matr[j*nrep+i]+=0.5*matrixelementr(i,j,sx,sy,qx,qy,rz,qs1,qs2,tr);
	}
      }
    }
    matr[i*nrep+i]+=ez;
    matr[i*nrep+i]-=hfldz*sec;
  }
}

//--------------------------------------------------------------------------------------
void fci::hamiltonian_vector_multiplication(dcmplex* vec1, dcmplex* vec2){
//--------------------------------------------------------------------------------------
  unsigned long i,j,t0,t1,sa,sb,sbb,b0,b1;
  int k,s0,s1,sx,sy,qx,qy,myrank,rz,qs1,qs2,tr;
  double ez,time1,time2;
  ofstream fout;
  char sector[10],name[100],xx[10],yy[10],ppx[10],ppy[10],zz[10],ssx[10],ssy[10],rr[10];

  //#pragma omp parallel for
  for(i=0;i<nrep;i++)
    vec2[i]=0;
  
  //time1=omp_get_wtime();
  //#pragma omp parallel for default(shared) private(i,j,sa,sb,sbb,b0,b1,t0,t1,myrank,k,s0,s1,sx,sy,qx,qy,qs1,qs2,rz,tr,ez) schedule(dynamic,1)
  for(i=0;i<nrep;i++){
    //myrank=omp_get_thread_num();
    sa=repr[i];
    ez=0;
    for(k=0;k<nb;k++){
      s0=bst[k][0];
      s1=bst[k][1];
      t0=(unsigned long)1<<s0;
      t1=(unsigned long)1<<s1;
      b0=sa&t0;
      b1=sa&t1;
      b0>>=s0;
      b1>>=s1;
      if(b0==b1)
	if(k>=2*ns)ez+=0.25*jcoup;
	else ez+=0.25;
      else{
	if(k>=2*ns)ez-=0.25*jcoup;
	else ez-=0.25;
	sb=(sa^t0)^t1;
	findrepresentative(sb,sbb,sx,sy,qx,qy,rz,qs1,qs2,tr);
	findstate(sbb,j);
	if(j!=nrep){
	  //#pragma omp critical
	  {
	    if(k>=2*ns)vec2[j]+=0.5*jcoup*matrixelement(i,j,sx,sy,qx,qy,rz,qs1,qs2,tr)*vec1[i];
	    else vec2[j]+=0.5*matrixelement(i,j,sx,sy,qx,qy,rz,qs1,qs2,tr)*vec1[i];
	  }
	}
      }
    }
    //#pragma omp critical
    {
      vec2[i]+=ez*vec1[i];
      vec2[i]-=hfldz*sec*vec1[i];
    }
  }
  //time2=omp_get_wtime();
  sprintf(sector,"%d",sec);
  sprintf(xx,"%d",kx);
  sprintf(yy,"%d",ky);
  sprintf(ppx,"%d",px);
  sprintf(ppy,"%d",py);
  sprintf(zz,"%d",qz);
  sprintf(ssx,"%d",sgm1);
  sprintf(ssy,"%d",sgm2);
  sprintf(rr,"%d",qr);
  strcpy(name,"out-");
  strcat(name,sector);
  strcat(name,"-");
  strcat(name,xx);
  strcat(name,"-");
  strcat(name,yy);
  strcat(name,"-");
  strcat(name,ppx);
  strcat(name,"-");
  strcat(name,ppy);
  strcat(name,"-");
  strcat(name,zz);
  strcat(name,"-");
  strcat(name,ssx);
  strcat(name,"-");
  strcat(name,ssy);
  strcat(name,"-");
  strcat(name,rr);
  strcat(name,".dat");
  fout.open(name,ios::app);
  fout<<"Hv time="<<time2-time1<<endl;
  fout.close();
}

//--------------------------------------------------------------------------------------
void fci::hamiltonian_vector_multiplication(double* vec1, double* vec2){
//--------------------------------------------------------------------------------------
  unsigned long i,j,sa,sb,sbb,b0,b1;
  int k,s0,s1,sx,sy,qx,qy,myrank,rz,qs1,qs2,tr;
  double ez;
  double tm1,tm2,tm3,timefindstate[100],timefindrepr[100],time1,time2;
  double **vecdup;
  ofstream fout;
  char sector[10],name[100],xx[10],yy[10],ppx[10],ppy[10],zz[10],ssx[10],ssy[10],rr[10];
  vecdup=new double*[psize];
  for(k=0;k<psize;k++)
    vecdup[k]=new double[nrep];
  for(i=0;i<psize;i++){
    timefindstate[i]=0;
    timefindrepr[i]=0;
  }
  //#pragma omp parallel for default(shared) private(k)
  for(i=0;i<nrep;i++){
    vec2[i]=0;
    for(k=0;k<psize;k++)
      vecdup[k][i]=0;
  }
  //time1=omp_get_wtime();
  //#pragma omp parallel for default(shared) private(i,j,sa,sb,sbb,b0,b1,myrank,k,s0,s1,sx,sy,qx,qy,qs1,qs2,rz,tr,ez,tm1,tm2,tm3) schedule(dynamic,1)
  //#pragma omp parallel for default(shared) private(i,j,sa,sb,sbb,b0,b1,myrank,k,s0,s1,sx,sy,qx,qy,qs1,qs2,rz,ez) schedule(dynamic,1)
  for(i=0;i<nrep;i++){
    //myrank=omp_get_thread_num();
    sa=repr[i];
    ez=0;
    for(k=0;k<nb;k++){
      s0=bst[k][0];
      s1=bst[k][1];
      b0=sa&bitmap[s0];
      b1=sa&bitmap[s1];
      b0>>=s0;
      b1>>=s1;
      if(b0==b1)
	if(k>=2*ns)ez+=0.25*jcoup;
	else ez+=0.25;
      else{
	if(k>=2*ns)ez-=0.25*jcoup;
	else ez-=0.25;
	sb=(sa^bitmap[s0])^bitmap[s1];
	//tm1=omp_get_wtime();
	findrepresentative(sb,sbb,sx,sy,qx,qy,rz,qs1,qs2,tr);
	//tm2=omp_get_wtime();
	findstate(sbb,j);
	//tm3=omp_get_wtime();
	timefindstate[myrank]+=tm3-tm2;
	timefindrepr[myrank]+=tm2-tm1;
	if(j!=nrep){
	  if(k>=2*ns)vecdup[myrank][j]+=0.5*jcoup*matrixelementr(i,j,sx,sy,qx,qy,rz,qs1,qs2,tr)*vec1[i];
	  else vecdup[myrank][j]+=0.5*matrixelementr(i,j,sx,sy,qx,qy,rz,qs1,qs2,tr)*vec1[i];
	}
      }
    }
    vecdup[myrank][i]+=(ez-hfldz*sec)*vec1[i];
  }
  //time2=omp_get_wtime();
  cout<<"aHv time="<<time2-time1<<endl;
  //#pragma omp parallel for default(shared) private(k)
  for(i=0;i<nrep;i++)
    for(k=0;k<psize;k++)
      vec2[i]+=vecdup[k][i];
  for(k=0;k<psize;k++)
    delete []vecdup[k];
  delete []vecdup;
  tm1=0;
  tm2=0;
  for(i=0;i<psize;i++){
    tm1+=timefindstate[i];
    tm2+=timefindrepr[i];
  }
  cout<<"time findstate="<<tm1<<"\ttime findrepr="<<tm2<<endl;
  sprintf(sector,"%d",sec);
  sprintf(xx,"%d",kx);
  sprintf(yy,"%d",ky);
  sprintf(ppx,"%d",px);
  sprintf(ppy,"%d",py);
  sprintf(zz,"%d",qz);
  sprintf(ssx,"%d",sgm1);
  sprintf(ssy,"%d",sgm2);
  sprintf(rr,"%d",qr);
  strcpy(name,"out-");
  strcat(name,sector);
  strcat(name,"-");
  strcat(name,xx);
  strcat(name,"-");
  strcat(name,yy);
  strcat(name,"-");
  strcat(name,ppx);
  strcat(name,"-");
  strcat(name,ppy);
  strcat(name,"-");
  strcat(name,zz);
  strcat(name,"-");
  strcat(name,ssx);
  strcat(name,"-");
  strcat(name,ssy);
  strcat(name,"-");
  strcat(name,rr);
  strcat(name,".dat");
  fout.open(name,ios::app);
  fout<<"time findstate="<<tm1<<"\ttime findrepr="<<tm2<<"\tHv time="<<time2-time1<<endl;
  fout.close();
}

//--------------------------------------------------------------------------------------
void fci::measure(){
//--------------------------------------------------------------------------------------
  int i;
  if(kx%(perix/2)==0&&ky%(periy/2)==0&&qr%2==0){ //use real hamiltonian matrix
    if(nrep>=nrep_ed_max){
      for(i=0;i<neig;i++)
	totspn[i]=measure_total_spin(eigvecr[i]);
    }
    else 
      for(i=0;i<nrep;i++)
	totspn[i]=measure_total_spin(&(matr[i*nrep]));
  }
  else{
    if(nrep>=nrep_ed_max){
      for(i=0;i<neig;i++)
	totspn[i]=measure_total_spin(eigvec[i]);
    }
    else
      for(i=0;i<nrep;i++)
	totspn[i]=measure_total_spin(&(mat[i*nrep]));
  }
}
//--------------------------------------------------------------------------------------
void fci::measure_bond_enr(double* vec1){
//--------------------------------------------------------------------------------------
  unsigned long i,j,sa,sb,sbb,b0,b1;
  int k,l,s0,s1,sx,sy,qx,qy,myrank,rz,qs1,qs2,tr;
  double ez;
  double tm1,tm2,tm3,tm0,tm4,**bondenr,time1,time2;
  double *normtmp,phase1;
  ofstream fout;
  char sector[10],name[100],xx[10],yy[10],ppx[10],ppy[10],zz[10],ssx[10],ssy[10],rr[10];
  normtmp=new double[psize];
  bondenr=new double*[4];
  bondenr[0]=new double[psize*4];
  for(k=1;k<4;k++)
    bondenr[k]=&(bondenr[0][k*psize]);
  for(k=0;k<psize;k++){
    bondenr[0][k]=0;
    bondenr[1][k]=0;
    bondenr[2][k]=0;
    bondenr[3][k]=0;
    normtmp[k]=0;
  }
  //time1=omp_get_wtime();
  //#pragma omp parallel for default(shared) private(i,j,sa,sb,sbb,b0,b1,myrank,k,l,s0,s1,sx,sy,qx,qy,qs1,qs2,rz,tr,ez) schedule(dynamic,1)
  for(i=0;i<nrep;i++){
    //myrank=omp_get_thread_num();
    sa=repr[i];
    normtmp[myrank]+=vec1[i]*vec1[i];
    for(k=0;k<ntrans;k++){
      for(l=0;l<4;l++){
	s0=mapptr[k][bst[l*ns][0]];
	s1=mapptr[k][bst[l*ns][1]];
	b0=sa&bitmap[s0];
	b1=sa&bitmap[s1];
	b0>>=s0;
	b1>>=s1;
	if(b0==b1){
	  if(abs(qz)==1)ez=0.5;
	  else ez=0.25;
	}
	else{
	  if(abs(qz)==1)ez=-0.5;
	  else ez=-0.25;
	  sb=(sa^bitmap[s0])^bitmap[s1];
	  findrepresentative(sb,sbb,sx,sy,qx,qy,rz,qs1,qs2,tr);
	  findstate(sbb,j);
	  if(j!=nrep){
	    if(abs(qz)==1){
	      bondenr[l][myrank]+=0.5*get_phase(i,j,k,sx,sy,qx,qy,rz,qs1,qs2,tr)*vec1[i]*vec1[j];
	      bondenr[l][myrank]+=qz*0.5*get_phase(i,j,k,sx,sy,qx,qy,(rz+1)%2,qs1,qs2,tr)*vec1[i]*vec1[j];
	    }
	    else
	      bondenr[l][myrank]+=0.5*get_phase(i,j,k,sx,sy,qx,qy,rz,qs1,qs2,tr)*vec1[i]*vec1[j];
	  }
	}
	bondenr[l][myrank]+=ez*vec1[i]*vec1[i];
      }
    }
  }
  //time2=omp_get_wtime();
  cout<<"measure bond energy time="<<time2-time1<<endl;
  tm0=0;
  tm1=0;
  tm2=0;
  tm3=0;
  tm4=0;
  for(i=0;i<psize;i++){
    tm0+=bondenr[0][i];
    tm1+=bondenr[1][i];
    tm2+=bondenr[2][i];
    tm3+=bondenr[3][i];
    tm4+=normtmp[i];
  }
  sprintf(sector,"%d",sec);
  sprintf(xx,"%d",kx);
  sprintf(yy,"%d",ky);
  sprintf(ppx,"%d",px);
  sprintf(ppy,"%d",py);
  sprintf(zz,"%d",qz);
  sprintf(ssx,"%d",sgm1);
  sprintf(ssy,"%d",sgm2);
  sprintf(rr,"%d",qr);
  strcpy(name,"bondenr-");
  strcat(name,sector);
  strcat(name,"-");
  strcat(name,xx);
  strcat(name,"-");
  strcat(name,yy);
  strcat(name,"-");
  strcat(name,ppx);
  strcat(name,"-");
  strcat(name,ppy);
  strcat(name,"-");
  strcat(name,zz);
  strcat(name,"-");
  strcat(name,ssx);
  strcat(name,"-");
  strcat(name,ssy);
  strcat(name,"-");
  strcat(name,rr);
  strcat(name,".dat");
  fout.open(name,ios::out);
  fout<<setprecision(12)<<tm0/(ns*ns)<<"\t"<<tm1/(ns*ns)<<"\t"<<tm2/(ns*ns)<<"\t"<<tm3/(ns*ns)<<"\t"<<(tm0+tm1+(tm2+tm3)*jcoup)/(double)(ns*ns)<<"\t"<<tm4<<endl;
  fout.close();
  delete []bondenr[0];
  delete []bondenr;
  delete []normtmp;
}

//--------------------------------------------------------------------------------------
double fci::measure_rotation_quantum_number(double* vec1){
//--------------------------------------------------------------------------------------
  unsigned long i,j,sa,sb,sbb;
  int sx,sy,qx,qy,myrank,rz,qs1,qs2,tr;
  double *totsdup,tots;
  totsdup=new double[psize];
  for(i=0;i<psize;i++)
    totsdup[i]=0;
  //#pragma omp parallel for default(shared) private(i,j,sa,sb,sbb,myrank,sx,sy,qx,qy,qs1,qs2,rz,tr) schedule(dynamic,1)
  for(i=0;i<nrep;i++){
    //myrank=omp_get_thread_num();
    sa=repr[i];
    rotaterepr(sa,sb);
    findrepresentative(sb,sbb,sx,sy,qx,qy,rz,qs1,qs2,tr);
    findstate(sbb,j);
    if(j!=nrep)
      totsdup[myrank]+=matrixelementr(i,j,sx,sy,qx,qy,rz,qs1,qs2,tr)*vec1[i]*vec1[j];
  }
  tots=0;
  for(i=0;i<psize;i++)
    tots+=totsdup[i];
  delete []totsdup;
  return tots;
}

//--------------------------------------------------------------------------------------
double fci::measure_total_spin(double* vec1){
//--------------------------------------------------------------------------------------
  unsigned long i,j,t0,t1,sa,sb,sbb,b0,b1;
  int k,s0,s1,sx,sy,qx,qy,myrank,rz,qs1,qs2,tr;
  double ez;
  double tm1,tm2,tm3,timefindstate[100],timefindrepr[100];
  double *totsdup,tots;
  ofstream fout;
  char sector[10],name[100],xx[10],yy[10],ppx[10],ppy[10],zz[10],ssx[10],ssy[10],rr[10];
  totsdup=new double[psize];
  for(i=0;i<psize;i++){
    timefindstate[i]=0;
    timefindrepr[i]=0;
    totsdup[i]=0;
  }

  //#pragma omp parallel for default(shared) private(i,j,sa,sb,sbb,b0,b1,t0,t1,myrank,k,s0,s1,sx,sy,qx,qy,qs1,qs2,rz,tr,ez,tm1,tm2,tm3) schedule(dynamic,1)
  for(i=0;i<nrep;i++){
    //myrank=omp_get_thread_num();
    sa=repr[i];
    ez=0;
    for(s0=0;s0<ns;s0++)
      for(s1=s0+1;s1<ns;s1++){
	t0=(unsigned long)1<<s0;
	t1=(unsigned long)1<<s1;
	b0=sa&t0;
	b1=sa&t1;
	b0>>=s0;
	b1>>=s1;
	if(b0==b1)
	  ez+=0.25;
	else{
	  ez-=0.25;
	  sb=(sa^t0)^t1;
	  //tm1=omp_get_wtime();
	  findrepresentative(sb,sbb,sx,sy,qx,qy,rz,qs1,qs2,tr);
	  //tm2=omp_get_wtime();
	  findstate(sbb,j);
	  //tm3=omp_get_wtime();
	  timefindstate[myrank]+=tm3-tm2;
	  timefindrepr[myrank]+=tm2-tm1;
	  if(j!=nrep)
	    totsdup[myrank]+=0.5*matrixelementr(i,j,sx,sy,qx,qy,rz,qs1,qs2,tr)*vec1[i]*vec1[j];
	}
      }
    totsdup[myrank]+=ez*vec1[i]*vec1[i];
  }
  
  tm1=0;
  tm2=0;
  tots=0;
  for(i=0;i<psize;i++){
    tm1+=timefindstate[i];
    tm2+=timefindrepr[i];
    tots+=totsdup[i];
  }
  tots*=2.;
  tots+=(double)ns*0.75;
  delete []totsdup;
  sprintf(sector,"%d",sec);
  sprintf(xx,"%d",kx);
  sprintf(yy,"%d",ky);
  sprintf(ppx,"%d",px);
  sprintf(ppy,"%d",py);
  sprintf(zz,"%d",qz);
  sprintf(ssx,"%d",sgm1);
  sprintf(ssy,"%d",sgm2);
  sprintf(rr,"%d",qr);
  strcpy(name,"out-");
  strcat(name,sector);
  strcat(name,"-");
  strcat(name,xx);
  strcat(name,"-");
  strcat(name,yy);
  strcat(name,"-");
  strcat(name,ppx);
  strcat(name,"-");
  strcat(name,ppy);
  strcat(name,"-");
  strcat(name,zz);
  strcat(name,"-");
  strcat(name,ssx);
  strcat(name,"-");
  strcat(name,ssy);
  strcat(name,"-");
  strcat(name,rr);
  strcat(name,".dat");
  fout.open(name,ios::app);
  fout<<"time findstate="<<tm1<<"\ttime findrepr="<<tm2<<"\t tots="<<tots<<endl;
  fout.close();
  return tots;
}

//--------------------------------------------------------------------------------------
double fci::measure_total_spin(dcmplex* vec1){
//--------------------------------------------------------------------------------------
  unsigned long i,j,t0,t1,sa,sb,sbb,b0,b1;
  int k,s0,s1,sx,sy,qx,qy,myrank,rz,qs1,qs2,tr;
  double ez;
  double tm1,tm2,tm3,timefindstate[100],timefindrepr[100];
  double *totsdup,tots,nor1,nor2,overlap;
  dcmplex *vtmp;
  totsdup=new double[psize];
  for(i=0;i<psize;i++){
    timefindstate[i]=0;
    timefindrepr[i]=0;
    totsdup[i]=0;
  }
  vtmp=new dcmplex[nrep];
  hamiltonian_vector_multiplication(vec1,vtmp);
  nor1=real(inner_prod(nrep,vtmp,vtmp));
  nor2=real(inner_prod(nrep,vtmp,vec1));
  overlap=fabs(nor2/sqrt(nor1));
  cout<<setprecision(12)<<j<<"\tcheck eigenvector eigval="<<sqrt(nor1)<<"\toverlap="<<overlap<<endl;  
  delete []vtmp;
  //#pragma omp parallel for default(shared) private(i,j,sa,sb,sbb,b0,b1,t0,t1,myrank,k,s0,s1,sx,sy,qx,qy,qs1,qs2,rz,tr,ez,tm1,tm2,tm3) schedule(dynamic,1)
  for(i=0;i<nrep;i++){
    //myrank=omp_get_thread_num();
    sa=repr[i];
    ez=0;
    for(s0=0;s0<ns;s0++)
      for(s1=s0+1;s1<ns;s1++){
	t0=(unsigned long)1<<s0;
	t1=(unsigned long)1<<s1;
	b0=sa&t0;
	b1=sa&t1;
	b0>>=s0;
	b1>>=s1;
	if(b0==b1)
	  ez+=0.25;
	else{
	  ez-=0.25;
	  sb=(sa^t0)^t1;
	  //tm1=omp_get_wtime();
	  findrepresentative(sb,sbb,sx,sy,qx,qy,rz,qs1,qs2,tr);
	  //tm2=omp_get_wtime();
	  findstate(sbb,j);
	  //tm3=omp_get_wtime();
	  timefindstate[myrank]+=tm3-tm2;
	  timefindrepr[myrank]+=tm2-tm1;
	  if(j!=nrep)
	    totsdup[myrank]+=real(0.5*matrixelement(i,j,sx,sy,qx,qy,rz,qs1,qs2,tr)*vec1[i]*conj(vec1[j]));
	}
      }
    totsdup[myrank]+=real(ez*vec1[i]*conj(vec1[i]));
  }
  
  tm1=0;
  tm2=0;
  tots=0;
  for(i=0;i<psize;i++){
    tm1+=timefindstate[i];
    tm2+=timefindrepr[i];
    tots+=totsdup[i];
  }
  tots*=2.;
  tots+=(double)ns*0.75;
  delete []totsdup;
  cout<<"time findstate="<<tm1<<"\ttime findrepr="<<tm2<<"\t tots="<<tots<<endl;
  return tots;
}

//--------------------------------------------------------------------------------------
void obtain_hermitian_matrix_eigenvector1(dcmplex *mtr, double *w, int mdim,int vflag){
//--------------------------------------------------------------------------------------
  char jobz,uplo;
  dcmplex *work;
  double *rwork;
  int info,lwork,i;

  lwork=-1;
  work=new dcmplex[2];
  rwork=new double[3*mdim];
  if(vflag==0)
    zheev_("N","U",&mdim,mtr,&mdim,w,work,&lwork,rwork,&info);
  else
    zheev_("V","U",&mdim,mtr,&mdim,w,work,&lwork,rwork,&info);
  lwork=(int)real(work[0]);
  delete []work;
  work=new dcmplex[lwork];
  if(vflag==0)
    zheev_("N","U",&mdim,mtr,&mdim,w,work,&lwork,rwork,&info);
  else
    zheev_("V","U",&mdim,mtr,&mdim,w,work,&lwork,rwork,&info);
  if(info!=0){
    cout<<"zheev info="<<info<<endl;
  }
  delete []work;
  delete []rwork;
  if(info!=0){
    /*
    for(i=0;i<mdim;i++){
      cout<<i<<"\t"<<w[i]<<endl;
    }
    */
    exit(0);
  }
}

//--------------------------------------------------------------------------------------
void obtain_symmetric_matrix_eigenvector1(double *mtr, double *w, int mdim, int flag){
//--------------------------------------------------------------------------------------
  //this function destroy the symmetric matrix mtr and output its eigenvectors
  int i,info,lwork;
  double *work;
  lwork=-1;
  work=new double[2];
  if(flag==0)
    dsyev_("N","U",&mdim,mtr,&mdim,w,work,&lwork,&info);
  else if(flag==1)
    dsyev_("V","U",&mdim,mtr,&mdim,w,work,&lwork,&info);
  lwork=(int)(work[0]);
  cout<<"lwork="<<lwork<<"\t"<<work[0]<<endl;
  delete []work;
  work=new double[lwork];
  if(flag==0)
    dsyev_("N","U",&mdim,mtr,&mdim,w,work,&lwork,&info);
  else if(flag==1)
    dsyev_("V","U",&mdim,mtr,&mdim,w,work,&lwork,&info);
  if(info!=0){
    cout<<"dsyev info="<<info<<endl;
  }
  delete []work;
  if(info!=0){
    exit(0);
  }
}

//--------------------------------------------------------------------------------------
void obtain_hermitian_matrix_eigenvector(dcmplex *mtr, double *w, int mdim, int dcut0){
//--------------------------------------------------------------------------------------
  //this function destroy the hermitian matrix mtr and output its eigenvectors
  char jobz,uplo,transa,transb;
  int n,lda,lwork,i,info,il,iu,m,*iwork,*ifail;
  dcmplex *work,alpha=1,beta=0,*z;
  double *rwork,vl,vu,abstol=0;

  il=mdim-dcut0+1;
  iu=mdim;
  lwork=-1;
  work=new dcmplex[2];
  rwork=new double[7*mdim];
  iwork=new int[5*mdim];
  ifail=new int[mdim];
  z=new dcmplex[mdim*dcut0];
  zheevx_("V","I","U",&mdim,mtr,&mdim,&vl,&vu,&il,&iu,&abstol,&m,w,z,&mdim,work,&lwork,rwork,iwork,ifail,&info);
  lwork=(int)real(work[0]);
  delete []work;
  lwork=lwork+mdim*100;
  work=new dcmplex[lwork];
  zheevx_("V","I","U",&mdim,mtr,&mdim,&vl,&vu,&il,&iu,&abstol,&m,w,z,&mdim,work,&lwork,rwork,iwork,ifail,&info);
  if(info!=0){
    cout<<"zheevx info="<<info<<endl;
    //exit(0);
  }
  for(i=0;i<mdim*dcut0;i++)
    mtr[i]=z[i];
  delete []work;
  delete []rwork;
  delete []iwork;
  delete []ifail;
  delete []z;
  if(info!=0){
    for(i=0;i<mdim;i++){
      cout<<i<<"\t"<<w[i]<<endl;
    }
    exit(0);
  }
}

//--------------------------------------------------------------------------------------
void fci::diagonalize(int vflag,int eflag){
//--------------------------------------------------------------------------------------
  unsigned long i;
  int j;
  ofstream fout;
  char base[100],name[100],sector[10],xx[10],yy[10],ppx[10],ppy[10],zz[10],ssx[10],ssy[10],rr[10];
  sprintf(sector,"%d",sec);
  sprintf(xx,"%d",kx);
  sprintf(yy,"%d",ky);
  sprintf(ppx,"%d",px);
  sprintf(ppy,"%d",py);
  sprintf(zz,"%d",qz);
  sprintf(ssx,"%d",sgm1);
  sprintf(ssy,"%d",sgm2);
  sprintf(rr,"%d",qr);
  strcpy(name,"eigval-");
  strcat(name,sector);
  strcat(name,"-");
  strcat(name,xx);
  strcat(name,"-");
  strcat(name,yy);
  strcat(name,"-");
  strcat(name,ppx);
  strcat(name,"-");
  strcat(name,ppy);
  strcat(name,"-");
  strcat(name,zz);
  strcat(name,"-");
  strcat(name,ssx);
  strcat(name,"-");
  strcat(name,ssy);
  strcat(name,"-");
  strcat(name,rr);
  strcat(name,".dat");
  if(kx%(perix/2)==0&&ky%(periy/2)==0&&qr%2==0) //use real hamiltonian matrix
    obtain_symmetric_matrix_eigenvector1(matr,eigval,nrep,vflag);
  else
    obtain_hermitian_matrix_eigenvector1(mat,eigval,nrep,vflag);
  fout.open(name,ios::out);
  if(fout.is_open()&&eflag==0){
    for(i=0;i<nrep;i++)
      fout<<i<<"\t"<<setprecision(12)<<eigval[i]<<endl;
    fout.close();
  }
}

//--------------------------------------------------------------------------------------
void fci::checkhermitian(){
//--------------------------------------------------------------------------------------
  unsigned long i,j,a,b;
  ofstream fout;
  fout.open("out.dat",ios::app);
  for(i=0;i<nrep;i++)
    for(j=i;j<nrep;j++){
      if(norm(mat[i*nrep+j]-conj(mat[j*nrep+i]))>=eps){
	fout<<setprecision(12)<<"complex fail\t"<<i<<"\t"<<j<<"\t"<<normfac[i]<<"\t"<<normfac[j]<<"\t"<<mat[j*nrep+i]<<"\t"<<mat[i*nrep+j]<<"\t"<<norm(mat[i*nrep+j]-conj(mat[j*nrep+i]))<<endl;
	exit(0);
      }
    }
  fout.close();
}

//--------------------------------------------------------------------------------------
void fci::checkhermitianr(){
//--------------------------------------------------------------------------------------
  unsigned long i,j,a,b;
  ofstream fout;
  fout.open("out.dat",ios::app);
  for(i=0;i<nrep;i++)
    for(j=i;j<nrep;j++){
      if(fabs(matr[i*nrep+j]-matr[j*nrep+i])>=eps){
	fout<<setprecision(12)<<"real fail\t"<<i<<"\t"<<j<<"\t"<<normfac[i]<<"\t"<<normfac[j]<<"\t"<<matr[j*nrep+i]<<"\t"<<matr[i*nrep+j]<<"\t"<<fabs(matr[i*nrep+j]-matr[j*nrep+i])<<endl;
	exit(0);
      }
    }
  fout.close();
}

//--------------------------------------------------------------------------------------
void fci::findstate(unsigned long rep, unsigned long &b){
//--------------------------------------------------------------------------------------
  unsigned long bmin,bmax,n1,n2;
  n1=lookuptable[rep%ntable];
  n2=lookuptable[(rep%ntable)+1];
  if(n1==n2||rep<repr[n1]||rep>repr[n2-1]){
    b=nrep;
    return;
  }
  bmin=n1;
  bmax=n2-1;
  if(repr[bmin]==rep){
    b=bmin;
    return;
  }
  if(repr[bmax]==rep){
    b=bmax;
    return;
  }
  do{
    b=bmin+(bmax-bmin)/2;
    if(bmin==bmax&&rep!=repr[b]){
      b=nrep;
      return;
    }
    if(rep<repr[b]) 
      bmax=b-1;
    else if(rep>repr[b])
      bmin=b+1;
    else return;
  }while(bmin<=bmax);
  b=nrep;
}

//--------------------------------------------------------------------------------------
bool fci::makefgfunctions(){
//--------------------------------------------------------------------------------------
  double kkx,kky,ppr;
  double pi2=2.*3.14159265358979;
  int i,j,k,l;
  dcmplex cnum;
  kkx=(double)kx*pi2/(double)perix;
  kky=(double)ky*pi2/(double)periy;
  ppr=(double)qr*pi2/(double)4;

  f4func=new dcmplex**[ly*2+1];
  for(j=0;j<ly*2+1;j++){
    f4func[j]=new dcmplex*[lx*2+1];
    for(i=0;i<lx*2+1;i++)
      f4func[j][i]=new dcmplex[4];
  }
  for(j=-ly;j<=ly;j++)
    for(i=-lx;i<=lx;i++)
      for(k=0;k<4;k++){
	cnum=std::complex<double>(cos(kkx*i+kky*j+ppr*k),-sin(kkx*i+kky*j+ppr*k));
	f4func[j][i][k]=cnum;
      }
  return true;
}
  
//--------------------------------------------------------------------------------------
dcmplex fci::matrixelement(unsigned long a, unsigned long b, int sx, int sy,int lpx, int lpy,int tz,int ls1,int ls2,int tr){
//--------------------------------------------------------------------------------------
  dcmplex mab;
  int fac;
  fac=0;
  if(sec==0&&qz==-1&&tz==1)fac++;
  if(sgm1==-1&&ls1==1)fac++;
  if(sgm2==-1&&ls2==1)fac++;
  if(px==-1&&lpx==1)fac++;
  if(py==-1&&lpy==1)fac++;
  mab=sqrt(normfac[b]/normfac[a]);
  mab*=f4func[ly+sy][lx+sx][tr];
  if((fac%2)==0)return mab;
  else return -mab;
}

//--------------------------------------------------------------------------------------
double fci::matrixelementr(unsigned long a, unsigned long b, int sx, int sy,int lpx, int lpy,int tz,int ls1,int ls2,int tr){
//--------------------------------------------------------------------------------------
  double mab;
  int fac;
  fac=0;
  if(sec==0&&qz==-1&&tz==1)fac+=1;
  if(sgm1==-1&&ls1==1)fac+=1;
  if(sgm2==-1&&ls2==1)fac+=1;
  if(px==-1&&lpx==1)fac++;
  if(py==-1&&lpy==1)fac++;
  mab=sqrt(normfac[b]/normfac[a]);
  mab*=real(f4func[ly+sy][lx+sx][tr]);
  if((fac%2)==0)return mab;
  else return -mab;
}

//--------------------------------------------------------------------------------------
double fci::get_phase(unsigned long a, unsigned long b, int i,int sx, int sy,int lpx, int lpy,int tz,int ls1,int ls2,int tr){
//--------------------------------------------------------------------------------------
  double mab;
  int fac,tx,ty,tpx,tpy,tsx,tsy,trr;
  fac=0;
  if(sec==0&&qz==-1&&tz==1)fac+=1;
  if(sgm1==-1&&ls1==1)fac+=1;
  if(sgm2==-1&&ls2==1)fac+=1;
  if(px==-1&&lpx==1)fac++;
  if(py==-1&&lpy==1)fac++;
  mab=sqrt(normfac[b]/normfac[a]);
  mab*=real(f4func[ly+sy][lx+sx][tr]);
  if((fac%2)==0)return mab;
  else return -mab;
}

//--------------------------------------------------------------------------------------
void fci::build_basis(int nup){
//--------------------------------------------------------------------------------------
//kx,ky are good quantum number
  unsigned long i,i0,maxcode;
  int j,n1,*nrp,myrank;
  double multi;
  bool pass;
  ofstream fout;
  nrp=new int[psize];
  for(i=0;i<psize;i++)
    nrp[i]=0;
  maxcode=0;
  for(j=0;j<ns;j++)maxcode=maxcode|((unsigned long)1<<j);
  //#pragma omp parallel for default(shared) private(i,i0,j,n1,myrank,multi,pass) schedule(dynamic,1024)
  for(i=0;i<=maxcode;i++){
    //myrank=omp_get_thread_num();
    n1=0;
    i0=i;
    for(j=0;j<ns;j++){
      if(i0%2==1) n1++;
      i0/=2;
    }
    if(n1!=nup)continue;
    if(sec==0)
      pass=checkrepresentative2(i,multi);
    else
      pass=checkrepresentative1(i,multi);
    if(pass&&fabs(multi)>eps)
      nrp[myrank]++;
  }
  nrep=0;
  for(i=0;i<psize;i++)
    nrep+=nrp[i];
  nbasis=nrep;
  cout<<"nrep="<<nrep<<endl;
  delete []nrp;
}

//--------------------------------------------------------------------------------------
void fci::build_basis1(int nup){
//--------------------------------------------------------------------------------------
  unsigned long i,i0,maxcode,**reprdup,*nrepdup;
  double **normdup;
  int j,n1,myrank;
  double multi;
  bool pass;
  ofstream fout;
  char sector[10],name[100],xx[10],yy[10],ppx[10],ppy[10],zz[10],ssx[10],ssy[10],rr[10];

  reprdup=new unsigned long*[psize];
  normdup=new double*[psize];
  nrepdup=new unsigned long[psize];
  for(j=0;j<psize;j++){
    reprdup[j]=new unsigned long[nrep];
    normdup[j]=new double[nrep];
    nrepdup[j]=0;
  }
  maxcode=0;
  for(j=0;j<ns;j++)maxcode=maxcode|((unsigned long)1<<j);
  nrep=0;
  //#pragma omp parallel for default(shared) private(i,i0,j,n1,myrank,multi,pass) schedule(dynamic,1024)
  for(i=0;i<=maxcode;i++){
    //myrank=omp_get_thread_num();
    n1=0;
    i0=i;
    for(j=0;j<ns;j++){
      if(i0%2==1) n1++;
      i0/=2;
    }
    if(n1!=nup)continue;
    if(sec==0)
      pass=checkrepresentative2(i,multi);
    else
      pass=checkrepresentative1(i,multi);
    if(pass&&fabs(multi)>eps){
      reprdup[myrank][nrepdup[myrank]]=i;
      normdup[myrank][nrepdup[myrank]]=multi;
      nrepdup[myrank]++;
    }
  }
  nrep=0;
  for(j=0;j<psize;j++){
    for(i=0;i<nrepdup[j];i++){
      repr[nrep]=reprdup[j][i];
      normfac[nrep]=normdup[j][i];
      nrep++;
    }
    delete []reprdup[j];
    delete []normdup[j];
  }
  delete []reprdup;
  delete []normdup;
  delete []nrepdup;
  cout<<"nrep="<<nrep<<endl;
  nbasis=nrep;
  make_lookup_table();
  sprintf(sector,"%d",sec);
  sprintf(xx,"%d",kx);
  sprintf(yy,"%d",ky);
  sprintf(ppx,"%d",px);
  sprintf(ppy,"%d",py);
  sprintf(zz,"%d",qz);
  sprintf(ssx,"%d",sgm1);
  sprintf(ssy,"%d",sgm2);
  sprintf(rr,"%d",qr);
  strcpy(name,"out-");
  strcat(name,sector);
  strcat(name,"-");
  strcat(name,xx);
  strcat(name,"-");
  strcat(name,yy);
  strcat(name,"-");
  strcat(name,ppx);
  strcat(name,"-");
  strcat(name,ppy);
  strcat(name,"-");
  strcat(name,zz);
  strcat(name,"-");
  strcat(name,ssx);
  strcat(name,"-");
  strcat(name,ssy);
  strcat(name,"-");
  strcat(name,rr);
  strcat(name,".dat");
  fout.open(name,ios::app);
  fout<<"nrep="<<nrep<<"\t"<<sec<<"\t"<<kx<<"\t"<<ky<<"\t"<<px<<"\t"<<py<<"\t"<<qz<<"\t"<<sgm1<<"\t"<<sgm2<<"\t"<<qr<<endl;
  fout.close();
}

//--------------------------------------------------------------------------------------
bool fci::checkrepresentative1(unsigned long s0, double& multi){
//--------------------------------------------------------------------------------------
  int i,j,k,l,m,n,nvisit,tx,ty,tpx,tpy,tsx,tsy,sgn,tr,rx,ry;
  bool *visit,zeroflag;
  double tmp;
  if(findrxry(s0,rx,ry,true)==false)return false;
  visit=new bool[ntransc];
  nvisit=0;
  for(i=0;i<ntransc;i++){
    if(transtep[0][i]<rx&&transtep[1][i]<ry){    
      visit[i]=true;
      nvisit+=1;
    }
    else
      visit[i]=false;
  }
  for(j=ns-1;j>=0;j--){
    zeroflag=false;
    for(i=0;i<ntransc;i++)
      if(visit[i]&&((s0&bitmap[mapptr[i][j]])==0)){
	zeroflag=true;
	break;
      }
    if(zeroflag&&nvisit>1)
      for(i=0;i<ntransc;i++)
	if(visit[i]&&((s0&bitmap[mapptr[i][j]])>0)){
	  if(i==0){
	    delete []visit;
	    return false;
	  }
	  visit[i]=false;
	  nvisit--;
	}
    if(nvisit==1)break;
  }
  multi=0;
  for(i=0;i<ntransc;i++)
    if(visit[i]){
      tx=transtep[0][i];
      ty=transtep[1][i];
      tpx=transtep[2][i];
      tpy=transtep[3][i];
      tsx=transtep[4][i];
      tsy=transtep[5][i];
      tr=transtep[6][i];
      tmp=real(f4func[ly+ty][lx+tx][tr]);
      sgn=0;
      if(tsx&&sgm1==-1)sgn++;
      if(tsy&&sgm2==-1)sgn++;
      if(px==-1&&tpx==1)sgn++;
      if(py==-1&&tpy==1)sgn++;
      if((sgn%2)==0) multi+=tmp;
      else multi-=tmp;
    }
  multi*=perix*periy/(rx*ry);
  delete []visit;
  return true;
}
/*
//--------------------------------------------------------------------------------------
bool fci::checkrepresentative1(unsigned long s0, double& multi){
//--------------------------------------------------------------------------------------
  int i,j,k,tx,ty,tz,rx,ry,rpx,rpy,rs1,rs2,tpx,tpy,ts1,ts2,facz,facs1,facs2;
  unsigned long st,stx,sty,spx,spy,stz,sts1;
  double fac;
  bool pass;
  pass=false;
  //calculate periodicity rx
  rx=lx;
  for(i=2;i<lx;i+=2){
    translatex(s0,i,st);
    if(st<s0)
      return pass;
    else if(st==s0){
      if(((kx*i)%perix)!=0)
	return pass;
      rx=i;
      break;      
    }
  }
  //calculate periodicity ry
  ry=ly;
  for(j=2;j<ly;j+=2){
    translatey(s0,j,st);
    if(st<s0)
      return pass;
    else if(st==s0){
      if(((ky*j)%periy)!=0)return pass;
      ry=j;
      break;
    }
  }
  if((2*kx==0||2*kx==(perix/2))&&(2*ky==0||2*ky==(periy/2))&&kx==ky&&kx==0){
    if(((px*ry)%periy)!=0)return pass;
    if(((py*rx)%perix)!=0)return pass;
    multi=0;
    for(tpx=0;tpx<2;tpx++)
      for(tpy=0;tpy<2;tpy++)
	for(tx=0;tx<rx;tx+=2)
	  for(ty=0;ty<ry;ty+=2)
	    for(ts1=0;ts1<2;ts1++)
	      for(ts2=0;ts2<2;ts2++){
		st=latticetransformations(s0,tx,ty,tpx,tpy,ts1,ts2,0);
		if(st<s0)return pass;
		else if(st==s0)
		  multi+=pow(sgm1,ts1)*pow(sgm2,ts2)*real(f4func[ly+ty][lx+tx][lx+tpy][ly+tpx]);
	      }
    //multi/=(rx*ry);
    multi*=(perix*periy)/(rx*ry);
    pass=true;
    return pass;
  }
}
//--------------------------------------------------------------------------------------
unsigned long fci::latticetransformations(unsigned long s0, int tx, int ty, int tpx, int tpy, int tsx, int tsy, int tz){
//--------------------------------------------------------------------------------------
  int i,*ptr;
  unsigned long sx;
  
  if(tx==0&&ty==0&&tpx==0&&tpy==0&&tsx==0&&tsy==0){
    if(tz==0)
      return s0;
    else
      return s0^bitallone;
  }
  ptr=latticemaps[tx][ty][tpx][tpy][tsx][tsy];
  sx=0;
  if(tz==0)
    for(i=0;i<ns;i++){
      if((s0&bitmap[ptr[i]])>0)
	sx|=((unsigned long)1<<i);
    }
  else
    for(i=0;i<ns;i++){
      if((s0&bitmap[ptr[i]])==0)
	sx|=((unsigned long)1<<i);
    }
  return sx;
}
*/

//--------------------------------------------------------------------------------------
void fci::translatex(unsigned long s0, int i0, unsigned long &sx){
//--------------------------------------------------------------------------------------
  int i;
  if(i0==0){
    sx=s0;
    return;
  }
  sx=0;
  for(i=0;i<ns;i++){
    if(s0%2==1)
      sx=sx|((unsigned long)1<<txmap[i0*ns+i]);
    s0/=2;
  }
}

//--------------------------------------------------------------------------------------
void fci::translatey(unsigned long s0, int i0, unsigned long &sy){
//--------------------------------------------------------------------------------------
  int i;
  if(i0==0){
    sy=s0;
    return;
  }
  sy=0;
  for(i=0;i<ns;i++){
    if(s0%2==1)
      sy=sy|((unsigned long)1<<tymap[i0*ns+i]);
    s0/=2;
  }
}

//--------------------------------------------------------------------------------------
void fci::rotaterepr(unsigned long s0, unsigned long &sx){
//--------------------------------------------------------------------------------------
  int i;
  sx=0;
  for(i=0;i<ns;i++){
    if(s0%2==1)
      sx=sx|((unsigned long)1<<rotmap[ns+i]);
    s0/=2;
  }
}

//--------------------------------------------------------------------------------------
bool fci::findrxry(unsigned long s0, int& rx, int& ry, bool flag){
//--------------------------------------------------------------------------------------
  int i,j,rx1,ry1;
  unsigned long s1,s2,r,st;
  rx=perix;
  ry=periy;
  r=s0;
  for(i=1;i<periy;i++){
    s2=r&row1;
    s2<<=(ns-lx);
    s1=(r>>lx);
    r=s1+s2;
    if(r==s0){
      if(((ky*i)%periy)!=0)return false;
      ry=i;
      break;
    }
    else if(flag&&r<s0)return false;
  }
  r=s0;
  for(i=1;i<perix;i++){
    s2=r&col1;
    s2<<=(lx-1);
    r&=col1rev;
    s1=(r>>1);
    r=s1+s2;
    if(r==s0){
      if(((kx*i)%perix)!=0)return false;
      rx=i;
      break;
    }
    else if(flag&&r<s0)return false;
  }

  //calculate periodicity rx
  rx1=lx;
  for(i=1;i<lx;i++){
    translatex(s0,i,st);
    //if(st<s0)
    //return pass;
    //else if(st==s0){
    //if(((kx*i)%perix)!=0)
    //return pass;
    if(st==s0){
      rx1=i;
      break;      
    }
  }
  //calculate periodicity ry
  ry1=ly;
  for(j=1;j<ly;j++){
    translatey(s0,j,st);
    //if(st<s0)
    //return pass;
    //else if(st==s0){
    //if(((ky*j)%periy)!=0)return pass;
    if(st==s0){
      ry1=j;
      break;
    }
  }
  if(rx!=rx1||ry!=ry1){
    cout<<s0<<"\twrong find rx ry\trx="<<rx<<"\t"<<rx1<<"\try="<<ry<<"\t"<<ry1<<endl;
    exit(0);
  }
  return true;

}

//--------------------------------------------------------------------------------------
bool fci::checkrepresentative2(unsigned long s0, double &multi){
//--------------------------------------------------------------------------------------
  int i,j,k,l,m,n,nvisit,tx,ty,tpx,tpy,tsx,tsy,sgn,rx,ry,tr;
  bool *visit,zeroflag;
  double tmp;
  unsigned long s1;
  if(findrxry(s0,rx,ry,true)==false)return false;
  visit=new bool[ntransc*2];
  nvisit=0;
  for(i=0;i<ntransc;i++){
    if(transtep[0][i]<rx&&transtep[1][i]<ry){    
      visit[i]=true;
      visit[i+ntransc]=true;
      nvisit+=2;
    }
    else{
      visit[i]=false;
      visit[i+ntransc]=false;
    }
  }
  for(j=ns-1;j>=0;j--){
    zeroflag=false;
    for(i=0;i<ntransc;i++)
      if(visit[i]||visit[i+ntransc]){
	s1=s0&bitmap[mapptr[i][j]];
	if(visit[i]&&s1==0){
	  zeroflag=true;
	  break;
	}
	else if(visit[i+ntransc]&&s1>0){
	  zeroflag=true;
	  break;
	}
      }
    if(zeroflag&&nvisit>1)
      for(i=0;i<ntransc;i++)
	if(visit[i]||visit[i+ntransc]){
	  s1=s0&bitmap[mapptr[i][j]];
	  if(s1>0&&visit[i]){
	    if(i==0){
	      delete []visit;
	      return false;
	    }
	    visit[i]=false;
	    nvisit--;
	  }
	  else if(s1==0&&visit[i+ntransc]){
	    visit[i+ntransc]=false;
	    nvisit--;
	  }
	}
    if(nvisit==1) break;
  }
  multi=0;
  for(i=0;i<ntransc;i++)
    if(visit[i]||visit[i+ntransc]){
      tx=transtep[0][i];
      ty=transtep[1][i];
      tpx=transtep[2][i];
      tpy=transtep[3][i];
      tsx=transtep[4][i];
      tsy=transtep[5][i];
      tr=transtep[6][i];
      tmp=real(f4func[ly+ty][lx+tx][tr]);
      sgn=0;
      if(tsx&&sgm1==-1)sgn++;
      if(tsy&&sgm2==-1)sgn++;
      if(px==-1&&tpx==1)sgn++;
      if(py==-1&&tpy==1)sgn++;
      if(visit[i+ntransc]&&(qz==-1))sgn++;
      if((sgn%2)==0) multi+=tmp;
      else multi-=tmp;
    }
  multi*=perix*periy/(rx*ry);
  //cout<<"s0="<<s0<<"\tmulti="<<multi<<endl;
  delete []visit;
  return true;
}

//--------------------------------------------------------------------------------------
void fci::findrepresentative(unsigned long s0, unsigned long &r, int &sx, int &sy ,int &lpx, int &lpy,int &rz,int &ls1, int &ls2,int& tr){
//--------------------------------------------------------------------------------------
  if(sec!=0){
    rz=0;
    findrepresentative1(s0,r,sx,sy,lpx,lpy,ls1,ls2,tr);
  }
  else if(sec==0)
    findrepresentative2(s0,r,sx,sy,lpx,lpy,ls1,ls2,rz,tr);
}
///*
//--------------------------------------------------------------------------------------
void fci::findrepresentative1(unsigned long s0, unsigned long& r, int& tx, int& ty, int& tpx, int& tpy, int& tsx, int& tsy, int& tr){
//--------------------------------------------------------------------------------------
  int i,j,k,l,m,n,nvisit;
  bool *visit,zeroflag;
  visit=new bool[ntrans];
  nvisit=ntrans;
  for(i=0;i<ntrans;i++)visit[i]=true;
  for(j=ns-1;j>=0;j--){
    zeroflag=false;
    for(i=0;i<ntrans;i++)
      if(visit[i]&&((s0&bitmap[mapptr[i][j]])==0)){
	zeroflag=true;
	break;
      }
    if(zeroflag&&nvisit>1)
      for(i=0;i<ntrans;i++)
	if(visit[i]&&((s0&bitmap[mapptr[i][j]])>0)){
	  visit[i]=false;
	  nvisit--;
	}
    if(nvisit==1)break;
  }
  for(i=0;i<ntrans;i++)
    if(visit[i]){
      tx=transtep[0][i];
      ty=transtep[1][i];
      tpx=transtep[2][i];
      tpy=transtep[3][i];
      tsx=transtep[4][i];
      tsy=transtep[5][i];
      tr=transtep[6][i];
      if(i==0){
	delete []visit;
	r=s0;
	return;
      }
      r=0;
      for(j=0;j<ns;j++)
	if((s0&bitmap[mapptr[i][j]])>0)
	  r|=bitmap[j];
      delete []visit;
      return;
    }
}
/*
//--------------------------------------------------------------------------------------
void fci::findrepresentative1(unsigned long s0, unsigned long &r, int &sx, int &sy ,int &lpx, int &lpy, int &ls1, int &ls2){
//--------------------------------------------------------------------------------------
  unsigned long st,stx,sty,spx,spy,stz,sts1;
  int tx,ty,tpx,tpy,tz,ts1,ts2;
  sx=0;
  sy=0;
  lpx=0;
  lpy=0;
  ls1=0;
  ls2=0;
  r=s0;
  if((2*kx==0||2*kx==(perix/2))&&(2*ky==0||2*ky==(periy/2))&&kx==ky&&kx==0){
    for(tpx=0;tpx<2;tpx++)
      for(tpy=0;tpy<2;tpy++)
	for(tx=0;tx<perix;tx+=2)
	  for(ty=0;ty<periy;ty+=2)
	    for(ts1=0;ts1<2;ts1++)
	      for(ts2=0;ts2<2;ts2++){
		st=latticetransformations(s0,tx,ty,tpx,tpy,ts1,ts2,0);
		if(st<r){
		  r=st;
		  sx=tx;
		  sy=ty;
		  lpx=tpx;
		  lpy=tpy;
		  ls1=ts1;
		  ls2=ts2;
		}
	      }
    return;    
  }
}
*/
//--------------------------------------------------------------------------------------
void fci::findrepresentative2(unsigned long s0, unsigned long& r, int& tx, int& ty, int& tpx, int& tpy, int& tsx, int& tsy, int& tz,int& tr){
//--------------------------------------------------------------------------------------
  int i,j,k,l,m,n,nvisit,rx,ry;
  bool *visit,zeroflag;
  unsigned long s1;
  /*
  if(findrxry(s0,rx,ry,false)==false){
    r=nrep;
    return;
  }
  nvisit=0;
  visit=new bool[ntrans*2];
  for(i=0;i<ntrans;i++){
    if(transtep[0][i]<rx&&transtep[1][i]<ry){    
      visit[i]=true;
      visit[i+ntrans]=true;
      nvisit+=2;
    }
    else{
      visit[i]=false;
      visit[i+ntrans]=false;
    }
  }
  */
  nvisit=2*ntrans;
  visit=new bool[ntrans*2];
  for(i=0;i<2*ntrans;i++)visit[i]=true;
  for(j=ns-1;j>=0;j--){
    zeroflag=false;
    for(i=0;i<ntrans;i++)
      if(visit[i]||visit[i+ntrans]){
	s1=s0&bitmap[mapptr[i][j]];
	if(visit[i]&&s1==0){
	  zeroflag=true;
	  break;
	}
	else if(visit[i+ntrans]&&s1>0){
	  zeroflag=true;
	  break;
	}
      }
    if(zeroflag&&nvisit>1)
      for(i=0;i<ntrans;i++)
	if(visit[i]||visit[i+ntrans]){
	  s1=s0&bitmap[mapptr[i][j]];
	  if(visit[i]&&s1>0){
	    visit[i]=false;
	    nvisit--;
	  }
	  else if(visit[i+ntrans]&&s1==0){
	    visit[i+ntrans]=false;
	    nvisit--;
	  }
	}
    if(nvisit==1) break;
  }
  for(i=0;i<ntrans;i++)
    if(visit[i]||visit[i+ntrans]){
      tx=transtep[0][i];
      ty=transtep[1][i];
      tpx=transtep[2][i];
      tpy=transtep[3][i];
      tsx=transtep[4][i];
      tsy=transtep[5][i];
      tr=transtep[6][i];
      if(i==0){
	if(visit[i]){
	  tz=0;
	  r=s0;
	}
	else{
	  tz=1;
	  r=s0^bitallone;
	}
	delete []visit;
	return;
      }
      r=0;
      for(j=0;j<ns;j++)
	if((s0&bitmap[mapptr[i][j]])>0)
	  r|=bitmap[j];
      if(visit[i])
	tz=0;
      else{
	tz=1;
	r^=bitallone;
      }
      delete []visit;
      return;
    }
}

//--------------------------------------------------------------------------------------
bool fci::build_latticemaps(){
//--------------------------------------------------------------------------------------
  int i,j,k,l,m,n,q,p,tx,ty,tpx,tpy,tsx,tsy,tr;
  dcmplex phase1,phase2;
  ntrans0=64*perix*periy;
  mapptr=new int*[ntrans0];
  transtep=new int*[7];
  transtep[0]=new int[ntrans0*7];
  for(i=1;i<7;i++)
    transtep[i]=&(transtep[0][ntrans0*i]);
  for(i=0;i<64;i++)
    bitmap[i]=(1UL<<i);
  bitallone=0;
  for(i=0;i<ns;i++)
    bitallone|=(1UL<<i);
  row1=0;
  for(i=0;i<lx;i++)
    row1+=bitmap[i];
  col1=0;
  for(i=0;i<ly;i++)
    col1+=bitmap[i*lx];
  col1rev=col1^bitallone;
  cout<<"bitallone="<<bitallone<<endl;
  latticemaps=new int*******[perix];
  for(i=0;i<perix;i+=1){
    latticemaps[i]=new int******[periy];
    for(j=0;j<periy;j+=1){
      latticemaps[i][j]=new int*****[2];
      for(k=0;k<2;k++){
	latticemaps[i][j][k]=new int****[2];
	for(l=0;l<2;l++){
	  latticemaps[i][j][k][l]=new int***[2];
	  for(m=0;m<2;m++){
	    latticemaps[i][j][k][l][m]=new int**[2];
	    for(n=0;n<2;n++){
	      latticemaps[i][j][k][l][m][n]=new int*[4];
	      for(p=0;p<4;p++){
		latticemaps[i][j][k][l][m][n][p]=new int[ns*2];
	      }
	    }
	  }
	}
      }
    }
  }
  for(tx=0;tx<perix;tx+=1)
    for(ty=0;ty<periy;ty+=1)
      for(tpx=0;tpx<2;tpx++)
	for(tpy=0;tpy<2;tpy++)
	  for(tsx=0;tsx<2;tsx++)
	    for(tsy=0;tsy<2;tsy++)
	      for(tr=0;tr<4;tr++)
		for(i=0;i<ns;i++){
		  k=sigma2map[i+tsy*ns];
		  j=sigma1map[k+tsx*ns];
		  k=tymap[j+ty*ns];
		  j=txmap[k+tx*ns];
		  k=pymap[j+tpy*ns];
		  j=pxmap[k+tpx*ns];
		  k=rotmap[j+tr*ns];
		  latticemaps[tx][ty][tpx][tpy][tsx][tsy][tr][k]=i;
		  latticemaps[tx][ty][tpx][tpy][tsx][tsy][tr][i+ns]=k;
		}
  cout<<"done latticemaps"<<endl;
  ntrans=0;
  if(kx!=0&&kx!=(perix/2)&&ky!=0&&ky!=(periy/2)&&kx!=ky&&(kx+ky)%perix==0){//6*6 case 2kx=-2ky, 2kx=2pi/3,4pi/3
    for(i=0;i<perix;i+=1)
      for(j=0;j<periy;j+=1)
	for(n=0;n<2;n++){
	  transtep[0][ntrans]=i;
	  transtep[1][ntrans]=j;
	  transtep[2][ntrans]=0;
	  transtep[3][ntrans]=0;
	  transtep[4][ntrans]=0;
	  transtep[5][ntrans]=n;
	  transtep[6][ntrans]=0;
	  mapptr[ntrans]=latticemaps[i][j][0][0][0][n][0];
	  ntrans++;
	}
  }
  else if(kx!=0&&kx!=(perix/2)&&ky!=0&&ky!=(periy/2)&&kx==ky){//6*6 case, kx!=0,pi, kx==ky
    for(i=0;i<perix;i+=1)
      for(j=0;j<periy;j+=1)
	for(m=0;m<2;m++){
	  transtep[0][ntrans]=i;
	  transtep[1][ntrans]=j;
	  transtep[2][ntrans]=0;
	  transtep[3][ntrans]=0;
	  transtep[4][ntrans]=m;
	  transtep[5][ntrans]=0;
	  transtep[6][ntrans]=0;
	  mapptr[ntrans]=latticemaps[i][j][0][0][m][0][0];
	  ntrans++;
	}
  }
  else if(kx!=0&&kx!=(perix/2)&&ky!=0&&ky!=(periy/2)&&kx!=ky&&(kx+ky)%perix!=0){//6*6 case, kx!=0,pi, kx==ky
    for(i=0;i<perix;i+=1)
      for(j=0;j<periy;j+=1){
	transtep[0][ntrans]=i;
	transtep[1][ntrans]=j;
	transtep[2][ntrans]=0;
	transtep[3][ntrans]=0;
	transtep[4][ntrans]=0;
	transtep[5][ntrans]=0;
	transtep[6][ntrans]=0;
	mapptr[ntrans]=latticemaps[i][j][0][0][0][0][0];
	ntrans++;
      }
  }
  else if((kx==0||kx==(perix/2))&&ky!=0&&ky!=(periy/2)){//use TxTyGx
    for(i=0;i<perix;i+=1)
      for(j=0;j<periy;j+=1)
	for(k=0;k<2;k++){
	  transtep[0][ntrans]=i;
	  transtep[1][ntrans]=j;
	  transtep[2][ntrans]=k;
	  transtep[3][ntrans]=0;
	  transtep[4][ntrans]=0;
	  transtep[5][ntrans]=0;
	  transtep[6][ntrans]=0;
	  mapptr[ntrans]=latticemaps[i][j][k][0][0][0][0];
	  ntrans++;
	}    
  }
  else if((ky==0||ky==(periy/2))&&kx!=0&&kx!=(perix/2)){//use TxTyGy
    for(i=0;i<perix;i+=1)
      for(j=0;j<periy;j+=1)
	for(l=0;l<2;l++){
	  transtep[0][ntrans]=i;
	  transtep[1][ntrans]=j;
	  transtep[2][ntrans]=0;
	  transtep[3][ntrans]=l;
	  transtep[4][ntrans]=0;
	  transtep[5][ntrans]=0;
	  transtep[6][ntrans]=0;
	  mapptr[ntrans]=latticemaps[i][j][0][l][0][0][0];
	  ntrans++;
	}    
  }
  else if(((ky==0||ky==(periy/2))||(kx==0||kx==(perix/2)))&&kx!=ky){
    for(i=0;i<perix;i+=1)
      for(j=0;j<periy;j+=1)
	for(k=0;k<2;k++)
	  for(l=0;l<2;l++){
	    transtep[0][ntrans]=i;
	    transtep[1][ntrans]=j;
	    transtep[2][ntrans]=k;
	    transtep[3][ntrans]=l;
	    transtep[4][ntrans]=0;
	    transtep[5][ntrans]=0;
	    transtep[6][ntrans]=0;
	    mapptr[ntrans]=latticemaps[i][j][k][l][0][0][0];
	    ntrans++;
	  }    
  }
  else if(kx==ky&&(kx==0||kx==perix/2)&&px==py){//4*4 case, pi,pi, use TxTySigma1
    for(n=0;n<2;n++)
      for(m=0;m<2;m++)
	//for(p=0;p<4;p++)
	  for(l=0;l<2;l++)
	    for(k=0;k<2;k++)
	      for(j=0;j<periy;j+=1)
		for(i=0;i<perix;i+=1){
		  transtep[0][ntrans]=i;
		  transtep[1][ntrans]=j;
		  transtep[2][ntrans]=k;
		  transtep[3][ntrans]=l;
		  transtep[4][ntrans]=m;
		  transtep[5][ntrans]=n;
		  //transtep[6][ntrans]=p;
		  transtep[6][ntrans]=0;
		  //mapptr[ntrans]=latticemaps[i][j][k][l][m][n][p];
		  mapptr[ntrans]=latticemaps[i][j][k][l][m][n][0];
		  ntrans++;
		}
  }
  else if(kx==ky&&(kx==0||kx==perix/2)&&px!=py){//4*4 case, pi,pi, use TxTySigma1
    //for(n=0;n<2;n++)
    //for(m=0;m<2;m++)
	for(l=0;l<2;l++)
	  for(k=0;k<2;k++)
	    for(j=0;j<periy;j+=1)
	      for(i=0;i<perix;i+=1){
		transtep[0][ntrans]=i;
		transtep[1][ntrans]=j;
		transtep[2][ntrans]=k;
		transtep[3][ntrans]=l;
		//transtep[4][ntrans]=m;
		//transtep[5][ntrans]=n;
		transtep[4][ntrans]=0;
		transtep[5][ntrans]=0;
		transtep[6][ntrans]=0;
		//mapptr[ntrans]=latticemaps[i][j][k][l][m][n][0];
		mapptr[ntrans]=latticemaps[i][j][k][l][0][0][0];
		ntrans++;
	      }
  }
  else if(0&&kx==ky&&kx==(perix/2)){//4*4 case, pi,pi, use TxTySigma1
    for(n=0;n<2;n++)
      for(m=0;m<2;m++)
	//for(l=0;l<2;l++)
	//for(k=0;k<2;k++)
	    for(j=0;j<periy;j+=1)
	      for(i=0;i<perix;i+=1){
		transtep[0][ntrans]=i;
		transtep[1][ntrans]=j;
		transtep[2][ntrans]=0;
		transtep[3][ntrans]=0;
		transtep[4][ntrans]=m;
		transtep[5][ntrans]=n;
		transtep[6][ntrans]=0;
		//mapptr[ntrans]=latticemaps[i][j][k][l][m][n][0];
		mapptr[ntrans]=latticemaps[i][j][0][0][m][n][0];
		ntrans++;
	      }
  }
  ntransc=ntrans;
  //if(kx==ky&&kx==0)
  //ntrans/=2;
  cout<<"kx="<<kx<<"\tky="<<ky<<"\tntrans="<<ntrans<<endl;
  bool *flag,check;
  flag=new bool[ntrans];
  for(i=0;i<ntrans;i++)flag[i]=true;
  n=0;
  for(i=0;i<ntrans;i++)
    if(flag[i]){
      tx=transtep[0][i];
      ty=transtep[1][i];
      tpx=transtep[2][i];
      tpy=transtep[3][i];
      tsx=transtep[4][i];
      tsy=transtep[5][i];
      tr=transtep[6][i];
      phase1=f4func[ly+ty][lx+tx][tr];
      if(tsx==1&&sgm1==-1)phase1*=-1;
      if(tsy==1&&sgm2==-1)phase1*=-1;
      if(px==-1&&tpx==1)phase1*=-1;
      if(py==-1&&tpy==1)phase1*=-1;
      for(j=i+1;j<ntrans;j++)
	if(flag[j]){
	  check=true;
	  for(k=0;k<ns;k++)
	    if(mapptr[i][k]!=mapptr[j][k]){
	      check=false;
	      break;
	    }
	  if(check){
	    tx=transtep[0][j];
	    ty=transtep[1][j];
	    tpx=transtep[2][j];
	    tpy=transtep[3][j];
	    tsx=transtep[4][j];
	    tsy=transtep[5][j];
	    tr=transtep[6][j];
	    phase2=f4func[ly+ty][lx+tx][tr];
	    if(tsx==1&&sgm1==-1)phase2*=-1;
	    if(tsy==1&&sgm2==-1)phase2*=-1;
	    if(px==-1&&tpx==1)phase2*=-1;
	    if(py==-1&&tpy==1)phase2*=-1;
	    if(norm(phase1-phase2)>1.e-12){
	      cout<<"the quantum number is not consistent"<<endl;
	      cout<<"phase1="<<phase1<<"\tphase2="<<phase2<<endl;
	      delete []flag;
	      return false;
	    }
	    flag[j]=false;
	    n++;
	  }
	}
    }
  k=0;
  for(i=0;i<ntrans;i++){
    if(flag[i]&&k<i){
      mapptr[k]=mapptr[i];
      for(j=0;j<7;j++)
	transtep[j][k]=transtep[j][i];
    }
    if(flag[i])k++;
  }
  if(k!=ntrans-n){
    cout<<"wrong counting ntrans"<<endl;
    exit(0);
  }
  ntrans=k;
  delete []flag;
  cout<<"kx="<<kx<<"\tky="<<ky<<"\tntrans="<<ntrans<<endl;
  tranbit=new unsigned long[ntrans];
  ntransite=new int[ntrans];
  for(i=0;i<ntrans;i++){
    ntransite[i]=0;
    tranbit[i]=0;
    for(j=0;j<ns;j++)
      if(mapptr[i][j]==j){
	ntransite[i]++;
	tranbit[i]|=bitmap[j];
      }
    ntransite[i]=ns-ntransite[i];
  }
  transite=new int*[ntrans];
  for(i=0;i<ntrans;i++){
    transite[i]=new int[ntransite[i]];
    k=0;
    for(j=0;j<ns;j++)
      if(mapptr[i][j]!=j){
	transite[i][k]=j;
	k++;
      }
  }
  //cout<<"output number of transition sites"<<endl;
  //for(i=0;i<ntrans;i++)
  //cout<<i<<"\t"<<ntransite[i]<<endl;
  return true;
}

//--------------------------------------------------------------------------------------
void fci::check_latticemaps(){
//--------------------------------------------------------------------------------------
  int i,j,k,l,m,n,q,tx,ty,tpx,tpy,tsx,tsy,tr,t1,t2,t3,t4,t5,t6,t7,rxry;
  bool check;
  rxry=perix*periy;
  l=0;
  tr=0;
  t7=0;
  for(tpx=0;tpx<2;tpx++)
    for(tpy=0;tpy<2;tpy++)
      for(tsx=0;tsx<2;tsx++)
	for(tsy=0;tsy<2;tsy++)
	  //for(tr=0;tr<4;tr++)
	    for(tx=0;tx<perix;tx+=1)
	      for(ty=0;ty<periy;ty+=1)
		for(t3=0;t3<2;t3++)
		  for(t4=0;t4<2;t4++)
		    for(t5=0;t5<2;t5++)
		      for(t6=0;t6<2;t6++)
			//for(t7=0;t7<4;t7++)
			  for(t1=0;t1<perix;t1+=1)
			    for(t2=0;t2<periy;t2+=1){
			      check=true;
			      for(i=0;i<ns;i++)
				if(latticemaps[tx][ty][tpx][tpy][tsx][tsy][tr][i]!=latticemaps[t1][t2][t3][t4][t5][t6][t7][i]){
				  check=false;
				  break;
				}
			      if(check&&!(tx==t1&&ty==t2&&tpx==t3&&tpy==t4&&tsx==t5&&tsy==t6&&tr==t7)){
				j=tx+ty*perix+tpx*rxry+tpy*2*rxry+tsx*4*rxry+tsy*8*rxry+tr*16*rxry;
				k=t1+t2*perix+t3*rxry+t4*2*rxry+t5*4*rxry+t6*8*rxry+t7*16*rxry;
				if(j<k){
				  l++;
				  cout<<j<<"\tfind identical1 \t"<<tx<<"\t"<<ty<<"\t"<<tpx<<"\t"<<tpy<<"\t"<<tsx<<"\t"<<tsy<<"\t"<<tr<<endl;
				  cout<<k<<"\tfind identical2 \t"<<t1<<"\t"<<t2<<"\t"<<t3<<"\t"<<t4<<"\t"<<t5<<"\t"<<t6<<"\t"<<t7<<endl;
				  cout<<endl;
				}
			      }
			    }
  cout<<"number of duplicate maps="<<l<<endl;
  exit(0);
}

//--------------------------------------------------------------------------------------
void fci::compute_specific_heat(double *ee, double *e2, double *nor, double* mz, double* m2, double tau, int ntau, double enr_gs ){
//--------------------------------------------------------------------------------------
  int i,j;
  double tmp,beta;
  for(j=0;j<ntau;j++){
    beta=1./((double)(j+1)*tau);
    ee[j]=0;
    e2[j]=0;
    mz[j]=0;
    m2[j]=0;
    nor[j]=0;
    for(i=0;i<nrep;i++){
      tmp=exp(-beta*(eigval[i]-enr_gs));
      mz[j]+=(double)sec*tmp;
      m2[j]+=(double)(sec*sec)*tmp;
      ee[j]+=eigval[i]*tmp;
      e2[j]+=eigval[i]*eigval[i]*tmp;
      nor[j]+=tmp;
    }
  }
}

//--------------------------------------------------------------------------------------
dcmplex inner_prod(unsigned long nrep, dcmplex* v1, dcmplex* v2){
//--------------------------------------------------------------------------------------
  unsigned long i;
  double prod;
  dcmplex pd;
  prod=0;
  //#pragma omp parallel for  reduction(+ : prod)
  for(i=0;i<nrep;i++)
    prod+=real(conj(v1[i])*v2[i]);
  pd=prod;
  return pd;
}

//--------------------------------------------------------------------------------------
double inner_prod(unsigned long nrep, double* v1, double* v2){
//--------------------------------------------------------------------------------------
  unsigned long i;
  double prod;
  prod=0;
  //#pragma omp parallel for reduction(+ : prod)
  for(i=0;i<nrep;i++)
    prod+=v1[i]*v2[i];
  return prod;
}

//--------------------------------------------------------------------------------------
void fci::initialize_lanczos(){
//--------------------------------------------------------------------------------------
  unsigned long i;
  double nor;
  if(kx%(perix/2)==0&&ky%(periy/2)==0&&qr%2==0){ //use real hamiltonian matrix
    nor=0;
    for(i=0;i<nrep;i++){
      ffr[0][i]=ran_()*2.-1.;
      //ffr[0][i]=1;
      nor+=ffr[0][i]*ffr[0][i];
    }
    nor=sqrt(nor);
    //#pragma omp parallel for
    for(i=0;i<nrep;i++)
      ffr[0][i]/=nor;
  }
  else{
    nor=0;
    for(i=0;i<nrep;i++){
      //real(ff[0][i])=ran_()*2.-1.;
      //imag(ff[0][i])=ran_()*2.-1.;
      ff[0][i]=std::complex<double>(1,1);
      nor+=norm(ff[0][i]);
    }
    nor=sqrt(nor);
    //#pragma omp parallel for
    for(i=0;i<nrep;i++)
      ff[0][i]/=nor;
  }
}

//--------------------------------------------------------------------------------------
int fci::lanczos1(){
//--------------------------------------------------------------------------------------
  unsigned long i;
  int m,j;
  dcmplex q1,q2;
  double eval[200],t1,t2,t3;
  bool check;
  cout<<"start lanczos1"<<endl;
  for(i=0;i<mlanc;i++){
    eigval[i]=0;
    eval[i]=0;
  }
  for(i=0;i<mlanc;i++)
    aal[i]=0;
  for(i=0;i<mlanc;i++)
    nnl[i]=0;
  nnl[0]=1;
  //t1=omp_get_wtime();  
  hamiltonian_vector_multiplication(ff[0],ff[1]);
  //t2=omp_get_wtime();  
  cout<<"Hv time="<<t2-t1<<endl;
  aal[0]=real(inner_prod(nrep,ff[0],ff[1]));
  //#pragma omp parallel for 
  for(i=0;i<nrep;i++)
    ff[1][i]=ff[1][i]-aal[0]*ff[0][i];
  nnl[1]=sqrt(real(inner_prod(nrep,ff[1],ff[1])));
  //#pragma omp parallel for 
  for(i=0;i<nrep;i++)
    ff[1][i]/=nnl[1];
  //t3=omp_get_wtime();  
  cout<<"m="<<1<<"\tvec process time\t"<<t3-t2<<endl;

  if(nnl[1]<1.e-8){
    diatridiag(2);
    return 2;
  }

  for(m=2;m<mlanc;m++){
    //t1=omp_get_wtime();  
    hamiltonian_vector_multiplication(ff[m-1],ff[m]);
    //t2=omp_get_wtime();  
    cout<<"m="<<m<<"\tHv time="<<t2-t1<<endl;
    aal[m-1]=real(inner_prod(nrep,ff[m-1],ff[m]));
    //#pragma omp parallel for 
    for(i=0;i<nrep;i++)
      ff[m][i]=ff[m][i]-aal[m-1]*ff[m-1][i]-nnl[m-1]*ff[m-2][i];
    nnl[m]=sqrt(real(inner_prod(nrep,ff[m],ff[m])));
    //#pragma omp parallel for 
    for(i=0;i<nrep;i++)
      ff[m][i]/=nnl[m];

    for(j=0;j<m;j++){
      q1=inner_prod(nrep,ff[j],ff[m]);
      if(norm(q1)<1.e-20)continue;
      q2=1./sqrt(1.-norm(q1));
      //#pragma omp parallel for 
      for(i=0;i<nrep;i++)
	ff[m][i]=q2*(ff[m][i]-q1*ff[j][i]);
    }
    //t3=omp_get_wtime();  
    cout<<"m="<<m<<"\tvec process time\t"<<t3-t2<<endl;
    if(nnl[m]<1.e-8){
      diatridiag(m);
      return m;
    }
    if((m%5)==0){
      diatridiag(m);
      check=true;
      if(m==5){
	check=false;
	for(j=0;j<neig;j++)
	  eval[j]=eigval[j];
      }
      else{
	for(j=0;j<2;j++){
	  if(fabs(eigval[j]-eval[j])>1.e-8)check=false;
	  eval[j]=eigval[j];
	}
      }
      if(check&&compute_eigenvector(m))return m;
    }
  }
  diatridiag(mlanc-1);
  check=compute_eigenvector(mlanc-1);
  return mlanc-1;
}

//--------------------------------------------------------------------------------------
int fci::lanczos1real(){
//--------------------------------------------------------------------------------------
  unsigned long i;
  int m,j;
  double q1,q2,t1,t2,t3;
  double eval[200];
  bool check;
  cout<<"start lanczos1real"<<endl;
  for(i=0;i<mlanc;i++){
    eigval[i]=0;
    eval[i]=0;
  }
  for(i=0;i<mlanc;i++)
    aal[i]=0;
  for(i=0;i<mlanc;i++)
    nnl[i]=0;
  nnl[0]=1;
  //t1=omp_get_wtime();  
  hamiltonian_vector_multiplication(ffr[0],ffr[1]);
  //t2=omp_get_wtime();  
  cout<<"Hv time="<<t2-t1<<endl;
  aal[0]=inner_prod(nrep,ffr[0],ffr[1]);
  //#pragma omp parallel for 
  for(i=0;i<nrep;i++)
    ffr[1][i]=ffr[1][i]-aal[0]*ffr[0][i];
  nnl[1]=sqrt(inner_prod(nrep,ffr[1],ffr[1]));
  //#pragma omp parallel for 
  for(i=0;i<nrep;i++)
    ffr[1][i]/=nnl[1];
  //t3=omp_get_wtime();  
  cout<<"m="<<1<<"\tvec process time\t"<<t3-t2<<endl;

  if(nnl[1]<1.e-8){
    diatridiag(2);
    return 2;
  }

  for(m=2;m<mlanc;m++){
    //t1=omp_get_wtime();  
    hamiltonian_vector_multiplication(ffr[m-1],ffr[m]);
    //t2=omp_get_wtime();  
    cout<<"m="<<m<<"\tHv time="<<t2-t1<<endl;
    aal[m-1]=inner_prod(nrep,ffr[m-1],ffr[m]);
    //#pragma omp parallel for 
    for(i=0;i<nrep;i++)
      ffr[m][i]=ffr[m][i]-aal[m-1]*ffr[m-1][i]-nnl[m-1]*ffr[m-2][i];
    nnl[m]=sqrt(inner_prod(nrep,ffr[m],ffr[m]));
    //#pragma omp parallel for 
    for(i=0;i<nrep;i++)
      ffr[m][i]/=nnl[m];
    for(j=0;j<m;j++){
      q1=inner_prod(nrep,ffr[j],ffr[m]);
      if(q1*q1<1.e-20)continue;
      q2=1./sqrt(1.-q1*q1);
      //#pragma omp parallel for 
      for(i=0;i<nrep;i++)
	ffr[m][i]=q2*(ffr[m][i]-q1*ffr[j][i]);
    }
    //t3=omp_get_wtime();  
    cout<<"m="<<m<<"\tvec process time\t"<<t3-t2<<endl;
    if(nnl[m]<1.e-8){
      diatridiag(m);
      return m;
    }
    if((m%5)==0){
      diatridiag(m);
      check=true;
      if(m==5){
	check=false;
	for(j=0;j<neig;j++)
	  eval[j]=eigval[j];
      }
      else{
	for(j=0;j<2;j++){
	  if(fabs(eigval[j]-eval[j])>1.e-8)check=false;
	  eval[j]=eigval[j];
	}
      }
      if(check&&compute_eigenvectorr(m))return m;
    }
  }
  diatridiag(mlanc-1);
  check=compute_eigenvectorr(mlanc-1);
  return mlanc-1;
}

//--------------------------------------------------------------------------------------
int fci::power1real(){
//--------------------------------------------------------------------------------------
  int i;
  unsigned long j;
  double nor1,nor2,overlap;
  for(i=0;i<mlanc;i++){
    hamiltonian_vector_multiplication(ffr[i%2],ffr[1-(i%2)]);
    nor1=inner_prod(nrep,ffr[1-(i%2)],ffr[1-(i%2)]);    
    nor2=inner_prod(nrep,ffr[1-(i%2)],ffr[i%2]);    
    nor1=sqrt(nor1);
    overlap=fabs(nor2/nor1);
    cout<<setprecision(12)<<"i="<<i<<"\teigval="<<nor1<<"\toverlap="<<overlap<<endl;
    if(fabs(1-overlap)<1.e-8)break;
    //#pragma omp parallel for
    for(j=0;j<nrep;j++)
      ffr[1-(i%2)][j]/=nor1;
  }
}

//--------------------------------------------------------------------------------------
void fci::diatridiag(int n){
//--------------------------------------------------------------------------------------
  double *d,*e,*work;
  int i,j,info;
  ofstream fout;
  char sector[10],base[100],name[100],xx[10],yy[10],ppx[10],ppy[10],zz[10],ssx[10],ssy[10],rr[10];
  d=new double[n];
  e=new double[n];
  work=new double[2*n];
  for(i=0;i<n;i++){
    d[i]=aal[i];
    e[i]=nnl[i+1];
  }
  dstev_("V",&n,d,e,lanczos_vec,&n,work,&info);
  for(i=0;i<n;i++){
    eigval[i]=d[i];
  }
  sprintf(sector,"%d",sec);
  sprintf(xx,"%d",kx);
  sprintf(yy,"%d",ky);
  sprintf(ppx,"%d",px);
  sprintf(ppy,"%d",py);
  sprintf(zz,"%d",qz);
  sprintf(ssx,"%d",sgm1);
  sprintf(ssy,"%d",sgm2);
  sprintf(rr,"%d",qr);
  strcpy(name,sector);
  strcat(name,"-");
  strcat(name,xx);
  strcat(name,"-");
  strcat(name,yy);
  strcat(name,"-");
  strcat(name,ppx);
  strcat(name,"-");
  strcat(name,ppy);
  strcat(name,"-");
  strcat(name,zz);
  strcat(name,"-");
  strcat(name,ssx);
  strcat(name,"-");
  strcat(name,ssy);
  strcat(name,"-");
  strcat(name,rr);
  strcpy(base,"eigval-");
  strcat(base,name);
  strcat(base,".dat");
  fout.open(base,ios::out);
  for(j=0;j<neig;j++)
    fout<<j<<"\t"<<setprecision(12)<<eigval[j]<<endl;
  fout.close();

  strcpy(base,"out-");
  strcat(base,name);
  strcat(base,".dat");
  fout.open(base,ios::app);
  fout<<n<<"\t";
  for(j=0;j<neig;j++)
    fout<<setprecision(12)<<eigval[j]<<"\t";
  fout<<endl;
  fout.close();
  delete []d;
  delete []e;
  delete []work;
}

//--------------------------------------------------------------------------------------
void fci::save_eigval(){
//--------------------------------------------------------------------------------------
  int i, j;
  ofstream fout;
  char base[100],name[100],sector[10],xx[10],yy[10],ppx[10],ppy[10],zz[10],ssx[10],ssy[10],rr[10];
  sprintf(sector,"%d",sec);
  sprintf(xx,"%d",kx);
  sprintf(yy,"%d",ky);
  sprintf(ppx,"%d",px);
  sprintf(ppy,"%d",py);
  sprintf(zz,"%d",qz);
  sprintf(ssx,"%d",sgm1);
  sprintf(ssy,"%d",sgm2);
  sprintf(rr,"%d",qr);
  strcpy(name,"eigval-");
  strcat(name,sector);
  strcat(name,"-");
  strcat(name,xx);
  strcat(name,"-");
  strcat(name,yy);
  strcat(name,"-");
  strcat(name,ppx);
  strcat(name,"-");
  strcat(name,ppy);
  strcat(name,"-");
  strcat(name,zz);
  strcat(name,"-");
  strcat(name,ssx);
  strcat(name,"-");
  strcat(name,ssy);
  strcat(name,"-");
  strcat(name,rr);
  strcpy(base,name);
  strcat(name,".dat");
  fout.open(name,ios::out);
  if(fout.is_open()){
    if(nrep>nrep_ed_max)
      for(i=0;i<neig;i++)
	//fout<<i<<"\t"<<setprecision(12)<<eigval[i]<<"\t"<<totspn[i]<<endl;
	fout<<i<<"\t"<<setprecision(12)<<eigval[i]<<endl;
    else
      for(i=0;i<nrep;i++)
	//fout<<i<<"\t"<<setprecision(12)<<eigval[i]<<"\t"<<totspn[i]<<endl;
	fout<<i<<"\t"<<setprecision(12)<<eigval[i]<<endl;
    fout.close();
  }
  strcpy(name,base);
  strcat(name,".bin");
  fout.open(name,ios::out | ios::binary);
  if(fout.is_open()){
    if(nrep>nrep_ed_max){
      if(kx%(perix/2)==0&&ky%(periy/2)==0&&qr%2==0)
	fout.write((char*)eigvecr[0],nrep*sizeof(double));    
      else
	fout.write((char*)eigvec[0],nrep*sizeof(dcmplex));
    }
    else{
      if(kx%(perix/2)==0&&ky%(periy/2)==0&&qr%2==0)
	fout.write((char*)matr,nrep*sizeof(double));    
      else
	fout.write((char*)mat,nrep*sizeof(dcmplex));
    }
    fout.close();
  }
}

//--------------------------------------------------------------------------------------
void fci::read_eigenvector(){
//--------------------------------------------------------------------------------------
  int i, j;
  ifstream fin;
  char base[100],name[100],sector[10],xx[10],yy[10],ppx[10],ppy[10],zz[10],ssx[10],ssy[10],rr[10];
  double r;
  sprintf(sector,"%d",sec);
  sprintf(xx,"%d",kx);
  sprintf(yy,"%d",ky);
  sprintf(ppx,"%d",px);
  sprintf(ppy,"%d",py);
  sprintf(zz,"%d",qz);
  sprintf(ssx,"%d",sgm1);
  sprintf(ssy,"%d",sgm2);
  sprintf(rr,"%d",qr);
  strcpy(name,"eigval-");
  strcat(name,sector);
  strcat(name,"-");
  strcat(name,xx);
  strcat(name,"-");
  strcat(name,yy);
  strcat(name,"-");
  strcat(name,ppx);
  strcat(name,"-");
  strcat(name,ppy);
  strcat(name,"-");
  strcat(name,zz);
  strcat(name,"-");
  strcat(name,ssx);
  strcat(name,"-");
  strcat(name,ssy);
  strcat(name,"-");
  strcat(name,rr);
  strcpy(base,name);
  strcpy(name,base);
  strcat(name,".bin");
  fin.open(name,ios::in | ios::binary);
  if(fin.is_open()){
    if(nrep>nrep_ed_max){
      if(kx%(perix/2)==0&&ky%(periy/2)==0&&qr%2==0){
	fin.read((char*)eigvecr[0],nrep*sizeof(double));
	measure_bond_enr(eigvecr[0]);
	r=measure_rotation_quantum_number(eigvecr[0]);
      }
      else{
	fin.read((char*)eigvec[0],nrep*sizeof(dcmplex));
	//measure_bond_enr(eigvec[0]);
      }
    }
    else{
      if(kx%(perix/2)==0&&ky%(periy/2)==0&&qr%2==0){
	fin.read((char*)matr,nrep*sizeof(double));
	measure_bond_enr(matr);
	r=measure_rotation_quantum_number(matr);
      }
      else{
	fin.read((char*)mat,nrep*sizeof(dcmplex));
	//measure_bond_enr(mat);
      }
    }
    fin.close();
  }
  else{
    cout<<"couldn't read eigenvector"<<endl;
  }
  cout<<"rotation quantum number is "<<r<<endl;
}

//--------------------------------------------------------------------------------------
bool fci::compute_eigenvector(int mlanc_curr){
//--------------------------------------------------------------------------------------
  unsigned long i;
  int j,k;
  double nor1,nor2,tots1,tots2,overlap[10];
  dcmplex *vectmp;
  bool check;
  ofstream fout;
  char sector[10],name[100],xx[10],yy[10],ppx[10],ppy[10],zz[10],ssx[10],ssy[10],rr[10];
  for(k=0;k<neig;k++)
    //#pragma omp parallel for 
    for(i=0;i<nrep;i++)
      eigvec[k][i]=0;
  for(k=0;k<neig;k++)
    for(j=0;j<mlanc_curr;j++)
      //#pragma omp parallel for 
      for(i=0;i<nrep;i++)
	eigvec[k][i]+=lanczos_vec[j+k*mlanc_curr]*ff[j][i];
  cout<<"done compute eigenvector"<<endl;
  vectmp=new dcmplex[nrep];
  check=true;
  sprintf(sector,"%d",sec);
  sprintf(xx,"%d",kx);
  sprintf(yy,"%d",ky);
  sprintf(ppx,"%d",px);
  sprintf(ppy,"%d",py);
  sprintf(zz,"%d",qz);
  sprintf(ssx,"%d",sgm1);
  sprintf(ssy,"%d",sgm2);
  sprintf(rr,"%d",qr);
  strcpy(name,"out-");
  strcat(name,sector);
  strcat(name,"-");
  strcat(name,xx);
  strcat(name,"-");
  strcat(name,yy);
  strcat(name,"-");
  strcat(name,ppx);
  strcat(name,"-");
  strcat(name,ppy);
  strcat(name,"-");
  strcat(name,zz);
  strcat(name,"-");
  strcat(name,ssx);
  strcat(name,"-");
  strcat(name,ssy);
  strcat(name,"-");
  strcat(name,rr);
  strcat(name,".dat");
  fout.open(name,ios::app);
  for(j=0;j<2;j++){
    hamiltonian_vector_multiplication(eigvec[j],vectmp);
    nor1=real(inner_prod(nrep,vectmp,vectmp));
    nor2=real(inner_prod(nrep,vectmp,eigvec[j]));
    overlap[j]=fabs(nor2/sqrt(nor1));
    //tots1=measure_total_spin(eigvec[j]);
    //cout<<setprecision(12)<<j<<"\tcheck eigenvector eigval="<<sqrt(nor1)<<"\toverlap="<<overlap[j]<<"\ttotal_spin="<<tots1<<endl;
    fout<<setprecision(12)<<j<<"\tcheck eigenvector eigval="<<sqrt(nor1)<<"\toverlap="<<overlap[j]<<endl;
    if(fabs(1-overlap[j])>1.e-7)check=false;
  }
  fout.close();
  delete []vectmp;
  return check;
}

//--------------------------------------------------------------------------------------
bool fci::compute_eigenvectorr(int mlanc_curr){
//--------------------------------------------------------------------------------------
  unsigned long i;
  int j,k;
  double *vectmp,nor1,nor2,tots1,overlap[10];
  bool check;
  ofstream fout;
  char sector[10],name[100],xx[10],yy[10],ppx[10],ppy[10],zz[10],ssx[10],ssy[10],rr[10];
  for(k=0;k<neig;k++)
    //#pragma omp parallel for 
    for(i=0;i<nrep;i++)
      eigvecr[k][i]=0;
  for(k=0;k<neig;k++)
    for(j=0;j<mlanc_curr;j++)
      //#pragma omp parallel for 
      for(i=0;i<nrep;i++)
	eigvecr[k][i]+=lanczos_vec[j+k*mlanc_curr]*ffr[j][i];
  cout<<"done compute eigenvector"<<endl;
  vectmp=new double[nrep];
  check=true;
  sprintf(sector,"%d",sec);
  sprintf(xx,"%d",kx);
  sprintf(yy,"%d",ky);
  sprintf(ppx,"%d",px);
  sprintf(ppy,"%d",py);
  sprintf(zz,"%d",qz);
  sprintf(ssx,"%d",sgm1);
  sprintf(ssy,"%d",sgm2);
  sprintf(rr,"%d",qr);
  strcpy(name,"out-");
  strcat(name,sector);
  strcat(name,"-");
  strcat(name,xx);
  strcat(name,"-");
  strcat(name,yy);
  strcat(name,"-");
  strcat(name,ppx);
  strcat(name,"-");
  strcat(name,ppy);
  strcat(name,"-");
  strcat(name,zz);
  strcat(name,"-");
  strcat(name,ssx);
  strcat(name,"-");
  strcat(name,ssy);
  strcat(name,"-");
  strcat(name,rr);
  strcat(name,".dat");
  fout.open(name,ios::app);
  for(j=0;j<2;j++){
    hamiltonian_vector_multiplication(eigvecr[j],vectmp);
    nor1=inner_prod(nrep,vectmp,vectmp);
    nor2=inner_prod(nrep,vectmp,eigvecr[j]);
    overlap[j]=fabs(nor2/sqrt(nor1));
    fout<<setprecision(12)<<j<<"\tcheck eigenvector eigval="<<sqrt(nor1)<<"\toverlap="<<overlap[j]<<endl;
    if(fabs(1-overlap[j])>1.e-7)check=false;
  }
  fout<<"done compute eigenvector"<<endl;
  fout.close();
  delete []vectmp;
  return check;
}

//--------------------------------------------------------------------------------------
void fci::make_lookup_table(){
//--------------------------------------------------------------------------------------
  unsigned long *reprdup,nt,i,j,n,n1;
  double *normdup;
  double t1,t2,t3;
  //t1=omp_get_wtime();
  nt=256;
  ntable=nt;
  lookuptable=new unsigned long[nt+1];
  reprdup=new unsigned long[nrep];
  normdup=new double[nrep];
  n=0;
  lookuptable[0]=0;
  for(i=0;i<nt;i++){
    n1=0;
    for(j=0;j<nrep;j++)
      if((repr[j]%nt)==i){
	reprdup[n+n1]=repr[j];
	normdup[n+n1]=normfac[j];
	n1++;
      }
    n+=n1;
    lookuptable[i+1]=n;
    //cout<<"i="<<i<<"\tlookuptable n1="<<n1<<endl;
  }
  if(lookuptable[nt]!=nrep){
    cout<<"make lookup table wrong"<<endl;
    exit(0);
  }
  //#pragma omp parallel for default(shared) schedule(dynamic,1)
  for(i=0;i<nrep;i++){
    repr[i]=reprdup[i];
    normfac[i]=normdup[i];
  }
  delete []reprdup;
  delete []normdup;
  //#pragma omp parallel for default(shared) private(n,n1)  schedule(dynamic,1)
  for(i=0;i<nt;i++){
    n=lookuptable[i];
    n1=lookuptable[i+1]-lookuptable[i];
    idsort(n1,&(repr[n]),&(normfac[n]));
  }
  //t2=omp_get_wtime();
  cout<<"resort nrep time="<<t2-t1<<endl;
}

//--------------------------------------------------------------------------------------
void idsort(unsigned long n, unsigned long *arr1, double *arr3){
//--------------------------------------------------------------------------------------
//bubble sorting method
  unsigned long i,j,igap,tmp;
  double dtmp;
  if(n<2)return;
  igap=n/2;
  do{
    for(i=igap;i<n;i++){
      j=i-igap;
    l50:
      if(arr1[j]>arr1[j+igap]){
	tmp=arr1[j];
	arr1[j]=arr1[j+igap];
	arr1[j+igap]=tmp;
	dtmp=arr3[j];
	arr3[j]=arr3[j+igap];
	arr3[j+igap]=dtmp;
      }
      else continue;
      if(j>=igap){
	j-=igap;
	goto l50;
      }
    }
    igap/=2;
  }while(igap>0);
  for(i=0;i<n-1;i++){
    if(arr1[i]>arr1[i+1])
      cout<<i<<"\t"<<arr1[i]<<endl;
  }
}
