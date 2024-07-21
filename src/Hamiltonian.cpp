#include "Hamiltonian.hpp"
namespace TightBinding
{
void Hamiltonian::PrintGeometryTotal(){
    int numberofatoms=0;
    for(auto &cell : systemcell.supercellMargin){
	numberofatoms+=cell->AtomNum;
    }
    for(auto &cell : systemcell.supercell){
	numberofatoms+=cell->AtomNum;
    }

    FILE *geo;
    geo=fopen("./Totalgeometry.data","w");
    fprintf(geo,"#BIND_OUTPUT version: 2.1\n");
    fprintf(geo,"#NUMBER OF ATOMS:\n");
    fprintf(geo,"\t%d\n",numberofatoms);
    fprintf(geo,"#ATOMIC POSITIONS\n");
    int o=0;
    for(auto &cell : systemcell.supercellMargin){
	for(int i=0;i<cell->AtomNum;++i){
	    o+=1;
	    Eigen::Vector3d pos=cell->Unit[i];
	    fprintf(geo,"   %d",o);
	    if(i%2){fprintf(geo,"   As");}else{fprintf(geo,"   Ga");}
	    fprintf(geo,"   %2.5f   %2.5f   %2.5f",pos(0)*5,pos(1)*5,pos(2)*5);
	    fprintf(geo,"\n");
	}
    }

    for(auto &cell : systemcell.supercell){
	for(int i=0;i<cell->AtomNum;++i){
	    o+=1;
	    Eigen::Vector3d pos=cell->Unit[i];
	    fprintf(geo,"   %d",o);
	    if(i%2){fprintf(geo,"   As");}else{fprintf(geo,"   Ga");}
	    //fprintf(geo,"   %2s",(cell->UnitTexture[i]).c_str());
	    fprintf(geo,"   %2.5f   %2.5f   %2.5f",pos(0)*5,pos(1)*5,pos(2)*5);
	    fprintf(geo,"\n");
	}
    }


    fclose(geo);
}

void Hamiltonian::PrintGeometry(std::vector<std::unique_ptr<Geometry>> & printCell){
    //center();
    int numberofatoms=0;
    double avgX=0;
    double avgY=0;
    double avgZ=0;
    for(auto &cell : printCell){
	for(int i=0;i<cell->AtomNum;++i){
	    avgX+=cell->Unit[i](0);
	    avgY+=cell->Unit[i](1);
	    avgZ+=cell->Unit[i](2);
	}
	numberofatoms+=cell->AtomNum;
    }
    FILE *geo;
    geo=fopen("./geometry.data","w");
    fprintf(geo,"#BIND_OUTPUT version: 2.1\n");
    fprintf(geo,"#NUMBER OF ATOMS:\n");
    fprintf(geo,"\t%d\n",numberofatoms);
    fprintf(geo,"#ATOMIC POSITIONS\n");
    int o=0;
    for(auto &cell : printCell){
	for(int i=0;i<cell->AtomNum;++i){
	    o+=1;
	    Eigen::Vector3d pos=cell->Unit[i];
	    pos(0)-=avgX/numberofatoms;
	    pos(1)-=avgY/numberofatoms;
	    pos(2)-=avgZ/numberofatoms;
	    fprintf(geo,"   %d",o);
	    //if(i%2){fprintf(geo,"   As");}else{fprintf(geo,"   Ga");}
/*
	    if(cell->UnitType[i]=="c"){
		fprintf(geo,"   %2s",((std::string)cell->texture).substr(0,2).c_str());
	    }else{
		fprintf(geo,"   %2s",((std::string)cell->texture).substr(2,2).c_str());
	    }
*/
	    fprintf(geo,"   %2s","Ga");
	    fprintf(geo,"   %2.5f   %2.5f   %2.5f",pos(0)*5,pos(1)*5,pos(2)*5);
	    fprintf(geo,"\n");
	}
    }

    fclose(geo);
}

    /*
       void Hamiltonian::FoldTheCol2(std::vector<double> &egv,std::vector<Eigen::MatrixXcd> &egvct,int cellme,int blpo,int AtomO,std::vector<int> vblpo,std::vector<int> vAtomO,const Eigen::Vector3d& k,std::vector<Eigen::Vector3d> chain,std::pair<double,double> EnergyCutoff,int newcolnum,int colnum,Eigen::Matrix3d &Ugm,Eigen::Matrix3d &ugm,Eigen::Matrix3d &uam,Eigen::MatrixXcd &yfitting,int dtag)
       {
       int UnitCellNumber=systemcell.UnitCellNumber;
       Eigen::Vector3d shift=systemcell.supercell[cellme]->Unit[0];
       assert(vblpo.size()==vAtomO.size());
       int blnum=vblpo.size();

       Eigen::VectorXcd retVxcd;
       retVxcd.resize(AtomO);
       retVxcd.setZero();
       int fitposCol=0;
       std::complex<double> ii(0.0,1.0);
       int fitposRow=0;

       int i=dtag;
       for(int j=0;j<blnum;++j)
       {
       Eigen::Vector3d ddd=chain[(i)*systemcell.UnitCellNumber]-shift;
       Eigen::Vector3d AtomPos=chain[(i)*systemcell.UnitCellNumber]-shift;
       Eigen::Vector3d dk(0,0,0);

       Eigen::Vector3d projAtomPos=uam*AtomPos;
       Eigen::Vector3d poshift((int)projAtomPos(0),(int)projAtomPos(1),(int)projAtomPos(2));

       int LpV,minDistance;
       std::vector<int> pushV2(systemcell.pushVec.begin()+1,systemcell.pushVec.end()-1);
       int bpushV2=pushV2[0];
       for(auto &pv : pushV2)
       {
       pv=pv-bpushV2;
       }

       for(int ini=0;ini<pushV2.size();++ini)
       {
       if(pushV2[ini]==poshift(0))
       {
       minDistance=ini;
       LpV=pushV2[minDistance];
       break;
       }
       if(pushV2[ini]>poshift(0))
       {
       minDistance=ini-1;
       LpV=pushV2[minDistance];
       break;
       }
       }
       int cNUM=systemcell.uvector.size()*systemcell.UnitCellNumber*orbNum;
    //std::cout<<"cNum="<<systemcell.uvector.size()*systemcell.UnitCellNumber*orbNum<<std::endl;
    dk=ugm*(AtomPos-uam.inverse()*poshift)+Ugm*poshift;//dk1
    for(int c=0;c<UnitCellNumber;++c)
    {
    Eigen::Vector3d posc=chain[(j)*UnitCellNumber+c]-shift;
    std::complex<double> expr=exp(((dk).dot(posc))*2*M_PI*ii);
    for(int ki=0;ki<orbNum;++ki)
    {
    int rowIndex=c*orbNum+ki;

    if(
    ( (fitposRow+rowIndex)>=(LpV)*cNUM
    && (fitposRow+rowIndex)<(LpV+systemcell.spNum[minDistance+1])*cNUM)
    || USE_RowTruncation==0
    )
    {
    yfitting(fitposRow+rowIndex,newcolnum)=egvct[i](rowIndex,colnum)*expr/sqrt(blnum);//5
    }
    }

    }
    fitposRow+=vAtomO[j];
}

}
*/
/*
   void Hamiltonian::PrepareHamiltonian(Eigen::MatrixXcd &blockMatrix,int tot)
   {

   int ZEROi = 0;
   int MONEi=-1; 
   int ONEi=1;    
   double ONEd=1.0;
   double ZEROd=0.0;

   int N=tot; 
   int M=tot; 
   int Nb=20; 
   int Mb=20; 
   int MM=-1,NZ=-1; 


   double MPIelapsed;
   double MPIt2;
   double MPIt1;

   int id,numprocs;
   blacs_pinfo_(&id, &numprocs);
//ctxt, myrow, mycol;
blacs_get_(&MONEi, &ZEROi, &ctxt);
bool mpiroot = (id == 0);
int info; 

descA=(int *)calloc(Nb,sizeof(int));
descZ=(int *)calloc(Mb,sizeof(int));


double *Ar=(double *)calloc(N*M,sizeof(double));
double *Ai=(double *)calloc(N*M,sizeof(double));

if (mpiroot) {
for(int i=0;i<N;++i)
{
for(int j=0;j<M;++j)
{
Ar[i+j*M]=blockMatrix(i,j).real();
Ai[i+j*M]=blockMatrix(i,j).imag();
}
}
}
//MPIt1=MPI_Wtime();
int procrows = sqrt(numprocs), proccols = sqrt(numprocs);
//double abstol=pdlamch(&ctxt,"U");
blacs_gridinit_(&ctxt, "R", &procrows, &proccols);  
blacs_gridinfo_(&ctxt, &procrows, &proccols, &myrow, &mycol );  
nrows = numroc_(&N, &Nb, &myrow, &ZEROi, &procrows);
int lda = std::max(1,nrows);
ncols = numroc_(&M, &Mb, &mycol, &ZEROi, &proccols);
descinit_(descA, &N, &M, &Nb, &Mb,&ZEROi,&ZEROi,&ctxt, &lda, &info);
descinit_(descZ, &N, &M, &Nb, &Mb,&ZEROi,&ZEROi,&ctxt, &lda, &info);

double *Ar_loc=(double *)calloc(nrows*ncols,sizeof(double));
double *Ai_loc=(double *)calloc(nrows*ncols,sizeof(double));


blacs_barrier_(&ctxt,"A");
sendr = 0, sendc = 0, recvr = 0, recvc = 0;
for (int r = 0; r < N; r += Nb, sendr=(sendr+1)%procrows) {
sendc = 0;
int nr = Nb;
if (N-r < Nb)
nr = N-r;
for (int c = 0; c < M; c += Mb, sendc=(sendc+1)%proccols) {
int nc = Mb;
if (M-c < Mb)
nc = M-c;
if (mpiroot) {
    dgesd2d_(&ctxt, &nr, &nc, Ai+N*c+r, &N, &sendr, &sendc);
}
if (myrow == sendr && mycol == sendc) {
    dgerv2d_(&ctxt, &nr, &nc, Ai_loc+nrows*recvc+recvr, &nrows,&ZEROi,&ZEROi);
    recvc = (recvc+nc)%ncols;
}
}
if (myrow == sendr)
    recvr = (recvr+nr)%nrows;
    }
blacs_barrier_(&ctxt,"A");
sendr = 0, sendc = 0, recvr = 0, recvc = 0;
for (int r = 0; r < N; r += Nb, sendr=(sendr+1)%procrows) {
    sendc = 0;
    int nr = Nb;
    if (N-r < Nb)
        nr = N-r;
    for (int c = 0; c < M; c += Mb, sendc=(sendc+1)%proccols) {
        int nc = Mb;
        if (M-c < Mb)
            nc = M-c;
        if (mpiroot) {
            dgesd2d_(&ctxt, &nr, &nc, Ar+N*c+r, &N, &sendr, &sendc);
        }
        if (myrow == sendr && mycol == sendc) {
            dgerv2d_(&ctxt, &nr, &nc, Ar_loc+nrows*recvc+recvr, &nrows,&ZEROi,&ZEROi);
            recvc = (recvc+nc)%ncols;
        }
    }
    if (myrow == sendr)
        recvr = (recvr+nr)%nrows;
}

blacs_barrier_(&ctxt,"A");

A_loc=(lapack_complex_double *)calloc(nrows*ncols,sizeof(lapack_complex_double));
Z_loc=(lapack_complex_double *)calloc(nrows*ncols,sizeof(lapack_complex_double));

blacs_barrier_(&ctxt,"A");
for(int i=0;i<nrows*ncols;++i)
{
    A_loc[i].real=Ar_loc[i];
    A_loc[i].imag=Ai_loc[i];
}
blacs_barrier_(&ctxt,"A");

free(Ai);
free(Ar);
free(Ai_loc);
free(Ar_loc);
}

void Hamiltonian::GatherHamiltonian(Eigen::VectorXd &WW,Eigen::MatrixXcd &ZZ,int tot)
{

    ZZ.resize(0,0);
    WW.resize(0);
    int id,numprocs;
    blacs_pinfo_(&id, &numprocs);
    int procrows = sqrt(numprocs), proccols = sqrt(numprocs);
    bool mpiroot = (id == 0);
    int N=tot; 
    int M=tot; 
    WW.resize(N);
    WW.setZero();
    ZZ.resize(N,M);
    ZZ.setZero();
    int Nb=20; 
    int Mb=20; 
    std::complex<double> ii(0.,1.);
    int ZEROi = 0;
    int MONEi=-1; 
    int ONEi=1;    
    double ONEd=1.0;
    double ZEROd=0.0;
    for(int i=0;i<N;++i)
    {
        WW[i]=W[i];
    }
    free(W);

    double *Zr_loc=(double *)calloc(nrows*ncols,sizeof(double));
    double *Zi_loc=(double *)calloc(nrows*ncols,sizeof(double));
    for(int i=0;i<nrows*ncols;++i)
    {
        Zr_loc[i]=Z_loc[i].real;
        Zi_loc[i]=Z_loc[i].imag;
    }
    free(Z_loc);

    double *Zr=(double *)calloc(N*M,sizeof(double));
    double *Zi=(double *)calloc(N*M,sizeof(double));
    //Gather matrix 
    blacs_barrier_(&ctxt,"A");
    sendr = 0;
    for (int r = 0; r < N; r += Nb, sendr=(sendr+1)%procrows) {
        sendc = 0;
        int nr = Nb;
        if (N-r < Nb)
            nr = N-r;
        for (int c = 0; c < M; c += Mb, sendc=(sendc+1)%proccols) {
            int nc = Mb;
            if (M-c < Mb)
                nc = M-c;
            if (myrow == sendr && mycol == sendc) {
                dgesd2d_(&ctxt, &nr, &nc, Zr_loc+nrows*recvc+recvr, &nrows,&ZEROi,&ZEROi);
                recvc = (recvc+nc)%ncols;
            }
            if (mpiroot) {
                dgerv2d_(&ctxt, &nr, &nc, Zr+N*c+r, &N, &sendr, &sendc);
            }
        }
        if (myrow == sendr)
            recvr = (recvr+nr)%nrows;
    }
    blacs_barrier_(&ctxt,"A");
    sendr = 0;
    for (int r = 0; r < N; r += Nb, sendr=(sendr+1)%procrows) {
        sendc = 0;
        int nr = Nb;
        if (N-r < Nb)
            nr = N-r;
        for (int c = 0; c < M; c += Mb, sendc=(sendc+1)%proccols) {
            int nc = Mb;
            if (M-c < Mb)
                nc = M-c;
            if (myrow == sendr && mycol == sendc) {
                dgesd2d_(&ctxt, &nr, &nc, Zi_loc+nrows*recvc+recvr, &nrows,&ZEROi,&ZEROi);
                recvc = (recvc+nc)%ncols;
            }
            if (mpiroot) {
                dgerv2d_(&ctxt, &nr, &nc, Zi+N*c+r, &N, &sendr, &sendc);
            }
        }

        if (myrow == sendr)
            recvr = (recvr+nr)%nrows;
    }
    blacs_barrier_(&ctxt,"A");
    //if(bcast_wave_func){
    MPI_Bcast(Zi,N*M,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Bcast(Zr,N*M,MPI_DOUBLE,0,MPI_COMM_WORLD);
    blacs_barrier_(&ctxt,"A");
    for(int i=0;i<N;++i)
    {
        for(int j=0;j<M;++j)
        {
            ZZ(i,j)=Zr[j*M+i]+Zi[j*M+i]*ii;
        }
    }
    blacs_barrier_(&ctxt,"A");
    //}




    free(Zi);
    free(Zr);
    free(Zi_loc);
    free(Zr_loc);

    blacs_gridexit_(&ctxt);
}
*/

inline std::complex<double> Hamiltonian::g0(const Eigen::Vector3d& k)
{
    const Eigen::Vector3d v = k * M_PI / 2.;

    return std::complex<double>(cos(v(0)) * cos(v(1)) * cos(v(2)), -sin(v(0)) * sin(v(1)) * sin(v(2)));
}

inline std::complex<double> Hamiltonian::g1(const Eigen::Vector3d& k)
{
    const Eigen::Vector3d v = k * M_PI / 2.;

    return std::complex<double>(-cos(v(0)) * sin(v(1)) * sin(v(2)), sin(v(0)) * cos(v(1)) * cos(v(2)));
}


inline std::complex<double> Hamiltonian::g2(const Eigen::Vector3d& k)
{
    const Eigen::Vector3d v = k * M_PI / 2.;

    return std::complex<double>(-sin(v(0)) * cos(v(1)) * sin(v(2)), cos(v(0)) * sin(v(1)) * cos(v(2)));
}


inline std::complex<double> Hamiltonian::g3(const Eigen::Vector3d& k)
{
    const Eigen::Vector3d v = k * M_PI / 2.;

    return std::complex<double>(-sin(v(0)) * sin(v(1)) * cos(v(2)), cos(v(0)) * cos(v(1)) * sin(v(2)));
}




/*
   void Hamiltonian::Method2(Eigen::MatrixXcd matTOT,Eigen::MatrixXcd &tt,Eigen::MatrixXcd &rr,Eigen::MatrixXcd &matII,std::vector<Eigen::MatrixXcd> vT,Eigen::Vector3d pl,Eigen::MatrixXcd T1,Eigen::MatrixXcd T2)
   {
   int bksize=systemcell.orbNum;
   int scsize=systemcell.supercell.size()*2.0;

   matTOT.block(0,0,bksize,bksize*scsize).setZero();
   matTOT.block(bksize*(scsize-1),0,bksize,bksize*scsize).setZero();

   Eigen::MatrixXcd DL=vT.back()*vT.front();
   Eigen::MatrixXcd DR=vT.back()*vT.front();
   csolver.compute(DL,Eigen::DecompositionOptions::ComputeEigenvectors);
   Eigen::VectorXcd tmpVX=csolver.eigenvalues();
   DL=csolver.eigenvectors();
//ArrangeD(tmpVX,DL,pl,0);
csolver.compute(DR,Eigen::DecompositionOptions::ComputeEigenvectors);
tmpVX=csolver.eigenvalues();
DR=csolver.eigenvectors();
//ArrangeD(tmpVX,DR,pl,0);


matTOT.block(0,0,bksize,bksize)=Eigen::MatrixXcd::Identity(bksize,bksize);
matTOT.block(0,bksize,bksize,bksize)=
-DL.topRightCorner(bksize,bksize)*DL.bottomRightCorner(bksize,bksize).inverse();
matTOT.block(bksize*(scsize-1),bksize*(scsize-2),bksize,bksize)=
-DR.bottomLeftCorner(bksize,bksize)*DL.topLeftCorner(bksize,bksize).inverse();
matTOT.block(bksize*(scsize-1),bksize*(scsize-1),bksize,bksize)=Eigen::MatrixXcd::Identity(bksize,bksize);

Eigen::MatrixXcd II;
II.resize(scsize*bksize,bksize);
II.setZero();
II.block(0,0,bksize,bksize)=
(DL.topLeftCorner(bksize,bksize)-DL.topRightCorner(bksize,bksize)*DL.bottomRightCorner(bksize,bksize).inverse()*DL.bottomLeftCorner(bksize,bksize))*
//Eigen::MatrixXcd::Ones(bksize,bksize);
Eigen::MatrixXcd::Identity(bksize,bksize);
II=matTOT.inverse()*II;

matII=II;
tt=DR.bottomLeftCorner(bksize,bksize).inverse()*II.block(bksize*(scsize-1),0,bksize,bksize);//tt=DR^-1*Cn;

rr=DL.topRightCorner(bksize,bksize).inverse()*(II.block(0,0,bksize,bksize)-DL.topLeftCorner(bksize,bksize));
//rr=DL.bottomRightCorner(bksize,bksize).inverse()*(II.block(bksize,0,bksize,bksize)-DL.bottomLeftCorner(bksize,bksize));

}
void Hamiltonian::Method4(Eigen::MatrixXcd matTOT,Eigen::MatrixXcd &tt,Eigen::MatrixXcd &rr,std::vector<Eigen::MatrixXcd> vT,Eigen::Vector3d pl,Eigen::MatrixXcd T1,Eigen::MatrixXcd T2)
{

int bksize=systemcell.orbNum;
int scsize=systemcell.supercell.size()*2.0;
Eigen::MatrixXcd minusDL=matTOT.block(0,bksize*(scsize-1),bksize,bksize);
Eigen::MatrixXcd middleDL=matTOT.block(0,0,bksize,bksize);
Eigen::MatrixXcd plusDL=matTOT.block(0,bksize,bksize,bksize);
Eigen::MatrixXcd DL;
DL.resize(bksize*2,bksize*2);
DL.setZero();
DL.block(0,0,bksize,bksize)=-minusDL.inverse()*middleDL;
DL.block(0,bksize,bksize,bksize)=-minusDL.inverse()*plusDL;
DL.block(bksize,0,bksize,bksize)=Eigen::MatrixXcd::Identity(bksize,bksize);
for(int i=0;i<vT.size();++i)
{
vT[i].resize(bksize*2,bksize*2);
vT[i].setZero();
Eigen::MatrixXcd minus=matTOT.block(bksize*(i+1),bksize*(i),bksize,bksize);
Eigen::MatrixXcd middle=matTOT.block(bksize*(i+1),bksize*(i+1),bksize,bksize);
Eigen::MatrixXcd plus=matTOT.block(bksize*(i+1),bksize*(i+2),bksize,bksize);
vT[i].block(0,0,bksize,bksize)=-minus.inverse()*middle;
vT[i].block(0,bksize,bksize,bksize)=-minus.inverse()*plus;
vT[i].block(bksize,0,bksize,bksize)=Eigen::MatrixXcd::Identity(bksize,bksize);
}
Eigen::MatrixXcd minusDR=matTOT.block(bksize*(scsize-1),bksize*(scsize-2),bksize,bksize);
Eigen::MatrixXcd middleDR=matTOT.block(bksize*(scsize-1),bksize*(scsize-1),bksize,bksize);
Eigen::MatrixXcd plusDR=matTOT.block(bksize*(scsize-1),0,bksize,bksize);
Eigen::MatrixXcd DR;
DR.resize(bksize*2,bksize*2);
DR.setZero();
DR.block(0,0,bksize,bksize)=-minusDR.inverse()*middleDR;
DR.block(0,bksize,bksize,bksize)=-minusDR.inverse()*plusDR;
DR.block(bksize,0,bksize,bksize)=Eigen::MatrixXcd::Identity(bksize,bksize);
matTOT.block(0,0,bksize,bksize*scsize).setZero();
matTOT.block(bksize*(scsize-1),0,bksize,bksize*scsize).setZero();

DL=DL*vT.front();
DR=vT.back()*DR;
csolver.compute(DL,Eigen::DecompositionOptions::ComputeEigenvectors);
Eigen::VectorXcd tmpVX=csolver.eigenvalues();
DL=csolver.eigenvectors();
//ArrangeD(tmpVX,DL,pl,0);
csolver.compute(DR,Eigen::DecompositionOptions::ComputeEigenvectors);
tmpVX=csolver.eigenvalues();
DR=csolver.eigenvectors();
//ArrangeD(tmpVX,DR,pl,0);


Eigen::MatrixXcd Ttot=Eigen::MatrixXcd::Identity(2*bksize,2*bksize);
for(int i=0;i<vT.size();++i)
{
    Ttot=Ttot*vT[i];
}
Eigen::MatrixXcd MTtot=Ttot;

Ttot=DL.inverse()*Ttot*DR;
tt=Ttot.topLeftCorner(bksize,bksize).inverse();
rr=Ttot.bottomLeftCorner(bksize,bksize)*tt;

Eigen::MatrixXcd II;
II.resize(bksize*2,bksize);
II.block(0,0,bksize,bksize)=rr;
II.block(bksize,0,bksize,bksize)=tt;
for(int i=0;i<bksize;++i)
{
    II.col(i)=II.col(i)/sqrt(II.col(i).squaredNorm());
}
tt=II.block(bksize,0,bksize,bksize);
rr=II.block(0,0,bksize,bksize);


std::vector<Eigen::MatrixXcd> vM;
for(int i=0;i<vT.size()-1;++i)
{
    //Eigen::MatrixXcd Mtmp=(vT[i]*vT[i+1]).inverse();
    Eigen::MatrixXcd Mtmp=vT[i]*vT[i+1];
    Eigen::MatrixXcd MMtmp=Mtmp.block(0,0,bksize,bksize*2);
    Mtmp.block(0,0,bksize,bksize*2)=Mtmp.block(bksize,0,bksize,bksize*2);
    Mtmp.block(bksize,0,bksize,bksize*2)=MMtmp;

    Eigen::MatrixXcd Mtmp2=Mtmp;
    Mtmp2.topLeftCorner(bksize,bksize)=Mtmp.bottomLeftCorner(bksize,bksize)*Mtmp.topLeftCorner(bksize,bksize).inverse();
    Mtmp2.topRightCorner(bksize,bksize)=Mtmp.bottomRightCorner(bksize,bksize)-Mtmp.bottomLeftCorner(bksize,bksize)*Mtmp.topLeftCorner(bksize,bksize).inverse()*Mtmp.topRightCorner(bksize,bksize);
    Mtmp2.bottomLeftCorner(bksize,bksize)=Mtmp.topLeftCorner(bksize,bksize).inverse();
    Mtmp2.bottomRightCorner(bksize,bksize)=-Mtmp.topRightCorner(bksize,bksize);
    vM.push_back(Mtmp2);
}

Eigen::MatrixXcd tmpS;
Eigen::MatrixXcd S=vM[0];
for(int i=1;i<vM.size();++i)
{
    Eigen::MatrixXcd Sm=S;
    Eigen::MatrixXcd Sp=vM[i];
    Eigen::MatrixXcd Sm11=Sm.topLeftCorner(bksize,bksize);
    Eigen::MatrixXcd Sm12=Sm.topRightCorner(bksize,bksize);
    Eigen::MatrixXcd Sm21=Sm.bottomLeftCorner(bksize,bksize);
    Eigen::MatrixXcd Sm22=Sm.bottomRightCorner(bksize,bksize);
    Eigen::MatrixXcd Sp11=Sp.topLeftCorner(bksize,bksize);
    Eigen::MatrixXcd Sp12=Sp.topRightCorner(bksize,bksize);
    Eigen::MatrixXcd Sp21=Sp.bottomLeftCorner(bksize,bksize);
    Eigen::MatrixXcd Sp22=Sp.bottomRightCorner(bksize,bksize);
    Eigen::MatrixXcd RsD=Sm.topRightCorner(bksize,bksize)*((Eigen::MatrixXcd::Identity(bksize,bksize)-Sp.topLeftCorner(bksize,bksize)*Sm.bottomRightCorner(bksize,bksize)).inverse());
    Eigen::MatrixXcd RsF=Sp.bottomLeftCorner(bksize,bksize)*((Eigen::MatrixXcd::Identity(bksize,bksize)-Sm.bottomRightCorner(bksize,bksize)*Sp.topLeftCorner(bksize,bksize)).inverse());

    S.topLeftCorner(bksize,bksize)=Sm11+RsD*Sp11*Sm21;
    S.topRightCorner(bksize,bksize)=RsD*Sp12;
    S.bottomLeftCorner(bksize,bksize)=RsF*Sm21;
    S.bottomRightCorner(bksize,bksize)=Sp22+RsF*Sm22*Sp12;

}

//Eigen::MatrixXcd tmpSS=S.block(0,0,bksize,bksize*2);
//S.block(0,0,bksize,bksize*2)=S.block(bksize,0,bksize,bksize*2);
//S.block(bksize,0,bksize,bksize*2)=tmpSS;


MTtot.topLeftCorner(bksize,bksize)=Ttot.bottomLeftCorner(bksize,bksize)*Ttot.topLeftCorner(bksize,bksize).inverse();
MTtot.topRightCorner(bksize,bksize)=Ttot.bottomRightCorner(bksize,bksize)-Ttot.bottomLeftCorner(bksize,bksize)*Ttot.topLeftCorner(bksize,bksize).inverse()*Ttot.topRightCorner(bksize,bksize);
MTtot.bottomLeftCorner(bksize,bksize)=Ttot.topLeftCorner(bksize,bksize).inverse();
MTtot.bottomRightCorner(bksize,bksize)=-Ttot.topRightCorner(bksize,bksize);



//print_matrix(MTtot,"MTtot");
//print_matrix(S,"S");
//print_matrix(S-MTtot,"test");

}

void Hamiltonian::Method3(Eigen::MatrixXcd matTOT,Eigen::MatrixXcd &tt,Eigen::MatrixXcd &rr,Eigen::MatrixXcd &matII,std::vector<Eigen::MatrixXcd> vT,Eigen::Vector3d pl,Eigen::MatrixXcd T1,Eigen::MatrixXcd T2)
{

    int bksize=systemcell.orbNum;
    int scsize=systemcell.supercell.size()*2.0;

    Eigen::MatrixXcd DL=vT.back()*vT.front();
    Eigen::MatrixXcd DR=vT.back()*vT.front();
    csolver.compute(DL,Eigen::DecompositionOptions::ComputeEigenvectors);
    Eigen::VectorXcd tmpVX=csolver.eigenvalues();
    DL=csolver.eigenvectors();
    //ArrangeD(tmpVX,DL,pl,0);
    csolver.compute(DR,Eigen::DecompositionOptions::ComputeEigenvectors);
    tmpVX=csolver.eigenvalues();
    DR=csolver.eigenvectors();
    //ArrangeD(tmpVX,DR,pl,0);

    Eigen::MatrixXcd Ttot=Eigen::MatrixXcd::Identity(2*bksize,2*bksize);
    for(int i=1;i<vT.size()-1;++i)
    {
        Ttot=Ttot*vT[i];
    }
    Ttot=DL.inverse()*Ttot*DR;
    tt=Ttot.topLeftCorner(bksize,bksize).inverse();
    rr=Ttot.bottomLeftCorner(bksize,bksize)*tt;
}
*/

void Hamiltonian::SetMatrix(const Eigen::Vector3d& k)
{
    MatrixSize=10;
    matrix.resize(MatrixSize,MatrixSize);
    matrix.setZero();
    const std::complex<double> g0v = g0(k);
    const std::complex<double> g1v = g1(k);
    const std::complex<double> g2v = g2(k);
    const std::complex<double> g3v = g3(k);

    const std::complex<double> g0vc = std::conj(g0v);
    const std::complex<double> g1vc = std::conj(g1v);
    const std::complex<double> g2vc = std::conj(g2v);
    const std::complex<double> g3vc = std::conj(g3v);

    const Material& m_material=material.materialsShort.at("InAs");
    // only the lower triangular of matrix is set because the diagonalization method only needs that

    // first column - 0 - for sa
    matrix(0, 0) = m_material.m_Esa;
    matrix(1, 0) = m_material.m_Vss * g0vc;
    matrix(2, 0) = 0;
    matrix(3, 0) = 0;
    matrix(4, 0) = 0;
    matrix(5, 0) = m_material.m_Vsapc * g1vc;
    matrix(6, 0) = m_material.m_Vsapc * g2vc;
    matrix(7, 0) = m_material.m_Vsapc * g3vc;
    matrix(8, 0) = 0;
    matrix(9, 0) = 0;

    // col 1 - for sc
    matrix(1, 1) = m_material.m_Esc;
    matrix(2, 1) = -m_material.m_Vscpa * g1v;
    matrix(3, 1) = -m_material.m_Vscpa * g2v;
    matrix(4, 1) = -m_material.m_Vscpa * g3v;
    matrix(5, 1) = 0;
    matrix(6, 1) = 0;
    matrix(7, 1) = 0;
    matrix(8, 1) = 0;
    matrix(9, 1) = 0;

    // col 2 - for pxa
    matrix(2, 2) = m_material.m_Epa;
    matrix(3, 2) = 0;
    matrix(4, 2) = 0;
    matrix(5, 2) = m_material.m_Vxx * g0vc;
    matrix(6, 2) = m_material.m_Vxy * g3vc;
    matrix(7, 2) = m_material.m_Vxy * g2vc;
    matrix(8, 2) = 0;
    matrix(9, 2) = -m_material.m_Vpasstarc * g1vc;

    // col 3 - for pya
    matrix(3, 3) = m_material.m_Epa;
    matrix(4, 3) = 0;
    matrix(5, 3) = m_material.m_Vxy * g3vc;
    matrix(6, 3) = m_material.m_Vxx * g0vc;
    matrix(7, 3) = m_material.m_Vxy * g1vc;
    matrix(8, 3) = 0;
    matrix(9, 3) = -m_material.m_Vpasstarc * g2vc;

    // col 4 - for pza
    matrix(4, 4) = m_material.m_Epa;
    matrix(5, 4) = m_material.m_Vxy * g2vc;
    matrix(6, 4) = m_material.m_Vxy * g1vc;
    matrix(7, 4) = m_material.m_Vxx * g0vc;
    matrix(8, 4) = 0;
    matrix(9, 4) = -m_material.m_Vpasstarc * g3vc;


    // col 5 - for pxc
    matrix(5, 5) = m_material.m_Epc;
    matrix(6, 5) = 0;
    matrix(7, 5) = 0;
    matrix(8, 5) = m_material.m_Vsstarapc * g1v;
    matrix(9, 5) = 0;

    // col 6 - for pyc
    matrix(6, 6) = m_material.m_Epc;
    matrix(7, 6) = 0;
    matrix(8, 6) = m_material.m_Vsstarapc * g2v;
    matrix(9, 6) = 0;

    // col 7 - for pzc
    matrix(7, 7) = m_material.m_Epc;
    matrix(8, 7) = m_material.m_Vsstarapc * g3v;
    matrix(9, 7) = 0;

    // col 8 - for s*a
    matrix(8, 8) = m_material.m_Esstara;
    matrix(9, 8) = 0. * g0vc; // it's set to 0

    // col 9 - for s*c
    matrix(9, 9) = m_material.m_Esstarc;
}





void Hamiltonian::printOutWaveFunctions(){
    int pos=0;
    FILE * wave_atom;
    std::fstream file_obj; 
    file_obj.open("./wave.dat",std::ios::out|std::ios::trunc);
    file_obj.close();

    file_obj.open("./wave.dat",std::ios::out|std::ios::app);
    for(auto &cell : systemcell.supercell){
        for(int i=0;i<cell->Unit.size();++i){
            file_obj<<pos<<"\t"<<cell->UnitWave_Functions[i].real()<<"\t"<<cell->UnitWave_Functions[i].imag()<<std::endl;
            //fprintf(wave_atom,"%d\t%1.5f\n",pos,cell->UnitWave_Functions[i].real());
            pos+=1;
        }
    }
    file_obj.close();
}

void Hamiltonian::WaveFunctions(){
    Eigen::MatrixXcd matSfull=gsolver.eigenvectors();
    int totalatom=0;

    for(auto &cell : systemcell.supercell){
        totalatom+=cell->AtomNum;
        cell->UnitWave_Functions.clear();
        cell->UnitWave_Functions.resize(cell->AtomNum);
    }
    //Eigen::VectorXcd vecS=matSfull.col(0);
    int pos=0;
    Eigen::VectorXcd vecS;
    //for(int at=0;at<totalatom*BandGapIndex/*matSfull.cols()*/;++at){
    for(int at=totalatom*BandGapIndex-1;at<totalatom*BandGapIndex;++at){
        //for(int at=totalatom*BandGapIndex;at<totalatom*BandGapIndex+1;++at){
        //for(int at=totalatom*BandGapIndex+1;at<totalatom*BandGapIndex+2;++at){
        //for(int at=totalatom*BandGapIndex+2;at<totalatom*BandGapIndex+3;++at){
        //for(int at=0;at<1;++at){
        pos=0;
        vecS=matSfull.col(at);
        for(auto &cell : systemcell.supercell){
            for(int i=0;i<cell->Unit.size();++i){
                std::complex<double> tr=vecS.segment(pos,orbNum).adjoint()*vecS.segment(pos,orbNum);
                cell->UnitWave_Functions[i]+=tr;
                pos+=orbNum;
            }
        }
    }

    }

    void Hamiltonian::GetUq(Eigen::Vector3d q){
        std::complex<double> ii(0.0,1.0);
        int dimX=MatrixSize;
        Uq.resize(dimX,dimX);
        Uq.setZero();
        std::unordered_map<std::string, Material> &ms=material.materials;
        if(USE_long_parameters==0)
        {
        ms=material.materialsShort;
        }
        for(auto &cell : systemcell.ListofMatrixElements_Type1){
            Eigen::Vector3d ed=cell->ed;
            Eigen::Vector3d dd=cell->d;
            Eigen::Vector3d dn=dd.normalized();
            //std::complex<double> efactor=(dd.norm()!=0)?pow(2*M_PI,-1.5)*exp((2*M_PI*(dd.dot(q)))*ii)/dd.norm():0;
            std::complex<double> efactor=exp(2.*M_PI*((dd).dot(q))*ii)/(dd.norm()+1e-6);
            //if(dd.norm()==0) efactor=0;
            int atom2=cell->Table.real();
            int atom1=cell->Table.imag();
            std::string leftType=(atom1%2==0?"c":"a");
            std::string rightType=(atom2%2==0?"c":"a");
            atom2*=orbNum;
            atom1*=orbNum;
            for(int io=0;io<orbNum;++io){
                for(int jo=0;jo<orbNum;++jo){
                    Eigen::Vector4i leftorbital=qlabels[jo];
                    Eigen::Vector4i rightorbital=qlabels[io];

                    if(leftorbital(3)==rightorbital(3))
                    {
                        int rows=atom2+io;
                        int cols=atom1+jo;
                        if(rows==cols)
                        {
                            Uq(rows,cols)=0.0;
                        }
                        else if(leftType!=rightType)
                        { 
                            double res=0;
                            if(rows>cols)
                            {
                                res=(SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->cellme)));
                            }else if(rows<cols)
                            {
                                res=(SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->SuperCell)));
                            }

                            Uq(rows,cols)+=efactor;

                        }
                        if(
                                atom2==atom1 
                                && leftorbital(1)==rightorbital(1) 
                                && leftorbital(1)==1
                                &&USE_SOC
                                && 0
                          )
                        {
                            if((leftorbital(2)+rightorbital(2))==0)
                            {
                                double lambda=ms.at(cell->cellme).param.at(leftType);
                                //matrix(rows,cols)+=ii*(0.5*rightorbital(2))*lambda;
                            }
                        }

                    }else{
                        int rows=atom2+io;
                        int cols=atom1+jo;

                        double lambda=1;
                        double m1=leftorbital(2);
                        double m2=rightorbital(2);
                        double leftsign=leftorbital(3);
                        double rightsign=rightorbital(3);
                        if(
                                atom2==atom1 
                                && leftorbital(1)==rightorbital(1) 
                                && leftorbital(1)==1
                                &&USE_SOC
                                && 1
                          )
                        {

                            assert(leftType==rightType);
                            lambda=ms.at(cell->cellme).param.at(leftType);
                            std::complex<double> soc;
                            std::complex<double> theta1=(1-(0==m1?1:0))*1.0/sqrt(2)*((m1>=0?1.0:0.0)-ii*((-m1)>=0?1.0:0.0))+0.5*(0==m1?1:0);
                            std::complex<double> theta2=(1-(0==m2?1:0))*1.0/sqrt(2)*((m2>=0?1.0:0.0)+ii*((-m2)>=0?1.0:0.0))+0.5*(0==m2?1:0);
                            std::complex<double> mu=theta1*theta2*(1+((-abs(m1*m2))>=0?1.0:0.0));
                            soc=mu*((abs(m1)==(abs(m2)+leftsign)?1.0:0.0)-pow(-1,(m1>=0?1:0)+(m2>=0?1:0))*(abs(m1)==(abs(m2)+rightsign)?1.0:0.0));

                            //matrix(rows,cols)+=lambda*soc;
                        }

                    }
/*
                    int rows=atom2+io;
                    int cols=atom1+jo;
                    if(rows==cols){
                        Uq(rows,cols)=7.64; 
                    }else if(atom1==atom2 && io!=jo){
                        Uq(rows,cols)=0; 
                    }else if(io==jo){

                        double res=0;
                        if(rows>cols)
                        {
                            res=(SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->cellme)));
                        }else if(rows<cols)
                        {
                            res=(SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->SuperCell)));
                        }

                        Uq(rows,cols)+=efactor; 
                    }
*/
                }
            }



        }

        //print_matrix(Uq,"Uq");
    }

    void Hamiltonian::ParseParameters(std::string inpfile)
    {
        orbNum=ParseBlockTitle<int>("Orbital",inpfile,0);
        qlabels.resize(orbNum);
        for(int i=0;i<orbNum;++i)
        {
            ParseBlockElement<int>("Orbital",inpfile,i,1,qlabels[i](0));
            ParseBlockElement<int>("Orbital",inpfile,i,2,qlabels[i](1));
            ParseBlockElement<int>("Orbital",inpfile,i,3,qlabels[i](2));
            ParseBlockElement<int>("Orbital",inpfile,i,4,qlabels[i](3));
        }
        if(USE_extra_k_point)
        {
            KextraNum=ParseBlockTitle<int>("Kextra",inpfile,0);
            KextraName.resize(KextraNum);
            KextraRange.resize(KextraNum);

            for(int i=0;i<KextraNum;++i)
            {
                ParseBlockElement<std::string>("Kextra",inpfile,i,0,KextraName[i]);
                ParseBlockElement<double>("Kextra",inpfile,i,1,KextraRange[i].first);
                ParseBlockElement<double>("Kextra",inpfile,i,2,KextraRange[i].second);
            }
        }
    }

    void Hamiltonian::get_transform_remote(Eigen::Vector3d kp,Eigen::MatrixXcd &transf,int displacement,int bandnumber,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements){
        Eigen::MatrixXcd RemoteMatrix;
        SetMatrix(kp,ListofMatrixElements,RemoteMatrix);
        preDiagonalize(gsoul,RemoteMatrix,false,matrix);

        transf.resize(MatrixSize,systemcell.supercell.size()*bandnumber);
        transf.setZero();
        for(auto &cell : systemcell.supercell){
            if(cell->type){
                int pp=cell->Tag2;
                int mAtom=cell->SecDim;
                int blpo=mposi(cell->Tag2);
                int posipp=0;
                std::for_each(systemcell.supercell.begin(),systemcell.supercell.end(),
                        [&posipp,pp]
                        (auto &scell)
                        {if(scell->Tag2<pp){posipp+=scell->AtomNum;}}
                        );
                transf.block(posipp*orbNum,pp*bandnumber,cell->AtomNum*orbNum,bandnumber)=gsoul[cell->Tag3]->eigenvectors().block(0,cell->AtomNum*BandGapIndex+displacement,cell->AtomNum*orbNum,bandnumber); 
            }
        }
    }


    void Hamiltonian::ReadMaterials()
    {
        //material.ReadMaterials("./MaterialData/Material.data","GaAs");
        //material.ReadMaterials("./MaterialData/Material.data","InAs");
        //material.ReadMaterials("./MaterialData/Material.data","AlAs");
        material.MixTwoMaterials("AlAs","GaAs",MixCoeff,"AlGaAs");
        material.MixTwoMaterialsShort("AlAs","GaAs",MixCoeff,"AlGaAs");
    }
    void Hamiltonian::destroy()
    {
        EigenVectorsBlocks.resize(0);
        EigenValuesBlocks.resize(0);
        matrix.resize(0,0);
        dMatrix.resize(0,0);
        matrix3Db.resize(0,0);
        matrixScatter.resize(0,0);
        gsoul.resize(0);
        gsoul3D.resize(0);
        systemcell.supercell.resize(0);
        systemcell.supercellMargin.resize(0);
        systemcell.spNum.resize(0);
        systemcell.fspNum.resize(0);
        systemcell.spStr.resize(0);

        systemcell.ListofMatrixElements_Type1.resize(0);
        systemcell.ListofMatrixElements_MixType1.resize(0);
        systemcell.ListofMatrixElements_Type2.resize(0);
        systemcell.ListofMatrixElements_Type3.resize(0);
        systemcell.ListofMatrixCoupling.resize(0);

    }
    Hamiltonian::~Hamiltonian()
    {
        matrix.resize(0,0);
        gsoul.resize(0);
        gsoul3D.resize(0);
        systemcell.supercell.resize(0);
        systemcell.supercellMargin.resize(0);
        systemcell.spNum.resize(0);
        systemcell.spStr.resize(0);

        systemcell.ListofMatrixElements_Type1.resize(0);
        systemcell.ListofMatrixElements_MixType1.resize(0);
        systemcell.ListofMatrixElements_Type2.resize(0);
        systemcell.ListofMatrixElements_Type3.resize(0);
        systemcell.ListofMatrixCoupling.resize(0);

    }

    inline int Hamiltonian::mposi(int ii)
    {
        int res=0;
        int types=systemcell.supercell[ii]->type;
        for(auto &cell : systemcell.supercell)
        {
            if(cell->Tag2<ii && cell->type==types)
            {
                res+=cell->SecDim;
            }
        }
        return res;
    }

    inline std::complex<double> Hamiltonian::MatrixCoupling11(const Eigen::MatrixXcd &lM ,int row,int col,const Eigen::MatrixXcd &rM,std::complex<int> & maplist)
    {
        const Eigen::VectorXcd & lV=lM.conjugate().col(row);
        const Eigen::VectorXcd & rV=rM.col(col);
        std::complex<double> result={0,0};
        for(auto &edge : thelist[maplist])
        {
            result+=lV(edge->Table1)*edge->value*rV(edge->Table2);
        } 
        return result;
    }
    void Hamiltonian::DestoryCouplingElements(Eigen::MatrixXcd &ConsideredMatrix)
    {
        for(auto &cell : systemcell.ListofMatrixCoupling)
        {

            int p=cell->TagPosition1;
            int pp=cell->TagPosition2;
            int atom2xo=cell->Table.real();
            int atom1xo=cell->Table.imag();
            int ppp=cell->cellme;
            int pTarget=cell->SuperCell;
            std::complex<int> maplist={cell->TagPosition1,cell->TagPosition2};

            for(int io=0;io<orbNum;++io)
            {
                for(int jo=0;jo<orbNum;++jo)
                {
                    int rows=atom2xo+io;
                    int cols=atom1xo+jo;
                    ConsideredMatrix(rows,cols)=0;
                }
            }
        }
    }
    void Hamiltonian::ReadCouplingElements2(Eigen::MatrixXcd &ConsideredMatrix)
    {
        for(auto it= thelist.begin();it!=thelist.end();++it)
        {
            it->second.clear();
        }
        thelist.clear();
        for(auto &cell : systemcell.ListofMatrixCoupling)
        {

            int p=cell->TagPosition1;
            int pp=cell->TagPosition2;
            int atom2xo=cell->Table.real();
            int atom1xo=cell->Table.imag();
            int ppp=cell->cellme;
            int pTarget=cell->SuperCell;
            std::complex<int> maplist={cell->TagPosition1,cell->TagPosition2};

            for(int io=0;io<orbNum;++io)
            {
                for(int jo=0;jo<orbNum;++jo)
                {
                    int rows=atom2xo+io;
                    int cols=atom1xo+jo;
                    thelist[maplist].push_back(std::make_unique<EdgeCoupling>(rows-ppp,cols-pTarget,ConsideredMatrix(rows,cols)));
                    //ConsideredMatrix(rows,cols)=0;
                }
            }
        }
    }
    void Hamiltonian::ReadCouplingElements3(Eigen::MatrixXcd &ConsideredMatrix)
    {
        for(auto it= thelist.begin();it!=thelist.end();++it)
        {
            it->second.clear();
        }
        thelist.clear();
        for(auto &cell : systemcell.ListofMatrixCoupling)
        {

            int p=cell->TagPosition1;
            int pp=cell->TagPosition2;
            int atom2xo=cell->Table.real();
            int atom1xo=cell->Table.imag();
            int ppp=cell->cellme;
            int pTarget=cell->SuperCell;
            std::complex<int> maplist={cell->TagPosition1,cell->TagPosition2};

            for(int io=0;io<orbNum;++io)
            {
                for(int jo=0;jo<orbNum;++jo)
                {
                    int rows=atom2xo+io;
                    int cols=atom1xo+jo;
                    std::complex<int> Index(rows,cols);
                    std::complex<double> element(0,0);
                    if(SparseMatrix.count(Index))
                    {
                        element={SparseMatrix[Index].real,SparseMatrix[Index].imag};
                    }
                    thelist[maplist].push_back(std::make_unique<EdgeCoupling>(rows-ppp,cols-pTarget,element));
                }
            }
        }
        SparseMatrix.clear();
    }

    void Hamiltonian::ReadCouplingElements(Eigen::MatrixXcd &ConsideredMatrix)
    {
        for(auto it= thelist.begin();it!=thelist.end();++it)
        {
            it->second.clear();
        }
        thelist.clear();
        for(auto &cell : systemcell.ListofMatrixCoupling)
        {

            int p=cell->TagPosition1;
            int pp=cell->TagPosition2;
            int atom2xo=cell->Table.real();
            int atom1xo=cell->Table.imag();
            int ppp=cell->cellme;
            int pTarget=cell->SuperCell;
            std::complex<int> maplist={cell->TagPosition1,cell->TagPosition2};

            for(int io=0;io<orbNum;++io)
            {
                for(int jo=0;jo<orbNum;++jo)
                {
                    int rows=atom2xo+io;
                    int cols=atom1xo+jo;
                    thelist[maplist].push_back(std::make_unique<EdgeCoupling>(rows-ppp,cols-pTarget,ConsideredMatrix(rows,cols)));
                    ConsideredMatrix(rows,cols)=0;
                    ConsideredMatrix(cols,rows)=0;
                }
            }
        }
    }
    void Hamiltonian::pzheevx3PreDiagonalize(std::vector<Eigen::MatrixXcd> &SubblockContainerOfEigenVectors,std::vector<Eigen::VectorXd> &SubblockContainerOfEigenValues,Eigen::MatrixXcd &ConsideredMatrix,bool onoff,Eigen::MatrixXcd &MatrixToBeTransformed,int flag)
    {
        int id,numprocs;
        std::complex<double> ii(0,1);
        blacs_pinfo_(&id, &numprocs);
        matrix_bottom=0;
        int size_of_blocks=0;

        for(auto &cell : systemcell.supercell)
        {
            if(cell->Tag3==size_of_blocks){size_of_blocks+=1;}
        }
        std::vector<Eigen::MatrixXcd> EigenVectors(size_of_blocks);
        std::vector<Eigen::VectorXd> EigenValues(size_of_blocks);
        //std::vector<Eigen::MatrixXcd> EigenVectors2(size_of_blocks);
        //std::vector<Eigen::VectorXd> EigenValues2(size_of_blocks);

        /*
           for(auto &cell : systemcell.supercell)
           {
           if(cell->Tag3==tmpo){
           tmpo+=1;
           solver.compute(ConsideredMatrix.block(cell->blpo,cell->blpo,cell->AtomO,cell->AtomO),Eigen::DecompositionOptions::ComputeEigenvectors);
        //EigenVectors[cell->Tag3]=(solver.eigenvectors());
        //EigenValues[cell->Tag3]=(solver.eigenvalues());
        }
        }
        */
        //Eigen::VectorXd tmp;
        Eigen::MatrixXd tmp2r;
        Eigen::MatrixXd tmp2i;
        int section=(size_of_blocks)/numprocs+1;

        //for(int i=0;i<size_of_blocks;++i)
        for(int i=std::min(id*section,size_of_blocks);i<std::min((id+1)*section,size_of_blocks);++i)
        {
            for(auto &cell : systemcell.supercell)
            {
                if(cell->Tag3==i && id==(i/section))
                {
                    solver.compute(ConsideredMatrix.block(cell->blpo,cell->blpo,cell->AtomO,cell->AtomO),Eigen::DecompositionOptions::ComputeEigenvectors);
                    EigenValues[cell->Tag3]=solver.eigenvalues();
                    EigenVectors[cell->Tag3]=solver.eigenvectors();
                    break;
                }
            }
        }

        for(int i=0;i<EigenValues.size();++i)
        {
            for(auto &cell : systemcell.supercell)
            {
                if(cell->Tag3==i)
                {
                    /*
                    //tmp.resize(cell->AtomO);
                    //if(id==i/section) tmp=EigenValues[i];
                    if(id==i/section) tmp=EigenValues[cell->Tag3];
                    MPI_Bcast(tmp.data(),cell->AtomO,MPI_DOUBLE,i/section,MPI_COMM_WORLD);
                    EigenValues[i]=tmp;
                    */
                    if(id!=i/section) EigenValues[cell->Tag3].resize(cell->AtomO);
                    MPI_Bcast(EigenValues[cell->Tag3].data(),cell->AtomO,MPI_DOUBLE,i/section,MPI_COMM_WORLD);

                    tmp2r.resize(cell->AtomO,cell->AtomO);
                    tmp2i.resize(cell->AtomO,cell->AtomO);
                    if(id==i/section) tmp2r=EigenVectors[cell->Tag3].real();
                    if(id==i/section) tmp2i=EigenVectors[cell->Tag3].imag();
                    MPI_Bcast(tmp2r.data(),cell->AtomO*cell->AtomO,MPI_DOUBLE,i/section,MPI_COMM_WORLD);
                    MPI_Bcast(tmp2i.data(),cell->AtomO*cell->AtomO,MPI_DOUBLE,i/section,MPI_COMM_WORLD);
                    EigenVectors[i]=tmp2r+tmp2i*ii;
                }
            }

        }

        SubblockContainerOfEigenVectors.clear();
        SubblockContainerOfEigenVectors.resize(size_of_blocks);
        SubblockContainerOfEigenValues.clear();
        SubblockContainerOfEigenValues.resize(size_of_blocks);
        int Atomtot=0;
        int tmpo=0;
        if(onoff) transform_matrix.resize(0,0);
        for(auto &cell : systemcell.supercell)
        {
            int ptag=cell->Tag2;
            int posip=0;
            if(onoff)
            {
                std::for_each(systemcell.supercell.begin(),systemcell.supercell.end(),
                        [&posip,ptag]
                        (auto &scell)
                        {if(scell->Tag2<ptag){posip+=scell->AtomNum;}}
                        );
            }
            if(cell->type==0)
            {
                cell->SecDim=cell->AtomO;
                cell->Startp=0;
                Atomtot+=cell->SecDim;
                matrix_bottom+=cell->Startp;
            }else if(cell->Tag3==/*SubblockContainerOfEigenValues.size()*/tmpo)
            {

                if(flag){
                    //Eigen::MatrixXcd Z,Z2;
                    //Eigen::VectorXd W,W2;
                    //Eigen::MatrixXcd blockMatrix=ConsideredMatrix.block(cell->blpo,cell->blpo,cell->AtomO,cell->AtomO);

                    //SubblockContainerOfEigenVectors.push_back(Eigen::MatrixXcd(0,0));
                    //SubblockContainerOfEigenValues.push_back(Eigen::VectorXd(0));
                    //Diagonal_parallel_pzheevx_all(SubblockContainerOfEigenValues.back(),SubblockContainerOfEigenVectors.back(),cell->blpo,0,cell->AtomO,1,1,MatrixSize,1);
                    tmpo+=1;
                    //SubblockContainerOfEigenVectors.push_back(Z2);
                    //SubblockContainerOfEigenValues.push_back(W2);
                }else{
                    //solver.compute(ConsideredMatrix.block(cell->blpo,cell->blpo,cell->AtomO,cell->AtomO),Eigen::DecompositionOptions::ComputeEigenvectors);
                    //SubblockContainerOfEigenVectors.push_back(solver.eigenvectors());
                    //SubblockContainerOfEigenValues.push_back(solver.eigenvalues());
                    SubblockContainerOfEigenVectors[cell->Tag3]=EigenVectors[cell->Tag3];
                    SubblockContainerOfEigenValues[cell->Tag3]=EigenValues[cell->Tag3];
                    tmpo+=1;
                }

                //SubblockContainerOfEigenVectors.push_back(Z.block(cell->blpo,cell->blpo,cell->AtomO,cell->AtomO));
                //SubblockContainerOfEigenValues.push_back(W.segment(cell->blpo,cell->AtomO));

                const Eigen::VectorXd &soup=SubblockContainerOfEigenValues[cell->Tag3];
                /*
                   const double *pfirst=soup.data();
                   const double *plast=soup.data()+soup.size();
                   int o=0,oo=0,oco=0;
                   for(;pfirst!=plast;++pfirst){
                   double eigvals=(*pfirst);
                   if(eigvals<leftCutoff){o+=1;}
                //else if(eigvals<=0){oo+=1;}
                //else if(eigvals<=1.2){oco+=1;}
                else if(eigvals<=rightCutoff){oo+=1;}
                else{break;}
                }
                cell->SecDim=oo;
                cell->Startp=o;
                */
                const auto upper=std::upper_bound(soup.data(),soup.data()+soup.size(),rightCutoff);
                const auto lower=std::lower_bound(soup.data(),upper,leftCutoff);
                cell->SecDim=std::distance(lower,upper);
                cell->Startp=std::distance(soup.data(),lower);

                Atomtot+=cell->SecDim;
                matrix_bottom+=cell->Startp;

            }else{
                for(auto &innercell : systemcell.supercell)
                {
                    if(innercell->Tag3==cell->Tag3)
                    {
                        cell->SecDim=innercell->SecDim;
                        cell->Startp=innercell->Startp;
                        break;
                    }
                }
                Atomtot+=cell->SecDim;
                matrix_bottom+=cell->Startp;

            }

            if(onoff)
            {
                int endcols=transform_matrix.cols();
                transform_matrix.conservativeResize(MatrixSize,endcols+cell->SecDim);
                transform_matrix.block(0,endcols,MatrixSize,cell->SecDim).setZero();
                transform_matrix.block(posip*orbNum,endcols,cell->AtomNum*orbNum,cell->SecDim)=SubblockContainerOfEigenVectors[cell->Tag3].block(0,cell->Startp,cell->AtomNum*orbNum,cell->SecDim); 
            }

        }

        if(onoff)
        {
            Eigen::MatrixXcd trMatrix2;
            trMatrix2.resize(MatrixSize,transform_X.cols()+transform_L.cols());
            trMatrix2.setZero();
            trMatrix2<< transform_X , transform_L;
            extendMat=trMatrix2.adjoint()*MatrixToBeTransformed*trMatrix2;
            extendSquareMat=trMatrix2.adjoint()*MatrixToBeTransformed*transform_matrix;

            transform_mix.resize(MatrixSize,transform_matrix.cols()+trMatrix2.cols());
            transform_mix.setZero();
            transform_mix<< transform_matrix , trMatrix2;
        }
        MatrixSizeReduced=Atomtot;
    }

    void Hamiltonian::pzheevxPreDiagonalize(const Eigen::Vector3d k,std::vector<Eigen::MatrixXcd> &SubblockContainerOfEigenVectors,std::vector<Eigen::VectorXd> &SubblockContainerOfEigenValues,Eigen::MatrixXcd &ConsideredMatrix,bool onoff,Eigen::MatrixXcd &MatrixToBeTransformed,int flag)
    {
        int id,numprocs;
        blacs_pinfo_(&id, &numprocs);
        matrix_bottom=0;
        SubblockContainerOfEigenVectors.clear();
        SubblockContainerOfEigenVectors.resize(0);
        SubblockContainerOfEigenValues.clear();
        SubblockContainerOfEigenValues.resize(0);
        int Atomtot=0;

        if(onoff) transform_matrix.resize(0,0);
        for(auto &cell : systemcell.supercell)
        {
            int ptag=cell->Tag2;
            int posip=0;
            if(onoff)
            {
                std::for_each(systemcell.supercell.begin(),systemcell.supercell.end(),
                        [&posip,ptag]
                        (auto &scell)
                        {if(scell->Tag2<ptag){posip+=scell->AtomNum;}}
                        );
            }
            if(cell->type==0)
            {
                cell->SecDim=cell->AtomO;
                cell->Startp=0;
                Atomtot+=cell->SecDim;
                matrix_bottom+=cell->Startp;
            }else if(cell->Tag3==SubblockContainerOfEigenValues.size())
            {

                if(flag){
                    //Eigen::MatrixXcd Z,Z2;
                    //Eigen::VectorXd W,W2;
                    //Eigen::MatrixXcd blockMatrix=ConsideredMatrix.block(cell->blpo,cell->blpo,cell->AtomO,cell->AtomO);

                    SubblockContainerOfEigenVectors.push_back(Eigen::MatrixXcd(0,0));
                    SubblockContainerOfEigenValues.push_back(Eigen::VectorXd(0));
                    Eigen::Matrix<lapack_complex_double ,-1,-1> blockMatrix;
                    SetMatrixLocNew(k,systemcell.ListofMatrixElements_Type1,blockMatrix,cell->blpo,cell->AtomO);
                    Diagonal_parallel_pzheevx(SubblockContainerOfEigenValues.back(),SubblockContainerOfEigenVectors.back(),0,0,cell->AtomO,1,1,cell->AtomO,1,blockMatrix);
                    //Diagonal_parallel_pzheevx_all(SubblockContainerOfEigenValues.back(),SubblockContainerOfEigenVectors.back(),cell->blpo,0,cell->AtomO,1,1,MatrixSize,1);
                    //SubblockContainerOfEigenVectors.push_back(Z2);
                    //SubblockContainerOfEigenValues.push_back(W2);
                }else{
                    solver.compute(ConsideredMatrix.block(cell->blpo,cell->blpo,cell->AtomO,cell->AtomO),Eigen::DecompositionOptions::ComputeEigenvectors);
                    SubblockContainerOfEigenVectors.push_back(solver.eigenvectors());
                    SubblockContainerOfEigenValues.push_back(solver.eigenvalues());
                }

                //SubblockContainerOfEigenVectors.push_back(Z.block(cell->blpo,cell->blpo,cell->AtomO,cell->AtomO));
                //SubblockContainerOfEigenValues.push_back(W.segment(cell->blpo,cell->AtomO));

                //const Eigen::VectorXd &soup=SubblockContainerOfEigenValues[cell->Tag3];
                /*
                   const double *pfirst=soup.data();
                   const double *plast=soup.data()+soup.size();
                   int o=0,oo=0,oco=0;
                   for(;pfirst!=plast;++pfirst){
                   double eigvals=(*pfirst);
                   if(eigvals<leftCutoff){o+=1;}
                //else if(eigvals<=0){oo+=1;}
                //else if(eigvals<=1.2){oco+=1;}
                else if(eigvals<=rightCutoff){oo+=1;}
                else{break;}
                }
                cell->SecDim=oo;
                cell->Startp=o;
                */
                /*
                   const auto upper=std::upper_bound(soup.data(),soup.data()+soup.size(),rightCutoff);
                   const auto lower=std::lower_bound(soup.data(),upper,leftCutoff);
                   cell->SecDim=std::distance(lower,upper);
                   cell->Startp=std::distance(soup.data(),lower);
                   */
                Filtrate_the_states(SubblockContainerOfEigenValues[cell->Tag3],cell,true);
                Atomtot+=cell->SecDim;
                matrix_bottom+=cell->Startp;

            }else{

                Filtrate_the_states(SubblockContainerOfEigenValues[cell->Tag3],cell,false);
                /*
                   for(auto &innercell : systemcell.supercell)
                   {
                   if(innercell->Tag3==cell->Tag3)
                   {
                   cell->SecDim=innercell->SecDim;
                   cell->Startp=innercell->Startp;
                   break;
                   }
                   }
                   */
                Atomtot+=cell->SecDim;
                matrix_bottom+=cell->Startp;

            }

            if(onoff)
            {
                int endcols=transform_matrix.cols();
                transform_matrix.conservativeResize(MatrixSize,endcols+cell->SecDim);
                transform_matrix.block(0,endcols,MatrixSize,cell->SecDim).setZero();
                transform_matrix.block(posip*orbNum,endcols,cell->AtomNum*orbNum,cell->SecDim)=SubblockContainerOfEigenVectors[cell->Tag3].block(0,cell->Startp,cell->AtomNum*orbNum,cell->SecDim); 
            }

        }

        if(onoff)
        {
            Eigen::MatrixXcd trMatrix2;
            trMatrix2.resize(MatrixSize,transform_X.cols()+transform_L.cols());
            trMatrix2.setZero();
            trMatrix2<< transform_X , transform_L;
            extendMat=trMatrix2.adjoint()*MatrixToBeTransformed*trMatrix2;
            extendSquareMat=trMatrix2.adjoint()*MatrixToBeTransformed*transform_matrix;

            transform_mix.resize(MatrixSize,transform_matrix.cols()+trMatrix2.cols());
            transform_mix.setZero();
            transform_mix<< transform_matrix , trMatrix2;
        }
        MatrixSizeReduced=Atomtot;
    }

    void Hamiltonian::Filtrate_the_states(const Eigen::VectorXd &soup,const std::unique_ptr<Geometry> &cell,bool newone,int &tmpStartp,int &tmpSecDim)
    {
        if(newone)
        {
            if(USE_OneRange_or_MoreRanges)
            {
                if(AutoMatic)
                {
                    if(USE_FixedSize_or_EnergyRange)
                    {
                        tmpSecDim=rightCutoff+leftCutoff;
                        tmpStartp=cell->AtomNum*BandGapIndex-leftCutoff;
                    }else{
                        const auto upper=std::upper_bound(soup.data(),soup.data()+soup.size(),rightCutoff);
                        const auto lower=std::lower_bound(soup.data(),upper,leftCutoff);
                        tmpSecDim=std::distance(lower,upper);
                        tmpStartp=std::distance(soup.data(),lower);
                    }
                }else{
                    if(USE_FixedSize_or_EnergyRange)
                    {
                        tmpSecDim=rightCutoff+leftCutoff;
                        tmpStartp=cell->AtomNum*BandGapIndex-leftCutoff;
                    }else{
                        const auto upper=std::upper_bound(soup.data(),soup.data()+soup.size(),rightCutoff);
                        const auto lower=std::lower_bound(soup.data(),upper,leftCutoff);
                        tmpSecDim=std::distance(lower,upper);
                        tmpStartp=std::distance(soup.data(),lower);
                    }
                }

            }else{
                if(AutoMatic){
                    if(USE_FixedSize_or_EnergyRange)
                    {

                    }else{

                    }
                }else{
                    if(USE_FixedSize_or_EnergyRange)
                    {
                        if(USE_decomposition)
                        {
                            tmpSecDim=cell->EnergyCutoff.first+cell->EnergyCutoff.second;
                            tmpStartp=(soup.size()/orbNum)*BandGapIndex-cell->EnergyCutoff.first;
                        }else{
                            tmpSecDim=cell->EnergyCutoff.first+cell->EnergyCutoff.second;
                            tmpStartp=cell->AtomNum*BandGapIndex-cell->EnergyCutoff.first;
                        }
                        //std::cout<<"new["<<cell->Tag2<<"]l="<<cell->EnergyCutoff.first<<"\tr="<<cell->EnergyCutoff.second<<std::endl;
                    }else{
                        const auto upper=std::upper_bound(soup.data(),soup.data()+soup.size(),cell->EnergyCutoff.second);
                        const auto lower=std::lower_bound(soup.data(),upper,cell->EnergyCutoff.first);
                        tmpSecDim=std::distance(lower,upper);
                        tmpStartp=std::distance(soup.data(),lower);
                    }
                }
            }
        }else{//oldone
            if(USE_OneRange_or_MoreRanges)
            {
                for(auto &innercell : systemcell.supercell)
                {
                    if(innercell->Tag3==cell->Tag3)
                    {
                        tmpSecDim=innercell->SecDim;
                        tmpStartp=innercell->Startp;
                        break;
                    }
                }

            }else{
                if(AutoMatic){
                    if(USE_FixedSize_or_EnergyRange)
                    {

                    }else{

                    }
                }else{
                    if(USE_FixedSize_or_EnergyRange)
                    {
                        if(USE_decomposition)
                        {
                            tmpSecDim=cell->EnergyCutoff.first+cell->EnergyCutoff.second;
                            tmpStartp=(soup.size()/orbNum)*BandGapIndex-cell->EnergyCutoff.first;
                        }else{
                            tmpSecDim=cell->EnergyCutoff.first+cell->EnergyCutoff.second;
                            tmpStartp=cell->AtomNum*BandGapIndex-cell->EnergyCutoff.first;
                        }
                        //std::cout<<"["<<cell->Tag3<<"]SecDim="<<cell->SecDim<<"\tStarp="<<cell->Startp<<std::endl;
                        //std::cout<<"old["<<cell->Tag2<<"]l="<<cell->EnergyCutoff.first<<"\tr="<<cell->EnergyCutoff.second<<std::endl;
                    }else{
                        if(soup.size())
                        {
                            const auto upper=std::upper_bound(soup.data(),soup.data()+soup.size(),cell->EnergyCutoff.second);
                            const auto lower=std::lower_bound(soup.data(),upper,cell->EnergyCutoff.first);
                            tmpSecDim=std::distance(lower,upper);
                            tmpStartp=std::distance(soup.data(),lower);
                        }else{
                            for(auto &innercell : systemcell.supercell)
                            {
                                if(innercell->Tag3==cell->Tag3)
                                {
                                    tmpSecDim=innercell->SecDim;
                                    tmpStartp=innercell->Startp;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }


    void Hamiltonian::Filtrate_the_states(const Eigen::VectorXd &soup,const std::unique_ptr<Geometry> &cell,bool newone)
    {
        if(newone)
        {
            //newone
            if(USE_OneRange_or_MoreRanges)
            {
                if(AutoMatic)
                {
                    if(USE_FixedSize_or_EnergyRange)
                    {
                        //FixedSize OneRange
                        cell->SecDim=rightCutoff+leftCutoff;
                        cell->Startp=cell->AtomNum*BandGapIndex-leftCutoff;
                    }else{
                        //EnergyRange OneRange
                        const auto upper=std::upper_bound(soup.data(),soup.data()+soup.size(),rightCutoff);
                        const auto lower=std::lower_bound(soup.data(),upper,leftCutoff);
                        cell->SecDim=std::distance(lower,upper);
                        cell->Startp=std::distance(soup.data(),lower);
                    }
                }else{
                    if(USE_FixedSize_or_EnergyRange)
                    {
                        //FixedSize OneRange
                        cell->SecDim=rightCutoff+leftCutoff;
                        cell->Startp=cell->AtomNum*BandGapIndex-leftCutoff;
                    }else{
                        //EnergyRange OneRange
                        const auto upper=std::upper_bound(soup.data(),soup.data()+soup.size(),rightCutoff);
                        const auto lower=std::lower_bound(soup.data(),upper,leftCutoff);
                        cell->SecDim=std::distance(lower,upper);
                        cell->Startp=std::distance(soup.data(),lower);
                    }
                }

            }else{
                if(AutoMatic){
                    if(USE_FixedSize_or_EnergyRange)
                    {

                    }else{

                    }
                }else{
                    if(USE_FixedSize_or_EnergyRange)
                    {
                        if(USE_decomposition)
                        {
                            //FixedSize MoreRange Decomposition
                            cell->SecDim=cell->EnergyCutoff.first+cell->EnergyCutoff.second;
                            cell->Startp=(soup.size()/orbNum)*BandGapIndex-cell->EnergyCutoff.first;
                        }else{
                            //FixedSize MoreRange No_Decomposition
                            cell->SecDim=cell->EnergyCutoff.first+cell->EnergyCutoff.second;
                            cell->Startp=cell->AtomNum*BandGapIndex-cell->EnergyCutoff.first;
                        }
                    }else{
                        //EnergyRange MoreRange Decomposition
                        const auto upper=std::upper_bound(soup.data(),soup.data()+soup.size(),cell->EnergyCutoff.second);
                        const auto lower=std::lower_bound(soup.data(),upper,cell->EnergyCutoff.first);
                        cell->SecDim=std::distance(lower,upper);
                        cell->Startp=std::distance(soup.data(),lower);
                    }
                }
            }
        }else{
            //oldone
            if(USE_OneRange_or_MoreRanges)
            {
                //OneRange
                for(auto &innercell : systemcell.supercell)
                {
                    if(innercell->Tag3==cell->Tag3)
                    {
                        cell->SecDim=innercell->SecDim;
                        cell->Startp=innercell->Startp;
                        break;
                    }
                }

            }else{
                if(AutoMatic){
                    if(USE_FixedSize_or_EnergyRange)
                    {

                    }else{

                    }
                }else{
                    if(USE_FixedSize_or_EnergyRange)
                    {
                        if(USE_decomposition)
                        {
                            //FixedSize MoreRange Decomposition
                            cell->SecDim=cell->EnergyCutoff.first+cell->EnergyCutoff.second;
                            cell->Startp=(soup.size()/orbNum)*BandGapIndex-cell->EnergyCutoff.first;
                        }else{
                            //FixedSize MoreRange No_Decomposition
                            cell->SecDim=cell->EnergyCutoff.first+cell->EnergyCutoff.second;
                            cell->Startp=cell->AtomNum*BandGapIndex-cell->EnergyCutoff.first;
                        }
                    }else{
                        if(soup.size())
                        {
                            //EnergyRange MoreRange
                            const auto upper=std::upper_bound(soup.data(),soup.data()+soup.size(),cell->EnergyCutoff.second);
                            const auto lower=std::lower_bound(soup.data(),upper,cell->EnergyCutoff.first);
                            cell->SecDim=std::distance(lower,upper);
                            cell->Startp=std::distance(soup.data(),lower);
                        }else{
                            for(auto &innercell : systemcell.supercell)
                            {
                                if(innercell->Tag3==cell->Tag3)
                                {
                                    cell->SecDim=innercell->SecDim;
                                    cell->Startp=innercell->Startp;
                                    break;
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    void Hamiltonian::pDiagonalize(std::vector<std::unique_ptr<Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd>>> &SubblockContainer,Eigen::MatrixXcd &ConsideredMatrix,bool onoff,Eigen::MatrixXcd &MatrixToBeTransformed)
    {
        matrix_bottom=0;
        SubblockContainer.clear();
        SubblockContainer.resize(0);
        int Atomtot=0;

        if(onoff) transform_matrix.resize(0,0);
        for(auto &cell : systemcell.supercell)
        {
            int ptag=cell->Tag2;
            int posip=0;
            if(onoff)
            {
                std::for_each(systemcell.supercell.begin(),systemcell.supercell.end(),
                        [&posip,ptag]
                        (auto &scell)
                        {if(scell->Tag2<ptag){posip+=scell->AtomNum;}}
                        );
            }
            if(cell->type==0)
            {
                cell->SecDim=cell->AtomO;
                cell->Startp=0;
                Atomtot+=cell->SecDim;
                matrix_bottom+=cell->Startp;
            }else if(cell->Tag3==SubblockContainer.size())
            {
                SubblockContainer.push_back(std::make_unique<Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd>>());
                SubblockContainer[cell->Tag3]->compute(ConsideredMatrix.block(cell->blpo,cell->blpo,cell->AtomO,cell->AtomO),Eigen::DecompositionOptions::ComputeEigenvectors);
                assert(SubblockContainer[cell->Tag3]->info() == Eigen::ComputationInfo::Success);

                Filtrate_the_states(SubblockContainer[cell->Tag3]->eigenvalues(),cell,true);

                Atomtot+=cell->SecDim;
                matrix_bottom+=cell->Startp;

            }else{
                Filtrate_the_states(SubblockContainer[cell->Tag3]->eigenvalues(),cell,false);
                /*
                   for(auto &innercell : systemcell.supercell)
                   {
                   if(innercell->Tag3==cell->Tag3)
                   {
                   cell->SecDim=innercell->SecDim;
                   cell->Startp=innercell->Startp;
                   break;
                   }
                   }
                   */
                Atomtot+=cell->SecDim;
                matrix_bottom+=cell->Startp;

            }

            if(onoff)
            {
                int endcols=transform_matrix.cols();
                transform_matrix.conservativeResize(MatrixSize,endcols+cell->SecDim);
                transform_matrix.block(0,endcols,MatrixSize,cell->SecDim).setZero();
                transform_matrix.block(posip*orbNum,endcols,cell->AtomNum*orbNum,cell->SecDim)=SubblockContainer[cell->Tag3]->eigenvectors().block(0,cell->Startp,cell->AtomNum*orbNum,cell->SecDim); 
            }

        }

        if(onoff)
        {
            Eigen::MatrixXcd trMatrix2;
            trMatrix2.resize(MatrixSize,transform_X.cols()+transform_L.cols());
            trMatrix2.setZero();
            trMatrix2<< transform_X , transform_L;
            extendMat=trMatrix2.adjoint()*MatrixToBeTransformed*trMatrix2;
            extendSquareMat=trMatrix2.adjoint()*MatrixToBeTransformed*transform_matrix;

            transform_mix.resize(MatrixSize,transform_matrix.cols()+trMatrix2.cols());
            transform_mix.setZero();
            transform_mix<< transform_matrix , trMatrix2;
        }
        MatrixSizeReduced=Atomtot;
    }

    void Hamiltonian::preDiagonalize(std::vector<std::unique_ptr<Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd>>> &SubblockContainer,Eigen::MatrixXcd &ConsideredMatrix,bool onoff,Eigen::MatrixXcd &MatrixToBeTransformed)
    {
        matrix_bottom=0;
        SubblockContainer.clear();
        SubblockContainer.resize(0);
        int Atomtot=0;

        if(onoff) transform_matrix.resize(0,0);
        for(auto &cell : systemcell.supercell)
        {
            int ptag=cell->Tag2;
            int posip=0;
            if(onoff)
            {
                std::for_each(systemcell.supercell.begin(),systemcell.supercell.end(),
                        [&posip,ptag]
                        (auto &scell)
                        {if(scell->Tag2<ptag){posip+=scell->AtomNum;}}
                        );
            }
            if(cell->type==0)
            {
                cell->SecDim=cell->AtomO;
                cell->Startp=0;
                Atomtot+=cell->SecDim;
                matrix_bottom+=cell->Startp;
            }else if(cell->Tag3==SubblockContainer.size())
            {
                SubblockContainer.push_back(std::make_unique<Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd>>());
                SubblockContainer[cell->Tag3]->compute(ConsideredMatrix.block(cell->blpo,cell->blpo,cell->AtomO,cell->AtomO),Eigen::MatrixXcd::Identity(cell->AtomO,cell->AtomO) ,Eigen::DecompositionOptions::ComputeEigenvectors|Eigen::Ax_lBx);
                assert(SubblockContainer[cell->Tag3]->info() == Eigen::ComputationInfo::Success);

                Filtrate_the_states(SubblockContainer[cell->Tag3]->eigenvalues(),cell,true);

                Atomtot+=cell->SecDim;
                matrix_bottom+=cell->Startp;

            }else{

                Filtrate_the_states(SubblockContainer[cell->Tag3]->eigenvalues(),cell,false);
                /*
                   for(auto &innercell : systemcell.supercell)
                   {
                   if(innercell->Tag3==cell->Tag3)
                   {
                   cell->SecDim=innercell->SecDim;
                   cell->Startp=innercell->Startp;
                   break;
                   }
                   }
                   */
                Atomtot+=cell->SecDim;
                matrix_bottom+=cell->Startp;

            }

            if(USE_extra_k_point)
            {

                for(int i=0;i<KextraNum;++i)
                {

                    //if((KextraRange[i].first+KextraRange[i].second)<=0)
                    if((KextraRange[i].first+KextraRange[i].second)<0 && KextraName[i]=="G")
                    {
                        matrix_bottom-=(KextraRange[i].second);
                    }

                }
            }

            if(onoff)
            {
                int endcols=transform_matrix.cols();
                transform_matrix.conservativeResize(MatrixSize,endcols+cell->SecDim);
                transform_matrix.block(0,endcols,MatrixSize,cell->SecDim).setZero();
                transform_matrix.block(posip*orbNum,endcols,cell->AtomNum*orbNum,cell->SecDim)=SubblockContainer[cell->Tag3]->eigenvectors().block(0,cell->Startp,cell->AtomNum*orbNum,cell->SecDim); 
            }

        }

        if(onoff)
        {
            Eigen::MatrixXcd trMatrix2;
            //trMatrix2.resize(MatrixSize,transform_X.cols()+transform_L.cols());
            //trMatrix2.setZero();
            //trMatrix2<< transform_X , transform_L;

            trMatrix2.resize(0,0);
            trMatrix2.setZero();
            for(int i=0;i<KextraNum;++i)
            {
                trMatrix2.conservativeResize(MatrixSize,trMatrix2.cols()+KextraTransform[i].cols());
                trMatrix2.block(0,trMatrix2.cols()-KextraTransform[i].cols(),MatrixSize,KextraTransform[i].cols())= KextraTransform[i];
            }

            extendMat=trMatrix2.adjoint()*MatrixToBeTransformed*trMatrix2;
            extendSquareMat=trMatrix2.adjoint()*MatrixToBeTransformed*transform_matrix;

            transform_mix.resize(MatrixSize,transform_matrix.cols()+trMatrix2.cols());
            transform_mix.setZero();
            transform_mix<< transform_matrix , trMatrix2;
        }
        MatrixSizeReduced=Atomtot;
    }


    void Hamiltonian::FillDiagonalBlocksByTransform(Eigen::MatrixXcd &ConsideredMatrix,std::vector<std::unique_ptr<Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd>>> &SubblockData)
    {
        for(auto &cell : systemcell.supercell)
        {
            if(cell->type)
            {
                int pp=cell->Tag2;
                int mAtom=cell->SecDim;
                int blpo=mposi(cell->Tag2);

                Eigen::MatrixXcd diag=
                    SubblockData[cell->Tag3]->eigenvectors().block(0,cell->Startp,cell->AtomO,mAtom).adjoint()
                    *ConsideredMatrix.block(cell->blpo,cell->blpo,cell->AtomO,cell->AtomO)
                    *SubblockData[cell->Tag3]->eigenvectors().block(0,cell->Startp,cell->AtomO,mAtom);

                ConsideredMatrix.block(cell->blpo,cell->blpo,cell->AtomO,cell->AtomO).setZero();
                ConsideredMatrix.block(blpo,blpo,mAtom,mAtom)=diag;
            }
        }
    }

    void Hamiltonian::SparseTransform(Eigen::MatrixXcd &ConsideredMatrix,std::vector<Eigen::MatrixXcd> &BlockData)
    {

        bool MethodFlag=0;//1 for sparse transform

        for(auto &cell : systemcell.supercell)
        {
            if(cell->type)
            {

                int pp=cell->Tag2;
                int mAtom=cell->SecDim;
                int blpo=mposi(cell->Tag2);
                Eigen::MatrixXcd diag;

                const Eigen::MatrixXcd & lM=BlockData[cell->Tag3].block(0,cell->Startp,cell->AtomO,mAtom);
                const Eigen::MatrixXcd & rM=BlockData[cell->Tag3].block(0,cell->Startp,cell->AtomO,mAtom);

                if(MethodFlag)
                {
                }else{
                    diag=
                        lM.adjoint()
                        *ConsideredMatrix.block(cell->blpo,cell->blpo,cell->AtomO,cell->AtomO)
                        *rM;
                }
                std::complex<double> ii(0,1);
                int ZEROi = 0;
                int MONEi=-1; 
                int ONEi=1;    
                double ONEd=1.0;
                double ZEROd=0.0;
                int N=mAtom; 
                int M=mAtom; 
                int Nb=1;
                int Mb=1;

                int id,numprocs;
                blacs_pinfo_(&id, &numprocs);
                int ctxt, myrow, mycol;
                blacs_get_(&MONEi, &ZEROi, &ctxt);
                int info; 

                int procrows = sqrt(numprocs), proccols = sqrt(numprocs);
                blacs_gridinit_(&ctxt, "R", &procrows, &proccols);  
                blacs_gridinfo_(&ctxt, &procrows, &proccols, &myrow, &mycol );  
                int nrows = numroc_(&N, &Nb, &myrow, &ZEROi, &procrows);
                int lda = std::max(1,nrows);
                int ncols = numroc_(&M, &Mb, &mycol, &ZEROi, &proccols);
                dMatrix2.resize(nrows,ncols);	
                dMatrix2.setConstant({0,0});


                for(int rows=0;rows<diag.rows();++rows)
                {
                    for(int cols=0;cols<diag.cols();++cols)
                    {
                        std::complex<double> element(0,0);
                        if(MethodFlag)
                        {
                            const Eigen::VectorXcd & lV=lM.conjugate().col(rows);
                            const Eigen::VectorXcd & rV=rM.col(cols);
                            for(auto &incell : SparseMatrix)
                            {
                                std::complex<double> result={incell.second.real,incell.second.imag};
                                element+=lV(incell.first.real())*result*rV(incell.first.imag());
                            }
                        }
                        int sendr=rows%procrows;
                        int sendc=cols%proccols;
                        int recvr;
                        int recvc;
                        if (myrow == sendr && mycol == sendc){
                            recvr=rows/procrows;
                            recvc=cols/proccols;
                        }
                        if(myrow==sendr && mycol == sendc){
                            if(MethodFlag)
                            {
                                dMatrix2(recvr,recvc).real=real(element);
                                dMatrix2(recvr,recvc).imag=imag(element);
                            }else{
                                dMatrix2(recvr,recvc).real=real(diag(rows,cols));
                                dMatrix2(recvr,recvc).imag=imag(diag(rows,cols));
                            }
                        }

                    }
                }
                /*
                   diag.resize(mAtom,mAtom);
                   diag.setZero();
                   for(int i=0;i<diag.rows();++i)
                   {
                   for(int j=0;j<diag.cols();++j)
                   {
                   const Eigen::VectorXcd & lV=lM.conjugate().col(i);
                   const Eigen::VectorXcd & rV=rM.col(j);
                   for(auto &incell : SparseMatrix)
                   {
                   std::complex<double> result={incell.second.real,incell.second.imag};
                   diag(i,j)+=lV(incell.first.real())*result*rV(incell.first.imag());
                   }
                   }
                   }
                   */
                //ConsideredMatrix.block(cell->blpo,cell->blpo,cell->AtomO,cell->AtomO).setZero();
                //ConsideredMatrix.block(blpo,blpo,mAtom,mAtom)=diag;
                diag.resize(0,0);
            }
        }
    }


    void Hamiltonian::FillDiagonalBlocksByTransform(Eigen::MatrixXcd &ConsideredMatrix,std::vector<Eigen::MatrixXcd> &BlockData)
    {
        for(auto &cell : systemcell.supercell)
        {
            if(cell->type)
            {

                int pp=cell->Tag2;
                int mAtom=cell->SecDim;
                int blpo=mposi(cell->Tag2);
                Eigen::MatrixXcd diag;
                transform_matrix=BlockData[cell->Tag3].block(0,cell->Startp,cell->AtomO,mAtom);
                diag=
                    BlockData[cell->Tag3].block(0,cell->Startp,cell->AtomO,mAtom).adjoint()
                    *ConsideredMatrix.block(cell->blpo,cell->blpo,cell->AtomO,cell->AtomO)
                    *BlockData[cell->Tag3].block(0,cell->Startp,cell->AtomO,mAtom);

                ConsideredMatrix.block(cell->blpo,cell->blpo,cell->AtomO,cell->AtomO).setZero();
                ConsideredMatrix.block(blpo,blpo,mAtom,mAtom)=diag;
            }
        }
    }
    void Hamiltonian::FillDiagonalBlocksByData(Eigen::MatrixXcd &ConsideredMatrix,std::vector<Eigen::VectorXd> &SubblockData)
    {
        for(auto &cell : systemcell.supercell)
        {
            if(cell->type)
            {
                int pp=cell->Tag2;
                int mAtom=cell->SecDim;
                int blpo=mposi(cell->Tag2);
                ConsideredMatrix.block(blpo,blpo,mAtom,mAtom)=SubblockData[cell->Tag3].segment(cell->Startp,mAtom).asDiagonal();
            }
        }
    }

    void Hamiltonian::pzheevxFillDiagonalBlocksByData(Eigen::MatrixXcd &ConsideredMatrix,std::vector<Eigen::VectorXd> &SubblockData)
    {

        int id,numprocs;
        int ZEROi = 0;
        int MONEi=-1; 
        int ctxt, myrow, mycol;
        blacs_pinfo_(&id, &numprocs);
        blacs_get_(&MONEi, &ZEROi, &ctxt);
        int procrows = sqrt(numprocs), proccols = sqrt(numprocs);
        blacs_gridinit_(&ctxt, "R", &procrows, &proccols);  
        blacs_gridinfo_(&ctxt, &procrows, &proccols, &myrow, &mycol );  

        for(auto &cell : systemcell.supercell)
        {
            if(cell->type)
            {
                int pp=cell->Tag2;
                int mAtom=cell->SecDim;
                int blpo=mposi(cell->Tag2);
                for(int i=0;i<mAtom;++i)
                {

                    int rows=blpo+i;
                    int cols=blpo+i;
                    int sendr=rows%procrows;
                    int sendc=cols%proccols;
                    int recvr;
                    int recvc;
                    if (myrow == sendr && mycol == sendc){
                        recvr=rows/procrows;
                        recvc=cols/proccols;
                        dMatrix2(recvr,recvc).real=SubblockData[cell->Tag3](cell->Startp+i);
                    }
                    //ConsideredMatrix(blpo+i,blpo+i)=SubblockData[cell->Tag3](cell->Startp+i);
                }
            }
        }

        blacs_gridexit_(&ctxt);
    }
    void Hamiltonian::gFillDiagonalBlocksByData(Eigen::MatrixXcd &ConsideredMatrix,std::vector<std::unique_ptr<Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd>>> &SubblockData)
    {
        for(auto &cell : systemcell.supercell)
        {
            if(cell->type)
            {
                int pp=cell->Tag2;
                int mAtom=cell->SecDim;
                int blpo=mposi(cell->Tag2);
                ConsideredMatrix.block(blpo,blpo,mAtom,mAtom)=SubblockData[cell->Tag3]->eigenvalues().segment(cell->Startp,mAtom).asDiagonal();
            }
        }
    }

    void Hamiltonian::FillDiagonalBlocksByData(Eigen::MatrixXcd &ConsideredMatrix,std::vector<std::unique_ptr<Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd>>> &SubblockData)
    {
        for(auto &cell : systemcell.supercell)
        {
            if(cell->type)
            {
                int pp=cell->Tag2;
                int mAtom=cell->SecDim;
                int blpo=mposi(cell->Tag2);
                ConsideredMatrix.block(blpo,blpo,mAtom,mAtom)=SubblockData[cell->Tag3]->eigenvalues().segment(cell->Startp,mAtom).asDiagonal();
            }
        }
    }
    void Hamiltonian::BlockCoupling(Eigen::MatrixXcd &ReducedMatrix,std::vector<Eigen::MatrixXcd> &BlockData)
    {
        std::complex<double> ii(0,1);
        Materialblock.clear();
        for(auto &cell : systemcell.ListofMatrixCoupling)
        {
            int p=cell->TagPosition1;
            int pp=cell->TagPosition2;
            int pmposi=mposi(p);
            int ppmposi=mposi(pp);
            std::complex<int> maplist={p,pp};
            std::complex<int> blocklist={cell->ListPosition1,cell->ListPosition2};
            const Eigen::MatrixXcd & lM=BlockData[cell->ListPosition1];
            const Eigen::MatrixXcd & rM=BlockData[cell->ListPosition2];
            std::string FrontStr;
            std::string BackStr;
            std::for_each(
                    systemcell.supercell.front()->UnitTexture.begin(),
                    systemcell.supercell.front()->UnitTexture.end(),
                    [&FrontStr]
                    (auto &scell)
                    {FrontStr+=scell;}
                    );
            std::for_each(
                    systemcell.supercell.back()->UnitTexture.begin(),
                    systemcell.supercell.back()->UnitTexture.end(),
                    [&BackStr]
                    (auto &scell)
                    {BackStr+=scell;}
                    );

            if((Materialblock.count(blocklist)!=1 || (FrontStr!=BackStr && (p==0||pp==0)))||USE_OneRange_or_MoreRanges==0 )
            {
                for(int kk=0;kk<systemcell.supercell[p]->SecDim;++kk)
                {
                    for(int jj=0;jj<systemcell.supercell[pp]->SecDim;++jj)
                    {
                        std::complex<double> elemt=MatrixCoupling11(lM,kk+systemcell.supercell[p]->Startp,jj+systemcell.supercell[pp]->Startp,rM,maplist);
                        ReducedMatrix(kk+pmposi,jj+ppmposi)=elemt;
                        ReducedMatrix(jj+ppmposi,kk+pmposi)=elemt.real()-elemt.imag()*ii;
                    }
                }


                if(abs(pp-p)==1)
                {
                    Materialblock[blocklist]=ReducedMatrix.block(pmposi,ppmposi,systemcell.supercell[p]->SecDim,systemcell.supercell[pp]->SecDim);
                }else{
                    Materialblock[blocklist]=ReducedMatrix.block(pmposi,ppmposi,systemcell.supercell[p]->SecDim,systemcell.supercell[pp]->SecDim).adjoint();
                }

            }else{

                if(abs(pp-p)==1)
                {

                    ReducedMatrix.block(pmposi,ppmposi,systemcell.supercell[p]->SecDim,systemcell.supercell[pp]->SecDim)=Materialblock[blocklist];
                    ReducedMatrix.block(ppmposi,pmposi,systemcell.supercell[pp]->SecDim,systemcell.supercell[p]->SecDim)=Materialblock[blocklist].adjoint();
                }else{

                    ReducedMatrix.block(pmposi,ppmposi,systemcell.supercell[p]->SecDim,systemcell.supercell[pp]->SecDim)=Materialblock[blocklist].adjoint();
                    ReducedMatrix.block(ppmposi,pmposi,systemcell.supercell[pp]->SecDim,systemcell.supercell[p]->SecDim)=Materialblock[blocklist];
                }

            }
        }


    }

    void Hamiltonian::pzheevxBlockCoupling(Eigen::MatrixXcd &ReducedMatrix,std::vector<Eigen::MatrixXcd> &BlockData)
    {

        int id,numprocs;
        int ZEROi = 0;
        int MONEi=-1; 
        int ctxt, myrow, mycol;
        blacs_pinfo_(&id, &numprocs);
        blacs_get_(&MONEi, &ZEROi, &ctxt);
        int procrows = sqrt(numprocs), proccols = sqrt(numprocs);
        blacs_gridinit_(&ctxt, "R", &procrows, &proccols);  
        blacs_gridinfo_(&ctxt, &procrows, &proccols, &myrow, &mycol );  

        std::complex<double> ii(0,1);
        Materialblock2.clear();
        for(auto &cell : systemcell.ListofMatrixCoupling)
        {
            int p=cell->TagPosition1;
            int pp=cell->TagPosition2;
            int pmposi=mposi(p);
            int ppmposi=mposi(pp);
            std::complex<int> maplist={p,pp};
            std::complex<int> blocklist={cell->ListPosition1,cell->ListPosition2};
            const Eigen::MatrixXcd & lM=BlockData[cell->ListPosition1];
            const Eigen::MatrixXcd & rM=BlockData[cell->ListPosition2];
            std::string FrontStr;
            std::string BackStr;
            std::for_each(
                    systemcell.supercell.front()->UnitTexture.begin(),
                    systemcell.supercell.front()->UnitTexture.end(),
                    [&FrontStr]
                    (auto &scell)
                    {FrontStr+=scell;}
                    );
            std::for_each(
                    systemcell.supercell.back()->UnitTexture.begin(),
                    systemcell.supercell.back()->UnitTexture.end(),
                    [&BackStr]
                    (auto &scell)
                    {BackStr+=scell;}
                    );

            Eigen::Matrix<lapack_complex_double,-1,-1> tmpMatrix;
            if(Materialblock2.count(blocklist)!=1 || (FrontStr!=BackStr && (p==0||pp==0)) )
                //if(1)
            {
                tmpMatrix.resize(systemcell.supercell[p]->SecDim,systemcell.supercell[pp]->SecDim);
                tmpMatrix.setConstant({0,0});
                for(int kk=0;kk<systemcell.supercell[p]->SecDim;++kk)
                {
                    for(int jj=0;jj<systemcell.supercell[pp]->SecDim;++jj)
                    {
                        std::complex<double> elemt=MatrixCoupling11(lM,kk+systemcell.supercell[p]->Startp,jj+systemcell.supercell[pp]->Startp,rM,maplist);
                        tmpMatrix(kk,jj).real=real(elemt);
                        tmpMatrix(kk,jj).imag=imag(elemt);
                        //tmpMatrix(jj,kk).real=real(std::conj(elemt));
                        //tmpMatrix(jj,kk).imag=imag(std::conj(elemt));
                        int rows=kk+pmposi;
                        int cols=jj+ppmposi;
                        int sendr=rows%procrows;
                        int sendc=cols%proccols;
                        int recvr=0;
                        int recvc=0;
                        if (myrow == sendr && mycol == sendc){
                            recvr=rows/procrows;
                            recvc=cols/proccols;
                            dMatrix2(recvr,recvc).real=real(elemt);
                            dMatrix2(recvr,recvc).imag=imag(elemt);
                        }
                        if (myrow == sendc && mycol == sendr){
                            recvr=rows/procrows;
                            recvc=cols/proccols;
                            dMatrix2(recvc,recvr).real=real(std::conj(elemt));
                            dMatrix2(recvc,recvr).imag=imag(std::conj(elemt));
                        }

                    }
                }


                if(abs(pp-p)==1)
                {
                    Materialblock2[blocklist]=tmpMatrix;
                }else{
                    Materialblock2[blocklist]=tmpMatrix.adjoint();
                }

            }else{

                if(abs(pp-p)==1)
                {
                    for(int kk=0;kk<systemcell.supercell[p]->SecDim;++kk)
                    {
                        for(int jj=0;jj<systemcell.supercell[pp]->SecDim;++jj)
                        {
                            std::complex<double> elemt=Materialblock2[blocklist](kk,jj).real+Materialblock2[blocklist](kk,jj).imag*ii;
                            int rows=kk+pmposi;
                            int cols=jj+ppmposi;
                            int sendr=rows%procrows;
                            int sendc=cols%proccols;
                            int recvr;
                            int recvc;
                            if (myrow == sendr && mycol == sendc){
                                recvr=rows/procrows;
                                recvc=cols/proccols;
                                dMatrix2(recvr,recvc).real=real(elemt);
                                dMatrix2(recvr,recvc).imag=imag(elemt);
                            }
                            if (myrow == sendc && mycol == sendr){
                                recvr=rows/procrows;
                                recvc=cols/proccols;
                                dMatrix2(recvc,recvr).real=real(std::conj(elemt));
                                dMatrix2(recvc,recvr).imag=imag(std::conj(elemt));
                            }

                        }
                    }

                }else{

                    for(int kk=0;kk<systemcell.supercell[p]->SecDim;++kk)
                    {
                        for(int jj=0;jj<systemcell.supercell[pp]->SecDim;++jj)
                        {
                            std::complex<double> elemt=Materialblock2[blocklist](jj,kk).real-Materialblock2[blocklist](jj,kk).imag*ii;
                            int rows=kk+pmposi;
                            int cols=jj+ppmposi;
                            int sendr=rows%procrows;
                            int sendc=cols%proccols;
                            int recvr;
                            int recvc;
                            if (myrow == sendr && mycol == sendc){
                                recvr=rows/procrows;
                                recvc=cols/proccols;
                                dMatrix2(recvr,recvc).real=real(elemt);
                                dMatrix2(recvr,recvc).imag=imag(elemt);
                            }
                            if (myrow == sendc && mycol == sendr){
                                recvr=rows/procrows;
                                recvc=cols/proccols;
                                dMatrix2(recvc,recvr).real=real(std::conj(elemt));
                                dMatrix2(recvc,recvr).imag=imag(std::conj(elemt));
                            }

                        }
                    }

                }

            }
        }

        blacs_gridexit_(&ctxt);
    }
    void Hamiltonian::gBlockCoupling(Eigen::MatrixXcd &ReducedMatrix,std::vector<std::unique_ptr<Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd>>> &BlockData)
    {
        std::complex<double> ii(0,1);
        Materialblock.clear();
        for(auto &cell : systemcell.ListofMatrixCoupling)
        {
            int p=cell->TagPosition1;
            int pp=cell->TagPosition2;
            int pmposi=mposi(p);
            int ppmposi=mposi(pp);
            std::complex<int> maplist={p,pp};
            std::complex<int> blocklist={cell->ListPosition1,cell->ListPosition2};
            const Eigen::MatrixXcd & lM=BlockData[cell->ListPosition1]->eigenvectors();
            const Eigen::MatrixXcd & rM=BlockData[cell->ListPosition2]->eigenvectors();
            std::string FrontStr;
            std::string BackStr;
            std::for_each(
                    systemcell.supercell.front()->UnitTexture.begin(),
                    systemcell.supercell.front()->UnitTexture.end(),
                    [&FrontStr]
                    (auto &scell)
                    {FrontStr+=scell;}
                    );
            std::for_each(
                    systemcell.supercell.back()->UnitTexture.begin(),
                    systemcell.supercell.back()->UnitTexture.end(),
                    [&BackStr]
                    (auto &scell)
                    {BackStr+=scell;}
                    );

            if((Materialblock.count(blocklist)!=1 || (FrontStr!=BackStr && (p==0||pp==0))) || USE_OneRange_or_MoreRanges==0 )
            {
                for(int kk=0;kk<systemcell.supercell[p]->SecDim;++kk)
                {
                    for(int jj=0;jj<systemcell.supercell[pp]->SecDim;++jj)
                    {
                        std::complex<double> elemt=MatrixCoupling11(lM,kk+systemcell.supercell[p]->Startp,jj+systemcell.supercell[pp]->Startp,rM,maplist);
                        ReducedMatrix(kk+pmposi,jj+ppmposi)=elemt;
                        ReducedMatrix(jj+ppmposi,kk+pmposi)=elemt.real()-elemt.imag()*ii;
                    }
                }


                if(abs(pp-p)==1)
                {
                    Materialblock[blocklist]=ReducedMatrix.block(pmposi,ppmposi,systemcell.supercell[p]->SecDim,systemcell.supercell[pp]->SecDim);
                }else{
                    Materialblock[blocklist]=ReducedMatrix.block(pmposi,ppmposi,systemcell.supercell[p]->SecDim,systemcell.supercell[pp]->SecDim).adjoint();
                }

            }else{

                if(abs(pp-p)==1)
                {

                    ReducedMatrix.block(pmposi,ppmposi,systemcell.supercell[p]->SecDim,systemcell.supercell[pp]->SecDim)=Materialblock[blocklist];
                    ReducedMatrix.block(ppmposi,pmposi,systemcell.supercell[pp]->SecDim,systemcell.supercell[p]->SecDim)=Materialblock[blocklist].adjoint();
                }else{

                    ReducedMatrix.block(pmposi,ppmposi,systemcell.supercell[p]->SecDim,systemcell.supercell[pp]->SecDim)=Materialblock[blocklist].adjoint();
                    ReducedMatrix.block(ppmposi,pmposi,systemcell.supercell[pp]->SecDim,systemcell.supercell[p]->SecDim)=Materialblock[blocklist];
                }

            }
        }

    }

    void Hamiltonian::BlockCoupling(Eigen::MatrixXcd &ReducedMatrix,std::vector<std::unique_ptr<Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd>>> &BlockData)
    {
        std::complex<double> ii(0,1);
        Materialblock.clear();
        for(auto &cell : systemcell.ListofMatrixCoupling)
        {
            int p=cell->TagPosition1;
            int pp=cell->TagPosition2;
            int pmposi=mposi(p);
            int ppmposi=mposi(pp);
            std::complex<int> maplist={p,pp};
            std::complex<int> blocklist={cell->ListPosition1,cell->ListPosition2};
            const Eigen::MatrixXcd & lM=BlockData[cell->ListPosition1]->eigenvectors();
            const Eigen::MatrixXcd & rM=BlockData[cell->ListPosition2]->eigenvectors();
            std::string FrontStr;
            std::string BackStr;
            std::for_each(
                    systemcell.supercell.front()->UnitTexture.begin(),
                    systemcell.supercell.front()->UnitTexture.end(),
                    [&FrontStr]
                    (auto &scell)
                    {FrontStr+=scell;}
                    );
            std::for_each(
                    systemcell.supercell.back()->UnitTexture.begin(),
                    systemcell.supercell.back()->UnitTexture.end(),
                    [&BackStr]
                    (auto &scell)
                    {BackStr+=scell;}
                    );

            if((Materialblock.count(blocklist)!=1 || (FrontStr!=BackStr && (p==0||pp==0)))||USE_OneRange_or_MoreRanges==0 )
            {
                for(int kk=0;kk<systemcell.supercell[p]->SecDim;++kk)
                {
                    for(int jj=0;jj<systemcell.supercell[pp]->SecDim;++jj)
                    {
                        std::complex<double> elemt=MatrixCoupling11(lM,kk+systemcell.supercell[p]->Startp,jj+systemcell.supercell[pp]->Startp,rM,maplist);
                        ReducedMatrix(kk+pmposi,jj+ppmposi)=elemt;
                        ReducedMatrix(jj+ppmposi,kk+pmposi)=elemt.real()-elemt.imag()*ii;
                    }
                }


                if(abs(pp-p)==1)
                {
                    Materialblock[blocklist]=ReducedMatrix.block(pmposi,ppmposi,systemcell.supercell[p]->SecDim,systemcell.supercell[pp]->SecDim);
                }else{
                    Materialblock[blocklist]=ReducedMatrix.block(pmposi,ppmposi,systemcell.supercell[p]->SecDim,systemcell.supercell[pp]->SecDim).adjoint();
                }

            }else{

                if(abs(pp-p)==1)
                {

                    ReducedMatrix.block(pmposi,ppmposi,systemcell.supercell[p]->SecDim,systemcell.supercell[pp]->SecDim)=Materialblock[blocklist];
                    ReducedMatrix.block(ppmposi,pmposi,systemcell.supercell[pp]->SecDim,systemcell.supercell[p]->SecDim)=Materialblock[blocklist].adjoint();
                }else{

                    ReducedMatrix.block(pmposi,ppmposi,systemcell.supercell[p]->SecDim,systemcell.supercell[pp]->SecDim)=Materialblock[blocklist].adjoint();
                    ReducedMatrix.block(ppmposi,pmposi,systemcell.supercell[pp]->SecDim,systemcell.supercell[p]->SecDim)=Materialblock[blocklist];
                }

            }
        }

    }

    void Hamiltonian::Append_matrix()
    {
        int KextraSize=0; 
        for(int i=0;i<KextraTransform.size();++i)
        {
            KextraSize+=KextraTransform[i].cols();
        }
        matrix.conservativeResize(MatrixSizeReduced+KextraSize,MatrixSizeReduced+KextraSize);
        matrix.block(MatrixSizeReduced,MatrixSizeReduced,KextraSize,KextraSize)=extendMat;
        matrix.block(MatrixSizeReduced,0,KextraSize,MatrixSizeReduced)=extendSquareMat;
        //matrix.block(0,MatrixSizeReduced,MatrixSizeReduced,transform_X.cols()+transform_L.cols())=extendSquareMat.adjoint();
        MatrixSizeReduced=matrix.rows();
        /*
           matrix.conservativeResize(MatrixSizeReduced+transform_X.cols()+transform_L.cols(),MatrixSizeReduced+transform_X.cols()+transform_L.cols());
           matrix.block(MatrixSizeReduced,MatrixSizeReduced,transform_X.cols()+transform_L.cols(),transform_X.cols()+transform_L.cols())=extendMat;
           matrix.block(MatrixSizeReduced,0,transform_X.cols()+transform_L.cols(),MatrixSizeReduced)=extendSquareMat;
        //matrix.block(0,MatrixSizeReduced,MatrixSizeReduced,transform_X.cols()+transform_L.cols())=extendSquareMat.adjoint();
        MatrixSizeReduced=matrix.rows();
        */
    }

    void Hamiltonian::gMatrixDiagonalize(bool onoff,Eigen::MatrixXcd OverlapMatrix)
    {
        int MR= OverlapMatrix.rows();
        if(onoff)
        {
            solver.compute(matrix.block(0,0,MR,MR),Eigen::DecompositionOptions::ComputeEigenvectors);

        }else{
            //gsolver.compute(matrix.block(0,0,MR,MR),OverlapMatrix, Eigen::DecompositionOptions::EigenvaluesOnly|Eigen::Ax_lBx);
            solver.compute(matrix.block(0,0,MR,MR),Eigen::DecompositionOptions::EigenvaluesOnly);
            //solver.compute(matrix.block(0,0,MR,MR), Eigen::DecompositionOptions::ComputeEigenvectors);
        }
        //gsolver.compute(matrix.block(0,0,MR,MR),OverlapMatrix,Eigen::DecompositionOptions::ComputeEigenvectors|Eigen::Ax_lBx);
        assert(solver.info() == Eigen::ComputationInfo::Success);
        /*
           for(int cs=0;cs<10;++cs)
           {
           std::cout<<"eigen["<<cs<<"]="<<solver.eigenvalues()(8+cs-6)<<std::endl;
           }
           */
        //std::cout<<"eigen="<<solver.eigenvalues()<<std::endl;
        if(onoff)
        {
            //WaveFunctions();
            //printOutWaveFunctions();
        }
    }

    void Hamiltonian::MatrixDiagonalize(bool onoff,Eigen::MatrixXcd OverlapMatrix)
    {
        int MR= OverlapMatrix.rows();
        if(onoff)
        {
            gsolver.compute(matrix.block(0,0,MR,MR),OverlapMatrix, Eigen::DecompositionOptions::ComputeEigenvectors|Eigen::Ax_lBx);

        }else{
            gsolver.compute(matrix.block(0,0,MR,MR),OverlapMatrix, Eigen::DecompositionOptions::EigenvaluesOnly|Eigen::Ax_lBx);
            //gsolver.compute(matrix.block(0,0,MR,MR),OverlapMatrix, Eigen::DecompositionOptions::ComputeEigenvectors|Eigen::Ax_lBx);
        }
        //gsolver.compute(matrix.block(0,0,MR,MR),OverlapMatrix,Eigen::DecompositionOptions::ComputeEigenvectors|Eigen::Ax_lBx);
        assert(gsolver.info() == Eigen::ComputationInfo::Success);
        if(onoff)
        {
            WaveFunctions();
            printOutWaveFunctions();
        }
    }




    void Hamiltonian::print_Hamiltonian()
    {

        Eigen::Matrix<std::complex<double>, -1, -1> ffmatrix;
        int Atom=0;
        for(int i=0;i<systemcell.supercell.size();i++)
        {
            Atom+=systemcell.supercell[i]->AtomNum;
        }
        ffmatrix.resize(Atom,Atom);	
        ffmatrix.setZero();
        FILE *fp;
        fp=fopen("./Pattern.data","w+");
        fprintf(fp, "\n matrix\n");
        for(int i=0;i<ffmatrix.rows();i++)
        {
            for(int j=0;j<ffmatrix.cols();j++)
            {
                int flag=in_the_list(i,j);
                std::string Sflag=Stringin_the_list(i,j);
                if(flag)
                {
                    //fprintf(fp,"{%2.0f,%2.0f}",i*1.,j*1.);
                    //fprintf(fp,"{%d}",flag);
                    fprintf(fp,"{%s}",Sflag.c_str());
                }else{
                    //fprintf(fp,"{     }");
                    fprintf(fp,"{  }");
                    //fprintf(fp,"{ }");
                }

            }fprintf(fp,"\n");
        }
        fclose(fp);
        ffmatrix.setZero();
    }


    std::string Hamiltonian::Stringin_the_list(int ii,int j)
    {
        //for(auto &cell : systemcell.ListofMatrixElements_Type2)
        for(auto &cell : systemcell.ListofMatrixElements_Type1)
            //for(auto &cell : systemcell.ListofMatrixElements_Type3)
            //if(ii==cell->Table.real() && j==cell->Table.imag())
            //for(auto &cell : list)
            //for(auto &cell : systemcell.ListofMatrixCoupling)
            //if(ii*systemcell.orbNum==cell->Table.real() && j*systemcell.orbNum==cell->Table.imag())
        {
            if(ii==cell->Table.real() && j==cell->Table.imag())
            {
                return cell->cellme.substr(0,1)+cell->SuperCell.substr(0,1);
                //return cell->Typeme+1;
                //return "1111";
            }
        }
        return "0000";
    }

    int Hamiltonian::in_the_list(int ii,int j)
    {
        //for(auto &cell : systemcell.ListofMatrixElements_Type2)
        for(auto &cell : systemcell.ListofMatrixElements_Type1)
            //for(auto &cell : systemcell.ListofMatrixElements_Type3)
            //if(ii==cell->Table.real() && j==cell->Table.imag())
            //for(auto &cell : list)
            //for(auto &cell : systemcell.ListofMatrixCoupling)
        {
            //if(ii*systemcell.orbNum==cell->Table.real() && j*systemcell.orbNum==cell->Table.imag())
            if(ii==cell->Table.real() && j==cell->Table.imag())
            {
                return cell->expflag+1;
                //return cell->Typeme+1;
                //return 1;
            }
        }
        return 0;
    }

    void Hamiltonian::printSupercell(std::string name)
    {
        int NUM=0;
        printf("==============%s============\n",name.c_str());
        for(auto &cell : systemcell.supercell)
        {
            printf("Tag3(%d)Tag2[(%d)]tag[%d]unit%d->%d \t",cell->Tag3,cell->Tag2,cell->tag,cell->AtomNum,cell->AtomNum/systemcell.UnitCellNumber);
            for(auto &utexture : cell->UnitTexture)
            {
                printf("%s ",utexture.c_str());
            }printf("\n");
            NUM+=cell->AtomNum;
        }printf("NumTotal=%d\n",NUM);
    }


    void Hamiltonian::print_rmatrix(Eigen::MatrixXd mat,std::string name)
    {
        fprintf(stderr, "\n %s %ld %ld\n",name.c_str(),mat.rows(),mat.cols());

        for(int i=0;i<mat.rows();i++)
        {
            for(int j=0;j<mat.cols();j++)
            {
                if(fabs(mat(i,j))<1e-5)
                {
                    fprintf(stderr,"{     }");
                    //fprintf(stderr,"{  }");
                }else{
                    fprintf(stderr,"{%2.3f}",mat(i,j));
                    //fprintf(stderr,"{%2.0f}",mat(i,j));
                }
            }fprintf(stderr,"\n");
        }
    }


    void Hamiltonian::print_matrix(Eigen::MatrixXcd mat,std::string name)
    {
        std::complex<double> trace=mat.trace();
        //fprintf(stderr, "\n %s %lutrace[%1.2f,%1.2f]\n",name.c_str(),mat.cols(),trace.real(),trace.imag());
        fprintf(stderr, "\n %s %ld %ld\n",name.c_str(),mat.rows(),mat.cols());

        for(int i=0;i<mat.rows();i++)
        {
            for(int j=0;j<mat.cols();j++)
            {
                //if(fabs(mat(i,j).imag()-0)<1e-10 && fabs(mat(i,j).real()-0)<1e-10)
                if(fabs(mat(i,j).imag()-0)<1e-3 && fabs(mat(i,j).real()-0)<1e-3)
                {
                    //fprintf(stderr,"{   }");
                    //fprintf(stderr,"{      }");
                    //fprintf(stderr,"{        }");
                    //fprintf(stderr,"{  }");
                    fprintf(stderr,"{ }");
                    //fprintf(stderr,"   ");
                }else{
                    //fprintf(stderr,"{%2.0f}",mat(i,j).real());
                    //fprintf(stderr,"{%2.3f}",mat(i,j).real());
                    //fprintf(stderr,"{%2.3f/%2.3f}",fabs(mat(i,j).real()),fabs(mat(i,j).imag()));
                    //fprintf(stderr,"{%2.1f/%2.1f}",(mat(i,j).real()),(mat(i,j).imag()));
                    //fprintf(stderr,"{%1.0f}",abs(mat(i,j).real()));
                    fprintf(stderr,"{%1.0f}",mat(i,j).real());
                    //fprintf(stderr," + ");
                }
            }fprintf(stderr,"\n");
        }

    }
    /*
       int Hamiltonian::FoldTheCol(std::vector<double> &egv,std::vector<Eigen::MatrixXcd> &egvct,int cellme,int blpo,int AtomO,std::vector<int> vblpo,std::vector<int> vAtomO,const Eigen::Vector3d& k,std::vector<Eigen::Vector3d> chain,std::pair<double,double> EnergyCutoff,int newcolnum,int colnum,Eigen::Matrix3d &Ugm,Eigen::Matrix3d &ugm,Eigen::Matrix3d &uam,Eigen::MatrixXcd &yfitting,int DataTag)
       {
       int retvalue=0;
       int UnitCellNumber=systemcell.UnitCellNumber;
       Eigen::Vector3d shift=systemcell.supercell[cellme]->Unit[0];
       assert(vblpo.size()==vAtomO.size());
       int blnum=vblpo.size();

       Eigen::VectorXcd retVxcd;
       retVxcd.resize(AtomO);
       retVxcd.setZero();
       int fitposCol=0;
       std::complex<double> ii(0.0,1.0);
       for(int i=0;i<blnum;++i)
       {
       int fitposRow=0;
       int bigtrans=i/systemcell.uvector.size();
       Eigen::Vector3d pshift=Eigen::Vector3d{bigtrans,0,0};
       Eigen::Vector3d inv_pshift=Eigen::Vector3d{bigtrans==0?0:(1.0/bigtrans),0,0};
       for(int j=0;j<blnum;++j)
       {
       Eigen::Vector3d ddd=chain[(i)*systemcell.UnitCellNumber]-shift-pshift;
       Eigen::Vector3d AtomPos=chain[(i)*systemcell.UnitCellNumber]-shift;
       Eigen::Vector3d dk(0,0,0);

       Eigen::Vector3d projAtomPos=uam*AtomPos;
       Eigen::Vector3d poshift((int)projAtomPos(0),(int)projAtomPos(1),(int)projAtomPos(2));
       int LpV,minDistance;
       std::vector<int> pushV2(systemcell.pushVec.begin()+1,systemcell.pushVec.end()-1);
       int bpushV2=pushV2[0];
       for(auto &pv : pushV2)
       {
       pv=pv-bpushV2;
       }

       for(int ini=0;ini<pushV2.size();++ini)
       {
       if(pushV2[ini]==poshift(0))
       {
       minDistance=ini;
       LpV=pushV2[minDistance];
       break;
       }
       if(pushV2[ini]>poshift(0))
       {
       minDistance=ini-1;
       LpV=pushV2[minDistance];
       break;
       }
       }
       int cNUM=80;
    //std::cout<<"cNum="<<systemcell.uvector.size()*systemcell.UnitCellNumber*orbNum<<std::endl;
    cNUM=systemcell.uvector.size()*systemcell.UnitCellNumber*orbNum;
    dk=ugm*(AtomPos-uam.inverse()*poshift)+Ugm*poshift;//dk1
    for(int cc=0;cc<UnitCellNumber;++cc)
    {
    for(int kj=0;kj<orbNum;++kj)
    {
    int colIndex=cc*orbNum+kj;

    for(int c=0;c<UnitCellNumber;++c)
    {
    Eigen::Vector3d posc=chain[(j)*UnitCellNumber+c]-shift;
    std::complex<double> expr=exp(((dk).dot(posc))*2*M_PI*ii);
    for(int ki=0;ki<orbNum;++ki)
    {
    int rowIndex=c*orbNum+ki;

    //fitting(fitposRow+rowIndex,fitposCol+colIndex)=egvct[i](rowIndex,colIndex)*expr/sqrt(blnum);//5
    if(
    (fitposCol+colIndex)==colnum
        && USE_RowTruncation==0
        )
        {
            //if(pshift(0)==0)
            {
                if(
                        (fitposCol+colIndex)==colnum
                        //&& (fitposRow+rowIndex)>=(LpV)*cNUM
                        //&& (fitposRow+rowIndex)<(LpV+systemcell.spNum[minDistance+1])*cNUM
                  )
                {
                    yfitting(fitposRow+rowIndex,newcolnum)=egvct[i](rowIndex,colIndex)*expr/sqrt(blnum);//5
                    retvalue=1;
                }
            }
        }else{

            if(
                    (fitposCol+colIndex)==colnum
                    && (fitposRow+rowIndex)>=(LpV)*cNUM
                    && (fitposRow+rowIndex)<(LpV+systemcell.spNum[minDistance+1])*cNUM
                    //&& (fitposRow+rowIndex)>=poshift(0)*cNUM
                    //&& (fitposRow+rowIndex)<(poshift(0)+1)*cNUM
              )
            {
                //retVxcd(fitposRow+rowIndex)=egvct[i](rowIndex,colIndex)*expr/sqrt(blnum);//5
                yfitting(fitposRow+rowIndex,newcolnum)=egvct[i](rowIndex,colIndex)*expr/sqrt(blnum);//5
                retvalue=1;
            }
        }

    }
    }

    }
    }
    fitposRow+=vAtomO[j];
    }
    fitposCol+=vAtomO[i];
    }


    return retvalue;
    }
    */
        /* 
           void Hamiltonian::FoldTheBlockHamiltonian2(int cellme,int blpo,int AtomO,std::vector<int> vblpo,std::vector<int> vAtomO,const Eigen::Vector3d& k,std::vector<Eigen::Vector3d> chain,Eigen::MatrixXcd &fitting)//without dtag
           {

           int UnitCellNumber=systemcell.UnitCellNumber;
           int NNN=systemcell.supercell[cellme]->AtomNum/(systemcell.UnitCellNumber*systemcell.uvector.size());

           systemcell.Unita1=systemcell.a1*systemcell.supercell[cellme]->AtomNum/(systemcell.UnitCellNumber*systemcell.uvector.size());


           systemcell.Unita2=systemcell.a2;
           systemcell.Unita3=systemcell.a3;
           double SVolumn=systemcell.Unita1.dot((systemcell.Unita2.cross(systemcell.Unita3)));
           systemcell.Ug1=systemcell.Unita2.cross(systemcell.Unita3)/SVolumn;
           systemcell.Ug2=systemcell.Unita3.cross(systemcell.Unita1)/SVolumn;
           systemcell.Ug3=systemcell.Unita1.cross(systemcell.Unita2)/SVolumn;

           SVolumn=systemcell.a1.dot((systemcell.a2.cross(systemcell.a3)));
           Eigen::Vector3d ua1=systemcell.a2.cross(systemcell.a3)/SVolumn;
           Eigen::Vector3d ua2=systemcell.a3.cross(systemcell.a1)/SVolumn;
           Eigen::Vector3d ua3=systemcell.a1.cross(systemcell.a2)/SVolumn;

        //Eigen::VectorXd m;
        Eigen::Matrix3d Ugm,ugm,uam;
        uam<<ua1,ua2,ua3;
        Ugm<<systemcell.Ug1,systemcell.Ug2,systemcell.Ug3;
        ugm<<systemcell.ug1,systemcell.ug2,systemcell.ug3;
        Eigen::Matrix3d uga,Uga;
        uga<<systemcell.u1,systemcell.u2,systemcell.u3;
        Uga<<systemcell.Unita1,systemcell.Unita2,systemcell.Unita3;
        //Eigen::Matrix3d A=Eigen::Map<Eigen::Matrix3d>(m.data());

        int printM=cellme;
        Eigen::Vector3d shift=systemcell.supercell[cellme]->Unit[0];
        assert(vblpo.size()==vAtomO.size());
        //std::vector<double> egv,egvTMP;
        Eigen::VectorXd egv;
        std::vector<Eigen::MatrixXcd> egvct;
        egvct.resize(0);
        int blnum=vblpo.size();
        //Eigen::MatrixXcd fitting;

        fitting.resize(AtomO,AtomO);
        fitting.setZero();

        std::complex<double> ii(0,1);
        int fitposCol=0;

        egv.resize(AtomO);
        egv.setZero();

        for(int i=0;i<blnum;++i)
        {
        int blen=vAtomO[i];
        gsolver.compute(matrixScatter.block(vblpo[i],vblpo[i],blen,blen),Eigen::MatrixXcd::Identity(blen,blen) ,Eigen::DecompositionOptions::ComputeEigenvectors|Eigen::Ax_lBx);
        egv.segment(vblpo[i]-vblpo[0],blen)=gsolver.eigenvalues();
        //egv.insert(egv.end(),gsolver.eigenvalues().data(),gsolver.eigenvalues().data()+gsolver.eigenvalues().size());
        egvct.push_back(gsolver.eigenvectors().transpose());
        int fitposRow=0;
        int bigtrans=i/systemcell.uvector.size();
        Eigen::Vector3d pshift=Eigen::Vector3d{bigtrans,0,0};
        Eigen::Vector3d inv_pshift=Eigen::Vector3d{bigtrans==0?0:(1.0/bigtrans),0,0};
        for(int j=0;j<blnum;++j)
        {
        Eigen::Vector3d ddd=chain[(i)*systemcell.UnitCellNumber]-shift-pshift;
        Eigen::Vector3d AtomPos=chain[(i)*systemcell.UnitCellNumber]-shift;
        for(int c=0;c<UnitCellNumber;++c)
        {
        Eigen::Vector3d dk(0,0,0);
        int IND=(j)*UnitCellNumber+c;
        Eigen::Vector3d posc=chain[(j)*UnitCellNumber+c]-shift;

        Eigen::Vector3d projAtomPos=uam*AtomPos;
    Eigen::Vector3d poshift((int)projAtomPos(0),(int)projAtomPos(1),(int)projAtomPos(2));

    dk=ugm*(AtomPos-uam.inverse()*poshift)+Ugm*poshift;//dk1
    //if((AtomPos-poshift).cross(systemcell.u1)==Eigen::Vector3d(0,0,0))

    std::complex<double> expr=exp(((dk).dot(posc))*2*M_PI*ii);

    for(int ki=0;ki<orbNum;++ki)
    {
        int rowIndex=c*orbNum+ki;
        fitting.block(fitposRow,fitposCol,vAtomO[j],vAtomO[i]).row(rowIndex)=gsolver.eigenvectors().row(rowIndex)*expr/sqrt(blnum);//5
    }
    }
    fitposRow+=vAtomO[j];
    }
    fitposCol+=vAtomO[i];
    }


    for(int i=0;i<egv.size();++i){
        for(int j=0;j<egv.size()-i-1;++j){
            //if(ord[j].imag()>ord[j+1].imag())
            if(egv(j)>egv(j+1))
            {
                double tmp=egv(j);
                egv(j)=egv(j+1);
                egv(j+1)=tmp;
                //std::complex<double> tmp2=ord[j];
                //ord[j]=ord[j+1];
                //ord[j+1]=tmp2;
                fitting.col(j).swap(fitting.col(j+1));
            }
        }
    }

    //return Startp;
    }
    */
        /*
           void Hamiltonian::FoldTheBlockHamiltonian(int cellme,int blpo,int AtomO,std::vector<int> vblpo,std::vector<int> vAtomO,const Eigen::Vector3d& k,std::vector<Eigen::Vector3d> chain,std::pair<double,double> EnergyCutoff,int &Startp,int &SecDim,Eigen::MatrixXcd &sfitting)
           {

           int UnitCellNumber=systemcell.UnitCellNumber;
           int NNN=systemcell.supercell[cellme]->AtomNum/(systemcell.UnitCellNumber*systemcell.uvector.size());

           systemcell.Unita1=systemcell.a1*systemcell.supercell[cellme]->AtomNum/(systemcell.UnitCellNumber*systemcell.uvector.size());


           systemcell.Unita2=systemcell.a2;
           systemcell.Unita3=systemcell.a3;
           double SVolumn=systemcell.Unita1.dot((systemcell.Unita2.cross(systemcell.Unita3)));
           systemcell.Ug1=systemcell.Unita2.cross(systemcell.Unita3)/SVolumn;
           systemcell.Ug2=systemcell.Unita3.cross(systemcell.Unita1)/SVolumn;
           systemcell.Ug3=systemcell.Unita1.cross(systemcell.Unita2)/SVolumn;

           SVolumn=systemcell.a1.dot((systemcell.a2.cross(systemcell.a3)));
           Eigen::Vector3d ua1=systemcell.a2.cross(systemcell.a3)/SVolumn;
           Eigen::Vector3d ua2=systemcell.a3.cross(systemcell.a1)/SVolumn;
           Eigen::Vector3d ua3=systemcell.a1.cross(systemcell.a2)/SVolumn;

        //Eigen::VectorXd m;
        Eigen::Matrix3d Ugm,ugm,uam;
        uam<<ua1,ua2,ua3;
        Ugm<<systemcell.Ug1,systemcell.Ug2,systemcell.Ug3;
        ugm<<systemcell.ug1,systemcell.ug2,systemcell.ug3;
        Eigen::Matrix3d uga,Uga;
        uga<<systemcell.u1,systemcell.u2,systemcell.u3;
        Uga<<systemcell.Unita1,systemcell.Unita2,systemcell.Unita3;
        //Eigen::Matrix3d A=Eigen::Map<Eigen::Matrix3d>(m.data());

        int printM=cellme;
        Eigen::Vector3d shift=systemcell.supercell[cellme]->Unit[0];
        assert(vblpo.size()==vAtomO.size());
        std::vector<double> egv,egvTMP;
        std::vector<Eigen::MatrixXcd> egvct;
        egv.resize(0);
        egvct.resize(0);
        int blnum=vblpo.size();
        Eigen::MatrixXcd fitting;
        fitting.resize(AtomO,AtomO);
        fitting.setZero();
        std::complex<double> ii(0,1);
        int fitposCol=0;
        for(int i=0;i<blnum;++i)
        {
        int blen=vAtomO[i];
        gsolver.compute(matrixScatter.block(vblpo[i],vblpo[i],blen,blen),Eigen::MatrixXcd::Identity(blen,blen) ,Eigen::DecompositionOptions::ComputeEigenvectors|Eigen::Ax_lBx);
        egv.insert(egv.end(),gsolver.eigenvalues().data(),gsolver.eigenvalues().data()+gsolver.eigenvalues().size());
        egvct.push_back(gsolver.eigenvectors().transpose());
        int fitposRow=0;
        int bigtrans=i/systemcell.uvector.size();
        Eigen::Vector3d pshift=Eigen::Vector3d{bigtrans,0,0};
        Eigen::Vector3d inv_pshift=Eigen::Vector3d{bigtrans==0?0:(1.0/bigtrans),0,0};
        for(int j=0;j<blnum;++j)
        {
        Eigen::Vector3d ddd=chain[(i)*systemcell.UnitCellNumber]-shift-pshift;
        Eigen::Vector3d AtomPos=chain[(i)*systemcell.UnitCellNumber]-shift;
        for(int c=0;c<UnitCellNumber;++c)
        {
        Eigen::Vector3d dk(0,0,0);
        int IND=(j)*UnitCellNumber+c;
        Eigen::Vector3d posc=chain[(j)*UnitCellNumber+c]-shift;

        Eigen::Vector3d projAtomPos=uam*AtomPos;
        Eigen::Vector3d poshift((int)projAtomPos(0),(int)projAtomPos(1),(int)projAtomPos(2));

        dk=ugm*(AtomPos-uam.inverse()*poshift)+Ugm*poshift;//dk1
        //if((AtomPos-poshift).cross(systemcell.u1)==Eigen::Vector3d(0,0,0))
        //if(AtomPos(1)==0.5 &&AtomPos(2)==0.5)
        //{
        //dk=Ugm*ugm*(AtomPos+(NNN-1)*(AtomPos-pshift));//dk2
        //}else{
        //dk=Ugm*ugm*(AtomPos);//dk2
        //}
        std::complex<double> expr=exp(((dk).dot(posc))*2*M_PI*ii);

    for(int ki=0;ki<orbNum;++ki)
    {
        int rowIndex=c*orbNum+ki;
        fitting.block(fitposRow,fitposCol,vAtomO[j],vAtomO[i]).row(rowIndex)=gsolver.eigenvectors().row(rowIndex)*expr;//5
    }
    }
    fitposRow+=vAtomO[j];
    }
    fitposCol+=vAtomO[i];
    }

    for(int i=0;i<fitting.cols();++i)
    {
        fitting.col(i)=fitting.col(i)/sqrt(fitting.col(i).squaredNorm());
    }//NORMALIZE

    std::vector<std::complex<double>> ord;
    sfitting=fitting;
    for(int i=0;i<egv.size();++i){
        ord.push_back({i*1.,egv[i]});
    }
    for(int i=0;i<egv.size();++i){
        for(int j=0;j<egv.size()-i-1;++j){
            if(ord[j].imag()>ord[j+1].imag()){
                //double tmp=egv[j];
                //egv[j]=egv[j+1];
                //egv[j+1]=tmp;
                std::complex<double> tmp2=ord[j];
                ord[j]=ord[j+1];
                ord[j+1]=tmp2;
                sfitting.col(j).swap(sfitting.col(j+1));
            }
        }
    }




    Eigen::MatrixXcd yfitting;
    if(USE_OneRange_or_MoreRanges)
    {
        if(USE_FixedSize_or_EnergyRange)
        {
            Startp=0;
            yfitting.resize(fitting.rows(),leftCutoff+rightCutoff);
            yfitting.setZero();
            for(int yi=0;yi<yfitting.cols();++yi)
            {
                yfitting.col(yi)=fitting.col(ord[AtomO*BandGapIndex/orbNum+yi-leftCutoff].real());
            }
            Startp=AtomO*BandGapIndex/orbNum-leftCutoff;
            SecDim=rightCutoff+leftCutoff;
        }else{
            yfitting.resize(fitting.rows(),fitting.cols());
            yfitting.setZero();
            int flag=0,yflag=0;
            Startp=0;
            for(int i=0;i<egv.size();++i)
            {
                if(egv[i]<leftCutoff) Startp+=1;
                if(egv[i]<leftCutoff || egv[i]>rightCutoff)
                {
                    egv.erase(egv.begin()+i);
                    i-=1;
                    yflag+=1;
                }else{
                    yfitting.col(flag)=fitting.col(yflag+i);
                    flag+=1;
                }
            }
            yfitting.conservativeResize(AtomO,flag);
            SecDim=yfitting.cols();
        }
        //gsoulFold_matrix.push_back(yfitting);
    }else{

        if(USE_FixedSize_or_EnergyRange)
        {
            Startp=AtomO*BandGapIndex/orbNum-EnergyCutoff.first;
            SecDim=EnergyCutoff.first+EnergyCutoff.second;
        }else{
            int flag=0,yflag=0;
            Startp=0;
            for(int i=0;i<egv.size();++i)
            {
                if(egv[i]<EnergyCutoff.first) Startp+=1;
                if(egv[i]<EnergyCutoff.first || egv[i]>EnergyCutoff.second)
                {
                    egv.erase(egv.begin()+i);
                    i-=1;
                    yflag+=1;
                }else{
                    flag+=1;
                }
            }
            SecDim=flag;
        }

        yfitting.resize(fitting.rows(),fitting.cols());
        yfitting.setZero();
        for(int yi=0;yi<yfitting.cols();++yi)
        {
            yfitting.col(yi)=fitting.col(ord[yi].real());
        }

        //gsoulFold_matrix.push_back(yfitting);
    }

    //return Startp;
    }
    */



        void Hamiltonian::update_pD()
        {
            mappD.clear();
            //for(auto &cell : tpfast)
            for(auto &cell : systemcell.ListofMatrixElements_Type1)
            {
                Eigen::Vector3d dn=cell->d.normalized();
                if(mappD.count(dn(2))!=1)
                {
                    update_params(dn);
                    mappD[dn(2)]=pD;
                }else{
                    //update_params(dn);
                    //printf("%d\n",mappD[dn(2)]==pD);
                }
            }
        }
    void Hamiltonian::update_params(Eigen::Vector3d dn)
    {

        pD.resize(3);
        for(int l=0;l<3;++l)
        {
            pD[l].resize(0);
            for(int m=0;m<=l;++m)
            {
                for(int pm=0;pm<=l;++pm)
                {
                    pD[l].push_back(std::complex<double>(func_D(l,m,pm,dn),func_D(l,m,-pm,dn)));
                }
            }
        }
    }
    void Hamiltonian::SetMatrix5(const Eigen::Vector3d kk,Eigen::Vector3d pl,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::MatrixXcd &matrix,int Wsplit)
    {
        std::unordered_map<std::string, Material> &ms=material.materials;
        if(USE_long_parameters==0)
        {
        ms=material.materialsShort;
        }
        int Atom=MatrixSize;
        matrix.resize(Atom,Atom);	
        matrix.setZero();
        for(auto &cell : ListofMatrixElements)
        {
            std::complex<double> ii(0.0,1.0);
            Eigen::Vector3d dd=cell->d;
            Eigen::Vector3d ed=cell->ed;
            Eigen::Vector3d edn=ed.normalized();
            Eigen::Vector3d dn=dd.normalized();
            //std::complex<double> efactor=exp(2*M_PI*((dd).dot(k))*ii);

            //if(dd.dot(systemcell.Unita1)>0) dd(2)=0;
            //dd(0)=0;
            std::complex<double> efactor=0;

            int atom2=cell->Table.real();
            int atom1=cell->Table.imag();
            std::string leftType=(atom1%2==0?"c":"a");
            //std::string leftType="a";
            //UnitType[atom2];

            std::string rightType=(atom2%2==0?"c":"a");
            //std::string rightType="c";
            atom2*=orbNum;
            atom1*=orbNum;
            //UnitType[atom1];
            

            if(Wsplit==999) efactor=exp(2*M_PI*((dd).dot(kk))*ii);
            else if(Wsplit==0 && dd.dot(pl)==0) efactor=exp(2*M_PI*((dd).dot(kk))*ii);
            else if(Wsplit==1 && dd.dot(pl)>0) efactor=exp(2*M_PI*((dd).dot(kk-kk.dot(pl)*pl))*ii);
            else if(Wsplit==-1 && dd.dot(pl)<0) efactor=exp(2*M_PI*((dd).dot(kk-kk.dot(pl)*pl))*ii);
            else if(Wsplit==-999) efactor=exp(2*M_PI*(dd.dot(kk-kk.dot(pl)*pl))*ii);

            for(int io=0;io<orbNum;++io)
            {
                for(int jo=0;jo<orbNum;++jo)
                {
                    //std::string left=labels[jo];
                    //std::string right=labels[io];
                    Eigen::Vector4i leftorbital=qlabels[jo];
                    Eigen::Vector4i rightorbital=qlabels[io];
                    std::string QWMaterial="GaAs";
                    std::string EBMaterial="AlAs";

                    if(leftorbital(3)==rightorbital(3))
                    {
                        int rows=atom2+io;
                        int cols=atom1+jo;
                        if(rows==cols && (abs(Wsplit)==999 || Wsplit==0))
                        {
                            if(cell->cellme==QWMaterial && cell->SuperCell==QWMaterial)
                            {
                                matrix(rows,cols)=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType)+BandOffSet;
                            }else{
                                matrix(rows,cols)=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType);
                            }
                        }
                        else if(leftType!=rightType)
                        { 
                            //double res=(SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->cellme))*0.5+SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->SuperCell))*0.5);
                            double res=0;
                            if(rows>cols)
                            {
                                res=(SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->cellme)));
                            }else if(rows<cols)
                            {
                                res=(SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->SuperCell)));
                            }
                            matrix(rows,cols)+=res*efactor;
                        }
                        if(
                                atom2==atom1 
                                && leftorbital(1)==rightorbital(1) 
                                && leftorbital(1)==1
                                &&USE_SOC
                                && 0
                          )
                        {
                            if((leftorbital(2)+rightorbital(2))==0)
                            {
                                double lambda=ms.at(cell->cellme).param.at(leftType);
                                matrix(rows,cols)+=ii*(0.5*rightorbital(2))*lambda;
                            }
                        }
                    }else{
                        int rows=atom2+io;
                        int cols=atom1+jo;

                        double lambda=1;
                        double m1=leftorbital(2);
                        double m2=rightorbital(2);
                        double leftsign=leftorbital(3);
                        double rightsign=rightorbital(3);
                        if(
                                atom2==atom1 
                                && leftorbital(1)==rightorbital(1) 
                                && leftorbital(1)==1
                                &&USE_SOC
                                && 1
                          )
                        {

                            assert(leftType==rightType);
                            lambda=ms.at(cell->cellme).param.at(leftType);
                            std::complex<double> soc;
                            std::complex<double> theta1=(1-(0==m1?1:0))*1.0/sqrt(2)*((m1>=0?1.0:0.0)-ii*((-m1)>=0?1.0:0.0))+0.5*(0==m1?1:0);
                            std::complex<double> theta2=(1-(0==m2?1:0))*1.0/sqrt(2)*((m2>=0?1.0:0.0)+ii*((-m2)>=0?1.0:0.0))+0.5*(0==m2?1:0);
                            std::complex<double> mu=theta1*theta2*(1+((-abs(m1*m2))>=0?1.0:0.0));
                            soc=mu*((abs(m1)==(abs(m2)+leftsign)?1.0:0.0)-pow(-1,(m1>=0?1:0)+(m2>=0?1:0))*(abs(m1)==(abs(m2)+rightsign)?1.0:0.0));

                            matrix(rows,cols)+=lambda*soc;
                        }

                    }
                }
            }
        }
    }

    void Hamiltonian::SetMatrixLocNew(const Eigen::Vector3d k,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::Matrix<lapack_complex_double ,-1,-1> &blockmatrix,int blpo,int AtomO)
    {

        std::complex<double> ii(0,1);
        int ZEROi = 0;
        int MONEi=-1; 
        int ONEi=1;    
        double ONEd=1.0;
        double ZEROd=0.0;
        int N=AtomO; 
        int M=AtomO; 
        int Nb=1;
        int Mb=1;

        int id,numprocs;
        blacs_pinfo_(&id, &numprocs);
        int ctxt, myrow, mycol;
        blacs_get_(&MONEi, &ZEROi, &ctxt);
        int info; 

        int procrows = sqrt(numprocs), proccols = sqrt(numprocs);
        blacs_gridinit_(&ctxt, "R", &procrows, &proccols);  
        blacs_gridinfo_(&ctxt, &procrows, &proccols, &myrow, &mycol );  
        int nrows = numroc_(&N, &Nb, &myrow, &ZEROi, &procrows);
        int lda = std::max(1,nrows);
        int ncols = numroc_(&M, &Mb, &mycol, &ZEROi, &proccols);
        blockmatrix.resize(nrows,ncols);	
        blockmatrix.setConstant({0,0});

        std::unordered_map<std::string, Material> &ms=material.materials;
        if(USE_long_parameters==0)
        {
        ms=material.materialsShort;
        }
        for(auto &cell : ListofMatrixElements)
        {
            Eigen::Vector3d dd=cell->d;
            Eigen::Vector3d ed=cell->ed;
            Eigen::Vector3d edn=ed.normalized();
            Eigen::Vector3d dn=dd.normalized();
            //std::cout<<"dd="<<dd<<std::endl;
            std::complex<double> efactor=exp(2.*M_PI*((dd).dot(k))*ii);
            //std::complex<double> efactor=exp(((dn*2.62).dot(k))*ii);




            int atom2=cell->Table.real();
            int atom1=cell->Table.imag();
            std::string leftType=(atom1%2==0?"c":"a");
            //std::string leftType="a";
            //UnitType[atom2];

            std::string rightType=(atom2%2==0?"c":"a");
            //std::string rightType="c";
            atom2*=orbNum;
            atom1*=orbNum;
            //UnitType[atom1];
            for(int io=0;io<orbNum;++io)
            {
                for(int jo=0;jo<orbNum;++jo)
                {
                    //std::string left=labels[jo];
                    //std::string right=labels[io];
                    Eigen::Vector4i leftorbital=qlabels[jo];
                    Eigen::Vector4i rightorbital=qlabels[io];
                    int rows=atom2+io-blpo;
                    int cols=atom1+jo-blpo;
                    if((rows>=0 && rows<AtomO) && (cols>=0 && cols<AtomO)) 
                    {

                        int sendr=rows%procrows;
                        int sendc=cols%proccols;
                        int recvr;
                        int recvc;
                        if (myrow == sendr && mycol == sendc){
                            //recvr=rows%nrows;
                            //recvc=cols%ncols;
                            recvr=rows/procrows;
                            recvc=cols/proccols;
                        }
                        std::string QWMaterial="GaAs";
                        std::string EBMaterial="AlAs";
                        if(rows==cols)
                        {
                            if((cell->cellme==QWMaterial && cell->SuperCell==EBMaterial) || (cell->cellme==EBMaterial && cell->SuperCell==QWMaterial))
                            {
                                if(myrow==sendr && mycol == sendc){

                                    blockmatrix(recvr,recvc).real=
                                        ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType)*0.5
                                        +ms.at(cell->SuperCell).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType)*0.5;

                                    //blockmatrix(recvr,recvc).real=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType);
                                }
                            }
                            else if(cell->cellme==QWMaterial && cell->SuperCell==QWMaterial)
                            {
                                if(myrow==sendr && mycol == sendc){
                                    blockmatrix(recvr,recvc).real=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType)+BandOffSet;
                                    //blockmatrix(recvr,recvc).real=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType)+BandOffSet;
                                }
                            }else{
                                if(myrow==sendr && mycol == sendc){
                                    blockmatrix(recvr,recvc).real=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType);

                                    //blockmatrix(recvr,recvc).real=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType);
                                }
                            }
                        }
                        else if(leftType!=rightType)
                        { 
                            if(myrow==sendr && mycol == sendc){
                                double res=0;

                                if(rows>cols)
                                {
                                    res=(SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->cellme)));
                                }else if(rows<cols)
                                {
                                    res=(SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->SuperCell)));
                                }
                                blockmatrix(recvr,recvc).real+=real(res*efactor);
                                blockmatrix(recvr,recvc).imag+=imag(res*efactor);

                            }
                        }
                    }
                }
            }
        }
        /*
           for(auto &cell : systemcell.ListofMatrixElements_MixType1)
           {

           int atom2=cell->Table.real();
           int atom1=cell->Table.imag();
           std::string leftType=(atom1%2==0?"c":"a");
           std::string rightType=(atom2%2==0?"c":"a");

           atom2*=orbNum;

           for(int io=0;io<orbNum;++io)
           {
           Eigen::Vector4i rightorbital=qlabels[io];
           Eigen::Vector4i leftorbital=qlabels[io];
           int rows=atom2+io-blpo;

           if((rows>=0 && rows<AtomO)) 
           {
           int sendr=rows%procrows;
           int sendc=rows%proccols;
           int recvr;
           int recvc;
           if (myrow == sendr && mycol == sendc){
        //recvr=rows%nrows;
        //recvc=cols%ncols;
        recvr=rows/procrows;
        recvc=rows/proccols;
        }
        if(myrow==sendr && mycol == sendc){

        blockmatrix(recvr,recvc).real=
        (ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType))*0.5
        +(ms.at(cell->SuperCell).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType))*0.5
        +BandOffSet*0.5;
        }

        }

        }
        }
        */
        blacs_gridexit_(&ctxt);
    }


    void Hamiltonian::SetMatrixLocWithoutShift(const Eigen::Vector3d k,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::Matrix<lapack_complex_double,-1,-1> &matrix)
    {

        std::complex<double> ii(0,1);
        int ZEROi = 0;
        int MONEi=-1; 
        int ONEi=1;    
        double ONEd=1.0;
        double ZEROd=0.0;
        int N=MatrixSize; 
        int M=MatrixSize; 
        int Nb=1;
        int Mb=1;

        int id,numprocs;
        blacs_pinfo_(&id, &numprocs);
        int ctxt, myrow, mycol;
        blacs_get_(&MONEi, &ZEROi, &ctxt);
        int info; 

        int procrows = sqrt(numprocs), proccols = sqrt(numprocs);
        blacs_gridinit_(&ctxt, "R", &procrows, &proccols);  
        blacs_gridinfo_(&ctxt, &procrows, &proccols, &myrow, &mycol );  
        int nrows = numroc_(&N, &Nb, &myrow, &ZEROi, &procrows);
        int lda = std::max(1,nrows);
        int ncols = numroc_(&M, &Mb, &mycol, &ZEROi, &proccols);
        matrix.resize(nrows,ncols);	
        matrix.setConstant({0,0});

        std::unordered_map<std::string, Material> &ms=material.materials;
        if(USE_long_parameters==0)
        {
            ms=material.materialsShort;
        }
        /*
           if(id==0)
           {
           for(auto &cell : ListofMatrixElements)
           {
           int atom2=cell->Table.real();
           int atom1=cell->Table.imag();
           atom2*=orbNum;
           atom1*=orbNum;
           int rows=atom2;
           int cols=atom1;
           int shifts=systemcell.supercell[cell->Typeme]->blpo;
           printf("typeme=%d[%d]\n",shifts,cols-shifts);
           }
           }
           */
        for(auto &cell : ListofMatrixElements)
        {
            Eigen::Vector3d dd=cell->d;
            Eigen::Vector3d ed=cell->ed;
            Eigen::Vector3d edn=ed.normalized();
            Eigen::Vector3d dn=dd.normalized();
            //std::cout<<"dd="<<dd<<std::endl;
            std::complex<double> efactor=exp(2.*M_PI*((dd).dot(k))*ii);
            //std::complex<double> efactor=exp(((dn*2.62).dot(k))*ii);




            int atom2=cell->Table.real();
            int atom1=cell->Table.imag();
            std::string leftType=(atom1%2==0?"c":"a");
            //std::string leftType="a";
            //UnitType[atom2];

            std::string rightType=(atom2%2==0?"c":"a");
            //std::string rightType="c";
            atom2*=orbNum;
            atom1*=orbNum;
            //UnitType[atom1];
            for(int io=0;io<orbNum;++io)
            {
                for(int jo=0;jo<orbNum;++jo)
                {
                    //std::string left=labels[jo];
                    //std::string right=labels[io];
                    Eigen::Vector4i leftorbital=qlabels[jo];
                    Eigen::Vector4i rightorbital=qlabels[io];
                    int rows=atom2+io;
                    int cols=atom1+jo;

                    //if(id==0) printf("newcols[%d]=%d{%d}\n",cell->Typeme,newcols,systemcell.supercell[cell->Typeme]->blpo);
                    int sendr=rows%procrows;
                    int sendc=cols%proccols;
                    int recvr;
                    int recvc;
                    if (myrow == sendr && mycol == sendc){
                        //recvr=rows%nrows;
                        //recvc=cols%ncols;
                        recvr=rows/procrows;
                        recvc=cols/proccols;
                    }
                    std::string QWMaterial="GaAs";
                    std::string EBMaterial="AlAs";
                    /*
                       if(rows==cols)
                       {
                       if((cell->cellme==QWMaterial && cell->SuperCell==EBMaterial) || (cell->cellme==EBMaterial && cell->SuperCell==QWMaterial))
                       {
                       if(myrow==sendr && mycol == sendc){
                       matrix(recvr,recvc).real=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType);
                       }
                       }
                       else if(cell->cellme==QWMaterial && cell->SuperCell==QWMaterial)
                       {
                       if(myrow==sendr && mycol == sendc){
                       matrix(recvr,recvc).real=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType)+BandOffSet;
                       }
                       }else{
                       if(myrow==sendr && mycol == sendc){
                       matrix(recvr,recvc).real=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType);
                       }
                       }
                       }
                       else if(leftType!=rightType)
                       { 
                       if(myrow==sendr && mycol == sendc){
                       double res=(SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->cellme))*0.5+SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->SuperCell))*0.5);
                       matrix(recvr,recvc).real+=real(res*efactor);
                       matrix(recvr,recvc).imag+=imag(res*efactor);
                       }
                       }
                       */
                    if(rows==cols)
                    {
                        if((cell->cellme==QWMaterial && cell->SuperCell==EBMaterial) || (cell->cellme==EBMaterial && cell->SuperCell==QWMaterial))
                        {
                            if(myrow==sendr && mycol == sendc){

                                matrix(recvr,recvc).real=
                                    ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType)*0.5
                                    +ms.at(cell->SuperCell).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType)*0.5;

                                //matrix(recvr,recvc).real=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType);
                            }
                        }
                        else if(cell->cellme==QWMaterial && cell->SuperCell==QWMaterial)
                        {
                            if(myrow==sendr && mycol == sendc){
                                matrix(recvr,recvc).real=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType)+BandOffSet;
                                //matrix(recvr,recvc).real=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType)+BandOffSet;
                            }
                        }else{
                            if(myrow==sendr && mycol == sendc){
                                matrix(recvr,recvc).real=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType);

                                //matrix(recvr,recvc).real=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType);
                            }
                        }
                    }
                    else if(leftType!=rightType)
                    { 
                        if(myrow==sendr && mycol == sendc){
                            double res=0;
                            /*
                               if((cell->cellme==QWMaterial && cell->SuperCell==EBMaterial)){
                               res=(SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->cellme)));
                               }else if((cell->cellme==EBMaterial && cell->SuperCell==QWMaterial)){
                               res=(SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->SuperCell)));
                               }else{
                               res=(SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->cellme)));
                               }
                               */
                            if(rows>cols)
                            {
                                res=(SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->cellme)));
                            }else if(rows<cols)
                            {
                                res=(SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->SuperCell)));
                            }
                            matrix(recvr,recvc).real+=real(res*efactor);
                            matrix(recvr,recvc).imag+=imag(res*efactor);


                            //double res=(SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->cellme))*0.5+SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->SuperCell))*0.5);
                            //matrix(recvr,recvc).real+=real(res*efactor);
                            //matrix(recvr,recvc).imag+=imag(res*efactor);
                        }
                    }

                }
            }
        }
        /*
           for(auto &cell : systemcell.ListofMatrixElements_MixType1)
           {

           int atom2=cell->Table.real();
           int atom1=cell->Table.imag();
           std::string leftType=(atom1%2==0?"c":"a");
           std::string rightType=(atom2%2==0?"c":"a");

           atom2*=orbNum;

           for(int io=0;io<orbNum;++io)
           {
           Eigen::Vector4i rightorbital=qlabels[io];
           Eigen::Vector4i leftorbital=qlabels[io];
           int rows=atom2+io;

           int sendr=rows%procrows;
           int sendc=rows%proccols;
           int recvr;
           int recvc;
           if (myrow == sendr && mycol == sendc){
        //recvr=rows%nrows;
        //recvc=cols%ncols;
        recvr=rows/procrows;
        recvc=rows/proccols;
        }
        if(myrow==sendr && mycol == sendc){

        matrix(recvr,recvc).real=
        (ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType))*0.5
        +(ms.at(cell->SuperCell).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType))*0.5
        +BandOffSet*0.5;
        }


        }
        }
        */
        blacs_gridexit_(&ctxt);
    }
    void Hamiltonian::SetMatrixLoc(const Eigen::Vector3d k,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::Matrix<lapack_complex_double ,-1,-1> &matrix)
    {

        std::complex<double> ii(0,1);
        int ZEROi = 0;
        int MONEi=-1; 
        int ONEi=1;    
        double ONEd=1.0;
        double ZEROd=0.0;
        int N=MatrixSize; 
        int M=MatrixSize; 
        int Nb=1;
        int Mb=1;

        int id,numprocs;
        blacs_pinfo_(&id, &numprocs);
        int ctxt, myrow, mycol;
        blacs_get_(&MONEi, &ZEROi, &ctxt);
        int info; 

        int procrows = sqrt(numprocs), proccols = sqrt(numprocs);
        blacs_gridinit_(&ctxt, "R", &procrows, &proccols);  
        blacs_gridinfo_(&ctxt, &procrows, &proccols, &myrow, &mycol );  
        int nrows = numroc_(&N, &Nb, &myrow, &ZEROi, &procrows);
        int lda = std::max(1,nrows);
        int ncols = numroc_(&M, &Mb, &mycol, &ZEROi, &proccols);
        matrix.resize(nrows,ncols);	
        matrix.setConstant({0,0});
        /*
        //std::unordered_map<std::string, Material> &ms=material.materialsShort;
        std::unordered_map<std::string, Material> &ms=material.materials;

        for(auto &cell : ListofMatrixElements)
        {
        Eigen::Vector3d dd=cell->d;
        Eigen::Vector3d ed=cell->ed;
        Eigen::Vector3d edn=ed.normalized();
        Eigen::Vector3d dn=dd.normalized();
        //std::cout<<"dd="<<dd<<std::endl;
        std::complex<double> efactor=exp(2.*M_PI*((dd).dot(k))*ii);
        //std::complex<double> efactor=exp(((dn*2.62).dot(k))*ii);




        int atom2=cell->Table.real();
        int atom1=cell->Table.imag();
        std::string leftType=(atom1%2==0?"c":"a");
        //std::string leftType="a";
        //UnitType[atom2];

        std::string rightType=(atom2%2==0?"c":"a");
        //std::string rightType="c";
        atom2*=orbNum;
        atom1*=orbNum;
        //UnitType[atom1];
        for(int io=0;io<orbNum;++io)
        {
        for(int jo=0;jo<orbNum;++jo)
        {
        //std::string left=labels[jo];
        //std::string right=labels[io];
        Eigen::Vector4i leftorbital=qlabels[jo];
        Eigen::Vector4i rightorbital=qlabels[io];
        int rows=atom2+io;
        int cols=atom1+jo;

        int shiftcols=cols-systemcell.supercell[cell->Typeme]->blpo;
        if(shiftcols>=0){
        //if(id==0) printf("newcols[%d]=%d{%d}\n",cell->Typeme,newcols,systemcell.supercell[cell->Typeme]->blpo);
        int sendr=rows%procrows;
        int sendc=cols%proccols;
        int recvr;
        int recvc;
        if (myrow == sendr && mycol == sendc){
        //recvr=rows%nrows;
        //recvc=cols%ncols;
        recvr=rows/procrows;
        recvc=shiftcols/proccols;
        }
        std::string QWMaterial="GaAs";
        std::string EBMaterial="AlAs";

        if(rows==cols)
        {
        if((cell->cellme==QWMaterial && cell->SuperCell==EBMaterial) || (cell->cellme==EBMaterial && cell->SuperCell==QWMaterial))
        {
        if(myrow==sendr && mycol == sendc){

        matrix(recvr,recvc).real=
        ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType)*0.5
        +ms.at(cell->SuperCell).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType)*0.5;

        //matrix(recvr,recvc).real=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType);
        }
        }
        else if(cell->cellme==QWMaterial && cell->SuperCell==QWMaterial)
        {
        if(myrow==sendr && mycol == sendc){
        matrix(recvr,recvc).real=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType)+BandOffSet;
        //matrix(recvr,recvc).real=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType)+BandOffSet;
    }
    }else{
        if(myrow==sendr && mycol == sendc){
            matrix(recvr,recvc).real=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType);

            //matrix(recvr,recvc).real=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType);
        }
    }
    }
        else if(leftType!=rightType)
        { 
            if(myrow==sendr && mycol == sendc){
                double res=0;

                if(rows>cols)
                {
                    res=(SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->cellme)));
                }else if(rows<cols)
                {
                    res=(SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->SuperCell)));
                }
                matrix(recvr,recvc).real+=real(res*efactor);
                matrix(recvr,recvc).imag+=imag(res*efactor);


                //double res=(SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->cellme))*0.5+SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->SuperCell))*0.5);
                //matrix(recvr,recvc).real+=real(res*efactor);
                //matrix(recvr,recvc).imag+=imag(res*efactor);
            }
        }

    }
    }
    }
    }
    */
        blacs_gridexit_(&ctxt);
    }




    double Hamiltonian::SKFormular(std::string leftType,std::string rightType,int n1,int l1,int m1,int n2,int l2,int m2,Eigen::Vector3d dn,const Material *matter)
    {
        int Lless=std::min(l1,l2);
        int Lmax=std::max(l1,l2);
        std::string MatSTR;
        if(Lless==l1)
        {
            MatSTR=std::to_string(n1)+std::to_string(l1)+leftType+std::to_string(n2)+std::to_string(l2)+rightType;
        }else{
            MatSTR=std::to_string(n2)+std::to_string(l2)+rightType+std::to_string(n1)+std::to_string(l1)+leftType;
        }
        //fprintf(stderr,"%s--over\n",MatSTR.c_str());
        if(dn(2)==1)
        {
            return (m1==m2)?pow(-1,(l1-l2+abs(l1-l2))*0.5)*matter->param.at(MatSTR+str_delta_m[m1]):0;
        }else{
            double sum=0;
            double parame;
            //if(matter->param.count(MatSTR)!=1){
            //parame=0;
            //}else{
            parame=matter->param.at(MatSTR);
            //}
            sum+=2*func_A(m1,dn)*func_A(m2,dn)*func_D(l1,abs(m1),0,dn)*func_D(l2,abs(m2),0,dn)*parame;
            //sum+=2*func_A(m1,dn)*func_A(m2,dn)*getfunc_D(l1,abs(m1),0,dn(2))*getfunc_D(l2,abs(m2),0,dn(2))*parame;
            //sum+=2*func_A(m1,dn)*func_A(m2,dn)*mappD[dn(2)][l1][abs(m1)*(l1+1)].real()*mappD[dn(2)][l2][abs(m2)*(l2+1)].real()*parame;
            for(int m_prime=1;m_prime<=Lless;++m_prime)
            {
                std::string MatSTRplus;
                MatSTRplus=MatSTR+str_delta_m[abs(m_prime)];
                //fprintf(stderr,"%s\n",MatSTRplus.c_str());
                double parame;
                //if(matter->param.count(MatSTRplus)!=1){
                //parame=0;
                //}else{
                parame=matter->param.at(MatSTRplus);
                //}
                sum+=(func_S(l1,m1,abs(m_prime),dn)*func_S(l2,m2,abs(m_prime),dn)+func_T(l1,m1,abs(m_prime),dn)*func_T(l2,m2,abs(m_prime),dn))*parame;
            }

            return std::pow(-1,(l1-l2+abs(l1-l2))/2)*sum;
        }

    }

    double Hamiltonian::func_A(int m,Eigen::Vector3d dn)
    {
        if(m==0) return 1/sqrt(2);
        std::complex<double> gamma={-dn(0)/sqrt(1-pow(dn(2),2)),dn(1)/sqrt(1-pow(dn(2),2))};
        double cos_m=std::pow(gamma,abs(m)).real();
        double sin_m=std::pow(gamma,abs(m)).imag();
        return pow(-1,abs(m))*(m>0?cos_m:-sin_m);
        //return pow(-1,abs(m))*(m>0?cos(abs(m)*asin(dn(1)/sqrt(1-pow(dn(2),2)))):-sin(abs(m)*asin(dn(1)/sqrt(1-pow(dn(2),2)))));
    }

    double Hamiltonian::func_B(int m,Eigen::Vector3d dn)
    {
        std::complex<double> gamma={-dn(0)/sqrt(1-pow(dn(2),2)),dn(1)/sqrt(1-pow(dn(2),2))};
        double cos_m=std::pow(gamma,abs(m)).real();
        double sin_m=std::pow(gamma,abs(m)).imag();
        return pow(-1,abs(m))*(m>0?sin_m:cos_m);
        //return pow(-1,abs(m))*(m>0?sin(abs(m)*asin(dn(1)/sqrt(1-pow(dn(2),2)))):cos(abs(m)*asin(dn(1)/sqrt(1-pow(dn(2),2)))));
    }


    double Hamiltonian::func_T(int l,int m,int m_prime,Eigen::Vector3d dn)
    {
        if(m==0) return 0;
        else return func_B(m,dn)*(pow(-1,abs(m_prime))*func_D(l,abs(m),abs(m_prime),dn)-func_D(l,abs(m),-abs(m_prime),dn));
        //else return func_B(m,dn)*(pow(-1,abs(m_prime))*getfunc_D(l,abs(m),abs(m_prime),dn(2))-getfunc_D(l,abs(m),-abs(m_prime),dn(2)));
        //else return func_B(m,dn)*(pow(-1,abs(m_prime))*mappD[dn(2)][l][abs(m)*(l+1)+abs(m_prime)].real()-mappD[dn(2)][l][abs(m)*(l+1)+abs(m_prime)].imag());
    }

    double Hamiltonian::func_S(int l,int m,int m_prime,Eigen::Vector3d dn)
    {
        return func_A(m,dn)*(pow(-1,abs(m_prime))*func_D(l,abs(m),abs(m_prime),dn)+func_D(l,abs(m),-abs(m_prime),dn));
        //return func_A(m,dn)*(pow(-1,abs(m_prime))*getfunc_D(l,abs(m),abs(m_prime),dn(2))+getfunc_D(l,abs(m),-abs(m_prime),dn(2)));
        //return func_A(m,dn)*(pow(-1,abs(m_prime))*mappD[dn(2)][l][abs(m)*(l+1)+abs(m_prime)].real()+mappD[dn(2)][l][abs(m)*(l+1)+abs(m_prime)].imag());
    }

    double Hamiltonian::factorials(int n)
    {
        double ret=1;
        for(int i=1;i<=n;++i)
        {
            ret*=i;
        }
        return ret;
    }

    double Hamiltonian::getfunc_D(int l,int m,int m_prime,double N)
    {
        int pos=m*(l+1)+abs(m_prime);
        if(m_prime<0)
        {
            return mappD[N][l][pos].imag();
        }else{
            return mappD[N][l][pos].real();
        }
    }

    double Hamiltonian::func_D(int l,int m,int m_prime,Eigen::Vector3d dn)
    {
        double N=dn(2);
        double ret;
        double tmp=0;

        if(l+m_prime >=0 &&l-m_prime>=0 &&l+m>=0 &&l-m>=0)
        {
            ret=std::pow((1+N)/2,l)*std::pow((1-N)*1./(1+N),(m*1./2-m_prime*1./2))*sqrt(factorials(l+m_prime)*factorials(l-m_prime)*factorials(l+m)*factorials(l-m));
        }else{
            ret=0; 
        }
        for(int t=0;t<=2*l+1;++t)
        {
            double dnom=0;
            if((l+m_prime-t)>=0 && (l-m-t)>=0 && t>=0 && (t+m-m_prime)>=0)
            {
                dnom=factorials(l+m_prime-t)*factorials(l-m-t)*factorials(t)*factorials(t+m-m_prime);
            }
            if(dnom!=0) tmp+=pow(-1.0,t)*std::pow((1-N)*1./(1+N),t)/dnom;

        }
        return ret*tmp; 
    }






    }
