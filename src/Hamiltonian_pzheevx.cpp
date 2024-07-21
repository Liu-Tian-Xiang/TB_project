
#include "Hamiltonian.hpp"

namespace TightBinding
{
    void Hamiltonian::PFTBH(int cellme,int blpo,int AtomO,std::vector<int> vblpo,std::vector<int> vAtomO,const Eigen::Vector3d& k,std::vector<Eigen::Vector3d> chain,Eigen::Matrix<lapack_complex_double,-1,-1> &sMatrix,int &Startp,int &SecDim,Eigen::Matrix3d &Ugm,Eigen::Matrix3d &ugm,Eigen::Matrix3d &uam,std::vector<Eigen::MatrixXcd> &egvct,Eigen::VectorXd &svalue,int &myrow,int &mycol)//without dtag
{

    int printM=cellme;
    Eigen::Vector3d shift=systemcell.supercell[cellme]->Unit[0];
    assert(vblpo.size()==vAtomO.size());
    int blnum=vblpo.size();

    std::complex<double> ii(0,1);
    int fitposCol=0;

    std::vector<int> LAtomO;
    std::vector<std::pair<double,double>> tmpEnergyCutoff;
    LAtomO.resize(0);
    tmpEnergyCutoff.resize(0);

    if(systemcell.supercell[0]->vblpo2.size())
    {
	LAtomO.push_back(systemcell.supercell[0]->vblpo2[0]);
    }else{
	LAtomO.push_back(systemcell.supercell[0]->vblpo.size()/systemcell.uvector.size());
    }

    tmpEnergyCutoff.push_back(systemcell.supercell[0]->EnergyCutoff);
    for(int iii=0;iii<systemcell.supercell[0]->vAtomO2.size();++iii)
    {
	LAtomO.push_back(systemcell.supercell[0]->vAtomO2[iii]);
	tmpEnergyCutoff.push_back(systemcell.supercell[0]->vEnergyCutoff[iii]);
    }


    Startp=0;
    int gi=0;
    for(int i1=0;i1<LAtomO.size();++i1)
    {
	for(int t=0;t<LAtomO[i1];++t)
	{
	    for(int i2=0;i2<systemcell.uvector.size();++i2)
	    {
		int blen=vAtomO[gi];
		Eigen::VectorXd tmpGsolver=svalue.segment(vblpo[gi]-vblpo[0],blen);
		for(int in=0;in<tmpGsolver.size();++in)
		{
		    if(tmpGsolver(in)<tmpEnergyCutoff[i1].first) Startp+=1;
		}
		gi+=1;
	    }
	}
    }




    SecDim=0;
    int i=0;
    gi=0;
    for(int i1=0;i1<LAtomO.size();++i1)
    {
	for(int t=0;t<LAtomO[i1];++t)
	{
	    for(int i2=0;i2<systemcell.uvector.size();++i2)
	    {
		int blen=vAtomO[gi];
		//gsolver.compute(matrixScatter.block(vblpo[gi],vblpo[gi],blen,blen),Eigen::MatrixXcd::Identity(blen,blen) ,Eigen::DecompositionOptions::ComputeEigenvectors|Eigen::Ax_lBx);
		//assert(gsolver.info() == Eigen::ComputationInfo::Success);
		Eigen::VectorXd tmpGsolver=svalue.segment(vblpo[gi]-vblpo[0],blen);
		//for(int in=0;in<gsolver.eigenvalues().size();++in)
		for(int in=0;in<tmpGsolver.size();++in)
		{

		    //if(gsolver.eigenvalues()(in)<tmpEnergyCutoff[i1].first) Startp+=1;
		    //if(gsolver.eigenvalues()(in)>tmpEnergyCutoff[i1].first && gsolver.eigenvalues()(in)<tmpEnergyCutoff[i1].second)
		    //if(tmpGsolver(in)<tmpEnergyCutoff[i1].first) Startp+=1;
		    if(tmpGsolver(in)>tmpEnergyCutoff[i1].first && tmpGsolver(in)<tmpEnergyCutoff[i1].second)
		    {
			int colnum=i;
			int Fi=colnum/(systemcell.UnitCellNumber*orbNum);
			int Fcc=colnum%(systemcell.UnitCellNumber*orbNum)/orbNum;
			int Fkj=colnum%(systemcell.UnitCellNumber*orbNum)%orbNum;
			//PFCol3(egvct[gi].col(Fcc*orbNum+Fkj),cellme,blpo,AtomO,vblpo,vAtomO, k, chain,Startp+SecDim,Fi,Fcc,Fkj,Ugm,ugm,uam,sMatrix,myrow,mycol);
			PFCol3(egvct[gi].col(Fcc*orbNum+Fkj),cellme,blpo,AtomO,vblpo,vAtomO, k, chain,SecDim,Fi,Fcc,Fkj,Ugm,ugm,uam,sMatrix,myrow,mycol);
			SecDim+=1;
		    }
		    i+=1;
		}
		gi+=1;
	    } 
	}
    }
    //fitting.conservativeResize(AtomO,SecDim);
}

    void Hamiltonian::tmpPFTBH(int cellme,int blpo,int AtomO,std::vector<int> vblpo,std::vector<int> vAtomO,const Eigen::Vector3d& k,std::vector<Eigen::Vector3d> chain,Eigen::Matrix<lapack_complex_double,-1,-1> &sMatrix,int &Startp,int &SecDim,Eigen::Matrix3d &Ugm,Eigen::Matrix3d &ugm,Eigen::Matrix3d &uam,int &myrow,int &mycol)//without dtag
{
    int printM=cellme;
    Eigen::Vector3d shift=systemcell.supercell[cellme]->Unit[0];
    assert(vblpo.size()==vAtomO.size());
    std::vector<Eigen::MatrixXcd> egvct;
    egvct.resize(0);
    int blnum=vblpo.size();

    std::complex<double> ii(0,1);
    int fitposCol=0;

    std::vector<int> LAtomO;
    std::vector<std::pair<double,double>> tmpEnergyCutoff;
    LAtomO.resize(0);
    tmpEnergyCutoff.resize(0);
    if(systemcell.supercell[0]->vblpo2.size())
    {
	LAtomO.push_back(systemcell.supercell[0]->vblpo2[0]);
    }else{
	LAtomO.push_back(systemcell.supercell[0]->vblpo.size()/systemcell.uvector.size());
    }
    tmpEnergyCutoff.push_back(systemcell.supercell[0]->EnergyCutoff);
    for(int iii=0;iii<systemcell.supercell[0]->vAtomO2.size();++iii)
    {
	LAtomO.push_back(systemcell.supercell[0]->vAtomO2[iii]);
	tmpEnergyCutoff.push_back(systemcell.supercell[0]->vEnergyCutoff[iii]);
    }

    Startp=0;
    for(int i1=0;i1<LAtomO.size();++i1)
    {
	int cellpos=systemcell.uvector.size()*LAtomO[i1]*vAtomO[0]*BandGapIndex/orbNum;
	Startp+=cellpos-tmpEnergyCutoff[i1].first;
    }

    SecDim=0;
    int i=0;
    int gi=0;
    int lai=0;
    for(int i1=0;i1<LAtomO.size();++i1)
    {
	int cellpos=systemcell.uvector.size()*LAtomO[i1]*vAtomO[0]*BandGapIndex/orbNum;
	//Startp+=cellpos-tmpEnergyCutoff[i1].first;

	Eigen::VectorXd tmpegv;
	tmpegv.resize(systemcell.uvector.size()*LAtomO[i1]*vAtomO[0]);
	tmpegv.setZero();
	std::vector<Eigen::MatrixXcd> egvct;
	egvct.resize(0);

	int flag=0;
	for(int t=0;t<LAtomO[i1];++t)
	{
	    for(int i2=0;i2<systemcell.uvector.size();++i2)
	    {
		int blen=vAtomO[gi];
		Eigen::MatrixXcd tmpMatrix;
		SetMatrixPrime(k,systemcell.ListofMatrixElements_Type3,tmpMatrix,vblpo[gi],blen);
/*
		gsolver.compute(tmpMatrix,Eigen::MatrixXcd::Identity(blen,blen) ,Eigen::DecompositionOptions::ComputeEigenvectors|Eigen::Ax_lBx);
		//gsolver.compute(matrixScatter.block(vblpo[gi],vblpo[gi],blen,blen),Eigen::MatrixXcd::Identity(blen,blen) ,Eigen::DecompositionOptions::ComputeEigenvectors|Eigen::Ax_lBx);
		assert(gsolver.info() == Eigen::ComputationInfo::Success);
		//Eigen::VectorXd gtmpegv=gsolver.eigenvalues();
		tmpegv.segment(flag*blen,blen)=gsolver.eigenvalues();
		egvct.push_back(gsolver.eigenvectors());
*/
		solver.compute(tmpMatrix,Eigen::DecompositionOptions::ComputeEigenvectors);
		assert(solver.info() == Eigen::ComputationInfo::Success);
		tmpegv.segment(flag*blen,blen)=solver.eigenvalues();
		egvct.push_back(solver.eigenvectors());

		gi+=1;
		flag+=1;
	    }
	}

	std::vector<std::complex<double>> ord;
	for(int oi=0;oi<tmpegv.size();++oi){
	    ord.push_back({oi*1.,tmpegv[oi]});
	}
	for(int oi=0;oi<tmpegv.size();++oi){
	    for(int oj=0;oj<tmpegv.size()-oi-1;++oj){
		if(ord[oj].imag()>ord[oj+1].imag())
		    //if(tmpegv(oj)>tmpegv(oj+1))
		{
		    //double tmp=tmpegv(oj);
		    //tmpegv(oj)=tmpegv(oj+1);
		    //tmpegv(oj+1)=tmp;
		    std::complex<double> tmp2=ord[oj];
		    ord[oj]=ord[oj+1];
		    ord[oj+1]=tmp2;
		    //fitting.col(j).swap(fitting.col(j+1));
		}
	    }
	}

	flag=0;
	for(int t=0;t<LAtomO[i1];++t)
	{
	    for(int i2=0;i2<systemcell.uvector.size();++i2)
	    {

		for(int in=0;in<egvct[flag].rows();++in)
		{
		    int blen=egvct[flag].rows();
		    if((flag*blen+in)>=(cellpos-tmpEnergyCutoff[i1].first) && (flag*blen+in)<(cellpos+tmpEnergyCutoff[i1].second))
		    {
			int fitpos=lai+ord[flag*blen+in].real();
			int colnum=fitpos;
			int Fi=colnum/(systemcell.UnitCellNumber*orbNum);
			int Fcc=colnum%(systemcell.UnitCellNumber*orbNum)/orbNum;
			int Fkj=colnum%(systemcell.UnitCellNumber*orbNum)%orbNum;
			//PFCol3(egvct[((int)(ord[flag*blen+in].real()))/blen].col(Fcc*orbNum+Fkj),cellme,blpo,AtomO,vblpo,vAtomO, k, chain,Startp+SecDim,Fi,Fcc,Fkj,Ugm,ugm,uam,sMatrix,myrow,mycol);
			PFCol3(egvct[((int)(ord[flag*blen+in].real()))/blen].col(Fcc*orbNum+Fkj),cellme,blpo,AtomO,vblpo,vAtomO, k, chain,SecDim,Fi,Fcc,Fkj,Ugm,ugm,uam,sMatrix,myrow,mycol);
			SecDim+=1;
		    }
		    i+=1;
		}
		flag+=1;
	    }
	}
	lai+=flag*vAtomO[0];
    }
    //fitting.conservativeResize(AtomO,SecDim);
}

    void Hamiltonian::PTransform(Eigen::Matrix<lapack_complex_double,-1,-1> &ConsideredMatrix,std::vector<Eigen::Matrix<lapack_complex_double,-1,-1>> &BlockData)
    {

                int ZEROi = 0;
                int MONEi=-1; 
                int ONEi=1;    
                double ONEd=1.0;
                double ZEROd=0.0;

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

        for(auto &cell : systemcell.supercell)
        {
            if(cell->type)
            {

                int pp=cell->Tag2;
                int mAtom=cell->SecDim;
                int blpo=0;//mposi(cell->Tag2);

                const Eigen::Matrix<lapack_complex_double,-1,-1> &lM=BlockData[cell->Tag3];//.block(0,cell->Startp,cell->AtomO,mAtom);
                const Eigen::Matrix<lapack_complex_double,-1,-1> &rM=BlockData[cell->Tag3];//.block(0,cell->Startp,cell->AtomO,mAtom);

                int N=cell->AtomO; 
                int M=cell->AtomO; 

                int nrows = numroc_(&N, &Nb, &myrow, &ZEROi, &procrows);
                int lda = std::max(1,nrows);
                int ncols = numroc_(&M, &Mb, &mycol, &ZEROi, &proccols);
                int descA[9];
                descinit_(descA, &N, &M, &Nb, &Mb,&ZEROi,&ZEROi,&ctxt, &lda, &info);

                Eigen::Matrix<lapack_complex_double,-1,-1> D_loc;
                D_loc.resize(nrows,ncols);
                D_loc.setConstant({0,0});
                lapack_complex_double ONEz={1,0};

                
                   int MN=cell->SecDim;
                   int posX=cell->Startp+1;
 
                   int sectionSZ=MN;
		   while(N%sectionSZ)
		   {
		       sectionSZ+=1;
		   }
		   Eigen::Matrix<lapack_complex_double,-1,-1> C_loc;
		   int tmp_nrows = numroc_(&sectionSZ, &Nb, &myrow, &ZEROi, &procrows);
		   int tmp_lda = std::max(1,tmp_nrows);
		   int tmp_ncols = numroc_(&sectionSZ, &Mb, &mycol, &ZEROi, &proccols);
		   int descC[9];
		   descinit_(descC, &sectionSZ, &sectionSZ, &Nb, &Mb,&ZEROi,&ZEROi,&ctxt, &tmp_lda, &info);
		   C_loc.resize(tmp_nrows,tmp_ncols);
                   int sectionNUM=N/sectionSZ;
                   for(int o=1;o<=sectionNUM;++o)
                   {int pmno=(o-1)*sectionSZ+1;
                       for(int oo=1;oo<=sectionNUM;++oo)
                       {int pmnoo=(oo-1)*sectionSZ+1;

                           if((pmno)<=MN && (pmnoo)<=MN && o<=oo)
                           {

                               for(int oo2=1;oo2<=sectionNUM;++oo2)
                               {int pmnoo2=(oo2-1)*sectionSZ+1;
                                   for(int oo1=1;oo1<=sectionNUM;++oo1)
                                   {int pmnoo1=(oo1-1)*sectionSZ+1;

                                       if((fabs(oo1-oo2)<=1 || (fabs(oo2-oo1)==sectionNUM-1 && (oo1==1 || oo1==sectionNUM))))
                                       {

					   C_loc.setConstant({0,0});
					   pzgemm_("c","n",&sectionSZ,&sectionSZ,&sectionSZ,&ONEz,lM.data(),&pmnoo1,&pmno,descA,ConsideredMatrix.data(),&pmnoo1,&pmnoo2,descA,&ONEz,C_loc.data(),&ONEi,&ONEi,descC);
                                           pzgemm_("n","n",&sectionSZ,&sectionSZ,&sectionSZ,&ONEz,C_loc.data(),&ONEi,&ONEi,descC,rM.data(),&pmnoo2,&pmnoo,descA,&ONEz,D_loc.data(),&pmno,&pmnoo,descA);
                                       }
                                   }
                               }

                           }

                       }
                   }
		   C_loc.resize(0,0);
                   ConsideredMatrix=D_loc;

                  


if(0)
{
                Eigen::Matrix<lapack_complex_double,-1,-1> C_loc;
                C_loc.resize(nrows,ncols);
                C_loc.setConstant({0,0});
                   int sectionSZ=512;//MN;
                   int sectionNUM=N/sectionSZ;
                   for(int o=1;o<=sectionNUM;++o)
		   {
			   int mno=sectionSZ;
			   int pmno=(o-1)*sectionSZ+1;
			   if((pmno)<=MN)
			   {
				   for(int oo=1;oo<=sectionNUM;++oo)
				   {
					   //if(o<=(oo+1))
					   {
						   int mnoo=sectionSZ;
						   int pmnoo=(oo-1)*sectionSZ+1;
						   for(int oo1=1;oo1<=sectionNUM;++oo1)
						   {
							   int mnoo1=sectionSZ;
							   int pmnoo1=(oo1-1)*sectionSZ+1;
							   if((fabs(oo-oo1)<=1 || (fabs(oo-oo1)==sectionNUM-1 && (oo==1 || oo==sectionNUM))))
							   {
								   pzgemm_("c","n",&mno,&mno,&mno,&ONEz,lM.data(),&pmnoo1,&pmno,descA,ConsideredMatrix.data(),&pmnoo1,&pmnoo,descA,&ONEz,C_loc.data(),&pmno,&pmnoo,descA);
							   }
						   }
					   } 
				   }
			   } 
		   }
                   //pzgemm_("c","n",&N,&N,&N,&ONEz,lM.data(),&ONEi,&ONEi,descA,ConsideredMatrix.data(),&ONEi,&ONEi,descA,&ONEz,C_loc.data(),&ONEi,&ONEi,descA);
                   ConsideredMatrix.setConstant({0,0});
/*
                   for(int o=1;o<=sectionNUM;++o)
                   {
                       int mno=sectionSZ;
                       int pmno=(o-1)*sectionSZ+1;
                       if((pmno)<=MN)
                       {
                           for(int oo=1;oo<=sectionNUM;++oo)
                           {
                               int mnoo=sectionSZ;
                               int pmnoo=(oo-1)*sectionSZ+1;

                                   if(o<=oo)
                                   {
                               for(int oo1=1;oo1<=sectionNUM;++oo1)
                               {
                                   int mnoo1=sectionSZ;
                                   int pmnoo1=(oo1-1)*sectionSZ+1;
                                   //pzgemm_("c","n",&mno,&mno,&mno,&ONEz,lM.data(),&pmnoo1,&pmno,descA,ConsideredMatrix.data(),&pmnoo1,&pmnoo,descA,&ONEz,C_loc.data(),&pmno,&pmnoo,descA);
                                       pzgemm_("n","n",&mno,&mno,&mno,&ONEz,C_loc.data(),&pmno,&pmnoo1,descA,rM.data(),&pmnoo1,&pmnoo,descA,&ONEz,ConsideredMatrix.data(),&pmno,&pmnoo,descA);
                               }
                                   }                   

                           }
                       }

                   }
*/
                   pzgemm_("n","n",&MN,&MN,&N,&ONEz,C_loc.data(),&ONEi,&ONEi,descA,rM.data(),&ONEi,&ONEi,descA,&ONEz,ConsideredMatrix.data(),&ONEi,&ONEi,descA);
		   C_loc.resize(0,0);
}
                   //blacs_barrier_(&ctxt,"A");
                //DataCom(blposX,blposY,MN,MN,NNN,MMM,Nb,Mb,C,C_loc,nrows,ncols,myrow,mycol,0,ctxt);
                //if(id==0) print_matrix(C,"test");
                
                
                //ConsideredMatrix.block(cell->blpo,cell->blpo,cell->AtomO,cell->AtomO).setZero();
                //ConsideredMatrix.block(blpo,blpo,mAtom,mAtom)=diag;
            }
        }

        blacs_gridexit_(&ctxt);
    }

    void Hamiltonian::PFCol3(const Eigen::VectorXcd &egvct,int cellme,int blpo,int AtomO,std::vector<int> vblpo,std::vector<int> vAtomO,const Eigen::Vector3d& k,std::vector<Eigen::Vector3d> chain,int newcolnum,int i,int cc,int kj,Eigen::Matrix3d &Ugm,Eigen::Matrix3d &ugm,Eigen::Matrix3d &uam,Eigen::Matrix<lapack_complex_double,-1,-1> &sMatrix,int &myrow,int &mycol)
    {

        int id,numprocs;
        blacs_pinfo_(&id, &numprocs);
        int procrows = sqrt(numprocs), proccols = sqrt(numprocs);

        //Eigen::MatrixXcd yfitting;
        int UnitCellNumber=systemcell.UnitCellNumber;
        Eigen::Vector3d shift=systemcell.supercell[cellme]->Unit[0];
        assert(vblpo.size()==vAtomO.size());
        int blnum=vblpo.size();

        std::complex<double> ii(0.0,1.0);
        int fitposRow=0;

        for(int j=0;j<blnum;++j)
        {
            Eigen::Vector3d ddd=chain[(i)*systemcell.UnitCellNumber]-shift;
            Eigen::Vector3d AtomPos=chain[(i)*systemcell.UnitCellNumber]-shift-Eigen::Vector3d(USE_Shift,0,0);
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

            dk=ugm*(AtomPos-uam.inverse()*poshift)+Ugm*poshift;//dk1
            int colIndex=cc*orbNum+kj;
            for(int c=0;c<UnitCellNumber;++c)
            {
                Eigen::Vector3d posc=chain[(j)*UnitCellNumber+c]-shift;
                std::complex<double> expr=exp(((dk).dot(posc))*2*M_PI*ii);
                for(int ki=0;ki<orbNum;++ki)
                {
                    int rowIndex=c*orbNum+ki;

                    if(
                            //(fitposCol+colIndex)==colnum
                            USE_RowTruncation==0
                      )
                    {
                        int rows=fitposRow+rowIndex;
                        int cols=newcolnum;
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
                        if(myrow==sendr && mycol == sendc){
                            sMatrix(recvr,recvc).real=real(egvct(rowIndex)*expr/sqrt(blnum));
                            sMatrix(recvr,recvc).imag=imag(egvct(rowIndex)*expr/sqrt(blnum));
                        }
                    }else{

                        if(
                                (fitposRow+rowIndex)>=(LpV)*cNUM
                                && (fitposRow+rowIndex)<(LpV+systemcell.spNum[minDistance+1])*cNUM
                          )
                        {

                            int rows=fitposRow+rowIndex;
                            int cols=newcolnum;
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
                            if(myrow==sendr && mycol == sendc){
                                sMatrix(recvr,recvc).real=real(egvct(rowIndex)*expr/sqrt(blnum));
                                sMatrix(recvr,recvc).imag=imag(egvct(rowIndex)*expr/sqrt(blnum));
                            }
                            //yfitting(fitposRow+rowIndex,newcolnum)=egvct(rowIndex)*expr/sqrt(blnum);//5
                            //yfitting(fitposRow+rowIndex,newcolnum)=egvct(rowIndex)*expr/sqrt(blnum);//5
                        }
                    }

                }
            }
            fitposRow+=vAtomO[j];
        }
    }

void Hamiltonian::PQuickFoldTheBlockHamiltonian(int cellme,int blpo,int AtomO,std::vector<int> vblpo,std::vector<int> vAtomO,const Eigen::Vector3d& k,std::vector<Eigen::Vector3d> chain,std::pair<double,double> EnergyCutoff,int &Startp,int &SecDim,int dtag,Eigen::VectorXd &svalue,Eigen::Matrix<lapack_complex_double,-1,-1> &sMatrix,int &myrow,int &mycol)
{
    int id,numprocs;
    blacs_pinfo_(&id, &numprocs);
    Eigen::MatrixXcd yfitting;
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

    int printM=cellme;
    Eigen::Vector3d shift=systemcell.supercell[cellme]->Unit[0];
    assert(vblpo.size()==vAtomO.size());


    svalue.resize(AtomO);
    svalue.setZero();
    std::vector<std::complex<double>> ord;
    std::complex<double> ii(0.0,1.0);

    int blnum=vblpo.size();
    std::vector<Eigen::MatrixXcd> egvct(blnum);
    int MethodParallel=1;

    if(MethodParallel)
    {    
	int section=(blnum)/numprocs+1;
	//for(int i=0;i<blnum;++i)
	for(int i=std::min(id*section,blnum);i<std::min((id+1)*section,blnum);++i)
	{
	    int blen=vAtomO[i];
	    Eigen::MatrixXcd tmpMatrix;
	    SetMatrixPrime(k,systemcell.ListofMatrixElements_Type3,tmpMatrix,vblpo[i],blen);
	    solver.compute(tmpMatrix,Eigen::DecompositionOptions::ComputeEigenvectors);
	    assert(solver.info() == Eigen::ComputationInfo::Success);

	    svalue.segment(vblpo[i]-vblpo[0],blen)=solver.eigenvalues();
	    egvct[i]=solver.eigenvectors();
	}

	//blacs_barrier_(&ctxt,"A");
	Eigen::MatrixXd tmp2r;
	Eigen::MatrixXd tmp2i;
	for(int i=0;i<egvct.size();++i)
	{
	    int blen=vAtomO[i];

	    //if(id!=i/section) egvct[i].resize(blen,blen);
	    MPI_Bcast(svalue.data()+vblpo[i]-vblpo[0],blen,MPI_DOUBLE,i/section,MPI_COMM_WORLD);

	    tmp2r.resize(blen,blen);
	    tmp2i.resize(blen,blen);
	    if(id==i/section) tmp2r=egvct[i].real();
	    if(id==i/section) tmp2i=egvct[i].imag();
	    MPI_Bcast(tmp2r.data(),blen*blen,MPI_DOUBLE,i/section,MPI_COMM_WORLD);
	    MPI_Bcast(tmp2i.data(),blen*blen,MPI_DOUBLE,i/section,MPI_COMM_WORLD);
	    egvct[i]=tmp2r+tmp2i*ii;
	}
    }else{    
	for(int i=0;i<blnum;++i)
	{
	    int blen=vAtomO[i];
	    Eigen::MatrixXcd tmpMatrix;
	    SetMatrixPrime(k,systemcell.ListofMatrixElements_Type3,tmpMatrix,vblpo[i],blen);
	    solver.compute(tmpMatrix,Eigen::DecompositionOptions::ComputeEigenvectors);
	    assert(solver.info() == Eigen::ComputationInfo::Success);


	    svalue.segment(vblpo[i]-vblpo[0],blen)=solver.eigenvalues();
	    egvct[i]=solver.eigenvectors();
	}
    }


    if(USE_FixedSize_or_EnergyRange)
    {
	//NoteUSE_FixedSize
	if(USE_OneRange_or_MoreRanges)
	{
	    //NoteUSE_OneRange FixedSize
	    if(USE_decomposition && fabs(dtag)!=999)
	    {


		//yfitting.resize(AtomO,leftCutoff+rightCutoff);
		yfitting.resize(AtomO,vAtomO[dtag]);
		yfitting.setZero();
		Startp=vAtomO[dtag]*BandGapIndex/orbNum-leftCutoff;
		SecDim=rightCutoff+leftCutoff;
		//for(int yi=0;yi<yfitting.cols();++yi)
		for(int yi=0;yi<SecDim;++yi)
		{
		    //yfitting.col(yi)=fitting.col(ord[AtomO*BandGapIndex/orbNum+yi-leftCutoff].real());
		    int colnum=dtag*vAtomO[dtag]+vAtomO[dtag]*BandGapIndex/orbNum+yi-leftCutoff;
		    int Fi=colnum/(systemcell.UnitCellNumber*orbNum);
		    int Fcc=colnum%(systemcell.UnitCellNumber*orbNum)/orbNum;
		    int Fkj=colnum%(systemcell.UnitCellNumber*orbNum)%orbNum;
		    PFCol3(egvct[Fi].col(Fcc*orbNum+Fkj),cellme,blpo,AtomO,vblpo,vAtomO, k, chain,Startp+yi,dtag,Fcc,Fkj,Ugm,ugm,uam,sMatrix,myrow,mycol);
		}
		svalue=svalue.segment(vblpo[dtag]-vblpo[0],vAtomO[dtag]);



	    }else{
	    //NoteUSE_OneRange FixedSize

		for(int i=0;i<svalue.size();++i){
		    ord.push_back({i*1.,svalue(i)});
		}
		for(int i=0;i<svalue.size();++i){
		    for(int j=0;j<svalue.size()-i-1;++j){
			if(ord[j].imag()>ord[j+1].imag()){
			    double tmp=svalue(j);
			    svalue(j)=svalue(j+1);
			    svalue(j+1)=tmp;
			    std::complex<double> tmp2=ord[j];
			    ord[j]=ord[j+1];
			    ord[j+1]=tmp2;
			    //fitting.col(j).swap(fitting.col(j+1));
			}
		    }
		}



		Startp=AtomO*BandGapIndex/orbNum-leftCutoff;
		SecDim=rightCutoff+leftCutoff;
		for(int yi=0;yi<SecDim;++yi)
		{
		    int colnum=ord[AtomO*BandGapIndex/orbNum+yi-leftCutoff].real();
		    int Fi=colnum/(systemcell.UnitCellNumber*orbNum);
		    int Fcc=colnum%(systemcell.UnitCellNumber*orbNum)/orbNum;
		    int Fkj=colnum%(systemcell.UnitCellNumber*orbNum)%orbNum;
		    //PFCol3(egvct[Fi].col(Fcc*orbNum+Fkj),cellme,blpo,AtomO,vblpo,vAtomO, k, chain,Startp+yi,Fi,Fcc,Fkj,Ugm,ugm,uam,sMatrix,myrow,mycol);
		    PFCol3(egvct[Fi].col(Fcc*orbNum+Fkj),cellme,blpo,AtomO,vblpo,vAtomO, k, chain,yi,Fi,Fcc,Fkj,Ugm,ugm,uam,sMatrix,myrow,mycol);
		}
	    }
	}else{
	    //NoteUSE_MoreRanges FixedSize
	    if(USE_decomposition && fabs(dtag)!=999)
	    {//place

		//std::cout<<"USE_MoreRanges FixedSize Tag="<<vblpo.size()<<std::endl;
		SecDim=EnergyCutoff.first+EnergyCutoff.second;
		Startp=(vAtomO[dtag]/orbNum)*BandGapIndex-EnergyCutoff.first;
		yfitting.resize(AtomO,vAtomO[dtag]);
		yfitting.setZero();
		for(int yi=0;yi<yfitting.cols();++yi)
		{
		    int colnum=dtag*egvct[dtag].cols()+yi;
		    int Fi=colnum/(systemcell.UnitCellNumber*orbNum);
		    int Fcc=colnum%(systemcell.UnitCellNumber*orbNum)/orbNum;
		    int Fkj=colnum%(systemcell.UnitCellNumber*orbNum)%orbNum;
		    PFCol3(egvct[Fi].col(Fcc*orbNum+Fkj),cellme,blpo,AtomO,vblpo,vAtomO, k, chain,yi,dtag,Fcc,Fkj,Ugm,ugm,uam,sMatrix,myrow,mycol);

		}
		svalue=svalue.segment(vblpo[dtag]-vblpo[0],vAtomO[dtag]);

	    }else{
		if(systemcell.supercell.size()!=1)
		{

		    for(int i=0;i<svalue.size();++i){
			ord.push_back({i*1.,svalue(i)});
		    }
		    for(int i=0;i<svalue.size();++i){
			for(int j=0;j<svalue.size()-i-1;++j){
			    if(ord[j].imag()>ord[j+1].imag()){
				double tmp=svalue(j);
				svalue(j)=svalue(j+1);
				svalue(j+1)=tmp;
				std::complex<double> tmp2=ord[j];
				ord[j]=ord[j+1];
				ord[j+1]=tmp2;
				//fitting.col(j).swap(fitting.col(j+1));
			    }
			}
		    }


		    SecDim=EnergyCutoff.first+EnergyCutoff.second;
		    Startp=(AtomO/orbNum)*BandGapIndex-EnergyCutoff.first;
		    //for(int yi=0;yi<AtomO;++yi)
		    for(int yi=Startp;yi<Startp+SecDim;++yi)
		    {
			int colnum=ord[yi].real();
			int Fi=colnum/(systemcell.UnitCellNumber*orbNum);
			int Fcc=colnum%(systemcell.UnitCellNumber*orbNum)/orbNum;
			int Fkj=colnum%(systemcell.UnitCellNumber*orbNum)%orbNum;
			//PFCol3(egvct[Fi].col(Fcc*orbNum+Fkj),cellme,blpo,AtomO,vblpo,vAtomO, k, chain,yi,Fi,Fcc,Fkj,Ugm,ugm,uam,sMatrix,myrow,mycol);
			PFCol3(egvct[Fi].col(Fcc*orbNum+Fkj),cellme,blpo,AtomO,vblpo,vAtomO, k, chain,yi-Startp,Fi,Fcc,Fkj,Ugm,ugm,uam,sMatrix,myrow,mycol);
		    }
		}else{
		    /*
		       if(MethodIn==1)
		       {
		       if(MethodParallel)
		       {    
		       int section=(blnum)/numprocs+1;
		    //for(int i=0;i<blnum;++i)
		    for(int i=std::min(id*section,blnum);i<std::min((id+1)*section,blnum);++i)
		    {
		    int blen=vAtomO[i];
		    Eigen::MatrixXcd tmpMatrix;
		    SetMatrixPrime(k,systemcell.ListofMatrixElements_Type3,tmpMatrix,vblpo[i],blen);
		    solver.compute(tmpMatrix,Eigen::DecompositionOptions::ComputeEigenvectors);
		    assert(solver.info() == Eigen::ComputationInfo::Success);

		    svalue.segment(vblpo[i]-vblpo[0],blen)=solver.eigenvalues();
		    egvct[i]=solver.eigenvectors();
		    }

		    //blacs_barrier_(&ctxt,"A");
		    Eigen::MatrixXd tmp2r;
		    Eigen::MatrixXd tmp2i;
		    for(int i=0;i<egvct.size();++i)
		    {
		    int blen=vAtomO[i];

		    //if(id!=i/section) egvct[i].resize(blen,blen);
		    MPI_Bcast(svalue.data()+vblpo[i]-vblpo[0],blen,MPI_DOUBLE,i/section,MPI_COMM_WORLD);

		    tmp2r.resize(blen,blen);
		    tmp2i.resize(blen,blen);
		    if(id==i/section) tmp2r=egvct[i].real();
		    if(id==i/section) tmp2i=egvct[i].imag();
		    MPI_Bcast(tmp2r.data(),blen*blen,MPI_DOUBLE,i/section,MPI_COMM_WORLD);
		    MPI_Bcast(tmp2i.data(),blen*blen,MPI_DOUBLE,i/section,MPI_COMM_WORLD);
		    egvct[i]=tmp2r+tmp2i*ii;
		    }
		    }else{    
		    for(int i=0;i<blnum;++i)
		    {
		    int blen=vAtomO[i];
		    Eigen::MatrixXcd tmpMatrix;
		    SetMatrixPrime(k,systemcell.ListofMatrixElements_Type3,tmpMatrix,vblpo[i],blen);
		    solver.compute(tmpMatrix,Eigen::DecompositionOptions::ComputeEigenvectors);
		    assert(solver.info() == Eigen::ComputationInfo::Success);


		    svalue.segment(vblpo[i]-vblpo[0],blen)=solver.eigenvalues();
		    egvct[i]=solver.eigenvectors();
		    }
		    }
		    for(int i=0;i<svalue.size();++i){
		    ord.push_back({i*1.,svalue(i)});
		    }
		    for(int i=0;i<svalue.size();++i){
		    for(int j=0;j<svalue.size()-i-1;++j){
		    if(ord[j].imag()>ord[j+1].imag()){
		    double tmp=svalue(j);
		    svalue(j)=svalue(j+1);
		    svalue(j+1)=tmp;
		    std::complex<double> tmp2=ord[j];
		    ord[j]=ord[j+1];
		    ord[j+1]=tmp2;
		    //fitting.col(j).swap(fitting.col(j+1));
		    }
		    }
		    }
		    }
		     */
		    //std::cout<<"tmpPFTBH - USE_MoreRanges FixedSize ONE"<<std::endl;
		    tmpPFTBH(cellme,blpo,AtomO,vblpo,vAtomO,k,chain,sMatrix,Startp,SecDim,Ugm,ugm,uam,myrow,mycol);
		}
	    }

	}
    }else{
	//NoteUSE_EnergyRange
	if(USE_OneRange_or_MoreRanges)
	{
	    //NoteUSE_OneRange EnergyRange
	    if(USE_decomposition &&fabs(dtag)!=999)
	    {


		//place
		yfitting.resize(AtomO,vAtomO[dtag]);
		yfitting.setZero();
		for(int yi=0;yi<yfitting.cols();++yi)
		{
		    int colnum=dtag*egvct[dtag].cols()+yi;
		    int Fi=colnum/(systemcell.UnitCellNumber*orbNum);
		    int Fcc=colnum%(systemcell.UnitCellNumber*orbNum)/orbNum;
		    int Fkj=colnum%(systemcell.UnitCellNumber*orbNum)%orbNum;
		    PFCol3(egvct[Fi].col(Fcc*orbNum+Fkj),cellme,blpo,AtomO,vblpo,vAtomO, k, chain,yi,dtag,Fcc,Fkj,Ugm,ugm,uam,sMatrix,myrow,mycol);
		}
		svalue=svalue.segment(vblpo[dtag]-vblpo[0],vAtomO[dtag]);
		const auto upper=std::upper_bound(svalue.data(),svalue.data()+svalue.size(),rightCutoff);
		const auto lower=std::lower_bound(svalue.data(),upper,leftCutoff);
		SecDim=std::distance(lower,upper);
		Startp=std::distance(svalue.data(),lower);


	    }else{

		for(int i=0;i<svalue.size();++i){
		    ord.push_back({i*1.,svalue(i)});
		}
		for(int i=0;i<svalue.size();++i){
		    for(int j=0;j<svalue.size()-i-1;++j){
			if(ord[j].imag()>ord[j+1].imag()){
			    double tmp=svalue(j);
			    svalue(j)=svalue(j+1);
			    svalue(j+1)=tmp;
			    std::complex<double> tmp2=ord[j];
			    ord[j]=ord[j+1];
			    ord[j+1]=tmp2;
			    //fitting.col(j).swap(fitting.col(j+1));
			}
		    }
		}



		const auto upper=std::upper_bound(svalue.data(),svalue.data()+svalue.size(),rightCutoff);
		const auto lower=std::lower_bound(svalue.data(),upper,leftCutoff);
		SecDim=std::distance(lower,upper);
		Startp=std::distance(svalue.data(),lower);
		//for(int yi=0;yi<AtomO;++yi)
		for(int yi=Startp;yi<Startp+SecDim;++yi)
		{
		    int colnum=ord[yi].real();
		    int Fi=colnum/(systemcell.UnitCellNumber*orbNum);
		    int Fcc=colnum%(systemcell.UnitCellNumber*orbNum)/orbNum;
		    int Fkj=colnum%(systemcell.UnitCellNumber*orbNum)%orbNum;
		    //PFCol3(egvct[Fi].col(Fcc*orbNum+Fkj),cellme,blpo,AtomO,vblpo,vAtomO, k, chain,yi,Fi,Fcc,Fkj,Ugm,ugm,uam,sMatrix,myrow,mycol);
		    PFCol3(egvct[Fi].col(Fcc*orbNum+Fkj),cellme,blpo,AtomO,vblpo,vAtomO, k, chain,yi-Startp,Fi,Fcc,Fkj,Ugm,ugm,uam,sMatrix,myrow,mycol);
		}
	    }
	}else{
	    //NoteUSE_MoreRanges EnergyRange
	    if(USE_decomposition && fabs(dtag)!=999)
	    {

		//place
		yfitting.resize(AtomO,vAtomO[dtag]);
		yfitting.setZero();
		for(int yi=0;yi<yfitting.cols();++yi)
		{
		    int colnum=dtag*egvct[dtag].cols()+yi;
		    int Fi=colnum/(systemcell.UnitCellNumber*orbNum);
		    int Fcc=colnum%(systemcell.UnitCellNumber*orbNum)/orbNum;
		    int Fkj=colnum%(systemcell.UnitCellNumber*orbNum)%orbNum;
		    PFCol3(egvct[Fi].col(Fcc*orbNum+Fkj),cellme,blpo,AtomO,vblpo,vAtomO, k, chain,yi,dtag,Fcc,Fkj,Ugm,ugm,uam,sMatrix,myrow,mycol);
		}
		svalue=svalue.segment(vblpo[dtag]-vblpo[0],vAtomO[dtag]);
		//std::cout<<"USE_MoreRanges EnergyRange decomposition "<<std::endl;
		const auto upper=std::upper_bound(svalue.data(),svalue.data()+svalue.size(),EnergyCutoff.second);
		const auto lower=std::lower_bound(svalue.data(),upper,EnergyCutoff.first);
		SecDim=std::distance(lower,upper);
		Startp=std::distance(svalue.data(),lower);

		for(int i=0;i<svalue.size();++i){
		    ord.push_back({i*1.,svalue(i)});
		}
		for(int i=0;i<svalue.size();++i){
		    for(int j=0;j<svalue.size()-i-1;++j){
			if(ord[j].imag()>ord[j+1].imag()){
			    double tmp=svalue(j);
			    svalue(j)=svalue(j+1);
			    svalue(j+1)=tmp;
			    std::complex<double> tmp2=ord[j];
			    ord[j]=ord[j+1];
			    ord[j+1]=tmp2;
			    //fitting.col(j).swap(fitting.col(j+1));
			}
		    }
		}




	    }else{
		if(systemcell.supercell.size()!=1)
		{

		    for(int i=0;i<svalue.size();++i){
			ord.push_back({i*1.,svalue(i)});
		    }
		    for(int i=0;i<svalue.size();++i){
			for(int j=0;j<svalue.size()-i-1;++j){
			    if(ord[j].imag()>ord[j+1].imag()){
				double tmp=svalue(j);
				svalue(j)=svalue(j+1);
				svalue(j+1)=tmp;
				std::complex<double> tmp2=ord[j];
				ord[j]=ord[j+1];
				ord[j+1]=tmp2;
				//fitting.col(j).swap(fitting.col(j+1));
			    }
			}
		    }


		    const auto upper=std::upper_bound(svalue.data(),svalue.data()+svalue.size(),EnergyCutoff.second);
		    const auto lower=std::lower_bound(svalue.data(),upper,EnergyCutoff.first);
		    SecDim=std::distance(lower,upper);
		    Startp=std::distance(svalue.data(),lower);
		    //for(int yi=0;yi<AtomO;++yi)
		    for(int yi=Startp;yi<Startp+SecDim;++yi)
		    {
			int colnum=ord[yi].real();
			int Fi=colnum/(systemcell.UnitCellNumber*orbNum);
			int Fcc=colnum%(systemcell.UnitCellNumber*orbNum)/orbNum;
			int Fkj=colnum%(systemcell.UnitCellNumber*orbNum)%orbNum;
			//PFCol3(egvct[Fi].col(Fcc*orbNum+Fkj),cellme,blpo,AtomO,vblpo,vAtomO, k, chain,yi,Fi,Fcc,Fkj,Ugm,ugm,uam,sMatrix,myrow,mycol);
			PFCol3(egvct[Fi].col(Fcc*orbNum+Fkj),cellme,blpo,AtomO,vblpo,vAtomO, k, chain,yi-Startp,Fi,Fcc,Fkj,Ugm,ugm,uam,sMatrix,myrow,mycol);
		    }
		}else{


		    //std::cout<<"PFTBH - USE_MoreRanges EnergyRange ONE NODEcomposition"<<std::endl;
		    PFTBH(cellme,blpo,AtomO,vblpo,vAtomO,k,chain,sMatrix,Startp,SecDim,Ugm,ugm,uam,egvct,svalue,myrow,mycol);

		    for(int i=0;i<svalue.size();++i){
			ord.push_back({i*1.,svalue(i)});
		    }
		    for(int i=0;i<svalue.size();++i){
			for(int j=0;j<svalue.size()-i-1;++j){
			    if(ord[j].imag()>ord[j+1].imag()){
				double tmp=svalue(j);
				svalue(j)=svalue(j+1);
				svalue(j+1)=tmp;
				std::complex<double> tmp2=ord[j];
				ord[j]=ord[j+1];
				ord[j+1]=tmp2;
				//fitting.col(j).swap(fitting.col(j+1));
			    }
			}
		    }
		}
	    }
	}
    }

    if(USE_RowTruncation)
    {
	for(int i=0;i<yfitting.cols();++i)
	{
	    yfitting.col(i)=yfitting.col(i)/sqrt(yfitting.col(i).squaredNorm());
	}
    }
}

    void Hamiltonian::PfoldDiagonalize(std::vector<Eigen::Matrix<lapack_complex_double,-1,-1>> &SubblockContainer,Eigen::MatrixXcd &ConsideredMatrix,bool onoff,Eigen::MatrixXcd &MatrixToBeTransformed,const Eigen::Vector3d& k,int dtag)
    {
        int ZEROi = 0;
        int MONEi=-1; 
        int ONEi=1;    
        double ONEd=1.0;
        double ZEROd=0.0;
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


        matrix_bottom=0;
        SubblockContainer.clear();
        SubblockContainer.resize(0);
        //gsoulFold_vector.clear();
        //gsoulFold_vector.resize(0);
        int Atomtot=0;

        std::vector<Eigen::VectorXd> EigValue;
        EigValue.resize(0);

        if(onoff) transform_matrix.resize(0,0);
        for(auto &cell : systemcell.supercell)
        {

            if(cell->type==0)
            {
                cell->SecDim=cell->AtomO;
                cell->Startp=0;
                Atomtot+=cell->SecDim;
                matrix_bottom+=cell->Startp;
            }else if(cell->Tag3==SubblockContainer.size())
            {


                std::unordered_set<std::string> celltex(cell->UnitTexture.begin(),cell->UnitTexture.end());

                //if(celltex.size()!=1/* && USE_direct_construction==0*/)
                if(0)
                {

                    gsolver.compute(ConsideredMatrix.block(cell->blpo,cell->blpo,cell->AtomO,cell->AtomO),Eigen::MatrixXcd::Identity(cell->AtomO,cell->AtomO) ,Eigen::DecompositionOptions::ComputeEigenvectors|Eigen::Ax_lBx);
                    assert(gsolver.info() == Eigen::ComputationInfo::Success);
                    EigValue.push_back(gsolver.eigenvalues());
                    Filtrate_the_states(gsolver.eigenvalues(),cell,true);

                    if(USE_OneRange_or_MoreRanges)
                    {
                        //SubblockContainer.push_back(gsolver.eigenvectors().block(0,cell->Startp,cell->AtomO,cell->SecDim));
                    }else{
                        //SubblockContainer.push_back(gsolver.eigenvectors());
                    }
                    Atomtot+=cell->SecDim;
                    matrix_bottom+=cell->Startp;
                }else{

                    int N=cell->AtomO; 
                    int M=cell->AtomO; 
                    int nrows = numroc_(&N, &Nb, &myrow, &ZEROi, &procrows);
                    int lda = std::max(1,nrows);
                    int ncols = numroc_(&M, &Mb, &mycol, &ZEROi, &proccols);
                    Eigen::Matrix<lapack_complex_double,-1,-1> sMatrix;
                    sMatrix.resize(nrows,ncols);	
                    sMatrix.setConstant({0,0});

                    Eigen::VectorXd svalue;

                    PQuickFoldTheBlockHamiltonian(cell->Tag2,cell->blpo,cell->AtomO,cell->vblpo,cell->vAtomO,k,cell->Unit,cell->EnergyCutoff,cell->Startp,cell->SecDim,dtag,svalue,sMatrix,myrow,mycol);
                    //print_matrix(sfitting.adjoint()*matrix*sfitting,"testing");

                    SubblockContainer.push_back(sMatrix);
                    EigValue.push_back(svalue);
                    Atomtot+=cell->SecDim;
                    matrix_bottom+=cell->Startp;

                }
                ///////////////////
            }else{
                //Filtrate_the_states(Eigen::VectorXd(),cell,false);
                Filtrate_the_states(EigValue[cell->Tag3],cell,false);

                Atomtot+=cell->SecDim;
                matrix_bottom+=cell->Startp;
            }
        }
        MatrixSizeReduced=Atomtot;

        blacs_gridexit_(&ctxt);
    }

    void Hamiltonian::Diagonal_parallel_pzheevx(Eigen::VectorXd &WW,Eigen::MatrixXcd &ZZ,int blposX,int blposY,int tot,int Nb,int Mb,int NN,int flag,Eigen::Matrix<lapack_complex_double,-1,-1> blockMatrix)
    {
        std::complex<double> ii(0.,1.);
        int ZEROi = 0;
        int MONEi=-1; 
        int ONEi=1;    
        double ONEd=1.0;
        double ZEROd=0.0;
        double ORFAC=-1;
        int MM=-1,NZ=-1; 

        int N=tot; 
        int M=tot; 
        int NNN=NN;
        int MMM=NN;

        WW.resize(N);
        WW.setZero();


        int id,numprocs;
        int ctxt,myrow,mycol;
        blacs_pinfo_(&id, &numprocs);
        blacs_get_(&MONEi, &ZEROi, &ctxt);
        int procrows = sqrt(numprocs), proccols = sqrt(numprocs);
        bool mpiroot = (id == 0);
        double abstol=0;//1e-25;//pdlamch_(&ctxt,"U");
        int info;

        int descA[9];
        int descZ[9];
        blacs_gridinit_(&ctxt, "r", &procrows, &proccols);  
        blacs_gridinfo_(&ctxt, &procrows, &proccols, &myrow, &mycol );  
        int nrows = numroc_(&NNN, &Nb, &myrow, &ZEROi, &procrows);
        int lda = std::max(1,nrows);
        int ncols = numroc_(&MMM, &Mb, &mycol, &ZEROi, &proccols);
        descinit_(descA, &NNN, &MMM, &Nb, &Mb,&ZEROi,&ZEROi,&ctxt, &lda, &info);
        assert(info==0);
        descinit_(descZ, &NNN, &MMM, &Nb, &Mb,&ZEROi,&ZEROi,&ctxt, &lda, &info);
        assert(info==0);
        int sendr = 0, sendc = 0, recvr = 0, recvc = 0;

        double LRWORKOPT;
        int LIWORKOPT; 
        int Eselect=N;
        int *IFAIL=(int*)calloc(N,sizeof(int));
        double *GAP=(double*)calloc(N,sizeof(double));
        int *ICLUSTER=(int*)calloc(N*M,sizeof(int));
        lapack_complex_double LWORKOPT;
        int xONEi=1+blposX;
        int yONEi=1+blposY;

        //lapack_complex_double *A_loc=(lapack_complex_double *)calloc(nrows*ncols,sizeof(lapack_complex_double));
        //lapack_complex_double *Z_loc=(lapack_complex_double *)calloc(nrows*ncols,sizeof(lapack_complex_double));
        Eigen::Matrix<lapack_complex_double,-1,-1> Z_loc;

        if(flag)
        {
            Z_loc.resize(nrows,ncols);
            Z_loc.setConstant({0,0});
        }

        blacs_barrier_(&ctxt,"A");

        if(flag)
            //if(1)
        {
            pzheevx_("v","a","u",&N, blockMatrix.data(), &xONEi , &yONEi, descA, &ZEROd, &ZEROd, &ONEi,&Eselect, &abstol, &MM, &NZ, WW.data(), &ORFAC, Z_loc.data(), &xONEi, &yONEi, descZ,&LWORKOPT,&MONEi,&LRWORKOPT,&MONEi, &LIWORKOPT, &MONEi, IFAIL, ICLUSTER, GAP, &info);
        }else{
            pzheevx_("n","a","u",&N, blockMatrix.data(), &xONEi , &yONEi, descA, &ZEROd, &ZEROd, &ONEi,&Eselect, &abstol, &MM, &NZ, WW.data(), &ORFAC, Z_loc.data(), &xONEi, &yONEi, descZ,&LWORKOPT,&MONEi,&LRWORKOPT,&MONEi, &LIWORKOPT, &MONEi, IFAIL, ICLUSTER, GAP, &info);
        }
        //pzheev_("v","u",&N, A_loc, &xONEi , &yONEi, descA, W,  Z_loc, &xONEi, &yONEi, descZ,&LWORKOPT,&MONEi,&LRWORKOPT,&MONEi, &info);
        //if(id==0) printf("info=%d\n",info);
        assert(info==0);

        blacs_barrier_(&ctxt,"A");
        LIWORKOPT=LIWORKOPT*scaleFactor;
        int LWORK=LWORKOPT.real*scaleFactor;
        int LRWORK=((int)LRWORKOPT)*scaleFactor;

        int *IWORK=(int *)calloc(LIWORKOPT,sizeof(int));
        double *RWORK=(double*)calloc(LRWORK,sizeof(double));
        lapack_complex_double *WORK=(lapack_complex_double*)calloc(LWORK,sizeof(lapack_complex_double));
        if(flag)
            //if(1)
        {
            pzheevx_("v","a","u",&N, blockMatrix.data(), &xONEi , &yONEi, descA, &ZEROd, &ZEROd, &ONEi,&Eselect, &abstol, &MM, &NZ, WW.data(), &ORFAC, Z_loc.data(), &xONEi, &yONEi, descZ, WORK,&LWORK,RWORK,&LRWORK, IWORK, &LIWORKOPT, IFAIL, ICLUSTER, GAP, &info );
        }else{
            pzheevx_("n","a","u",&N, blockMatrix.data(), &xONEi , &yONEi, descA, &ZEROd, &ZEROd, &ONEi,&Eselect, &abstol, &MM, &NZ, WW.data(), &ORFAC, Z_loc.data(), &xONEi, &yONEi, descZ, WORK,&LWORK,RWORK,&LRWORK, IWORK, &LIWORKOPT, IFAIL, ICLUSTER, GAP, &info );
        }
        //pzheev_("v","u",&N, A_loc, &xONEi , &yONEi, descA,  W,  Z_loc, &xONEi, &yONEi, descZ, WORK,&LWORK,RWORK,&LRWORK,  &info );
        //if(id==0) printf("info=%d\n",info);
        assert(info==0);
        blacs_barrier_(&ctxt,"A");

        free(WORK);
        free(RWORK);
        free(IWORK);
        free(IFAIL);
        free(GAP);
        free(ICLUSTER);

        if(flag)
        {
            DataCom(blposX,blposY,N,M,NNN,MMM,Nb,Mb,ZZ,Z_loc,nrows,ncols,myrow,mycol,0,ctxt);
        }
        Z_loc.resize(0,0);

        blacs_gridexit_(&ctxt);


    }

    void Hamiltonian::Diagonal_parallel_pzheevx_all(Eigen::VectorXd &WW,Eigen::MatrixXcd &ZZ,int blposX,int blposY,int tot,int Nb,int Mb,int NN,int flag)
    {
        std::complex<double> ii(0.,1.);
        int ZEROi = 0;
        int MONEi=-1; 
        int ONEi=1;    
        double ONEd=1.0;
        double ZEROd=0.0;
        double ORFAC=-1;
        int MM=-1,NZ=-1; 

        int N=tot; 
        int M=tot; 
        int NNN=NN;
        int MMM=NN;

        WW.resize(N);
        WW.setZero();


        int id,numprocs;
        int ctxt,myrow,mycol;
        blacs_pinfo_(&id, &numprocs);
        blacs_get_(&MONEi, &ZEROi, &ctxt);
        int procrows = sqrt(numprocs), proccols = sqrt(numprocs);
        bool mpiroot = (id == 0);
        double abstol=0;//1e-25;//pdlamch_(&ctxt,"U");
        int info;

        int descA[9];
        int descZ[9];
        blacs_gridinit_(&ctxt, "r", &procrows, &proccols);  
        blacs_gridinfo_(&ctxt, &procrows, &proccols, &myrow, &mycol );  
        int nrows = numroc_(&NNN, &Nb, &myrow, &ZEROi, &procrows);
        int lda = std::max(1,nrows);
        int ncols = numroc_(&MMM, &Mb, &mycol, &ZEROi, &proccols);
        descinit_(descA, &NNN, &MMM, &Nb, &Mb,&ZEROi,&ZEROi,&ctxt, &lda, &info);
        assert(info==0);
        descinit_(descZ, &NNN, &MMM, &Nb, &Mb,&ZEROi,&ZEROi,&ctxt, &lda, &info);
        assert(info==0);
        int sendr = 0, sendc = 0, recvr = 0, recvc = 0;

        double LRWORKOPT;
        int LIWORKOPT; 
        int Eselect=N;
        int *IFAIL=(int*)calloc(N,sizeof(int));
        double *GAP=(double*)calloc(N,sizeof(double));
        int *ICLUSTER=(int*)calloc(N*M,sizeof(int));
        lapack_complex_double LWORKOPT;
        int xONEi=1+blposX;
        int yONEi=1+blposY;

        //lapack_complex_double *A_loc=(lapack_complex_double *)calloc(nrows*ncols,sizeof(lapack_complex_double));
        //lapack_complex_double *Z_loc=(lapack_complex_double *)calloc(nrows*ncols,sizeof(lapack_complex_double));
        Eigen::Matrix<lapack_complex_double,-1,-1> Z_loc;

        if(flag)
        {
            Z_loc.resize(nrows,ncols);
            Z_loc.setConstant({0,0});
        }
        /*
           for(int i=0;i<dMatrix2.rows();++i)
           {
           for(int j=0;j<dMatrix2.cols();++j)
           {
        //A_loc[i+j*nrows].real=dMatrix(i,j).real();
        //A_loc[i+j*nrows].imag=dMatrix(i,j).imag();
        A_loc[i+j*nrows].real=dMatrix2(i,j).real;
        A_loc[i+j*nrows].imag=dMatrix2(i,j).imag;
        }
        }
        */
        blacs_barrier_(&ctxt,"A");

        if(flag)
            //if(1)
        {
            pzheevx_("v","a","u",&N, dMatrix2.data(), &xONEi , &yONEi, descA, &ZEROd, &ZEROd, &ONEi,&Eselect, &abstol, &MM, &NZ, WW.data(), &ORFAC, Z_loc.data(), &xONEi, &yONEi, descZ,&LWORKOPT,&MONEi,&LRWORKOPT,&MONEi, &LIWORKOPT, &MONEi, IFAIL, ICLUSTER, GAP, &info);
        }else{
            pzheevx_("n","a","u",&N, dMatrix2.data(), &xONEi , &yONEi, descA, &ZEROd, &ZEROd, &ONEi,&Eselect, &abstol, &MM, &NZ, WW.data(), &ORFAC, Z_loc.data(), &xONEi, &yONEi, descZ,&LWORKOPT,&MONEi,&LRWORKOPT,&MONEi, &LIWORKOPT, &MONEi, IFAIL, ICLUSTER, GAP, &info);
        }
        //pzheev_("v","u",&N, A_loc, &xONEi , &yONEi, descA, W,  Z_loc, &xONEi, &yONEi, descZ,&LWORKOPT,&MONEi,&LRWORKOPT,&MONEi, &info);
        //if(id==0) printf("info=%d\n",info);
        assert(info==0);

        blacs_barrier_(&ctxt,"A");
        LIWORKOPT=LIWORKOPT*scaleFactor;
        int LWORK=LWORKOPT.real*scaleFactor;
        int LRWORK=((int)LRWORKOPT)*scaleFactor;

        int *IWORK=(int *)calloc(LIWORKOPT,sizeof(int));
        double *RWORK=(double*)calloc(LRWORK,sizeof(double));
        lapack_complex_double *WORK=(lapack_complex_double*)calloc(LWORK,sizeof(lapack_complex_double));
        /*
           for(int i=0;i<nrows;++i)
           {
           for(int j=0;j<i;++j)
           {
        //lapack_complex_double tmp=A_loc[i+j*nrows];
        //A_loc[i+j*nrows]=A_loc[i*ncols+j];
        //A_loc[i*ncols+j]=tmp;
        //A_loc[i+j*nrows].real=0;
        //A_loc[i*ncols+j].real=0;
        }
        }
        for(int i=0;i<nrows;++i)
        {
        for(int j=0;j<i;++j)
        {
        //lapack_complex_double tmp=A_loc[i+j*nrows];
        //A_loc[i+j*nrows]=A_loc[i*ncols+j];
        //A_loc[i*ncols+j]=tmp;
        //A_loc[i+j*nrows].real=0;
        //	A_loc[i*ncols+j].real=0;
        }
        }
        */

        if(flag)
            //if(1)
        {
            pzheevx_("v","a","u",&N, dMatrix2.data(), &xONEi , &yONEi, descA, &ZEROd, &ZEROd, &ONEi,&Eselect, &abstol, &MM, &NZ, WW.data(), &ORFAC, Z_loc.data(), &xONEi, &yONEi, descZ, WORK,&LWORK,RWORK,&LRWORK, IWORK, &LIWORKOPT, IFAIL, ICLUSTER, GAP, &info );
        }else{
            pzheevx_("n","a","u",&N, dMatrix2.data(), &xONEi , &yONEi, descA, &ZEROd, &ZEROd, &ONEi,&Eselect, &abstol, &MM, &NZ, WW.data(), &ORFAC, Z_loc.data(), &xONEi, &yONEi, descZ, WORK,&LWORK,RWORK,&LRWORK, IWORK, &LIWORKOPT, IFAIL, ICLUSTER, GAP, &info );
        }
        //pzheev_("v","u",&N, A_loc, &xONEi , &yONEi, descA,  W,  Z_loc, &xONEi, &yONEi, descZ, WORK,&LWORK,RWORK,&LRWORK,  &info );
        //if(id==0) printf("info=%d\n",info);
        assert(info==0);
        blacs_barrier_(&ctxt,"A");

        free(WORK);
        free(RWORK);
        free(IWORK);
        free(IFAIL);
        free(GAP);
        free(ICLUSTER);

        if(flag)
        {
            DataCom(blposX,blposY,N,M,NNN,MMM,Nb,Mb,ZZ,Z_loc,nrows,ncols,myrow,mycol,0,ctxt);
        }

        Z_loc.resize(0,0);

        blacs_gridexit_(&ctxt);


    }

    void Hamiltonian::Diagonal_parallel_pzheevx_test(Eigen::VectorXd &WW,Eigen::MatrixXcd &ZZ,int blposX,int blposY,int tot,int Nb,int Mb,int NN,int flag)
    {
        std::complex<double> ii(0.,1.);
        int ZEROi = 0;
        int MONEi=-1; 
        int ONEi=1;    
        double ONEd=1.0;
        double ZEROd=0.0;
        double ORFAC=-1;
        int MM=-1,NZ=-1; 

        int N=tot; 
        int M=tot; 
        int NNN=NN;
        int MMM=NN;

        WW.resize(N);
        WW.setZero();


        int id,numprocs;
        int ctxt,myrow,mycol;
        blacs_pinfo_(&id, &numprocs);
        blacs_get_(&MONEi, &ZEROi, &ctxt);
        int procrows = sqrt(numprocs), proccols = sqrt(numprocs);
        bool mpiroot = (id == 0);
        double abstol=0;//1e-25;//pdlamch_(&ctxt,"U");
        int info;

        int descA[9];
        int descZ[9];
        blacs_gridinit_(&ctxt, "r", &procrows, &proccols);  
        blacs_gridinfo_(&ctxt, &procrows, &proccols, &myrow, &mycol );  
        int nrows = numroc_(&NNN, &Nb, &myrow, &ZEROi, &procrows);
        int lda = std::max(1,nrows);
        int ncols = numroc_(&MMM, &Mb, &mycol, &ZEROi, &proccols);
        descinit_(descA, &NNN, &MMM, &Nb, &Mb,&ZEROi,&ZEROi,&ctxt, &lda, &info);
        assert(info==0);
        descinit_(descZ, &NNN, &MMM, &Nb, &Mb,&ZEROi,&ZEROi,&ctxt, &lda, &info);
        assert(info==0);
        int sendr = 0, sendc = 0, recvr = 0, recvc = 0;

        double LRWORKOPT;
        int LIWORKOPT; 
        int Eselect=N;
        int *IFAIL=(int*)calloc(N,sizeof(int));
        double *GAP=(double*)calloc(N,sizeof(double));
        int *ICLUSTER=(int*)calloc(N*M,sizeof(int));
        lapack_complex_double LWORKOPT;
        int xONEi=1+blposX;
        int yONEi=1+blposY;

        //lapack_complex_double *A_loc=(lapack_complex_double *)calloc(nrows*ncols,sizeof(lapack_complex_double));
        //lapack_complex_double *Z_loc=(lapack_complex_double *)calloc(nrows*ncols,sizeof(lapack_complex_double));
        Eigen::Matrix<lapack_complex_double,-1,-1> Z_loc;

        if(flag)
        {
            Z_loc.resize(nrows,ncols);
            Z_loc.setConstant({0,0});
        }
        /*
           for(int i=0;i<dMatrix2.rows();++i)
           {
           for(int j=0;j<dMatrix2.cols();++j)
           {
        //A_loc[i+j*nrows].real=dMatrix(i,j).real();
        //A_loc[i+j*nrows].imag=dMatrix(i,j).imag();
        A_loc[i+j*nrows].real=dMatrix2(i,j).real;
        A_loc[i+j*nrows].imag=dMatrix2(i,j).imag;
        }
        }
        */
        blacs_barrier_(&ctxt,"A");

        if(flag)
            //if(1)
        {
            pzheevx_("v","a","u",&N, dMatrix2.data(), &xONEi , &yONEi, descA, &ZEROd, &ZEROd, &ONEi,&Eselect, &abstol, &MM, &NZ, WW.data(), &ORFAC, Z_loc.data(), &xONEi, &yONEi, descZ,&LWORKOPT,&MONEi,&LRWORKOPT,&MONEi, &LIWORKOPT, &MONEi, IFAIL, ICLUSTER, GAP, &info);
        }else{
            pzheevx_("n","a","u",&N, dMatrix2.data(), &xONEi , &yONEi, descA, &ZEROd, &ZEROd, &ONEi,&Eselect, &abstol, &MM, &NZ, WW.data(), &ORFAC, Z_loc.data(), &xONEi, &yONEi, descZ,&LWORKOPT,&MONEi,&LRWORKOPT,&MONEi, &LIWORKOPT, &MONEi, IFAIL, ICLUSTER, GAP, &info);
        }
        //pzheev_("v","u",&N, A_loc, &xONEi , &yONEi, descA, W,  Z_loc, &xONEi, &yONEi, descZ,&LWORKOPT,&MONEi,&LRWORKOPT,&MONEi, &info);
        //if(id==0) printf("info=%d\n",info);
        assert(info==0);

        blacs_barrier_(&ctxt,"A");
        LIWORKOPT=LIWORKOPT*scaleFactor;
        int LWORK=LWORKOPT.real*scaleFactor;
        int LRWORK=((int)LRWORKOPT)*scaleFactor;

        int *IWORK=(int *)calloc(LIWORKOPT,sizeof(int));
        double *RWORK=(double*)calloc(LRWORK,sizeof(double));
        lapack_complex_double *WORK=(lapack_complex_double*)calloc(LWORK,sizeof(lapack_complex_double));
        /*
           for(int i=0;i<nrows;++i)
           {
           for(int j=0;j<i;++j)
           {
        //lapack_complex_double tmp=A_loc[i+j*nrows];
        //A_loc[i+j*nrows]=A_loc[i*ncols+j];
        //A_loc[i*ncols+j]=tmp;
        //A_loc[i+j*nrows].real=0;
        //A_loc[i*ncols+j].real=0;
        }
        }
        for(int i=0;i<nrows;++i)
        {
        for(int j=0;j<i;++j)
        {
        //lapack_complex_double tmp=A_loc[i+j*nrows];
        //A_loc[i+j*nrows]=A_loc[i*ncols+j];
        //A_loc[i*ncols+j]=tmp;
        //A_loc[i+j*nrows].real=0;
        //	A_loc[i*ncols+j].real=0;
        }
        }
        */

        if(flag)
            //if(1)
        {
            pzheevx_("v","a","u",&N, dMatrix2.data(), &xONEi , &yONEi, descA, &ZEROd, &ZEROd, &ONEi,&Eselect, &abstol, &MM, &NZ, WW.data(), &ORFAC, Z_loc.data(), &xONEi, &yONEi, descZ, WORK,&LWORK,RWORK,&LRWORK, IWORK, &LIWORKOPT, IFAIL, ICLUSTER, GAP, &info );
        }else{
            pzheevx_("n","a","u",&N, dMatrix2.data(), &xONEi , &yONEi, descA, &ZEROd, &ZEROd, &ONEi,&Eselect, &abstol, &MM, &NZ, WW.data(), &ORFAC, Z_loc.data(), &xONEi, &yONEi, descZ, WORK,&LWORK,RWORK,&LRWORK, IWORK, &LIWORKOPT, IFAIL, ICLUSTER, GAP, &info );
        }
        //pzheev_("v","u",&N, A_loc, &xONEi , &yONEi, descA,  W,  Z_loc, &xONEi, &yONEi, descZ, WORK,&LWORK,RWORK,&LRWORK,  &info );
        //if(id==0) printf("info=%d\n",info);
        assert(info==0);
        blacs_barrier_(&ctxt,"A");

        free(WORK);
        free(RWORK);
        free(IWORK);
        free(IFAIL);
        free(GAP);
        free(ICLUSTER);

        if(flag)
        {
            DataCom(blposX,blposY,N,M,NNN,MMM,Nb,Mb,ZZ,Z_loc,nrows,ncols,myrow,mycol,0,ctxt);
        }
        
        Z_loc.resize(0,0);

        blacs_gridexit_(&ctxt);


    }





}
