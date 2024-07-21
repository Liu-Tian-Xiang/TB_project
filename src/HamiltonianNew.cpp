
#include "Hamiltonian.hpp"
namespace TightBinding
{
    void Hamiltonian::MtoS(Eigen::MatrixXcd &M,Eigen::MatrixXcd &S,int flag){
        S.resize(M.rows(),M.cols());

        Eigen::MatrixXcd M11;
        Eigen::MatrixXcd M12;
        Eigen::MatrixXcd M21;
        Eigen::MatrixXcd M22;

        if(flag)
        {
            //M11=M.topRightCorner(M.rows()/2,M.cols()/2);
            //M12=M.topLeftCorner(M.rows()/2,M.cols()/2);
            //M21=M.bottomRightCorner(M.rows()/2,M.cols()/2);
            //M22=M.bottomLeftCorner(M.rows()/2,M.cols()/2);
        }else{
            M11=M.topLeftCorner(M.rows()/2,M.cols()/2);
            M12=M.topRightCorner(M.rows()/2,M.cols()/2);
            M21=M.bottomLeftCorner(M.rows()/2,M.cols()/2);
            M22=M.bottomRightCorner(M.rows()/2,M.cols()/2);
        }
        S.topLeftCorner(S.rows()/2,S.cols()/2)=M21*(M11.inverse());
        S.topRightCorner(S.rows()/2,S.cols()/2)=M22-M21*M11.inverse()*M12;
        S.bottomLeftCorner(S.rows()/2,S.cols()/2)=M11.inverse();
        S.bottomRightCorner(S.rows()/2,S.cols()/2)=-M11.inverse()*M12;
    }
    void Hamiltonian::SxS(Eigen::MatrixXcd &Sa,Eigen::MatrixXcd &Sb){

        Eigen::MatrixXcd Sa11=Sa.topLeftCorner(Sa.rows()/2,Sa.cols()/2);
        Eigen::MatrixXcd Sa12=Sa.topRightCorner(Sa.rows()/2,Sa.cols()/2);
        Eigen::MatrixXcd Sa21=Sa.bottomLeftCorner(Sa.rows()/2,Sa.cols()/2);
        Eigen::MatrixXcd Sa22=Sa.bottomRightCorner(Sa.rows()/2,Sa.cols()/2);

        Eigen::MatrixXcd Sb11=Sb.topLeftCorner(Sb.rows()/2,Sb.cols()/2);
        Eigen::MatrixXcd Sb12=Sb.topRightCorner(Sb.rows()/2,Sb.cols()/2);
        Eigen::MatrixXcd Sb21=Sb.bottomLeftCorner(Sb.rows()/2,Sb.cols()/2);
        Eigen::MatrixXcd Sb22=Sb.bottomRightCorner(Sb.rows()/2,Sb.cols()/2);

        Eigen::MatrixXcd Id=Eigen::MatrixXcd::Identity(Sa.rows()/2,Sa.cols()/2);

        Eigen::MatrixXcd D=Sa12*((Id-Sb11*Sa22).inverse());
        Eigen::MatrixXcd F=Sb21*((Id-Sa22*Sb11).inverse());

        Sb.topLeftCorner(Sb.rows()/2,Sb.cols()/2)=Sa11+D*Sb11*Sa21;
        Sb.topRightCorner(Sb.rows()/2,Sb.cols()/2)=D*Sb12;
        Sb.bottomLeftCorner(Sb.rows()/2,Sb.cols()/2)=F*Sa21;
        Sb.bottomRightCorner(Sb.rows()/2,Sb.cols()/2)=Sb22+F*Sa22*Sb12;
    }
    void Hamiltonian::readEigenVectors(Eigen::MatrixXcd & eigvects,const char* filename,int buffsize){
        int cols=0,rows=0;
        std::ifstream infile;
        infile.open(filename);
        while (!infile.eof())
        {
            std::string line;
            std::getline(infile, line);

            int temp_cols = 0;
            std::stringstream stream(line);
            std::complex<double> bufftmp;
            while(stream >> bufftmp && !stream.eof()){
                //printf("(%1.4f,%1.4f)\t",bufftmp.real(),bufftmp.imag());
                eigvects(rows,temp_cols++)=bufftmp;
            }
            rows++;
        }
        infile.close();
    }

    void Hamiltonian::readEigenValues(Eigen::VectorXd& eigvals,const char* filename){
        int cols=0,rows=0;
        std::ifstream infile;
        infile.open(filename);
        while (!infile.eof())
        {
            std::string line;
            getline(infile, line);

            int temp_cols = 0;
            std::stringstream stream(line);
            while(!stream.eof())
                stream >> eigvals[cols*rows+temp_cols++];

            if (temp_cols == 0)
                continue;

            if (cols == 0)
                cols = temp_cols;

            rows++;
        }
        infile.close();
    }
    void Hamiltonian::store_precomputed_data(const Eigen::Vector3d& kp,int data_index,const Eigen::VectorXd& SAVEDeigenvalues,const Eigen::MatrixXcd& SAVEDeigenvectors)
    {
        //for(int data_index=0;data_index<soul.size();++data_index)
        {
            std::fstream file_obj; 
            char openfile[200];
            sprintf(openfile,"./Data_%1.2f-%1.2f/eigenvalues_%1.5f-%1.5f-%1.5f_%d.dat",leftCutoff,rightCutoff,kp(0),kp(1),kp(2),data_index);
            Eigen::IOFormat CommaInitFmt(Eigen::StreamPrecision, /*Eigen::DontAlignCols*/0, ",", ",", "", "", "", "");
            //Eigen::IOFormat VectFmt(Eigen::StreamPrecision, 0, ",", ",", "", "", "),", ",(");
            file_obj.open(openfile,std::ios::out|std::ios::trunc);
            //file_obj<<soul[data_index]->eigenvalues()/*.format(CommaInitFmt)*/;
            file_obj<<SAVEDeigenvalues;
            file_obj.close();

            char openfile2[200];
            sprintf(openfile2,"./Data_%1.2f-%1.2f/eigenvectors_%1.5f-%1.5f-%1.5f_%d.dat",leftCutoff,rightCutoff,kp(0),kp(1),kp(2),data_index);
            file_obj.open(openfile2,std::ios::out|std::ios::trunc);
            //const Eigen::MatrixXcd &mcd=soul[data_index]->eigenvectors();
            //file_obj<<mcd;
            file_obj<<SAVEDeigenvectors;
            file_obj.close();
        }
    }

    std::string Hamiltonian::store_transform_matrix(const Eigen::Vector3d& kp){

        std::fstream file_obj; 
        char openfile[200];
        sprintf(openfile,"./Data_%1.2f-%1.2f/KPOINT_%1.5f-%1.5f-%1.5f_.dat",leftCutoff,rightCutoff,kp(0),kp(1),kp(2));
        file_obj.open(openfile,std::ios::out|std::ios::trunc);
        file_obj<<transform_mix;
        file_obj.close();
        return std::string(openfile);
    }

    void Hamiltonian::store_materialblocks(std::string name,Eigen::MatrixXcd mat){
        std::fstream file_obj; 
        file_obj.open(name,std::ios::out|std::ios::trunc);
        file_obj<<mat;

/*
        for(int i=0;i<mat.rows();++i){
            for(int j=0;j<mat.cols();++j){
                file_obj<<mat(i,j).real()<<"+";
                file_obj<<mat(i,j).imag()<<"i";
                if(j!=mat.cols()-1) {file_obj<<"\t";}else{file_obj<<std::endl;}
            }
        }
*/
        file_obj.close();
    }

    void Hamiltonian::read_transform_matrix(Eigen::MatrixXcd & eigvects,const char* filename){
        int cols=0,rows=0;
        std::ifstream infile;
        infile.open(filename);
        while (!infile.eof())
        {
            std::string line;
            std::getline(infile, line);

            int temp_cols = 0;
            std::stringstream stream(line);
            std::complex<double> bufftmp;
            while(stream >> bufftmp && !stream.eof()){
                //printf("(%1.4f,%1.4f)\t",bufftmp.real(),bufftmp.imag());
                eigvects(rows,temp_cols++)=bufftmp;
            }
            rows++;
        }
        infile.close();
    }

    void Hamiltonian::GetVK(Eigen::MatrixXcd &vk,Eigen::MatrixXcd &kHvect,Eigen::VectorXcd &kkv,Eigen::Vector3d npl,std::complex<double> Eomega)
    {
        std::complex<double> ii(0.0,1.0);
        Eigen::Vector3d pl={0.5,0.,0.};
        Eigen::MatrixXcd ckHvect2=kHvect;
        ckHvect2.setZero();
        Eigen::MatrixXcd bk=kkv.asDiagonal();
        Eigen::MatrixXcd kval=kkv.asDiagonal();
        Eigen::MatrixXcd kabs=kkv.rowwise().norm().asDiagonal();
        vk=kkv.asDiagonal();vk.setZero();
        for(int i=0;i<kval.rows();++i)
        {
            kval(i,i)=(ii*std::log(kkv(i)))*((pl*2.).norm()/M_PI);
        }
        Eigen::MatrixXcd kWmatPlus,kWmatMinus;
        SetMatrixKW5(0,npl,systemcell.ListofMatrixElements_Type3,kWmatMinus,-1);
        SetMatrixKW5(0,npl,systemcell.ListofMatrixElements_Type3,kWmatPlus,1);
        kWmatMinus=kWmatMinus.block(0,0,kkv.size(),kkv.size());
        kWmatPlus=kWmatPlus.block(0,0,kkv.size(),kkv.size());

        Eigen::MatrixXcd kmat;

        for(int i=0;i<kval.rows();++i)
        {
            std::complex<double> kz=kval(i,i);
            std::complex<double> PeFactor=exp(2*M_PI*(0.25*kz)*ii);
            std::complex<double> MeFactor=exp(-2*M_PI*(0.25*kz)*ii);
            ckHvect2.col(i)=kHvect.block(0,i,kHvect.rows(),1);
            ckHvect2.block(0,i,kHvect.rows()/2,1)*=1;
            ckHvect2.block(kHvect.rows()/2,i,kHvect.rows()/2,1)/=PeFactor;          
            kmat=(kWmatMinus*(2*M_PI*-ii*0.25)*MeFactor+kWmatPlus*(2*M_PI*ii*0.25)*PeFactor);

            std::complex<double> tmp=(ckHvect2.col(i).adjoint()*kmat*ckHvect2.col(i));
            std::complex<double> tmp_overlap=(ckHvect2.col(i).adjoint()*ckHvect2.col(i));

            vk(i,i)=tmp/tmp_overlap;
            //Vk(i,i)=tmp;

            //G0+=-M_PI*ii*CcolStar.segment(0,5)*Ccol.segment(0,5).adjoint()/Vk(i,i);
            //G0+=-M_PI*ii*cvec*cvecStar.adjoint()*(cvecStar.adjoint()*kmat*cvec).diagonal().asDiagonal()*cStarc.inverse();
            //print_matrix_in_details(kval,"kval");
            //print_matrix_in_details(kabs,"kabs");
        }

    }
    void Hamiltonian::SetMatrixKW5(std::complex<double> kz,Eigen::Vector3d pl,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::MatrixXcd &Kmatrix5,int Wsplit)
    {
        std::unordered_map<std::string, Material> &ms=material.materials;
        if(USE_long_parameters==0)
        {
            ms=material.materialsShort;
        }
        int Atom=MatrixSize;
        Kmatrix5.resize(Atom,Atom);	
        Kmatrix5.setZero();
        for(auto &cell : ListofMatrixElements)
        {
            std::complex<double> ii(0.0,1.0);
            Eigen::Vector3d dd=cell->d;
            Eigen::Vector3d ed=cell->ed;
            Eigen::Vector3d edn=ed.normalized();
            Eigen::Vector3d dn=dd.normalized();

            int atom2=cell->Table.real();
            int atom1=cell->Table.imag();
            std::string leftType=(atom1%2==0?"c":"a");

            std::string rightType=(atom2%2==0?"c":"a");
            atom2*=orbNum;
            atom1*=orbNum;

            Eigen::Vector3d d_1_=dd.dot(pl)*pl;
            Eigen::Vector3d d11=dd-d_1_;

            //std::complex<double> efactor=pow(2*M_PI*ii*dd.dot(pl),Wsplit)*exp(2*M_PI*(dd.dot(pl)*kz)*ii);
            std::complex<double> efactor=0;

            if(dd.dot(pl)>0 && Wsplit==1)
            {
                //std::cout<<rows<<">"<<cols<<"---"<<dd.dot(pl)<<std::endl;
                //efactor=pow(2*M_PI*ii*dd.dot(pl),Wsplit)*exp(2*M_PI*(dd.dot(pl)*kz)*ii);
                efactor=1;
            }else if(dd.dot(pl)<0 && Wsplit==-1){
                //std::cout<<rows<<"<"<<cols<<"---"<<dd.dot(pl)<<std::endl;
                //efactor=pow(2*M_PI*ii*dd.dot(pl),Wsplit)*exp(2*M_PI*(dd.dot(pl)*kz)*ii);
                efactor=1;
            }else if(dd.dot(pl)==0 && Wsplit==0){
                efactor=1;
            }

            for(int io=0;io<orbNum;++io)
            {
                for(int jo=0;jo<orbNum;++jo)
                {
                    Eigen::Vector4i leftorbital=qlabels[jo];
                    Eigen::Vector4i rightorbital=qlabels[io];
                    int rows=atom2+io;
                    int cols=atom1+jo;
                    std::string QWMaterial="GaAs";
                    std::string EBMaterial="AlAs";

                    if(leftorbital(3)==rightorbital(3))
                    {
                        int rows=atom2+io;
                        int cols=atom1+jo;
                        if(rows==cols && (Wsplit==0))
                        {
                            if(cell->cellme==QWMaterial && cell->SuperCell==QWMaterial)
                            {
                                Kmatrix5(rows,cols)=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType)+BandOffSet;
                            }else{
                                Kmatrix5(rows,cols)=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType);
                            }

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

                            Kmatrix5(rows,cols)+=res*efactor;

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
                                Kmatrix5(rows,cols)+=ii*(0.5*rightorbital(2))*lambda;
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

                            Kmatrix5(rows,cols)+=lambda*soc;
                        }

                    }

                }
            }
        }


    }

    
    void Hamiltonian::CaseSurfaceDOS(std::complex<double> Eomega,Eigen::Vector3d k)
    {

        std::complex<double> ii(0.0,1.0);
        int bksize=systemcell.orbNum;
        int scsize=systemcell.supercell.size()*2.0;
        int hfM=bksize*2;
        Eigen::MatrixXcd tt2,rr1,rr2,rr3,tt3,tt1;
        Eigen::MatrixXcd tt4,rr4,tt,rr,Tac,Tca;
        Eigen::MatrixXcd J0,JN;
        Eigen::MatrixXcd matTOT;
        Eigen::MatrixXcd matII;
        std::vector<Eigen::MatrixXcd> vT,vT2,vW,vJ,vG;

        Eigen::Vector3d kk={0.,0.,0.};
        Eigen::Vector3d pl={1,1,1};
        //Eigen::Vector3d pl={1.0,0.,0.};
        //Eigen::Vector3d pl={1,1,0.};
        Eigen::Vector3d npl=pl.normalized();

        Eigen::MatrixXcd G0;
        //GenH(Eomega,matTOT,bksize,kk,npl);
        //eGenVT(matTOT,vT,vT2,vW,vG,bksize,scsize);

        Eigen::MatrixXcd kH;//=vT[1]*vT[0];

        if(pl==Eigen::Vector3d(1,1,0))
        //if(0)
        {
            GenH2(Eomega,matTOT,bksize*2,scsize/2,kk,npl);
            eGenVT(matTOT,vT,vT2,vJ,vG,bksize*2,scsize/2);
            kH=vT.front(); 
        }else{
            GenH(Eomega,matTOT,vW,bksize,kk,npl);
            eGenVT(matTOT,vT,vT2,vJ,vG,bksize,scsize);
            kH=vT.back()*vT.front();
        }


        csolver.compute(kH,Eigen::DecompositionOptions::ComputeEigenvectors);
        Eigen::VectorXcd kkv=csolver.eigenvalues();
        Eigen::MatrixXcd kHvect=csolver.eigenvectors();
        csolver.compute(kH.adjoint(),Eigen::DecompositionOptions::ComputeEigenvectors);
        Eigen::VectorXcd Skkv=csolver.eigenvalues();
        Eigen::MatrixXcd SkHvect=csolver.eigenvectors();
        //SkHvect=kHvect.inverse().adjoint();

        //Dreverse(kkv,kHvect);
        Eigen::MatrixXcd tmpkHvect=kHvect;
        Eigen::MatrixXcd ckHvect=kHvect;
        Eigen::MatrixXcd kHvect2=kHvect;
        Eigen::MatrixXcd kHvectStar=kHvect;
        Eigen::MatrixXcd kHvectLong(kHvect.rows()*2,kHvect.cols());
        Eigen::MatrixXcd kHvectLongStar(kHvect.rows()*2,kHvect.cols());
        kHvect2.setZero();
        tmpkHvect.setZero();
        ckHvect.setZero();
        kHvectStar.setZero();
        kHvectLong.setZero();
        kHvectLongStar.setZero();
        Eigen::MatrixXcd ckHvect2=kHvect2;
        Eigen::MatrixXcd longCkHvect;
        Eigen::MatrixXcd longCkHvect2;
        longCkHvect.resize(kHvect2.rows()*2,kHvect2.cols());
        longCkHvect2.resize(kHvect2.rows()*2,kHvect2.cols());
        longCkHvect.setZero();
        longCkHvect2.setZero();
        //ArrangeDRE(kkv,kHvect);
        Eigen::MatrixXcd bk=kkv.asDiagonal();
        Eigen::MatrixXcd kval=kkv.asDiagonal();
        Eigen::MatrixXcd kabs=kkv.rowwise().norm().asDiagonal();
        Eigen::MatrixXcd Vk=kkv.asDiagonal();Vk.setZero();
        //for(int i=0;i<kH.rows();++i)
        {
            //kval(i,i)=(ii*std::log(kkv(i)))*((pl*2.).norm()/M_PI);
            
            //kval(i,i)=(ii*std::log(kkv(i)))*((npl).norm()/M_PI);
        }
        G0=Vk;
        Eigen::MatrixXcd kWmatPlus,kWmatMinus;
        SetMatrixKW5(0,npl,systemcell.ListofMatrixElements_Type3,kWmatMinus,-1);
        SetMatrixKW5(0,npl,systemcell.ListofMatrixElements_Type3,kWmatPlus,1);
        kWmatMinus=kWmatMinus.block(0,0,10,10);
        kWmatPlus=kWmatPlus.block(0,0,10,10);
/*
        int enlargefactor=2;
        vW[0].resize(bksize*2*enlargefactor,bksize*2*enlargefactor);
        vW[0].setZero();
        vW[0].block(0,bksize*enlargefactor,bksize*enlargefactor,bksize*enlargefactor)=kWmatPlus;
        vW[0].block(bksize*enlargefactor,0,bksize*enlargefactor,bksize*enlargefactor)=kWmatMinus;
*/ 

        Eigen::MatrixXcd kmat;
        for(int i=0;i<kH.rows();++i)
        {
            //std::complex<double> kz=(ii*std::log(kkv(i)))*((pl*2.).norm()/M_PI);
            //std::complex<double> PeFactor=exp(2*M_PI*(0.25*kz)*ii);
            //std::complex<double> MeFactor=exp(-2*M_PI*(0.25*kz)*ii);
            std::complex<double> PeFactor=1.0/sqrt(kkv(i));
            std::complex<double> MeFactor=sqrt(kkv(i));
            std::complex<double> tmp=0;
            std::complex<double> tmp_overlap=1;
            if(1 /*&& abs(abs(kkv(i))-1)<1e-3*/)
            {
                //kmat=(kWmatMinus*(-2*M_PI*ii*0.25)*MeFactor+kWmatPlus*(2*M_PI*ii*0.25)*PeFactor);
                //Eigen::MatrixXcd IvW=vW[0];
                //IvW.block(0,bksize*enlargefactor,bksize*enlargefactor,bksize*enlargefactor)*=PeFactor;
                //IvW.block(bksize*enlargefactor,0,bksize*enlargefactor,bksize*enlargefactor)*=MeFactor;
if(0)       
{

                ckHvect2.col(i)=kHvect.block(0,i,kHvect2.rows(),1);
                ckHvect2.block(0,i,kHvect2.rows()/2,1)*=1;
                ckHvect2.block(kHvect2.rows()/2,i,kHvect2.rows()/2,1)*=MeFactor;          
/*
                ckHvect.col(i)=SkHvect.block(0,i,SkHvect.rows(),1);
                ckHvect.block(0,i,SkHvect.rows()/2,1)*=1;
                ckHvect.block(SkHvect.rows()/2,i,SkHvect.rows()/2,1)*=MeFactor;          
*/

                longCkHvect.block(0,i,ckHvect2.col(i).size(),1)=ckHvect2.col(i);
                longCkHvect2.block(0,i,ckHvect2.col(i).size(),1)=ckHvect2.col(i);
                longCkHvect.block(ckHvect2.col(i).size(),i,ckHvect2.col(i).size(),1)=ckHvect2.col(i)*(ii)*std::conj(MeFactor);
                longCkHvect2.block(ckHvect2.col(i).size(),i,ckHvect2.col(i).size(),1)=ckHvect2.col(i)*(ii)/MeFactor;
                //std::complex<double> Itmp=(ckHvect2.col(i).adjoint()*kmat*ckHvect2.col(i));
                //std::complex<double> longItmp=(longCkHvect.adjoint()*vW[0]*longCkHvect2)(i,i);
                //std::complex<double> Itmp_overlap=(ckHvect2.col(i).adjoint()*ckHvect2.col(i));
                //tmp=longItmp;
                //tmp_overlap=Itmp_overlap;
}else{
    /*
                std::complex<double> SMeFactor=1;
                ckHvect2.col(i)=kHvect.block(0,i,kHvect2.rows(),1);
                //ckHvect2.block(0,i,kHvect2.rows()/2,1)*=1;
                //ckHvect2.block(kHvect2.rows()/2,i,kHvect2.rows()/2,1)*=ii;
                ckHvect.col(i)=SkHvect.block(0,i,SkHvect.rows(),1);
                //ckHvect.block(0,i,SkHvect.rows()/2,1)*=1;
                //ckHvect.block(SkHvect.rows()/2,i,SkHvect.rows()/2,1)*=ii;          
                longCkHvect.block(0,i,ckHvect2.col(i).size(),1)=ckHvect.col(i);
                longCkHvect.block(ckHvect2.col(i).size(),i,ckHvect2.col(i).size(),1)=ckHvect.col(i)*ii;
                longCkHvect2.block(0,i,ckHvect2.col(i).size(),1)=ckHvect2.col(i);
                longCkHvect2.block(ckHvect2.col(i).size(),i,ckHvect2.col(i).size(),1)=ckHvect2.col(i)*ii;
               */ 
                //std::complex<double> longItmp=(longCkHvect.adjoint()*vW[0]*longCkHvect2)(i,i);
                //std::complex<double> longItmp=(ckHvect2.col(i).adjoint()*vW[0].topRightCorner(vW[0].rows()/2,vW[0].cols()/2)*ckHvect2.col(i));
                //std::complex<double> longItmp=(ckHvect2.col(i).adjoint()*vW[0].bottomLeftCorner(vW[0].rows()/2,vW[0].cols()/2)*ckHvect2.col(i));
                //std::complex<double> Itmp_overlap=(ckHvect2.col(i).adjoint()*ckHvect2.col(i));

                //tmp=(ckHvect2.adjoint()*vW[0].bottomLeftCorner(vW[0].rows()/2,vW[0].cols()/2)*ckHvect2)(i,i);
                //tmp=longItmp;
                //tmp_overlap=Itmp_overlap;
}
                
            }
            if(0)
            {
                ckHvect2.col(i)=kHvect.block(0,i,kHvect2.rows(),1);
                ckHvect2.block(0,i,kHvect2.rows()/2,1)*=1;
                ckHvect2.block(kHvect2.rows()/2,i,kHvect2.rows()/2,1)/=PeFactor;          
                kmat=(kWmatMinus*(-2*M_PI*ii*0.25)*MeFactor+kWmatPlus*(2*M_PI*ii*0.25)*PeFactor);
                std::complex<double> Itmp=(ckHvect2.col(i).adjoint()*kmat*ckHvect2.col(i));
                std::complex<double> Itmp_overlap=(ckHvect2.col(i).adjoint()*ckHvect2.col(i));
                tmp=Itmp;
                tmp_overlap=Itmp_overlap;
            }
            if(0)
            {
                Eigen::MatrixXcd trix;
                std::complex<double> kz=(ii*std::log(kkv(i)))*((pl*2.).norm()/M_PI);
                SetMatrixA(k,kz,npl,systemcell.ListofMatrixElements_Type3,trix);
                trix=trix.block(0,0,10,10);
                print_matrix(trix*ckHvect2.col(i)-Eomega*ckHvect2.col(i),"HC-eC");
                print_matrix(kH*kHvect.col(i)-kkv(i)*kHvect.col(i),"TC-kC");
            }
            
            //print_matrix(kmat,"kkmat");
            //print_matrix(vW[0],"vW[0]");
            

            //Vk(i,i)=tmp/tmp_overlap;

            //Vk(i,i)=tmp;

            //G0+=-M_PI*ii*CcolStar.segment(0,5)*Ccol.segment(0,5).adjoint()/Vk(i,i);
            //G0+=-M_PI*ii*cvec*cvecStar.adjoint()*(cvecStar.adjoint()*kmat*cvec).diagonal().asDiagonal()*cStarc.inverse();
            //print_matrix_in_details(kval,"kval");
            //print_matrix_in_details(kabs,"kabs");
        }
        //Vk=((longCkHvect.adjoint()*vW[0]*longCkHvect2).diagonal()).asDiagonal();
        //Vk=(kHvect.inverse()*(vW[0].topRightCorner(vW[0].rows()/2,vW[0].cols()/2)-vW[0].bottomLeftCorner(vW[0].rows()/2,vW[0].cols()/2))*kHvect)*ii;
        
        Eigen::MatrixXcd tmpvW=vW[0].topRightCorner(vW[0].rows()/2,vW[0].cols()/2);
        Eigen::MatrixXcd MtmpvW=tmpvW.bottomLeftCorner(tmpvW.rows()/2,tmpvW.cols()/2)-tmpvW.topRightCorner(tmpvW.rows()/2,tmpvW.cols()/2);
        
        Vk=(kHvect.inverse()*(tmpvW-tmpvW.adjoint())*kHvect)*ii;
        Eigen::MatrixXcd X=kHvect*(Vk.diagonal().asDiagonal())*kHvect.inverse();

        Vk=kHvect.inverse()*X*kHvect;
        print_matrix_in_details(Vk,"Vk");
        print_matrix_in_details(kHvect.inverse()*kHvect,"kH*kH");
        print_matrix_in_details(kHvect.inverse()*kH*kHvect,"kH*kH*kH");

        //print_matrix(longCkHvect.adjoint()*vW[0]*longCkHvect2-Vk,"test3");
        //print_matrix(Vk,"Vk");

        //Vk=(ckHvect2.transpose()*ckHvect2);
        //print_matrix(kH*kHvect-kHvect*bk,"Tuv=k*uv");
        //print_matrix(Vk,"Vk");

        //MatrixO(kHvect2);
        //MatrixO(kHvectStar);
        int Vsize=Vk.rows();
        Eigen::MatrixXcd Vk1=Vk.block(0,0,Vsize/2,Vsize/2).reverse();
        Eigen::MatrixXcd Vk2=Vk.block(Vsize/2,Vsize/2,Vsize/2,Vsize/2);
        Eigen::MatrixXcd bk1=bk.block(0,0,Vsize/2,Vsize/2).reverse();
        Eigen::MatrixXcd bk2=bk.block(Vsize/2,Vsize/2,Vsize/2,Vsize/2);
        std::complex<double> Vtrace=0;
/*
        Eigen::MatrixXcd M=kHvect;
        M.block(0,M.cols()/2,M.rows(),M.cols()/2)=-kH*M.block(0,M.cols()/2,M.rows(),M.cols()/2);

        Eigen::MatrixXcd Inject=Eigen::MatrixXcd::Identity(kH.rows(),kH.rows()/2);
        //Inject.block(kH.rows()/2,0,kH.rows()/2,1).setZero();
        //Inject.block(0,0,kH.rows()/2,1).setZero();
        print_matrix(Inject,"inject");

        Eigen::MatrixXcd mtr=M.inverse()*kH*kHvect*Inject;
        print_matrix(mtr,"mtr");
        print_matrix(mtr.adjoint()*mtr,"mtr");
        for(int i=0;i<kH.rows()/2;++i)
        {
            //Vtrace+=std::conj(bk(i,i))*bk(i,i)*Vk(i,i);
            //Vtrace+=std::conj(tr(i))*tr(i)*Vk(i,i);
            Vtrace+=std::conj(mtr(i,0))*mtr(i,0)*Vk(i,i);
            mtr(i,i)=std::conj(mtr(i,0))*mtr(i,0)*Vk(i,i)/Vk(i,i);
        }

        std::cout<<std::endl;
        std::cout<<std::endl;
        std::cout<<std::endl;
        std::cout<<std::fixed<<std::setprecision(5)
            <<"                                   vtrace="<<Vtrace<<std::endl;
        std::cout<<std::endl;
        std::cout<<std::endl;
        std::cout<<std::endl;
*/
        //print_matrix(bk1,"bk");
        //print_matrix(bk2,"|bk|");

        //print_matrix(Vk1,"Vk1");
        //print_matrix(Vk2,"Vk2");
        //print_matrix(bk*Vk,"DOT");

        //print_matrix_in_details(Vk.block(0,0,5,5)+Vk.block(5,5,5,5),"3");
        //print_matrix_in_details(Vk.block(0,0,5,5)-Vk.block(5,5,5,5),"4");
        //print_matrix_in_details(Vk,"Vk");
        
        //G0=-M_PI*ii*(kHvect2.block(0,0,Vsize/2,Vsize/2)*kHvectStar.block(0,0,Vsize/2,Vsize/2).adjoint()-kHvect2.block(Vsize/2,Vsize/2,Vsize/2,Vsize/2).adjoint()*kHvectStar.block(Vsize/2,Vsize/2,Vsize/2,Vsize/2));
        //G0=-M_PI*ii*(kHvect2.block(0,0,Vsize,Vsize/2)*Vk1*kHvect2.block(0,0,Vsize,Vsize/2).adjoint());
        //G0=-M_PI*ii*(kHvectStar.block(0,0,Vsize/2,Vsize/2)*Vk1*kHvect2.block(0,0,Vsize/2,Vsize/2).adjoint()-kHvectStar.block(Vsize/2,Vsize/2,Vsize/2,Vsize/2)*Vk2*kHvect2.block(Vsize/2,Vsize/2,Vsize/2,Vsize/2).adjoint());
        
        //G0=-M_PI*ii*Vk1.inverse();
        G0=-M_PI*ii*Vk1;
        //G0=Vk1;
        Eigen::MatrixXcd eII;


        std::fstream file_obj; 
        file_obj.open("./Cdata.txt",std::ios::out|std::ios::app);
        file_obj<<Eomega.real();
        file_obj<<"\t"<<abs((G0.trace()).imag());
        //file_obj<<"\t"<<((G0.trace()).imag());
        //file_obj<<"\t"<<-((G0.trace()).imag());
        //file_obj<<"\t"<<abs((G0.trace()));

        //file_obj<<"\t"<<-M_PI*(G0.trace()).imag();
        //file_obj<<"\t"<<abs(G0.determinant().imag());

        //file_obj<<"\t"<<abs(-M_PI*(vG[0].trace()).imag());
        //file_obj<<"\t"<<abs(vG[0].determinant().real());

        file_obj<<std::endl;
        file_obj.close();
    }

}
