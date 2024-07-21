
#include "Hamiltonian.hpp"
namespace TightBinding
{

    void Hamiltonian::SetMatrixSparse2(const Eigen::Vector3d k,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::MatrixXcd &matrix)
    {
        std::unordered_map<std::string, Material> &ms=material.materials;
        if(USE_long_parameters==0)
        {
        ms=material.materialsShort;
        }
        SparseMatrix.clear();
        std::string QWMaterial="GaAs";
        std::string EBMaterial="AlAs";
        for(auto &cell : ListofMatrixElements)
        {
            std::complex<double> ii(0.0,1.0);
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
                    if(leftorbital(3)==rightorbital(3))
                    {
                        int rows=atom2+io;
                        int cols=atom1+jo;
                        std::complex<int> Index(rows,cols);
                        if(rows==cols)
                        {
                            double res=0;
                            if(cell->cellme==QWMaterial && cell->SuperCell==QWMaterial)
                            {
                                res=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType)+BandOffSet;
                            }else{
                                res=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType);
                            }
                            if(SparseMatrix.count(Index)){
                                SparseMatrix[Index].real+=res;
                                SparseMatrix[Index].imag=0;
                            }else{
                                SparseMatrix[Index].real=res;
                                SparseMatrix[Index].imag=0;
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
                            //matrix(rows,cols)+=res*efactor;

                            if(SparseMatrix.count(Index)){
                                SparseMatrix[Index].real+=real(res*efactor);
                                SparseMatrix[Index].imag+=imag(res*efactor);
                            }else{
                                SparseMatrix[Index].real=real(res*efactor);
                                SparseMatrix[Index].imag=imag(res*efactor);
                            }

                        }
                        if(
                                atom2==atom1 
                                && leftorbital(1)==rightorbital(1) 
                                && leftorbital(1)==1
                                //&&USE_SOC
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
                        std::complex<int> Index(rows,cols);

                        double lambda=1;
                        double m1=leftorbital(2);
                        double m2=rightorbital(2);
                        double leftsign=leftorbital(3);
                        double rightsign=rightorbital(3);
                        if(
                                atom2==atom1 
                                && leftorbital(1)==rightorbital(1) 
                                && leftorbital(1)==1
                                //&&USE_SOC
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


    void Hamiltonian::SetMatrixSparse(const Eigen::Vector3d k,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::MatrixXcd &matrix)
    {
        std::unordered_map<std::string, Material> &ms=material.materials;
        if(USE_long_parameters==0)
        {
        ms=material.materialsShort;
        }
        SparseMatrix.clear();
        for(auto &cell : ListofMatrixElements)
        {
            std::complex<double> ii(0.0,1.0);
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
                    std::complex<int> Index(rows,cols);
                    double BandOffSet=0.;
                    std::string QWMaterial="GaAs";
                    std::string EBMaterial="AlAs";
                    //double BandOffSet=0.46*0.5;
                    //double BandOffSet=0.46;
                    if(rows==cols)
                    {
                        if((cell->cellme==QWMaterial && cell->SuperCell==EBMaterial) || (cell->cellme==EBMaterial && cell->SuperCell==QWMaterial))
                        {
                            //matrix(rows,cols)=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType)+BandOffSet*0.5;
                        }
                        else if(cell->cellme==QWMaterial && cell->SuperCell==QWMaterial)
                        {
                            //matrix(rows,cols)=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType)+BandOffSet;
                        }else{
                            //matrix(rows,cols)=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType);
                        }
                    }
                    else if(leftType!=rightType && rows>cols)
                    { 
                        double res=0;
                        res=(SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->cellme)));

                        /*
                           if((cell->cellme==QWMaterial && cell->SuperCell==EBMaterial)){
                           res=(SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->cellme)));
                           }else if((cell->cellme==EBMaterial && cell->SuperCell==QWMaterial)){
                           res=(SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->SuperCell)));
                           }else{
                           res=(SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->cellme)));
                           }                       
                           */
                        //matrix(rows,cols)+=res*efactor;
                        if(SparseMatrix.count(Index)){
                            SparseMatrix[Index].real+=real(res*efactor);
                            SparseMatrix[Index].imag+=imag(res*efactor);
                        }else{
                            SparseMatrix[Index].real=real(res*efactor);
                            SparseMatrix[Index].imag=imag(res*efactor);
                        }

                    }
                }
            }
        }
    }


    void Hamiltonian::ShiftMatrix(Eigen::MatrixXcd &Omatrix,Eigen::MatrixXcd &Smatrix)
    {
        Smatrix=Omatrix;
        Smatrix.setZero();
        for(auto &cell : systemcell.supercell)
        {
            Eigen::MatrixXcd blockMatrix=Omatrix.block(cell->blpo,cell->blpo,cell->AtomO,cell->AtomO);
            Smatrix.block(cell->blpo,0,cell->AtomO,cell->AtomO)=blockMatrix;
        }
    }

    void Hamiltonian::MatrixO(Eigen::MatrixXcd &matrix)
    {
        for(int i=1;i<matrix.cols();++i)
        {
            Eigen::VectorXcd cl=matrix.col(i);
            for(int j=0;j<i;++j)
            {
                Eigen::VectorXcd il=matrix.col(j);
                std::complex<double> pj=il.adjoint()*cl;
                std::complex<double> jj=il.adjoint()*il;
                cl=cl-(pj/jj)*il;
                //cl=cl.normalized();
            }
            matrix.col(i)=cl.normalized();
        }
    }
    void Hamiltonian::SetMatrixA(const Eigen::VectorXcd k,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::MatrixXcd &matrix)
    {
        std::unordered_map<std::string, Material> &ms=material.materials;
        if(USE_long_parameters==0)
        {
        ms=material.materialsShort;
        }
        int Atom=MatrixSize;
        matrix.resize(Atom,Atom);	
        matrix.setZero();
        std::string QWMaterial="GaAs";
        std::string EBMaterial="AlAs";

        for(auto &cell : ListofMatrixElements)
        {
            /*
               if(cell->cellme=="GaAs") cell->cellme="AlAs";
               else if(cell->cellme=="AlAs") cell->cellme="GaAs";
               if(cell->SuperCell=="GaAs") cell->SuperCell="AlAs";
               else if(cell->SuperCell=="AlAs") cell->SuperCell="GaAs";
               */
            std::complex<double> ii(0.0,1.0);
            Eigen::Vector3d dd=cell->d;
            Eigen::Vector3d ed=cell->ed;
            Eigen::Vector3d edn=ed.normalized();
            Eigen::Vector3d dn=dd.normalized();
            //std::cout<<"dd="<<dd<<std::endl;
            //double phaseFactor=(2*M_PI*dd).dot(k);
            //std::complex<double> efactor=exp(phaseFactor*ii);
            std::complex<double> efactor=exp(2.*M_PI*(dd(0)*k(0)+dd(1)*k(1)+dd(2)*k(2))*ii);




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
                    if(leftorbital(3)==rightorbital(3))
                    {
                        int rows=atom2+io;
                        int cols=atom1+jo;
                        if(rows==cols)
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
                        matrix(rows,rows)=
                            (ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType))*0.5
                            +(ms.at(cell->SuperCell).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType))*0.5
                            +BandOffSet*0.5;

            }
        }
*/
        //print_matrix(matrix,"");
    }

    void Hamiltonian::SetMatrixA(Eigen::Vector3d kk,std::complex<double> kz,Eigen::Vector3d pl,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::MatrixXcd &SmatrixA)
    {
        std::unordered_map<std::string, Material> &ms=material.materials;
        if(USE_long_parameters==0)
        {
        ms=material.materialsShort;
        }
        int Atom=MatrixSize;
        SmatrixA.resize(Atom,Atom);	
        SmatrixA.setZero();
        for(auto &cell : ListofMatrixElements)
        {
            std::complex<double> ii(0.0,1.0);
            Eigen::Vector3d dd=cell->d;
            Eigen::Vector3d ed=cell->ed;
            Eigen::Vector3d edn=ed.normalized();
            Eigen::Vector3d dn=dd.normalized();

            Eigen::Vector3d d_1_=dd.dot(pl)*pl;
            Eigen::Vector3d d11=dd-d_1_;
            std::complex<double> efactor=exp(2*M_PI*(dd.dot(pl)*kz)*ii);



            //std::complex<double> efactor=exp(2*M_PI*(dd.dot(kk))*ii);



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
                    double BandOffSet=0.;
                    std::string QWMaterial="GaAs";
                    std::string EBMaterial="AlAs";
                    //double BandOffSet=0.46*0.5;
                    //double BandOffSet=0.35;
                    if(rows==cols)
                    {
                        if((cell->cellme==QWMaterial && cell->SuperCell==EBMaterial) || (cell->cellme==EBMaterial && cell->SuperCell==QWMaterial))
                        {
                            SmatrixA(rows,cols)=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType)+BandOffSet*0.5;
                        }
                        else if(cell->cellme==QWMaterial && cell->SuperCell==QWMaterial)
                        {
                            SmatrixA(rows,cols)=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType)+BandOffSet;
                        }else{
                            SmatrixA(rows,cols)=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType);
                        }
                    }
                    else if(leftType!=rightType)
                    { 
                        double res=(SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->cellme))*0.5+SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->SuperCell))*0.5);
                        SmatrixA(rows,cols)+=res*efactor;

                    }
                }
            }
        }
    }

    

    void Hamiltonian::SetMatrixK5(const Eigen::MatrixXcd &Zmatrix,std::complex<double> kz,Eigen::Vector3d pl,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::MatrixXcd &Kmatrix5,int Wsplit)
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

            std::complex<double> efactor=pow(2*M_PI*ii*dd.dot(pl),Wsplit)*exp(2*M_PI*(dd.dot(pl)*kz)*ii);
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
                        if(rows==cols)
                        {
                        
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

        //Kmatrix5=Zmatrix.adjoint()*transform_matrix.adjoint()*Kmatrix5*transform_matrix*Zmatrix;
        //Kmatrix5=transform_matrix.adjoint()*Kmatrix5*transform_matrix;
        Kmatrix5=Zmatrix.adjoint()*Kmatrix5*Zmatrix;

    }

    void Hamiltonian::GenH2(std::complex<double> Eomega,Eigen::MatrixXcd &emat,int bksize,int scsize,Eigen::Vector3d kk,Eigen::Vector3d npl)
    {

        std::complex<double> ii(0.0,1.0);
        Eigen::MatrixXcd matTOT,mat0,matPlus,matMinus;

        SetMatrix5(kk,npl,systemcell.ListofMatrixElements_Type3,mat0,0);
        SetMatrix5(kk,npl,systemcell.ListofMatrixElements_Type3,matPlus,1);
        SetMatrix5(kk,npl,systemcell.ListofMatrixElements_Type3,matMinus,-1);
        

        ReadCouplingElements(mat0);
        ReadCouplingElements(matPlus);
        ReadCouplingElements(matMinus);
        matTOT.resize(mat0.rows()+2*bksize,mat0.cols()+2*bksize);
        matTOT.setZero();
        matTOT.block(bksize,bksize,mat0.rows(),mat0.cols())=mat0-Eigen::MatrixXcd::Identity(mat0.rows(),mat0.cols())*(Eomega);
        for(int i=0;i<scsize;++i)
        {
            matTOT.block(bksize*(i+1),bksize*i,bksize,bksize)=matMinus.block(bksize*(i),bksize*(i),bksize,bksize);//Eigen::MatrixXcd::Identity(bksize,bksize);
            matTOT.block(bksize*(i+1),bksize*(i+2),bksize,bksize)=matPlus.block(bksize*(i),bksize*(i),bksize,bksize);//-Eigen::MatrixXcd::Identity(bksize,bksize);
        }
        emat=matTOT;
        //Eigen::MatrixXcd tmp=matTOT.block(bksize,bksize,mat0.rows(),mat0.cols());
        //print_matrix(mat0,"mat0");
        //print_matrix(matPlus,"matPlus");
        //print_matrix(matMinus,"matMinus");
        //print_matrix(matTOT,"matTOT");
        //print_matrix(tmp,"tmp");
        //print_matrix(tmp.adjoint()-tmp,"Ttmp");

        //print_matrix(matTOT,"emat");

    }

    void Hamiltonian::GenH(std::complex<double> Eomega,Eigen::MatrixXcd &emat,std::vector<Eigen::MatrixXcd> &vW,int bksize,Eigen::Vector3d kk,Eigen::Vector3d npl)
    {

        std::complex<double> ii(0.0,1.0);
        Eigen::MatrixXcd matTOT,mat0,matPlus,matMinus;
        //Eigen::Vector3d kk={0.,0.,0.};
        if(0){
            SetMatrix5(kk,npl,systemcell.ListofMatrixElements_Type1,matTOT,-999);
            matTOT-=Eigen::MatrixXcd::Identity(MatrixSize,MatrixSize)*Eomega;
        }else{
            //SetMatrix5(kk,npl,systemcell.ListofMatrixElements_Type3,mat0,0);
            //SetMatrix5(kk,npl,systemcell.ListofMatrixElements_Type3,matPlus,1);
            //SetMatrix5(kk,npl,systemcell.ListofMatrixElements_Type3,matMinus,-1);
            SetMatrixKW5(0,npl,systemcell.ListofMatrixElements_Type3,matMinus,-1);
            SetMatrixKW5(0,npl,systemcell.ListofMatrixElements_Type3,matPlus,1);
            SetMatrixKW5(0,npl,systemcell.ListofMatrixElements_Type3,mat0,0);

            //print_matrix(matPlus,"matPlus");
            //print_matrix(matMinus,"matMinus");

            ReadCouplingElements(mat0);
            ReadCouplingElements(matPlus);
            ReadCouplingElements(matMinus);
            vW.resize(0);
            int enlargefactor=2;
            for(auto &cell : systemcell.supercell)
            {
                vW.push_back(Eigen::MatrixXcd());
                vW.back().resize(bksize*2*enlargefactor,bksize*2*enlargefactor);
                vW.back().setZero();
                vW.back().block(0,bksize*enlargefactor,bksize*enlargefactor,bksize*enlargefactor)=matPlus.block(cell->blpo,cell->blpo,bksize*enlargefactor,bksize*enlargefactor);
                vW.back().block(bksize*enlargefactor,0,bksize*enlargefactor,bksize*enlargefactor)=matMinus.block(cell->blpo,cell->blpo,bksize*enlargefactor,bksize*enlargefactor);
            }

            //print_matrix(matPlus,"matPlus");
            //print_matrix(matMinus,"matMinus");
            //print_matrix(vW[0],"vW");


            for(auto &cell : systemcell.supercell)
            {
                matPlus.block(cell->blpo+bksize,(cell->blpo+2*bksize)%MatrixSize,cell->AtomO/2.,cell->AtomO/2.)=matPlus.block(cell->blpo+bksize,cell->blpo,cell->AtomO/2.,cell->AtomO/2.);
                matPlus.block(cell->blpo+bksize,cell->blpo,cell->AtomO/2.,cell->AtomO/2.).setZero();
            }

            for(auto &cell : systemcell.supercell)
            {
                int pos=((cell->blpo-bksize)<0)?(MatrixSize+(cell->blpo-bksize)):(cell->blpo-bksize);
                matMinus.block(cell->blpo,pos,cell->AtomO/2.,cell->AtomO/2.)=matMinus.block(cell->blpo,cell->blpo+bksize,cell->AtomO/2.,cell->AtomO/2.);
                matMinus.block(cell->blpo,cell->blpo+bksize,cell->AtomO/2.,cell->AtomO/2.).setZero();
            }


            matTOT=(matMinus+matMinus.adjoint()+matPlus.adjoint()+matPlus)*0.5;

            matTOT=matTOT+mat0-Eigen::MatrixXcd::Identity(mat0.rows(),mat0.cols())*(Eomega);
        }


        emat.resize(matTOT.rows()+2*bksize,matTOT.cols()+2*bksize);
        emat.setZero();
        emat.block(bksize,bksize,matTOT.rows(),matTOT.cols())=matTOT;
        emat.block(bksize,0,bksize,bksize)=emat.block(bksize,matTOT.cols(),bksize,bksize);
        emat.block(bksize,matTOT.cols(),bksize,bksize).setZero();

        emat.block(matTOT.rows(),matTOT.cols()+bksize,bksize,bksize)=emat.block(matTOT.rows(),bksize,bksize,bksize);
        emat.block(matTOT.rows(),bksize,bksize,bksize).setZero();


        //print_matrix(emat,"emat");

    }
    void Hamiltonian::ArrangeT(Eigen::VectorXcd &values,Eigen::MatrixXcd &vectors)
    {
        Eigen::MatrixXcd tmpV;
        int dim=values.size()/2;
        tmpV.resize(dim*2,dim*2);
        for(int i=0;i<dim;++i)
            //for(int i=dim-1;i<=0;--i)
        {
            double Rek=std::log(values(i)).real();
            double Imk=std::log(values(i)).imag();
            tmpV.col(2*i)=vectors.col(i);
            tmpV.col(2*i+1)=vectors.col(dim*2-i-1);
            /*
               if(Imk>0){
               tmpV.col(2*i)=vectors.col(i);
               tmpV.col(2*i+1)=vectors.col(dim*2-i-1);
               }else{
               tmpV.col(2*i)=vectors.col(dim*2-i-1);
               tmpV.col(2*i+1)=vectors.col(i);
               }
               */
        }
        /*
           for(int i=0;i<dim;++i)
           {
           for(int j=0;j<dim-i-1;++j)
           {
           double nm1=abs(values(j));
           double nm2=abs(values(j+1));

           if(nm1>nm2)
           {
           Eigen::VectorXcd tmpVXcd=tmpV.col(j);
           tmpV.col(j)=tmpV.col(j+1);
           tmpV.col(j+1)=tmpVXcd;
           }
           }
           }

           for(int i=0;i<dim;++i)
           {
           for(int j=0;j<dim-i-1;++j)
           {
           double nm1=abs(values(j+dim));
           double nm2=abs(values(j+1+dim));

           if(nm1>nm2)
           {
           Eigen::VectorXcd tmpVXcd=tmpV.col(j+dim);
           tmpV.col(j+dim)=tmpV.col(j+1+dim);
           tmpV.col(j+1+dim)=tmpVXcd;
           }
           }
           }
           */
        vectors=tmpV;
    }

    void Hamiltonian::ArrangeDDD(Eigen::VectorXcd &values,Eigen::MatrixXcd &vectors)
    {
        std::complex<double> ii(0.0,1.0);
        Eigen::VectorXcd DLe=values;
        for(int i=0;i<DLe.size();++i)
        {
            DLe(i)=-ii*std::log(DLe(i));
        }
        int bksize=systemcell.orbNum;
        Eigen::MatrixXcd tmpII;
        tmpII.resize(2*bksize,1);
        tmpII.setZero();
        int flag_1=0;

        for(int i=0;i<values.size();++i)
        {
            double Re=DLe(i).real();
            double Im=DLe(i).imag();
            double val=std::abs(values(i));
            if     ((Re>=0) && fabs(Im-0)<1e-7 )  {tmpII(i,0)=2;flag_1+=1;}
            else if((Re<0) && fabs(Im-0)<1e-7 )   {tmpII(i,0)=-1;}
            else if((val-1)<1e-5)     {tmpII(i,0)=-2;}
            else if((val-1)>1e-5)     {tmpII(i,0)=1;}
        }

        for(int i=0;i<values.size();++i)
        {
            for(int j=0;j<values.size()-i-1;++j)
            {

                //if(tmpII(j,0).real()>tmpII(j+1,0).real())
                //if(DLe(j).real()>DLe(j+1).real())
                if(DLe(j).imag()>DLe(j+1).imag())
                {
                    /*
                       std::complex<double> tmpI=tmpII(j,0);
                       tmpII(j,0)=tmpII(j+1,0);
                       tmpII(j+1,0)=tmpI;
                       */
                    std::complex<double> tmpI=DLe(j);
                    DLe(j)=DLe(j+1);
                    DLe(j+1)=tmpI;

                    Eigen::VectorXcd tmpMXcd=vectors.col(j);
                    vectors.col(j)=vectors.col(j+1);
                    vectors.col(j+1)=tmpMXcd;
                }
            }
        }
    }

    int Hamiltonian::ArrangeDD(Eigen::VectorXcd &values,Eigen::MatrixXcd &vectors,Eigen::MatrixXcd &II,std::vector<Eigen::MatrixXcd> &vW)
    {

        std::complex<double> ii(0.0,1.0);
        Eigen::VectorXcd vL=-(vectors.adjoint()*vW.front()*vectors).diagonal();
        Eigen::VectorXcd DLe=values;
        for(int i=0;i<DLe.size();++i)
        {
            double Re=(std::log(DLe(i))/ii).real();
            double Im=(std::log(DLe(i))/ii).imag();
            DLe(i)=std::complex<double>(Re,Im);
        }
        int bksize=systemcell.orbNum;
        Eigen::MatrixXcd tmpII;
        tmpII.resize(2*bksize,1);
        tmpII.setZero();
        int flag_1=0;
        for(int i=0;i<values.size();++i)
        {

            double Rev=values(i).real();
            double Imv=values(i).imag();
            double vI=vL(i).imag();
            double Re=DLe(i).real();
            double Im=DLe(i).imag();
            double val=std::abs(values(i));
            if     ((Re>=0) && fabs(Im-0)<1e-5 )  {tmpII(i,0)=2;flag_1+=1;}
            else if((Re<0) && fabs(Im-0)<1e-5 )   {tmpII(i,0)=-1;}
            else if(Im<-1e-7)     {tmpII(i,0)=-2;}
            else if(Im>1e-7)     {tmpII(i,0)=1;flag_1+=1;}

            /* 
               double Rev=values(i).real();
               double Imv=values(i).imag();
               if(std::fabs(Rev-0)<1e-10 || std::fabs(Imv-0)<1e-10)
               {
               if(std::fabs(Rev-0)<1e-10)
               {
               if(Imv>0) {tmpII(i,0)=1;}
               else if(Imv<=0) {tmpII(i,0)=-1;}
               }else if(std::fabs(Imv-0)<1e-10){
               if(Rev>=1) {tmpII(i,0)=1;flag_1+=1;}
               else if(Rev<1) {tmpII(i,0)=-1;}
               }
               }else{
               double Re=(std::log(values(i))/ii).real();
               double Im=(std::log(values(i))/ii).imag();
               if     (Im==0 && Re<0)   {tmpII(i,0)=-1;}
               else if(Im==0 && Re>0)   {tmpII(i,0)=1;flag_1+=1;}
               else if(Im==0 && Re==0)  {tmpII(i,0)=1;flag_1+=1;}
               else if(Im>0 && Re==0)   {tmpII(i,0)=1;}
               else if(Im>0 && Re<0)    {tmpII(i,0)=-1;}
               else if(Im>0 && Re>0)    {tmpII(i,0)=1;}
               else if(Im<0 && Re==0)   {tmpII(i,0)=1;}
               else if(Im<0 && Re>0)    {tmpII(i,0)=1;}
               else if(Im<0 && Re<0)    {tmpII(i,0)=-1;}
               }
               */
        }

        for(int i=0;i<values.size();++i)
        {
            for(int j=0;j<values.size()-i-1;++j)
            {

                if(tmpII(j,0).real()<tmpII(j+1,0).real())
                {
                    std::complex<double> tmpI=tmpII(j,0);
                    tmpII(j,0)=tmpII(j+1,0);
                    tmpII(j+1,0)=tmpI;

                    Eigen::VectorXcd tmpMXcd=vectors.col(j);
                    vectors.col(j)=vectors.col(j+1);
                    vectors.col(j+1)=tmpMXcd;
                }
            }
        }

        //std::cout<<tmpII<<std::endl;
        //Dreverse(values,vectors);
        return flag_1;
    }

    void Hamiltonian::ArrangeDIM(Eigen::VectorXcd &values,Eigen::MatrixXcd &vectors,Eigen::MatrixXcd &II)
    {

        std::complex<double> ii(0.0,1.0);
        int bksize=systemcell.orbNum;
        Eigen::MatrixXcd tmpII;
        tmpII.resize(2*bksize,1);
        int flag_1=0;
        for(int i=0;i<values.size();++i)
        {

            double Rev=values(i).real();
            double Imv=values(i).imag();
            if(std::fabs(Rev-0)<1e-10 || std::fabs(Imv-0)<1e-10)
            {
                if(std::fabs(Rev-0)<1e-10)
                {
                    if(Imv>0) tmpII(i,0)=1;
                    else if(Imv<=0) tmpII(i,0)=-1;
                }else if(std::fabs(Imv-0)<1e-10){
                    if(Rev>=1) tmpII(i,0)=1;
                    else if(Rev<1) tmpII(i,0)=-1;
                }
            }else{
                double Re=(std::log(values(i))/ii).real();
                double Im=(std::log(values(i))/ii).imag();
                if     (Im==0 && Re<0)   {tmpII(i,0)=-1;}
                else if(Im==0 && Re>0)   {tmpII(i,0)=-1;}
                else if(Im==0 && Re==0)  {tmpII(i,0)=-1;}
                else if(Im>0 && Re==0)   {tmpII(i,0)=1;}
                else if(Im>0 && Re<0)    {tmpII(i,0)=1;}
                else if(Im>0 && Re>0)    {tmpII(i,0)=1;}
                else if(Im<0 && Re==0)   {tmpII(i,0)=-1;}
                else if(Im<0 && Re>0)    {tmpII(i,0)=-1;}
                else if(Im<0 && Re<0)    {tmpII(i,0)=-1;}
            }
        }

        for(int i=0;i<values.size();++i)
        {
            for(int j=0;j<values.size()-i-1;++j)
            {

                if(tmpII(j,0).real()>tmpII(j+1,0).real())
                {
                    std::complex<double> tmpI=tmpII(j,0);
                    tmpII(j,0)=tmpII(j+1,0);
                    tmpII(j+1,0)=tmpI;

                    Eigen::VectorXcd tmpMXcd=vectors.col(j);
                    vectors.col(j)=vectors.col(j+1);
                    vectors.col(j+1)=tmpMXcd;
                }
            }
        }

    }

    int Hamiltonian::ArrangeDRE(Eigen::VectorXcd &values,Eigen::MatrixXcd &vectors)
    {
        std::complex<double> ii(0.0,1.0);
        int M=values.size()/2;

        for(int i=0;i<M;++i)
        {
            for(int j=0;j<M-i-1;++j)
            {

                std::complex<double> Rek1=-ii*std::log(values(j));
                std::complex<double> Rek2=-ii*std::log(values(j+1));
                if(Rek1.imag()>=Rek2.imag())
                {
                    std::complex<double> tmpVXcd=values(j);
                    values(j)=values(j+1);
                    values(j+1)=tmpVXcd;
                    Eigen::VectorXcd tmpMXcd=vectors.col(j);
                    vectors.col(j)=vectors.col(j+1);
                    vectors.col(j+1)=tmpMXcd;
                }
            }
        }

        int flag=0;
        return flag;
    }

    void Hamiltonian::Dreverse(Eigen::VectorXcd &values,Eigen::MatrixXcd &vectors)
    {

        int dim=values.size();
        Eigen::MatrixXcd tmpV;
        tmpV.resize(dim,dim);
        tmpV.setZero();
        for(int i=0;i<dim;++i)
        {
            tmpV.col(i)=vectors.col(dim-i-1);
        }
        vectors=tmpV;
    }
    void Hamiltonian::Dswap(Eigen::VectorXcd &values,Eigen::MatrixXcd &vectors)
    {
        int valSize=values.size();
        Eigen::MatrixXcd tmp=vectors.block(0,valSize*0.5,valSize,valSize*0.5);
        vectors.block(0,valSize*0.5,valSize,valSize*0.5)=vectors.block(0,0,valSize,valSize*0.5);
        vectors.block(0,0,valSize,valSize*0.5)=tmp;

    }
    void Hamiltonian::print_matrix_in_details(Eigen::MatrixXcd mat,std::string name)
    {
        std::complex<double> trace=mat.trace();
        //fprintf(stderr, "\n %s %lutrace[%1.2f,%1.2f]\n",name.c_str(),mat.cols(),trace.real(),trace.imag());
        fprintf(stderr, "\n %s %ld %ld\n",name.c_str(),mat.rows(),mat.cols());

        for(int i=0;i<mat.rows();i++)
        {
            for(int j=0;j<mat.cols();j++)
            {
                if(fabs(mat(i,j).imag()-0)<1e-5 && fabs(mat(i,j).real()-0)<1e-5)
                {
                    //fprintf(stderr,"{   }");
                    //fprintf(stderr,"{      }");
                    fprintf(stderr,"{            }");
                    //fprintf(stderr,"{  }");
                    //fprintf(stderr,"{ }");
                }else{
                    //fprintf(stderr,"{%2.0f}",mat(i,j).real());
                    //fprintf(stderr,"{%2.3f}",mat(i,j).real());
                    //fprintf(stderr,"{%2.3f/%2.3f}",fabs(mat(i,j).real()),fabs(mat(i,j).imag()));
                    fprintf(stderr,"{%2.3f/%2.3f}",(mat(i,j).real()),(mat(i,j).imag()));
                    //fprintf(stderr,"{%1.0f}",abs(mat(i,j).real()));
                    //fprintf(stderr,"{%1.0f}",mat(i,j).real());
                }
            }fprintf(stderr,"\n");
        }

    }
    void Hamiltonian::SetMatrixK5(std::complex<double> kz,Eigen::Vector3d pl,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::MatrixXcd &Kmatrix5,int Wsplit)
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

            std::complex<double> efactor=pow(2*M_PI*ii*dd.dot(pl),Wsplit)*exp(2*M_PI*(dd.dot(pl)*kz)*ii);
            //std::complex<double> efactor=1;
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
                        if(rows==cols)
                        {

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
    void Hamiltonian::eGenVT(Eigen::MatrixXcd &matTOT,std::vector<Eigen::MatrixXcd> &vT,std::vector<Eigen::MatrixXcd> &vT2,std::vector<Eigen::MatrixXcd> &vJ,std::vector<Eigen::MatrixXcd> &vG,int bksize,int scsize)
    {
        std::complex<double> ii(0.0,1.0);
        vT.resize(0);
        vJ.resize(0);
        for(int i=0;i<scsize;++i)
        {
            Eigen::MatrixXcd vTi,vJi;
            vTi.resize(bksize*2,bksize*2);
            vTi.setZero();
            vJi.resize(bksize*2,bksize*2);
            vJi.setZero();
            Eigen::MatrixXcd minus=matTOT.block(bksize*(i+1),bksize*(i),bksize,bksize);
            Eigen::MatrixXcd middle=matTOT.block(bksize*(i+1),bksize*(i+1),bksize,bksize);
            Eigen::MatrixXcd plus=matTOT.block(bksize*(i+1),bksize*(i+2),bksize,bksize);

            vTi.block(0,0,bksize,bksize)=-minus.inverse()*middle;
            vTi.block(0,bksize,bksize,bksize)=-minus.inverse()*plus;
            vTi.block(bksize,0,bksize,bksize)=Eigen::MatrixXcd::Identity(bksize,bksize);
            
            vT.push_back(vTi);

            vJi.block(0,bksize,bksize,bksize)=plus;
            vJi.block(bksize,0,bksize,bksize)=minus;

            vJ.push_back(vJi);

            vTi.block(0,0,bksize,bksize)=-plus.inverse()*minus;
            vTi.block(0,bksize,bksize,bksize)=-plus.inverse()*middle;
            vTi.block(bksize,0,bksize,bksize)=Eigen::MatrixXcd::Identity(bksize,bksize);
            vT2.push_back(vTi);

            //matTOT.block(bksize*(i+1),bksize*(i),bksize,bksize)=plus.adjoint();
            //matTOT.block(bksize*(i+1),bksize*(i+2),bksize,bksize)=minus.adjoint();

/*
            if(i>0){
                Eigen::MatrixXcd T=vT[i-1]*vT[i];
                csolver.compute(T,Eigen::DecompositionOptions::ComputeEigenvectors);
                Eigen::VectorXcd Tval=csolver.eigenvalues();
                Eigen::MatrixXcd Tvec=csolver.eigenvectors();
                //ArrangeDRE(Tval,Tvec);
                //Dreverse(Tval,Tvec);
                Eigen::MatrixXcd S1=Tvec.block(0,bksize,bksize,bksize);
                Eigen::MatrixXcd S2=Tvec.block(bksize,bksize,bksize,bksize);
                //T=(vT[i]*vT[i-1]).inverse();
                //csolver.compute(T,Eigen::DecompositionOptions::ComputeEigenvectors);
                //Tval=csolver.eigenvalues();
                //Tvec=csolver.eigenvectors();
                //ArrangeDRE(Tval,Tvec);
                //Dreverse(Tval,Tvec);
                Eigen::MatrixXcd S3=Tvec.block(0,0,bksize,bksize);
                Eigen::MatrixXcd S4=Tvec.block(bksize,0,bksize,bksize);
                //Eigen::MatrixXcd G0=-(middle+minus*S2*S1.inverse()+minus.adjoint()*S3*S4.inverse()).inverse();
                //Eigen::MatrixXcd G0=-(middle+plus*S2*S1.inverse()).inverse();
                Eigen::MatrixXcd G0=-(middle+plus*S1*S2.inverse()).inverse();
                vG.push_back(G0);
            }
*/
            //print_matrix(G0-(-middle-minus*G0*plus).inverse(),"h00");

            //vTi.block(0,0,bksize,bksize)=-plus.inverse()*middle;
            //vTi.block(0,bksize,bksize,bksize)=-plus.inverse()*minus;
            //vTi.block(bksize,0,bksize,bksize)=Eigen::MatrixXcd::Identity(bksize,bksize);

        }

/*
        for(int i=0;i<vT.size()-1;++i)
        {

            vW.push_back(matTOT.block(i*bksize+bksize,i*bksize+bksize,bksize*2,bksize*2));
            //vW.push_back(matTOT.block(i*bksize,i*bksize,bksize*2,bksize*2));

            vW[i].block(0,0,bksize,bksize).setZero();
            vW[i].block(bksize,bksize,bksize,bksize).setZero();

            //vW[i].block(bksize,0,bksize,bksize)=-vW[i].block(bksize,0,bksize,bksize);
            //vW[i]=ii*vW[i];

        }
*/
        /*
           for(int i=0;i<vT.size()-2;++i)
           {
           vW.push_back(matTOT.block(i*bksize+bksize+bksize,i*bksize+bksize+bksize,bksize*2,bksize*2));
        //vW.push_back(matTOT.block(i*bksize,i*bksize,bksize*2,bksize*2));
        vW[i].block(0,0,bksize,bksize).setZero();
        vW[i].block(bksize,bksize,bksize,bksize).setZero();
        vW[i].block(bksize,0,bksize,bksize)=-vW[i].block(bksize,0,bksize,bksize);
        }
        */
        //print_matrix_in_details(vW.front(),"test");

    }

    void Hamiltonian::Test_vW(Eigen::MatrixXcd matII,std::vector<Eigen::MatrixXcd> &vW,std::vector<Eigen::MatrixXcd> &vT,Eigen::MatrixXcd &J0,Eigen::MatrixXcd &JN,int bksize)
    {
        int col=matII.cols();
        J0.resize(0,0);
        JN.resize(0,0);
        if(col!=0){
            for(int i=0;i<vW.size();++i)
            {
                Eigen::MatrixXcd CC,tmpvW;
                CC=matII.block(i*(bksize)+bksize,0,2*bksize,col);
                tmpvW=CC.adjoint()*vW[i]*CC;
                //print_matrix_in_details(tmpvW,"vW");
                //vW[i]=tmpvW.diagonal().normalized().asDiagonal();
                vW[i]=tmpvW.diagonal().asDiagonal();
                //vW[i]*=-1.;
            }
            J0=vW.front();
            JN=vW.back();
        }
    }

    void Hamiltonian::NormTR(Eigen::MatrixXcd &tt,Eigen::MatrixXcd &rr,int bksize)
    {
        Eigen::MatrixXcd II;
        II.resize(rr.rows()+tt.rows(),rr.cols());
        II.block(0,0,rr.rows(),rr.cols())=rr;
        II.block(rr.rows(),0,tt.rows(),tt.cols())=tt;
        for(int i=0;i<II.cols();++i)
        {
            II.col(i)=II.col(i)/sqrt(II.col(i).squaredNorm());
        }
        rr=II.block(0,0,rr.rows(),rr.cols());
        tt=II.block(rr.rows(),0,tt.rows(),tt.cols());
    }
    void Hamiltonian::ComplexBands(std::complex<double> Eomega,Eigen::Vector3d k)
    {

        std::complex<double> ii(0.0,1.0);
        int bksize=systemcell.orbNum;
        int scsize=systemcell.supercell.size()*2.0;
        Eigen::MatrixXcd matTOT,minus,middle,plus;
        std::vector<Eigen::MatrixXcd> vT,vT2,vW,vG;

        Eigen::Vector3d kk={0.,0.,0.};
        Eigen::Vector3d pl={1,0,0};
        Eigen::Vector3d npl=pl.normalized();

        Eigen::MatrixXcd kH;
        Eigen::MatrixXcd vTi,vTj;
        Eigen::MatrixXcd kTi,kTj;

        if(pl==Eigen::Vector3d(1,1,0))
        {
            GenH2(Eomega,matTOT,bksize*2,scsize/2,kk,npl);
            eGenVT(matTOT,vT,vT2,vW,vG,bksize*2,scsize/2);
            kH=vT.front(); 
        }else{
            GenH(Eomega,matTOT,vW,bksize,kk,npl);
            eGenVT(matTOT,vT,vT2,vW,vG,bksize,scsize);
            kH=vT.front()*vT.back();
        }

        csolver.compute(kH,Eigen::DecompositionOptions::ComputeEigenvectors);
        Eigen::VectorXcd kkv=csolver.eigenvalues();
        Eigen::MatrixXcd kval=kkv.asDiagonal();
        for(int i=0;i<kval.rows();++i)
        {
            kval(i,i)=(-ii*std::log(kval(i,i)))*((pl*2.).norm()/M_PI);
        }


        std::fstream file_obj; 
        file_obj.open("./Cdata.txt",std::ios::out|std::ios::app);
        file_obj<<Eomega.real();
        for(int i=0;i<kval.rows();++i)
        {
            double Rek=kval(i,i).real();
            double Imk=kval(i,i).imag();
            double newImk=-std::log(csolver.eigenvalues()(i)).real();
            std::complex<double> kfit=std::exp(ii*(Rek+ii*Imk)*pl.norm()*2.);
            //std::cout<<std::fixed<<std::setprecision(5)<<abs(kfit-csolver.eigenvalues()(i))<<std::endl;
            //if(Imk>0 && abs(abs(csolver.eigenvalues()(i))-1)<1e-5) file_obj<<"\t"<<Imk; else file_obj<<"\t"<<0;
            //if(Rek<0 && abs(abs(csolver.eigenvalues()(i))-1)<1e-5) file_obj<<"\t"<<Rek; else file_obj<<"\t"<<0;

            if(Imk>0) file_obj<<"\t"<<Imk; else file_obj<<"\t"<<0;
            if(Rek<0) file_obj<<"\t"<<Rek; else file_obj<<"\t"<<0;

            //file_obj<<"\t"<<Imk;
            //if(abs(kabs(i,i).real()-1)<1e-5) file_obj<<"\t"<<Rek; else file_obj<<"\t"<<-sqrt(3)*0.5;
            //file_obj<<"\t"<<Rek;

        }
        file_obj<<std::endl;
        file_obj.close();

    }

   
void Hamiltonian::CaseComplexBands5(std::complex<double> Eomega,Eigen::Vector3d k)
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
    std::vector<Eigen::MatrixXcd> vT,vT2,vW,vG;

    Eigen::Vector3d kk={0.,0.,0.};
    Eigen::Vector3d pl={0.25,0.25,0.25};
    //Eigen::Vector3d pl={0.25,0.25,0.};
    //Eigen::Vector3d pl={1,0.,0.};
    Eigen::Vector3d npl=pl.normalized();

    Eigen::MatrixXcd G0;
    GenH(Eomega,matTOT,vW,bksize,kk,npl);
    eGenVT(matTOT,vT,vT2,vW,vG,bksize,scsize);
    Eigen::MatrixXcd minus=matTOT.block(bksize*(1),bksize*(0),bksize,bksize).adjoint();
    Eigen::MatrixXcd middle=matTOT.block(bksize*(1),bksize*(1),bksize,bksize);
    Eigen::MatrixXcd plus=matTOT.block(bksize*(1),bksize*(2),bksize,bksize).adjoint();
    Eigen::MatrixXcd kH;
    Eigen::MatrixXcd vTi,vTj;
    Eigen::MatrixXcd kTi,kTj;
    vTi.resize(bksize*2,bksize*2);
    vTi.setZero();
    kTi.resize(bksize*2,bksize*2);
    kTi.setZero();
    kTi.block(0,0,bksize,bksize)=minus.inverse();
    vTi.block(0,0,bksize,bksize)=-minus.inverse()*middle;
    vTi.block(0,bksize,bksize,bksize)=-minus.inverse()*plus;
    vTi.block(bksize,0,bksize,bksize)=Eigen::MatrixXcd::Identity(bksize,bksize);

    minus=matTOT.block(bksize*(2),bksize*(1),bksize,bksize).adjoint();
    middle=matTOT.block(bksize*(2),bksize*(2),bksize,bksize);
    plus=matTOT.block(bksize*(2),bksize*(3),bksize,bksize).adjoint();
    vTj.resize(bksize*2,bksize*2);
    vTj.setZero();
    kTj.resize(bksize*2,bksize*2);
    kTj.setZero();
    kTj.block(0,0,bksize,bksize)=minus.inverse();
    vTj.block(0,0,bksize,bksize)=-minus.inverse()*middle;
    vTj.block(0,bksize,bksize,bksize)=-minus.inverse()*plus;
    vTj.block(bksize,0,bksize,bksize)=Eigen::MatrixXcd::Identity(bksize,bksize);

    kH=vTi*vTj;
    Eigen::MatrixXcd kT=kTi*vTj+vTi*kTj;
    //Eigen::MatrixXcd kH=vT[0]*vT[1];

    csolver.compute(kH,Eigen::DecompositionOptions::ComputeEigenvectors);
    Eigen::VectorXcd kkv=csolver.eigenvalues();
    Eigen::MatrixXcd kval=kkv.asDiagonal();

    for(int i=0;i<kval.rows();++i)
    {
        kval(i,i)=(-ii*std::log(kval(i,i)))*((pl*2.).norm()/M_PI);
    }
    Eigen::MatrixXcd kabs=kkv.rowwise().norm().asDiagonal();
    Eigen::MatrixXcd kHvect=csolver.eigenvectors();
    //print_matrix_in_details(kHvect.inverse()*kH*kHvect,"before");
    //MatrixO(kHvect);
    //print_matrix_in_details(kHvect.inverse()*kH*kHvect,"after");
    //ArrangeDRE(kkv,kHvect);
    Eigen::MatrixXcd Vk2=(kHvect.inverse()*kT*kHvect).diagonal().asDiagonal();
    for(int i=0;i<Vk2.rows();++i)
    {
        Vk2(i,i)=Vk2(i,i)/(ii*kkv(i)/((pl*2.).norm()/M_PI));

        if(
                (Vk2(i,i).real()<0 && abs(Vk2(i,i).real())>1e-5 && abs(Vk2(i,i).imag())<1e-5)
          )
        {
            Vk2(i,i)=-Vk2(i,i);//.real()+Vk2(i,i).imag()*ii;
        }else if(abs(Vk2(i,i).imag())>1e-5 && Vk2(i,i).imag()<0 /*&& abs(Vk2(i,i).real())<1e-5*/){
            //std::cout<<std::fixed<<std::setprecision(5)<<"V="<<Vk2(i,i)<<std::endl;
            Vk2(i,i)=-Vk2(i,i);//.real()+Vk2(i,i).imag()*ii;
        }
    }
    int Vsize=Vk2.rows();
    Eigen::MatrixXcd CC=kHvect;
    Eigen::MatrixXcd S2=CC.block(0,0,5,5);
    Eigen::MatrixXcd S1=CC.block(5,0,5,5);
    Eigen::MatrixXcd S4=CC.block(0,5,5,5);
    Eigen::MatrixXcd S3=CC.block(5,5,5,5);

    minus=matTOT.block(bksize*(1),bksize*(0),bksize,bksize);
    middle=matTOT.block(bksize*(1),bksize*(1),bksize,bksize);
    plus=matTOT.block(bksize*(1),bksize*(2),bksize,bksize);
    G0=-M_PI*ii*Vk2;
    //print_matrix_in_details(G0,"");
    /*
       double results=0;
       std::complex<double> Cres=1;
       for(int i=0;i<Vk2.rows();++i)
       {
       results+=abs(Vk2(i,i).real());
       Cres*=Vk2(i,i);
       }
       */
    //G0=-(middle+minus*S1*S2.inverse()/*+plus*S3*S4.inverse()*/).inverse();

    std::fstream file_obj; 
    file_obj.open("./Cdata.txt",std::ios::out|std::ios::app);
    file_obj<<Eomega.real();
    file_obj<<"\t"<<-M_PI*abs(G0.trace().imag());
    //file_obj<<"\t"<<abs((G0).determinant().real());

    //file_obj<<"\t"<<abs(-M_PI*(vG[0].trace()).imag());
    //file_obj<<"\t"<<abs(vG[0].determinant().real());
    for(int i=0;i<kval.rows();++i)
    {
        double Rek=kval(i,i).real();
        double Imk=kval(i,i).imag();
        double newImk=-std::log(csolver.eigenvalues()(i)).real();
        std::complex<double> kfit=std::exp(ii*(Rek+ii*Imk)*pl.norm()*2.);
        //std::cout<<std::fixed<<std::setprecision(5)<<abs(kfit-csolver.eigenvalues()(i))<<std::endl;
        //if(Imk>0 && abs(abs(csolver.eigenvalues()(i))-1)<1e-5) file_obj<<"\t"<<Imk; else file_obj<<"\t"<<0;
        //if(Rek<0 && abs(abs(csolver.eigenvalues()(i))-1)<1e-5) file_obj<<"\t"<<Rek; else file_obj<<"\t"<<0;

        //if(Imk>0) file_obj<<"\t"<<Imk; else file_obj<<"\t"<<0;
        //if(Rek<0) file_obj<<"\t"<<Rek; else file_obj<<"\t"<<0;

        //file_obj<<"\t"<<Imk;
        //if(abs(kabs(i,i).real()-1)<1e-5) file_obj<<"\t"<<Rek; else file_obj<<"\t"<<-sqrt(3)*0.5;
        //file_obj<<"\t"<<Rek;

    }
    file_obj<<std::endl;
    file_obj.close();
}


    void Hamiltonian::Case5(std::complex<double> Eomega)
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
        std::vector<Eigen::MatrixXcd> vT,vT2,vW,vG,vJ;
        Eigen::Vector3d kk={0.,0.,0.};
        //Eigen::Vector3d pl={0.25,0.25,0.25};
        //pl=pl*sqrt(2)*8.;
        


        //Eigen::Vector3d pl={1,1,1};
        Eigen::Vector3d pl={1,0,0};
        Eigen::Vector3d npl=pl.normalized();

        Eigen::MatrixXcd D1vk;
        Eigen::MatrixXcd D2vk;
        if(pl==Eigen::Vector3d(1,1,0))
        {
            GenH2(Eomega,matTOT,bksize*2,scsize/2,kk,npl);
            eGenVT(matTOT,vT,vT2,vJ,vG,bksize*2,scsize/2);
            Method1(npl,matTOT,tt1,rr1,matII,vT,Tac,Tca,Eomega,bksize*2,scsize/2,vW,D1vk,D2vk);
            //Test_vW(matII,vW,vT,J0,JN,bksize*2);//TTT
        }else{
            GenH(Eomega,matTOT,vW,bksize,kk,npl);
            eGenVT(matTOT,vT,vT2,vJ,vG,bksize,scsize);
            Method1(npl,matTOT,tt1,rr1,matII,vT,Tac,Tca,Eomega,bksize,scsize,vW,D1vk,D2vk);
            //Test_vW(matII,vW,vT,J0,JN,bksize);//TTT
        }

        tt=tt1;
        rr=rr1;
        //NormTR(tt,rr,bksize);
        std::complex<double> ww={2,2};
        std::fstream file_obj; 
        file_obj.open("./Cdata.txt",std::ios::out|std::ios::app);
        file_obj<<Eomega.real();
        double Rrets=0,Trets=0;
        //print_matrix(D1vk,"D1vk");
        for(int j=0;j<tt.cols();++j)
        {
            for(int i=0;i<tt.rows();++i)
            {
                //if(abs(D2vk(i,j).imag())<1e-1 && abs(D1vk(j,j).imag())<1e-1)
                {
                    Trets+=(std::conj(tt(i,j))*tt(i,j)*abs(D2vk(i,j))/abs(D1vk(j,j))).real();
                }
            }
        }
        for(int j=0;j<rr.cols();++j)
        {
            for(int i=0;i<rr.rows();++i)
            {
                //if(abs(D1vk(rr.rows()/2+i,rr.rows()/2+j).imag())<1e-1 && abs(D1vk(j,j).imag())<1e-1)
                {
                    Rrets+=(std::conj(rr(i,j))*rr(i,j)*abs(D1vk(rr.rows()/2+i,rr.rows()/2+j))/abs(D1vk(j,j))).real();
                }
            }
        }

        Eigen::MatrixXcd G0=-M_PI*ii*D1vk.block(0,0,D1vk.rows()/2,D1vk.cols()/2);
        //file_obj<<"\t"<<abs((G0.trace()).imag());
        file_obj<<"\t"<<Rrets;
        file_obj<<"\t"<<Trets;
        //file_obj<<"\t"<<abs(Trets)+abs(Rrets);
        file_obj<<std::endl;
        file_obj.close();
    }
/*
    void Hamiltonian::SmallMethod1(Eigen::MatrixXcd matTOT,Eigen::MatrixXcd &tt,Eigen::MatrixXcd &rr,Eigen::MatrixXcd &matII,std::vector<Eigen::MatrixXcd> vT,Eigen::MatrixXcd T1,Eigen::MatrixXcd T2,std::complex<double> Eomega,int bksize,int scsize,std::vector<Eigen::MatrixXcd> &vW)
    {

        Eigen::MatrixXcd D1=vT.front()*vT.back();
        csolver.compute(D1,Eigen::DecompositionOptions::ComputeEigenvectors);
        Eigen::VectorXcd tmpD1=csolver.eigenvalues();
        Eigen::MatrixXcd DL=csolver.eigenvectors();
        Eigen::MatrixXcd eII;

        int flag=ArrangeDD(tmpD1,DL,eII,vW);

        int small_dim=matTOT.rows()-2.*bksize;
        eII.resize(small_dim,flag);
        eII.setZero();
        eII.block(0,0,flag,flag)=Eigen::MatrixXcd::Identity(flag,flag);
        //int flag=ArrangeDRE(tmpD1,DL);
        DL=DL.inverse();//DR=DR.inverse();
        Eigen::MatrixXcd DR=DL;

        matTOT.block(0+bksize,0+bksize,bksize,bksize)=DL.topLeftCorner(bksize,bksize);
        matTOT.block(0+bksize,bksize+bksize,bksize,bksize)=DL.topRightCorner(bksize,bksize);
        matTOT.block(matTOT.rows()-bksize-bksize,small_dim-bksize,bksize,bksize)=DR.bottomLeftCorner(bksize,bksize);
        matTOT.block(matTOT.rows()-bksize-bksize,matTOT.cols()-bksize-bksize,bksize,bksize)=DR.bottomRightCorner(bksize,bksize);
        //print_matrix(matTOT,"mattot");

        if(eII.cols()!=0){
            matII=matTOT.block(bksize,bksize,small_dim,small_dim).inverse()*eII;
            tt=DR.block(0,0,bksize,2.*bksize)*matII.block(small_dim-2*bksize,0,2.*bksize,matII.cols());
            rr=DL.block(bksize,0,bksize,2.*bksize)*matII.block(0,0,2.*bksize,matII.cols());
        }

    }
*/
    void Hamiltonian::Method1(Eigen::Vector3d npl,Eigen::MatrixXcd matTOT,Eigen::MatrixXcd &tt,Eigen::MatrixXcd &rr,Eigen::MatrixXcd &matII,std::vector<Eigen::MatrixXcd> vT,Eigen::MatrixXcd T1,Eigen::MatrixXcd T2,std::complex<double> Eomega,int bksize,int scsize,std::vector<Eigen::MatrixXcd> &vW,Eigen::MatrixXcd &D1vk,Eigen::MatrixXcd &D2vk)
    {

        std::complex<double> ii(0.0,1.0);
        Eigen::MatrixXcd eII,eII2,b;
        //Eigen::MatrixXcd DD=vT.back()*vT.front();
        //Eigen::MatrixXcd DD=(vT.front()*vT.back()).inverse();
        Eigen::MatrixXcd DD=(vT.front()*vT.back());

        csolver.compute(DD,Eigen::DecompositionOptions::ComputeEigenvectors);
        Eigen::VectorXcd DDvalue=csolver.eigenvalues();
        Eigen::MatrixXcd DDvector=csolver.eigenvectors();
        Eigen::MatrixXcd DDvk;

        Eigen::MatrixXcd tmpvW=vW[0].topRightCorner(vW[0].rows()/2,vW[0].cols()/2);
        DDvk=(DDvector.inverse()*(tmpvW-tmpvW.adjoint())*DDvector)*ii;
        //GetVK(DDvk,DDvector,DDvalue,npl,Eomega);


        Eigen::MatrixXcd D1=(vT.front()*vT.back()).inverse();
        Eigen::MatrixXcd D2=(vT.front()*vT.back()).inverse();

        csolver.compute(D1,Eigen::DecompositionOptions::ComputeEigenvectors);
        Eigen::VectorXcd D1value=csolver.eigenvalues();
        Eigen::MatrixXcd D1vector=csolver.eigenvectors();
        csolver.compute(D2,Eigen::DecompositionOptions::ComputeEigenvectors);
        Eigen::VectorXcd D2value=csolver.eigenvalues();
        Eigen::MatrixXcd D2vector=csolver.eigenvectors();

        //GetVK(D1vk,D1vector,D1value,npl,Eomega);
        //GetVK(D2vk,D2vector,D2value,npl,Eomega);
        D1vk=DDvk;
        D2vk=DDvk;
/*
        Eigen::MatrixXcd M=D1vector;
        for(int i=0;i<vT.size();++i)
        {
            M=vT[i].inverse()*M;
        }
        M=D2vector.inverse()*M;
        //Eigen::MatrixXcd MM=D1vector.inverse()*D2vector;//M.inverse();
        Eigen::MatrixXcd MM=M.inverse();
*/ 

        Eigen::MatrixXcd MM=D2vector;
        for(int i=vT.size()-1;i>=0;--i){
            MM=vT[i]*MM;
        }
        MM=D1vector.inverse()*MM;
        Eigen::MatrixXcd S0;
        MtoS(MM,S0,0);



        Eigen::MatrixXcd theS,newS;
        Eigen::MatrixXcd M01=D2vector;
        MtoS(M01,theS,0);
        for(int i=vT.size()-1;i>=0;--i){
            M01=vT[i];
            MtoS(M01,newS,0);
            SxS(newS,theS);
        }
        M01=D1vector.inverse();
        MtoS(M01,newS,0);
        SxS(newS,theS);
        print_matrix(theS-S0,"theS");

        int colInj=D1vector.cols()/2;
        Eigen::MatrixXcd Inj=
            //Eigen::MatrixXcd::Zero(M.rows()/2,colInj);
            //Eigen::MatrixXcd::Ones(M.rows()/2,colInj);
            Eigen::MatrixXcd::Identity(D1vector.rows()/2,colInj);

        //Inj(3,3)=1;
        //Inj(4,4)=1;
/*
        Inj(0,4)=0;
        Inj(1,4)=0;
        Inj(2,4)=0;
        Inj(3,4)=0;
        Inj(4,4)=1;
*/

   
    
        tt=MM.block(0,0,MM.rows()/2,MM.cols()/2).inverse()*Inj;
        rr=MM.block(MM.rows()/2,0,MM.rows()/2,MM.cols()/2)*tt;
       

        Eigen::MatrixXcd t01=M01.block(0,0,MM.rows()/2,MM.cols()/2).inverse()*Inj;
        Eigen::MatrixXcd r01=M01.block(MM.rows()/2,0,MM.rows()/2,MM.cols()/2)*t01;
        /*
        Eigen::MatrixXcd rt01;
        rt01.resize(MM.rows(),MM.cols()/2);
        rt01.block(0,0,MM.cols()/2,MM.cols()/2)=r01;
        rt01.block(MM.cols()/2,0,MM.cols()/2,MM.cols()/2)=t01;
*/
        Eigen::MatrixXcd t0,Ir,rt,I0;


/*
        t0.resize(MM.rows(),MM.cols()/2);
        Ir.resize(MM.rows(),MM.cols()/2);
        Ir.block(0,0,MM.cols()/2,MM.cols()/2)=Inj;
        Ir.block(MM.cols()/2,0,MM.cols()/2,MM.cols()/2)=r01;
        t0.block(0,0,MM.cols()/2,MM.cols()/2)=t01;
        t0.block(MM.cols()/2,0,MM.cols()/2,MM.cols()/2).setZero();

        print_matrix(Ir-M01*t0,"I'r'-S01_t0");


        print_matrix(S0,"S2");


        Ir.resize(MM.rows(),MM.cols()/2);
        t0.resize(MM.rows(),MM.cols()/2);
        Ir.block(0,0,MM.cols()/2,MM.cols()/2)=Inj;
        Ir.block(MM.cols()/2,0,MM.cols()/2,MM.cols()/2)=rr;
        t0.block(0,0,MM.cols()/2,MM.cols()/2)=tt;
        t0.block(MM.cols()/2,0,MM.cols()/2,MM.cols()/2).setZero();

        print_matrix(Ir-M02*M01*t0,"Ir-St0");

 
        rt.resize(MM.rows(),MM.cols()/2);
        I0.resize(MM.rows(),MM.cols()/2);
        rt.block(0,0,MM.cols()/2,MM.cols()/2)=r01;
        rt.block(MM.cols()/2,0,MM.cols()/2,MM.cols()/2)=t01;
        I0.block(0,0,MM.cols()/2,MM.cols()/2)=Eigen::MatrixXcd::Identity(MM.cols()/2,MM.cols()/2);
        I0.block(MM.cols()/2,0,MM.cols()/2,MM.cols()/2).setZero();
        print_matrix(rt-S01*I0,"r01t01-S01_I0");
 
        rt.resize(MM.rows(),MM.cols()/2);
        I0.resize(MM.rows(),MM.cols()/2);
        rt.block(0,0,MM.cols()/2,MM.cols()/2)=rr;
        rt.block(MM.cols()/2,0,MM.cols()/2,MM.cols()/2)=tt;
        I0.block(0,0,MM.cols()/2,MM.cols()/2)=Eigen::MatrixXcd::Identity(MM.cols()/2,MM.cols()/2);
        I0.block(MM.cols()/2,0,MM.cols()/2,MM.cols()/2).setZero();
        print_matrix(rt-S0*I0,"rt-S0_I0");

        //print_matrix(matTOT,"matTOT");
 */       
/*
        for(int i=0;i<vJ.size();++i)
        {
            Eigen::MatrixXcd mxC=D1vector;
            for(int j=0;j<i;++j)
            {
                mxC=vT[j]*mxC;
            }
            Eigen::MatrixXcd mxC2=mxC;
            mxC2.block(mxC2.rows()/2,0,mxC2.rows()/2,mxC2.cols()/2)=-mxC2.block(mxC2.rows()/2,0,mxC2.rows()/2,mxC2.cols()/2);
            //print_matrix(vJ[i],"vJ");
            print_matrix((mxC.adjoint()*vJ[i]*mxC2).real(),"vJ");
        }
*/
        //print_matrix(M*Ir-tzero,"test");
       /* 
        csolver.compute(M,Eigen::DecompositionOptions::ComputeEigenvectors);
        Eigen::MatrixXcd Mvalue=csolver.eigenvalues().asDiagonal();
        std::cout<<Mvalue.diagonal()<<std::endl;
        std::cout<<Mvalue.diagonal().rowwise().norm()<<std::endl;
        Eigen::MatrixXcd Mvector=csolver.eigenvectors();
       
        Eigen::MatrixXcd TEST2=(Mvalue*Mvector.inverse())*(Ir)-Mvector.inverse()*tzero;

        print_matrix(Mvector.inverse()*Ir,"Ir");
        print_matrix(Mvector.inverse()*tzero,"tzero");
        print_matrix(TEST2,"test3");
*/
        //print_matrix(D2vk*Ir-D1vk*tzero,"test2");
        //print_matrix(tzero,"tzero");
        //Eigen::MatrixXcd TT=tt.adjoint()*tt*D2vk.block(0,0,D2vk.rows()/2,D2vk.cols()/2)*Icoef.inverse();
        //Eigen::MatrixXcd RR=rr.adjoint()*rr*D1vk.block(D1vk.rows()/2,D1vk.cols()/2,D1vk.rows()/2,D1vk.cols()/2)*Icoef.inverse();
        //tt=TT;
        //rr=TT;


        int flag=D1value.size()/2;



/*
        Eigen::MatrixXcd matTOT2=matTOT;
        matTOT2.setZero();

        for(int i=0;i<scsize/2;++i)
        {
            //matTOT2.block((i+1)*bksize*2,(i+1)*bksize*2,bksize*2,bksize*2)=-D1;
            matTOT2.block((i+1)*bksize*2,(i+1)*bksize*2,bksize*2,bksize*2)=-vT[i*2]*vT[i*2+1];
            matTOT2.block(i*bksize*2,(i+1)*bksize*2,bksize*2,bksize*2)=Eigen::MatrixXcd::Identity(bksize*2,bksize*2);
        }
        //matTOT2.block(0,bksize,bksize*2,bksize)=-D1*D.block(0,bksize,bksize*2,bksize);
        //matTOT2.block(bksize*2*scsize/2,0,bksize*2,bksize)=D.block(0,0,bksize*2,bksize);
        
        matTOT2.block(0,bksize,bksize*2,bksize)=-D1vector*D1.block(0,0,bksize*2,bksize);
        matTOT2.block(bksize*2*scsize/2,0,bksize*2,bksize)=D1.block(0,bksize,bksize*2,bksize);
        
        eII2.resize(matTOT.rows(),flag);
        eII2.setZero();

        //eII2.block(0,0,bksize*2,flag)=D1*D.block(0,bksize,bksize*2,flag);
        eII2.block(0,0,bksize*2,flag)=D1*D1vector.block(0,0,bksize*2,flag);

        eII2=matTOT2.inverse()*eII2;
*/
        
        //int flag=ArrangeDRE(D1value,DL);
        
        //Eigen::MatrixXcd DL=D.adjoint();
        //Eigen::MatrixXcd DR=D.inverse();

        //print_matrix(matTOT,"matTOT");
       /* 
        matTOT.block(0,0,bksize,bksize)=DL.topLeftCorner(bksize,bksize);
        matTOT.block(0,bksize,bksize,bksize)=DL.topRightCorner(bksize,bksize);
        matTOT.block(matTOT.rows()-bksize,matTOT.cols()-2*bksize,bksize,bksize)=DR.bottomLeftCorner(bksize,bksize);
        matTOT.block(matTOT.rows()-bksize,matTOT.cols()-bksize,bksize,bksize)=DR.bottomRightCorner(bksize,bksize);
*/
        //matTOT=matTOT.block(bksize,bksize,matTOT.rows()-2.*bksize,matTOT.cols()-2.*bksize);
        
        Eigen::MatrixXcd DL=D1vector.inverse();
        Eigen::MatrixXcd DR=D2vector.inverse();
        matTOT.block(0,0,bksize,bksize)=DL.topLeftCorner(bksize,bksize);
        matTOT.block(0,bksize,bksize,bksize)=DL.topRightCorner(bksize,bksize);
        matTOT.block(matTOT.rows()-bksize,matTOT.cols()-2*bksize,bksize,bksize)=DR.bottomLeftCorner(bksize,bksize);
        matTOT.block(matTOT.rows()-bksize,matTOT.cols()-bksize,bksize,bksize)=DR.bottomRightCorner(bksize,bksize);
        //print_matrix(matTOT,"matTOT--after");

        eII.resize(matTOT.rows(),flag);
        eII.setZero();
        eII.block(0,0,flag,flag)=Eigen::MatrixXcd::Identity(flag,flag);

        if(eII.cols()!=0){
            matII=matTOT.inverse()*eII;

            Eigen::MatrixXcd ttt=DR.block(0,0,bksize,2.*bksize)*matII.block(matII.rows()-2*bksize,0,2.*bksize,matII.cols());
            Eigen::MatrixXcd rrr=DL.block(bksize,0,bksize,2.*bksize)*matII.block(0,0,2.*bksize,matII.cols());

            //print_matrix(tt-ttt,"t--t");
            //print_matrix(rr-rrr,"r--r");

            rr=rrr;
            tt=ttt;

            //rr=DR.block(0,0,bksize,2.*bksize)*matII.block(0,0,2.*bksize,matII.cols());
            //tt=DL.block(bksize,0,bksize,2.*bksize)*matII.block(matII.rows()-2*bksize,0,2.*bksize,matII.cols());
            
        }
        //print_matrix(rr*rr.adjoint(),"rr");
        //print_matrix(tt*tt.adjoint(),"tt");
    }
    void Hamiltonian::DataCom(int blposX,int blposY,int bN,int bM,int N,int M,int Nb,int Mb,Eigen::MatrixXcd &blockMatrix,Eigen::Matrix<lapack_complex_double,-1,-1> &mat_loc,int &nrows,int &ncols,int &myrow,int &mycol,double flag,int &ctxt)
    {

        int id,numprocs;
        blacs_pinfo_(&id, &numprocs);
        int procrows = sqrt(numprocs), proccols = sqrt(numprocs);
        bool mpiroot = (id == 0);

        std::complex<double> ii(0,1);
        int ZEROi = 0;
        int MONEi=-1; 
        int ONEi=1;    
        double ONEd=1.0;
        double ZEROd=0.0;
        if(flag){
            mat_loc.resize(nrows,ncols);
            mat_loc.setConstant({0,0});
            int sendr = 0, sendc = 0, recvr = 0, recvc = 0;
            for (int r = 0; r < N; r += Nb ) {
                int nr = std::min(Nb,N-r);
                //sendc = 0;
                for (int c = 0; c < M; c += Mb) {
                    int nc = std::min(Mb,M-c);

                    if (mpiroot){
                        double tmp=blockMatrix(r,c).imag();
                        dgesd2d_(&ctxt, &nr, &nc, &tmp, &N, &sendr, &sendc);
                    }
                    if (myrow == sendr && mycol == sendc){
                        double tmp=0;
                        dgerv2d_(&ctxt, &nr, &nc, &tmp, &nrows,&ZEROi,&ZEROi);
                        mat_loc(recvr,recvc).imag+=tmp;
                    }
                    if (mpiroot){
                        double tmp=blockMatrix(r,c).real();
                        dgesd2d_(&ctxt, &nr, &nc, &tmp, &N, &sendr, &sendc);
                    }
                    if (myrow == sendr && mycol == sendc){
                        double tmp=0;
                        dgerv2d_(&ctxt, &nr, &nc, &tmp, &nrows,&ZEROi,&ZEROi);
                        mat_loc(recvr,recvc).real+=tmp;
                    }
                    if (myrow == sendr && mycol == sendc){
                        recvc = (recvc+nc)%ncols;
                    }

                    //sendc=(sendc+1)%proccols;
                    sendc=(c+1)%proccols;
                }
                if (myrow == sendr)
                {
                    recvr = (recvr+nr)%nrows;
                }
                //sendr=(sendr+1)%procrows;
                sendr=(r+1)%procrows;
            }
        }else{
            blockMatrix.resize(bN,bM);
            //blockMatrix.resize(N,M);
            blockMatrix.setZero();
            int sendr = 0, sendc = 0, recvr = 0, recvc = 0;
            //for (int r = blposX; r < bN; r += Nb ) {
            for (int r = 0; r < N; r += Nb ) {
                int nr = std::min(Nb,N-r);
                sendc = 0;
                //for (int c = blposY; c < bM; c += Mb) {
                for (int c = 0; c < M; c += Mb) {
                    int nc = std::min(Mb,M-c);

                    ///rrrr00
                    if((r-blposX)>=0 && (r-blposX)<bN && (c-blposY)>=0&& (c-blposY)<bM){
                        for(int pr=0;pr<procrows;++pr){
                            for(int pc=0;pc<proccols;++pc){

                                if (myrow == sendr && mycol == sendc){
                                    dgesd2d_(&ctxt, &nr, &nc, &mat_loc(recvr,recvc).imag, &nrows,&pr,&pc);
                                }
                                if (myrow == pr && mycol == pc){
                                    double tmp=0;
                                    dgerv2d_(&ctxt, &nr, &nc, &tmp, &N, &sendr, &sendc);
                                    blockMatrix(r-blposX,c-blposY)+=tmp*ii;
                                    //blockMatrix(r,c)+=tmp*ii;
                                }

                                if (myrow == sendr && mycol == sendc){
                                    dgesd2d_(&ctxt, &nr, &nc, &mat_loc(recvr,recvc).real, &nrows,&pr,&pc);
                                }
                                if (myrow == pr && mycol == pc){
                                    double tmp=0;
                                    dgerv2d_(&ctxt, &nr, &nc, &tmp, &N, &sendr, &sendc);
                                    blockMatrix(r-blposX,c-blposY)+=tmp;
                                    //blockMatrix(r,c)+=tmp;
                                }

                            }
                        }
                    } 

                    if (myrow == sendr && mycol == sendc){
                        recvc = (recvc+nc)%ncols;
                    }

                    sendc=(sendc+1)%proccols;
                    //sendc=(c+1)%proccols;
                }
                if (myrow == sendr)
                {
                    recvr = (recvr+nr)%nrows;
                }
                sendr=(sendr+1)%procrows;
                //sendr=(r+1)%procrows;
            }
            }


        }
        void Hamiltonian::ScatterHamiltonian(Eigen::MatrixXcd &blockMatrix,Eigen::Matrix<lapack_complex_double,-1,-1> &mat_loc,int tot,int Nb,int Mb)
        {

            std::complex<double> ii(0,1);
            int ZEROi = 0;
            int MONEi=-1; 
            int ONEi=1;    
            double ONEd=1.0;
            double ZEROd=0.0;

            int N=tot; 
            int M=tot; 

            int id,numprocs;
            blacs_pinfo_(&id, &numprocs);
            int ctxt, myrow, mycol;
            blacs_get_(&MONEi, &ZEROi, &ctxt);
            bool mpiroot = (id == 0);
            int info; 

            int procrows = sqrt(numprocs), proccols = sqrt(numprocs);
            blacs_gridinit_(&ctxt, "R", &procrows, &proccols);  
            blacs_gridinfo_(&ctxt, &procrows, &proccols, &myrow, &mycol );  
            int nrows = numroc_(&N, &Nb, &myrow, &ZEROi, &procrows);
            int lda = std::max(1,nrows);
            int ncols = numroc_(&M, &Mb, &mycol, &ZEROi, &proccols);


            mat_loc.resize(nrows,ncols);
            mat_loc.setConstant({0,0});
            DataCom(0,0,N,M,N,M,Nb,Mb,blockMatrix,mat_loc,nrows,ncols,myrow,mycol,1,ctxt);

            /*
               int sendr = 0, sendc = 0, recvr = 0, recvc = 0;
               for (int r = 0; r < N; r += Nb ) {
               int nr = std::min(Nb,N-r);
            //sendc = 0;
            for (int c = 0; c < M; c += Mb) {
            int nc = std::min(Mb,M-c);

            if (mpiroot){
            double tmp=blockMatrix(r,c).imag();
            dgesd2d_(&ctxt, &nr, &nc, &tmp, &N, &sendr, &sendc);
            }
            if (myrow == sendr && mycol == sendc){
            double tmp=0;
            dgerv2d_(&ctxt, &nr, &nc, &tmp, &nrows,&ZEROi,&ZEROi);
            mat_loc(recvr,recvc)+=tmp*ii;
            }
            if (mpiroot){
            double tmp=blockMatrix(r,c).real();
            dgesd2d_(&ctxt, &nr, &nc, &tmp, &N, &sendr, &sendc);
            }
            if (myrow == sendr && mycol == sendc){
            double tmp=0;
            dgerv2d_(&ctxt, &nr, &nc, &tmp, &nrows,&ZEROi,&ZEROi);
            mat_loc(recvr,recvc)+=tmp;
            }
            if (myrow == sendr && mycol == sendc){
            recvc = (recvc+nc)%ncols;
            }

            sendc=(sendc+1)%proccols;
            }
            if (myrow == sendr)
            {
            recvr = (recvr+nr)%nrows;
            }
            sendr=(sendr+1)%procrows;
            }

*/
            /*
               blacs_barrier_(&ctxt,"A");
               if(id==0)
               print_matrix(mat_loc,"A_loc0");
               blacs_barrier_(&ctxt,"A");
               if(id==1)
               print_matrix(mat_loc,"A_loc1");
               blacs_barrier_(&ctxt,"A");
               if(id==2)
               print_matrix(mat_loc,"A_loc2");
               blacs_barrier_(&ctxt,"A");
               if(id==3)
               print_matrix(mat_loc,"A_loc3");
               blacs_barrier_(&ctxt,"A");
               */
            blacs_gridexit_(&ctxt);

        }


        }
