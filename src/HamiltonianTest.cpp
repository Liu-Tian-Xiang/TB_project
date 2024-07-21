
#include "Hamiltonian.hpp"
namespace TightBinding
{
/*
    void Hamiltonian::NewSetMatrix(const Eigen::Vector3d k,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::MatrixXcd &matrix)
    {
        //std::unordered_map<std::string, Material> &ms=material.materialsShort;
        std::unordered_map<std::string, Material> &ms=material.materials;
        int Atom=MatrixSize;
        matrix.resize(Atom,Atom);	
        matrix.setZero();
        for(auto &cell : ListofMatrixElements)
        {
            
            //if(cell->expflag==0)
            {
                std::complex<double> ii(0.0,1.0);
                Eigen::Vector3d dd=cell->d;
                Eigen::Vector3d ed=cell->ed;
                Eigen::Vector3d edn=ed.normalized();
                Eigen::Vector3d dn=dd.normalized();
                //std::cout<<"dd="<<dd<<std::endl;
                //double phaseFactor=(2*M_PI*dd).dot(k);
                //std::complex<double> efactor=exp(phaseFactor*ii);
                std::complex<double> efactor=exp(2.*M_PI*((dd).dot(k))*ii);




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
                            std::string QWMaterial="GaAs";
                            std::string EBMaterial="AlAs";



                            if(rows==cols)
                            {
                                if(cell->cellme == cell->SuperCell)
                                {
                                    if(cell->cellme==QWMaterial)
                                    {
                                        matrix(rows,cols)= ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType)+BandOffSet;
                                    }else{
                                        matrix(rows,cols)= ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType);
                                    }
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
                                matrix(rows,cols)+=res*efactor;
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
        for(auto &cell : ListofMatrixElements)
        {
            
            //if(cell->expflag==0)
            {
                std::complex<double> ii(0.0,1.0);

                int atom2=cell->Table.real();
                int atom1=cell->Table.imag();
                std::string leftType=(atom1%2==0?"c":"a");

                std::string rightType=(atom2%2==0?"c":"a");
                atom2*=orbNum;
                atom1*=orbNum;
                for(int io=0;io<orbNum;++io)
                {
                    for(int jo=0;jo<orbNum;++jo)
                    {
                        Eigen::Vector4i leftorbital=qlabels[jo];
                        Eigen::Vector4i rightorbital=qlabels[io];
                        if(leftorbital(3)==rightorbital(3))
                        {
                            int rows=atom2+io;
                            int cols=atom1+jo;
                            std::string QWMaterial="GaAs";
                            std::string EBMaterial="AlAs";

                            if(cell->cellme!=cell->SuperCell){
                                if(rows>cols)
                                {
                                    matrix(cols,cols)=
                                        +ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType)*0.5
                                        +ms.at(cell->SuperCell).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType)*0.5+BandOffSet*0.5;
                                }
                            }
                        }
                    }
                }
            }

        }
    }
*/

    void Hamiltonian::SetMatrix(const Eigen::Vector3d k,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::MatrixXcd &matrix)
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
            std::complex<double> efactor=exp(2.*M_PI*((dd).dot(k))*ii);




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
}
