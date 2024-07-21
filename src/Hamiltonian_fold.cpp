
#include "Hamiltonian.hpp"
namespace TightBinding
{
    void Hamiltonian::vecSetMatrixFold(const Eigen::Vector3d k,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::MatrixXcd &matrix,std::string mat)
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
	    std::string ccellme=mat;
	    std::string cSuperCell=mat;

	    std::complex<double> ii(0.0,1.0);
	    Eigen::Vector3d dd=cell->d;
	    Eigen::Vector3d ed=cell->ed;
	    Eigen::Vector3d edn=ed.normalized();
	    Eigen::Vector3d dn=dd.normalized();


	    int atom2=cell->Table.real();
	    int atom1=cell->Table.imag();

	    int sign=cell->expflag;
	    Eigen::Vector3d dk={0,0,0};
	    //Eigen::Vector3d kk={0.3,0.2,0.2};
	    //Eigen::Vector3d pshift=Eigen::Vector3d{(sign)/systemcell.uvector.size(),0,0};
	    //Eigen::Vector3d inv_pshift=Eigen::Vector3d{((sign)/systemcell.uvector.size())==0?0:(1.0/((sign)/systemcell.uvector.size())),0,0};
	    Eigen::Vector3d pshift((sign)/systemcell.uvector.size(),0,0);
	    Eigen::Vector3d inv_pshift(((sign)/systemcell.uvector.size())==0?0:(1.0/((sign)/systemcell.uvector.size())),0,0);


	    Eigen::Vector3d ddd(0,0,0);
	    Eigen::Vector3d AtomPos(0,0,0);
	    systemcell.Unita1=systemcell.a1;
	    systemcell.Unita2=systemcell.a2;
	    systemcell.Unita3=systemcell.a3;
	    if(cell->Typeme<systemcell.supercell.size())
	    {
		ddd =systemcell.supercell[cell->Typeme]->Unit[sign*2]-systemcell.supercell[cell->Typeme]->Unit[0]-pshift;
		AtomPos =systemcell.supercell[cell->Typeme]->Unit[sign*2]-systemcell.supercell[cell->Typeme]->Unit[0]-Eigen::Vector3d(USE_Shift,0,0);
		//ddd-=supercell[0]->Unit[0];

		systemcell.Unita1=systemcell.a1*systemcell.supercell[cell->Typeme]->AtomNum/(systemcell.UnitCellNumber*systemcell.uvector.size());
	    }
	    int NNN=systemcell.supercell[cell->Typeme]->AtomNum/(systemcell.UnitCellNumber*systemcell.uvector.size());

	    double SVolumn=systemcell.Unita1.dot((systemcell.Unita2.cross(systemcell.Unita3)));
	    systemcell.Ug1=systemcell.Unita2.cross(systemcell.Unita3)/SVolumn;
	    systemcell.Ug2=systemcell.Unita3.cross(systemcell.Unita1)/SVolumn;
	    systemcell.Ug3=systemcell.Unita1.cross(systemcell.Unita2)/SVolumn;


	    SVolumn=systemcell.a1.dot((systemcell.a2.cross(systemcell.a3)));
	    Eigen::Vector3d ua1=systemcell.a2.cross(systemcell.a3)/SVolumn;
	    Eigen::Vector3d ua2=systemcell.a3.cross(systemcell.a1)/SVolumn;
	    Eigen::Vector3d ua3=systemcell.a1.cross(systemcell.a2)/SVolumn;


	    Eigen::Matrix3d Ugm,ugm,uam;
	    uam<<ua1,ua2,ua3;
	    Ugm<<systemcell.Ug1,systemcell.Ug2,systemcell.Ug3;
	    ugm<<systemcell.ug1,systemcell.ug2,systemcell.ug3;
	    Eigen::Matrix3d uga,Uga;
	    uga<<systemcell.u1,systemcell.u2,systemcell.u3;
	    Uga<<systemcell.Unita1,systemcell.Unita2,systemcell.Unita3;

	    Eigen::Vector3d projAtomPos=uam*AtomPos;
	    Eigen::Vector3d poshift((int)projAtomPos(0),(int)projAtomPos(1),(int)projAtomPos(2));

	    dk=ugm*(AtomPos-uam.inverse()*poshift)+Ugm*poshift;//dk2
	    /*
	       if(AtomPos(1)==0.5 &&AtomPos(2)==0.5)
	    //if((AtomPos-poshift).cross(systemcell.u1)==Eigen::Vector3d(0,0,0))
	    {
	    dk=Ugm*ugm*(AtomPos+(NNN-1)*(AtomPos-pshift));//dk1
	    }else{
	    dk=Ugm*ugm*(AtomPos);//dk1
	    }
	     */
	    std::string leftType=(atom1%2==0?"c":"a");
	    std::string rightType=(atom2%2==0?"c":"a");
	    std::complex<double> efactor=exp(2*M_PI*(dd.dot(k+dk))*ii);
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
			if(rows==cols)
			{
			    if((ccellme==QWMaterial && cSuperCell==EBMaterial) || (ccellme==EBMaterial && cSuperCell==QWMaterial))
			    {
				matrix(rows,cols)=ms.at(ccellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType);
			    }
			    else if(ccellme==QWMaterial && cSuperCell==QWMaterial)
			    {
				matrix(rows,cols)=ms.at(ccellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType);//+BandOffSet;
			    }else{
				matrix(rows,cols)=ms.at(ccellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType);
			    }

			}
			else if(leftType!=rightType)
			{ 

			    std::complex<double> res=(SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(ccellme))*0.5+SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cSuperCell))*0.5);
			    res*=efactor;
			    matrix(rows,cols)+=res;
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
				double lambda=ms.at(ccellme).param.at(leftType);
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
			    lambda=ms.at(ccellme).param.at(leftType);
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

    void Hamiltonian::SetMatrixPrime(const Eigen::Vector3d k,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::MatrixXcd &matrix,int vblpoi,int blen)
    {

        std::unordered_map<std::string, Material> &ms=material.materials;
        if(USE_long_parameters==0)
        {
            ms=material.materialsShort;
        }
        //int Atom=MatrixSize;
	//matrix.resize(Atom,Atom);	
	matrix.resize(blen,blen);	
	matrix.setZero();
	std::string ccellme;//="AlGaAs";
	std::string cSuperCell;//="AlGaAs";

	for(auto &cell : ListofMatrixElements)
	{

	    if(USE_MaterialMix)
	    {
		std::string ccellme=cell->cellme;
		std::string cSuperCell=cell->SuperCell;
		cell->cellme="AlGaAs";
		cell->SuperCell="AlGaAs";
		//if(cell->cellme=="GaAs") cell->cellme="GaAs";
		//else if(cell->cellme=="AlAs") cell->cellme="AlAs";
		//if(cell->SuperCell=="GaAs") cell->SuperCell="GaAs";
		//else if(cell->SuperCell=="AlAs") cell->SuperCell="AlAs";
	    }else{
		std::string ccellme=cell->cellme;
		std::string cSuperCell=cell->SuperCell;
		//if(cell->cellme=="GaAs") cell->cellme="AlAs";
		//else if(cell->cellme=="AlAs") cell->cellme="GaAs";
		//if(cell->SuperCell=="GaAs") cell->SuperCell="AlAs";
		//else if(cell->SuperCell=="AlAs") cell->SuperCell="GaAs";
	    }
	    std::complex<double> ii(0.0,1.0);
	    Eigen::Vector3d dd=cell->d;
	    Eigen::Vector3d ed=cell->ed;
	    Eigen::Vector3d edn=ed.normalized();
	    Eigen::Vector3d dn=dd.normalized();


	    int atom2=cell->Table.real();
	    int atom1=cell->Table.imag();

	    int sign=cell->expflag;
	    Eigen::Vector3d dk={0,0,0};
	    //Eigen::Vector3d kk={0.3,0.2,0.2};
	    //Eigen::Vector3d pshift=Eigen::Vector3d{(sign)/systemcell.uvector.size(),0,0};
	    //Eigen::Vector3d inv_pshift=Eigen::Vector3d{((sign)/systemcell.uvector.size())==0?0:(1.0/((sign)/systemcell.uvector.size())),0,0};
	    Eigen::Vector3d pshift((sign)/systemcell.uvector.size(),0,0);
	    Eigen::Vector3d inv_pshift(((sign)/systemcell.uvector.size())==0?0:(1.0/((sign)/systemcell.uvector.size())),0,0);


	    Eigen::Vector3d ddd(0,0,0);
	    Eigen::Vector3d AtomPos(0,0,0);
	    systemcell.Unita1=systemcell.a1;
	    systemcell.Unita2=systemcell.a2;
	    systemcell.Unita3=systemcell.a3;
	    if(cell->Typeme<systemcell.supercell.size())
	    {
		ddd =systemcell.supercell[cell->Typeme]->Unit[sign*2]-systemcell.supercell[cell->Typeme]->Unit[0]-pshift;
		AtomPos =systemcell.supercell[cell->Typeme]->Unit[sign*2]-systemcell.supercell[cell->Typeme]->Unit[0]-Eigen::Vector3d(USE_Shift,0,0);
		//ddd-=supercell[0]->Unit[0];

		systemcell.Unita1=systemcell.a1*systemcell.supercell[cell->Typeme]->AtomNum/(systemcell.UnitCellNumber*systemcell.uvector.size());
	    }
	    int NNN=systemcell.supercell[cell->Typeme]->AtomNum/(systemcell.UnitCellNumber*systemcell.uvector.size());

	    double SVolumn=systemcell.Unita1.dot((systemcell.Unita2.cross(systemcell.Unita3)));
	    systemcell.Ug1=systemcell.Unita2.cross(systemcell.Unita3)/SVolumn;
	    systemcell.Ug2=systemcell.Unita3.cross(systemcell.Unita1)/SVolumn;
	    systemcell.Ug3=systemcell.Unita1.cross(systemcell.Unita2)/SVolumn;


	    SVolumn=systemcell.a1.dot((systemcell.a2.cross(systemcell.a3)));
	    Eigen::Vector3d ua1=systemcell.a2.cross(systemcell.a3)/SVolumn;
	    Eigen::Vector3d ua2=systemcell.a3.cross(systemcell.a1)/SVolumn;
	    Eigen::Vector3d ua3=systemcell.a1.cross(systemcell.a2)/SVolumn;


	    Eigen::Matrix3d Ugm,ugm,uam;
	    uam<<ua1,ua2,ua3;
	    Ugm<<systemcell.Ug1,systemcell.Ug2,systemcell.Ug3;
	    ugm<<systemcell.ug1,systemcell.ug2,systemcell.ug3;
	    Eigen::Matrix3d uga,Uga;
	    uga<<systemcell.u1,systemcell.u2,systemcell.u3;
	    Uga<<systemcell.Unita1,systemcell.Unita2,systemcell.Unita3;

	    Eigen::Vector3d projAtomPos=uam*AtomPos;
	    Eigen::Vector3d poshift((int)projAtomPos(0),(int)projAtomPos(1),(int)projAtomPos(2));

	    dk=ugm*(AtomPos-uam.inverse()*poshift)+Ugm*poshift;//dk2
	    /*
	       if(AtomPos(1)==0.5 &&AtomPos(2)==0.5)
	    //if((AtomPos-poshift).cross(systemcell.u1)==Eigen::Vector3d(0,0,0))
	    {
	    dk=Ugm*ugm*(AtomPos+(NNN-1)*(AtomPos-pshift));//dk1
	    }else{
	    dk=Ugm*ugm*(AtomPos);//dk1
	    }
	     */
	    std::string leftType=(atom1%2==0?"c":"a");
	    std::string rightType=(atom2%2==0?"c":"a");
	    std::complex<double> efactor=exp(2*M_PI*(dd.dot(k+dk))*ii);
	    atom2*=orbNum;
	    atom1*=orbNum;
	    for(int io=0;io<orbNum;++io)
	    {
		for(int jo=0;jo<orbNum;++jo)
		{

		    Eigen::Vector4i leftorbital=qlabels[jo];
		    Eigen::Vector4i rightorbital=qlabels[io];
			int rows=atom2+io-vblpoi;
			int cols=atom1+jo-vblpoi;
                    if((rows>=0 && rows<blen) && (cols>=0 && cols<blen)) 
{
		    if(leftorbital(3)==rightorbital(3))
		    {
			std::string QWMaterial="GaAs";
			std::string EBMaterial="AlAs";
			if(rows==cols)
			{
			    if((ccellme==QWMaterial && cSuperCell==EBMaterial) || (ccellme==EBMaterial && cSuperCell==QWMaterial))
			    {
				matrix(rows,cols)=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType);
			    }
			    else if(ccellme==QWMaterial && cSuperCell==QWMaterial)
			    {
				matrix(rows,cols)=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType);//+BandOffSet;
			    }else{
				matrix(rows,cols)=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType);
			    }

			}
			else if(leftType!=rightType)
			{ 

			    std::complex<double> res=(SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->cellme))*0.5+SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->SuperCell))*0.5);
			    res*=efactor;
			    matrix(rows,cols)+=res;
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

    }


    void Hamiltonian::SetMatrixFold(const Eigen::Vector3d k,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::MatrixXcd &matrix)
    {

        std::unordered_map<std::string, Material> &ms=material.materials;
        if(USE_long_parameters==0)
        {
            ms=material.materialsShort;
        }
        int Atom=MatrixSize;
	matrix.resize(Atom,Atom);	
	matrix.setZero();
	std::string ccellme;//="AlGaAs";
	std::string cSuperCell;//="AlGaAs";

	for(auto &cell : ListofMatrixElements)
	{

	    if(USE_MaterialMix)
	    {
		std::string ccellme=cell->cellme;
		std::string cSuperCell=cell->SuperCell;
		cell->cellme="AlGaAs";
		cell->SuperCell="AlGaAs";
		//if(cell->cellme=="GaAs") cell->cellme="GaAs";
		//else if(cell->cellme=="AlAs") cell->cellme="AlAs";
		//if(cell->SuperCell=="GaAs") cell->SuperCell="GaAs";
		//else if(cell->SuperCell=="AlAs") cell->SuperCell="AlAs";
	    }else{
		std::string ccellme=cell->cellme;
		std::string cSuperCell=cell->SuperCell;
		//if(cell->cellme=="GaAs") cell->cellme="AlAs";
		//else if(cell->cellme=="AlAs") cell->cellme="GaAs";
		//if(cell->SuperCell=="GaAs") cell->SuperCell="AlAs";
		//else if(cell->SuperCell=="AlAs") cell->SuperCell="GaAs";
	    }
	    std::complex<double> ii(0.0,1.0);
	    Eigen::Vector3d dd=cell->d;
	    Eigen::Vector3d ed=cell->ed;
	    Eigen::Vector3d edn=ed.normalized();
	    Eigen::Vector3d dn=dd.normalized();


	    int atom2=cell->Table.real();
	    int atom1=cell->Table.imag();

	    int sign=cell->expflag;
	    Eigen::Vector3d dk={0,0,0};
	    //Eigen::Vector3d kk={0.3,0.2,0.2};
	    //Eigen::Vector3d pshift=Eigen::Vector3d{(sign)/systemcell.uvector.size(),0,0};
	    //Eigen::Vector3d inv_pshift=Eigen::Vector3d{((sign)/systemcell.uvector.size())==0?0:(1.0/((sign)/systemcell.uvector.size())),0,0};
	    Eigen::Vector3d pshift((sign)/systemcell.uvector.size(),0,0);
	    Eigen::Vector3d inv_pshift(((sign)/systemcell.uvector.size())==0?0:(1.0/((sign)/systemcell.uvector.size())),0,0);


	    Eigen::Vector3d ddd(0,0,0);
	    Eigen::Vector3d AtomPos(0,0,0);
	    systemcell.Unita1=systemcell.a1;
	    systemcell.Unita2=systemcell.a2;
	    systemcell.Unita3=systemcell.a3;
	    if(cell->Typeme<systemcell.supercell.size())
	    {
		ddd =systemcell.supercell[cell->Typeme]->Unit[sign*2]-systemcell.supercell[cell->Typeme]->Unit[0]-pshift;
		AtomPos =systemcell.supercell[cell->Typeme]->Unit[sign*2]-systemcell.supercell[cell->Typeme]->Unit[0]-Eigen::Vector3d(USE_Shift,0,0);
		//ddd-=supercell[0]->Unit[0];

		systemcell.Unita1=systemcell.a1*systemcell.supercell[cell->Typeme]->AtomNum/(systemcell.UnitCellNumber*systemcell.uvector.size());
	    }
	    int NNN=systemcell.supercell[cell->Typeme]->AtomNum/(systemcell.UnitCellNumber*systemcell.uvector.size());

	    double SVolumn=systemcell.Unita1.dot((systemcell.Unita2.cross(systemcell.Unita3)));
	    systemcell.Ug1=systemcell.Unita2.cross(systemcell.Unita3)/SVolumn;
	    systemcell.Ug2=systemcell.Unita3.cross(systemcell.Unita1)/SVolumn;
	    systemcell.Ug3=systemcell.Unita1.cross(systemcell.Unita2)/SVolumn;


	    SVolumn=systemcell.a1.dot((systemcell.a2.cross(systemcell.a3)));
	    Eigen::Vector3d ua1=systemcell.a2.cross(systemcell.a3)/SVolumn;
	    Eigen::Vector3d ua2=systemcell.a3.cross(systemcell.a1)/SVolumn;
	    Eigen::Vector3d ua3=systemcell.a1.cross(systemcell.a2)/SVolumn;


	    Eigen::Matrix3d Ugm,ugm,uam;
	    uam<<ua1,ua2,ua3;
	    Ugm<<systemcell.Ug1,systemcell.Ug2,systemcell.Ug3;
	    ugm<<systemcell.ug1,systemcell.ug2,systemcell.ug3;
	    Eigen::Matrix3d uga,Uga;
	    uga<<systemcell.u1,systemcell.u2,systemcell.u3;
	    Uga<<systemcell.Unita1,systemcell.Unita2,systemcell.Unita3;

	    Eigen::Vector3d projAtomPos=uam*AtomPos;
	    Eigen::Vector3d poshift((int)projAtomPos(0),(int)projAtomPos(1),(int)projAtomPos(2));

	    dk=ugm*(AtomPos-uam.inverse()*poshift)+Ugm*poshift;//dk2
	    /*
	       if(AtomPos(1)==0.5 &&AtomPos(2)==0.5)
	    //if((AtomPos-poshift).cross(systemcell.u1)==Eigen::Vector3d(0,0,0))
	    {
	    dk=Ugm*ugm*(AtomPos+(NNN-1)*(AtomPos-pshift));//dk1
	    }else{
	    dk=Ugm*ugm*(AtomPos);//dk1
	    }
	     */
	    std::string leftType=(atom1%2==0?"c":"a");
	    std::string rightType=(atom2%2==0?"c":"a");
	    std::complex<double> efactor=exp(2*M_PI*(dd.dot(k+dk))*ii);
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
			if(rows==cols)
			{
			    if((ccellme==QWMaterial && cSuperCell==EBMaterial) || (ccellme==EBMaterial && cSuperCell==QWMaterial))
			    {
				matrix(rows,cols)=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType);
			    }
			    else if(ccellme==QWMaterial && cSuperCell==QWMaterial)
			    {
				matrix(rows,cols)=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType);//+BandOffSet;
			    }else{
				matrix(rows,cols)=ms.at(cell->cellme).param.at(std::to_string(leftorbital(0))+std::to_string(leftorbital(1))+leftType+std::to_string(rightorbital(0))+std::to_string(rightorbital(1))+rightType);
			    }

			}
			else if(leftType!=rightType)
			{ 

			    std::complex<double> res=(SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->cellme))*0.5+SKFormular(leftType,rightType,leftorbital(0),leftorbital(1),leftorbital(2),rightorbital(0),rightorbital(1),rightorbital(2),dn,&ms.at(cell->SuperCell))*0.5);
			    res*=efactor;
			    matrix(rows,cols)+=res;
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
    void Hamiltonian::foldDiagonalize(std::vector<Eigen::MatrixXcd> &SubblockContainer,Eigen::MatrixXcd &ConsideredMatrix,bool onoff,Eigen::MatrixXcd &MatrixToBeTransformed,const Eigen::Vector3d& k,int dtag)
    {
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
			SubblockContainer.push_back(gsolver.eigenvectors().block(0,cell->Startp,cell->AtomO,cell->SecDim));
		    }else{
			SubblockContainer.push_back(gsolver.eigenvectors());
		    }
		    Atomtot+=cell->SecDim;
		    matrix_bottom+=cell->Startp;
		}else{

		    Eigen::MatrixXcd sfitting;
		    Eigen::VectorXd svalue;
		    QuickFoldTheBlockHamiltonian(cell->Tag2,cell->blpo,cell->AtomO,cell->vblpo,cell->vAtomO,k,cell->Unit,cell->EnergyCutoff,cell->Startp,cell->SecDim,dtag,svalue,sfitting);
		    //print_matrix(sfitting.adjoint()*matrix*sfitting,"testing");

		    SubblockContainer.push_back(sfitting);
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
    }
    /*
       void Hamiltonian::FoldTheBlockHamiltonian2(int cellme,int blpo,int AtomO,std::vector<int> vblpo,std::vector<int> vAtomO,const Eigen::Vector3d& k,std::vector<Eigen::Vector3d> chain,Eigen::MatrixXcd &fitting,Eigen::VectorXd &egv,int dtag)//with dtag
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
    //Eigen::VectorXd egv;
    //std::vector<Eigen::MatrixXcd> egvct;
    //egvct.resize(0);
    int blnum=vblpo.size();
    //Eigen::MatrixXcd fitting;

    if(USE_decomposition)
    {
    fitting.resize(AtomO,vAtomO[dtag]);
    }else{
    fitting.resize(AtomO,AtomO);
    }
    fitting.setZero();

    std::complex<double> ii(0,1);
    int fitposCol=0;

    //for(int i=0;i<blnum;++i)
    {
    int i=dtag;
    int blen=vAtomO[i];
    gsolver.compute(matrixScatter.block(vblpo[i],vblpo[i],blen,blen),Eigen::MatrixXcd::Identity(blen,blen) ,Eigen::DecompositionOptions::ComputeEigenvectors|Eigen::Ax_lBx);
    //egv.segment(vblpo[i]-vblpo[0],blen)=gsolver.eigenvalues();

    egv=gsolver.eigenvalues();

    //egv.insert(egv.end(),gsolver.eigenvalues().data(),gsolver.eigenvalues().data()+gsolver.eigenvalues().size());
    //egvct.push_back(gsolver.eigenvectors().transpose());
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
	Eigen::Vector3d posc=chain[(j)*UnitCellNumber+c]-shift;

	Eigen::Vector3d projAtomPos=uam*AtomPos;
	Eigen::Vector3d poshift((int)projAtomPos(0),(int)projAtomPos(1),(int)projAtomPos(2));

	dk=ugm*(AtomPos-uam.inverse()*poshift)+Ugm*poshift;//dk1
	//int IND=(j)*UnitCellNumber+c;
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
//fitposCol+=vAtomO[i];
}

//std::vector<std::complex<double>> ord;
//sfitting=fitting;


//return Startp;
}
    */
void Hamiltonian::QuickFoldTheBlockHamiltonian(int cellme,int blpo,int AtomO,std::vector<int> vblpo,std::vector<int> vAtomO,const Eigen::Vector3d& k,std::vector<Eigen::Vector3d> chain,std::pair<double,double> EnergyCutoff,int &Startp,int &SecDim,int dtag,Eigen::VectorXd &svalue,Eigen::MatrixXcd &yfitting)
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

    int printM=cellme;
    Eigen::Vector3d shift=systemcell.supercell[cellme]->Unit[0];
    assert(vblpo.size()==vAtomO.size());
    std::vector<Eigen::MatrixXcd> egvct;
    egvct.resize(0);

    Eigen::VectorXd svalue2;
    svalue2.resize(AtomO);
    svalue2.setZero();

    svalue.resize(AtomO);
    svalue.setZero();
    std::vector<std::complex<double>> ord;
    std::complex<double> ii(0.0,1.0);

    int blnum=vblpo.size();
    //std::vector<double> egv;
    //egv.resize(0);
    for(int i=0;i<blnum;++i)
    {
	int blen=vAtomO[i];
	if(fabs(dtag)!=999)
	{
	    gsolver.compute(matrixScatter.block(vblpo[i],vblpo[i],blen,blen),Eigen::MatrixXcd::Identity(blen,blen) ,Eigen::DecompositionOptions::ComputeEigenvectors|Eigen::Ax_lBx);
	    assert(gsolver.info() == Eigen::ComputationInfo::Success);
	}else if(dtag==999){
	    gsolver.compute(vecMatrixFold[0].block(vblpo[i],vblpo[i],blen,blen),Eigen::MatrixXcd::Identity(blen,blen) ,Eigen::DecompositionOptions::ComputeEigenvectors|Eigen::Ax_lBx);
	    assert(gsolver.info() == Eigen::ComputationInfo::Success);
	}else if(dtag==-999){
	    gsolver.compute(vecMatrixFold[1].block(vblpo[i],vblpo[i],blen,blen),Eigen::MatrixXcd::Identity(blen,blen) ,Eigen::DecompositionOptions::ComputeEigenvectors|Eigen::Ax_lBx);
	    assert(gsolver.info() == Eigen::ComputationInfo::Success);
	}
	//gsolver.compute(matrixScatter.block(vblpo[i],vblpo[i],blen,blen),Eigen::MatrixXcd::Identity(blen,blen) ,Eigen::DecompositionOptions::ComputeEigenvectors|Eigen::Ax_lBx);
	//assert(gsolver.info() == Eigen::ComputationInfo::Success);

	//egv.insert(egv.end(),gsolver.eigenvalues().data(),gsolver.eigenvalues().data()+gsolver.eigenvalues().size());
	svalue.segment(vblpo[i]-vblpo[0],blen)=gsolver.eigenvalues();
	//ord.insert(ord.end(),gsolver.eigenvalues().data(),gsolver.eigenvalues().data()+gsolver.eigenvalues().size());
	egvct.push_back(gsolver.eigenvectors());
    }
    svalue2=svalue;
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
		    FCol3(egvct[Fi].col(Fcc*orbNum+Fkj),cellme,blpo,AtomO,vblpo,vAtomO, k, chain,Startp+yi,dtag,Fcc,Fkj,Ugm,ugm,uam,yfitting);
		}
	    }else{
		//yfitting.resize(AtomO,leftCutoff+rightCutoff);
		yfitting.resize(AtomO,AtomO);
		yfitting.setZero();
		Startp=AtomO*BandGapIndex/orbNum-leftCutoff;
		SecDim=rightCutoff+leftCutoff;
		//for(int yi=0;yi<yfitting.cols();++yi)
		for(int yi=0;yi<SecDim;++yi)
		{
		    //yfitting.col(yi)=fitting.col(ord[AtomO*BandGapIndex/orbNum+yi-leftCutoff].real());
		    int colnum=ord[AtomO*BandGapIndex/orbNum+yi-leftCutoff].real();
		    int Fi=colnum/(systemcell.UnitCellNumber*orbNum);
		    int Fcc=colnum%(systemcell.UnitCellNumber*orbNum)/orbNum;
		    int Fkj=colnum%(systemcell.UnitCellNumber*orbNum)%orbNum;
		    FCol3(egvct[Fi].col(Fcc*orbNum+Fkj),cellme,blpo,AtomO,vblpo,vAtomO, k, chain,Startp+yi,Fi,Fcc,Fkj,Ugm,ugm,uam,yfitting);
		}
	    }
	}else{
	    //NoteUSE_MoreRanges FixedSize
	    if(USE_decomposition && fabs(dtag)!=999)
	    {//place
		std::cout<<"USE_MoreRanges FixedSize Tag="<<vblpo.size()<<std::endl;
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
		    FCol3(egvct[Fi].col(Fcc*orbNum+Fkj),cellme,blpo,AtomO,vblpo,vAtomO, k, chain,yi,dtag,Fcc,Fkj,Ugm,ugm,uam,yfitting);

		}
		svalue=svalue2.segment(vblpo[dtag]-vblpo[0],vAtomO[dtag]);
	    }else{
		if(systemcell.supercell.size()!=1)
		{
		    SecDim=EnergyCutoff.first+EnergyCutoff.second;
		    Startp=(AtomO/orbNum)*BandGapIndex-EnergyCutoff.first;
		    yfitting.resize(AtomO,AtomO);
		    yfitting.setZero();
		    for(int yi=0;yi<yfitting.cols();++yi)
		    {
			int colnum=ord[yi].real();
			int Fi=colnum/(systemcell.UnitCellNumber*orbNum);
			int Fcc=colnum%(systemcell.UnitCellNumber*orbNum)/orbNum;
			int Fkj=colnum%(systemcell.UnitCellNumber*orbNum)%orbNum;
			FCol3(egvct[Fi].col(Fcc*orbNum+Fkj),cellme,blpo,AtomO,vblpo,vAtomO, k, chain,yi,Fi,Fcc,Fkj,Ugm,ugm,uam,yfitting);
		    }
		}else{
		    std::cout<<"tmpFTBH - USE_MoreRanges FixedSize ONE"<<std::endl;
		    tmpFTBH(cellme,blpo,AtomO,vblpo,vAtomO,k,chain,yfitting,Startp,SecDim,Ugm,ugm,uam);
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
		    FCol3(egvct[Fi].col(Fcc*orbNum+Fkj),cellme,blpo,AtomO,vblpo,vAtomO, k, chain,yi,dtag,Fcc,Fkj,Ugm,ugm,uam,yfitting);
		}
		svalue=svalue2.segment(vblpo[dtag]-vblpo[0],vAtomO[dtag]);
		const auto upper=std::upper_bound(svalue.data(),svalue.data()+svalue.size(),rightCutoff);
		const auto lower=std::lower_bound(svalue.data(),upper,leftCutoff);
		SecDim=std::distance(lower,upper);
		Startp=std::distance(svalue.data(),lower);
	    }else{
		const auto upper=std::upper_bound(svalue.data(),svalue.data()+svalue.size(),rightCutoff);
		const auto lower=std::lower_bound(svalue.data(),upper,leftCutoff);
		SecDim=std::distance(lower,upper);
		Startp=std::distance(svalue.data(),lower);
		yfitting.resize(AtomO,AtomO);
		yfitting.setZero();
		for(int yi=0;yi<yfitting.cols();++yi)
		{
		    int colnum=ord[yi].real();
		    int Fi=colnum/(systemcell.UnitCellNumber*orbNum);
		    int Fcc=colnum%(systemcell.UnitCellNumber*orbNum)/orbNum;
		    int Fkj=colnum%(systemcell.UnitCellNumber*orbNum)%orbNum;
		    FCol3(egvct[Fi].col(Fcc*orbNum+Fkj),cellme,blpo,AtomO,vblpo,vAtomO, k, chain,yi,Fi,Fcc,Fkj,Ugm,ugm,uam,yfitting);
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
		    FCol3(egvct[Fi].col(Fcc*orbNum+Fkj),cellme,blpo,AtomO,vblpo,vAtomO, k, chain,yi,dtag,Fcc,Fkj,Ugm,ugm,uam,yfitting);
		}
		svalue=svalue2.segment(vblpo[dtag]-vblpo[0],vAtomO[dtag]);
		std::cout<<"USE_MoreRanges EnergyRange decomposition "<<std::endl;
		const auto upper=std::upper_bound(svalue.data(),svalue.data()+svalue.size(),EnergyCutoff.second);
		const auto lower=std::lower_bound(svalue.data(),upper,EnergyCutoff.first);
		SecDim=std::distance(lower,upper);
		Startp=std::distance(svalue.data(),lower);
	    }else{
		if(systemcell.supercell.size()!=1)
		{
		    const auto upper=std::upper_bound(svalue.data(),svalue.data()+svalue.size(),EnergyCutoff.second);
		    const auto lower=std::lower_bound(svalue.data(),upper,EnergyCutoff.first);
		    SecDim=std::distance(lower,upper);
		    Startp=std::distance(svalue.data(),lower);
		    yfitting.resize(AtomO,AtomO);
		    yfitting.setZero();
		    for(int yi=0;yi<yfitting.cols();++yi)
		    {
			int colnum=ord[yi].real();
			int Fi=colnum/(systemcell.UnitCellNumber*orbNum);
			int Fcc=colnum%(systemcell.UnitCellNumber*orbNum)/orbNum;
			int Fkj=colnum%(systemcell.UnitCellNumber*orbNum)%orbNum;
			FCol3(egvct[Fi].col(Fcc*orbNum+Fkj),cellme,blpo,AtomO,vblpo,vAtomO, k, chain,yi,Fi,Fcc,Fkj,Ugm,ugm,uam,yfitting);
		    }
		}else{
		    std::cout<<"FTBH - USE_MoreRanges EnergyRange ONE NODEcomposition"<<std::endl;
		    FTBH(cellme,blpo,AtomO,vblpo,vAtomO,k,chain,yfitting,Startp,SecDim,Ugm,ugm,uam,egvct,svalue2);
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
void Hamiltonian::tmpFTBH(int cellme,int blpo,int AtomO,std::vector<int> vblpo,std::vector<int> vAtomO,const Eigen::Vector3d& k,std::vector<Eigen::Vector3d> chain,Eigen::MatrixXcd &fitting,int &Startp,int &SecDim,Eigen::Matrix3d &Ugm,Eigen::Matrix3d &ugm,Eigen::Matrix3d &uam)//without dtag
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

    fitting.resize(AtomO,AtomO);
    fitting.setZero();
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
		gsolver.compute(matrixScatter.block(vblpo[gi],vblpo[gi],blen,blen),Eigen::MatrixXcd::Identity(blen,blen) ,Eigen::DecompositionOptions::ComputeEigenvectors|Eigen::Ax_lBx);
		assert(gsolver.info() == Eigen::ComputationInfo::Success);
		//Eigen::VectorXd gtmpegv=gsolver.eigenvalues();
		tmpegv.segment(flag*blen,blen)=gsolver.eigenvalues();
		egvct.push_back(gsolver.eigenvectors());

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
			FCol3(egvct[((int)(ord[flag*blen+in].real()))/blen].col(Fcc*orbNum+Fkj),cellme,blpo,AtomO,vblpo,vAtomO, k, chain,Startp+SecDim,Fi,Fcc,Fkj,Ugm,ugm,uam,fitting);
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
void Hamiltonian::FTBH(int cellme,int blpo,int AtomO,std::vector<int> vblpo,std::vector<int> vAtomO,const Eigen::Vector3d& k,std::vector<Eigen::Vector3d> chain,Eigen::MatrixXcd &fitting,int &Startp,int &SecDim,Eigen::Matrix3d &Ugm,Eigen::Matrix3d &ugm,Eigen::Matrix3d &uam,std::vector<Eigen::MatrixXcd> &egvct,Eigen::VectorXd &svalue2)//without dtag
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

    fitting.resize(AtomO,AtomO);
    fitting.setZero();

    Startp=0;
    int gi=0;
    for(int i1=0;i1<LAtomO.size();++i1)
    {
	for(int t=0;t<LAtomO[i1];++t)
	{
	    for(int i2=0;i2<systemcell.uvector.size();++i2)
	    {
		int blen=vAtomO[gi];
		Eigen::VectorXd tmpGsolver=svalue2.segment(vblpo[gi]-vblpo[0],blen);
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
		Eigen::VectorXd tmpGsolver=svalue2.segment(vblpo[gi]-vblpo[0],blen);
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
			FCol3(egvct[gi].col(Fcc*orbNum+Fkj),cellme,blpo,AtomO,vblpo,vAtomO, k, chain,Startp+SecDim,Fi,Fcc,Fkj,Ugm,ugm,uam,fitting);
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
/*
   int Hamiltonian::FoldTheCol3(std::vector<double> &egv,std::vector<Eigen::MatrixXcd> &egvct,int cellme,int blpo,int AtomO,std::vector<int> vblpo,std::vector<int> vAtomO,const Eigen::Vector3d& k,std::vector<Eigen::Vector3d> chain,std::pair<double,double> EnergyCutoff,int newcolnum,int colnum,Eigen::Matrix3d &Ugm,Eigen::Matrix3d &ugm,Eigen::Matrix3d &uam,Eigen::MatrixXcd &yfitting,int DataTag)
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
//for(int i=0;i<blnum;++i)
{
int i=colnum/(UnitCellNumber*orbNum);
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
//for(int cc=0;cc<UnitCellNumber;++cc)
{
int cc=colnum%(UnitCellNumber*orbNum)/orbNum;
//for(int kj=0;kj<orbNum;++kj)
{
int kj=colnum%(UnitCellNumber*orbNum)%orbNum;
int colIndex=cc*orbNum+kj;

{

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
	yfitting(fitposRow+rowIndex,newcolnum)=egvct[i](rowIndex,colIndex)*expr/sqrt(blnum);//5
	retvalue=1;
    }else{

	if(
		(fitposRow+rowIndex)>=(LpV)*cNUM
		&& (fitposRow+rowIndex)<(LpV+systemcell.spNum[minDistance+1])*cNUM
	  )
	{
	    yfitting(fitposRow+rowIndex,newcolnum)=egvct[i](rowIndex,colIndex)*expr/sqrt(blnum);//5
	    retvalue=1;
	}
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
void Hamiltonian::FCol3(const Eigen::VectorXcd &egvct,int cellme,int blpo,int AtomO,std::vector<int> vblpo,std::vector<int> vAtomO,const Eigen::Vector3d& k,std::vector<Eigen::Vector3d> chain,int newcolnum,int i,int cc,int kj,Eigen::Matrix3d &Ugm,Eigen::Matrix3d &ugm,Eigen::Matrix3d &uam,Eigen::MatrixXcd &yfitting)
{
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
		    yfitting(fitposRow+rowIndex,newcolnum)=egvct(rowIndex)*expr/sqrt(blnum);//5
		}else{

		    if(
			    (fitposRow+rowIndex)>=(LpV)*cNUM
			    && (fitposRow+rowIndex)<(LpV+systemcell.spNum[minDistance+1])*cNUM
		      )
		    {
			yfitting(fitposRow+rowIndex,newcolnum)=egvct(rowIndex)*expr/sqrt(blnum);//5
		    }
		}

	    }
	}
	fitposRow+=vAtomO[j];
    }
}
/*
   int Hamiltonian::FoldTheCol3(const Eigen::MatrixXcd &egvct,int cellme,int blpo,int AtomO,std::vector<int> vblpo,std::vector<int> vAtomO,const Eigen::Vector3d& k,std::vector<Eigen::Vector3d> chain,std::pair<double,double> EnergyCutoff,int newcolnum,int colnum,Eigen::Matrix3d &Ugm,Eigen::Matrix3d &ugm,Eigen::Matrix3d &uam,Eigen::MatrixXcd &yfitting,int DataTag)
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
//for(int i=0;i<blnum;++i)
{
int i=colnum/(UnitCellNumber*orbNum);
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
//for(int cc=0;cc<UnitCellNumber;++cc)
{
int cc=colnum%(UnitCellNumber*orbNum)/orbNum;
//for(int kj=0;kj<orbNum;++kj)
{
int kj=colnum%(UnitCellNumber*orbNum)%orbNum;
int colIndex=cc*orbNum+kj;

{

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
	yfitting(fitposRow+rowIndex,newcolnum)=egvct(rowIndex,colIndex)*expr/sqrt(blnum);//5
	retvalue=1;
    }else{

	if(
		(fitposRow+rowIndex)>=(LpV)*cNUM
		&& (fitposRow+rowIndex)<(LpV+systemcell.spNum[minDistance+1])*cNUM
	  )
	{
	    yfitting(fitposRow+rowIndex,newcolnum)=egvct(rowIndex,colIndex)*expr/sqrt(blnum);//5
	    retvalue=1;
	}
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

}
