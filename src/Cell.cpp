#include "Cell.hpp"


namespace TightBinding
{
    void Cell::ProcessCell()
    {
        for(int i=0;i<spJoint.size();++i)
        {
            int jsize=(spJoint[i].imag()==0?spJoint[i].imag()+1:spJoint[i].imag());
            for(int j=0;j<(uvector.size()*jsize-1);++j)
            {
                JointSupercell(i);
            }
        }

        for(int i=0;i<TypeCollection.size();++i)
        {
            int jsize=(TypeCollection[i].imag()==0?TypeCollection[i].imag()+1:TypeCollection[i].imag());
            for(int j=0;j<(jsize-1);++j)
            {
                //JointSupercell(i);
                JointSupercell2(i);
            }
        }
    }
    
        void Cell::ProcessParameters(int paramd,int InAsLength,double leftCutoff,double rightCutoff)
        {
            dest=0;
            if(paramd==1) dest=paramd+2+3;
            else if(pow(2,DividedFactor)==paramd) dest=paramd+2;
            else dest=paramd+2+2;
            spNum.clear();
            spNum.resize(dest);
            int spbn=pow(2,DividedFactor)/paramd;
            modInAsLength=InAsLength%spbn;
            floorInAsLength=InAsLength/spbn;
            if(spbn==1){floorInAsLength--;}
            for(int i=1;i<dest-1;++i)
            {
                spNum[i]=spbn;
            }
            if(paramd==1)
            {
                spNum[1]/=2;
                spNum[1]-=InAsLength;
                spNum[2]=InAsLength;
                spNum[3]=InAsLength;
                spNum[4]/=2;
                spNum[4]-=InAsLength;
            }
            else if(pow(2,DividedFactor)!=paramd)
            {
                int halfd=dest/2;
                if(spNum[halfd-2-floorInAsLength]>modInAsLength)
                {
                    spNum[halfd-2-floorInAsLength]-=modInAsLength;
                }
                spNum[halfd-1-floorInAsLength]=modInAsLength;
                spNum[halfd+floorInAsLength]=modInAsLength;
                if(spNum[halfd+1+floorInAsLength]>modInAsLength)
                {
                    spNum[halfd+1+floorInAsLength]-=modInAsLength;
                }
            }


            for(auto cell=spNum.begin()+1;cell!=spNum.end()-1;++cell)
            {
                if((*cell)==0)
                {
                    spNum.erase(cell);
                    --cell;
                }
            }

            dest=spNum.size();
            spNum[0]=spNum[dest-2];
            spNum[dest-1]=spNum[1];


            spStr.resize(0);
            int correctfloor=(modInAsLength==0 && spbn!=1)?(floorInAsLength-1):(floorInAsLength);

            for(int i=0;i<dest-2;++i)
            {
                std::string texture;
                if(i+1>=spNum.size()/2-1-correctfloor && i+1<=spNum.size()/2+correctfloor)
                {
                    texture="GaAs";
                }else{
                    texture="AlAs";
                }
                spStr.push_back(texture);
            }//stacking


            pushVec.resize(dest);
            for(int i=0;i<dest;++i)
            {
                int tmp=0;
                for(int j=0;j<i;++j)
                {
                    tmp+=spNum[j];
                }
                pushVec[i]=tmp;
            }

            std::vector<int> JointValue;
            JointValue.resize(0);
            spJoint.resize(0);
            std::vector<std::string> spStr_tmp;
	    vecEnergyCutoff.resize(0);
            spStr_tmp.resize(0);
            for(int i=0;i<dest-2;++i)
            {
                int tmpint=spNum[i+1];
                spJoint.push_back(std::complex<int>{accumulate(JointValue.begin(),JointValue.end(),0),tmpint});
                JointValue.push_back(tmpint);
                std::string tmpspStr=spStr[i];
                for(int s=0;s<tmpint*uvector.size();++s)
                {
                    spStr_tmp.push_back(tmpspStr);
		    vecEnergyCutoff.push_back(std::pair<double,double>(leftCutoff,rightCutoff));
                }
            }
            JointValue.resize(0);

            BlockStr=spStr;
            spStr=spStr_tmp;
            FifthLine.resize(0);
            if(paramd==1)
            {
                for(int l=0;l<BlockStr.size();++l)
                {
                    FifthLine.push_back(0);
                }
            }
            else if(paramd==pow(2,DividedFactor))
            {
                for(int l=0;l<BlockStr.size();++l)
                {
                    FifthLine.push_back(l);
                }
            }
            else
            {
                int tmp=0;
                for(int l=1;l<spNum.size();++l)
                {
                    int criteria=(pow(2,DividedFactor)/paramd);
                    if(spNum[l]==criteria)
                    {
                        tmp=0;
                        FifthLine.push_back(l);
                    }
                    else
                    {
                        tmp+=spNum[l];
                        if(spNum[l]<criteria && tmp<=criteria)
                        {
                            FifthLine.push_back(spNum.size()+1);
                        }
                        else
                        {
                            FifthLine.push_back(spNum.size()+2);
                        }
                    }
                }
            }

            TypeCollection.resize(0);
            std::vector<int> tmpfoo2;
            for(int i=0;i<dest-2;++i)
            {
                int tmpint=FifthLine[i];
                tmpfoo2.push_back(FifthLine[i]);
                int count=0;
                for(int j=tmpfoo2.size()-1;j>=0;--j)
                {  
                    if(tmpfoo2[i]!=tmpfoo2[j]) break;
                    else count++;
                }
                if(TypeCollection.size()==0 | TypeCollection.back().real()!=tmpint) TypeCollection.push_back(std::complex<int>{tmpint,--count});
                if(TypeCollection.back().real()==tmpint) TypeCollection.back()=std::complex<int>{tmpint,count};
            }
            tmpfoo2.resize(0);

            for(int i=0;i<TypeCollection.size();++i)
            {
                int pos=0;
                for(int j=i-1;j>=0;--j)
                {
                    pos+=TypeCollection[j].imag()?TypeCollection[j].imag():1;
                }            
                TypeCollection[i]=std::complex<int>{pos,TypeCollection[i].imag()};
            }

        }

        void Cell::GenSuperCellRequared()
        {
            for(auto &pcell : supercell)
            {
                int ptag=pcell->Tag2;
                int posip=0;
                std::for_each(supercell.begin(),supercell.end(),
                        [&posip,ptag]
                        (auto &scell)
                        {if(scell->Tag2<ptag){posip+=scell->AtomNum;}}
                        );
                pcell->blpo=posip*orbNum;
                pcell->AtomO=pcell->AtomNum*orbNum;
            }
        }
/*
	void Cell::OldSelectedMasterGen(int onoff)
       {
           
            auto Is_Neighbor=[dt=dis](Eigen::Vector3d &vec1,Eigen::Vector3d &vec2)
            {
                Eigen::Vector3d di=vec1-vec2;
                double distance=di.squaredNorm();
                return (fabs(distance-dt)<1e-10 || distance==0)?true:false;
            };
            if(onoff & LoMEType1) ListofMatrixElements_Type1.clear();
            if(onoff & LoMEType2) ListofMatrixElements_Type2.clear();
            if(onoff & LoMEType3)
            {
                ListofMatrixElements_Type3.clear();
                for(auto &pcell : supercell)
                {
                    pcell->vblpo.resize(0);
                    pcell->vAtomO.resize(0);
                }
            }
            if(onoff & LoMC) ListofMatrixCoupling.clear();

            int usize=uvector.size();
            //fspNum=spNum={1,3,1,1,3};
            for(auto &pcell : supercell)
            {
                int p=pcell->Tag2;

                int ptag=pcell->tag;
                int posip=0;//posi(pcell->tag);
                std::for_each(supercell.begin(),supercell.end(),
                        [&posip,ptag]
                        (auto &scell)
                        {if(scell->tag<ptag){posip+=scell->AtomNum;}}
                        );

                int posip2=0;
                std::for_each(supercell.begin(),supercell.end(),
                        [&posip2,p]
                        (auto &scell)
                        {if(scell->Tag2<p){posip2+=scell->AtomNum;}}
                        );

                if(onoff & LoMEType3)
                {
                    pcell->vblpo.push_back(posip*orbNum);
                    pcell->vAtomO.push_back(pcell->AtomNum*orbNum);
                }
                for(auto &ppcell : supercell)
                {

                    int pp=ppcell->Tag2;

                    int pptag=ppcell->tag;
                    int posipp=0;
                    std::for_each(supercell.begin(),supercell.end(),
                            [&posipp,pptag]
                            (auto &scell)
                            {if(scell->tag<pptag){posipp+=scell->AtomNum;}}
                            );

                    int posipp2=0;
                    std::for_each(supercell.begin(),supercell.end(),
                        [&posipp2,pp]
                        (auto &scell)
                        {if(scell->Tag2<pp){posipp2+=scell->AtomNum;}}
                        );

                    for(int k=0;k<pcell->AtomNum;++k)
                    {
                        for(int j=0;j<ppcell->AtomNum;++j)
                        {
                            //modify
                            if(Is_Neighbor(ppcell->Unit[j],pcell->Unit[k]))
                            {

                                std::string pstexture=pcell->UnitTexture[k/UnitCellNumber];
                                std::string ppstexture=ppcell->UnitTexture[j/UnitCellNumber];
                                if(onoff & LoMC)
                                {
                                    if(
                                            (k+posip>j+posipp)
                                            && (p!=pp)
                                      )
                                    {
                                        ListofMatrixCoupling.push_back(std::make_unique<GenTheList>(std::complex<int>{(k+posip)*orbNum,(j+posipp)*orbNum},pcell->Tag3,ppcell->Tag3,p,pp,posipp2*orbNum,posip2*orbNum));
                                    }
                                }
                                if(onoff & LoMEType2)
                                {//gmarkonoff

                                    std::string pstexture_l=supercell[posip2/UnitCellNumber]->UnitTexture[0];
                                    std::string pstexture_r=supercell[posip2/UnitCellNumber+fspNum[p]*usize-1]->UnitTexture[0];
                                    std::string ppstexture_l=supercell[posipp2/UnitCellNumber]->UnitTexture[0];
                                    std::string ppstexture_r=supercell[posipp2/UnitCellNumber+fspNum[pp]*usize-1]->UnitTexture[0]; 

                                    ListofMatrixElements_Type2.push_back(std::make_unique<Neighbor>(std::complex<int>{k+posip,j+posipp},ppcell->Unit[j]-pcell->Unit[k],ppstexture,pstexture,ppcell->Unit[0]-pcell->Unit[0],0,0));//half
                                    if(fspNum.size()==1)
                                    {

                                    }
                                    else if((pp-p)==1)
                                    {

                                            Eigen::Vector3d vec=ppcell->Unit[j]-pcell->Unit[k];
                                            double DirectFlag=vec.dot(a1);
                                            if(DirectFlag>0)
                                            {
                                                ListofMatrixElements_Type2.push_back(std::make_unique<Neighbor>(std::complex<int>{k+posip,j+posipp-fspNum[p]*usize*UnitCellNumber},ppcell->Unit[j]-pcell->Unit[k],pstexture_l,pstexture_r,ppcell->Unit[0]-pcell->Unit[0],0,0));
                                            }else{
                                                ListofMatrixElements_Type2.push_back(std::make_unique<Neighbor>(std::complex<int>{k+posip+fspNum[p]*usize*UnitCellNumber,j+posipp},ppcell->Unit[j]-pcell->Unit[k],ppstexture_r,ppstexture_l,ppcell->Unit[0]-pcell->Unit[0],0,0));
                                            }

                                        }
                                        else if((pp-p)==-1)
                                        {
                                            Eigen::Vector3d vec=ppcell->Unit[j]-pcell->Unit[k];
                                            double DirectFlag=vec.dot(a1);
                                            if(DirectFlag>0)
                                            {
                                                ListofMatrixElements_Type2.push_back(std::make_unique<Neighbor>(std::complex<int>{k+posip,j+posipp+fspNum[pp]*usize*UnitCellNumber},ppcell->Unit[j]-pcell->Unit[k],pstexture_l,pstexture_r,ppcell->Unit[0]-pcell->Unit[0],0,0));
                                            }else{
                                                ListofMatrixElements_Type2.push_back(std::make_unique<Neighbor>(std::complex<int>{k+posip-fspNum[pp]*usize*UnitCellNumber,j+posipp},ppcell->Unit[j]-pcell->Unit[k],ppstexture_r,ppstexture_l,ppcell->Unit[0]-pcell->Unit[0],0,0));
                                            }

                                        }

                                }
                                if(onoff & LoMEType1)
                                {
                                    ListofMatrixElements_Type1.push_back(std::make_unique<Neighbor>(std::complex<int>{k+posip,j+posipp},ppcell->Unit[j]-pcell->Unit[k],ppstexture,pstexture,ppcell->Unit[0]-pcell->Unit[0],0,0));//half
                                }
                                if(onoff & LoMEType3)//gonoffmark
                                {

                                    if(pp!=p)
                                    {
                                        ListofMatrixElements_Type3.push_back(std::make_unique<Neighbor>(std::complex<int>{k+posip,j+posipp},ppcell->Unit[j]-pcell->Unit[k],ppstexture,pstexture,ppcell->Unit[0]-pcell->Unit[0],0,0));//half
                                    }


                                   ListofMatrixElements_Type3.push_back(std::make_unique<Neighbor>(std::complex<int>{k+posip,j+posip},ppcell->Unit[j]-pcell->Unit[k],pstexture,pstexture,ppcell->Unit[0]-pcell->Unit[0],(ptag-posip2/UnitCellNumber),p));//half

                                }

                            }
                        }
                    }
                }

                for(auto &ppcell : supercellMargin)
                {
                    int pp=ppcell->Tag2;
                    int pptag=ppcell->tag;
                    int posipp=0;//posi(ppcell->tag);
                    std::for_each(supercell.begin(),supercell.end(),
                            [&posipp,pptag]
                            (auto &scell)
                            {if(scell->tag<pptag){posipp+=scell->AtomNum;}}
                            );
                    
                    int posipp2=0;
                    std::for_each(supercell.begin(),supercell.end(),
                        [&posipp2,pp]
                        (auto &scell)
                        {if(scell->Tag2<pp){posipp2+=scell->AtomNum;}}
                        );

                    for(int k=0;k<pcell->AtomNum;++k)
                    {
                        for(int j=0;j<ppcell->AtomNum;++j)
                        {
                            if(Is_Neighbor(ppcell->Unit[j],pcell->Unit[k]))
                            {
                                std::string pstexture=pcell->UnitTexture[k/UnitCellNumber];
                                std::string ppstexture=supercell[pptag]->UnitTexture[j/UnitCellNumber];

                                
                                if(onoff & LoMC)
                                {
                                    if(
                                            (k+posip>j+posipp)
                                            && (p!=pp)
                                      )
                                    {
                                        ListofMatrixCoupling.push_back(std::make_unique<GenTheList>(std::complex<int>{(k+posip)*orbNum,(j+posipp)*orbNum},pcell->Tag3,ppcell->Tag3,p,pp,posipp2*orbNum,posip2*orbNum));
                                    }
                                }
                                if(onoff & LoMEType2)//gmarkonoffg
                                {
                                    std::string pstexture_l=supercell[posip2/UnitCellNumber]->UnitTexture[0];
                                    std::string pstexture_r=supercell[posip2/UnitCellNumber+fspNum[p]*usize-1]->UnitTexture[0];
                                    std::string ppstexture_l=supercell[posipp2/UnitCellNumber]->UnitTexture[0];
                                    std::string ppstexture_r=supercell[posipp2/UnitCellNumber+fspNum[pp]*usize-1]->UnitTexture[0]; 

                                    ListofMatrixElements_Type2.push_back(std::make_unique<Neighbor>(std::complex<int>{k+posip,j+posipp},ppcell->Unit[j]-pcell->Unit[k],ppstexture,pstexture,ppcell->Unit[0]-pcell->Unit[0],0,0));//half
                                    if(fspNum.size()==1)
                                    {

                                    }
                                    else if((pp-p)==1)
                                    {

                                           Eigen::Vector3d vec=ppcell->Unit[j]-pcell->Unit[k];
                                            double DirectFlag=vec.dot(a1);
                                            if(DirectFlag>0)
                                            {
                                                ListofMatrixElements_Type2.push_back(std::make_unique<Neighbor>(std::complex<int>{k+posip,j+posipp-fspNum[p]*usize*UnitCellNumber},ppcell->Unit[j]-pcell->Unit[k],pstexture_l,pstexture_r,ppcell->Unit[0]-pcell->Unit[0],0,0));
                                            }else{
                                                ListofMatrixElements_Type2.push_back(std::make_unique<Neighbor>(std::complex<int>{k+posip+fspNum[p]*usize*UnitCellNumber,j+posipp},ppcell->Unit[j]-pcell->Unit[k],ppstexture_r,ppstexture_l,ppcell->Unit[0]-pcell->Unit[0],0,0));
                                            }

                                        }
                                        else if((pp-p)==-1)
                                        {

                                            Eigen::Vector3d vec=ppcell->Unit[j]-pcell->Unit[k];
                                            double DirectFlag=vec.dot(a1);
                                            if(DirectFlag>0)
                                            {
                                                ListofMatrixElements_Type2.push_back(std::make_unique<Neighbor>(std::complex<int>{k+posip,j+posipp+fspNum[pp]*usize*UnitCellNumber},ppcell->Unit[j]-pcell->Unit[k],pstexture_l,pstexture_r,ppcell->Unit[0]-pcell->Unit[0],0,0));
                                            }else{
                                                ListofMatrixElements_Type2.push_back(std::make_unique<Neighbor>(std::complex<int>{k+posip-fspNum[pp]*usize*UnitCellNumber,j+posipp},ppcell->Unit[j]-pcell->Unit[k],ppstexture_r,ppstexture_l,ppcell->Unit[0]-pcell->Unit[0],0,0));
                                            }

                                        }
                                        else if((pp-p)==(fspNum.size()-1) &&spNum.size()-2 >2)
                                        {
 
                                           Eigen::Vector3d vec=ppcell->Unit[j]-pcell->Unit[k];
                                            double DirectFlag=vec.dot(a1);
                                            if(DirectFlag<0)
                                            {
                                                ListofMatrixElements_Type2.push_back(std::make_unique<Neighbor>(std::complex<int>{k+posip+posipp2,j+posipp},ppcell->Unit[j]-pcell->Unit[k],ppstexture_r,ppstexture_l,ppcell->Unit[0]-pcell->Unit[0],0,0));
                                            }

                                       }
                                        else if((pp-p)==-(fspNum.size()-1) &&spNum.size()-2 >2)
                                        {

                                            Eigen::Vector3d vec=ppcell->Unit[j]-pcell->Unit[k];
                                            double DirectFlag=vec.dot(a1);
                                            if(DirectFlag>0)
                                            {
                                                ListofMatrixElements_Type2.push_back(std::make_unique<Neighbor>(std::complex<int>{k+posip,j+posipp+posip2},ppcell->Unit[j]-pcell->Unit[k],pstexture_l,pstexture_r,ppcell->Unit[0]-pcell->Unit[0],0,0));
                                            }
                                        }

                                }
                                if(onoff & LoMEType1)
                                {

                                    ListofMatrixElements_Type1.push_back(std::make_unique<Neighbor>(std::complex<int>{k+posip,j+posipp},ppcell->Unit[j]-pcell->Unit[k],ppstexture,pstexture,ppcell->Unit[0]-pcell->Unit[0],0,0));
                                }
                                if(onoff & LoMEType3)//gonoffmark
                                {

                                    if(pp!=p)
                                    {
                                        ListofMatrixElements_Type3.push_back(std::make_unique<Neighbor>(std::complex<int>{k+posip,j+posipp},ppcell->Unit[j]-pcell->Unit[k],ppstexture,pstexture,ppcell->Unit[0]-pcell->Unit[0],0,0));
                                    }
                                    ListofMatrixElements_Type3.push_back(std::make_unique<Neighbor>(std::complex<int>{k+posip,j+posip},ppcell->Unit[j]-pcell->Unit[k],pstexture,pstexture,ppcell->Unit[0]-pcell->Unit[0],(ptag-posip2/UnitCellNumber),p));//half

                                }
                            }
                        }
                    }
                }

            }

            fspNum.resize(0);
        }
*/
        void Cell::SelectedMasterGen(int onoff)
        {

            int id,numprocs;
            MPI_Comm_rank(MPI_COMM_WORLD,&id);

            auto Is_Neighbor=[dt=dis](Eigen::Vector3d &vec1,Eigen::Vector3d &vec2)
            {
                Eigen::Vector3d di=vec1-vec2;
                double distance=di.squaredNorm();
                return (fabs(distance-dt)<1e-10 || distance==0)?true:false;
            };
            if(onoff & LoMEType1) 
            {
                ListofMatrixElements_Type1.clear();
                ListofMatrixElements_MixType1.clear();
            }
            if(onoff & LoMEType2) ListofMatrixElements_Type2.clear();
            if(onoff & LoMEType3)
            {
                ListofMatrixElements_Type3.clear();
                for(auto &pcell : supercell)
                {
                    pcell->vblpo2.resize(0);
                    pcell->vAtomO2.resize(0);
                    pcell->vblpo.resize(0);
                    pcell->vAtomO.resize(0);
                    pcell->vEnergyCutoff.resize(0);
                }
            }
            if(onoff & LoMC) ListofMatrixCoupling.clear();

            int usize=uvector.size();
            //fspNum=spNum={1,3,1,1,3};
            int NeighborFlag=0;
            for(auto &pcell : supercell)
            {
                int p=pcell->Tag2;

                int ptag=pcell->tag;
                int posip=0;//posi(pcell->tag);
                std::for_each(supercell.begin(),supercell.end(),
                        [&posip,ptag]
                        (auto &scell)
                        {if(scell->tag<ptag){posip+=scell->AtomNum;}}
                        );

                int posip2=0;
                std::for_each(supercell.begin(),supercell.end(),
                        [&posip2,p]
                        (auto &scell)
                        {if(scell->Tag2<p){posip2+=scell->AtomNum;}}
                        );

                if(onoff & LoMEType3)
                {
                    pcell->vblpo.push_back(posip*orbNum);
                    pcell->vAtomO.push_back(pcell->AtomNum*orbNum);
                }

                int pdis=std::distance(supercell.begin(),std::find(supercell.begin(),supercell.end(),pcell));
                for(int k=0;k<pcell->AtomNum;++k)
                {
                    NeighborFlag=0;
                    if(NeighborFlag<4)
                    {

                        for(int ppp=-1;ppp<=1;++ppp)
                        {
                            int SubcellNumber=uvector.size();
                            Eigen::Vector3d mshift(0,0,0);
                            int cshift=pdis/SubcellNumber+ppp;
                            if(ppp==0)
                            {
                            }else if(cshift==-1)
                            {
                                cshift=supercell.size()/SubcellNumber-1;
                                mshift=-Unita1;
                            }else if(cshift==supercell.size()/SubcellNumber)
                            {
                                cshift=0;
                                mshift=Unita1;
                            }else{
                            }
                            for(auto &UCMargin : UnitCellMargin)
                            {
                        
                            for(int unitcells=0;unitcells<uvector.size();++unitcells)
                            {

                                auto &ppcell=*(supercell.begin()+cshift*SubcellNumber+unitcells);

                                int pp=ppcell->Tag2;
                                int pptag=ppcell->tag;
                                int posipp=0;//posi(ppcell->tag);
                                std::for_each(supercell.begin(),supercell.end(),
                                        [&posipp,pptag]
                                        (auto &scell)
                                        {if(scell->tag<pptag){posipp+=scell->AtomNum;}}
                                        );

                                int posipp2=0;
                                std::for_each(supercell.begin(),supercell.end(),
                                        [&posipp2,pp]
                                        (auto &scell)
                                        {if(scell->Tag2<pp){posipp2+=scell->AtomNum;}}
                                        );
                                //if((ppcell->type==0 || fabs(pptag-ptag)<=1))
                                if((fabs(pptag/SubcellNumber-ptag/SubcellNumber)==supercell.size()/SubcellNumber-1 || fabs(pptag/SubcellNumber-ptag/SubcellNumber)<=1))
                                {
                                    for(int j=0;j<ppcell->AtomNum;++j)
                                    {
                                        Eigen::Vector3d Mshift=mshift+ppcell->Unit[j]+UCMargin;
                                        if(Is_Neighbor(Mshift,pcell->Unit[k]))
                                        {
                                            NeighborFlag+=1;
                                            std::string pstexture=pcell->UnitTexture[k/UnitCellNumber];
                                            std::string ppstexture=supercell[pptag]->UnitTexture[j/UnitCellNumber];


                                            if((onoff & LoMC ))
                                            {
                                                if(
                                                        (k+posip>j+posipp)
                                                        && (p!=pp)
                                                  )
                                                {
                                                    /*if(USE_3Dpbc || pcell->type*ppcell->type)*/ ListofMatrixCoupling.push_back(std::make_unique<GenTheList>(std::complex<int>{(k+posip)*orbNum,(j+posipp)*orbNum},pcell->Tag3,ppcell->Tag3,p,pp,posipp2*orbNum,posip2*orbNum));
                                                }
                                            }
                                            if(onoff & LoMEType2)//gmarkonoffg
                                            {
                                                std::string pstexture_l=supercell[posip2/UnitCellNumber]->UnitTexture[0];
                                                std::string pstexture_r=supercell[posip2/UnitCellNumber+fspNum[p]*usize-1]->UnitTexture[0];
                                                std::string ppstexture_l=supercell[posipp2/UnitCellNumber]->UnitTexture[0];
                                                std::string ppstexture_r=supercell[posipp2/UnitCellNumber+fspNum[pp]*usize-1]->UnitTexture[0]; 

                                                ListofMatrixElements_Type2.push_back(std::make_unique<Neighbor>(std::complex<int>{k+posip,j+posipp},Mshift-pcell->Unit[k],ppstexture,pstexture,ppcell->Unit[0]-pcell->Unit[0],0,0));//half
                                                if(fspNum.size()==1)
                                                {

                                                }
                                                else if((pp-p)==1)
                                                {

                                                    Eigen::Vector3d vec=Mshift-pcell->Unit[k];
                                                    double DirectFlag=vec.dot(a1);
                                                    if(DirectFlag>0)
                                                    {
                                                        ListofMatrixElements_Type2.push_back(std::make_unique<Neighbor>(std::complex<int>{k+posip,j+posipp-fspNum[p]*usize*UnitCellNumber},Mshift-pcell->Unit[k],pstexture_l,pstexture_r,ppcell->Unit[0]-pcell->Unit[0],0,0));
                                                    }else{
                                                        ListofMatrixElements_Type2.push_back(std::make_unique<Neighbor>(std::complex<int>{k+posip+fspNum[p]*usize*UnitCellNumber,j+posipp},Mshift-pcell->Unit[k],ppstexture_r,ppstexture_l,ppcell->Unit[0]-pcell->Unit[0],0,0));
                                                    }

                                                }
                                                else if((pp-p)==-1)
                                                {

                                                    Eigen::Vector3d vec=Mshift-pcell->Unit[k];
                                                    double DirectFlag=vec.dot(a1);
                                                    if(DirectFlag>0)
                                                    {
                                                        ListofMatrixElements_Type2.push_back(std::make_unique<Neighbor>(std::complex<int>{k+posip,j+posipp+fspNum[pp]*usize*UnitCellNumber},Mshift-pcell->Unit[k],pstexture_l,pstexture_r,ppcell->Unit[0]-pcell->Unit[0],0,0));
                                                    }else{
                                                        ListofMatrixElements_Type2.push_back(std::make_unique<Neighbor>(std::complex<int>{k+posip-fspNum[pp]*usize*UnitCellNumber,j+posipp},Mshift-pcell->Unit[k],ppstexture_r,ppstexture_l,ppcell->Unit[0]-pcell->Unit[0],0,0));
                                                    }

                                                }
                                                else if((pp-p)==(fspNum.size()-1) &&spNum.size()-2 >2)
                                                {

                                                    Eigen::Vector3d vec=Mshift-pcell->Unit[k];
                                                    double DirectFlag=vec.dot(a1);
                                                    if(DirectFlag<0)
                                                    {
                                                        ListofMatrixElements_Type2.push_back(std::make_unique<Neighbor>(std::complex<int>{k+posip+posipp2,j+posipp},Mshift-pcell->Unit[k],ppstexture_r,ppstexture_l,ppcell->Unit[0]-pcell->Unit[0],0,0));
                                                    }

                                                }
                                                else if((pp-p)==-(fspNum.size()-1) &&spNum.size()-2 >2)
                                                {

                                                    Eigen::Vector3d vec=Mshift-pcell->Unit[k];
                                                    double DirectFlag=vec.dot(a1);
                                                    if(DirectFlag>0)
                                                    {
                                                        ListofMatrixElements_Type2.push_back(std::make_unique<Neighbor>(std::complex<int>{k+posip,j+posipp+posip2},Mshift-pcell->Unit[k],pstexture_l,pstexture_r,ppcell->Unit[0]-pcell->Unit[0],0,0));
                                                    }
                                                }

                                            }
                                            if((onoff & LoMEType1) )
                                            {

                                                /*if(USE_3Dpbc || pcell->type*ppcell->type)*/ ListofMatrixElements_Type1.push_back(std::make_unique<Neighbor>(std::complex<int>{k+posip,j+posipp},Mshift-pcell->Unit[k],ppstexture,pstexture,ppcell->Unit[0]-pcell->Unit[0],0,p));

                                                if(ppstexture!=pstexture)
                                                {
                                                    if(k+posip>j+posipp)
                                                    {
                                                        /*if(USE_3Dpbc || pcell->type*ppcell->type)*/ ListofMatrixElements_MixType1.push_back(std::make_unique<Neighbor>(std::complex<int>{j+posipp,j+posipp},Mshift-pcell->Unit[k],ppstexture,pstexture,ppcell->Unit[0]-pcell->Unit[0],1,p));
                                                    }
                                                }

                                            }
                                            if(onoff & LoMEType3)//gonoffmark
                                            {

                                                if(pp!=p)
                                                {
                                                    ListofMatrixElements_Type3.push_back(std::make_unique<Neighbor>(std::complex<int>{k+posip,j+posipp},Mshift-pcell->Unit[k],ppstexture,pstexture,ppcell->Unit[0]-pcell->Unit[0],0,0));
                                                }
                                                ListofMatrixElements_Type3.push_back(std::make_unique<Neighbor>(std::complex<int>{k+posip,j+posip},Mshift-pcell->Unit[k],pstexture,pstexture,ppcell->Unit[0]-pcell->Unit[0],(ptag-posip2/UnitCellNumber),p));//half

                                            }
                                        }
                                    }
                                }
                            }
                            }
                        }

                    }else{
                        break;
                    }
                }

            }

            fspNum.resize(0);
        }
       void Cell::JointSupercell2(int side2)
       {
           //supercell[side2]->vblpo2.push_back(supercell[side2]->AtomNum/(uvector.size()*UnitCellNumber));
           //supercell[side2]->vAtomO2.push_back(supercell[side2+1]->AtomNum/(uvector.size()*UnitCellNumber));

	   if((supercell[side2]->UnitTexture.back()!=supercell[side2+1]->UnitTexture.front()) || supercell[side2]->vblpo2.size()==0)
	   {
	       supercell[side2]->vblpo2.push_back(supercell[side2]->AtomNum/(uvector.size()*UnitCellNumber));
	       supercell[side2]->vAtomO2.push_back(supercell[side2+1]->AtomNum/(uvector.size()*UnitCellNumber));
	       supercell[side2]->vEnergyCutoff.push_back(supercell[side2+1]->EnergyCutoff);
	   }else{
	       supercell[side2]->vAtomO2.back()+=supercell[side2+1]->AtomNum/(uvector.size()*UnitCellNumber);
	       //supercell[side2]->vEnergyCutoff.back().first+=supercell[side2+1]->EnergyCutoff.first;
	       //supercell[side2]->vEnergyCutoff.back().second+=supercell[side2+1]->EnergyCutoff.second;
	       supercell[side2]->vEnergyCutoff.push_back(supercell[side2+1]->EnergyCutoff);
	   }

           supercell[side2]->UnitTexture.insert(supercell[side2]->UnitTexture.end(),supercell[side2+1]->UnitTexture.begin(),supercell[side2+1]->UnitTexture.end());
           supercell[side2]->Unit.insert(supercell[side2]->Unit.end(),supercell[side2+1]->Unit.begin(),supercell[side2+1]->Unit.end());

           supercell[side2]->AtomNum+=supercell[side2+1]->AtomNum;
           supercell[side2]->vblpo.insert(supercell[side2]->vblpo.end(),supercell[side2+1]->vblpo.begin(),supercell[side2+1]->vblpo.end());
           supercell[side2]->vAtomO.insert(supercell[side2]->vAtomO.end(),supercell[side2+1]->vAtomO.begin(),supercell[side2+1]->vAtomO.end());

           //supercell[side2]->vblpo2.push_back(supercell[side2+1]->blpo);
           //supercell[side2]->vAtomO2.push_back(supercell[side2+1]->AtomO);
           
           //supercell[side2]->vEnergyCutoff.push_back(supercell[side2+1]->EnergyCutoff);
           
           //supercell[side2]->vblpo2.insert(supercell[side2]->vblpo2.end(),supercell[side2+1]->vblpo.begin(),supercell[side2+1]->vblpo.end());
           //supercell[side2]->vAtomO2.insert(supercell[side2]->vAtomO2.end(),supercell[side2+1]->vAtomO.begin(),supercell[side2+1]->vAtomO.end());

           supercell.erase(supercell.begin()+side2+1);
       }

        void Cell::JointSupercell(int side2)
        {
            supercell[side2]->UnitTexture.insert(supercell[side2]->UnitTexture.end(),supercell[side2+1]->UnitTexture.begin(),supercell[side2+1]->UnitTexture.end());
            supercell[side2]->Unit.insert(supercell[side2]->Unit.end(),supercell[side2+1]->Unit.begin(),supercell[side2+1]->Unit.end());
            supercell[side2]->AtomNum+=supercell[side2+1]->AtomNum;
            supercell[side2]->vblpo.insert(supercell[side2]->vblpo.end(),supercell[side2+1]->vblpo.begin(),supercell[side2+1]->vblpo.end());
            supercell[side2]->vAtomO.insert(supercell[side2]->vAtomO.end(),supercell[side2+1]->vAtomO.begin(),supercell[side2+1]->vAtomO.end());
            supercell.erase(supercell.begin()+side2+1);
        }

        void Cell::AddBoundaryConditions(int numSize)
        {

            supercellMargin.clear();
            supercellMargin.resize(0);

	    UnitCellMargin.resize(0);
	    UnitCellMargin.push_back(Eigen::Vector3d(0,0,0));
	    UnitCellMargin.push_back(Unita2);
	    UnitCellMargin.push_back(-Unita2);

if(USE_3Dpbc==1)
{
	    UnitCellMargin.push_back(Unita3);
	    UnitCellMargin.push_back(-Unita3);
	    UnitCellMargin.push_back(Unita2+Unita3);
	    UnitCellMargin.push_back(-Unita2+Unita3);
	    UnitCellMargin.push_back(Unita2-Unita3);
	    UnitCellMargin.push_back(-Unita2-Unita3);
}

/*
*/

/*
            for(auto &cell : supercell)
            {
                supercellMargin.push_back(std::make_unique<Geometry>(Unita2,cell->Unit,cell->tag,cell->Tag3,cell->Tag2,1)); 
                supercellMargin.push_back(std::make_unique<Geometry>(-Unita2,cell->Unit,cell->tag,cell->Tag3,cell->Tag2,1)); 
                supercellMargin.push_back(std::make_unique<Geometry>(Unita3,cell->Unit,cell->tag,cell->Tag3,cell->Tag2,1)); 
                supercellMargin.push_back(std::make_unique<Geometry>(-Unita3,cell->Unit,cell->tag,cell->Tag3,cell->Tag2,1)); 

                supercellMargin.push_back(std::make_unique<Geometry>(Unita2+Unita3,cell->Unit,cell->tag,cell->Tag3,cell->Tag2,1)); 
                supercellMargin.push_back(std::make_unique<Geometry>(-Unita2+Unita3,cell->Unit,cell->tag,cell->Tag3,cell->Tag2,1)); 
                supercellMargin.push_back(std::make_unique<Geometry>(Unita2-Unita3,cell->Unit,cell->tag,cell->Tag3,cell->Tag2,1)); 
                supercellMargin.push_back(std::make_unique<Geometry>(-Unita2-Unita3,cell->Unit,cell->tag,cell->Tag3,cell->Tag2,1)); 

            }

            for(int num=0;num<numSize;num++)
            {
                supercellMargin.push_back(std::make_unique<Geometry>(Unita1,supercell[num]->Unit,supercell[num]->tag,supercell[num]->Tag3,supercell[num]->Tag2,0)); 
                supercellMargin.push_back(std::make_unique<Geometry>(Unita1+Unita2,supercell[num]->Unit,supercell[num]->tag,supercell[num]->Tag3,supercell[num]->Tag2,0)); 
                supercellMargin.push_back(std::make_unique<Geometry>(Unita1-Unita2,supercell[num]->Unit,supercell[num]->tag,supercell[num]->Tag3,supercell[num]->Tag2,0)); 
                supercellMargin.push_back(std::make_unique<Geometry>(Unita1+Unita3,supercell[num]->Unit,supercell[num]->tag,supercell[num]->Tag3,supercell[num]->Tag2,0)); 
                supercellMargin.push_back(std::make_unique<Geometry>(Unita1-Unita3,supercell[num]->Unit,supercell[num]->tag,supercell[num]->Tag3,supercell[num]->Tag2,0)); 

                supercellMargin.push_back(std::make_unique<Geometry>(Unita1+Unita2+Unita3,supercell[num]->Unit,supercell[num]->tag,supercell[num]->Tag3,supercell[num]->Tag2,0)); 
                supercellMargin.push_back(std::make_unique<Geometry>(Unita1+Unita2-Unita3,supercell[num]->Unit,supercell[num]->tag,supercell[num]->Tag3,supercell[num]->Tag2,0)); 
                supercellMargin.push_back(std::make_unique<Geometry>(Unita1-Unita2+Unita3,supercell[num]->Unit,supercell[num]->tag,supercell[num]->Tag3,supercell[num]->Tag2,0)); 
                supercellMargin.push_back(std::make_unique<Geometry>(Unita1-Unita2-Unita3,supercell[num]->Unit,supercell[num]->tag,supercell[num]->Tag3,supercell[num]->Tag2,0)); 
            }
            ////
            int backIt=supercell.size()-1;
            for(int num=0;num<numSize;++num)
            {
                supercellMargin.push_back(std::make_unique<Geometry>(-Unita1,supercell[backIt-num]->Unit,supercell[backIt-num]->tag,supercell[backIt-num]->Tag3,supercell[backIt-num]->Tag2,0)); 
                supercellMargin.push_back(std::make_unique<Geometry>(-Unita1+Unita2,supercell[backIt-num]->Unit,supercell[backIt-num]->tag,supercell[backIt-num]->Tag3,supercell[backIt-num]->Tag2,0)); 
                supercellMargin.push_back(std::make_unique<Geometry>(-Unita1-Unita2,supercell[backIt-num]->Unit,supercell[backIt-num]->tag,supercell[backIt-num]->Tag3,supercell[backIt-num]->Tag2,0)); 
                supercellMargin.push_back(std::make_unique<Geometry>(-Unita1+Unita3,supercell[backIt-num]->Unit,supercell[backIt-num]->tag,supercell[backIt-num]->Tag3,supercell[backIt-num]->Tag2,0)); 
                supercellMargin.push_back(std::make_unique<Geometry>(-Unita1-Unita3,supercell[backIt-num]->Unit,supercell[backIt-num]->tag,supercell[backIt-num]->Tag3,supercell[backIt-num]->Tag2,0)); 

                supercellMargin.push_back(std::make_unique<Geometry>(-Unita1+Unita2+Unita3,supercell[backIt-num]->Unit,supercell[backIt-num]->tag,supercell[backIt-num]->Tag3,supercell[backIt-num]->Tag2,0)); 
                supercellMargin.push_back(std::make_unique<Geometry>(-Unita1+Unita2-Unita3,supercell[backIt-num]->Unit,supercell[backIt-num]->tag,supercell[backIt-num]->Tag3,supercell[backIt-num]->Tag2,0)); 
                supercellMargin.push_back(std::make_unique<Geometry>(-Unita1-Unita2+Unita3,supercell[backIt-num]->Unit,supercell[backIt-num]->tag,supercell[backIt-num]->Tag3,supercell[backIt-num]->Tag2,0)); 
                supercellMargin.push_back(std::make_unique<Geometry>(-Unita1-Unita2-Unita3,supercell[backIt-num]->Unit,supercell[backIt-num]->tag,supercell[backIt-num]->Tag3,supercell[backIt-num]->Tag2,0)); 
            }
*/

        }

        int Cell::GetMatrixSize()
        {
            int MatrixSize=0;
            for(auto &cell:supercell)
            {
                MatrixSize+=cell->AtomNum*orbNum;
            }
            return MatrixSize; 
        }
/*
        void Cell::StackingStructure3D()
        {

            std::vector<std::string> cellString;
            std::vector<std::string> ReducedcellString;
            cellString.resize(0);
            ReducedcellString.resize(0);
            supercell.clear();
            int count=0;

            fspNum.resize(0);
            std::vector<std::string> tmpStrfoo;
            tmpStrfoo.resize(0);
            int ct2=0;
            for(int i=0;i<TypeCollection.size();++i)
            {
                if(((ct2)==TypeCollection[i].real())&&TypeCollection[i].imag()!=0)
                {
                    fspNum.push_back(accumulate(spNum.begin()+TypeCollection[i].real()+1,spNum.begin()+TypeCollection[i].real()+TypeCollection[i].imag()+1,0));
                    for(int j=0;j<(TypeCollection[i].imag());++j)
                    {
                        std::string strfortmp("");
                        for(int istr=0;istr<TypeCollection[i].imag();++istr)
                        {
                            int key=TypeCollection[i].real()+istr;
                            for(int r=0;r<spNum[key+1];++r)
                            {
                                strfortmp+=BlockStr[key]+"-";
                            }
                        }
                        tmpStrfoo.push_back(strfortmp);
                        ct2+=1;
                    }
                }else{
                    std::string strfortmp("");
                    for(int r=0;r<spNum[ct2+1];++r)
                    {
                        strfortmp+=BlockStr[ct2]+"-";
                    }
                    tmpStrfoo.push_back(strfortmp);
                    fspNum.push_back(spNum[ct2+1]);
                    ct2+=1;
                }
            }

            int Tag2=0;
            int TypeStatus=FifthLine[0];
            for(int i=0;i<spNum.size()-2;++i)
            {
                std::string ctex=tmpStrfoo[i];

                if(std::find(cellString.begin(),cellString.end(),ctex)==cellString.end())
                {
                    ReducedcellString.push_back(ctex);
                }
                cellString.push_back(ctex);
                int Tag3=std::distance(ReducedcellString.begin(),std::find(ReducedcellString.begin(),ReducedcellString.end(),ctex));
                if(TypeStatus!=FifthLine[i]){Tag2+=1;TypeStatus=FifthLine[i];}

                for(int in=0;in<spNum[i+1];++in)
                {
                    for(int ii=0;ii<uvector.size();++ii)
                    {
                        int Tag1=uvector.size()*count+ii;
                        supercell.push_back(std::make_unique<Geometry>(a1*count+uvector[ii],UnitCell,Tag1,Tag3,Tag2,spStr[Tag1],std::pair<double,double>(16.,64.))); 
                    }
                    count+=1;
                }
            }
            cellString.resize(0);
            ReducedcellString.resize(0);
            spStr.resize(0);
            FifthLine.resize(0);
            tmpStrfoo.resize(0);
            BlockStr.resize(0);
            AddBoundaryConditions(uvector.size());
        
        }
*/
        void Cell::StackingStructure()
        {

            std::vector<std::string> cellString;
            std::vector<std::string> ReducedcellString;
            cellString.resize(0);
            ReducedcellString.resize(0);
            supercell.clear();
            int count=0;

            fspNum.resize(0);
            std::vector<std::string> tmpStrfoo;
            tmpStrfoo.resize(0);
            int ct2=0;
            for(int i=0;i<TypeCollection.size();++i)
            {
                if(((ct2)==TypeCollection[i].real())&&TypeCollection[i].imag()!=0)
                {
                    fspNum.push_back(accumulate(spNum.begin()+TypeCollection[i].real()+1,spNum.begin()+TypeCollection[i].real()+TypeCollection[i].imag()+1,0));
                    for(int j=0;j<(TypeCollection[i].imag());++j)
                    {
                        std::string strfortmp("");
                        for(int istr=0;istr<TypeCollection[i].imag();++istr)
                        {
                            int key=TypeCollection[i].real()+istr;
                            for(int r=0;r<spNum[key+1];++r)
                            {
                                strfortmp+=BlockStr[key]+"-";
                            }
                        }
                        tmpStrfoo.push_back(strfortmp);
                        ct2+=1;
                    }
                }else{
                    std::string strfortmp("");
                    for(int r=0;r<spNum[ct2+1];++r)
                    {
                        strfortmp+=BlockStr[ct2]+"-";
                    }
                    tmpStrfoo.push_back(strfortmp);
                    fspNum.push_back(spNum[ct2+1]);
                    ct2+=1;
                }
            }

            int Tag2=0;
            int TypeStatus=FifthLine[0];
            for(int i=0;i<spNum.size()-2;++i)
            {
                std::string ctex=tmpStrfoo[i];

                if(std::find(cellString.begin(),cellString.end(),ctex)==cellString.end())
                {
                    ReducedcellString.push_back(ctex);
                }
                cellString.push_back(ctex);
                int Tag3=std::distance(ReducedcellString.begin(),std::find(ReducedcellString.begin(),ReducedcellString.end(),ctex));
                if(TypeStatus!=FifthLine[i]){Tag2+=1;TypeStatus=FifthLine[i];}

                for(int in=0;in<spNum[i+1];++in)
                {
                    for(int ii=0;ii<uvector.size();++ii)
                    {
                        int Tag1=uvector.size()*count+ii;
                        if(vecEnergyCutoff.size()==0)
                        {
                        supercell.push_back(std::make_unique<Geometry>(a1*count+uvector[ii],UnitCell,Tag1,Tag3,Tag2,spStr[Tag1],std::pair<double,double>{0,0})); 
                        }else{
                        supercell.push_back(std::make_unique<Geometry>(a1*count+uvector[ii],UnitCell,Tag1,Tag3,Tag2,spStr[Tag1],vecEnergyCutoff[Tag1])); 
                        }
                    }
                    count+=1;
                }
            }
            cellString.resize(0);
            ReducedcellString.resize(0);
            spStr.resize(0);
            //FifthLine.resize(0);
            tmpStrfoo.resize(0);
            BlockStr.resize(0);
            //if(USE_3Dpbc){
                AddBoundaryConditions(uvector.size());
            //}else{
                //AddBoundaryConditions(0);
            //}
            /*
               std::cout<<"Unita1=("<<Unita1(0)<<","<<Unita1(1)<<","<<Unita1(2)<<")"<<std::endl;
               std::cout<<"Unita2=("<<Unita2(0)<<","<<Unita2(1)<<","<<Unita2(2)<<")"<<std::endl;
               std::cout<<"Unita3=("<<Unita3(0)<<","<<Unita3(1)<<","<<Unita3(2)<<")"<<std::endl;
               Eigen::Vector3d rUnita1=RotateAround(Eigen::Vector3d(1,0,0),Unita1,-M_PI*0.25);
               Eigen::Vector3d rUnita2=RotateAround(Eigen::Vector3d(1,0,0),Unita2,-M_PI*0.25)*sqrt(2);
               Eigen::Vector3d rUnita3=RotateAround(Eigen::Vector3d(1,0,0),Unita3,-M_PI*0.25)*sqrt(2);

               std::cout<<"Rnita1=("<<rUnita1(0)<<","<<rUnita1(1)<<","<<rUnita1(2)<<")"<<std::endl;
               std::cout<<"Rnita2=("<<rUnita2(0)<<","<<rUnita2(1)<<","<<rUnita2(2)<<")"<<std::endl;
               std::cout<<"Rnita3=("<<rUnita3(0)<<","<<rUnita3(1)<<","<<rUnita3(2)<<")"<<std::endl;
               */
        }

        void Cell::ParseCell(std::string inpfile)
        {

            ParseBlockElement<double>("CellTranslatedVector",inpfile,0,1,a1(0));
            ParseBlockElement<double>("CellTranslatedVector",inpfile,0,2,a1(1));
            ParseBlockElement<double>("CellTranslatedVector",inpfile,0,3,a1(2));

            ParseBlockElement<double>("CellTranslatedVector",inpfile,1,1,a2(0));
            ParseBlockElement<double>("CellTranslatedVector",inpfile,1,2,a2(1));
            ParseBlockElement<double>("CellTranslatedVector",inpfile,1,3,a2(2));

            ParseBlockElement<double>("CellTranslatedVector",inpfile,2,1,a3(0));
            ParseBlockElement<double>("CellTranslatedVector",inpfile,2,2,a3(1));
            ParseBlockElement<double>("CellTranslatedVector",inpfile,2,3,a3(2));


            uvector.resize(0);
            uvector.push_back(Eigen::Vector3d(0,0,0));
            for(int i=0;i<ParseBlockTitle("CellConstructedVector",inpfile,0);++i)
            {
                Eigen::Vector3d tmp;
                ParseBlockElement<double>("CellConstructedVector",inpfile,i,1,tmp(0));
                ParseBlockElement<double>("CellConstructedVector",inpfile,i,2,tmp(1));
                ParseBlockElement<double>("CellConstructedVector",inpfile,i,3,tmp(2));
                uvector.push_back(tmp);
            }

            ParseBlockElement<double>("CellConstructedVector",inpfile,0,1,u1(0));
            ParseBlockElement<double>("CellConstructedVector",inpfile,0,2,u1(1));
            ParseBlockElement<double>("CellConstructedVector",inpfile,0,3,u1(2));

            ParseBlockElement<double>("CellConstructedVector",inpfile,1,1,u2(0));
            ParseBlockElement<double>("CellConstructedVector",inpfile,1,2,u2(1));
            ParseBlockElement<double>("CellConstructedVector",inpfile,1,3,u2(2));

            ParseBlockElement<double>("CellConstructedVector",inpfile,2,1,u3(0));
            ParseBlockElement<double>("CellConstructedVector",inpfile,2,2,u3(1));
            ParseBlockElement<double>("CellConstructedVector",inpfile,2,3,u3(2));


            //UnitCellNumber=ParseValue<int>("AtomNum",-1);
            orbNum=ParseBlockTitle<int>("Orbital",inpfile,0);

            double PVolumn=u1.dot((u2.cross(u3)));
            ug1=u2.cross(u3)/PVolumn;
            ug2=u3.cross(u1)/PVolumn;
            ug3=u1.cross(u2)/PVolumn;


            char coor[12]="Coordinates";
            UnitCell.resize(UnitCellNumber);
            for(int i=0;i<UnitCellNumber;++i)
            {
                ParseBlockElement<double>("Coordinates",inpfile,i,1,UnitCell[i](0));
                ParseBlockElement<double>("Coordinates",inpfile,i,2,UnitCell[i](1));
                ParseBlockElement<double>("Coordinates",inpfile,i,3,UnitCell[i](2));
            }

            Unita1=a1*pow(2,DividedFactor);
            Unita2=a2;
            Unita3=a3;

            double SVolumn=Unita1.dot((Unita2.cross(Unita3)));
            Unitg1=Unita2.cross(Unita3)/SVolumn;
            Unitg2=Unita3.cross(Unita1)/SVolumn;
            Unitg3=Unita1.cross(Unita2)/SVolumn;
        }
/*
        void Cell::ParseSegments_Z()
        {
            Zdest=ParseValue<int>("ZSegmentsLength",-888,0)+2;
            ZspNum.resize(Zdest); 
            for(int i=0;i<Zdest;++i)
            {
                int foo=i;
                if(i==0) foo=Zdest-2;
                ParseBlockElement<int>("ZSegments",0,(foo-1)%(Zdest-2),ZspNum[i]);
            }
            

            ZTypeCollection.resize(0);
            ZFifthLine.resize(0);
            std::vector<int> tmpfoo;
            for(int i=0;i<Zdest-2;++i)
            {
                int tmpint;
                ParseBlockElement<int>("ZSegments",2,i,tmpint);
                tmpfoo.push_back(tmpint);
                ZFifthLine.push_back(tmpint);
                int count=0;
                for(int j=tmpfoo.size()-1;j>=0;--j)
                {  
                    if(tmpfoo[i]!=tmpfoo[j]) break;
                    else count++;
                }
                if(ZTypeCollection.size()==0 | ZTypeCollection.back().real()!=tmpint) ZTypeCollection.push_back(std::complex<int>{tmpint,--count});
                if(ZTypeCollection.back().real()==tmpint) ZTypeCollection.back()=std::complex<int>{tmpint,count};
            }
            tmpfoo.resize(0);

            for(int i=0;i<ZTypeCollection.size();++i)
            {
                int pos=0;
                for(int j=i-1;j>=0;--j)
                {
                    pos+=ZTypeCollection[j].imag()?ZTypeCollection[j].imag():1;
                }            
                ZTypeCollection[i]=std::complex<int>{pos,ZTypeCollection[i].imag()};
            }



            std::vector<int> JointValue;
            ZspStr.resize(0); 
            ZBlockStr.resize(0);
            ZspJoint.resize(0);
            for(int i=0;i<Zdest-2;++i)
            {
                int tmpint;
                ParseBlockElement<int>("ZSegments",0,i,tmpint);
                ZspJoint.push_back(std::complex<int>{accumulate(JointValue.begin(),JointValue.end(),0),tmpint});
                JointValue.push_back(tmpint);
                std::string tmpspStr;
                ParseBlockElement<std::string>("ZSegments",1,i,tmpspStr);
                ZBlockStr.push_back(tmpspStr);
                for(int s=0;s<tmpint*uvector.size();++s)
                {
                    ZspStr.push_back(tmpspStr);
                }
            }
            JointValue.resize(0);       
        }

        void Cell::ParseSegments_Y()
        {
            Ydest=ParseValue<int>("YSegmentsLength",-888,0)+2;
            YspNum.resize(Ydest); 
            for(int i=0;i<Ydest;++i)
            {
                int foo=i;
                if(i==0) foo=Ydest-2;
                ParseBlockElement<int>("YSegments",0,(foo-1)%(Ydest-2),YspNum[i]);
            }
            

            YTypeCollection.resize(0);
            YFifthLine.resize(0);
            std::vector<int> tmpfoo;
            for(int i=0;i<Ydest-2;++i)
            {
                int tmpint;
                ParseBlockElement<int>("YSegments",2,i,tmpint);
                tmpfoo.push_back(tmpint);
                YFifthLine.push_back(tmpint);
                int count=0;
                for(int j=tmpfoo.size()-1;j>=0;--j)
                {  
                    if(tmpfoo[i]!=tmpfoo[j]) break;
                    else count++;
                }
                if(YTypeCollection.size()==0 | YTypeCollection.back().real()!=tmpint) YTypeCollection.push_back(std::complex<int>{tmpint,--count});
                if(YTypeCollection.back().real()==tmpint) YTypeCollection.back()=std::complex<int>{tmpint,count};
            }
            tmpfoo.resize(0);

            for(int i=0;i<YTypeCollection.size();++i)
            {
                int pos=0;
                for(int j=i-1;j>=0;--j)
                {
                    pos+=YTypeCollection[j].imag()?YTypeCollection[j].imag():1;
                }            
                YTypeCollection[i]=std::complex<int>{pos,YTypeCollection[i].imag()};
            }



            std::vector<int> JointValue;
            YspStr.resize(0); 
            YBlockStr.resize(0);
            YspJoint.resize(0);
            for(int i=0;i<Ydest-2;++i)
            {
                int tmpint;
                ParseBlockElement<int>("YSegments",0,i,tmpint);
                YspJoint.push_back(std::complex<int>{accumulate(JointValue.begin(),JointValue.end(),0),tmpint});
                JointValue.push_back(tmpint);
                std::string tmpspStr;
                ParseBlockElement<std::string>("YSegments",1,i,tmpspStr);
                YBlockStr.push_back(tmpspStr);
                for(int s=0;s<tmpint*uvector.size();++s)
                {
                    YspStr.push_back(tmpspStr);
                }
            }
            JointValue.resize(0);       
        }
*/

        void Cell::ParseSegments_X(std::string inpfile)
        {
            //dest=ParseValue<int>("XSegmentsLength",-888,0)+2;
	    dest=ParseBlockTitle("XSegments",inpfile,-999)+2;
            spNum.resize(dest); 
            for(int i=0;i<dest;++i)
            {
                int foo=i;
                if(i==0) foo=dest-2;
                ParseBlockElement<int>("XSegments",inpfile,(foo-1)%(dest-2),2,spNum[i]);
            }
            pushVec.resize(dest);
            for(int i=0;i<dest;++i)
            {
                int tmp=0;
                for(int j=0;j<i;++j)
                {
                    tmp+=spNum[j];
                }
                pushVec[i]=tmp;
            }

            TypeCollection.resize(0);
            FifthLine.resize(0);
            std::vector<int> tmpfoo;
            for(int i=0;i<dest-2;++i)
            {
                int tmpint;
                ParseBlockElement<int>("XSegments",inpfile,i,0,tmpint);
                tmpfoo.push_back(tmpint);
                FifthLine.push_back(tmpint);
                int count=0;
                for(int j=tmpfoo.size()-1;j>=0;--j)
                {  
                    if(tmpfoo[i]!=tmpfoo[j]) break;
                    else count++;
                }
                if(TypeCollection.size()==0 | TypeCollection.back().real()!=tmpint) TypeCollection.push_back(std::complex<int>{tmpint,--count});
                if(TypeCollection.back().real()==tmpint) TypeCollection.back()=std::complex<int>{tmpint,count};
            }
            tmpfoo.resize(0);

            for(int i=0;i<TypeCollection.size();++i)
            {
                int pos=0;
                for(int j=i-1;j>=0;--j)
                {
                    pos+=TypeCollection[j].imag()?TypeCollection[j].imag():1;
                }            
                TypeCollection[i]=std::complex<int>{pos,TypeCollection[i].imag()};
            }



            std::vector<int> JointValue;
            spStr.resize(0); 
	    vecEnergyCutoff.resize(0);
            BlockStr.resize(0);
            spJoint.resize(0);
            for(int i=0;i<dest-2;++i)
            {
                int tmpint;
                ParseBlockElement<int>("XSegments",inpfile,i,2,tmpint);
                spJoint.push_back(std::complex<int>{accumulate(JointValue.begin(),JointValue.end(),0),tmpint});
                JointValue.push_back(tmpint);
                std::string tmpspStr;
                ParseBlockElement<std::string>("XSegments",inpfile,i,1,tmpspStr);
                BlockStr.push_back(tmpspStr);

		double tmpEnergyCutoffFirst,tmpEnergyCutoffSecond;
                ParseBlockElement<double>("XSegments",inpfile,i,3,tmpEnergyCutoffFirst);
                ParseBlockElement<double>("XSegments",inpfile,i,4,tmpEnergyCutoffSecond);
                for(int s=0;s<tmpint*uvector.size();++s)
                {
                    spStr.push_back(tmpspStr);
		    vecEnergyCutoff.push_back(std::pair<double,double>(tmpEnergyCutoffFirst,tmpEnergyCutoffSecond));
                }
            }
	    //std::cout<<"vec_size="<<vecEnergyCutoff.size()<<std::endl;
            JointValue.resize(0);       
        }

        void Cell::ParseCellExtra(std::string inpfile)
        {
            ParseSegments_X(inpfile);
            //ParseSegments_Y();
            //ParseSegments_Z();

            double PVolumn=a1.dot((a2.cross(a3)));
            g1=a2.cross(a3)/PVolumn;
            g2=a3.cross(a1)/PVolumn;
            g3=a1.cross(a2)/PVolumn;
            //Unita1=a1*(pushVec[dest-1]-pushVec[1]);//systemLength;
            Unita1=a1*accumulate(spNum.begin()+1,spNum.end()-1,0);//systemLength;
            Unita2=a2;
            Unita3=a3;
            

            double SVolumn=Unita1.dot((Unita2.cross(Unita3)));
            Unitg1=Unita2.cross(Unita3)/SVolumn;
            Unitg2=Unita3.cross(Unita1)/SVolumn;
            Unitg3=Unita1.cross(Unita2)/SVolumn;

        }
    }

