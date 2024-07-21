#include "Util.hpp"
#ifndef TB_Cell_H
#define TB_Cell_H

namespace TightBinding
{
    class EdgeCoupling
    {
        public:
            EdgeCoupling(){};
            EdgeCoupling(int Tab1,int Tab2,std::complex<double> max):Table1(Tab1),Table2(Tab2),value(max){};
            int Table1,Table2;
            std::complex<double> value={0.0,0.0};
    };

    class GenTheList
    {
	public:
	    GenTheList(){};
        GenTheList(std::complex<int> Tab,int lp1,int lp2,int tp1,int tp2,int o,int cellmein):Table(Tab),ListPosition1(lp1),ListPosition2(lp2),TagPosition1(tp1),TagPosition2(tp2),SuperCell(o),cellme(cellmein){};
	    std::complex<int> Table;
	    int ListPosition1,ListPosition2;
	    int TagPosition1,TagPosition2;
	    int cellme;
	    int SuperCell=0;
    };


    class Neighbor
    {
	public:
	    Neighbor(){};
        Neighbor(std::complex<int> Tab,Eigen::Vector3d dvec,std::string o,std::string oo,Eigen::Vector3d ad,int expfl,int Type):Table(Tab),d(dvec),SuperCell(o),cellme(oo),ed(ad),expflag(expfl),Typeme(Type){};

	    std::complex<int> Table;
        int expflag;
        int Typeme;

	    Eigen::Vector3d d;
	    Eigen::Vector3d ed;
	    std::string cellme;
	    std::string SuperCell;
    }; 
    class Geometry
    {
        private:
        public:
            int Tag3;
            int Tag2;
            int splitSize;
            int type=0;
            int tag;
            int blpo;
            int AtomO;
            std::vector<int> vblpo;
            std::vector<int> vAtomO;

            std::vector<int> vblpo2;
            std::vector<int> vAtomO2;
            int AtomNum=0;
            int gatherAtomNum=0;
            int SecDim=0;
            int Startp=0;
            Geometry(){};
            ~Geometry(){Unit.resize(0);};
            std::vector<Eigen::Vector3d> Unit;
            std::vector<std::string> UnitTexture;
            std::vector<std::complex<double>> UnitWave_Functions;
            std::string stexture="GaAs";
            std::vector<Eigen::Vector3d> UnitGather;
	    std::pair<double,double> EnergyCutoff;
        std::vector<std::pair<double,double>> vEnergyCutoff;
////////////
            Geometry(Eigen::Vector3d vec,std::vector<Eigen::Vector3d> &unit_trans,int tag_trans,int Tag3_trans,int Tag2_trans,int type_trans):tag(tag_trans),Tag3(Tag3_trans),Tag2(Tag2_trans),type(type_trans),AtomNum(unit_trans.size()),Unit(unit_trans){
                UnitTexture.resize(0);
                for(auto &cell : Unit){
                    cell+=vec;
                }
            }

            Geometry(Eigen::Vector3d vec,std::vector<Eigen::Vector3d> &UnitCell,int tag1,int tag3,int tag2,std::string texture,std::pair<double,double> pairDouble):tag(tag1),Tag3(tag3),Tag2(tag2),type(1),AtomNum(UnitCell.size()),Unit(UnitCell),EnergyCutoff(pairDouble){
                UnitTexture.resize(0);
                UnitTexture.push_back(texture);
                for(auto &cell : Unit){
                    cell+=vec;
                }
            }



    };



    class Cell 
    {
        private:
        public:

	    int modInAsLength=0;
        int DividedFactor;
        int dest=0;
        int Ydest=0;
        int Zdest=0;
        int USE_3Dpbc;
	    int floorInAsLength=0;
	    double dis=0.0;
        Eigen::Vector3d Unita1;
	    Eigen::Vector3d Unita2;
	    Eigen::Vector3d Unita3;
        Eigen::Vector3d Unitg1;
	    Eigen::Vector3d Unitg2;
	    Eigen::Vector3d Unitg3;

        Eigen::Vector3d g1;
        Eigen::Vector3d g2;
	    Eigen::Vector3d g3;

        Eigen::Vector3d ug1;
        Eigen::Vector3d ug2;
	    Eigen::Vector3d ug3;
        Eigen::Vector3d Ug1;
        Eigen::Vector3d Ug2;
	    Eigen::Vector3d Ug3;
	    std::vector<int> pushVec;
	    std::vector<int> spNum;
	    std::vector<int> YspNum;
	    std::vector<int> ZspNum;
	    std::vector<int> fspNum;
	    std::vector<int> FifthLine;
	    std::vector<int> YFifthLine;
	    std::vector<int> ZFifthLine;
	    int UnitCellNumber=0;
	    std::vector<Eigen::Vector3d> UnitCell;
	    std::vector<Eigen::Vector3d> UnitCellMargin;

	    std::vector<std::string> spStr;
	    std::vector<std::string> YspStr;
	    std::vector<std::string> ZspStr;
	    std::vector<std::string> BlockStr;
	    std::vector<std::string> YBlockStr;
	    std::vector<std::string> ZBlockStr;
	    std::vector<std::complex<int>> spJoint;
	    std::vector<std::complex<int>> TypeCollection;
	    std::vector<std::complex<int>> YspJoint;
	    std::vector<std::complex<int>> YTypeCollection;
	    std::vector<std::complex<int>> ZspJoint;
	    std::vector<std::complex<int>> ZTypeCollection;
	    int orbNum=0;
	    std::vector<std::pair<double,double>> vecEnergyCutoff;

        Eigen::Vector3d a1;
	    Eigen::Vector3d a2;
	    Eigen::Vector3d a3;

        Eigen::Vector3d u1;
	    Eigen::Vector3d u2;
	    Eigen::Vector3d u3;

        std::vector<Eigen::Vector3d> uvector;


        Cell(){};
	    ~Cell(){};
	    std::vector<std::unique_ptr<Geometry>> supercell;
	    std::vector<std::unique_ptr<Geometry>> supercellMargin;

	    std::vector<std::unique_ptr<Neighbor>> ListofMatrixElements_Type1;
	    std::vector<std::unique_ptr<Neighbor>> ListofMatrixElements_MixType1;
	    std::vector<std::unique_ptr<Neighbor>> ListofMatrixElements_Type2;
        std::vector<std::unique_ptr<Neighbor>> ListofMatrixElements_Type3;
	    std::vector<std::unique_ptr<GenTheList>> ListofMatrixCoupling;

        void CellOptimize();
	    void ParseCell(std::string inpfile);
	    void ParseCellExtra(std::string inpfile);
        void ParseSegments_X(std::string inpfile);
        void ParseSegments_Z();
        void ParseSegments_Y();
	    void AddBoundaryConditions(int num);
	    void AddBoundaryConditions2(int num);
    void AddBoundaryConditions();
    void AddBoundaryConditionsS(int numSize);
        void StackingStructure();
        void StackingStructure3D();
	    void JointSupercell(int side2);
	    void JointSupercell2(int side2);
        int GetMatrixSize();
	    void GenTpfast3D(int onoff,std::vector<std::unique_ptr<Neighbor>> &mlist);
	    void MasterGen(int onoff,std::vector<std::unique_ptr<Neighbor>> &mlist);
	    void SelectedMasterGen(int onoff);
	    void OldSelectedMasterGen(int onoff);
        void GenTpfast(std::vector<std::unique_ptr<Neighbor>> &mlist);
        void StackingFolding(std::vector<std::unique_ptr<Neighbor>> &tpfastFolding);
        void GenTransformLists(std::vector<std::unique_ptr<GenTheList>> &list,int onoff);
        void GenSuperCellRequared();
        inline bool Is_Neighbor(Eigen::Vector3d &vec1,Eigen::Vector3d &vec2);
        void Initialize(int thrd,int midlen,int f1,int f2);
        inline int posi(int ii);
        void ProcessParameters(int paramd,int midlen,double leftCutoff,double rightCutoff);
        void ProcessSuperCell(int paramd);
        void ProcessCell();
    };



}
#endif
