#include "Material.hpp"
#include "Cell.hpp"

#ifndef TB_Hamiltonian_H
#define TB_Hamiltonian_H


namespace TightBinding
{
    template<typename T> struct less {};
    template<typename T>
	struct less<std::complex<T> >
	{
	    bool operator()(std::complex<T> const& a, std::complex<T> const& b) const
	    {
		return std::array<T,2>{a.real(),a.imag()} < std::array<T,2>{b.real(),b.imag()};
	    }
	};

    class Hamiltonian
    {
	private:
	public:
	    double BandOffSet;
	    double BandGapIndex;
	    double scaleFactor;
	    //int ctxt;
	    //std::vector<std::unique_ptr<GenTheList>> list;
	    std::vector<std::vector<std::complex<double>>> pD;
	    std::map<double,std::vector<std::vector<std::complex<double>>>> mappD;

	    int orbNum=0;
	    int UnitCellNumber=0;
	    int KextraNum=0;
	    std::vector<std::string> KextraName;
	    std::vector<std::pair<double,double>> KextraRange;

	    Eigen::ComplexEigenSolver<Eigen::MatrixXcd> csolver;
	    Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> solver;
	    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> gsolver;

	    std::map<std::complex<int>,std::vector<std::unique_ptr<EdgeCoupling>>,less<std::complex<int>>> thelist;
	    std::map<std::complex<int>,Eigen::MatrixXcd,less<std::complex<int>>> Materialblock;
	    std::map<std::complex<int>,Eigen::Matrix<lapack_complex_double,-1,-1>,less<std::complex<int>>> Materialblock2;
	    std::map<std::complex<int>,lapack_complex_double,less<std::complex<int>>> SparseMatrix;
	    Eigen::MatrixXcd dMatrix;
	    Eigen::Matrix<lapack_complex_double,-1,-1> dMatrix2;
	    Eigen::MatrixXcd matrix3Db;
	    Eigen::MatrixXcd matrixScatter;
        std::vector<Eigen::MatrixXcd> vecMatrixFold;
	    Eigen::MatrixXcd transform_X;
	    std::vector<Eigen::MatrixXcd> KextraTransform;
	    Eigen::MatrixXcd transform_L;
	    Eigen::MatrixXcd transform_mix;
	    Eigen::MatrixXcd transform_matrix;
	    Eigen::MatrixXcd Uq;
	    Eigen::MatrixXcd extendMat;
	    Eigen::MatrixXcd extendSquareMat;
	    Eigen::MatrixXcd matrix,matrixfull;
	    std::vector<Eigen::Vector4i> qlabels;
	    Materials material;
	    Cell systemcell;

	    std::vector<std::unique_ptr<Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd>>> gsoul;
	    std::vector<std::unique_ptr<Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd>>> soul;
	    std::vector<Eigen::MatrixXcd> EigenVectorsBlocks;
	    std::vector<Eigen::Matrix<lapack_complex_double,-1,-1>> EigenVectorsBlocks2;
	    std::vector<Eigen::VectorXd> EigenValuesBlocks;
std::vector<Eigen::MatrixXcd> PureFold_matrix2;
        std::vector<Eigen::VectorXd> PureFold_vector2;
        std::vector<Eigen::MatrixXcd> PureFold_matrix;
        std::vector<Eigen::VectorXd> PureFold_vector;
	    std::vector<Eigen::MatrixXcd> gsoulFold_matrix;
	    std::vector<Eigen::Matrix<lapack_complex_double,-1,-1>> PgsoulFold_matrix;
	    std::vector<Eigen::VectorXd> gsoulFold_vector;
	    std::vector<std::unique_ptr<Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd>>> gsoul3D;

        void MtoS(Eigen::MatrixXcd &M,Eigen::MatrixXcd &S,int flag);
        void SxS(Eigen::MatrixXcd &Sa,Eigen::MatrixXcd &Sb);
	    void PrepareHamiltonian(Eigen::MatrixXcd &blockMatrix,int tot);
	    void GatherHamiltonian(Eigen::VectorXd &WW,Eigen::MatrixXcd &ZZ,int tot);
	    void FoldTheBlockHamiltonian(int cellme,int blpo,int AtomO,std::vector<int> vblpo,std::vector<int> vAtomO,const Eigen::Vector3d& k,std::vector<Eigen::Vector3d> chain,std::pair<double,double> EnergyCutoff,int &Startp,int &SecDim,Eigen::MatrixXcd &sfitting);
	    void PrintGeometryTotal();
	    void PrintGeometry(std::vector<std::unique_ptr<Geometry>> & printCell);
	    void FoldTheBlockHamiltonian2(int cellme,int blpo,int AtomO,std::vector<int> vblpo,std::vector<int> vAtomO,const Eigen::Vector3d& k,std::vector<Eigen::Vector3d> chain,Eigen::MatrixXcd &sfitting,Eigen::VectorXd &vfitting,int dtag);
	    void FoldTheBlockHamiltonian2(int cellme,int blpo,int AtomO,std::vector<int> vblpo,std::vector<int> vAtomO,const Eigen::Vector3d& k,std::vector<Eigen::Vector3d> chain,Eigen::MatrixXcd &sfitting);
        void GetVK(Eigen::MatrixXcd &vk,Eigen::MatrixXcd &kHvect,Eigen::VectorXcd &kkv,Eigen::Vector3d npl,std::complex<double> Eomega);
	    void PFTBH(int cellme,int blpo,int AtomO,std::vector<int> vblpo,std::vector<int> vAtomO,const Eigen::Vector3d& k,std::vector<Eigen::Vector3d> chain,Eigen::Matrix<lapack_complex_double,-1,-1> &sMatrix,int &Startp,int &SecDim,Eigen::Matrix3d &Ugm,Eigen::Matrix3d &ugm,Eigen::Matrix3d &uam,std::vector<Eigen::MatrixXcd> &egvct,Eigen::VectorXd &svalue2,int &myrow,int &mycol);
	    void tmpPFTBH(int cellme,int blpo,int AtomO,std::vector<int> vblpo,std::vector<int> vAtomO,const Eigen::Vector3d& k,std::vector<Eigen::Vector3d> chain,Eigen::Matrix<lapack_complex_double,-1,-1> &sMatrix,int &Startp,int &SecDim,Eigen::Matrix3d &Ugm,Eigen::Matrix3d &ugm,Eigen::Matrix3d &uam,int &myrow,int &mycol);
	    void FTBH(int cellme,int blpo,int AtomO,std::vector<int> vblpo,std::vector<int> vAtomO,const Eigen::Vector3d& k,std::vector<Eigen::Vector3d> chain,Eigen::MatrixXcd &sfitting,int &Startp,int &SecDim,Eigen::Matrix3d &Ugm,Eigen::Matrix3d &ugm,Eigen::Matrix3d &uam,std::vector<Eigen::MatrixXcd> &egvct,Eigen::VectorXd &svalue2);
	    void tmpFTBH(int cellme,int blpo,int AtomO,std::vector<int> vblpo,std::vector<int> vAtomO,const Eigen::Vector3d& k,std::vector<Eigen::Vector3d> chain,Eigen::MatrixXcd &sfitting,int &Startp,int &SecDim,Eigen::Matrix3d &Ugm,Eigen::Matrix3d &ugm,Eigen::Matrix3d &uam);
	    int FoldTheCol(std::vector<double> &egv,std::vector<Eigen::MatrixXcd> &egvct,int cellme,int blpo,int AtomO,std::vector<int> vblpo,std::vector<int> vAtomO,const Eigen::Vector3d& k,std::vector<Eigen::Vector3d> chain,std::pair<double,double> EnergyCutoff,int newcolnum,int colnum,Eigen::Matrix3d &Ugm,Eigen::Matrix3d &ugm,Eigen::Matrix3d &uam,Eigen::MatrixXcd &yfitting,int dtag);
	    int FoldTheCol3(std::vector<double> &egv,std::vector<Eigen::MatrixXcd> &egvct,int cellme,int blpo,int AtomO,std::vector<int> vblpo,std::vector<int> vAtomO,const Eigen::Vector3d& k,std::vector<Eigen::Vector3d> chain,std::pair<double,double> EnergyCutoff,int newcolnum,int colnum,Eigen::Matrix3d &Ugm,Eigen::Matrix3d &ugm,Eigen::Matrix3d &uam,Eigen::MatrixXcd &yfitting,int dtag);
	    int FoldTheCol3(const Eigen::MatrixXcd &egvct,int cellme,int blpo,int AtomO,std::vector<int> vblpo,std::vector<int> vAtomO,const Eigen::Vector3d& k,std::vector<Eigen::Vector3d> chain,std::pair<double,double> EnergyCutoff,int newcolnum,int colnum,Eigen::Matrix3d &Ugm,Eigen::Matrix3d &ugm,Eigen::Matrix3d &uam,Eigen::MatrixXcd &yfitting,int dtag);
	    void FCol3(const Eigen::VectorXcd &egvct,int cellme,int blpo,int AtomO,std::vector<int> vblpo,std::vector<int> vAtomO,const Eigen::Vector3d& k,std::vector<Eigen::Vector3d> chain,int newcolnum,int i,int cc,int kj,Eigen::Matrix3d &Ugm,Eigen::Matrix3d &ugm,Eigen::Matrix3d &uam,Eigen::MatrixXcd &yfitting);
	    void PFCol3(const Eigen::VectorXcd &egvct,int cellme,int blpo,int AtomO,std::vector<int> vblpo,std::vector<int> vAtomO,const Eigen::Vector3d& k,std::vector<Eigen::Vector3d> chain,int newcolnum,int i,int cc,int kj,Eigen::Matrix3d &Ugm,Eigen::Matrix3d &ugm,Eigen::Matrix3d &uam,Eigen::Matrix<lapack_complex_double,-1,-1> &sMatrix,int &myrow,int &mycol);
	    void FoldTheCol2(std::vector<double> &egv,std::vector<Eigen::MatrixXcd> &egvct,int cellme,int blpo,int AtomO,std::vector<int> vblpo,std::vector<int> vAtomO,const Eigen::Vector3d& k,std::vector<Eigen::Vector3d> chain,std::pair<double,double> EnergyCutoff,int newcolnum,int colnum,Eigen::Matrix3d &Ugm,Eigen::Matrix3d &ugm,Eigen::Matrix3d &uam,Eigen::MatrixXcd &yfitting,int dtag);
	    void QuickFoldTheBlockHamiltonian(int cellme,int blpo,int AtomO,std::vector<int> vblpo,std::vector<int> vAtomO,const Eigen::Vector3d& k,std::vector<Eigen::Vector3d> chain,std::pair<double,double> EnergyCutoff,int &Startp,int &SecDim,int dtag,Eigen::VectorXd &svalue,Eigen::MatrixXcd &sfitting);
	    void PQuickFoldTheBlockHamiltonian(int cellme,int blpo,int AtomO,std::vector<int> vblpo,std::vector<int> vAtomO,const Eigen::Vector3d& k,std::vector<Eigen::Vector3d> chain,std::pair<double,double> EnergyCutoff,int &Startp,int &SecDim,int dtag,Eigen::VectorXd &svalue,Eigen::Matrix<lapack_complex_double,-1,-1> &sfitting,int &myrow,int &mycol);

	    //data
	    std::unordered_map<int,std::string> str_delta_m;
	    double rightCutoff;
	    double leftCutoff;
	    double MixCoeff;
	    int USE_FixedSize_or_EnergyRange;
	    int USE_SOC;
	    int USE_long_parameters;
	    double USE_Shift;
	    int AutoMatic;
	    int USE_extra_k_point;
	    int USE_OneRange_or_MoreRanges;
	    int USE_direct_construction;
	    int USE_RowTruncation;
	    int USE_decomposition;
	    int USE_ColTruncation;
	    int USE_MaterialMix;
	    void destroy();
	    void ParseParameters(std::string inpfile);
	    void print_Hamiltonian();
	    int in_the_list(int ii,int j);
	    std::string Stringin_the_list(int ii,int j);
	    int MatrixSize;
	    int matrix_bottom;
	    int MatrixSizeReduced;
	    void pDiagonalize(std::vector<std::unique_ptr<Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd>>> &SubblockContainer,Eigen::MatrixXcd &ConsideredMatrix,bool onoff,Eigen::MatrixXcd &MatrixToBeTransformed);
	    void Filtrate_the_states(const Eigen::VectorXd &soup,const std::unique_ptr<Geometry> &cell,bool newone);
	    void Filtrate_the_states(const Eigen::VectorXd &soup,const std::unique_ptr<Geometry> &cell,bool newone,int &Startp,int &SecDim);

	    void ScatterHamiltonian(Eigen::MatrixXcd &blockMatrix,Eigen::Matrix<lapack_complex_double,-1,-1> &mat_loc,int tot,int Nb,int Mb);
	    void GetUq(Eigen::Vector3d q);
	    void DataCom(int blposX,int blposY,int bN,int bM,int N,int M,int Nb,int Mb,Eigen::MatrixXcd &blockMatrix,Eigen::Matrix<lapack_complex_double,-1,-1> &mat_loc,int &nrows,int &ncols,int &myrow,int &mycol,double flag,int &ctxt);
	    void Initialize(int paramd,int midlen,int f1,int f2);
	    ///function
	    double SKFormular(std::string leftType,std::string rightType,int n1,int l1,int m1,int n2,int l2,int m2,Eigen::Vector3d dn,const Material *matter);
	    void ReadMaterials();
        void read_transform_matrix(Eigen::MatrixXcd & eigvects,const char* filename);
        void readEigenValues(Eigen::VectorXd& eigvals,const char* filename);
        void readEigenVectors(Eigen::MatrixXcd & eigvects,const char* filename,int buffsize);
        void store_precomputed_data(const Eigen::Vector3d& kp,int data_index,const Eigen::VectorXd& SAVEDeigenvalues,const Eigen::MatrixXcd& SAVEDeigenvectors);
        void store_materialblocks(std::string name,Eigen::MatrixXcd mat);
        std::string store_transform_matrix(const Eigen::Vector3d& kp);
	    void ArrangeT(Eigen::VectorXcd &values,Eigen::MatrixXcd &vectors);
	    int ArrangeDRE(Eigen::VectorXcd &values,Eigen::MatrixXcd &vectors);
	    void Diagonal_parallel_pzheevx(int tot);
	    void Diagonal_parallel_pzheevx(Eigen::VectorXd &WW,Eigen::MatrixXcd &ZZ,int blposX,int blposY,int tot,int Nb,int Mb,int NN,int flag,Eigen::Matrix<lapack_complex_double,-1,-1> blockMatrix);
	    void Diagonal_parallel_pzheevx_all(Eigen::VectorXd &WW,Eigen::MatrixXcd &ZZ,int blposX,int blposY,int tot,int Nb,int Mb,int NN,int flag);
	    void Diagonal_parallel_pzheevx_test(Eigen::VectorXd &WW,Eigen::MatrixXcd &ZZ,int blposX,int blposY,int tot,int Nb,int Mb,int NN,int flag);
	    void Diagonal_parallel_pzheevx_blocks(int blpos,int tot);
	    int Diagonal_parallel_pdsyevx(double *H,int m,int flag,int tot);
	    void ArrangeDIM(Eigen::VectorXcd &values,Eigen::MatrixXcd &vectors,Eigen::MatrixXcd &II);
	    int  ArrangeDD(Eigen::VectorXcd &values,Eigen::MatrixXcd &vectors,Eigen::MatrixXcd &II,std::vector<Eigen::MatrixXcd> &vW);
	    void ArrangeDDD(Eigen::VectorXcd &values,Eigen::MatrixXcd &vectors);
	    void Dswap(Eigen::VectorXcd &values,Eigen::MatrixXcd &vectors);
	    void Dreverse(Eigen::VectorXcd &values,Eigen::MatrixXcd &vectors);
	    void WaveFunctions();
	    void printOutWaveFunctions();
	    void NormTR(Eigen::MatrixXcd &tt,Eigen::MatrixXcd &rr,int bksize);
	    void GenVT(Eigen::MatrixXcd matTOT,std::vector<Eigen::MatrixXcd> &vT,std::vector<Eigen::MatrixXcd> &vW);
	    void GenH(std::complex<double> Eomega,Eigen::MatrixXcd &emat,std::vector<Eigen::MatrixXcd> &vW,int bksize,Eigen::Vector3d kk,Eigen::Vector3d npl);
	    void GenH2(std::complex<double> Eomega,Eigen::MatrixXcd &emat,int bksize,int scsize,Eigen::Vector3d kk,Eigen::Vector3d npl);
	    void eGenVT(Eigen::MatrixXcd &matTOT,std::vector<Eigen::MatrixXcd> &vT,std::vector<Eigen::MatrixXcd> &vT2,std::vector<Eigen::MatrixXcd> &vW,std::vector<Eigen::MatrixXcd> &vG,int bksize,int scsize);
	    void Test_vW(Eigen::MatrixXcd matTOT,std::vector<Eigen::MatrixXcd> &vW,std::vector<Eigen::MatrixXcd> &vT,Eigen::MatrixXcd &J0,Eigen::MatrixXcd &JN,int bksize);
	    void FillDiagonalBlocksByTransform(Eigen::MatrixXcd &ConsideredMatrix,std::vector<std::unique_ptr<Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd>>> &SubblockData);
	    void preDiagonalize(std::vector<std::unique_ptr<Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd>>> &SubblockContainer,Eigen::MatrixXcd &ConsideredMatrix,bool onoff,Eigen::MatrixXcd &MatrixToBeTransformed);
	    void pzheevxPreDiagonalize(const Eigen::Vector3d k,std::vector<Eigen::MatrixXcd> &SubblockContainerOfEigenVectors,std::vector<Eigen::VectorXd> &SubblockContainerOfEigenValues,Eigen::MatrixXcd &ConsideredMatrix,bool onoff,Eigen::MatrixXcd &MatrixToBeTransformed,int flag);
	    void pzheevx3PreDiagonalize(std::vector<Eigen::MatrixXcd> &SubblockContainerOfEigenVectors,std::vector<Eigen::VectorXd> &SubblockContainerOfEigenValues,Eigen::MatrixXcd &ConsideredMatrix,bool onoff,Eigen::MatrixXcd &MatrixToBeTransformed,int flag);
	    void lCase5(std::complex<double> Eomega);
	    void Case5(std::complex<double> Eomega);
	    void ShiftMatrix(Eigen::MatrixXcd &Omatrix,Eigen::MatrixXcd &Smatrix);
	    void CaseComplexBands5(std::complex<double> Eomega,Eigen::Vector3d k);
	    void CaseSurfaceDOS(std::complex<double> Eomega,Eigen::Vector3d k);
	    void ComplexBands(std::complex<double> Eomega,Eigen::Vector3d k);
	    void gFillDiagonalBlocksByData(Eigen::MatrixXcd &ConsideredMatrix,std::vector<std::unique_ptr<Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd>>> &SubblockData);
	    void FillDiagonalBlocksByData(Eigen::MatrixXcd &ConsideredMatrix,std::vector<std::unique_ptr<Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd>>> &SubblockData);
	    void pzheevxFillDiagonalBlocksByData(Eigen::MatrixXcd &ConsideredMatrix,std::vector<Eigen::VectorXd> &SubblockData);
	    void FillDiagonalBlocksByData(Eigen::MatrixXcd &ConsideredMatrix,std::vector<Eigen::VectorXd> &SubblockData);
	    void Method1(Eigen::Vector3d npl,Eigen::MatrixXcd matTOT,Eigen::MatrixXcd &tt,Eigen::MatrixXcd &rr,Eigen::MatrixXcd &matII,std::vector<Eigen::MatrixXcd> vT,Eigen::MatrixXcd T,Eigen::MatrixXcd T2,std::complex<double> energy,int bksize,int scsize,std::vector<Eigen::MatrixXcd> &vW,Eigen::MatrixXcd &D1vk,Eigen::MatrixXcd &D2vk);
	    void SmallMethod1(Eigen::MatrixXcd matTOT,Eigen::MatrixXcd &tt,Eigen::MatrixXcd &rr,Eigen::MatrixXcd &matII,std::vector<Eigen::MatrixXcd> vT,Eigen::MatrixXcd T,Eigen::MatrixXcd T2,std::complex<double> energy,int bksize,int scsize,std::vector<Eigen::MatrixXcd> &vW);
	    void Method2(Eigen::MatrixXcd matTOT,Eigen::MatrixXcd &tt,Eigen::MatrixXcd &rr,Eigen::MatrixXcd &matII,std::vector<Eigen::MatrixXcd> vT,Eigen::Vector3d pl,Eigen::MatrixXcd T,Eigen::MatrixXcd T2);
	    void Method3(Eigen::MatrixXcd matTOT,Eigen::MatrixXcd &tt,Eigen::MatrixXcd &rr,Eigen::MatrixXcd &matII,std::vector<Eigen::MatrixXcd> vT,Eigen::Vector3d pl,Eigen::MatrixXcd T,Eigen::MatrixXcd T2);
	    void Method4(Eigen::MatrixXcd matTOT,Eigen::MatrixXcd &tt,Eigen::MatrixXcd &rr,std::vector<Eigen::MatrixXcd> vT,Eigen::Vector3d pl,Eigen::MatrixXcd T,Eigen::MatrixXcd T2);
	    double Diagonalize_reduced_H(bool io);
	    void BlockCoupling(Eigen::MatrixXcd &ReducedMatrix,std::vector<std::unique_ptr<Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd>>> &BlockData);
	    void BlockCoupling(Eigen::MatrixXcd &ReducedMatrix,std::vector<Eigen::MatrixXcd> &BlockData);
	    void gBlockCoupling(Eigen::MatrixXcd &ReducedMatrix,std::vector<std::unique_ptr<Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd>>> &BlockData);
	    void pzheevxBlockCoupling(Eigen::MatrixXcd &ReducedMatrix,std::vector<Eigen::MatrixXcd> &BlockData);
	    double blockCoupling3D3D(const Eigen::Vector3d& k);
	    void ReadCouplingElements(Eigen::MatrixXcd &ConsideredMatrix);
	    void ReadCouplingElements3(Eigen::MatrixXcd &ConsideredMatrix);
	    void DestoryCouplingElements(Eigen::MatrixXcd &ConsideredMatrix);
	    void ReadCouplingElements2(Eigen::MatrixXcd &ConsideredMatrix);
	    Hamiltonian(){
		str_delta_m[0]="sigma";
		str_delta_m[1]="pi";
		str_delta_m[2]="delta";
	    };
	    ~Hamiltonian();
	    inline int mposi(int i); 
	    inline int mposi(int i,int flag); 
	    inline static std::complex<double> g0(const Eigen::Vector3d& k);
	    inline static std::complex<double> g1(const Eigen::Vector3d& k);
	    inline static std::complex<double> g2(const Eigen::Vector3d& k);
	    inline static std::complex<double> g3(const Eigen::Vector3d& k);


	    void print_matrix(Eigen::MatrixXcd mat,std::string name);
	    void print_matrix_in_details(Eigen::MatrixXcd mat,std::string name);
	    void printSupercell(std::string name);
	    void print_rmatrix(Eigen::MatrixXd mat,std::string name);
	    void update_params(Eigen::Vector3d dn);
	    void update_pD();
	    void SetMatrixSparse(const Eigen::Vector3d k,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::MatrixXcd &matrix);
	    void SetMatrixSparse2(const Eigen::Vector3d k,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::MatrixXcd &matrix);
	    void SetMatrix(const Eigen::Vector3d k,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::MatrixXcd &matrix);
	    void NewSetMatrix(const Eigen::Vector3d k,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::MatrixXcd &matrix);
	    void SetMatrixFold(const Eigen::Vector3d k,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::MatrixXcd &matrix);
	    void SetMatrixPrime(const Eigen::Vector3d k,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::MatrixXcd &matrix,int vblpoi,int blen);
	    void vecSetMatrixFold(const Eigen::Vector3d k,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::MatrixXcd &matrix,std::string material);
	    void SetMatrixLoc(const Eigen::Vector3d k,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::Matrix<lapack_complex_double,-1,-1> &matrix);
	    void SetMatrixLocNew(const Eigen::Vector3d k,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::Matrix<lapack_complex_double ,-1,-1> &matrix,int blpo,int AtomO);
	    void SetMatrixLocWithoutShift(const Eigen::Vector3d k,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::Matrix<lapack_complex_double,-1,-1> &matrix);
	    void SetMatrixA(Eigen::Vector3d kk,std::complex<double> kz,Eigen::Vector3d pl,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::MatrixXcd &matrix);
        void SetMatrixA(const Eigen::VectorXcd k,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::MatrixXcd &matrix);
	    void SetMatrix5(const Eigen::Vector3d kk,Eigen::Vector3d pl,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::MatrixXcd &matrix,int Wsplit);
	    void SetMatrixK5(const Eigen::MatrixXcd &Zmatrix,std::complex<double> kz,Eigen::Vector3d pl,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::MatrixXcd &Kmatrix5,int Wsplit);
	    void SetMatrixK5(std::complex<double> kz,Eigen::Vector3d pl,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::MatrixXcd &Kmatrix5,int Wsplit);
	    void SetMatrixKW5(std::complex<double> kz,Eigen::Vector3d pl,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements,Eigen::MatrixXcd &Kmatrix5,int Wsplit);
	    void SetMatrix(const Eigen::Vector3d& k);
	    double func_A(int m,Eigen::Vector3d dn);
	    double func_B(int m,Eigen::Vector3d dn);
	    double func_T(int l,int m,int m_prime,Eigen::Vector3d dn);
	    double func_S(int l,int m,int m_prime,Eigen::Vector3d dn);
	    void MatrixO(Eigen::MatrixXcd &matrix);
	    double factorials(int n);
	    double func_D(int l,int m,int m_prime,Eigen::Vector3d dn);
	    double getfunc_D(int l,int m,int m_prime,double N);
	    double SetMatrix_folding(const Eigen::Vector3d k);
	    void FillDiagonalBlocksByTransform(Eigen::MatrixXcd &ReducedMatrix,std::vector<Eigen::MatrixXcd> &BlockData);
	    void SparseTransform(Eigen::MatrixXcd &ReducedMatrix,std::vector<Eigen::MatrixXcd> &BlockData);
	    void PTransform(Eigen::Matrix<lapack_complex_double,-1,-1> &ReducedMatrix,std::vector<Eigen::Matrix<lapack_complex_double,-1,-1>> &BlockData);
	    inline std::complex<double> MatrixCoupling11(const Eigen::MatrixXcd &lM ,int row,int col,const Eigen::MatrixXcd &rM,std::complex<int> & maplist);

	    void GenEdgeFast3D(const Eigen::Vector3d& k);
	    void MatrixDiagonalize(bool onoff,Eigen::MatrixXcd OverlapMatrix);
	    void gMatrixDiagonalize(bool onoff,Eigen::MatrixXcd OverlapMatrix);
	    void MatrixDiagonalize2(bool io,Eigen::MatrixXcd &OverlapMatrix);
	    void foldDiagonalize(std::vector<Eigen::MatrixXcd> &SubblockContainer,Eigen::MatrixXcd &ConsideredMatrix,bool onoff,Eigen::MatrixXcd &MatrixToBeTransformed,const Eigen::Vector3d& k,int dtag);
	    void PfoldDiagonalize(std::vector<Eigen::Matrix<lapack_complex_double,-1,-1>> &SubblockContainer,Eigen::MatrixXcd &ConsideredMatrix,bool onoff,Eigen::MatrixXcd &MatrixToBeTransformed,const Eigen::Vector3d& k,int dtag);
	    const Eigen::VectorXd& eigenvalues() const { return gsolver.eigenvalues(); }
	    const Eigen::MatrixXcd& eigenvectors() const { return gsolver.eigenvectors(); }
	    void get_transform_remote(Eigen::Vector3d kp,Eigen::MatrixXcd &transf,int bandnumber,int displacement,std::vector<std::unique_ptr<Neighbor>> &ListofMatrixElements);
	    void Square_matrix(std::vector<std::unique_ptr<Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd>>> &SubblockContainer,Eigen::MatrixXcd &ConsideredMatrix);
	    void Append_matrix();
    };

}

#endif
