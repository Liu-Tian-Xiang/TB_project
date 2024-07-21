#include "Hamiltonian.hpp"
#include "SymmetryPoints.hpp"

#ifndef TB_BANDSTRUCTURE_H
#define TB_BANDSTRUCTURE_H


namespace TightBinding
{


    struct EnumClassHash
    {
	template <typename T>
	    std::size_t operator()(T t) const
	    {
		return static_cast<std::size_t>(t);
	    }
    };


    class Timer
    {
	private:
	    // Type aliases to make accessing nested type easier
	    using clock_t = std::chrono::high_resolution_clock;
	    using second_t = std::chrono::duration<double, std::ratio<1> >;

	    std::chrono::time_point<clock_t> m_beg;

	public:
	    Timer() : m_beg(clock_t::now())
	{
	}

	    void reset()
	    {
		m_beg = clock_t::now();
	    }

	    double elapsed() const
	    {
		return std::chrono::duration_cast<second_t>(clock_t::now() - m_beg).count();
	    }
    };

    class BandStructure
    {
	private:
	    enum TableOFactions
	    {
		Running_5 = 5,//no reduction
		Running_0 = 0,//no reduction
		Running_Epsilon = 3333,//no reduction
		Running_transform2D2D = 4,
		Running_transform3D2D = 6,
		Running_transform3D3D = 9,
		Running_transform2D2Dek = 41,
		Running_transform3D2Dek = 61,
		Running_transform3D3Dek = 91,
		Running_bigtransform2D2Dek = 42,
		Running_bigtransform3D2Dek = 62,
		Running_bigtransform3D3Dek = 92,
		Running_2D2D = 2222,
		Running_pzheevx = 222,
		Running_3D2D = 32,
		Running_3D3D = 33,
		Running_2D2Dek = 221,//with remote k points extension
		Running_3D2Dek = 321,//folding
		Running_3D3Dek = 331,//folding
		Folding = 11,//folding
		apply_method,
		apply_method5,
		method_one_transform,
		apply_method_mpi,
		apply_method3,
		apply_method4,
		apply_method6,
		apply_method7,
		apply_method8,
		apply_method9,
		apply_method10,
		apply_method11,
		do_nothing,
		Init,
		parse,
		AutoMaticInit,
		CustomInit,
		still,
		band_design,
		read_elements,
		read_elements3,
		reconstruct_Hamiltonian_mpi,
		blocks_diagonalization4,
		reconstruct_Hamiltonian2,
		reconstruct_Hamiltonian4,
		reconstruct_Hamiltonian5,
		reconstruct_Hamiltonian6,
		reconstruct_Hamiltonian,
		reconstruct_Hamiltonian3,
		blocks_diagonalization,
		blocks_diagonalization8,
		blocks_diagonalization9,
		blocks_diagonalization3,
		diagonalize,
		diagonalize_pzheevx,
		diagonalize_Eigen,
		diagonalize_Mix,
		output_data,
		output_data_pzheevx,
		pre_calculation,
		transform_Hamiltonian,
		set_Hamiltonian,
		set_Hamiltonian_Type1,
		set_Hamiltonian_Type2,
		set_Hamiltonian_Type3,
		set_Hamiltonian_TypeVec,
		set_Hamiltonian_sparse,
		set_Hamiltonian_sparse2,
		set_Hamiltonian_distributed,
		set_Hamiltonian_distributed_with_shift,
		design_running_pattern,
		design_running_pattern_single,
		create_structure,
		find_neighbors,
		reshape_structure,
		band_design_pzheevx,
		band_design2,
		MeshBZ,
		COMPUTE_K_POINTS,
		COMPUTE_Epsilon,
		COMPUTE_transform2D2D,
		COMPUTE_transform3D2D,
		COMPUTE_transform3D3D,
		COMPUTE_transform2D2Dek,
		COMPUTE_transform3D2Dek,
		COMPUTE_transform3D3Dek,
		COMPUTE_bigtransform2D2Dek,
		COMPUTE_bigtransform3D2Dek,
		COMPUTE_bigtransform3D3Dek,
		COMPUTE_2D2Dek,
		COMPUTE_3D2Dek,
		COMPUTE_3D3Dek,
		COMPUTE_2D2D,
		COMPUTE_pzheevx,
		COMPUTE_3D2D,
		COMPUTE_3D3D,
		COMPUTE_no_reduction,
		COMPUTE_5,
		COMPUTE_Folding,
		With_ek3D,
		With_ek2D,
		Output_and_Recycle,
		DirectFinish,
		RedirectFinish,
		OutputEpsilon,
		Output_and_Recycle_Full,
		End
    };
	public:

    std::unordered_map<TableOFactions,std::string,EnumClassHash> TableOFactionsName=
	//std::unordered_map<int,std::string,EnumClassHash> TableOFactionsName=
    {
	{parse,"parse"},
	{Init,"Init"},
	{CustomInit,"CustomInit"},
	{AutoMaticInit,"AutoMaticInit"},
	{still,"still"},
	{band_design,"band_design"},
	{read_elements,"read_elements"},
	{read_elements3,"read_elements3"},
	{reconstruct_Hamiltonian_mpi,"reconstruct_Hamiltonian_mpi"},
	{blocks_diagonalization4,"blocks_diagonalization4"},
	{reconstruct_Hamiltonian2,"reconstruct_Hamiltonian2"},
	{reconstruct_Hamiltonian4,"reconstruct_Hamiltonian4"},
	{reconstruct_Hamiltonian5,"reconstruct_Hamiltonian5"},
	{reconstruct_Hamiltonian6,"reconstruct_Hamiltonian6"},
	{reconstruct_Hamiltonian,"reconstruct_Hamiltonian"},
	{reconstruct_Hamiltonian3,"reconstruct_Hamiltonian3"},
	{blocks_diagonalization,"blocks_diagonalization"},
	{blocks_diagonalization8,"blocks_diagonalization8"},
	{blocks_diagonalization9,"blocks_diagonalization9"},
	{blocks_diagonalization3,"blocks_diagonalization3"},
	{diagonalize,"diagonalize"},
	{diagonalize_pzheevx,"diagonalize_pzheevx"},
	{diagonalize_Eigen,"diagonalize_Eigen"},
	{diagonalize_Mix,"diagonalize_Mix"},
	{output_data,"output_data"},
	{output_data_pzheevx,"output_data_pzheevx"},
	{pre_calculation,"pre_calculation"},
	{transform_Hamiltonian,"transform_Hamiltonian"},
	{set_Hamiltonian,"set_Hamiltonian"},
	{set_Hamiltonian_Type1,"set_Hamiltonian_Type1"},
	{set_Hamiltonian_Type2,"set_Hamiltonian_Type2"},
	{set_Hamiltonian_Type3,"set_Hamiltonian_Type3"},
	{set_Hamiltonian_TypeVec,"set_Hamiltonian_TypeVec"},
	{set_Hamiltonian_sparse,"set_Hamiltonian_sparse"},
	{set_Hamiltonian_sparse2,"set_Hamiltonian_sparse2"},
	{set_Hamiltonian_distributed,"set_Hamiltonian_distributed"},
	{set_Hamiltonian_distributed_with_shift,"set_Hamiltonian_distributed_with_shift"},
	{design_running_pattern,"design_running_pattern"},
	{design_running_pattern_single,"design_running_pattern_single"},
	{create_structure,"create_structure"},
	{find_neighbors,"find_neighbors"},
	{reshape_structure,"reshape_structure"},
	{band_design_pzheevx,"band_design_pzheevx"},
	{band_design2,"band_design2"},
	{MeshBZ,"MeshBZ"},
	{COMPUTE_K_POINTS,"COMPUTE_K_POINTS"},
	{COMPUTE_Epsilon,"COMPUTE_Epsilon"},
	{COMPUTE_transform2D2D,"COMPUTE_transform2D2D"},
	{COMPUTE_transform3D2D,"COMPUTE_transform3D2D"},
	{COMPUTE_transform3D3D,"COMPUTE_transform3D3D"},
	{COMPUTE_transform2D2Dek,"COMPUTE_transform2D2Dek"},
	{COMPUTE_transform3D2Dek,"COMPUTE_transform3D2Dek"},
	{COMPUTE_transform3D3Dek,"COMPUTE_transform3D3Dek"},
	{COMPUTE_bigtransform2D2Dek,"COMPUTE_bigtransform2D2Dek"},
	{COMPUTE_bigtransform3D2Dek,"COMPUTE_bigtransform3D2Dek"},
	{COMPUTE_bigtransform3D3Dek,"COMPUTE_bigtransform3D3Dek"},
	{COMPUTE_2D2Dek,"COMPUTE_2D2Dek"},
	{COMPUTE_3D2Dek,"COMPUTE_3D2Dek"},
	{COMPUTE_3D3Dek,"COMPUTE_3D3Dek"},
	{COMPUTE_2D2D,"COMPUTE_2D2D"},
	{COMPUTE_pzheevx,"COMPUTE_pzheevx"},
	{COMPUTE_3D2D,"COMPUTE_3D2D"},
	{COMPUTE_3D3D,"COMPUTE_3D3D"},
	{COMPUTE_no_reduction,"COMPUTE_no_reduction"},
	{COMPUTE_5,"COMPUTE_5"},
	{COMPUTE_Folding,"COMPUTE_Folding"},
	{With_ek3D,"With_ek3D"},
	{With_ek2D,"With_ek2D"},
	{Output_and_Recycle,"Output_and_Recycle"},
	{DirectFinish,"DirectFinish"},
	{RedirectFinish,"RedirectFinish"},
	{OutputEpsilon,"OutputEpsilon"},
	{Output_and_Recycle_Full,"Output_and_Recycle_Full"},
	{End,"End"}
    };
    bool FirstLoop=true;
    Eigen::VectorXcd NumMass,Mmass;
    Timer timing;
    Timer loopTiming;
    Hamiltonian hamiltonian;
    SymmetryPoints symmetryPoints;

    Eigen::MatrixXcd Xmatrix;
    Eigen::MatrixXcd eigenvectors;
    Eigen::VectorXd eigenvalues;
    double mesh_factor;
    double integrate_factor;
    int bandbegin;
    int bandend;
    int range_py;
    int range_output;
    int DataAppending;
    int appending;
    int LengthBegin;
    int LengthEnd;
    int kpnumber;
    int RunningFlag=0;
    int kppos;
    int bandgapid;
    int AtomNumbers;
    int thrd;
    int midlen;
    int DataTag;
    int tmpName;
    int USE_parallel_K_points;
    int USE_parallel_output;
    int USE_CalcDOS;
    int USE_Infiniteloop;
    int USE_bulk_states_construction;
    int USE_parallel_in_preCalc;
    int USE_pzheevx_in_blocks_diag;
    int USE_Eigen_in_blocks_diag;
    int USE_blocks_transform;
    int USE_blocks_transform_bulk_states;
    int USE_transform_matrix;
    int USE_distributed;
    int USE_pzheevx_OR_Eigen_in_preCalc;

    int USE_LoMEType1;
    int USE_LoMEType2;
    int USE_LoMEType3;
    int USE_LoMC;


    double bandgap;
    double conductionedge,valenceedge;
    int conductionpos,valencepos;
    std::string kppath;
    std::string dataVersion;
    std::string inpfile;
    std::vector<int> EventTable;
    std::vector<int> MatrixBottom;
    std::vector<unsigned int> symmetryPointsPositions;
    std::vector<double> ConductionBand;
    std::vector<double> ValenceBand;
    std::vector<double> BandValue;
    std::complex<double> sumX;

    //TableOFactionsName[Init]="Init";
    std::unordered_map<std::string,double &> parseElementsDouble=
    {
	{"down",hamiltonian.leftCutoff},
	{"up",hamiltonian.rightCutoff},
	//{"MixCoeff",hamiltonian.MixCoeff},
	//{"USE_Shift",hamiltonian.USE_Shift},
	{"dis",hamiltonian.systemcell.dis}
	//{"BandOffSet",hamiltonian.BandOffSet},
	//{"scaleFactor",hamiltonian.scaleFactor}
    };
    std::unordered_map<std::string,double &> parseUSE_double=
    {
	{"MixCoeff",hamiltonian.MixCoeff},
	{"USE_Shift",hamiltonian.USE_Shift},
	{"BandOffSet",hamiltonian.BandOffSet},
	{"scaleFactor",hamiltonian.scaleFactor}
    };

    std::unordered_map<std::string,int &> parseUSE_int=
    {
	{"AutoMatic",hamiltonian.AutoMatic},
	{"USE_FixedSize_or_EnergyRange",hamiltonian.USE_FixedSize_or_EnergyRange},
	{"USE_SOC",hamiltonian.USE_SOC},
	{"USE_long_parameters",hamiltonian.USE_long_parameters},
	{"USE_CalcDOS",USE_CalcDOS},
	{"USE_Infiniteloop",USE_Infiniteloop},
	{"USE_MaterialMix",hamiltonian.USE_MaterialMix},
	{"USE_OneRange_or_MoreRanges",hamiltonian.USE_OneRange_or_MoreRanges},
	{"USE_parallel_in_preCalc",USE_parallel_in_preCalc},
	{"USE_pzheevx_OR_Eigen_in_preCalc",USE_pzheevx_OR_Eigen_in_preCalc},
	{"USE_parallel_output",USE_parallel_output},
	{"USE_bulk_states_construction",USE_bulk_states_construction},
	{"USE_direct_construction",hamiltonian.USE_direct_construction},
	{"USE_decomposition",hamiltonian.USE_decomposition},
	{"USE_RowTruncation",hamiltonian.USE_RowTruncation},
	{"USE_ColTruncation",hamiltonian.USE_ColTruncation},
	{"USE_parallel_K_points",USE_parallel_K_points},
	{"USE_pzheevx_in_blocks_diag",USE_pzheevx_in_blocks_diag},
	{"USE_Eigen_in_blocks_diag",USE_Eigen_in_blocks_diag},
	{"USE_blocks_transform",USE_blocks_transform},
	{"USE_blocks_transform_bulk_states",USE_blocks_transform_bulk_states},
	{"USE_transform_matrix",USE_transform_matrix},
	{"USE_distributed",USE_distributed},
	{"USE_extra_k_point",hamiltonian.USE_extra_k_point},
	{"USE_LoMEType1",USE_LoMEType1},
	{"USE_LoMEType2",USE_LoMEType2},
	{"USE_LoMEType3",USE_LoMEType3},
	{"USE_LoMC",USE_LoMC},
	{"USE_3Dpbc",hamiltonian.systemcell.USE_3Dpbc}
    };
    std::unordered_map<std::string,int &> parseElements=
    {
	{"bandbegin",bandbegin},
	{"bandend",bandend},
	{"LengthBegin",LengthBegin},
	{"LengthEnd",LengthEnd},
	{"kpnumber",kpnumber},
	{"DataTag",DataTag},
	//{"AutoMatic",hamiltonian.AutoMatic},
	//{"USE_FixedSize_or_EnergyRange",hamiltonian.USE_FixedSize_or_EnergyRange},
	//{"USE_SOC",hamiltonian.USE_SOC},
	//{"USE_OneRange_or_MoreRanges",hamiltonian.USE_OneRange_or_MoreRanges},
	//{"USE_parallel_in_preCalc",USE_parallel_in_preCalc},
	//{"USE_pzheevx_OR_Eigen_in_preCalc",USE_pzheevx_OR_Eigen_in_preCalc},
	//{"USE_parallel_output",USE_parallel_output},
	//{"USE_bulk_states_construction",USE_bulk_states_construction},
	//{"USE_direct_construction",hamiltonian.USE_direct_construction},
	//{"USE_decomposition",hamiltonian.USE_decomposition},
	//{"USE_RowTruncation",hamiltonian.USE_RowTruncation},
	//{"USE_ColTruncation",hamiltonian.USE_ColTruncation},
	//{"USE_parallel_K_points",USE_parallel_K_points},
	//{"USE_pzheevx_in_blocks_diag",USE_pzheevx_in_blocks_diag},
	//{"USE_Eigen_in_blocks_diag",USE_Eigen_in_blocks_diag},
	//{"USE_blocks_transform",USE_blocks_transform},
	//{"USE_blocks_transform_bulk_states",USE_blocks_transform_bulk_states},
	//{"USE_transform_matrix",USE_transform_matrix},
	//{"USE_distributed",USE_distributed},
	//{"USE_extra_k_point",hamiltonian.USE_extra_k_point},
	//{"USE_LoMEType1",USE_LoMEType1},
	//{"USE_LoMEType2",USE_LoMEType2},
	//{"USE_LoMEType3",USE_LoMEType3},
	//{"USE_LoMC",USE_LoMC},
	//{"USE_3Dpbc",hamiltonian.systemcell.USE_3Dpbc},
	{"DividedFactor",hamiltonian.systemcell.DividedFactor},
	{"AtomNum",hamiltonian.systemcell.UnitCellNumber}
    };
    std::unordered_map<std::string,std::string &> parseString=
    {
	{"kppath",kppath},
	{"dataVersion",dataVersion}
    };
    bool flagthrd;
    std::vector<std::vector<double>> results;

    double kkparam;
    void ChooseOptions(int &in);
    void alter_run_modes(int mode,std::vector<int> &alter_table);
    auto loop_sub_run(std::vector<int> &runlists,int level) -> void;
    void BandOutput(const Eigen::VectorXd &eigenvals);
    void BandOutputW(const Eigen::VectorXd &WW,int opTag1,int opTag2);
    BandStructure(){kkparam=(pow(h_bar*iA,2)/(m_e*e_c));outputflag=0;computeflag=0; bandbegin=0; bandend=0; LengthBegin=0; LengthEnd=0;};
    void Output(int No,int interval,int fileNo,std::string dataVersion,int range_ouput,int appending);
    void NewOutput(int No,int interval,int fileNo);
    void CalculateMass(const Eigen::MatrixXcd &VectM,const Eigen::VectorXd &eigenvals,Eigen::VectorXcd &Mmass);
    void CalculateMassByNumerical(const Eigen::Vector3d zz,const Eigen::Vector3d gr,double ShiftValue);
    void CalculateMassByNumerical(const Eigen::Vector3d zz,const Eigen::Vector3d gr,Eigen::VectorXcd &fNumMass);
    void FindBandEdge(const std::complex<double> kz,Eigen::Vector3d pl,const Eigen::MatrixXcd &VectM,const Eigen::VectorXd &eigenvals,int ShiftValue);
	    void FindBandEdge(const Eigen::VectorXd &eigenvals,Eigen::VectorXcd &Mmass,int ShiftValue);
    void PrintMass(const Eigen::VectorXd &eigenvals,Eigen::VectorXcd &Mmass,Eigen::VectorXcd &Mmass2,int leftSf,int rightSf);
    void GenerateChi();
    void SimpleOutput(int No,int interval,int fileNo,std::string dataVersion,int range_ouput,int appending);
    void MeshBrillouinZone();
    void OutPutStatus(int flag);
    void H_recycle();
    auto H_recycle(std::vector<int>&table,int &level,int opTag1,int opTag2) -> void;
    void H_recyclePlus3(std::vector<int>&table,int &level);
    void System_design();
    void Band_Init_k_path();
    void AutoMaticBandInit();
    void CustomBandInit();
    void BandInitEpsilon();
    void Initialize(std::vector<std::string> path,unsigned int nrPoints = 600);
    void H_store_results();
    void Generate_mesh(int q);
    void RandomKpointList(int q);
    void H_loop_k_points(int backwordNum);
    void H_loop_k_points(int backwordNum,int &level);
    void H_clean();
    void ComputePrepareMPI();
    void ComputePrepareDFTB();
    void ParseParameters();
    unsigned int GetPointsNumber() const { return static_cast<unsigned int>(kpoints.size()); }
    const std::vector<std::string>& GetPath() const { return m_path; }
    void define_run_modes(int calc_mode);
    auto set_run_modes(int calc_mode,std::vector<int> &table,int flag) -> void;
    void loop_run();
    auto master_loop_run(int mode,std::vector<int> &table,int level,int flag) -> void;
    void MakeDecision();
    void MakeDecision(std::vector<int>&table,int &level);
    void Output_Epsilon();

	private:
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> k1solver;
    Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> k2solver;
    int CurrentMode;
    int intpar;
    int nrPoints;
    int interval;
    int startPoint;
    int cendPoint;
    int outputflag;
    int computeflag;
    std::vector<std::string> m_path;
    std::vector<Eigen::Vector3d> kpoints;
    std::vector<Eigen::Vector3d> kpoints_mesh;
    std::vector<int> kpoints_degen;
    std::vector<Eigen::Vector3d> qpoints;
    int Monkhorstpack(int r,int q);
    };
}

#endif
