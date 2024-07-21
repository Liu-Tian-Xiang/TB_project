#include "BandStructure.hpp"
namespace po = boost::program_options;

Eigen::Vector3d RotateAround(const Eigen::Vector3d &Axis,const Eigen::Vector3d &ts,double angle)
{
    Eigen::Vector3d result;
	const Eigen::Vector3d u = Axis.normalized();
	double sinAngle = sin(angle);
	double cosAngle = cos(angle);
    double oneMinus = 1. - cosAngle;

	result(0) = (cosAngle + oneMinus * u(0) * u(0)) * ts(0) + (oneMinus * u(0) * u(1) - u(2) * sinAngle) * ts(1) + (oneMinus * u(0) * u(2) + u(1) * sinAngle) * ts(2);
	result(1) = (oneMinus * u(0) * u(1) + u(2) * sinAngle) * ts(0) + (cosAngle + oneMinus * u(1) * u(1)) * ts(1) + (oneMinus * u(1) * u(2) - u(0) * sinAngle) * ts(2);
	result(2) = (oneMinus * u(0) * u(2) - u(1) * sinAngle) * ts(0) + (oneMinus * u(1) * u(2) + u(0) * sinAngle) * ts(1) + (cosAngle + oneMinus * u(2) * u(2)) * ts(2);
    return result;
}
auto init_options(TightBinding::BandStructure &bandStructure) -> po::options_description
{
    po::options_description description;
    description.add_options()
        ( "help", "Produce the help message.\n" );
    description.add_options()
        ( "appending", "appending mode for output\n" );
    description.add_options()
        ( "run"
          , po::value<std::string>()//->required()
          , "Running Modes:"
          "\n22 --> 2D eigenvector transform 2D block Hamiltonian"
          "\n32 --> 3D eigenvector transform 2D block Hamiltonian"
          "\n3  --> Almost the same as 32 but including the folding algorithm"
          "\n33 --> 3D eigenvector transform 3D block Hamiltonian\n"
        );
    description.add_options()
        ( "compare"
          , po::value<std::string>()//->required()
          , "Compare Running Modes:"
          " ./ecode --compare 22,33,32"
        );
    description.add_options()
        ( "tmpName"
          , po::value<int>()//->required()
          , "Output tmpName:"
          " ./ecode --tmpName 02"
        );
/*
    description.add_options()
        ( "MixCoeff"
          , po::value<double>()//->required()
          , "Set MixCoeff:"
          " ./ecode --MixCoeff 0.2"
        );
    description.add_options()
        ( "USE_Shift"
          , po::value<double>()//->required()
          , "Set Shift:"
          " ./ecode --USE_Shift 0.2"
        );
    description.add_options()
        ( "BandOffSet"
          , po::value<double>()//->required()
          , "Set BandOffSet:"
          " ./ecode --BandOffSet 0.2"
        );
*/
    description.add_options()//Intparse
        ( "dataTag"
          , po::value<int>()//->required()
          , "Output dataTag:"
          " ./ecode --dataTag 02"
        );
    description.add_options()
        ( "dataFlag"
          , po::value<std::string>()//->required()
          , "Output dataFlag:"
          " ./ecode --dataFlag 02"
        );
    description.add_options()
        ( "inpfile"
          , po::value<std::string>()//->required()
          , "Output inpfile:"
          " ./ecode --inpfile 02"
        );
    description.add_options()
        ( "range"
          , po::value<std::string>()//->required()
          , "Set energy range:"
          " ./ecode --range -2,10"
        );
    description.add_options()
        ( "StructureSettings"
          , po::value<std::string>()//->required()
          , "Setting structure index:"
          " ./ecode --StructureSettings 0(bandbegin),3(bandend),4(DividedFactor)"
        );
    description.add_options()
        ( "LengthSettings"
          , po::value<std::string>()//->required()
          , "Setting quantum well wide index:"
          " ./ecode --LengthSettings 0(LengthBegin),3(LengthEnd),4(DividedFactor)"
        );
        //;
/*
    std::vector<std::string> USE_String=
    {
        {"AutoMatic"},
        {"USE_FixedSize_or_EnergyRange"},
        {"USE_SOC"},
        {"USE_MaterialMix"},
        {"USE_OneRange_or_MoreRanges"},
        {"USE_parallel_in_preCalc"},
        {"USE_pzheevx_OR_Eigen_in_preCalc"},
        {"USE_parallel_output"},
        {"USE_bulk_states_construction"},
        {"USE_direct_construction"},
        {"USE_decomposition"},
        {"USE_RowTruncation"},
        {"USE_ColTruncation"},
        {"USE_parallel_K_points"},
        {"USE_pzheevx_in_blocks_diag"},
        {"USE_Eigen_in_blocks_diag"},
        {"USE_blocks_transform"},
        {"USE_blocks_transform_bulk_states"},
        {"USE_transform_matrix"},
        {"USE_distributed"},
        {"USE_extra_k_point"},
        {"USE_LoMEType1"},
        {"USE_LoMEType2"},
        {"USE_LoMEType3"},
        {"USE_LoMC"},
        {"USE_3Dpbc"},
    };
    for(auto &str : USE_String)
    {
        description.add_options()
            ( str.c_str()
              , po::value<int>()//->required()
              , str.c_str()
            );
    }
*/

    for(auto &str : bandStructure.parseUSE_int)
    {
        description.add_options()
            ( str.first.c_str()
              , po::value<int>()//->required()
              , str.first.c_str()
            );
    }

    for(auto &str : bandStructure.parseUSE_double)
    {
        description.add_options()
            ( str.first.c_str()
              , po::value<double>()//->required()
              , str.first.c_str()
            );
    }

    return description;
}


template<class _Help, class _Run>
auto process_command_line( int argc, char** argv
        , _Help&& help
        , _Run&& run 
        , _Run&& compare
        ) -> void
{
    TightBinding::BandStructure bandStructure;
    auto const description = init_options(bandStructure);
    po::variables_map vm;

    po::store( po::command_line_parser(argc, argv)
            .options(description)
            .run()
            , vm );
    if(vm.count("help")) {
        po::notify(vm);
        help(description);
        return;
    }
    if(vm.count("compare")) {
        po::notify(vm);
        //compare(vm);
        return;
    }
    if(vm.count("run")) {
        po::notify(vm);
        run(vm,bandStructure);
        return;
    }
}

template <class _R>
auto parse_run_mode(std::string const& str) -> std::vector<_R>
{

    auto const parse_number = [](auto const& s) {
        return boost::lexical_cast<_R>(boost::trim_copy(s));
    };

    std::vector<std::string> tokens=split<char const >(str,',');


    std::vector<_R> qs;
    qs.reserve(tokens.size());
    copy( tokens | boost::adaptors::transformed(std::cref(parse_number)) 
            , std::back_inserter(qs) );

    return qs;

}

void treate_vm(TightBinding::BandStructure &bandStructure,boost::program_options::variables_map const& vm)
{
    int id;
    MPI_Comm_rank(MPI_COMM_WORLD,&id);

   std::fstream file_obj; 
    
    ///vm treating
    bandStructure.DataAppending=0;
    bandStructure.FirstLoop=true;
    if(vm.count("appending")){
	bandStructure.DataAppending=1;
    }
    if(vm.count("tmpName")){
	bandStructure.tmpName=vm["tmpName"].as<int>();
    }else{
	bandStructure.tmpName=0;
    }
    if(id==0)
    {
	file_obj.open(("./out"+std::to_string(bandStructure.tmpName)+".txt").c_str(),std::ios::out|std::ios::trunc);
        file_obj.close();
        file_obj.open(strdup("./out.code"),std::ios::out|std::ios::trunc);
        file_obj.close();
    }

/*
    if(vm.count("MixCoeff")){
	bandStructure.hamiltonian.MixCoeff=vm["MixCoeff"].as<double>();
	bandStructure.parseElementsDouble.erase("MixCoeff");
	if(id==0)
	{
	    std::fstream file_obj; 
	    file_obj.open("./out"+std::to_string(bandStructure.tmpName)+".txt",std::ios::out|std::ios::app);
	    file_obj<<"MixCoeff"<<" = "<<bandStructure.hamiltonian.MixCoeff<<std::endl;
	    file_obj.close();
	}
    }

    if(vm.count("USE_Shift")){
        bandStructure.hamiltonian.USE_Shift=vm["USE_Shift"].as<double>();
        bandStructure.parseElementsDouble.erase("USE_Shift");
        if(id==0)
        {
            std::fstream file_obj; 
            file_obj.open("./out"+std::to_string(bandStructure.tmpName)+".txt",std::ios::out|std::ios::app);
            file_obj<<"USE_Shift"<<" = "<<bandStructure.hamiltonian.USE_Shift<<std::endl;
            file_obj.close();
        }
    }

    if(vm.count("BandOffSet")){
	bandStructure.hamiltonian.BandOffSet=vm["BandOffSet"].as<double>();
	bandStructure.parseElementsDouble.erase("BandOffSet");
	if(id==0)
	{
	    std::fstream file_obj; 
	    file_obj.open("./out"+std::to_string(bandStructure.tmpName)+".txt",std::ios::out|std::ios::app);
	    file_obj<<"BandOffSet"<<" = "<<bandStructure.hamiltonian.BandOffSet<<std::endl;
	    file_obj.close();
	}
    }
*/
    if(vm.count("dataTag")){
	bandStructure.DataTag=vm["dataTag"].as<int>();
	bandStructure.parseElements.erase("DataTag");
	if(id==0)
	{
	    std::fstream file_obj; 
	    file_obj.open("./out"+std::to_string(bandStructure.tmpName)+".txt",std::ios::out|std::ios::app);
	    file_obj<<"DataTag"<<" = "<<bandStructure.DataTag<<std::endl;
	    file_obj.close();
	}
    }

    std::vector<std::string> vecStr;
    vecStr.resize(0);
    for(auto &USEmap : bandStructure.parseUSE_int)
    {
        if(vm.count(USEmap.first)){
            USEmap.second=vm[USEmap.first].as<int>();
            //std::cout<<USEmap.first<<"="<<USEmap.second<<std::endl;
            if(id==0)
            {
                std::fstream file_obj; 
                file_obj.open("./out"+std::to_string(bandStructure.tmpName)+".txt",std::ios::out|std::ios::app);
                file_obj<<USEmap.first<<" = "<<USEmap.second<<std::endl;
                file_obj.close();
            }
            vecStr.push_back(USEmap.first);
	}else{
	    USEmap.second=0;
	}
    }
    for(auto &erStr:vecStr)
    {
        bandStructure.parseUSE_int.erase(erStr);
    }
    vecStr.resize(0);
/////double
    for(auto &USEmap : bandStructure.parseUSE_double)
    {
        if(vm.count(USEmap.first)){
            USEmap.second=vm[USEmap.first].as<double>();
            if(id==0)
            {
                std::fstream file_obj; 
                file_obj.open("./out"+std::to_string(bandStructure.tmpName)+".txt",std::ios::out|std::ios::app);
                file_obj<<USEmap.first<<" = "<<USEmap.second<<std::endl;
                file_obj.close();
            }
            vecStr.push_back(USEmap.first);
	}else{
            USEmap.second=0;
	}
    }
    for(auto &erStr:vecStr)
    {
        bandStructure.parseUSE_double.erase(erStr);
    }
    vecStr.resize(0);


    if(vm.count("inpfile")){
        bandStructure.inpfile=vm["inpfile"].as<std::string>();
        if(id==0)
        {
            std::fstream file_obj; 
            file_obj.open("./out"+std::to_string(bandStructure.tmpName)+".txt",std::ios::out|std::ios::app);
            file_obj<<"inpfile"<<" = "<<bandStructure.inpfile<<std::endl;
            file_obj.close();
        }
    }else{
        bandStructure.inpfile="inp";
        if(id==0)
        {
            std::fstream file_obj; 
            file_obj.open("./out"+std::to_string(bandStructure.tmpName)+".txt",std::ios::out|std::ios::app);
            file_obj<<"inpfile"<<" = "<<bandStructure.inpfile<<std::endl;
            file_obj.close();
        }
    }

    if(vm.count("dataFlag")){
        bandStructure.dataVersion=vm["dataFlag"].as<std::string>();
        bandStructure.parseString.erase("dataVersion");
    if(id==0)
    {
        std::fstream file_obj; 
        file_obj.open("./out"+std::to_string(bandStructure.tmpName)+".txt",std::ios::out|std::ios::app);
        file_obj<<"dataVersion"<<" = "<<bandStructure.dataVersion<<std::endl;
        file_obj.close();
    }
    }
    if(vm.count("range")){
        auto Erange = parse_run_mode<double>(vm["range"].as<std::string>());
        if(Erange.size()>0) {bandStructure.hamiltonian.leftCutoff=Erange[0];bandStructure.parseElementsDouble.erase("down");}
        if(Erange.size()>1) {bandStructure.hamiltonian.rightCutoff=Erange[1];bandStructure.parseElementsDouble.erase("up");}
    if(id==0)
    {
        std::fstream file_obj; 
        file_obj.open("./out"+std::to_string(bandStructure.tmpName)+".txt",std::ios::out|std::ios::app);
        if(Erange.size()>0) file_obj<<"down"<<" = "<<Erange[0]<<std::endl;
        if(Erange.size()>1) file_obj<<"up"<<" = "<<Erange[1]<<std::endl;
        file_obj.close();
    }
    }
    if(vm.count("StructureSettings")){
        auto Ssettings = 
            parse_run_mode<double>(vm["StructureSettings"].as<std::string>());

        if(Ssettings.size()>0) {bandStructure.bandbegin=Ssettings[0];bandStructure.parseElements.erase("bandbegin");}
        if(Ssettings.size()>1) {bandStructure.bandend=Ssettings[1];bandStructure.parseElements.erase("bandend");}
        if(Ssettings.size()>2) {bandStructure.hamiltonian.systemcell.DividedFactor=Ssettings[2];bandStructure.parseElements.erase("DividedFactor");}

	if(bandStructure.bandbegin>pow(2,bandStructure.hamiltonian.systemcell.DividedFactor-1)){bandStructure.bandbegin=0;}
	if(bandStructure.bandend>pow(2,bandStructure.hamiltonian.systemcell.DividedFactor-1)){bandStructure.bandend=0;}
        
    if(id==0)
    {
        std::fstream file_obj; 
        file_obj.open("./out"+std::to_string(bandStructure.tmpName)+".txt",std::ios::out|std::ios::app);
        if(Ssettings.size()>0) file_obj<<"bandbegin"<<" = "<<Ssettings[0]<<std::endl;
        if(Ssettings.size()>1) file_obj<<"bandend"<<" = "<<Ssettings[1]<<std::endl;
        if(Ssettings.size()>2) file_obj<<"DividedFactor"<<" = "<<Ssettings[2]<<std::endl;
        file_obj.close();
    }
    }
    if(vm.count("LengthSettings")){
        auto Lsettings = 
            parse_run_mode<double>(vm["LengthSettings"].as<std::string>());

        if(Lsettings.size()>0) {bandStructure.LengthBegin=Lsettings[0];bandStructure.parseElements.erase("LengthBegin");}
        if(Lsettings.size()>1) {bandStructure.LengthEnd=Lsettings[1];bandStructure.parseElements.erase("LengthEnd");}
        if(Lsettings.size()>2) {bandStructure.hamiltonian.systemcell.DividedFactor=Lsettings[2];bandStructure.parseElements.erase("DividedFactor");}

	if(bandStructure.LengthBegin>pow(2,bandStructure.hamiltonian.systemcell.DividedFactor-1)){bandStructure.LengthBegin=0;}
	if(bandStructure.LengthEnd>pow(2,bandStructure.hamiltonian.systemcell.DividedFactor-1)){bandStructure.LengthEnd=0;}
        
    if(id==0)
    {
        std::fstream file_obj; 
        file_obj.open("./out"+std::to_string(bandStructure.tmpName)+".txt",std::ios::out|std::ios::app);
        if(Lsettings.size()>0) file_obj<<"LengthBegin"<<" = "<<Lsettings[0]<<std::endl;
        if(Lsettings.size()>1) file_obj<<"LengthEnd"<<" = "<<Lsettings[1]<<std::endl;
        if(Lsettings.size()>2) file_obj<<"DividedFactor"<<" = "<<Lsettings[2]<<std::endl;
        file_obj.close();
    }
    }


    bandStructure.ParseParameters();

}
/*
template <class _C>
auto compare(boost::program_options::variables_map const& vm) -> void
{
    using namespace boost;
    auto RunModes = 
        parse_run_mode<_C>(vm["compare"].as<std::string>());

    int id;
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD,&id);
    

    TightBinding::BandStructure bandStructure;
    treate_vm(bandStructure,vm);
    std::fstream file_obj; 
    if(id==0)
    {
        file_obj.open(strdup("./out.code"),std::ios::out|std::ios::app);
        file_obj<<"RunModes =";
        for(int i=0;i<RunModes.size();++i)
        {
            file_obj<<" "<<RunModes[i];
        }
        file_obj<<std::endl;
        file_obj.close();
    }
    for(int i=0;i<RunModes.size();++i)
    {
        file_obj.open(("./BandData/"+bandStructure.dataVersion+"0/bandGap"+std::to_string(RunModes[i])+".data").c_str(),std::ios::out|std::ios::trunc);
        file_obj.close();
        file_obj.open(("./BandData/"+bandStructure.dataVersion+"0/TimeSpend"+std::to_string(RunModes[i])+".data").c_str(),std::ios::out|std::ios::trunc);
        file_obj.close();

        bandStructure.define_run_modes(RunModes[i]);
        bandStructure.loop_run();
    }

    MPI_Finalize();
}
*/
template <class _C>
auto run(boost::program_options::variables_map const& vm,TightBinding::BandStructure &bandStructure) -> void
{
    using namespace boost;
    //using _R = typename _C::value_type;

    auto RunMode = boost::lexical_cast<_C>(vm["run"].as<std::string>());

    int id;
    MPI_Init(NULL,NULL);
    MPI_Comm_rank(MPI_COMM_WORLD,&id);
    
    //TightBinding::BandStructure bandStructure;

    std::fstream file_obj; 
    treate_vm(bandStructure,vm);
    if(id==0)
    {
        file_obj.open(strdup("./out.code"),std::ios::out|std::ios::app);
        file_obj<<"RunModes =";
        file_obj<<" "<<RunMode;
        file_obj<<std::endl;
        file_obj.close();
    }

    
    bandStructure.set_run_modes(RunMode,bandStructure.EventTable,1);
    bandStructure.master_loop_run(1,bandStructure.EventTable,bandStructure.EventTable.size()-1,1);

    MPI_Finalize();
}

auto dispatch(boost::program_options::variables_map const& vm,TightBinding::BandStructure &bandStructure) -> void
{
    run<int>(vm,bandStructure);
}

auto dispatchCompare(boost::program_options::variables_map const& vm,TightBinding::BandStructure &bandStructure) -> void
{
    //compare<int>(vm);
}


int main(int argc,char *argv[])
{
    process_command_line
        ( argc, argv
          , [](auto desc) { std::cout << desc << '\n'; }
          , &dispatch
          , &dispatchCompare
        );
    return 0;
}

