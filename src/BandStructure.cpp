#include "BandStructure.hpp"
namespace TightBinding
{

    void BandStructure::GenerateChi()
    {
        int numprocs;
        MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
        unsigned int ic=startPoint+computeflag;

        //Eigen::Vector3d q=kpoints_mesh[i];///qqq
        //double omega=0;
        
        Eigen::Vector3d q=Eigen::Vector3d(0,0,0);///qqq
        double omega=kpoints_mesh[outputflag](0);

        Eigen::Vector3d k2=kpoints[ic];
        Eigen::Vector3d k1=k2-q;
        Eigen::MatrixXcd MeshEigenvectors;
        Eigen::VectorXd MeshEigenvalues;

        eigenvectors=hamiltonian.transform_matrix*eigenvectors;

        //set 1
        if(k1!=k2)
        {
            hamiltonian.SetMatrix(k1,hamiltonian.systemcell.ListofMatrixElements_Type1,hamiltonian.matrix);
            k1solver.compute(hamiltonian.matrix,Eigen::MatrixXcd::Identity(hamiltonian.matrix.rows(),hamiltonian.matrix.cols()), Eigen::DecompositionOptions::ComputeEigenvectors |Eigen::Ax_lBx);
            assert(k1solver.info() == Eigen::ComputationInfo::Success);
            MeshEigenvectors=k1solver.eigenvectors();
            MeshEigenvalues=k1solver.eigenvalues();
        }else{
            MeshEigenvectors=eigenvectors;
            MeshEigenvalues=eigenvalues;
        }

        
        //hamiltonian.store_precomputed_data(k2,ic,eigenvalues,eigenvectors);
        //hamiltonian.store_precomputed_data(k2,-ic,MeshEigenvalues,MeshEigenvectors);
        
        //set 2
        //hamiltonian.SetMatrix(k2,hamiltonian.systemcell.ListofMatrixElements_Type1,hamiltonian.matrix);
        //k2solver.compute(hamiltonian.matrix,Eigen::MatrixXcd::Identity(hamiltonian.matrix.rows(),hamiltonian.matrix.cols()), Eigen::DecompositionOptions::ComputeEigenvectors |Eigen::Ax_lBx);
        //assert(k2solver.info() == Eigen::ComputationInfo::Success);

        std::complex<double> ii(0.0,1.0);
        int dimX=hamiltonian.MatrixSize;
        Eigen::MatrixXcd G;
        G.resize(hamiltonian.MatrixSizeReduced,hamiltonian.MatrixSizeReduced);
        G.setZero();
        //for(int i=0;i<dimX;++i)
        for(int i=0;i<hamiltonian.MatrixSizeReduced;++i)
        {
            //for(int j=0;j<dimX;++j)
            for(int j=0;j<hamiltonian.MatrixSizeReduced;++j)
            {
                double Ei=MeshEigenvalues(i);//k1solver.eigenvalues()(i);
                double Ej=eigenvalues(j);
                //G(i,j)=((i<AtomNumbers*hamiltonian.BandGapIndex?1.0:0.0)-(j<AtomNumbers*hamiltonian.BandGapIndex?1.0:0.0))/(Ei-Ej-omega-0.05*ii);
                //if(abs(omega-Ei)<10 && abs(omega-Ej)<10)
                //if(abs(Ej-Ei)<10)
                //if((Ei<5 && Ei>-1) && (Ej<5 && Ej>-1) || (abs(omega-Ei)<3 && abs(omega-Ej)<3))
		//if((abs(omega-Ei)<=3 && abs(omega-Ej<=3)) || (abs(Ei-Ej)<6))
		if((abs(Ei-Ej)<7))
		//if(1)
        /*
		if(
		(MatrixBottom[ic]+i)<AtomNumbers*hamiltonian.BandGapIndex*1.62 && 
		(MatrixBottom[ic]+j)<AtomNumbers*hamiltonian.BandGapIndex*1.62 //&& 
		//(MatrixBottom[ic]+i)>AtomNumbers*hamiltonian.BandGapIndex/6.0 && 
		//(MatrixBottom[ic]+j)>AtomNumbers*hamiltonian.BandGapIndex/6.0 
		  )
*/
		{
                G(i,j)=
                    (
                     (((MatrixBottom[ic]+i)<AtomNumbers*hamiltonian.BandGapIndex)?1.0:0.0)
                     -(((MatrixBottom[ic]+j)<AtomNumbers*hamiltonian.BandGapIndex)?1.0:0.0)
                     )/(Ei-Ej-omega-0.08*ii);
        //std::cout<<"iF="<<integrate_factor<<std::endl;
		}
                //integrate_factor+=1;
            }
        }
        for(int alpha=0;alpha<dimX;++alpha)
        {
            for(int alprime=0;alprime<dimX;++alprime)
            {
                Eigen::VectorXcd A;
                Eigen::VectorXcd B;
                A.resize(MeshEigenvectors.cols());
                A.setZero();
                B.resize(eigenvectors.cols());
                B.setZero();
                //for(int i=0;i<dimX;++i)
                for(int i=0;i<hamiltonian.MatrixSizeReduced;++i)
                {

                    A(i)= MeshEigenvectors(alpha,i)*
                        MeshEigenvectors.conjugate()(alprime,i);
                    B(i)= eigenvectors(alpha,i)*
                        eigenvectors.conjugate()(alprime,i);

                    /*
                       A(i)= k1solver.eigenvectors()(alpha,i)*
                       k1solver.eigenvectors().adjoint()(alprime,i);
                       B(i)= k2solver.eigenvectors()(alpha,i)*
                       k2solver.eigenvectors().adjoint()(alprime,i);
                       */
                }
                //auto const tmp=(G.transpose()*A).transpose()*B.conjugate();
                //Xmatrix(alpha,alprime)+=tmp(0);
                Xmatrix(alpha,alprime)+=(A.transpose()*G*B)(0);
            }
        }


    }

    void BandStructure::MeshBrillouinZone()
    {
        int id,numprocs;
        MPI_Comm_rank(MPI_COMM_WORLD,&id);
        MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
        AtomNumbers=0;
        std::for_each(hamiltonian.systemcell.supercell.begin(),hamiltonian.systemcell.supercell.end(),
                [&atom=AtomNumbers]
                (auto &scell)
                {atom+=scell->AtomNum;}
                );

        Generate_mesh(8);
        //kpoints=kpoints_mesh;
        int dimX=hamiltonian.MatrixSize;
        Xmatrix.resize(dimX,dimX);
        Xmatrix.setZero();
        nrPoints = kpoints.size();//GetPointsNumber();
        interval=nrPoints;
        startPoint = 0;
        cendPoint=nrPoints;
        MatrixBottom.clear();
        if(id==0 && FirstLoop)
        {

            std::fstream file_obj; 
            file_obj.open(("./out"+std::to_string(tmpName)+".txt").c_str(),std::ios::out|std::ios::app);
            file_obj.close();
            //char dirname[40];    
            //sprintf(dirname,"./Data_%1.2f-%1.2f",hamiltonian.leftCutoff,hamiltonian.rightCutoff);
            //mkdir(dirname,S_IRWXU|S_IRWXG|S_IRWXO);

            rename(("./out"+std::to_string(tmpName)+".txt").c_str(),("./BandData/"+dataVersion+"0/out"+std::to_string(DataTag)+".txt").c_str());//DOS
        }

        if(FirstLoop) FirstLoop=false;
    }

    void BandStructure::ChooseOptions(int &in)
    {

        if(0)
        {

        }else{
            if(USE_distributed)
            {
                if(USE_blocks_transform)
                {
                    in = apply_method3;
                }else{

                }
                if(USE_bulk_states_construction)
                {
                    in = apply_method9;
                }
            }else{
                if(USE_blocks_transform)
                {
                    in = apply_method;
                    if(USE_blocks_transform_bulk_states)
                    {
                        if(hamiltonian.USE_extra_k_point){
                            in = apply_method7;
                        }else{
                            in = apply_method5;
                        }
                    }
                    if(USE_bulk_states_construction)
                    {
                        in = apply_method8;
                    }
                    if(USE_parallel_in_preCalc)
                    {
                        in = apply_method6;
                    }

                }
                if(USE_transform_matrix)
                {
                    in = method_one_transform;
                }
            }

        }
    }
    void BandStructure::BandOutputW(const Eigen::VectorXd &WW,int opTag1,int opTag2)
    {
        int id,numprocs;
        MPI_Comm_rank(MPI_COMM_WORLD,&id);
        MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
        int i=startPoint+computeflag;
        const Eigen::VectorXd &eigenvals=WW;


        if(conductionpos>0 && conductionpos<eigenvals.size())
        {
            conductionedge=eigenvals(conductionpos);
            ConductionBand.push_back(conductionedge);
        }else{
            std::cout<<"Energy range(conduction) too narrow -999!"<<std::endl;
            conductionedge=-999;
            ConductionBand.push_back(-999);
        }

        if(valencepos>0 && valencepos<eigenvals.size())
        {
            valenceedge=eigenvals(valencepos);
            ValenceBand.push_back(valenceedge);
        }else{
            std::cout<<"Energy range(valence) too narrow -999!"<<std::endl;
            valenceedge=-999;
            ValenceBand.push_back(-999);
        }
        bandgap=conductionedge-valenceedge;
        BandValue.push_back(bandgap);



        char name[40];    

        if(id==0 || USE_parallel_output)
        {
            std::fstream file_obj; 
            file_obj.open(("./BandData/"+dataVersion+std::to_string(opTag1)+"_"+std::to_string(opTag2)+"/M"+std::to_string(CurrentMode)+"band"+std::to_string(id)+".data").c_str(),std::ios::out|std::ios::app);
            file_obj<<i<<"\t";
            for(int ii=-hamiltonian.matrix_bottom;ii<=valencepos;++ii)
            {
                if(ii>=0 && ii<eigenvals.size())
                {
                    file_obj<<std::fixed<<std::setprecision(10)<<eigenvals(ii)<<"\t";
                }else{
                    if(hamiltonian.leftCutoff*hamiltonian.rightCutoff<0 && USE_CalcDOS!=1)
                    {
                        file_obj<<std::fixed<<std::setprecision(10)<<hamiltonian.leftCutoff<<"\t";
                        //file_obj<<std::fixed<<std::setprecision(10)<<-999<<"\t";
                    }else{
                        file_obj<<std::fixed<<std::setprecision(10)<<-999<<"\t";
                    }
                }
            }
            for(int ii=conductionpos;ii<hamiltonian.MatrixSize-hamiltonian.matrix_bottom;++ii)
            {
                if(ii>=0 && ii<eigenvals.size())
                {
                    file_obj<<std::fixed<<std::setprecision(10)<<eigenvals(ii)<<"\t";
                }else{
                    if(hamiltonian.leftCutoff*hamiltonian.rightCutoff<0&& USE_CalcDOS!=1)
                    {
                        file_obj<<std::fixed<<std::setprecision(10)<<hamiltonian.rightCutoff<<"\t";
                        //file_obj<<std::fixed<<std::setprecision(10)<<999<<"\t";
                    }else{
                        file_obj<<std::fixed<<std::setprecision(10)<<999<<"\t";
                    }
                }
            }
            file_obj<<std::endl;
            file_obj.close();
        }
    }

    auto BandStructure::set_run_modes(int calc_mode,std::vector<int> &table,int flag) -> void
    {
        if(flag) CurrentMode=calc_mode;
        switch(calc_mode)
        {
            case method_one_transform:
                table.clear();
                table.push_back(End);
                table.push_back(transform_Hamiltonian);
                table.push_back(blocks_diagonalization);
                break;
            case apply_method6:
                table.clear();
                table.push_back(End);
                table.push_back(reconstruct_Hamiltonian3);
                table.push_back(read_elements);
                table.push_back(blocks_diagonalization3);
                break;
            case apply_method7:
                table.clear();
                table.push_back(End);
                table.push_back(reconstruct_Hamiltonian4);
                table.push_back(read_elements);
                table.push_back(blocks_diagonalization);
                table.push_back(With_ek3D);
                break;
            case apply_method9:
                table.clear();
                table.push_back(End);
                table.push_back(reconstruct_Hamiltonian6);
                table.push_back(blocks_diagonalization9);
                break;
            case apply_method8:
                table.clear();
                table.push_back(End);
                table.push_back(reconstruct_Hamiltonian5);
                table.push_back(read_elements);
                table.push_back(blocks_diagonalization8);
                break;
            case apply_method5:
                table.clear();
                table.push_back(End);
                table.push_back(reconstruct_Hamiltonian2);
                table.push_back(read_elements);
                table.push_back(blocks_diagonalization);
                break;
            case apply_method:
                table.clear();
                table.push_back(End);
                table.push_back(reconstruct_Hamiltonian);
                table.push_back(read_elements);
                table.push_back(blocks_diagonalization);
                break;
            case apply_method3:
                table.clear();
                table.push_back(End);
                table.push_back(reconstruct_Hamiltonian_mpi);
                table.push_back(read_elements3);
                table.push_back(blocks_diagonalization3);
                break;
            case do_nothing:
                table.clear();
                table.push_back(End);
                table.push_back(still);
                break;
            case Running_pzheevx:
                //Running mode 222
                table.clear();
                table.push_back(End);
                table.push_back(Output_and_Recycle);
                table.push_back(COMPUTE_K_POINTS);
                table.push_back(output_data);
                table.push_back(diagonalize);
                table.push_back(pre_calculation);
                table.push_back(set_Hamiltonian);
                table.push_back(design_running_pattern);
                table.push_back(reshape_structure);
                table.push_back(find_neighbors);
                table.push_back(create_structure);
                table.push_back(Init);
                break;
            case Running_2D2D:
                //Running mode 22
                table.clear();
                table.push_back(End);
                table.push_back(Output_and_Recycle);
                table.push_back(COMPUTE_K_POINTS);
                table.push_back(output_data);
                table.push_back(diagonalize);
                table.push_back(pre_calculation);
                table.push_back(set_Hamiltonian);
                table.push_back(design_running_pattern);
                table.push_back(reshape_structure);
                table.push_back(find_neighbors);
                table.push_back(create_structure);
                table.push_back(Init);
                break;
            case Running_Epsilon:
                //Running mode 22
                table.clear();
                table.push_back(End);

                table.push_back(OutputEpsilon);//Output_Epsilon
                table.push_back(COMPUTE_K_POINTS);
                table.push_back(COMPUTE_Epsilon);//GenerateChi
                table.push_back(diagonalize);
                table.push_back(pre_calculation);
                table.push_back(set_Hamiltonian);
                table.push_back(MeshBZ);//MeshBrillouinZone
                table.push_back(reshape_structure);
                table.push_back(find_neighbors);
                table.push_back(create_structure);

                table.push_back(Init);//Band_Init_k_path
                break;
            case Running_5:
                //Running mode 0
                table.clear();
                table.push_back(End);
                table.push_back(COMPUTE_K_POINTS);
                table.push_back(COMPUTE_5);
                table.push_back(set_Hamiltonian);
                table.push_back(design_running_pattern);
                table.push_back(reshape_structure);
                table.push_back(find_neighbors);
                table.push_back(create_structure);
                table.push_back(Init);
                break;
                //default:
        }
        //RunningFlag=table.size()-1;

    }

    void BandStructure::MakeDecision(std::vector<int>&table,int &level)
    {

        int id,numprocs;
        MPI_Comm_rank(MPI_COMM_WORLD,&id);
        if(
                table[level] == COMPUTE_Epsilon
          )
        {
            //int posNum=std::distance(std::find(table.begin(),table.end(),table[level]),std::find(table.begin(),table.end(),COMPUTE_Epsilon))+1;
            
            //H_loop_k_points(1,level);//+1
        }

        if(
                table[level] == COMPUTE_K_POINTS
          )
        {
            int posNum=std::distance(std::find(table.begin(),table.end(),table[level]),std::find(table.begin(),table.end(),set_Hamiltonian))+1;
            H_loop_k_points(posNum,level);//+1
        }
        if(table[level] == OutputEpsilon)
        {
            //int posNum=std::distance(std::find(table.begin(),table.end(),table[level]),std::find(table.begin(),table.end(),Init_Epsilon))+1;
            //H_loop_k_points(posNum,level);
            H_recyclePlus3(table,level);//+3
        }

        if(table[level] == Output_and_Recycle)
        {
            if(hamiltonian.AutoMatic<=0)
            {
                NewOutput(id,startPoint,DataTag);
                H_recycle(table,level,DataTag,midlen);//+3
            }else{
                NewOutput(id,startPoint,thrd);
                H_recycle(table,level,thrd,midlen);//+3
            }
        }

        if(table[level] == pre_calculation)
        {
            std::vector<int> alter_table;
            alter_table.resize(0);
            int in=do_nothing;
            ChooseOptions(in);
            set_run_modes(in,alter_table,0);
            master_loop_run(1,alter_table,alter_table.size()-1,0);
        }

        if(table[level] == diagonalize)
        {
            std::vector<int> alter_table;
            alter_table.resize(0);
            alter_table.clear();
            alter_table.push_back(End);

            if(USE_distributed)
            {
                if(USE_bulk_states_construction)
                {
                    //alter_table.push_back(diagonalize_Eigen);
                    alter_table.push_back(diagonalize_pzheevx);
                }else{
                    alter_table.push_back(diagonalize_pzheevx);
                }
            }else if(USE_distributed==0){
                if(hamiltonian.USE_extra_k_point){
                    alter_table.push_back(diagonalize_Mix);
                }else{
                    alter_table.push_back(diagonalize_Eigen);
                }
            }

            master_loop_run(1,alter_table,alter_table.size()-1,0);
        }


        if(table[level] == set_Hamiltonian)
        {
            std::vector<int> alter_table;
            alter_table.resize(0);
            alter_table.clear();
            alter_table.push_back(End);

            if(USE_distributed && USE_blocks_transform)
            {
                if(USE_bulk_states_construction)
                {
                    alter_table.push_back(set_Hamiltonian_distributed);
                    //alter_table.push_back(set_Hamiltonian_sparse2);
                }else{
                    alter_table.push_back(set_Hamiltonian_distributed_with_shift);
                    alter_table.push_back(set_Hamiltonian_sparse);
                }
            }else if(USE_distributed){
                alter_table.push_back(set_Hamiltonian_distributed);
            }

            if(USE_LoMEType1)
            {
                /*
                   if(USE_distributed==1 && USE_bulk_states_construction==0){

                   }else{
                   alter_table.push_back(set_Hamiltonian_Type1);
                   }
                   */
                if(USE_distributed==1){

                }else{
                    alter_table.push_back(set_Hamiltonian_Type1);
                }
            }
            if(USE_LoMEType2)
            {
                alter_table.push_back(set_Hamiltonian_Type2);
            }
            if(USE_LoMEType3)
            {
                if(USE_distributed==1){
                    //alter_table.push_back(set_Hamiltonian_Type3);
                }else{
                    alter_table.push_back(set_Hamiltonian_Type3);
                }
            }
            if(hamiltonian.USE_MaterialMix)
            {
                //alter_table.push_back(set_Hamiltonian_TypeVec);
            }

            master_loop_run(1,alter_table,alter_table.size()-1,0);
        }

        if(table[level] == Init)
        {
            std::vector<int> alter_table;
            alter_table.resize(0);
            alter_table.clear();
            alter_table.push_back(End);

            if(hamiltonian.AutoMatic!=0)
            {
                alter_table.push_back(AutoMaticInit);
            }else{
                alter_table.push_back(CustomInit);
            }

            master_loop_run(1,alter_table,alter_table.size()-1,0);
        }


        level--;
    }


    auto BandStructure::master_loop_run(int mode,std::vector<int> &table,int level,int flag) -> void
    {
        int id,numprocs;
        MPI_Comm_rank(MPI_COMM_WORLD,&id);
        MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

        if(flag) timing.reset();
        while(table[level]!=End)
        {
            loopTiming.reset();
            switch(table[level])
            {
                case OutputEpsilon:
                    Output_Epsilon();
                    hamiltonian.systemcell.spNum.resize(0);
                    hamiltonian.systemcell.spStr.resize(0);
                    hamiltonian.EigenVectorsBlocks.resize(0);
                    hamiltonian.EigenValuesBlocks.resize(0);
                    hamiltonian.matrix.resize(0,0);
                    hamiltonian.dMatrix.resize(0,0);
                    hamiltonian.matrix3Db.resize(0,0);
                    hamiltonian.matrixScatter.resize(0,0);
                    hamiltonian.gsoul.resize(0);
                    hamiltonian.gsoul3D.resize(0);
                    hamiltonian.systemcell.supercellMargin.resize(0);
                    hamiltonian.systemcell.spNum.resize(0);
                    hamiltonian.systemcell.fspNum.resize(0);
                    hamiltonian.systemcell.spStr.resize(0);

                    //hamiltonian.systemcell.ListofMatrixElements_Type2.resize(0);
                    //hamiltonian.systemcell.ListofMatrixElements_Type3.resize(0);
                    hamiltonian.systemcell.ListofMatrixCoupling.resize(0);

                    MakeDecision(table,level);
                    break;
                case COMPUTE_Epsilon:
                    GenerateChi();
                    MakeDecision(table,level);
                    //RunningFlag--;
                    break;
                case MeshBZ:
                    MeshBrillouinZone();
                    MakeDecision(table,level);
                    break;
                case Init:
                    MakeDecision(table,level);
                    break;
                case CustomInit:
                    CustomBandInit();
                    MakeDecision(table,level);
                    break;
                case AutoMaticInit:
                    AutoMaticBandInit();
                    MakeDecision(table,level);
                    break;
                case With_ek3D:
                    hamiltonian.KextraTransform.resize(hamiltonian.KextraNum);
                    for(int i=0;i<hamiltonian.KextraNum;++i)
                    {
                        hamiltonian.KextraTransform[i].resize(0,0);
                    }
                    for(int i=0;i<hamiltonian.KextraNum;++i)
                    {
                        //hamiltonian.get_transform_remote(symmetryPoints.symmetryPoints[hamiltonian.KextraName[i]].position,hamiltonian.KextraTransform[i],hamiltonian.KextraRange[i].first,hamiltonian.KextraRange[i].second,hamiltonian.systemcell.ListofMatrixElements_Type1);
                        hamiltonian.get_transform_remote(symmetryPoints.symmetryPoints[hamiltonian.KextraName[i]].position,hamiltonian.KextraTransform[i],hamiltonian.KextraRange[i].first,hamiltonian.KextraRange[i].second,hamiltonian.systemcell.ListofMatrixElements_Type2);
                    }
                    MakeDecision(table,level);
                    break;
                case pre_calculation:
                    MakeDecision(table,level);

                    //if(hamiltonian.USE_decomposition)
                    //{
                    //conductionpos=hamiltonian.systemcell.UnitCellNumber*hamiltonian.BandGapIndex-hamiltonian.matrix_bottom;
                    //valencepos=conductionpos-1;
                    //}else{
                    conductionpos=AtomNumbers*hamiltonian.BandGapIndex-hamiltonian.matrix_bottom;
                    valencepos=conductionpos-1;
                    //}

                    //std::cout<<"numCorrect="<<conductionpos<<"\tworing="<<hamiltonian.MatrixSize/hamiltonian.orbNum*hamiltonian.BandGapIndex-hamiltonian.matrix_bottom<<std::endl;
                    //RunningFlag--;
                    break;
                case blocks_diagonalization3:
                    MakeDecision(table,level);
                    hamiltonian.matrix_bottom=0;
                    if(USE_parallel_in_preCalc) hamiltonian.pzheevx3PreDiagonalize(hamiltonian.EigenVectorsBlocks,hamiltonian.EigenValuesBlocks,hamiltonian.matrix,false,hamiltonian.matrix,0);
                    if(USE_parallel_in_preCalc==0) hamiltonian.pzheevxPreDiagonalize(kpoints[startPoint+computeflag],hamiltonian.EigenVectorsBlocks,hamiltonian.EigenValuesBlocks,hamiltonian.matrix,false,hamiltonian.matrix,USE_pzheevx_OR_Eigen_in_preCalc);
                    MatrixBottom.push_back(hamiltonian.matrix_bottom);
                    //RunningFlag--;
                    break;
                case blocks_diagonalization9:
                    MakeDecision(table,level);
                    hamiltonian.matrix_bottom=0;
                    hamiltonian.PfoldDiagonalize(hamiltonian.PgsoulFold_matrix,hamiltonian.matrix3Db,hamiltonian.USE_extra_k_point,hamiltonian.matrix,kpoints[startPoint+computeflag],tmpName);
                    //hamiltonian.foldDiagonalize(hamiltonian.gsoulFold_matrix,hamiltonian.matrix3Db,hamiltonian.USE_extra_k_point,hamiltonian.matrix,kpoints[startPoint+computeflag],tmpName);
                    MatrixBottom.push_back(hamiltonian.matrix_bottom);
                    //RunningFlag--;
                    break;
                case blocks_diagonalization8:

                    //hamiltonian.foldDiagonalize(hamiltonian.PureFold_matrix2,hamiltonian.matrix3Db,hamiltonian.USE_extra_k_point,hamiltonian.matrix,kpoints[startPoint+computeflag],-999);
                    //hamiltonian.foldDiagonalize(hamiltonian.PureFold_matrix,hamiltonian.matrix3Db,hamiltonian.USE_extra_k_point,hamiltonian.matrix,kpoints[startPoint+computeflag],999);
                    MakeDecision(table,level);
                    hamiltonian.matrix_bottom=0;
                    hamiltonian.foldDiagonalize(hamiltonian.gsoulFold_matrix,hamiltonian.matrix3Db,hamiltonian.USE_extra_k_point,hamiltonian.matrix,kpoints[startPoint+computeflag],tmpName);
                    MatrixBottom.push_back(hamiltonian.matrix_bottom);
                    //RunningFlag--;
                    break;
                case blocks_diagonalization:
                    MakeDecision(table,level);
                    hamiltonian.matrix_bottom=0;
                    //if(USE_blocks_transform) hamiltonian.pDiagonalize(hamiltonian.soul,hamiltonian.matrix,USE_transform_matrix,hamiltonian.matrix);
                    //if(USE_transform_matrix) hamiltonian.pDiagonalize(hamiltonian.soul,hamiltonian.matrix,true,hamiltonian.matrix);
                    if(USE_blocks_transform_bulk_states==0) hamiltonian.pDiagonalize(hamiltonian.soul,hamiltonian.matrix,USE_transform_matrix,hamiltonian.matrix);
                    if(USE_blocks_transform_bulk_states) hamiltonian.preDiagonalize(hamiltonian.gsoul3D,hamiltonian.matrix3Db,hamiltonian.USE_extra_k_point,hamiltonian.matrix);
                    if(USE_bulk_states_construction) hamiltonian.foldDiagonalize(hamiltonian.gsoulFold_matrix,hamiltonian.matrix3Db,hamiltonian.USE_extra_k_point,hamiltonian.matrix,kpoints[startPoint+computeflag],tmpName);
                    MatrixBottom.push_back(hamiltonian.matrix_bottom);
                    //RunningFlag--;
                    break;
                case reconstruct_Hamiltonian3:
                    MakeDecision(table,level);
                    hamiltonian.FillDiagonalBlocksByData(hamiltonian.matrix,hamiltonian.EigenValuesBlocks);
                    hamiltonian.BlockCoupling(hamiltonian.matrix,hamiltonian.EigenVectorsBlocks);
                    //RunningFlag--;
                    break;
                case reconstruct_Hamiltonian_mpi:
                    MakeDecision(table,level);
                    hamiltonian.pzheevxFillDiagonalBlocksByData(hamiltonian.matrix,hamiltonian.EigenValuesBlocks);
                    hamiltonian.pzheevxBlockCoupling(hamiltonian.matrix,hamiltonian.EigenVectorsBlocks);
                    //RunningFlag--;
                    break;
                case transform_Hamiltonian:
                    MakeDecision(table,level);
                    hamiltonian.matrix=hamiltonian.transform_matrix.adjoint()*hamiltonian.matrix*hamiltonian.transform_matrix;
                    //RunningFlag--;
                    break;
                case reconstruct_Hamiltonian4:
                    MakeDecision(table,level);
                    hamiltonian.DestoryCouplingElements(hamiltonian.matrix);
                    hamiltonian.FillDiagonalBlocksByTransform(hamiltonian.matrix,hamiltonian.gsoul3D);
                    //hamiltonian.FillDiagonalBlocksByData(hamiltonian.matrix,hamiltonian.gsoul3D);
                    hamiltonian.BlockCoupling(hamiltonian.matrix,hamiltonian.gsoul3D);
                    hamiltonian.Append_matrix();
                    //RunningFlag--;
                    break;
                case reconstruct_Hamiltonian6:
                    MakeDecision(table,level);
                    //hamiltonian.SparseTransform(hamiltonian.matrix,hamiltonian.gsoulFold_matrix);
                    hamiltonian.PTransform(hamiltonian.dMatrix2,hamiltonian.PgsoulFold_matrix);
                    break;
                case reconstruct_Hamiltonian5:
                    MakeDecision(table,level);
                    hamiltonian.FillDiagonalBlocksByTransform(hamiltonian.matrix,hamiltonian.gsoulFold_matrix);
                    hamiltonian.BlockCoupling(hamiltonian.matrix,hamiltonian.gsoulFold_matrix);
                    //hamiltonian.FillDiagonalBlocksByTransform(hamiltonian.matrix,hamiltonian.PureFold_matrix);
                    break;
                case reconstruct_Hamiltonian2:
                    MakeDecision(table,level);
                    hamiltonian.FillDiagonalBlocksByTransform(hamiltonian.matrix,hamiltonian.gsoul3D);
                    //hamiltonian.FillDiagonalBlocksByData(hamiltonian.matrix,hamiltonian.gsoul3D);
                    hamiltonian.BlockCoupling(hamiltonian.matrix,hamiltonian.gsoul3D);
                    //RunningFlag--;
                    break;
                case reconstruct_Hamiltonian:
                    MakeDecision(table,level);
                    hamiltonian.gFillDiagonalBlocksByData(hamiltonian.matrix,hamiltonian.soul);
                    hamiltonian.gBlockCoupling(hamiltonian.matrix,hamiltonian.soul);
                    //RunningFlag--;
                    break;
                case read_elements3:
                    MakeDecision(table,level);
                    hamiltonian.ReadCouplingElements3(hamiltonian.matrix);
                    //hamiltonian.dMatrix.setZero();
                    hamiltonian.dMatrix2.setConstant({0,0});
                    //RunningFlag--;
                    break;
                case read_elements:
                    MakeDecision(table,level);
                    //if(hamiltonian.USE_extra_k_point){
                    //hamiltonian.ReadCouplingElements2(hamiltonian.matrix);
                    //}else{
                    hamiltonian.ReadCouplingElements(hamiltonian.matrix);
                    //}
                    if(USE_blocks_transform_bulk_states==0) hamiltonian.matrix.setZero();
                    //RunningFlag--;
                    break;
                case still: MakeDecision(table,level);
                            if(USE_distributed==0)
                            {
                                hamiltonian.transform_matrix=Eigen::MatrixXcd::Identity(hamiltonian.MatrixSizeReduced,hamiltonian.MatrixSizeReduced);
                                hamiltonian.transform_mix=Eigen::MatrixXcd::Identity(hamiltonian.MatrixSizeReduced,hamiltonian.MatrixSizeReduced);
                            }
                            hamiltonian.matrix_bottom=0;
                            MatrixBottom.push_back(hamiltonian.matrix_bottom);
                            break;
                case output_data:
                            MakeDecision(table,level);
                            if(hamiltonian.AutoMatic<=0)
                            {
                                BandOutputW(eigenvalues,DataTag,midlen);
                            }else{
                                BandOutputW(eigenvalues,thrd,midlen);
                            }
                            eigenvectors.resize(0,0);
                            eigenvalues.resize(0);
                            //RunningFlag--;
                            break;
                case diagonalize_pzheevx:
                            MakeDecision(table,level);
                            //hamiltonian.print_matrix(hamiltonian.dMatrix2.data(),"dM");

                            if(USE_bulk_states_construction)
                            {
                                //hamiltonian.Diagonal_parallel_pzheevx_all(eigenvalues,eigenvectors,hamiltonian.systemcell.supercell[0]->Startp,hamiltonian.systemcell.supercell[0]->Startp,hamiltonian.MatrixSizeReduced,1,1,hamiltonian.MatrixSize,0);
                                hamiltonian.Diagonal_parallel_pzheevx_all(eigenvalues,eigenvectors,0,0,hamiltonian.MatrixSizeReduced,1,1,hamiltonian.MatrixSize,0);
                                //hamiltonian.Diagonal_parallel_pzheevx_all(eigenvalues,eigenvectors,0,0,hamiltonian.MatrixSizeReduced,1,1,hamiltonian.MatrixSizeReduced,0);
                            }else{
                                hamiltonian.Diagonal_parallel_pzheevx_all(eigenvalues,eigenvectors,0,0,hamiltonian.MatrixSizeReduced,1,1,hamiltonian.MatrixSize,0);
                                //hamiltonian.Diagonal_parallel_pzheevx_test(eigenvalues,eigenvectors,0,0,hamiltonian.MatrixSizeReduced,1,1,hamiltonian.MatrixSize,0);
                            }
                            //std::cout<<"egv1111="<<eigenvalues<<std::endl;
                            //hamiltonian.gMatrixDiagonalize(false,Eigen::MatrixXcd::Identity(hamiltonian.MatrixSizeReduced,hamiltonian.MatrixSizeReduced));
                            //eigenvalues=hamiltonian.solver.eigenvalues();
                            //std::cout<<"egv1111="<<eigenvalues<<std::endl;

                            break;
                case diagonalize:
                            MakeDecision(table,level);
                            break;
                case diagonalize_Eigen:
                            MakeDecision(table,level);
                            hamiltonian.gMatrixDiagonalize(true,Eigen::MatrixXcd::Identity(hamiltonian.MatrixSizeReduced,hamiltonian.MatrixSizeReduced));
                            if(true)
                            {
                                eigenvectors=hamiltonian.solver.eigenvectors();
                                //hamiltonian.print_matrix(hamiltonian.matrix.block(0,0,hamiltonian.MatrixSizeReduced,hamiltonian.MatrixSizeReduced),"Matrix");
                            }
                            eigenvalues=hamiltonian.solver.eigenvalues();
                            //FindBandEdge(std::complex<double>(0),Eigen::Vector3d(0,0,1),eigenvectors,eigenvalues,6);
                            //CalculateMassByNumerical(Eigen::Vector3d(0,0,1e-4),Eigen::Vector3d(0,0,0),6);

                            //RunningFlag--;
                            break;
                case diagonalize_Mix:
                            MakeDecision(table,level);
                            hamiltonian.MatrixDiagonalize(false,hamiltonian.transform_mix.adjoint()*hamiltonian.transform_mix);
                            //std::cout<<"matrixSize="<<hamiltonian.matrix.rows()<<"\tReducedSize="<<hamiltonian.MatrixSizeReduced<<"\tmix_size="<<(hamiltonian.transform_mix.adjoint()*hamiltonian.transform_mix).rows()<<std::endl;
                            //hamiltonian.MatrixDiagonalize(false,Eigen::MatrixXcd::Identity(hamiltonian.MatrixSizeReduced,hamiltonian.MatrixSizeReduced));

                            if(false)
                            {
                                eigenvectors=hamiltonian.gsolver.eigenvectors();
                            }
                            eigenvalues=hamiltonian.gsolver.eigenvalues();

                            //hamiltonian.gMatrixDiagonalize(false,Eigen::MatrixXcd::Identity(hamiltonian.MatrixSizeReduced,hamiltonian.MatrixSizeReduced));
                            //eigenvectors=hamiltonian.solver.eigenvectors();
                            //eigenvalues=hamiltonian.solver.eigenvalues();
                            //RunningFlag--;
                            break;
                case set_Hamiltonian_distributed:
                            MakeDecision(table,level);
                            hamiltonian.SetMatrixLocWithoutShift(kpoints[startPoint+computeflag],hamiltonian.systemcell.ListofMatrixElements_Type1,hamiltonian.dMatrix2);
                            //OutPutStatus(0);
                            //RunningFlag--;
                            break;
                case set_Hamiltonian_distributed_with_shift:
                            MakeDecision(table,level);
                            hamiltonian.SetMatrixLoc(kpoints[startPoint+computeflag],hamiltonian.systemcell.ListofMatrixElements_Type1,hamiltonian.dMatrix2);
                            //hamiltonian.SetMatrixLocWithoutShift(kpoints[startPoint+computeflag],hamiltonian.systemcell.ListofMatrixElements_Type1,hamiltonian.dMatrix2);
                            break;
                case set_Hamiltonian_sparse2:
                            MakeDecision(table,level);
                            hamiltonian.SetMatrixSparse2(kpoints[startPoint+computeflag],hamiltonian.systemcell.ListofMatrixElements_Type1,hamiltonian.matrix);
                            break;
                case set_Hamiltonian_sparse:
                            MakeDecision(table,level);
                            hamiltonian.SetMatrixSparse(kpoints[startPoint+computeflag],hamiltonian.systemcell.ListofMatrixElements_Type1,hamiltonian.matrix);
                            break;
                case COMPUTE_5:
                            if(1)
                            {
                                std::complex<double> ii(0,1);
                                hamiltonian.Case5(kpoints[startPoint+computeflag](0)+1e-4*ii);
                                //hamiltonian.CaseComplexBands5(kpoints[startPoint+computeflag](0)+1e-7*ii,kpoints[startPoint+computeflag]);
                                //hamiltonian.CaseSurfaceDOS(kpoints[startPoint+computeflag](0)+1e-3*ii,kpoints[startPoint+computeflag]);
                                //hamiltonian.ComplexBands(kpoints[startPoint+computeflag](0)+1e-7*ii,kpoints[startPoint+computeflag]);
                            }
                            MakeDecision(table,level);
                            break;
                case set_Hamiltonian_Type1:
                            MakeDecision(table,level);
                            hamiltonian.SetMatrix(kpoints[startPoint+computeflag],hamiltonian.systemcell.ListofMatrixElements_Type1,hamiltonian.matrix);
                            break;
                case set_Hamiltonian_Type2:
                            MakeDecision(table,level);
                            hamiltonian.SetMatrix(kpoints[startPoint+computeflag],hamiltonian.systemcell.ListofMatrixElements_Type2,hamiltonian.matrix3Db);
                            break;
                case set_Hamiltonian_Type3:
                            MakeDecision(table,level);
                            hamiltonian.SetMatrixFold(kpoints[startPoint+computeflag],hamiltonian.systemcell.ListofMatrixElements_Type3,hamiltonian.matrixScatter);
                            break;
                case set_Hamiltonian_TypeVec:
                            MakeDecision(table,level);
                            hamiltonian.vecMatrixFold.resize(0);
                            hamiltonian.vecMatrixFold.push_back(Eigen::MatrixXcd());
                            hamiltonian.vecSetMatrixFold(kpoints[startPoint+computeflag],hamiltonian.systemcell.ListofMatrixElements_Type3,hamiltonian.vecMatrixFold.back(),"AlAs");
                            hamiltonian.vecMatrixFold.push_back(Eigen::MatrixXcd());
                            hamiltonian.vecSetMatrixFold(kpoints[startPoint+computeflag],hamiltonian.systemcell.ListofMatrixElements_Type3,hamiltonian.vecMatrixFold.back(),"GaAs");
                            break;
                case set_Hamiltonian:
                            MakeDecision(table,level);
                            OutPutStatus(USE_parallel_K_points);
                            //RunningFlag--;
                            break;
                case design_running_pattern_single:
                            MakeDecision(table,level);
                            //ComputePrepareMPI();
                            //RunningFlag--;
                            break;
                case design_running_pattern:
                            MakeDecision(table,level);
                            ComputePrepareDFTB();
                            //RunningFlag--;
                            break;
                case reshape_structure:
                            MakeDecision(table,level);
                            if(USE_blocks_transform)
                            {
                                hamiltonian.systemcell.ProcessCell();
                                //hamiltonian.systemcell.GenSuperCellRequared();
                            }
                            hamiltonian.systemcell.GenSuperCellRequared();
                            //hamiltonian.print_Hamiltonian();
                            //RunningFlag--;
                            break;
                case find_neighbors:
                            MakeDecision(table,level);
                            /*
                               hamiltonian.systemcell.SelectedMasterGen(
                               (USE_LoMEType1?LoMEType1:0)
                               | (USE_LoMEType2?LoMEType2:0)
                               | (USE_LoMEType3?LoMEType3:0)
                               | (USE_LoMC?LoMC:0)
                               );
                               */

                            hamiltonian.systemcell.SelectedMasterGen(
                                    (USE_LoMEType1?LoMEType1:0)
                                    | (USE_LoMEType2?LoMEType2:0)
                                    | (USE_LoMEType3?LoMEType3:0)
                                    | (USE_LoMC?LoMC:0)
                                    );

                            //RunningFlag--;
                            break;
                case create_structure:
                            MakeDecision(table,level);
                            hamiltonian.systemcell.StackingStructure();
                            hamiltonian.MatrixSize=hamiltonian.systemcell.GetMatrixSize();
                            hamiltonian.MatrixSizeReduced=hamiltonian.MatrixSize;
                            //RunningFlag--;
                            break;
                case COMPUTE_K_POINTS:
                            MakeDecision(table,level);
                            break;
                case Output_and_Recycle:
                            H_clean();
                            MakeDecision(table,level);
                            //NewOutput(id,startPoint,DataTag);
                            break;
                            //default:
                            //MakeDecision(table,level);//RunningFlag--;
            }

            if(id==0)
            {
                std::fstream file_obj; 
                file_obj.open("./log"+std::to_string(tmpName),std::ios::out|std::ios::app);
                file_obj<<std::fixed<<std::setprecision(4)<<loopTiming.elapsed()<<"   "<<TableOFactionsName[static_cast<TableOFactions>(table[level+1])]<<std::endl;
                file_obj.close();
            }

            loopTiming.reset();
        }
    }
    /*
       void BandStructure::PrintMass(const Eigen::VectorXd &eigenvals,Eigen::VectorXcd &Mmass,Eigen::VectorXcd &Mmass2,int leftSf,int rightSf)
       {
       if((startPoint+computeflag)==kppos)
       {
       std::cout<<std::endl;
       std::cout<<"====="<<std::endl;
       std::cout<<"("<<conductionedge<<","<<valenceedge<<")BandGap="<<(conductionedge-valenceedge)<<std::endl;


       for(int ii=-leftSf;ii<rightSf;++ii)
       {
       int pos=AtomNumbers*hamiltonian.BandGapIndex-hamiltonian.matrix_bottom+ii;
    //pos+=1;
    if(ii==-1) std::cout<<std::endl;
    if(pos>=0 && pos<Mmass.size())
    {
    std::cout<<std::fixed<<std::setprecision(20)
    <<"["<<(pos)<<"]\t"
    <<(1.0/(Mmass(pos).real()));
    }
    if(pos>=0 && pos<eigenvals.size())
    {
    std::cout<<"\t("<<eigenvals(pos)<<")";
    }
    if(pos>=0 && pos<Mmass2.size())
    {
    double ratio=Mmass(pos).real()/Mmass2(pos).real();
    ratio=(fabs(ratio-1)<1e-3)?1:ratio;
    std::cout<<"\t"<<1.0/(Mmass2(pos).real())<<"\t"
    <<ratio;
    }
    std::cout<<std::endl;
    if(ii==0) std::cout<<std::endl;
    }

    std::cout<<"====="<<std::endl;
    std::cout<<std::endl;

    }
    }
    */

    void BandStructure::FindBandEdge(const std::complex<double> kz,Eigen::Vector3d pl,const Eigen::MatrixXcd &VectM,const Eigen::VectorXd &eigenvals,int ShiftValue)
    {
        int fastmode=0;
        int printRange=20;
        conductionpos=-99;
        valencepos=-99;

        Eigen::MatrixXcd Kmatrix,KKmatrix;
        hamiltonian.SetMatrixK5(VectM,kz,pl.normalized(),hamiltonian.systemcell.ListofMatrixElements_Type1,Kmatrix,1);
        hamiltonian.SetMatrixK5(VectM,kz,pl.normalized(),hamiltonian.systemcell.ListofMatrixElements_Type1,KKmatrix,2);

        double hhbar=1.05457162825*9.10938188;
        std::cout<<std::endl;

        for(int i=0;i<ShiftValue;++i)
        {
            int pos=AtomNumbers*hamiltonian.BandGapIndex-hamiltonian.matrix_bottom+i;
            if(pos>=0 && pos<eigenvals.size())
            {

                std::complex<double> Mass=KKmatrix(pos,pos);
                for(int inm=0;inm<Kmatrix.rows();++inm)
                {
                    if(abs(eigenvals(pos)-eigenvals(inm))>1e-4)
                    {
                        Mass+=2.*pow(std::abs(Kmatrix(inm,pos)),2)/(eigenvals(pos)-eigenvals(inm));
                    }
                }
                double conductionmass=hhbar/Mass.real();
                if(i<printRange)
                {
                    std::cout
                        <<std::fixed<<std::setprecision(20)
                        <<"["<<(pos)<<"]\t"
                        <<(conductionmass)<<"\t("
                        <<eigenvals(pos)
                        <<std::endl;
                }
                if(fastmode){
                    if((conductionmass-1e-3)>0)
                    {
                        conductionpos=pos;
                        break;
                    }
                }else{
                    //if((conductionmass-1e-3)>0 && conductionpos==-99)
                    if(conductionpos==-99)
                    {
                        conductionpos=pos;
                    }
                }
            }
        }

        std::cout<<std::endl;

        for(int i=1;i<ShiftValue;++i)
        {
            int pos=AtomNumbers*hamiltonian.BandGapIndex-hamiltonian.matrix_bottom-i;
            if(pos>=0 && pos<eigenvals.size())
            {

                std::complex<double> Mass=KKmatrix(pos,pos);
                for(int inm=0;inm<Kmatrix.rows();++inm)
                {
                    if(abs(eigenvals(pos)-eigenvals(inm))>1e-4)
                    {
                        Mass+=2.*pow(std::abs(Kmatrix(inm,pos)),2)/(eigenvals(pos)-eigenvals(inm));
                    }
                }

                double valencemass=hhbar/Mass.real();
                if(i<printRange)
                {
                    std::cout<<std::fixed<<std::setprecision(20)
                        <<"["<<(pos)<<"]\t"
                        <<(valencemass)<<"\t("
                        <<eigenvals(pos)
                        <<std::endl;
                }
                if(fastmode){
                    if((valencemass+1e-3)<0)
                    {
                        valencepos=pos;
                        break;
                    }
                }else{
                    if((valencemass+1e-3)<0 && valencepos == -99)
                    {
                        valencepos=pos;
                    }
                }
            }
        }
        std::cout
            <<std::fixed<<std::setprecision(5)
            <<"("<<eigenvals(conductionpos)<<","<<eigenvals(valencepos)<<")BandGap="<<(eigenvals(conductionpos)-eigenvals(valencepos))<<std::endl;
        std::cout<<std::endl;


    }
    void BandStructure::CalculateMassByNumerical(const Eigen::Vector3d zz,const Eigen::Vector3d gr,Eigen::VectorXcd &fNumMass)
    {
        Eigen::MatrixXcd _m,m_,m;
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> _Msolver;
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> Msolver;
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> M_solver;

        hamiltonian.SetMatrix(kpoints[startPoint+computeflag]+gr-zz,hamiltonian.systemcell.ListofMatrixElements_Type1,_m);
        hamiltonian.SetMatrix(kpoints[startPoint+computeflag]+gr,hamiltonian.systemcell.ListofMatrixElements_Type1,m);
        hamiltonian.SetMatrix(kpoints[startPoint+computeflag]+gr+zz,hamiltonian.systemcell.ListofMatrixElements_Type1,m_);

        _Msolver.compute(_m,Eigen::MatrixXcd::Identity(_m.rows(),_m.cols()), Eigen::DecompositionOptions::ComputeEigenvectors|Eigen::Ax_lBx);
        Msolver.compute(m,Eigen::MatrixXcd::Identity(m.rows(),m.cols()), Eigen::DecompositionOptions::ComputeEigenvectors|Eigen::Ax_lBx);
        M_solver.compute(m_,Eigen::MatrixXcd::Identity(m_.rows(),m_.cols()), Eigen::DecompositionOptions::ComputeEigenvectors|Eigen::Ax_lBx);
        NumMass=(M_solver.eigenvalues()+_Msolver.eigenvalues()-2.*Msolver.eigenvalues());
        NumMass=NumMass/zz.squaredNorm();
    }

    void BandStructure::CalculateMassByNumerical(const Eigen::Vector3d zz,const Eigen::Vector3d gr,double ShiftValue)
    {
        Eigen::MatrixXcd _m,m_,m;
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> _Msolver;
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> Msolver;
        Eigen::GeneralizedSelfAdjointEigenSolver<Eigen::MatrixXcd> M_solver;

        hamiltonian.SetMatrix(kpoints[startPoint+computeflag]+gr-zz,hamiltonian.systemcell.ListofMatrixElements_Type1,_m);
        //hamiltonian.preDiagonalize(hamiltonian.gsoul,_m,true,_m);
        int _mb=0;//hamiltonian.matrix_bottom;
        //_m=hamiltonian.transform_matrix.adjoint()*_m*hamiltonian.transform_matrix;
        /*
           hamiltonian.ReadCouplingElements(_m);
           _m.setZero();
           hamiltonian.FillDiagonalBlocksByData(_m,hamiltonian.gsoul);
           hamiltonian.BlockCoupling(_m,hamiltonian.gsoul);
           */
        hamiltonian.SetMatrix(kpoints[startPoint+computeflag]+gr,hamiltonian.systemcell.ListofMatrixElements_Type1,m);
        //hamiltonian.preDiagonalize(hamiltonian.gsoul,m,true,m);
        int mb=0;//hamiltonian.matrix_bottom;
        //m=hamiltonian.transform_matrix.adjoint()*m*hamiltonian.transform_matrix;
        /*
           hamiltonian.ReadCouplingElements(m);
           m.setZero();
           hamiltonian.FillDiagonalBlocksByData(m,hamiltonian.gsoul);
           hamiltonian.BlockCoupling(m,hamiltonian.gsoul);
           */

        hamiltonian.SetMatrix(kpoints[startPoint+computeflag]+gr+zz,hamiltonian.systemcell.ListofMatrixElements_Type1,m_);
        //hamiltonian.preDiagonalize(hamiltonian.gsoul,m_,true,m_);
        int mb_=0;//hamiltonian.matrix_bottom;
        //m_=hamiltonian.transform_matrix.adjoint()*m_*hamiltonian.transform_matrix;
        /*
           hamiltonian.ReadCouplingElements(m_);
           m_.setZero();
           hamiltonian.FillDiagonalBlocksByData(m_,hamiltonian.gsoul);
           hamiltonian.BlockCoupling(m_,hamiltonian.gsoul);
           */


        _Msolver.compute(_m,Eigen::MatrixXcd::Identity(_m.rows(),_m.cols()), Eigen::DecompositionOptions::ComputeEigenvectors|Eigen::Ax_lBx);
        Msolver.compute(m,Eigen::MatrixXcd::Identity(m.rows(),m.cols()), Eigen::DecompositionOptions::ComputeEigenvectors|Eigen::Ax_lBx);
        M_solver.compute(m_,Eigen::MatrixXcd::Identity(m_.rows(),m_.cols()), Eigen::DecompositionOptions::ComputeEigenvectors|Eigen::Ax_lBx);

        const Eigen::VectorXd &eigenvals=Msolver.eigenvalues();
        const Eigen::VectorXd &reigenvals=M_solver.eigenvalues();
        const Eigen::VectorXd &leigenvals=_Msolver.eigenvalues();

        //NumMass=(M_solver.eigenvalues()+_Msolver.eigenvalues()-2.*Msolver.eigenvalues())/zz.squaredNorm();

        int fastmode=0;
        int printRange=10;
        conductionpos=-99;
        valencepos=-99;

        std::cout<<std::endl;
        //std::cout<<"Index="<<DataTag<<std::endl;
        //std::cout<<std::fixed<<std::setprecision(20)<<"leig="<<leigenvals<<std::endl;
        //std::cout<<std::fixed<<std::setprecision(20)<<"eig="<<eigenvals<<std::endl;
        //std::cout<<std::fixed<<std::setprecision(20)<<"reig="<<reigenvals<<std::endl;

        double hhbar=1.05457162825*9.10938188;
        //double hhbar=0.106;
        //double hhbar=0.156;
        //double hhbar=0.109;
        //double hhbar=pow(hhh/(m_e*2*M_PI),1)*1e3;
        //double hhbar=1./(m_e*1e31);
        for(int i=0;i<ShiftValue;++i)
        {
            int pos=AtomNumbers*hamiltonian.BandGapIndex+i;
            if(
                    (pos-mb>=0 && pos-mb<eigenvals.size())
                    && (pos-mb_>=0 && pos-mb_<reigenvals.size())
                    && (pos-_mb>=0 && pos-_mb<leigenvals.size())
              )
            {

                double NumMassPos=(M_solver.eigenvalues()(pos-mb_)+_Msolver.eigenvalues()(pos-_mb)-2.*Msolver.eigenvalues()(pos-mb))/zz.squaredNorm();
                std::complex<double> Mass=NumMassPos;
                //std::complex<double> Mass=NumMass(pos);

                double conductionmass=hhbar/NumMassPos;
                if(i<printRange)
                {
                    std::cout
                        <<std::fixed<<std::setprecision(20)
                        <<"["<<(pos-mb)<<"]\t"
                        //<<(NumMassPos)<<"\t("
                        <<(conductionmass)<<"\t("
                        <<eigenvals(pos-mb)
                        <<std::endl;
                }
                if(fastmode){
                    if((conductionmass-1e-3)>0)
                    {
                        conductionpos=pos-mb;
                        break;
                    }
                }else{
                    if((conductionmass-1e-3)>0 && conductionpos==-99)
                    {
                        conductionpos=pos-mb;
                    }
                }
            }
        }

        std::cout<<std::endl;

        for(int i=1;i<ShiftValue;++i)
        {
            int pos=AtomNumbers*hamiltonian.BandGapIndex-i;
            if(
                    (pos-mb>=0 && pos-mb<eigenvals.size())
                    && (pos-mb_>=0 && pos-mb_<reigenvals.size())
                    && (pos-_mb>=0 && pos-_mb<leigenvals.size())
              )
            {

                double NumMassPos=(M_solver.eigenvalues()(pos-mb_)+_Msolver.eigenvalues()(pos-_mb)-2.*Msolver.eigenvalues()(pos-mb))/zz.squaredNorm();
                std::complex<double> Mass=NumMassPos;
                //std::complex<double> Mass=NumMass(pos);

                double valencemass=hhbar/NumMassPos;
                if(i<printRange)
                {
                    std::cout<<std::fixed<<std::setprecision(20)
                        <<"["<<(pos-mb)<<"]\t"
                        //<<(NumMassPos)<<"\t("
                        <<(valencemass)<<"\t("
                        <<eigenvals(pos-mb)
                        <<std::endl;
                }
                if(fastmode){
                    if((valencemass+1e-3)<0)
                    {
                        valencepos=pos-mb;
                        break;
                    }
                }else{
                    if((valencemass+1e-3)<0 && valencepos == -99)
                    {
                        valencepos=pos-mb;
                    }
                }
            }
        }
        std::cout
            <<std::fixed<<std::setprecision(5)
            <<"("<<eigenvals(conductionpos)<<","<<eigenvals(valencepos)<<")BandGap="<<(eigenvals(conductionpos)-eigenvals(valencepos))<<std::endl;
        std::cout<<std::endl;

    }
    void BandStructure::CalculateMass(const Eigen::MatrixXcd &VectM,const Eigen::VectorXd &eigenvals,Eigen::VectorXcd &Mmass)
    {

        Eigen::MatrixXcd EEmatrix,Ematrix,matrix1,matrix_1,Kmatrix,KKmatrix,Mmatrix;
        Mmass.resize(hamiltonian.MatrixSizeReduced);
        Mmass.setZero();

        Eigen::Vector3d pl={1.,1.,1.};
        std::complex<double> kz=0.0;
        hamiltonian.SetMatrixK5(VectM,kz*0.,pl.normalized(),hamiltonian.systemcell.ListofMatrixElements_Type1,Kmatrix,1);
        hamiltonian.SetMatrixK5(VectM,kz*0.,pl.normalized(),hamiltonian.systemcell.ListofMatrixElements_Type1,KKmatrix,2);

        for(int mi=0;mi<Mmass.size();++mi)
            //for(int mi=(AtomNumbers*hamiltonian.BandGapIndex-hamiltonian.matrix_bottom)/2;mi<(AtomNumbers*hamiltonian.BandGapIndex-hamiltonian.matrix_bottom)+5;++mi)
        {
            Mmass(mi)=KKmatrix(mi,mi);
            for(int inm=0;inm<Kmatrix.rows();++inm)
                //for(int inm=(AtomNumbers*hamiltonian.BandGapIndex-hamiltonian.matrix_bottom-pushValue);inm<(AtomNumbers*hamiltonian.BandGapIndex-hamiltonian.matrix_bottom+pushValue);++inm)
            {
                if(abs(eigenvals(mi)-eigenvals(inm))>1e-4)
                    //if(eigenvals(mi)!=eigenvals(inm))
                {
                    Mmass(mi)+=2.*pow(std::abs(Kmatrix(inm,mi)),2)/(eigenvals(mi)-eigenvals(inm));
                    //Mmass(mi)+=1.732*(Kmatrix(mi,inm)*Kmatrix(inm,mi))/(eigenvals(mi)-eigenvals(inm));
                    //Mmass(mi)+=2.*(Kmatrix(mi,inm)*Kmatrix(inm,mi))/(eigenvals(mi)-eigenvals(inm));
                }
            }
        }

        //FindBandEdge(eigenvals,Mmass,AtomNumbers*2-hamiltonian.matrix_bottom);
        //end== 
    }

    auto BandStructure::ParseParameters() -> void
    {
        int id,numprocs;
        MPI_Comm_rank(MPI_COMM_WORLD,&id);
        MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

        hamiltonian.ReadMaterials();
        double BackupValue=0;
        for(auto &pInt : parseUSE_int)
        {
            pInt.second=ParseValue<int>(pInt.first,inpfile,std::to_string(tmpName),BackupValue,id);
        }
        for(auto &pInt : parseElements)
        {
            pInt.second=ParseValue<int>(pInt.first,inpfile,std::to_string(tmpName),BackupValue,id);
        }
        for(auto &pInt : parseUSE_double)
        {
            pInt.second=ParseValue<double>(pInt.first,inpfile,std::to_string(tmpName),BackupValue,id);
        }
        for(auto &pInt : parseElementsDouble)
        {
            pInt.second=ParseValue<double>(pInt.first,inpfile,std::to_string(tmpName),BackupValue,id);
        }
        for(auto &pInt : parseString)
        {
            pInt.second=ParseValue<std::string>(pInt.first,inpfile,std::to_string(tmpName),"",id);
        }
        if(hamiltonian.USE_SOC) 
        {
            hamiltonian.BandGapIndex=4;
        }else{
            hamiltonian.BandGapIndex=2;
        }
        std::fstream file_obj; 
        if(id==0)
        {
            mkdir(strdup("./BandData"),S_IRWXU|S_IRWXG|S_IRWXO);
            mkdir(("./BandData/"+dataVersion+"0").c_str(),S_IRWXU|S_IRWXG|S_IRWXO);

            file_obj.open(("./out"+std::to_string(tmpName)+".txt").c_str(),std::ios::out|std::ios::app);

            if(USE_distributed)
            {
                file_obj<<"numprocs = "<<1<<std::endl;
            }else{
                file_obj<<"numprocs = "<<numprocs<<std::endl;
            }
            file_obj<<"InAsLength = "<<LengthBegin<<std::endl;
            file_obj<<"dataVersion = \""<<dataVersion<<"\""<<std::endl;
            file_obj.close();
        }
    }

    auto BandStructure::Initialize(std::vector<std::string> path, unsigned int kpnumber) -> void
    {
        //reset
        bandgapid=0;
        bandgap=0;
        ConductionBand.clear();
        ValenceBand.clear();
        BandValue.clear();
        kpoints.clear();
        results.clear();
        kpoints.reserve(kpnumber);
        //results.reserve(nrPoints);
        //calculate k points
        m_path.swap(path);

        kpoints = symmetryPoints.GeneratePoints(m_path, kpnumber, symmetryPointsPositions);///key
        kpoints_mesh = symmetryPoints.GeneratePoints(m_path, kpnumber, symmetryPointsPositions);

        if(USE_CalcDOS) RandomKpointList(1);//here DOS
    }

    auto BandStructure::AutoMaticBandInit() -> void
    {
        int id,numprocs;
        MPI_Comm_rank(MPI_COMM_WORLD,&id);
        MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
        ///Parse
        hamiltonian.systemcell.ParseCell(inpfile);
        hamiltonian.ParseParameters(inpfile);

        hamiltonian.MatrixSize=0;
        bandgapid=0;


        int bandlen=bandend-bandbegin+1;
        int lengthlen=LengthEnd-LengthBegin+1;
        midlen=LengthBegin+outputflag/bandlen;
        thrd=bandbegin+outputflag%bandlen;

        if(id==0)
        {
            char dirname[40];    
            if(hamiltonian.AutoMatic<=0)
            {
                sprintf(dirname,"./BandData/%s%d_%d",dataVersion.c_str(),DataTag,midlen);
            }else{
                sprintf(dirname,"./BandData/%s%d_%d",dataVersion.c_str(),thrd,midlen);
            }
            mkdir(dirname,S_IRWXU|S_IRWXG|S_IRWXO);
        }

        hamiltonian.systemcell.ProcessParameters(pow(2,thrd),midlen,hamiltonian.leftCutoff,hamiltonian.rightCutoff);

        symmetryPoints.GenSymmetryPoints(hamiltonian.systemcell.Unitg1,hamiltonian.systemcell.Unitg2,hamiltonian.systemcell.Unitg3,hamiltonian.systemcell.Unitg1,hamiltonian.systemcell.Unitg2,hamiltonian.systemcell.Unitg3);

        Initialize(split<std::string>(kppath,"->"),kpnumber);
        kppos=symmetryPointsPositions[std::distance(m_path.begin(),std::find(m_path.begin(),m_path.end(),"G"))];

        if(id==0 && FirstLoop)
        {
            std::fstream file_obj; 
            file_obj.open("./Chi"+std::to_string(tmpName)+".dat",std::ios::out|std::ios::trunc);
            file_obj.close();
            file_obj.open("./log"+std::to_string(tmpName),std::ios::out|std::ios::trunc);
            file_obj.close();



            /*
               if(DataAppending==0)
               {
               for(int wl=LengthBegin;wl<=LengthEnd;++wl){
               file_obj.open(("./BandData/"+dataVersion+"0/TimeSpend"+std::to_string(CurrentMode)+"_"+std::to_string(wl)+".data").c_str(),std::ios::out|std::ios::trunc);
               file_obj.close();
               file_obj.open(("./BandData/"+dataVersion+"0/ReducedHamiltonian"+std::to_string(CurrentMode)+"_"+std::to_string(wl)+".data").c_str(),std::ios::out|std::ios::trunc);
               file_obj.close();
               }
               }
               */
            file_obj.open(("./out"+std::to_string(tmpName)+".txt").c_str(),std::ios::out|std::ios::app);
            file_obj<<"kvalue =";
            for(int i=0;i<m_path.size();++i)
            {
                file_obj<<" "<<symmetryPointsPositions[i];
            }
            file_obj<<std::endl;
            file_obj<<"pnumber"<<" = "<<bandend<<std::endl;
            file_obj<<"kpos"<<" = "<<kppos<<std::endl;
            file_obj.close();
        }

        if(FirstLoop)
        {
            if(DataAppending==0)
            {
                std::fstream file_obj; 
                for(int fl=bandbegin;fl<=bandend;++fl){
                    for(int wl=LengthBegin;wl<=LengthEnd;++wl){
                        mkdir(("./BandData/"+dataVersion+std::to_string(fl)+"_"+std::to_string(wl)).c_str(),S_IRWXU|S_IRWXG|S_IRWXO);
                        file_obj.open(("./BandData/"+dataVersion+std::to_string(fl)+"_"+std::to_string(wl)+"/M"+std::to_string(CurrentMode)+"band"+std::to_string(id)+".data").c_str(),std::ios::out|std::ios::trunc);
                        file_obj.close();
                    }
                }
            }
        }


    }


    auto BandStructure::CustomBandInit() -> void
    {
        int id,numprocs;
        MPI_Comm_rank(MPI_COMM_WORLD,&id);
        MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
        ///Parse
        hamiltonian.systemcell.ParseCell(inpfile);
        hamiltonian.ParseParameters(inpfile);

        hamiltonian.MatrixSize=0;
        bandgapid=0;

        //midlen=LengthEnd;

        int bandlen=bandend-bandbegin+1;
        int lengthlen=LengthEnd-LengthBegin+1;
        midlen=LengthBegin+outputflag/bandlen;
        thrd=bandbegin+outputflag%bandlen;

        bandbegin=0;
        bandend=0;
        LengthBegin=0;
        LengthEnd=0;


        if(id==0)
        {
            char dirname[40];    
            sprintf(dirname,"./BandData/%s%d_%d",dataVersion.c_str(),DataTag,midlen);
            mkdir(dirname,S_IRWXU|S_IRWXG|S_IRWXO);
        }

        hamiltonian.systemcell.ParseCellExtra(inpfile);


        symmetryPoints.GenSymmetryPoints(hamiltonian.systemcell.Unitg1,hamiltonian.systemcell.Unitg2,hamiltonian.systemcell.Unitg3,hamiltonian.systemcell.Unitg1,hamiltonian.systemcell.Unitg2,hamiltonian.systemcell.Unitg3);

        Initialize(split<std::string>(kppath,"->"),kpnumber);
        kppos=symmetryPointsPositions[std::distance(m_path.begin(),std::find(m_path.begin(),m_path.end(),"G"))];

        if(id==0 && FirstLoop)
        {
            std::fstream file_obj; 
            file_obj.open("./Chi"+std::to_string(tmpName)+".dat",std::ios::out|std::ios::trunc);
            file_obj.close();
            file_obj.open("./log"+std::to_string(tmpName),std::ios::out|std::ios::trunc);
            file_obj.close();
            /*
               if(DataAppending==0)
               {
               for(int wl=LengthBegin;wl<=LengthEnd;++wl){
               file_obj.open(("./BandData/"+dataVersion+"0/TimeSpend"+std::to_string(CurrentMode)+"_"+std::to_string(wl)+".data").c_str(),std::ios::out|std::ios::trunc);
               file_obj.close();
               file_obj.open(("./BandData/"+dataVersion+"0/ReducedHamiltonian"+std::to_string(CurrentMode)+"_"+std::to_string(wl)+".data").c_str(),std::ios::out|std::ios::trunc);
               file_obj.close();
               }
               }
               */
            file_obj.open(("./out"+std::to_string(tmpName)+".txt").c_str(),std::ios::out|std::ios::app);
            file_obj<<"kvalue =";
            for(int i=0;i<m_path.size();++i)
            {
                file_obj<<" "<<symmetryPointsPositions[i];
            }
            file_obj<<std::endl;
            file_obj<<"pnumber"<<" = "<<bandend<<std::endl;
            file_obj<<"kpos"<<" = "<<kppos<<std::endl;
            file_obj.close();
        }

        if(FirstLoop)
        {
            if(DataAppending==0)
            {
                std::fstream file_obj; 
                file_obj.open("./Cdata.txt",std::ios::out|std::ios::trunc);
                file_obj.close();
                for(int fl=bandbegin;fl<=bandend;++fl){
                    for(int wl=LengthBegin;wl<=LengthEnd;++wl){
                        mkdir(("./BandData/"+dataVersion+std::to_string(fl)+"_"+std::to_string(wl)).c_str(),S_IRWXU|S_IRWXG|S_IRWXO);
                        file_obj.open(("./BandData/"+dataVersion+std::to_string(fl)+"_"+std::to_string(wl)+"/M"+std::to_string(CurrentMode)+"band"+std::to_string(id)+".data").c_str(),std::ios::out|std::ios::trunc);
                        file_obj.close();
                    }
                }
            }
        }


    }


    auto BandStructure::H_recyclePlus3(std::vector<int>&table,int &level) -> void
    {

        int id,numprocs;
        MPI_Comm_rank(MPI_COMM_WORLD,&id);
        MPI_Comm_size(MPI_COMM_WORLD,&numprocs);


        outputflag+=1;
        unsigned int lim_outputflag=kpoints_mesh.size();//qqq
        if(outputflag<lim_outputflag)
        {
            std::cout<<"outputflag="<<outputflag<<std::endl;
            //level=table.size()-2;//loop from begin

            int posNum=std::distance(std::find(table.begin(),table.end(),table[level]),std::find(table.begin(),table.end(),MeshBZ))+1;
            level+=posNum;
            timing.reset();
        }else{
            outputflag=0;
            timing.reset();
        }
    }


    auto BandStructure::H_recycle(std::vector<int>&table,int &level,int opTag1,int opTag2) -> void
    {

        int id,numprocs;
        MPI_Comm_rank(MPI_COMM_WORLD,&id);
        MPI_Comm_size(MPI_COMM_WORLD,&numprocs);


        int bandlen=bandend-bandbegin+1;
        int lengthlen=LengthEnd-LengthBegin+1;
        outputflag+=1;
        //if(USE_Infiniteloop && USE_CalcDOS) outputflag-=1;
        if(outputflag<lengthlen*bandlen)
        {
            //RunningFlag+=3;//loop in the middle
            //if(USE_Infiniteloop && USE_CalcDOS)   
            //{
            //int posNum=std::distance(std::find(table.begin(),table.end(),table[level]),std::find(table.begin(),table.end(),design_running_pattern))+1;
            //level+=posNum;//loop from begin
            level=table.size();//loop from begin
            //}else{
            //level=table.size();//loop from begin
            //}
            //elaps.push_back(timing.elapsed());

            if(id==0)
            {
                int DF=0;
                if(bandbegin!=bandend && LengthBegin==LengthEnd)
                {
                    DF=hamiltonian.systemcell.DividedFactor;
                }else if(bandbegin==bandend && LengthBegin!=LengthEnd)
                {
                    DF=pow(2.,hamiltonian.systemcell.DividedFactor-1);
                }else if(bandbegin!=bandend && LengthBegin!=LengthEnd)
                {
                    DF=hamiltonian.systemcell.DividedFactor;
                }

                std::fstream file_obj; 
                file_obj.open(("./BandData/"+dataVersion+"0/TimeSpend"+std::to_string(CurrentMode)+"-"+std::to_string(opTag1)+"_"+std::to_string(opTag2)+".data").c_str(),std::ios::out|std::ios::trunc);
                //file_obj<<(DF-opTag1)<<"\t"<<timing.elapsed()<<std::endl;
                file_obj<<(opTag1)<<"\t"<<timing.elapsed()<<std::endl;
                file_obj.close();
                file_obj.open(("./BandData/"+dataVersion+"0/ReducedHamiltonian"+std::to_string(CurrentMode)+"-"+std::to_string(opTag1)+"_"+std::to_string(opTag2)+".data").c_str(),std::ios::out|std::ios::trunc);
                file_obj<<(opTag1)<<"\t"<<hamiltonian.MatrixSizeReduced<<std::endl;
                file_obj.close();

                //file_obj.open(("./out"+std::to_string(tmpName)+".txt").c_str(),std::ios::out|std::ios::app);
                //file_obj<<(DF-opTag1)<<"\t"<<timing.elapsed()<<std::endl;
                //file_obj<<"TimeSpend"<<(opTag1)<<"_"<<(opTag2)<<" = "<<timing.elapsed()<<std::endl;
                //file_obj<<"ReducedHamiltonian"<<(opTag1)<<"_"<<(opTag2)<<" = "<<hamiltonian.MatrixSizeReduced<<std::endl;
                //file_obj.close();

            }

            timing.reset();
        }else{

            if(id==0)
            {
                int DF=0;
                if(bandbegin!=bandend && LengthBegin==LengthEnd)
                {
                    DF=hamiltonian.systemcell.DividedFactor;
                }else if(bandbegin==bandend && LengthBegin!=LengthEnd)
                {
                    DF=pow(2.,hamiltonian.systemcell.DividedFactor-1);
                }else if(bandbegin!=bandend && LengthBegin!=LengthEnd)
                {
                    DF=hamiltonian.systemcell.DividedFactor;
                }
                std::fstream file_obj; 
                file_obj.open(("./BandData/"+dataVersion+"0/TimeSpend"+std::to_string(CurrentMode)+"-"+std::to_string(opTag1)+"_"+std::to_string(opTag2)+".data").c_str(),std::ios::out|std::ios::trunc);
                //file_obj<<(DF-opTag1)<<"\t"<<timing.elapsed()<<std::endl;
                file_obj<<(opTag1)<<"\t"<<timing.elapsed()<<std::endl;
                file_obj.close();
                file_obj.open(("./BandData/"+dataVersion+"0/ReducedHamiltonian"+std::to_string(CurrentMode)+"-"+std::to_string(opTag1)+"_"+std::to_string(opTag2)+".data").c_str(),std::ios::out|std::ios::trunc);
                file_obj<<(opTag1)<<"\t"<<hamiltonian.MatrixSizeReduced<<std::endl;
                file_obj.close();

                //file_obj.open(("./out"+std::to_string(tmpName)+".txt").c_str(),std::ios::out|std::ios::app);
                //file_obj<<(DF-DataTag)<<"\t"<<timing.elapsed()<<std::endl;
                //file_obj<<"TimeSpend"<<(DataTag)<<"_"<<(opTag2)<<" = "<<timing.elapsed()<<std::endl;
                //file_obj<<"ReducedHamiltonian"<<(DataTag)<<"_"<<(opTag2)<<" = "<<hamiltonian.MatrixSizeReduced<<std::endl;
                //file_obj.close();
            }
            outputflag=0;

            timing.reset();

        }
    }


    auto BandStructure::NewOutput(int No,int startpoint,int fileNo) -> void{

        int id;
        MPI_Comm_rank(MPI_COMM_WORLD,&id);

        auto _minC=std::min_element(ConductionBand.begin(),ConductionBand.end());
        auto _maxV=std::max_element(ValenceBand.begin(),ValenceBand.end());
        double minC_in_id=(*_minC);
        double maxV_in_id=(*_maxV);
        double minC,maxV;
        MPI_Allreduce(&minC_in_id,&minC,1,MPI_DOUBLE,MPI_MIN,MPI_COMM_WORLD);
        MPI_Allreduce(&maxV_in_id,&maxV,1,MPI_DOUBLE,MPI_MAX,MPI_COMM_WORLD);
        double BandGapSearched=minC-maxV;
        /*
           char name[40];    
           if(id==0)
           {
           char dirname[40];    
           sprintf(dirname,"./BandData/%s%d",dataVersion.c_str(),fileNo);
           mkdir(dirname,S_IRWXU|S_IRWXG|S_IRWXO);
           }
           */
        auto maxPos=std::max_element(MatrixBottom.begin(),MatrixBottom.end());
        int maxValue_in_id=(*maxPos);
        int maxValue;
        MPI_Allreduce(&maxValue_in_id,&maxValue,1,MPI_INT,MPI_MAX,MPI_COMM_WORLD);

        //auto _minR=std::min_element(results.begin(),results.end(),[](std::vector<double> &i,std::vector<double> &j){return i.size()<j.size();});
        //int minResults_in_id=(*_minR).size();
        //int minResults;
        //MPI_Allreduce(&minResults_in_id,&minResults,1,MPI_INT,MPI_MIN,MPI_COMM_WORLD);

        if(id==0 && FirstLoop)
        {
            std::fstream file_obj; 
            file_obj.open(("./out"+std::to_string(tmpName)+".txt").c_str(),std::ios::out|std::ios::app);
            file_obj<<"GapLevel"<<fileNo<<" = "<<AtomNumbers*hamiltonian.BandGapIndex<<std::endl;
            file_obj.close();
            rename(("./out"+std::to_string(tmpName)+".txt").c_str(),("./BandData/"+dataVersion+"0/out"+std::to_string(fileNo)+".txt").c_str());
        }
        if(FirstLoop) FirstLoop=false;
        MatrixBottom.clear();

    }
    /*
       auto BandStructure::ComputePrepareMPI() -> void
       {
       int id,numprocs;
       MPI_Comm_rank(MPI_COMM_WORLD,&id);
       MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

       AtomNumbers=0;
       std::for_each(hamiltonian.systemcell.supercell.begin(),hamiltonian.systemcell.supercell.end(),
       [&atom=AtomNumbers]
       (auto &scell)
       {atom+=scell->AtomNum;}
       );
       if(hamiltonian.USE_decomposition)
       {
       AtomNumbers=hamiltonian.systemcell.UnitCellNumber;
       }

       nrPoints = GetPointsNumber();
       std::cout<<"correct="<<nrPoints<<"\twrong="<<qpoints.size()<<std::endl;

       if(USE_parallel_K_points)
       {
       interval=(nrPoints%numprocs)?(nrPoints/numprocs+1):(nrPoints/numprocs);
       startPoint = id*interval;
       cendPoint=std::min(startPoint + interval, nrPoints);
       }else{
       interval=nrPoints;
       startPoint = 0;
       cendPoint=nrPoints;
       }
       MatrixBottom.clear();
       if(cendPoint<=kppos) bandgapid=1;
       int bandid=0;
       MPI_Allreduce(&bandgapid,&bandid,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
       if(numprocs==1) bandid=0;
       if(id==0 && FirstLoop)
       {
       std::fstream file_obj; 
       file_obj.open(("./BandData/"+dataVersion+"0/out"+std::to_string(DataTag)+".txt").c_str(),std::ios::out|std::ios::app);
       file_obj<<"kpos_id"<<" = "<<bandid<<std::endl;
       file_obj.close();
       }
       if(FirstLoop) FirstLoop=false;
       }
       */


    auto BandStructure::ComputePrepareDFTB() -> void
    {
        int id,numprocs;
        MPI_Comm_rank(MPI_COMM_WORLD,&id);
        MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
        AtomNumbers=0;
        /*
           if(USE_CalcDOS) 
           {
           bandgapid=0;
           bandgap=0;
           ConductionBand.clear();
           ValenceBand.clear();
           BandValue.clear();
           kpoints.clear();
           results.clear();
           RandomKpointList(1);//here DOS
           }
           */
        std::for_each(hamiltonian.systemcell.supercell.begin(),hamiltonian.systemcell.supercell.end(),
                [&atom=AtomNumbers]
                (auto &scell)
                {atom+=scell->AtomNum;}
                );
        if(hamiltonian.USE_decomposition)
        {
            AtomNumbers=hamiltonian.systemcell.UnitCellNumber;
        }

        nrPoints = GetPointsNumber();
        if(USE_parallel_K_points)
        {
            interval=(nrPoints%numprocs)?(nrPoints/numprocs+1):(nrPoints/numprocs);
            startPoint = id*interval;
            cendPoint=std::min(startPoint + interval, nrPoints);
        }else{
            interval=nrPoints;
            startPoint = 0;
            cendPoint=nrPoints;
        }
        MatrixBottom.clear();

        if(cendPoint<=kppos) bandgapid=1;
        int bandid=0;
        if(USE_parallel_K_points)
        {
            MPI_Allreduce(&bandgapid,&bandid,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        }else{
            bandid=0;
        }
        if(numprocs==1) bandid=0;

        if(id==0 && FirstLoop)
        {

            std::fstream file_obj; 
            //file_obj.open(("./BandData/"+dataVersion+"0/out"+std::to_string(DataTag)+".txt").c_str(),std::ios::out|std::ios::app);
            file_obj.open(("./out"+std::to_string(tmpName)+".txt").c_str(),std::ios::out|std::ios::app);
            file_obj<<"kpos_id"<<" = "<<bandid<<std::endl;
            file_obj.close();

            if(USE_CalcDOS && USE_Infiniteloop) rename(("./out"+std::to_string(tmpName)+".txt").c_str(),("./BandData/"+dataVersion+"0/out"+std::to_string(DataTag)+".txt").c_str());//DOS
        }
        if(FirstLoop && USE_CalcDOS && USE_Infiniteloop) FirstLoop=false;//DOS

    }

    int BandStructure::Monkhorstpack(int r,int q){
        return 2*r-q-1;
    }

    void BandStructure::RandomKpointList(int q){
        int id,numprocs;
        MPI_Comm_rank(MPI_COMM_WORLD,&id);
        MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

        mesh_factor=1.0/(2*q);
        kpoints.resize(0);
        int o=0;
        /*
           std::random_device rd;
           std::default_random_engine engDevice(rd());
           */
        int seed = 0;
        //if(id==0) seed=distr(engDevice);
        if(id==0) seed=std::chrono::system_clock::now().time_since_epoch().count();
        MPI_Bcast(&seed,1,MPI_INT,0,MPI_COMM_WORLD);
        //srand(seed);
        std::mt19937 eng (seed);
        std::uniform_real_distribution<double> distr(0,1);
        //std::chrono::system_clock::now().time_since_epoch().count();

        for(int kz=1;kz<=q;++kz){
            for(int ky=1;ky<=q;++ky){
                for(int kx=1;kx<=q;++kx){
                    //Eigen::Vector3d tmp=hamiltonian.systemcell.Unitg1*Monkhorstpack(kx,q)+hamiltonian.systemcell.Unitg2*Monkhorstpack(ky,q)+hamiltonian.systemcell.Unitg3*Monkhorstpack(kz,q);
                    //tmp=tmp*mesh_factor;

                    Eigen::Vector3d tmp=hamiltonian.systemcell.Unitg1*distr(eng)+hamiltonian.systemcell.Unitg2*distr(eng)+hamiltonian.systemcell.Unitg3*distr(eng);
                    //Eigen::Vector3d tmp=hamiltonian.systemcell.Unitg1*distr(engDevice)+hamiltonian.systemcell.Unitg2*distr(engDevice)+hamiltonian.systemcell.Unitg3*distr(engDevice);
                    kpoints.push_back(tmp);
                    o+=1;
                }
            }
        }
        integrate_factor=hamiltonian.systemcell.Unitg1.dot(hamiltonian.systemcell.Unitg2.cross(hamiltonian.systemcell.Unitg3))/pow(q*1.,3.);
    }


    void BandStructure::Generate_mesh(int q){
        mesh_factor=1.0/(2*q);
        //kpoints_mesh.resize(pow(q,3));
        //kpoints_degen.resize(pow(q,3));
        kpoints.resize(0);
        kpoints_degen.resize(pow(q,3),0);
        int o=0;
        for(int kz=1;kz<=q;++kz){
            for(int ky=1;ky<=q;++ky){
                for(int kx=1;kx<=q;++kx){
                    Eigen::Vector3d tmp=hamiltonian.systemcell.Unitg1*Monkhorstpack(kx,q)+hamiltonian.systemcell.Unitg2*Monkhorstpack(ky,q)+hamiltonian.systemcell.Unitg3*Monkhorstpack(kz,q);
                    tmp=tmp*mesh_factor;
                    kpoints.push_back(tmp);
                    //printf("[%d]={%2.1f,%2.1f,%2.1f}\n",o+1,tmp.X,tmp.Y,tmp.Z);
                    o+=1;
                }
            }
        }
        integrate_factor=hamiltonian.systemcell.Unitg1.dot(hamiltonian.systemcell.Unitg2.cross(hamiltonian.systemcell.Unitg3))/pow(q*1.,3.);
        //integrate_factor=0;
        //std::cout<<"iF="<<integrate_factor<<std::endl;
    }

    auto BandStructure::H_loop_k_points(int backwordNum,int &level) -> void
    {

        int id,numprocs;
        MPI_Comm_rank(MPI_COMM_WORLD,&id);
        MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

        computeflag+=1;

        if(USE_Infiniteloop && USE_CalcDOS)
        {
            computeflag-=1;
            RandomKpointList(1);
        }
        if(computeflag<cendPoint-startPoint)
        {
            level+=backwordNum;
        }else{
            computeflag=0;
        }

    }

    auto BandStructure::OutPutStatus(int flag) -> void
    {

        int id,numprocs;
        MPI_Comm_rank(MPI_COMM_WORLD,&id);
        MPI_Comm_size(MPI_COMM_WORLD,&numprocs);

        //if(id==0 || flag) fprintf(stderr,"    H_setMatrix:i=%d  {%1.3f,%1.3f,%1.3f} %d\n",startPoint+computeflag,kpoints[startPoint+computeflag](0),kpoints[startPoint+computeflag](1),kpoints[startPoint+computeflag](2),hamiltonian.MatrixSize);
        /*
           if(id==0 && numprocs!=1) 
           {
           double ScaleFactor=0.4;
           int Pos=startPoint+computeflag+1;
           std::cout <<"\r["<<CurrentMode<<"]Calculating:("<<thrd+1<<"/"<<hamiltonian.systemcell.DividedFactor+1<<")["<<hamiltonian.MatrixSize<<"] [";
           for(int i=0;i<nrPoints*ScaleFactor;++i)
           {
           int iflag=i/numprocs;
           if(iflag<Pos*ScaleFactor)
           {
           std::cout <<"|";
           }else{
           std::cout <<" ";
           }
           }
           std::cout <<"]-"<<Pos*100./interval<<"%    ";
           }
           else if(numprocs==1)
           {
           double ScaleFactor=10;
           int Pos=outputflag;
           std::cout <<"\r["<<CurrentMode<<"]Calculating:["<<hamiltonian.MatrixSize<<"] [";
           for(int i=0;i<hamiltonian.systemcell.DividedFactor*ScaleFactor;++i)
           {
           if(i<Pos*ScaleFactor)
           {
           std::cout <<"|";
           }else{
           std::cout <<" ";
           }
           }
           std::cout <<"]-"<<(Pos)*100./hamiltonian.systemcell.DividedFactor<<"%    ";
           }
           */
    }

    void BandStructure::Output_Epsilon()
    {
        double dimX=Xmatrix.cols();
        //Eigen::Vector3d q=kpoints_mesh[outputflag];//qqq
        Eigen::Vector3d q=Eigen::Vector3d(0,0,0);//qqq
        hamiltonian.GetUq(q);
        /*
           for(int i=0;i<Xmatrix.rows();++i){
           for(int j=0;j<Xmatrix.cols();++j){
           Xmatrix(i,j)*=hamiltonian.Uq(i,j);
           }
           }
           */
        //Xmatrix=Eigen::MatrixXcd::Identity(dimX,dimX);//-/*hamiltonian.Uq*/Xmatrix/*4.*M_PI/(pow(q.normalized()+0.1,2.))*/;
        Xmatrix=Eigen::MatrixXcd::Identity(dimX,dimX)-hamiltonian.Uq*Xmatrix/*4.*M_PI/(pow(q.normalized()+0.1,2.))*/;
        //sumX=Xmatrix.sum()*hamiltonian.systemcell.Unitg1.dot(hamiltonian.systemcell.Unitg2.cross(hamiltonian.systemcell.Unitg3))/integrate_factor;
        sumX=Xmatrix.sum()*integrate_factor;
        //fprintf(Chi,"%1.2f\t%1.3f\t\t%1.3f\t%1.2f\t%1.2f\n",q.normalized(),sumX.real(),sumX.imag(),hamiltonian.g1*hamiltonian.g2.CrossProduct(hamiltonian.g3),hamiltonian.Unita1*hamiltonian.Unita2.CrossProduct(hamiltonian.Unita3));
        std::fstream file_obj; 
        file_obj.open("./Chi"+std::to_string(tmpName)+".dat",std::ios::out|std::ios::app);
        //if(q.norm()<=1.) fprintf(stderr,"%1.2f\t%1.3f\t\t%1.3f\t%1.2f\t%1.2f\n",q.norm(),sumX.real(),sumX.imag(),hamiltonian.systemcell.g1.dot(hamiltonian.systemcell.g2.cross(hamiltonian.systemcell.g3)),hamiltonian.systemcell.Unita1.dot(hamiltonian.systemcell.Unita2.cross(hamiltonian.systemcell.Unita3)));
        if(q.norm()<=1.)
        {
            //file_obj<<q.norm()<<"\t"<<sumX.real()<<"\t"<<sumX.imag()<<"\t"<<hamiltonian.systemcell.g1.dot(hamiltonian.systemcell.g2.cross(hamiltonian.systemcell.g3))<<"\t"<<hamiltonian.systemcell.Unita1.dot(hamiltonian.systemcell.Unita2.cross(hamiltonian.systemcell.Unita3))<<std::endl;
            file_obj<<kpoints_mesh[outputflag](0)<<"\t"<<sumX.real()<<"\t"<<sumX.imag()<<"\t"<<hamiltonian.systemcell.g1.dot(hamiltonian.systemcell.g2.cross(hamiltonian.systemcell.g3))<<"\t"<<hamiltonian.systemcell.Unita1.dot(hamiltonian.systemcell.Unita2.cross(hamiltonian.systemcell.Unita3))<<std::endl;///qqq
        }
        file_obj.close();

    }
    auto BandStructure::H_clean() -> void
    {
        int id,numprocs;
        MPI_Comm_rank(MPI_COMM_WORLD,&id);
        MPI_Comm_size(MPI_COMM_WORLD,&numprocs);


        hamiltonian.destroy();   

        hamiltonian.systemcell.spNum.resize(0);
        hamiltonian.systemcell.spStr.resize(0);
        //int bandid=0;
        //MPI_Allreduce(&bandgapid,&bandid,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD);
        //if(numprocs==1) bandid=0;
        //MPI_Bcast(&bandgap,1,MPI_DOUBLE,bandid,MPI_COMM_WORLD);
    }

}
