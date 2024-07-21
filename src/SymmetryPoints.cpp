#include "SymmetryPoints.hpp"

namespace TightBinding
{

    void SymmetryPoints::GenSymmetryPoints(Eigen::Vector3d &Ug1,Eigen::Vector3d &Ug2,Eigen::Vector3d &Ug3,Eigen::Vector3d &pg1,Eigen::Vector3d &pg2,Eigen::Vector3d &pg3)
    {
/*
        std::cout<<"pg1={"<<pg1(0)<<","<<pg1(1)<<","<<pg1(2)<<"}"<<std::endl;
        std::cout<<"pg2={"<<pg2(0)<<","<<pg2(1)<<","<<pg2(2)<<"}"<<std::endl;
        std::cout<<"pg3={"<<pg3(0)<<","<<pg3(1)<<","<<pg3(2)<<"}"<<std::endl;

        std::cout<<"Ug1={"<<Ug1(0)<<","<<Ug1(1)<<","<<Ug1(2)<<"}"<<std::endl;
        std::cout<<"Ug2={"<<Ug2(0)<<","<<Ug2(1)<<","<<Ug2(2)<<"}"<<std::endl;
        std::cout<<"Ug3={"<<Ug3(0)<<","<<Ug3(1)<<","<<Ug3(2)<<"}"<<std::endl;
*/
        symmetryPoints["E1"]=SymmetryPoint("E1", 
                //Eigen::Vector3d(-15,0.0,0.0)
                //Eigen::Vector3d(0,0.0,0.0)
                Eigen::Vector3d(1.,0.0,0.0)
                );
        symmetryPoints["E2"]=SymmetryPoint("E2", 
                //Eigen::Vector3d(10,0.0,0.0)
                //Eigen::Vector3d(6.5,0.0,0.0)
                Eigen::Vector3d(4,0.0,0.0)
                );
        symmetryPoints["XX"]=SymmetryPoint("XX", 
                Eigen::Vector3d(1.0,0.0,0.0).normalized()
                );
	symmetryPoints["YY"]=SymmetryPoint("YY", 
                Eigen::Vector3d(1,1,0.0).normalized()
                );
	symmetryPoints["ZZ"]=SymmetryPoint("ZZ", 
                Eigen::Vector3d(1.,1.,1.).normalized()
                );

        symmetryPoints["G"]=SymmetryPoint("G", 
                    0.*Ug1
                    +0.0*Ug2
                    +0.*Ug3
                );

        symmetryPoints["R"]=SymmetryPoint("R", 
                0.5*Ug1
                +0.5*Ug2
                +0.*Ug3
                );
        symmetryPoints["M"]=SymmetryPoint("M", 
                0.*Ug1
                +0.5*Ug2
                +0.5*Ug3
                );
        symmetryPoints["Z"]=SymmetryPoint("Z", 
                0.5*Ug1
                +0.*Ug2
                +0.*Ug3
                );
        symmetryPoints["Xbar"]=SymmetryPoint("Xbar", 
                0.*Ug1
                +0.5*Ug2
                +0.*Ug3
                );
        symmetryPoints["A"]=SymmetryPoint("A", 
                0.5*Ug1
                +0.5*Ug2
                +0.5*Ug3
                );

        ///////////////////////////////////
        
        symmetryPoints["(1,0,0)"]=SymmetryPoint("(1,0,0)", 
                //pgI
                pg3
                );

        Eigen::Vector3d pgI(1,0,0);
        symmetryPoints["(0,1,1)/2"]=SymmetryPoint("(0,1,1)/2", 
                //pgI
                0.0*pg1
                +0.5*pg2
                +0.5*pg3
                );

        pgI=Eigen::Vector3d(0.5,0.5,0.5);
        symmetryPoints["(1,1,1)/2"]=SymmetryPoint("(1,1,1)/2", 
                //pgI
                0.5*pg1
                +0.5*pg2
                +0.5*pg3
                );

        symmetryPoints["X"]=SymmetryPoint("X", 
                0.0*pg1
                +0.5*pg2
                +0.5*pg3
                );

        symmetryPoints["L"]=SymmetryPoint("L", 
                0.5*pg1
                +0.5*pg2
                +0.5*pg3
                );

        symmetryPoints["W"]=SymmetryPoint("W", 
                0.25*pg1
                +0.5*pg2
                +0.75*pg3
                );

        symmetryPoints["U"]=SymmetryPoint("U", 
                (1./4)*pg1
                +(5./8)*pg2
                +(5./8)*pg3
                );

        symmetryPoints["K"]=SymmetryPoint("K", 
                (3./8)*pg1
                +(3./8)*pg2
                +(3./4)*pg3
                );

        symmetryPoints["X'"]=SymmetryPoint("X'", 
                (1./2)*pg1
                +(1./2)*pg2
                +(1.)*pg3
                );


    }

    SymmetryPoints::SymmetryPoints()
    {

        //symmetryPoints["-(0.25,0.25,0)"] = SymmetryPoint("-IK", Eigen::Vector3d(-0.25, -0.25, 0));
        //symmetryPoints["(0.25,0.25,0)"] = SymmetryPoint("IK", Eigen::Vector3d(0.25, 0.25, 0));
        //symmetryPoints["-I"] = SymmetryPoint("-I", Eigen::Vector3d(-0.25, 0, 0));
        //symmetryPoints["I"] = SymmetryPoint("I", Eigen::Vector3d(0.25, 0, 0));

        symmetryPoints["L"] = SymmetryPoint("L", Eigen::Vector3d(0.5, 0.5, 0.5));
        symmetryPoints["G"] = SymmetryPoint("G", Eigen::Vector3d(0., 0., 0.));
        symmetryPoints["X"] = SymmetryPoint("X", Eigen::Vector3d(1., 0., 0.));
        symmetryPoints["W"] = SymmetryPoint("W", Eigen::Vector3d(1., 0.5, 0.));
        symmetryPoints["K"] = SymmetryPoint("K", Eigen::Vector3d(0.75, 0.75, 0.));
        symmetryPoints["U"] = SymmetryPoint("U", Eigen::Vector3d(1., 0.25, 0.25));
        symmetryPoints["UK"] = SymmetryPoint("UK", Eigen::Vector3d(1., 1., 0.));

        symmetryPoints["A"] = SymmetryPoint("A", Eigen::Vector3d(0.5,0.5,0.5));
        symmetryPoints["R"] = SymmetryPoint("R", 0.5*Eigen::Vector3d(0.,1.,1.));
        symmetryPoints["M"] = SymmetryPoint("M", 0.5*Eigen::Vector3d(1.,1.,0.));
        symmetryPoints["Xbar"] = SymmetryPoint("Xbar", 0.5*Eigen::Vector3d(0.,1.,0.));
        symmetryPoints["Z"] = SymmetryPoint("Z", 0.5*Eigen::Vector3d(0.,0.,1));


    }




    std::vector<Eigen::Vector3d> SymmetryPoints::GeneratePoints(const std::vector<std::string>& path, int kpnumber, std::vector<unsigned int>& symmetryPointsPositions)
    {
        std::vector<Eigen::Vector3d> result;

        symmetryPointsPositions.clear();
        symmetryPointsPositions.reserve(path.size());

        //if (kpnumber <= path.size() * 2 + 1) return result;
        if (kpnumber < 1) return result;

        result.reserve(kpnumber);

        double length = 0;

        for (int i = 1; i < path.size(); ++i)
        {
            const Eigen::Vector3d dif = symmetryPoints[path[i]].position - symmetryPoints[path[i - 1]].position;
            length += dif.norm();
        }

        const double stepSize = length / (kpnumber - 1);

        for (int i = 1; i < path.size(); ++i)
        {
            const Eigen::Vector3d startPos = symmetryPoints[path[i - 1]].position;
            const Eigen::Vector3d dif = symmetryPoints[path[i]].position - startPos;
            const double difLength = dif.norm();

            const Eigen::Vector3d stepVec = dif / difLength * stepSize;

            symmetryPointsPositions.push_back(static_cast<unsigned int>(result.size()));

            for (Eigen::Vector3d pos = startPos; (pos - startPos).norm() < difLength; pos += stepVec)
                result.push_back(pos);

            if(path.size() == (i+1))
            {
                symmetryPointsPositions.push_back(static_cast<unsigned int>(result.size()));
                result.push_back(symmetryPoints[path[i]].position);
            }

        }
        return std::move(result);
    }

}
