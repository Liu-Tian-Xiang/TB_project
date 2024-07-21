#include "Util.hpp"
#ifndef TB_SYM_H
#define TB_SYM_H

namespace TightBinding
{

	class SymmetryPoint
	{
	public:
		SymmetryPoint() {}

		SymmetryPoint(const std::string& Name, const Eigen::Vector3d pos) : name(Name), position(pos) {}

		SymmetryPoint(const SymmetryPoint& sym)
			: name(sym.name), position(sym.position)
		{
		}

		SymmetryPoint& operator=(const SymmetryPoint& sym)
		{
			name = sym.name;
			position = sym.position;

			return *this;
		}

		std::string name;
        Eigen::Vector3d position;
	};


	class SymmetryPoints
	{
	public:
		SymmetryPoints();

        void GenSymmetryPoints(Eigen::Vector3d &Ug1,Eigen::Vector3d &Ug2,Eigen::Vector3d &Ug3,Eigen::Vector3d &g1,Eigen::Vector3d &g2,Eigen::Vector3d &g3);
		std::map<std::string, SymmetryPoint> symmetryPoints;
		std::vector<Eigen::Vector3d> GeneratePoints(const std::vector<std::string>& path, int nrPoints, std::vector<unsigned int>& symmetryPointsPositions);
	};

}

#endif

