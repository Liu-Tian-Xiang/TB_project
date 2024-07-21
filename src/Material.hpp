#include "Util.hpp"
#ifndef TB_Mat_H
#define TB_Mat_H

namespace TightBinding
{

    class Material
    {
	public:
	    //Material(){}
		Material() : m_Esa(0), m_Epa(0), m_Esc(0), m_Epc(0), m_Esstara(0), m_Esstarc(0), m_Vss(0), m_Vxx(0), m_Vxy(0), m_Vsapc(0), m_Vscpa(0), m_Vsstarapc(0), m_Vpasstarc(0) {}
		Material(const std::string& Name, double Esa, double Epa, double Esc, double Epc, double Esstara, double Esstarc, double Vss, double Vxx, double Vxy, double Vsapc, double Vscpa, double Vsstarapc, double Vpasstarc,double lambdaA,double lambdaC);
	    Material(const std::string& Name, double Esa, double Esc, double Epa, double Epc, double Ed, double Es_, double SS, double S_S_,double S_aSc,double SaS_c,double  SaPcSIGMA,double  ScPaSIGMA,double  S_aPcSIGMA,double  S_cPaSIGMA ,double SaDcSIGMA,double ScDaSIGMA,double S_aDcSIGMA,double S_cDaSIGMA,double PPSIGMA,double PPPI,double PaDcSIGMA,double PcDaSIGMA,double PaDcPI,double PcDaPI,double DDSIGMA,double DDPI,double DDDELTA ,double lambdaA,double lambdaC);
	    Material(const std::string &NameStr):name(NameStr){}
	    std::string name;
		double m_Esa;
		double m_Epa;
		double m_Esc;
		double m_Epc;
		double m_Esstara;
		double m_Esstarc;
		double m_Vsapc;
		double m_Vscpa;
		double m_Vsstarapc;
		double m_Vpasstarc;
		double m_Vss;

		double m_Vxx;
		double m_Vxy;
        std::unordered_map<std::string,double> param;
    };


    class Materials
    {
	public:

	    Materials();//{};

	    void MixTwoMaterials(const std::string &mats1,const std::string &mats2,double coef_mats1,const std::string &NewMtl);
	    void MixTwoMaterialsShort(const std::string &mats1,const std::string &mats2,double coef_mats1,const std::string &NewMtl);
	    void ReadMaterials(const std::string &fileName,const std::string &mtl);
		std::unordered_map<std::string, Material> materialsShort;
	    std::unordered_map<std::string, Material> materials;
    };


}
#endif
