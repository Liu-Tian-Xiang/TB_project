#include "Material.hpp"


namespace TightBinding
{
	Material::Material(const std::string& Name, double Esa, double Epa, double Esc, double Epc, double Esstara, double Esstarc, double Vss, double Vxx, double Vxy, double Vsapc, double Vscpa, double Vsstarapc, double Vpasstarc,double lambdaA,double lambdaC)
		: name(Name), m_Esa(Esa), m_Epa(Epa), m_Esc(Esc), m_Epc(Epc), m_Esstara(Esstara), m_Esstarc(Esstarc), m_Vss(Vss), m_Vxx(Vxx), m_Vxy(Vxy), m_Vsapc(Vsapc), m_Vscpa(Vscpa), m_Vsstarapc(Vsstarapc), m_Vpasstarc(Vpasstarc)
    {

        param[ "00a00a" ]= Esa;
        param[ "00c00c" ]= Esc;
        param[ "01a01a" ]= Epa;
        param[ "01c01c" ]= Epc;

        param[ "10a10a" ]= Esstara;
        param[ "10c10c" ]= Esstarc;
        param[ "00a00c" ]= Vss*0.25;
        param[ "00c00a" ]= Vss*0.25;

        param[ "10a10c" ]= 0;
        param[ "10c10a" ]= 0;

        param[ "00a01c" ]= Vsapc*(sqrt(3))*0.25;
        param[ "00c01a" ]= Vscpa*(sqrt(3))*0.25;
        param[ "10a01c" ]= Vsstarapc*(sqrt(3))*0.25;
        param[ "10c01a" ]= Vpasstarc*(sqrt(3))*0.25;

        param[ "01a01c" ]=(Vxx+2*Vxy)*0.25;
        param[ "01c01a" ]=(Vxx+2*Vxy)*0.25;
        param[ "01a01cpi" ]= (Vxx-Vxy)*0.25;
        param[ "01c01api" ]= (Vxx-Vxy)*0.25;
        param[ "10a00c" ]= 0;
        param[ "00c10a" ]= 0;
        param[ "00a10c" ]= 0;
        param[ "10c00a" ]= 0;
       
        param[ "a" ]= lambdaA;
        param[ "c" ]= lambdaC;
    }      

    Material::Material(const std::string& Name, double Esa, double Esc, double Epa, double Epc, double Ed, double Es_, double SS, double S_S_,double S_aSc,double SaS_c,double  SaPcSIGMA,double  ScPaSIGMA,double  S_aPcSIGMA,double  S_cPaSIGMA ,double SaDcSIGMA,double ScDaSIGMA,double S_aDcSIGMA,double S_cDaSIGMA,double PPSIGMA,double PPPI,double PaDcSIGMA,double PcDaSIGMA,double PaDcPI,double PcDaPI,double DDSIGMA,double DDPI,double DDDELTA ,double lambdaA,double lambdaC)
	: name(Name)
    {

	param[ "00a00a" ]= Esa;
	param[ "00c00c" ]= Esc;
	param[ "01a01a" ]= Epa;
	param[ "01c01c" ]= Epc;
	param[ "10a10a" ]= Es_;
	param[ "10c10c" ]= Es_;
	param[ "02a02a" ]= Ed;
	param[ "02c02c" ]= Ed;

	param[ "00a00c" ]= SS;
	param[ "00c00a" ]= SS;
	param[ "00a01c" ]= SaPcSIGMA;
	param[ "00c01a" ]= ScPaSIGMA;
	param[ "01a01c" ]=PPSIGMA;
	param[ "01c01a" ]=PPSIGMA;
	param[ "01a01cpi" ]= PPPI;
	param[ "01c01api" ]= PPPI;
	param[ "00a02c" ]= SaDcSIGMA;
	param[ "00c02a" ]= ScDaSIGMA;
    param[ "01a02c" ]= PaDcSIGMA;
    param[ "01c02a" ]= PcDaSIGMA;
    param[ "01a02cpi" ]= PaDcPI;
    param[ "01c02api" ]= PcDaPI;
    param[ "02a02c" ]= DDSIGMA;
    param[ "02c02a" ]= DDSIGMA;
    param[ "02a02cpi" ]= DDPI;
    param[ "02c02api" ]= DDPI;
    param[ "02a02cdelta" ]= DDDELTA;
    param[ "02c02adelta" ]= DDDELTA;
   	param[ "10a10c" ]= S_S_;
	param[ "10c10a" ]= S_S_;
	param[ "10a00c" ]= S_aSc;
	param[ "00c10a" ]= S_aSc;
	param[ "00a10c" ]= SaS_c;
	param[ "10c00a" ]= SaS_c;
	param[ "10a01c" ]= S_aPcSIGMA;
	param[ "10c01a" ]= S_cPaSIGMA;
 	param[ "10a02c" ]= S_aDcSIGMA;
 	param[ "10c02a" ]= S_cDaSIGMA;
 	param[ "a" ]= lambdaA;
 	param[ "c" ]= lambdaC;

    }

    Materials::Materials()
    {
        materialsShort["C"] = Material("C",  -4.545,  3.84, -4.545, 3.84, 11.37, 11.37, -22.7250, 3.8400, 11.6700, 15.2206, 15.2206, 8.2109, 8.2109,0,0);
		materialsShort["Si"] = Material("Si", -4.2, 1.715, -4.2, 1.715, 6.685, 6.685, -8.3, 1.715, 4.575, 5.7292, 5.7292, 5.3749, 5.3749,0,0);
		materialsShort["Ge"] = Material("Ge", -5.88, 1.61, -5.88, 1.61, 6.39, 6.39, -6.78, 1.61, 4.9, 5.4649, 5.4649, 5.2191, 5.2191,0,0);
		materialsShort["Sn"] = Material("Sn", -5.67, 1.33, -5.67, 1.33, 5.9, 5.9, -5.67, 1.33, 4.08, 4.5116, 4.5116, 5.8939, 5.8939,0,0);
		materialsShort["SiC"] = Material("SiC", -8.4537, 2.1234, -4.8463, 4.3466, 9.6534, 9.3166, -12.4197, 3.038, 5.9216, 9.49, 9.2007, 8.7138, 4.4051,0,0);
		materialsShort["AlP"] = Material("AlP", -7.8466, 1.3169, -1.2534, 4.2831, 8.7069, 7.4231, -7.4535, 2.3749, 4.8378, 5.2451, 5.2775, 5.2508, 4.6388,0.067,0.024);
		//materialsShort["AlAs"] = Material("AlAs", -7.5273, 0.9833, -1.1627, 3.5867, 7.4833, 6.7267, -6.6642, 1.878, 4.2919, 5.1106, 5.4965, 4.5216, 4.995,0.421,0.024);
		materialsShort["AlAs"] = Material("AlAs", -8.2663110, 0.3442887, -1.6298230, 2.9476890, 6.8442390, 6.0876890, -6.6642, 1.878, 3.8600, 5.6000, 7.6000, 4.2200, 8.3000,0.421,0.024);//modified by Boykin et.al band shift.
		//materialsShort["AlAs"] = Material("AlAs", -3.21537, -0.09711, -9.52462, 4.97139, 12.05550, 3.99445, -8.84261, -0.01434, 4.25949, 2.42476, 13.20317, 5.83246, 4.60075,0.29145,0.03152);//modified by Klimeck et.al band shift.
		//materialsShort["AlAs"] = Material("AlAs", -8.266310, 0.344290, -1.782020, 2.947690, 6.845424, 6.087690, -6.6642, 1.878, 3.8600, 5.6000, 7.6000, 4.2200, 8.3000,0.421,0.024);//modified by Boykin et.al without band shift.
		//materialsShort["AlSb"] = Material("AlSb", -6.1714, 0.9807, -2.0716, 3.0163, 6.7607, 6.1543, -5.6448, 1.7199, 3.6648, 4.9121, 4.2137, 4.3662, 3.0739,0,0);
		materialsShort["AlSb"] = Material("AlSb", -5.24996, 1.10214, -1.65016, 3.13774, 6.88214, 6.27574, -5.6648, 1.7199, 3.6648, 5.5000, 6.2137, 5.3000, 5.2739,0,0);//Timothy Type
		materialsShort["GaP"] = Material("GaP", -8.1124, 1.125, -2.1976, 4.115, 8.515, 7.185, -7.4709, 2.1516, 5.1369, 4.2771, 6.319, 4.6541, 5.095,0.067,0.174);
		//materialsShort["GaAs"] = Material("GaAs", -8.3431, 1.0414, -2.6569, 3.6686, 8.5914, 6.7386, -6.4513, 1.9546, 5.0779, 4.48, 5.7839, 4.8422, 4.8077,0.421,0.174);
		materialsShort["GaAs"] = Material("GaAs", -8.39000, 1.07475, -2.65405, 3.55475, 8.57475, 6.70475, -6.4513, 1.9546, 4.7700, 4.6800, 7.7000, 4.8500, 6.9000,0.421,0.174);//modified by Timothy et.al.
		//materialsShort["GaAs"] = Material("GaAs", -3.53284, 0.27772, -8.11499, 4.57341, 12.33930, 4.31241, -6.87653, 1.33572, 5.07596, 2.85929, 11.09774, 6.31619, 5.02335,0.32703,0.12000);//modified by Klimeck et.al.
		materialsShort["GaSb"] = Material("GaSb", -7.3207, 0.8554, -3.8993, 2.9146, 6.6354, 5.9846, -6.1567, 1.5789, 4.1285, 4.9601, 4.6675, 4.9895, 4.218,0.973,0.179);
		materialsShort["InP"] = Material("InP", -8.5274, 0.8735, -1.4826, 4.0465, 8.2635, 7.0665, -5.3614, 1.8801, 4.2324, 2.2265, 5.5825, 3.4623, 4.4814,0.067,0.392);
		//materialsShort["InAs"] = Material("InAs", -9.5381, 0.9099, -2.7219, 3.7201, 7.4099, 6.7401, -5.6052, 1.8398, 4.4693, 3.0354, 5.4389, 3.3744, 3.9097,0.421,0.392);
		materialsShort["InAs"] = Material("InAs", -9.6081, 0.9099, -2.5519, 3.7201, 7.4099, 6.7401, -5.4052, 1.8398, 4.4693, 3.3054, 5.4389, 3.3744, 3.9097,0.421,0.392);//Timothy type
		materialsShort["InSb"] = Material("InSb", -8.0157, 0.6738, -3.4643, 2.9162, 6.453, 5.9362, -5.5193, 1.4018, 3.8761, 3.788, 4.59, 3.5666, 3.4048,0.973,0.392);
		materialsShort["ZnSe"] = Material("ZnSe", -11.8383, 1.5072, 0.0183, 5.9928, 7.5872, 8.9928, -6.2163, 3.0054, 5.9942, 3.498, 6.3191, 2.5891, 3.9533,0.48,0.074);
		materialsShort["ZnTe"] = Material("ZnTe", -9.8150, 1.4834, 0.935, 5.2666, 7.0834, 8.2666, -6.5765, 2.7951, 5.467, 5.9827, 5.8199, 1.3196, 0.,0,0);
	// complete with all compounds params from Table 1 from P. Vogl et al.
	materials["Si"] = Material("Si", -2.0196,-2.0196,4.5448,4.5448,14.1836,19.6748,-1.9413,-3.3081,-1.6933,-1.6933,2.7836,2.7836,2.8428,2.8428,-2.7998,-2.7998,-0.7003,-0.7003,4.1068,-1.5934,-2.1073,-2.1073,1.9977,1.9977,-1.2327,2.5145,-2.4734,0.0195,0.0195);
	materials["C"] = Material("C", -1.0458,-1.0458,7.0850,7.0850,27.9267,38.2661,-4.3882,-2.6737,-2.3899,-2.3899,5.4951,5.4951,5.1709,5.1709,-2.7655,-2.7655,-2.3034,-2.3034,7.5480,-2.6363,-2.1621,-2.1621,3.9281,3.9281,-4.1813,4.9779,-3.9884,0,0);
	materials["AlP"] = Material("AlP", -5.3355,0.9573,3.3471,6.3392,14.1717,20.5963,-1.7403,-3.6444,-1.6448,-1.4307,2.6146,2.7804,2.0623,2.3361,-2.5253,-2.1687,-0.7810,-0.7211,4.0355,-1.3077,-1.6750,-1.8239,1.8760,2.1848,-1.3479,2.3750,-1.8464,0.0196,0.0073);
	materials["GaP"] = Material("GaP", -5.3379,-0.4005,3.3453,6.3844,14.0431,20.3952,-1.7049,-3.5704,-1.6034,-1.6358,2.8074,2.9800,2.3886,2.1482,-2.7840,-2.3143,-0.6426,-0.6589,4.1988,-1.4340,-1.7911,-1.8106,1.8574,2.1308,-1.2268,2.2752,-2.0124,0.0301,0.0408);
	materials["InP"] = Material("InP", -5.3321,0.3339,3.3447,6.4965,12.7756,18.8738,-1.4010,-3.6898,-1.8450,-1.2867,2.1660,2.6440,2.5652,2.0521,-2.5559,-2.2192,-0.7912,-0.8166,4.0203,-1.2807,-1.9239,-1.8851,1.5679,1.7763,-1.2482,2.1487,-1.6857,0.0228,0.1124);
	materials["AlAs"] = Material("AlAs", -5.9819,0.9574,3.5826,6.3386,13.0570,19.5133,-1.7292,-3.6094,-1.6167,-1.2688,2.5175,2.7435,2.1190,2.1989,-2.5535,-2.3869,-0.8064,-0.7442,4.2460,-1.3398,-1.7240,-1.7601,1.7776,2.0928,-1.2175,2.1693,-1.7540,0.1721,0.0072);
	materials["GaAs"] = Material("GaAs", -5.9819,-0.4028,3.5820,6.3853,13.1023,19.4220,-1.6187,-3.6761,-1.9927,-1.5648,2.4912,2.9382,2.1835,2.2086,-2.7333,-2.4095,-0.6906,-0.6486,4.4094,-1.4572,-1.7811,-1.8002,1.7821,2.0709,-1.1409,2.2030,-1.9770,0.1824,0.0408);
	materials["InAs"] = Material("InAs", -5.9801,0.3333,3.5813,6.4939,12.1954,17.8411,-1.4789,-3.8514,-2.1320,-1.2219,2.3159,2.8006,2.6467,1.9012,-2.5828,-2.4499,-0.8497,-0.8371,4.1188,-1.3687,-2.1222,-2.0584,1.5462,1.7106,-1.2009,2.1820,-1.7788,0.1763,0.1248);
	materials["AlSb"] = Material("AlSb", -4.9565,0.9521,4.0739,6.3386,11.4691,16.4173,-1.6179,-3.3145,-1.6983,-1.2097,2.5918,2.9334,2.4649,1.8889,-2.7920,-2.0008,-0.7307,-0.7878,4.1042,-1.5273,-1.9819,-1.9726,2.1292,1.8364,-1.1395,2.1206,-1.7260,0.3912,0.0079);
	
    }
    void Materials::MixTwoMaterialsShort(const std::string &mtl1,const std::string &mtl2,double coef_mtl1,const std::string &NewMtl)
    {

        if(materialsShort.count(NewMtl)!=1)
        {
            materialsShort[NewMtl]=Material(NewMtl);
        }else{
            return;
        }
        for(auto &mtp1 : materialsShort[mtl1].param)
        {
            materialsShort[NewMtl].param.insert(std::pair<std::string,double>(mtp1.first ,(mtp1.second*coef_mtl1+materialsShort[mtl2].param[mtp1.first]*(1.0-coef_mtl1))));
        }
    }

    void Materials::MixTwoMaterials(const std::string &mtl1,const std::string &mtl2,double coef_mtl1,const std::string &NewMtl)
    {

        if(materials.count(NewMtl)!=1)
        {
            materials[NewMtl]=Material(NewMtl);
        }else{
            return;
        }
        for(auto &mtp1 : materials[mtl1].param)
        {
            materials[NewMtl].param.insert(std::pair<std::string,double>(mtp1.first ,(mtp1.second*coef_mtl1+materials[mtl2].param[mtp1.first]*(1.0-coef_mtl1))));
        }
    }
    void Materials::ReadMaterials(const std::string &fileName,const std::string &mtl)
    {

        if(materials.count(mtl)!=1)
        {
            materials[mtl]=Material(mtl);
        }else{
            return;
        }

        std::ifstream newfile;
        newfile.open(fileName,std::ios::in);
        std::string tp;
        while(getline(newfile,tp)&&newfile.is_open())
        {
            std::stringstream convert;
            convert<<tp;

            std::string NameOfMat;
            convert>>NameOfMat;
            if(mtl==NameOfMat)
            {
                //std::cout<<NameOfMat<<std::endl;
                convert>>materials[mtl].param["00a00a"];
                convert>>materials[mtl].param["00c00c"];
                convert>>materials[mtl].param["01a01a"];
                convert>>materials[mtl].param["01c01c"];
                convert>>materials[mtl].param["02a02a"];
                convert>>materials[mtl].param["02c02c"];
                convert>>materials[mtl].param["10a10a"];
                convert>>materials[mtl].param["10c10c"];
                convert>>materials[mtl].param["00a00c"];
                convert>>materials[mtl].param["00c00a"];
                convert>>materials[mtl].param["10a10c"];
                convert>>materials[mtl].param["10c10a"];
                convert>>materials[mtl].param["10a00c"];
                convert>>materials[mtl].param["00c10a"];
                convert>>materials[mtl].param["00a10c"];
                convert>>materials[mtl].param["10c00a"];
                convert>>materials[mtl].param["00a01c"];
                convert>>materials[mtl].param["00c01a"];
                convert>>materials[mtl].param["10a01c"];
                convert>>materials[mtl].param["10c01a"];
                convert>>materials[mtl].param["00a02c"];
                convert>>materials[mtl].param["00c02a"];
                convert>>materials[mtl].param["10a02c"];
                convert>>materials[mtl].param["10c02a"];
                convert>>materials[mtl].param["01a01c"];
                convert>>materials[mtl].param["01c01a"];

                convert>>materials[mtl].param["01a01cpi"];
                convert>>materials[mtl].param["01c01api"];
                convert>>materials[mtl].param["01a02c"];
                convert>>materials[mtl].param["01c02a"];
                convert>>materials[mtl].param["01a02cpi"];
                convert>>materials[mtl].param["01c02api"];
                convert>>materials[mtl].param["02a02c"];
                convert>>materials[mtl].param["02c02a"];
                convert>>materials[mtl].param["02a02cpi"];
                convert>>materials[mtl].param["02c02api"];
                convert>>materials[mtl].param["02a02cdelta"];
                convert>>materials[mtl].param["02c02adelta"];
                convert>>materials[mtl].param["a"];
                convert>>materials[mtl].param["c"];

                newfile.close();
            }
        }

    }


}
