#ifndef Eigen_H
#define Eigen_H
#include "eigen-3.4.0/Eigen/Eigen"
#endif

#ifndef TB_MODULE_H
#define TB_MODULE_H

#define LoMEType1 1
#define LoMEType2 2
#define LoMEType3 4
#define LoMC 8
#define RED     "\x1b[31m"
#define GREEN   "\x1b[32m"
#define YELLOW  "\x1b[33m"
#define BLUE    "\x1b[34m"
#define MAGENTA "\x1b[35m"
#define CYAN    "\x1b[36m"
#define RESET   "\x1b[0m"


#define iA (1e10)/* m */
#define e_c (1.602176487e-19)/* A s */
#define h_bar (1.05457162825e-34) /* kg m^2 / s */
//#define h_bar (1.05457e-34) /* kg m^2 / s */
#define hhh (6.62606896e-34) /* kg m^2 / s */
#define m_e (9.10938188e-31) /* kg */

#include "mkl_lapacke.h"
#include "mkl_lapack.h"
#include "mkl_scalapack.h"
#include <mkl_pblas.h>
#include <mkl_blacs.h>
#include <map>
#include <array>
#include <ctime>
#include <ratio>
#include <random>
#include <mpi.h>
#include <string>
#include <memory>
#include <unistd.h>
#include <vector>
#include <chrono>
#include <fstream>
#include <complex>
#include <numeric>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <sys/time.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <unordered_map>
#include <unordered_set>
#include <iterator>
#include <regex>
#include <boost/program_options.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/range/adaptors.hpp>
#include <boost/range/algorithm.hpp>

extern "C" {
    int numroc_(const int*, const int*, const int*, const int*, const int*);
    void descinit_(int *, const int *, const int *, const int *, const int *, const int *, const int *,const int *, const int *,int *);
}

//std::vector<std::string> split(std::string stringToBeSplitted, std::string delimeter);
Eigen::Vector3d RotateAround(const Eigen::Vector3d &Axis,const Eigen::Vector3d &ts,double angle);

template<class _D>
std::vector<std::string> split(std::string stringToBeSplitted, _D TypeDdelimeter)
{
    std::stringstream convert;
    convert<<TypeDdelimeter;std::string delimeter;
    convert>>delimeter;

    std::vector<std::string> splittedString;
    int startIndex = 0;
    int  endIndex = 0;
    while( (endIndex = stringToBeSplitted.find(delimeter, startIndex)) < stringToBeSplitted.size() )
    {

        std::string val = stringToBeSplitted.substr(startIndex, endIndex - startIndex);
        splittedString.push_back(val);
        startIndex = endIndex + delimeter.size();

    }
    if(startIndex < stringToBeSplitted.size())
    {
        std::string val = stringToBeSplitted.substr(startIndex);
        splittedString.push_back(val);
    }
    return splittedString;
}


template<class _R,class _L>
_R ParseValue(std::string key,std::string inpfile,std::string outtag,_R bakup, _L signal){
    std::ifstream newfile;
    newfile.open(inpfile,std::ios::in);
    std::string tp;
    while(newfile.is_open()&&getline(newfile,tp)){
        boost::trim(tp);
        if(tp.find("#")!=std::string::npos) tp.erase(tp.find("#"),std::string::npos);
        if(tp.find("=")!=std::string::npos){
            std::vector<std::string> vLine;
            vLine=split<std::string>(tp,"=");
            for(auto&vline:vLine) boost::trim(vline);
            if(key==vLine[0] && vLine.size()>1){
                std::stringstream convert;
                convert<<vLine[1];
                _R ret;
                convert>>ret;

                if(signal==0){
                    std::fstream file_obj; 
                    file_obj.open("./out"+outtag+".txt",std::ios::out|std::ios::app);
                    file_obj<<vLine[0]<<" = "<<ret<<std::endl;
                    file_obj.close();
                }

                newfile.close();
                return ret;
            }
        }
    }
    newfile.close();
    return bakup;
}

template<class _R>
_R ParseBlockTitle(std::string key,std::string inpfile,_R NotFound){
    std::ifstream newfile;
    newfile.open(inpfile,std::ios::in);
    std::string tp;
    while(newfile.is_open()&&getline(newfile,tp)){
        boost::trim(tp);
        if(tp.find("#")!=std::string::npos) tp.erase(tp.find("#"),std::string::npos);
        if(tp.find("%")!=std::string::npos && tp!="%"){
            std::vector<std::string> vLine;
            vLine=split<std::string>(tp,"%");
            for(auto&vline:vLine) boost::trim(vline);
            //for(auto&vline:vLine) std::cout<<"vLine[]="<<vline<<std::endl;
            if(key==vLine[1] && vLine.size()>2){
                std::stringstream convert;
                convert<<vLine[2];
                _R ret;
                convert>>ret;
                newfile.close();
                return ret;
            }

        }
    }
    newfile.close();
    return NotFound;
}

template<class _L>
void ParseBlockElement(std::string key,std::string inpfile,int r,int c, _L &Input){
    std::ifstream newfile;
    newfile.open(inpfile,std::ios::in);
    std::string tp;
    while(newfile.is_open()&&getline(newfile,tp)){
        boost::trim(tp);
        if(tp.find("#")!=std::string::npos) tp.erase(tp.find("#"),std::string::npos);
        if(tp.find("%")!=std::string::npos && tp!="%"){
            std::vector<std::string> vLine;
            vLine=split<std::string>(tp,"%");
            for(auto&vline:vLine) boost::trim(vline);
            //for(auto&vline:vLine) std::cout<<"vLine[]="<<vline<<std::endl;
            if(key==vLine[1]){
                for(int l=0;l<=r;++l) getline(newfile,tp);
                if(tp.find("#")!=std::string::npos) tp.erase(tp.find("#"),std::string::npos);
                std::vector<std::string> sLine;
                sLine=split<std::string>(tp,"|");
                for(auto&sline:sLine) boost::trim(sline);
                std::stringstream convert;
                convert<<sLine[c];_L ret;
                convert>>ret;
                Input=ret;
                newfile.close();
            }

        }
    }
    newfile.close();
}

#endif
