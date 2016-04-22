//
//  runCase.h
//  CPanel
//
//  Created by Chris Satterwhite on 10/13/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__runCase__
#define __CPanel__runCase__

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <Eigen/Dense>
#include <Eigen/IterativeLinearSolvers>
#include <boost/filesystem/operations.hpp>
#include "geometry.h"
#include "VTUfile.h"
#include "bodyStreamline.h"
#include "cpFile.h"
#include "inputParams.h"
//#include "streamline.h"
#include "particle.h"
#include "edge.h" //VPP



class cpCase
{
    geometry *geom;

    inputParams* params;
    double Vmag;
    double mach;
    double PG; // Prandtl-Glauert Correction - (1-M^2)^(1/2)
    double alpha;
    double beta;
    int timeStep = 0; //VPP
    bool vortPartFlag; //VPP
    double dt;
    
    
    Eigen::Vector3d Vinf;
    Eigen::Matrix3d transform;
    
    
    std::vector<bodyPanel*>* bPanels;
    std::vector<wakePanel*>* wPanels;
//    std::vector<wakePanel*>* wake2panels; //2BW
    std::vector<particle*> particles;

    Eigen::VectorXd sigmas;
    Eigen::VectorXd wake2Doublets;
    Eigen::VectorXd particleStrengths;
    Eigen::MatrixXd D;
    
    double CL_trefftz;
    double CD_trefftz;
    Eigen::Vector3d Fbody;
    Eigen::Vector3d Fwind;
    Eigen::Vector3d CM; //[roll,pitch,yaw]
    Eigen::VectorXd spanLoc;
    Eigen::VectorXd Cl;
    Eigen::VectorXd Cd;
    
    Eigen::Vector3d dF_dAlpha;
    Eigen::Vector3d dF_dBeta;
    Eigen::Vector3d dM_dAlpha;
    Eigen::Vector3d dM_dBeta;
    
    std::vector<bodyStreamline*> bStreamlines;
    Eigen::Vector3d windToBody(double V,double alpha,double beta);
    
    Eigen::Vector3d bodyToWind(const Eigen::Vector3d &vec);
    void setSourceStrengths();
    void setVPSourceStrengths();
    bool solveMatrixEq();
    bool solveVPmatrixEq(); //VPP
    void compVelocity();
    void trefftzPlaneAnalysis();
    void createStreamlines();
    void stabilityDerivatives();
    void writeVTU(std::string filename);
    void writeFiles();
    void writeBodyData(boost::filesystem::path path, const Eigen::MatrixXd &nodeMat);
    void writeWakeData(boost::filesystem::path path, const Eigen::MatrixXd &nodeMat);
//    void writeBuffWake2Data(boost::filesystem::path path, const Eigen::MatrixXd &nodeMat); //2BW
    void writeParticleData(boost::filesystem::path path);
    void writeSpanwiseData(boost::filesystem::path path);
    void writeBodyStreamlines(boost::filesystem::path path);
//    void collapsePanels(); //BW2
    void collapseBufferWake();
    void convectParticles();
    void vortexStretching();


//    void convectBufferWake(); //VPP

    
public:
    cpCase(geometry *geom, double V, double alpha, double beta, double mach, inputParams* inParams) : geom(geom), Vmag(V), alpha(alpha), beta(beta), mach(mach), params(inParams)
    {
        Vinf = windToBody(V,alpha,beta);
        bPanels = geom->getBodyPanels();
        wPanels = geom->getWakePanels();
//        wake2panels = geom->getWake2Panels(); //2BW
        PG = sqrt(1-pow(mach,2));
        vortPartFlag = inParams->vortPartFlag; //VPP
//        dt = inParams->timeK*inParams->cref/Vinf.norm(); //VPP
        dt = geom->dt;
    }
    
    virtual ~cpCase();
    
    void run(bool printFlag, bool surfStreamFlag, bool stabDerivFlag, bool vortPartFlag);
    
    double getMach() {return mach;}
    double getV() {return Vmag;}
    double getAlpha() {return alpha;}
    double getBeta() {return beta;}
    double getTimeStep() {return timeStep;} //VPP
    double getCL() {return CL_trefftz;}
    double getCD() {return CD_trefftz;}
    Eigen::Vector3d getMoment() {return CM;}
    Eigen::Vector3d getBodyForces() {return Fbody;}
    Eigen::Vector3d getWindForces() {return Fwind;}
    Eigen::Vector3d get_dF_dAlpha() {return dF_dAlpha;}
    Eigen::Vector3d get_dF_dBeta() {return dF_dBeta;}
    Eigen::Vector3d get_dM_dAlpha() {return dM_dAlpha;}
    Eigen::Vector3d get_dM_dBeta() {return dM_dBeta;}
    
};
#endif /* defined(__CPanel__runCase__) */
