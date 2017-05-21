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
#include <vector>
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
#include "particleOctree.h"
#include "vortexFil.h"
#include "particleFMM.h"
#include "octreeFile.h"



class cpCase
{
protected:
    geometry *geom;

    inputParams* params;
    double Vmag;
    double mach;
    double PG; // Prandtl-Glauert Correction - (1-M^2)^(1/2)
    double alpha;
    double beta;
    int timeStep = 0;
    bool vortPartFlag; //VPP
    bool manualStepsSet;
    double numSteps = 1000; // Run for a LOT of steps before convergence criteria kills solution
    double dt;
    std::vector<double> CL; //VPP
    std::vector<double> CM_x; //VPP
    
    bool highAccuracy;
    bool accelerate;
    bool unsteady;

    Eigen::MatrixXd bodyKin;
    
    Eigen::Matrix3d transform;
    
    
    std::vector<bodyPanel*>* bPanels;
    std::vector<wakePanel*>* wPanels;
    std::vector<wakePanel*>* w2panels; //2BW
    std::vector<particle*> particles;
    std::vector<vortexFil*> filaments;

    Eigen::VectorXd sigmas;
    Eigen::VectorXd wake2Doublets;
    Eigen::MatrixXd D;
    
    double CL_trefftz = 0;
    double CD_trefftz = 0;
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
    Eigen::Vector3d Vinf(Eigen::Vector3d POI);
    Eigen::Vector3d windToBody(double V,double alpha,double beta);
    
    Eigen::Vector3d bodyToWind(const Eigen::Vector3d &vec);
    void setSourceStrengths();
    bool solveMatrixEq();
    bool solveVPmatrixEq(); //VPP
    bool solutionConvergence();
    Eigen::Vector3d VinfPlusVecPot(Eigen::Vector3d POI);
    void compVelocity();
    void trefftzPlaneAnalysis();
    void createStreamlines();
    void stabilityDerivatives();
    void writeVTU(std::string filename);
    void writeFiles();
    void writeBodyData(boost::filesystem::path path, const Eigen::MatrixXd &nodeMat);
    void writeWakeData(boost::filesystem::path path, const Eigen::MatrixXd &nodeMat);
    void writeBuffWake2Data(boost::filesystem::path path, const Eigen::MatrixXd &nodeMat); //2BW
    void writeParticleData(boost::filesystem::path path);
    void writeFilamentData(boost::filesystem::path path);
    void writeSpanwiseData(boost::filesystem::path path);
    void writeBodyStreamlines(boost::filesystem::path path);
//    void collapsePanels(); //BW2
    void collapseBufferWake();
    void collapseWakeForEachEdge();
    bool edgeIsUsed(edge* thisEdge, std::vector<edge*> pEdges);
    Eigen::Vector3d edgeStrength(wakePanel* pan, edge* curEdge, int edgeNum);
    Eigen::Vector3d seedPos(wakePanel* pan, int edgeNum);


    void convectParticles();
    void vortexStretching();
    void particleStrengthUpdate();
    void particleStrengthUpdateGaussian();
    
    Eigen::Vector3d velocityInflFromEverything(Eigen::Vector3d POI);
    Eigen::Vector3d velocityInflFromEverything(particle* part);
    
    particleOctree partOctree;
    particleFMM FMM;

    void convectBufferWake(); //VPP
    void readBodyKinFile();
    
    
    Eigen::Vector3d rungeKuttaStepper(Eigen::Vector3d POI);
    
    
    struct meshDat {
        std::vector<Eigen::Vector3d> velocity;
        std::vector<double> vorticity;
        std::vector<double> pressure;
        std::vector<Eigen::Vector3d> cellCenter;
    } volMeshDat;
    
    std::vector<Eigen::VectorXi> cells;
    Eigen::MatrixXd pts;

    void createVolMesh();
    void writeVolMeshData(boost::filesystem::path path, Eigen::MatrixXd &nodeMat, std::vector<Eigen::VectorXi> cells);

    double changeCDnm1, changeCLnm1;
    

public:
    cpCase(geometry *geom, double V, double alpha, double beta, double mach, inputParams* inParams) : geom(geom), Vmag(V), alpha(alpha), beta(beta), mach(mach), params(inParams)
    {
//        Vinf = windToBody(V,alpha,beta); 
        bPanels = geom->getBodyPanels();
        wPanels = geom->getWakePanels();
        w2panels = geom->getBufferWake2Panels(); //2BW
        PG = sqrt(1-pow(mach,2));
        
        vortPartFlag = inParams->vortPartFlag; //VPP
        dt = geom->getDt();
        if(inParams->numSteps != 0){
            manualStepsSet = true;
            numSteps = inParams->numSteps;
        }
        accelerate = inParams->accel;
        highAccuracy  = inParams->high_accuracy;
        
        unsteady = inParams->unsteady;
        
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
    
    struct partWakeMaxDims{};
    partWakeMaxDims partMaxDims();
    
};
#endif /* defined(__CPanel__runCase__) */
