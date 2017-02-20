//
//  cpCaseVP.h
//  CPanel
//
//  Created by Connor Sousa on 1/31/17.
//  Copyright (c) 2017 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__cpCaseVP__
#define __CPanel__cpCaseVP__

#include <iostream>
#include "cpCase.h"
#include "particle.h"
#include "vortexFil.h"
#include "particleFMM.h"

#include "octreeFile.h"



class cpCaseVP : public cpCase{
    
    std::vector<particle*> particles;
    std::vector<vortexFil*> filaments;
    std::vector<wakePanel*>* w2panels; // Second row of buffer wake
    
    
    bool solveMatrixEq();
    bool solveVPmatrixEq();
    void convectBufferWake();
    bool edgeIsUsed(edge* thisEdge, std::vector<edge*> pEdges);
    void collapseBufferWake();
    void setSourceStrengthsVP();
    void compVelocity();
    void trefftzPlaneAnalysis();
    void particleStrengthUpdate();
    void convectParticles();
    void writeFiles();
    
    bool solutionConvergence();
    double changeCDnm1, changeCLnm1;
    
    Eigen::Vector3d Vinf(Eigen::Vector3d POI);
    Eigen::Vector3d VinfPlusVecPot(Eigen::Vector3d POI);
    Eigen::Vector3d edgeStrength(wakePanel* pan, edge* curEdge, int edgeNum);
    Eigen::Vector3d rungeKuttaStepper(Eigen::Vector3d POI);
    
    Eigen::Vector3d velocityInflFromEverything(Eigen::Vector3d POI);
    Eigen::Vector3d velocityInflFromEverything(particle* part);
    
    void writeFilesVP();
    void writeBodyDataVP(boost::filesystem::path path,const Eigen::MatrixXd &nodeMat);
    void writeWakeDataVP(boost::filesystem::path path, const Eigen::MatrixXd &nodeMat);
    void writeBuffWake2Data(boost::filesystem::path path, const Eigen::MatrixXd &nodeMat);
    void writeFilamentData(boost::filesystem::path path);
    void writeParticleData(boost::filesystem::path path);
    
    void readBodyKinFile();
    
    Eigen::VectorXd wake2Doublets;
    int timestep = 1;
    
    double dt;
    int numSteps;
    bool accelerate;
    bool unsteady;
    
    
    particleOctree partOctree;
    particleFMM FMM;
    Eigen::MatrixXd bodyKin;

    
    
    
public:
    cpCaseVP( geometry *geom, double V, double alpha, double beta, double mach, inputParams* inParams ) : cpCase( geom, V, alpha, beta, mach, inParams )
    {
        bPanels = geom->getBodyPanels();
        wPanels = geom->getWakePanels();
        w2panels = geom->getBufferWake2Panels();
        dt = geom->getDt();
        PG = sqrt(1-pow(mach,2));
        
        numSteps = inParams->numSteps;
        accelerate = inParams->accelerateCode;
        unsteady = inParams->unsteady;
    }
    
    
    void run(bool printFlag, bool surfStreamFlag, bool stabDerivFlag);
  
};





#endif /* defined(__CPanel__cpCaseVP__) */
