//
//  inputParams.h
//  CPanel
//
//  Created by Chris Satterwhite on 1/16/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__inputParams__
#define __CPanel__inputParams__

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <Eigen/Dense>
#include "cpFile.h"

struct inputParams
{
    cpFile* inputFile;
    cpFile* geomFile;
    bool normFlag; //Normals included in input file
    double Sref;
    double bref;
    double cref;
    double timeStep;
    double numSteps;
    bool stepsSetMaunally;
    double streamSpacing;
    bool accelerateCode = false;
    bool unsteady = false;
    std::string bodyKinFileLoc;

    Eigen::Vector3d cg;
    Eigen::VectorXd velocities,alphas,betas,machs;

    bool surfStreamFlag = false;
    bool stabDerivFlag = false;
    bool writeCoeffFlag = false;
    bool vortexParticles = false;
    bool volMeshFlag = false;

    std::vector<double> volMeshBounds;
    std::vector<int> volMeshRes;


    inputParams(cpFile* inFile) : inputFile(inFile) {}

    ~inputParams()
    {
        delete geomFile;
    }

    bool set();
    void print(std::ostream &stream);

private:
    bool checkGeomFile();
    void makeWorkingDir();
    void writeInputFile();
    void printVec(Eigen::VectorXd &vec,std::ostream &stream);
};


#endif /* defined(__CPanel__inputParams__) */
