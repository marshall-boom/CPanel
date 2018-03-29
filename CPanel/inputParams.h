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
    int numSteps;
    bool stepsSetManually;
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

    using volMeshBounds_type = std::vector<double>;
    using volMeshBounds_index_type = volMeshBounds_type::size_type;
    using volMeshRes_type = std::vector<int>;
    using volMeshRes_index_type = volMeshRes_type::size_type;

    volMeshBounds_type volMeshBounds;
    volMeshRes_type volMeshRes;


    inputParams(cpFile* inFile)
      : inputFile(inFile), geomFile(nullptr), normFlag(false), Sref(1), bref(1), cref(1),
		timeStep(0), numSteps(0), stepsSetManually(false), streamSpacing(0) {}

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
