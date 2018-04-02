/*******************************************************************************
 * Copyright (c) 2015 Chris Satterwhite
 * Copyright (c) 2018 David D. Marshall <ddmarsha@calpoly.edu>
 *
 * This program and the accompanying materials are made
 * available under the terms of the Eclipse Public License 2.0
 * which is available at https://www.eclipse.org/legal/epl-2.0/
 *
 * See LICENSE.md file in the project root for full license information.
 *
 * SPDX-License-Identifier: EPL-2.0
 *
 * Contributors:
 *    Chris Satterwhite - initial code and implementation
 *    Connor Sousa - Vortex particle & unsteady implementation
 *    David D. Marshall - misc. changes
 ******************************************************************************/

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
    using volMeshRes_type = std::vector<size_t>;
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
