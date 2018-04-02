/*******************************************************************************
 * Copyright (c) 2015 Chris Satterwhite
 * Copyright (c) 2016 Connor Sousa
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
 *    Connor Sousa - particle code and implementation
 *    David D. Marshall - misc. improvements
 ******************************************************************************/

#ifndef __CPanel__influenceTests__
#define __CPanel__influenceTests__

#include <stdio.h>
#include <Eigen/Dense>
#include <random>
#include "src/cpptest.h"
#include "bodyPanel.h"
#include "cpNode.h"
#include "edge.h"
#include "geometry.h"
#include <cmath>
#include "DoubVelData.h"
#include "particle.h"
#include "particleOctree.h"
#include "particleFMM.h"
#include "octreeFile.h"


class influenceTests{
    
    double maxCoefficient(Eigen::MatrixXd mat);
    void barnesHutTest(std::vector<particle*> testParts);
    void interactionListTest(std::vector<particle*> testParts);
    std::vector<particle*> randomParticleGenerator(int numParts);
    void BarnesHutSpeedTest(std::vector<particle*> testParts);
    void speedTestGraphComparison();
    void accuracyTest();
    void BarnesHutCarpetPlotData(std::vector<particle*> testParts);


    // Random number generator using code from here: https://isocpp.org/files/papers/n3551.pdf
    double pick_a_number( double from, double upto );
    std::vector<int> linspace(int a, int b, int n);
    std::vector<double> linspace(double a, double b, int n);

    

public:
    void panelVConstTest(Eigen::MatrixXd testingPoints, Eigen::MatrixXd velocityData, bodyPanel* tPan);
    void panelVPntTest(Eigen::MatrixXd testingPoints, Eigen::MatrixXd velocityData, bodyPanel* tPan);
    void velocityComparer(DoubVelData* dat);
    void potentialComparer(bodyPanel* pan);
    
    void FMMtests();
};


//
//class influenceTests : public Test::Suite
//{
//    Eigen::MatrixXd testNodes;
//    bodyPanel* testPan;
//    
//    bool isEqual(double var1, double var2, int decimalPrecision);
//    
//    void createTestPanel();
//    void testPntClose();
//    void testPntFar();
////    void testCollocationPnt();
//public:
//    influenceTests();
//    
//    ~influenceTests() {delete testPan;}
//};
//
#endif /* defined(__CPanel__influenceTests__) */
