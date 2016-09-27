//
//  influenceTests.h
//  CPanel
//
//  Created by Chris Satterwhite on 1/24/15.
//  Recreated by Connor Sousa on 8/22/16.
//
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

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
    // Random number generator using code from here: https://isocpp.org/files/papers/n3551.pdf
    double pick_a_number( double from, double upto );
    std::vector<int> linspace(int a, int b, int n);

    

public:
    void panelVConstTest(Eigen::MatrixXd testingPoints, Eigen::MatrixXd velocityData, bodyPanel* tPan);
    void panelVPntTest(Eigen::MatrixXd testingPoints, Eigen::MatrixXd velocityData, bodyPanel* tPan);
    void velocityComparer(DoubVelData* dat);
    
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
