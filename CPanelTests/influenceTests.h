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
#include "src/cpptest.h"
#include "bodyPanel.h"
#include "cpNode.h"
#include "edge.h"
#include "geometry.h"
#include <cmath>
#include "DoubVelData.h"


class influenceTests{
    
    double maxCoefficient(Eigen::MatrixXd mat);
    
public:
    void panelVConstTest(Eigen::MatrixXd testingPoints, Eigen::MatrixXd velocityData, bodyPanel* tPan);
    void panelVPntTest(Eigen::MatrixXd testingPoints, Eigen::MatrixXd velocityData, bodyPanel* tPan);
    void velocityComparer(DoubVelData* dat);

    
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
