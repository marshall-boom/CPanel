//
//  DoubVelData.h
//  CPanel - Unstructured Panel Code
//
//  Created by Connor Sousa on 8/30/16.
//  Copyright (c) 2016 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel___Unstructured_Panel_Code__DoubVelData__
#define __CPanel___Unstructured_Panel_Code__DoubVelData__

#include <iostream>
#include <Eigen/Dense>
#include "influenceTests.h"
#include "geometry.h"
#include "cpFile.h"
#include "inputParams.h"
#include "DoubVelData.h"

struct DoubVelData {
    Eigen::Vector3d colocationPnt, n1, n2, n3, n4;
    double panelSideLength, mu, sigma;
    int numVertPts, numDiagPts, numMedPts;
    Eigen::MatrixXd vertSurveyPts; Eigen::MatrixXd diagSurveyPts; Eigen::MatrixXd MedSurveyPts;
    Eigen::MatrixXd vertVelPnt; Eigen::MatrixXd vertVelConst; Eigen::MatrixXd diagVelPnt;
    Eigen::MatrixXd diagVelConst; Eigen::MatrixXd medVelPnt; Eigen::MatrixXd medVelConst;
    
    std::vector<cpNode*> nodes;
    std::vector<edge*> pEdges;
    bodyPanel* testPan;
    
};



#endif /* defined(__CPanel___Unstructured_Panel_Code__DoubVelData__) */
