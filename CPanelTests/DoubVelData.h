/*******************************************************************************
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
 *    Connor Sousa - initial code and implementation
 *    David D. Marshall - porting to GoogleTest
 ******************************************************************************/

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
