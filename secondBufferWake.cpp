//
//  secondBufferWake.cpp
//  CPanel - Unstructured Panel Code
//
//  Created by Connor Sousa on 2/29/16.
//  Copyright (c) 2016 Chris Satterwhite. All rights reserved.
//

#include "secondBufferWake.h"

secondBufferWake::secondBufferWake(std::vector<cpNode*> nodes, std::vector<edge*> pEdges, Eigen::Vector3d bezNorm, wake* parentWake, wakePanel* parentPan, int surfID) : wakePanel(nodes,pEdges,bezNorm,parentWake,surfID), parentPanel(parentPan)
{
    for (int i=0; i<pEdges.size(); i++)
    {
        pEdges[i]->addSBWPan(this);
    }
}

