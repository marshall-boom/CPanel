//
//  wakeLine.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 10/26/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "wakeLine.h"

wakeLine::wakeLine(bodyPanel* uupper, bodyPanel* llower, Eigen::Vector3d nnormal)
  : upper(uupper), lower(llower), normal(nnormal)
{
    setDimensions();
}

double wakeLine::getStrength()
{
    return upper->getMu()-lower->getMu();
}

void wakeLine::setDimensions()
{
    std::vector<edge*> es = upper->getEdges();
    std::vector<cpNode*> TEnodes;
    for (size_t i=0; i<es.size(); i++)
    {
        if (es[i]->isTE())
        {
            TEnodes = es[i]->getNodes();
        }
    }

    p1 = TEnodes[0]->getPnt();
    p2 = TEnodes[1]->getPnt();
    pMid = (p1+p2)/2;
}
