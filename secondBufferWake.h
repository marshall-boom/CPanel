//
//  secondBufferWake.h
//  CPanel - Unstructured Panel Code
//
//  Created by Connor Sousa on 2/29/16.
//  Copyright (c) 2016 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel___Unstructured_Panel_Code__secondBufferWake__
#define __CPanel___Unstructured_Panel_Code__secondBufferWake__

#include <iostream>
#include "wakePanel.h"


class secondBufferWake : public wakePanel{
    
    wakePanel* parentPanel;
    
public:
    secondBufferWake(std::vector<cpNode*> nodes, std::vector<edge*> pEdges, Eigen::Vector3d bezNorm, wake* parentWake, wakePanel* parentPanel, int surfID);
    
    
};



#endif /* defined(__CPanel___Unstructured_Panel_Code__secondBufferWake__) */
