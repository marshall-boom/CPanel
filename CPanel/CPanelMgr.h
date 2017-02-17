//
//  CPanelMgr.h
//  CPanel
//
//  Created by Chris Satterwhite on 1/16/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__CPanelMgr__
#define __CPanel__CPanelMgr__

#include <stdio.h>
#include "cpCase.h"
#include "cpCaseVP.h"
#include "inputParams.h"
#include "geometry.h"

class caseMgr
{
    inputParams* p;
    geometry* geom;
    Eigen::VectorXi outSpacing;
 
    std::vector<cpCase*> cases;
    std::vector<cpCaseVP*> casesVP;
    
    void setCases();
    void runCases();
    void writeSummary();
    void writeCase(int caseNumber, cpCase* c, std::ofstream &outStream);
public:
    caseMgr(inputParams* p,geometry* geom) : p(p), geom(geom)
    {
        setCases();
        runCases();
        writeSummary();
    }
};


#endif /* defined(__CPanel__CPanelMgr__) */
