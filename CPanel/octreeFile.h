//
//  octreeFile.h
//  CPanel
//
//  Created by Chris Satterwhite on 11/15/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__octreeFile__
#define __CPanel__octreeFile__

#include <stdio.h>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include "panelOctree.h"
#include "particleOctree.h"

class octreeFile
{
    void writeFile(std::string filename,panelOctree* oct);
    void writeFile(std::string filename,particleOctree* oct);
    
public:
    octreeFile(std::string filename,panelOctree* oct)
    {
        writeFile(filename,oct);
    }
    
    octreeFile(std::string filename,particleOctree* oct)
    {
        writeFile(filename,oct);
    }
};
#endif /* defined(__CPanel__octreeFile__) */
