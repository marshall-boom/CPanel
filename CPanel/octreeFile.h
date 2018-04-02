/*******************************************************************************
 * Copyright (c) 2014 Chris Satterwhite
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
 *    Connor Sousa - Vortex particle implementation
 *    David D. Marshall - misc. changes
 ******************************************************************************/

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
