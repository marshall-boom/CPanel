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
 *    David D. Marshall - misc. improvements
 ******************************************************************************/

#include <iostream>
#include "octree.h"
#include <vector>

class thing
{
public:
    std::array<double,3> center;
    
    thing()
    {
        center[0] = 0;
        center[1] = 0;
        center[2] = 0;
    }
};

class testOctree : public octree<thing>
{
public:
    std::array<double,3> findRefPoint(const thing* obj)
    {
        return (obj->center);
    }
    
    testOctree() : octree() {}
};

int main()
{
    std::vector<thing*> data;
    thing* obj;
    int nX = 10;
    int nY = 10;
    int nZ = 10;
    for (int i=0; i<nX; i++)
    {
        for (int j=0; j<nY; j++)
        {
            for (int k=0; k<nZ; k++)
            {
                obj = new thing;
                obj->center[0] = i;
                obj->center[1] = j;
                obj->center[2] = k;
                data.push_back(obj);
            }
        }
    }
    testOctree myOctree;
    myOctree.addData(data);
    
}

