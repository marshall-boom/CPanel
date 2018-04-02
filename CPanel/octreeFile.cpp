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

#include "octreeFile.h"

void octreeFile::writeFile(std::string filename,panelOctree* oct)
{
    std::vector<Eigen::Vector3d> centers;
    std::vector<Eigen::Vector3d> halfDs;
    std::vector<node<panel>*> nodes;
    nodes = oct->getNodes();
    
    for (size_t i=0; i<nodes.size(); i++)
    {
        centers.push_back(nodes[i]->getOrigin());
        halfDs.push_back(nodes[i]->getHalfDimension());
    }
    assert(centers.size() == halfDs.size());
    std::ofstream fid;
    fid.open(filename);
    if (fid.is_open())
    {
        fid << (nodes.size()+1) << "\n";
        for (size_t i=0; i<centers.size(); i++)
        {
            fid << centers[i](0) << "\t" << centers[i](1) << "\t" << centers[i](2) << "\n";
        }
        for (size_t i=0; i<halfDs.size(); i++)
        {
            fid << halfDs[i](0) << "\t" << halfDs[i](1) << "\t" << halfDs[i](2) << "\n";
        }
        fid.close();
    }
    else
    {
        std::cout << "Octree file could not be opened" << std::endl;
    }
}

void octreeFile::writeFile(std::string filename,particleOctree* oct)
{
    std::vector<Eigen::Vector3d> centers;
    std::vector<Eigen::Vector3d> halfDs;
    std::vector<node<particle>*> nodes;
    nodes = oct->getNodes();
    
    for (size_t i=0; i<nodes.size(); i++)
    {
        centers.push_back(nodes[i]->getOrigin());
        halfDs.push_back(nodes[i]->getHalfDimension());
    }
    assert(centers.size() == halfDs.size());
    std::ofstream fid;
    fid.open(filename);
    if (fid.is_open())
    {
        fid << (nodes.size()+1) << "\n";
        for (size_t i=0; i<centers.size(); i++)
        {
            fid << centers[i](0) << "\t" << centers[i](1) << "\t" << centers[i](2) << "\n";
        }
        for (size_t i=0; i<halfDs.size(); i++)
        {
            fid << halfDs[i](0) << "\t" << halfDs[i](1) << "\t" << halfDs[i](2) << "\n";
        }
        fid.close();
    }
    else
    {
        std::cout << "Octree file could not be opened" << std::endl;
    }
}
