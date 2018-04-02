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
 *    David D. Marshall - misc. changes
 ******************************************************************************/

#ifndef __CPanel___Unstructured_Panel_Code__particleFMM__
#define __CPanel___Unstructured_Panel_Code__particleFMM__

#include <iostream>
#include "particleOctree.h"

class particleFMM
{
    // Might structure something like: build the multipole expansion, then have separate functions that either: calc velocity for a specific point or step through all of the particles and step


    particleOctree* partTree;
    size_t levels;

    void computeMultExp();

    Eigen::Vector3d findExpPos(node<particle>* thisNode);
    Eigen::Vector3d findExpPos(std::vector<particle*> parts);

    Eigen::Vector3d findExpStrength(node<particle>* thisNode);
    Eigen::Vector3d findExpStrength(std::vector<particle*> parts);

    double findExpRadius(node<particle>* thisNode);
    double findExpRadius(std::vector<particle*> parts);

    void printParts();
    void printMultExps();

    std::vector<node<particle>*> getChildren(std::vector<node<particle>*> nodeVec);


public:
    particleFMM() : partTree(nullptr), levels(0) {};

    void build(particleOctree* tree);
    Eigen::Vector3d barnesHutVel(Eigen::Vector3d POI);
    Eigen::Vector3d barnesHutVel(particle* part);
    Eigen::Vector3d barnesHutStretch(particle* part);
    Eigen::Vector3d barnesHutDiff(particle* part);
};


#endif /* defined(__CPanel___Unstructured_Panel_Code__particleFMM__) */
