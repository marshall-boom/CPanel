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
 *    David D. Marshall - misc. changes
 ******************************************************************************/

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
