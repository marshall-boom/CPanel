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

#ifndef __CPanel___Unstructured_Panel_Code__particleOctree__
#define __CPanel___Unstructured_Panel_Code__particleOctree__

#include <iostream>
#include "octree.h"
#include "particle.h"

class particle;

class particleOctree : public octree<particle>
{

public:
    particleOctree() : octree() {}
    
    Eigen::Vector3d findRefPoint(const particle &obj);
};



#endif /* defined(__CPanel___Unstructured_Panel_Code__particleOctree__) */
