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

#ifndef CPanel_member_h
#define CPanel_member_h

#include "Eigen/Dense"

template<typename type>
class member
{
    type* pObj;
    Eigen::Vector3d refPoint;

public:
    member(type* object, Eigen::Vector3d referencePoint) : pObj(object), refPoint(referencePoint) {}

    Eigen::Vector3d getRefPoint() const {return refPoint;}

    type *getObj() {return pObj;}
};


#endif
