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

#include "convexHull.h"

convexHull::convexHull(Eigen::MatrixXd points, bool bound) : boundary(bound)
{
    Eigen::Vector3d basePnt = points.row(0);
    double minTheta = M_PI;
    double maxTheta = 0;
    
    for (int i=1; i<points.rows(); i++)
    {
        if ((points(i,1) == basePnt(1) && points(i,0) < basePnt(0)) || (points(i,1) < basePnt(1)))
        {
            basePnt = points.row(i);
        }
    }
    for (int i=0; i<points.rows(); i++)
    {
        member* m = new member(points.row(i),basePnt);
        if (points(i,0) != basePnt(0) || points(i,1) != basePnt(1) || points(i,2) != basePnt(2))
        {
            if (m->theta < minTheta)
            {
                minTheta = m->theta;
            }
            if (m->theta > maxTheta)
            {
                maxTheta = m->theta;
            }
        }
        if (points(i,0) == basePnt(0) && points(i,1) == basePnt(1) && points(i,2) == basePnt(2))
        {
            members.insert(members.begin(),m);
        }
        else
        {
            members.push_back(m);
        }
    }
    std::sort(members.begin()+1,members.end(),compareTheta());
    size_t minBegin = 0;
    bool minBeginSet = false;
    size_t minEnd = 0;
    bool minEndSet = false;
    members_index_type maxBegin = members.size();
    for (members_index_type i=1; i<members.size()-1; i++)
    {
        if (members[i]->theta == minTheta && !minBeginSet)
        {
            minBegin = i;
        }
        else if (members[i]->theta != minTheta && minBeginSet && minEndSet)
        {
            minEnd = i;
        }
        else if (members[i]->theta == maxTheta && maxBegin == members.size())
        {
            maxBegin = i;
        }
    }
    std::sort(members.begin()+static_cast<members_type::difference_type>(minBegin),members.begin()+static_cast<members_type::difference_type>(minEnd),compareDistAscending());
    std::sort(members.begin()+static_cast<members_type::difference_type>(maxBegin),members.end(),compareDistDescending());
    computeHull();
}


convexHull::member::member(const Eigen::Vector3d &point, const Eigen::Vector3d &basePnt)
{
    x = point(0);
    y = point(1);
    z = point(2);
    theta = atan2((y-basePnt(1)),(x-basePnt(0)));
    d = sqrt(pow(y-basePnt(1),2)+pow(x-basePnt(0),2));
}

void convexHull::computeHull()
{
    hull.push_back(members.front());
    for (members_index_type i=0; i<members.size()-1; i++)
    {
        Eigen::Vector3d v1,v2;
        if (i==members.size()-2)
        {
            v1 = makeVector(members[i],members[i+1]);
            v2 = makeVector(members[i],members[0]);
        }
        else
        {
            v1 = makeVector(members[i],members[i+1]);
            v2 = makeVector(members[i],members[i+2]);
        }
        
        double crossZ = v1.cross(v2)(2);
        
        if ((crossZ > 0) || (crossZ == 0 && !boundary))
        {
            hull.push_back(members[i+1]);
        }
    }
}

Eigen::Vector3d convexHull::makeVector(member* p1, member* p2)
{
    Eigen::Vector3d vec;
    vec(0) = p2->x-p1->x;
    vec(1) = p2->y-p1->y;
    vec(2) = 0; //Only computing 2D complex hull
    return vec;
}

bool convexHull::compareNodes(std::vector<Eigen::Vector3d> nodesLocal)
{
    if (hull.size() > nodesLocal.size())
    {
        return false;
    }
    
    Eigen::Vector3d pMember;
    bool breakFlag = false;
    
    for (members_index_type i=0; i<hull.size(); i++)
    {
        pMember(0) = hull[i]->x;
        pMember(1) = hull[i]->y;
        pMember(2) = hull[i]->z;
        for (members_index_type j=0; j<nodesLocal.size(); j++)
        {
            if (pMember == nodesLocal[j])
            {
                breakFlag = true;
                break;
            }
        }
        if (breakFlag)
        {
            breakFlag = false;
            continue;
        }
        else
        {
            return false;
        }
    }
    return true;
}
