//
//  particleOctree.cpp
//  CPanel - Unstructured Panel Code
//
//  Created by Connor Sousa on 5/9/16.
//  Copyright (c) 2016 Connor Sousa. All rights reserved.
//

#include "particleOctree.h"

Eigen::Vector3d particleOctree::findRefPoint(const particle &obj){
    return obj.pos;
}


