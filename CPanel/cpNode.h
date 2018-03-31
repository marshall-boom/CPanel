//
//  cpNode.h
//  CPanel - Unstructured Panel Code
//
//  Created by Chris Satterwhite on 2/2/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel___Unstructured_Panel_Code__cpNode__
#define __CPanel___Unstructured_Panel_Code__cpNode__

#include <stdio.h>
#include <Eigen/Dense>
#include <vector>

class edge;
class bodyPanel;


class cpNode
{
	using edges_type = std::vector<edge *>;
	using edges_index_type = edges_type::size_type;
	using bodyPanels_type = std::vector<bodyPanel *>;
	using bodyPanels_index_type = bodyPanels_type::size_type;

    Eigen::Vector3d pnt;
    size_t index;
    edges_type edges;
    bodyPanels_type bodyPans;
    double c_w = 1.3; //vpp
    
    bool TEnode;
    
public:
    cpNode(Eigen::Vector3d pnt,size_t iindex);
        
    Eigen::Vector3d operator-=(const cpNode &rhs);
    
    Eigen::Vector3d operator+=(const cpNode &rhs);
    
    void addEdge(edge* e);
    void addBodyPanel(bodyPanel* p);
    
    edge* getTE(edge* exception);
    
    void setTE();
    void setIndex(size_t i);
    
    void setPnt(Eigen::Vector3d pos){pnt = pos;}
    Eigen::Vector3d getPnt() const {return pnt;}
    
    size_t getIndex() const {return index;}
    
    std::vector<edge*> getEdges() {return edges;}
    
    std::vector<bodyPanel*> getBodyPans() {return bodyPans;}
    
    edge* getOtherTrailEdge(edge* current);
    std::vector<edge*> getTrailingEdges();
    
    Eigen::Vector3d nodeWakeProjAngle();
    Eigen::Vector3d firstProjNode(double dt, double inputV);
    Eigen::Vector3d secProjNode(double dt, double inputV);
    
    bool isTE() {return TEnode;}

    
};

#endif /* defined(__CPanel___Unstructured_Panel_Code__cpNode__) */
