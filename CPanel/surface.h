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

#ifndef __CPanel__surface__
#define __CPanel__surface__

#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include "panel.h"
#include "wakePanel.h"
#include "bodyPanel.h"

class geometry;
class bodyStreamline;

class surface
{
protected:
	using panels_type = std::vector<bodyPanel *>;
	using panels_index_type = panels_type::size_type;

    geometry* geom;
    panels_type panels;
    int surfID;
    bool TEflag; //Surface has sharp trailing edges
    bool LSflag; //Surface is a lifting surface
    std::vector<edge*> trailingEdges;
    Eigen::Vector3d rearStagnationPnt(const Eigen::Vector3d &Vinf, double PG, bodyPanel* &p);
    std::vector<edge*> getTrailingEdges();
    
public:
    surface(const int &surfaceID,geometry* geom);
    
    virtual ~surface();
    
//    surface(const surface& copy) : surfID(copy.surfID)
//    {
//        for (int i=0; i<copy.panels.size(); i++)
//        {
//            panels[i] = new bodyPanel(*copy.panels[i]);
//        }
//    }
    
    virtual void addPanel(bodyPanel* bPan);
    
    void setTEflag() {TEflag = true;}
    void setLSflag() {LSflag = true;}

    std::vector<bodyPanel*> getPanels() const {return panels;}
    int getID() const {return surfID;}
    std::vector<std::pair<Eigen::Vector3d,bodyPanel*>> getStreamlineStartPnts(const Eigen::Vector3d &Vinf,double PG);
    bool sharpTE() {return TEflag;}
    bool isLiftingSurf() {return LSflag;}

    
};

#endif /* defined(__CPanel__surface__) */
