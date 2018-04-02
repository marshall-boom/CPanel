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

#ifndef __CPanel__liftingSurf__
#define __CPanel__liftingSurf__

#include <stdio.h>
#include "surface.h"
#include "wake.h"

class liftingSurf : public surface
{
    wake* wakeSurf;
    
public:
    liftingSurf(int surfID,geometry* geom);
    
    ~liftingSurf()
    {
        delete wakeSurf;
    }
    
    void addPanel(bodyPanel* bPan);
    void addPanel(wakePanel* wPan);
    wake* getWake() {return wakeSurf;}
    std::vector<panel*> getAllPanels();
    std::vector<wakePanel*> getWakePanels() {return wakeSurf->getPanels();}
};

#endif /* defined(__CPanel__liftingSurf__) */
