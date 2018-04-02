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

#include "liftingSurf.h"

liftingSurf::liftingSurf(int ssurfID,geometry* ggeom) : surface(ssurfID,ggeom), wakeSurf(nullptr)
{
//    wakeSurf = new wake;
}

void liftingSurf::addPanel(bodyPanel* bPan)
{
    panels.push_back(bPan);
}
void liftingSurf::addPanel(wakePanel* wPan)
{
    wakeSurf->addPanel(wPan);
}
std::vector<panel*> liftingSurf::getAllPanels()
{
    std::vector<panel*> allPans;
    std::vector<wakePanel*> wakePans = wakeSurf->getPanels();
    allPans.insert(allPans.end(),panels.begin(),panels.end());
    allPans.insert(allPans.end(),wakePans.begin(),wakePans.end());
    return allPans;
}

