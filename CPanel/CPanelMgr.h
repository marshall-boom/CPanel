/*******************************************************************************
 * Copyright (c) 2015 Chris Satterwhite
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

#ifndef __CPanel__CPanelMgr__
#define __CPanel__CPanelMgr__

#include <stdio.h>
#include "cpCase.h"
#include "cpCaseVP.h"
#include "inputParams.h"
#include "geometry.h"

class caseMgr
{
	using cases_type = std::vector<cpCase *>;
	using cases_index_type = cases_type::size_type;
	using casesVP_type = std::vector<cpCaseVP *>;
	using casesVP_index_type = casesVP_type::size_type;

    inputParams* p;
    geometry* geom;
    Eigen::VectorXi outSpacing;
    
    cases_type cases;
    casesVP_type casesVP;
    
    void setCases();
    void runCases();
    void writeSummary();
    void writeCase(size_t caseNumber, cpCase* c, std::ofstream &outStream);
public:
    caseMgr(inputParams* pp,geometry* ggeom) : p(pp), geom(ggeom)
    {
        setCases();
        runCases();
        writeSummary();
    }
};


#endif /* defined(__CPanel__CPanelMgr__) */
