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

#ifndef __CPanel___Unstructured_Panel_Code__cpFile__
#define __CPanel___Unstructured_Panel_Code__cpFile__

#include <string>

struct cpFile
{
    std::string file;
    std::string path;
    std::string name;
    std::string ext;

    cpFile(std::string filename);
    cpFile(const char* filename);

    void changePath(std::string newPath);

private:
    void parsefile();
};

#endif /* defined(__CPanel___Unstructured_Panel_Code__cpFile__) */
