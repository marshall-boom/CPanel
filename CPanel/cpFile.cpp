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

#include "cpFile.h"

#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

cpFile::cpFile(std::string filename) : file(filename)
{
    parsefile();
}

cpFile::cpFile(const char* filename)
{
    file = filename;
    parsefile();
}

void cpFile::parsefile()
{
    std::size_t pathEnd = file.find_last_of("/")+1;
    std::size_t nameEnd = file.find_last_of(".");
    bool relPath = false;
    if (file.substr(0,1) != "/")
    {
        relPath = true;
    }
    // Get extension
    ext = file.substr(nameEnd,file.size()-nameEnd);
    name = file.substr(pathEnd,nameEnd-pathEnd);

    if (!relPath)
    {
        path = file.substr(0,pathEnd);
    }
    else
    {
        std::string relativePath = file.substr(0,pathEnd);
        std::string currentPath = boost::filesystem::current_path().string();

        std::stringstream ss;
        ss << currentPath << "/" << relativePath;
        path = ss.str();
        std::stringstream ss2;
        ss2 << path << name << ext;
        file = ss2.str();
    }
}

void cpFile::changePath(std::string newPath)
{
    path = newPath;
    std::stringstream temp;
    temp << newPath << "/" << name << ext;
    file = temp.str();
}
