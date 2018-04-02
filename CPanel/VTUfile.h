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

#ifndef __CPanel__VTUfile__
#define __CPanel__VTUfile__

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <Eigen/Dense>
#include <typeinfo>


struct cellDataArray
{
    std::string name;
    Eigen::MatrixXd data;
    
    cellDataArray(std::string nname) : name(nname) {}
};

struct pntDataArray
{
    std::string name;
    Eigen::MatrixXd data;
    
    pntDataArray(std::string nname) : name(nname) {}
};

struct piece
{
    Eigen::MatrixXd pnts;
    Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> connectivity;
    std::vector<cellDataArray> cellData;
    std::vector<pntDataArray> pntData;
};

class VTUfile
{
	using pieces_type = std::vector<piece>;
	using pieces_index_type = pieces_type::size_type;

	pieces_type pieces;
    std::string name;
    
    void printDoubleArray(std::ofstream &f,std::string name,Eigen::MatrixXd array);
    
    void printIntArray(std::ofstream &f,std::string name,Eigen::MatrixXi array);
    
    void printSizeTArray(std::ofstream &f,std::string nname,Eigen::Matrix<size_t, Eigen::Dynamic, Eigen::Dynamic> array);

    void write();
    
public:
    VTUfile(std::string nname, piece p) : name(nname)
    {
        pieces.push_back(p);
        write();
    }
    
    VTUfile(std::string nname, std::vector<piece> ppieces) : pieces(ppieces), name(nname)
    {
        write();
    }
    
    std::string getName() {return name;}
};

#endif /* defined(__CPanel__VTUfile__) */
