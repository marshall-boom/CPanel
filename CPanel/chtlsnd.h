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

#ifndef __CPanel__chtlsnd__
#define __CPanel__chtlsnd__

#include <stdio.h>
#include <Eigen/Dense>
#include <Eigen/QR>
#include <iostream>

class chtlsnd
{
    Eigen::MatrixXd F;
    Eigen::MatrixXd G;
    Eigen::MatrixXd H;
    
    Eigen::MatrixXi derivSequence(int q, int N);
    Eigen::MatrixXi sortBySum(Eigen::MatrixXi m);
    Eigen::MatrixXi insertRow(const Eigen::MatrixXi &m, const Eigen::MatrixXi &insert, int row);
    
    
public:
    chtlsnd(const Eigen::Matrix<double,1,Eigen::Dynamic> &X0, const Eigen::MatrixXd &Xf, int order, const Eigen::MatrixXd &Xb, const Eigen::MatrixXd &Vb, Eigen::VectorXd V0);
    static int factorial(int i);
    Eigen::MatrixXd getF() {return F;}
    Eigen::MatrixXd getG() {return G;} 
    Eigen::MatrixXd getH() {return H;}
};

#endif /* defined(__CPanel__chtlsnd__) */
