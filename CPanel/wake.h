//
//  wake.h
//  CPanel
//
//  Created by Chris Satterwhite on 10/15/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#ifndef __CPanel__wake__
#define __CPanel__wake__

#include <stdio.h>
#include <vector>
#include <cmath>
#include <Eigen/Dense>
#include "particle.h"

//#include "wakePanel.h"
//#include "wakeLine.h"

class wakePanel;
class wakeLine;
class geometry;
class particle;

class wake
{
	using wakePanels_type = std::vector<wakePanel *>;
	using wakePanels_index_type = wakePanels_type::size_type;
	using wakeLines_type = std::vector<wakeLine *>;
	using wakeLines_index_type = wakeLines_type::size_type;

    size_t ID;
    geometry* geom;
    std::vector<wakePanel*> wpanels;
    std::vector<wakePanel*> TEpanels;
    std::vector<wakeLine*> wakeLines;
    std::vector<wakePanel*> vortexSheets; //Not currently in use
    double x0,xf,z0,zf,yMin,yMax;
    Eigen::Vector3d normal;
    
    double CL;
    double CD;
    Eigen::VectorXd yLoc;
    Eigen::VectorXd Cl;
    Eigen::VectorXd Cd;
    
    
    void setWakeDimensions();
    wakeLine* findWakeLine(double y);
    double Vradial(Eigen::Vector3d pWake);
//    Eigen::Vector3d Vradial2(Eigen::Vector3d pWake);
    Eigen::Vector3d pntInWake(double x, double y);
//    Eigen::Vector3d pntVel(Eigen::Matrix<double,1,3> POI, Eigen::MatrixXd pntCloud,Eigen::Matrix<bool,Eigen::Dynamic,1> upperFlag, Eigen::Vector3d Vinf);
    
    double particlePntInWakeY(Eigen::Vector3d Spts , particle* SptsP1 , particle* SptsP2);
    double dPhiWeighted(Eigen::Vector3d pt , particle* P1 , particle* P2);
    double stretchFactor( particle* P1 , particle* P2 );


    
public:
    wake(size_t wakeID, geometry* ggeom)
      : ID(wakeID), geom(ggeom), x0(0), xf(0), z0(0), zf(0), yMin(0), yMax(0), CL(0), CD(0) {}
    
    ~wake();
    
//    wake(const wake& copy);
    
    bool isSameWake(wake* other);
    void mergeWake(wake* other);
    void addPanel(wakePanel* wPan);
    void addTEPanel(wakePanel* p);
    void addWakeLine(wakeLine* wl);
    
    std::vector<wakePanel*> getPanels() const {return wpanels;}
    
    void trefftzPlane(double Vinf,double Sref);
    void trefftzPlaneVP(double Vinf,double Sref, std::vector<particle*>* particles, int numSimSteps);
    Eigen::Vector3d lambVectorInt(Eigen::VectorXd &yLoc);
    
    
    
    std::vector<wakeLine*> getWakeLines() {return wakeLines;}
    std::vector<wakePanel*> getTrailingEdgePanels() {return TEpanels;}
    
    double wakeStrength(double y);
    
//    std::vector<wakePanel*> getVortexSheets() {return vortexSheets;}
    double getYMin() {return yMin;}
    double getYMax() {return yMax;}
    double getX0() {return x0;}
    double getXf() {return xf;}
    double getZ0() {return z0;}
    double getZf() {return zf;}
    
    double getCL() {return CL;}
    double getCD() {return CD;}
    
    Eigen::VectorXd getSpanwiseCl() {return Cl;}
    Eigen::VectorXd getSpanwiseCd() {return Cd;}
    Eigen::VectorXd getSpanwisePnts() {return yLoc;}
    
};

#endif /* defined(__CPanel__wake__) */
