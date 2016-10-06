//
//  main.cpp
//  CPanelTests
//
//  Created by Chris Satterwhite on 1/22/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#include <iostream>
#include "geometryTests.h"
#include "influenceTests.h"
#include "geometry.h"
#include "cpFile.h"
#include "inputParams.h"
#include "wakePanel.h"
#include "DoubVelData.h"
#include "surface.h"
#include <Eigen/Dense>
#include <cmath>
#include <string>
#include <fstream>



DoubVelData* readDataFile(std::string fid);
void createTestPanel(DoubVelData*, geometry* dGeom);


int main(int argc, const char * argv[])
{
    // Connor is re-writing all of these tests becuase they were made before CPanel re-write so they don't work anymore
    //
    // An test data file was written from Calebretta's velocity influence functions.
    // Input file format is under CPanelDevelopment->PanelVInfluenceTest->OutputFileFormat
    
    
    cpFile inFile(argv[1]);
    inputParams inData(&inFile);
    inData.set();
    geometry geom(&inData);
    
    
    int numFiles = 7;
    
    for (int i=0; i<numFiles; i++) {
        std::cout << "\n\nFile Number " << i+1 << ":" << std::endl;
        
//        std::string fid = "/Users/C_Man/Desktop/CPanelDevelopment/PanelVInfluenceTest/doubletVelTest" + std::to_string(i+1) + ".txt";
//        std::string fid = "/Users/C_Man/Desktop/CPanelDevelopment/PanelVInfluenceTest/sourceVelTest" + std::to_string(i+1) + ".txt";
        std::string fid = "/Users/C_Man/Desktop/CPanelDevelopment/PanelVInfluenceTest/sourceOutput.txt";
        
        DoubVelData *dat;
        dat = readDataFile(fid);
        createTestPanel(dat, &geom);
        
        GeomTests *geomTest = new GeomTests(dat->testPan);
        geomTest->conductTests(dat->testPan, &geom);
        
        
        influenceTests* test = new influenceTests();
        
        test->FMMtests();
        test->velocityComparer(dat);
        
    }
    
}




DoubVelData* readDataFile(std::string fid)
{

    DoubVelData* d = new DoubVelData();

    std::ifstream velData;
    velData.open(fid);
    
    if(!velData.is_open())
    {
        std::cout << "Error: File not opened" << std::endl;
        exit(EXIT_FAILURE);
    }
    
    velData >> d->colocationPnt.x() >> d->colocationPnt.y() >> d->colocationPnt.z();
    velData >> d->panelSideLength;
    velData >> d->n1.x() >> d->n1.y() >> d->n1.z();
    velData >> d->n2.x() >> d->n2.y() >> d->n2.z();
    velData >> d->n3.x() >> d->n3.y() >> d->n3.z();
    velData >> d->n4.x() >> d->n4.y() >> d->n4.z();
    velData >> d->mu;
    velData >> d->sigma;

    // Read in vertical velocity survey points
    double dummy;
    velData >> dummy; // Assigning to a double then converting to int. so cpp recognizes it properly
    d->numVertPts = (int) dummy;
    d->vertSurveyPts.resize(d->numVertPts,3);

    for(int i=0; i<d->vertSurveyPts.rows(); i++)
    {
        velData >> d->vertSurveyPts(i,0) >> d->vertSurveyPts(i,1) >> d->vertSurveyPts(i,2);
    }
    
    // Read in diagonal velocity survey points
    velData >> dummy; // Assigning to a double then converting to int. so cpp recognizes it properly
    d->numDiagPts = (int) dummy;
    d->diagSurveyPts.resize(d->numDiagPts,3);
    for(int i=0; i<d->diagSurveyPts.rows(); i++)
    {
        velData >> d->diagSurveyPts(i,0) >> d->diagSurveyPts(i,1) >> d->diagSurveyPts(i,2);
    }
    
    // Read in median velocity survey points
    velData >> dummy; // Assigning to a double then converting to int. so cpp recognizes it properly
    d->numMedPts = (int) dummy;
    d->MedSurveyPts.resize(d->numMedPts,3);
    for(int i=0; i<d->MedSurveyPts.rows(); i++)
    {
        velData >> d->MedSurveyPts(i,0) >> d->MedSurveyPts(i,1) >> d->MedSurveyPts(i,2);
    }
    
    // Vertical vel. survey: point doublet then const. strength doub
    d->vertVelPnt.resize(d->numVertPts,3);
    for(int i=0; i<d->vertVelPnt.rows(); i++)
    {
        velData >> d->vertVelPnt(i,0) >> d->vertVelPnt(i,1) >> d->vertVelPnt(i,2);
    }
    d->vertVelConst.resize(d->numVertPts,3);
    for(int i=0; i<d->vertVelConst.rows(); i++)
    {
        velData >> d->vertVelConst(i,0) >> d->vertVelConst(i,1) >> d->vertVelConst(i,2);
    }
    
    // Diagonal vel. survey: point doublet then const. strength doub
    d->diagVelPnt.resize(d->numDiagPts,3);
    for(int i=0; i<d->diagVelPnt.rows(); i++)
    {
        velData >> d->diagVelPnt(i,0) >> d->diagVelPnt(i,1) >> d->diagVelPnt(i,2);
    }
    d->diagVelConst.resize(d->numDiagPts,3);
    for(int i=0; i<d->diagVelConst.rows(); i++)
    {
        velData >> d->diagVelConst(i,0) >> d->diagVelConst(i,1) >> d->diagVelConst(i,2);
    }
    
    // Median vel. survey: point doublet then const. strength doub
    d->medVelPnt.resize(d->numMedPts,3);
    for(int i=0; i<d->medVelPnt.rows(); i++)
    {
        velData >> d->medVelPnt(i,0) >> d->medVelPnt(i,1) >> d->medVelPnt(i,2);
    }
    
    d->medVelConst.resize(d->numMedPts,3);
    for(int i=0; i<d->medVelConst.rows(); i++)
    {
        velData >> d->medVelConst(i,0) >> d->medVelConst(i,1) >> d->medVelConst(i,2);
    }
    
    velData.close();
    
    return d;
}


void createTestPanel(DoubVelData* d,geometry* dGeom){
    
    // Make the test panel
    // Will use wake panel for the doublet verification
    // Will use body panel for source verification
    
    
    // Create Nodes
    d->nodes.push_back(new cpNode(d->n1, 0));
    d->nodes.push_back(new cpNode(d->n2, 0));
    d->nodes.push_back(new cpNode(d->n3, 0));
    d->nodes.push_back(new cpNode(d->n4, 0));
    
    // Create Edges
    for(int  i=0; i<d->nodes.size(); i++)
    {
        if(i == d->nodes.size()-1)
        {
            d->pEdges.push_back(new edge(d->nodes[i], d->nodes[0], dGeom));
        }
        else
        {
            d->pEdges.push_back(new edge(d->nodes[i], d->nodes[i+1], dGeom));
        }
    }
    
    // Dummy normal vector (if zero, panel constructor automatically calculates it)
    Eigen::Vector3d normVec = Eigen::Vector3d::Zero();
    wake* dummyWake;
    
    const int surfID = 0;
    surface* dummySurf = new surface(surfID, dGeom);
    
    // Create test panel
//    d->testPan = new wakePanel(d->nodes, d->pEdges, normVec, dummyWake, 0);
    d->testPan = new bodyPanel(d->nodes, d->pEdges, normVec, dummySurf, surfID);
    d->testPan->setMu(d->mu);
    d->testPan->setSigma(Eigen::Vector3d::Zero(), d->sigma); //looking at code, this should work.
    
    
    
}







//    Test::Suite tests;
////    tests.add(std::auto_ptr<Test::Suite>(new GeomTests));
//    tests.add(std::auto_ptr<Test::Suite>(new influenceTests));
//
//    Test::TextOutput output(Test::TextOutput::Verbose);
//    return tests.run(output);


