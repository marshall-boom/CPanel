//
//  influenceTests.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 1/24/15.
//  Copyright (c) 2015 Chris Satterwhite. All rights reserved.
//

#include "influenceTests.h"

void influenceTests::panelVConstTest(Eigen::MatrixXd testingPoints, Eigen::MatrixXd velocityData, bodyPanel* tPan)
{
    // Test panel vertical velocity survey
    Eigen::MatrixXd cpanelVelSurv(testingPoints.rows(),3), survDifference(testingPoints.rows(),3);
    for(int i=0; i<cpanelVelSurv.rows(); i++)
    {
        cpanelVelSurv.row(i) = tPan->panelV(testingPoints.row(i))*4*M_PI;
//        Eigen::Vector3d survpnt = testingPoints.row(i);
//        Eigen::Vector3d pntVel = tPan->pntVInf(testingPoints.row(i))*4*M_PI;
//        Eigen::Vector3d constVel = tPan->panelV(testingPoints.row(i))*4*M_PI;
//        Eigen::Vector3d expected = velocityData.row(i);
        
        for(int j=0; j<testingPoints.cols(); j++)
        {
            survDifference(i,j) = std::abs(cpanelVelSurv(i,j)-velocityData(i,j))*2/(velocityData(i,j)+cpanelVelSurv(i,j))*100;
        }
    }
    
//    std::cout << "\n\ntesting pts " << testingPoints << std::endl;
//    std::cout << "\n\nvelocity dat" << velocityData  << std::endl;
//    std::cout << "\n\ncpanel Vel  " << cpanelVelSurv << std::endl;
//    
//    std::cout << "\n\nsurv difference" << survDifference << std::endl;
//    
    std::cout << maxCoefficient(survDifference);

    

}

void influenceTests::panelVPntTest(Eigen::MatrixXd testingPoints, Eigen::MatrixXd velocityData, bodyPanel* tPan)
{
    // Test panel vertical velocity survey
    Eigen::MatrixXd cpanelVelSurv(testingPoints.rows(),3), survDifference(testingPoints.rows(),3);
    for(int i=0; i<cpanelVelSurv.rows(); i++)
    {
        cpanelVelSurv.row(i) = tPan->pntVInf(testingPoints.row(i))*4*M_PI;
        Eigen::Vector3d survpnt = testingPoints.row(i);
        Eigen::Vector3d pntVel = tPan->pntVInf(testingPoints.row(i))*4*M_PI;
        Eigen::Vector3d constVel = tPan->panelV(testingPoints.row(i))*4*M_PI;
        Eigen::Vector3d expected = velocityData.row(i);
        
        for(int j=0; j<testingPoints.cols(); j++)
        {
            survDifference(i,j) = std::abs(cpanelVelSurv(i,j)-velocityData(i,j))*2/(velocityData(i,j)+cpanelVelSurv(i,j))*100;
        }
    }
//    std::cout << "\n\ntesting pts " << testingPoints << std::endl;
//    std::cout << "\n\nvelocity dat" << velocityData  << std::endl;
//    std::cout << "\n\ncpanel Vel  " << cpanelVelSurv << std::endl;
//
//    std::cout << "\n\nsurv difference" << survDifference << std::endl;
//

    std::cout << maxCoefficient(survDifference);
    
    
}


double influenceTests::maxCoefficient(Eigen::MatrixXd mat)
{
    double maxCoef = 0;
    for(int i=0; i< mat.rows(); i++){
        for (int j=0; j<mat.cols(); j++) {
            if (isnormal(mat(i,j))) {
                if(std::abs(mat(i,j)) > maxCoef){
                    maxCoef = std::abs(mat(i,j));
                }
            }
        }
    }
    
    return maxCoef;
}


void influenceTests::velocityComparer(DoubVelData* dat)
{
    
    

influenceTests* test = new influenceTests();

// Constant doublet vertical survey
std::cout << "\nVertical const doublet velocity max difference in percent %: " << std::flush;
test->panelVConstTest(dat->vertSurveyPts, dat->vertVelConst, dat->testPan);

// Constant doublet diagonal survey
std::cout << "\nDiagonal const doublet velocity max difference in percent %: " << std::flush;
test->panelVConstTest(dat->diagSurveyPts, dat->diagVelConst, dat->testPan);

// Constant doublet median survey
std::cout << "\nMedian const doublet velocity max difference in percent %: " << std::flush;
test->panelVConstTest(dat->MedSurveyPts, dat->medVelConst, dat->testPan);

// Point doublet vertical survey
std::cout << "\nVertical pnt doublet velocity max difference in percent %: " << std::flush;
test->panelVPntTest(dat->vertSurveyPts, dat->vertVelPnt, dat->testPan);

// Point doublet diagonal survey
std::cout << "\nDiagonal pnt doublet velocity max difference in percent %: " << std::flush;
test->panelVPntTest(dat->diagSurveyPts, dat->diagVelPnt, dat->testPan);

// Point doublet median survey
std::cout << "\nMedian pnt doublet velocity max difference in percent %: " << std::flush;
test->panelVPntTest(dat->MedSurveyPts, dat->medVelPnt, dat->testPan);




}


std::vector<particle*> influenceTests::randomParticleGenerator(int numParts){
    // Function generates random particles with position and strength within 0 and 1 for each component
    std::vector<particle*> particles;
    
    double radius = .5/numParts; // ~ avg. dist btw parts
    for (int i=0; i<numParts; i++) {
        
        // Create random position and strength
        Eigen::Vector3d pos, strength;
        for (int j=0; j<pos.size(); j++)
        {
            pos(j) = pick_a_number(0,1);
            strength(j) = pick_a_number(0,1);
        }
        
        particle* p = new particle(pos,strength,radius, {0,0,0}, {0,0,0});
        particles.push_back(p);
    }
    
    return particles;
}


void influenceTests::FMMtests(){
    
    std::vector<particle*> testParts = randomParticleGenerator(100);
    
//    BarnesHutSpeedTest(testParts);
//    barnesHutTest(testParts);
//    interactionListTest(testParts);
    speedTestGraphComparison();
}

void influenceTests::barnesHutTest(std::vector<particle*> testParts){
    // This test creates an octree with random particles and then computes the velocity influnce by both the brute force method and then the Barnes Hut algorithm
    
    // Test params
    Eigen::Vector3d POI = {0,0,0};

    
    // Manually sum up particle influences
    std::cout << "Brute force method..." << std::endl;
    Eigen::Vector3d allPartsVel = Eigen::Vector3d::Zero();
    
    clock_t partBegin = clock();
    for (int i=0; i<testParts.size(); i++)
    {
        allPartsVel+=testParts[i]->partVelInfl(POI);
    }
    clock_t partEnd = clock();
    std::cout << "Individual Particle Calc Time = " <<  (partEnd-partBegin) << " clocks" << std::endl;

    
    
    // Barnes-Hut
    std::cout << "Barnes-Hut method:" << std::endl;
    
    std::cout << "Building octree..." << std::endl;
    clock_t octBuild_start = clock();
        particleOctree partTree;
        partTree.setMaxMembers(10);
        partTree.addData(testParts);
    clock_t octBuild_end = clock();
    std::cout << "Treecode build time = " <<  (octBuild_end-octBuild_start) << " clocks" << std::endl;

    // Writing octree output (not inlcuded in time)
    std::string file_name = "/Users/C_Man/Desktop/CPanelDevelopment/FMM/testTree.txt";
    octreeFile* oct;
    oct = new octreeFile(file_name,&partTree);
    
    
    std::cout << "Calculating expansions..." << std::endl;
    clock_t expBuild_start = clock();
        particleFMM FMM;
        FMM.build(&partTree);
    clock_t expBuild_end = clock();
    std::cout << "Particle Expansions build time = " <<  (expBuild_end-expBuild_start) << " clocks" << std::endl;

    
    std::cout << "Velocity using Barnes-Hut..." << std::endl;
    clock_t bhBegin = clock();
        Eigen::Vector3d barnesHutVel = FMM.barnesHut(POI);
    clock_t bhEnd = clock();
    std::cout << "Barnes Hut time = " <<  (bhEnd-bhBegin) << " clocks" << std::endl;


    std::cout << "Brute force method vel = {"<<allPartsVel.x()<<", "<<allPartsVel.y()<<","<<allPartsVel.z()<<"}"<<std::endl;
    std::cout << "Barnes-Hut velocity    = {"<<barnesHutVel.x()<<", "<<barnesHutVel.y()<<","<<barnesHutVel.z()<<"}"<<std::endl;

    double percDiffX = (allPartsVel.x()-barnesHutVel.x())/allPartsVel.x()*100;
    double percDiffY = (allPartsVel.y()-barnesHutVel.y())/allPartsVel.y()*100;
    double percDiffZ = (allPartsVel.z()-barnesHutVel.z())/allPartsVel.z()*100;

    std::cout << "% diff = {" << percDiffX<<"%, "<<percDiffY<<"%, "<<percDiffZ<<"%} "<<std::endl;
    
    std::cout <<  std::endl;
    
}

void influenceTests::BarnesHutSpeedTest(std::vector<particle*> testParts){
    // Function goes through each interaction like what would happen to advance a cloud of particles.
    std::cout << "------------Number of test particles = " << testParts.size() << "--------------" << std::endl;
    
    
    // Brute Force
//    std::cout << "Brute force method..." << std::endl;
    clock_t bruteForceBegin = clock();
    std::vector<Eigen::Vector3d> velOnVec;
    for (int i=0; i<testParts.size(); i++)
    {
        Eigen::Vector3d velOn = Eigen::Vector3d::Zero();
        for (int j=0; j<testParts.size(); j++)
        {
            velOn += testParts[j]->partVelInfl(testParts[i]->pos);
        }
        velOnVec.push_back(velOn);
    }
    clock_t bfTime = (clock()-bruteForceBegin);
    std::cout << "Brute force method = " <<  bfTime << " clocks" << std::endl;

    
    
    // Barnes Hut
//    clock_t barnesHutBegin = clock();
    
    // Build octree
    particleOctree partTree;
    partTree.setMaxMembers(1);
    partTree.addData(testParts);
    
    // Build Expansions
    particleFMM FMM;
    FMM.build(&partTree);
    
    std::vector<Eigen::Vector3d> velOnVecBH;
    clock_t barnesHutBegin = clock();
    for (int i=0; i<testParts.size(); i++)
    {
//        clock_t oneBH = clock();
        velOnVecBH.push_back(FMM.barnesHut(testParts[i]->pos));
//        clock_t oneBHend = clock();
//        std::cout << "one BH calc: " << oneBHend-oneBH << std::endl;
    }
    clock_t bhTime = (clock()-barnesHutBegin);
    std::cout << "Barnes Hut method  = " << bhTime  << " clocks" << std::endl;

    // Converting clock time to double so can do math with it
    std::string bhTimeString = std::to_string(bhTime);
    std::string::size_type sz;     // alias of size_t
    double bhTimeD = std::stod (bhTimeString, &sz);
    
    std::string bfTimeString = std::to_string(bfTime);
    double bfTimeD = std::stod (bfTimeString, &sz);
    
    std::cout << "Speed up  = " << (bfTimeD/bhTimeD)-1 << " times faster" << std::endl;
}

void influenceTests::speedTestGraphComparison()
{
    int a = 500;
    int b = 10000;
    int n = 20;
    std::vector<int> numParts = linspace(a, b, n);
    
    for (int i=0; i<numParts.size(); i++)
    {
        std::vector<particle*> testParts = randomParticleGenerator(numParts[i]);
        BarnesHutSpeedTest(testParts);
    }

}


std::vector<int> influenceTests::linspace(int a, int b, int n)
{
    std::vector<int> linSpaced;
    
    int step = (b-a)/(n-1);
    
    for (int i=0; i<n; i++) {
        linSpaced.push_back(a+step*i);
    }
    return linSpaced;
}

std::default_random_engine & global_urng( )
{
    static std::default_random_engine u{};
    return u;
}

double influenceTests::pick_a_number( double from, double upto )
{
    static std::uniform_real_distribution<> d{};
    using parm_t = decltype(d)::param_type;
    
    return d( global_urng(), parm_t{from, upto} );
}





























//influenceTests::influenceTests()
//{
//    createTestPanel();
//        TEST_ADD(influenceTests::testPntClose)
//        TEST_ADD(influenceTests::testPntFar)
////        TEST_ADD(influenceTests::testCollocationPnt)
//    //    TEST_ADD(influenceTests::testDoublet)
//    //    TEST_ADD(influenceTests::testSource)
//}
//
//void influenceTests::createTestPanel()  //Connor: I think this is old and doesn't work with new panels
//{
//    testNodes.resize(3,3);
//    testNodes(0,0) = 1;
//    testNodes(0,1) = 0;
//    testNodes(0,2) = 0;
//    testNodes(1,0) = -1;
//    testNodes(1,1) = 1;
//    testNodes(1,2) = 0;
//    testNodes(2,0) = -1;
//    testNodes(2,1) = -1;
//    testNodes(2,2) = 0;
//    Eigen::Vector3i verts;
//    verts << 0,1,2;
//    Eigen::Vector3d bezNorm;
//    bezNorm << 0,0,1;
//    std::vector<edge*> edges;
//    testPan = new bodyPanel(&testNodes,edges,bezNorm,1,false);
//}
//
//
//void influenceTests::testPntClose()
//{
//    Eigen::Vector3d POI;
//    POI << 2,2,2; // |POI|/(side length) ~ 1.2 so panel formulation will be used;
//
//    double phiSrc,phiDub;
//    Eigen::Vector3d vSrc,vDub;
//    phiSrc = 0;
//    phiDub = 0;
//    vSrc = Eigen::Vector3d::Zero();
//    vDub = Eigen::Vector3d::Zero();
//    testPan->panelPhiInf(POI, phiSrc, phiDub);
//    testPan->panelVInf(POI, vSrc, vDub);
//    TEST_ASSERT(isEqual(phiSrc,0.043411,5))
//    TEST_ASSERT(isEqual(vSrc(0),0.007363,5))
//    TEST_ASSERT(isEqual(vSrc(1),0.006451,5))
//    TEST_ASSERT(isEqual(vSrc(2),0.006645,5))
//
//    TEST_ASSERT(isEqual(phiDub,0.006645,5))
//    TEST_ASSERT(isEqual(vDub(0),-0.003381,5))
//    TEST_ASSERT(isEqual(vDub(1),-0.003026,5))
//    TEST_ASSERT(isEqual(vDub(2),0.000158,5))
//}
//
//void influenceTests::testPntFar()
//{
//    Eigen::Vector3d POI;
//    POI << 7,7,7; // |POI|/(side length) ~ 5.5 so far field approximation will be used;
//
//    double phiSrc,phiDub;
//    Eigen::Vector3d vSrc,vDub;
//    phiSrc = 0;
//    phiDub = 0;
//    vSrc = Eigen::Vector3d::Zero();
//    vDub = Eigen::Vector3d::Zero();
//    testPan->panelPhiInf(POI, phiSrc, phiDub);
//    testPan->panelVInf(POI, vSrc, vDub);
//    TEST_ASSERT(isEqual(phiSrc,0.012919,5))
//    TEST_ASSERT(isEqual(vSrc(0),0.000625,5))
//    TEST_ASSERT(isEqual(vSrc(1),0.000596,5))
//    TEST_ASSERT(isEqual(vSrc(2),0.000596,5))
//    
//    TEST_ASSERT(isEqual(phiDub,0.0005958,5))
//    TEST_ASSERT(isEqual(vDub(0),-0.000086,6))
//    TEST_ASSERT(isEqual(vDub(1),-0.000082,6))
//    TEST_ASSERT(isEqual(vDub(2),0.0000027,6))
//}

//void influenceTests::testCollocationPnt()
//{
//    Eigen::Vector3d POI;
//    POI << -1.0/3,0,0;
//
//    double phiSource;
//    Eigen::Vector3d vSource;
//    phiSource = p.sourcePhi(1,POI,vertsLocal);
//    vSource = p.sourceV(1,POI,vertsLocal);
//    TEST_ASSERT(phiSource == 0)
//    TEST_ASSERT(vSource(0) == 0)
//    TEST_ASSERT(vSource(1) == 0)
//    TEST_ASSERT(vSource(2) == 0.5)
//
//    double phiDoublet;
//    Eigen::Vector3d vDoublet;
//    phiDoublet = p.doubletPhi(1,POI,vertsLocal);
//    vDoublet = p.doubletV(1,POI,vertsLocal);
//    TEST_ASSERT(phiDoublet == -0.5)
//    TEST_ASSERT(vDoublet(0) == 0)
//    TEST_ASSERT(vDoublet(1) == 0)
//    TEST_ASSERT(isEqual(vDoublet(2),1.432394487,9))
//}

//bool influenceTests::isEqual(double var1, double var2, int decimalPrecision)
//{
//    double eps = pow(10,-decimalPrecision);
//    return (std::abs(var1-var2)<eps);
//    
//}

//void influenceTests::testDoublet()
//{
//    std::string path = "/Users/Chris/Desktop/Thesis/Code/Geometry and Solution Files/";
//    std::string filename = "doubletCheck.txt";
//    std::ifstream fid;
//    fid.open(path+filename);
//    if (fid.is_open())
//    {
//        Eigen::MatrixXd nodes;
//        nodes.resize(3,3);
//        Eigen::Vector3i verts;
//        verts << 0,1,2;
//        Eigen::Vector3d col;
//        Eigen::MatrixXd pnts(9,3);
//        Eigen::VectorXd pot(9);
//        Eigen::MatrixXd vel(9,3);
//
//        for (int i=0; i<3; i++)
//        {
//            fid >> nodes(i,0) >> nodes(i,1) >> nodes(i,2);
//        }
//        fid >> col(0) >> col(1) >> col(2);
//        for (int i=0; i<9; i++)
//        {
//            fid >> pnts(i,0) >> pnts(i,1) >> pnts(i,2);
//        }
//        for (int i=0; i<9; i++)
//        {
//            fid >> pot(i);
//        }
//        for (int i=0; i<9; i++)
//        {
//            fid >> vel(i,0) >> vel(i,1) >> vel(i,2);
//        }
//        vel *= -1; //Different sign convention
//        Eigen::VectorXd phiCalc(9);
//        Eigen::MatrixXd Vcalc(9,3);
//        bodyPanel pan(verts,&nodes,1);
//        Eigen::Matrix3d localVerts = pan.getLocalVerts();
//
//        for (int i=0; i<9; i++)
//        {
//            phiCalc(i) = pan.doubletPhi(1,pnts.row(i),localVerts);
//            TEST_ASSERT(isEqual(phiCalc(i),pot(i),5))
//            Vcalc.row(i) = pan.doubletV(true, pnts.row(i));
//            for (int j=0; j<3; j++)
//            {
//                std::cout << Vcalc(i,j) << "," << vel(i,j) << std::endl;
//                TEST_ASSERT(isEqual(Vcalc(i,j),vel(i,j),5))
//            }
//        }
//    }
//    else
//    {
//        std::cout << "ERROR : File not found" << std::endl;
//        exit(EXIT_FAILURE);
//    }
//
//}
//
//void influenceTests::testSource()
//{
//    std::string path = "/Users/Chris/Desktop/Thesis/Code/Geometry and Solution Files/";
//    std::string filename = "sourceCheck.txt";
//    std::ifstream fid;
//    fid.open(path+filename);
//    if (fid.is_open())
//    {
//        Eigen::MatrixXd nodes;
//        nodes.resize(3,3);
//        Eigen::Vector3i verts;
//        verts << 0,1,2;
//        Eigen::Vector3d col;
//        Eigen::MatrixXd pnts(9,3);
//        Eigen::VectorXd pot(9);
//        Eigen::MatrixXd vel(9,3);
//
//        for (int i=0; i<3; i++)
//        {
//            fid >> nodes(i,0) >> nodes(i,1) >> nodes(i,2);
//        }
//        fid >> col(0) >> col(1) >> col(2);
//        for (int i=0; i<9; i++)
//        {
//            fid >> pnts(i,0) >> pnts(i,1) >> pnts(i,2);
//        }
//        for (int i=0; i<9; i++)
//        {
//            fid >> pot(i);
//        }
//        for (int i=0; i<9; i++)
//        {
//            fid >> vel(i,0) >> vel(i,1) >> vel(i,2);
//        }
//        pot *= -1; //Different sign convention;
//        Eigen::VectorXd phiCalc(9);
//        Eigen::MatrixXd Vcalc(9,3);
//        bodyPanel pan(verts,&nodes,1);
//        Eigen::Matrix3d localVerts = pan.getLocalVerts();
//
//        for (int i=0; i<9; i++)
//        {
//            phiCalc(i) = pan.sourcePhi(1,pnts.row(i),localVerts);
//            TEST_ASSERT(isEqual(phiCalc(i),pot(i),5))
//            Vcalc.row(i) = pan.sourceV(true, pnts.row(i));
//            for (int j=0; j<3; j++)
//            {
//                TEST_ASSERT(isEqual(-Vcalc(i,j),vel(i,j),5))
//            }
//        }
//    }
//    else
//    {
//        std::cout << "ERROR : File not found" << std::endl;
//        exit(EXIT_FAILURE);
//    }
//
//}
