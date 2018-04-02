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
 *    David D. Marshall - porting to GoogleTest
 ******************************************************************************/

#include "geometryTests.h"



void GeomTests::conductTests(bodyPanel* bPan, geometry* geom)
{
    nearFilamentTest();
//    wPanelVectorTest(geom);
//    wPanelEdgesInOrderTest(geom);
}


void GeomTests::nearFilamentTest()
{
    // test will output
    Eigen::Vector3d p1 = {-1,0,0};
    Eigen::Vector3d p2 = {1,0,0};
    double minX = -1.3;
    double maxX = 1.3;
    double minY = -0.2;
    double maxY = 0.2;
    const int nX = 60;
    const int nY = 20;
    std::vector<double> xVec = linspace(minX, maxX, nX);
    std::vector<double> yVec = linspace(minY, maxY, nY);
    
    
    
//    Eigen::Matrix<bool, nX, nY> isNearField;
    std::cout << "\n\nclose all; clear all; clc; format compact; hold on;" << std::endl;
    std::cout << "plot(" << p1.x() <<","<< p1.y() <<",'bs');" << std::endl;
    std::cout << "plot(" << p2.x() <<","<< p2.y() <<",'bs');" << std::endl;
    std::cout << "plot(["<<p1.x()<<","<<p2.x()<<"],["<<p1.y()<<","<<p2.y()<<"],'k','lineWidth',2);" << std::endl;
//    std::cout << "plot([" << std::to_string(minX) << " "<<std::to_string(maxX)<<"],[-.1 -.1])" << std::endl;
//    std::cout << "plot([" << std::to_string(minX) << " "<<std::to_string(maxX)<<"],[.1 .1])" << std::endl;
    
    for (int i=0; i<nX; i++) {
        for (int j=0; j<nY; j++){
            Eigen::Vector3d testPnt = {xVec[i],yVec[j],0};
            
            bool isNear = testPanel->nearFilamentCheck(p1, p2, testPnt);
            if(isNear){
                std::cout << "plot(" << testPnt.x() <<","<< testPnt.y() <<",'r*');" << std::endl;
            }else{
                std::cout << "plot(" << testPnt.x() <<","<< testPnt.y() <<",'b*');" << std::endl;
            }
        }
    }
    std::cout << "set(gcf,'color','white');" << std::endl;
    std::cout << "axis equal tight;" << std::endl;
    std::cout << "title('Vortex Filament Near Field Test \delta = 5%')" << std::endl;
    
}

void GeomTests::wPanelVectorTest(geometry* geom){

    std::vector<wakePanel*> wPans = (*geom->getWakePanels());
    
    // Plots the corner points for pts in order. node zero has smallest x and largest y value
    
    std::cout << "\n\nclose all; clear all; clc; format compact; hold on;" << std::endl;
    for (int i=0; i<wPans.size(); i++) {
        std::vector<cpNode*> pNodes = wPans[i]->pointsInOrder();
        Eigen::Vector3d node = pNodes[0]->getPnt();
        std::cout << "plot3("<<node.x() <<","<< node.y() <<","<< node.z() <<",'rs');"<<std::endl;
        node = pNodes[1]->getPnt();
        std::cout << "plot3("<<node.x() <<","<< node.y() <<","<< node.z() <<",'bx');"<<std::endl;
        node = pNodes[2]->getPnt();
        std::cout << "plot3("<<node.x() <<","<< node.y() <<","<< node.z() <<",'m+');"<<std::endl;
        node = pNodes[3]->getPnt();
        std::cout << "plot3("<<node.x() <<","<< node.y() <<","<< node.z() <<",'ko');"<<std::endl;
        
    }
    
    std::cout << "axis equal; set(gcf,'color','white');" << std::endl;
}

void GeomTests::wPanelEdgesInOrderTest(geometry* geom){
    // Plot all of the midpoints of the panel edges in the correct order
    std::vector<wakePanel*> wPans = (*geom->getWakePanels());
    
    std::cout << "\n\nclose all; clear all; clc; format compact; hold on;" << std::endl;
    for (int i=0; i<wPans.size(); i++) {
        std::vector<edge*> pEdges = wPans[i]->edgesInOrder();
        Eigen::Vector3d node = pEdges[0]->getMidPoint();
        std::cout << "plot3("<<node.x() <<","<< node.y() <<","<< node.z() <<",'rs');"<<std::endl;
        node = pEdges[1]->getMidPoint();
        std::cout << "plot3("<<node.x() <<","<< node.y() <<","<< node.z() <<",'bx');"<<std::endl;
        node = pEdges[2]->getMidPoint();
        std::cout << "plot3("<<node.x() <<","<< node.y() <<","<< node.z() <<",'m+');"<<std::endl;
        node = pEdges[3]->getMidPoint();
        std::cout << "plot3("<<node.x() <<","<< node.y() <<","<< node.z() <<",'ko');"<<std::endl;
        
    }
    std::cout << "axis equal; set(gcf,'color','white');" << std::endl;
    
}


std::vector<double> GeomTests::linspace(double a, double b, int n)
{
    std::vector<double> linSpaced;
    
    double step = (b-a)/(n-1);
    
    for (int i=0; i<n; i++) {
        linSpaced.push_back(a+step*i);
    }
    return linSpaced;
}





//GeomTests::GeomTests()
//{
//    makeTestTriFile();
//    TEST_ADD(GeomTests::test_readTri);
//    TEST_ADD(GeomTests::test_neighborSearch);
//    TEST_ADD(GeomTests::test_createOctree);
//}
//
//void GeomTests::test_readTri()
//{
//    // Tests createSurvaces method in geometry class, as well as addPanel method in the surface class.
//    // Create Geometry
//    Eigen::Vector3d cg = Eigen::Vector3d::Zero();
//    testGeom = new geometry(testTriFile,false,1,1,1,cg);
//    
//    // 5 body panels should be created in a nonlifting surface with surfID of 2;
//    std::vector<surface*> NLsurfs = testGeom->getNonLiftingSurfs();
//    TEST_ASSERT(NLsurfs.size() == 1);
//    TEST_ASSERT(NLsurfs[0]->getID() == 2);
//    TEST_ASSERT(NLsurfs[0]->getPanels().size() == 5);
//    
//    // 6 body panels should be created in a lifting surface with surfID of 1.  4 corresponding wake panels should have been added to its wake
//    std::vector<liftingSurf*> Lsurfs = testGeom->getLiftingSurfs();
//    TEST_ASSERT(Lsurfs.size() == 1);
//    TEST_ASSERT(Lsurfs[0]->getID() == 1);
//    TEST_ASSERT(Lsurfs[0]->getPanels().size() == 6); // Just Wing Panels
//    TEST_ASSERT(Lsurfs[0]->getWakePanels().size() == 4);
//    TEST_ASSERT(Lsurfs[0]->getAllPanels().size() == 10); // Wing+Wake Panels
//}
//
//void GeomTests::test_neighborSearch()
//{
//    std::vector<bodyPanel*> *bPanels = testGeom->getBodyPanels();
//    Eigen::VectorXi nNeighbs(11);
//    nNeighbs << 2,1,2,3,1,2,2,2,3,1,1;
//    for (int i=0; i<(*bPanels).size(); i++)
//    {
//        TEST_ASSERT_MSG((*bPanels)[i]->getNeighbors().size() == nNeighbs(i), "Incorrect Number of Neighboring Panels Set" );
//    }
//    
//    std::vector<wakePanel*> *wPanels = testGeom->getWakePanels();
//    Eigen::VectorXi TEpans(4);
//    TEpans << 0,1,0,1;
//    for (int i=0; i<wPanels->size(); i++)
//    {
//        TEST_ASSERT_MSG((*wPanels)[i]->isTEpanel() == TEpans(i), "Trailing Edge Panels Not Set Properly");
//    }
//    
//    TEST_ASSERT_MSG((*wPanels)[1]->getUpper() == (*bPanels)[0], "Upper Parent Panel not set properly")
//    
//    TEST_ASSERT_MSG((*wPanels)[1]->getLower() == (*bPanels)[1], "Lower Parent Panel not set properly")
//    
//    TEST_ASSERT_MSG((*wPanels)[3]->getUpper() == (*bPanels)[3], "Upper Parent Panel not set properly at wing-body joint")
//    
//    TEST_ASSERT_MSG((*wPanels)[3]->getLower() == (*bPanels)[4], "Lower Parent Panel not set properly at wing-body joint")
//}
//
//void GeomTests::test_createOctree()
//{
//    TEST_ASSERT(testGeom->getOctree()->getMembers().size() == nTris);
//}
//
//void GeomTests::makeTestTriFile()
//{
//    // Construct tri file to test
//    std::ofstream fid;
//    testTriFile = "testfile.tri";
//    fid.open(testTriFile);
//    if (fid)
//    {
//        // See Visual of testfile.tri in geometryTests.h
//        // Nodes with z=-1 (2,4,11) are for lower panel tris.
//        
//        nTris = 15;
//        nNodes = 19;
//        double eps = pow(10,-7); // Addresses floating point error in wake nodes intended to be coincident with body nodes seen in .tri files output by OpenVSP.
//
//        fid << nNodes << "\t" << nTris << "\n";
//        
//        fid << 0 << "\t" << -3 << "\t" << 0 << "\n"
//        << 0 << "\t" << -3 << "\t" << -1 << "\n"
//        << 0 << "\t" << -2 << "\t" << 0 << "\n"
//        << 0 << "\t" << -2 << "\t" << -1 << "\n"
//        << 0 << "\t" << -1 << "\t" << 0 << "\n"
//        << 0 << "\t" << 0 << "\t" << 0 << "\n"
//        << 1 << "\t" << -3 << "\t" << 0 << "\n"
//        << 1 << "\t" << -2 << "\t" << 0 << "\n"
//        << 1 << "\t" << -1 << "\t" << 0 << "\n"
//        << 1 << "\t" << 0 << "\t" << 0 << "\n"
//        << 1 << "\t" << 0 << "\t" << -1 << "\n"
//        << 2 << "\t" << -1 << "\t" << 0 << "\n"
//        << 2 << "\t" << 0 << "\t" << 0 << "\n"
//        << 1 << "\t" << -3 << "\t" << eps << "\n"
//        << 1 << "\t" << -2 << "\t" << eps << "\n"
//        << 1 << "\t" << -1 << "\t" << eps << "\n"
//        << 2 << "\t" << -3 << "\t" << 0 << "\n"
//        << 2 << "\t" << -2 << "\t" << 0 << "\n"
//        << 2 << "\t" << -1 << "\t" << 0 << "\n";
//        
//        fid << 1 << "\t" << 7 << "\t" << 8 << "\n"
//        << 2 << "\t" << 8 << "\t" << 7 << "\n"
//        << 1 << "\t" << 8 << "\t" << 3 << "\n"
//        << 3 << "\t" << 8 << "\t" << 9 << "\n"
//        << 4 << "\t" << 9 << "\t" << 8 << "\n"
//        << 3 << "\t" << 9 << "\t" << 5 << "\n"
//        << 5 << "\t" << 9 << "\t" << 6 << "\n"
//        << 6 << "\t" << 9 << "\t" << 10 << "\n"
//        << 10 << "\t" << 9 << "\t" << 12 << "\n"
//        << 11 << "\t" << 12 << "\t" << 9 << "\n"
//        << 10 << "\t" << 12 << "\t" << 13 << "\n"
//        << 14 << "\t" << 17 << "\t" << 18 << "\n"
//        << 14 << "\t" << 18 << "\t" << 15 << "\n"
//        << 15 << "\t" << 18 << "\t" << 19 << "\n"
//        << 15 << "\t" << 19 << "\t" << 16 << "\n";
//        
//        fid << 1 << "\n" << 1 << "\n" << 1 << "\n" << 1 << "\n" << 1 << "\n" << 1 << "\n" << 2 << "\n" << 2 << "\n" << 2 << "\n" << 2 << "\n" << 2 << "\n" << 10001 << "\n" << 10001 << "\n" << 10001 << "\n" << 10001;
//
//        fid.close();
//    }
//}
