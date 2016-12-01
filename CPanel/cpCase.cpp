//
//  runCase.cpp
//  CPanel
//
//  Created by Chris Satterwhite on 10/13/14.
//  Copyright (c) 2014 Chris Satterwhite. All rights reserved.
//

#include "cpCase.h"

cpCase::~cpCase()
{
    for (int i=0; i<bStreamlines.size(); i++)
    {
        delete bStreamlines[i];
    }
}

void cpCase::run(bool printFlag, bool surfStreamFlag, bool stabDerivFlag, bool vortPartFlag)
{
    if(unsteady){
        readBodyKinFile();
    }
    
    
//    int nNodes = geom->getNodes().size();
//    std::vector<cpNode*> nodes = geom->getNodes();
//    double maxPt = 0;
//    for (int i=0; i<nNodes; i++) {
//        if(nodes[i]->getPnt().x() > maxPt){
//            maxPt = nodes[i]->getPnt().x();
//        }
//    }
//    std::cout << "maxPtx = " << maxPt << std::endl;
    
    
    bool converged;
    std::string check = "\u2713";
    setSourceStrengths();
    
    converged = solveMatrixEq();
    

    
    if(vortPartFlag){
        
        converged = false;
        bool matrixConversion = false;
        
        if(unsteady){
//            compVelocity(); // This doesn't need to be here unless solving Cl/Cd etc
        }
        writeFiles(); // Might take this out to or at least make it an option?

        
        // Initially convect buffer wake
        timeStep++;
        convectBufferWake();
        converged = solveVPmatrixEq();
        
        if (unsteady) {
//            compVelocity();
        }
        writeFiles();

        
        for(int i=0; i<(numSteps-2); i++){

            timeStep++;
            
            std::cout << "Time step " << timeStep << "/" << numSteps << ".  Flow time = " << timeStep*dt;
            std::cout << ".  " << particles.size() << " particles" << std::endl;

            if(highAccuracy){
                collapseWakeForEachEdge(); // Need to modify this too...
            }else{
                collapseBufferWake();

            }
            
            convectBufferWake();
            
//            particleStrengthUpdate();
            particleStrengthUpdateGaussian();

            
            if(accelerate){
                
                partOctree.removeData();
                partOctree.setMaxMembers(10); //Barnes Hut
                
                if(highAccuracy){
                    partOctree.setMaxTheta(0.25);
                }
                partOctree.addData(particles);
                
                FMM.build(&partOctree);
            }
            
            
            setSourceStrengths();
            
            matrixConversion = solveVPmatrixEq();
            
            // Check for simulation convergence
            if (!matrixConversion) {
                std::cout << "*** Warning : Solution did not converge ***" << std::endl;
                std::exit(0);
            }else{
                bool conv = solutionConvergence();
                if (conv) {
//                    break;
                }
            }
            
            if(unsteady){
//                compVelocity();
            }
            
            writeFiles();

            convectParticles();
            
        }
        
        timeStep++; // Because the last timestep incrementor TAKE OUT WHEN RE-ARRANGE LOOP
    }
    
    
    if (printFlag)
    {
        std::cout << std::setw(17) << std::left << check << std::flush;
    }
    
    compVelocity();
    
    
    if (printFlag)
    {
        std::cout << std::setw(17) << std::left << check << std::flush;
    }
    
    trefftzPlaneAnalysis();
    
    if (printFlag)
    {
        std::cout << std::setw(18) << std::left << check << std::flush;
    }

    if (surfStreamFlag)
    {
        createStreamlines();
        if (printFlag)
        {
            std::cout << std::setw(16) << std::left << check << std::flush;
        }
    }
    else
    {
        if (printFlag)
        {
            std::cout << std::setw(16) << std::left << "X" << std::flush;
        }
    }
    
    if (stabDerivFlag)
    {
        stabilityDerivatives();
        if (printFlag)
        {
            std::cout << std::setw(23) << std::left << check << std::endl;
        }
    }
    else
    {
        if (printFlag)
        {
            std::cout << std::setw(23) << std::left << "X" << std::endl;
        }
    }
    
        
    if (!converged && printFlag)
    {
        std::cout << "*** Warning : Solution did not converge ***" << std::endl;
    }
    
    if (printFlag){
        writeFiles();
    }
    
}


Eigen::Vector3d cpCase::Vinf(Eigen::Vector3d POI)
{
    if (!unsteady)
    {
        return windToBody( Vmag , alpha , beta );
    }
    else
    {
        Eigen::Vector3d localVel;
        
        // U = U3 + (-q*z + r*y)
        localVel.x() = bodyKin(timeStep, 0) - bodyKin(timeStep, 4)*POI.z() + bodyKin(timeStep, 5)*POI.y();

        // V = V3 + (-r*x + p*z)
        localVel.y() = bodyKin(timeStep, 1) - bodyKin(timeStep, 5)*POI.x() + bodyKin(timeStep, 3)*POI.z();
        
        // W = W3 + (-p*y + q*x)
        localVel.z() = bodyKin(timeStep, 2) - bodyKin(timeStep, 3)*POI.y() + bodyKin(timeStep, 4)*POI.x();
    
        return localVel;
    }
}



Eigen::Vector3d cpCase::windToBody(double V, double alpha, double beta)
{
    alpha *= M_PI/180;
    beta *= M_PI/180;
    
    Eigen::Matrix3d T;
    T(0,0) = cos(alpha)*cos(beta);
    T(0,1) = cos(alpha)*sin(beta);
    T(0,2) = -sin(alpha);
    T(1,0) = -sin(beta);
    T(1,1) = cos(beta);
    T(1,2) = 0;
    T(2,0) = sin(alpha)*cos(beta);
    T(2,1) = sin(alpha)*sin(beta);
    T(2,2) = cos(alpha);
    Eigen::Vector3d Vt;
    Eigen::Vector3d Vel;
    Vt << V,0,0;
    transform = T;

    Vel = transform*Vt;
    
    return Vel;
}

Eigen::Vector3d cpCase::bodyToWind(const Eigen::Vector3d &vec)
{
    Eigen::Matrix3d T = transform.transpose();
    return T*vec;
}

void cpCase::setSourceStrengths()
{
    sigmas.resize(bPanels->size());
    for (int i=0; i<bPanels->size(); i++)
    {
        Eigen::Vector3d sumVelInfl = Eigen::Vector3d::Zero();
        
        if (accelerate && timeStep > 0){
            sumVelInfl += FMM.barnesHutVel((*bPanels)[i]->getCenter());
        }
        else
        {
            for(int j=0; j<particles.size(); j++)
            {
                sumVelInfl += particles[j]->partVelInflGaussian((*bPanels)[i]->getCenter());
            }
        }
        
        for(int j=0; j<filaments.size(); j++)
        {
            sumVelInfl += filaments[j]->velInfl( (*bPanels)[i]->getCenter() );
        }
        
        (*bPanels)[i]->setSigma( Vinf( (*bPanels)[i]->getCenter() ) + sumVelInfl , 0 ); //function called VinfLocal(POI), and this func won't change
        sigmas(i) = (*bPanels)[i]->getSigma();
    }
}

bool cpCase::solveMatrixEq()
{
    bool converged = true;
    
    Eigen::MatrixXd* A = geom->getA();
    Eigen::MatrixXd* B = geom->getB();
    Eigen::VectorXd RHS = -(*B)*sigmas;
    Eigen::VectorXd doubletStrengths(bPanels->size());

    
    Eigen::BiCGSTAB<Eigen::MatrixXd> res;
    res.compute((*A));
    doubletStrengths = res.solve(RHS);
    if (res.error() > pow(10,-10))
    {
        converged = false;
    }
    
    for (int i=0; i<bPanels->size(); i++)
    {
        (*bPanels)[i]->setMu(doubletStrengths(i));
        (*bPanels)[i]->setPotential(Vinf((*bPanels)[i]->getCenter()));
    }
    
    for (int i=0; i<wPanels->size(); i++)
    {
        (*wPanels)[i]->setMu();
        (*wPanels)[i]->setPotential(Vinf((*wPanels)[i]->getCenter()));
    }
    
    return converged;
}

bool cpCase::solveVPmatrixEq()
{
    bool converged = true;
    
    Eigen::MatrixXd* A = geom->getA();
    Eigen::MatrixXd* B = geom->getB();
    Eigen::MatrixXd* C = geom->getC();
    Eigen::VectorXd RHS = -(*B)*(sigmas) - (*C)*(wake2Doublets);
    Eigen::VectorXd doubletStrengths(bPanels->size());
    
    
    Eigen::BiCGSTAB<Eigen::MatrixXd> res;
    res.compute((*A));
    doubletStrengths = res.solve(RHS);
    if (res.error() > pow(10,-10))
    {
        converged = false;
    }

    for (int i=0; i<bPanels->size(); i++)
    {
        (*bPanels)[i]->setMu(doubletStrengths(i));
        (*bPanels)[i]->setPotential(VinfPlusVecPot((*bPanels)[i]->getCenter()));

    }

    // VinfPlusVecPot
    for (int i=0; i<wPanels->size(); i++)
    {
        (*wPanels)[i]->setPrevStrength((*wPanels)[i]->getMu());
        (*wPanels)[i]->setMu();
        
        (*wPanels)[i]->setPotential(VinfPlusVecPot((*wPanels)[i]->getCenter()));
        
        (*w2panels)[i]->setPotential(VinfPlusVecPot((*w2panels)[i]->getCenter())); // Included in this loop because there are the same number of w2pans as w1pans

    }
    return converged;
}

bool cpCase::solutionConvergence(){
    
    if(unsteady | (timeStep < 5)){
        return false;
    }
    
    // Check for simulation convergence. Will use the Trefftz plane because it is so cheap. Medium size mesh takes <1% of time step time.
    
    double strengthDiff = 0;
    for(int i=0; i<(*w2panels).size(); i++)
    {
        double pan2mu = (*w2panels)[i]->getMu();
        double pan1mu = (*w2panels)[i]->getBufferParent()->getMu();
        
        strengthDiff += std::abs(pan2mu - pan1mu);
    }
    
    strengthDiff /= (*w2panels).size();
    
    // Store previous time step values
    double CDnm1 = CD_trefftz;
    double CLnm1 = CL_trefftz;
    
    if(timeStep > 2){
        // Need particles to create the S curve
        trefftzPlaneAnalysis();
    }
    
    std::cout << timeStep << " " << CL_trefftz << "  " << CD_trefftz << std::endl;
    
    double changeCD = std::abs((CD_trefftz - CDnm1)/CDnm1);
    double changeCL = std::abs((CL_trefftz - CLnm1)/CLnm1);
    
//    std::cout << timeStep<< "        " << CL_trefftz << "   " <<changeCL << "        " << CD_trefftz << "       " <<  changeCD  << std::endl;
    
    double val = 0.0001;
    
    if((changeCL < val) && (changeCD < val) && (changeCDnm1 < val) && (changeCDnm1 < val)){
        return true;
    }
    
    return false;
    
}




Eigen::Vector3d cpCase::VinfPlusVecPot(Eigen::Vector3d POI)
{
    Eigen::Vector3d vInfluence = Vinf(POI);
    
    // Particles
    if(accelerate && timeStep > 1)
    {
        vInfluence += FMM.barnesHutVel(POI);
    }else
    {
        for (int i=0; i<particles.size(); i++) {
            vInfluence += particles[i]->partVelInflGaussian(POI);
        }
    }
    
    // Filaments
    for (int i=0; i<filaments.size(); i++)
    {
        vInfluence += filaments[i]->velInfl(POI);
    }
    
    return vInfluence;
}

void cpCase::compVelocity()
{
    //  Velocity Survey with known doublet and source strengths
    CM.setZero();
    Eigen::Vector3d moment;
    Fbody = Eigen::Vector3d::Zero();
    bodyPanel* p;
    
    std::vector<double> panLift;
    bool liftDist = false;
    
    for (int i=0; i<bPanels->size(); i++)
    {
        p = (*bPanels)[i];
        p->computeVelocity(PG,Vinf(p->getCenter()));
        if (unsteady) {
            p->computeCp( Vinf((*bPanels)[i]->getCenter()).norm() , dt ); // Katz 13.169
        }else{
            p->computeCp( Vmag );
        }
        Fbody += -p->getCp()*p->getArea()*p->getBezNormal()/params->Sref;
        moment = p->computeMoments(params->cg);
        CM(0) += moment(0)/(params->Sref*params->bref);
        CM(1) += moment(1)/(params->Sref*params->cref);
        CM(2) += moment(2)/(params->Sref*params->bref);
    
//        if(liftDist)
//        {
//            panLift.push_back( (-p->getCp()*p->getArea()*p->getBezNormal()/params->cref ).z() );
//        }
    
    }
    Fwind = bodyToWind(Fbody);
    CL.push_back(Fbody.z());
    CM_x.push_back(CM(1));

    
//    if(liftDist)
//    {
//        std::cout << "panY = [";
//        for (int i = 0; i<(*bPanels).size(); i++) {
//            std::cout << (*bPanels)[i]->getCenter().y() << ", ";
//        }
//        std::cout << "];" << std::endl;
//        
//        std::cout << "panLift = [";
//        for (int i = 0; i<(*bPanels).size(); i++) {
//            std::cout << panLift[i] << ", ";
//        }
//        std::cout << "];" << std::endl;
//    }
    
//    std::cout << "CL = [";
//    for (int i=0; i<CL.size(); i++) {
//        std::cout << CL[i] << ", ";
//    }
//    std::cout << "]" << std::endl;
//    
//    std::cout << "CM_pitch = [";
//    for (int i=0; i<CM_x.size(); i++) {
//        std::cout << CM_x[i] << ", ";
//    }
//    std::cout << "]" << std::endl;
    
}

void cpCase::trefftzPlaneAnalysis()
{
    std::vector<wake*> wakes = geom->getWakes();
    CL_trefftz = 0;
    CD_trefftz = 0;
    double CD_trefftzVel = 0;

    
    for (int i=0; i<wakes.size(); i++)
    {
        if(vortPartFlag){
            wakes[i]->trefftzPlaneVP(Vmag,params->Sref, &particles, timeStep);
        } else{
            wakes[i]->trefftzPlane(Vmag,params->Sref);
        }
        CL_trefftz += wakes[i]->getCL()/PG;
        CD_trefftz += wakes[i]->getCD()/pow(PG,2);
    }
    double CD_T = CD_trefftz;
}

void cpCase::createStreamlines()
{
    // Gather TE Panels
    std::vector<surface*> surfs = geom->getSurfaces();
    std::vector<std::pair<Eigen::Vector3d,bodyPanel*>> streamPnts;
    bodyStreamline* s;
    
    for (int i=0; i<surfs.size(); i++)
    {
        streamPnts = surfs[i]->getStreamlineStartPnts(Vinf({0,0,0}),PG); // will need to modify this in order to find the stagnation points for a moving body...
        for (int j=0; j<streamPnts.size(); j++)
        {
            s = new bodyStreamline(std::get<0>(streamPnts[j]),std::get<1>(streamPnts[j]),Vinf({0,0,0}),PG,geom,3,false);
            bStreamlines.push_back(s);
        }
    }
    
//    // Off Body Streamline (Temporary)
//    std::vector<streamline*> sLines;
//    streamline* sline;
//    Eigen::Vector3d start;
//    int ny = 1;
//    int nz = 6;
//    double dz = 0.1;
//    for (int i=0; i<ny; i++)
//    {
//        for (int j=0; j<nz; j++)
//        {
////            start << -1,-0.5*params->bref+params->bref*i/(ny-1),-0.35+j*dz;
//            start << -1,-2.85,-0.45+j*dz;
//            sline = new streamline(start,20,0.0001,Vinf,PG,geom);
//            sLines.push_back(sline);
//        }
//    }
//    
//    piece p;
//    std::vector<piece> pieces;
//    std::vector<pntDataArray> data;
//    pntDataArray vel("Velocity");
//    Eigen::MatrixXi con;
//    Eigen::MatrixXd pntMat;
//    std::vector<Eigen::Vector3d> pnts,velocities;
//    
//    for (int i=0; i<sLines.size(); i++)
//    {
//        pnts = sLines[i]->getPnts();
//        velocities = sLines[i]->getVelocities();
//        vel.data.resize(velocities.size(),3);
//        pntMat.resize(pnts.size(),3);
//        con.resize(pnts.size()-1,2);
//        for (int j=0; j<pnts.size(); j++)
//        {
//            pntMat.row(j) = pnts[j];
//            vel.data.row(j) = velocities[j];
//            if (j<con.rows())
//            {
//                con(j,0) = j;
//                con(j,1) = j+1;
//            }
//        }
//        data.push_back(vel);
//        p.pnts = pntMat;
//        p.connectivity = con;
//        p.pntData = data;
//        
//        pieces.push_back(p);
//        data.clear();
//    }
//    
//    boost::filesystem::path path = boost::filesystem::current_path();
//    
//    std::string fname = path.string()+"/OffBodyStreamlines.vtu";
//    VTUfile wakeFile(fname,pieces);
//    
//    for (int i=0; i<sLines.size(); i++)
//    {
//        delete sLines[i];
//    }
}

void cpCase::stabilityDerivatives()
{
    double delta = 0.5;
    double dRad = delta*M_PI/180;
    cpCase dA(geom,Vmag,alpha+delta,beta,mach,params);
    cpCase dB(geom,Vmag,alpha,beta+delta,mach,params);
    
    dA.run(false,false,false,false);
    dB.run(false,false,false,false);
    
    Eigen::Vector3d FA = dA.getWindForces();
    FA(2) = dA.getCL();
    FA(0) = dA.getCD();
    
    Eigen::Vector3d FB = dB.getWindForces();
    FB(2) = dB.getCL();
    FB(0) = dB.getCD();
    
    Eigen::Vector3d F = Fwind;
    F(2) = CL_trefftz;
    F(0) = CD_trefftz;
    
    // Finish calculating stability derivatives
    dF_dAlpha = (FA-F)/dRad;
    dF_dBeta = (FB-F)/dRad;
    
    dM_dAlpha = (dA.getMoment()-CM)/dRad;
    dM_dBeta = (dB.getMoment()-CM)/dRad;
    
    
}

void cpCase::writeFiles()
{
    std::stringstream caseLabel;
    caseLabel << "/V" << Vmag << "_Mach" << mach << "_alpha" << alpha << "_beta" << beta;
    boost::filesystem::path subdir = boost::filesystem::current_path().string()+caseLabel.str();
    if (!boost::filesystem::exists(subdir))
    {
        boost::filesystem::create_directories(subdir);
    }
    Eigen::MatrixXd nodeMat = geom->getNodePnts();
    writeBodyData(subdir,nodeMat);
    if (geom->getWakes().size() > 0)
    {
        writeWakeData(subdir,nodeMat);
        writeBuffWake2Data(subdir,nodeMat);
        writeSpanwiseData(subdir);
    }
    
    if (false) {
        createVolMesh();
        writeVolMeshData(subdir, pts, cells);
    }
    
    if(vortPartFlag){
//        if(timeStep > 0)
        {
            writeParticleData(subdir);
            writeFilamentData(subdir);
        }
    }
    
    if (params->surfStreamFlag)
    {
        writeBodyStreamlines(subdir);
    }

}

void cpCase::writeBodyData(boost::filesystem::path path,const Eigen::MatrixXd &nodeMat)
{
    std::vector<cellDataArray> data;
    cellDataArray mu("Doublet Strengths"),sigma("Source Strengths"),pot("Velocity Potential"),V("Velocity"),Cp("Cp"),bN("bezNormals"),x("xPosition"),y("yPosition"),z("zPostition");
    Eigen::MatrixXi con(bPanels->size(),3);
    mu.data.resize(bPanels->size(),1);
    sigma.data.resize(bPanels->size(),1);
    pot.data.resize(bPanels->size(),1);
    V.data.resize(bPanels->size(),3);
    Cp.data.resize(bPanels->size(),1);
    bN.data.resize(bPanels->size(),3);
    x.data.resize(bPanels->size(),1);
    y.data.resize(bPanels->size(),1);
    z.data.resize(bPanels->size(),1);
    
    for (int i=0; i<bPanels->size(); i++)
    {
        mu.data(i,0) = (*bPanels)[i]->getMu();
        sigma.data(i,0) = (*bPanels)[i]->getSigma();
        pot.data(i,0) = (*bPanels)[i]->getPotential();
        V.data.row(i) = (*bPanels)[i]->getGlobalV();
        Cp.data(i,0) = (*bPanels)[i]->getCp();
        con.row(i) = (*bPanels)[i]->getVerts();
        bN.data.row(i) = (*bPanels)[i]->getBezNormal();
        x.data(i,0) = (*bPanels)[i]->getCenter().x();
        y.data(i,0) = (*bPanels)[i]->getCenter().y();
        z.data(i,0) = (*bPanels)[i]->getCenter().z();
    }
    
    data.push_back(mu);
    data.push_back(sigma);
    data.push_back(pot);
    data.push_back(V);
    data.push_back(Cp);
    data.push_back(bN);
    data.push_back(x);
    data.push_back(y);
    data.push_back(z);
    
    piece body;
    body.pnts = nodeMat;
    body.connectivity = con;
    body.cellData = data;
    
    std::string fname = path.string()+"/surfaceData-" + std::to_string(timeStep)+".vtu";
    VTUfile bodyFile(fname,body);
}

void cpCase::writeWakeData(boost::filesystem::path path, const Eigen::MatrixXd &nodeMat)
{
    std::vector<cellDataArray> data;
    cellDataArray mu("Doublet Strengths"),pot("Velocity Potential");
    Eigen::MatrixXi con;
    if(vortPartFlag){ //VPP
        con.resize(wPanels->size(),4);
    }else{
        con.resize(wPanels->size(),3);
    }
    mu.data.resize(wPanels->size(),1);
    pot.data.resize(wPanels->size(),1);
    for (int i=0; i<wPanels->size(); i++)
    {
        mu.data(i,0) = (*wPanels)[i]->getMu();
        pot.data(i,0) = (*wPanels)[i]->getPotential();
        con.row(i) = (*wPanels)[i]->getVerts();
    }
    data.push_back(mu);
    data.push_back(pot);
    
    piece wake;
    wake.pnts = nodeMat;
    wake.connectivity = con;
    wake.cellData = data;
    std::string fname = path.string()+"/wakeData-"+std::to_string(timeStep)+".vtu";
    VTUfile wakeFile(fname,wake);
}

void cpCase::writeBuffWake2Data(boost::filesystem::path path, const Eigen::MatrixXd &nodeMat)
{
    std::vector<cellDataArray> data;
    cellDataArray mu("Doublet Strengths"),pot("Velocity Potential");
    Eigen::MatrixXi con;
    con.resize(w2panels->size(),4);
    mu.data.resize(w2panels->size(),1);
    pot.data.resize(w2panels->size(),1);
    for (int i=0; i<w2panels->size(); i++)
    {
        mu.data(i,0) = (*w2panels)[i]->getMu();
        pot.data(i,0) = (*w2panels)[i]->getPotential();
        con.row(i) = (*w2panels)[i]->getVerts();
    }
    data.push_back(mu);
    data.push_back(pot);
    
    piece wake;
    wake.pnts = nodeMat;
    wake.connectivity = con;
    wake.cellData = data;
    std::string fname = path.string()+"/bufferWake2Data-"+std::to_string(timeStep)+".vtu";
    VTUfile wakeFile(fname,wake);
}

void cpCase::writeFilamentData(boost::filesystem::path path)
{
    std::vector<cellDataArray> data;
    cellDataArray mu("Gamma");
    Eigen::MatrixXi con;
    
    Eigen::MatrixXd nodeMat(2*filaments.size(),3);
    for(int i=0; i<filaments.size(); i++){
        nodeMat.row(2*i) = filaments[i]->getP1();
        nodeMat.row(2*i+1) = filaments[i]->getP2();
    }
    

    con.resize(filaments.size(),2); // 2 pnts per filament
    
    mu.data.resize(filaments.size(),1);

    int nodeCounter = 0; // Used to make filament end points
    for (int i=0; i<filaments.size(); i++)
    {
        mu.data(i,0) = filaments[i]->getStrength();
        con(i,0) = nodeCounter;
        con(i,1) = nodeCounter+1;
        nodeCounter = nodeCounter+2;
    }
    
    data.push_back(mu);
    
    piece fils;
    fils.pnts = nodeMat;
    fils.connectivity = con;
    fils.cellData = data;
    std::string fname = path.string()+"/filaments-"+std::to_string(timeStep)+".vtu";
    VTUfile filFile(fname,fils);
}

void cpCase::writeParticleData(boost::filesystem::path path)
{

    Eigen::MatrixXd partMat(particles.size(),3);
    for (int i=0; i<particles.size(); i++)
    {
        partMat.row(i) = particles[i]->pos;
    }
    
    std::vector<cellDataArray> data;
    cellDataArray strength("Strength"), shedTime("Time Step Shed");
    Eigen::MatrixXi con(particles.size(),1);
    Eigen::MatrixXi shed(particles.size(),1);
    
    strength.data.resize(particles.size(),3);
    shedTime.data.resize(particles.size(),1);
    for (int i=0; i<particles.size(); i++)
    {
        strength.data.row(i) = particles[i]->strength;
        shedTime.data(i,0) = particles[i]->shedTime;
        con(i) = i;
    }
    data.push_back(strength);
    data.push_back(shedTime);
    
    piece parts;
    parts.pnts = partMat;
    parts.connectivity = con;
    parts.cellData = data;
    
    std::string fname = path.string()+"/particleData-"+std::to_string(timeStep)+".vtu";
    VTUfile partFile(fname,parts);
}


void cpCase::writeSpanwiseData(boost::filesystem::path path)
{
    std::vector<wake*> wakes = geom->getWakes();
    for (int i=0; i<wakes.size(); i++)
    {
        Eigen::VectorXd spanLoc,Cl,Cd;
        spanLoc = 2*wakes[i]->getSpanwisePnts()/params->bref;
        Cl = wakes[i]->getSpanwiseCl()/PG;
        Cd = wakes[i]->getSpanwiseCd()/pow(PG,2);
        
        std::stringstream ss;
        ss << path.string() << "/spanwiseData_Wake" << i+1 << ".csv";
        std::string fname = ss.str();
        std::ofstream fout;
        fout.open(fname);
        if (fout)
        {
            fout << "2y/b,Cl,Cdi" << std::endl;
            for (int i=0; i<spanLoc.size(); i++)
            {
                fout << spanLoc(i) << "," << Cl(i) << "," << Cd(i) << std::endl;
            }
        }
        fout.close();
    }
}

void cpCase::writeBodyStreamlines(boost::filesystem::path path)
{
    piece p;
    std::vector<piece> pieces;
    std::vector<pntDataArray> data;
    pntDataArray vel("Velocity");
    Eigen::MatrixXi con;
    Eigen::MatrixXd pntMat;
    std::vector<Eigen::Vector3d> pnts,velocities;
    
    for (int i=0; i<bStreamlines.size(); i++)
    {
        pnts = bStreamlines[i]->getPnts();
        velocities = bStreamlines[i]->getVelocities();
        vel.data.resize(velocities.size(),3);
        pntMat.resize(pnts.size(),3);
        con.resize(pnts.size()-1,2);
        for (int j=0; j<pnts.size(); j++)
        {
            pntMat.row(j) = pnts[j];
            vel.data.row(j) = velocities[j];
            if (j<con.rows())
            {
                con(j,0) = j;
                con(j,1) = j+1;
            }
        }
        data.push_back(vel);
        p.pnts = pntMat;
        p.connectivity = con;
        p.pntData = data;
        
        pieces.push_back(p);
        data.clear();
    }
    
    std::string fname = path.string()+"/streamlines.vtu";
    VTUfile wakeFile(fname,pieces);
}


void cpCase::writeVolMeshData(boost::filesystem::path path, Eigen::MatrixXd &nodeMat, std::vector<Eigen::VectorXi> cells)
{
    int nCells = (int)cells.size();
    
    std::vector<cellDataArray> data;
    cellDataArray vel("Velocity"),vortMag("Vorticity Magnitude"), Cp("Cp");
    Eigen::MatrixXi con;
    con.resize(nCells,8);
    vel.data.resize(nCells,3);
    vortMag.data.resize(nCells,1);
    Cp.data.resize(nCells,1);
    for (int i=0; i<nCells; i++)
    {
        vel.data.row(i) = volMeshDat.velocity[i];
        vortMag.data(i,0) = volMeshDat.vorticity[i];
        Cp.data(i,0) = volMeshDat.pressure[i];
        con.row(i) = cells[i];
    }
    data.push_back(vel);
    data.push_back(vortMag);
    data.push_back(Cp);
    
    piece volMesh;
    volMesh.pnts = nodeMat;
    volMesh.connectivity = con;
    volMesh.cellData = data;
    std::string fname = path.string()+"/volumeMesh-"+std::to_string(timeStep)+".vtu";
    VTUfile wakeFile(fname,volMesh);
}

void cpCase::convectBufferWake()
{
    wake2Doublets.resize((*w2panels).size());
    for (int i=0; i<(*w2panels).size(); i++)
    {
        (*w2panels)[i]->setPrevStrength((*w2panels)[i]->getMu());
        double parentMu = (*w2panels)[i]->getBufferParent()->getMu();
        (*w2panels)[i]->setMu(parentMu);
        wake2Doublets[i] = parentMu;
    }
    
}



void cpCase::collapseBufferWake(){
    
    std::vector<edge*> usedEdges;
    
    for(int i=0; i<(*w2panels).size(); i++)
    {
        std::vector<edge*> pEdges = (*w2panels)[i]->edgesInOrder();
        
        Eigen::Vector3d strength = Eigen::Vector3d::Zero();
        for (int j=1; j<4; j++)
        {
            if (!edgeIsUsed(pEdges[j],usedEdges))
            {
                usedEdges.push_back(pEdges[j]);
                strength += edgeStrength((*w2panels)[i], pEdges[j], j); // Don't need to pass in pEdges...
            }
        }
        
        Eigen::Vector3d pos = rungeKuttaStepper((*w2panels)[i]->getCenter());
        Eigen::Vector3d ptVel = Vinf(pos);
        double radius = (*w2panels)[i]->getPartRadius(ptVel,dt); // VinfLocal
        
        particle* p = new particle(pos, strength, radius, {0,0,0}, {0,0,0}, timeStep); // Last two are previous pos and strength values used for advanced time stepper.
        p->parentPanel = (*w2panels)[i];
        particles.push_back(p);
        
    }
    
    // Create filament
    if(filaments.size() == 0)
    {
        for(int i=0; i<(*w2panels).size(); i++)
        {
            vortexFil* fil;
            Eigen::Vector3d p1,p2;

            p1 = (*w2panels)[i]->pointsInOrder()[2]->getPnt();
            p2 = (*w2panels)[i]->pointsInOrder()[3]->getPnt();
            
            fil = new vortexFil(p1, p2,-(*w2panels)[i]->getMu(), (*w2panels)[i]); // Negative strength is because filament is actually the upstream edge being convected which is oriented the opposite direction as downstream edge

            filaments.push_back(fil);
            (*w2panels)[i]->setVortFil(fil);
        }
    }
    else
    {
        for(int i=0; i<(*w2panels).size(); i++){
            filaments[i]->setStrength(-(*w2panels)[i]->getMu()); // Negative strength is because filament is actually the upstream edge being convected which is oriented the opposite direction as downstream edge
        }
        
    }
}

void cpCase::collapseWakeForEachEdge(){
    // Go through each panel. Ignore the trialing edge as it has no circulation. After creating a particle at each edge, put a pointer to the edge in a vector. Before making any particles, make sure the edge has not yet been used. Edge collapse follows
    
    std::vector<edge*> usedEdges;
    for(int i=0; i<(*w2panels).size(); i++)
    {
        std::vector<edge*> pEdges = (*w2panels)[i]->edgesInOrder();
        for (int j=1; j<4; j++)
        {
            if (!edgeIsUsed(pEdges[j],usedEdges))
            {
                usedEdges.push_back(pEdges[j]);
                Eigen::Vector3d pos = seedPos((*w2panels)[i], j);
                
                Eigen::Vector3d strength = edgeStrength((*w2panels)[i], pEdges[j], j);
                double radius = (*w2panels)[i]->getPartRadius(Vinf((*w2panels)[i]->getCenter()),dt); // VinfLocal
                
                particle* p = new particle(pos, strength, radius, {0,0,0}, {0,0,0}, timeStep); //last two are previous pos and strength values used for advanced time stepper.
                particles.push_back(p);
            }
        }
    }

    // Create filament
    if(filaments.size() == 0 )
    {
        for(int i=0; i<(*w2panels).size(); i++)
        {
            vortexFil* fil;
            Eigen::Vector3d p1,p2;
            
            p1 = (*w2panels)[i]->pointsInOrder()[2]->getPnt();
            p2 = (*w2panels)[i]->pointsInOrder()[3]->getPnt();
            
            fil = new vortexFil(p1, p2,-(*w2panels)[i]->getMu(), (*w2panels)[i]); // Negative strength is because filament is actually the upstream edge being convected which is oriented the opposite direction as downstream edge
            
            filaments.push_back(fil);
            (*w2panels)[i]->setVortFil(fil);
        }
    }
    else
    {
        for(int i=0; i<(*w2panels).size(); i++)
        {
            filaments[i]->setStrength(-(*w2panels)[i]->getMu());//-(*w2panels)[i]->getMu());
        }
        
    }
    
}

bool cpCase::edgeIsUsed(edge* thisEdge, std::vector<edge*> pEdges){
    
    for(int i=0; i<pEdges.size(); i++){
        if(thisEdge == pEdges[i]){
            return true;
        }
    }
    return false;
}

Eigen::Vector3d cpCase::edgeStrength(wakePanel* pan, edge* curEdge, int edgeNum){
    //wait, why can't this go in the panel class? Ithink it can
    
    
    Eigen::Vector3d strength;
    std::vector<cpNode*> ptsIO = pan->pointsInOrder();
    
    if(edgeNum == 0){
        std::cout << "Don't try to collapse the upstream edge" << std::endl;
        std::exit(0);
        // Equation for part strength comes from Martin eq. 5.25
//        Eigen::Vector3d Rj = pan->pointsInOrder()[0]->getPnt();
//        Eigen::Vector3d Ri = pan->pointsInOrder()[1]->getPnt();
//        strength = (pan->getMu())*(Ri-Rj);
    }
    if(edgeNum == 2)
    {
        Eigen::Vector3d Rj = pan->pointsInOrder()[2]->getPnt();
        Eigen::Vector3d Ri = pan->pointsInOrder()[3]->getPnt();
        strength = (pan->getMu()-pan->getPrevStrength())*(Ri-Rj);
        
    }
    else if(edgeNum == 1)
    {
        wakePanel* otherPan = curEdge->getOtherWakePan(pan);
        Eigen::Vector3d Rj = ptsIO[1]->getPnt();
        Eigen::Vector3d Ri = ptsIO[2]->getPnt();
        
        if(otherPan) // Panel has neighbor
        {
            strength = (pan->getMu()-otherPan->getMu())*(Ri-Rj);
        }else{
            strength = pan->getMu()*(Ri-Rj);
        }
    }
    else // Is edge 3
    {
        wakePanel* otherPan = curEdge->getOtherWakePan(pan);
        Eigen::Vector3d Rj = ptsIO[3]->getPnt();
        Eigen::Vector3d Ri = ptsIO[0]->getPnt();
        
        if(otherPan)
        {
            strength = (pan->getMu()-otherPan->getMu())*(Ri-Rj);
        }else{
            strength = pan->getMu()*(Ri-Rj);
        }
    }
    return strength;
}

Eigen::Vector3d cpCase::rungeKuttaStepper(Eigen::Vector3d POI){
    
    // RK4 algorithm makes symmeterized influences VERY impractical
    
    Eigen::Vector3d k1, k2, k3, k4;
    
    k1 = velocityInflFromEverything(POI);
    k2 = velocityInflFromEverything(POI+k1*dt/2);
    k3 = velocityInflFromEverything(POI+k2*dt/2);
    k4 = velocityInflFromEverything(POI+k3*dt);
    
    return POI + dt*(k1/6 + k2/3 + k3/3 + k4/6);

}

Eigen::Vector3d cpCase::seedPos(wakePanel* pan, int edgeNum){
    std::vector<edge*> edgesIO = pan->edgesInOrder();
    
    Eigen::Vector3d partStart = edgesIO[edgeNum]->getMidPoint();

    return rungeKuttaStepper(partStart);
}

Eigen::Vector3d cpCase::velocityInflFromEverything(Eigen::Vector3d POI){
    
    // Freestream influence
    Eigen::Vector3d velOnPart = Vinf(POI);
//    std::cout<< "\n\n" << velOnPart << "\n\n\n" << std::endl;

    // Particle influence
    if(accelerate && timeStep > 2){
        velOnPart += FMM.barnesHutVel(POI);
    }
    else{
        for(int j=0;j<particles.size();j++)
        {
            velOnPart += particles[j]->partVelInfl(POI);
        }
    }

//    std::cout<< "\n\n" << velOnPart << "\n\n\n" << std::endl;
    // Body panel influence
    for(int j=0;j<(*bPanels).size();j++){
        velOnPart += (*bPanels)[j]->panelV(POI);
    }
//    std::cout<< "\n\n" << velOnPart << "\n\n\n" << std::endl;

    // Buffer wake influence
    for(int j=0;j<(*wPanels).size();j++){
        velOnPart += (*wPanels)[j]->panelV(POI);
    }
//    std::cout<< "\n\n" << velOnPart << "\n\n\n" << std::endl;

    for(int j=0;j<(*w2panels).size();j++){
        velOnPart += (*w2panels)[j]->panelV(POI);
    }

//    std::cout<< "\n\n" << velOnPart << "\n\n\n" << std::endl;

    // Vortex Filament influence
    for(int i=0; i<filaments.size(); i++){
        velOnPart +=filaments[i]->velInfl(POI);
    }
//    std::cout<< "\n\n" << velOnPart << "\n\n\n" << std::endl;

    return velOnPart;
}

Eigen::Vector3d cpCase::velocityInflFromEverything(particle* part)
{
    Eigen::Vector3d pos = part->pos;
    
    // Freestream influence
    Eigen::Vector3d velOnPart = Vinf(part->pos);
    
    // If particle reaches very far field condition
    if ((pos-Eigen::Vector3d::Zero()).norm() > 20*params->cref)
    {
        return velOnPart;
    }
    
    // Particle influence
    if(accelerate && timeStep > 2){
        velOnPart += FMM.barnesHutVel(part);
    }
    else{
        for(int j=0;j<particles.size();j++)
        {
            velOnPart += particles[j]->partVelInfl(pos);
        }
    }
    
    // If particle reaches first farfield condition
    if((pos-Eigen::Vector3d::Zero()).norm() > 12*params->cref)
    {
        return velOnPart;
    }
    
    // Body panel influence
    for(int j=0;j<(*bPanels).size();j++){
        velOnPart += (*bPanels)[j]->panelV(pos);
    }
    
    // Buffer wake influence
    for(int j=0;j<(*wPanels).size();j++){
        velOnPart += (*wPanels)[j]->panelV(pos);
    }
    for(int j=0;j<(*w2panels).size();j++){
        velOnPart += (*w2panels)[j]->panelV(pos);
    }
    
    // Vortex Filament influence
    for(int i=0; i<filaments.size(); i++){
        velOnPart +=filaments[i]->velInfl(pos);
    }
    
    return velOnPart;
}

void cpCase::convectParticles(){
    std::vector<Eigen::Vector3d> newPartPositions;
    
    for(int i=0;i<particles.size();i++){
        
        Eigen::Vector3d newPos;
        
        if(highAccuracy)
        {
            newPartPositions.push_back(rungeKuttaStepper(particles[i]->pos));

        }
        else
        {
            // Adams bashforth scheme
            Eigen::Vector3d downwash = Eigen::Vector3d::Zero(); // Used for rotor hover estimation
//            downwash.z() = -0.1137;
            Eigen::Vector3d velOnPart = velocityInflFromEverything(particles[i]) + downwash;
            
            
            if(particles[i]->getPrevVelInfl().isZero())
            {
                newPos = particles[i]->pos + dt*velOnPart;
            }
            else{
                newPos = particles[i]->pos + dt*(1.5*velOnPart - 0.5*particles[i]->getPrevVelInfl());
            }

            particles[i]->setPrevVelInfl(velOnPart);
            newPartPositions.push_back(newPos);
        }
    }
    
    for(int i=0;i<particles.size();i++){
        particles[i]->setPos(newPartPositions[i]);

    }
}



void cpCase::particleStrengthUpdate(){
    // This function uses the combined vortex stretching and diffusion equation used by Wincklemans for a regularized vortex core with high algebraic smoothing. (he refers to it as the strength update equation.)
    
    std::vector<Eigen::Vector3d> stretchDiffVec; // Creating a vector of diffusion values because the strength change needs to be set after all particle influences have been calculated
    for(int i=0; i<particles.size(); i++)
    {
        Eigen::Vector3d dAlpha = Eigen::Vector3d::Zero();
        
        // FarField2 condition
        if((particles[i]->pos-Eigen::Vector3d::Zero()).norm() > 20*params->cref)
        {
            stretchDiffVec.push_back(dAlpha);
        }
        // FarField1 condition
        else if((particles[i]->pos-Eigen::Vector3d::Zero()).norm() > 12*params->cref)
        {
            // Stretching from particles
            for(int j=0; j<particles.size(); j++)
            {
                dAlpha += particles[i]->partStrengthUpdate(particles[j]);
            }
            stretchDiffVec.push_back(dAlpha);
        }
        // Treat normally
        else{
            
            // Stretching from particles
            for(int j=0; j<particles.size(); j++)
            {
                dAlpha += particles[i]->partStrengthUpdate(particles[j]);
            }
            
            // Stretching from body panels
            for(int j=0; j<(*bPanels).size(); j++)
            {
                dAlpha += (*bPanels)[j]->partStretching(particles[i]);
            }
            
            // Stretching from wake panels
            for(int j=0;j<(*wPanels).size(); j++)
            {
                dAlpha += (*wPanels)[j]->partStretching(particles[i]);
            }
            
            stretchDiffVec.push_back(dAlpha);
        }
    }
    
    // No need for Kutta accuracy
    for(int i=0;i<particles.size();i++){
        Eigen::Vector3d newStrength;
        if(particles[i]->getprevStrengthUpdate().isZero())
        {
            newStrength = particles[i]->strength + stretchDiffVec[i]*dt;
        }else{
            //Adams bashforth
            newStrength = particles[i]->strength + dt*(1.5*stretchDiffVec[i] - 0.5*particles[i]->getprevStrengthUpdate());
        }
        particles[i]->setprevStrengthUpdate(stretchDiffVec[i]);
        particles[i]->setStrength(newStrength);
    }
    
}

void cpCase::particleStrengthUpdateGaussian(){
    // This function uses the update equations found in 'Vortex Methods for DNS of...' by Plouhmhans. It uses the particle strength exchange for the viscous diffusion
    
    // Can use this strength exchange function and call it from the normal stretching function instead of the winklemans function?
    
    std::vector<Eigen::Vector3d> stretchDiffVec; // Creating vector values because the strength change needs to be set after all particle influences have been calculated
    for(int i=0; i<particles.size(); i++)
    {
        Eigen::Vector3d dAlpha_diff = Eigen::Vector3d::Zero();
        Eigen::Vector3d dAlpha_stretch = Eigen::Vector3d::Zero();
        
        // FarField2 condition
        if((particles[i]->pos-Eigen::Vector3d::Zero()).norm() > 20*params->cref)
        {
            stretchDiffVec.push_back(Eigen::Vector3d::Zero());
        }
        // FarField1 condition
        else if((particles[i]->pos-Eigen::Vector3d::Zero()).norm() > 12*params->cref)
        {
            if(unsteady)
            {
                dAlpha_stretch = FMM.barnesHutStretch(particles[i]);
                dAlpha_diff = FMM.barnesHutDiff(particles[i]);
            }
            else
            {
                // Stretching from particles
                for(int j=0; j<particles.size(); j++)
                {
                    if(i!=j) // Kroneger Delta Function
                    {
                        dAlpha_diff += particles[i]->viscousDiffusionGaussian(particles[j]);
                        dAlpha_stretch += particles[i]->vortexStretchingGaussian(particles[j]);
                    }
                }
            }
            stretchDiffVec.push_back( dAlpha_diff + dAlpha_stretch );
        }
        else // Treat normally
        {
            // Influence from Particles
            if(accelerate && timeStep > 2){
                dAlpha_stretch += FMM.barnesHutStretch(particles[i]);
                dAlpha_diff += FMM.barnesHutDiff(particles[i]);
            }else
            {
                for(int j=0; j<particles.size(); j++)
                {
                    if(i!=j){ // Kroneger Delta Function
                        dAlpha_diff += particles[i]->viscousDiffusionGaussian(particles[j]);
                        dAlpha_stretch += particles[i]->vortexStretchingGaussian(particles[j]);
                    }
                }
            }
            
            // Stretching from body panels
            for(int j=0; j<(*bPanels).size(); j++)
            {
                dAlpha_diff += (*bPanels)[j]->partStretching(particles[i]);
            }
            
            // Stretching from wake panels
            for(int j=0;j<(*wPanels).size(); j++)
            {
                dAlpha_diff += (*wPanels)[j]->partStretching(particles[i]);
            }
            
            stretchDiffVec.push_back( dAlpha_diff + dAlpha_stretch );
        }
    }
        for(int i=0;i<particles.size();i++)
        {
            Eigen::Vector3d newStrength;
            if(particles[i]->getprevStrengthUpdate().isZero())
            {
                newStrength = particles[i]->strength + stretchDiffVec[i]*dt;
            }else
            {   // Adams bashforth
                newStrength = particles[i]->strength + dt*(1.5*stretchDiffVec[i] - 0.5*particles[i]->getprevStrengthUpdate());
            }
            particles[i]->setprevStrengthUpdate(stretchDiffVec[i]);
            particles[i]->setStrength(newStrength);
        }
}



void cpCase::readBodyKinFile(){
    std::ifstream fid;
    fid.open(params->bodyKinFileLoc);

    bodyKin.resize(numSteps, 6); // U, V, W, p, q, r
    
    for (int i=0; i<numSteps; i++)
    {
        for(int j=0; j < 6; j++)
        {
            fid >> bodyKin(i,j);
        }
    }
        
    fid.close();
}
















void cpCase::createVolMesh(){
    
    // Clear Past Timestep Mesh
    cells.clear();
    volMeshDat.cellCenter.clear();
    volMeshDat.velocity.clear();
    volMeshDat.vorticity.clear();
    
    // limits of mesh
    double x0, xf, y0, yf, z0, zf;
    x0 = -0.5;
    xf = 6.5;
    y0 = -8;
    yf = 8;
    z0 = -0.5;
    zf = 1.25;
    
    // resolution
    int nX, nY, nZ;
    nX = 500;
    nY = 1;
    nZ = 105;
    int nCells = nX * nY * nZ;
    
    double hx, hy, hz;
    hx = (xf-x0)/nX;
    hy = (yf-y0)/nY;
    hz = (zf-z0)/nZ;
    
    // Add one for number of points
    int nXp = nX+1;
    int nYp = nY+1;
    int nZp = nZ+1;
    int numPts = nXp * nYp * nZp;
    
    
    
    // creating pnt vector
    pts.resize(numPts, 3);
    int count = 0;
    for (int k=0; k<nZp; k++) {
        for (int j=0; j<nYp; j++) {
            for (int i=0; i<nXp; i++) {
                pts(count,0) = x0 + hx*i;
                pts(count,1) = y0 + hy*j;
                pts(count,2) = z0 + hz*k;
                count++;
            }
        }
    }
    
    
    // Create cells
    for (int i=0; i<nX; i++) {
        for (int j=0; j<nY; j++) {
            for (int k=0; k<nZ; k++) {
                Eigen::VectorXi cell_pts(8);
                cell_pts[0] = (k*(nXp*nYp) + j*(nXp) + i);
                cell_pts[1] = (k*(nXp*nYp) + j*(nXp) + i) + 1;
                
                cell_pts[2] = (k*(nXp*nYp) + (j+1)*(nXp) + i);
                cell_pts[3] = (k*(nXp*nYp) + (j+1)*(nXp) + i) + 1;
                
                cell_pts[4] = ((k+1)*(nXp*nYp) + j*(nXp) + i);
                cell_pts[5] = ((k+1)*(nXp*nYp) + j*(nXp) + i) + 1;
                
                cell_pts[6] = ((k+1)*(nXp*nYp) + (j+1)*(nXp) + i);
                cell_pts[7] = ((k+1)*(nXp*nYp) + (j+1)*(nXp) + i) + 1;
                
                cells.push_back(cell_pts);
                
                // Find cell center
                Eigen::MatrixXd cellCorners(8,3);
                for (int m=0; m<cellCorners.rows(); m++) {
                    cellCorners.row(m) = pts.row(cell_pts[m]);
                }
                
                Eigen::Vector3d center;
                for (int m=0; m<cellCorners.cols(); m++) {
                    center(m) = cellCorners.col(m).mean();
                }
                volMeshDat.cellCenter.push_back(center);
            }
        }
    }
    
    for (int i=0; i<nCells; i++) {
        // Find cell center
        Eigen::Vector3d velInCell = velocityInflFromEverything(volMeshDat.cellCenter[i]);
        volMeshDat.velocity.push_back(velInCell);
    }
    
    for (int i=0; i<nCells; i++) {
        // Incompressible Bernoulli equation
        double Cp = pow( volMeshDat.velocity[i].norm() / Vmag , 2 );
        volMeshDat.pressure.push_back(Cp);
    }
    
    for (int i=0; i<nCells; i++) {
        volMeshDat.vorticity.push_back(volMeshDat.cellCenter[i].norm());
    }
    
    
}




