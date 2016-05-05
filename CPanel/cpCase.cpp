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
    bool converged;
    std::string check = "\u2713";
    setSourceStrengths();
    converged = solveMatrixEq();

    
    
    if(vortPartFlag){
        std::cout << "Writing timestep " << timeStep << " files..." << std::endl;
        
        compVelocity();
        writeFiles();
        timeStep++;
        
        findSeedPoints(); //find seed points and radii once so don't have to every time step
        findSeedRadii();
        
//        std::cout << "time step = " << dt << std::endl;
        for(int i=0; i<numSteps; i++){
            std::cout << "Time step " << timeStep << ". Flow time = " << timeStep*dt << std::endl;
            std::cout << "Tracking " << particles.size() << " particles..." << std::endl;
            convectParticles();
//            vortexStretching();
            viscousDiffusion();
            collapseBufferWake();
            if(panStrengthChange.norm() < 0.02 && manualStepsSet == false && timeStep > 2){
                std::cout << "Convergence criteria met" << std::endl;
                break;
            }
            setVPSourceStrengths();
            std::cout << "Solving singularity strengths..." << std::endl;
            converged = solveVPmatrixEq();
//            if(i%5 == 0){
//                std::cout << "Computing velocity..." << std::endl;
//                compVelocity();
//            }
            writeFiles();
            timeStep++;
            std::cout << "part strength change norm = " << panStrengthChange.norm() << std::endl;
        }
        
        std::cout << "CL=[";
        for(int i=0; i<CL.size(); i++){
            std::cout << CL[i] << ", ";
        }
        std::cout << "];" << std::endl;
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
        std::cout << "*** Warning : Solution did not converge ***" << std::endl;;
    }
    
    if (printFlag){
        if(vortPartFlag == false){
            writeFiles();
        }
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
        (*bPanels)[i]->setSigma(Vinf,0);
        sigmas(i) = (*bPanels)[i]->getSigma();
    }
}

void cpCase::setVPSourceStrengths()
{
    sigmas.resize(bPanels->size());
    for (int i=0; i<bPanels->size(); i++)
    {
        (*bPanels)[i]->setVPSigma(Vinf,0,particles);
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
        (*bPanels)[i]->setPotential(Vinf);
    }
    
    for (int i=0; i<wPanels->size(); i++)
    {
        (*wPanels)[i]->setMu();
        (*wPanels)[i]->setPotential(Vinf);
    }
    
    if(vortPartFlag){
        panStrengthChange.resize(wPanels->size());
        for(int i=0; i<wPanels->size(); i++){
            prevPanStrength.push_back((*wPanels)[i]->getMu());
        }
    }

    
//    if(vortPartFlag){ //2BW
//        // Set second bufferwake strength equal to first
//        wake2Doublets.resize(wPanels->size());
//        
//        for (int i=0; i<wPanels->size(); i++){ //VPP
//            (*wake2panels)[i]->manuallySetMu((*wPanels)[i]->getMu());
//            (*wake2panels)[i]->setPotential(Vinf);
//            wake2Doublets(i) = (*wPanels)[i]->getMu();
//        }
//    }
    
    
    return converged;
}

bool cpCase::solveVPmatrixEq()
{
    bool converged = true;
    
    Eigen::MatrixXd* A = geom->getA();
    Eigen::MatrixXd* B = geom->getB();
    Eigen::VectorXd RHS = -(*B)*(sigmas);
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
        (*bPanels)[i]->setPotential(Vinf);
    }
    for (int i=0; i<wPanels->size(); i++)
    {
        prevPanStrength[i] = (*wPanels)[i]->getMu();
        (*wPanels)[i]->setMu();
        (*wPanels)[i]->setPotential(Vinf);
        panStrengthChange[i] = (*wPanels)[i]->getMu()-prevPanStrength[i];
    }
    
    return converged;
}


void cpCase::compVelocity()
{
    //  Velocity Survey with known doublet and source strengths
    CM.setZero();
    Eigen::Vector3d moment;
    Fbody = Eigen::Vector3d::Zero();
    bodyPanel* p;
    for (int i=0; i<bPanels->size(); i++)
    {
        p = (*bPanels)[i];
        p->computeVelocity(PG,Vinf);
        p->computeCp(Vmag);
        Fbody += -p->getCp()*p->getArea()*p->getBezNormal()/params->Sref;
        moment = p->computeMoments(params->cg);
        CM(0) += moment(0)/(params->Sref*params->bref);
        CM(1) += moment(1)/(params->Sref*params->cref);
        CM(2) += moment(2)/(params->Sref*params->bref);
    }
    Fwind = bodyToWind(Fbody);
    std::cout << "CL = " << Fbody.z() << std::endl;
    CL.push_back(Fbody.z());
}

void cpCase::trefftzPlaneAnalysis()
{
    std::vector<wake*> wakes = geom->getWakes();
    CL_trefftz = 0;
    CD_trefftz = 0;
    for (int i=0; i<wakes.size(); i++)
    {
        wakes[i]->trefftzPlane(Vmag,params->Sref);
        CL_trefftz += wakes[i]->getCL()/PG;
        CD_trefftz += wakes[i]->getCD()/pow(PG,2);
    }
}

void cpCase::createStreamlines()
{
    // Gather TE Panels
    std::vector<surface*> surfs = geom->getSurfaces();
    std::vector<std::pair<Eigen::Vector3d,bodyPanel*>> streamPnts;
    bodyStreamline* s;
    
    for (int i=0; i<surfs.size(); i++)
    {
        streamPnts = surfs[i]->getStreamlineStartPnts(Vinf,PG);
        for (int j=0; j<streamPnts.size(); j++)
        {
            s = new bodyStreamline(std::get<0>(streamPnts[j]),std::get<1>(streamPnts[j]),Vinf,PG,geom,3,false);
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
        writeSpanwiseData(subdir);
    }
    
    if(vortPartFlag){
        if(timeStep > 0){
//          writeBuffWake2Data(subdir,nodeMat);
            writeParticleData(subdir);
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
    cellDataArray mu("Doublet Strengths"),pot("Velocity Potential"),V("Velocity"),Cp("Cp"),bN("bezNormals");
    Eigen::MatrixXi con(bPanels->size(),3);
    mu.data.resize(bPanels->size(),1);
    pot.data.resize(bPanels->size(),1);
    V.data.resize(bPanels->size(),3);
    Cp.data.resize(bPanels->size(),1);
    bN.data.resize(bPanels->size(),3);
    for (int i=0; i<bPanels->size(); i++)
    {
        mu.data(i,0) = (*bPanels)[i]->getMu();
        pot.data(i,0) = (*bPanels)[i]->getPotential();
        V.data.row(i) = (*bPanels)[i]->getGlobalV();
        Cp.data(i,0) = (*bPanels)[i]->getCp();
        con.row(i) = (*bPanels)[i]->getVerts();
        bN.data.row(i) = (*bPanels)[i]->getBezNormal();
    }
    
    data.push_back(mu);
    data.push_back(pot);
    data.push_back(V);
    data.push_back(Cp);
    data.push_back(bN);
    
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

//void cpCase::writeBuffWake2Data(boost::filesystem::path path, const Eigen::MatrixXd &nodeMat)
//{
//    std::vector<cellDataArray> data;
//    cellDataArray mu("Doublet Strengths"),pot("Velocity Potential");
//    Eigen::MatrixXi con(wake2panels->size(),4);
//    
//    mu.data.resize(wake2panels->size(),1);
//    pot.data.resize(wake2panels->size(),1);
//    for (int i=0; i<wake2panels->size(); i++)
//    {
//        mu.data(i,0) = (*wake2panels)[i]->getMu();
//        pot.data(i,0) = (*wake2panels)[i]->getPotential();
//        con.row(i) = (*wake2panels)[i]->getVerts();
//    }
//    
//    data.push_back(mu);
//    data.push_back(pot);
//    
//    piece wake;
//    wake.pnts = nodeMat;
//    wake.connectivity = con;
//    wake.cellData = data;
//    
//    std::string fname = path.string()+"/buffer2Data-"+std::to_string(timeStep)+".vtu";
//    VTUfile wakeFile(fname,wake);
//}

void cpCase::writeParticleData(boost::filesystem::path path)
{

    Eigen::MatrixXd partMat(particles.size(),3);
    for (int i=0; i<particles.size(); i++)
    {
        partMat.row(i) = particles[i]->getPos();
    }
    
    std::vector<cellDataArray> data;
    cellDataArray strength("Strength");
    Eigen::MatrixXi con(particles.size(),1);
    
    strength.data.resize(particles.size(),3);
    for (int i=0; i<particles.size(); i++)
    {
        strength.data.row(i) = particles[i]->getStrength();
        con(i) = i;
    }
    data.push_back(strength);
    
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

void cpCase::findSeedPoints(){
    
    for(int i=0; i<(*wPanels).size(); i++){
        seedPts.push_back((*wPanels)[i]->partSeedPt(Vinf,dt));
    }
}

void cpCase::findSeedRadii(){
    
    for(int i=0; i<(*wPanels).size(); i++){
        seedRadii.push_back((*wPanels)[i]->getPartRadius(Vinf, dt));
    }
}


void cpCase::collapseBufferWake(){
    particle* p = nullptr;
    
    for(int i=0;i<(*wPanels).size();i++){
        
        Eigen::Vector3d pos = seedPts[i];   // (*wPanels)[i]->partSeedPt(Vinf,dt);
        double radius = seedRadii[i]; //(*wPanels)[i]->getPartRadius(Vinf,dt);
        Eigen::Vector3d strength;
        if(timeStep == 1){
//            strength = (*wPanels)[i]->panToPartStrengthT1();
            strength = (*wPanels)[i]->panToPartStrength();
        }else{
            strength = (*wPanels)[i]->panToPartStrength();
        }
        
        p = new particle(pos, strength, radius);
        particles.push_back(p);
    }
    
//    for(int i=0; i<(*wPanels).size(); i++){ // In a separate loop because the top loop calls neighboring panels...
//        (*wPanels)[i]->manuallySetMu(-1e10);
//    }
}


void cpCase::convectParticles(){
    
    for(int i=0;i<particles.size();i++){
        
        // Freestream influence
        Eigen::Vector3d FSinfl = Vinf; // Initialize with freestream
        
        // Body panel influence
        Eigen::Vector3d bPanInfl = Eigen::Vector3d::Zero();
        for(int j=0;j<(*bPanels).size();j++){
            bPanInfl += (*bPanels)[j]->panelV(particles[i]->getPos());
        }
        
        // Buffer wake influence
        Eigen::Vector3d wPanInfl = Eigen::Vector3d::Zero();
        for(int j=0;j<(*wPanels).size();j++){
            wPanInfl += (*wPanels)[j]->panelV(particles[i]->getPos());
        }
        
        // Particle influence
        Eigen::Vector3d partInfl = Eigen::Vector3d::Zero();
        for(int j=0;j<particles.size();j++){
            if(i != j){ // should be able to take this out with regularization. t
                partInfl += particles[j]->partVelInflGaussian(particles[i]->getPos());
//                partInfl += particles[j]->partVelInfl(particles[i]->getPos());
            }
        }
        
        // Track Particles
        Eigen::Vector3d velOnPart = FSinfl + bPanInfl + wPanInfl + partInfl;
        Eigen::Vector3d newPos = particles[i]->getPos()+velOnPart*dt;
        particles[i]->setPos(newPos);

    }
}

//void cpCase::convectBufferWake(){ //VPP
//    
//    for(int i=0; i<(*wPanels).size(); i++){
//        (*wake2panels)[i]->manuallySetMu((*wPanels)[i]->getMu());
//        (*wPanels)[i]->manuallySetMu(-100000000); // Can take this out later, just want to make sure it doesn't affect anything.
//    }
//    
//}

void cpCase::vortexStretching(){
    // due to particles
    std::vector<Eigen::Matrix3d> stretchingVec;
    std::vector<Eigen::Vector3d> HAstretchVec;
    for(int i=0; i<particles.size(); i++){
        Eigen::Matrix3d stretching;
        stretching = Eigen::Matrix3d::Zero();
        
        Eigen::Vector3d HAstretch;
        HAstretch = Eigen::Vector3d::Zero();
        for(int j=0; j<particles.size(); j++){
            if(i !=j){
                HAstretch += particles[i]->partStretching(particles[j]);
                stretching += particles[i]->partStretchingGaussian(particles[j]);
            }
        }
        
        stretchingVec.push_back(stretching);
        HAstretchVec.push_back(HAstretch);
    }

    
    for(int i=0;i<particles.size();i++){
        Eigen::Vector3d stren = particles[i]->getStrength(); //WWWWW
        Eigen::Vector3d gStretch = stretchingVec[i]*stren;
        
        std::cout << "\n High Alg: " << HAstretchVec[i].x() << " , " << HAstretchVec[i].y() << ", " << HAstretchVec[i].z() << std::endl;
        std::cout << " Gaussian: " << gStretch.x() << " , " << gStretch.y() << ", " << gStretch.z() << std::endl;

        
//        std::cout << "      "  " pre: " << particles[i]->getStrength().x() << " , " << particles[i]->getStrength().y() << ", " << particles[i]->getStrength().z() << std::endl;

//        particles[i]->setStrength(stren+stretchingVec[i]);
//        std::cout << "      "  " post: " << particles[i]->getStrength().x() << " , " << particles[i]->getStrength().y() << ", " << particles[i]->getStrength().z() << std::endl;

    }
    
    
    // due to panels
    //...
    std::cout << "fix sigma for vel infl" << std::endl;
    

}

void cpCase::viscousDiffusion(){
    
    std::vector<Eigen::Vector3d> diffVec;
    for(int i=0; i<particles.size(); i++){
        Eigen::Vector3d diffusion;
        diffusion = Eigen::Vector3d::Zero();
        for(int j=0; j<particles.size(); j++){
                diffusion += particles[i]->partDiffusion(particles[j]);
        }
        diffVec.push_back(diffusion);
    }
    
    for(int i=0;i<particles.size();i++){
//        Eigen::Vector3d stren = particles[i]->getStrength();
        particles[i]->setStrength(particles[i]->getStrength()+diffVec[i]);
        
    }
    
}
