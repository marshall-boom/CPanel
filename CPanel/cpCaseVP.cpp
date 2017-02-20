//
//  cpCaseVP.cpp
//  CPanel
//
//  Created by Connor Sousa on 1/31/17.
//  Copyright (c) 2017 Chris Satterwhite. All rights reserved.
//

#include "cpCaseVP.h"

void cpCaseVP::run(bool printFlag, bool surfStreamFlag, bool stabDerivFlag){
    
    if(unsteady){
        readBodyKinFile();
    }
    
    bool converged = false;
    std::string check = "\u2713";
    setSourceStrengths();
    
    converged = solveMatrixEq();
    writeFilesVP();
    
    bool matrixConversion = false;
    
    if(unsteady){
//            compVelocity(); // This doesn't need to be here unless solving Cl/Cd etc
    }
    
    
    // Initially convect buffer wake and solve
    timestep++;
    convectBufferWake(); ///
    converged = solveVPmatrixEq();
    if (unsteady) { /* compVelocity(); */ }
    writeFilesVP();
    
    // For rest of timesteps
    while (timestep < numSteps){
        
        timestep++;
        
        std::cout << "Time step " << timestep << "/" << numSteps << ".  Flow time = " << timestep*dt;
        std::cout << ".  " << particles.size() << " particles" << std::endl;
        
        convectParticles();
        
        // Turn second row of wake panels into particles
        collapseBufferWake();
        
        // Advance first row of wake panels to second row
        convectBufferWake();
        
        // Influence of body and particles onto particles
        particleStrengthUpdate();
        
        
        // SOMETHING SOMETHING OCTREEE------------------
        if(accelerate){
            
            partOctree.removeData();
            partOctree.setMaxMembers(10); //Barnes Hut
            
//            if(highAccuracy){
//                partOctree.setMaxTheta(0.25);
//            }
            partOctree.addData(particles);
            
            FMM.build(&partOctree);
        }
        
        
        // Influence of particles onto body
        setSourceStrengthsVP();
        
        // Solve system of equations
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
        
        // If forces are desired in middle of sim
        if(unsteady){
            compVelocity();
        }
        
//        writeFilesVP();
//        timestep++;
//        // Move Geometry
//        std::vector<double> bodyKin { 1, 0, 0, 0, 0, 0 };
//        geom->moveGeom(bodyKin);
        
        writeFilesVP();
        
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
        writeFilesVP();
    }
    
    
}

void cpCaseVP::convectBufferWake()
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

void cpCaseVP::setSourceStrengthsVP()
{
    sigmas.resize(bPanels->size());
    
    for (int i=0; i<bPanels->size(); i++)
    {
        Eigen::Vector3d sumVelInfl = Eigen::Vector3d::Zero();
        
        // Influence from particles
        if (accelerate && timestep > 0){ //??
            sumVelInfl += FMM.barnesHutVel((*bPanels)[i]->getCenter());
        }
        else
        {
            for(int j=0; j<particles.size(); j++)
            {
                sumVelInfl += particles[j]->partVelInfl((*bPanels)[i]->getCenter());
            }
        }
        
        // Influence from vortex filaments
        for(int j=0; j<filaments.size(); j++)
        {
            sumVelInfl += filaments[j]->velInfl( (*bPanels)[i]->getCenter() );
        }
        
        (*bPanels)[i]->setSigma( Vinf( (*bPanels)[i]->getCenter() ) + sumVelInfl , 0 );
        sigmas(i) = (*bPanels)[i]->getSigma();
    }
}

bool cpCaseVP::solutionConvergence(){
    
    if(unsteady | (timestep < 5)){
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
    
    if(timestep > 2){
        // Need particles to create the S curve
        trefftzPlaneAnalysis();
    }
    
    //    std::cout << timeStep << " " << CL_trefftz << "  " << CD_trefftz << std::endl;
    
    double changeCD = std::abs((CD_trefftz - CDnm1)/CDnm1);
    double changeCL = std::abs((CL_trefftz - CLnm1)/CLnm1);
    
    //    std::cout << timeStep<< "        " << CL_trefftz << "   " <<changeCL << "        " << CD_trefftz << "       " <<  changeCD  << std::endl;
    
    double val = 0.0001;
    
    if((changeCL < val) && (changeCD < val) && (changeCDnm1 < val) && (changeCDnm1 < val)){
        return true;
    }
    
    return false;
    
}



bool cpCaseVP::solveMatrixEq()
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

bool cpCaseVP::solveVPmatrixEq()
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
    
    for (int i=0; i<wPanels->size(); i++)
    {
        (*wPanels)[i]->setPrevStrength((*wPanels)[i]->getMu());
        (*wPanels)[i]->setMu();
        
        (*wPanels)[i]->setPotential(VinfPlusVecPot((*wPanels)[i]->getCenter()));
        
        (*w2panels)[i]->setPotential(VinfPlusVecPot((*w2panels)[i]->getCenter())); // Included in this loop because there are always the same number of w2pans as w1pans
        
    }
    return converged;
}


bool cpCaseVP::edgeIsUsed(edge* thisEdge, std::vector<edge*> pEdges){
    
    for(int i=0; i<pEdges.size(); i++){
        if(thisEdge == pEdges[i]){
            return true;
        }
    }
    
    return false;
}

void cpCaseVP::collapseBufferWake(){
    
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
        
        particle* p = new particle(pos, strength, radius, {0,0,0}, {0,0,0}, timestep); // Last two are previous pos and strength values used for advanced time stepper.
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

void cpCaseVP::compVelocity()
{
    //  Velocity Survey with known doublet and source strengths
    CM.setZero();
    Eigen::Vector3d moment;
    Fbody = Eigen::Vector3d::Zero();
    bodyPanel* p;
    
    //    std::vector<double> panLift;
    //    bool liftDist = false;
    
    for (int i=0; i<bPanels->size(); i++)
    {
        p = (*bPanels)[i];
        p->computeVelocity(PG,Vinf(p->getCenter()));
        if (unsteady) {
            p->computeCp( Vinf((*bPanels)[i]->getCenter()).norm() , dt ); // Katz 13.169 ///
        }else{
            p->computeCp( Vmag );
        }
        Fbody += -p->getCp()*p->getArea()*p->getBezNormal()/params->Sref;
        moment = p->computeMoments(params->cg);
        CM(0) += moment(0)/(params->Sref*params->bref);
        CM(1) += moment(1)/(params->Sref*params->cref);
        CM(2) += moment(2)/(params->Sref*params->bref);
        
    }
    Fwind = bodyToWind(Fbody);
    //    CL.push_back(Fbody.z());
    //    CM_x.push_back(CM(1));
    
    std::cout << Fbody.z() << std::endl;
    
    //    std::cout << std::endl;
}

void cpCaseVP::trefftzPlaneAnalysis(){
    
    std::vector<wake*> wakes = geom->getWakes();
    CL_trefftz = 0;
    CD_trefftz = 0;
    for (int i=0; i<wakes.size(); i++)
    {
        wakes[i]->trefftzPlaneVP(Vmag,params->Sref, &particles, timestep);
        CL_trefftz += wakes[i]->getCL()/PG;
        CD_trefftz += wakes[i]->getCD()/pow(PG,2);
    }
}

void cpCaseVP::particleStrengthUpdate(){
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
                dAlpha_stretch = FMM.barnesHutStretch(particles[i]); ///
                dAlpha_diff = FMM.barnesHutDiff(particles[i]);
            }
            else
            {
                // Stretching from particles
                for(int j=0; j<particles.size(); j++)
                {
                    if(i!=j) // Kroneger Delta Function
                    {
                        dAlpha_diff += particles[i]->viscousDiffusion(particles[j]);
                        dAlpha_stretch += particles[i]->vortexStretching(particles[j]);
                    }
                }
            }
            stretchDiffVec.push_back( dAlpha_diff + dAlpha_stretch );
        }
        else // Treat normally
        {
            // Influence from Particles
            if(accelerate && timestep > 3){
                dAlpha_stretch += FMM.barnesHutStretch(particles[i]);
                dAlpha_diff += FMM.barnesHutDiff(particles[i]);
            }else
            {
                for(int j=0; j<particles.size(); j++)
                {
                    if(i!=j){ // Kroneger Delta Function
                        dAlpha_diff += particles[i]->viscousDiffusion(particles[j]);
                        dAlpha_stretch += particles[i]->vortexStretching(particles[j]);
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

void cpCaseVP::convectParticles()
{
    std::vector<Eigen::Vector3d> newPartPositions;
    
    for(int i=0;i<particles.size();i++){
        
        Eigen::Vector3d newPos;
        
        // Adams bashforth scheme
        Eigen::Vector3d velOnPart = velocityInflFromEverything(particles[i]);
        
        if(particles[i]->getPrevVelInfl().isZero())
        {
            newPos = particles[i]->pos + dt*velOnPart;
        }
        else{
            newPos = particles[i]->pos + dt*(1.5*velOnPart - 0.5*particles[i]->getPrevVelInfl());
        }
        
        particles[i]->setPrevVelInfl(velOnPart);
        newPartPositions.push_back(newPos);
        //        }
    }
    
    for(int i=0;i<particles.size();i++){
        particles[i]->setPos(newPartPositions[i]);
        
    }
}

double bodyKin(int dum1, int dum2){return 3.4;}; ///

Eigen::Vector3d cpCaseVP::Vinf(Eigen::Vector3d POI)
{
    if (!unsteady)
    {
        return windToBody( Vmag , alpha , beta );
    }
    else
    {
        Eigen::Vector3d localVel;
        
        // U = U3 + (-q*z + r*y)
        localVel.x() = bodyKin(timestep, 0) - bodyKin(timestep, 4)*POI.z() + bodyKin(timestep, 5)*POI.y();
        
        // V = V3 + (-r*x + p*z)
        localVel.y() = bodyKin(timestep, 1) - bodyKin(timestep, 5)*POI.x() + bodyKin(timestep, 3)*POI.z();
        
        // W = W3 + (-p*y + q*x)
        localVel.z() = bodyKin(timestep, 2) - bodyKin(timestep, 3)*POI.y() + bodyKin(timestep, 4)*POI.x();
        
        return localVel;
    }
}

Eigen::Vector3d cpCaseVP::VinfPlusVecPot(Eigen::Vector3d POI)
{
    Eigen::Vector3d vInfluence = Vinf(POI);
    
    // Particles
    if(accelerate && timestep > 2)
    {
        vInfluence += FMM.barnesHutVel(POI);
    }else
    {
        for (int i=0; i<particles.size(); i++) {
            vInfluence += particles[i]->partVelInfl(POI);
        }
    }
    
    
    // Filaments
    for (int i=0; i<filaments.size(); i++)
    {
        vInfluence += filaments[i]->velInfl(POI);
    }
    
    return vInfluence;
}

Eigen::Vector3d cpCaseVP::edgeStrength(wakePanel* pan, edge* curEdge, int edgeNum){
    //wait, why can't this go in the panel class?? Ithink it can
    
    
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

Eigen::Vector3d cpCaseVP::rungeKuttaStepper( Eigen::Vector3d POI ){
    
    // RK4 algorithm makes symmeterized influences VERY impractical
    
    Eigen::Vector3d k1, k2, k3, k4;
    
    k1 = velocityInflFromEverything(POI);
    k2 = velocityInflFromEverything(POI+k1*dt/2);
    k3 = velocityInflFromEverything(POI+k2*dt/2);
    k4 = velocityInflFromEverything(POI+k3*dt);
    
    return POI + dt*(k1/6 + k2/3 + k3/3 + k4/6);
    
}

Eigen::Vector3d cpCaseVP::velocityInflFromEverything( Eigen::Vector3d POI )
{
    // Freestream influence
    Eigen::Vector3d velOnPart = Vinf(POI);
    //    std::cout << velOnPart.x() << " , " << velOnPart.y() << " , " << velOnPart.z() << std::endl;
    
    
    // Particle influence
    if(accelerate && timestep > 3){ //??
        velOnPart += FMM.barnesHutVel(POI);
    }
    else{
        for(int j=0;j<particles.size();j++)
        {
            velOnPart += particles[j]->partVelInfl(POI);
        }
    }
    
    
    // Body panel influence
    for(int j=0;j<(*bPanels).size();j++){
        velOnPart += (*bPanels)[j]->panelV(POI);
    }
    
    
    // Buffer wake influence
    for(int j=0;j<(*wPanels).size();j++){
        velOnPart += (*wPanels)[j]->panelV(POI);
    }
    
    for(int j=0;j<(*w2panels).size();j++){
        velOnPart += (*w2panels)[j]->panelV(POI);
    }
    
    // Vortex Filament influence
    for(int i=0; i<filaments.size(); i++){
        velOnPart +=filaments[i]->velInfl(POI);
    }
    
    //    std::cout << velOnPart.x() << " , " << velOnPart.y() << " , " << velOnPart.z() << std::endl;
    
    return velOnPart;
}

Eigen::Vector3d cpCaseVP::velocityInflFromEverything(particle* part)
{
    // Function is overloaded for symmeterized velocity influence for particles
    Eigen::Vector3d pos = part->pos;
    
    // Freestream influence
    Eigen::Vector3d velOnPart = Vinf(pos);
    
    // If particle reaches very far field condition
    if ((pos-Eigen::Vector3d::Zero()).norm() > 20*params->cref)
    {
        return velOnPart;
    }
    
    // Particle influence
    if(accelerate && timestep > 2){
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
    
//    std::cout << velOnPart.x() << " , " << velOnPart.y() << " , " << velOnPart.z() << std::endl;
    
    return velOnPart;
}

void cpCaseVP::writeFilesVP()
{
    std::stringstream caseLabel;
    caseLabel << "/V" << Vmag << "_Mach" << mach << "_alpha" << alpha << "_beta" << beta;
    boost::filesystem::path subdir = boost::filesystem::current_path().string()+caseLabel.str();
    if (!boost::filesystem::exists(subdir))
    {
        boost::filesystem::create_directories(subdir);
    }
    Eigen::MatrixXd nodeMat = geom->getNodePnts();
    
    writeBodyDataVP(subdir,nodeMat);
    writeWakeDataVP(subdir,nodeMat);
    writeBuffWake2Data(subdir,nodeMat);
    writeSpanwiseData(subdir);
    writeParticleData(subdir);
    writeFilamentData(subdir);
    
    if (false) {
//        createVolMesh();
//        writeVolMeshData(subdir, pts, cells);
    }
    
    if (params->surfStreamFlag)
    {
        writeBodyStreamlines(subdir);
    }
    
}

void cpCaseVP::writeBodyDataVP(boost::filesystem::path path,const Eigen::MatrixXd &nodeMat)
{
    std::vector<cellDataArray> data;
    cellDataArray mu("Doublet Strengths"),sigma("Source Strengths"),pot("Velocity Potential"),V("Velocity"),Cp("Cp"),bN("bezNormals");
    Eigen::MatrixXi con(bPanels->size(),3);
    mu.data.resize(bPanels->size(),1);
    sigma.data.resize(bPanels->size(),1);
    pot.data.resize(bPanels->size(),1);
    V.data.resize(bPanels->size(),3);
    Cp.data.resize(bPanels->size(),1);
    bN.data.resize(bPanels->size(),3);
    for (int i=0; i<bPanels->size(); i++)
    {
        mu.data(i,0) = (*bPanels)[i]->getMu();
        sigma.data(i,0) = (*bPanels)[i]->getSigma();
        pot.data(i,0) = (*bPanels)[i]->getPotential();
        V.data.row(i) = (*bPanels)[i]->getGlobalV();
        Cp.data(i,0) = (*bPanels)[i]->getCp();
        con.row(i) = (*bPanels)[i]->getVerts();
        bN.data.row(i) = (*bPanels)[i]->getBezNormal();
    }
    
    data.push_back(mu);
    data.push_back(sigma);
    data.push_back(pot);
    data.push_back(V);
    data.push_back(Cp);
    data.push_back(bN);
    
    piece body;
    body.pnts = nodeMat;
    body.connectivity = con;
    body.cellData = data;
    
    std::string fname = path.string()+"/surfaceData-"+std::to_string(timestep)+".vtu";
    VTUfile bodyFile(fname,body);
}

void cpCaseVP::writeWakeDataVP(boost::filesystem::path path, const Eigen::MatrixXd &nodeMat)
{
    std::vector<cellDataArray> data;
    cellDataArray mu("Doublet Strengths"),pot("Velocity Potential");
    Eigen::MatrixXi con(wPanels->size(),(*wPanels)[0]->getVerts().size()); // Assumes wake won't mix tris and quads
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
    
    std::string fname = path.string()+"/wakeData-"+std::to_string(timestep)+".vtu";
    VTUfile wakeFile(fname,wake);
}


void cpCaseVP::writeBuffWake2Data(boost::filesystem::path path, const Eigen::MatrixXd &nodeMat)
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
    std::string fname = path.string()+"/bufferWake2Data-"+std::to_string(timestep)+".vtu";
    VTUfile wakeFile(fname,wake);
}

void cpCaseVP::writeFilamentData(boost::filesystem::path path)
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
    std::string fname = path.string()+"/filaments-"+std::to_string(timestep)+".vtu";
    VTUfile filFile(fname,fils);
}

void cpCaseVP::writeParticleData(boost::filesystem::path path)
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
    
    std::string fname = path.string()+"/particleData-"+std::to_string(timestep)+".vtu";
    VTUfile partFile(fname,parts);
}


void cpCaseVP::readBodyKinFile(){
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

//void cpCase::particleStrengthUpdate(){
//    // This function uses the combined vortex stretching and diffusion equation used by Wincklemans for a regularized vortex core with high algebraic smoothing. (he refers to it as the strength update equation.)
//
//    std::vector<Eigen::Vector3d> stretchDiffVec; // Creating a vector of diffusion values because the strength change needs to be set after all particle influences have been calculated
//    for(int i=0; i<particles.size(); i++)
//    {
//        Eigen::Vector3d dAlpha = Eigen::Vector3d::Zero();
//
//        // FarField2 condition
//        if((particles[i]->pos-Eigen::Vector3d::Zero()).norm() > 20*params->cref)
//        {
//            stretchDiffVec.push_back(dAlpha);
//        }
//        // FarField1 condition
//        else if((particles[i]->pos-Eigen::Vector3d::Zero()).norm() > 12*params->cref)
//        {
//            // Stretching from particles
//            for(int j=0; j<particles.size(); j++)
//            {
//                dAlpha += particles[i]->partStrengthUpdate(particles[j]);
//            }
//            stretchDiffVec.push_back(dAlpha);
//        }
//        // Treat normally
//        else{
//
//            // Stretching from particles
//            for(int j=0; j<particles.size(); j++)
//            {
//                dAlpha += particles[i]->partStrengthUpdate(particles[j]);
//            }
//
//            // Stretching from body panels
//            for(int j=0; j<(*bPanels).size(); j++)
//            {
//                dAlpha += (*bPanels)[j]->partStretching(particles[i]);
//            }
//
//            // Stretching from wake panels
//            for(int j=0;j<(*wPanels).size(); j++)
//            {
//                dAlpha += (*wPanels)[j]->partStretching(particles[i]);
//            }
//
//            stretchDiffVec.push_back(dAlpha);
//        }
//    }
//
//    // No need for Kutta accuracy
//    for(int i=0;i<particles.size();i++){
//        Eigen::Vector3d newStrength;
//        if(particles[i]->getprevStrengthUpdate().isZero())
//        {
//            newStrength = particles[i]->strength + stretchDiffVec[i]*dt;
//        }else{
//            //Adams bashforth
//            newStrength = particles[i]->strength + dt*(1.5*stretchDiffVec[i] - 0.5*particles[i]->getprevStrengthUpdate());
//        }
//        particles[i]->setprevStrengthUpdate(stretchDiffVec[i]);
//        particles[i]->setStrength(newStrength);
//    }
//
//}

