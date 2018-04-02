/*******************************************************************************
 * Copyright (c) 2017 Connor Sousa
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
 *    Connor Sousa - initial code and implementation
 *    David D. Marshall - porting to GoogleTest
 ******************************************************************************/

#include "cpCaseVP.h"

void cpCaseVP::run(bool printFlag, bool surfStreamFlag, bool stabDerivFlag){
    
    if(unsteady){
        readBodyKinFile();
        solnMat.resize(numSteps, 12);
    }
    
    bool matrix_convergence = false;
    std::string check = "\u2713";
    std::setprecision(5);
    if(printFlag) std::cout << std::setw(11) << " Time step " << std::to_string(timeStep) + "/" + std::to_string(numSteps) + ".  " << "Flow time = " << std::setw(5) << std::to_string(timeStep*dt) + ". " << std::endl;
    
    // First time step
    setSourceStrengths();
    matrix_convergence = solveMatrixEq();
    if (unsteady) compVelocity();
    writeFilesVP();
    //    moveGeometry();
    
    
    // Initially convect buffer wake and solve
    timeStep++;
    if(printFlag) std::cout << std::setw(11) << " Time step " << std::to_string(timeStep) + "/" + std::to_string(numSteps) + ".  " << "Flow time = " << std::setw(5) << std::to_string(timeStep*dt) + ". " << std::endl;
    
    convectBufferWake();
    matrix_convergence = solveVPmatrixEq();
    
    if (unsteady) compVelocity();
    writeFilesVP();
    //    moveGeometry();
    
    
    // For rest of timesteps
    while (timeStep < numSteps){
        
        timeStep++;
        if(printFlag) std::cout << std::setw(11) << " Time step " << (std::to_string(timeStep) + "/" + std::to_string(numSteps) + ".  ") << "Flow time = " << std::setw(8) << std::to_string(timeStep*dt) + ".  " << particles.size() << " particles" << std::endl;
        
        convectParticles();
        
        // Turn second row of wake panels into particles
        collapseBufferWake();
        
        // Move Geometry
        //        moveGeometry();
        
        // Advance first row of wake panels to second row
        convectBufferWake();
        
        // Influence of body and particles onto particles
        particleStrengthUpdate();
        
        // Build/update particle octree
        if(accelerate){
            partOctree.removeData();
            partOctree.setMaxMembers(10); //Barnes Hut
            partOctree.addData(particles);
            
            FMM.build(&partOctree);
        }
        
        // Influence of particles onto body
        setSourceStrengthsVP();
        
        // Solve system of equations
        matrix_convergence = solveVPmatrixEq();
        
        // Check for simulation convergence
        if(!matrix_convergence){
            std::cout << "*** Warning : Matrix solution did not converge ***" << std::endl;
            std::exit(0);
        }else{
            bool soln_conv = solutionConvergence();
            if (soln_conv) {
                if(printFlag) std::cout << "Converged..." << std::endl;
                break;
            }
        }
        
        if(unsteady) compVelocity();
        
        writeFilesVP();
        
    }
    
    if(printFlag) std::cout << "\n (\u2713 - Complete, X - Not Requested)\n"<< std::setw(17) << std::left << " Solve System" << std::setw(17) << std::left << "Surface Data" << std::setw(18) << std::left << "Trefftz Plane" <<  std::setw(14) << std::left << "Streamlines" << std::setw(24) << std::left << "Stability Derivatives" << std::setw(23) << std::left << "Volume Mesh" << std::endl << std::setw(19) << std::left << " " + check << std::flush;
    
    if(!unsteady) compVelocity(); // calculated each time step in unsteady
    if(printFlag) std::cout << std::setw(19) << std::left << check << std::flush;
    
    trefftzPlaneAnalysisVP();
    if(printFlag) std::cout << std::setw(20) << std::left << check << std::flush;
    
    
    if (surfStreamFlag)
    {
        createStreamlines();
        if(printFlag) std::cout << std::setw(14) << std::left << check << std::flush;
    }
    else
    {
        if(printFlag) std::cout << std::setw(14) << std::left << "X" << std::flush;
    }
    
    if (stabDerivFlag)
    {
        stabilityDerivativesVP();
        if(printFlag) std::cout << std::setw(24) << std::left << check << std::flush;
    }
    else
    {
        if(printFlag) std::cout << std::setw(24) << std::left << "X" << std::flush;
    }
    
    
    if (params->volMeshFlag) {
        createVolMesh();
        populateVolMesh();
        if (printFlag) std::cout << std::setw(23) << std::left << check << std::endl;
    }
    else
    {
        if (printFlag) std::cout << std::setw(23) << std::left << "X" << std::endl;
    }
    
    
    if (!matrix_convergence && printFlag)
    {
        std::cout << "*** Warning : Solution did not converge ***" << std::endl;
    }
    
    if (printFlag) writeFilesVP();
    
    
}

void cpCaseVP::convectBufferWake(){
    wake2Doublets.resize((*w2panels).size());
    for (int i=0; i<(*w2panels).size(); i++)
    {
        (*w2panels)[i]->setPrevStrength((*w2panels)[i]->getMu());
        double parentMu = (*w2panels)[i]->getBufferParent()->getMu();
        (*w2panels)[i]->setMu(parentMu);
        wake2Doublets[i] = parentMu;
    }
    
}


void cpCaseVP::setSourceStrengthsVP(){
    sigmas.resize(bPanels->size());
    
    for (int i=0; i<bPanels->size(); i++)
    {
        Eigen::Vector3d sumVelInfl = Eigen::Vector3d::Zero();
        
        // Influence from particles
        if (accelerate){
            sumVelInfl += FMM.barnesHutVel((*bPanels)[i]->getCenter());
        }
        else
        {
            for(int j=0; j<particles.size(); j++)
            {
                sumVelInfl += particles[j]->partVelInflGaussian((*bPanels)[i]->getCenter());
            }
        }
        
        // Influence from vortex filaments
        for(int j=0; j<filaments.size(); j++)
        {
            sumVelInfl += filaments[j]->velInfl( (*bPanels)[i]->getCenter() );
        }
        
        (*bPanels)[i]->setSigma( Vinf((*bPanels)[i]->getCenter()) + sumVelInfl , 0 );
        
        sigmas(i) = (*bPanels)[i]->getSigma();
    }
}

bool cpCaseVP::solutionConvergence(){
    
    if(unsteady | (timeStep < 5) | params->stepsSetMaunally){
        return false;
    }
    
    // Check for simulation convergence. Will use the Trefftz plane because it is so cheap. Medium size mesh takes <1% of time step time.
    
    // Store previous time step values
    double prevCD = CD_trefftz;
    double prevCL = CL_trefftz;
    
    trefftzPlaneAnalysisVP();
    
    double changeCD = std::abs((CD_trefftz - prevCD)/prevCD);
    double changeCL = std::abs((CL_trefftz - prevCL)/prevCL);
    
    double val = 0.0001; // 0.01%
    
    if((changeCL < val) && (changeCD < val)){
        return true;
    }
    
    return false;
    
}



bool cpCaseVP::solveMatrixEq(){
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

bool cpCaseVP::solveVPmatrixEq(){
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
                strength += (*w2panels)[i]->edgeStrength( pEdges[j], j ); // Don't need to pass in pEdges...
            }
        }
        
        Eigen::Vector3d pos = rungeKuttaStepper((*w2panels)[i]->getCenter());
        //        Eigen::Vector3d pos = (*w2panels)[i]->getCenter(); // For moving geometry
        Eigen::Vector3d ptVel = Vinf(pos);
        double radius = (*w2panels)[i]->getPartRadius(ptVel,dt); // VinfLocal
        
        particle* p = new particle(pos, strength, radius, {0,0,0}, {0,0,0}, timeStep); // Zeroes are previous pos and strength values used for Adams-Bashforth scheme
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

void cpCaseVP::compVelocity(){
    //  Velocity Survey with known doublet and source strengths
    CM.setZero();
    Eigen::Vector3d moment;
    Fbody = Eigen::Vector3d::Zero();
    bodyPanel* p;
    
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
    
    if(unsteady){
        if(timeStep > 3) trefftzPlaneAnalysisVP(); // No convergence check for unsteady sims
        
        //   flow_time | CL_tr | CDi_tr | CN | CA | CY | CL | CD | CY | Cm | Cl | Cn |
        solnMat(timeStep-1,0) = timeStep*dt;
        solnMat(timeStep-1,1) = CL_trefftz;
        solnMat(timeStep-1,2) = CD_trefftz;
        
        solnMat(timeStep-1,3) = Fbody.x();
        solnMat(timeStep-1,4) = Fbody.y();
        solnMat(timeStep-1,5) = Fbody.z();
        
        solnMat(timeStep-1,6) = Fwind.x();
        solnMat(timeStep-1,7) = Fwind.y();
        solnMat(timeStep-1,8) = Fwind.z();
        
        solnMat(timeStep-1,9)  = CM.x();
        solnMat(timeStep-1,10) = CM.y();
        solnMat(timeStep-1,11) = CM.z();
    }
    
    
}

void cpCaseVP::trefftzPlaneAnalysisVP(){
    
    std::vector<wake*> wakes = geom->getWakes();
    CL_trefftz = 0;
    CD_trefftz = 0;
    for (int i=0; i<wakes.size(); i++)
    {
        wakes[i]->trefftzPlaneVP(Vmag,params->Sref, &particles, timeStep);
        CL_trefftz += wakes[i]->getCL()/PG;
        CD_trefftz += wakes[i]->getCD()/pow(PG,2);
    }
}

void cpCaseVP::stabilityDerivativesVP()
{
    double delta = 0.5;
    double dRad = delta*M_PI/180;
    cpCaseVP dA(geom,Vmag,alpha+delta,beta,mach,params);
    cpCaseVP dB(geom,Vmag,alpha,beta+delta,mach,params);
    
    dA.run(false,false,false);
    dB.run(false,false,false);
    
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


void cpCaseVP::particleStrengthUpdate(){
    // This function uses the update equations found in 'Vortex Methods for DNS of...' by Plouhmhans. It uses the particle strength exchange for the viscous diffusion
    
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
            if(accelerate)
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
            if(accelerate && timeStep > 3){
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

void cpCaseVP::convectParticles(){
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


Eigen::Vector3d cpCaseVP::Vinf(Eigen::Vector3d POI){
    if (!unsteady)
    {
        return windToBody( Vmag , alpha , beta );
    }
    else
    {
        Eigen::Vector3d localVel;
        
        // U = U3 + (-q*z + r*y)
        localVel.x() = bodyKin(timeStep-1, 0) - bodyKin(timeStep-1, 4)*POI.z() + bodyKin(timeStep-1, 5)*POI.y();
        
        // V = V3 + (-r*x + p*z)
        localVel.y() = bodyKin(timeStep-1, 1) - bodyKin(timeStep-1, 5)*POI.x() + bodyKin(timeStep-1, 3)*POI.z();
        
        // W = W3 + (-p*y + q*x)
        localVel.z() = bodyKin(timeStep-1, 2) - bodyKin(timeStep-1, 3)*POI.y() + bodyKin(timeStep-1, 4)*POI.x();
        
        return localVel;
    }
}

Eigen::Vector3d cpCaseVP::VinfPlusVecPot(Eigen::Vector3d POI){
    Eigen::Vector3d vInfluence = Vinf(POI);
    
    // Particles
    if(accelerate && timeStep > 2)
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



Eigen::Vector3d cpCaseVP::rungeKuttaStepper( Eigen::Vector3d POI ){
    
    // RK4 algorithm makes symmeterized influences VERY impractical
    
    Eigen::Vector3d k1, k2, k3, k4;
    
    k1 = velocityInflFromEverything(POI);
    k2 = velocityInflFromEverything(POI+k1*dt/2);
    k3 = velocityInflFromEverything(POI+k2*dt/2);
    k4 = velocityInflFromEverything(POI+k3*dt);
    
    return POI + dt*(k1/6 + k2/3 + k3/3 + k4/6);
    
}

Eigen::Vector3d cpCaseVP::velocityInflFromEverything( Eigen::Vector3d POI ){
    // Freestream influence
    Eigen::Vector3d velOnPart = Vinf(POI);
    
    // Particle influence
    if(accelerate && timeStep > 3){
        velOnPart += FMM.barnesHutVel(POI);
    }
    else{
        for(int j=0;j<particles.size();j++)
        {
            velOnPart += particles[j]->partVelInflGaussian(POI);
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

Eigen::Vector3d cpCaseVP::velocityInflFromEverything(particle* part){
    
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
    if(accelerate && timeStep > 2){
        velOnPart += FMM.barnesHutVel(part);
    }
    else{
        for(int j=0;j<particles.size();j++)
        {
            velOnPart += particles[j]->partVelInflGaussian(pos);
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

void cpCaseVP::writeFilesVP(){
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
    
    if (params->volMeshFlag && cells.size() > 0) //
    {
        writeVolMeshData(subdir, pts, cells);
    }
    
    if (params->surfStreamFlag)
    {
        writeBodyStreamlines(subdir);
    }
    
}

void cpCaseVP::writeBodyDataVP(boost::filesystem::path path,const Eigen::MatrixXd &nodeMat){
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
    
    std::string fname = path.string()+"/surfaceData-"+std::to_string(timeStep)+".vtu";
    VTUfile bodyFile(fname,body);
}

void cpCaseVP::writeWakeDataVP(boost::filesystem::path path, const Eigen::MatrixXd &nodeMat){
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
    
    std::string fname = path.string()+"/wakeData-"+std::to_string(timeStep)+".vtu";
    VTUfile wakeFile(fname,wake);
}


void cpCaseVP::writeBuffWake2Data(boost::filesystem::path path, const Eigen::MatrixXd &nodeMat){
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

void cpCaseVP::writeFilamentData(boost::filesystem::path path){
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

void cpCaseVP::writeParticleData(boost::filesystem::path path){
    
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


void cpCaseVP::populateVolMesh(){
    
    // Clear Past Timestep Mesh
    volMeshDat.velocity.clear();
    volMeshDat.coef_press.clear();
    
    for (int i=0; i<cells.size(); i++) {
        Eigen::Vector3d velInCell = velocityInflFromEverything(volMeshDat.cellCenter[i]);
        volMeshDat.velocity.push_back(velInCell);
    }
    
    for (int i=0; i<cells.size(); i++) {
        // Incompressible Bernoulli equation
        double Cp = pow( volMeshDat.velocity[i].norm() / Vmag , 2 );
        volMeshDat.coef_press.push_back(Cp);
    }
    
}




//void cpCaseVP::moveGeometry(){
//    // Steps to continue development include (but are not limited to)
//        // In the 'velocityInflFromEverything' functions, initialize with zero velocity instead of freestream
//        // Uncomment all 'moveGeometry' calls
//        // The 'collapseBufferWake' function must seed particle at the panel center instead of tracking it with freestream.
//        // Instead of the 'bodyMovement' vector, you'll probably want to use Eigen::MatrixXd and just slice with .row
//
//    std::vector<double> bodyMovement = { -99.619, 0, -8.7156, 0, 0, 0 };
//
//    geom->moveGeom(bodyMovement);
//
//    // Filaments aren't part of geom... but they should be.
//    for (int i=0; i<filaments.size(); i++) {
//        filaments[i]->moveFilament(bodyMovement, dt);
//    }
//
//
//}





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




































////
////  cpCaseVP.cpp
////  CPanel_BEM
////
////  Created by Connor Sousa on 5/11/17.
////  Copyright (c) 2017 Chris Satterwhite. All rights reserved.
////
//
//#include "cpCaseVP.h"
//
//void cpCaseVP::run(bool printFlag, bool surfStreamFlag, bool stabDerivFlag){
//    
//    if(unsteady){
//        readBodyKinFile();
////        solnMat.resize(numSteps, 12);
//    }
//    
//    bool matrix_convergence = false;
//    std::string check = "\u2713";
//    std::setprecision(5);
//    if(printFlag) std::cout << std::setw(11) << " Time step " << std::to_string(timeStep) + "/" + std::to_string(numSteps) + ".  " << "Flow time = " << std::setw(5) << std::to_string(timeStep*dt) + ". " << std::endl;
//    
//    // First time step
//    setSourceStrengths();
//    matrix_convergence = solveMatrixEq();
//    if (unsteady) compVelocity();
////    writeFilesVP();111
//    writeFiles();
//    
//    //    moveGeometry();
//    
//    
//    // Initially convect buffer wake and solve
//    timeStep++;
//    if(printFlag) std::cout << std::setw(11) << " Time step " << std::to_string(timeStep) + "/" + std::to_string(numSteps) + ".  " << "Flow time = " << std::setw(5) << std::to_string(timeStep*dt) + ". " << std::endl;
//    
//    convectBufferWake();
//    matrix_convergence = solveVPmatrixEq();
//    
//    if (unsteady) compVelocity();
////    writeFilesVP();111
//    writeFiles();
//    //    moveGeometry();
//    
//    
//    // For rest of timesteps
//    while (timeStep < numSteps){
//        
//        timeStep++;
//        if(printFlag) std::cout << std::setw(11) << " Time step " << (std::to_string(timeStep) + "/" + std::to_string(numSteps) + ".  ") << "Flow time = " << std::setw(8) << std::to_string(timeStep*dt) + ".  " << particles.size() << " particles" << std::endl;
//        
//        convectParticles();
//        
//        // Turn second row of wake panels into particles
//        collapseBufferWake();
//        
//        // Move Geometry
//        //        moveGeometry();
//        
//        // Advance first row of wake panels to second row
//        convectBufferWake();
//        
//        // Influence of body and particles onto particles
////        particleStrengthUpdate();111
//        particleStrengthUpdateGaussian();
//        
//        // Build/update particle octree
//        if(accelerate){
//            partOctree.removeData();
//            partOctree.setMaxMembers(10); //Barnes Hut
//            partOctree.addData(particles);
//            
//            FMM.build(&partOctree);
//        }
//        
//        // Influence of particles onto body
////        setSourceStrengthsVP();111
//        setSourceStrengths();
//        
//        // Solve system of equations
//        matrix_convergence = solveVPmatrixEq();
//        
//        // Check for simulation convergence
//        if(!matrix_convergence){
//            std::cout << "*** Warning : Matrix solution did not converge ***" << std::endl;
//            std::exit(0);
//        }else{
//            bool soln_conv = solutionConvergence();
//            if (soln_conv) {
//                if(printFlag) std::cout << "Converged..." << std::endl;
//                break;
//            }
//        }
//        
//        if(unsteady) compVelocity();
//        
////        writeFilesVP();111
//        writeFiles();
//        
//    }
//    
//    if(printFlag) std::cout << "\n (\u2713 - Complete, X - Not Requested)\n"<< std::setw(17) << std::left << " Solve System" << std::setw(17) << std::left << "Surface Data" << std::setw(18) << std::left << "Trefftz Plane" <<  std::setw(14) << std::left << "Streamlines" << std::setw(24) << std::left << "Stability Derivatives" << std::setw(23) << std::left << "Volume Mesh" << std::endl << std::setw(19) << std::left << " " + check << std::flush;
//    
//    if(!unsteady) compVelocity(); // calculated each time step in unsteady
//    if(printFlag) std::cout << std::setw(19) << std::left << check << std::flush;
//    
//    trefftzPlaneAnalysisVP();
//    if(printFlag) std::cout << std::setw(20) << std::left << check << std::flush;
//    
//    
//    if (surfStreamFlag)
//    {
//        createStreamlines();
//        if(printFlag) std::cout << std::setw(14) << std::left << check << std::flush;
//    }
//    else
//    {
//        if(printFlag) std::cout << std::setw(14) << std::left << "X" << std::flush;
//    }
//    
////    if (stabDerivFlag)
////    {
////        stabilityDerivativesVP();
////        if(printFlag) std::cout << std::setw(24) << std::left << check << std::flush;
////    }
////    else
////    {
////        if(printFlag) std::cout << std::setw(24) << std::left << "X" << std::flush;
////    }
//    
//    
////    if (params->volMeshFlag) {
////        createVolMesh();
////        populateVolMesh();
////        if (printFlag) std::cout << std::setw(23) << std::left << check << std::endl;
////    }
////    else
////    {
////        if (printFlag) std::cout << std::setw(23) << std::left << "X" << std::endl;
////    }
//    
//    
//    if (!matrix_convergence && printFlag)
//    {
//        std::cout << "*** Warning : Solution did not converge ***" << std::endl;
//    }
//    
////    if (printFlag) writeFilesVP();111
//    if (printFlag) writeFiles();
//    
//    
//}
//
//Eigen::Vector3d cpCaseVP::Vinf(Eigen::Vector3d POI)
//{
//    if (!unsteady)
//    {
//        return windToBody( Vmag , alpha , beta );
//    }
//    else
//    {
//        Eigen::Vector3d localVel;
//        
//        // U = U3 + (-q*z + r*y)
//        localVel.x() = bodyKin(timeStep, 0) - bodyKin(timeStep, 4)*POI.z() + bodyKin(timeStep, 5)*POI.y();
//        
//        // V = V3 + (-r*x + p*z)
//        localVel.y() = bodyKin(timeStep, 1) - bodyKin(timeStep, 5)*POI.x() + bodyKin(timeStep, 3)*POI.z();
//        
//        // W = W3 + (-p*y + q*x)
//        localVel.z() = bodyKin(timeStep, 2) - bodyKin(timeStep, 3)*POI.y() + bodyKin(timeStep, 4)*POI.x();
//        
//        return localVel;
//    }
//}
//
//
//
//
//
//
//
//void cpCaseVP::setSourceStrengths()
//{
//    sigmas.resize(bPanels->size());
//    for (int i=0; i<bPanels->size(); i++)
//    {
//        Eigen::Vector3d sumVelInfl = Eigen::Vector3d::Zero();
//        
//        if (accelerate && timeStep > 0){
//            sumVelInfl += FMM.barnesHutVel((*bPanels)[i]->getCenter());
//        }
//        else
//        {
//            for(int j=0; j<particles.size(); j++)
//            {
//                sumVelInfl += particles[j]->partVelInflGaussian((*bPanels)[i]->getCenter());
//            }
//        }
//        
//        for(int j=0; j<filaments.size(); j++)
//        {
//            sumVelInfl += filaments[j]->velInfl( (*bPanels)[i]->getCenter() );
//        }
//        
//        (*bPanels)[i]->setSigma( Vinf( (*bPanels)[i]->getCenter() ) + sumVelInfl , 0 ); //function called VinfLocal(POI), and this func won't change
//        sigmas(i) = (*bPanels)[i]->getSigma();
//    }
//}
//
//bool cpCaseVP::solveMatrixEq()
//{
//    bool converged = true;
//    
//    Eigen::MatrixXd* A = geom->getA();
//    Eigen::MatrixXd* B = geom->getB();
//    Eigen::VectorXd RHS = -(*B)*sigmas;
//    Eigen::VectorXd doubletStrengths(bPanels->size());
//    
//    
//    Eigen::BiCGSTAB<Eigen::MatrixXd> res;
//    res.compute((*A));
//    doubletStrengths = res.solve(RHS);
//    if (res.error() > pow(10,-10))
//    {
//        converged = false;
//    }
//    
//    for (int i=0; i<bPanels->size(); i++)
//    {
//        (*bPanels)[i]->setMu(doubletStrengths(i));
//        (*bPanels)[i]->setPotential(Vinf((*bPanels)[i]->getCenter()));
//    }
//    
//    for (int i=0; i<wPanels->size(); i++)
//    {
//        (*wPanels)[i]->setMu();
//        (*wPanels)[i]->setPotential(Vinf((*wPanels)[i]->getCenter()));
//    }
//    
//    return converged;
//}
//
//bool cpCaseVP::solveVPmatrixEq()
//{
//    bool converged = true;
//    
//    Eigen::MatrixXd* A = geom->getA();
//    Eigen::MatrixXd* B = geom->getB();
//    Eigen::MatrixXd* C = geom->getC();
//    Eigen::VectorXd RHS = -(*B)*(sigmas) - (*C)*(wake2Doublets);
//    Eigen::VectorXd doubletStrengths(bPanels->size());
//    
//    
//    Eigen::BiCGSTAB<Eigen::MatrixXd> res;
//    res.compute((*A));
//    doubletStrengths = res.solve(RHS);
//    if (res.error() > pow(10,-10))
//    {
//        converged = false;
//    }
//    
//    for (int i=0; i<bPanels->size(); i++)
//    {
//        (*bPanels)[i]->setMu(doubletStrengths(i));
//        (*bPanels)[i]->setPotential(VinfPlusVecPot((*bPanels)[i]->getCenter()));
//        
//    }
//    
//    // VinfPlusVecPot
//    for (int i=0; i<wPanels->size(); i++)
//    {
//        (*wPanels)[i]->setPrevStrength((*wPanels)[i]->getMu());
//        (*wPanels)[i]->setMu();
//        
//        (*wPanels)[i]->setPotential(VinfPlusVecPot((*wPanels)[i]->getCenter()));
//        
//        (*w2panels)[i]->setPotential(VinfPlusVecPot((*w2panels)[i]->getCenter())); // Included in this loop because there are the same number of w2pans as w1pans
//        
//    }
//    return converged;
//}
//bool cpCaseVP::solutionConvergence(){
//    
//    if(unsteady | (timeStep < 5) | params->stepsSetMaunally){
//        return false;
//    }
//    
//    // Check for simulation convergence. Will use the Trefftz plane because it is so cheap. Medium size mesh takes <1% of time step time.
//    
//    // Store previous time step values
//    double prevCD = CD_trefftz;
//    double prevCL = CL_trefftz;
//    
//    trefftzPlaneAnalysisVP();
//    
//    double changeCD = std::abs((CD_trefftz - prevCD)/prevCD);
//    double changeCL = std::abs((CL_trefftz - prevCL)/prevCL);
//    
//    double val = 0.0001; // 0.01%
//    
//    if((changeCL < val) && (changeCD < val)){
//        return true;
//    }
//    
//    return false;
//    
//}
//
//
//
//
//
//Eigen::Vector3d cpCaseVP::VinfPlusVecPot(Eigen::Vector3d POI)
//{
//    Eigen::Vector3d vInfluence = Vinf(POI);
//    
//    // Particles
//    if(accelerate && timeStep > 2)
//    {
//        vInfluence += FMM.barnesHutVel(POI);
//    }else
//    {
//        for (int i=0; i<particles.size(); i++) {
//            vInfluence += particles[i]->partVelInflGaussian(POI);
//        }
//    }
//    
//    // Filaments
//    for (int i=0; i<filaments.size(); i++)
//    {
//        vInfluence += filaments[i]->velInfl(POI);
//    }
//    
//    return vInfluence;
//}
//
//void cpCaseVP::compVelocity(){
//    //  Velocity Survey with known doublet and source strengths
//    CM.setZero();
//    Eigen::Vector3d moment;
//    Fbody = Eigen::Vector3d::Zero();
//    bodyPanel* p;
//    
//    for (int i=0; i<bPanels->size(); i++)
//    {
//        p = (*bPanels)[i];
//        p->computeVelocity(PG,Vinf(p->getCenter()));
//        if (unsteady) {
//            p->computeCp( Vinf((*bPanels)[i]->getCenter()).norm() , dt ); // Katz 13.169 ///
//        }else{
//            p->computeCp( Vmag );
//        }
//        Fbody += -p->getCp()*p->getArea()*p->getBezNormal()/params->Sref;
//        moment = p->computeMoments(params->cg);
//        CM(0) += moment(0)/(params->Sref*params->bref);
//        CM(1) += moment(1)/(params->Sref*params->cref);
//        CM(2) += moment(2)/(params->Sref*params->bref);
//    }
//    Fwind = bodyToWind(Fbody);
//    
////    if(unsteady){
////        if(timestep > 3) trefftzPlaneAnalysisVP(); // No convergence check for unsteady sims
////        
////        //   flow_time | CL_tr | CDi_tr | CN | CA | CY | CL | CD | CY | Cm | Cl | Cn |
////        solnMat(timestep-1,0) = timestep*dt;
////        solnMat(timestep-1,1) = CL_trefftz;
////        solnMat(timestep-1,2) = CD_trefftz;
////        
////        solnMat(timestep-1,3) = Fbody.x();
////        solnMat(timestep-1,4) = Fbody.y();
////        solnMat(timestep-1,5) = Fbody.z();
////        
////        solnMat(timestep-1,6) = Fwind.x();
////        solnMat(timestep-1,7) = Fwind.y();
////        solnMat(timestep-1,8) = Fwind.z();
////        
////        solnMat(timestep-1,9)  = CM.x();
////        solnMat(timestep-1,10) = CM.y();
////        solnMat(timestep-1,11) = CM.z();
////    }
//    
//    
//}
//
////void cpCaseVP::trefftzPlaneAnalysis()
////{
////    std::vector<wake*> wakes = geom->getWakes();
////    CL_trefftz = 0;
////    CD_trefftz = 0;
////    double CD_trefftzVel = 0;
////    
////    
////    for (int i=0; i<wakes.size(); i++)
////    {
////        if(vortPartFlag){
////            wakes[i]->trefftzPlaneVP(Vmag,params->Sref, &particles, timeStep);
////        } else{
////            wakes[i]->trefftzPlane(Vmag,params->Sref);
////        }
////        CL_trefftz += wakes[i]->getCL()/PG;
////        CD_trefftz += wakes[i]->getCD()/pow(PG,2);
////    }
////    double CD_T = CD_trefftz;
////}
//
//void cpCaseVP::trefftzPlaneAnalysisVP(){
//    
//    std::vector<wake*> wakes = geom->getWakes();
//    CL_trefftz = 0;
//    CD_trefftz = 0;
//    for (int i=0; i<wakes.size(); i++)
//    {
//        wakes[i]->trefftzPlaneVP(Vmag,params->Sref, &particles, timeStep);
//        CL_trefftz += wakes[i]->getCL()/PG;
//        CD_trefftz += wakes[i]->getCD()/pow(PG,2);
//    }
//}
//
//
//void cpCaseVP::stabilityDerivatives()
//{
//    double delta = 0.5;
//    double dRad = delta*M_PI/180;
//    cpCase dA(geom,Vmag,alpha+delta,beta,mach,params);
//    cpCase dB(geom,Vmag,alpha,beta+delta,mach,params);
//    
//    dA.run(false,false,false,false);
//    dB.run(false,false,false,false);
//    
//    Eigen::Vector3d FA = dA.getWindForces();
//    FA(2) = dA.getCL();
//    FA(0) = dA.getCD();
//    
//    Eigen::Vector3d FB = dB.getWindForces();
//    FB(2) = dB.getCL();
//    FB(0) = dB.getCD();
//    
//    Eigen::Vector3d F = Fwind;
//    F(2) = CL_trefftz;
//    F(0) = CD_trefftz;
//    
//    // Finish calculating stability derivatives
//    dF_dAlpha = (FA-F)/dRad;
//    dF_dBeta = (FB-F)/dRad;
//    
//    dM_dAlpha = (dA.getMoment()-CM)/dRad;
//    dM_dBeta = (dB.getMoment()-CM)/dRad;
//    
//    
//}
//
//void cpCaseVP::writeFiles()
//{
//    std::stringstream caseLabel;
//    caseLabel << "/V" << Vmag << "_Mach" << mach << "_alpha" << alpha << "_beta" << beta;
//    boost::filesystem::path subdir = boost::filesystem::current_path().string()+caseLabel.str();
//    if (!boost::filesystem::exists(subdir))
//    {
//        boost::filesystem::create_directories(subdir);
//    }
//    Eigen::MatrixXd nodeMat = geom->getNodePnts();
//    writeBodyData(subdir,nodeMat);
//    if (geom->getWakes().size() > 0)
//    {
//        writeWakeData(subdir,nodeMat);
//        writeBuffWake2Data(subdir,nodeMat);
//        writeSpanwiseData(subdir);
//    }
//    
//    if (false) {
//        createVolMesh();
//        writeVolMeshData(subdir, pts, cells);
//    }
//    
//    if(vortPartFlag){
//        //        if(timeStep > 0)
//        {
//            writeParticleData(subdir);
//            writeFilamentData(subdir);
//        }
//    }
//    
//    if (params->surfStreamFlag)
//    {
//        writeBodyStreamlines(subdir);
//    }
//    
//}
//
//void cpCaseVP::writeBodyData(boost::filesystem::path path,const Eigen::MatrixXd &nodeMat)
//{
//    std::vector<cellDataArray> data;
//    cellDataArray mu("Doublet Strengths"),sigma("Source Strengths"),pot("Velocity Potential"),V("Velocity"),Cp("Cp"),bN("bezNormals"),x("xPosition"),y("yPosition"),z("zPostition");
//    Eigen::MatrixXi con(bPanels->size(),3);
//    mu.data.resize(bPanels->size(),1);
//    sigma.data.resize(bPanels->size(),1);
//    pot.data.resize(bPanels->size(),1);
//    V.data.resize(bPanels->size(),3);
//    Cp.data.resize(bPanels->size(),1);
//    bN.data.resize(bPanels->size(),3);
//    x.data.resize(bPanels->size(),1);
//    y.data.resize(bPanels->size(),1);
//    z.data.resize(bPanels->size(),1);
//    
//    for (int i=0; i<bPanels->size(); i++)
//    {
//        mu.data(i,0) = (*bPanels)[i]->getMu();
//        sigma.data(i,0) = (*bPanels)[i]->getSigma();
//        pot.data(i,0) = (*bPanels)[i]->getPotential();
//        V.data.row(i) = (*bPanels)[i]->getGlobalV();
//        Cp.data(i,0) = (*bPanels)[i]->getCp();
//        con.row(i) = (*bPanels)[i]->getVerts();
//        bN.data.row(i) = (*bPanels)[i]->getBezNormal();
//        x.data(i,0) = (*bPanels)[i]->getCenter().x();
//        y.data(i,0) = (*bPanels)[i]->getCenter().y();
//        z.data(i,0) = (*bPanels)[i]->getCenter().z();
//    }
//    
//    data.push_back(mu);
//    data.push_back(sigma);
//    data.push_back(pot);
//    data.push_back(V);
//    data.push_back(Cp);
//    data.push_back(bN);
//    data.push_back(x);
//    data.push_back(y);
//    data.push_back(z);
//    
//    piece body;
//    body.pnts = nodeMat;
//    body.connectivity = con;
//    body.cellData = data;
//    
//    std::string fname = path.string()+"/surfaceData-" + std::to_string(timeStep)+".vtu";
//    VTUfile bodyFile(fname,body);
//}
//
//void cpCaseVP::writeWakeData(boost::filesystem::path path, const Eigen::MatrixXd &nodeMat)
//{
//    std::vector<cellDataArray> data;
//    cellDataArray mu("Doublet Strengths"),pot("Velocity Potential");
//    Eigen::MatrixXi con;
//    if(vortPartFlag){ //VPP
//        con.resize(wPanels->size(),4);
//    }else{
//        con.resize(wPanels->size(),3);
//    }
//    mu.data.resize(wPanels->size(),1);
//    pot.data.resize(wPanels->size(),1);
//    for (int i=0; i<wPanels->size(); i++)
//    {
//        mu.data(i,0) = (*wPanels)[i]->getMu();
//        pot.data(i,0) = (*wPanels)[i]->getPotential();
//        con.row(i) = (*wPanels)[i]->getVerts();
//    }
//    data.push_back(mu);
//    data.push_back(pot);
//    
//    piece wake;
//    wake.pnts = nodeMat;
//    wake.connectivity = con;
//    wake.cellData = data;
//    std::string fname = path.string()+"/wakeData-"+std::to_string(timeStep)+".vtu";
//    VTUfile wakeFile(fname,wake);
//}
//
//void cpCaseVP::writeBuffWake2Data(boost::filesystem::path path, const Eigen::MatrixXd &nodeMat)
//{
//    std::vector<cellDataArray> data;
//    cellDataArray mu("Doublet Strengths"),pot("Velocity Potential");
//    Eigen::MatrixXi con;
//    con.resize(w2panels->size(),4);
//    mu.data.resize(w2panels->size(),1);
//    pot.data.resize(w2panels->size(),1);
//    for (int i=0; i<w2panels->size(); i++)
//    {
//        mu.data(i,0) = (*w2panels)[i]->getMu();
//        pot.data(i,0) = (*w2panels)[i]->getPotential();
//        con.row(i) = (*w2panels)[i]->getVerts();
//    }
//    data.push_back(mu);
//    data.push_back(pot);
//    
//    piece wake;
//    wake.pnts = nodeMat;
//    wake.connectivity = con;
//    wake.cellData = data;
//    std::string fname = path.string()+"/bufferWake2Data-"+std::to_string(timeStep)+".vtu";
//    VTUfile wakeFile(fname,wake);
//}
//
//void cpCaseVP::writeFilamentData(boost::filesystem::path path)
//{
//    std::vector<cellDataArray> data;
//    cellDataArray mu("Gamma");
//    Eigen::MatrixXi con;
//    
//    Eigen::MatrixXd nodeMat(2*filaments.size(),3);
//    for(int i=0; i<filaments.size(); i++){
//        nodeMat.row(2*i) = filaments[i]->getP1();
//        nodeMat.row(2*i+1) = filaments[i]->getP2();
//    }
//    
//    
//    con.resize(filaments.size(),2); // 2 pnts per filament
//    
//    mu.data.resize(filaments.size(),1);
//    
//    int nodeCounter = 0; // Used to make filament end points
//    for (int i=0; i<filaments.size(); i++)
//    {
//        mu.data(i,0) = filaments[i]->getStrength();
//        con(i,0) = nodeCounter;
//        con(i,1) = nodeCounter+1;
//        nodeCounter = nodeCounter+2;
//    }
//    
//    data.push_back(mu);
//    
//    piece fils;
//    fils.pnts = nodeMat;
//    fils.connectivity = con;
//    fils.cellData = data;
//    std::string fname = path.string()+"/filaments-"+std::to_string(timeStep)+".vtu";
//    VTUfile filFile(fname,fils);
//}
//
//void cpCaseVP::writeParticleData(boost::filesystem::path path)
//{
//    
//    Eigen::MatrixXd partMat(particles.size(),3);
//    for (int i=0; i<particles.size(); i++)
//    {
//        partMat.row(i) = particles[i]->pos;
//    }
//    
//    std::vector<cellDataArray> data;
//    cellDataArray strength("Strength"), shedTime("Time Step Shed");
//    Eigen::MatrixXi con(particles.size(),1);
//    Eigen::MatrixXi shed(particles.size(),1);
//    
//    strength.data.resize(particles.size(),3);
//    shedTime.data.resize(particles.size(),1);
//    for (int i=0; i<particles.size(); i++)
//    {
//        strength.data.row(i) = particles[i]->strength;
//        shedTime.data(i,0) = particles[i]->shedTime;
//        con(i) = i;
//    }
//    data.push_back(strength);
//    data.push_back(shedTime);
//    
//    piece parts;
//    parts.pnts = partMat;
//    parts.connectivity = con;
//    parts.cellData = data;
//    
//    std::string fname = path.string()+"/particleData-"+std::to_string(timeStep)+".vtu";
//    VTUfile partFile(fname,parts);
//}
//
//
//void cpCaseVP::writeSpanwiseData(boost::filesystem::path path)
//{
//    std::vector<wake*> wakes = geom->getWakes();
//    for (int i=0; i<wakes.size(); i++)
//    {
//        Eigen::VectorXd spanLoc,Cl,Cd;
//        spanLoc = 2*wakes[i]->getSpanwisePnts()/params->bref;
//        Cl = wakes[i]->getSpanwiseCl()/PG;
//        Cd = wakes[i]->getSpanwiseCd()/pow(PG,2);
//        
//        std::stringstream ss;
//        ss << path.string() << "/spanwiseData_Wake" << i+1 << ".csv";
//        std::string fname = ss.str();
//        std::ofstream fout;
//        fout.open(fname);
//        if (fout)
//        {
//            fout << "2y/b,Cl,Cdi" << std::endl;
//            for (int i=0; i<spanLoc.size(); i++)
//            {
//                fout << spanLoc(i) << "," << Cl(i) << "," << Cd(i) << std::endl;
//            }
//        }
//        fout.close();
//    }
//}
//
//void cpCaseVP::writeBodyStreamlines(boost::filesystem::path path)
//{
//    piece p;
//    std::vector<piece> pieces;
//    std::vector<pntDataArray> data;
//    pntDataArray vel("Velocity");
//    Eigen::MatrixXi con;
//    Eigen::MatrixXd pntMat;
//    std::vector<Eigen::Vector3d> pnts,velocities;
//    
//    for (int i=0; i<bStreamlines.size(); i++)
//    {
//        pnts = bStreamlines[i]->getPnts();
//        velocities = bStreamlines[i]->getVelocities();
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
//    std::string fname = path.string()+"/streamlines.vtu";
//    VTUfile wakeFile(fname,pieces);
//}
//
//
//void cpCaseVP::writeVolMeshData(boost::filesystem::path path, Eigen::MatrixXd &nodeMat, std::vector<Eigen::VectorXi> cells)
//{
//    int nCells = (int)cells.size();
//    
//    std::vector<cellDataArray> data;
//    cellDataArray vel("Velocity"),vortMag("Vorticity Magnitude"), Cp("Cp");
//    Eigen::MatrixXi con;
//    con.resize(nCells,8);
//    vel.data.resize(nCells,3);
//    vortMag.data.resize(nCells,1);
//    Cp.data.resize(nCells,1);
//    for (int i=0; i<nCells; i++)
//    {
//        vel.data.row(i) = volMeshDat.velocity[i];
//        vortMag.data(i,0) = volMeshDat.vorticity[i];
//        Cp.data(i,0) = volMeshDat.pressure[i];
//        con.row(i) = cells[i];
//    }
//    data.push_back(vel);
//    data.push_back(vortMag);
//    data.push_back(Cp);
//    
//    piece volMesh;
//    volMesh.pnts = nodeMat;
//    volMesh.connectivity = con;
//    volMesh.cellData = data;
//    std::string fname = path.string()+"/volumeMesh-"+std::to_string(timeStep)+".vtu";
//    VTUfile wakeFile(fname,volMesh);
//}
//
//void cpCaseVP::convectBufferWake(){
//    wake2Doublets.resize((*w2panels).size());
//    for (int i=0; i<(*w2panels).size(); i++)
//    {
//        (*w2panels)[i]->setPrevStrength((*w2panels)[i]->getMu());
//        double parentMu = (*w2panels)[i]->getBufferParent()->getMu();
//        (*w2panels)[i]->setMu(parentMu);
//        wake2Doublets[i] = parentMu;
//    }
//    
//}
//
//
//void cpCaseVP::collapseBufferWake(){
//    
//    std::vector<edge*> usedEdges;
//    
//    for(int i=0; i<(*w2panels).size(); i++)
//    {
//        std::vector<edge*> pEdges = (*w2panels)[i]->edgesInOrder();
//        
//        Eigen::Vector3d strength = Eigen::Vector3d::Zero();
//        for (int j=1; j<4; j++)
//        {
//            if (!edgeIsUsed(pEdges[j],usedEdges))
//            {
//                usedEdges.push_back(pEdges[j]);
//                strength += (*w2panels)[i]->edgeStrength( pEdges[j], j ); // Don't need to pass in pEdges...
//            }
//        }
//        
//        Eigen::Vector3d pos = rungeKuttaStepper((*w2panels)[i]->getCenter());
//        //        Eigen::Vector3d pos = (*w2panels)[i]->getCenter(); // For moving geometry
//        Eigen::Vector3d ptVel = Vinf(pos);
//        double radius = (*w2panels)[i]->getPartRadius(ptVel,dt); // VinfLocal
//        
//        particle* p = new particle(pos, strength, radius, {0,0,0}, {0,0,0}, timeStep); // Zeroes are previous pos and strength values used for Adams-Bashforth scheme
//        p->parentPanel = (*w2panels)[i];
//        particles.push_back(p);
//        
//    }
//    
//    // Create filament
//    if(filaments.size() == 0)
//    {
//        for(int i=0; i<(*w2panels).size(); i++)
//        {
//            vortexFil* fil;
//            Eigen::Vector3d p1,p2;
//            
//            p1 = (*w2panels)[i]->pointsInOrder()[2]->getPnt();
//            p2 = (*w2panels)[i]->pointsInOrder()[3]->getPnt();
//            
//            fil = new vortexFil(p1, p2,-(*w2panels)[i]->getMu(), (*w2panels)[i]); // Negative strength is because filament is actually the upstream edge being convected which is oriented the opposite direction as downstream edge
//            
//            filaments.push_back(fil);
//            (*w2panels)[i]->setVortFil(fil);
//        }
//    }
//    else
//    {
//        for(int i=0; i<(*w2panels).size(); i++){
//            filaments[i]->setStrength(-(*w2panels)[i]->getMu()); // Negative strength is because filament is actually the upstream edge being convected which is oriented the opposite direction as downstream edge
//        }
//        
//    }
//}
//
//
//
//
//bool cpCaseVP::edgeIsUsed(edge* thisEdge, std::vector<edge*> pEdges){
//    
//    for(int i=0; i<pEdges.size(); i++){
//        if(thisEdge == pEdges[i]){
//            return true;
//        }
//    }
//    return false;
//}
//
//
//
//Eigen::Vector3d cpCaseVP::rungeKuttaStepper(Eigen::Vector3d POI){
//    
//    // RK4 algorithm makes symmeterized influences VERY impractical
//    
//    Eigen::Vector3d k1, k2, k3, k4;
//    
//    k1 = velocityInflFromEverything(POI);
//    k2 = velocityInflFromEverything(POI+k1*dt/2);
//    k3 = velocityInflFromEverything(POI+k2*dt/2);
//    k4 = velocityInflFromEverything(POI+k3*dt);
//    
//    return POI + dt*(k1/6 + k2/3 + k3/3 + k4/6);
//    
//}
//
//
//Eigen::Vector3d cpCaseVP::velocityInflFromEverything( Eigen::Vector3d POI ){
//    // Freestream influence
//    Eigen::Vector3d velOnPart = Vinf(POI);
//    
//    // Particle influence
//    if(accelerate && timeStep > 3){
//        velOnPart += FMM.barnesHutVel(POI);
//    }
//    else{
//        for(int j=0;j<particles.size();j++)
//        {
//            velOnPart += particles[j]->partVelInflGaussian(POI);
//        }
//    }
//    
//    
//    // Body panel influence
//    for(int j=0;j<(*bPanels).size();j++){
//        velOnPart += (*bPanels)[j]->panelV(POI);
//    }
//    
//    
//    // Buffer wake influence
//    for(int j=0;j<(*wPanels).size();j++){
//        velOnPart += (*wPanels)[j]->panelV(POI);
//    }
//    
//    for(int j=0;j<(*w2panels).size();j++){
//        velOnPart += (*w2panels)[j]->panelV(POI);
//    }
//    
//    // Vortex Filament influence
//    for(int i=0; i<filaments.size(); i++){
//        velOnPart +=filaments[i]->velInfl(POI);
//    }
//    
//    //    std::cout << velOnPart.x() << " , " << velOnPart.y() << " , " << velOnPart.z() << std::endl;
//    
//    return velOnPart;
//}
//
//Eigen::Vector3d cpCaseVP::velocityInflFromEverything(particle* part){
//    
//    // Function is overloaded for symmeterized velocity influence for particles
//    Eigen::Vector3d pos = part->pos;
//    
//    // Freestream influence
//    Eigen::Vector3d velOnPart = Vinf(pos);
//    
//    // If particle reaches very far field condition
//    if ((pos-Eigen::Vector3d::Zero()).norm() > 20*params->cref)
//    {
//        return velOnPart;
//    }
//    
//    // Particle influence
//    if(accelerate && timeStep > 2){
//        velOnPart += FMM.barnesHutVel(part);
//    }
//    else{
//        for(int j=0;j<particles.size();j++)
//        {
//            velOnPart += particles[j]->partVelInflGaussian(pos);
//        }
//    }
//    
//    // If particle reaches first farfield condition
//    if((pos-Eigen::Vector3d::Zero()).norm() > 12*params->cref)
//    {
//        return velOnPart;
//    }
//    
//    
//    // Body panel influence
//    for(int j=0;j<(*bPanels).size();j++){
//        velOnPart += (*bPanels)[j]->panelV(pos);
//    }
//    
//    
//    // Buffer wake influence
//    for(int j=0;j<(*wPanels).size();j++){
//        velOnPart += (*wPanels)[j]->panelV(pos);
//    }
//    for(int j=0;j<(*w2panels).size();j++){
//        velOnPart += (*w2panels)[j]->panelV(pos);
//    }
//    
//    
//    // Vortex Filament influence
//    for(int i=0; i<filaments.size(); i++){
//        velOnPart +=filaments[i]->velInfl(pos);
//    }
//    
//    return velOnPart;
//}
//
//void cpCaseVP::convectParticles(){
//    std::vector<Eigen::Vector3d> newPartPositions;
//    
//    for(int i=0;i<particles.size();i++){
//        
//        Eigen::Vector3d newPos;
//        
//        // Adams bashforth scheme
//        Eigen::Vector3d velOnPart = velocityInflFromEverything(particles[i]);
//        
//        if(particles[i]->getPrevVelInfl().isZero())
//        {
//            newPos = particles[i]->pos + dt*velOnPart;
//        }
//        else{
//            newPos = particles[i]->pos + dt*(1.5*velOnPart - 0.5*particles[i]->getPrevVelInfl());
//        }
//        
//        particles[i]->setPrevVelInfl(velOnPart);
//        newPartPositions.push_back(newPos);
//        //        }
//    }
//    
//    for(int i=0;i<particles.size();i++){
//        particles[i]->setPos(newPartPositions[i]);
//        
//    }
//}
//
//
//
//void cpCaseVP::particleStrengthUpdate(){
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
//
//void cpCaseVP::particleStrengthUpdateGaussian(){
//    // This function uses the update equations found in 'Vortex Methods for DNS of...' by Plouhmhans. It uses the particle strength exchange for the viscous diffusion
//    
//    // Can use this strength exchange function and call it from the normal stretching function instead of the winklemans function?
//    
//    std::vector<Eigen::Vector3d> stretchDiffVec; // Creating vector values because the strength change needs to be set after all particle influences have been calculated
//    for(int i=0; i<particles.size(); i++)
//    {
//        Eigen::Vector3d dAlpha_diff = Eigen::Vector3d::Zero();
//        Eigen::Vector3d dAlpha_stretch = Eigen::Vector3d::Zero();
//        
//        // FarField2 condition
//        if((particles[i]->pos-Eigen::Vector3d::Zero()).norm() > 20*params->cref)
//        {
//            stretchDiffVec.push_back(Eigen::Vector3d::Zero());
//        }
//        // FarField1 condition
//        else if((particles[i]->pos-Eigen::Vector3d::Zero()).norm() > 12*params->cref)
//        {
//            if(unsteady)
//            {
//                dAlpha_stretch = FMM.barnesHutStretch(particles[i]);
//                dAlpha_diff = FMM.barnesHutDiff(particles[i]);
//            }
//            else
//            {
//                // Stretching from particles
//                for(int j=0; j<particles.size(); j++)
//                {
//                    if(i!=j) // Kroneger Delta Function
//                    {
//                        dAlpha_diff += particles[i]->viscousDiffusionGaussian(particles[j]);
//                        dAlpha_stretch += particles[i]->vortexStretchingGaussian(particles[j]);
//                    }
//                }
//            }
//            stretchDiffVec.push_back( dAlpha_diff + dAlpha_stretch );
//        }
//        else // Treat normally
//        {
//            // Influence from Particles
//            if(accelerate && timeStep > 2){
//                dAlpha_stretch += FMM.barnesHutStretch(particles[i]);
//                dAlpha_diff += FMM.barnesHutDiff(particles[i]);
//            }else
//            {
//                for(int j=0; j<particles.size(); j++)
//                {
//                    if(i!=j){ // Kroneger Delta Function
//                        dAlpha_diff += particles[i]->viscousDiffusionGaussian(particles[j]);
//                        dAlpha_stretch += particles[i]->vortexStretchingGaussian(particles[j]);
//                    }
//                }
//            }
//            
//            // Stretching from body panels
//            for(int j=0; j<(*bPanels).size(); j++)
//            {
//                dAlpha_diff += (*bPanels)[j]->partStretching(particles[i]);
//            }
//            
//            // Stretching from wake panels
//            for(int j=0;j<(*wPanels).size(); j++)
//            {
//                dAlpha_diff += (*wPanels)[j]->partStretching(particles[i]);
//            }
//            
//            stretchDiffVec.push_back( dAlpha_diff + dAlpha_stretch );
//        }
//    }
//    for(int i=0;i<particles.size();i++)
//    {
//        Eigen::Vector3d newStrength;
//        if(particles[i]->getprevStrengthUpdate().isZero())
//        {
//            newStrength = particles[i]->strength + stretchDiffVec[i]*dt;
//        }else
//        {   // Adams bashforth
//            newStrength = particles[i]->strength + dt*(1.5*stretchDiffVec[i] - 0.5*particles[i]->getprevStrengthUpdate());
//        }
//        particles[i]->setprevStrengthUpdate(stretchDiffVec[i]);
//        particles[i]->setStrength(newStrength);
//    }
//}
//
//
//
//void cpCaseVP::readBodyKinFile(){
//    std::ifstream fid;
//    fid.open(params->bodyKinFileLoc);
//    
//    bodyKin.resize(numSteps, 6); // U, V, W, p, q, r
//    
//    for (int i=0; i<numSteps; i++)
//    {
//        for(int j=0; j < 6; j++)
//        {
//            fid >> bodyKin(i,j);
//        }
//    }
//    
//    fid.close();
//}
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//
//void cpCaseVP::createVolMesh(){
//    
//    // Clear Past Timestep Mesh
//    cells.clear();
//    volMeshDat.cellCenter.clear();
//    volMeshDat.velocity.clear();
//    volMeshDat.vorticity.clear();
//    
//    // limits of mesh
//    double x0, xf, y0, yf, z0, zf;
//    x0 = -0.5;
//    xf = 6.5;
//    y0 = -8;
//    yf = 8;
//    z0 = -0.5;
//    zf = 1.25;
//    
//    // resolution
//    int nX, nY, nZ;
//    nX = 500;
//    nY = 1;
//    nZ = 105;
//    int nCells = nX * nY * nZ;
//    
//    double hx, hy, hz;
//    hx = (xf-x0)/nX;
//    hy = (yf-y0)/nY;
//    hz = (zf-z0)/nZ;
//    
//    // Add one for number of points
//    int nXp = nX+1;
//    int nYp = nY+1;
//    int nZp = nZ+1;
//    int numPts = nXp * nYp * nZp;
//    
//    
//    
//    // creating pnt vector
//    pts.resize(numPts, 3);
//    int count = 0;
//    for (int k=0; k<nZp; k++) {
//        for (int j=0; j<nYp; j++) {
//            for (int i=0; i<nXp; i++) {
//                pts(count,0) = x0 + hx*i;
//                pts(count,1) = y0 + hy*j;
//                pts(count,2) = z0 + hz*k;
//                count++;
//            }
//        }
//    }
//    
//    
//    // Create cells
//    for (int i=0; i<nX; i++) {
//        for (int j=0; j<nY; j++) {
//            for (int k=0; k<nZ; k++) {
//                Eigen::VectorXi cell_pts(8);
//                cell_pts[0] = (k*(nXp*nYp) + j*(nXp) + i);
//                cell_pts[1] = (k*(nXp*nYp) + j*(nXp) + i) + 1;
//                
//                cell_pts[2] = (k*(nXp*nYp) + (j+1)*(nXp) + i);
//                cell_pts[3] = (k*(nXp*nYp) + (j+1)*(nXp) + i) + 1;
//                
//                cell_pts[4] = ((k+1)*(nXp*nYp) + j*(nXp) + i);
//                cell_pts[5] = ((k+1)*(nXp*nYp) + j*(nXp) + i) + 1;
//                
//                cell_pts[6] = ((k+1)*(nXp*nYp) + (j+1)*(nXp) + i);
//                cell_pts[7] = ((k+1)*(nXp*nYp) + (j+1)*(nXp) + i) + 1;
//                
//                cells.push_back(cell_pts);
//                
//                // Find cell center
//                Eigen::MatrixXd cellCorners(8,3);
//                for (int m=0; m<cellCorners.rows(); m++) {
//                    cellCorners.row(m) = pts.row(cell_pts[m]);
//                }
//                
//                Eigen::Vector3d center;
//                for (int m=0; m<cellCorners.cols(); m++) {
//                    center(m) = cellCorners.col(m).mean();
//                }
//                volMeshDat.cellCenter.push_back(center);
//            }
//        }
//    }
//    
//    for (int i=0; i<nCells; i++) {
//        // Find cell center
//        Eigen::Vector3d velInCell = velocityInflFromEverything(volMeshDat.cellCenter[i]);
//        volMeshDat.velocity.push_back(velInCell);
//    }
//    
//    for (int i=0; i<nCells; i++) {
//        // Incompressible Bernoulli equation
//        double Cp = pow( volMeshDat.velocity[i].norm() / Vmag , 2 );
//        volMeshDat.pressure.push_back(Cp);
//    }
//    
//    for (int i=0; i<nCells; i++) {
//        volMeshDat.vorticity.push_back(volMeshDat.cellCenter[i].norm());
//    }
//    
//    
//}
//
//
//
//
