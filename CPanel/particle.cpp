//
//  particle.cpp
//  CPanel - Unstructured Panel Code
//
//  Created by Connor Sousa on 3/5/16.
//  Copyright (c) 2016 Chris Satterwhite. All rights reserved.
//

#include "particle.h"
#include <Eigen/Geometry>
#include "math.h"


particle::particle(Eigen::Vector3d pos, Eigen::Vector3d strength, double radius, Eigen::Vector3d previousVelInfl, Eigen::Vector3d previousStrengthUpdate, int shedTime) : pos(pos), strength(strength), radius(radius), previousVelInfl(previousVelInfl), previousStrengthUpdate(previousStrengthUpdate), shedTime(shedTime) {};



Eigen::Vector3d particle::partVelInfl(particle* part){
    // Velocity influence with Gaussian smoothing
    
    if(this == part){return Eigen::Vector3d::Zero();} // Particle doesn't influence itself
    
    // 'part' is the influenced particle
    double sigma = std::sqrt(pow(coreOverlap*this->radius,2) + pow(coreOverlap*part->radius,2))/std::sqrt(2);
    
    Eigen::Vector3d dist = part->pos-this->pos;
    double rho = dist.norm()/sigma;
    
    double K = (1/(4*M_PI*rho)*erf(rho/pow(2,0.5))-1/(pow(2*M_PI,1.5))*exp(-0.5*rho*rho))/(rho*rho);
    
    return -1/pow(sigma,3)*K*dist.cross(strength);
};

Eigen::Vector3d particle::partVelInfl(const Eigen::Vector3d &POI){
    // This overloaded function is to calculate the velocity at a point, instead of on a particle. The only difference is that the smoothing radius is no longer symmeterized
    
    if((this->pos - POI).norm() == 0){return Eigen::Vector3d::Zero();}
    
    double sigma = coreOverlap*radius;
    
    Eigen::Vector3d dist = POI-pos;
    double rho = dist.norm()/sigma;
    
    double K = (1/(4*M_PI*rho)*erf(rho/pow(2,0.5))-1/(pow(2*M_PI,1.5))*exp(-0.5*rho*rho))/(rho*rho);
    
    return -1/pow(sigma,3)*K*dist.cross(strength);
    
    
};


Eigen::Vector3d particle::vortexStretching(particle* part)
{
    
    if(this == part){return Eigen::Vector3d::Zero();} // Particle doesn't influence itself
    
    Eigen::Vector3d Xi = this->pos;
    Eigen::Vector3d Xj = part->pos;
    
    if ((Xi-Xj).norm() > 5*this->radius)
    {
        return Eigen::Vector3d::Zero();
    }
    
    double sigij = std::sqrt(pow(coreOverlap*this->radius,2) + pow(coreOverlap*part->radius,2))/std::sqrt(2);
    
    double rho = (Xi-Xj).norm()/sigij;
    double xi = 1/(pow(2*M_PI,1.5))*exp(-rho*rho/2);
    
    double G = 1/(4*M_PI*rho)*erf(rho/std::sqrt(2));
    double K = (G-xi)/(rho*rho);
    double F = (3*K-xi)/(rho*rho);
    
    Eigen::Matrix3d inflMat;
    
    for(int k=0; k<3; k++){ // Building matrix found in He and Zhao Eq. 15
        for(int l=0; l<3; l++){
            if(k==l){
                inflMat(k,l) = K - 1/(sigij*sigij)*F*(Xi(k)-Xj(k))*(Xi(l)-Xj(l));
            }else{
                inflMat(k,l) =    -1/(sigij*sigij)*F*(Xi(k)-Xj(k))*(Xi(l)-Xj(l));
            }
        }
    }
    
    Eigen::Matrix3d alphaTilda;
    Eigen::Vector3d a = part->strength;
    alphaTilda << 0,-a.z(),a.y(),a.z(),0,-a.x(),-a.y(),a.x(),0; // cross product matrix
    
    Eigen::Matrix3d velGradient = 1/pow(sigij,3)*alphaTilda*inflMat;
    
    return velGradient*this->strength;
}



Eigen::Vector3d particle::viscousDiffusion(particle* part)
{
    // From Ploumhans: 'Vortex Methods for 3D Bluff...'
    
    if(this == part){return Eigen::Vector3d::Zero();} // Part doesn't influence itself
    
    double nu = 1.983e-5; //Change based off of input...
    
    Eigen::Vector3d Xi = this->pos;
    Eigen::Vector3d Xj = part->pos;
    
    double sigij = std::sqrt(pow(coreOverlap*this->radius,2) + pow(coreOverlap*part->radius,2))/std::sqrt(2);
    double rho = (Xi-Xj).norm()/sigij;
    
    double Vi = 4*M_PI/3*pow(this->radius,3); //Not 100% sure this is the best way to find Volume. Could put it in constructor? Also, should it be spherical, or cubic. Maybe more info in Calebretta...
    double Vj = 4*M_PI/3*pow(part->radius,3);
    
    double xi = 1/(pow(2*M_PI,1.5)*pow(rho,3))*exp(-rho*rho/2);
    
    return 2*nu/(sigij*sigij)*xi*(Vi*part->strength-Vj*this->strength);
}




//=======================================//

// This function is much faster than using the Gaussian influence (above). However, it is not technically as stable. This influencing function could be used for almost all simulations. However, I wanted to make sure the solver was robust in case CPanel is extended to a very parallel, high resultion solver.


//Eigen::Vector3d particle::partVelInfl(const Eigen::Vector3d &POI){
//
//    double sigma = coreOverlap*radius;
//    Eigen::Vector3d dist = POI-pos;
//
//    // ** High algebraic smoothing ** //
//    return -1/(4*M_PI)*(pow(dist.norm(),2) + 2.5*pow(sigma,2))/(pow(pow(dist.norm(),2) + pow(sigma,2),2.5)) * dist.cross(strength);
//};


// This funtion is the strength update function from Winklemans which includes stretching and diffusion (diffusion through the value of nu)

//Eigen::Vector3d particle::partStrengthUpdate(particle* part){
//    // Transpose scheme from Winklemans (appendix F)
//    double sigma = coreOverlap*radius; // coreOverlap of 1.3;
//
//    double volP = 4*M_PI/3*pow(this->radius,3);
//    double volQ = 4*M_PI/3*pow(part->radius,3);
//
//    Eigen::Vector3d dist = part->pos - this->pos;
//    double nu = 1.983e-5;
//
//    Eigen::Vector3d alphaP = this->strength;
//    Eigen::Vector3d alphaQ = part->strength;
//
//
//    Eigen::Vector3d firstTerm = (dist.norm()*dist.norm() + 2.5*sigma*sigma)/(pow(dist.norm()*dist.norm() + sigma*sigma,2.5)) * alphaP.cross(alphaQ);
//
//    Eigen::Vector3d secTerm = 3*(dist.norm()*dist.norm()+3.5*sigma*sigma)/(pow(dist.norm()*dist.norm()+sigma*sigma,3.5))*(alphaP.dot(dist.cross(alphaQ)))*dist;
//
//    Eigen::Vector3d thirdTerm = 105*nu*pow(sigma,4)/(pow(dist.norm()*dist.norm()+sigma*sigma,4.5))*(volP*alphaQ - volQ*alphaP);
//
//    return -(1/(4*M_PI))*(firstTerm + secTerm + thirdTerm);
//}






