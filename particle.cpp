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


particle::particle(Eigen::Vector3d pos, Eigen::Vector3d strength, double radius) : pos(pos), strength(strength), radius(radius) {};

Eigen::Vector3d particle::partVelInfl(const Eigen::Vector3d &POI){
    
    
    double coreOverlap = 1.3; // The amount of core overlap between neighboring particles. Wincklemans used 1.3. More info: calebretta pp. 47
    
    // Research more into what's going on here
    double sigma = coreOverlap*radius;
    
    Eigen::Vector3d dist = POI-pos;
    
    // ** Regularized Influence ** //
    return -1/(4*M_PI)*(pow(dist.norm(),2) + 2.5*pow(sigma,2))/(pow(pow(dist.norm(),2) + pow(sigma,2),2.5)) * dist.cross(this->strength);
    
    // ** Singular Influence ** //
//    return -1/(4*M_PI)*(1/pow(dist.norm(),3))*dist.cross(strength);
    
};

Eigen::Vector3d particle::partVelInflGaussian(const Eigen::Vector3d &POI){
    
    double coreOverlap = 1; // The amount of core overlap between neighboring particles. Wincklemans used 1.3. More info: calebretta pp. 47
    
    double sigma = std::sqrt(pow(coreOverlap*this->radius,2) + pow(coreOverlap*this->radius,2))/std::sqrt(2);

    
    Eigen::Vector3d dist = POI-pos;
    double rho = dist.norm()/sigma;
    double eta = 1/(pow(2*M_PI,1.5)*pow(sigma,3))*exp(-0.5*rho*rho);
    
    double K = (1/(4*M_PI*rho)*erf(rho/pow(2,0.5))-1/(pow(2*M_PI,1.5))*exp(-0.5*rho*rho))/(rho*rho);
    
    return -1/pow(sigma,3)*K*dist.cross(strength);
    
};

Eigen::Vector3d particle::partStretching(particle* part){
    
//    std::cout << "Pre stretching: " << this->strength.x() << ", " << this->strength.y() << ", " << this->strength.z() << std::endl;
    double sigma = 1.3*radius; // coreOverlap of 1.3;
    
    double volP = 4*M_PI/3*pow(this->radius,3);
    double volQ = 4*M_PI/3*pow(this->radius,3);

    Eigen::Vector3d dist = part->pos - this->pos;
    double nu = 1.983e-5;
    
    Eigen::Vector3d alphaQ = part->getStrength();
    Eigen::Vector3d alphaP = this->getStrength();

    Eigen::Vector3d strengthChange = 1/(4*M_PI)*( 1.5*(dist.norm()*dist.norm() + 3.5*sigma*sigma)/(pow(dist.norm()*dist.norm() + sigma*sigma,3.5))*alphaP.dot(dist)*dist.cross(alphaQ) +
                                                 (alphaP.dot(dist.cross(alphaQ)))*dist+
                                                 
                                                 105*nu*pow(sigma,4)/(pow(dist.norm()*dist.norm()+sigma*sigma,4.5))*(volP*alphaQ - volQ*alphaP));
    
//    
//    Eigen::Vector3d firstTerm = 1.5*(dist.norm()*dist.norm() + 3.5*sigma*sigma)/(pow(dist.norm()*dist.norm() + sigma*sigma,3.5))*alphaP.dot(dist)*dist.cross(alphaQ);
//    
//    Eigen::Vector3d secTerm =  (alphaP.dot(dist.cross(alphaQ)))*dist;
//    
//    Eigen::Vector3d thirdTerm = 105*nu*pow(sigma,4)/(pow(dist.norm()*dist.norm()+sigma*sigma,4.5))*(volP*alphaQ - volQ*alphaP);

    
//    std::cout << "firstTerm: " << firstTerm.x() << ", " << firstTerm.y() << ", " << firstTerm.z() << std::endl;
//    std::cout << "secTerm: " << secTerm.x() << ", " << secTerm.y() << ", " << secTerm.z() << std::endl;
//    std::cout << "thirdTerm: " << thirdTerm.x() << ", " << thirdTerm.y() << ", " << thirdTerm.z() << "\n" <<std::endl;
    
//    std::cout << "firstTerm: " << firstTerm.norm() << std::endl;
//    std::cout << "secTerm  : " << secTerm.norm() << "\n" << std::endl;
//    std::cout << "thirdTerm: " << thirdTerm.norm() << "\n" <<std::endl;
    
    return strengthChange;
}

Eigen::Matrix3d particle::partStretchingGaussian(particle* part){
//    double nu = 1.983e-5;
    Eigen::Vector3d Xi = this->getPos();
    Eigen::Vector3d Xj = part->getPos();
    
    double coreOverlap = 1; // The amount of core overlap between neighboring particles. Wincklemans used 1.3. More info: calebretta pp. 47
    
//    std::cout << "sigma = " << this->radius;
    double sigij = std::sqrt(pow(coreOverlap*this->radius,2) + pow(coreOverlap*part->radius,2))/std::sqrt(2);
//    std::cout << ".  sigIJ = " << sigij << std::endl;
    
    double rho = (Xi-Xj).norm()/sigij;
    
    double K = (1/(4*M_PI*rho)*erf(rho/std::sqrt(2))-1/(pow(2*M_PI,1.5))*exp(-0.5*rho*rho))/(rho*rho);

    double xi = 1/(pow(2*M_PI,1.5))*exp(-rho*rho/2);

    double F = (3*K-xi)/(rho*rho);
    
    Eigen::Matrix3d inflMat;
    
    for(int k=0; k<3; k++){ // Building matrix found in He and Zhao Eq. 15
        for(int l=0; l<3; l++){
            if(k==l){
                inflMat(k,l) = K - 1/(sigij*sigij)*F*(Xi(k)-Xj(k))*(Xi(l)-Xj(l));
            }else{
                inflMat(k,l) = -1/(sigij*sigij)*F*(Xi(k)-Xj(k))*(Xi(l)-Xj(l));
            }
        }
    }
    
    Eigen::Matrix3d alphaTilda;
    Eigen::Vector3d a = part->getStrength();
    alphaTilda << 0,-a.z(),a.y(),a.z(),0,-a.x(),-a.y(),a.x(),0;
    
    Eigen::Matrix3d gradient = 1/pow(sigij,3)*alphaTilda*inflMat;
    
    return gradient;
}


