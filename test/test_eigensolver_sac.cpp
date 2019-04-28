/******************************************************************************
 * Author:   Laurent Kneip                                                    *
 * Contact:  kneip.laurent@gmail.com                                          *
 * License:  Copyright (c) 2013 Laurent Kneip, ANU. All rights reserved.      *
 *                                                                            *
 * Redistribution and use in source and binary forms, with or without         *
 * modification, are permitted provided that the following conditions         *
 * are met:                                                                   *
 * * Redistributions of source code must retain the above copyright           *
 *   notice, this list of conditions and the following disclaimer.            *
 * * Redistributions in binary form must reproduce the above copyright        *
 *   notice, this list of conditions and the following disclaimer in the      *
 *   documentation and/or other materials provided with the distribution.     *
 * * Neither the name of ANU nor the names of its contributors may be         *
 *   used to endorse or promote products derived from this software without   *
 *   specific prior written permission.                                       *
 *                                                                            *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"*
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE  *
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE *
 * ARE DISCLAIMED. IN NO EVENT SHALL ANU OR THE CONTRIBUTORS BE LIABLE        *
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL *
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR *
 * SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER *
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT         *
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY  *
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF     *
 * SUCH DAMAGE.                                                               *
 ******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <iomanip>
#include <opengv/relative_pose/methods.hpp>

#include "EigensolverSacProblem.hpp"
#include "EigensolverAdapter.hpp"
#include <opengv/sac/MultiRansac.hpp>
#include <opengv/sac/Lmeds.hpp>
#include <sstream>
#include <fstream>

#include "random_generators.hpp"
#include "experiment_helpers.hpp"
#include "time_measurement.hpp"
#include <opengv/math/cayley.hpp>
#include <limits.h>
#include <Eigen/Eigen>



using namespace std;
using namespace Eigen;
using namespace opengv;

int main(int argc, char **argv) {
  // initialize random seed
  initializeRandomSeed();

  //set experiment parameters
  double noise = 10.0;
  double outlierFraction = 0.4;
  size_t numberPoints = 200;

  //generate a random pose for viewpoint 1
  translation_t position1 = Eigen::Vector3d::Zero();
  rotation_t rotation1 = Eigen::Matrix3d::Identity();

  //generate a random pose for viewpoint 2
//  translation_t position2 = generateRandomDirectionTranslation(0.0005);
//  rotation_t rotation2 = generateRandomRotation(0.5);

  translation_t position2;
  rotation_t rotation2;
  double theta = 9 * M_PI / 180.0;
  rotation2 << cos(theta), sin(theta), 0,
          -sin(theta), cos(theta), 0,
          0, 0, 1;
  position2 << 1 - cos(theta), sin(theta), 0;

  position2 = position2;
  std::cout << "rotation2: " << rotation2 << std::endl;
  std::cout << "position2: " << position2.transpose() << std::endl;

  //create a fake central camera
  translations_t camOffsets;
  rotations_t camRotations;
//    generateRandomCameraSystem(4, camOffsets, camRotations);
  camOffsets.resize(4); camRotations.resize(4);
  camOffsets[0]<<-0.2, 0, 0;
  camOffsets[1]<< 0.2, 0, 0;
  camOffsets[2]<< 0, 0.3, 0;
  camOffsets[3]<< 0, -0.3, 0;

  camRotations[0]<<0, 0, -1, 1, 0, 0, 0, -1, 0;
  camRotations[1]<<0, 0, 1, -1, 0, 0, 0, -1, 0;
  camRotations[2]<<1, 0, 0, 0, 0, 1, 0, -1, 0;
  camRotations[3]<<-1, 0, 0, 0, 0, -1, 0, -1, 0;


  //derive correspondences based on random point-cloud
  std::vector<std::shared_ptr<bearingVectors_t> > multiBearingVectors1;
  std::vector<std::shared_ptr<bearingVectors_t> > multiBearingVectors2;
  std::vector<std::shared_ptr<Eigen::MatrixXd> > gt;

  generateMulti2D2DCorrespondences(
          position1, rotation1, position2, rotation2, camOffsets, camRotations,
          numberPoints, noise, outlierFraction,
          multiBearingVectors1, multiBearingVectors2, gt );

  //derive correspondences based on random point-cloud
/*    bearingVectors_t bearingVectors1;
    bearingVectors_t bearingVectors2;
    std::vector<int> camCorrespondences1; //unused in the central case
    std::vector<int> camCorrespondences2; //unused in the central case
    Eigen::MatrixXd gt2(3, numberPoints);
    generateRandom2D2DCorrespondences(
            position1, rotation1, position2, rotation2,
            camOffsets, camRotations, numberPoints, noise, outlierFraction,
            bearingVectors1, bearingVectors2,
            camCorrespondences1, camCorrespondences2, gt2);*/


  /*  for (int i = 0; i < multiBearingVectors1.size(); i++) {
        for (int j = 0; j < (*multiBearingVectors1[i]).size(); j++) {
            (*multiBearingVectors1[i])[j] = camRotations[i] * (*multiBearingVectors1[i])[j];
            (*multiBearingVectors2[i])[j] = camRotations[i] * (*multiBearingVectors2[i])[j];
        }
    }*/


  std::cout << "bv size: " << multiBearingVectors1.size() << std::endl;

  //Extract the relative pose
  translation_t position;
  rotation_t rotation;
  extractRelativePose(
          position1, position2, rotation1, rotation2, position, rotation);

  //create a central relative adapter
  polyview::demo::EigensolverAdapter adapter(
          multiBearingVectors1,
          multiBearingVectors2,
          camOffsets,
          camRotations,
          position,
          rotation);

  rotation_t R_perturbed;
  double theta_P = 0.5 * M_PI / 180.0;
  R_perturbed << cos(theta_P), sin(theta_P), 0,
          -sin(theta_P), cos(theta_P), 0,
          0, 0, 1;

//  R_perturbed = camRotations[0].transpose() * R_perturbed * camRotations[0];

  std::cout << R_perturbed << std::endl;

  adapter.setR12(R_perturbed);

  //print experiment characteristics
//  printExperimentCharacteristics( position, R_perturbed, noise, outlierFraction );

  /*   rotation_t eigensolver_rotation;
     eigensolverOutput_t output;
     output.rotation = adapter.getR12(); //transferring the initial value
     eigensolver_rotation = relative_pose::eigensolver(adapter, output);

     double z = opengv::math::rot2cayley(eigensolver_rotation)[2];
     double pureError = asin((2 * z) / (1 + pow(z, 2))) / -theta;
     //print results
     std::cout << "results from eigensystem based rotation solver:" << std::endl;
     std::cout << pureError << std::endl << std::endl;*/

  //Create an EigensolverSacProblem and Ransac
  //The number of samples can be configured
  opengv::sac::MultiRansac<
          polyview::demos::EigensolverSacProblem> ransac;
  std::shared_ptr<
          polyview::demos::EigensolverSacProblem> eigenproblem_ptr(
          new polyview::demos::EigensolverSacProblem (adapter, 12));
  ransac.sac_model_ = eigenproblem_ptr;
  ransac.threshold_ = 2.0*(1.0 - cos(atan(sqrt(2.0)*0.5/800.0)));
  ransac.max_iterations_ = 1000;



  //Run the experiment
  struct timeval tic;
  struct timeval toc;
  gettimeofday(&tic, 0);
  ransac.computeModel();
  gettimeofday(&toc, 0);
  double ransac_time = TIMETODOUBLE(timeval_minus(toc, tic));

  //do final polishing of the model over all inliers
  polyview::demos::EigensolverSacProblem::model_t optimizedModel;
  eigenproblem_ptr->optimizeModelCoefficients(
          ransac.inliers_,
          ransac.model_coefficients_,
          optimizedModel);

  double z2 = opengv::math::rot2cayley(ransac.model_coefficients_.block<3,3>(0,0))[2];
  double ransacError = asin((2 * z2) / (1 + pow(z2, 2))) / -theta;

  double z3 = opengv::math::rot2cayley(optimizedModel.block<3,3>(0,0))[2];
  double optimizedError = asin((2 * z3) / (1 + pow(z3, 2))) / -theta;

  //print the results
  std::cout << "the ransac results is: " << std::endl;
  std::cout << ransac.model_coefficients_ << std::endl << std::endl;
  std::cout << ransacError << std::endl << std::endl;
  std::cout << "Ransac needed " << ransac.iterations_ << std::endl;
  for(size_t i = 0; i < ransac.inliers_.size(); i++){
    std::cout << "inliers_" <<ransac.inliers_.size()<<": ";
    for(size_t j = 0; j < ransac.inliers_[i].size(); j++)
      std::cout << ransac.inliers_[i][j] << " ";
    std::cout << std::endl;
  }



/*    std::cout << "the optimized result is: " << std::endl;
    std::cout << optimizedError << std::endl << std::endl;*/

/*  std::cout << ransac.model_coefficients_.rotation << std::endl << std::endl;
  std::cout << "Ransac needed " << ransac.iterations_ << " iterations and ";
  std::cout << ransac_time << " seconds" << std::endl << std::endl;
  std::cout << "the number of inliers is: " << ransac.inliers_.size();
  std::cout << std::endl << std::endl;
  std::cout << "the found inliers are: " << std::endl;
  for(size_t i = 0; i < ransac.inliers_.size(); i++)
    std::cout << ransac.inliers_[i] << " ";
  std::cout << std::endl << std::endl;
  std::cout << "the optimized result is: " << std::endl;
  std::cout << optimizedModel.rotation << std::endl;*/

  // Create Lmeds
/*    sac::Lmeds<polyview::demos::EigensolverSacProblem> lmeds;
    lmeds.sac_model_ = eigenproblem_ptr;
    lmeds.threshold_ = 1.0;
    lmeds.max_iterations_ = 2000;

    //Run the experiment
    gettimeofday(&tic, 0);
    lmeds.computeModel();
    gettimeofday(&toc, 0);
    double lmeds_time = TIMETODOUBLE(timeval_minus(toc, tic));

    double z4 = opengv::math::rot2cayley(lmeds.model_coefficients_.rotation)[2];
    double lmedError = asin((2 * z4) / (1 + pow(z4, 2))) / -theta;
    std::cout << "the lmeds results is: " << std::endl;
    std::cout << lmedError << std::endl << std::endl;*/

  //print the results
/*  std::cout << "the lmeds results is: " << std::endl;
  std::cout << lmeds.model_coefficients_.rotation << std::endl << std::endl;
  std::cout << "lmeds needed " << lmeds.iterations_ << " iterations and ";
  std::cout << lmeds_time << " seconds" << std::endl << std::endl;
  std::cout << "the number of inliers is: " << lmeds.inliers_.size();
  std::cout << std::endl << std::endl;
  std::cout << "the found inliers are: " << std::endl;
  for(size_t i = 0; i < lmeds.inliers_.size(); i++)
    std::cout << lmeds.inliers_[i] << " ";
  std::cout << std::endl << std::endl;*/
}
