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
#include <limits.h>
#include <Eigen/Eigen>
#include <opengv/relative_pose/methods.hpp>
#include <opengv/relative_pose/NoncentralRelativeMultiAdapter.hpp>
#include <opengv/sac/MultiRansac.hpp>
#include <opengv/sac_problems/relative_pose/MultiNoncentralRelativePoseSacProblem.hpp>
#include <sstream>
#include <fstream>

#include "random_generators.hpp"
#include "experiment_helpers.hpp"
#include "time_measurement.hpp"
#include <opengv/math/cayley.hpp>

using namespace std;
using namespace Eigen;
using namespace opengv;

int main(int argc, char **argv) {

    for (int level = 0; level < 11; level++) {
        double noise = 2.0;
        std::cout<<"noise level: "<<noise<<std::endl;
        for (size_t loop = 0; loop < 1000; loop++) {
            // initialize random seed
            initializeRandomSeed();

            //set experiment parameters
            double outlierFraction = level * 0.03;
            size_t numberPoints = 100;

            //generate a random pose for viewpoint 1
            translation_t position1 = Eigen::Vector3d::Zero();
            rotation_t rotation1 = Eigen::Matrix3d::Identity();

            //generate a random pose for viewpoint 2

            translation_t position2 = Eigen::Vector3d::Zero();
            rotation_t rotation2 = Eigen::Matrix3d::Identity();
            double random_t = (((double) rand()) / ((double) RAND_MAX)) * 0.5 + 0.5;
            for (int i = 0; i < 6; i++) {
                double theta = (0.75 + i * 0.05) * random_t;
                theta = theta / 180 * M_PI;
                rotation_t rel_r;
                translation_t rel_t;
                rel_r << cos(theta), -sin(theta), 0,
                        sin(theta), cos(theta), 0,
                        0, 0, 1;
                rel_t << -sin(theta), cos(theta), 0;
                rel_t = 0.1 * rel_t / rel_t.norm();
                position2 = position2 + rotation2 * rel_t;
                rotation2 = rotation2 * rel_r;
            }
            //    std::cout << rotation2 << std::endl;

            double angle = opengv::math::rot2cayley(rotation2)[2];
            angle = asin((2 * angle) / (1 + pow(angle, 2)));

            //  std::cout << "rotation2: " << rotation2 << std::endl;
            //    std::cout << "position2: " << position2.transpose() << std::endl;
            //    std::cout << "direction: " << (position2 / position2.norm()).transpose() << std::endl;

            //create a fake central camera
            translations_t camOffsets;
            rotations_t camRotations;

            bool use_random_offset = false;

            if (use_random_offset) {
                generateRandomCameraSystem(4, camOffsets, camRotations);
            } else {
                camOffsets.resize(4);
                camRotations.resize(4);
                camOffsets[0] << -0.6, 0, 0;
                camOffsets[1] << 0.6, 0, 0;
                camOffsets[2] << 0, -1.0, 0;
                camOffsets[3] << 0, 1.0, 0;

                camRotations[0] << 0, 0, -1, 1, 0, 0, 0, -1, 0;
                camRotations[1] << 0, 0, 1, -1, 0, 0, 0, -1, 0;
                camRotations[2] << 1, 0, 0, 0, 0, 1, 0, -1, 0;
                camRotations[3] << -1, 0, 0, 0, 0, -1, 0, -1, 0;

            }

            //derive correspondences based on random point-cloud
            std::vector<std::shared_ptr<bearingVectors_t> > multiBearingVectors1;
            std::vector<std::shared_ptr<bearingVectors_t> > multiBearingVectors2;
            std::vector<std::shared_ptr<Eigen::MatrixXd> > gt;

            generateMulti2D2DCorrespondences(
                    position1, rotation1, position2, rotation2, camOffsets, camRotations,
                    numberPoints, noise, outlierFraction,
                    multiBearingVectors1, multiBearingVectors2, gt);

            //Extract the relative pose
            translation_t position;
            rotation_t rotation;
            extractRelativePose(
                    position1, position2, rotation1, rotation2, position, rotation, false);

            //create a non-central relative multi-adapter
            relative_pose::NoncentralRelativeMultiAdapter adapter(
                    multiBearingVectors1,
                    multiBearingVectors2,
                    camOffsets,
                    camRotations);

            struct timeval tic;
            struct timeval toc;

            //Create a MultiNoncentralRelativePoseSacProblem and Ransac
            gettimeofday( &tic, 0 );
            opengv::sac::MultiRansac<
                    sac_problems::relative_pose::MultiNoncentralRelativePoseSacProblem> ransac;
            std::shared_ptr<
                    sac_problems::relative_pose::MultiNoncentralRelativePoseSacProblem> relposeproblem_ptr(
                    new sac_problems::relative_pose::MultiNoncentralRelativePoseSacProblem(
                            adapter,
                            sac_problems::relative_pose::MultiNoncentralRelativePoseSacProblem::GE));
            ransac.sac_model_ = relposeproblem_ptr;
            ransac.threshold_ = 2.0 * (1.0 - cos(atan(sqrt(2.0) * 3.0 / 800.0)));
        //    ransac.threshold_ = 100;
            ransac.max_iterations_ = 1000;

            //Run the experiment
            ransac.computeModel();
            gettimeofday( &toc, 0 );
            double eight_time = TIMETODOUBLE(timeval_minus(toc,tic));
            size_t eight_iter = ransac.iterations_;
            double eight_inliers = 0;
            for (size_t i = 0; i < ransac.inliers_.size(); i++)
                for (int j = 0; j < ransac.inliers_[i].size(); j++) {
                    if (ransac.inliers_[i][j] > numberPoints * outlierFraction) {
                        eight_inliers = eight_inliers + 1;
                    }
                }

            eight_inliers = eight_inliers / (numberPoints * 4 * (1 - outlierFraction));

            /*sac_problems::relative_pose::MultiNoncentralRelativePoseSacProblem::model_t optimizedModel;
            relposeproblem_ptr->optimizeModelCoefficients(
                    ransac.inliers_,
                    ransac.model_coefficients_,
                    optimizedModel);

            double z3 = opengv::math::rot2cayley(rotation2.transpose() * optimizedModel.block<3, 3>(0, 0))[2];
            double optimizedError = asin((2 * z3) / (1 + pow(z3, 2)));

            //    std::cout << optimizedModel << std::endl << std::endl;
            //    std::cout << abs(optimizedError) << std::endl << std::endl;

            std::ofstream output1;

            std::string basePath("/home/ifwang/Documents/3DV/experiment/simulation");
            std::string FolderPath("GE");

            char t_error_Path[200];
            sprintf(t_error_Path, "%s/%s/noise_%.1f.txt", basePath.c_str(), FolderPath.c_str(), noise);

            output1.precision(6);
            output1.open(t_error_Path, std::ios::app);
            output1 << abs(optimizedError) << std::endl;*/

            std::ofstream output1,output2,output3;
            std::string basePath("/home/ifwang/Documents/3DV/experiment/simulation_iter");
            std::string FolderPath1("GE");

            char time_eight[200];
            char iter_eight[200];
            char inlier_eight[200];

            sprintf(time_eight, "%s/%s/time_%.2f.txt", basePath.c_str(),FolderPath1.c_str(),outlierFraction);
            sprintf(iter_eight, "%s/%s/iter_%.2f.txt", basePath.c_str(),FolderPath1.c_str(),outlierFraction);
            sprintf(inlier_eight, "%s/%s/inlier_%.2f.txt", basePath.c_str(),FolderPath1.c_str(),outlierFraction);

            output1.precision(6);
            output1.open(time_eight, std::ios::app);
            output1 << eight_time << std::endl;

            output2.precision(6);
            output2.open(iter_eight, std::ios::app);
            output2 << eight_iter << std::endl;

            output3.precision(6);
            output3.open(inlier_eight, std::ios::app);
            output3 << eight_inliers << std::endl;

        }
    }
    std::cout<<"done"<<std::endl;

/*  //print the results
  std::cout << "the ransac threshold is: " << ransac.threshold_ << std::endl;
  std::cout << "the ransac results is: " << std::endl;
  std::cout << optimizedModel << std::endl << std::endl;
  std::cout << "Ransac needed " << ransac.iterations_ << " iterations and ";
  std::cout << ransac_time << " seconds" << std::endl << std::endl;
  size_t numberInliers = 0;
  for(size_t i = 0; i < ransac.inliers_.size(); i++)
    numberInliers += ransac.inliers_[i].size();
  std::cout << "the number of inliers is: " << numberInliers;
  std::cout << std::endl << std::endl;
  std::cout << "the found inliers are: " << std::endl;
  for(size_t i = 0; i < ransac.inliers_.size(); i++)
  {
    for(size_t j = 0; j < ransac.inliers_[i].size(); j++)
      std::cout << ransac.inliers_[i][j] << " ";
  }
  std::cout << std::endl << std::endl;*/
}
