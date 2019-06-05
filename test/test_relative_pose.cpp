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
#include <opengv/relative_pose/CentralRelativeAdapter.hpp>
#include <sstream>
#include <fstream>

#include "random_generators.hpp"
#include "experiment_helpers.hpp"
#include "time_measurement.hpp"
#include <opencv2/opencv.hpp>
#include <opencv2/core/eigen.hpp>
#include <opencv2/calib3d.hpp>
#include <opengv/math/cayley.hpp>

using namespace std;
using namespace Eigen;
using namespace opengv;

int main(int argc, char **argv) {
    // initialize random seed
    initializeRandomSeed();

    for (int level = 0; level < 11; level++) {
        double noise = level * 0.5;
        std::cout << "noise level: " << noise << std::endl;
        for (size_t loop = 0; loop < 1000; loop++) {
            //set experiment parameters
            double outlierFraction = 0.0;
            size_t numberPoints = 100;

            //generate a random pose for viewpoint 1
            translation_t position1 = Eigen::Vector3d::Zero();
            rotation_t rotation1 = Eigen::Matrix3d::Identity();

            //generate a random pose for viewpoint 2
            //  translation_t position2 = generateRandomDirectionTranslation(0.0005);
            //  rotation_t rotation2 = generateRandomRotation(0.5);

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

            double angle = opengv::math::rot2cayley(rotation2)[2];
            angle = asin((2 * angle) / (1 + pow(angle, 2)));

            //create a fake central camera
            translations_t camOffsets;
            rotations_t camRotations;


            bool use_random_offset = false;

            if (use_random_offset) {
                generateCentralCameraSystem(camOffsets, camRotations);
            } else {
                camOffsets.resize(1);
                camRotations.resize(1);
                size_t mode = 1;
                switch (mode) {
                    case 1:
                        camOffsets[0] << 0, 0, 0;
                        camRotations[0] << 1, 0, 0, 0, 0, 1, 0, -1, 0;
                        break;
                    case 2:
                        camOffsets[0] << 0.2, 0, 0;
                        camRotations[0] << 0, 0, 1, -1, 0, 0, 0, -1, 0;
                        break;
                    case 3:
                        camOffsets[0] << 0, 0.3, 0;
                        camRotations[0] << 1, 0, 0, 0, 0, 1, 0, -1, 0;
                        break;
                    case 4:
                        camOffsets[0] << 0, -0.3, 0;
                        camRotations[0] << -1, 0, 0, 0, 0, -1, 0, -1, 0;
                        break;
                }
            }

            //derive correspondences based on random point-cloud
            bearingVectors_t bearingVectors1;
            bearingVectors_t bearingVectors2;
            std::vector<int> camCorrespondences1; //unused in the central case
            std::vector<int> camCorrespondences2; //unused in the central case
            Eigen::MatrixXd gt(3, numberPoints);
            generateRandom2D2DCorrespondences(
                    position1, rotation1, position2, rotation2,
                    camOffsets, camRotations, numberPoints, noise, outlierFraction,
                    bearingVectors1, bearingVectors2,
                    camCorrespondences1, camCorrespondences2, gt);

            //Extract the relative pose
            translation_t position;
            rotation_t rotation;
            extractRelativePose(
                    position1, position2, rotation1, rotation2, position, rotation);

            //print experiment characteristics
//            printExperimentCharacteristics(position, rotation, noise, outlierFraction);

            //compute and print the essential-matrix
        //    printEssentialMatrix(position, rotation);

            //create a central relative adapter
            relative_pose::CentralRelativeAdapter adapter(
                    bearingVectors1,
                    bearingVectors2,
                    rotation);

            //timer
            struct timeval tic;
            struct timeval toc;
            size_t iterations = 50;

            //running experiments
        //    std::cout << "running eightpt" << std::endl;
            essential_t eightpt_essential;
            gettimeofday(&tic, 0);
            for (size_t i = 0; i < iterations; i++)
                eightpt_essential = relative_pose::eightpt(adapter);
            gettimeofday(&toc, 0);
            double eightpt_time = TIMETODOUBLE(timeval_minus(toc, tic)) / iterations;


        //    std::cout << "results from eight-point algorithm:" << std::endl;
        //    std::cout << eightpt_essential << std::endl << std::endl;

            cv::Mat_<double> essential;
            cv::eigen2cv(eightpt_essential, essential);
            cv::Mat_<double> R1, R2;
            cv::Mat_<double> t;
            cv::decomposeEssentialMat(essential, R1, R2, t);

            Eigen::Matrix3d R1_e, R2_e;
            Eigen::Vector3d t_e;
            cv::cv2eigen(R1, R1_e);
            cv::cv2eigen(R2, R2_e);
            cv::cv2eigen(t, t_e);

            Eigen::Vector3d z1 = opengv::math::rot2cayley(R1_e);
            Eigen::Vector3d z2 = opengv::math::rot2cayley(R2_e);
            Eigen::Vector3d z;
            if (z1.norm() > 1)
                z = z2;
            else
                z = z1;

            rotation_t final_rot = opengv::math::cayley2rot(z);
            final_rot = camRotations[0] * final_rot * camRotations[0].transpose();
//            std::cout << "correct R : " <<std::endl<< final_rot << std::endl;
//
//            std::cout << "R1 : " <<std::endl<< camRotations[0] * R1_e * camRotations[0].transpose() << std::endl;
//            std::cout << "R2 : " <<std::endl<< camRotations[0] * R2_e * camRotations[0].transpose() << std::endl;

            double z3 = opengv::math::rot2cayley(rotation2.transpose() * final_rot)[2];
            double optimizedError = asin((2 * z3) / (1 + pow(z3, 2)));

            std::ofstream output1;

            std::string basePath("/home/ifwang/Documents/3DV/experiment/simulation");
            std::string FolderPath("eightpt");

            char t_error_Path[200];
            sprintf(t_error_Path, "%s/%s/noise_%.1f.txt", basePath.c_str(), FolderPath.c_str(), noise);

            output1.precision(6);
            output1.open(t_error_Path, std::ios::app);
            output1 << abs(optimizedError) << std::endl;
        }
    }
}
