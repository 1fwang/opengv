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
#include <opengv/relative_pose/NoncentralRelativeMultiAdapter.hpp>
#include <opengv/sac_problems/relative_pose/MultiNoncentralRelativePoseSacProblem.hpp>
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

    for (int level = 0; level < 11; level++) {
        double noise = 2.0;
        std::cout<<"noise level: "<<noise<<std::endl;
        for (size_t loop = 0; loop < 1; loop++) {
            // initialize random seed
            initializeRandomSeed();

            //set experiment parameters
            double outlierFraction = 0.0;
            size_t numberPoints = 20;

            //generate a random pose for viewpoint 1
            translation_t position1 = Eigen::Vector3d::Zero();
            rotation_t rotation1 = Eigen::Matrix3d::Identity();

            //generate a random pose for viewpoint 2

            translation_t position2 = Eigen::Vector3d::Zero();
            rotation_t rotation2 = Eigen::Matrix3d::Identity();
            double random_t = (((double) rand()) / ((double) RAND_MAX)) * 0.5 + 0.5;
            for (int i = 0; i < 6; i++) {
            //    double theta = (0.75 + i * 0.05) * random_t;
                double theta = 0.2 + double(i * level) / 25 ;
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
                    position1, position2, rotation1, rotation2, position, rotation);

            //create a central relative adapter
            polyview::demo::EigensolverAdapter adapter(
                    multiBearingVectors1,
                    multiBearingVectors2,
                    camOffsets,
                    camRotations,
                    position,
                    rotation);

            relative_pose::NoncentralRelativeMultiAdapter adapter2(
                    multiBearingVectors1,
                    multiBearingVectors2,
                    camOffsets,
                    camRotations);

/*            rotation_t R_perturbed;
            double theta_P = 0.9 * random_t * 6 * M_PI / 180.0;
            R_perturbed << cos(theta_P), sin(theta_P), 0,
                    -sin(theta_P), cos(theta_P), 0,
                    0, 0, 1;

            adapter.setR12(R_perturbed);*/

            //Create an EigensolverSacProblem and Ransac
            //The number of samples can be configured

            struct timeval tic;
            struct timeval toc;

            gettimeofday( &tic, 0 );
            opengv::sac::MultiRansac<
                    polyview::demos::EigensolverSacProblem> ransac_eig;
            std::shared_ptr<
                    polyview::demos::EigensolverSacProblem> eigenproblem_ptr(
                    new polyview::demos::EigensolverSacProblem(adapter, 9));
            ransac_eig.sac_model_ = eigenproblem_ptr;
            ransac_eig.threshold_ = 2.0 * (1.0 - cos(atan(sqrt(2.0) * 0.5 / 800.0)));
        //    ransac_eig.threshold_ = 100;
            ransac_eig.max_iterations_ = 1000;

            ransac_eig.computeModel();
            gettimeofday( &toc, 0 );
            double our_time = TIMETODOUBLE(timeval_minus(toc,tic));
            size_t our_iter = ransac_eig.iterations_;
            double our_inliers = 0;
            for (size_t i = 0; i < ransac_eig.inliers_.size(); i++)
                for (int j = 0; j < ransac_eig.inliers_[i].size(); j++) {
                    if (ransac_eig.inliers_[i][j] > numberPoints * outlierFraction) {
                        our_inliers = our_inliers + 1;
                    }
                }

            our_inliers = our_inliers / (numberPoints * 4 * (1 - outlierFraction));

            polyview::demos::EigensolverSacProblem::model_t optimizedModel_eig;
            eigenproblem_ptr->optimizeModelCoefficients(
                    ransac_eig.inliers_,
                    ransac_eig.model_coefficients_,
                    optimizedModel_eig);


            double z = opengv::math::rot2cayley(rotation2.transpose() * ransac_eig.model_coefficients_.block<3, 3>(0, 0))[2];
            double optimizedError = asin((2 * z) / (1 + pow(z, 2)));

            std::cout<<"optimizedError: "<<optimizedError<<std::endl;
/*            std::ofstream output1;

            std::string basePath("/home/ifwang/Documents/3DV/experiment/simulation_minimal");
            std::string FolderPath("our");

            char t_error_Path[200];
            sprintf(t_error_Path, "%s/%s/noise_%.1f.txt", basePath.c_str(), FolderPath.c_str(), double(level) / 10);

            output1.precision(6);
            output1.open(t_error_Path, std::ios::app);
            output1 << fabs(optimizedError) << std::endl;*/


            /*gettimeofday( &tic, 0 );
            opengv::sac::MultiRansac<
                    sac_problems::relative_pose::MultiNoncentralRelativePoseSacProblem> ransac;
            std::shared_ptr<
                    sac_problems::relative_pose::MultiNoncentralRelativePoseSacProblem> relposeproblem_ptr(
                    new sac_problems::relative_pose::MultiNoncentralRelativePoseSacProblem(
                            adapter2,
                            sac_problems::relative_pose::MultiNoncentralRelativePoseSacProblem::SEVENTEENPT));
            ransac.sac_model_ = relposeproblem_ptr;
            ransac.threshold_ = 2.0 * (1.0 - cos(atan(sqrt(2.0) * 3.0 / 800.0)));
            //ransac.threshold_ = 100;
            ransac.max_iterations_ = 1000;

            //Run the experiment
            ransac.computeModel();
            gettimeofday( &toc, 0 );
            double seventeen_time = TIMETODOUBLE(timeval_minus(toc,tic));
            size_t seventeen_iter = ransac.iterations_;
            double seventeen_inliers = 0;
            for (size_t i = 0; i < ransac.inliers_.size(); i++)
                for (int j = 0; j < ransac.inliers_[i].size(); j++) {
                    if (ransac.inliers_[i][j] > numberPoints * outlierFraction) {
                        seventeen_inliers = seventeen_inliers + 1;
                    }
                }
            seventeen_inliers = seventeen_inliers / (numberPoints * 4 * (1 - outlierFraction));*/

/*            for(size_t i = 0; i < ransac_eig.inliers_.size(); i++)
                std::cout<<"numberInliers of cam "<<i <<" has "<< ransac_eig.inliers_[i].size()<<std::endl;*/

            //do final polishing of the model over all inliers

 /*           polyview::demos::EigensolverSacProblem::model_t optimizedModel_eig;
            eigenproblem_ptr->optimizeModelCoefficients(
                    ransac_eig.inliers_,
                    ransac_eig.model_coefficients_,
                    optimizedModel_eig);


            sac_problems::relative_pose::MultiNoncentralRelativePoseSacProblem::model_t optimizedModel;
            relposeproblem_ptr->optimizeModelCoefficients(
                    ransac.inliers_,
                    ransac.model_coefficients_,
                    optimizedModel);

            double z_eig = opengv::math::rot2cayley(rotation2.transpose() * optimizedModel_eig.block<3, 3>(0, 0))[2];
            double optimizedError_eig = asin((2 * z_eig) / (1 + pow(z_eig, 2)));


            double z = opengv::math::rot2cayley(rotation2.transpose() * optimizedModel.block<3, 3>(0, 0))[2];
            double optimizedError = asin((2 * z) / (1 + pow(z, 2)));


            std::ofstream output1,output2,output3;

            std::string basePath("/home/ifwang/Documents/3DV/experiment/simulation_iter");
            std::string FolderPath1("our");
            std::string FolderPath2("SEVENTEENPT");

            char error_Path[200];
            char error_Path_eig[200];
            char time[200];
            sprintf(error_Path_eig, "%s/%s/noise_%.1f.txt", basePath.c_str(),FolderPath1.c_str(),noise);
            sprintf(error_Path, "%s/%s/noise_%.1f.txt", basePath.c_str(),FolderPath2.c_str(),noise);
            sprintf(time, "%s/%s/time_%.1f.txt", basePath.c_str(),FolderPath1.c_str(),noise);

            output1.precision(6);
            output2.precision(6);
            output3.precision(6);
            output1.open(error_Path_eig, std::ios::app);
            output1 << abs(optimizedError_eig) << std::endl;

            output2.open(error_Path, std::ios::app);
            output2 << abs(optimizedError) << std::endl;

            output3.open(time, std::ios::app);
            output3 << eig_time << std::endl;*/

            /*std::ofstream output1,output2,output3,output4,output5,output6;
            std::string basePath("/home/ifwang/Documents/3DV/experiment/simulation_iter");
            std::string FolderPath1("our");
            std::string FolderPath2("SEVENTEENPT");

            char time_our[200];
            char iter_our[200];
            char inlier_our[200];
            char time_seventeen[200];
            char iter_seventeen[200];
            char inlier_seventeen[200];

            sprintf(time_our, "%s/%s/time_%.2f.txt", basePath.c_str(),FolderPath1.c_str(),outlierFraction);
            sprintf(iter_our, "%s/%s/iter_%.2f.txt", basePath.c_str(),FolderPath1.c_str(),outlierFraction);
            sprintf(inlier_our, "%s/%s/inlier_%.2f.txt", basePath.c_str(),FolderPath1.c_str(),outlierFraction);

            sprintf(time_seventeen, "%s/%s/time_%.2f.txt", basePath.c_str(),FolderPath2.c_str(),outlierFraction);
            sprintf(iter_seventeen, "%s/%s/iter_%.2f.txt", basePath.c_str(),FolderPath2.c_str(),outlierFraction);
            sprintf(inlier_seventeen, "%s/%s/inlier_%.2f.txt", basePath.c_str(),FolderPath2.c_str(),outlierFraction);

            output1.precision(6);
            output1.open(time_our, std::ios::app);
            output1 << our_time << std::endl;

            output2.precision(6);
            output2.open(iter_our, std::ios::app);
            output2 << our_iter << std::endl;

            output3.precision(6);
            output3.open(inlier_our, std::ios::app);
            output3 << our_inliers << std::endl;

            output4.precision(6);
            output4.open(time_seventeen, std::ios::app);
            output4 << seventeen_time << std::endl;

            output5.precision(6);
            output5.open(iter_seventeen, std::ios::app);
            output5 << seventeen_iter << std::endl;

            output6.precision(6);
            output6.open(inlier_seventeen, std::ios::app);
            output6 << seventeen_inliers << std::endl;*/

        }
    }
    std::cout<<"done"<<std::endl;
}
