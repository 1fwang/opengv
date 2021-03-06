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

#include "EigensolverSacProblem.hpp"
#include <opengv/relative_pose/modules/eigensolver/modules.hpp>
#include <opengv/relative_pose/methods.hpp>
#include <opengv/math/cayley.hpp>
#include <opengv/triangulation/methods.hpp>

#include <opengv/OptimizationFunctor.hpp>
#include <unsupported/Eigen/NonLinearOptimization>
#include <unsupported/Eigen/NumericalDiff>
#include <iostream>
#include <fstream>
#include <sstream>

double
polyview::demos::EigensolverSacProblem::
getSmallestEVwithJacobian(
        const EigensolverSacProblem::adapter_t &adapter,
        const std::vector<std::vector<int> > &indices,
        const size_t &camID,
        const opengv::cayley_t &cayley,
        Eigen::Matrix<double, 1, 3> &jacobian) {
    Eigen::Matrix3d M = Eigen::Matrix3d::Zero();
    //  std::cout<<"CamID: "<<camID<<std::endl;
    for (size_t i = 0; i < indices[camID].size(); i++) {
        opengv::bearingVector_t f1 =
                adapter.getCamRotation(camID) * adapter.getBearingVector1(camID, indices[camID][i]);
        opengv::bearingVector_t f2 =
                adapter.getCamRotation(camID) * adapter.getBearingVector2(camID, indices[camID][i]);
        opengv::rotation_t R = opengv::math::cayley2rot(cayley);
        Eigen::Vector3d n = f1.cross(R * f2);
        //    std::cout<<"n norm: "<<n.norm()<<std::endl;
        /*n = n / n.norm();*/
        M = M + n * n.transpose();
    }
    Eigen::EigenSolver<Eigen::Matrix3d> Eig(M, true);
    Eigen::Matrix<std::complex<double>, 3, 1> D_complex = Eig.eigenvalues();
    Eigen::Matrix<std::complex<double>, 3, 3> V_complex = Eig.eigenvectors();


    opengv::eigenvalues_t D;
    opengv::eigenvectors_t V;
    for (size_t i = 0; i < 3; i++) {
        D[i] = D_complex[i].real();
        for (size_t j = 0; j < 3; j++)
            V(i, j) = V_complex(i, j).real();
    }

    int index;
    double temp;
    index = 0;
    temp = D[0];

    for (int k = 0; k < 3; k++) {
        if (D[k] < temp && D[k] != 0) {
            temp = D[k];
            index = k;
        }
    }

    //std::cout<<"smallest Eigenvalue: "<<temp<<std::endl;

    opengv::translation_t t_1 = V.col(index);
    jacobian(0, 2) = t_1[2];
    //t_1[2] = 0;
    double distance = 0;
    for (size_t i = 0; i < indices[camID].size(); i++) {
        opengv::bearingVector_t f1 =
                adapter.getCamRotation(camID) * adapter.getBearingVector1(camID, indices[camID][i]);
        opengv::bearingVector_t f2 =
                adapter.getCamRotation(camID) * adapter.getBearingVector2(camID, indices[camID][i]);
        opengv::rotation_t R = opengv::math::cayley2rot(cayley);
        Eigen::Vector3d n = f1.cross(R * f2);
        //    std::cout<<"n norm: "<<n.norm()<<std::endl;
        n = n / n.norm();
        distance = distance + pow(n.transpose() * t_1, 2);
    }
    //std::cout<<"distance: "<<distance<<std::endl;

    //Retrieve the smallest Eigenvalue by the following closed form solution.
    //Plus Jacobian.
//    M = opengv::relative_pose::modules::eigensolver::composeMwithJacobians(
//            xxF,yyF,zzF,xyF,yzF,zxF,cayley,M_jac1,M_jac2,M_jac3);
    return distance;
}

namespace polyview {
    namespace demos {
        struct Eigensolver_step : opengv::OptimizationFunctor<double> {
            const std::vector<Eigen::Matrix3d> &_xxF;
            const std::vector<Eigen::Matrix3d> &_yyF;
            const std::vector<Eigen::Matrix3d> &_zzF;
            const std::vector<Eigen::Matrix3d> &_xyF;
            const std::vector<Eigen::Matrix3d> &_yzF;
            const std::vector<Eigen::Matrix3d> &_zxF;
            const EigensolverSacProblem::adapter_t &_adapter;
            const std::vector<std::vector<int> > &_indices;
            const int _rotation_type;

            Eigensolver_step(
                    const std::vector<Eigen::Matrix3d> &xxF,
                    const std::vector<Eigen::Matrix3d> &yyF,
                    const std::vector<Eigen::Matrix3d> &zzF,
                    const std::vector<Eigen::Matrix3d> &xyF,
                    const std::vector<Eigen::Matrix3d> &yzF,
                    const std::vector<Eigen::Matrix3d> &zxF,
                    const EigensolverSacProblem::adapter_t &adapter,
                    const std::vector<std::vector<int> > &indices,
                    const int rotation_type) :
                    opengv::OptimizationFunctor<double>(1, 1),
                    _xxF(xxF), _yyF(yyF), _zzF(zzF), _xyF(xyF), _yzF(yzF), _zxF(zxF), _adapter(adapter),
                    _indices(indices), _rotation_type(rotation_type) {}

            int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const {
                opengv::cayley_t cayley;
                switch (_rotation_type) {
                    case 1:
                        cayley << x, 0, 0;
                        break;
                    case 2:
                        cayley << 0, x, 0;
                        break;
                    case 3:
                        cayley << 0, 0, x;
                        break;
                    default:
                        std::cout << "Not A planar Motion!" << std::endl;
                }

                //    std::ofstream output;
                //    std::string folder("/home/ifwang/datasets/energy.txt");
                //    output.precision(6);
                //    output.open(folder, std::ios::app);
                fvec[0] = 0;
                for (int i = 0; i < _xxF.size(); i++) {
                    Eigen::Matrix<double, 1, 3> jacobian;
                    double eigenvalue;
                    //    eigenvalue = polyview::demos::EigensolverSacProblem::getSmallestEVwithJacobian(
                    //            _xxF[i], _yyF[i], _zzF[i], _xyF[i], _yzF[i], _zxF[i], _adapter, _indices, i, cayley, jacobian);

                    eigenvalue = opengv::relative_pose::modules::eigensolver::getSmallestEVwithJacobian(
                            _xxF[i], _yyF[i], _zzF[i], _xyF[i], _yzF[i], _zxF[i], cayley, jacobian);
                    fvec[0] = fvec[0] + pow(eigenvalue, 2);
                    //    fvec[i] = eigenvalue;

                }

                return 0;
            }
        };

    }
}

namespace polyview {
    namespace demos {
        struct Eigensolver_optimize : opengv::OptimizationFunctor<double> {
            const std::vector<Eigen::Matrix3d> &_xxF;
            const std::vector<Eigen::Matrix3d> &_yyF;
            const std::vector<Eigen::Matrix3d> &_zzF;
            const std::vector<Eigen::Matrix3d> &_xyF;
            const std::vector<Eigen::Matrix3d> &_yzF;
            const std::vector<Eigen::Matrix3d> &_zxF;
            const EigensolverSacProblem::adapter_t &_adapter;
            const std::vector<std::vector<int> > &_indices;
            const int _rotation_type;

            Eigensolver_optimize(
                    const std::vector<Eigen::Matrix3d> &xxF,
                    const std::vector<Eigen::Matrix3d> &yyF,
                    const std::vector<Eigen::Matrix3d> &zzF,
                    const std::vector<Eigen::Matrix3d> &xyF,
                    const std::vector<Eigen::Matrix3d> &yzF,
                    const std::vector<Eigen::Matrix3d> &zxF,
                    const EigensolverSacProblem::adapter_t &adapter,
                    const std::vector<std::vector<int> > &indices,
                    const int rotation_type) :
                    opengv::OptimizationFunctor<double>(1, 1),
                    _xxF(xxF), _yyF(yyF), _zzF(zzF), _xyF(xyF), _yzF(yzF), _zxF(zxF), _adapter(adapter), _indices(indices), _rotation_type(rotation_type) {}

            int operator()(const Eigen::VectorXd &x, Eigen::VectorXd &fvec) const {
                opengv::cayley_t cayley;
                switch (_rotation_type) {
                    case 1:
                        cayley << x, 0, 0;
                        break;
                    case 2:
                        cayley << 0, x, 0;
                        break;
                    case 3:
                        cayley << 0, 0, x;
                        break;
                    default:
                        std::cout << "Not A planar Motion!" << std::endl;
                }

                fvec[0] = 0;
//                fvec[1] = 0;
                //    std::cout<<"number of cams: "<<_adapter.getNumberPairs()<<std::endl;
                for (int i = 0; i < _adapter.getNumberPairs(); i++) {
                    Eigen::Matrix<double, 1, 3> jacobian;
                    double eigenvalue,eigenvalue2;
                    eigenvalue = polyview::demos::EigensolverSacProblem::getSmallestEVwithJacobian(
                            _adapter, _indices, i, cayley, jacobian);
                    fvec[0] = fvec[0] + pow(eigenvalue, 2);

/*                    eigenvalue2 = opengv::relative_pose::modules::eigensolver::getSmallestEVwithJacobian(
                            _xxF[i], _yyF[i], _zzF[i], _xyF[i], _yzF[i], _zxF[i], cayley, jacobian);
                    fvec[1] = fvec[1] + pow(eigenvalue2, 2);*/

                }
                return 0;
            }
        };

    }
}


bool
polyview::demos::EigensolverSacProblem::computeModelCoefficients(
        const std::vector<std::vector<int> > &indices,
        model_t &outModel) const {
    double maxVariation = 0.1; //****** 0.01 ******//

    //randomize the starting point a bit
    opengv::rotation_t rotation = _adapter.getR12();

    std::vector<Eigen::Matrix3d> xxF, yyF, zzF, xyF, yzF, zxF;

    //Initial summation terms
    for (size_t k = 0; k < indices.size(); k++) {
        xxF.push_back(Eigen::Matrix3d::Zero());
        yyF.push_back(Eigen::Matrix3d::Zero());
        zzF.push_back(Eigen::Matrix3d::Zero());
        xyF.push_back(Eigen::Matrix3d::Zero());
        yzF.push_back(Eigen::Matrix3d::Zero());
        zxF.push_back(Eigen::Matrix3d::Zero());

        //Fill summation terms
        for (size_t i = 0; i < indices[k].size(); i++) {
            opengv::bearingVector_t f1 = _adapter.getCamRotation(k) * _adapter.getBearingVector1(k, indices[k][i]);
            opengv::bearingVector_t f2 = _adapter.getCamRotation(k) * _adapter.getBearingVector2(k, indices[k][i]);
            Eigen::Matrix3d F = f2 * f2.transpose();

            double weight = 1.0;

            xxF[k] = xxF[k] + weight * f1[0] * f1[0] * F;
            yyF[k] = yyF[k] + weight * f1[1] * f1[1] * F;
            zzF[k] = zzF[k] + weight * f1[2] * f1[2] * F;
            xyF[k] = xyF[k] + weight * f1[0] * f1[1] * F;
            yzF[k] = yzF[k] + weight * f1[1] * f1[2] * F;
            zxF[k] = zxF[k] + weight * f1[2] * f1[0] * F;
        }
    }

    //Do minimization
    const int n = 3;
    int rotation_type = -1;
    Eigen::VectorXd x(n);
    Eigen::VectorXd x_planar(1);

    x = opengv::math::rot2cayley(rotation);
//    std::cout<<"x"<<x<<std::endl;

    if (x[0] != 0 && x[1] == 0 && x[2] == 0) {
        rotation_type = 1;
        x_planar << x[0];
    } else if (x[0] == 0 && x[1] != 0 && x[2] == 0) {
        rotation_type = 2;
        x_planar << x[1];
    } else if (x[0] == 0 && x[1] == 0) {
        rotation_type = 3;
        x_planar << x[2];
    }

    rotation_type = 3;
    //PLANAR MOTION
    Eigensolver_step functor(xxF, yyF, zzF, xyF, yzF, zxF, _adapter, indices, rotation_type);
    Eigen::NumericalDiff<Eigensolver_step> numDiff(functor);
    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<Eigensolver_step>> lm(numDiff);

    lm.resetParameters();
    lm.parameters.ftol = 0.00005;
    lm.parameters.xtol = 1.E1 * NumTraits<double>::epsilon();
    lm.parameters.maxfev = 1000;
    lm.minimize(x_planar);

    /*if (indices[0].size() > 3){
          std::ofstream output1,output2;
          std::string folder1("/home/ifwang/datasets/energy_test/eigenvalue.txt");
          std::string folder2("/home/ifwang/datasets/energy_test/jacobian.txt");
          output1.precision(6);
          output2.precision(6);
          output1.open(folder1, std::ios::trunc);
          output2.open(folder2, std::ios::trunc);

          for (int j = -1000; j < 1000; j++) {
              double eigen_energy = 0;
              double jacobian_energy = 0;
              for (int i = 0; i < xxF.size(); i++) {
                  opengv::cayley_t cayley;
                  cayley << 0, 0, 0.0001 * j;
                  double eigenvalue,eigenvalue2;
                  Eigen::Matrix<double, 1, 3> jacobian;
                  eigenvalue = opengv::relative_pose::modules::eigensolver::getSmallestEVwithJacobian(
                          xxF[i], yyF[i], zzF[i], xyF[i], yzF[i], zxF[i], cayley, jacobian);
                  eigenvalue2 = polyview::demos::EigensolverSacProblem::getSmallestEVwithJacobian(
                          _adapter, indices, i, cayley, jacobian);
                  output1 << eigenvalue << " ";
                  output2 << eigenvalue2 << " ";
                  eigen_energy = eigen_energy + pow(eigenvalue,2);
                  jacobian_energy = jacobian_energy + pow(eigenvalue2,2);
              }
              output1 << eigen_energy << " ";
              output2 << jacobian_energy << " ";
              output1 << std::endl;
              output2 << std::endl;
          }
      }*/


    switch (rotation_type) {
        case 1:
            x[0] = x_planar[0];
            break;
        case 2:
            x[1] = x_planar[0];
            break;
        case 3:
            x[2] = x_planar[0];
            break;
        default:
            std::cout << "Not A planar Motion!" << std::endl;
    }

    opengv::cayley_t cayley_computed = x;
    //   std::cout<<"x optimized"<<x<<std::endl;
    opengv::rotation_t R = opengv::math::cayley2rot(cayley_computed);
    opengv::translation_t t;

    Eigen::MatrixXd LargeMatrix = Eigen::MatrixXd::Zero(indices.size() * 3, indices.size() + 3);
    Eigen::VectorXd Entries = Eigen::VectorXd::Zero(indices.size() * 3);

//    LargeMatrix.block(0,indices.size(), indices.size() * 3,3) = - Eigen::MatrixXd::Ones( indices.size() * 3, 3);
    //GET absolute R and t with scale
    for (int i = 0; i < indices.size(); i++) {

/*        opengv::rotation_t rel_R = _adapter.getCamRotation(i).transpose() * R * _adapter.getCamRotation(i);
        opengv::cayley_t cayley_temp = opengv::math::rot2cayley(rel_R);*/

        Eigen::Matrix3d M = opengv::relative_pose::modules::eigensolver::composeM(xxF[i], yyF[i], zzF[i], xyF[i],
                                                                                  yzF[i], zxF[i], cayley_computed);
        Eigen::EigenSolver<Eigen::Matrix3d> Eig(M, true);
        Eigen::Matrix<std::complex<double>, 3, 1> D_complex = Eig.eigenvalues();
        Eigen::Matrix<std::complex<double>, 3, 3> V_complex = Eig.eigenvectors();
        opengv::eigenvalues_t D;
        opengv::eigenvectors_t V;
        for (size_t k = 0; k < 3; k++) {
            D[k] = D_complex[k].real();
            for (size_t j = 0; j < 3; j++)
                V(k, j) = V_complex(k, j).real();
        }

        int index;
        double temp;
        index = 0;
        temp = D[0];

        for (int k = 0; k < 3; k++) {
            if (D[k] < temp && D[k] != 0) {
                temp = D[k];
                index = k;
            }
        }
        //    double translationMagnitude = sqrt(pow(D[1], 2) + pow(D[2], 2));
        opengv::translation_t t = V.col(index);

        //Correct the translation
        //    std::cout<<"indices size()"<<indices[i].size()<<std::endl;

        size_t pick;
        for (size_t index = 0; index < indices.size(); index++) {
            if (indices[index].size() != 0) {
                pick = index;
                break;
            }
        }

        opengv::bearingVector_t f1 = _adapter.getCamRotation(pick) * _adapter.getBearingVector1(i, indices[pick][0]);
        opengv::bearingVector_t f2 = _adapter.getCamRotation(pick) * _adapter.getBearingVector2(i, indices[pick][0]);
        f2 = R * f2;
        Eigen::Vector3d opticalFlow = f1 - f2;
        if (opticalFlow.dot(t) < 0.0)
            t = -t;

        t[2] = 0;
//        std::cout<<"t"<<i<<": "<<(t/t.norm()).transpose()<<std::endl;
//        std::cout<<"R"<<i<<": "<<cayley_computed.transpose()<<std::endl;
        //    std::cout<<"D"<<i<<": "<<V<<std::endl;
        //    std::cout<<"V"<<i<<": "<<D<<std::endl;
        Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);

        LargeMatrix.block(i * 3, i, 3, 1) = t;
        LargeMatrix.block(i * 3, indices.size(), 3, 3) = -I;
        Entries.segment(i * 3, 3) = (R - I) * _adapter.getCamOffset(i);
    }

    /*   std::cout<<"LargeMatrix"<<std::endl;
       std::cout<<LargeMatrix<<std::endl;
       std::cout<<"Entries"<<std::endl;
       std::cout<<Entries<<std::endl;
       std::cout<<std::endl;*/

    Eigen::VectorXd solution = LargeMatrix.colPivHouseholderQr().solve(Entries);
    t = solution.segment(indices.size(), 3);
//    std::cout<<"centre T:"<<t.transpose()<<std::endl;

    outModel.col(3) = t;
    outModel.block<3, 3>(0, 0) = R;

    //Correct the translation
/*    opengv::bearingVector_t f1 = _adapter.getBearingVector1(i, indices[i][0]);
    opengv::bearingVector_t f2 = _adapter.getBearingVector2(i, indices[i][0]);
    f2 = outModel.block<3,3>(0,0) * f2;
    Eigen::Vector3d opticalFlow = f1 - f2;
    if (opticalFlow.dot(outModel.col(3)) < 0.0)
        outModel.col(3) = -outModel.col(3);*/

    return true;
}

void
polyview::demos::EigensolverSacProblem::getSelectedDistancesToModel(
        const model_t & model,
        const std::vector<std::vector<int> > & indices,
        std::vector<std::vector<double> > & scores) const
{
    opengv::translation_t init_t = _adapter.gett12();
    opengv::rotation_t init_R = _adapter.getR12();

    opengv::translation_t translation = model.col(3);
    opengv::rotation_t rotation = model.block<3,3>(0,0);
//    std::cout<<"translation: "<<std::endl<<translation.transpose()<<std::endl;
//    std::cout<<"rotation: "<<std::endl<<rotation<<std::endl;

    Eigen::Matrix<double,4,1> p_hom;
    p_hom[3] = 1.0;

    for( size_t camIndex = 0; camIndex < indices.size(); camIndex++ )
    {
        opengv::translation_t cam1Offset = _adapter.getCamOffset(camIndex);
        opengv::rotation_t cam1Rotation = _adapter.getCamRotation(camIndex);
        opengv::translation_t cam2Offset = _adapter.getCamOffset(camIndex);
        opengv::rotation_t cam2Rotation = _adapter.getCamRotation(camIndex);

        opengv::translation_t directTranslation =
                cam1Rotation.transpose() *
                ((translation - cam1Offset) + rotation * cam2Offset);
        opengv::rotation_t directRotation =
                cam1Rotation.transpose() * rotation * cam2Rotation;

        _adapter.sett12(directTranslation);
        _adapter.setR12(directRotation);

    /*    opengv::translation_t cam1Offset = _adapter.getCamOffset(camIndex);
        opengv::translation_t cam2Offset = _adapter.getCamOffset(camIndex);

        opengv::translation_t directTranslation = translation - cam1Offset + rotation * cam2Offset;
        opengv::rotation_t directRotation = rotation;

        _adapter.sett12(directTranslation);
        _adapter.setR12(directRotation);*/

        opengv::transformation_t inverseSolution;
        inverseSolution.block<3,3>(0,0) = directRotation.transpose();
        inverseSolution.col(3) =
                -inverseSolution.block<3,3>(0,0)*directTranslation;

        for(
                size_t correspondenceIndex = 0;
                correspondenceIndex < indices[camIndex].size();
                correspondenceIndex++ )
        {
            p_hom.block<3,1>(0,0) =
                    opengv::triangulation::triangulate2(
                            _adapter,
                            _adapter.convertMultiIndex(
                                    camIndex, indices[camIndex][correspondenceIndex] ));

            opengv::bearingVector_t reprojection1 = p_hom.block<3,1>(0,0);
            opengv::bearingVector_t reprojection2 = inverseSolution * p_hom;
            reprojection1 = reprojection1 / reprojection1.norm();
            reprojection2 = reprojection2 / reprojection2.norm();
            opengv::bearingVector_t f1 =
                    _adapter.getBearingVector1(camIndex,correspondenceIndex);
            opengv::bearingVector_t f2 =
                    _adapter.getBearingVector2(camIndex,correspondenceIndex);

            //bearing-vector based outlier criterium (select threshold accordingly):
            //1-(f1'*f2) = 1-cos(alpha) \in [0:2]
            double reprojError1 = 1.0 - (f1.transpose() * reprojection1);
            double reprojError2 = 1.0 - (f2.transpose() * reprojection2);
            scores[camIndex].push_back(reprojError1 + reprojError2);
      //      if (correspondenceIndex < 30)
      //          std::cout<<"correspondenceIndex ID: "<<correspondenceIndex<<" reprojection error: "<<reprojError1 + reprojError2<<std::endl;
      //      std::cout << "camIndex: " << camIndex <<" bv index "<<correspondenceIndex<< " reprojection error: " << reprojError1 + reprojError2 << std::endl;
        }
   //     std::cout<<"reprojection error: "<<indices[camIndex].size()<<std::endl;
    }
    _adapter.sett12(init_t);
    _adapter.setR12(init_R);
}

void
polyview::demos::EigensolverSacProblem::optimizeModelCoefficients(
        const std::vector<std::vector<int> > &inliers,
        const model_t &model,
        model_t &optimized_model) {
    optimized_model = model; //todo: include non-linear optimization of model

    opengv::translation_t translation = model.col(3);
    opengv::rotation_t rotation = model.block<3, 3>(0, 0);

/*    _adapter.sett12(translation);
    _adapter.setR12(rotation);*/

/*      optimized_model = opengv::relative_pose::optimize_nonlinear(
              _adapter,_adapter.convertMultiIndices(inliers));*/

    polyview::demos::EigensolverSacProblem::computeModelCoefficients(inliers, optimized_model);

//    translation = optimized_model.col(3);
//    rotation = optimized_model.block<3, 3>(0, 0);

//    _adapter.sett12(translation);
//    _adapter.setR12(rotation);


/*    std::vector<Eigen::Matrix3d> xxF, yyF, zzF, xyF, yzF, zxF;

    //Initial summation terms
    for (size_t k = 0; k < inliers.size(); k++) {
        xxF.push_back(Eigen::Matrix3d::Zero());
        yyF.push_back(Eigen::Matrix3d::Zero());
        zzF.push_back(Eigen::Matrix3d::Zero());
        xyF.push_back(Eigen::Matrix3d::Zero());
        yzF.push_back(Eigen::Matrix3d::Zero());
        zxF.push_back(Eigen::Matrix3d::Zero());

        //Fill summation terms
        for (size_t i = 0; i < inliers[k].size(); i++) {
            opengv::bearingVector_t f1 = _adapter.getCamRotation(k) * _adapter.getBearingVector1(k, inliers[k][i]);
            opengv::bearingVector_t f2 = _adapter.getCamRotation(k) * _adapter.getBearingVector2(k, inliers[k][i]);
            Eigen::Matrix3d F = f2 * f2.transpose();

            double weight = 1.0;

            xxF[k] = xxF[k] + weight * f1[0] * f1[0] * F;
            yyF[k] = yyF[k] + weight * f1[1] * f1[1] * F;
            zzF[k] = zzF[k] + weight * f1[2] * f1[2] * F;
            xyF[k] = xyF[k] + weight * f1[0] * f1[1] * F;
            yzF[k] = yzF[k] + weight * f1[1] * f1[2] * F;
            zxF[k] = zxF[k] + weight * f1[2] * f1[0] * F;
        }
    }

    //randomize the starting point a bit
    rotation = _adapter.getR12();

    //Do minimization
    const int n = 3;
    int rotation_type = -1;
    Eigen::VectorXd x(n);
    Eigen::VectorXd x_planar(1);

    x = opengv::math::rot2cayley(rotation);
//    std::cout<<"x"<<x<<std::endl;

    if (x[0] != 0 && x[1] == 0 && x[2] == 0) {
        rotation_type = 1;
        x_planar << x[0];
    } else if (x[0] == 0 && x[1] != 0 && x[2] == 0) {
        rotation_type = 2;
        x_planar << x[1];
    } else if (x[0] == 0 && x[1] == 0) {
        rotation_type = 3;
        x_planar << x[2];
    }

    rotation_type = 3;
    //PLANAR MOTION
    Eigensolver_optimize functor(xxF, yyF, zzF, xyF, yzF, zxF, _adapter, inliers, rotation_type);
    Eigen::NumericalDiff<Eigensolver_optimize> numDiff(functor);
    Eigen::LevenbergMarquardt<Eigen::NumericalDiff<Eigensolver_optimize>> lm(numDiff);

    lm.resetParameters();
    lm.parameters.ftol = 0.00005;
    lm.parameters.xtol = 1.E1 * NumTraits<double>::epsilon();
    lm.parameters.maxfev = 100;
    lm.minimize(x_planar);

    if ((inliers[0].size() + inliers[1].size() + inliers[2].size() + inliers[3].size()) > 0) {
        std::ofstream output1,output2;
        std::string folder1("/home/ifwang/datasets/energy_test/eigenvalue.txt");
        std::string folder2("/home/ifwang/datasets/energy_test/jacobian.txt");
        output1.precision(6);
        output2.precision(6);
        output1.open(folder1, std::ios::trunc);
        output2.open(folder2, std::ios::trunc);

        for (double theta = -M_PI/2; theta < M_PI/2; theta = theta + 0.001) {
            double eigen_energy = 0;
            double jacobian_energy = 0;
            for (int i = 0; i < xxF.size(); i++) {

                opengv::cayley_t cayley;
                opengv::rotation_t rel_r;
                rel_r << cos(theta), -sin(theta), 0,
                        sin(theta), cos(theta), 0,
                        0, 0, 1;
                cayley = opengv::math::rot2cayley(rel_r);

                double eigenvalue,eigenvalue2;
                Eigen::Matrix<double, 1, 3> jacobian;
                eigenvalue = opengv::relative_pose::modules::eigensolver::getSmallestEVwithJacobian(
                        xxF[i], yyF[i], zzF[i], xyF[i], yzF[i], zxF[i], cayley, jacobian);
                eigenvalue2 = polyview::demos::EigensolverSacProblem::getSmallestEVwithJacobian(
                        _adapter, inliers, i, cayley, jacobian);
                output1 << pow(eigenvalue,2) << " ";
                output2 << pow(eigenvalue2,1) << " ";
                eigen_energy = eigen_energy + pow(eigenvalue,2);
                jacobian_energy = jacobian_energy + pow(eigenvalue2,1);
            }
            output1 << eigen_energy << " ";
            output2 << jacobian_energy << " ";
            output1 << std::endl;
            output2 << std::endl;
        }
    }

    switch (rotation_type) {
        case 1:
            x[0] = x_planar[0];
            break;
        case 2:
            x[1] = x_planar[0];
            break;
        case 3:
            x[2] = x_planar[0];
            break;
        default:
            std::cout << "Not A planar Motion!" << std::endl;
    }

    opengv::cayley_t cayley_computed = x;
    //   std::cout<<"x optimized"<<x<<std::endl;
    opengv::rotation_t R = opengv::math::cayley2rot(cayley_computed);
    opengv::translation_t t;

    Eigen::MatrixXd LargeMatrix = Eigen::MatrixXd::Zero(inliers.size() * 3, inliers.size() + 3);
    Eigen::VectorXd Entries = Eigen::VectorXd::Zero(inliers.size() * 3);

//    LargeMatrix.block(0,indices.size(), indices.size() * 3,3) = - Eigen::MatrixXd::Ones( indices.size() * 3, 3);
    //GET absolute R and t with scale
    for (int i = 0; i < inliers.size(); i++) {

        Eigen::Matrix3d M = opengv::relative_pose::modules::eigensolver::composeM(xxF[i], yyF[i], zzF[i], xyF[i],
                                                                                  yzF[i], zxF[i], cayley_computed);
        Eigen::EigenSolver<Eigen::Matrix3d> Eig(M, true);
        Eigen::Matrix<std::complex<double>, 3, 1> D_complex = Eig.eigenvalues();
        Eigen::Matrix<std::complex<double>, 3, 3> V_complex = Eig.eigenvectors();
        opengv::eigenvalues_t D;
        opengv::eigenvectors_t V;
        for (size_t k = 0; k < 3; k++) {
            D[k] = D_complex[k].real();
            for (size_t j = 0; j < 3; j++)
                V(k, j) = V_complex(k, j).real();
        }

        int index;
        double temp;
        index = 0;
        temp = D[0];

        for (int k = 0; k < 3; k++) {
            if (D[k] < temp && D[k] != 0) {
                temp = D[k];
                index = k;
            }
        }
        //    double translationMagnitude = sqrt(pow(D[1], 2) + pow(D[2], 2));
        opengv::translation_t t = V.col(index);

        //Correct the translation
        //    std::cout<<"indices size()"<<indices[i].size()<<std::endl;

        size_t pick;
        for (size_t index = 0; index < inliers.size(); index++) {
            if (inliers[index].size() != 0) {
                pick = index;
                break;
            }
        }

        opengv::bearingVector_t f1 = _adapter.getCamRotation(pick) * _adapter.getBearingVector1(i, inliers[pick][0]);
        opengv::bearingVector_t f2 = _adapter.getCamRotation(pick) * _adapter.getBearingVector2(i, inliers[pick][0]);
        f2 = R * f2;
        Eigen::Vector3d opticalFlow = f1 - f2;
        if (opticalFlow.dot(t) < 0.0)
            t = -t;

        Eigen::Matrix3d I = Eigen::Matrix3d::Identity(3, 3);

        LargeMatrix.block(i * 3, i, 3, 1) = t;
        LargeMatrix.block(i * 3, inliers.size(), 3, 3) = -I;
        Entries.segment(i * 3, 3) = (R - I) * _adapter.getCamOffset(i);
    }

    Eigen::VectorXd solution = LargeMatrix.colPivHouseholderQr().solve(Entries);
    t = solution.segment(inliers.size(), 3);
    //  std::cout<<"centre T:"<<t.transpose()<<std::endl;

    optimized_model.col(3) = t;
    optimized_model.block<3, 3>(0, 0) = R;*/
}


/*std::vector<int>
polyview::demos::EigensolverSacProblem::getSampleSize() const {
    std::vector<int> sampleSizes;
    for( size_t i = 0; i < _adapter.getNumberPairs(); i++ )
        sampleSizes.push_back(_sampleSize);
    return sampleSizes;
}*/

std::vector<int>
polyview::demos::EigensolverSacProblem::getSampleSizes() const {
    std::vector<int> sampleSizes;
    for (size_t i = 0; i < _adapter.getNumberPairs(); i++)
        sampleSizes.push_back(0);

    int sampleSize = _sampleSize;

    //set it to a random cam index where to start
    size_t binIndex = floor(((double) _adapter.getNumberPairs()) *
                            ((double) rand()) / ((double) RAND_MAX));
    for (int i = 0; i < sampleSize; i++) {
        sampleSizes[binIndex]++;
        binIndex++;
        if (binIndex >= sampleSizes.size())
            binIndex = 0;
    }
/*    sampleSizes[0] = 3;
    sampleSizes[1] = 3;
    sampleSizes[2] = 3;
    sampleSizes[3] = 3;*/
//       std::cout<<"sampleSizes: "<<sampleSizes[0]<<" "<<sampleSizes[1]<<" "<<sampleSizes[2]<<" "<<sampleSizes[3]<<std::endl;
    return sampleSizes;
}