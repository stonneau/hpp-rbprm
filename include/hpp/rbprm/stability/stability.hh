//
// Copyright (c) 2014 CNRS
// Authors: Steve Tonneau (steve.tonneau@laas.fr)
//
// This file is part of hpp-rbprm.
// hpp-rbprm is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// hpp-rbprm is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Lesser Public License for more details.  You should have
// received a copy of the GNU Lesser General Public License along with
// hpp-core  If not, see
// <http://www.gnu.org/licenses/>.

#ifndef HPP_RBPRM_STABILITY_HH
# define HPP_RBPRM_STABILITY_HH

#include <hpp/model/device.hh>
#include <hpp/rbprm/rbprm-state.hh>
#include <hpp/rbprm/rbprm-fullbody.hh>
#include <robust-equilibrium-lib/static_equilibrium.hh>

#include <map>
#include <memory>

namespace hpp {

  namespace rbprm {
  namespace stability{


    typedef Eigen::Matrix <model::value_type, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MatrixXX;
    typedef Eigen::Matrix <model::value_type, Eigen::Dynamic, 1>                               VectorX;

    /// Using the polytope computation of the gravito inertial wrench cone, performs
    /// a static equilibrium test on the robot.
    ///
    /// \param fullbody The considered robot for static equilibrium
    /// \param state The current State of the robots, in terms of contact creation
    /// \return Whether the configuration is statically balanced
    double IsStable(const RbPrmFullBodyPtr_t fullbody, State& state, const robust_equilibrium::StaticEquilibriumAlgorithm = robust_equilibrium::STATIC_EQUILIBRIUM_ALGORITHM_DLP);


    /// Using the polytope computation of the gravito inertial wrench cone,
    /// returns the CWC of the robot at a given state
    ///
    /// \param fullbody The considered robot for static equilibrium
    /// \param state The current State of the robots, in terms of contact creation
    std::pair<MatrixXX, VectorX> ComputeCentroidalCone(const RbPrmFullBodyPtr_t fullbody, State& state, const core::value_type friction = 0.5);

    /* --- GCR criterion --- */
    /// Data structure to define a plane corresponding to the following equation : ax + by + cz + d = 0
    struct Plane
    {
      double a;
      double b;
      double c;
      double d;
      Plane() : a(0), b(0), c(1), d(0) {}
      Plane(double aa, double bb, double cc, double dd) : a(aa), b(bb), c(cc), d(dd) {}
      Plane(const Plane & pe) : a(pe.a), b(pe.b), c(pe.c), d(pe.d) {}
      Plane & operator=(const Plane & pe);
    };

    /// Data structure to store 2-dimensional informations (2D vectors)
    struct Vec2D
    {
      double x;
      double y;
      Vec2D() : x(0), y(0) {}
      Vec2D(double xx, double yy) : x(xx), y(yy) {}
      Vec2D(const Vec2D & c2D) : x(c2D.x), y(c2D.y) {}
      Vec2D & operator=(const Vec2D & c);
    };

    /// Data structure to store the parameters of the GCR criterion
    struct GCRParam
    {
      double edge_; // supposed optimal size for the edges of the support polygon
      double groundThreshold_; // distance from the ground plane under which we will consider a contact as a ground contact (and perform its orthogonal projection in the ground plane)
      double step1Threshold_; // maximum error accepted on the distance between the CoM projection and the support polygon center
      double step2Threshold_; // maximum error accepted on the distance of each limb to the support polygon center
      GCRParam() : edge_(-1.0), groundThreshold_(-1.0), step1Threshold_(-1.0), step2Threshold_(-1.0) {} // negative values in order to say that the attributes have to be settled
      GCRParam(double e, double g, double s1, double s2) : edge_(e), groundThreshold_(g), step1Threshold_(s1), step2Threshold_(s2) {}
      GCRParam(const GCRParam & p) : edge_(p.edge_), groundThreshold_(p.groundThreshold_), step1Threshold_(p.step1Threshold_), step2Threshold_(p.step2Threshold_) {}
      GCRParam & operator=(const GCRParam & p);
    };

    /// Performs the 2-steps Ground Contacts Relevance criterion validation
    ///
    /// \param state The considered state of the robot
    /// \param comPosition The center of mass position of the robot
    /// \param config The parameters of the criterion validation
    /// \return true if steps 1 and 2 are successful, false otherwise
    bool GCRValidation(const hpp::rbprm::State & state, const fcl::Vec3f & comPosition, const GCRParam & config);

    /// Performs the 2-steps evaluation in order to quantify the quality of a state according to the GCR criterion
    ///
    /// \param state The considered state of the robot
    /// \param comPosition The center of mass position of the robot
    /// \param config The parameters of the criterion evaluation (only config.edge and config.groundThreshold will be used for the evaluation)
    /// \return The cost of the considered state (the lower the cost is, the better the state is --> cost to minimize)
    double GCREvaluation(const hpp::rbprm::State & state, const fcl::Vec3f & comPosition, const GCRParam & config);

    /* --- End GCR criterion --- */

  } // namespace stability
} // namespace rbprm
} // namespace hpp
#endif // HPP_RBPRM_STABILITY_HH
