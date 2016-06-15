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

#ifndef HPP_RBPRM_PATH_INTERPOLATION_HH
# define HPP_RBPRM_PATH_INTERPOLATION_HH

# include <hpp/rbprm/config.hh>
# include <hpp/rbprm/rbprm-fullbody.hh>
# include <hpp/core/path-vector.hh>
# include <hpp/model/device.hh>

# include <vector>

namespace hpp {
  namespace rbprm {
    namespace interpolation {
    HPP_PREDEF_CLASS(RbPrmInterpolation);

    /// Interpolation class for transforming a path computed by RB-PRM into
    /// a discrete sequence of balanced contact configurations.
    ///
    class RbPrmInterpolation;
    typedef boost::shared_ptr <RbPrmInterpolation> RbPrmInterpolationPtr_t;

    class HPP_RBPRM_DLLAPI RbPrmInterpolation
    {
    public:
        /// Creates a smart pointer to the Interpolation class
        ///
        /// \param path the path returned by RB-PRM computation
        /// \param robot the FullBody instance considered for extending the part
        /// \param start the start full body configuration of the problem
        /// \param end the end full body configuration of the problem
        /// \return a pointer to the created RbPrmInterpolation instance
        static RbPrmInterpolationPtr_t create (const RbPrmFullBodyPtr_t robot, const State& start, const State& end,
                                               const core::PathVectorConstPtr_t path = core::PathVectorConstPtr_t());

    public:
        ~RbPrmInterpolation();

        /// Transforms the path computed by RB-PRM into
        /// a discrete sequence of balanced contact configurations.
        ///
        /// \param affordances the set of 3D objects to consider for contact creation.
        /// \param affFilters a vector of strings determining which affordance
        ///  types are to be used in generating contacts for each limb.
        /// \param timeStep the discretization step of the path.
        /// \param robustnessTreshold minimum value of the static equilibrium robustness criterion required to accept the configuration (0 by default).
        /// \return a pointer to the created RbPrmInterpolation instance
        std::vector<State> Interpolate(const affMap_t& affordances,
																	  	 const std::map<std::string, std::vector<std::string> >& affFilters,
                                       const double timeStep = 0.01, const double robustnessTreshold=0.);

        /// Transforms a discrete sequence of configurations into
        /// a discrete sequence of balanced contact configurations.
        ///
        /// \param affordances the set of 3D objects to consider for contact creation.
        /// \param affFilters a vector of strings determining which affordance
        ///  types are to be used in generating contacts for each limb.
        /// \param configs
        /// \param robustnessTreshold minimum value of the static equilibrium robustness criterion required to accept the configuration (0 by default).
        /// \return a pointer to the created RbPrmInterpolation instance
        std::vector<State> Interpolate(const affMap_t& affordances,
																			 const std::map<std::string, std::vector<std::string> >& affFilters,
                                       const std::vector<model::Configuration_t>& configs, const double robustnessTreshold=0.);

    public:
        const core::PathVectorConstPtr_t path_;
        const State start_;
        const State end_;

    private:
        RbPrmFullBodyPtr_t robot_;

    protected:
      RbPrmInterpolation (const core::PathVectorConstPtr_t path,
				const RbPrmFullBodyPtr_t robot,const State& start, const State& end);

      ///
      /// \brief Initialization.
      ///
      void init (const RbPrmInterpolationWkPtr_t& weakPtr);

    private:
      RbPrmInterpolationWkPtr_t weakPtr_;
    }; // class RbPrmLimb
    } // namespace interpolation
  } // namespace rbprm
} // namespace hpp

#endif // HPP_RBPRM_PATH_INTERPOLATION_HH
