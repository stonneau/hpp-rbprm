// Copyright (c) 2014, LAAS-CNRS
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
// hpp-rbprm. If not, see <http://www.gnu.org/licenses/>.

#include <hpp/rbprm/interpolation/rbprm-path-interpolation.hh>

#ifdef PROFILE
    #include "hpp/rbprm/rbprm-profiler.hh"
#endif

namespace hpp {
  namespace rbprm {
    namespace interpolation {

    RbPrmInterpolationPtr_t RbPrmInterpolation::create (const hpp::rbprm::RbPrmFullBodyPtr_t robot,
                                                        const hpp::rbprm::State &start, const hpp::rbprm::State &end,
                                                        const core::PathVectorConstPtr_t path)
    {
        RbPrmInterpolation* rbprmDevice = new RbPrmInterpolation(path, robot, start, end);
        RbPrmInterpolationPtr_t res (rbprmDevice);
        res->init (res);
        return res;
    }

    RbPrmInterpolation::~RbPrmInterpolation()
    {
        // NOTHING
    }

    // ========================================================================

    namespace
    {
    core::Configuration_t configPosition(core::ConfigurationIn_t previous, const core::PathVectorConstPtr_t path, double i)
    {
        core::Configuration_t configuration = previous;
        const core::Configuration_t configPosition = path->operator ()(std::min(i, path->timeRange().second));
        configuration.head(configPosition.rows()) = configPosition;
        return configuration;
    }
    }

    std::vector<State> RbPrmInterpolation::Interpolate(const affMap_t& affordances,
			const std::map<std::string, std::vector<std::string> >& affFilters,
			const double timeStep, const double robustnessTreshold)
    {
        if(!path_) throw std::runtime_error ("Cannot interpolate; no path given to interpolator ");
        std::vector<model::Configuration_t> configs;
        const core::interval_t& range = path_->timeRange();
        configs.push_back(start_.configuration_);
        for(double i = range.first + timeStep; i< range.second; i+= timeStep)
        {
            configs.push_back(configPosition(configs.back(),path_,i));
        }
        return Interpolate(affordances, affFilters, configs, robustnessTreshold);
    }

    std::vector<State> RbPrmInterpolation::Interpolate(const affMap_t& affordances,
			const std::map<std::string, std::vector<std::string> >& affFilters,
			const std::vector<model::Configuration_t>& configs,
			const double robustnessTreshold)
    {
        int nbFailures = 0;
        std::vector<State> states;
        states.push_back(this->start_);
        std::size_t nbRecontacts = 0;
        bool allowFailure = true;
#ifdef PROFILE
    RbPrmProfiler& watch = getRbPrmProfiler();
    watch.reset_all();
    watch.start("complete generation");
#endif
        for(std::vector<model::Configuration_t>::const_iterator cit = configs.begin()+1; cit != configs.end(); ++cit)
        {
            const State& previous = states.back();
            core::Configuration_t configuration = *cit;
            Eigen::Vector3d dir = configuration.head<3>() - previous.configuration_.head<3>();
            fcl::Vec3f direction(dir[0], dir[1], dir[2]);
            bool nonZero(false);
            direction.normalize(&nonZero);
            if(!nonZero) direction = fcl::Vec3f(0,0,1.);
            // TODO Direction 6d
            bool sameAsPrevious(true);
            bool multipleBreaks(false);
            State newState = ComputeContacts(previous, robot_,configuration,
						affordances,affFilters,direction,
							sameAsPrevious, multipleBreaks,allowFailure,robustnessTreshold);
            if(allowFailure && multipleBreaks)
            {
                ++ nbFailures;
                ++cit;
if (nbFailures > 1)
{
#ifdef PROFILE
    watch.stop("complete generation");
    watch.add_to_count("planner failed", 1);
    std::ofstream fout;
    fout.open("log.txt", std::fstream::out | std::fstream::app);
    std::ostream* fp = &fout;
    watch.report_count(*fp);
    fout.close();
#endif
    return states;
}
            }
            if(multipleBreaks && !allowFailure)
            {
                ++nbRecontacts;
                cit--;
            }
            else
            {
                nbRecontacts = 0;
            }
            if(sameAsPrevious)
            {
                states.pop_back();
            }
            newState.nbContacts = newState.contactNormals_.size();
            states.push_back(newState);
            allowFailure = nbRecontacts > robot_->GetLimbs().size() + 6;
        }
        states.push_back(this->end_);
#ifdef PROFILE
        watch.add_to_count("planner succeeded", 1);
        watch.stop("complete generation");
        std::ofstream fout;
        fout.open("log.txt", std::fstream::out | std::fstream::app);
        std::ostream* fp = &fout;
        watch.report_all_and_count(2,*fp);
        fout.close();
#endif
        return states;
    }

    void RbPrmInterpolation::init(const RbPrmInterpolationWkPtr_t& weakPtr)
    {
        weakPtr_ = weakPtr;
    }

    RbPrmInterpolation::RbPrmInterpolation (const core::PathVectorConstPtr_t path, const hpp::rbprm::RbPrmFullBodyPtr_t robot, const hpp::rbprm::State &start, const hpp::rbprm::State &end)
        : path_(path)
        , start_(start)
        , end_(end)
        , robot_(robot)
    {
        // TODO
    }
    } // interpolation
  } // rbprm
} //hpp
