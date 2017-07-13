//
// Copyright (c) 2017 CNRS
// Authors: Steve Tonneau
//
// This file is part of hpp-rbprm
// hpp-core is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// hpp-core is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// General Lesser Public License for more details.  You should have
// received a copy of the GNU Lesser General Public License along with
// hpp-core  If not, see
// <http://www.gnu.org/licenses/>.

#include <hpp/rbprm/contact_generation/contact_generation.hh>
#include <hpp/rbprm/contact_generation/algorithm.hh>
#include <hpp/rbprm/stability/stability.hh>
#include <hpp/rbprm/tools.hh>
# include <hpp/rbprm/rbprm-state.hh>
# include <hpp/rbprm/rbprm-fullbody.hh>



namespace hpp {
namespace rbprm {
namespace contact{

ContactReport::ContactReport()
    : projection::ProjectionReport()
    , contactMaintained_(false)
    , multipleBreaks_(false)
    , contactCreated_(false)
    , repositionedInPlace_(false)
{
    // NOTHING
}

ContactReport::ContactReport(const projection::ProjectionReport& parent)
    : projection::ProjectionReport(parent)
    , contactMaintained_(false)
    , multipleBreaks_(false)
    , contactCreated_(false)
    , repositionedInPlace_(false)
{
    // NOTHING
}

ContactReport generateContactReport(const projection::ProjectionReport& parent, const ContactGenHelper& helper, bool repositionedInPlace=false)
{
    const State& previous = helper.previousState_;
    const State& result = parent.result_;
    ContactReport report(parent) ;
    report.contactCreated_ = (result.fixedContacts(previous).size() == previous.nbContacts);
    report.multipleBreaks_ = (result.contactBreaks(previous).size() > helper.maxContactBreaks_);
    report.repositionedInPlace_ = repositionedInPlace;
    report.contactMaintained_ = !repositionedInPlace && !(result.contactCreations(previous).size() > 0);
    return report;
}

projection::ProjectionReport genContactFromOneMaintainCombinatorial(ContactGenHelper& helper)
{
    // retrieve the first feasible result of maintain combinatorial...
    projection::ProjectionReport rep = contact::maintain_contacts(helper);
    if(rep.success_)
    {
        // ... if found, then try to generate feasible contact for this combinatorial.
        helper.workingState_ = rep.result_;
        return gen_contacts(helper);
    }
    return rep;
}

// if contact generation failed, tries to reposition the contacts without moving the root
ContactReport handleFailure(ContactGenHelper& helper)
{
    helper.workingState_ = helper.previousState_;
    projection::ProjectionReport rep = repositionContacts(helper);

    std::vector<std::string> breaks = rep.result_.contactBreaks(helper.previousState_);
    std::vector<std::string> creations = rep.result_.contactCreations(helper.previousState_);
    return generateContactReport(rep,helper,true);
}

ContactReport oneStep(ContactGenHelper& helper)
{
    projection::ProjectionReport rep;
    do
        rep = genContactFromOneMaintainCombinatorial(helper);
    while(!rep.success_ && !helper.candidates_.empty());
    if(!rep.success_) // TODO only possible in quasi static
        return handleFailure(helper);
    return generateContactReport(rep,helper);
}



bool ContactExistsWithinGroup(const hpp::rbprm::RbPrmLimbPtr_t& limb,
                              const hpp::rbprm::RbPrmFullBody::T_LimbGroup& limbGroups,
                              const State& current)
{
    const std::vector<std::string>& group = limbGroups.at(limb->limb_->name());
    for(std::vector<std::string>::const_iterator cit = group.begin();
        cit != group.end(); ++cit)
    {
        if(current.contactPositions_.find(*cit) != current.contactPositions_.end())
            return true;
    }
    return false;
}

ContactComputationStatus ComputeStableContact(const hpp::rbprm::RbPrmFullBodyPtr_t& body,
                          State& current,
                          core::CollisionValidationPtr_t validation,
                          const std::string& limbId,
                          const hpp::rbprm::RbPrmLimbPtr_t& limb,
                          model::ConfigurationIn_t rbconfiguration,
                          model::ConfigurationOut_t configuration,
                          const affMap_t& affordances,
                          const std::map<std::string, std::vector<std::string> >& affFilters,
                          const fcl::Vec3f& direction,
                          fcl::Vec3f& position, fcl::Vec3f& normal, const double robustnessTreshold,
                          bool contactIfFails = true, bool stableForOneContact = true,
                          const sampling::heuristic evaluate = 0)
{
    contact::ContactGenHelper contactGenHelper(body,current,current.configuration_,affordances,affFilters,robustnessTreshold,1,1,false,true,
                                      direction,fcl::Vec3f(0,0,0),contactIfFails,stableForOneContact);

    hpp::rbprm::projection::ProjectionReport rep = contact::generate_contact(contactGenHelper,limbId,evaluate);
    current = rep.result_;
    configuration = rep.result_.configuration_;
    if(rep.status_ != NO_CONTACT)
    {
        position = rep.result_.contactPositions_[limbId];
        normal = rep.result_.contactNormals_[limbId];
    }
    if(rep.status_ == STABLE_CONTACT)
        current.stable = true;
    return rep.status_;
}

hpp::rbprm::State ComputeContacts(const hpp::rbprm::RbPrmFullBodyPtr_t& body,
        model::ConfigurationIn_t configuration, const affMap_t& affordances,
  const std::map<std::string, std::vector<std::string> >& affFilters, const fcl::Vec3f& direction,
        const double robustnessTreshold)
{
    const T_Limb& limbs = body->GetLimbs();
    const rbprm::RbPrmFullBody::T_LimbGroup& limbGroups = body->GetGroups();
    const std::map<std::string, core::CollisionValidationPtr_t>& limbcollisionValidations = body->GetLimbCollisionValidation();
    State result;
    // save old configuration
    core::ConfigurationIn_t save = body->device_->currentConfiguration();
    result.configuration_ = configuration;
    body->device_->currentConfiguration(configuration);
    body->device_->computeForwardKinematics();
    for(T_Limb::const_iterator lit = limbs.begin(); lit != limbs.end(); ++lit)
    {
        if(!ContactExistsWithinGroup(lit->second, limbGroups ,result))
        {
            fcl::Vec3f normal, position;
            ComputeStableContact(body,result,
                                limbcollisionValidations.at(lit->first), lit->first,
                                lit->second, configuration, result.configuration_, affordances,affFilters,
                                direction, position, normal, robustnessTreshold, true, true);
        }
        result.nbContacts = result.contactNormals_.size();
    }
    // reload previous configuration
    body->device_->currentConfiguration(save);
    return result;
}

hpp::rbprm::contact::ContactReport ComputeContacts(const hpp::rbprm::State& previous,
        const hpp::rbprm::RbPrmFullBodyPtr_t& body,
        model::ConfigurationIn_t configuration, const affMap_t& affordances,
        const std::map<std::string, std::vector<std::string> >& affFilters,
        const fcl::Vec3f& direction, const double robustnessTreshold)
{
    // save old configuration
    core::ConfigurationIn_t save = body->device_->currentConfiguration();
    model::Device::Computation_t flag = body->device_->computationFlag ();
    model::Device::Computation_t newflag = static_cast <model::Device::Computation_t> (model::Device::JOINT_POSITION);
    // load new root position
    body->device_->controlComputation (newflag);
    body->device_->currentConfiguration(configuration);
    body->device_->computeForwardKinematics ();
    // try to maintain previous contacts
    contact::ContactGenHelper cHelper(body,previous,configuration,affordances,affFilters,robustnessTreshold,1,1,false,
                                      true,direction,fcl::Vec3f(0,0,0),false,false);
    contact::ContactReport rep = contact::oneStep(cHelper);

    // copy extra dofs
    if(rep.success_)
    {
        const model::size_type& extraDim = body->device_->extraConfigSpace().dimension();
        rep.result_.configuration_.tail(extraDim) = configuration.tail(extraDim);
    }
    body->device_->currentConfiguration(save);
    body->device_->controlComputation (flag);
    return rep;
}

} // namespace projection
} // namespace rbprm
} // namespace hpp
