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
#include <hpp/rbprm/stability/stability.hh>
#include <hpp/rbprm/tools.hh>
#ifdef PROFILE
    #include "hpp/rbprm/rbprm-profiler.hh"
#endif


namespace hpp {
namespace rbprm {
namespace contact{


ContactGenHelper::ContactGenHelper(RbPrmFullBodyPtr_t fb, const State& ps, model::ConfigurationIn_t configuration,
                                    const hpp::rbprm::affMap_t &affordances, const std::map<std::string, std::vector<std::string> > &affFilters,
                                    const double robustnessTreshold,
                                    const std::size_t maxContactBreaks, const std::size_t maxContactCreations,
                                    const bool checkStabilityMaintain, const bool checkStabilityGenerate,
                                    const fcl::Vec3f& direction,
                                    const fcl::Vec3f& acceleration,
                                    const bool contactIfFails,
                                    const bool stableForOneContact)
: fullBody_(fb)
, previousState_(ps)
, checkStabilityMaintain_(checkStabilityMaintain)
, contactIfFails_(contactIfFails)
, stableForOneContact_(stableForOneContact)
, acceleration_(acceleration)
, direction_(direction)
, robustnessTreshold_(robustnessTreshold)
, maxContactBreaks_(maxContactBreaks)
, maxContactCreations_(maxContactCreations)
, affordances_(affordances)
, affFilters_(affFilters)
, workingState_(previousState_)
, checkStabilityGenerate_(checkStabilityGenerate)
{
    workingState_.configuration_ = configuration;
    workingState_.stable = false;
}

typedef std::vector<T_State > T_DepthState;

bool push_if_new(T_State& states, const State currentState)
{
    for(CIT_State cit = states.begin(); cit != states.end(); ++cit)
    {
        if(currentState.contactOrder_== cit->contactOrder_)
            return false;
    }
    states.push_back(currentState);
    return true;
}

void maintain_contacts_combinatorial_rec(const hpp::rbprm::State& currentState, const std::size_t  depth,
                                         const std::size_t maxBrokenContacts, T_DepthState& res)
{
    if (!push_if_new(res[depth], currentState) || depth>=maxBrokenContacts) return;
    std::queue<std::string> contactOrder = currentState.contactOrder_;
    int size = contactOrder.size(); int i = 0;
    while(!contactOrder.empty() && size != i)
    {
        hpp::rbprm::State copyState = currentState;
std::vector<std::string> fixed = currentState.fixedContacts(currentState);
        const std::string contactRemoved = contactOrder.front();
//if(!
//((std::find(fixed.begin(), fixed.end(),std::string("hrp2_rleg_rom")) == fixed.end() && contactRemoved == std::string("hrp2_lleg_rom")) ||
//(std::find(fixed.begin(), fixed.end(),std::string("hrp2_lleg_rom")) == fixed.end() && contactRemoved == std::string("hrp2_rleg_rom"))))
{
        copyState.RemoveContact(contactRemoved);
        maintain_contacts_combinatorial_rec(copyState, depth+1, maxBrokenContacts, res);
}
/*else
{
 std::cout << "avoided both leg removed"    << std::endl;
 contactOrder.push(contactRemoved);
}*/
++i;
contactOrder.pop();
    }
}

Q_State flatten(const T_DepthState& depthstates)
{
    Q_State res;
    for(T_DepthState::const_iterator cit = depthstates.begin(); cit != depthstates.end(); ++cit)
    {
        for(CIT_State ccit = cit->begin(); ccit != cit->end(); ++ccit)
            res.push(*ccit);
    }
    return res;
}

Q_State maintain_contacts_combinatorial(const hpp::rbprm::State& currentState, const std::size_t maxBrokenContacts)
{
    T_DepthState res(maxBrokenContacts+1);
    maintain_contacts_combinatorial_rec(currentState, 0, maxBrokenContacts,res);
    return flatten(res);
}

using namespace projection;

bool maintain_contacts_stability_rec(hpp::rbprm::RbPrmFullBodyPtr_t fullBody,
                        model::ConfigurationIn_t targetRootConfiguration,
                        Q_State& candidates,const std::size_t contactLength,
                        const fcl::Vec3f& acceleration, const double robustness,
                        ProjectionReport& currentRep)
{
    if(stability::IsStable(fullBody,currentRep.result_, acceleration) > robustness)
    {
        currentRep.result_.stable = true;
        return true;
    }
    currentRep.result_.stable = false;
    if(!candidates.empty())
    {
        State cState = candidates.front();
        candidates.pop();
         // removed more contacts, cannot be stable if previous state was not
        if(cState.contactOrder_.size() < contactLength) return false;
        ProjectionReport rep = projectToRootConfiguration(fullBody,targetRootConfiguration,cState);
        Q_State copy_candidates = candidates;
        if(maintain_contacts_stability_rec(fullBody,targetRootConfiguration,copy_candidates,contactLength,acceleration, robustness, rep))
        {
            currentRep = rep;
            candidates = copy_candidates;
            return true;
        }
    }
    return false;
}


hpp::model::ObjectVector_t getAffObjectsForLimb(const std::string& limb,
    const affMap_t& affordances, const std::map<std::string, std::vector<std::string> >& affFilters)
{
    model::ObjectVector_t affs;
    std::vector<std::string> affTypes;
    bool settingFound = false;
    for (std::map<std::string, std::vector<std::string> >::const_iterator fIt =
        affFilters.begin (); fIt != affFilters.end (); ++fIt)
    {
        std::size_t found = fIt->first.find(limb);
        if (found != std::string::npos)
        {
            affTypes = fIt->second;
            settingFound = true;
            break;
        }
    }
    if (!settingFound)
    {
        // TODO: Keep warning or delete it?
        std::cout << "No affordance filter setting found for limb " << limb
            << ". Has such filter been set?" << std::endl;
        // Use all AFF OBJECTS as default if no filter setting exists
        for (affMap_t::const_iterator affordanceIt = affordances.begin ();
            affordanceIt != affordances.end (); ++affordanceIt)
        {
            std::copy (affordanceIt->second.begin (), affordanceIt->second.end (), std::back_inserter (affs));
        }
    }
    else
    {
        for (std::vector<std::string>::const_iterator affTypeIt = affTypes.begin ();
            affTypeIt != affTypes.end (); ++affTypeIt)
        {
            affMap_t::const_iterator affIt = affordances.find(*affTypeIt);
            std::copy (affIt->second.begin (), affIt->second.end (), std::back_inserter (affs));
        }
    }
    if (affs.empty())
        throw std::runtime_error ("No aff objects found for limb " + limb);
    return affs;
}

ProjectionReport maintain_contacts_stability(ContactGenHelper &contactGenHelper, ProjectionReport& currentRep)
{
    const std::size_t contactLength(currentRep.result_.contactOrder_.size());
//contactGenHelper.candidates_.pop(); // TODO REMOVE (TEST)
    maintain_contacts_stability_rec(contactGenHelper.fullBody_,
                                    contactGenHelper.workingState_.configuration_,
                                    contactGenHelper.candidates_,
                                    contactLength, contactGenHelper.acceleration_,
                                    contactGenHelper.robustnessTreshold_, currentRep);
    return currentRep;
}


std::vector<std::string> extractEffectorsName(const rbprm::T_Limb& limbs)
{
    std::vector<std::string> res;
    for(rbprm::T_Limb::const_iterator cit = limbs.begin(); cit != limbs.end(); ++cit)
        res.push_back(cit->first);
    return res;
}

std::vector<std::string> sortLimbs(const State& currentState, const std::vector<std::string>& freeLimbs)
{
    std::vector<std::string> res1;
    std::vector<std::string> res2;
    // first add unused limbs
    // then undraw from contact order
    std::queue<std::string> order = currentState.contactOrder_;
    while(!order.empty())
    {
        const std::string l = order.front(); order.pop();
        if(std::find(freeLimbs.begin(), freeLimbs.end(), l) != freeLimbs.end())
            tools::insertIfNew(res1, l);
    }
    for(std::vector<std::string>::const_iterator cit = freeLimbs.begin(); cit!= freeLimbs.end(); ++cit)
    {
        if(std::find(res1.begin(), res1.end(), *cit) == res1.end())
        {
            res2.push_back(*cit);
        }
    }
    res2.insert(res2.end(), res1.begin(), res1.end());
    //res2.insert(res2.begin(), res1.begin(), res1.end());
    return res2;
}

ProjectionReport genColFree(ContactGenHelper &contactGenHelper, ProjectionReport& currentRep)
{
    ProjectionReport res = currentRep;
    // identify broken limbs and find collision free configurations for each one of them.
    std::vector<std::string> effNames(extractEffectorsName(contactGenHelper.fullBody_->GetLimbs()));
    std::vector<std::string> freeLimbs = rbprm::freeEffectors(currentRep.result_,effNames.begin(), effNames.end() );
    freeLimbs = sortLimbs(contactGenHelper.workingState_, freeLimbs);
    for(std::vector<std::string>::const_iterator cit = freeLimbs.begin(); cit != freeLimbs.end() && res.success_; ++cit)
        res = projection::setCollisionFree(contactGenHelper.fullBody_,contactGenHelper.fullBody_->GetLimbCollisionValidation().at(*cit),*cit,res.result_);

    return res;
}

void stringCombinatorialRec(std::vector<std::vector<std::string> >& res, const std::vector<std::string>& candidates, const std::size_t depth)
{
    if(depth == 0) return;
    std::vector<std::vector<std::string> > newStates;
    for(std::vector<std::vector<std::string> >::iterator it = res.begin(); it != res.end(); ++it)
    {
        for(std::vector<std::string>::const_iterator canditates_it = candidates.begin(); canditates_it != candidates.end(); ++canditates_it)
        {
            std::vector<std::string> contacts = *it;
            if(tools::insertIfNew(contacts,*canditates_it))
            {
                newStates.push_back(contacts);
            }
        }
    }
    stringCombinatorialRec(newStates, candidates, depth-1);
    res.insert(res.end(),newStates.begin(),newStates.end());
}


std::vector<std::vector<std::string> > stringCombinatorial(const std::vector<std::string>& candidates, const std::size_t maxDepth)
{
    std::vector<std::vector<std::string> > res;
    std::vector<std::string> tmp;
    res.push_back(tmp);
    stringCombinatorialRec(res, candidates, maxDepth);
    return res;
}

void gen_contacts_combinatorial_rec(const std::vector<std::string>& freeEffectors, const State& previous, T_ContactState& res, const std::size_t maxCreatedContacts)
{
    std::vector<std::vector<std::string> > allNewStates = stringCombinatorial(freeEffectors, maxCreatedContacts);
    for(std::vector<std::vector<std::string> >::const_iterator cit = allNewStates.begin(); cit!=allNewStates.end();++cit)
    {
        ContactState contactState; contactState.first = previous; contactState.second = *cit;
        res.push(contactState);
    }
}

T_ContactState gen_contacts_combinatorial(const std::vector<std::string>& freeEffectors, const State& previous, const std::size_t maxCreatedContacts)
{
    T_ContactState res;;
    gen_contacts_combinatorial_rec(freeEffectors, previous, res, maxCreatedContacts);
    return res;
}

T_ContactState gen_contacts_combinatorial(ContactGenHelper& contactGenHelper)
{
    State& cState = contactGenHelper.workingState_;
    std::vector<std::string> effNames(extractEffectorsName(contactGenHelper.fullBody_->GetLimbs()));
    const std::vector<std::string> freeLimbs = rbprm::freeEffectors(cState,effNames.begin(), effNames.end() );
    return gen_contacts_combinatorial(freeLimbs, cState, contactGenHelper.maxContactCreations_);
}


ProjectionReport maintain_contacts(ContactGenHelper &contactGenHelper)
{
    ProjectionReport rep;
    Q_State& candidates = contactGenHelper.candidates_;
    if(candidates.empty())
        candidates = maintain_contacts_combinatorial(contactGenHelper.workingState_,contactGenHelper.maxContactBreaks_);
    else
        candidates.pop(); // first candidate already treated.
    while(!candidates.empty() && !rep.success_)
    {
        //retrieve latest state
        State cState = candidates.front();
        candidates.pop();
        rep = projectToRootConfiguration(contactGenHelper.fullBody_,contactGenHelper.workingState_.configuration_,cState);
        if(rep.success_)
            rep = genColFree(contactGenHelper, rep);
        if(rep.success_)
        {
            //collision validation
            hpp::core::ValidationReportPtr_t valRep (new hpp::core::CollisionValidationReport);
            rep.success_ = contactGenHelper.fullBody_->GetCollisionValidation()->validate(rep.result_.configuration_, valRep);
        }
    }
    if(rep.success_ && contactGenHelper.checkStabilityMaintain_)
        return maintain_contacts_stability(contactGenHelper, rep);
    return rep;
}


sampling::T_OctreeReport CollideOctree(const ContactGenHelper &contactGenHelper, const std::string& limbName,
                                                    RbPrmLimbPtr_t limb, const sampling::heuristic evaluate)
{
    fcl::Transform3f transform = limb->octreeRoot(); // get root transform from configuration
    hpp::model::ObjectVector_t affordances = getAffObjectsForLimb (limbName,contactGenHelper.affordances_, contactGenHelper.affFilters_);

    //#pragma omp parallel for
    // request samples which collide with each of the collision objects
    sampling::heuristic eval =  evaluate == 0 ? limb->evaluate_ : evaluate;
    std::size_t i (0);
    if (affordances.empty ())
      throw std::runtime_error ("No aff objects found!!!");

    std::vector<sampling::T_OctreeReport> reports(affordances.size());
    for(model::ObjectVector_t::const_iterator oit = affordances.begin();
        oit != affordances.end(); ++oit, ++i)
    {
        if(eval)
            sampling::GetCandidates(limb->sampleContainer_, transform, *oit, contactGenHelper.direction_, reports[i], eval);
        else
            sampling::GetCandidates(limb->sampleContainer_, transform, *oit, contactGenHelper.direction_, reports[i]);
    }
    sampling::T_OctreeReport finalSet;
    // order samples according to EFORT
    for(std::vector<sampling::T_OctreeReport>::const_iterator cit = reports.begin();
        cit != reports.end(); ++cit)
    {
        finalSet.insert(cit->begin(), cit->end());
    }
    return finalSet;
}

hpp::rbprm::State findValidCandidate(const ContactGenHelper &contactGenHelper, const std::string& limbId,
                        RbPrmLimbPtr_t limb, core::CollisionValidationPtr_t validation, bool& found_sample,
                                     bool& unstableContact, const sampling::heuristic evaluate = 0)
{
    State current = contactGenHelper.workingState_;
    current.stable = false;
    sampling::T_OctreeReport finalSet = CollideOctree(contactGenHelper, limbId, limb, evaluate);
    core::Configuration_t moreRobust, configuration;
    configuration = current.configuration_;
    double maxRob = -std::numeric_limits<double>::max();
    sampling::T_OctreeReport::const_iterator it = finalSet.begin();
    fcl::Vec3f position, normal;
    fcl::Matrix3f rotation;
    ProjectionReport rep ;
    for(;!found_sample && it!=finalSet.end(); ++it)
    {
        const sampling::OctreeReport& bestReport = *it;
        /*ProjectionReport */rep = projectSampleToObstacle(contactGenHelper.fullBody_, limbId, limb, bestReport, validation, configuration, current);
        if(rep.success_)
        {
            double robustness = stability::IsStable(contactGenHelper.fullBody_,rep.result_, contactGenHelper.acceleration_);
            if(    !contactGenHelper.checkStabilityGenerate_
                || (rep.result_.nbContacts == 1 && !contactGenHelper.stableForOneContact_)
                || robustness>=contactGenHelper.robustnessTreshold_)
            {
                maxRob = std::max(robustness, maxRob);
                position = limb->effector_->currentTransformation().getTranslation();
                rotation = limb->effector_->currentTransformation().getRotation();
                normal = rep.result_.contactNormals_.at(limbId);
                found_sample = true;
            }
            // if no stable candidate is found, select best contact
            // anyway
            else if((robustness > maxRob) && contactGenHelper.contactIfFails_)
            {
                moreRobust = configuration;
                maxRob = robustness;
                position = limb->effector_->currentTransformation().getTranslation();
                rotation = limb->effector_->currentTransformation().getRotation();
                normal = rep.result_.contactNormals_.at(limbId);
                unstableContact = true;
            }
        }
    }
    if(found_sample || unstableContact)
    {
        current.contacts_[limbId] = true;
        current.contactNormals_[limbId] = normal;
        current.contactPositions_[limbId] = position;
        current.contactRotation_[limbId] = rotation;
        current.contactOrder_.push(limbId);
    }
    if(found_sample)
    {
        current.configuration_ = configuration;
        current.stable = true;
    }
    if(unstableContact)
    {
        current.configuration_ = moreRobust;
        current.stable = false;
    }
    return current;
}

ProjectionReport generate_contact(const ContactGenHelper &contactGenHelper, const std::string& limbName,
                                  const sampling::heuristic evaluate)
{
    ProjectionReport rep;

    RbPrmLimbPtr_t limb = contactGenHelper.fullBody_->GetLimbs().at(limbName);
    core::CollisionValidationPtr_t validation = contactGenHelper.fullBody_->GetLimbCollisionValidation().at(limbName);
    limb->limb_->robot()->currentConfiguration(contactGenHelper.workingState_.configuration_);
    limb->limb_->robot()->computeForwardKinematics ();

    // pick first sample which is collision free
    bool found_sample(false);
    bool unstableContact(false); //set to true in case no stable contact is found
    rep.result_ = findValidCandidate(contactGenHelper,limbName,limb, validation, found_sample,unstableContact, evaluate);
    if(found_sample)
    {
        rep.status_ = STABLE_CONTACT;
        rep.success_ = true;
#ifdef PROFILE
  RbPrmProfiler& watch = getRbPrmProfiler();
  watch.add_to_count("contact", 1);
#endif
    }
    else if(unstableContact)
    {
        rep.status_ = UNSTABLE_CONTACT;
        rep.success_ = !contactGenHelper.checkStabilityGenerate_;
#ifdef PROFILE
  RbPrmProfiler& watch = getRbPrmProfiler();
  watch.add_to_count("unstable contact", 1);
#endif
    }
    else
    {
        rep =  setCollisionFree(contactGenHelper.fullBody_,validation,limbName,rep.result_);
        rep.status_ = NO_CONTACT;
        rep.success_ = false;
#ifdef PROFILE
  RbPrmProfiler& watch = getRbPrmProfiler();
  watch.add_to_count("no contact", 1);
#endif
    }
    return rep;
}

ProjectionReport gen_contacts(ContactGenHelper &contactGenHelper)
{
    ProjectionReport rep;
    T_ContactState candidates = gen_contacts_combinatorial(contactGenHelper);
    std::vector <std::string> reqLimbs(contactGenHelper.fullBody_->getRequiredLimbs());
    while(!candidates.empty() && !rep.success_)
    {
        //retrieve latest state
        ContactState cState = candidates.front();
        candidates.pop();
        bool checkStability(contactGenHelper.checkStabilityGenerate_);
        contactGenHelper.checkStabilityGenerate_ = false; // stability not mandatory before last contact is created
        if(cState.second.empty() && contactGenHelper.workingState_.stable)
        {
            bool reqLimValid(true);
            // check if all required contacts are active
            for(unsigned int i = 0 ; reqLimValid && (i < reqLimbs.size()) ; ++i)
            {
                if(contactGenHelper.workingState_.contactPositions_.count(reqLimbs[i]) == 0)
                    reqLimValid = false;
            }            

            // if all required contacts are active
            if(reqLimValid)
            {
                rep.result_ = contactGenHelper.workingState_;
                rep.status_ = NO_CONTACT;
                rep.success_ = true;
                return rep;
            }
        }
        for(std::vector<std::string>::const_iterator cit = cState.second.begin();
            cit != cState.second.end(); ++cit)
        {
            if(cit+1 == cState.second.end())
                contactGenHelper.checkStabilityGenerate_ = checkStability;
            rep = generate_contact(contactGenHelper,*cit);
            if(rep.success_)
            {
                contactGenHelper.workingState_ = rep.result_;
            }
            //else
            //    break;
        }
    }
    return rep;
}

projection::ProjectionReport repositionContacts(ContactGenHelper& helper)
{
    ProjectionReport resultReport;
    State result = helper.workingState_;
    result.stable = false;
    State previous = result;
    // replace existing contacts
    // start with older contact created
    std::stack<std::string> poppedContacts;
    std::queue<std::string> oldOrder = result.contactOrder_;
    std::queue<std::string> newOrder;
    std::string nContactName ="";
    core::Configuration_t savedConfig = helper.previousState_.configuration_;
    core::Configuration_t config = savedConfig;
    while(!result.stable &&  !oldOrder.empty())
    {
        std::string previousContactName = oldOrder.front();
        std::string groupName = helper.fullBody_->GetLimbs().at(previousContactName)->limb_->name();
        const std::vector<std::string>& group = helper.fullBody_->GetGroups().at(groupName);
        oldOrder.pop();
        core::ConfigurationIn_t save = helper.fullBody_->device_->currentConfiguration();
        bool notFound(true);
        for(std::vector<std::string>::const_iterator cit = group.begin();
            notFound && cit != group.end(); ++cit)
        {
            result.RemoveContact(*cit);
            helper.workingState_ = result;
            projection::ProjectionReport rep = contact::generate_contact(helper,*cit);
            if(rep.status_ == STABLE_CONTACT)
            {
                nContactName = *cit;
                notFound = false;
                result = rep.result_;
            }
            else
            {
                result = previous;
                config = savedConfig;
            }
        }
        if(notFound)
        {
            config = savedConfig;
            result.configuration_ = savedConfig;
            poppedContacts.push(previousContactName);
            helper.fullBody_->device_->currentConfiguration(save);
        }
    }
    while(!poppedContacts.empty())
    {
        newOrder.push(poppedContacts.top());
        poppedContacts.pop();
    }
    while(!oldOrder.empty())
    {
        newOrder.push(oldOrder.front());
        oldOrder.pop();
    }
    if(result.stable)
    {
        newOrder.push(nContactName);
        resultReport.status_ = STABLE_CONTACT;
        resultReport.success_ = true;
    }
    result.contactOrder_ = newOrder;
    resultReport.result_ = result;
    return resultReport;
}

/* ZMP criterion */

Vec2D & Vec2D::operator=(const Vec2D & c)
{
    if(this != &c)
    {
        this->x = c.x;
        this->y = c.y;
    }
    return *this;
}
double Vec2D::operator[](int idx) const
{
    idx = idx % 2;
    if(idx == 0)
        return this->x;
    else
        return this->y;
}
double & Vec2D::operator[](int idx)
{
    idx = idx % 2;
    if(idx == 0)
        return this->x;
    else
        return this->y;
}
bool operator==(const Vec2D & v1, const Vec2D & v2)
{
    //return ((v1.x == v2.x) && (v1.y == v2.y));
    return ((std::abs(v1.x - v2.x) < 1e-9) && (std::abs(v1.y - v2.y) < 1e-9));
}
bool operator!=(const Vec2D & v1, const Vec2D & v2)
{
    return !(v1 == v2);
}
std::ostream & operator<<(std::ostream & out, const Vec2D & v)
{
    out << "(" << v.x << ", " << v.y << ")";
    return out;
}
Plane & Plane::operator=(const Plane & pe)
{
    if(this != &pe)
    {
        this->a = pe.a;
        this->b = pe.b;
        this->c = pe.c;
        this->d = pe.d;
    }
    return *this;
}
fcl::Vec3f orthogonalProjection(const fcl::Vec3f & point, const Plane & plane)
{
    double k(((plane.a * point[0]) + (plane.b * point[1]) + (plane.c * point[2]) + plane.d) / (std::pow(plane.a, 2) + std::pow(plane.b, 2) + std::pow(plane.c, 2)));
    return fcl::Vec3f(point[0] - k*plane.a, point[1] - k*plane.b, point[2] - k*plane.c);
}
std::vector <Vec2D> compute_support_polygon(const std::map <std::string, fcl::Vec3f> & contactPositions)
{
    Plane h_plane(0, 0, 1, 0); // horizontal plane
    std::vector <Vec2D> res;
    for(std::map<std::string, fcl::Vec3f>::const_iterator cit = contactPositions.begin(); cit != contactPositions.end(); ++cit)
    {
        //fcl::Vec3f proj(orthogonalProjection(cit->second, h_plane));
        //Vec2D vertex_2D(proj[0], proj[1]);
        //res.push_back(vertex_2D);
        res.push_back(Vec2D(cit->second[0], cit->second[1])); // because the plane is horizontal, we just have to remove the z (vertical) component
    }
    return res;
}
double computeAngle(const Vec2D & center, const Vec2D & end1, const Vec2D & end2)
{
    // build vector1
    Vec2D vector1(end1.x - center.x, end1.y - center.y);
    double norm1(std::sqrt(std::pow(vector1.x, 2) + std::pow(vector1.y, 2)));

    // build vector2
    Vec2D vector2(end2.x - center.x, end2.y - center.y);
    double norm2(std::sqrt(std::pow(vector2.x, 2) + std::pow(vector2.y, 2)));

    // find the angle between the two vectors using their scalar product
    double sp((vector1.x*vector2.x) + (vector1.y*vector2.y));
    return std::acos(sp/(norm1*norm2));
}
void scanningProcess(const Vec2D & basePoint, std::vector <Vec2D> & subset, double & angle, const Vec2D & currentPoint, bool higher, bool direction)
{
    // higher == true --> currentPoint is above basePoint
    // higher == false --> currentPoint is below basePoint
    // direction == true --> scan to the right
    // direction == false --> scan to the left
    int higher_val = higher ? 1 : -1;
    int direction_val = direction ? 1 : -1;

    if(subset.size() == 1)
    {
        // init
        angle = computeAngle(basePoint, Vec2D(basePoint.x + direction_val, basePoint.y), currentPoint);
        subset.push_back(currentPoint);
    }
    else if((higher_val*currentPoint.y) >= (higher_val*subset.back().y))
    {
        double opening(computeAngle(subset.back(), Vec2D(subset.back().x + direction_val, subset.back().y), currentPoint));
        if(opening <= angle)
        {
            subset.push_back(currentPoint);
            angle = opening;
        }
        else
        {
            subset.pop_back();
            opening = computeAngle(subset.back(), Vec2D(subset.back().x + direction_val, subset.back().y), currentPoint);
            bool convex(false);
            if(subset.size() == 1)
                convex = true;
            while(!convex)
            {
                Vec2D base(subset[subset.size() - 2]);
                angle = computeAngle(base, Vec2D(base.x + direction_val, base.y), subset.back());
                if(angle < opening)
                {
                    subset.pop_back();
                    opening = computeAngle(subset.back(), Vec2D(subset.back().x + direction_val, subset.back().y), currentPoint);
                }
                else
                    convex = true;
                if(subset.size() == 1)
                    convex = true;
            }
            subset.push_back(currentPoint);
            angle = opening;
        }
    }
}
std::vector <Vec2D> convexHull(std::vector <Vec2D> set)
{
    std::vector <Vec2D> res_tmp, res;
    
    if(!set.empty())
    {
        /* sort the input set by x increasing */
        std::vector <Vec2D> sortedSet;
        while(!set.empty())
        {
            double index(0);
            double min(set[index].x);
            for(unsigned int i = 0; i < set.size(); ++i)
            {
                if(set[i].x < min)
                {
                    min = set[i].x;
                    index = i;
                }
            }
            sortedSet.push_back(set[index]);
            set.erase(set.begin()+index);
        }

        /* first scanning, to the right */
        Vec2D basePoint(sortedSet[0]);
        std::vector <Vec2D> tr_upper_set; tr_upper_set.push_back(basePoint);
        std::vector <Vec2D> tr_lower_set; tr_lower_set.push_back(basePoint);
        double upperAngle, lowerAngle;
        for(unsigned int i = 1; i < sortedSet.size(); ++i)
        {
            if(sortedSet[i].y >= basePoint.y) // if the point is upper than basePoint
            {
                scanningProcess(basePoint, tr_upper_set, upperAngle, sortedSet[i], true, true);
            }
            else // if the point is lower than basePoint
            {
                scanningProcess(basePoint, tr_lower_set, lowerAngle, sortedSet[i], false, true);
            }
        }

        /* second scanning, to the left */
        basePoint = sortedSet.back();
        std::vector <Vec2D> tl_upper_set; tl_upper_set.push_back(basePoint);
        std::vector <Vec2D> tl_lower_set; tl_lower_set.push_back(basePoint);
        for(int i = (int)sortedSet.size() - 2; i >= 0; --i)
        {
            if(sortedSet[i].y >= basePoint.y) // if the point is upper than basePoint
            {
                scanningProcess(basePoint, tl_upper_set, upperAngle, sortedSet[i], true, false);
            }
            else // if the point is lower than basePoint
            {
                scanningProcess(basePoint, tl_lower_set, lowerAngle, sortedSet[i], false, false);
            }
        }

        /* merge the four subsets without keeping the duplicates (subsets boundaries, ...) */
        for(unsigned int i = 0; i < tr_upper_set.size(); ++i)
        {
            res_tmp.push_back(tr_upper_set[i]);
        }
        for(int i = (int)tl_upper_set.size() - 1; i >= 0; --i)
        {
            res_tmp.push_back(tl_upper_set[i]);
        }
        for(unsigned int i = 0; i < tl_lower_set.size(); ++i)
        {
            res_tmp.push_back(tl_lower_set[i]);
        }
        for(int i = (int)tr_lower_set.size() - 1; i >= 0; --i)
        {
            res_tmp.push_back(tr_lower_set[i]);
        }
        for(unsigned int i = 0; i < res_tmp.size(); ++i)
        {
            if(!contains(res, res_tmp[i]))
                res.push_back(res_tmp[i]);
        }
    }

    return res;
}
double orientedAngle2D(const Vec2D & center, const Vec2D & base, const Vec2D & goal)
{
    // get the two vectors
    Vec2D vbase(base.x - center.x, base.y - center.y);
    Vec2D vgoal(goal.x - center.x, goal.y - center.y);

    // normalize the vectors
    double normBase(std::sqrt(std::pow(vbase.x, 2) + std::pow(vbase.y, 2)));
    double normGoal(std::sqrt(std::pow(vgoal.x, 2) + std::pow(vgoal.y, 2)));
    for(int i = 0; i < 2; ++i)
    {
        vbase[i] /= normBase;
        vgoal[i] /= normGoal;
    }

    // update the new norms
    normBase = std::sqrt(std::pow(vbase.x, 2) + std::pow(vbase.y, 2));
    normGoal = std::sqrt(std::pow(vgoal.x, 2) + std::pow(vgoal.y, 2));

    // calculate the angle between the two vectors
    double dotProduct((vbase.x * vgoal.x) + (vbase.y * vgoal.y));
    double angle(std::acos(dotProduct/(normBase * normGoal)));

    // calculate the orientation of the angle
    if( ((vbase.x >= 0) && (vgoal.x > 0)) || ((vbase.x > 0) && (vgoal.x >= 0)) ) // both vectors in the right side
    {
        if(vbase.y > vgoal.y)
            angle = -angle;
    }
    else if( ((vbase.x <= 0) && (vgoal.x < 0)) || ((vbase.x < 0) && (vgoal.x <= 0)) ) // both vectors in the left side
    {
        if(vbase.y < vgoal.y)
            angle = -angle;
    }
    else if( ((vbase.y >= 0) && (vgoal.y > 0)) || ((vbase.y > 0) && (vgoal.y >= 0)) ) // both vectors in the top side
    {
        if(vbase.x < vgoal.x)
            angle = -angle;
    }
    else if( ((vbase.y <= 0) && (vgoal.y < 0)) || ((vbase.y < 0) && (vgoal.y <= 0)) ) // both vectors in the bottom side
    {
        if(vbase.x > vgoal.x)
            angle = -angle;
    }
    else if( ((vbase.x == 0) && (vgoal.x == 0)) || ((vbase.y == 0) && (vgoal.y == 0)) ) // both vectors along the same axis
    {
        // Unknown direction, we cannot determine whether the angle goes from base to goal by the first or the second side of the straight line
    }
    else // vectors in opposite quadrants
    {
        if( ((vbase.x > 0) && (vbase.y > 0)) || ((vgoal.x < 0) && (vgoal.y < 0)) ) // first and third quadrants (posX-posY, negX-negY)
        {
            vbase.x = std::abs(vbase.x); vbase.y = std::abs(vbase.y);
			vgoal.x = std::abs(vgoal.x); vgoal.y = std::abs(vgoal.y);
            if(vbase.y < vgoal.y)
				angle = -angle;
        }
        else // second and fourth quadrants (negX-posY, posX-negY)
        {
            vbase.x = std::abs(vbase.x); vbase.y = std::abs(vbase.y);
			vgoal.x = std::abs(vgoal.x); vgoal.y = std::abs(vgoal.y);
            if(vbase.y > vgoal.y)
				angle = -angle;
        }
    }
    // return the result
    return angle;
}
bool isInside(const Vec2D & point, const std::vector <Vec2D> & polygon)
{
    // get the winding number (sumAngles == windingNumber*2*pi)
    double sumAngles(orientedAngle2D(point, polygon.back(), polygon.front())); // sumAngles is the sum of all subtended angles by each polygon edge from the considered point
    Vec2D base, goal;
    for(unsigned int i = 0; i < polygon.size() - 1; ++i)
    {
        base = polygon[i];
        goal = polygon[i+1];
        sumAngles += orientedAngle2D(point, base, goal);
    }

    double pi(3.14159265);
    if(std::abs(sumAngles) < 1e-5) // if sumAngles == 0 --> point is outside polygon
    {
        return false;
    }
    else if(std::abs(sumAngles - 2*pi) < 1e-5) // if sumAngles == 2*pi --> point is inside polygon
    {
        return true;
    }
    else // Impossible case
    {
        throw std::string("The polygon contains loops, please ensure that the polygon vertices are in order");
    }
}
double euclideanDist(const Vec2D & v1, const Vec2D & v2)
{
    return std::sqrt(std::pow(v1.x - v2.x, 2) + std::pow(v1.y - v2.y, 2));
}
Vec2D weightedCentroidConvex2D(const std::vector <Vec2D> & convexPolygon)
{
    if(convexPolygon.empty())
        throw std::string("Impossible to find the weighted centroid of nothing (the specified convex polygon has no vertices)");

    Vec2D res;
    if(convexPolygon.size() == 1)
        res = convexPolygon[0];
    else if(convexPolygon.size() == 2)
    {
        double resX((convexPolygon[0].x + convexPolygon[1].x) / 2.0);
        double resY((convexPolygon[0].y + convexPolygon[1].y) / 2.0);
        res = Vec2D(resX, resY);
    }
    else
    {
        // get the longest edge and define the minimum admissible threshold for counting a vertex as a single point
        double maxDist(euclideanDist(convexPolygon.back(), convexPolygon.front()));
        for(unsigned int i = 0; i < convexPolygon.size() - 1; ++i)
        {
            double dist(euclideanDist(convexPolygon[i], convexPolygon[i+1]));
            if(dist > maxDist)
                maxDist = dist;
        }
        double threshold(maxDist/10.0);

        // shift the list until starting from a lonely (to the rear) point
        std::vector <Vec2D> shifted(convexPolygon);
        while(euclideanDist(shifted.back(), shifted.front()) <= threshold)
        {
            shifted.push_back(shifted.front());
            shifted.erase(shifted.begin());
        }

        // look over the shifted set
        std::vector <Vec2D> finalSet, localSubset;
        bool subsetOngoing = false;
        shifted.push_back(shifted.front());

        for(unsigned int i = 0; i < shifted.size() - 1; ++i)
        {
            if(euclideanDist(shifted[i], shifted[i+1]) > threshold)
            {
                if(!subsetOngoing)
                    finalSet.push_back(shifted[i]);
                else
                {
                    localSubset.push_back(shifted[i]);
                    double moyX(0.0), moyY(0.0);
                    for(unsigned int j = 0; j < localSubset.size(); ++j)
                    {
                        moyX += localSubset[j].x;
                        moyY += localSubset[j].y;
                    }
                    moyX /= localSubset.size();
                    moyY /= localSubset.size();
                    finalSet.push_back(Vec2D(moyX, moyY));
                    localSubset.clear();
                    subsetOngoing = false;
                }
            }
            else
            {
                localSubset.push_back(shifted[i]);
                if(!subsetOngoing)
                    subsetOngoing = true;
            }
        }

        double resX(0.0), resY(0.0);
        for(unsigned int i = 0; i < finalSet.size(); ++i)
        {
            resX += finalSet[i].x;
            resY += finalSet[i].y;
        }
        resX /= finalSet.size();
        resY /= finalSet.size();
        res = Vec2D(resX, resY);
    }
    return res;
}
Vec2D computeZMP(const fcl::Vec3f & comPos, const fcl::Vec3f & comAccel, double g)
{

}
bool isValidZMP(const hpp::rbprm::State & state, const fcl::Vec3f & comPos, const fcl::Vec3f & comAccel, double g)
{

}
double evaluateZMP(const hpp::rbprm::State & state, const fcl::Vec3f & comPos, const fcl::Vec3f & comAccel, double g)
{

}

} // namespace projection
} // namespace rbprm
} // namespace hpp
