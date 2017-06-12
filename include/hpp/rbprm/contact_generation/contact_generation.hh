/// Copyright (c) 2017 CNRS
/// Authors: stonneau
///
///
// This file is part of hpp-rbprm.
// hpp-rbprm is free software: you can redistribute it
// and/or modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation, either version
// 3 of the License, or (at your option) any later version.
//
// hpp-wholebody-step-planner is distributed in the hope that it will be
// useful, but WITHOUT ANY WARRANTY; without even the implied warranty
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Lesser Public License for more details. You should have
// received a copy of the GNU Lesser General Public License along with
// hpp-wholebody-step-planner. If not, see
// <http://www.gnu.org/licenses/>.

#ifndef HPP_RBPRM_CONTACT_GENERATION_HH
# define HPP_RBPRM_CONTACT_GENERATION_HH

# include <hpp/rbprm/rbprm-state.hh>
# include <hpp/rbprm/rbprm-fullbody.hh>
# include <hpp/rbprm/projection/projection.hh>

# include <queue>

namespace hpp {
namespace rbprm {
namespace contact{

typedef std::queue<hpp::rbprm::State> Q_State;
typedef std::pair <hpp::rbprm::State, std::vector<std::string> > ContactState;
typedef std::queue<ContactState> T_ContactState;

struct ContactGenHelper
{
     ContactGenHelper(RbPrmFullBodyPtr_t fb, const State& ps,
                      model::ConfigurationIn_t configuration,
                      const affMap_t& affordances,
                      const std::map<std::string, std::vector<std::string> >& affFilters,
                      const double robustnessTreshold = 0,
                      const std::size_t maxContactBreaks = 1,
                      const std::size_t maxContactCreations = 1,
                      const bool checkStabilityMaintain = false,
                      const bool checkStabilityGenerate = true,
                      const fcl::Vec3f& direction = fcl::Vec3f(0,0,1),
                      const fcl::Vec3f& acceleration = fcl::Vec3f(0,0,0),
                      const bool contactIfFails = false,
                      const bool stableForOneContact = false);
    ~ContactGenHelper(){}
    hpp::rbprm::RbPrmFullBodyPtr_t fullBody_;
    const hpp::rbprm::State previousState_;
    const bool checkStabilityMaintain_;
    bool contactIfFails_;
    const bool stableForOneContact_;
    const fcl::Vec3f acceleration_;
    const fcl::Vec3f direction_;
    const double robustnessTreshold_;
    const std::size_t maxContactBreaks_;
    const std::size_t maxContactCreations_;
    const affMap_t& affordances_;
    const std::map<std::string, std::vector<std::string> >& affFilters_;
    hpp::rbprm::State workingState_;
    bool checkStabilityGenerate_;
    Q_State candidates_;
};


hpp::model::ObjectVector_t HPP_RBPRM_DLLAPI getAffObjectsForLimb(const std::string& limb,
    const affMap_t& affordances, const std::map<std::string, std::vector<std::string> >& affFilters);

/// Generates all potentially valid cases of valid contact maintenance
/// given a previous configuration.
/// \param fullBody target Robot
/// \param target, desired root position
/// \param currentState current state of the robot (configuration and contacts)
/// \param maxBrokenContacts max number of contacts that can be broken in the process
/// \return a queue of contact states candidates for maintenance, ordered by number of contacts broken
/// and priority in the list wrt the contact order
Q_State maintain_contacts_combinatorial(const hpp::rbprm::State& currentState, const std::size_t maxBrokenContacts=1);

/// Generates all potentially valid cases of valid contact creation by removing the top state of the priority
/// stack
/// \param freeEffectors list of free candidates
/// \param previous current state of the robot
/// \param maxCreatedContacts max number of contacts that can be created in the process
/// \return a QUEUE of contact states candidates for maintenance, ordered by number of contacts broken
/// and priority in the list wrt the contact order
T_ContactState gen_contacts_combinatorial(const std::vector<std::string>& freeEffectors, const State& previous, const std::size_t maxCreatedContacts);

/// Given a combinatorial of possible contacts, generate
/// the first "valid" configuration, that is the first kinematic
/// configuration that removes the minimum number of contacts and
/// is collision free.
/// \param ContactGenHelper parametrization of the planner
/// \return the best candidate wrt the priority in the list and the contact order
projection::ProjectionReport maintain_contacts(ContactGenHelper& contactGenHelper);

/// Given a current state and an effector, tries to generate a contact configuration.
/// \param ContactGenHelper parametrization of the planner
/// \param limb the limb to create a contact with
/// \return the best candidate wrt the priority in the list and the contact order
projection::ProjectionReport generate_contact(const ContactGenHelper& contactGenHelper, const std::string& limb,
                                              const sampling::heuristic evaluate = 0);

/// Given a combinatorial of possible contacts, generate
/// the first "valid" contact configuration, that is the first contact
/// configuration is contact and equilibrium.
/// \param ContactGenHelper parametrization of the planner
/// \return the best candidate wrt the priority in the list and the contact order
projection::ProjectionReport gen_contacts(ContactGenHelper& contactGenHelper);

/// Tries to reposition one contact of a given state into
/// a new one, more balanced
/// \param ContactGenHelper parametrization of the planner
/// \return the best candidate wrt the priority in the list and the contact order
projection::ProjectionReport repositionContacts(ContactGenHelper& contactGenHelper);

/* ZMP criterion */

/// Data structure to store 2-dimensional informations (2D vectors)
struct Vec2D
{
    double x;
    double y;
	Vec2D() : x(0), y(0) {}
	Vec2D(double xx, double yy) : x(xx), y(yy) {}
	Vec2D(const Vec2D & c2D) : x(c2D.x), y(c2D.y) {}
    Vec2D & operator=(const Vec2D & c);
    double operator[](int idx) const;
    double & operator[](int idx);
};
bool operator==(const Vec2D & v1, const Vec2D & v2);
bool operator!=(const Vec2D & v1, const Vec2D & v2);
std::ostream & operator<<(std::ostream & out, const Vec2D & v);

/// Function to verify the existence of an element in a std::vector
template <typename T>
bool contains(const std::vector <T> & vect, const T & val)
{
    bool found(false);
    for(unsigned int i = 0; !found && (i < vect.size()); ++i)
    {
        if(vect[i] == val)
            found = true;
    }
    return found;
}

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

/// Compute the position of the Zero-Moment-Point (ZMP)
///
/// \param comPos The position of the center of mass of the robot
/// \param comAccel The center of mass acceleration
/// \param g The gravity acceleration, default value: -9.80665 m.s-1
/// \return The position of the ZMP
Vec2D computeZMP(const fcl::Vec3f & comPos, const fcl::Vec3f & comAccel, double g = -9.80665);

/// Perform the validation of a contact state according to the ZMP criterion
/// The ZMP must be inside the convex hull of the support polygon (relevant only for coplanar contacts on the ground horizontal plane)
/// The support polygon is formed by the orthogonal projection of each contact in the ground horizontal plane.
///
/// \param state The considered state of the robot
/// \param comPos The position of the center of mass of the robot
/// \param comAccel The center of mass acceleration
/// \param g The gravity acceleration, default value: -9.80665 m.s-1
/// \return true if the contact state is valid, false otherwise
bool isValidZMP(const hpp::rbprm::State & state, const fcl::Vec3f & comPos, const fcl::Vec3f & comAccel, double g = -9.80665);

/// Evaluate the cost of a contact state according to the ZMP criterion
///
/// \param state The considered state of the robot
/// \param comPos The position of the center of mass of the robot
/// \param comAccel The center of mass acceleration
/// \param g The gravity acceleration, default value: -9.80665 m.s-1
/// \return The cost of the considered contact state (the lower the cost is, the better the contact state is --> cost to minimize)
double evaluateZMP(const hpp::rbprm::State & state, const fcl::Vec3f & comPos, const fcl::Vec3f & comAccel, double g = -9.80665);


    } // namespace projection
  } // namespace rbprm
} // namespace hpp
#endif // HPP_RBPRM_CONTACT_GENERATION_HH
