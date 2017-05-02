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

#include <hpp/rbprm/stability/stability.hh>
#include <hpp/rbprm/stability/support.hh>
#include <hpp/model/device.hh>
#include <hpp/model/joint.hh>
#include <hpp/model/center-of-mass-computation.hh>
#include <hpp/rbprm/tools.hh>

#include <robust-equilibrium-lib/static_equilibrium.hh>

#include <Eigen/Dense>

#include <vector>
#include <map>
#include <string>

#ifdef PROFILE
    #include "hpp/rbprm/rbprm-profiler.hh"
#endif

using namespace hpp;
using namespace hpp::core;
using namespace hpp::model;
using namespace hpp::rbprm;
using namespace robust_equilibrium;

namespace hpp {
namespace rbprm {
namespace stability{

    void computeRectangleContact(const std::string& name, const RbPrmLimbPtr_t limb, const State& state, Ref_matrix43 p)
    {
        const double& lx = limb->x_, ly = limb->y_;
        const fcl::Vec3f& position = state.contactPositions_.at(name);
        //create rotation matrix from normal
        Eigen::Matrix3d R;
        p << lx,  ly, 0,
             lx, -ly, 0,
            -lx, -ly, 0,
            -lx,  ly, 0;
        if(limb->contactType_ == _3_DOF)
        {
            //create rotation matrix from normal
            const fcl::Vec3f& normal = state.contactNormals_.at(name);
            const fcl::Vec3f z_current = limb->effector_->currentTransformation().getRotation() * limb->normal_;
            const fcl::Matrix3f alignRotation = tools::GetRotationMatrix(z_current,normal);
            const fcl::Matrix3f rotation = alignRotation * limb->effector_->currentTransformation().getRotation();
            const fcl::Vec3f offset = rotation * limb->offset_;
            Eigen::Vector3d z,x,y;
            for(int i =0; i<3; ++i) z[i] = normal[i];
            x = z.cross(Eigen::Vector3d(0,-1,0));
            if(x.norm() < 10e-6)
            {
                y = z.cross(fcl::Vec3f(1,0,0));
                y.normalize();
                x = y.cross(z);
            }
            else
            {
                x.normalize();
                y = z.cross(x);
            }
            R.block<3,1>(0,0) = x;
            R.block<3,1>(0,1) = y;
            R.block<3,1>(0,2) = z;

            for(std::size_t i =0; i<4; ++i)
            {
                p.row(i) = position + (R*(p.row(i).transpose())) + offset;
            }
        }
        else
        {
            fcl::Vec3f z_axis(0,0,1);
            fcl::Matrix3f rotationLocal = tools::GetRotationMatrix(z_axis, limb->normal_);
            rotationLocal.inverse();
            fcl::Transform3f roWorld;
            roWorld.setRotation(state.contactRotation_.at(name));
            roWorld.setTranslation(position);
            for(std::size_t i =0; i<4; ++i)
            {
                fcl::Vec3f pLocal = rotationLocal*(p.row(i).transpose()) + limb->offset_;
                p.row(i) = (roWorld * pLocal).getTranslation();
            }
        }
    }

    void computePointContact(const std::string& name, const RbPrmLimbPtr_t limb, const State& state, Ref_vector3 p)
    {
        const fcl::Vec3f& position = state.contactPositions_.at(name);
        //create rotation matrix from normal
        const fcl::Vec3f& normal = state.contactNormals_.at(name);
        const fcl::Vec3f z_current = limb->effector_->currentTransformation().getRotation() * limb->normal_;
        const fcl::Matrix3f alignRotation = tools::GetRotationMatrix(z_current,normal);
        const fcl::Matrix3f rotation = alignRotation * limb->effector_->currentTransformation().getRotation();
        const fcl::Vec3f offset = rotation * limb->offset_;
        p = position + offset;
    }

    StaticEquilibrium initLibrary(const RbPrmFullBodyPtr_t fullbody)
    {
        return StaticEquilibrium(fullbody->device_->name(), fullbody->device_->mass(),4,SOLVER_LP_QPOASES,true,10,false);
    }

    const std::size_t numContactPoints(const RbPrmLimbPtr_t& limb)
    {
        if(limb->contactType_ == _3_DOF)
            return 1;
        else
            return 4;
    }

    const std::vector<std::size_t> numContactPoints(const T_Limb& limbs, const std::vector<std::string>& contacts, std::size_t& totalNumContacts)
    {
        std::size_t n;
        std::vector<std::size_t> res;
        for(std::vector<std::string>::const_iterator cit = contacts.begin();
            cit != contacts.end(); ++cit)
        {
            n = numContactPoints(limbs.at(*cit));
            totalNumContacts+= n;
            res.push_back(n);
        }
        return res;
    }

    robust_equilibrium::Vector3 setupLibrary(const RbPrmFullBodyPtr_t fullbody, State& state, StaticEquilibrium& sEq, StaticEquilibriumAlgorithm alg,
                                             const core::value_type friction = 0.3)
    {
        hpp::model::ConfigurationIn_t save = fullbody->device_->currentConfiguration();
        std::vector<std::string> contacts;
        for(std::map<std::string,bool>::const_iterator cit = state.contacts_.begin();
            cit!=state.contacts_.end(); ++ cit)
        {
            if(cit->second) contacts.push_back(cit->first);
        }
        fullbody->device_->currentConfiguration(state.configuration_);
        fullbody->device_->computeForwardKinematics();
        const T_Limb limbs = fullbody->GetLimbs();
        std::size_t nbContactPoints(0);
        std::vector<std::size_t> contactPointsInc = numContactPoints(limbs, contacts,nbContactPoints);
        robust_equilibrium::MatrixX3 normals  (nbContactPoints,3);
        robust_equilibrium::MatrixX3 positions(nbContactPoints,3);
        std::size_t currentIndex(0), c(0);
        for(std::vector<std::size_t>::const_iterator cit = contactPointsInc.begin();
            cit != contactPointsInc.end(); ++cit, ++c)
        {
            const RbPrmLimbPtr_t limb =limbs.at(contacts[c]);
            const fcl::Vec3f& n = state.contactNormals_.at(contacts[c]);
            Vector3 normal(n[0],n[1],n[2]);
            const std::size_t& inc = *cit;
            if(inc > 1)
                computeRectangleContact(contacts[c], limb,state,positions.middleRows<4>(currentIndex));
            else
                computePointContact(contacts[c], limb,state,positions.middleRows<1>(currentIndex,inc));
            for(int i =0; i < inc; ++i)
            {
                normals.middleRows<1>(currentIndex+i) = normal;
            }
            currentIndex += inc;
        }
        robust_equilibrium::Vector3 com;
/*model::CenterOfMassComputationPtr_t comcptr = model::CenterOfMassComputation::create(fullbody->device_);
comcptr->add(fullbody->device_->getJointByName("romeo/base_joint_xyz"));
comcptr->computeMass();
comcptr->compute();
const fcl::Vec3f comfcl = comcptr->com();*/
        const fcl::Vec3f comfcl = fullbody->device_->positionCenterOfMass();
        state.com_ = comfcl;
        for(int i=0; i< 3; ++i) com(i)=comfcl[i];
        fullbody->device_->currentConfiguration(save);
        sEq.setNewContacts(positions,normals,friction,alg);
        return com;
    }

    std::pair<MatrixXX, VectorX> ComputeCentroidalCone(const RbPrmFullBodyPtr_t fullbody, State& state, const hpp::core::value_type friction)
    {
        std::pair<MatrixXX, VectorX> res;
        MatrixXX& H = res.first;
        VectorX& h = res.second;
#ifdef PROFILE
        RbPrmProfiler& watch = getRbPrmProfiler();
        watch.start("test balance");
#endif
        StaticEquilibrium staticEquilibrium(initLibrary(fullbody));
        setupLibrary(fullbody,state,staticEquilibrium,STATIC_EQUILIBRIUM_ALGORITHM_PP, friction);
#ifdef PROFILE
    watch.stop("test balance");
#endif
        LP_status status = LP_STATUS_OPTIMAL;
        if(status != LP_STATUS_OPTIMAL)
        {
            std::cout << "error " << std::endl;
        }
        else
        {
            status = staticEquilibrium.getPolytopeInequalities(H,h);
            if(status != LP_STATUS_OPTIMAL)
            {
                std::cout << "error " << std::endl;
                H = Eigen::MatrixXd::Zero(6,6);
                h = Eigen::MatrixXd::Zero(6,1);
            }
        }
        return res;
    }


    double IsStable(const RbPrmFullBodyPtr_t fullbody, State& state, const robust_equilibrium::StaticEquilibriumAlgorithm algorithm)
    {
#ifdef PROFILE
    RbPrmProfiler& watch = getRbPrmProfiler();
    watch.start("test balance");
#endif
        StaticEquilibrium staticEquilibrium(initLibrary(fullbody));
        robust_equilibrium::Vector3 com = setupLibrary(fullbody,state,staticEquilibrium,algorithm);
        double res;LP_status status;
        if(algorithm == STATIC_EQUILIBRIUM_ALGORITHM_PP)
        {
            bool isStable(false);
            status = staticEquilibrium.checkRobustEquilibrium(com,isStable);
            res = isStable? 1. : -1.;
        }
        else // STATIC_EQUILIBRIUM_ALGORITHM_DLP
        {
           status = staticEquilibrium.computeEquilibriumRobustness(com,res);
        }
#ifdef PROFILE
    watch.stop("test balance");
#endif
        if(status != LP_STATUS_OPTIMAL)
        {
            if(status == LP_STATUS_INFEASIBLE || status == LP_STATUS_UNBOUNDED)
                //return 1.1; // completely arbitrary: TODO
                return -1.1; // completely arbitrary: TODO
            return -std::numeric_limits<double>::max();
        }
        return res ;
    }

    /* --- GCR criterion --- */
    // We suppose that x-axis = roll axis (forward oriented), y-axis = pitch axis (sideways oriented), z-axis = yaw axis (upwards oriented).

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
    Vec2D & Vec2D::operator=(const Vec2D & c)
    {
        if(this != &c)
        {
            this->x = c.x;
            this->y = c.y;
        }
        return *this;
    }
    Plane compute_ground_horizontal_plane(const std::map<std::string, fcl::Vec3f> & contactsPositions)
    {
        // find the lowest contact through which the ground plane will pass
        std::map<std::string, fcl::Vec3f>::const_iterator cit = contactsPositions.begin();
        fcl::Vec3f lowest(cit->second);
        for(++cit; cit != contactsPositions.end(); ++cit)
        {
            if(cit->second[2] < lowest[2])
                lowest = cit->second;
        }
        Plane res;
        res.a = 0; res.b = 0; res.c = 1; // definition of the normal vector of the ground plane
        res.d = -res.c*lowest[2];
        return res;
    }
    fcl::Vec3f orthogonal_projection(const fcl::Vec3f & point, const Plane & groundPlane)
    {
        double k(((groundPlane.a * point[0]) + (groundPlane.b * point[1]) + (groundPlane.c * point[2]) + groundPlane.d) /
                (std::pow(groundPlane.a, 2) + std::pow(groundPlane.b, 2) + std::pow(groundPlane.c, 2)));
        return fcl::Vec3f(point[0] - k*groundPlane.a, point[1] - k*groundPlane.b, point[2] - k*groundPlane.c);
    }
    double distance_to_plane(const fcl::Vec3f & contact, const Plane & groundPlane)
    {
        fcl::Vec3f proj(orthogonal_projection(contact, groundPlane));
        return std::sqrt(std::pow(contact[0] - proj[0], 2) + std::pow(contact[1] - proj[1], 2) + std::pow(contact[2] - proj[2], 2));
    }
    std::vector <fcl::Vec3f> compute_ground_contacts(const std::map<std::string, fcl::Vec3f> & contactsPositions, const Plane & groundPlane, double threshold)
    {
        std::vector <fcl::Vec3f> res;
        for(std::map<std::string, fcl::Vec3f>::const_iterator cit = contactsPositions.begin(); cit != contactsPositions.end(); ++cit)
        {
            if(distance_to_plane(cit->second, groundPlane) <= threshold)
                res.push_back(cit->second);
        }
        return res;
    }
    std::vector <Vec2D> compute_support_polygon(const std::vector <fcl::Vec3f> & groundContacts, const Plane & groundPlane)
    {
        std::vector <Vec2D> res;
        for(unsigned int i = 0; i < groundContacts.size(); ++i)
        {
            fcl::Vec3f proj(orthogonal_projection(groundContacts[i], groundPlane));
            Vec2D vertex_2D(proj[0], proj[1]);
            res.push_back(vertex_2D);
        }
        return res;
    }
    bool step1(const Vec2D & supportPolygonCenter, const Vec2D & comPos, double threshold)
    {
        return (std::sqrt(std::pow(supportPolygonCenter.x - comPos.x, 2) + std::pow(supportPolygonCenter.y - comPos.y, 2)) <= threshold);
    }
    double step1Eval(const Vec2D & supportPolygonCenter, const Vec2D & comPos)
    {
        return std::sqrt(std::pow(supportPolygonCenter.x - comPos.x, 2) + std::pow(supportPolygonCenter.y - comPos.y, 2));
    }
    double epsilon(double edge, std::size_t nbLegs)
    {
        double res(0);
        if(nbLegs > 1)
        {
            double pi(3.14159265);
            res = edge/(2*std::sin(pi/nbLegs));
        }
        return res;
    }
    bool step2(const std::vector <Vec2D> & supportPolygon, const Vec2D & supportPolygonCenter, double edge, double threshold)
    {
        bool is_valid(true);
        for(unsigned int i = 0; is_valid && (i < supportPolygon.size()); ++i)
        {
            double dist(std::sqrt(std::pow(supportPolygon[i].x - supportPolygonCenter.x, 2) + std::pow(supportPolygon[i].y - supportPolygonCenter.y, 2)));
            if(std::abs(dist - epsilon(edge, supportPolygon.size())) > threshold)
                is_valid = false;
        }
        return is_valid;
    }
    double step2Eval(const std::vector <Vec2D> & supportPolygon, const Vec2D & supportPolygonCenter, double edge)
    {
        double dist(std::sqrt(std::pow(supportPolygon[0].x - supportPolygonCenter.x, 2) + std::pow(supportPolygon[0].y - supportPolygonCenter.y, 2)));
        double epsi(epsilon(edge, supportPolygon.size()));
        double cost(std::abs(dist - epsi));
        // find the worst contact cost
        for(unsigned int i = 1; i < supportPolygon.size(); ++i)
        {
            dist = std::sqrt(std::pow(supportPolygon[i].x - supportPolygonCenter.x, 2) + std::pow(supportPolygon[i].y - supportPolygonCenter.y, 2));
            double eval(std::abs(dist - epsi));
            if(eval > cost)
                cost = eval;
        }
        return cost;
    }
    bool GCRValidation(const hpp::rbprm::State & state, const fcl::Vec3f & comPosition, const GCRParam & config)
    {
        Plane ground(compute_ground_horizontal_plane(state.contactPositions_));
        std::vector <fcl::Vec3f> groundContacts(compute_ground_contacts(state.contactPositions_, ground, config.groundThreshold_));
        std::vector <Vec2D> supportPolygon(compute_support_polygon(groundContacts, ground));
        
        Vec2D polyCenter;
        for(unsigned int i = 0; i < supportPolygon.size(); ++i)
        {
            polyCenter.x += supportPolygon[i].x;
            polyCenter.y += supportPolygon[i].y;
        }
        polyCenter.x /= supportPolygon.size();
        polyCenter.y /= supportPolygon.size();

        fcl::Vec3f comPosProj(orthogonal_projection(comPosition, ground));
        Vec2D comPos2D(comPosProj[0], comPosProj[1]);

        bool step1Valid(step1(polyCenter, comPos2D, config.step1Threshold_));
        bool step2Valid(step2(supportPolygon, polyCenter, config.edge_, config.step2Threshold_));

        return (step1Valid && step2Valid);
    }
    double GCREvaluation(const hpp::rbprm::State & state, const fcl::Vec3f & comPosition, const GCRParam & config)
    {
        Plane ground(compute_ground_horizontal_plane(state.contactPositions_));
        std::vector <fcl::Vec3f> groundContacts(compute_ground_contacts(state.contactPositions_, ground, config.groundThreshold_));
        std::vector <Vec2D> supportPolygon(compute_support_polygon(groundContacts, ground));
        
        Vec2D polyCenter;
        for(unsigned int i = 0; i < supportPolygon.size(); ++i)
        {
            polyCenter.x += supportPolygon[i].x;
            polyCenter.y += supportPolygon[i].y;
        }
        polyCenter.x /= supportPolygon.size();
        polyCenter.y /= supportPolygon.size();

        fcl::Vec3f comPosProj(orthogonal_projection(comPosition, ground));
        Vec2D comPos2D(comPosProj[0], comPosProj[1]);

        double cost1(step1Eval(polyCenter, comPos2D));
        double cost2(step2Eval(supportPolygon, polyCenter, config.edge_));

        return (cost1 + cost2);
    }
    /* --- End GCR criterion --- */
}
}
}
