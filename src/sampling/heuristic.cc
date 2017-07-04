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

#include <hpp/rbprm/sampling/heuristic.hh>
#include <time.h>

#include <Eigen/Eigen>

using namespace hpp;
using namespace hpp::model;
using namespace hpp::rbprm;
using namespace hpp::rbprm::sampling;

double ZMPHeuristic(const sampling::Sample & sample, const Eigen::Vector3d & /*direction*/, const Eigen::Vector3d & /*normal*/, const HeuristicParam & params)
{
    fcl::Vec3f effectorPosition = transform(sample.effectorPosition_, params.tfWorldRoot_.getTranslation(), params.tfWorldRoot_.getRotation());

    std::map <std::string, fcl::Vec3f> contacts;
    contacts.insert(params.contactPositions_.begin(), params.contactPositions_.end());
    contacts.insert(std::make_pair(params.sampleLimbName_, effectorPosition));

    Vec2D zmp;
    double g(params.g_);
    double result;

    if(params.lightVersion_)
    {
        double x_zmp(params.comPosition_[0] - (params.comPosition_[2]/g)*params.comAcceleration_[0]);
        double y_zmp(params.comPosition_[1] - (params.comPosition_[2]/g)*params.comAcceleration_[1]);
        zmp = Vec2D(x_zmp, y_zmp);
    }
    else
    {
        double zAccel(g + params.comAcceleration_[2]);
        double epsi(1e-9);
        // if the z-forces are in balance
        if(std::abs(zAccel) <= epsi) // zAccel == 0
        {
            if((std::abs(params.comAcceleration_[0]) > epsi) || (std::abs(params.comAcceleration_[1]) > epsi)) // (params.comAcceleration_[0] != 0) || (params.comAcceleration_[1] != 0)
                result = std::numeric_limits<double>::max();
            else
                result = Vec2D::euclideanDist(Vec2D(params.comPosition_[0], params.comPosition_[1]), weightedCentroidConvex2D(convexHull(computeSupportPolygon(contacts))));
            return -result; // '-' because minimize a value is equivalent to maximimze its opposite
        }
        double x_zmp(params.comPosition_[0] - (params.comPosition_[2]/zAccel)*params.comAcceleration_[0]);
        double y_zmp(params.comPosition_[1] - (params.comPosition_[2]/zAccel)*params.comAcceleration_[1]);
        zmp = Vec2D(x_zmp, y_zmp);
    }
    try
    {
        Vec2D wcentroid(weightedCentroidConvex2D(convexHull(computeSupportPolygon(contacts))));
        result = Vec2D::euclideanDist(zmp, wcentroid);
    }
    catch(std::string s)
    {
        std::cout << s << std::endl;
        result = std::numeric_limits<double>::max();
    }
    return -result; // '-' because minimize a value is equivalent to maximimze its opposite
}
double StraightCentroidHeuristic(const sampling::Sample & sample, const Eigen::Vector3d & direction, const Eigen::Vector3d & /*normal*/, const HeuristicParam & params)
{
    // Try to keep the centroid vector (formed by the previous centroid and the current centroid) aligned with the direction vector of movement
    double res;
    if(!params.previousContactPositions_.empty())
    {
        Vec2D base_p(weightedCentroidConvex2D(convexHull(computeSupportPolygon(params.previousContactPositions_))));
        fcl::Vec3f effectorPosition = transform(sample.effectorPosition_, params.tfWorldRoot_.getTranslation(), params.tfWorldRoot_.getRotation());
        std::map <std::string, fcl::Vec3f> contacts;
        contacts.insert(params.contactPositions_.begin(), params.contactPositions_.end());
        contacts.insert(std::make_pair(params.sampleLimbName_, effectorPosition));
        Vec2D end1_p(weightedCentroidConvex2D(convexHull(computeSupportPolygon(contacts))));
        Vec2D end2_p(base_p.x + direction[0], base_p.y + direction[1]);
        double coeff(1.0/(1e-9 + computeAngle(base_p, end1_p, end2_p))); // The more we move along the direction, the more the coeff is important (because the ideal angle is 0)
        res = coeff * Vec2D::euclideanDist(base_p, end1_p); // Try to increase the path length as far as possible in the direction of movement
    }
    else
    {
        fcl::Vec3f effectorPosition = transform(sample.effectorPosition_, params.tfWorldRoot_.getTranslation(), params.tfWorldRoot_.getRotation());
        std::map <std::string, fcl::Vec3f> contacts;
        contacts.insert(params.contactPositions_.begin(), params.contactPositions_.end());
        Vec2D base_p(weightedCentroidConvex2D(convexHull(computeSupportPolygon(contacts))));
        contacts.insert(std::make_pair(params.sampleLimbName_, effectorPosition));
        Vec2D wcen(weightedCentroidConvex2D(convexHull(computeSupportPolygon(contacts))));
        Vec2D end_p(base_p.x + direction[0], base_p.y + direction[1]);

        Line2D * line = StraightLine2DManager::straightLineFromPoints2D(base_p, end_p);
        double dist(StraightLine2DManager::distanceToStraightLine2D(wcen, line));
        delete line;

        double coeff(1.0/(1e-9 + dist)); // The more we are close of the direction line (based in 0), the more the coeff is important
        double resX( (end_p.x >= 0) ? wcen.x : -wcen.x );
        double resY( (end_p.y >= 0) ? wcen.y : -wcen.y );
        res = coeff * (resX + resY);
    }
    return res;
}
double CombinedCentroidZMPHeuristic(const sampling::Sample & sample, const Eigen::Vector3d & direction, const Eigen::Vector3d & normal, const HeuristicParam & params)
{
    double res;
    if((std::abs(params.comAcceleration_[0]) <= 1e-3) && (std::abs(params.comAcceleration_[1]) <= 1e-3))
        res = StraightCentroidHeuristic(sample, direction, normal, params);
    else
        res = ZMPHeuristic(sample, direction, normal, params);
    return res;
}

double EFORTHeuristic(const sampling::Sample& sample,
                      const Eigen::Vector3d& direction, const Eigen::Vector3d& normal, const HeuristicParam & /*params*/)
{
    double EFORT = -direction.transpose() * sample.jacobianProduct_.block<3,3>(0,0) * (-direction);
    return EFORT * Eigen::Vector3d::UnitZ().dot(normal);
}

double EFORTNormalHeuristic(const sampling::Sample& sample,
                      const Eigen::Vector3d& direction, const Eigen::Vector3d& normal, const HeuristicParam & /*params*/)
{
    double EFORT = -direction.transpose() * sample.jacobianProduct_.block<3,3>(0,0) * (-direction);
    return EFORT * direction.dot(normal);
}

double ManipulabilityHeuristic(const sampling::Sample& sample,
                               const Eigen::Vector3d& /*direction*/, const Eigen::Vector3d& normal, const HeuristicParam & /*params*/)
{
    if(Eigen::Vector3d::UnitZ().dot(normal) < 0.7) return -1;
    return sample.staticValue_ * 10000 * Eigen::Vector3d::UnitZ().dot(normal) * 100000  +  ((double)rand()) / ((double)(RAND_MAX));
}

double RandomHeuristic(const sampling::Sample& /*sample*/,
                       const Eigen::Vector3d& /*direction*/, const Eigen::Vector3d& /*normal*/, const HeuristicParam & /*params*/)
{
    return ((double)rand()) / ((double)(RAND_MAX));
}


double ForwardHeuristic(const sampling::Sample& sample,
                      const Eigen::Vector3d& direction, const Eigen::Vector3d& normal, const HeuristicParam & /*params*/)
{
    return sample.staticValue_ * 10000 * Eigen::Vector3d::UnitZ().dot(normal) * 100  + sample.effectorPosition_.dot(fcl::Vec3f(direction(0),direction(1),direction(2))) + ((double)rand()) / ((double)(RAND_MAX));
}



double BackwardHeuristic(const sampling::Sample& sample,
                      const Eigen::Vector3d& direction, const Eigen::Vector3d& normal, const HeuristicParam & /*params*/)
{
    return sample.staticValue_ * 10000 * Eigen::Vector3d::UnitZ().dot(normal) * 100  - sample.effectorPosition_.dot(fcl::Vec3f(direction(0),direction(1),direction(2))) + ((double)rand()) / ((double)(RAND_MAX));
}

double StaticHeuristic(const sampling::Sample& sample,
                      const Eigen::Vector3d& /*direction*/, const Eigen::Vector3d& /*normal*/, const HeuristicParam & /*params*/)
{
    /*hppDout(info,"sample : ");
    hppDout(info,"sample : "<<&sample);
    hppDout(info,"id = "<<sample.id_);
    hppDout(info,"length = "<<sample.length_);
    hppDout(info,"startRank = "<<sample.startRank_);
    hppDout(info,"effectorPosition = "<<sample.effectorPosition_);
    hppDout(info,"configuration = "<<sample.configuration_);
    hppDout(info,"staticValue = "<<sample.staticValue_);
    */
    return sample.staticValue_;

}


double DistanceToLimitHeuristic(const sampling::Sample& sample,
                      const Eigen::Vector3d& /*direction*/, const Eigen::Vector3d& /*normal*/, const HeuristicParam & /*params*/)
{
    return sample.configuration_.norm();
}

HeuristicFactory::HeuristicFactory()
{
    unsigned int seed =  (unsigned int) (time(NULL)) ;
    //seed = 1485441926; // prepare_jump
    // seed = 1486147856; // stairs (18)
    //seed = 1486392757; // sideWall HyQ
    // seed = 1486721923; //hrp2 downSLope
   // seed = 1487082431; // hrp2 +0.15 z axis
    //seed = 1488532591; // downSLope (good contact but not interp) (yaml 1)
    // seed = 1488545915 ; //downSlope (need test)
    // seed = 1488550692 ; // downslope 2
    //seed = 1491571994; // straight walk 2 m static
    //seed = 1491580336 ; // straight walk dynamic 0.5
    //seed = 1492176551; // walk 0.3 !!
    //seed = 1493208163 ; // stairs static contacts
    std::cout<<"seed = "<<seed<<std::endl;
    srand ( seed);
    hppDout(notice,"SEED for heuristic = "<<seed);
    heuristics_.insert(std::make_pair("static", &StaticHeuristic));
    heuristics_.insert(std::make_pair("EFORT", &EFORTHeuristic));
    heuristics_.insert(std::make_pair("EFORT_Normal", &EFORTNormalHeuristic));
    heuristics_.insert(std::make_pair("manipulability", &ManipulabilityHeuristic));
    heuristics_.insert(std::make_pair("random", &RandomHeuristic));
    heuristics_.insert(std::make_pair("forward", &ForwardHeuristic));
    heuristics_.insert(std::make_pair("backward", &BackwardHeuristic));
    heuristics_.insert(std::make_pair("jointlimits", &DistanceToLimitHeuristic));
    heuristics_.insert(std::make_pair("ZMP", &ZMPHeuristic));
    heuristics_.insert(std::make_pair("straightCentroid", &StraightCentroidHeuristic));
    heuristics_.insert(std::make_pair("combined", &CombinedCentroidZMPHeuristic));
}

HeuristicFactory::~HeuristicFactory(){}

bool HeuristicFactory::AddHeuristic(const std::string& name, const heuristic func)
{
    if(heuristics_.find(name) != heuristics_.end())
        return false;
    heuristics_.insert(std::make_pair(name,func));
    return true;
}
