#ifndef HPP_HEURISTIC_TOOLS_HH
#define HPP_HEURISTIC_TOOLS_HH

#include <hpp/model/device.hh> // trick to get the includes of fcl, ...
#include <map>

namespace hpp{
namespace rbprm{
namespace sampling{
    
    /// Defines a parameters set for the ZMP-based heuristic
    struct HeuristicParam
    {
        std::map<std::string, fcl::Vec3f> contactPositions_; // to get the others contacts (without the considered sample)
        std::map<std::string, fcl::Vec3f> previousContactPositions_; // to get the contacts positions at the previous state
        fcl::Vec3f comAcceleration_; // The CoM acceleration
        fcl::Vec3f comPosition_; // The CoM position
        std::string sampleLimbName_; // The name of the considered sample
        fcl::Transform3f tfWorldRoot_; // The transform between the world coordinate system and the root of the robot
        bool lightVersion_; // To true if we don't want to consider z-CoM accelerations
        double g_; // The gravity acceleration

        HeuristicParam() : lightVersion_(false), g_(-9.80665) {}
        HeuristicParam(const std::map<std::string, fcl::Vec3f> & cp, const std::map<std::string, fcl::Vec3f> & ocp, const fcl::Vec3f & comAcc, const fcl::Vec3f & comPos,
                          const std::string & sln, const fcl::Transform3f & tf, bool lv = false);
        HeuristicParam(const HeuristicParam & zhp);

        HeuristicParam & operator=(const HeuristicParam & zhp);
    };

    /// Computes the transform of a point
    ///
    /// \param p The considered point
    /// \param tr The translation to apply
    /// \param ro The rotation to apply
    /// \return The transformed point
    fcl::Vec3f transform(const fcl::Vec3f & p, const fcl::Vec3f & tr, const fcl::Matrix3f & ro);

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

    /// Computes the support polygon
    ///
    /// \param contactPositions The map of the contact positions
    /// \return The support polygon (orthogonal projection of the contact positions in the ground plane)
    std::vector <Vec2D> computeSupportPolygon(const std::map <std::string, fcl::Vec3f> & contactPositions);

    /// Computes the convex hull of a set
    ///
    /// \param set The set we want to get the convex hull
    /// \return The convex hull of the specified set
    std::vector <Vec2D> convexHull(std::vector <Vec2D> set);

    /// Computes the weighted centroid of a convex polygon
    /// This is the "real (visual) center" of a polygon (an approximation of it in the worst case)
    ///
    /// \param convexPolygon The convex polygon to whom we want to find the weighted centroid
    /// \return The weighted centroid of the specified convex polygon
    Vec2D weightedCentroidConvex2D(const std::vector <Vec2D> & convexPolygon);

} // namespace sampling
} // namespace rbprm
} // namespace hpp

#endif