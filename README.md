Implementation of RB-PRM planner using hpp

BUG:
root mobile, transform root

TODO:

So far only handles free flyer robots (first joint is 3d translation, second is 
SO(3) )

Init value of HRP2 ankle really weird?

Configure tolerance individually for each effector (translation, rotation)

only samples until effector position 

Assumes 6 DOF abstraction

Add cone to state to avoid stability recomputation

optimize stability request (escpecially for the ocllision free case)

Point offset depends on rotation with normal when matching contact

Handle not only position but also rotation in effector contact point offset.
(especially when initializing default start and goal states in corba server)

Check rotation matrices given as a target for ik are the closest for current position

Tests for fullbody contact generation

Tests for interpolation path

Implement RBPRM objects for HRP2,
and an automatic generation program.

Implement and test stability criterion.

Uniformize use of Vec3f vs Eigen

avoid always taking collision objects into parameters

Adding documentation.
