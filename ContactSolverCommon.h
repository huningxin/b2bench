#ifndef CONTACT_SOLVER_COMMON_H
#define CONTACT_SOLVER_COMMON_H

struct b2ManifoldPoint
{
	b2Vec2 localPoint;		///< usage depends on manifold type
	float32 normalImpulse;	///< the non-penetration impulse
	float32 tangentImpulse;	///< the friction impulse
	uint32 id;			///< uniquely identifies a contact point between two shapes
};

struct b2Manifold
{
	enum Type
	{
		e_circles,
		e_faceA,
		e_faceB
	};

	b2ManifoldPoint points[b2_maxManifoldPoints];	///< the points of contact
	b2Vec2 localNormal;								///< not use for Type::e_points
	b2Vec2 localPoint;								///< usage depends on manifold type
	Type type;
	int32 pointCount;								///< the number of manifold points
};

struct b2ContactPositionConstraint
{
	b2Vec2 localPoints[b2_maxManifoldPoints];
	b2Vec2 localNormal;
	b2Vec2 localPoint;
	int32 indexA;
	int32 indexB;
	float32 invMassA, invMassB;
	b2Vec2 localCenterA, localCenterB;
	float32 invIA, invIB;
	b2Manifold::Type type;
	float32 radiusA, radiusB;
	int32 pointCount;
};

struct b2ContactPositionConstraintSIMD
{
	float32* localPoints_x;
	float32* localPoints_y;
	float32* localNormal_x;
	float32* localNormal_y;
	float32* localPoint_x;
	float32* localPoint_y;
	float32* invMassA;
	float32* invMassB;
	float32* localCenterA_x;
	float32* localCenterA_y;
	float32* localCenterB_x;
	float32* localCenterB_y;
	float32* invIA;
	float32* invIB;
	float32* radiusA;
	float32* radiusB;
	float32* positionA_x;
	float32* positionA_y;
	float32* positionA_a;
	float32* positionB_x;
	float32* positionB_y;
	float32* positionB_a;
	b2Manifold::Type type;
	int32 pointCount;
};

/// This is an internal structure.
struct b2Position
{
	b2Vec2 c;
	float32 a;
};

#endif