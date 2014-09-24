#ifndef CONTACT_SOLVER_VECTOR_H
#define CONTACT_SOLVER_VECTOR_H

#include "b2Math.h"

typedef float float32x4  __attribute__((vector_size(16), aligned(1)));

struct b2ContactPositionConstraintVector
{
	float32x4 localPoints[b2_maxManifoldPoints];
	float32x4 localNormal;
	float32x4 localPoint;
	int32 indexA;
	int32 indexB;
	float32 invMassA, invMassB;
	float32x4 localCenterA, localCenterB;
	float32 invIA, invIB;
	b2Manifold::Type type;
	float32 radiusA, radiusB;
	int32 pointCount;
};

/// This is an internal structure.
struct b2PositionVector
{
	float32x4 c;
	float32 a;
};

inline float32x4 b2RotateVector(const float32x4& q, const float32x4& v)
{
	float32x4 tmp = {0, 1, 0, 0};
	return __builtin_shufflevector(q, q, 1, 0, 0, 0) * __builtin_shufflevector(v, v, 0, 0, 0, 0)
	       + q * __builtin_shufflevector(v, v, 1, 1, 0, 0) * tmp;
}

inline float32x4 b2MulTransformVector(const float32x4& q, const float32x4& p, const float32x4& v)
{
	return b2RotateVector(q, v) + p;
}

inline float32 b2NormalizeVector(float32x4& v)
{
	float32x4 tmp = v*v;
	float32 length = b2Sqrt(v[0] + v[1]);
	if (length < b2_epsilon)
	{
		return 0.0f;
	}
	float32 invLength = 1.0f / length;
	float32x4 invLengthx2 = {invLength, invLength, 0, 0};
	v *= invLengthx2;

	return length;
}

inline float32 b2DotVector(const float32x4& a, const float32x4& b)
{
	float32x4 tmp = a * b;
	return tmp[0] + tmp[1];
}

inline float32 b2CrossVector(const float32x4& a, const float32x4& b)
{
	float32x4 tmp = a * __builtin_shufflevector(b, b, 1, 0, 0, 0);
	return tmp[0] - tmp[1];
}

struct b2PositionSolverManifoldVector
{
	void Initialize(b2ContactPositionConstraintVector* pc,
		            const float32x4& xfA_q, const float32x4& xfA_p,
		            const float32x4& xfB_q, const float32x4& xfB_p,
		            int32 index)
	{
		b2Assert(pc->pointCount > 0);

		switch (pc->type)
		{
		case b2Manifold::e_circles:
			{
				float32x4 pointA = b2MulTransformVector(xfA_q, xfA_p, pc->localPoint);
				float32x4 pointB = b2MulTransformVector(xfB_q, xfB_p, pc->localPoints[0]);
				normal = pointB - pointA;
				b2NormalizeVector(normal);
				float32x4 tmp = {0.5f, 0.5f, 0, 0};
				point = tmp * (pointA + pointB);
				separation = b2DotVector(pointB - pointA, normal) - pc->radiusA - pc->radiusB;
			}
			break;

		case b2Manifold::e_faceA:
			{
				normal = b2RotateVector(xfA_q, pc->localNormal);
				float32x4 planePoint = b2MulTransformVector(xfA_q, xfA_p, pc->localPoint);

				float32x4 clipPoint = b2MulTransformVector(xfB_q, xfB_p, pc->localPoints[index]);
				separation = b2DotVector(clipPoint - planePoint, normal) - pc->radiusA - pc->radiusB;
				point = clipPoint;
			}
			break;

		case b2Manifold::e_faceB:
			{
				normal = b2RotateVector(xfB_q, pc->localNormal);
				float32x4 planePoint = b2MulTransformVector(xfB_q, xfB_p, pc->localPoint);

				float32x4 clipPoint = b2MulTransformVector(xfA_q, xfB_p, pc->localPoints[index]);
				separation = b2DotVector(clipPoint - planePoint, normal) - pc->radiusA - pc->radiusB;
				point = clipPoint;

				// Ensure normal points from A to B
				normal = -normal;
			}
			break;
		}
	}

	float32x4 normal;
	float32x4 point;
	float32 separation;
};

// Sequential solver.
bool SolvePositionConstraintsVector(
	int32 m_count, b2ContactPositionConstraintVector* m_positionConstraints, b2PositionVector* m_positions)
{
	float32 minSeparation = 0.0f;

	for (int32 i = 0; i < m_count; ++i)
	{
		b2ContactPositionConstraintVector* pc = m_positionConstraints + i;

		int32 indexA = pc->indexA;
		int32 indexB = pc->indexB;
		float32x4 localCenterA = pc->localCenterA;
		float32 mA = pc->invMassA;
		float32 iA = pc->invIA;
		float32x4 localCenterB = pc->localCenterB;
		float32 mB = pc->invMassB;
		float32 iB = pc->invIB;
		int32 pointCount = pc->pointCount;

		float32x4 cA = m_positions[indexA].c;
		float32 aA = m_positions[indexA].a;

		float32x4 cB = m_positions[indexB].c;
		float32 aB = m_positions[indexB].a;

		// Solve normal constraints
		for (int32 j = 0; j < pointCount; ++j)
		{
			float32x4 xfA_q = {sinf(aA), cosf(aA), 0, 0};
			float32x4 xfB_q = {sinf(aB), cosf(aB), 0, 0};
			float32x4 xfA_p = cA - b2RotateVector(xfA_q, localCenterA);
			float32x4 xfB_p = cB - b2RotateVector(xfB_q, localCenterB);

			b2PositionSolverManifoldVector psm;
			psm.Initialize(pc, xfA_q, xfA_p, xfB_q, xfB_p, j);
			float32x4 normal = psm.normal;

			float32x4 point = psm.point;
			float32 separation = psm.separation;

			float32x4 rA = point - cA;
			float32x4 rB = point - cB;

			// Track max constraint error.
			minSeparation = b2Min(minSeparation, separation);

			// Prevent large corrections and allow slop.
			float32 C = b2Clamp(b2_baumgarte * (separation + b2_linearSlop), -b2_maxLinearCorrection, 0.0f);

			// Compute the effective mass.
			float32 rnA = b2CrossVector(rA, normal);
			float32 rnB = b2CrossVector(rB, normal);
			float32 K = mA + mB + iA * rnA * rnA + iB * rnB * rnB;

			// Compute normal impulse
			float32 impulse = K > 0.0f ? - C / K : 0.0f;

			float32x4 impulse4 = {impulse, impulse, 0, 0};
			float32x4 P = impulse4 * normal;

			float32x4 mA4 = {mA, mA, 0, 0};
			cA -= mA4 * P;
			aA -= iA * b2CrossVector(rA, P);

			float32x4 mB4 = {mB, mB, 0, 0};
			cB += mB4 * P;
			aB += iB * b2CrossVector(rB, P);
		}

		m_positions[indexA].c = cA;
		m_positions[indexA].a = aA;

		m_positions[indexB].c = cB;
		m_positions[indexB].a = aB;
	}

	// We can't expect minSpeparation >= -b2_linearSlop because we don't
	// push the separation above -b2_linearSlop.
	return minSeparation >= -3.0f * b2_linearSlop;
}

#endif
