#ifndef CONTACT_SOLVER_SIMD_H
#define CONTACT_SOLVER_SIMD_H

#include "b2Math.h"
#include "ContactSolverCommon.h"

#include <xmmintrin.h>

#define B2_DEBUG_SIMD 0

float32 SolvePositionConstraintsSIMD(int32 count, b2ContactPositionConstraintSIMD* pc)
{
	static const __m128 SIGNMASK = _mm_castsi128_ps(_mm_set1_epi32(0x80000000));

	// pc_array contains pcs with same type and pointsCount
	b2Manifold::Type type = pc->type;
	int32 pointCount = pc->pointCount;

	__m128 minSeparation4 = _mm_set_ps1(0.0f);

	for (int32 i = 0; i < count; i += 4)
	{
		// SIMD load
		//b2Vec2 localCenterA = pc->localCenterA;
		__m128 localCenterA_x4 = _mm_loadu_ps(pc->localCenterA_x + i);
		__m128 localCenterA_y4 = _mm_loadu_ps(pc->localCenterA_y + i);
		//float32 mA = pc->invMassA;
		__m128 mA4 = _mm_loadu_ps(pc->invMassA + i);	
		//float32 iA = pc->invIA;
		__m128 iA4 = _mm_loadu_ps(pc->invIA + i);
		//b2Vec2 localCenterB = pc->localCenterB;
		__m128 localCenterB_x4 = _mm_loadu_ps(pc->localCenterB_x + i);
		__m128 localCenterB_y4 = _mm_loadu_ps(pc->localCenterB_y + i);
		//float32 mB = pc->invMassB;
		__m128 mB4 = _mm_loadu_ps(pc->invMassB + i);
		//float32 iB = pc->invIB;
		__m128 iB4 = _mm_loadu_ps(pc->invIB + i);
		//b2Vec2 cA = m_positions[indexA].c;
		__m128 cA_x4 = _mm_loadu_ps(pc->positionA_x + i);
		__m128 cA_y4 = _mm_loadu_ps(pc->positionA_y + i);
		//float32 aA = m_positions[indexA].a;
		__m128 aA4 = _mm_loadu_ps(pc->positionA_a + i);
		//b2Vec2 cB = m_positions[indexB].c;
		__m128 cB_x4 = _mm_loadu_ps(pc->positionB_x + i);
		__m128 cB_y4 = _mm_loadu_ps(pc->positionB_y + i);
		//float32 aB = m_positions[indexB].a;
		__m128 aB4 = _mm_loadu_ps(pc->positionB_a + i);

		__m128 localNormal_x4 = _mm_loadu_ps(pc->localNormal_x + i);
		__m128 localNormal_y4 = _mm_loadu_ps(pc->localNormal_y + i);
		__m128 localPoint_x4 = _mm_loadu_ps(pc->localPoint_x + i);
		__m128 localPoint_y4 = _mm_loadu_ps(pc->localPoint_y + i);
		__m128 radiusA4 = _mm_loadu_ps(pc->radiusA + i);
		__m128 radiusB4 = _mm_loadu_ps(pc->radiusB + i);

#if B2_DEBUG_SIMD == 1
		printf("input - cA.x = %f cA.y = %f aA = %f cB.x = %f cB.y = %f aB = %f\n", cA_x4[0], cA_y4[0], aA4[0], cB_x4[0], cB_y4[0], aB4[0]);
		printf("input - cA.x = %f cA.y = %f aA = %f cB.x = %f cB.y = %f aB = %f\n", cA_x4[1], cA_y4[1], aA4[1], cB_x4[1], cB_y4[1], aB4[1]);
		printf("input - cA.x = %f cA.y = %f aA = %f cB.x = %f cB.y = %f aB = %f\n", cA_x4[2], cA_y4[2], aA4[2], cB_x4[2], cB_y4[2], aB4[2]);
		printf("input - cA.x = %f cA.y = %f aA = %f cB.x = %f cB.y = %f aB = %f\n", cA_x4[3], cA_y4[3], aA4[3], cB_x4[3], cB_y4[3], aB4[3]);
#endif

		for (int32 j = 0; j < pointCount; ++j)
		{
			// xfA.q.Set(aA);
			__m128 xfA_q_s4, xfA_q_c4;
			xfA_q_s4[0] = sinf(aA4[0]);
			xfA_q_s4[1] = sinf(aA4[1]);
			xfA_q_s4[2] = sinf(aA4[2]);
			xfA_q_s4[3] = sinf(aA4[3]);
			xfA_q_c4[0] = cosf(aA4[0]);
			xfA_q_c4[1] = cosf(aA4[1]);
			xfA_q_c4[2] = cosf(aA4[2]);
			xfA_q_c4[3] = cosf(aA4[3]);

			// xfB.q.Set(aB);
			__m128 xfB_q_s4, xfB_q_c4;
			xfB_q_s4[0] = sinf(aB4[0]);
			xfB_q_s4[1] = sinf(aB4[1]);
			xfB_q_s4[2] = sinf(aB4[2]);
			xfB_q_s4[3] = sinf(aB4[3]);
			xfB_q_c4[0] = cosf(aB4[0]);
			xfB_q_c4[1] = cosf(aB4[1]);
			xfB_q_c4[2] = cosf(aB4[2]);
			xfB_q_c4[3] = cosf(aB4[3]);

			// xfA.p = cA - b2Mul(xfA.q, localCenterA);
			__m128 xfA_p_x4, xfA_p_y4;
			xfA_p_x4 = _mm_sub_ps(cA_x4, _mm_sub_ps(_mm_mul_ps(xfA_q_c4, localCenterA_x4),
												    _mm_mul_ps(xfA_q_s4, localCenterA_y4)));
			xfA_p_y4 = _mm_sub_ps(cA_y4, _mm_add_ps(_mm_mul_ps(xfA_q_s4, localCenterA_x4),
												    _mm_mul_ps(xfA_q_c4, localCenterA_y4)));

			// xfB.p = cB - b2Mul(xfB.q, localCenterB);
			__m128 xfB_p_x4, xfB_p_y4;
			xfB_p_x4 = _mm_sub_ps(cB_x4, _mm_sub_ps(_mm_mul_ps(xfB_q_c4, localCenterB_x4),
												    _mm_mul_ps(xfB_q_s4, localCenterB_y4)));
			xfB_p_y4 = _mm_sub_ps(cB_y4, _mm_add_ps(_mm_mul_ps(xfB_q_s4, localCenterB_x4),
												    _mm_mul_ps(xfB_q_c4, localCenterB_y4)));

#if B2_DEBUG_SIMD == 1
			printf("xfA.p.x = %f xfA.p.y = %f xfB.p.x = %f xfB.p.y = %f\n", xfA_p_x4[0], xfA_p_y4[0], xfB_p_x4[0], xfB_p_y4[0]);
			printf("xfA.p.x = %f xfA.p.y = %f xfB.p.x = %f xfB.p.y = %f\n", xfA_p_x4[1], xfA_p_y4[1], xfB_p_x4[1], xfB_p_y4[1]);
			printf("xfA.p.x = %f xfA.p.y = %f xfB.p.x = %f xfB.p.y = %f\n", xfA_p_x4[2], xfA_p_y4[2], xfB_p_x4[2], xfB_p_y4[2]);
			printf("xfA.p.x = %f xfA.p.y = %f xfB.p.x = %f xfB.p.y = %f\n", xfA_p_x4[3], xfA_p_y4[3], xfB_p_x4[3], xfB_p_y4[3]);
#endif

			// b2PositionSolverManifold psm;
			// psm.Initialize(pc, xfA, xfB, j);
			//b2Vec2 normal = psm.normal;

			//b2Vec2 point = psm.point;
			//float32 separation = psm.separation;
			__m128 normal_x4, normal_y4;
			__m128 point_x4, point_y4;
			__m128 separation4;
			if (type == b2Manifold::e_faceA) {
				// normal = b2Mul(xfA.q, pc->localNormal);
				normal_x4 = _mm_sub_ps(_mm_mul_ps(xfA_q_c4, localNormal_x4),
									   _mm_mul_ps(xfA_q_s4, localNormal_y4));
				normal_y4 = _mm_add_ps(_mm_mul_ps(xfA_q_s4, localNormal_x4),
									   _mm_mul_ps(xfA_q_c4, localNormal_y4));

#if B2_DEBUG_SIMD == 1
				printf("normal.x = %f normal.y = %f\n", normal_x4[0], normal_y4[0]);
				printf("normal.x = %f normal.y = %f\n", normal_x4[1], normal_y4[1]);
				printf("normal.x = %f normal.y = %f\n", normal_x4[2], normal_y4[2]);
				printf("normal.x = %f normal.y = %f\n", normal_x4[3], normal_y4[3]);
#endif

				// b2Vec2 planePoint = b2Mul(xfA, pc->localPoint);
				__m128 planePoint_x4, planePoint_y4;
				planePoint_x4 = _mm_add_ps(_mm_sub_ps(_mm_mul_ps(xfA_q_c4, localPoint_x4),
													  _mm_mul_ps(xfA_q_s4, localPoint_y4)),
										   xfA_p_x4);
				planePoint_y4 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(xfA_q_s4, localPoint_x4),
													  _mm_mul_ps(xfA_q_c4, localPoint_y4)),
										   xfA_p_y4);

				// b2Vec2 clipPoint = b2Mul(xfB, pc->localPoints[index]);
				__m128 clipPoint_x4, clipPoint_y4;
				__m128 localPoints_x4 = _mm_loadu_ps(pc->localPoints_x + i + j);
				__m128 localPoints_y4 = _mm_loadu_ps(pc->localPoints_y + i + j);
				
				clipPoint_x4 = _mm_add_ps(_mm_sub_ps(_mm_mul_ps(xfB_q_c4, localPoints_x4),
													 _mm_mul_ps(xfB_q_s4, localPoints_y4)),
										  xfB_p_x4);
				clipPoint_y4 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(xfB_q_s4, localPoints_x4),
													 _mm_mul_ps(xfB_q_c4, localPoints_y4)),
										  xfB_p_y4);

				// separation = b2Dot(clipPoint - planePoint, normal) - pc->radiusA - pc->radiusB;
				__m128 temp_x4 = _mm_sub_ps(clipPoint_x4, planePoint_x4);
				__m128 temp_y4 = _mm_sub_ps(clipPoint_y4, planePoint_y4);
				separation4 = _mm_sub_ps(_mm_sub_ps(_mm_add_ps(_mm_mul_ps(temp_x4, normal_x4),
															   _mm_mul_ps(temp_y4, normal_y4)),
												    radiusA4),
										 radiusB4);
				// point = clipPoint;
				point_x4 = clipPoint_x4;
				point_y4 = clipPoint_y4;
			} else if (type == b2Manifold::e_faceB) {
				//normal = b2Mul(xfB.q, pc->localNormal);
				normal_x4 = _mm_sub_ps(_mm_mul_ps(xfB_q_c4, localNormal_x4),
									   _mm_mul_ps(xfB_q_s4, localNormal_y4));
				normal_y4 = _mm_add_ps(_mm_mul_ps(xfB_q_s4, localNormal_x4),
									   _mm_mul_ps(xfB_q_c4, localNormal_y4));

				//b2Vec2 planePoint = b2Mul(xfB, pc->localPoint);
				__m128 planePoint_x4, planePoint_y4;
				planePoint_x4 = _mm_add_ps(_mm_sub_ps(_mm_mul_ps(xfB_q_c4, localPoint_x4),
													  _mm_mul_ps(xfB_q_s4, localPoint_y4)),
										   xfB_p_x4);
				planePoint_y4 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(xfB_q_s4, localPoint_x4),
													  _mm_mul_ps(xfB_q_c4, localPoint_y4)),
										   xfB_p_y4);
				//b2Vec2 clipPoint = b2Mul(xfA, pc->localPoints[index]);
				__m128 clipPoint_x4, clipPoint_y4;
				__m128 localPoints_x4 = _mm_loadu_ps(pc->localPoints_x+ i + j);
				__m128 localPoints_y4 = _mm_loadu_ps(pc->localPoints_y + i + j);
				
				clipPoint_x4 = _mm_add_ps(_mm_sub_ps(_mm_mul_ps(xfA_q_c4, localPoints_x4),
													 _mm_mul_ps(xfA_q_s4, localPoints_y4)),
										  xfB_p_x4);
				clipPoint_y4 = _mm_add_ps(_mm_add_ps(_mm_mul_ps(xfA_q_s4, localPoints_x4),
													 _mm_mul_ps(xfA_q_c4, localPoints_y4)),
										  xfB_p_y4);
				//separation = b2Dot(clipPoint - planePoint, normal) - pc->radiusA - pc->radiusB;
				__m128 temp_x4 = _mm_sub_ps(clipPoint_x4, planePoint_x4);
				__m128 temp_y4 = _mm_sub_ps(clipPoint_y4, planePoint_y4);
				separation4 = _mm_sub_ps(_mm_sub_ps(_mm_add_ps(_mm_mul_ps(temp_x4, normal_x4),
															   _mm_mul_ps(temp_y4, normal_y4)),
												    radiusA4),
										 radiusB4);
				//point = clipPoint;
				point_x4 = clipPoint_x4;
				point_y4 = clipPoint_y4;

				//normal = -normal;
				normal_x4 = _mm_xor_ps(normal_x4, SIGNMASK);
				normal_y4 = _mm_xor_ps(normal_y4, SIGNMASK);
			} else {
				b2Assert(false);
			}

			// b2Vec2 rA = point - cA;
			__m128 rA_x4 = _mm_sub_ps(point_x4, cA_x4);
			__m128 rA_y4 = _mm_sub_ps(point_y4, cA_y4);

			// b2Vec2 rB = point - cB;
			__m128 rB_x4 = _mm_sub_ps(point_x4, cB_x4);
			__m128 rB_y4 = _mm_sub_ps(point_y4, cB_y4);

			// minSeparation = b2Min(minSeparation, separation);
			minSeparation4 = _mm_min_ps(minSeparation4, separation4);

			// float32 C = b2Clamp(b2_baumgarte * (separation + b2_linearSlop), -b2_maxLinearCorrection, 0.0f);
			__m128 lowx4 = _mm_set_ps1(-b2_maxLinearCorrection);
			__m128 zero4 = _mm_set_ps1(0.0f);
			__m128 b2_baumgarte4 = _mm_set_ps1(b2_baumgarte);
			__m128 b2_linearSlop4 = _mm_set_ps1(b2_linearSlop);
			__m128 C4 = _mm_max_ps(lowx4, _mm_min_ps(_mm_mul_ps(b2_baumgarte4, _mm_add_ps(separation4, b2_linearSlop4)), zero4));

			// float32 rnA = b2Cross(rA, normal);
			__m128 rnA4 = _mm_sub_ps(_mm_mul_ps(rA_x4, normal_y4),
				                     _mm_mul_ps(rA_y4, normal_x4));
			// float32 rnB = b2Cross(rB, normal);
			__m128 rnB4 = _mm_sub_ps(_mm_mul_ps(rB_x4, normal_y4),
				                     _mm_mul_ps(rB_y4, normal_x4));
			// float32 K = mA + mB + iA * rnA * rnA + iB * rnB * rnB;
			__m128 K4 = _mm_add_ps(_mm_add_ps(_mm_add_ps(mA4, mB4),
				                              _mm_mul_ps(_mm_mul_ps(rnA4, rnA4), iA4)),
								   _mm_mul_ps(_mm_mul_ps(rnB4, rnB4), iB4));

			// float32 impulse = K > 0.0f ? - C / K : 0.0f;
			__m128 trueValue = _mm_div_ps(C4, K4);
			// negate it.
			trueValue = _mm_xor_ps(trueValue, SIGNMASK);

			 __m128 mask = _mm_cmpgt_ps(K4, zero4);
			__m128 impulse4 = _mm_or_ps(_mm_and_ps(mask, trueValue),
										_mm_andnot_ps(mask, zero4));

			// b2Vec2 P = impulse * normal;
			__m128 P_x4 = _mm_mul_ps(normal_x4, impulse4);
			__m128 P_y4 = _mm_mul_ps(normal_y4, impulse4);

#if B2_DEBUG_SIMD == 1
			printf("P.x = %f P.y = %f\n", P_x4[0], P_y4[0]);
			printf("P.x = %f P.y = %f\n", P_x4[1], P_y4[1]);
			printf("P.x = %f P.y = %f\n", P_x4[2], P_y4[2]);
			printf("P.x = %f P.y = %f\n", P_x4[3], P_y4[3]);
#endif

			// cA -= mA * P;
			cA_x4 = _mm_sub_ps(cA_x4, _mm_mul_ps(mA4, P_x4));
			cA_y4 = _mm_sub_ps(cA_y4, _mm_mul_ps(mA4, P_y4));

			// aA -= iA * b2Cross(rA, P);
			aA4 = _mm_sub_ps(aA4, _mm_mul_ps(iA4, _mm_sub_ps(_mm_mul_ps(rA_x4, P_y4), _mm_mul_ps(rA_y4, P_x4))));

			// cB += mB * P;
			cB_x4 = _mm_add_ps(cB_x4, _mm_mul_ps(mB4, P_x4));
			cB_y4 = _mm_add_ps(cB_y4, _mm_mul_ps(mB4, P_y4));

			// aB += iB * b2Cross(rB, P);
			aB4 = _mm_add_ps(aB4, _mm_mul_ps(iB4, _mm_sub_ps(_mm_mul_ps(rB_x4, P_y4), _mm_mul_ps(rB_y4, P_x4))));
		}

		// SIMD store
		_mm_storeu_ps(pc->positionA_x + i, cA_x4);
		_mm_storeu_ps(pc->positionA_y + i, cA_y4);
		_mm_storeu_ps(pc->positionA_a + i, aA4);
		_mm_storeu_ps(pc->positionB_x + i, cB_x4);
		_mm_storeu_ps(pc->positionB_y + i, cB_y4);
		_mm_storeu_ps(pc->positionB_a + i, aB4);

#if B2_DEBUG_SIMD == 1
		printf("output - cA.x = %f cA.y = %f aA = %f cB.x = %f cB.y = %f aB = %f\n",  cA_x4[0], cA_y4[0], aA4[0], cB_x4[0], cB_y4[0], aB4[0]);
		printf("output - cA.x = %f cA.y = %f aA = %f cB.x = %f cB.y = %f aB = %f\n",  cA_x4[1], cA_y4[1], aA4[1], cB_x4[1], cB_y4[1], aB4[1]);
		printf("output - cA.x = %f cA.y = %f aA = %f cB.x = %f cB.y = %f aB = %f\n",  cA_x4[2], cA_y4[2], aA4[2], cB_x4[2], cB_y4[2], aB4[2]);
		printf("output - cA.x = %f cA.y = %f aA = %f cB.x = %f cB.y = %f aB = %f\n",  cA_x4[3], cA_y4[3], aA4[3], cB_x4[3], cB_y4[3], aB4[3]);
#endif
	}

	float32 minSeparation = 0.0f;
	minSeparation = b2Min(minSeparation, minSeparation4[0]);
	minSeparation = b2Min(minSeparation,  minSeparation4[1]);
	minSeparation = b2Min(minSeparation,  minSeparation4[2]);
	minSeparation = b2Min(minSeparation,  minSeparation4[3]);
	return minSeparation;
}

#endif