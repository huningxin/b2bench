#include <cstdio>
#include <time.h>
#include <math.h>
#include <stdlib.h>

#include "ContactSolverScalar.h"
#include "ContactSolverSIMD.h"

#define DEBUG 0

#define WARMUP 2
#define FRAMES 64
#define ITERATIONS 1

#define COUNT 20000

typedef struct {
  float mean;
  float pc_5th;
  float pc_95th;
} result_t;

using namespace std;

const int e_count = 40;

// Simple nearest-rank %ile (on sorted array). We should have enough samples to make this reasonable.
float percentile(clock_t times[FRAMES], float pc) {
  int rank = (int)((pc * FRAMES) / 100);
  return times[rank];
}

int _cmp(const void *a, const void *b) {
	return *(clock_t*)a - *(clock_t*)b;
}

result_t measure(clock_t times[FRAMES]) {
  float values[FRAMES];
  result_t r;

	float total = 0;
	for (int i = 0; i < FRAMES; ++i) {
		values[i] = (float)times[i] / CLOCKS_PER_SEC * 1000;
	 	total += values[i];
	}
  r.mean = total / FRAMES;

	qsort(times, FRAMES, sizeof(clock_t), _cmp);
  r.pc_5th = percentile(times, 5) / CLOCKS_PER_SEC * 1000;
  r.pc_95th = percentile(times, 95) / CLOCKS_PER_SEC * 1000;
  return r;
}

result_t bench_scalar() {
	int32 count = COUNT;
	b2ContactPositionConstraint pc_array[COUNT];
	b2Position positions[COUNT];

	for (int32 i = 0; i < COUNT; ++i) {
		b2ContactPositionConstraint* pc = pc_array + i;
		pc->indexA = i;
		pc->indexB = i;
		pc->pointCount = 1;
		pc->type = b2Manifold::e_faceA;
		pc->localPoints[0] = b2Vec2(i, i);
		pc->localNormal = b2Vec2(i, i);
		pc->localPoint = b2Vec2(i, i);
		pc->invMassA = i;
		pc->invMassB = i;
		pc->localCenterA = b2Vec2(i, i);
		pc->localCenterB = b2Vec2(i, i);
		pc->invIA = i;
		pc->invIB = i;
		pc->radiusA = i;
		pc->radiusB = i;

		b2Position* p = positions + i;
		p->c = b2Vec2(i, i);
		p->a = i;
	}

	printf("warming...\n");
	for (int32 i = 0; i < WARMUP; ++i) {
		SolvePositionConstraints(count, pc_array, positions);
  	}

  	printf("benching...\n");
	clock_t times[FRAMES]; 
	for (int32 i = 0; i < FRAMES; ++i) {
		clock_t start = clock();
		for(int32 j = 0; j < ITERATIONS; ++j)
			SolvePositionConstraints(count, pc_array, positions);
		clock_t end = clock();
		times[i] = end - start;
		printf("%f\n", (float)times[i]);
	}

	float32 checksum = 0.0;
	for (int32 i = 0; i < COUNT; ++i) {
		b2Position* p = positions + i;
		checksum += p->c.x;
		checksum += p->c.y;
		checksum += p->a;
	}
	printf("checksum: %f\n", checksum);

  return measure(times);
}

result_t bench_simd() {
	int32 count = COUNT;
	b2ContactPositionConstraintSIMD pc;
	pc.type = b2Manifold::e_faceA;
	pc.pointCount = 1;
	pc.localPoints_x = new float32[COUNT];
	pc.localPoints_y = new float32[COUNT];
	pc.localNormal_x = new float32[COUNT];
	pc.localNormal_y = new float32[COUNT];
	pc.localPoint_x = new float32[COUNT];
	pc.localPoint_y = new float32[COUNT];
	pc.invMassA = new float32[COUNT];
	pc.invMassB = new float32[COUNT];
	pc.localCenterA_x = new float32[COUNT];
	pc.localCenterA_y = new float32[COUNT];
	pc.localCenterB_x = new float32[COUNT];
	pc.localCenterB_y = new float32[COUNT];
	pc.invIA = new float32[COUNT];
	pc.invIB = new float32[COUNT];
	pc.radiusA = new float32[COUNT];
	pc.radiusB = new float32[COUNT];
	pc.positionA_x = new float32[COUNT];
	pc.positionA_y = new float32[COUNT];
	pc.positionA_a = new float32[COUNT];
	pc.positionB_x = new float32[COUNT];
	pc.positionB_y = new float32[COUNT];
	pc.positionB_a = new float32[COUNT];

	for (int32 i = 0; i < COUNT; ++i) {
		pc.localPoints_x[i] = i;
		pc.localPoints_y[i] = i;
		pc.localNormal_x[i] = i;
		pc.localNormal_y[i] = i;
		pc.localPoint_x[i] = i;
		pc.localPoint_y[i] = i;
		pc.invMassA[i] = i;
		pc.invMassB[i] = i;
		pc.localCenterA_x[i] = i;
		pc.localCenterA_y[i] = i;
		pc.localCenterB_x[i] = i;
		pc.localCenterB_y[i] = i;
		pc.invIA[i] = i;
		pc.invIB[i] = i;
		pc.radiusA[i] = i;
		pc.radiusB[i] = i;
		pc.positionA_x[i] = i;
		pc.positionA_y[i] = i;
		pc.positionA_a[i] = i;
		pc.positionB_x[i] = i;
		pc.positionB_y[i] = i;
		pc.positionB_a[i] = i;
	}

	printf("warming...\n");
	for (int32 i = 0; i < WARMUP; ++i) {
		SolvePositionConstraintsSIMD(count, &pc);
  	}

  	printf("benching...\n");
	clock_t times[FRAMES]; 
	for (int32 i = 0; i < FRAMES; ++i) {
		clock_t start = clock();
		for(int32 j = 0; j < ITERATIONS; ++j)
			SolvePositionConstraintsSIMD(count, &pc);
		clock_t end = clock();
		times[i] = end - start;
		printf("%f\n", (float)times[i]);
	}

	float32 checksum = 0.0;
	for (int32 i = 0; i < COUNT; ++i) {
		checksum += pc.positionA_x[i];
		checksum += pc.positionA_y[i];
		checksum += pc.positionA_a[i];
		checksum += pc.positionB_x[i];
		checksum += pc.positionB_y[i];
		checksum += pc.positionB_a[i];
	}
	printf("checksum: %f\n", checksum);

  return measure(times);
}

int main(int argc, char** argv) {
  result_t result = bench_scalar();
  printf("Scalar benchmark complete.\n  ms/frame: %f 5th %%ile: %f 95th %%ile: %f\n", result.mean, result.pc_5th, result.pc_95th);

  result_t result_simd = bench_simd();
  printf("SIMD benchmark complete.\n  ms/frame: %f 5th %%ile: %f 95th %%ile: %f\n", result_simd.mean, result_simd.pc_5th, result_simd.pc_95th);
  return 0;
}