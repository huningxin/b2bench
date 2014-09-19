#include <cstdio>
#include <time.h>
#include <math.h>
#include <stdlib.h>

#include "ContactSolverScalar.h"
#include "ContactSolverSIMD.h"

#define DEBUG 0

#define WARMUP 2
#define FRAMES 64

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
	b2Position positions[2];

	for (int32 i = 0; i < COUNT; ++i) {
		b2ContactPositionConstraint* pc = pc_array + i;
		pc->indexA = 0;
		pc->indexB = 0;
		pc->pointCount = 2;
		pc->type = b2Manifold::e_faceA;
	}

	printf("warming...\n");
	for (int32 i = 0; i < WARMUP; ++i) {
		SolvePositionConstraints(count, pc_array, positions);
  	}

  	printf("benching...\n");
	clock_t times[FRAMES]; 
	for (int32 i = 0; i < FRAMES; ++i) {
		clock_t start = clock();
		for(int32 j = 0; j < 10; ++j)
		SolvePositionConstraints(count, pc_array, positions);
		clock_t end = clock();
		times[i] = end - start;
		printf("%f\n", (float)times[i] / CLOCKS_PER_SEC * 1000);
	}

  return measure(times);
}

result_t bench_simd() {
	int32 count = COUNT;
	b2ContactPositionConstraintSIMD pc;
	pc.type = b2Manifold::e_faceA;
	pc.pointCount = 2;
	pc.localPoints_x = new float32[COUNT * pc.pointCount];
	pc.localPoints_y = new float32[COUNT * pc.pointCount];
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

	printf("warming...\n");
	for (int32 i = 0; i < WARMUP; ++i) {
		SolvePositionConstraintsSIMD(count, &pc);
  	}

  	printf("benching...\n");
	clock_t times[FRAMES]; 
	for (int32 i = 0; i < FRAMES; ++i) {
		clock_t start = clock();
		for(int32 j = 0; j < 10; ++j)
		SolvePositionConstraintsSIMD(count, &pc);
		clock_t end = clock();
		times[i] = end - start;
		printf("%f\n", (float)times[i]);
	}

  return measure(times);
}

int main(int argc, char** argv) {
  result_t result = bench_scalar();
  printf("Scalar benchmark complete.\n  ms/frame: %f 5th %%ile: %f 95th %%ile: %f\n", result.mean, result.pc_5th, result.pc_95th);

  result_t result_simd = bench_simd();
  printf("SIMD benchmark complete.\n  ms/frame: %f 5th %%ile: %f 95th %%ile: %f\n", result_simd.mean, result_simd.pc_5th, result_simd.pc_95th);
  return 0;
}