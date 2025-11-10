#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <cstdint>

float tdiff(struct timeval *start, struct timeval *end) {
  return (end->tv_sec - start->tv_sec) + 1e-6 * (end->tv_usec - start->tv_usec);
}


uint64_t s[4];

static inline uint64_t rotl(const uint64_t x, int k) {
  return (x << k) | (x >> (64 - k));
}

// unsigned long long seed = 100;

static uint64_t splitmix64(uint64_t *state) {
  uint64_t z = (*state += 0x9e3779b97f4a7c15ULL);
  z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
  z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
  return z ^ (z >> 31);
}

void initRNG(uint64_t seed) {
  uint64_t sm_state = seed;
  s[0] = splitmix64(&sm_state);
  s[1] = splitmix64(&sm_state);
  s[2] = splitmix64(&sm_state);
  s[3] = splitmix64(&sm_state);
}

uint64_t randomU64() {
  const uint64_t result = rotl(s[0] + s[3], 23) + s[0];
  const uint64_t t = s[1] << 17;
  
  s[2] ^= s[0];
  s[3] ^= s[1];
  s[1] ^= s[2];
  s[0] ^= s[3];
  s[2] ^= t;
  s[3] = rotl(s[3], 45);
  
  return result;
}

double randomDouble() {
  uint64_t x = randomU64();
  return (x >> 11) * 0x1.0p-53;
}

int L;          // Lattice size (L x L)
double T;       // Temperature
double J = 1.0; // Coupling constant
int **lattice;

double expCache[9];

void initializeCache() {
  for (int i = -4; i <= 4; i++) {
    int idx = i + 4;
    double dE = 2.0 * J * i;
    expCache[idx] = exp(-dE / T);
  }
}

void initializeLattice() {
  lattice = (int **)malloc(sizeof(int *) * L);
  for (int i = 0; i < L; i++) {
    lattice[i] = (int *)malloc(sizeof(int) * L);
    for (int j = 0; j < L; j++) {
      lattice[i][j] = (randomDouble() < 0.5) ? -1 : 1;
    }
  }
}

double calculateTotalEnergy() {
  double energy = 0.0;

  for (int i = 0; i < L; i++) {
    for (int j = 0; j < L; j++) {
      int spin = lattice[i][j];

      int up = lattice[(i - 1 + L) % L][j];
      int down = lattice[(i + 1) % L][j];
      int left = lattice[i][(j - 1 + L) % L];
      int right = lattice[i][(j + 1) % L];

      energy += -J * spin * (up + down + left + right);
    }
  }
  return 0.5 * energy;
}

double calculateMagnetization() {
  double mag = 0.0;
  for (int i = 0; i < L; i++) {
    for (int j = 0; j < L; j++) {
      mag += lattice[i][j];
    }
  }
  return mag / (L * L);
}

double calculateLocalEnergy(int i, int j) {
  int spin = lattice[i][j];
  int up    = lattice[(i - 1 + L) % L][j];
  int down  = lattice[(i + 1) % L][j];
  int left  = lattice[i][(j - 1 + L) % L];
  int right = lattice[i][(j + 1) % L];
  
  return 2 * J * spin * (up + down + left + right);
}
void metropolisHastingsStep() {
  int i = (int)(randomDouble() * L);
  int j = (int)(randomDouble() * L);

  lattice[i][j] *= -1;
  double dE = calculateLocalEnergy(i, j);

  if (dE <= 0.0) {
    return;
  }
  double prob = expCache[(int)dE+4];
  if (randomDouble() >= prob) {
    lattice[i][j] *= -1;
  }
}

void saveLatticeImage(const char *png_filename) {
  char ppm_filename[256];
  snprintf(ppm_filename, sizeof(ppm_filename), "temp_%s.ppm", png_filename);

  FILE *f = fopen(ppm_filename, "wb");
  if (!f) {
    printf("Error: Could not create temporary file %s\n", ppm_filename);
    return;
  }

  fprintf(f, "P6\n");
  fprintf(f, "%d %d\n", L, L);
  fprintf(f, "255\n");

  for (int i = 0; i < L; i++) {
    for (int j = 0; j < L; j++) {
      unsigned char r, g, b;
      if (lattice[i][j] == 1) {
        r = 255;
        g = 255;
        b = 255;
      } else {
        r = 0;
        g = 50;
        b = 200;
      }
      fwrite(&r, 1, 1, f);
      fwrite(&g, 1, 1, f);
      fwrite(&b, 1, 1, f);
    }
  }

  fclose(f);

  char cmd[512];
  snprintf(cmd, sizeof(cmd), "convert %s %s 2>/dev/null", ppm_filename,
           png_filename);
  int result = system(cmd);

  if (result == 0) {
    printf("Saved visualization to %s\n", png_filename);
    remove(ppm_filename);
  } else {
    rename(ppm_filename, png_filename);
    printf("Saved visualization to %s (install ImageMagick for PNG)\n",
           png_filename);
  }
}

void sanityCheck(double energy, double mag_per_spin, const char *stage) {
  double energy_per_spin = energy / (L * L);
  double Tc = 2.0 * J / log(1.0 + sqrt(2.0));

  printf("Sanity check [%s]:\n", stage);

  // 1. Energy per spin
  if (energy_per_spin < -2.0 * J - 0.01 || energy_per_spin > 2.0 * J + 0.01) {
    printf("  [ERROR] Energy per spin (%.4f) outside expected bounds "
           "[%.2f, %.2f]\n",
           energy_per_spin, -2.0 * J, 2.0 * J);
  } else {
    printf("  [OK] Energy per spin = %.4f (within bounds [%.2f, %.2f])\n",
           energy_per_spin, -2.0 * J, 2.0 * J);
  }

  // 2. Magnetization per spin
  if (fabs(mag_per_spin) > 1.01) {
    printf("  [ERROR] Magnetization per spin (%.4f) outside physical bounds "
           "[-1, 1]\n",
           mag_per_spin);
  } else {
    printf("  [OK] Magnetization per spin = %.4f (within bounds [-1, 1])\n",
           mag_per_spin);
  }

  printf("\n");
}

void freeLattice() {
  for (int i = 0; i < L; i++) {
    free(lattice[i]);
  }
  free(lattice);
}

int main(int argc, const char **argv) {
  if (argc < 4) {
    printf("Usage: %s <lattice_size> <temperature> <steps>\n", argv[0]);
    printf("Example: %s 100 2.269 10000000\n", argv[0]);
    printf("\n2D Ising Model\n");
    printf("Critical temperature: Tc = 2J/ln(1+√2) ≈ 2.26918531421\n");
    return 1;
  }
  initRNG(100); // for xoshiro init

  L = atoi(argv[1]);
  T = atof(argv[2]);
  int steps = atoi(argv[3]);

  printf("2D Ising Model\n");
  printf("=================================================\n");
  printf("Lattice size: %d x %d (%d spins)\n", L, L, L * L);
  printf("Temperature: T = %.4f (Tc ≈ 2.269)\n", T);
  printf("Coupling constant: J = %.2f\n", J);
  printf("Number of Metropolis-Hastings steps: %d\n", steps);
  printf("=================================================\n\n");

  initializeLattice();

  initializeCache();

  double initial_energy = calculateTotalEnergy();
  double initial_mag = calculateMagnetization();
  printf("Initial energy: %.4f\n", initial_energy);
  printf("Initial magnetization: %.4f\n\n", initial_mag);

  sanityCheck(initial_energy, initial_mag, "Initial state");

  saveLatticeImage("initial_state.png");

  struct timeval start, end;
  gettimeofday(&start, NULL);

  for (int step = 0; step < steps; step++) {
    metropolisHastingsStep();
  }

  gettimeofday(&end, NULL);

  double final_energy = calculateTotalEnergy();
  double final_mag = calculateMagnetization();

  printf("\nFinal energy: %.4f\n", final_energy);
  printf("Final magnetization: %.4f\n\n", final_mag);

  sanityCheck(final_energy, final_mag, "Final state");

  printf("Total time: %0.6f seconds\n\n", tdiff(&start, &end));

  saveLatticeImage("final_state.png");

  freeLattice();
  return 0;
}
