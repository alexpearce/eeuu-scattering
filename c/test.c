#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Fine structure
float alpha;
// Weak mixing angle / Weinberg angle
float sin_2_theta_w;
float theta_w;
// e- coupling
float g_e;
// Z coupling
float g_z;
// Take electron to be massless
float m_e;
// Muon rest mass
float m_u;
// Z0 rest mass
float m_z;
// Z0 width
float gamma_z0;
// e- vectorial coupling
float C_e_v;
// e- axial coupling
float C_e_a;
// mu vectorial coupling
float C_u_v;
// mu axial coupling
float C_u_a;

// collider energy
float sqrt_s;
float s;
// energy variables
float epsilon;
float epsilon_2;
float zeta;

// Updates the global energy-dependant variables
// float energy in GeV
void set_collider_to(float energy) {
  sqrt_s = energy;
  s = sqrt_s * sqrt_s;
  epsilon = m_u / sqrt_s;
  epsilon_2 = epsilon * epsilon;
  zeta = m_z / sqrt_s;
}

double gamma_gamma(double cos_theta) {
  double epsilon_factor = 4.0*epsilon_2;
  
  // Already using alpha globally, so just to be sure call this alph
  double alph = (pow(g_e, 4.0) * sqrt(1 - epsilon_factor)) / (64.0*M_PI*M_PI*s);
  
  double cos_theta_2 = cos_theta * cos_theta;
  double sin_theta_2 = 1 - cos_theta_2;
  
  return alph*(1 + cos_theta_2 + (epsilon_factor*sin_theta_2));
}

double z_z(double cos_theta) {
  double cos_theta_2 = cos_theta * cos_theta;
  double sin_theta_2 = 1 - cos_theta_2;
  
  double epsilon_factor = 1 - (4.0*epsilon_2);
  double zeta_factor    = 1 - (zeta*zeta);
  
  double alph  = (pow(g_z, 4) * sqrt(epsilon_factor)) / (1024.0*M_PI*M_PI*s);
         alph /= (zeta_factor*zeta_factor) + (zeta*zeta*gamma_z0*gamma_z0 / s);
  double beta  = ((C_e_v*C_e_v) + (C_e_a*C_e_a))*(C_u_v)*(C_u_v);
  double gamma = ((C_e_v*C_e_v) + (C_e_a*C_e_a))*(C_u_a)*(C_u_a) * epsilon_factor;
  double delta = 8.0 * C_e_v * C_e_a * C_u_v * C_u_a * sqrt(epsilon_factor);

  
  return alph*(beta*(1 + cos_theta_2 + (4.0*epsilon_2*sin_theta_2)) + gamma*(1 + cos_theta_2) + delta*cos_theta);
}

double gamma_z(double cos_theta) {
  double cos_theta_2 = cos_theta * cos_theta;
  double sin_theta_2 = 1 - cos_theta_2;
  
  double epsilon_factor = sqrt(1 - epsilon_2);
  double zeta_factor    = 1 - (zeta*zeta);
  
  double alph  = (2.0 * g_e*g_e * g_z*g_z * epsilon_factor * zeta_factor) / (256.0*M_PI*M_PI*s);
         alph /= (zeta_factor*zeta_factor) + (zeta*zeta * gamma_z0*gamma_z0  / s);
  double beta  = C_e_v * C_u_v;
  double delta = 2.0 * C_e_a * C_u_a * epsilon_factor;
  
  return alph * (beta*(1 + cos_theta_2 + (4.0*epsilon_2*sin_theta_2)) + (delta*cos_theta));
}

// Numerically integrate a function f between a and b using N sampled points
double monte_carlo(double (*f)(double), double a, double b, int N) {
  // To call f, do (*f)(value)
  double difference = b - a;
  double prefactor  = difference / (double)N;
  double sum        = 0.0;
  
  int i;
  double x;
  for (i = 1; i <= N; i++) {
    // Random number in the interval [a, b], a > b
    x = (difference*drand48()) + a;
    sum += (*f)(x);
  }
  
  return prefactor*sum;
}

double trapezium(double (*f)(double), double a, double b, int N) {
  // Strip width
  double h         = (b - a)/(double)N;
  double prefactor = (h/2.0) * ((*f)(a) + (*f)(b));
  double sum       = 0.0;

  int i;
  for (i = 1; i <= N; i++) {
    sum += (*f)(a + ((double)i*h));
  }
  
  return prefactor + (h*sum);
}

double * trapezium_cross_section(double (*f)(double), double *ptr, double a, double b, double step_size, int strips) {
  int intervals = (b - a)/step_size;
  //double cross_section[intervals];
  
  int count = 0;
  int i;
  for (i = 0; i <= intervals; i++) {
    set_collider_to(a + (i * step_size));
    *ptr = trapezium((*f), -1.0, 1.0, strips);
    ptr++;
    count++;
  }
  
  return ptr;
}

// Seed the drand48 number generator
void seed_random(void) {
  unsigned int seed;
  FILE* urandom = fopen("/dev/urandom", "r");
  fread(&seed, sizeof(int), 1, urandom);
  fclose(urandom);
  srand48(seed);
}

double * ptr_fun(double *ptr) {
  double * ptr2;
  ptr2 = (ptr + 1);
  *ptr2 = 1.2;
  return ptr;
}

int main(void) {
  
  /* Constants */

  alpha         = 1.0 / 128.0;
  sin_2_theta_w = 0.23152;
  theta_w       = asin(sqrt(sin_2_theta_w));
  g_e           = sqrt(4*M_PI*alpha);
  g_z           = g_e / (cos(theta_w) * sin(theta_w));
  m_e           = 0.0;
  m_u           = 0.105658;
  m_z           = 91.1876;
  gamma_z0      = 2.5;
  C_e_v         = -0.5 + (2*sin_2_theta_w);
  C_e_a         = -0.5;
  C_u_v         = C_e_v;
  C_u_a         = C_e_a;
  
  /* End Constants */
  
  // Initialize the collider at 3GeV
  set_collider_to(3);
  
  // Sanity checks
  //printf("Trapezium:   %.15e\n", trapezium(gamma_gamma, -1.0, 1.0, 1000000));
  //printf("Monte Carlo: %.15e\n", monte_carlo(gamma_gamma, -1.0, 1.0, 10000));
  
  seed_random();
  
  // Set ranges and step sizes
  float a = 3;
  float b = 10;
  float step_size = 0.1;
  int N = 1000;
  
  // How many steps will we take? This many.
  int intervals = (b - a)/step_size;
  
  // Create array of size `intervals`
  double cross_section[intervals];
  // Pointer to the array
  double *ptr;
  // Point to first element of the array
  ptr = &cross_section[0];
  
  trapezium_cross_section(gamma_gamma, ptr, a, b, step_size, N);
  
  // int i;
  // for (i = 0; i < 6; i++) {
  //   printf("cross_section[%i] = %.15e\n", i, cross_section[i]);
  // }
  
  return 0;
}