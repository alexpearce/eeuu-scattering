#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/*
 * To compile:
 *   gcc -Wall test.c -o test && ./test
 *
 */

// Data structure that we use to return two 'arrays'
typedef struct {
  double *root_s;
  double *sigma;
} Data;

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

Data trapezium_cross_section(double (*f)(double), double *cross_section_ptr, double a, double b, double step_size, int strips) {
  int intervals = (b - a)/step_size;
    
  int i;
  double energies[intervals];
  double *energies_ptr = &energies[0];
  
  double new_energy = a;
  for (i = 0; i < intervals; i++) {
    // Set the energy
    set_collider_to(new_energy);
    *energies_ptr = new_energy;
    
    // Integrate
    *cross_section_ptr = trapezium((*f), -1.0, 1.0, strips);
    
    // Increase the energy for the next loop
    new_energy += step_size;
    
    cross_section_ptr++;
    energies_ptr++;
  }
  
  // Data struct
  Data rtn;
  
  // 'Reset' the pointers back to the beginning of their arrays
  rtn.root_s = energies_ptr - intervals;
  rtn.sigma  = cross_section_ptr - intervals;
  
  return rtn;
}

Data monte_carlo_cross_section(double (*f)(double), double *ptr, double a, double b, double step_size, int strips) {
  int intervals = (b - a)/step_size;
    
  int i;
  double energies[intervals];
  double new_energy = a;
  for (i = 0; i <= intervals; i++) {
    set_collider_to(new_energy);
    energies[i] = new_energy;
    // Increase the energy
    new_energy += step_size;
    
    *ptr = monte_carlo((*f), -1.0, 1.0, strips);
    ptr++;
  }
  
  Data rtn;
  rtn.root_s = &energies[0];
  rtn.sigma = ptr - (intervals - 1);
  
  return rtn;
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

// Writes an array of doubles to a file
// Expects `size` to be the array size
void write_to_file(Data dataset, int size) {
  // Pointers to the arrays
  double *root_s_ptr = dataset.root_s;
  double *sigma_ptr  = dataset.sigma;
  
  // File pointer
  FILE *file;
  
  // Open the file for writing
  file = fopen("data.csv", "w+");
  
  int i;
  for (i = 0; i < size; i++) {
    // Print a CSV line "x, y"
    fprintf(file, "%.15f, %.15f\n", *root_s_ptr, *sigma_ptr);
    root_s_ptr++;
    sigma_ptr++;
  }
  // Close the file
  fclose(file);
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
  double a = 3.0;
  double b = 200.0;
  double step_size = 0.1;
  int N = 1000;
  
  // How many steps will we take? This many.
  int intervals = (b - a)/step_size;
  
  // Create array of size `intervals`
  double cross_section[intervals];
  // Pointer to the array
  double *ptr;
  // Point to first element of the array
  ptr = &cross_section[0];
  
  Data dataset;
  dataset = trapezium_cross_section(z_z, &cross_section[0], a, b, step_size, N);
  
  // http://stackoverflow.com/a/4162948/596068
  int arr_length = sizeof(cross_section) / sizeof(cross_section[0]);
  
  write_to_file(dataset, arr_length);
  
  // int i;
  // for (i = 0; i < 6; i++) {
  //   printf("cross_section[%i] = %.15e\n", i, cross_section[i]);
  // }
  
  return 0;
}