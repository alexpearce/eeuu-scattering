#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

/*
 * To compile:
 *   gcc -Wall file.c -o file
 *
 */

// Data structure that we use to return two 'arrays'
typedef struct {
  float *root_s;
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

Data trapezium_cross_section(double (*f)(double), double *cross_section_ptr, float *energy_ptr, int intervals, int strips) {
  int i;
  for (i = 0; i < intervals; i++) {
    // Set the energy
    set_collider_to(*energy_ptr);
    // Integrate
    *cross_section_ptr = trapezium((*f), -1.0, 1.0, strips);
    // Next array element
    cross_section_ptr++;
    energy_ptr++;
  }
  
  // Data struct
  Data rtn;
  
  // 'Reset' the pointers back to the beginning of their arrays
  rtn.root_s = energy_ptr - intervals;
  rtn.sigma  = cross_section_ptr - intervals;
  
  return rtn;
}

Data monte_carlo_cross_section(double (*f)(double), double *cross_section_ptr, float *energy_ptr, int intervals, int strips) {
  int i;
  for (i = 0; i < intervals; i++) {
    // Set the energy
    set_collider_to(*energy_ptr);
    // Integrate
    *cross_section_ptr = monte_carlo((*f), -1.0, 1.0, strips);
    // Next array element
    cross_section_ptr++;
    energy_ptr++;
  }
  
  // Data struct
  Data rtn;
  
  // 'Reset' the pointers back to the beginning of their arrays
  rtn.root_s = energy_ptr - intervals;
  rtn.sigma  = cross_section_ptr - intervals;
  
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
  float *root_s_ptr = dataset.root_s;
  double *sigma_ptr  = dataset.sigma;
  
  // File pointer
  FILE *file;
  
  // Open the file for writing
  file = fopen("data.csv", "w");
  
  int i;
  for (i = 0; i < size; i++) {
    // Print a CSV line "x, y"
    // fprintf(file, "%.3f, %.15f\n", *root_s_ptr, *sigma_ptr);
    fprintf(file, "%.2f, %.15f\n", *root_s_ptr, *sigma_ptr);
    root_s_ptr++;
    sigma_ptr++;
  }
  // Close the file
  fclose(file);
}

// Calculate the total cross section
void cross_section(int export_to_file) {
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
  
  // Set ranges and step sizes
  int a = 3;
  int b = 200;
  float step_size = 0.1;
  int N = 1000;
  
  // How many steps will we take? This many.
  int intervals = (float)(b - a)/step_size;
  
  // Set up the energies
  float root_s[intervals];
  root_s[0] = a;
  int i;
  for (i = 1; i < intervals; i++) {
    // Set the energy
    root_s[i] = a + i*step_size;
  }
  
  // Create array of size `intervals`
  double g_g_sigma[intervals], z_z_sigma[intervals], g_z_sigma[intervals], total_sigma[intervals];
  
  Data g_g_data, z_z_data, g_z_data, total_data;
  // g_g_data = trapezium_cross_section(gamma_gamma, &g_g_sigma[0], &root_s[0], intervals, N);
  // z_z_data = trapezium_cross_section(gamma_z, &z_z_sigma[0], &root_s[0], intervals, N);
  // g_z_data = trapezium_cross_section(z_z, &g_z_sigma[0], &root_s[0], intervals, N);
  g_g_data = monte_carlo_cross_section(gamma_gamma, &g_g_sigma[0], &root_s[0], intervals, N);
  z_z_data = monte_carlo_cross_section(gamma_z, &z_z_sigma[0], &root_s[0], intervals, N);
  g_z_data = monte_carlo_cross_section(z_z, &g_z_sigma[0], &root_s[0], intervals, N);
  
  // http://stackoverflow.com/a/4162948/596068
  int arr_length = sizeof(root_s) / sizeof(root_s[0]);
  
  // Prepare the total data
  // i.e. sum the individual cross sections
  double *g_g_sigma_ptr = g_g_data.sigma;
  double *z_z_sigma_ptr = z_z_data.sigma;
  double *g_z_sigma_ptr = g_z_data.sigma;
  
  // We have an iterator just up there ^
  for (i = 0; i < arr_length; i++) {
    total_sigma[i] = *g_g_sigma_ptr + *z_z_sigma_ptr + *g_z_sigma_ptr;
    
    g_g_sigma_ptr++;
    z_z_sigma_ptr++;
    g_z_sigma_ptr++;
  }
  
  total_data.sigma = &total_sigma[0];
  total_data.root_s = &root_s[0];
  
  if (export_to_file == 1) {
    write_to_file(total_data, arr_length);
  }
}

int main(int argc, char *argv[]) {
  
  int i, should_time, export_to_file = 0;
  int times = 1;
  for (i = 0; i < argc; i++) {
    if (strcmp(argv[i], "--export") == 0) {
      // Export the cross section to file
      export_to_file = 1;
    } else if (strcmp(argv[i], "--time") == 0) {
      if (atoi(argv[i+1]) == 0) {
        // TODO this causes a segfault if i+1 isn't an argument
        printf("--time argument must have an integer value.\n");
      } else {
        // Time the cross section function N times
        should_time = 1;
        // Number of times to repeat the timing
        times = atoi(argv[i+1]);
      }
    }
  }
  
  // Initialize the collider at 3GeV
  set_collider_to(3);
  
  // Sanity checks
  //printf("Trapezium:   %.15e\n", trapezium(gamma_gamma, -1.0, 1.0, 1000000));
  //printf("Monte Carlo: %.15e\n", monte_carlo(gamma_gamma, -1.0, 1.0, 10000));
  
  seed_random();
  
  if (should_time) {
    double start, end;
    int j;
    start = clock();
    for (j = 0; j < times; j++) {
      // Don't export to file
      cross_section(0);
    }
    end = clock();
    
    float diff = (end - start)/CLOCKS_PER_SEC;
    printf("Total time for %i trials was %.3f seconds.\n", times, diff);
    printf("The average for one trial was %.3f seconds.\n", diff/(float)times);
  } else {
    cross_section(export_to_file);
  }
  
  return 0;
}