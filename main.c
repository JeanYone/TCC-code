// Libraries that provide necessary data ypes and functions
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>
// External Library used for writing images in PNG format
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"
// Defines LATYPE as int and represent the data type used for each point in the magnetic lattice.
#define LATYPE int // 
typedef unsigned long long u64;
typedef long long i64;
// A data type that stores other variables within it, so when it is "called," all the variables contained within it are read.
typedef struct {
    double energy, magnetization;
} energy_magnetization_pack;
// Same as above
typedef struct {
    double avg_E, avg_M, specific_heat, magnetic_sus; //ඞඞඞඞඞඞඞඞඞඞඞ
    double avg_E2, avg_M2;
} ising_data;
// a função intrinseca de numeros randomicos rand()
// nesse caso o numero aleátorio gerado esta entre 0 e 1
double rrand() {
    return (double)rand() / (double)RAND_MAX;
}
// initializes a random lattice by allocating memory for the lattice and 
// assigning random spin values (+1 or -1) to each lattice point. 
LATYPE* init_random_lattice(u64 n) {
    LATYPE *lattice = (LATYPE*)calloc(n * n, sizeof(LATYPE));
    for (u64 i = 0; i < n * n; ++i)
        lattice[i] = 2 * (rand() % 2) - 1;

    return lattice;
}
// unction takes an existing lattice and randomizes its spin values.
void randomize_lattice(LATYPE *lattice, u64 n) {
    for (u64 i = 0; i < n * n; ++i)
        lattice[i] = 2 * (rand() % 2) - 1;
}
// ensures that lattice indices wrap around the boundaries of the lattice, allowing periodic boundary conditions.
i64 bounds(i64 i, u64 n) {
    if (0 <= i && (u64)i < n) return i;
    return ((i % n) + n) % n;
}
// calculates the change in energy
double dE_to_flip(LATYPE S, double J, LATYPE sum_neigh, double H) {
    return 2.0 * S * (J * sum_neigh + H);
}
// calculates the sum of the nearest neighbor spins for a given lattice point.
LATYPE neighbor_sum(LATYPE *lattice, i64 i, i64 j, int n) {
    LATYPE left  = lattice[i * n + bounds(j - 1, n)];
    LATYPE right = lattice[i * n + bounds(j + 1, n)];
    LATYPE down  = lattice[bounds(i - 1, n) * n + j];
    LATYPE up    = lattice[bounds(i + 1, n) * n + j];
    return left + right + down + up;
}
// calculates the total energy and magnetization of the lattice by iterating
// over all lattice points and summing the corresponding contributions.
energy_magnetization_pack lattice_energy_magnetization(LATYPE *lattice, u64 n, double J, double H) {
    energy_magnetization_pack to_return = {0};
    for (i64 row = 0; (u64)row < n; ++row) {
        for (i64 col = 0; (u64)col < n; ++col) {
            to_return.energy += -H * lattice[row * n + col] - J * neighbor_sum(lattice, row, col, n) * lattice[row * n + col] / 4.0;
            to_return.magnetization += lattice[row * n + col];
        }
    }
    return to_return;
}
// Performs a single step of the Metropolis algorithm for thermalization. It randomly selects a lattice point, calculates the change in energy 
// if the spin were flipped, and decides whether to flip the spin based on the Metropolis criterion.
void termal_step(LATYPE *lattice, u64 n, double J, double H, double T) {

    for (u64 row = 0; row < n; ++row) {
        for (u64 col = 0; col < n; ++col) {
            u64 row_random = rand() % n;
            u64 col_random = rand() % n;

            double dE = 2.0 * lattice[row_random * n + col_random] * (J * neighbor_sum(lattice, row_random, col_random, n) + H);
            if (dE < 0)
                lattice[row_random * n + col_random] *= -1;
            else if (rrand() < exp(-dE / T))
                lattice[row_random * n + col_random] *= -1;
        }
    }
}
// Performs the simulation of the Ising model using the Metropolis algorithm.
// It takes a lattice, its size, the number of simulation steps, and the physical parameters (interaction energy J, external magnetic field H, and temperature T).
// It executes thermalization steps to reach thermal equilibrium and then performs additional steps to calculate statistical averages of energy and magnetization. The specific heat and magnetic susceptibility are also calculated using these averages
ising_data ising_metropolis(LATYPE *lattice, u64 n, u64 steps, double J, double H, double T) {
    ising_data data = {0};

    bool valid_lattice = lattice;
    if (!valid_lattice) {
        fprintf(stderr, "Warning: lattice provided is NULL. Will be created a lattice and freed afterwards. You will have no access to that lattice");
        lattice = init_random_lattice(n);
    }

    double energy        = 0.0;
    double magnetization = 0.0;


    for (u64 i = 0; i < steps; ++i) {
        termal_step(lattice, n, J, H, T);
        energy_magnetization_pack E_M = lattice_energy_magnetization(lattice, n, J, H);
        energy        = E_M.energy;
        magnetization = E_M.magnetization;

        data.avg_E  += energy;
        data.avg_M  += magnetization;
        data.avg_E2 += energy * energy;
        data.avg_M2 += magnetization * magnetization;
    }

    double avg_fac  = 1.0 / (steps * n * n);
    double avg_fac2 = 1.0 / (steps * steps * n * n);


    data.specific_heat = (avg_fac * data.avg_E2 - avg_fac2 * data.avg_E * data.avg_E) / (T * T);
    data.magnetic_sus  = (avg_fac * data.avg_M2 - avg_fac2 * data.avg_M * data.avg_M) / T;

    data.avg_E *= avg_fac;
    data.avg_M *= avg_fac;

    if (!valid_lattice)
        free(lattice);

    return data;
}
// exports the lattice as a PNG image. It converts the spin values to pixel colors (black or white) 
// and uses the STB library to write the image file.
void export_lattice(double T, LATYPE *lattice, int n) {
    unsigned int char_size = snprintf(NULL, 0, "./imgs/%e.png", T) + 1;
    char *str = (char*)calloc(char_size, 1);
    snprintf(str, char_size, "./imgs/%e.png", T);

    uint32_t *lattice_img = (uint32_t*)calloc(n * n, sizeof(uint32_t));
    for (int i = 0; i < n * n; ++i)
        lattice_img[i] = (int)lattice[i] == 1? 0xFFFFFFFF: 0xFF000000;

    stbi_write_png(str, n, n, 4, lattice_img, sizeof(uint32_t) * n);

    free(str);
    free(lattice_img);
}

//Compilar com:
//gcc -O3 -Wall -Wextra -pedantic -o main main.c -lm

// is the entry point of the program. It sets the simulation parameters, initializes the lattice, iterates over temperatures, performs the simulation for each temperature, and writes the results to a data file.
// Optionally, it exports the lattice as PNG images.
int main(void) {
    int n_temperatures = 100;
    u64 steps          = 1000000;
    u64 eq_steps       = 200000;
    u64 n              = 32;
    double Tmin        = 1.0;
    double Tmax        = 4.0;
    double dT          = (Tmax - Tmin) / n_temperatures;
    double trans_start = 2.1; //inicio da faixa de transicão
    double trans_end   = 2.8; //final da faixa de transicão
    double fac_trans   = 30.0;//fator de numero de pontos na região de transicão (esse valor x quantos teriam). Caso não queira, so colocar =1.0
    
    double J = 1.0;
    double H = 0.0;

    bool print_lattice = true; //desativar torna o programa bem mais rapido

    FILE *out_test = fopen("./data.dat", "w");

    double T = Tmin;
    int i = 0;
    double dT_var = dT;
    LATYPE *lattice = init_random_lattice(n);
    while (T <= Tmax) {
        srand(time(NULL));
        //randomize_lattice(lattice, n);

        if (i % (n_temperatures / 10) == 0)
            printf("T = %.5e\n", T);
        fflush(out_test);
        fflush(stdout);


        ising_data data = ising_metropolis(lattice, n, eq_steps, J, H, T); //for EQ
        data = ising_metropolis(lattice, n, steps, J, H, T);

        if (print_lattice)
            export_lattice(T, lattice, n);

        fprintf(out_test, "%e\t%e\t%e\t%e\t%e\n", T, data.avg_E, data.avg_M, data.specific_heat, data.magnetic_sus);
        ++i;

        dT_var = dT;
        if (trans_start < T && T < trans_end)
            dT_var = dT / fac_trans;

        T += dT_var;
    }
    free(lattice);
    fclose(out_test);
    return 0;
}
