#include "blake3.h"
#include "test.h"
#include "bench.h"
#include "common.h"
#include "sample_z_small.h"
#include "sample_z_large.h"

// Define constants used in the program
#define R       (HEIGHT+1)  // Define R as HEIGHT incremented by 1
#define V       (HEIGHT+3)  // Define V as HEIGHT incremented by 3
#define TAU     1000        // Define TAU as 1000
#define NTI     130         // Define NTI as 130

/**
 * @brief Computes a hash of the given matrices using BLAKE3.
 *
 * @param[out] h           - The resulting hash of size BLAKE3_OUT_LEN.
 * @param[in] A            - The first matrix of size [R][V].
 * @param[in] t            - The second matrix of size [TAU][V].
 * @param[in] W            - The third matrix of size [R][NTI].
 */
static void pibnd_hash(uint8_t h[BLAKE3_OUT_LEN], params::poly_q A[R][V],
                       params::poly_q t[TAU][V], params::poly_q W[R][NTI]) {

    blake3_hasher hasher;  // Declare a BLAKE3 hasher object

    // Initialize the BLAKE3 hasher
    blake3_hasher_init(&hasher);

    // Hash the public key (A matrix)
    // Iterate through each element of the A matrix and update the hasher
    for (size_t i = 0; i < R; i++) {
        for (int j = 0; j < V; j++) {
            blake3_hasher_update(&hasher, (const uint8_t *)A[i][j].data(), 16 * DEGREE);
        }
    }

    // Hash the t matrix
    // Iterate through each element of the t matrix and update the hasher
    for (size_t i = 0; i < TAU; i++) {
        for (int j = 0; j < V; j++) {
            blake3_hasher_update(&hasher, (const uint8_t *)t[i][j].data(), 16 * DEGREE);
        }
    }

    // Hash the W matrix
    // Iterate through each element of the W matrix and update the hasher
    for (size_t i = 0; i < R; i++) {
        for (int j = 0; j < NTI; j++) {
            blake3_hasher_update(&hasher, (const uint8_t *)W[i][j].data(), 16 * DEGREE);
        }
    }

    // Finalize the hash computation and store the result in h
    blake3_hasher_finalize(&hasher, h, BLAKE3_OUT_LEN);
}


/**
 * @brief Performs rejection sampling for the given matrices.
 *
 * @param[in] Z    - The first matrix of size [V][NTI].
 * @param[in] SC   - The second matrix of size [V][NTI].
 * @param[in] s2   - A scalar value.
 * @return int     - Returns 1 if rejection sampling is successful, 0 otherwise.
 */
static int pibnd_rej_sampling(params::poly_q Z[V][NTI], params::poly_q SC[V][NTI], uint64_t s2) {
    // Declare arrays to store coefficients of polynomials
    array<mpz_t, params::poly_q::degree> coeffs0, coeffs1;

    // Declare a polynomial variable and multiple large integer (mpz_t) variables
    params::poly_q t;
    mpz_t dot, norm, qDivBy2, tmp;

    // Declare variables for calculations
    double r, M = 3.0;
    int64_t seed;
    mpf_t u;  // Floating-point number with arbitrary precision
    uint8_t buf[8];  // Buffer to store random bytes
    gmp_randstate_t state;  // Random state for GMP library
    int result;

    // Initialize the arbitrary precision floating-point number
    mpf_init(u);

    // Initialize the random state with Mersenne Twister algorithm
    gmp_randinit_mt(state);

    // Initialize the mpz_t variables
    mpz_inits(dot, norm, qDivBy2, tmp, nullptr);

    // Initialize the coefficient arrays with appropriate size
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_init2(coeffs0[i], (params::poly_q::bits_in_moduli_product() << 2));
        mpz_init2(coeffs1[i], (params::poly_q::bits_in_moduli_product() << 2));
    }

    // Compute qDivBy2 as half of the moduli product
    mpz_fdiv_q_2exp(qDivBy2, params::poly_q::moduli_product(), 1);

    // Initialize norm and dot to zero
    mpz_set_ui(norm, 0);
    mpz_set_ui(dot, 0);

    // Iterate through matrices to compute dot and norm
    for (int i = 0; i < V - 1; i++) {
        for (int j = 0; j < NTI; j++) {
            t = Z[i][j];
            t.invntt_pow_invphi();
            t.poly2mpz(coeffs0);
            t = SC[i][j];
            t.invntt_pow_invphi();
            t.poly2mpz(coeffs1);
            for (size_t i = 0; i < params::poly_q::degree; i++) {
                util::center(coeffs0[i], coeffs0[i], params::poly_q::moduli_product(), qDivBy2);
                util::center(coeffs1[i], coeffs1[i], params::poly_q::moduli_product(), qDivBy2);
                mpz_mul(tmp, coeffs0[i], coeffs1[i]);
                mpz_add(dot, dot, tmp);
                mpz_mul(tmp, coeffs1[i], coeffs1[i]);
                mpz_add(norm, norm, tmp);
            }
        }
    }

    // Get random bytes into buf
    getrandom(buf, sizeof(buf), 0);

    // Copy the random bytes into seed
    memcpy(&seed, buf, sizeof(buf));

    // Seed the GMP random state with the obtained seed
    gmp_randseed_ui(state, seed);

    // Generate a random floating-point number
    mpf_urandomb(u, state, mpf_get_default_prec());

    // Check if dot is negative and update result
    result = mpz_get_d(dot) < 0;

    // Compute r based on dot and norm
    r = -2.0 * mpz_get_d(dot) + mpz_get_d(norm);
    r = r / (2.0 * s2);
    r = exp(r) / M;

    // Update result based on the comparison of u and 1/M
    result |= mpf_get_d(u) > (1/M);

    // Clear the memory of the arbitrary precision variables
    mpf_clear(u);
    mpz_clears(dot, norm, qDivBy2, tmp, nullptr);

    // Clear the memory of the coefficient arrays
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_clear(coeffs0[i]);
        mpz_clear(coeffs1[i]);
    }

    // Return the result
    return result;
}


/**
 * Test if the l2-norm is within bounds (4 * sigma * sqrt(N)).
 *
 * @param[in] r 			- the polynomial to compute the l2-norm.
 * @return the computed norm.
 */
/**
 * Test if the l2-norm of a polynomial is within bounds (4 * sigma * sqrt(N)).
 *
 * @param[in] r         - the polynomial to compute the l2-norm.
 * @param[in] sigma_sqr - the squared sigma value.
 * @return true if the norm is within bounds, false otherwise.
 */
bool pibnd_test_norm1(params::poly_q r, uint64_t sigma_sqr) {
    // Declare an array to store coefficients of the polynomial.
    array < mpz_t, params::poly_q::degree > coeffs;

    // Declare variables for the norm, half of the modulus, and a temporary variable.
    mpz_t norm, qDivBy2, tmp;

    // Unused polynomial variable.
    params::poly_q t;

    // Initialize the GMP variables.
    mpz_inits(norm, qDivBy2, tmp, nullptr);

    // Initialize the coefficients array with appropriate bit sizes.
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], (params::poly_q::bits_in_moduli_product() << 2));
    }

    // Convert the polynomial to its coefficient representation.
    r.poly2mpz(coeffs);

    // Compute qDivBy2 as half of the modulus.
    mpz_fdiv_q_2exp(qDivBy2, params::poly_q::moduli_product(), 1);

    // Initialize the norm to 0.
    mpz_set_ui(norm, 0);

    // Compute the squared l2-norm of the polynomial.
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        // Center the coefficient around 0.
        util::center(coeffs[i], coeffs[i], params::poly_q::moduli_product(), qDivBy2);

        // Square the coefficient.
        mpz_mul(tmp, coeffs[i], coeffs[i]);

        // Add the squared coefficient to the norm.
        mpz_add(norm, norm, tmp);
    }

    // Compute the bound as (sigma * sqrt(2N))^2 = sigma^2 * 2 * N.
    uint64_t bound = 2 * sigma_sqr * params::poly_q::degree;

    // Check if the computed norm is less than or equal to the bound.
    int result = mpz_cmp_ui(norm, bound) <= 0;

    // Clear the GMP variables to free memory.
    mpz_clears(norm, qDivBy2, tmp, nullptr);
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_clear(coeffs[i]);
    }

    // Return the result of the comparison.
    return result;
}


/**
 * Test if the l2-norm of a polynomial is within a specific bound defined by BOUND_B.
 *
 * @param[in] r - the polynomial to compute the l2-norm.
 * @return true if the norm is outside the bounds, false otherwise.
 */
bool pibnd_test_norm2(params::poly_q r) {
    // Declare an array to store coefficients of the polynomial.
    array < mpz_t, params::poly_q::degree > coeffs;

    // Declare variables for the norm, half of the modulus, a temporary variable, and the modulus q.
    mpz_t norm, qDivBy2, tmp, q;

    // Unused polynomial variable.
    params::poly_q t;

    // Initialize the GMP variables.
    mpz_inits(norm, qDivBy2, tmp, q, nullptr);

    // Initialize the coefficients array with appropriate bit sizes.
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], (params::poly_q::bits_in_moduli_product() << 2));
    }

    // Convert the polynomial to its coefficient representation.
    r.poly2mpz(coeffs);

    // Compute qDivBy2 as half of the modulus.
    mpz_fdiv_q_2exp(qDivBy2, params::poly_q::moduli_product(), 1);

    // Initialize the norm to 0.
    mpz_set_ui(norm, 0);

    // Set the modulus q from the predefined PRIMEQ constant.
    mpz_set_str(q, PRIMEQ, 10);

    // Compute the squared l2-norm of the polynomial.
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        // Reduce the coefficient modulo q.
        mpz_mod(coeffs[i], coeffs[i], q);

        // Center the coefficient around 0.
        util::center(coeffs[i], coeffs[i], params::poly_q::moduli_product(), qDivBy2);

        // Square the coefficient.
        mpz_mul(tmp, coeffs[i], coeffs[i]);

        // Add the squared coefficient to the norm.
        mpz_add(norm, norm, tmp);
    }

    // Set the bound from the predefined BOUND_B constant and square it.
    mpz_set_str(tmp, BOUND_B, 10);
    mpz_mul(tmp, tmp, tmp);

    // Check if the computed norm is less than or equal to the squared bound.
    int result = mpz_cmp(norm, tmp) <= 0;

    // Clear the GMP variables to free memory.
    mpz_clears(norm, qDivBy2, tmp, q, nullptr);
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_clear(coeffs[i]);
    }

    // Return the negation of the result (true if outside the bounds, false if within).
    return !result;
}


// Sample a challenge.
/**
 * Sample a challenge polynomial.
 *
 * @param[out] f - The polynomial to be sampled.
 */
void pibnd_sample_chall(params::poly_q & f) {
    // Sample a polynomial from the zero distribution.
    f = nfl::ZO_dist();

    // Transform the polynomial to NTT domain using the power of phi.
    f.ntt_pow_phi();
}


/**
 * The prover function for the BND protocol.
 *
 * @param[out] h - The hash output.
 * @param[out] Z - The polynomial matrix.
 * @param[in] A - The public key polynomial matrix.
 * @param[in] t - The polynomial matrix.
 * @param[in] s - The secret polynomial matrix.
 */
static void pibnd_prover(uint8_t h[BLAKE3_OUT_LEN], params::poly_q Z[V][NTI],
                         params::poly_q A[R][V], params::poly_q t[TAU][V],
                         params::poly_q s[TAU][V]) {

    // Declare polynomial matrices W, C, and SC.
    params::poly_q W[R][NTI], C[TAU][NTI], SC[V][NTI];

    // Declare an array to store coefficients.
    std::array < mpz_t, params::poly_q::degree > coeffs;

    // Declare variables for mathematical operations.
    mpz_t qDivBy2;
    int64_t coeff;

    // Initialize qDivBy2 and compute half of the moduli product.
    mpz_init(qDivBy2);
    mpz_fdiv_q_2exp(qDivBy2, params::poly_q::moduli_product(), 1);

    // Initialize the coefficients array.
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], (params::poly_q::bits_in_moduli_product() << 2));
    }

    // Main loop for the prover.
    do {
        // Sample polynomial matrix Y from Gaussian distribution.
        for (int i = 0; i < V; i++) {
            for (int j = 0; j < NTI; j++) {
                for (size_t k = 0; k < params::poly_q::degree; k++) {
                    // Determine which sigma to use based on the index.
                    if (i < V - 1) {
                        coeff = sample_z(0.0, SIGMA_B1);
                    } else {
                        coeff = sample_z((__float128) 0.0, (__float128) SIGMA_B2);
                    }
                    mpz_set_si(coeffs[k], coeff);
                }
                // Convert coefficients to polynomial and transform to NTT domain.
                Z[i][j].mpz2poly(coeffs);
                Z[i][j].ntt_pow_phi();
            }
        }

        // Compute W as the product of A and Y.
        for (int i = 0; i < R; i++) {
            for (int j = 0; j < NTI; j++) {
                W[i][j] = 0;
                for (int k = 0; k < V; k++) {
                    W[i][j] = W[i][j] + A[i][k] * Z[k][j];
                }
            }
        }

        // Compute the hash of A, t, and W.
        pibnd_hash(h, A, t, W);

        // Seed the random number generator with the hash.
        nfl::fastrandombytes_seed(h);

        // Sample the challenge matrix C.
        for (int i = 0; i < TAU; i++) {
            for (int j = 0; j < NTI; j++) {
                C[i][j] = nfl::ZO_dist();
                C[i][j].ntt_pow_phi();
            }
        }

        // Reseed the random number generator.
        nfl::fastrandombytes_reseed();

        // Compute Z as the sum of Y and SC.
        for (int i = 0; i < V; i++) {
            for (int j = 0; j < NTI; j++) {
                SC[i][j] = 0;
                for (int k = 0; k < TAU; k++) {
                    SC[i][j] = SC[i][j] + s[k][i] * C[k][j];
                }
                Z[i][j] = Z[i][j] + SC[i][j];
            }
        }
        // Continue until the rejection sampling condition is met.
    } while (pibnd_rej_sampling(Z, SC, SIGMA_B1 * SIGMA_B1) == 1);

    // Clear the allocated memory for coefficients and qDivBy2.
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_clear(coeffs[i]);
    }
    mpz_clear(qDivBy2);
}


/**
 * The verifier function for the BND protocol.
 *
 * @param[in] h1 - The hash input.
 * @param[in] Z - The polynomial matrix.
 * @param[in] A - The public key polynomial matrix.
 * @param[in] t - The polynomial matrix.
 * @return int - Returns 1 if verification is successful, 0 otherwise.
 */
int pibnd_verifier(uint8_t h1[BLAKE3_OUT_LEN], params::poly_q Z[V][NTI],
                   params::poly_q A[R][V], params::poly_q t[TAU][V]) {

    // Declare polynomial matrices W and C.
    params::poly_q W[R][NTI], C[TAU][NTI];

    // Declare a hash buffer for comparison.
    uint8_t h2[BLAKE3_OUT_LEN];

    // Result variable to store the outcome of the verification.
    int result;

    // Seed the random number generator with the provided hash.
    nfl::fastrandombytes_seed(h1);

    // Sample the challenge matrix C.
    for (int i = 0; i < TAU; i++) {
        for (int j = 0; j < NTI; j++) {
            C[i][j] = nfl::ZO_dist();
            C[i][j].ntt_pow_phi();
        }
    }

    // Compute W as the product of A and Z.
    for (int i = 0; i < R; i++) {
        for (int j = 0; j < NTI; j++) {
            W[i][j] = 0;
            for (int k = 0; k < V; k++) {
                W[i][j] = W[i][j] + A[i][k] * Z[k][j];
            }
        }
    }

    // Update W by subtracting the product of t and C.
    for (int i = 0; i < R; i++) {
        for (int j = 0; j < NTI; j++) {
            for (int k = 0; k < TAU; k++) {
                W[i][j] = W[i][j] - t[k][i] * C[k][j];
            }
        }
    }

    // Compute the hash of A, t, and W.
    pibnd_hash(h2, A, t, W);

    // Compare the computed hash with the provided hash.
    result = memcmp(h1, h2, BLAKE3_OUT_LEN) == 0;

    // Check the norms of the polynomials in Z.
    for (int i = 0; i < V; i++) {
        for (int j = 0; j < NTI; j++) {
            Z[i][j].invntt_pow_invphi();
            if (i < V - 1) {
                // Check the norm using the first test function.
                result &= pibnd_test_norm1(Z[i][j], SIGMA_B1 * SIGMA_B1);
            } else {
                // Check the norm using the second test function.
                result &= pibnd_test_norm2(Z[i][j]);
            }
        }
    }
    return result;
}


/**
 * Test function for the BND protocol.
 * This function tests the consistency of the BND proof by creating instances of the protocol,
 * generating a proof using the prover function, and then verifying the proof using the verifier function.
 */
static void test() {
    // Declare polynomial matrices and other necessary variables.
    params::poly_q A[R][V], s[TAU][V], t[TAU][V], Z[V][NTI];
    uint8_t h1[BLAKE3_OUT_LEN];
    std::array < mpz_t, params::poly_q::degree > coeffs;
    gmp_randstate_t prng;
    mpz_t q;

    // Initialize the prime number q and the coefficients array.
    mpz_init(q);
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], (params::poly_q::bits_in_moduli_product() << 2));
    }

    // Initialize the random number generator.
    gmp_randinit_default(prng);
    mpz_set_str(q, PRIMEQ, 10);

    // Create instances of the polynomial matrix A with random coefficients.
    for (int i = 0; i < R; i++) {
        for (int j = 0; j < V; j++) {
            for (size_t k = 0; k < params::poly_q::degree; k++) {
                mpz_urandomb(coeffs[k], prng, LEVEL);
                mpz_mod(coeffs[k], coeffs[k], q);
            }
            A[i][j].mpz2poly(coeffs);
            A[i][j].ntt_pow_phi();
        }
    }

    // Create TAU relations where each t_i is computed as the product of A and s_i.
    for (int i = 0; i < TAU; i++) {
        for (int j = 0; j < V; j++) {
            pibnd_sample_chall(s[i][j]);
        }
        for (int j = 0; j < R; j++) {
            t[i][j] = 0;
            for (int k = 0; k < V; k++) {
                t[i][j] = t[i][j] + A[j][k] * s[i][k];
            }
        }
    }

    // Test the consistency of the BND proof.
    TEST_ONCE("BND proof is consistent") {
        // Generate a proof using the prover function.
        pibnd_prover(h1, Z, A, t, s);
        // Verify the proof using the verifier function and assert the result.
        TEST_ASSERT(pibnd_verifier(h1, Z, A, t) == 1, end);
    } TEST_END;

    end:

    // Clear memory and cleanup.
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_clear(coeffs[i]);
    }
    mpz_clear(q);
    gmp_randclear(prng);
    return;
}


/**
 * Benchmark function for the BND protocol.
 * This function benchmarks the performance of the BND prover and verifier functions.
 * It initializes instances of the protocol and then measures the time taken by the prover and verifier functions.
 */
static void bench() {
    // Declare polynomial matrices and other necessary variables.
    params::poly_q A[R][V], s[TAU][V], t[TAU][V], Z[V][NTI];
    uint8_t h1[BLAKE3_OUT_LEN];
    std::array < mpz_t, params::poly_q::degree > coeffs;
    gmp_randstate_t prng;
    mpz_t q;

    // Initialize the prime number q and the coefficients array.
    mpz_init(q);
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], (params::poly_q::bits_in_moduli_product() << 2));
    }

    // Initialize the random number generator.
    gmp_randinit_default(prng);
    mpz_set_str(q, PRIMEQ, 10);

    // Create instances of the polynomial matrix A with random coefficients.
    for (int i = 0; i < R; i++) {
        for (int j = 0; j < V; j++) {
            for (size_t k = 0; k < params::poly_q::degree; k++) {
                mpz_urandomb(coeffs[k], prng, LEVEL);
                mpz_mod(coeffs[k], coeffs[k], q);
            }
            A[i][j].mpz2poly(coeffs);
            A[i][j].ntt_pow_phi();
        }
    }

    // Create TAU relations where each t_i is computed as the product of A and s_i.
    for (int i = 0; i < TAU; i++) {
        for (int j = 0; j < V; j++) {
            pibnd_sample_chall(s[i][j]);
        }
        for (int j = 0; j < R; j++) {
            t[i][j] = 0;
            for (int k = 0; k < V; k++) {
                t[i][j] = t[i][j] + A[j][k] * s[i][k];
            }
        }
    }

    // Benchmark the prover function.
    BENCH_SMALL("BND prover (N relations)", pibnd_prover(h1, Z, A, t, s));
    // Benchmark the verifier function.
    BENCH_SMALL("BND verifier (N relations)", pibnd_verifier(h1, Z, A, t));

    // Clear memory and cleanup.
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_clear(coeffs[i]);
    }
    mpz_clear(q);
    gmp_randclear(prng);
    return;
}

/**
 * Main function to run tests and benchmarks for the lattice-based BND proof.
 * @param argc: Number of command-line arguments.
 * @param argv: Array of command-line arguments.
 * @return: Returns 0 on successful execution.
 */
int main(int argc, char *argv[]) {
    // Print the header for the test section.
    printf("\n** Tests for lattice-based BND proof:\n\n");
    // Run the test function.
    test();

    // Print the header for the benchmark section.
    printf("\n** Benchmarks for lattice-based BND proof:\n\n");
    // Run the benchmark function.
    bench();

    return 0;
}

