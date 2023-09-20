#include <math.h>
#include <stdlib.h>
#include <iostream>

#include <flint/flint.h>
#include <flint/fmpz_mod_poly.h>

#include "blake3.h"
#include "common.h"
#include "test.h"
#include "bench.h"
#include "assert.h"
#include "sample_z_small.h"

/*============================================================================*/
/* Private definitions                                                        */
/*============================================================================*/

// Define the number of messages
#define MSGS        1000

/**
 * @brief Computes a linear hash using the BLAKE3 hashing algorithm.
 *
 * @param beta Output parameter to store the computed hash.
 * @param key Commitment key used in the hashing process.
 * @param x First commitment value.
 * @param y Second commitment value.
 * @param alpha Array of two polynomial values from the linear relation.
 * @param u Polynomial value used in the hashing process.
 * @param t Polynomial value used in the hashing process.
 * @param _t Another polynomial value used in the hashing process.
 */
static void lin_hash(params::poly_q & beta, comkey_t & key, commit_t x,
                     commit_t y, params::poly_q alpha[2], params::poly_q & u,
                     params::poly_q t, params::poly_q _t) {
    // Define a hash array with a length equal to the output length of BLAKE3
    uint8_t hash[BLAKE3_OUT_LEN];
    // Initialize a BLAKE3 hasher object
    blake3_hasher hasher;

    // Start the BLAKE3 hashing process
    blake3_hasher_init(&hasher);

    /* Hash public key. */
    // Iterate through the HEIGHT of the key and hash the A1 component
    for (size_t i = 0; i < HEIGHT; i++) {
        for (int j = 0; j < WIDTH - HEIGHT; j++) {
            blake3_hasher_update(&hasher, (const uint8_t *)key.A1[i][j].data(),
                                 16 * DEGREE);
        }
    }
    // Hash the A2 component of the key
    for (size_t j = 0; j < WIDTH; j++) {
        blake3_hasher_update(&hasher, (const uint8_t *)key.A2[0][j].data(), 16 * DEGREE);
    }

    /* Hash alpha, beta from linear relation. */
    // Hash the two polynomial values from the linear relation
    for (size_t i = 0; i < 2; i++) {
        blake3_hasher_update(&hasher, (const uint8_t *)alpha[i].data(), 16 * DEGREE);
    }

    // Hash the c1 component of the first and second commitment values
    blake3_hasher_update(&hasher, (const uint8_t *)x.c1.data(), 16 * DEGREE);
    blake3_hasher_update(&hasher, (const uint8_t *)y.c1.data(), 16 * DEGREE);
    // Hash the c2 component of the first and second commitment values
    for (size_t i = 0; i < x.c2.size(); i++) {
        blake3_hasher_update(&hasher, (const uint8_t *)x.c2[i].data(), 16 * DEGREE);
        blake3_hasher_update(&hasher, (const uint8_t *)y.c2[i].data(), 16 * DEGREE);
    }

    // Hash the polynomial values u, t, and _t
    blake3_hasher_update(&hasher, (const uint8_t *)u.data(), 16 * DEGREE);
    blake3_hasher_update(&hasher, (const uint8_t *)t.data(), 16 * DEGREE);
    blake3_hasher_update(&hasher, (const uint8_t *)_t.data(), 16 * DEGREE);

    // Finalize the hashing process and store the result in the hash array
    blake3_hasher_finalize(&hasher, hash, BLAKE3_OUT_LEN);

    /* Sample challenge from RNG seeded with hash. */
    // Seed the fast random bytes generator with the computed hash
    nfl::fastrandombytes_seed(hash);
    // Sample a challenge using the bdlop_sample_chal function
    bdlop_sample_chal(beta);
    // Reseed the fast random bytes generator
    nfl::fastrandombytes_reseed();
}


/**
 * @brief Computes the multiplicative inverse of a polynomial in a given field.
 *
 * @param inv Output parameter to store the computed inverse polynomial.
 * @param p Input polynomial for which the inverse is to be computed.
 */
static void poly_inverse(params::poly_q & inv, params::poly_q p) {
    // Declare an array to store coefficients of the polynomial
    std::array < mpz_t, params::poly_q::degree > coeffs;
    // Declare a variable to store the modulus of the field
    fmpz_t q;
    // Declare variables to store the polynomial and the irreducible polynomial
    fmpz_mod_poly_t poly, irred;
    // Declare a context variable for modular arithmetic operations
    fmpz_mod_ctx_t ctx_q;

    // Initialize the modulus variable
    fmpz_init(q);
    // Initialize the coefficients array with the appropriate bit size
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], (params::poly_q::bits_in_moduli_product() << 2));
    }

    // Set the modulus value to the product of moduli for the polynomial
    fmpz_set_mpz(q, params::poly_q::moduli_product());
    // Initialize the context for modular arithmetic with the modulus
    fmpz_mod_ctx_init(ctx_q, q);
    // Initialize the polynomial and irreducible polynomial variables in the context
    fmpz_mod_poly_init(poly, ctx_q);
    fmpz_mod_poly_init(irred, ctx_q);

    // Convert the polynomial to its coefficient representation
    p.poly2mpz(coeffs);
    // Define the irreducible polynomial for the field
    fmpz_mod_poly_set_coeff_ui(irred, params::poly_q::degree, 1, ctx_q);
    fmpz_mod_poly_set_coeff_ui(irred, 0, 1, ctx_q);

    // Set the polynomial coefficients from the array
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        fmpz_mod_poly_set_coeff_mpz(poly, i, coeffs[i], ctx_q);
    }
    // Compute the multiplicative inverse of the polynomial modulo the irreducible polynomial
    fmpz_mod_poly_invmod(poly, poly, irred, ctx_q);

    // Retrieve the coefficients of the inverse polynomial
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        fmpz_mod_poly_get_coeff_mpz(coeffs[i], poly, i, ctx_q);
    }

    // Convert the coefficient representation back to the polynomial form
    inv.mpz2poly(coeffs);

    // Clear the memory allocated for the modulus
    fmpz_clear(q);
    // Clear the memory allocated for the coefficients array
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_clear(coeffs[i]);
    }
}


/**
 * @brief Computes the multiplicative inverses of a set of polynomials simultaneously.
 *
 * This function computes the multiplicative inverses of a set of polynomials using
 * a simultaneous inversion algorithm. The inverses are computed in the NTT domain.
 *
 * @param inv Output array to store the computed inverse polynomials.
 * @param m Input array of polynomials for which the inverses are to be computed.
 */
static void simul_inverse(params::poly_q inv[MSGS], params::poly_q m[MSGS]) {
    // Declare a polynomial 'u' to store intermediate results
    params::poly_q u;
    // Declare an array 't' to store intermediate polynomials
    params::poly_q t[MSGS];

    // Initialize the first polynomial in the 'inv' and 't' arrays
    inv[0] = m[0];
    t[0] = m[0];

    // Compute the product of the polynomials up to the i-th polynomial
    for (size_t i = 1; i < MSGS; i++) {
        t[i] = m[i];
        inv[i] = inv[i - 1] * m[i];
    }

    // Set 'u' to the product of all the polynomials
    u = inv[MSGS - 1];
    // Convert 'u' from NTT domain to coefficient domain
    u.invntt_pow_invphi();
    // Compute the multiplicative inverse of 'u'
    poly_inverse(u, u);
    // Convert 'u' back to NTT domain
    u.ntt_pow_phi();

    // Compute the inverses for each polynomial using the simultaneous inversion algorithm
    for (size_t i = MSGS - 1; i > 0; i--) {
        inv[i] = u * inv[i - 1];
        u = u * t[i];
    }
    // Set the first inverse polynomial
    inv[0] = u;
}


/**
 * @brief Performs rejection sampling based on the given polynomials.
 *
 * This function performs rejection sampling using the provided polynomials `z` and `v`.
 * It computes the dot product and norm of the polynomials and then uses these values
 * to decide whether to accept or reject the sample.
 *
 * @param z Input array of polynomials.
 * @param v Input array of polynomials.
 * @param s2 A scalar value used in the rejection sampling computation.
 * @return Returns 1 if the sample is rejected, 0 otherwise.
 */
static int rej_sampling(params::poly_q z[WIDTH], params::poly_q v[WIDTH], uint64_t s2) {
    // Declare arrays to store coefficients of the polynomials in coefficient domain
    array<mpz_t, params::poly_q::degree> coeffs0, coeffs1;
    // Declare a polynomial for intermediate computations
    params::poly_q t;
    // Declare variables to store the dot product, norm, half of the modulus, and a temporary value
    mpz_t dot, norm, qDivBy2, tmp;
    // Declare variables for the rejection sampling computation
    double r, M = 1.75;
    // Seed for random number generation
    int64_t seed;
    // Floating-point variable for random number generation
    mpf_t u;
    // Buffer to store random bytes
    uint8_t buf[8];
    // Random state for GMP's random number generation
    gmp_randstate_t state;
    // Result of the rejection sampling (1 if rejected, 0 otherwise)
    int result;

    // Initialize the floating-point variable
    mpf_init(u);
    // Initialize GMP's random number generator with the Mersenne Twister algorithm
    gmp_randinit_mt(state);
    // Initialize the mpz variables
    mpz_inits(dot, norm, qDivBy2, tmp, nullptr);
    // Initialize the coefficient arrays with the appropriate bit size
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_init2(coeffs0[i], (params::poly_q::bits_in_moduli_product() << 2));
        mpz_init2(coeffs1[i], (params::poly_q::bits_in_moduli_product() << 2));
    }

    // Compute half of the modulus
    mpz_fdiv_q_2exp(qDivBy2, params::poly_q::moduli_product(), 1);
    // Initialize the dot product and norm to zero
    mpz_set_ui(norm, 0);
    mpz_set_ui(dot, 0);
    // Compute the dot product and norm using the provided polynomials
    for (int i = 0; i < WIDTH; i++) {
        t = z[i];
        t.invntt_pow_invphi();
        t.poly2mpz(coeffs0);
        t = v[i];
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

    // Generate a random seed using the system's random number generator
    ssize_t grret = getrandom(buf, sizeof(buf), 0);
    if (grret == -1) {
        throw "Could not generate a random seed. Check if something is wrong with system prng.";
    }
    // Convert the random bytes to an integer seed
    memcpy(&seed, buf, sizeof(buf));
    // Seed GMP's random number generator with the generated seed
    gmp_randseed_ui(state, seed);
    // Generate a random floating-point number
    mpf_urandomb(u, state, mpf_get_default_prec());

    // Compute the rejection threshold
    r = -2.0 * mpz_get_d(dot) + mpz_get_d(norm);
    r = r / (2.0 * s2);
    r = exp(r) / M;
    // Decide whether to accept or reject the sample
    result = mpf_get_d(u) > r;

    // Clear the allocated memory for the variables
    mpf_clear(u);
    mpz_clears(dot, norm, qDivBy2, tmp, nullptr);
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_clear(coeffs0[i]);
        mpz_clear(coeffs1[i]);
    }

    // Return the result of the rejection sampling
    return result;
}


/**
 * @brief Prover function for the linear proof system.
 *
 * This function implements the prover's algorithm for a linear proof system.
 * The prover samples polynomials, computes commitments, and then uses rejection sampling
 * to ensure that the samples are valid.
 *
 * @param y Output array of polynomials.
 * @param _y Output array of polynomials.
 * @param t Polynomial result of the prover's computations.
 * @param _t Polynomial result of the prover's computations.
 * @param u Polynomial result of the prover's computations.
 * @param x Commitment from the verifier.
 * @param _x Commitment from the verifier.
 * @param alpha Array of polynomials from the linear relation.
 * @param key Commitment key.
 * @param r Random polynomials.
 * @param _r Random polynomials.
 */
static void lin_prover(params::poly_q y[WIDTH], params::poly_q _y[WIDTH],
                       params::poly_q& t, params::poly_q& _t, params::poly_q& u,
                       commit_t x, commit_t _x, params::poly_q alpha[2],
                       comkey_t & key, vector<params::poly_q> r, vector<params::poly_q> _r) {

    // Declare variables for the challenge, temporary polynomials, and coefficient arrays
    params::poly_q beta, tmp[WIDTH], _tmp[WIDTH];
    // TODO: Array below is old
    // array<mpz_t, params::poly_q::degree> coeffs;
    array<mpz_t, params::poly_q::degree> coeffs1;
    array<mpz_t, params::poly_q::degree> coeffs2;

    // Declare a variable to store half of the modulus
    mpz_t qDivBy2;
    // Variables to store the results of the rejection sampling
    int rej0, rej1;

    // Initialize the variable to store half of the modulus
    mpz_init(qDivBy2);
    // Initialize the coefficient arrays with the appropriate bit size
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_init2(coeffs1[i], (params::poly_q::bits_in_moduli_product() << 2));
        mpz_init2(coeffs2[i], (params::poly_q::bits_in_moduli_product() << 2));
    }
    // Compute half of the modulus
    mpz_fdiv_q_2exp(qDivBy2, params::poly_q::moduli_product(), 1);

    do {
        // Prover samples y and y' from a Gaussian distribution
        for (int i = 0; i < WIDTH; i++) {
            for (size_t k = 0; k < params::poly_q::degree; k++) {
                int coeff1 = sample_z(0.0, SIGMA_C);
                int coeff2 = sample_z(0.0, SIGMA_C);
                mpz_set_si(coeffs1[k], coeff1);
                mpz_set_si(coeffs2[k], coeff2);
            }
            y[i].mpz2poly(coeffs1);
            y[i].ntt_pow_phi();
            _y[i].mpz2poly(coeffs2);
            _y[i].ntt_pow_phi();
            /* TODO: This is the old one
             * for (size_t k = 0; k < params::poly_q::degree; k++) {
                int64_t coeff = sample_z(0.0, SIGMA_C);
                mpz_set_si(coeffs[k], coeff);


            }
            y[i].mpz2poly(coeffs);
            y[i].ntt_pow_phi();
            for (size_t k = 0; k < params::poly_q::degree; k++) {
                int64_t coeff = sample_z(0.0, SIGMA_C);
                mpz_set_si(coeffs[k], coeff);
            }
            _y[i].mpz2poly(coeffs);
            _y[i].ntt_pow_phi();*/
        }

        // Compute t and _t using the sampled y and y' and the commitment key
        t = y[0];
        _t = _y[0];
        for (int i = 0; i < HEIGHT; i++) {
            for (int j = 0; j < WIDTH - HEIGHT; j++) {
                t = t + key.A1[i][j] * y[j + HEIGHT];
                _t = _t + key.A1[i][j] * _y[j + HEIGHT];
            }
        }

        // Compute u using the sampled y and y', alpha, and the commitment key
        u = 0;
        for (int i = 0; i < WIDTH; i++) {
            u = u + alpha[0] * (key.A2[0][i] * y[i]) - (key.A2[0][i] * _y[i]);
        }

        // Sample the challenge using the linear hash function
        lin_hash(beta, key, x, _x, alpha, u, t, _t);

        // Update y and _y using the challenge and random polynomials
        for (int i = 0; i < WIDTH; i++) {
            tmp[i] = beta * r[i];
            _tmp[i] = beta * _r[i];
            y[i] = y[i] + tmp[i];
            _y[i] = _y[i] + _tmp[i];
        }
        // Perform rejection sampling on y and _y
        rej0 = rej_sampling(y, tmp, SIGMA_C * SIGMA_C);
        rej1 = rej_sampling(_y, _tmp, SIGMA_C * SIGMA_C);
    } while (rej0 || rej1);  // Repeat until both samples are accepted

    // Clear the allocated memory for the variables
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_clear(coeffs1[i]);
        mpz_clear(coeffs2[i]);
    }
    mpz_clear(qDivBy2);
}


/**
 * @brief Verifier function for the linear proof system.
 *
 * This function implements the verifier's algorithm for a linear proof system.
 * The verifier samples a challenge, checks the norms of the provided polynomials,
 * computes commitments, and then verifies the correctness of the prover's computations.
 *
 * @param z Array of polynomials from the prover.
 * @param _z Array of polynomials from the prover.
 * @param t Polynomial result of the prover's computations.
 * @param _t Polynomial result of the prover's computations.
 * @param u Polynomial result of the prover's computations.
 * @param x Commitment from the prover.
 * @param _x Commitment from the prover.
 * @param alpha Array of polynomials from the linear relation.
 * @param key Commitment key.
 * @return Returns 1 if the verification is successful, otherwise 0.
 */
static int lin_verifier(params::poly_q z[WIDTH], params::poly_q _z[WIDTH],
                        params::poly_q t, params::poly_q _t, params::poly_q u,
                        commit_t x, commit_t _x, params::poly_q alpha[2], comkey_t & key) {

    // Declare variables for the challenge, temporary polynomials, and a zero polynomial
    params::poly_q beta, v, _v, tmp, zero = 0;
    // Initialize the result to 1 (true)
    int result = 1;

    // Sample the challenge using the linear hash function
    lin_hash(beta, key, x, _x, alpha, u, t, _t);

    // Verifier checks the norm of each polynomial in z and _z
    for (int i = 0; i < WIDTH; i++) {
        v = z[i];
        v.invntt_pow_invphi();
        result &= bdlop_test_norm(v, 4 * SIGMA_C * SIGMA_C);
        v = _z[i];
        v.invntt_pow_invphi();
        result &= bdlop_test_norm(v, 4 * SIGMA_C * SIGMA_C);
    }

    // Verifier computes A1z and A1z' using the commitment key
    v = z[0];
    _v = _z[0];
    for (int i = 0; i < HEIGHT; i++) {
        for (int j = 0; j < WIDTH - HEIGHT; j++) {
            v = v + key.A1[i][j] * z[j + HEIGHT];
            _v = _v + key.A1[i][j] * _z[j + HEIGHT];
        }
    }

    // Verifier checks the correctness of t and _t
    tmp = t + beta * x.c1 - v;
    tmp.invntt_pow_invphi();
    result &= (tmp == zero);
    tmp = _t + beta * _x.c1 - _v;
    tmp.invntt_pow_invphi();
    result &= (tmp == zero);

    // Verifier computes the polynomial v using the commitment key, z, _z, and alpha
    v = 0;
    for (int i = 0; i < WIDTH; i++) {
        v = v + alpha[0] * (key.A2[0][i] * z[i]) - (key.A2[0][i] * _z[i]);
    }
    t = (alpha[0] * x.c2[0] + alpha[1] - _x.c2[0]) * beta + u;

    // Convert t and v back to the coefficient representation
    t.invntt_pow_invphi();
    v.invntt_pow_invphi();

    // Check the correctness of t and v
    result &= ((t - v) == 0);
    // Return the result of the verification
    return result;
}

/**
 * @brief Computes a hash value based on commitments, polynomials, and other parameters.
 *
 * This function uses the BLAKE3 cryptographic hash function to compute a hash value
 * based on the provided input parameters. The hash value is then used to sample a challenge.
 *
 * @param beta Output polynomial challenge.
 * @param c Array of commitments.
 * @param d Array of commitments.
 * @param _ms Array of polynomials.
 * @param rho Array of polynomials of size SIZE.
 */
void shuffle_hash(params::poly_q & beta, commit_t c[MSGS], commit_t d[MSGS],
                  params::poly_q _ms[MSGS], params::poly_q rho[SIZE]) {

    // Declare a buffer to store the hash output
    uint8_t hash[BLAKE3_OUT_LEN];
    // Initialize a BLAKE3 hasher object
    blake3_hasher hasher;
    // Initialize the hasher
    blake3_hasher_init(&hasher);

    // (Redundant initialization, can be removed)
    blake3_hasher_init(&hasher);

    // Hash each polynomial in the _ms array and the c1 components of the c and d commitments
    for (int i = 0; i < MSGS; i++) {
        blake3_hasher_update(&hasher, (const uint8_t *)_ms[i].data(), 16 * DEGREE);
        blake3_hasher_update(&hasher, (const uint8_t *)c[i].c1.data(), 16 * DEGREE);
        blake3_hasher_update(&hasher, (const uint8_t *)d[i].c1.data(), 16 * DEGREE);

        // Hash the c2 components of the c and d commitments
        for (size_t j = 0; j < c[j].c2.size(); j++) {
            blake3_hasher_update(&hasher, (const uint8_t *)c[i].c2[j].data(), 16 * DEGREE);
            blake3_hasher_update(&hasher, (const uint8_t *)d[i].c2[j].data(), 16 * DEGREE);
        }
    }

    // Hash each polynomial in the rho array
    for (int i = 0; i < SIZE; i++) {
        blake3_hasher_update(&hasher, (const uint8_t *)rho[i].data(), 16 * DEGREE);
    }

    // Finalize the hash computation and store the result in the hash buffer
    blake3_hasher_finalize(&hasher, hash, BLAKE3_OUT_LEN);

    // Seed the random number generator with the computed hash
    nfl::fastrandombytes_seed(hash);
    // Sample a uniform polynomial challenge
    beta = nfl::uniform();
    // Reseed the random number generator for future use
    nfl::fastrandombytes_reseed();
}


/**
 * @brief Executes the shuffle proof generation.
 *
 * This function implements the prover's side of a shuffle proof protocol. It generates
 * proofs for the shuffling of polynomial commitments.
 *
 * @param y Array of polynomials.
 * @param _y Array of polynomials.
 * @param t Array of polynomials.
 * @param _t Array of polynomials.
 * @param u Array of polynomials.
 * @param d Array of commitments.
 * @param s Array of polynomials.
 * @param c Array of commitments.
 * @param ms Array of polynomials.
 * @param _ms Array of polynomials.
 * @param r Array of vectors of polynomials.
 * @param rho Array of polynomials.
 * @param key Commitment key.
 */
static void shuffle_prover(params::poly_q y[MSGS][WIDTH],
                           params::poly_q _y[MSGS][WIDTH], params::poly_q t[MSGS],
                           params::poly_q _t[MSGS], params::poly_q u[MSGS], commit_t d[MSGS],
                           params::poly_q s[MSGS], commit_t c[MSGS], params::poly_q ms[MSGS],
                           params::poly_q _ms[MSGS], vector < params::poly_q > r[MSGS],
                           params::poly_q rho[MSGS], comkey_t & key) {

    // Initialize temporary variables
    vector < params::poly_q > t0(1);
    vector < params::poly_q > _r[MSGS];
    params::poly_q alpha[2], beta, theta[MSGS], inv[MSGS];

    // Prover samples theta_i and computes commitments D_i
    for (size_t i = 0; i < MSGS - 1; i++) {
        theta[i] = nfl::ZO_dist();  // Sample from zero distribution
        theta[i].ntt_pow_phi();  // Apply NTT
        if (i == 0) {
            t0[0] = theta[0] * _ms[0];
        } else {
            t0[0] = theta[i - 1] * ms[i] + theta[i] * _ms[i];
        }
        t0[0].invntt_pow_invphi();  // Apply inverse NTT
        _r[i].resize(WIDTH);
        bdlop_sample_rand(_r[i]);  // Sample random values
        bdlop_commit(d[i], t0, key, _r[i]);  // Compute commitment
    }
    t0[0] = theta[MSGS - 2] * ms[MSGS - 1];
    t0[0].invntt_pow_invphi();
    _r[MSGS - 1].resize(WIDTH);
    bdlop_sample_rand(_r[MSGS - 1]);
    bdlop_commit(d[MSGS - 1], t0, key, _r[MSGS - 1]);

    // Compute the hash for the shuffle
    shuffle_hash(beta, c, d, _ms, rho);

    // Check relationship and compute inverse
    simul_inverse(inv, _ms);
    for (size_t i = 0; i < MSGS - 1; i++) {
        if (i == 0) {
            s[0] = theta[0] * _ms[0] - beta * ms[0];
        } else {
            s[i] = theta[i - 1] * ms[i] + theta[i] * _ms[i] - s[i - 1] * ms[i];
        }
        s[i] = s[i] * inv[i];
    }

    // Run multiple linear proof instances for each commitment
    for (size_t l = 0; l < MSGS; l++) {
        if (l < MSGS - 1) {
            t0[0] = s[l] * _ms[l];
        } else {
            if (MSGS & 1) {
                params::poly_q zero = 0;
                t0[0] = zero - beta * _ms[l];
            } else {
                t0[0] = beta * _ms[l];
            }
        }

        if (l == 0) {
            alpha[0] = beta;
        } else {
            alpha[0] = s[l - 1];
        }
        alpha[1] = t0[0];
        lin_prover(y[l], _y[l], t[l], _t[l], u[l], c[l], d[l], alpha, key, r[l],
                   _r[l]);
    }
}


/**
 * @brief Executes the shuffle verification.
 *
 * This function implements the verifier's side of a shuffle proof protocol. It verifies
 * proofs for the shuffling of polynomial commitments.
 *
 * @param y Array of polynomials.
 * @param _y Array of polynomials.
 * @param t Array of polynomials.
 * @param _t Array of polynomials.
 * @param u Array of polynomials.
 * @param d Array of commitments.
 * @param s Array of polynomials.
 * @param c Array of commitments.
 * @param _ms Array of polynomials.
 * @param rho Array of polynomials.
 * @param key Commitment key.
 * @return Returns 1 if the verification is successful, otherwise 0.
 */
static int shuffle_verifier(params::poly_q y[MSGS][WIDTH],
                            params::poly_q _y[MSGS][WIDTH], params::poly_q t[MSGS],
                            params::poly_q _t[MSGS], params::poly_q u[MSGS], commit_t d[MSGS],
                            params::poly_q s[MSGS], commit_t c[MSGS], params::poly_q _ms[MSGS],
                            params::poly_q rho[SIZE], comkey_t & key) {

    // Initialize variables
    params::poly_q alpha[2], beta;
    vector < params::poly_q > t0(1);
    int result = 1;  // Initialize result to true (1)

    // Compute the hash for the shuffle
    shuffle_hash(beta, c, d, _ms, rho);

    // Iterate over all messages
    for (size_t l = 0; l < MSGS; l++) {
        // Compute t0 based on the current message
        if (l < MSGS - 1) {
            t0[0] = s[l] * _ms[l];
        } else {
            if (MSGS & 1) {
                params::poly_q zero = 0;
                t0[0] = zero - beta * _ms[l];
            } else {
                t0[0] = beta * _ms[l];
            }
        }

        // Set alpha values based on the current message
        if (l == 0) {
            alpha[0] = beta;
        } else {
            alpha[0] = s[l - 1];
        }
        alpha[1] = t0[0];

        // Verify the linear proof for the current message and update the result
        result &=
                lin_verifier(y[l], _y[l], t[l], _t[l], u[l], c[l], d[l], alpha,
                             key);
    }

    return result;  // Return the final result
}


/**
 * @brief Executes the run function for shuffling polynomial commitments.
 *
 * This function extends commitments, adjusts the key, and then runs the shuffle prover
 * and verifier functions.
 *
 * @param com Array of commitments.
 * @param m 2D vector of polynomials.
 * @param _m 2D vector of polynomials.
 * @param key Commitment key.
 * @param r 2D vector of polynomials.
 * @return Returns 1 if the shuffle verification is successful, otherwise 0.
 */
static int run(commit_t com[MSGS], vector < vector < params::poly_q >> m,
               vector < vector < params::poly_q >> _m, comkey_t & key,
               vector < params::poly_q > r[MSGS]) {

    // Initialize variables
    params::poly_q ms[MSGS], _ms[MSGS];
    commit_t d[MSGS], cs[MSGS];
    vector < params::poly_q > t0(1);
    params::poly_q one, t1, rho[SIZE], s[MSGS];
    params::poly_q y[MSGS][WIDTH], _y[MSGS][WIDTH], t[MSGS], _t[MSGS], u[MSGS];
    comkey_t _key;

    // Extend commitments and adjust key
    rho[0] = 1;
    rho[0].ntt_pow_phi();
    for (size_t j = 1; j < SIZE; j++) {
        rho[j] = nfl::uniform();  // Sample uniformly random values for rho
    }

    // Adjust commitments based on the new rho values
    for (size_t i = 0; i < MSGS; i++) {
        ms[i] = m[i][0];
        ms[i].ntt_pow_phi();
        _ms[i] = _m[i][0];
        _ms[i].ntt_pow_phi();
        cs[i].c1 = com[i].c1;
        cs[i].c2.resize(1);
        cs[i].c2[0] = com[i].c2[0] * rho[0];
        for (size_t j = 1; j < SIZE; j++) {
            cs[i].c2[0] = cs[i].c2[0] + com[i].c2[j] * rho[j];
            t1 = m[i][j];
            t1.ntt_pow_phi();
            ms[i] = ms[i] + t1 * rho[j];
            t1 = _m[i][j];
            t1.ntt_pow_phi();
            _ms[i] = _ms[i] + t1 * rho[j];
        }
    }

    // Adjust the key based on the new rho values
    for (size_t i = 0; i < HEIGHT; i++) {
        for (size_t j = HEIGHT; j < WIDTH; j++) {
            _key.A1[i][j - HEIGHT] = key.A1[i][j - HEIGHT];
        }
    }
    _key.A2[0][0] = 0;
    _key.A2[0][1] = rho[0];
    for (size_t j = 2; j < WIDTH; j++) {
        _key.A2[0][j] = key.A2[0][j];
    }
    for (size_t i = 1; i < SIZE; i++) {
        for (size_t j = 2; j < WIDTH; j++) {
            _key.A2[0][j] = _key.A2[0][j] + rho[i] * key.A2[i][j];
        }
    }

    // Run the shuffle prover with the adjusted commitments and key
    shuffle_prover(y, _y, t, _t, u, d, s, cs, ms, _ms, r, rho, _key);

    // Run the shuffle verifier with the adjusted commitments and key
    return shuffle_verifier(y, _y, t, _t, u, d, s, cs, _ms, rho, _key);
}


#ifdef MAIN

/**
 * @brief Test function for polynomial commitments.
 *
 * This function generates a commitment key, commits to random messages,
 * shuffles the messages, and then tests the correctness of polynomial inverse
 * and shuffle proof.
 */
static void test() {
    // Initialize commitment key, commitments, and message vectors
    comkey_t key;
    commit_t com[MSGS];
    vector < vector < params::poly_q >> m(MSGS), _m(MSGS);
    vector < params::poly_q > r[MSGS];

    // Generate commitment key
    bdlop_keygen(key);

    // Commit to random messages
    for (int i = 0; i < MSGS; i++) {
        m[i].resize(SIZE);
        for (int j = 0; j < SIZE; j++) {
            m[i][j] = nfl::ZO_dist();  // Sample a random polynomial
        }
        r[i].resize(WIDTH);
        bdlop_sample_rand(r[i]);  // Sample random values for r
        bdlop_commit(com[i], m[i], key, r[i]);  // Commit to the message
    }

    // Prover shuffles messages (using a circular shift for simplicity)
    for (int i = 0; i < MSGS; i++) {
        _m[i].resize(SIZE);
        for (int j = 0; j < SIZE; j++) {
            _m[i][j] = m[(i + 1) % MSGS][j];  // Circular shift
        }
    }

    // Test the correctness of polynomial inverse
    TEST_ONCE("polynomial inverse is correct") {
        params::poly_q alpha[2] = { nfl::uniform(), nfl::uniform() };

        poly_inverse(alpha[1], alpha[0]);
        alpha[0].ntt_pow_phi();
        alpha[1].ntt_pow_phi();
        alpha[0] = alpha[0] * alpha[1];
        alpha[0] = alpha[0] * alpha[1];
        alpha[0].invntt_pow_invphi();
        alpha[1].invntt_pow_invphi();
        TEST_ASSERT(alpha[0] == alpha[1], end);  // Check if the results are equal
    } TEST_END;

    // Test the consistency of the shuffle proof
    TEST_ONCE("shuffle proof is consistent") {
        TEST_ASSERT(run(com, m, _m, key, r) == 1, end);  // Check if the shuffle proof is valid
    } TEST_END;

    end:
    return;
}

/**
 * @brief Microbenchmark function for polynomial operations.
 *
 * This function benchmarks polynomial addition, multiplication, and inverse operations.
 */
static void microbench() {
    params::poly_q alpha[2] = { nfl::uniform(), nfl::uniform() };

    alpha[0].ntt_pow_phi();
    alpha[1].ntt_pow_phi();

    // Benchmark polynomial addition
    BENCH_BEGIN("Polynomial addition") {
            BENCH_ADD(alpha[0] = alpha[0] + alpha[1]);
        } BENCH_END;

    // Benchmark polynomial multiplication
    BENCH_BEGIN("Polynomial multiplication") {
            BENCH_ADD(alpha[0] = alpha[0] * alpha[1]);
        } BENCH_END;

    alpha[0].invntt_pow_invphi();
    // Benchmark polynomial inverse
    BENCH_BEGIN("Polynomial inverse") {
            BENCH_ADD(poly_inverse(alpha[1], alpha[0]));
        } BENCH_END;
}

/**
 * @brief Benchmarking function for various cryptographic operations.
 *
 * This function benchmarks the linear hash, linear proof, and linear verifier operations.
 * It also benchmarks the shuffle-proof for N messages.
 */
static void bench() {
    // Initialize commitment key, commitments, and message vectors
    comkey_t key;
    commit_t com[MSGS];
    vector < vector < params::poly_q >> m(MSGS), _m(MSGS);
    vector < params::poly_q > r[MSGS];
    params::poly_q y[WIDTH], _y[WIDTH], t, _t, u, alpha[2], beta;

    // Generate commitment key
    bdlop_keygen(key);

    // Commit to random messages
    for (int i = 0; i < MSGS; i++) {
        m[i].resize(SIZE);
        for (int j = 0; j < SIZE; j++) {
            m[i][j] = nfl::ZO_dist();  // Sample a random polynomial
        }
        r[i].resize(WIDTH);
        bdlop_sample_rand(r[i]);  // Sample random values for r
        bdlop_commit(com[i], m[i], key, r[i]);  // Commit to the message
    }

    // Prover shuffles messages (using a circular shift for simplicity)
    for (int i = 0; i < MSGS; i++) {
        _m[i].resize(SIZE);
        for (int j = 0; j < SIZE; j++) {
            _m[i][j] = m[(i + 1) % MSGS][j];  // Circular shift
        }
    }

    // Sample random values for alpha and transform them to NTT domain
    alpha[0] = nfl::ZO_dist();
    alpha[1] = nfl::ZO_dist();
    alpha[0].ntt_pow_phi();
    alpha[1].ntt_pow_phi();

    // Sample a challenge value for beta
    bdlop_sample_chal(beta);

    // Benchmark the linear hash operation
    BENCH_BEGIN("linear hash") {
            BENCH_ADD(lin_hash(beta, key, com[0], com[1], alpha, u, t, _t));
        } BENCH_END;

    // Benchmark the linear proof operation
    BENCH_BEGIN("linear proof") {
            BENCH_ADD(lin_prover(y, _y, t, _t, u, com[0], com[1], alpha, key, r[0], r[0]));
        } BENCH_END;

    // Benchmark the linear verifier operation
    BENCH_BEGIN("linear verifier") {
            BENCH_ADD(lin_verifier(y, _y, t, _t, u, com[0], com[1], alpha, key));
        } BENCH_END;

    // Benchmark the shuffle-proof for N messages
    BENCH_SMALL("shuffle-proof (N messages)", run(com, m, _m, key, r));
}

/**
 * @brief Main function for testing and benchmarking.
 *
 * @param argc Number of command-line arguments.
 * @param argv Array of command-line arguments.
 * @return int Returns 0 on successful execution.
 */
int main(int argc, char *argv[]) {
    // Display the header for tests
    printf("\n** Tests for lattice-based shuffle proof:\n\n");
    test();  // Run tests

    // Display the header for microbenchmarks
    printf("\n** Microbenchmarks for polynomial arithmetic:\n\n");
    microbench();  // Run microbenchmarks

    // Display the header for benchmarks
    printf("\n** Benchmarks for lattice-based shuffle proof:\n\n");
    bench();  // Run benchmarks

    return 0;  // Return 0 indicating successful execution
}

#endif
