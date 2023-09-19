#include <math.h>
#include <stdlib.h>

#include <flint/flint.h>
#include <flint/fmpz_mod_poly.h>

#include "common.h"
#include "test.h"
#include "bench.h"
#include "assert.h"
#include "blake3.h"

#define ETA         325
#define R           (HEIGHT+2)
#define V           (WIDTH+3)
#define TAU         1000

/* Had to move those to global to avoid overflowing the stack. */
params::poly_q A[R][V], s[TAU][V], t[TAU][R];
params::poly_big H0[V], _H[V], H[TAU][3];

/**
 * This function computes a hash of the given commitment `com` using the BLAKE3 algorithm.
 * It then uses the hash to seed a random number generator, from which it generates
 * random numbers `x`, `beta0`, and `beta` modulo `q`.
 *
 * @param x: Output random number.
 * @param beta0: Output random number.
 * @param beta: 2D array to store output random numbers.
 * @param q: Modulus for random number generation.
 * @param com: Commitment structure containing data to be hashed.
 */
static void pismall_hash(fmpz_t x, fmpz_t beta0, fmpz_t beta[TAU][3], fmpz_t q, commit_t &com) {
    // Define a buffer to store the BLAKE3 hash output.
    uint8_t hash[BLAKE3_OUT_LEN];

    // Initialize a BLAKE3 hasher.
    blake3_hasher hasher;

    // Define a random state for FLINT random number generation.
    flint_rand_t rand;

    // Seed values for the random number generator.
    ulong seed[2];

    // Initialize the random state.
    flint_randinit(rand);

    // Initialize the BLAKE3 hasher.
    blake3_hasher_init(&hasher);

    // Update the hasher with the data from `com.c1`.
    blake3_hasher_update(&hasher, (const uint8_t *) com.c1.data(), 16 * DEGREE);

    // Update the hasher with the data from each element of `com.c2`.
    for (size_t i = 0; i < com.c2.size(); i++) {
        blake3_hasher_update(&hasher, (const uint8_t *) com.c2[i].data(), 16 * DEGREE);
    }

    // Finalize the hash computation and store the result in `hash`.
    blake3_hasher_finalize(&hasher, hash, BLAKE3_OUT_LEN);

    // Split the hash into two parts and use them as seeds for the random number generator.
    memcpy(&seed[0], hash, sizeof(ulong));
    memcpy(&seed[1], hash + BLAKE3_OUT_LEN / 2, sizeof(ulong));

    // Seed the random number generator with the obtained seeds.
    flint_randseed(rand, seed[0], seed[1]);

    // Generate a random number `x` modulo `q`.
    fmpz_randm(x, rand, q);

    // Generate a random number `beta0` modulo `q`.
    fmpz_randm(beta0, rand, q);

    // Generate random numbers for the `beta` 2D array modulo `q`.
    for (size_t i = 0; i < TAU; i++) {
        for (size_t j = 0; j < 3; j++) {
            fmpz_randm(beta[i][j], rand, q);
        }
    }
}


/**
 * This function converts a polynomial represented in the `fmpz_mod_poly_t` format
 * to a polynomial in the `params::poly_q` format. It also applies certain transformations
 * to the polynomial in the `params::poly_q` format.
 *
 * @param out: Output polynomial in the `params::poly_q` format.
 * @param in: Input polynomial in the `fmpz_mod_poly_t` format.
 * @param ctx: Context for modular arithmetic operations.
 */
static void poly_to(params::poly_q &out, fmpz_mod_poly_t &in, const fmpz_mod_ctx_t ctx) {
    // Define an array to store the coefficients of the polynomial.
    // The size of the array is determined by the degree of the polynomial in the `params::poly_q` format.
    array<mpz_t, params::poly_q::degree> coeffs;

    // Initialize each coefficient in the array with a size that is four times the number of bits
    // in the product of the moduli for the `params::poly_q` polynomial.
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], (params::poly_q::bits_in_moduli_product() << 2));
    }

    // Extract the coefficients from the `fmpz_mod_poly_t` polynomial and store them in the `coeffs` array.
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        fmpz_mod_poly_get_coeff_mpz(coeffs[i], in, i, ctx);
    }

    // Convert the coefficients from the `coeffs` array to the `params::poly_q` polynomial format.
    out.mpz2poly(coeffs);

    // Apply the NTT (Number Theoretic Transform) power of phi transformation to the `params::poly_q` polynomial.
    out.ntt_pow_phi();

    // Clear the memory allocated for the coefficients in the `coeffs` array.
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_clear(coeffs[i]);
    }
}


/**
 * This function encodes two polynomials (`in0` and `in1`) and an array of integers (`in`)
 * into a single polynomial in the `params::poly_big` format. After encoding, certain transformations
 * are applied to the resulting polynomial.
 *
 * @param out: Output polynomial in the `params::poly_big` format.
 * @param in0: First input polynomial in the `fmpz_mod_poly_t` format.
 * @param in1: Second input polynomial in the `fmpz_mod_poly_t` format.
 * @param in: Array of integers to be encoded.
 * @param ctx: Context for modular arithmetic operations.
 */
static void poly_encode(params::poly_big &out, fmpz_mod_poly_t in0,
                        fmpz_mod_poly_t in1, fmpz_t in[ETA], const fmpz_mod_ctx_t ctx) {
    // Define an array to store the coefficients of the resulting polynomial.
    // The size of the array is determined by the degree of the polynomial in the `params::poly_big` format.
    array<mpz_t, params::poly_big::degree> coeffs;
    size_t i;

    // Initialize each coefficient in the array with a size that is four times the number of bits
    // in the product of the moduli for the `params::poly_big` polynomial.
    for (i = 0; i < params::poly_big::degree; i++) {
        mpz_init2(coeffs[i], (params::poly_big::bits_in_moduli_product() << 2));
    }

    // Extract the coefficients from the first input polynomial (`in0`) and store them in the beginning of the `coeffs` array.
    for (i = 0; i < params::poly_q::degree; i++) {
        fmpz_mod_poly_get_coeff_mpz(coeffs[i], in0, i, ctx);
    }

    // Continue extracting coefficients from the second input polynomial (`in1`) and store them in the next segment of the `coeffs` array.
    for (; i < 2 * params::poly_q::degree; i++) {
        fmpz_mod_poly_get_coeff_mpz(coeffs[i], in1, i - params::poly_q::degree, ctx);
    }

    // Continue by encoding the integers from the `in` array into the final segment of the `coeffs` array.
    for (; i < 2 * params::poly_q::degree + ETA; i++) {
        fmpz_get_mpz(coeffs[i], in[i - 2 * params::poly_q::degree]);
    }

    // Convert the coefficients from the `coeffs` array to the `params::poly_big` polynomial format.
    out.mpz2poly(coeffs);

    // Apply the NTT (Number Theoretic Transform) power of phi transformation to the `params::poly_big` polynomial.
    out.ntt_pow_phi();

    // Clear the memory allocated for the coefficients in the `coeffs` array.
    for (size_t i = 0; i < params::poly_big::degree; i++) {
        mpz_clear(coeffs[i]);
    }
}


/**
 * This function converts a polynomial from the `params::poly_q` format to the `fmpz_mod_poly_t` format.
 * The conversion involves certain transformations and coefficient manipulations.
 *
 * @param out: Output polynomial in the `fmpz_mod_poly_t` format.
 * @param in: Input polynomial in the `params::poly_q` format.
 * @param ctx: Context for modular arithmetic operations.
 */
static void poly_from(fmpz_mod_poly_t &out, params::poly_q &in,
                      const fmpz_mod_ctx_t ctx) {
    // Define an array to store the coefficients of the input polynomial.
    // The size of the array is determined by the degree of the polynomial in the `params::poly_q` format.
    array<mpz_t, params::poly_q::degree> coeffs;

    // Apply the inverse NTT (Number Theoretic Transform) power of inverse phi transformation to the input polynomial.
    in.invntt_pow_invphi();

    // Initialize each coefficient in the array with a size that is four times the number of bits
    // in the product of the moduli for the `params::poly_q` polynomial.
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], (params::poly_q::bits_in_moduli_product() << 2));
    }

    // Convert the coefficients from the `params::poly_q` polynomial format to the `coeffs` array.
    in.poly2mpz(coeffs);

    // Initialize the output polynomial to zero and ensure it has the appropriate length.
    fmpz_mod_poly_zero(out, ctx);
    fmpz_mod_poly_fit_length(out, params::poly_q::degree, ctx);

    // Set the coefficients of the output polynomial using the values in the `coeffs` array.
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        fmpz_mod_poly_set_coeff_mpz(out, i, coeffs[i], ctx);
    }

    // Apply the NTT (Number Theoretic Transform) power of phi transformation to the input polynomial.
    // This step essentially undoes the earlier inverse transformation.
    in.ntt_pow_phi();

    // Clear the memory allocated for the coefficients in the `coeffs` array.
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_clear(coeffs[i]);
    }
}

/**
 * This function sets up the polynomial interpolation for a small number of points.
 * It computes the Lagrange basis polynomials for given points.
 *
 * @param lag: Array of polynomials to store the computed Lagrange basis polynomials.
 * @param a: Array of evaluation points.
 * @param q: Modulus for arithmetic operations.
 * @param prng: Pseudo-random number generator.
 * @param ctx: Context for modular arithmetic operations.
 */
static void pismall_setup(fmpz_mod_poly_t lag[], fmpz_t a[TAU], fmpz_t &q,
                          flint_rand_t prng, const fmpz_mod_ctx_t ctx) {
    // Define an array to store the coefficients of the polynomial.
    array<mpz_t, params::poly_q::degree> coeffs;
    fmpz_t t, u;  // Temporary variables for arithmetic operations.
    fmpz xs[TAU]; // Temporary array to store evaluation points.

    // Initialize the temporary variables.
    fmpz_init(t);
    fmpz_init(u);

    // Initialize each coefficient in the array with a size that is four times the number of bits
    // in the product of the moduli for the `params::poly_q` polynomial.
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], (params::poly_q::bits_in_moduli_product() << 2));
    }

    // Initialize and set the values of the temporary evaluation points array.
    for (size_t i = 0; i < TAU; i++) {
        fmpz_init(&xs[i]);
        fmpz_set(&xs[i], a[i]);
    }

    // Compute the zeroth Lagrange basis polynomial, l_0(X) = Product(X - a_i) for all i.
    fmpz_mod_poly_product_roots_fmpz_vec(lag[0], xs, TAU, ctx);

    // Compute the remaining Lagrange basis polynomials.
    for (size_t i = 1; i <= TAU; i++) {
        // Swap the i-th evaluation point with the last one.
        fmpz_set(t, &xs[i - 1]);
        fmpz_set(&xs[i - 1], &xs[TAU - 1]);

        // Compute the polynomial product without the i-th evaluation point.
        fmpz_mod_poly_product_roots_fmpz_vec(lag[i], xs, TAU - 1, ctx);

        // Restore the swapped evaluation point.
        fmpz_set(&xs[i - 1], t);

        // Compute the denominator of the Lagrange basis polynomial.
        fmpz_set_ui(u, 1);
        for (size_t j = 0; j < TAU; j++) {
            if (i - 1 != j) {
                fmpz_sub(t, a[i - 1], a[j]);
                fmpz_mul(u, u, t);
            }
        }

        // Compute the inverse of the denominator modulo q.
        fmpz_invmod(u, u, q);

        // Multiply the Lagrange basis polynomial by the inverse of the denominator.
        fmpz_mod_poly_scalar_mul_fmpz(lag[i], lag[i], u, ctx);
    }

    // Clear the memory allocated for the temporary variables.
    fmpz_clear(t);
    fmpz_clear(u);
    for (size_t i = 0; i < TAU; i++) {
        fmpz_clear(&xs[i]);
    }
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_clear(coeffs[i]);
    }
}

/**
 * This function implements the prover
 *
 * @param com: Commitment structure.
 * @param x: Evaluation point.
 * @param f: Array of polynomials.
 * @param rf: Array of random values associated with f.
 * @param h: 2D array of polynomials.
 * @param rh: Array of random values associated with h.
 * @param rd: Vector of random values for the commitment.
 * @param key: Commitment key.
 * @param lag: Array of Lagrange basis polynomials.
 * @param prng: Pseudo-random number generator.
 * @param ctx: Context for modular arithmetic operations.
 * @return: Returns 1 upon successful completion.
 */
static int pismall_prover(commit_t &com, fmpz_t x, fmpz_mod_poly_t f[V],
                          fmpz_t rf[ETA], fmpz_mod_poly_t h[2][V], fmpz_t rh[ETA],
                          vector<params::poly_q> rd, comkey_t &key,
                          fmpz_mod_poly_t lag[TAU + 1], flint_rand_t prng,
                          const fmpz_mod_ctx_t ctx) {
    // Define arrays to store coefficients and other polynomial-related values.
    array<mpz_t, params::poly_q::degree> coeffs, coeffs0, v[3][TAU][V];
    fmpz_mod_poly_t poly, zero; // Temporary polynomial variables.
    fmpz_t t, u, q, y[TAU], beta0, beta[TAU][3], r0[ETA], r[TAU][3][ETA];
    fmpz_mod_ctx_t ctx_q; // Context for modular arithmetic in q.
    vector<params::poly_q> d; // Vector to store polynomials.
    params::poly_q s0[V]; // Array to store polynomials.

    // Initialize temporary variables and context.
    fmpz_init(t);
    fmpz_init(u);
    fmpz_init(q);
    fmpz_set_mpz(q, params::poly_q::moduli_product());
    fmpz_mod_ctx_init(ctx_q, q);
    fmpz_set_str(q, PRIMEQ, 10);
    fmpz_mod_poly_init(poly, ctx);
    fmpz_mod_poly_init(zero, ctx);

    // Initialize coefficient arrays.
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_init(coeffs[i]);
        mpz_init(coeffs0[i]);
    }

    // Initialize y, beta, and r arrays.
    for (size_t i = 0; i < TAU; i++) {
        fmpz_init(y[i]);
        for (size_t j = 0; j < 3; j++) {
            for (size_t k = 0; k < V; k++) {
                for (size_t l = 0; l < params::poly_q::degree; l++) {
                    mpz_init(v[j][i][k][l]);
                }
            }
        }
    }
    fmpz_init(beta0);
    for (size_t i = 0; i < TAU; i++) {
        for (size_t j = 0; j < 3; j++) {
            fmpz_init(beta[i][j]);
            for (size_t k = 0; k < ETA; k++) {
                fmpz_init(r[i][j][k]);
                fmpz_randm(r[i][j][k], prng, q);
            }
        }
    }
    for (size_t i = 0; i < ETA; i++) {
        fmpz_init(r0[i]);
        fmpz_randm(r0[i], prng, q);
        fmpz_randm(rh[i], prng, q);
    }

    // Randomly initialize h polynomials.
    for (size_t i = 0; i < V; i++) {
        fmpz_mod_poly_randtest(h[0][i], prng, params::poly_q::degree, ctx);
        fmpz_mod_poly_randtest(h[1][i], prng, params::poly_q::degree, ctx);
    }

    /* Generate s_0 = Z_q^vN, h = Z_q^2vn. */
    for (size_t i = 0; i < V; i++) {
        fmpz_mod_poly_randtest(poly, prng, params::poly_q::degree, ctx);
        poly_to(s0[i], poly, ctx);
    }

    /* Compute d = -A * s_0 and convert back from NTT due to commitment. */
    d.resize(R);
    for (int j = 0; j < R; j++) {
        d[j] = 0;
        for (int k = 0; k < V; k++) {
            d[j] = d[j] - A[j][k] * s0[k];
        }
        d[j].invntt_pow_invphi();
    }
    bdlop_commit(com, d, key, rd);

    // Compute the explicit polynomial operations for v_i,j.
    // This section involves a lot of polynomial arithmetic and transformations.
    // The operations are performed in the NTT domain and then transformed back.
    // The exact nature of these operations depends on the specifics of the protocol.
    /* Compute f \circ (f + 1) \circ (f - 1) explicitly by computing the v_i,j. */
    // Iterate over all polynomials indexed by 'k'
    for (size_t k = 0; k < V; k++) {
        // Convert s0[k] from NTT (Number Theoretic Transform) form to standard polynomial form
        s0[k].invntt_pow_invphi();
        // Convert s0[k] polynomial to its coefficient representation in 'coeffs0'
        s0[k].poly2mpz(coeffs0);
        // Convert s0[k] back to NTT form
        s0[k].ntt_pow_phi();

        // Iterate over all polynomials indexed by 'i'
        for (size_t i = 0; i < TAU; i++) {
            // Convert s[i][k] from NTT form to standard polynomial form
            s[i][k].invntt_pow_invphi();
            // Convert s[i][k] polynomial to its coefficient representation in 'coeffs'
            s[i][k].poly2mpz(coeffs);
            // Convert s[i][k] back to NTT form
            s[i][k].ntt_pow_phi();

            // Iterate over all coefficients of the polynomial
            for (size_t l = 0; l < params::poly_q::degree; l++) {
                // The following operations seem to be computing specific polynomial transformations
                // or evaluations based on the coefficients of s[i][k] and s0[k]. The exact nature
                // of these operations would depend on the specifics of the protocol or algorithm.

                // Compute v[0][i][k][l] based on coeffs[l]
                fmpz_set_mpz(t, coeffs[l]);
                fmpz_mod_sub_ui(u, t, 1, ctx_q);
                fmpz_mod_mul(t, t, u, ctx_q);
                fmpz_mod_add_ui(u, u, 2, ctx_q);
                fmpz_mod_mul(t, t, u, ctx_q);
                fmpz_get_mpz(v[0][i][k][l], t);

                // Compute v[1][i][k][l] based on coeffs[l] and coeffs0[l]
                fmpz_set_mpz(t, coeffs[l]);
                fmpz_set_mpz(u, coeffs0[l]);
                fmpz_mod_sub_ui(t, t, 1, ctx_q);
                fmpz_mod_mul(u, u, t, ctx_q);
                fmpz_mod_add_ui(t, t, 1, ctx_q);
                fmpz_mod_mul(u, u, t, ctx_q);
                fmpz_get_mpz(v[1][i][k][l], u);

                // Update v[1][i][k][l] based on coeffs[l] and coeffs0[l]
                fmpz_set_mpz(t, coeffs[l]);
                fmpz_set_mpz(u, coeffs0[l]);
                fmpz_mod_sub_ui(u, u, 1, ctx_q);
                fmpz_mod_mul(u, u, t, ctx_q);
                fmpz_mod_add_ui(t, t, 1, ctx_q);
                fmpz_mod_mul(u, u, t, ctx_q);
                fmpz_set_mpz(t, v[1][i][k][l]);
                fmpz_mod_add_fmpz(u, u, t, ctx_q);
                fmpz_get_mpz(v[1][i][k][l], u);

                // Update v[1][i][k][l] again based on coeffs[l] and coeffs0[l]
                fmpz_set_mpz(t, coeffs[l]);
                fmpz_set_mpz(u, coeffs0[l]);
                fmpz_mod_add_ui(u, u, 1, ctx_q);
                fmpz_mod_mul(u, u, t, ctx_q);
                fmpz_mod_sub_ui(t, t, 1, ctx_q);
                fmpz_mod_mul(u, u, t, ctx_q);
                fmpz_set_mpz(t, v[1][i][k][l]);
                fmpz_mod_add_fmpz(u, u, t, ctx_q);
                fmpz_get_mpz(v[1][i][k][l], u);

                // Compute v[2][i][k][l] based on coeffs[l] and coeffs0[l]
                fmpz_set_mpz(t, coeffs[l]);
                fmpz_set_mpz(u, coeffs0[l]);
                fmpz_mod_sub_ui(u, u, 1, ctx_q);
                fmpz_mod_mul(t, t, u, ctx_q);
                fmpz_mod_add_ui(u, u, 2, ctx_q);
                fmpz_mod_mul(t, t, u, ctx_q);
                fmpz_get_mpz(v[2][i][k][l], t);

                // Update v[2][i][k][l] based on coeffs[l] and coeffs0[l]
                fmpz_set_mpz(u, coeffs0[l]);
                fmpz_mod_sub_ui(t, u, 1, ctx_q);
                fmpz_mod_mul(u, u, t, ctx_q);
                fmpz_set_mpz(t, coeffs[l]);
                fmpz_mod_add_ui(t, t, 1, ctx_q);
                fmpz_mod_mul(u, u, t, ctx_q);
                fmpz_set_mpz(t, v[1][i][k][l]);
                fmpz_mod_add_fmpz(u, u, t, ctx_q);
                fmpz_get_mpz(v[2][i][k][l], u);

                // Update v[2][i][k][l] again based on coeffs[l] and coeffs0[l]
                fmpz_set_mpz(u, coeffs0[l]);
                fmpz_mod_add_ui(t, u, 1, ctx_q);
                fmpz_mod_mul(u, u, t, ctx_q);
                fmpz_set_mpz(t, coeffs[l]);
                fmpz_mod_sub_ui(t, t, 1, ctx_q);
                fmpz_mod_mul(u, u, t, ctx_q);
                fmpz_set_mpz(t, v[1][i][k][l]);
                fmpz_mod_add_fmpz(u, u, t, ctx_q);
                fmpz_get_mpz(v[2][i][k][l], u);
            }
        }
    }


    // The goal of this block is to encode certain polynomials as larger NTT (Number Theoretic Transform) forms.

// Iterate over all polynomials indexed by 'k'
    for (size_t k = 0; k < V; k++) {
        // Convert the polynomial s0[k] from its NTT representation to a standard polynomial representation
        poly_from(poly, s0[k], ctx);

        // Encode the polynomial 'poly' (which is s0[k] in standard form) and store the result in H0[k]
        poly_encode(H0[k], poly, zero, r0, ctx);

        // Encode the polynomials h[0][k] and h[1][k] and store the result in _H[k]
        poly_encode(_H[k], h[0][k], h[1][k], rh, ctx);

        // Iterate over all polynomials indexed by 'i'
        for (size_t i = 0; i < TAU; i++) {
            // Iterate over the three possible values of 'j'
            for (size_t j = 0; j < 3; j++) {
                // Set the coefficients of the polynomial 'poly' using values from the v array
                for (size_t l = 0; l < params::poly_q::degree; l++) {
                    fmpz_mod_poly_set_coeff_mpz(poly, l, v[j][i][k][l], ctx_q);
                }

                // If j is 0, convert the polynomial s[i][j] from its NTT representation to a standard polynomial representation
                // and then encode it with 'poly' and store the result in H[i][j]
                if (j == 0) {
                    poly_from(zero, s[i][j], ctx_q);
                    poly_encode(H[i][j], zero, poly, r[i][j], ctx_q);
                }

                // Zero out the polynomial 'zero'
                fmpz_mod_poly_zero(zero, ctx);

                // Encode the zero polynomial with 'poly' and store the result in H[i][j]
                poly_encode(H[i][j], zero, poly, r[i][j], ctx_q);
            }
        }
    }


    pismall_hash(x, beta0, beta, q, com);

    /* Compute the polynomial f as the product of s_0 and l_0 evaluated at x. */
    /* Compute f = s_0 * l_0(x). */
    // Evaluate the polynomial l_0 at the point x and store the result in y[0]
    fmpz_mod_poly_evaluate_fmpz(y[0], lag[0], x, ctx_q);

    // Iterate over all polynomials indexed by 'i'
    for (int i = 0; i < V; i++) {
        // Convert the polynomial s0[i] from its NTT representation to a standard polynomial representation
        poly_from(poly, s0[i], ctx_q);

        // Multiply the polynomial 'poly' (which is s0[i] in standard form) by the scalar y[0] and store the result in f[i]
        fmpz_mod_poly_scalar_mul_fmpz(f[i], poly, y[0], ctx_q);
    }

    /* Compute the polynomial f as the sum of products of each s_i and l_i evaluated at x. */
    /* Compute remaining f = f(x) = \Sum s_i * l_i(x). */
    // Start iterating from the second polynomial (since the first one was handled in the previous block)
    for (size_t i = 1; i <= TAU; i++) {
        // Evaluate the polynomial l_i at the point x and store the result in y[i - 1]
        fmpz_mod_poly_evaluate_fmpz(y[i - 1], lag[i], x, ctx_q);

        // Iterate over all polynomials indexed by 'j'
        for (int j = 0; j < V; j++) {
            // Convert the polynomial s[i - 1][j] from its NTT representation to a standard polynomial representation
            poly_from(poly, s[i - 1][j], ctx_q);

            // Multiply the polynomial 'poly' by the scalar y[i - 1] (which is the value of l_i evaluated at x)
            fmpz_mod_poly_scalar_mul_fmpz(poly, poly, y[i - 1], ctx_q);

            // Add the resulting polynomial to the existing f[j] to accumulate the sum
            fmpz_mod_poly_add(f[j], f[j], poly, ctx_q);
        }
    }


    /* Compute _rf = r0 * l_0(x) + Sum r_i,j * l_i(x) * l_0(x)^j */
    /* and _rh = r0 * beta_0 + Sum r_i,j * beta_i,j */
    // Iterate over each element indexed by 'k' up to ETA
    for (size_t k = 0; k < ETA; k++) {
        // Multiply r0[k] by y[0] and store the result in rf[k]
        fmpz_mod_mul(rf[k], r0[k], y[0], ctx);

        // Multiply r0[k] by beta0 and store the result in 't'
        fmpz_mod_mul(t, r0[k], beta0, ctx);

        // Add the value of 't' to rh[k]
        fmpz_mod_add(rh[k], rh[k], t, ctx_q);

        // Iterate over each polynomial up to TAU
        for (size_t i = 1; i <= TAU; i++) {
            // For each polynomial, iterate over 3 coefficients
            for (size_t j = 0; j < 3; j++) {
                // Multiply r[i - 1][j][k] by y[i - 1] and store the result in 't'
                fmpz_mod_mul(t, r[i - 1][j][k], y[i - 1], ctx_q);

                // Multiply 't' by y[i] 'j' times
                for (size_t l = 0; l < j; l++) {
                    fmpz_mod_mul(t, t, y[i], ctx_q);
                }

                // Add the value of 't' to rf[k]
                fmpz_mod_add(rf[k], rf[k], t, ctx_q);

                // Multiply r[i - 1][j][k] by beta[i - 1][j] and store the result in 't'
                fmpz_mod_mul(t, r[i - 1][j][k], beta[i - 1][j], ctx_q);

                // Add the value of 't' to rh[k]
                fmpz_mod_add(rh[k], rh[k], t, ctx_q);
            }
        }
    }


    /* Compute _h = h + s_0 * beta_0 + \sum beta_i,j * (\delta_i * si, v_i,j) */
    // Iterate over each polynomial indexed by 'k' up to V
    for (size_t k = 0; k < V; k++) {
        // Convert the polynomial s0[k] from its representation to a more standard form
        poly_from(poly, s0[k], ctx);

        // Multiply the polynomial 'poly' by scalar 'beta0'
        fmpz_mod_poly_scalar_mul_fmpz(poly, poly, beta0, ctx);

        // Add the resulting polynomial to h[0][k]
        fmpz_mod_poly_add(h[0][k], h[0][k], poly, ctx);

        // Iterate over each polynomial up to TAU
        for (size_t i = 1; i <= TAU; i++) {
            // For each polynomial, iterate over 3 coefficients
            for (size_t j = 0; j < 3; j++) {
                // If the coefficient index is 0
                if (j == 0) {
                    // Convert the polynomial s[i][k] from its representation to a more standard form
                    poly_from(poly, s[i][k], ctx_q);

                    // Multiply the polynomial 'poly' by scalar 'beta[i - 1][j]'
                    fmpz_mod_poly_scalar_mul_fmpz(poly, poly, beta[i - 1][j], ctx_q);

                    // Add the resulting polynomial to h[0][k]
                    fmpz_mod_poly_add(h[0][k], h[0][k], poly, ctx);
                }

                // Iterate over each coefficient of the polynomial up to its degree
                for (size_t l = 0; l < params::poly_q::degree; l++) {
                    // Set the l-th coefficient of 'poly' to the value in v[j][i - 1][k][l]
                    fmpz_mod_poly_set_coeff_mpz(poly, l, v[j][i - 1][k][l], ctx_q);
                }

                // Multiply the polynomial 'poly' by scalar 'beta[i - 1][j]'
                fmpz_mod_poly_scalar_mul_fmpz(poly, poly, beta[i - 1][j], ctx_q);

                // Add the resulting polynomial to h[1][k]
                fmpz_mod_poly_add(h[1][k], h[1][k], poly, ctx);
            }
        }
    }


    // Clear memory for all initialized variables.
    fmpz_clear(t);
    fmpz_clear(u);
    fmpz_clear(q);
    fmpz_mod_poly_clear(poly, ctx);
    fmpz_mod_poly_clear(zero, ctx);
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_clear(coeffs[i]);
        mpz_clear(coeffs0[i]);
    }
    for (size_t i = 0; i < TAU; i++) {
        fmpz_clear(y[i]);
        for (size_t j = 0; j < 3; j++) {
            for (size_t k = 0; k < V; k++) {
                for (size_t l = 0; l < params::poly_q::degree; l++) {
                    mpz_clear(v[j][i][k][l]);
                }
            }
        }
    }
    fmpz_clear(beta0);
    for (size_t i = 0; i < TAU; i++) {
        for (size_t j = 0; j < 3; j++) {
            fmpz_clear(beta[i][j]);
            for (size_t k = 0; k < ETA; k++) {
                fmpz_clear(r[i][j][k]);
            }
        }
    }
    for (size_t i = 0; i < ETA; i++) {
        fmpz_clear(r0[i]);
    }
    return 1;
}

// Function to verify a commitment
static int pismall_verifier(commit_t &com, fmpz_mod_poly_t f[V], fmpz_t rf[ETA],
                            vector<params::poly_q> rd, comkey_t &key,
                            fmpz_mod_poly_t lag[TAU + 1], flint_rand_t prng,
                            const fmpz_mod_ctx_t ctx) {
    // Declare arrays and variables for polynomial coefficients, polynomials, and other parameters
    array<mpz_t, params::poly_q::degree> coeffs;
    fmpz_mod_poly_t poly, r[R], _d[R];
    fmpz_mod_ctx_t ctx_q;
    fmpz_t q, y, x, beta0, beta[TAU][3];
    params::poly_q one = 1;
    vector<params::poly_q> m;

    // Initialize variables for the prime modulus, hash values, and beta values
    fmpz_init(q);
    fmpz_init(x);
    fmpz_init(y);
    fmpz_init(beta0);
    for (size_t i = 0; i < TAU; i++) {
        for (size_t j = 0; j < 3; j++) {
            fmpz_init(beta[i][j]);
        }
    }

    // Set the prime modulus value
    fmpz_set_str(q, PRIMEQ, 10);

    // Compute the hash of the commitment
    pismall_hash(x, beta0, beta, q, com);

    // Set the modulus for polynomial arithmetic and initialize the polynomial context
    fmpz_set_mpz(q, params::poly_q::moduli_product());
    fmpz_mod_ctx_init(ctx_q, q);
    fmpz_mod_poly_init(poly, ctx_q);

    // Initialize the coefficients array and polynomials r and _d
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], (params::poly_q::bits_in_moduli_product() << 2));
    }
    for (int i = 0; i < R; i++) {
        fmpz_mod_poly_init(r[i], ctx_q);
        fmpz_mod_poly_zero(r[i], ctx_q);
        fmpz_mod_poly_init(_d[i], ctx);
    }

    // Evaluate the lagrange polynomials at x and compute the sum of the products of the polynomials t and the evaluated lagrange polynomials
    for (int i = 1; i <= TAU; i++) {
        fmpz_mod_poly_evaluate_fmpz(y, lag[i], x, ctx_q);
        for (int j = 0; j < R; j++) {
            poly_from(poly, t[i - 1][j], ctx_q);
            fmpz_mod_poly_scalar_mul_fmpz(poly, poly, y, ctx_q);
            fmpz_mod_poly_add(r[j], r[j], poly, ctx_q);
        }
    }

    // Compute d = A * _f by iterating over each polynomial and subtracting the product of A and f from r
    for (int j = 0; j < V; j++) {
        for (int i = 0; i < R; i++) {
            poly_to(one, f[j], ctx_q);
            one = one * A[i][j];
            poly_from(poly, one, ctx_q);
            fmpz_mod_poly_sub(r[i], r[i], poly, ctx_q);
        }
    }

    // Multiply each polynomial in r by the inverse of the evaluated lagrange polynomial at 0
    fmpz_mod_poly_evaluate_fmpz(y, lag[0], x, ctx_q);
    fmpz_invmod(y, y, q);
    m.resize(R);
    for (int i = 0; i < R; i++) {
        fmpz_mod_poly_scalar_mul_fmpz(r[i], r[i], y, ctx_q);
        poly_to(m[i], r[i], ctx_q);
        m[i].invntt_pow_invphi();
    }
    one = 1;
    one.ntt_pow_phi();

    // Verify the commitment using the bdlop_open function
    int result = bdlop_open(com, m, key, rd, one);

    // Clear all initialized variables to free memory
    fmpz_clear(q);
    fmpz_clear(x);
    fmpz_clear(y);
    fmpz_mod_poly_clear(poly, ctx);
    for (int i = 0; i < R; i++) {
        fmpz_mod_poly_clear(r[i], ctx);
        fmpz_mod_poly_clear(_d[i], ctx);
    }
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_clear(coeffs[i]);
    }
    fmpz_clear(beta0);
    for (size_t i = 0; i < TAU; i++) {
        for (size_t j = 0; j < 3; j++) {
            fmpz_clear(beta[i][j]);
        }
    }

    // Return the result of the verification
    return result;
}


// Function to test
static void test(flint_rand_t rand) {
    // Declare variables for the prime modulus, random values, polynomials, and other parameters
    fmpz_t q, x, a[TAU], rf[ETA], rh[ETA];
    fmpz_mod_poly_t poly;
    fmpz_mod_ctx_t ctx, ctx_q;
    comkey_t key;
    commit_t com;
    vector<params::poly_q> rd;
    fmpz_mod_poly_t f[V], lag[TAU + 1], h[2][V];

    // Initialize the prime modulus and random values
    fmpz_init(q);
    fmpz_init(x);

    // Set the modulus for polynomial arithmetic and initialize the polynomial context
    fmpz_set_mpz(q, params::poly_q::moduli_product());
    fmpz_mod_ctx_init(ctx_q, q);

    // Set the prime modulus value and initialize the polynomial context
    fmpz_set_str(q, PRIMEQ, 10);
    fmpz_mod_ctx_init(ctx, q);
    fmpz_mod_poly_init(poly, ctx);

    // Initialize polynomials f and h
    for (int i = 0; i < V; i++) {
        fmpz_mod_poly_init(f[i], ctx);
        fmpz_mod_poly_init(h[0][i], ctx);
        fmpz_mod_poly_init(h[1][i], ctx);
    }

    // Initialize and set random values for a
    for (size_t i = 0; i < TAU; i++) {
        fmpz_init(a[i]);
        fmpz_randm(a[i], rand, q);
    }

    // Generate random polynomials and store them in A
    for (int i = 0; i < R; i++) {
        for (int j = 0; j < V; j++) {
            fmpz_mod_poly_randtest(poly, rand, params::poly_q::degree, ctx);
            poly_to(A[i][j], poly, ctx);
        }
    }

    // Initialize rf and rh
    for (size_t i = 0; i < ETA; i++) {
        fmpz_init(rf[i]);
        fmpz_init(rh[i]);
    }

    /* Create a total of TAU relations t_i = A * s_i */
    for (int i = 0; i < TAU; i++) {
        fmpz_mod_poly_init(lag[i], ctx);
        for (int j = 0; j < V; j++) {
            s[i][j] = nfl::hwt_dist{
                    2 * DEGREE / 3};
            s[i][j].ntt_pow_phi();
        }
        for (int j = 0; j < R; j++) {
            t[i][j] = 0;
            for (int k = 0; k < V; k++) {
                t[i][j] = t[i][j] + A[j][k] * s[i][k];
            }
        }
    }
    fmpz_mod_poly_init(lag[TAU], ctx);

    // Resize rd to R
    rd.resize(R);

    // Test the conversion correctness
    TEST_BEGIN("conversion is correct") {
        // Check if the conversion between polynomial and coefficient representation is correct
        for (int i = 0; i < R; i++) {
            for (int j = 0; j < V; j++) {
                poly_from(f[j], A[i][j], ctx);
                poly_to(rd[0], f[j], ctx);
                TEST_ASSERT(rd[0] == A[i][j], end);
            }
        }

        // Check if polynomial multiplication is consistent with the coefficient representation
        fmpz_mod_poly_zero(poly, ctx_q);
        fmpz_mod_poly_set_coeff_ui(poly, params::poly_q::degree, 1, ctx_q);
        fmpz_mod_poly_set_coeff_ui(poly, 0, 1, ctx_q);

        for (int i = 0; i < R; i++) {
            for (int j = 0; j < V; j++) {
                rd[0] = A[i][j] * A[i][j];
                poly_from(f[j], A[i][j], ctx_q);
                fmpz_mod_poly_mulmod(f[j], f[j], f[j], poly, ctx_q);
                poly_to(rd[1], f[j], ctx_q);
                TEST_ASSERT(rd[0] == rd[1], end);
            }
        }
    }
    TEST_END;

    // Setup
    pismall_setup(lag, a, q, rand, ctx);

    // Generate the commitment key and sample random values for rd
    bdlop_keygen(key);
    bdlop_sample_rand(rd);

    // Test if the AEX proof is consistent
    TEST_ONCE("AEX proof is consistent")
    {
        pismall_prover(com, x, f, rf, h, rh, rd, key, lag, rand, ctx);
        TEST_ASSERT(pismall_verifier(com, f, rf, rd, key, lag, rand, ctx) == 1,
                    end);
    }
    TEST_END;

    // Print benchmark results for the lattice-based AEX proof
    printf("\n** Benchmarks for lattice-based AEX proof:\n\n");
    BENCH_SMALL("pismall_setup", pismall_setup(lag, a, q, rand, ctx));
    BENCH_SMALL("pismall_prover", pismall_prover(com, x, f, rf, h, rh, rd, key, lag,
                                                 rand, ctx));
    BENCH_SMALL("pismall_verifier", pismall_verifier(com, f, rf, rd, key, lag, rand,
                                                     ctx));

    // Cleanup: Clear all initialized variables to free memory
    end:
    for (int i = 0; i < V; i++) {
        fmpz_mod_poly_clear(f[i], ctx);
        fmpz_mod_poly_clear(h[0][i], ctx);
        fmpz_mod_poly_clear(h[1][i], ctx);
    }
    for (int i = 0; i <= TAU; i++) {
        fmpz_mod_poly_clear(lag[i], ctx);
    }
    for (size_t i = 0; i < TAU; i++) {
        fmpz_clear(a[i]);
    }
    for (size_t i = 0; i < ETA; i++) {
        fmpz_clear(rf[i]);
        fmpz_clear(rh[i]);
    }
    fmpz_mod_poly_clear(poly, ctx);
    fmpz_mod_ctx_clear(ctx);
    fmpz_clear(q);
    fmpz_clear(x);
}


int main() {
    flint_rand_t rand;
    flint_randinit(rand);

    printf("\n** Tests for lattice-based AEX proof:\n\n");
    test(rand);
    flint_randclear(rand);
}
