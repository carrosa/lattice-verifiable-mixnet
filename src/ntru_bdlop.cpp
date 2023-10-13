#include "test.h"
#include "bench.h"
#include "common.h"
#include <sys/random.h>
#include <array>

/*============================================================================*/
/* Public definitions                                                         */
/*============================================================================*/

using namespace std;

/**
 * Test if the l2-norm is within bounds (4 * sigma * sqrt(N)).
 *
 * @param[in] r 			- the polynomial to compute the l2-norm.
 * @return the computed norm.
 */
// This function tests if the norm of a polynomial `r` is less than a certain bound.
// The bound is determined by the parameter `sigma_sqr`.
bool ntru_bdlop_test_norm(ntru_params::poly_q r, uint64_t sigma_sqr) {
    // Declare an array to store the coefficients of the polynomial `r`.
    array<mpz_t, ntru_params::poly_q::degree> coeffs;

    // Declare variables to store the norm of the polynomial, half of the moduli product, and a temporary variable.
    mpz_t norm, qDivBy2, tmp;

    // Initialize the variables `norm`, `qDivBy2`, and `tmp`.
    mpz_inits(norm, qDivBy2, tmp, nullptr);

    // Initialize the coefficients of the polynomial with a size that is four times the bits in the moduli product.
    for (size_t i = 0; i < ntru_params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], (ntru_params::poly_q::bits_in_moduli_product() << 2));
    }

    // Convert the polynomial `r` to its coefficient representation and store it in `coeffs`.
    r.poly2mpz(coeffs);

    // Compute half of the moduli product and store it in `qDivBy2`.
    mpz_fdiv_q_2exp(qDivBy2, ntru_params::poly_q::moduli_product(), 1);

    // Initialize the norm to zero.
    mpz_set_ui(norm, 0);

    // For each coefficient of the polynomial:
    for (size_t i = 0; i < ntru_params::poly_q::degree; i++) {
        // Center the coefficient around zero using the moduli product and `qDivBy2`.
        util::center(coeffs[i], coeffs[i], ntru_params::poly_q::moduli_product(), qDivBy2);

        // Compute the square of the coefficient and store it in `tmp`.
        mpz_mul(tmp, coeffs[i], coeffs[i]);

        // Add the squared coefficient to the norm.
        mpz_add(norm, norm, tmp);
    }

    // Compute the bound as (4 * sigma * sqrt(N))^2 = 16 * sigma^2 * N.
    uint64_t bound = 16 * sigma_sqr * ntru_params::poly_q::degree;
//    std::cout << "\n\n\nBOUND:\n" << bound << "\n";
//    gmp_printf("\n\nNORM:\n %Zd\n", norm);
    // Compare the computed norm with the bound. If the norm is less than the bound, `result` will be true.
    int result = mpz_cmp_ui(norm, bound) < 0;
//    std::cout << "\n\nRESULT:\n" << result << "\n\n";

    // Clear the memory used by the variables `norm`, `qDivBy2`, and `tmp`.
    mpz_clears(norm, qDivBy2, tmp, nullptr);

    // Clear the memory used by the coefficients of the polynomial.
    for (size_t i = 0; i < ntru_params::poly_q::degree; i++) {
        mpz_clear(coeffs[i]);
    }

    // Return the result of the comparison.
    return result;
}


// This function samples random polynomials and stores them in the vector `r`.
// Each polynomial is sampled from a zero distribution and then transformed using NTT
void ntru_bdlop_sample_rand(vector<ntru_params::poly_q> &r) {
    // Iterate over each element of the vector `r`.
    for (size_t i = 0; i < r.size(); i++) {
        // Sample a random polynomial from a zero distribution.
        r[i] = nfl::ZO_dist();

        // Transform the sampled polynomial using NTT (Number Theoretic Transform).
        // The `ntt_pow_phi()` function likely computes the polynomial in the NTT domain.
        r[i].ntt_pow_phi();
    }
}


// This function samples a challenge polynomial `f`.
// The challenge is constructed by subtracting two polynomials, `c0` and `c1`,
// each of which is sampled from a Hamming weight distribution with a specified number of non-zero coefficients.
void ntru_bdlop_sample_chal(ntru_params::poly_q &f) {
    // Declare two polynomials `c0` and `c1`.
    ntru_params::poly_q c0, c1;

    // Sample polynomial `c0` from a Hamming weight distribution.
    // The `nfl::hwt_dist {NONZERO}` likely samples a polynomial with a fixed number of non-zero coefficients,
    // where the number is specified by the `NONZERO` constant.
    c0 = nfl::hwt_dist{NTRU_NONZERO};

    // Similarly, sample polynomial `c1` from the same Hamming weight distribution.
    c1 = nfl::hwt_dist{NTRU_NONZERO};

    // Construct the challenge polynomial `f` by subtracting `c1` from `c0`.
    f = c0 - c1;

    // Transform the challenge polynomial `f` into the NTT (Number Theoretic Transform) domain.
    // The `ntt_pow_phi()` function computes the polynomial in the NTT domain.
    f.ntt_pow_phi();
}


// This function generates a key pair and stores it in the `key` structure.
void ntru_bdlop_keygen(ntru_comkey_t &key) {
    // Initialize a polynomial `one` with the value 1.
    ntru_params::poly_q one = 1;

    // Transform the polynomial `one` into the NTT (Number Theoretic Transform) domain.
    // The `ntt_pow_phi()` function computes the polynomial in the NTT domain.
    one.ntt_pow_phi();

    // Populate the `A1` matrix of the key:
    // The outer loop iterates over the rows of the matrix, defined by the constant `NTRU_HEIGHT`.
    for (size_t i = 0; i < NTRU_HEIGHT; i++) {
        // The inner loop iterates over the columns of the matrix, defined by the difference between `WIDTH` and `NTRU_HEIGHT`.
        for (int j = 0; j < WIDTH - NTRU_HEIGHT; j++) {
            // Each element of the matrix is sampled from a uniform distribution.
            key.A1[i][j] = nfl::uniform();
        }
    }

    // Populate the `A2` matrix of the key:
    // The outer loop iterates over the rows of the matrix, defined by the constant `NTRU_SIZE`.
    // The first inner loop sets the initial columns (up to `NTRU_HEIGHT + NTRU_SIZE`) of the matrix to 0.
    for (size_t j = 0; j < NTRU_HEIGHT + NTRU_SIZE; j++) {
        key.A2[j] = 0;
    }
    // Set the diagonal element (offset by `NTRU_HEIGHT`) to the polynomial `one`.
    key.A2[NTRU_HEIGHT] = one;

    // The second inner loop populates the remaining columns (from `NTRU_SIZE + NTRU_HEIGHT` to `WIDTH`)
    // with values sampled from a uniform distribution.
    for (size_t j = NTRU_SIZE + NTRU_HEIGHT; j < WIDTH; j++) {
        key.A2[j] = nfl::uniform();
    }
}

// This function computes a commitment `com` to a message `m` using a commitment key `key` and randomness `r`.
void ntru_bdlop_commit(ntru_commit_t &com, ntru_params::poly_q &m, ntru_comkey_t &key,
                       vector<ntru_params::poly_q> r) {
    // Declare a temporary polynomial `_m`.
    ntru_params::poly_q _m;

    // Initialize the first component of the commitment `com.c1` with the first element of the randomness vector `r`.
    com.c1 = r[0];

    // Compute the first component of the commitment:
    // The outer loop iterates over the rows of the matrix `A1` in the commitment key, defined by the constant `NTRU_HEIGHT`.
    for (size_t i = 0; i < NTRU_HEIGHT; i++) {
        // The inner loop iterates over the elements of the randomness vector `r`, excluding the first `NTRU_HEIGHT` elements.
        for (size_t j = 0; j < r.size() - NTRU_HEIGHT; j++) {
            // Update the first component of the commitment by adding the product of the matrix element `A1[i][j]` and the randomness element `r[j + NTRU_HEIGHT]`.
            com.c1 = com.c1 + key.A1[i][j] * r[j + NTRU_HEIGHT];
        }
    }

    com.c2 = 0;

    // The inner loop iterates over the elements of the randomness vector `r`.
    for (size_t j = 0; j < r.size(); j++) {
        // Update the current element of the second component of the commitment by adding the product of the matrix element `A2[i][j]` and the randomness element `r[j]`.
        com.c2 = com.c2 + key.A2[j] * r[j];
    }

    // Transform the current element of the message `m` into the NTT (Number Theoretic Transform) domain.
    _m = m;
    _m.ntt_pow_phi();

    // Update the current element of the second component of the commitment by adding the transformed message element `_m`.
    com.c2 = com.c2 + _m;
}


// This function attempts to open a commitment on a message `m` using the provided randomness `r`, commitment key `key`, and factor `f`.
// It returns an integer indicating the success (true) or failure (false) of the opening process.
int ntru_bdlop_open(ntru_commit_t &com, ntru_params::poly_q m, ntru_comkey_t &key,
                    vector<ntru_params::poly_q> r, ntru_params::poly_q &f) {
    // Declare temporary polynomials for computations.
    ntru_params::poly_q c1, _c1, c2, _c2, _m;
    int result = true;

    // Compute the first component of the commitment using the matrix `A1` from the commitment key and a subset of the randomness vector `r`.
    c1 = r[0];
    for (size_t i = 0; i < NTRU_HEIGHT; i++) {
        for (size_t j = 0; j < r.size() - NTRU_HEIGHT; j++) {
            c1 = c1 + key.A1[i][j] * r[j + NTRU_HEIGHT];
        }
    }

    // Compute the second component of the commitment using the matrix `A2` from the commitment key, the entire randomness vector `r`, and the message `m`.
    // The message is also multiplied by the factor `f` before being added to the commitment.
    c2 = 0;
    for (size_t j = 0; j < r.size(); j++) {
        c2 = c2 + key.A2[j] * r[j];
    }
    _m = m;
    _m.ntt_pow_phi();
    c2 = c2 + f * _m;

    // Multiply the first component of the provided commitment by the factor `f`.
    _c1 = f * com.c1;

    // Check if the computed first component matches the provided commitment's first component and if the sizes of the second components match.
    if (_c1 != c1) {
        result = false;
        cout << "ERROR: Commit opening failed test for c1" << endl;
    } else {
        // For each element of the message, multiply the corresponding element of the provided commitment's second component by the factor `f` and check if it matches the computed second component.
        _c2 = f * com.c2;
        if (_c2 != c2) {
            cout << "ERROR: Commit opening failed test for c2" << endl;
            result = false;
        }

        // Check the norm of each element of the randomness vector `r` against a threshold (likely related to the security of the commitment scheme).
        for (size_t i = 0; i < r.size(); i++) {
            c1 = r[i];
            c1.invntt_pow_invphi();
            if (!ntru_bdlop_test_norm(c1, 16 * NTRU_SIGMA_C)) {
                cout << "ERROR: Commit opening failed norm test" << endl;
                result = false;
                break;
            }
        }
    }

    // Return the result of the opening process.
    return result;
}

#ifdef MAIN

// Test function for single message commitments.
static void test1() {
    ntru_comkey_t key;  // Commitment key.
    ntru_commit_t com, _com;  // Commitments.

    // Generate a commitment key.
    ntru_bdlop_keygen(key);

    // Define randomness and message vectors.
    vector<ntru_params::poly_q> r(WIDTH), s(WIDTH);
    ntru_params::poly_q f;
    ntru_params::poly_q m = nfl::uniform();  // Single message.

    TEST_BEGIN("commitment for single messages can be generated and opened") {
        // Sample random values.
        ntru_bdlop_sample_rand(r);

        // Generate a commitment for the message.
        ntru_bdlop_commit(com, m, key, r);

        // Set factor to 1 and transform it.
        f = 1;
        f.ntt_pow_phi();

        // Test if the commitment can be opened.
        TEST_ASSERT(ntru_bdlop_open(com, m, key, r, f) == 1, end);

        // Sample a challenge.
        ntru_bdlop_sample_chal(f);

        // Multiply randomness by the challenge.
        for (size_t j = 0; j < r.size(); j++) {
            s[j] = f * r[j];
        }

        // Test if the commitment can be opened with the modified randomness.
        TEST_ASSERT(ntru_bdlop_open(com, m, key, s, f) == 1, end);
    }
    TEST_END;

    TEST_BEGIN("commitments for single messages are linearly homomorphic") {
        // Test linearity.
        ntru_params::poly_q rho = nfl::uniform();
        for (size_t j = 0; j < r.size(); j++) {
            r[j] = 0;
        }
        ntru_bdlop_commit(_com, rho, key, r);
        com.c1 = com.c1 - _com.c1;
        com.c2 = com.c2 - _com.c2;
        m = m - rho;

        // Test if the commitment can be opened after linear operations.
        TEST_ASSERT(ntru_bdlop_open(com, m, key, s, f) == 1, end);
    }
    TEST_END;

    end:
    return;
}

// Test function for multiple message commitments.
static void test2() {
    ntru_comkey_t key, _key;  // Commitment keys.
    ntru_commit_t c, com, _com;  // Commitments.

    // Generate a commitment key.
    ntru_bdlop_keygen(key);

    // Define randomness vectors and factors.
    vector<ntru_params::poly_q> r(WIDTH), s(WIDTH);
    ntru_params::poly_q f, one;
//    ntru_params::poly_q t;

    // Define messages.
    ntru_params::poly_q m = nfl::uniform();

    TEST_BEGIN("commitment for multiple messages can be generated and opened") {
        // Sample random values.
        ntru_bdlop_sample_rand(r);

        // Generate a commitment for the messages.
        ntru_bdlop_commit(com, m, key, r);

        // Set factor to 1 and transform it.
        f = 1;
        f.ntt_pow_phi();

        // Test if the commitment can be opened.
        TEST_ASSERT(ntru_bdlop_open(com, m, key, r, f) == 1, end);

        // Sample a challenge.
        ntru_bdlop_sample_chal(f);

        // Multiply randomness by the challenge.
        for (size_t j = 0; j < r.size(); j++) {
            s[j] = f * r[j];
        }

        // Test if the commitment can be opened with the modified randomness.
        TEST_ASSERT(ntru_bdlop_open(com, m, key, s, f) == 1, end);
    }
    TEST_END;

    TEST_ONCE("commitments do not open for the wrong keys")
    {
        // Sample random values.
        ntru_bdlop_sample_rand(r);

        // Generate a commitment for the messages.
        ntru_bdlop_commit(com, m, key, r);

        // Sample a challenge.
        ntru_bdlop_sample_chal(f);

        // Multiply randomness by the challenge.
        for (size_t j = 0; j < r.size(); j++) {
            s[j] = f * r[j];
        }

        // Generate a new commitment key.
        ntru_bdlop_keygen(key);

        // Test if the commitment can be opened with the wrong key.
        TEST_ASSERT(ntru_bdlop_open(com, m, key, s, f) == 0, end);
    }
    TEST_END;

    end:
    return;
}

// Benchmarking function for lattice-based commitments.
static void bench() {
    ntru_comkey_t key;  // Commitment key.
    ntru_commit_t com;  // Commitment.
    ntru_params::poly_q f;  // Factor.
    vector<ntru_params::poly_q> r(WIDTH), s(WIDTH);  // Randomness vectors.
    ntru_params::poly_q m;  // Messages.
    ntru_params::poly_p _m;  // Message in another domain.

    ntru_params::poly_q pk, sk, c;  // Public key and secret key and cipher text for NTRU encryption.

    // Benchmark key generation.
    BENCH_SMALL("bdlp_keygen", ntru_bdlop_keygen(key));

    // Benchmark challenge sampling.
    BENCH_BEGIN("ntru_bdlop_sample_chal")
        {
            BENCH_ADD(ntru_bdlop_sample_chal(f));
        }
    BENCH_END;

    // Benchmark randomness sampling.
    BENCH_BEGIN("ntru_bdlop_sample_rand")
        {
            BENCH_ADD(ntru_bdlop_sample_rand(r));
        }
    BENCH_END;

    // Benchmark commitment generation for messages.
    BENCH_BEGIN("ntru_bdlop_commit")
        {
            m = nfl::uniform();
            BENCH_ADD(ntru_bdlop_commit(com, m, key, r));
        }
    BENCH_END;

    // Benchmark commitment opening.
    BENCH_BEGIN("ntru_bdlop_open")
        {
            m = nfl::uniform();
            ntru_bdlop_commit(com, m, key, r);
            for (size_t j = 0; j < r.size(); j++) {
                s[j] = f * r[j];
            }
            BENCH_ADD(ntru_bdlop_open(com, m, key, s, f));
        }
    BENCH_END;
}

// Main function.
int main(int argc, char *arv[]) {
    printf("\n** Tests for lattice-based commitments:\n\n");
    test1();
    test2();

    printf("\n** Benchmarks for lattice-based commitments:\n\n");
    bench();

    return 0;
}

#endif
