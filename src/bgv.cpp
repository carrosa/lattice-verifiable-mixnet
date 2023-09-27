#include "test.h"
#include "bench.h"
#include "common.h"
#include <sys/random.h>
//#include "sample_z_small.h"
//#include "flint/flint.h"
//#include "flint/fmpz_poly.h"
//#include "flint/fmpz_mod_poly.h"

#include <math.h>
#include <stdlib.h>

#include <flint/flint.h>
#include <flint/fmpz_mod_poly.h>

#include "blake3.h"
#include "common.h"
#include "test.h"
#include "bench.h"
#include "assert.h"
#include "sample_z_small.h"


// Function to sample a message for the GHL encryption scheme.
void ghl_sample_message(params::poly_q &m) {
    // Initialize an array to store polynomial coefficients.
    std::array<mpz_t, DEGREE> coeffs;
    uint64_t buf;

    // Calculate the number of bits required for the moduli product.
    size_t bits_in_moduli_product = params::poly_p::bits_in_moduli_product();

    // Initialize the coefficients with the appropriate bit size.
    for (size_t i = 0; i < params::poly_p::degree; i++) {
        mpz_init2(coeffs[i], bits_in_moduli_product << 2);
    }

    // Sample random coefficients for the polynomial.
    for (size_t j = 0; j < params::poly_p::degree / 2; j += 32) {
        ssize_t grret = getrandom(&buf, sizeof(buf), 0);
        if (grret == -1) {
            throw "Could not generate a random seed. Check if something is wrong with system prng.";
        }
        for (size_t k = 0; k < 64; k += 2) {
            mpz_set_ui(coeffs[j + k / 2], (buf >> k) % PRIMEP);
        }
    }

    // Convert the coefficients to a polynomial.
    m.mpz2poly(coeffs);

    // Clear the memory used by the coefficients.
    for (size_t i = 0; i < params::poly_p::degree; i++) {
        mpz_clear(coeffs[i]);
    }
}

// Function to sample a message for the BGV encryption scheme.
void bgv_sample_message(params::poly_p &m) {
    // Similar to ghl_sample_message but for BGV scheme.
    std::array<mpz_t, DEGREE> coeffs;
    uint64_t buf;

    size_t bits_in_moduli_product = params::poly_p::bits_in_moduli_product();
    for (size_t i = 0; i < params::poly_p::degree; i++) {
        mpz_init2(coeffs[i], bits_in_moduli_product << 2);
    }

    for (size_t j = 0; j < params::poly_p::degree; j += 32) {
        ssize_t grret = getrandom(&buf, sizeof(buf), 0);
        if (grret == -1) {
            throw "Could not generate a random seed. Check if something is wrong with system prng.";
        }
        for (size_t k = 0; k < 64; k += 2) {
            mpz_set_ui(coeffs[j + k / 2], (buf >> k) % PRIMEP);
        }
    }
    m.mpz2poly(coeffs);

    for (size_t i = 0; i < params::poly_p::degree; i++) {
        mpz_clear(coeffs[i]);
    }
}

// Function to reduce the coefficients of a polynomial modulo PRIMEP.
void bgv_rdc(params::poly_p &m) {
    std::array<mpz_t, DEGREE> coeffs;

    size_t bits_in_moduli_product = params::poly_p::bits_in_moduli_product();
    for (size_t i = 0; i < params::poly_p::degree; i++) {
        mpz_init2(coeffs[i], bits_in_moduli_product << 2);
    }

    // Convert polynomial to coefficients.
    m.poly2mpz(coeffs);

    // Reduce each coefficient modulo PRIMEP.
    for (size_t j = 0; j < params::poly_p::degree; j++) {
        mpz_mod_ui(coeffs[j], coeffs[j], PRIMEP);
    }

    // Convert the coefficients back to a polynomial.
    m.mpz2poly(coeffs);

    for (size_t i = 0; i < params::poly_p::degree; i++) {
        mpz_clear(coeffs[i]);
    }
}

// Function to generate a key pair for the BGV encryption scheme.
void bgv_keygen(bgvkey_t &pk, params::poly_q &sk) {
    // Have to change to work with NTRU
    params::poly_q e = nfl::ZO_dist();  // Sample error polynomial.
    pk.a = nfl::uniform();  // Sample uniform polynomial for public key.

    sk = nfl::ZO_dist();  // Sample secret key.
    sk.ntt_pow_phi();  // Transform secret key.
    e.ntt_pow_phi();  // Transform error polynomial.

    pk.b = pk.a * sk;
    for (size_t i = 0; i < PRIMEP; i++) {
        pk.b = pk.b + e;
    }
}

// TODO - MOVE BELOW TO ntru.cpp

/**
 * INFO - THIS FUNCTION IS COPIED FROM shuffle.cpp
 * @brief Computes the multiplicative inverse of a polynomial in a given field.
 *
 * @param inv Output parameter to store the computed inverse polynomial.
 * @param p Input polynomial for which the inverse is to be computed.
 */
void poly_inverse(params::poly_q &inv, params::poly_q p) {
    // Declare an array to store coefficients of the polynomial
    std::array<mpz_t, params::poly_q::degree> coeffs;
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

// Function to genereate a key pair for the NTRU encryption scheme
void ntru_keygen(params::poly_q &pk, params::poly_q &sk) {
    // Sample f and g
    params::poly_q f;
    params::poly_q g;
    array<mpz_t, params::poly_q::degree> coeffs_f;
    array<mpz_t, params::poly_q::degree> coeffs_g;
    for (size_t k = 0; k < params::poly_q::degree; k++) {
        int64_t coeff_f;
        unsigned long int f_mod_p;
        do {
            coeff_f = sample_z(0.0, SIGMA_NTRU);
        } while ((k == 0 && coeff_f % PRIMEP != 1) || (k >= 1 && coeff_f % PRIMEP != 0));
        int64_t coeff_g = sample_z(0.0, SIGMA_NTRU);
        mpz_set_si(coeffs_f[k], coeff_f);
        mpz_set_si(coeffs_g[k], coeff_g);
    }
    // Multiply
    f.mpz2poly(coeffs_f);
    g.mpz2poly(coeffs_g);
    poly_inverse(f, f);

    f.ntt_pow_phi();
    g.ntt_pow_phi();

    pk = f * g;
    sk = f;
}

// TODO - MOVE ABOVE TO ntru.cpp


// Function to split a secret key into multiple shares.
void bgv_keyshare(params::poly_q s[], size_t shares, params::poly_q &sk) {
    params::poly_q t = sk;
    for (size_t i = 1; i < shares; i++) {
        s[i] = nfl::uniform();  // Sample a random polynomial.
        s[i].ntt_pow_phi();  // Transform the polynomial.
        t = t - s[i];  // Subtract the share from the total.
    }
    s[0] = t;  // The first share is the remaining value.
}

void ntru_keyshare(params::poly_q s[], size_t shares, params::poly_q &sk)  {
    // Exactly the same as BGV
    // TODO
    params::poly_q t = sk;
    for (size_t i = 1; i < shares; i++) {
        s[i] = nfl::uniform();
        s[i].ntt_pow_phi();
        t = t-s[i];
    }
    s[0] = t;
}

void ntru_encrypt(params::poly_q &c, params::poly_q &pk, params::poly_p &m) {
    // Initialize s, e
    params::poly_q s = nfl::ZO_dist();
    params::poly_q e = nfl::ZO_dist();

    // NTT transformations
    s.ntt_pow_phi();
    e.ntt_pow_phi();

    // Transform m to be in Rq
    params::poly_q m_q;
    array<mpz_t, params::poly_p::degree> coeffs_m = m.poly2mpz();
    m_q.mpz2poly(coeffs_m);
    m_q.ntt_pow_phi();

    // pk is already in NTT domain, so we don't have to transform that
    c = m_q;
    for (size_t i = 0; i < PRIMEP; i++) {
        c = c + pk * s + e;
    }
}

// Function to encrypt a message using the BGV encryption scheme.
void bgv_encrypt(bgvenc_t &c, bgvkey_t &pk, params::poly_p &m) {
    // Change this to ntru
    std::array<mpz_t, DEGREE> coeffs;
    params::poly_q e1 = nfl::ZO_dist();  // Sample first error polynomial.
    params::poly_q e2 = nfl::ZO_dist();  // Sample second error polynomial.
    params::poly_q r = nfl::ZO_dist();  // Sample random polynomial.

    e1.ntt_pow_phi();  // Transform first error polynomial.
    e2.ntt_pow_phi();  // Transform second error polynomial.
    r.ntt_pow_phi();  // Transform random polynomial.

    c.u = pk.a * r;
    c.v = pk.b * r;
    for (size_t i = 0; i < PRIMEP; i++) {
        c.u = c.u + e1;
        c.v = c.v + e2;
    }

    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], params::poly_q::bits_in_moduli_product() << 2);
    }

    // Convert the message polynomial to coefficients.
    m.poly2mpz(coeffs);

    // Convert the coefficients to a polynomial and transform it.
    r.mpz2poly(coeffs);
    r.ntt_pow_phi();

    // Add the transformed polynomial to the second component of the ciphertext.
    c.v = c.v + r;

    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_clear(coeffs[i]);
    }
}

// Function to add two ciphertexts together.
void bgv_add(bgvenc_t &c, bgvenc_t &d, bgvenc_t &e) {
    c.u = d.u + e.u;  // Add the first components.
    c.v = d.v + e.v;  // Add the second components.
}

void ntru_add(params::poly_q &c, params::poly_q &c1, params::poly_q &c2) {
    c = c1 + c2;
}

void ntru_decrypt() {
    // TODO
}

// Function to decrypt a ciphertext using the BGV encryption scheme.
void bgv_decrypt(params::poly_p &m, bgvenc_t &c, params::poly_q &sk) {
    // Change to ntru
    std::array<mpz_t, DEGREE> coeffs;
    params::poly_q t = c.v - sk * c.u;  // Compute the decryption polynomial.
    mpz_t qDivBy2;  // Variable to store half the moduli product.

    mpz_init(qDivBy2);
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], params::poly_q::bits_in_moduli_product() << 2);
    }

    // Compute half the moduli product.
    mpz_fdiv_q_2exp(qDivBy2, params::poly_q::moduli_product(), 1);

    // Transform the decryption polynomial back to coefficient space.
    t.invntt_pow_invphi();

    // Convert the polynomial to coefficients.
    t.poly2mpz(coeffs);

    // Center and reduce each coefficient modulo PRIMEP.
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        util::center(coeffs[i], coeffs[i], params::poly_q::moduli_product(),
                     qDivBy2);
        mpz_mod_ui(coeffs[i], coeffs[i], PRIMEP);
    }

    // Convert the coefficients back to a polynomial.
    m.mpz2poly(coeffs);

    mpz_clear(qDivBy2);
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_clear(coeffs[i]);
    }
}

void ntru_distdec() {
    // TODO
}

// Function to perform distributed decryption of a ciphertext.
void bgv_distdec(params::poly_q &tj, bgvenc_t &c, params::poly_q &sj) {
    // Change to ntru
    std::array<mpz_t, DEGREE> coeffs;
    params::poly_q mj, Ej;  // Intermediate polynomials.
    mpz_t qDivBy2, bound;  // Variables to store half the moduli product and a bound.

    mpz_init(qDivBy2);
    mpz_init(bound);
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], params::poly_q::bits_in_moduli_product() << 2);
    }

    Ej = nfl::uniform();
    Ej.poly2mpz(coeffs);
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        util::center(coeffs[i], coeffs[i], params::poly_q::moduli_product(),
                     qDivBy2);
        mpz_mod(coeffs[i], coeffs[i], bound);
    }
    Ej.mpz2poly(coeffs);
    Ej.ntt_pow_phi();
    mj = sj * c.u;
    tj = mj;
    for (size_t i = 0; i < PRIMEP; i++) {
        tj = tj + Ej;
    }
}

void ntru_comb() {

}

// Function to combine the results of distributed decryption.
void bgv_comb(params::poly_p &m, bgvenc_t &c, params::poly_q t[],
              size_t shares) {
    std::array<mpz_t, DEGREE> coeffs;
    params::poly_q v;  // Intermediate polynomial.
    mpz_t qDivBy2;  // Variable to store half the moduli product.

    mpz_init(qDivBy2);
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], params::poly_q::bits_in_moduli_product() << 2);
    }

    // Compute half the moduli product.
    mpz_fdiv_q_2exp(qDivBy2, params::poly_q::moduli_product(), 1);

    // Subtract the first share from the second component of the ciphertext.
    v = c.v - t[0];

    // Subtract the remaining shares.
    for (size_t i = 1; i < shares; i++) {
        v = v - t[i];
    }

    // Transform the polynomial back to coefficient space.
    v.invntt_pow_invphi();

    // Convert the polynomial to coefficients.
    v.poly2mpz(coeffs);

    // Center and reduce each coefficient modulo PRIMEP.
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        util::center(coeffs[i], coeffs[i], params::poly_q::moduli_product(),
                     qDivBy2);
        mpz_mod_ui(coeffs[i], coeffs[i], PRIMEP);
    }

    // Convert the coefficients back to a polynomial.
    m.mpz2poly(coeffs);

    mpz_clear(qDivBy2);
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_clear(coeffs[i]);
    }
}

#ifdef MAIN
// Function to test the BGV encryption scheme.
static void test() {
    bgvkey_t pk;  // Public key.
    params::poly_q sk, s[PARTIES], t[PARTIES], acc, mq;  // Secret key, shares, intermediate results, accumulator, and another polynomial.
    params::poly_p m, _m;  // Message polynomials.
    bgvenc_t c1, c2;  // Ciphertexts.

    bgv_keygen(pk, sk);  // Generate a key pair.
    bgv_sample_message(m);  // Sample a message.

    // Test that BGV encryption is consistent.
    TEST_BEGIN("BGV encryption is consistent") {
        bgv_encrypt(c1, pk, m);  // Encrypt the message.
        bgv_decrypt(_m, c1, sk);  // Decrypt the ciphertext.
        TEST_ASSERT(m - _m == 0, end);  // Check that the decrypted message matches the original.

        bgv_sample_message(m);  // Sample another message.
        bgv_decrypt(_m, c1, sk);  // Decrypt the previous ciphertext.
        TEST_ASSERT(m - _m != 0, end);  // Check that the decrypted message does not match the new message.

        bgv_encrypt(c1, pk, m);  // Encrypt the new message.
        bgv_keygen(pk, sk);  // Generate a new key pair.
        bgv_decrypt(_m, c1, sk);  // Decrypt the ciphertext with the new secret key.
        TEST_ASSERT(m - _m != 0, end);  // Check that the decryption is not successful with the wrong key.
    } TEST_END;

    // Test that BGV encryption is additively homomorphic.
    TEST_BEGIN("BGV encryption is additively homomorphic") {
        bgv_sample_message(m);  // Sample a message.
        bgv_sample_message(_m);  // Sample another message.
        bgv_encrypt(c1, pk, m);  // Encrypt the first message.
        bgv_encrypt(c2, pk, _m);  // Encrypt the second message.
        bgv_add(c1, c1, c2);  // Add the two ciphertexts.
        m = m + _m;  // Add the two messages.
        bgv_rdc(m);  // Reduce the coefficients of the result.
        bgv_decrypt(_m, c1, sk);  // Decrypt the combined ciphertext.
        TEST_ASSERT(m - _m == 0, end);  // Check that the decrypted message matches the sum of the original messages.
    } TEST_END;

    // Test that BGV distributed decryption is consistent.
    TEST_BEGIN("BGV distributed decryption is consistent") {
        bgv_keygen(pk, sk);  // Generate a key pair.
        bgv_encrypt(c1, pk, m);  // Encrypt a message.
        bgv_keyshare(s, PARTIES, sk);  // Split the secret key into shares.

        // Combine the shares to check that they match the original secret key.
        acc = s[0];
        for (size_t j = 1; j < PARTIES; j++) {
            acc = acc + s[j];
        }
        TEST_ASSERT(sk - acc == 0, end);  // Check that the combined shares match the original secret key.

        // Perform distributed decryption using each share.
        for (size_t j = 0; j < PARTIES; j++) {
            bgv_distdec(t[j], c1, s[j]);
        }

        // Combine the results of distributed decryption.
        bgv_comb(_m, c1, t, PARTIES);
        TEST_ASSERT(m - _m == 0, end);  // Check that the combined result matches the original message.
    } TEST_END;

  end:
    return;
}

static void ntru_test() {
    params::poly_q sk, pk;
    params::poly_p m, _m;
    params::poly_q c1, c2;

    ntru_keygen(pk, sk);
}

// Function to benchmark the BGV encryption scheme.
static void bench() {
    bgvkey_t pk;  // Public key.
    params::poly_q sk, s[PARTIES], t[PARTIES], acc;  // Secret key, shares, intermediate results, and accumulator.
    params::poly_p m, _m;  // Message polynomials.
    bgvenc_t c;  // Ciphertext.

    // Benchmark key generation.
    BENCH_SMALL("bgv_keygen", bgv_keygen(pk, sk));

    // Benchmark message sampling.
    BENCH_BEGIN("bgv_sample_message") {
        BENCH_ADD(bgv_sample_message(m));
    } BENCH_END;

    // Benchmark encryption.
    BENCH_BEGIN("bgv_encrypt") {
        BENCH_ADD(bgv_encrypt(c, pk, m));
    } BENCH_END;

    // Benchmark decryption.
    BENCH_BEGIN("bgv_decrypt") {
        bgv_encrypt(c, pk, m);
        BENCH_ADD(bgv_decrypt(_m, c, sk));
    } BENCH_END;

    // Benchmark ciphertext addition.
    BENCH_BEGIN("bgv_add") {
        BENCH_ADD(bgv_add(c, c, c));
    } BENCH_END;

    // Benchmark distributed decryption.
    bgv_keyshare(t, PARTIES, sk);
    BENCH_BEGIN("bgv_distdec") {
        bgv_encrypt(c, pk, m);
        BENCH_ADD(bgv_distdec(t[0], c, s[0]));
    } BENCH_END;

    // Benchmark the combination of distributed decryption results.
    for (size_t i = 1; i < PARTIES; i++) {
        bgv_distdec(t[i], c, s[i]);
    }
    BENCH_BEGIN("bgv_comb") {
        BENCH_ADD(bgv_comb(_m, c, t, PARTIES));
    } BENCH_END;
}

// Main function to run tests and benchmarks.
int main(int argc, char *arv[]) {
    //printf("\n** Tests for BGV encryption:\n\n");
    //test();

    printf("\n** Tests for NTRU encryption:\n\n");
    ntru_test();

    //printf("\n** Benchmarks for BGV encryption:\n\n");
    //bench();

    return 0;
}
#endif
