#include "test.h"
#include "bench.h"
#include "common.h"
#include <sys/random.h>

// Function to sample a message for the GHL encryption scheme.
void ghl_sample_message(params::poly_q & m) {
    // Initialize an array to store polynomial coefficients.
    std::array < mpz_t, DEGREE > coeffs;
    uint64_t buf;

    // Calculate the number of bits required for the moduli product.
    size_t bits_in_moduli_product = params::poly_p::bits_in_moduli_product();

    // Initialize the coefficients with the appropriate bit size.
    for (size_t i = 0; i < params::poly_p::degree; i++) {
        mpz_init2(coeffs[i], bits_in_moduli_product << 2);
    }

    // Sample random coefficients for the polynomial.
    for (size_t j = 0; j < params::poly_p::degree / 2; j += 32) {
        getrandom(&buf, sizeof(buf), 0);
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
void bgv_sample_message(params::poly_p & m) {
    // Similar to ghl_sample_message but for BGV scheme.
    std::array < mpz_t, DEGREE > coeffs;
    uint64_t buf;

    size_t bits_in_moduli_product = params::poly_p::bits_in_moduli_product();
    for (size_t i = 0; i < params::poly_p::degree; i++) {
        mpz_init2(coeffs[i], bits_in_moduli_product << 2);
    }

    for (size_t j = 0; j < params::poly_p::degree; j += 32) {
        getrandom(&buf, sizeof(buf), 0);
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
void bgv_rdc(params::poly_p & m) {
    std::array < mpz_t, DEGREE > coeffs;

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
void bgv_keygen(bgvkey_t & pk, params::poly_q & sk) {
    params::poly_q e = nfl::ZO_dist();  // Sample error polynomial.
    pk.a = nfl::uniform();  // Sample uniform polynomial for public key.

    sk = nfl::ZO_dist();  // Sample secret key.
    sk.ntt_pow_phi();  // Transform secret key.
    e.ntt_pow_phi();  // Transform error polynomial.

    // Compute the public key's second component.
    pk.b = pk.a * sk + (e + e + e);
}

// Function to split a secret key into multiple shares.
void bgv_keyshare(params::poly_q s[], size_t shares, params::poly_q & sk) {
    params::poly_q t = sk;
    for (size_t i = 1; i < shares; i++) {
        s[i] = nfl::uniform();  // Sample a random polynomial.
        s[i].ntt_pow_phi();  // Transform the polynomial.
        t = t - s[i];  // Subtract the share from the total.
    }
    s[0] = t;  // The first share is the remaining value.
}

// Function to encrypt a message using the BGV encryption scheme.
void bgv_encrypt(bgvenc_t & c, bgvkey_t & pk, params::poly_p & m) {
    std::array < mpz_t, DEGREE > coeffs;
    params::poly_q e1 = nfl::ZO_dist();  // Sample first error polynomial.
    params::poly_q e2 = nfl::ZO_dist();  // Sample second error polynomial.
    params::poly_q r = nfl::ZO_dist();  // Sample random polynomial.

    e1.ntt_pow_phi();  // Transform first error polynomial.
    e2.ntt_pow_phi();  // Transform second error polynomial.
    r.ntt_pow_phi();  // Transform random polynomial.

    // Compute the first component of the ciphertext.
    c.u = pk.a * r + (e1 + e1 + e1);

    // Compute the second component of the ciphertext.
    c.v = pk.b * r + (e2 + e2 + e2);

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
void bgv_add(bgvenc_t & c, bgvenc_t & d, bgvenc_t & e) {
    c.u = d.u + e.u;  // Add the first components.
    c.v = d.v + e.v;  // Add the second components.
}

// Function to decrypt a ciphertext using the BGV encryption scheme.
void bgv_decrypt(params::poly_p & m, bgvenc_t & c, params::poly_q & sk) {
    std::array < mpz_t, DEGREE > coeffs;
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

// Function to perform distributed decryption of a ciphertext.
void bgv_distdec(params::poly_q & tj, bgvenc_t & c, params::poly_q & sj) {
    std::array < mpz_t, DEGREE > coeffs;
    params::poly_q mj, Ej;  // Intermediate polynomials.
    mpz_t qDivBy2, bound;  // Variables to store half the moduli product and a bound.

    mpz_init(qDivBy2);
    mpz_init(bound);
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        mpz_init2(coeffs[i], params::poly_q::bits_in_moduli_product() << 2);
    }

    // Compute half the moduli product.
    mpz_fdiv_q_2exp(qDivBy2, params::poly_q::moduli_product(), 1);
    mpz_set_str(bound, BOUND_D, 10);  // Set the bound.

    // Sample a random polynomial.
    Ej = nfl::uniform();

    // Convert the polynomial to coefficients and center them.
    Ej.poly2mpz(coeffs);
    for (size_t i = 0; i < params::poly_q::degree; i++) {
        util::center(coeffs[i], coeffs[i], params::poly_q::moduli_product(),
                     qDivBy2);
        mpz_mod(coeffs[i], coeffs[i], bound);
    }

    // Convert the coefficients back to a polynomial and transform it.
    Ej.mpz2poly(coeffs);
    Ej.ntt_pow_phi();

    // Compute the intermediate polynomial for distributed decryption.
    mj = sj * c.u;
    tj = mj + (Ej + Ej + Ej);
}

// Function to combine the results of distributed decryption.
void bgv_comb(params::poly_p & m, bgvenc_t & c, params::poly_q t[],
              size_t shares) {
    std::array < mpz_t, DEGREE > coeffs;
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
	printf("\n** Tests for BGV encryption:\n\n");
	test();

	printf("\n** Benchmarks for BGV encryption:\n\n");
	bench();

	return 0;
}
#endif
