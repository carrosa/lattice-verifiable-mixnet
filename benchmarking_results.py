NTRU = {
    "ntru_keygen": 172402073,
    "ntru_sample_message": 236925,
    "ntru_enc": 596368,
    "ntru_dec": 628323,
    "ntru_distdec": 807781,
    "ntru_comb": 494008,
}

BDLOP = {
    "bdlop_keygen": 40208,
    "bdlop_sample_chal": 55356,
    "bdlop_sample_rand": 141505,
    "bdlop_commit": 69991,
    "bdlop_open": 872940,
}

PISMALL = {
    "pismall_setup": 1687136462,
    "pismall_prover": 67074041350,
    "pismall_verifier": 3408190269,
}

SHUFFLE = {
    "shuffle_linear_hash": 152496,
    "shuffle_linear_proof": 18287390,
    "shuffle_linear_verifier": 2122465,
    "shuffle_proof": 21493411889,
}

PROCESSOR = {
    "op_modes": [32, 64],
    "addr_size": {"physical": 39, "virtual": 48},
    "byte_order": "Little Endian",
    "model_name": "11th Gen Inter(R) Core(TM) i5-1145G7 @ 2.60GHz",
    "hertz": 2600000000,
    # "hertz": 1500000000,
    "turbo": False,
}


def in_ms(value, freq):
    return "{:.20f}".format(value / (freq / 1000))


def print_times():
    for key, value in NTRU.items():
        print(key, in_ms(value, PROCESSOR["hertz"]), "ms")

    for key, value in BDLOP.items():
        print(key, in_ms(value, PROCESSOR["hertz"]), "ms")

    for key, value in PISMALL.items():
        print(key, in_ms(value, PROCESSOR["hertz"]), "ms")

    for key, value in SHUFFLE.items():
        print(key, in_ms(value, PROCESSOR["hertz"]), "ms")


print_times()
