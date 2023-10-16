NTRU = {
    "ntru_keygen": 172402073,
    "ntru_sample_message": 236925,
    "ntru_enc": 596368,
    "ntru_dec": 628323,
    "ntru_distdec": 807781,
    "ntru_comb": 494008,
}

BDLOP = {
    "bdlop_keygen": 117564,
    "bdlop_sample_chal": 163698,
    "bdlop_sample_rand": 413910,
    "bdlop_commit": 204523,
    "bdlop_open": 2551079,
}

PISMALL = {
    "pismall_setup": 4570786660,
    "pismall_prover": 183613860009,
    "pismall_verifier": 9895237671,
}

SHUFFLE = {
    "shuffle_linear_hash": 430810,
    "shuffle_linear_proof": 51794475,
    "shuffle_linear_verifier": 5842662,
    "shuffle_proof_and_verify": 60409852302,
    "shuffle_proof": 54084003493,
}

PROCESSOR = {
    "op_modes": [32, 64],
    "addr_size": {"physical": 39, "virtual": 48},
    "byte_order": "Little Endian",
    "model_name": "11th Gen Inter(R) Core(TM) i5-1145G7 @ 2.60GHz",
    # "hertz": 2600000000,
    # "hertz": 1500000000,
    "hertz": 3600000000,
    "turbo": False,
}


def in_ms(value, freq):
    return "{:.3f}".format(value / (freq / 1000))


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
