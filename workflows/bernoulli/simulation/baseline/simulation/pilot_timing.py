"""Timing-only pilot: 20 draws per size, with no saved scientific outcomes."""
import json
import statistics
import time

from engine import backend_info, certify, generate_network, make_rng, MIXTURE_IDS


def main():
    result = {"purpose": "timing only; scientific outcomes not recorded",
              "backend": backend_info(), "mixture": "B", "replicates_per_size": 20,
              "pilot_master_seed": 2026101399, "timings": []}
    for n in (100, 1600, 25600):
        generation, certificate = [], []
        for replicate in range(20):
            rng = make_rng(2026101399, MIXTURE_IDS["B"], n, replicate)
            start = time.perf_counter()
            parents, rules, degrees, _ = generate_network(n, "B", rng)
            middle = time.perf_counter()
            certify(parents, rules, degrees, backend="cpp")
            end = time.perf_counter()
            generation.append(middle-start)
            certificate.append(end-middle)
        result["timings"].append({"n": n,
            "generation_median_s": statistics.median(generation),
            "certificate_median_s": statistics.median(certificate),
            "generation_total_s": sum(generation),
            "certificate_total_s": sum(certificate)})
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
