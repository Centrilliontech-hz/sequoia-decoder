import os
import json
import csv
import gzip
import subprocess
import argparse
import sys
from collections import defaultdict

# Add the current directory to sys.path so we can import test_generator
sys.path.append(os.path.dirname(os.path.abspath(__file__)))
try:
    from test_generator import TestGenerator
except ImportError:
    # Fallback if the script is run from a different location
    sys.path.append(
        os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts")
    )
    from test_generator import TestGenerator


def run_compare():
    parser = argparse.ArgumentParser(
        description="Compare sequoia-decoder and PostMaster52"
    )
    parser.add_argument(
        "--n", type=int, default=100, help="Number of reads per scenario"
    )
    parser.add_argument(
        "--decoder",
        default="./target/release/sequoia-decoder",
        help="Path to sequoia-decoder",
    )
    parser.add_argument(
        "--postmaster",
        default="./../sequoia-workflow/sequoia_scripts/bin/PostMaster52",
        help="Path to PostMaster52",
    )
    parser.add_argument(
        "--outdir", default="tests/compare_results", help="Output directory"
    )
    parser.add_argument(
        "--threads", type=int, default=32, help="Number of threads to use"
    )
    args = parser.parse_args()

    # 1. Setup paths
    abs_outdir = os.path.abspath(args.outdir)
    if not os.path.exists(abs_outdir):
        os.makedirs(abs_outdir)

    data_dir = os.path.join(abs_outdir, "data")
    if not os.path.exists(data_dir):
        os.makedirs(data_dir)

    script_dir = os.path.dirname(os.path.abspath(__file__))
    codebook_path = os.path.join(script_dir, "../../src/codebook.csv")
    if not os.path.exists(codebook_path):
        # Alternative path check
        codebook_path = os.path.abspath("src/codebook.csv")

    # Check if tools exist
    if not os.path.exists(args.decoder):
        print(f"[!] Error: sequoia-decoder not found at {args.decoder}")
        sys.exit(1)
    if not os.path.exists(args.postmaster):
        print(f"[!] Error: PostMaster52 not found at {args.postmaster}")
        sys.exit(1)

    # 2. Generate Test Data
    print("[*] Generating test data...")
    gen = TestGenerator(codebook_path)
    fastq_prefix = os.path.join(data_dir, "test_HD")
    r1, r2 = gen.generate_test_set(
        "HD", n_per_scenario=args.n, output_prefix=fastq_prefix
    )

    # 3. Run sequoia-decoder
    print("[*] Running sequoia-decoder...")
    sequoia_out = os.path.join(abs_outdir, "sequoia_run")
    if not os.path.exists(sequoia_out):
        os.makedirs(sequoia_out)

    cmd_sequoia = [
        args.decoder,
        "-i",
        r1,
        r2,
        "-o",
        sequoia_out,
        "--anchor",
        "HD",
        "-j",
        str(args.threads),
        "--fuzzy-sep",
        "--strategy",
        "min-total-edit",
        "--output-mode",
        "merge",
    ]
    print(f"    Command: {' '.join(cmd_sequoia)}")
    subprocess.run(cmd_sequoia, check=True)

    # 4. Run PostMaster52
    print("[*] Running PostMaster52...")
    pm_out = os.path.join(abs_outdir, "postmaster_run")
    if not os.path.exists(pm_out):
        os.makedirs(pm_out)

    # Create parameters.json for PM52 based on the user's template
    pm_params = {
        "BatchSize": 500000,
        "BottomAdapt": "",
        "CodeFormat": "Complement",
        "FastqFileType": "PairEnd",
        "InputFastqFiles": [os.path.abspath(r1), os.path.abspath(r2)],
        "MAX_Separator_Edit_Distance": 1,
        "MaxEditDistance": 3,
        "NumberOfThreads": args.threads,
        "OutputRead2": True,
        "PolyTCaptureSequence": True,
        "Read1UMILength": 9,
        "SeparatorSeq": "GGG",
        "ShowProgress": True,
        "ZipcodeType": "Yosemite2017",
        "TopAdaptor": "GGTTTTTT",
        "OutputDir": os.path.abspath(pm_out) + "/",
    }

    params_path = os.path.join(abs_outdir, "pm_params.json")
    with open(params_path, "w") as f:
        json.dump(pm_params, f, indent=4)

    cmd_pm = [os.path.abspath(args.postmaster), params_path]
    print(f"    Command: {' '.join(cmd_pm)}")
    subprocess.run(cmd_pm, check=True)

    # 5. Evaluate and Compare
    print("[*] Evaluating results...")
    sequoia_csv = os.path.join(sequoia_out, "decoded_zipcodes.csv.gz")
    pm_csv = os.path.join(pm_out, "decoded_zipcodes.csv.gz")

    res_sequoia, _, _ = gen.evaluate(sequoia_csv)
    res_pm, _, _ = gen.evaluate(pm_csv)

    # 6. Generate Comparison Report
    report_path = os.path.join(abs_outdir, "comparison_report.md")
    print(f"[*] Generating comparison report at {report_path}...")

    with open(report_path, "w") as f:
        f.write("# Decoder Comparison Report (HD Mode)\n\n")
        f.write(
            "This report compares the performance of `sequoia-decoder` and `PostMaster52` using generated test data.\n\n"
        )

        f.write("## Accuracy Summary\n\n")
        f.write(
            "| Scenario | Decoder | Total | Decoded | Decode % | Correct % | Mean Dist |\n"
        )
        f.write("| :--- | :--- | :--- | :--- | :--- | :--- | :--- |\n")

        scenarios = sorted(res_sequoia.keys())
        for scenario in scenarios:
            s_res = res_sequoia[scenario]
            p_res = res_pm[scenario]

            s_dec_rate = s_res["decoded"] / s_res["total"]
            s_corr_rate = (
                s_res["correct"] / s_res["decoded"] if s_res["decoded"] > 0 else 0
            )
            s_mean_dist = (
                (sum(s_res["x_dists"]) + sum(s_res["y_dists"]))
                / (len(s_res["x_dists"]) + len(s_res["y_dists"]))
                if (len(s_res["x_dists"]) + len(s_res["y_dists"])) > 0
                else 0
            )

            p_dec_rate = p_res["decoded"] / p_res["total"]
            p_corr_rate = (
                p_res["correct"] / p_res["decoded"] if p_res["decoded"] > 0 else 0
            )
            p_mean_dist = (
                (sum(p_res["x_dists"]) + sum(p_res["y_dists"]))
                / (len(p_res["x_dists"]) + len(p_res["y_dists"]))
                if (len(p_res["x_dists"]) + len(p_res["y_dists"])) > 0
                else 0
            )

            f.write(
                f"| {scenario} | Sequoia | {s_res['total']} | {s_res['decoded']} | {s_dec_rate:.1%} | {s_corr_rate:.1%} | {s_mean_dist:.2f} |\n"
            )
            f.write(
                f"| | PM52 | {p_res['total']} | {p_res['decoded']} | {p_dec_rate:.1%} | {p_corr_rate:.1%} | {p_mean_dist:.2f} |\n"
            )
            f.write("| --- | --- | --- | --- | --- | --- | --- |\n")

    print(f"[+] Comparison report saved to: {report_path}")


if __name__ == "__main__":
    run_compare()
