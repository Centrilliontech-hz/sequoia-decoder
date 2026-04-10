import gzip
import random
import os
import csv
import argparse
import subprocess
import multiprocessing
from collections import defaultdict


# =============================================================================
# 1. DNA Utils Mixin
# =============================================================================
class DNAUtils:
    BASES = ["A", "C", "G", "T"]
    BASES_DICT = {
        "N": ["A", "C", "G", "T"],
        "H": ["A", "C", "T"],
        "D": ["A", "G", "T"],
        "V": ["A", "C", "G"],
        "B": ["C", "G", "T"],
    }

    @staticmethod
    def random_dna(n):
        return "".join(random.choice(DNAUtils.BASES) for _ in range(n))

    @staticmethod
    def random_by_string(s):
        # only >3 n can be generated
        s = str(s).upper()
        assert all(
            i in DNAUtils.BASES_DICT for i in s
        ), f"Invalid character in string: {s}"
        return "".join(random.choice(DNAUtils.BASES_DICT[i]) for i in s)

    @staticmethod
    def mutate_dna(seq, n=1, allow_ins=True, allow_del=True, op_weights=None):
        if not seq:
            if n > 0:
                return "".join(random.choice(DNAUtils.BASES) for _ in range(n))
            return seq

        if op_weights is None:
            op_weights = {"sub": 1, "ins": 1, "del": 1}

        # 根据开关过滤操作并准备权重
        current_weights = op_weights.copy()
        if not allow_ins:
            current_weights["ins"] = 0
        if not allow_del:
            current_weights["del"] = 0

        ops = list(current_weights.keys())
        weights = list(current_weights.values())

        # 如果权重全为 0，则强制回退到替换操作
        if sum(weights) == 0:
            ops = ["sub"]
            weights = [1]

        res = list(seq)
        for _ in range(n):
            # 使用权重进行随机选择
            op = random.choices(ops, weights=weights, k=1)[0]

            if op == "sub" and len(res) > 0:
                idx = random.randint(0, len(res) - 1)
                orig = res[idx]
                # 确保替换后的碱基一定不同
                res[idx] = random.choice([b for b in DNAUtils.BASES if b != orig])
            elif op == "ins":
                idx = random.randint(0, len(res))
                res.insert(idx, random.choice(DNAUtils.BASES))
            elif op == "del" and len(res) > 1:
                # 保持至少1个碱基
                idx = random.randint(0, len(res) - 1)
                res.pop(idx)
            else:
                # 如果操作不满足条件（如对单碱基进行 del），回退到 sub
                if len(res) > 0:
                    idx = random.randint(0, len(res) - 1)
                    orig = res[idx]
                    res[idx] = random.choice([b for b in DNAUtils.BASES if b != orig])
                else:
                    res.append(random.choice(DNAUtils.BASES))

        # 彻底解决“突变回退”问题
        final_seq = "".join(res)
        if n > 0 and final_seq == seq:
            return DNAUtils.mutate_dna(
                seq, 1, allow_ins=allow_ins, allow_del=allow_del, op_weights=op_weights
            )

        return final_seq

    @staticmethod
    def rc(seq):
        complement = {
            "A": "T",
            "C": "G",
            "G": "C",
            "T": "A",
            "N": "N",
            "D": "H",
            "H": "D",
            "R": "Y",
            "Y": "R",
            "M": "K",
            "K": "M",
            "W": "W",
            "S": "S",
            "B": "V",
            "V": "B",
        }
        return "".join(complement.get(base, base) for base in reversed(seq))

    @staticmethod
    def apply_error_distribution(seq, allow_ins=True, allow_del=True, op_weights=None):
        if not seq:
            return seq
        r = random.random()
        if r < 0.5:
            return seq  # 0 errors
        elif r < 0.8:
            return DNAUtils.mutate_dna(
                seq, 1, allow_ins, allow_del, op_weights
            )  # 1 error
        elif r < 0.95:
            return DNAUtils.mutate_dna(
                seq, 2, allow_ins, allow_del, op_weights
            )  # 2 errors
        else:
            return DNAUtils.mutate_dna(
                seq, 3, allow_ins, allow_del, op_weights
            )  # 3 errors


# =============================================================================
# 2. Anchor Mode Mixins
# =============================================================================
class ModeHDMixin:
    @staticmethod
    def build_hd(umi, x, y, flags):
        sep = "GGG"
        if flags.get("sep_error"):
            sep = DNAUtils.mutate_dna(sep, 1)
        tail = "GG" + "T" * 20
        if flags.get("adaptor_error"):
            tail = DNAUtils.mutate_dna("GG", 1) + "T" * 20

        if flags.get("structural_error"):
            # Swap UMI and X as a structural error example
            return x + umi + sep + y + tail
        return umi + x + sep + y + tail


class ModeHDCMixin:
    @staticmethod
    def build_hdc(umi, x, y, flags):
        sep = "GGG"
        if flags.get("sep_error"):
            sep = DNAUtils.mutate_dna(sep, 1)
        tail = "CTGTCTCTTATACACAT" + "T" * 10
        if flags.get("adaptor_error"):
            tail = DNAUtils.mutate_dna("CTGTCTCTTATACACAT", 2) + "T" * 10

        if flags.get("structural_error"):
            # Swap X and Y as a structural error example
            return umi + y + sep + x + tail
        return umi + x + sep + y + tail


class ModeHDCv3Mixin:
    @staticmethod
    def build_hdcv3(umi, x, y, flags):
        nnnhh = (
            DNAUtils.random_dna(3)
            + random.choice(["A", "T", "C"])
            + random.choice(["A", "T", "C"])
        )
        hhnnn = (
            random.choice(["A", "T", "C"])
            + random.choice(["A", "T", "C"])
            + DNAUtils.random_dna(3)
        )
        sep = "GGG"
        if flags.get("sep_error"):
            sep = DNAUtils.mutate_dna(sep, 1)

        prefix = "GGTGCCTGGTGTATCCATCCAA"
        if flags.get("structural_error"):
            prefix = DNAUtils.mutate_dna(prefix, 3)

        suffix = "CTGTCTCTTATACACA"
        if flags.get("adaptor_error"):
            suffix = DNAUtils.mutate_dna(suffix, 2)

        return prefix + nnnhh + "GG" + x + sep + y + "GG" + hhnnn + suffix


class ModeHDCTCRMixin:
    @staticmethod
    def build_hdc_tcr(umi, x, y, flags):
        sep = "CCC"
        if flags.get("sep_error"):
            sep = DNAUtils.mutate_dna(sep, 1)
        suffix = "AGATCGGAAGAGCGTCGTGTAG"
        if flags.get("adaptor_error"):
            suffix = DNAUtils.mutate_dna(suffix, 3)

        # {ZipcodeY‘:15 或 16}CCC{ZipcodeX':15或16}{UMI:N9}AGATCGGAAGAGCGTCGTGTAG
        if flags.get("structural_error"):
            # Swap X and Y
            return DNAUtils.rc(x) + sep + DNAUtils.rc(y) + DNAUtils.rc(umi) + suffix
        return DNAUtils.rc(y) + sep + DNAUtils.rc(x) + DNAUtils.rc(umi) + suffix


class ModeHDCv3TCRMixin:
    @staticmethod
    def build_hdcv3_tcr(umi, x, y, flags):
        sep = "CCC"
        if flags.get("sep_error"):
            sep = DNAUtils.mutate_dna(sep, 1)
        nnnddcc = (
            DNAUtils.random_dna(3)
            + random.choice(["A", "G", "T"])
            + random.choice(["A", "G", "T"])
            + "CC"
        )
        ccddnnn = (
            "CC"
            + random.choice(["A", "G", "T"])
            + random.choice(["A", "G", "T"])
            + DNAUtils.random_dna(3)
        )

        suffix = "TTGGATGGATACACCAGGCACC"
        if flags.get("adaptor_error"):
            suffix = DNAUtils.mutate_dna(suffix, 3)

        r1 = nnnddcc + DNAUtils.rc(y) + sep + DNAUtils.rc(x) + ccddnnn + suffix
        if flags.get("structural_error"):
            r1 = r1[::-1]  # reverse the whole thing
        return r1


# =============================================================================
# 3. Scenario Perturbation Mixins
# =============================================================================
class ScenarioPerturbationMixin:
    @staticmethod
    def apply_shift(r1):
        shift = random.choice([-2, -1, 1, 2])
        if shift > 0:
            return DNAUtils.random_dna(shift) + r1
        elif shift < 0:
            return r1[abs(shift) :]
        return r1

    @staticmethod
    def apply_interference(r1):
        # Inject random noise at a random position
        if len(r1) < 10:
            return r1
        pos = random.randint(0, len(r1) - 5)
        return r1[:pos] + DNAUtils.random_dna(3) + r1[pos + 3 :]


# =============================================================================
# 4. Main Test Generator Class
# =============================================================================
class TestGenerator(
    DNAUtils,
    ModeHDMixin,
    ModeHDCMixin,
    ModeHDCv3Mixin,
    ModeHDCTCRMixin,
    ModeHDCv3TCRMixin,
    ScenarioPerturbationMixin,
):

    def __init__(self, codebook_path):
        self.codebook = []
        with open(codebook_path, "r") as f:
            for line in f:
                if line.strip() and not line.startswith("#"):
                    parts = line.strip().split(",")
                    if len(parts) == 2:
                        self.codebook.append((parts[0], parts[1]))

        self.ground_truth = {}  # seq_id -> (expected_x, expected_y, umi, scenario)

    @staticmethod
    def generate_r1(mode, umi, x_seq, y_seq, scenario):
        flags = {
            "sep_error": scenario == "sep_error",
            "adaptor_error": scenario == "adaptor_error",
            "structural_error": scenario == "structural_error",
        }

        # 1. Build Base Structure based on Mode
        if mode == "HD":
            r1 = TestGenerator.build_hd(umi, x_seq, y_seq, flags)
        elif mode == "HDC":
            r1 = TestGenerator.build_hdc(umi, x_seq, y_seq, flags)
        elif mode == "HDCv3":
            r1 = TestGenerator.build_hdcv3(umi, x_seq, y_seq, flags)
        elif mode == "HDC-TCR":
            r1 = TestGenerator.build_hdc_tcr(umi, x_seq, y_seq, flags)
        elif mode == "HDCv3-TCR":
            r1 = TestGenerator.build_hdcv3_tcr(umi, x_seq, y_seq, flags)
        else:
            raise ValueError(f"Unsupported mode: {mode}")

        # 2. Apply Scenario Perturbations
        if scenario == "interference":
            r1 = TestGenerator.apply_interference(r1)
        elif scenario == "shift":
            r1 = TestGenerator.apply_shift(r1)

        return r1

    def generate_test_set(self, mode, n_per_scenario=100, output_prefix="test"):
        r1_path = f"{output_prefix}.R1.fastq.gz"
        r2_path = f"{output_prefix}.R2.fastq.gz"

        scenarios = [
            "normal",
            "ambiguity",
            "interference",
            "shift",
            "sep_error",
            "adaptor_error",
            "structural_error",
        ]

        # Prepare arguments for parallel processing
        tasks = []
        for idx, scenario in enumerate(scenarios):
            tasks.append(
                (
                    scenario,
                    n_per_scenario,
                    idx * n_per_scenario,  # offset for sequence ID
                    mode,
                    self.codebook,
                )
            )

        print(
            f"[*] Parallel generating {len(scenarios)} scenarios using {multiprocessing.cpu_count()} cores..."
        )

        # Use multiprocessing to generate scenarios in parallel
        # compresslevel=1 is much faster than default 9
        with multiprocessing.Pool() as pool:
            results = pool.map(self._generate_scenario_worker, tasks)

        with gzip.open(r1_path, "wt", compresslevel=1) as f1, gzip.open(
            r2_path, "wt", compresslevel=1
        ) as f2:

            for r1_data, r2_data, gt_data in results:
                f1.write(r1_data)
                f2.write(r2_data)
                self.ground_truth.update(gt_data)

        return r1_path, r2_path

    @staticmethod
    def _generate_scenario_worker(args):
        """Worker function for parallel generation."""
        scenario, n_per_scenario, start_count, mode, codebook = args
        # Re-seed random for each process to ensure different sequences
        random.seed()

        r1_lines = []
        r2_lines = []
        gt_part = {}

        n_codes = len(codebook)

        for i in range(n_per_scenario):
            count = start_count + i + 1
            seq_id = f"read_{count}_{scenario}"

            if i < n_per_scenario // 2:
                idx_x = i % n_codes
                idx_y = random.randint(0, n_codes - 1)
            else:
                idx_x = random.randint(0, n_codes - 1)
                idx_y = i % n_codes

            code_x, seq_x = codebook[idx_x]
            code_y, seq_y = codebook[idx_y]

            if scenario == "ambiguity":
                seq_x_gen = DNAUtils.mutate_dna(seq_x, 1)
                seq_y_gen = DNAUtils.mutate_dna(seq_y, 1)
            elif scenario == "normal":
                seq_x_gen = seq_x
                seq_y_gen = seq_y
            else:
                seq_x_gen = DNAUtils.apply_error_distribution(seq_x)
                seq_y_gen = DNAUtils.apply_error_distribution(seq_y)

            umi = DNAUtils.random_dna(9)
            r1_seq = TestGenerator.generate_r1(
                mode, umi, seq_x_gen, seq_y_gen, scenario
            )
            r2_seq = DNAUtils.random_dna(100)

            r1_lines.append(f"@{seq_id}\n{r1_seq}\n+\n{'I'*len(r1_seq)}\n")
            r2_lines.append(f"@{seq_id}\n{r2_seq}\n+\n{'I'*len(r2_seq)}\n")

            expected_x = code_x
            expected_y = code_y
            if scenario in ["structural_error"]:
                expected_x = None
                expected_y = None
            gt_part[seq_id] = (expected_x, expected_y, umi, scenario, r1_seq)

        return "".join(r1_lines), "".join(r2_lines), gt_part

    def evaluate(self, decoder_output_csv, error_csv_path=None):
        """Evaluate decoding results against ground truth."""
        print("[*] Evaluating results...")

        # Scenarios result summary
        results = defaultdict(
            lambda: {
                "total": 0,
                "decoded": 0,
                "x_correct": 0,
                "y_correct": 0,
                "x_dists": [],
                "y_dists": [],
                "correct": 0,  # Both X and Y correct
                "dist_sum": 0,
                "dist_max": 0,
            }
        )
        # 分别统计 X 和 Y 的情况
        # total: 总数, decoded: 成功解码数, errors: 误判数(解码但错误), dist_sum: 误判的累积距离
        scenario_x_stats = defaultdict(
            lambda: defaultdict(
                lambda: {"total": 0, "decoded": 0, "errors": 0, "dist_sum": 0}
            )
        )
        scenario_y_stats = defaultdict(
            lambda: defaultdict(
                lambda: {"total": 0, "decoded": 0, "errors": 0, "dist_sum": 0}
            )
        )

        codebook_map = {code: seq for code, seq in self.codebook}

        # Load decoded results
        decoded = {}
        if os.path.exists(decoder_output_csv):
            with gzip.open(decoder_output_csv, "rt") as f:
                reader = csv.DictReader(f)
                for row in reader:
                    decoded[row["SeqID"]] = {
                        "code1": row["X-Pos"],
                        "code2": row["Y-Pos"],
                        "umi": row["UMI"],
                    }

        error_records = []

        for seq_id, gt_info in self.ground_truth.items():
            # Support both old (4 elements) and new (5 elements) ground truth format for compatibility
            if len(gt_info) == 5:
                exp_x, exp_y, exp_umi, scenario, r1_seq = gt_info
            else:
                exp_x, exp_y, exp_umi, scenario = gt_info
                r1_seq = "N/A"

            res = results[scenario]
            res["total"] += 1

            if exp_x is not None:
                scenario_x_stats[scenario][exp_x]["total"] += 1
            if exp_y is not None:
                scenario_y_stats[scenario][exp_y]["total"] += 1

            got = decoded.get(seq_id)

            def get_code_num(c):
                if not c or c == "None":
                    return None
                return int(c.split("→")[1]) if "→" in c else int(c)

            exp_x_num = get_code_num(exp_x)
            exp_y_num = get_code_num(exp_y)

            got_x_num = None
            got_y_num = None
            if got:
                res["decoded"] += 1
                got_x_num = get_code_num(got["code1"])
                got_y_num = get_code_num(got["code2"])

            # 处理 X 端
            dist_x = 0
            if exp_x_num is not None:
                if got:
                    scenario_x_stats[scenario][exp_x]["decoded"] += 1
                    if got_x_num is not None:
                        dist_x = abs(got_x_num - exp_x_num)
                        res["x_dists"].append(dist_x)
                        if dist_x > 0:
                            scenario_x_stats[scenario][exp_x]["errors"] += 1
                            scenario_x_stats[scenario][exp_x]["dist_sum"] += dist_x
                        else:
                            res["x_correct"] += 1
                else:
                    # Not decoded but expected
                    dist_x = 999

            # 处理 Y 端
            dist_y = 0
            if exp_y_num is not None:
                if got:
                    scenario_y_stats[scenario][exp_y]["decoded"] += 1
                    if got_y_num is not None:
                        dist_y = abs(got_y_num - exp_y_num)
                        res["y_dists"].append(dist_y)
                        if dist_y > 0:
                            scenario_y_stats[scenario][exp_y]["errors"] += 1
                            scenario_y_stats[scenario][exp_y]["dist_sum"] += dist_y
                        else:
                            res["y_correct"] += 1
                else:
                    # Not decoded but expected
                    dist_y = 999

            # 记录错误详情
            # 仅记录显著错误（差距大于 4），忽略微小的歧义误差
            if dist_x > 4 or dist_y > 4:
                error_records.append(
                    {
                        "SeqID": seq_id,
                        "Scenario": scenario,
                        "Truth_X_Code": exp_x,
                        "Truth_X_Seq": codebook_map.get(exp_x, "N/A"),
                        "Truth_Y_Code": exp_y,
                        "Truth_Y_Seq": codebook_map.get(exp_y, "N/A"),
                        "Decoded_X_Code": got["code1"] if got else "None",
                        "Decoded_X_Seq": (
                            codebook_map.get(got["code1"], "N/A")
                            if got and got["code1"] != "None"
                            else "N/A"
                        ),
                        "Decoded_Y_Code": got["code2"] if got else "None",
                        "Decoded_Y_Seq": (
                            codebook_map.get(got["code2"], "N/A")
                            if got and got["code2"] != "None"
                            else "N/A"
                        ),
                        "Dist_X": dist_x,
                        "Dist_Y": dist_y,
                        "Sequence": r1_seq,
                    }
                )

            # 整体准确率统计
            if got and exp_x_num is not None and exp_y_num is not None:
                if dist_x == 0 and dist_y == 0:
                    res["correct"] += 1

                dist = dist_x + dist_y
                res["dist_sum"] += dist
                res["dist_max"] = max(res["dist_max"], dist)

        # Write error records to CSV
        if error_csv_path and error_records:
            print(
                f"[*] Saving {len(error_records)} error details to {error_csv_path}..."
            )
            keys = error_records[0].keys()
            with open(error_csv_path, "w", newline="") as f:
                dict_writer = csv.DictWriter(f, fieldnames=keys)
                dict_writer.writeheader()
                dict_writer.writerows(error_records)

        # Calculate Top 10 error codes PER scenario for X and Y based on AVG DISTANCE
        def get_top_10(stats_dict):
            top_10_map = {}
            for scenario, stats_map in stats_dict.items():
                code_metrics = []
                for code, stats in stats_map.items():
                    if stats["total"] > 0:
                        decode_rate = stats["decoded"] / stats["total"]
                        avg_dist = (
                            stats["dist_sum"] / stats["errors"]
                            if stats["errors"] > 0
                            else 0
                        )
                        # 排序权重：优先看平均距离，其次看错误数
                        if avg_dist > 0 or decode_rate < 1.0:
                            code_metrics.append(
                                (
                                    code,
                                    codebook_map.get(code, "N/A"),
                                    decode_rate,
                                    avg_dist,
                                    stats["errors"],
                                    stats["decoded"],
                                )
                            )
                # 按平均距离降序排列
                top_10_map[scenario] = sorted(
                    code_metrics, key=lambda x: x[3], reverse=True
                )[:10]
            return top_10_map

        scenario_top_10_x = get_top_10(scenario_x_stats)
        scenario_top_10_y = get_top_10(scenario_y_stats)

        return results, scenario_top_10_x, scenario_top_10_y


# =============================================================================
# 5. Execution Entry Point
# =============================================================================
def main():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--mode", choices=["HD", "HDC", "HDCv3", "HDC-TCR", "HDCv3-TCR"], default="HD"
    )
    parser.add_argument("--n", type=int, default=100)
    parser.add_argument("--decoder", default="./target/release/sequoia-decoder")
    parser.add_argument(
        "--outdir", default="tests", help="Base directory for all test outputs"
    )
    args = parser.parse_args()

    # Setup directories
    data_dir = os.path.join(args.outdir, "data")
    report_dir = os.path.join(args.outdir, "reports")
    for d in [data_dir, report_dir]:
        if not os.path.exists(d):
            os.makedirs(d)

    # Locate codebook
    script_dir = os.path.dirname(os.path.abspath(__file__))
    codebook_path = os.path.join(script_dir, "../../src/codebook.csv")
    if not os.path.exists(codebook_path):
        codebook_path = os.path.join(script_dir, "src/codebook.csv")

    gen = TestGenerator(codebook_path)

    print(f"[*] Mode: {args.mode}")
    print(f"[*] Generating {args.n} reads per scenario...")
    fastq_prefix = os.path.join(data_dir, f"test_{args.mode}")
    r1, r2 = gen.generate_test_set(
        args.mode, n_per_scenario=args.n, output_prefix=fastq_prefix
    )

    run_dir = os.path.join(args.outdir, f"run_{args.mode}")
    if not os.path.exists(run_dir):
        os.makedirs(run_dir)

    print(f"[*] Running decoder...")
    cmd = [args.decoder, "-i", r1, r2, "-o", run_dir, "--anchor", args.mode]
    cmd += ["--sep", "CCC"] if "TCR" in args.mode else ["--sep", "GGG"]
    subprocess.run(cmd, check=True)

    print(f"[*] Evaluating results...")
    decoded_csv = os.path.join(run_dir, "decoded_zipcodes.csv.gz")
    error_csv = os.path.join(report_dir, f"errors_{args.mode}.csv")
    results, scenario_top_10_x, scenario_top_10_y = gen.evaluate(
        decoded_csv, error_csv_path=error_csv
    )

    # Generate Report
    report_lines = [
        f"# Test Report for {args.mode}",
        "\n## Accuracy Summary",
        "\n| Scenario | Decode % | X-Corr % | Y-Corr % | X-Mean | Y-Mean | X-Q90 | Y-Q90 | X-Max | Y-Max |",
        "| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |",
    ]

    def get_q90(data):
        if not data:
            return 0
        sorted_data = sorted(data)
        idx = int(len(sorted_data) * 0.9)
        return sorted_data[min(idx, len(sorted_data) - 1)]

    for scenario, res in results.items():
        dec_rate = res["decoded"] / res["total"]
        x_corr = res["x_correct"] / res["decoded"] if res["decoded"] > 0 else 0
        y_corr = res["y_correct"] / res["decoded"] if res["decoded"] > 0 else 0

        x_mean = sum(res["x_dists"]) / len(res["x_dists"]) if res["x_dists"] else 0
        y_mean = sum(res["y_dists"]) / len(res["y_dists"]) if res["y_dists"] else 0

        x_q90 = get_q90(res["x_dists"])
        y_q90 = get_q90(res["y_dists"])

        x_max = max(res["x_dists"]) if res["x_dists"] else 0
        y_max = max(res["y_dists"]) if res["y_dists"] else 0

        report_lines.append(
            f"| {scenario} | {dec_rate:.1%} | {x_corr:.1%} | {y_corr:.1%} | {x_mean:.2f} | {y_mean:.2f} | {x_q90} | {y_q90} | {x_max} | {y_max} |"
        )

    report_lines.append(
        "\n## Code Dependency Analysis (Top 10 Abnormal Codes per Scenario)"
    )
    for scenario in [
        "normal",
        "ambiguity",
        "interference",
        "shift",
        "sep_error",
        "adaptor_error",
        "structural_error",
    ]:
        top_x = scenario_top_10_x.get(scenario, [])
        top_y = scenario_top_10_y.get(scenario, [])

        if top_x or top_y:
            report_lines.append(f"\n### Scenario: {scenario}")
            # 并列表格头
            report_lines.append(
                "\n| X-Code | X-Seq | X-Decode % | X-Err-Dist | Y-Code | Y-Seq | Y-Decode % | Y-Err-Dist |"
            )
            report_lines.append(
                "| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |"
            )

            # 取两者中较长的长度
            max_len = max(len(top_x), len(top_y))
            for i in range(max_len):
                row = []
                # X 部分
                if i < len(top_x):
                    code, seq, dec_rate, avg_dist, errs, dec_total = top_x[i]
                    row.extend(
                        [
                            code,
                            f"`{seq}`",
                            f"{dec_rate:.1%}",
                            f"{avg_dist:.1f} ({errs}/{dec_total})",
                        ]
                    )
                else:
                    row.extend(["-", "-", "-", "-"])

                # Y 部分
                if i < len(top_y):
                    code, seq, dec_rate, avg_dist, errs, dec_total = top_y[i]
                    row.extend(
                        [
                            code,
                            f"`{seq}`",
                            f"{dec_rate:.1%}",
                            f"{avg_dist:.1f} ({errs}/{dec_total})",
                        ]
                    )
                else:
                    row.extend(["-", "-", "-", "-"])

                report_lines.append(f"| {' | '.join(row)} |")

    report_path = os.path.join(report_dir, f"report_{args.mode}.md")
    with open(report_path, "w") as f:
        f.write("\n".join(report_lines))

    print(f"[+] Report saved: {report_path}")
    if os.path.exists(error_csv):
        print(f"[+] Error details saved: {error_csv}")


if __name__ == "__main__":
    main()
