import subprocess
import os
import shutil
import argparse
from datetime import datetime


def run_command(cmd):
    print(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True)


def main():
    parser = argparse.ArgumentParser(
        description="Run all decoder test modes and merge reports."
    )
    parser.add_argument(
        "-n", type=int, default=50, help="Number of sequences per scenario"
    )
    parser.add_argument(
        "--keep-data", action="store_true", help="Keep temporary fastq and run data"
    )
    args = parser.parse_args()

    modes = ["HD", "HDC", "HDCv3", "HDC-TCR", "HDCv3-TCR"]
    base_dir = "tests"
    temp_run_dir = os.path.join(base_dir, "batch_run_tmp")
    report_dir = os.path.join(base_dir, "reports")
    script_path = os.path.join(base_dir, "scripts/test_generator.py")

    # Ensure directories exist
    if os.path.exists(temp_run_dir):
        shutil.rmtree(temp_run_dir)
    os.makedirs(temp_run_dir)
    os.makedirs(report_dir, exist_ok=True)

    final_report_path = "full_test_report.md"
    with open(final_report_path, "w") as f_out:
        f_out.write(f"# Full Decoder Test Report\n")
        f_out.write(f"Generated at: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")

    for mode in modes:
        print(f"\n>>> Testing Mode: {mode}")
        try:
            # Run the generator for this mode
            cmd = [
                "python3",
                script_path,
                "--mode",
                mode,
                "--n",
                str(args.n),
                "--outdir",
                temp_run_dir,
            ]
            run_command(cmd)

            # Append mode report to final report
            mode_report_path = os.path.join(temp_run_dir, f"reports/report_{mode}.md")
            if os.path.exists(mode_report_path):
                with open(mode_report_path, "r") as f_in:
                    content = f_in.read()
                    with open(final_report_path, "a") as f_out:
                        f_out.write(content + "\n\n---\n\n")

                # Copy to permanent reports dir as well
                shutil.copy(
                    mode_report_path, os.path.join(report_dir, f"report_{mode}.md")
                )
        except Exception as e:
            print(f"Error testing mode {mode}: {e}")

    # Cleanup
    if not args.keep_data:
        print("\nCleaning up temporary files...")
        shutil.rmtree(temp_run_dir)
    else:
        print(f"\nTemporary files kept in: {temp_run_dir}")

    print(f"\nAll tests completed. Final report: {final_report_path}")


if __name__ == "__main__":
    main()
