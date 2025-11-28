import os
import shutil
import subprocess


ORDERS = [
    "Artiodactyla",
    "Carnivora",
    "Chiroptera",
    "Primates",
    "Rodentia",
]

def check_mafft():
    if shutil.which("mafft") is None:
        raise FileNotFoundError(
            "MAFFT not found. Install with:\n"
            "  conda install -c bioconda mafft"
        )

def run_mafft_dbd(infile, outfile, threads=None):
    """
    Run MAFFT L-INS-i specifically for DBD sequences.
    """
    check_mafft()
    if threads is None:
        try:
            threads = max(1, os.cpu_count() or 1)
        except Exception:
            threads = 1

    cmd = [
        "mafft",
        "--localpair",
        "--maxiterate", "1000",
        "--reorder",
        "--adjustdirection",
        "--thread", str(threads),
        infile
    ]

    print(f"\nRunning MAFFT on {infile}")
    print("Command:", " ".join(cmd))

    with open(outfile, "w") as outfh:
        proc = subprocess.run(cmd, stdout=outfh, stderr=subprocess.PIPE, text=True)

    if proc.returncode != 0:
        print("MAFFT failed:")
        print(proc.stderr)
        raise RuntimeError("MAFFT returned non-zero exit code.")

    print(f"MAFFT complete → {outfile}")

def trim_alignment(infile, outfile):
    """
    Trim protein alignments with trimAl automated mode.
    """
    if shutil.which("trimal") is None:
        raise FileNotFoundError(
            "trimAl not found. Install with:\n"
            "  conda install -c bioconda trimal"
        )

    cmd = ["trimal", "-in", infile, "-out", outfile, "-automated1"]
    print(f"Trimming {infile} with trimAl...")

    with open(outfile, "w") as outfh:
        proc = subprocess.run(cmd, stdout=outfh, stderr=subprocess.PIPE, text=True)

    if proc.returncode != 0:
        print("trimAl failed:")
        print(proc.stderr)
        raise RuntimeError("trimAl returned non-zero exit code.")

    print(f"Trimmed alignment → {outfile}")

def main():
    print("\n=== FOXO3 DBD MSA PIPELINE ===")

    for order in ORDERS:
        input_fasta = f"{order}_foxo3_seq.fasta"
        aligned_out = f"{order}_DBD_aligned.fasta"
        trimmed_out = f"{order}_DBD_trimmed.fasta"

        if not os.path.exists(input_fasta):
            print(f"\n[SKIP] {input_fasta} not found.")
            continue

        print(f"\n=== Processing: {order} ===")

        # Step 1: MAFFT
        run_mafft_dbd(input_fasta, aligned_out)

        # Step 2: trimAl
        trim_alignment(aligned_out, trimmed_out)

        # Remove alignment intermediate
        try:
            os.remove(aligned_out)
        except FileNotFoundError:
            pass

    print("\nAll done.\n")


if __name__ == "__main__":
    main()
