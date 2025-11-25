import os
import shutil
import subprocess
import sys


def check_mafft():
    if shutil.which("mafft") is None:
        msg = (
            "MAFFT not found in PATH.\n"
            "Install it (recommended):\n"
            "  conda install -c bioconda mafft\n"
            "or on Ubuntu/Debian:\n"
            "  sudo apt-get install mafft\n"
            "or download from https://mafft.cbrc.jp/\n"
        )
        raise FileNotFoundError(msg)

def run_mafft(infile, outfile, threads=None):
    check_mafft()
    if threads is None:
        try:
            threads = max(1, os.cpu_count() or 1)
        except Exception:
            threads = 1

    # sensible options:
    # --auto : automatically chooses algorithm depending on input size/complexity
    # --reorder : reorder sequences in the output to maximize alignment quality
    # --adjustdirection : auto-detect and reverse-complement sequences if needed (mostly for nucleotides; safe for proteins)
    # --thread N : use N threads
    cmd = [
        "mafft",
        "--auto",
        "--reorder",
        "--adjustdirection",
        f"--thread", str(threads),
        infile
    ]

    print("Running MAFFT (this may take seconds -> minutes depending on sequences)...")
    print("Command:", " ".join(cmd))
    with open(outfile, "w") as outfh:
        proc = subprocess.run(cmd, stdout=outfh, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        print("MAFFT failed. stderr:")
        print(proc.stderr)
        raise RuntimeError("MAFFT returned non-zero exit code.")
    print(f"MAFFT finished. Alignment saved to: {outfile}")

def trim_alignment(infile, outfile):
    """
    Trim protein alignment using trimAl automated mode.
    """
    if shutil.which("trimal") is None:
        raise FileNotFoundError("trimAl not found. Install via conda: conda install -c bioconda trimal")
    cmd = ["trimal", "-in", infile, "-out", outfile, "-automated1"]
    print("Trimming alignment with trimAl...")
    with open(outfile, "w") as outfh:
        proc = subprocess.run(cmd, stdout=outfh, stderr=subprocess.PIPE, text=True)
    if proc.returncode != 0:
        print("trimAl failed:", proc.stderr)
        raise RuntimeError("trimAl returned non-zero exit code.")
    print(f"Trimmed alignment saved to: {outfile}")


'''
def fasta_to_phylip(fasta_in=ALIGNED_FASTA, phylip_out=PHY_FILE):
    """
    Convert aligned FASTA -> Phylip (sequential)
    """
    try:
        from Bio import SeqIO
    except Exception:
        raise ImportError("Biopython is required for FASTA->PHYLIP conversion. Install with `pip install biopython` or `conda install biopython`")

    records = list(SeqIO.parse(fasta_in, "fasta"))
    if not records:
        raise RuntimeError("No sequences found in aligned FASTA")

    # Check sequence lengths (PHYLIP needs all sequences same length)
    lens = {len(r.seq) for r in records}
    if len(lens) != 1:
        raise RuntimeError("Aligned sequences have differing lengths — alignment failed or truncated.")

    # Write in sequential PHYLIP format
    SeqIO.write(records, phylip_out, "phylip-sequential")
    print(f"Wrote PHYLIP file: {phylip_out}")
'''

def main():
    if len(sys.argv) < 2:
        print("Usage: python MSA.py [test|full|anything]")
        sys.exit(1)

    prefix = sys.argv[1].lower()

    # Final input sequences from your extraction script
    name = "test" if prefix.startswith("t") else "full"

    FULL_INPUT_FASTA = f"{name}_foxo3_seq.fasta"
    DBD_INPUT_FASTA = f"{name}_foxo3_DBD.fasta"
    IDR_INPUT_FASTA = f"{name}_foxo3_IDR.fasta"

    # Intermediate alignments (MAFFT)
    FULL_ALIGNED_FASTA = f"{name}_full_seq_aligned.fasta"
    DBD_ALIGNED_FASTA  = f"{name}_DBD_aligned.fasta"
    IDR_ALIGNED_FASTA  = f"{name}_IDR_aligned.fasta"

    # Final trimmed outputs (trimAl)
    FINAL_FULL = f"{name}_full_final.fasta"
    FINAL_DBD  = f"{name}_DBD_final.fasta"
    FINAL_IDR  = f"{name}_IDR_final.fasta"

    # ------------------------------------------------------------------------------
    # MODE 1 — TEST MODE
    # Produce EVERYTHING: raw MAFFT outputs, plus trimmed versions
    # ------------------------------------------------------------------------------
    if prefix.startswith("t"):
        print("\n[TEST MODE] Producing all intermediate + final files.\n")

        run_mafft(FULL_INPUT_FASTA, FULL_ALIGNED_FASTA)
        run_mafft(DBD_INPUT_FASTA,  DBD_ALIGNED_FASTA)
        run_mafft(IDR_INPUT_FASTA,  IDR_ALIGNED_FASTA)

        trim_alignment(FULL_ALIGNED_FASTA, FINAL_FULL)
        trim_alignment(DBD_ALIGNED_FASTA,  FINAL_DBD)
        trim_alignment(IDR_ALIGNED_FASTA,  FINAL_IDR)

        print("\n[TEST MODE COMPLETE] All files written.\n")
        return

    # ------------------------------------------------------------------------------
    # MODE 2 — FULL RUN
    # Create ONLY the final trimmed outputs (no intermediate files retained)
    # ------------------------------------------------------------------------------

    print("\n[FULL MODE] Producing ONLY final trimmed alignments.\n")

    # create temporary alignment files
    tmp_full = FULL_ALIGNED_FASTA
    tmp_dbd  = DBD_ALIGNED_FASTA
    tmp_idr  = IDR_ALIGNED_FASTA

    run_mafft(FULL_INPUT_FASTA, tmp_full)
    run_mafft(DBD_INPUT_FASTA,  tmp_dbd)
    run_mafft(IDR_INPUT_FASTA,  tmp_idr)

    trim_alignment(tmp_full, FINAL_FULL)
    trim_alignment(tmp_dbd,  FINAL_DBD)
    trim_alignment(tmp_idr,  FINAL_IDR)

    # delete intermediate alignments
    for f in [tmp_full, tmp_dbd, tmp_idr]:
        try:
            os.remove(f)
        except FileNotFoundError:
            pass

    print("\n[FULL MODE COMPLETE] Only final trimmed files kept.\n")

    '''
    try:
        fasta_to_phylip()
    except Exception as e:
        print("PHYLIP conversion skipped or failed:", e)
        print("You can still use the aligned FASTA directly in many tools (MAFFT output).")
    '''

if __name__ == "__main__":
    main()


