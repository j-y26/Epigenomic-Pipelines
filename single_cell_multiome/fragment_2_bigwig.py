# Convert fragment files to BigWig format using bedtools, samtools, and bamCoverage.

# Usage: python fragment_2_bigwig.py -i/--input <input_bed> \
#                                    -g/--genome <genome_sizes> \
#                                    -o/--output <output_bigwig> \
#                                    [--binsize <bin_size>] \
#                                    [--normalize <CPM|RPKM|BPM|RPGC|None>] \
#                                    [--genome_size <effective_genome_size>] \
#                                    [--threads <num_threads>] \
#                                    [--barcodes <barcode_file>] \
#                                    [--keep-temp]


import argparse
import gzip
import os
import shutil
import subprocess
import tempfile
from pathlib import Path
from tqdm import tqdm


def parse_args():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(description="Convert fragment file to BigWig.")
    parser.add_argument("-i", "--input", required=True, help="Input BED or fragment.tsv(.gz) file")
    parser.add_argument("-g", "--genome", required=True, help="Genome chromosome sizes file, available at UCSC")
    parser.add_argument("-o", "--output", required=True, help="Output BigWig file")
    parser.add_argument("--binsize", type=int, default=10, help="Bin size (default: 10)")
    parser.add_argument("--normalize", default="RPGC", choices=["CPM", "RPKM", "BPM", "RPGC", "None"],
                        help="Normalization method (default: CPM)")
    parser.add_argument("--genome_size", help="Effective genome size for normalization (optional) See deepTools documentation for details.")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads (default: 4)")
    parser.add_argument("--barcodes", help="Optional file with barcodes to retain (1 per line)")
    parser.add_argument("--keep-temp", action="store_true", help="Keep intermediate files")
    return parser.parse_args()


def check_dependencies():
    for tool in ["bedtools", "samtools", "bamCoverage"]:
        if shutil.which(tool) is None:
            raise RuntimeError(f"Error: {tool} not found in PATH. Please install it.")


def read_barcodes(barcode_file):
    with open(barcode_file) as f:
        return set(line.strip() for line in f if line.strip())


def is_gzipped(filepath):
    return filepath.endswith(".gz")


def prepare_bed(input_file, output_bed, barcodes=None):
    barcodes_set = read_barcodes(barcodes) if barcodes else None
    open_func = gzip.open if is_gzipped(input_file) else open

    print("Reading input and writing BED...")
    with open_func(input_file, "rt") as infile, open(output_bed, "w") as out:
        for line in tqdm(infile, desc="Filtering fragments"):
            parts = line.strip().split("\t")
            if len(parts) < 3:
                continue
            if barcodes_set and parts[-1] not in barcodes_set:
                continue
            out.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\n")


def run_command(command, description):
    print(f"[Running] {description}: {' '.join(command)}")
    result = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if result.returncode != 0:
        print(result.stderr.decode())
        raise RuntimeError(f"{description} failed.")
    return result


def main():
    args = parse_args()
    check_dependencies()

    with tempfile.TemporaryDirectory() as tmpdir if not args.keep_temp else Path("bed2bigwig_temp") as tmpdir:
        tmpdir = Path(tmpdir)
        tmpdir.mkdir(exist_ok=True)

        bed_file = tmpdir / "fragments.bed"
        sorted_bed = tmpdir / "fragments.sorted.bed"
        bam_file = tmpdir / "fragments.bam"
        sorted_bam = tmpdir / "fragments.sorted.bam"

        prepare_bed(args.input, bed_file, args.barcodes)

        print("Sorting BED...")
        run_command(["sort", "-k1,1", "-k2,2n", str(bed_file)], "Sort BED")
        with open(sorted_bed, "w") as f_out:
            subprocess.run(["sort", "-k1,1", "-k2,2n", str(bed_file)], stdout=f_out)

        print("Converting BED to BAM...")
        with open(bam_file, "wb") as out_bam:
            subprocess.run(["bedtools", "bedtobam", "-i", str(sorted_bed), "-g", args.genome],
                           stdout=out_bam)

        print("Sorting BAM...")
        run_command(["samtools", "sort", "-@", str(args.threads), "-o", str(sorted_bam), str(bam_file)],
                    "Sort BAM")

        print("Indexing BAM...")
        run_command(["samtools", "index", str(sorted_bam)], "Index BAM")

        print("Running bamCoverage...")
        cmd = [
            "bamCoverage",
            "-b", str(sorted_bam),
            "-o", args.output,
            "--binSize", str(args.binsize),
            "--numberOfProcessors", str(args.threads),
        ]
        if args.normalize != "None":
            cmd += ["--normalizeUsing", args.normalize]
        run_command(cmd, "bamCoverage")

        print(f"\n✅ BigWig file written to: {args.output}")
        if args.keep_temp:
            print(f"🗂 Temporary files kept at: {tmpdir}")


if __name__ == "__main__":
    main()
