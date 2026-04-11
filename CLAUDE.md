# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Build & Run

```bash
# Build
make build

# Install to /usr/bin (requires sudo)
make bai

# Run
bam2ab1 --bam <file.bam> --ref <file.fasta> [-o output.ab1]

# Tests
cargo test
```

## Architecture

**bam2ab1** converts BAM alignment files to AB1 trace files. The pipeline:

1. **Input**: Aligned BAM + reference FASTA
2. **Pileup generation** (`pileup_counter/`): Builds a multiple sequence alignment matrix and computes base ratios at each position
3. **Transformation** (`ab1/`): Converts pileup data to AB1 format using Gaussian peak modeling
4. **Output**: AB1 files (one per contig/chunk)

### Key Components

- `bam2ab1.rs` - CLI entry point; handles BAM/FASTA parsing, chunked processing
- `pileup_counter/plp_from_records.rs` - MSA construction, max insertion computation per reference position
- `ab1/plp2ab1.rs` - Gaussian pulse modeling (`Gaussian::compute`), peak width calculations
- `ab1/data_process.rs` - Traits (`TPlp2Ab1`) for transformation strategies
- `ab1/ab1_file.rs` - ABIF file format parsing/writing

### Transformation Strategies

- `Plp2Ab1WithDeletionShrink` - Adjusts peak positions for deletion regions
- `Plp2Ab1WithInsIdentifier` - Identifies insertion regions as `N` or `acgt`
- `Plp2Ab1Normal` - Standard transformation

### Parameters

- `--chunkSize` / `--ovlpSize` - For long sequences, split into overlapping chunks
- `--insIdent` - Insertion region marker (`N` or `acgt`)
