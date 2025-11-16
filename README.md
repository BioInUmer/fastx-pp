<img width="200" height="200" alt="image" src="https://github.com/user-attachments/assets/6711c199-1e23-4574-bdff-a61b60d0f7f5" />

# fastx-pp

[![Python](https://img.shields.io/badge/Python-3.x-blue.svg)](https://www.python.org/)
[![License](https://img.shields.io/badge/License-MIT-green.svg)](LICENSE)
[![GitHub](https://img.shields.io/badge/GitHub-fastx--pp-black.svg)](https://github.com/BioInUmer/fastx-pp)

A command-line tool for preprocessing genomic sequence files (FASTA/FASTQ) with reverse-complement, trimming, and adaptor removal capabilities.

---

## 🎴 Features

- **Reverse-complement**: Generate the reverse complement of DNA sequences
- **Trim**: Remove a specified number of bases from both ends of sequences
- **Adaptor-removal**: Remove adaptor sequences from the start of reads
- **Format support**: Works with both FASTA and FASTQ files
- **Quality preservation**: Maintains quality scores in FASTQ files during all operations
- **Comprehensive statistics**: Detailed processing summaries with base composition analysis

## Documentation

> 📄 For detailed information about each function, implementation details, and comprehensive usage instructions, please refer to the **Report.pdf** document included in this repository. **Note: RUScript.sh = samurai.sh**

## Installation

### Prerequisites

- Python 3.x

### Clone the Repository

```bash
git clone https://github.com/BioInUmer/fastx-pp.git
cd fastx-pp
```

No additional dependencies are required beyond the Python standard library.

---

## ▶︎ Usage

### Basic Syntax

```bash
python3 fastx_pp.py --input <input_file> --output <output_file> --operation <operation> [options]
```

### Operations

#### 1. Reverse-complement

Generate the reverse complement of sequences:

```bash
python3 fastx_pp.py --input reads.fasta --output reads_rc.fasta --operation rc
```

#### 2. Trim

Remove bases from both ends of sequences:

```bash
python3 fastx_pp.py --input reads.fastq --output reads_trimmed.fastq --operation trim --trim-left 5 --trim-right 10
```

**Required arguments:**
- `--trim-left <num>`: Number of bases to remove from the 5' end
- `--trim-right <num>`: Number of bases to remove from the 3' end

#### 3. Adaptor-removal

Remove adaptor sequences from the beginning of reads:

```bash
python3 fastx_pp.py --input reads.fastq --output reads_clean.fastq --operation adaptor-removal --adaptor AGATCGGAAGAG
```

**Required arguments:**
- `--adaptor <sequence>`: Adaptor sequence to remove (only A, T, C, G allowed)

### Arguments

| Argument | Description | Required |
|----------|-------------|----------|
| `--input` | Input FASTA/FASTQ file | Yes |
| `--output` | Output FASTA/FASTQ file | Yes |
| `--operation` | Operation to perform: `rc`, `trim`, or `adaptor-removal` | Yes |
| `--trim-left` | Bases to trim from 5' end (for `trim` operation) | Conditional |
| `--trim-right` | Bases to trim from 3' end (for `trim` operation) | Conditional |
| `--adaptor` | Adaptor sequence to remove (for `adaptor-removal` operation) | Conditional |

### Important Notes

- Input and output files must have matching extensions (`.fasta` or `.fastq`)
- The tool automatically detects file format based on the first character (`@` for FASTQ, `>` for FASTA)
- For trimming operations, the total trim length must not exceed sequence length
- Adaptor sequences are case-insensitive
- The tool will prompt before overwriting existing output files

## Examples

### Example 1: Reverse-complement FASTA file

```bash
python3 fastx_pp.py --input sequences.fasta --output sequences_rc.fasta --operation rc
```

**Output:**
```
 File 'sequences.fasta' was successfully reversed-complemented ✅
 Check the output file → sequences_rc.fasta

=================================================================
                            SUMMARY                             
-----------------------------------------------------------------
 Total reads processed: 1.000
 Total bases processed: 50.000
 ↳(25% A, 25% C, 25% G, 25% T) | (0% N)
=================================================================
```

### Example 2: Trim FASTQ file

```bash
python3 fastx_pp.py --input reads.fastq --output reads_trimmed.fastq --operation trim --trim-left 3 --trim-right 7
```

**Output:**
```
 File 'reads.fastq' was successfully hard-trimmed ✅
 Check the output file → reads_trimmed.fastq

=================================================================
                            SUMMARY                             
-----------------------------------------------------------------
 Total reads processed: 500
 Total bases processed: 75.000
 ↳(23% A, 27% C, 26% G, 24% T) | (0% N)
·································································
 Total bases trimmed: 5.000
 ↳(22% A, 28% C, 25% G, 25% T) | (0% N)
=================================================================
```

### Example 3: Remove adaptor sequences

```bash
python3 fastx_pp.py --input reads.fastq --output reads_clean.fastq --operation adaptor-removal --adaptor AGATCGGAAGAG
```

**Output:**
```
 File 'reads.fastq' was successfully processed ✅
 Check the output file → reads_clean.fastq

=================================================================
                            SUMMARY                             
-----------------------------------------------------------------
 Total reads processed: 1.000
 Total bases processed: 150.000
 ↳(24% A, 26% C, 25% G, 25% T) | (0% N)
·································································
 Adaptor: AGATCGGAAGAG
 Total adaptors found & removed: 856
=================================================================
```

---

## File Structure

```
fastx-pp/
├── fastx_pp.py       # Main script
├── README.md         # This file
├── Report.pdf        # Detailed documentation
├── sample.fasta      # Sample FASTA file
├── sample.fastq      # Sample FASTQ file
└── LICENSE           # License file
```

---

## Error Handling

The tool includes comprehensive error handling for:
- Invalid file paths or permissions
- Incorrect file formats
- Invalid command-line arguments
- Trimming lengths exceeding sequence length
- Invalid nucleotide sequences

## Performance

fastx-pp is optimized for efficiency and can process large sequence files with minimal memory footprint by reading files line-by-line rather than loading entire files into memory.

---

## Contributing

Contributions, issues, and feature requests are welcome! Feel free to check the [issues page](https://github.com/BioInUmer/fastx-pp/issues).

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Citation

If you use fastx-pp in your research, please cite:

```
Kousar, M.U.H. (2025). fastx-pp: A FASTA/FASTQ Preprocessor. 
GitHub repository: https://github.com/BioInUmer/fastx-pp
```

## Acknowledgments

- Developed for bioinformatics sequence preprocessing workflows
- Designed with simplicity and efficiency in mind
- Inspired by common bioinformatics preprocessing needs

---

For questions or support, please open an issue on the [GitHub repository](https://github.com/BioInUmer/fastx-pp/issues).
