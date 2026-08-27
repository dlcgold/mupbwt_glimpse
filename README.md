# mupbwt_glimpse

## Build

```shell
make DNANEXUS=1     # builds with system boost/htslib
make                # downloads and builds htslib 1.16 + boost 1.73.0 before building
```

## Quick start

```shell
# 1. simulate low-coverage BAMs
mkdir -p unmasked_10k/reads100_150_default
./GLIMPSE2_simulate_bams_static --input-vcf target_unmasked_100_msprime_sim.bcf \
  --read-length 150 --contig 1 -O unmasked_10k/reads100_150_default --thread 7

# 2. build the BAM list phase expects (one absolute path per line)
realpath unmasked_10k/reads100_150_default/*.bam > unmasked_10k/reads100_150_default/all.txt

# 3. split + index the reference panel
./split_reference/bin/GLIMPSE2_split_reference \
  --reference 10k/ref_panel_10k_msprime_sim.bcf \
  --input-region 1:77-9999973 --output-region 1:77-9999973 \
  --output 10k/ref_panel_10k --mupbwt

# 4. phase / impute
./phase/bin/GLIMPSE2_phase \
  --bam-list unmasked_10k/reads100_150_default/all.txt \
  --reference 10k/ref_panel_10k_1_77_9999973.bin \
  --output unmasked_10k/imputed.bcf \
  --main 15 --burnin 5 --thread 5 --mupbwt

# 5. concordance — needs unmasked_10k/con.txt (created by hand, see below)
./concordance/bin/GLIMPSE2_concordance --input unmasked_10k/con.txt \
  --output unmasked_10k/concordance \
  --bins 0.00000 0.00100 0.00200 0.00500 0.01000 0.05000 0.10000 0.20000 0.50000 \
  --thread 4 --gt-val

python conc.py unmasked_10k/concordance.rsquare.grp.txt.gz unmasked_10k/plot.pdf
```

Create `unmasked_10k/con.txt` by hand — one whitespace-separated line per chromosome,
`<chr> <allele-freq bcf> <truth bcf> <imputed bcf>`:

```
1  unmasked_10k/frq.bcf  unmasked_10k/target_unmasked_100_msprime_sim.bcf  unmasked_10k/imputed.bcf
```

## mu-PBWT flags

### split_reference

use `--mupbwt` to build mu-PBWT alongside stock:

|        | without `--mupbwt` (default) | with `--mupbwt`                                         |
| ------ | ---------------------------- | ------------------------------------------------------- |
| Builds | stock PBWT only              | mu-PBWT index and stock required information            |
| Output | `<output>_<region>.bin` only | `<output>_<region>.bin` **and** `<output>_<region>.ser` |

### phase

Add `--mupbwt` to phase using the mu-PBWT. Optional flags:

| Flag                         | Default | Effect                                                                         |
| ---------------------------- | ------- | ------------------------------------------------------------------------------ |
| `--mupbwt-min-thr`           | 1000000 | divisor for the short match-length threshold                                   |
| `--mupbwt-medium-thr`        | 100000  | divisor for the medium match-length threshold (must stay < `--mupbwt-min-thr`) |
| `--mupbwt-max`               | 16      | multiplier for the max exact-match steps (bigger = more memory)                |
| `--mupbwt-chunk`             | 50      | site-scan chunk size (bigger = more memory)                                    |
| `--mupbwt-depth`             | 10      | selection depth (must be <= `--max-depth`, bigger = more memory)               |
| `--mupbwt-persistence`       | true    | carry over decayed matches from the previous iteration                         |
| `--mupbwt-persistence-decay` | 0.9     | per-iteration decay applied to carried-over match scores                       |
| `--mupbwt-persistence-floor` | 0.02    | drop a carried-over match once its decayed score falls below this              |

Use more memory should results in better R^2 results.

## Docker

A `Dockerfile` at the repo root builds all five tools (`chunk`, `split_reference`, `phase`, `ligate`,
`concordance`) plus `GLIMPSE2_simulate_bams_static` into a self-contained image. No local dependencies needed besides Docker (or Podman).

```shell
docker build --load -t glimpse2-mupbwt .
```

Mount your working directory (here the current directory) to `/data`; every path below is relative
to `/data`. Run the steps in order — the BAM list and `con.txt` are created between runs.

```shell
# 1. simulate low-coverage BAMs
docker run --rm -v "$PWD:/data" -w /data glimpse2-mupbwt -c "
mkdir -p unmasked_10k/reads100_150_default && \
/app/GLIMPSE2_simulate_bams_static --input-vcf target_unmasked_100_msprime_sim.bcf \
  --read-length 150 --contig 1 -O unmasked_10k/reads100_150_default --thread 7"

# 2. build the BAM list phase expects (one path per line)
docker run --rm -v "$PWD:/data" -w /data glimpse2-mupbwt -c "
realpath unmasked_10k/reads100_150_default/*.bam > unmasked_10k/reads100_150_default/all.txt"

# 3. split + index the reference panel
docker run --rm -v "$PWD:/data" -w /data glimpse2-mupbwt -c "
/app/split_reference/bin/GLIMPSE2_split_reference \
  --reference 10k/ref_panel_10k_msprime_sim.bcf \
  --input-region 1:77-9999973 --output-region 1:77-9999973 \
  --output 10k/ref_panel_10k --mupbwt"

# 4. phase / impute
docker run --rm -v "$PWD:/data" -w /data glimpse2-mupbwt -c "
/app/phase/bin/GLIMPSE2_phase \
  --bam-list unmasked_10k/reads100_150_default/all.txt \
  --reference 10k/ref_panel_10k_1_77_9999973.bin \
  --output unmasked_10k/imputed.bcf \
  --main 15 --burnin 5 --thread 5 --mupbwt"

# 5. concordance — create unmasked_10k/con.txt by hand first (see "Quick start"):
#    1  unmasked_10k/frq.bcf  unmasked_10k/target_unmasked_100_msprime_sim.bcf  unmasked_10k/imputed.bcf
docker run --rm -v "$PWD:/data" -w /data glimpse2-mupbwt -c "
/app/concordance/bin/GLIMPSE2_concordance --input unmasked_10k/con.txt \
  --output unmasked_10k/concordance \
  --bins 0.00000 0.00100 0.00200 0.00500 0.01000 0.05000 0.10000 0.20000 0.50000 \
  --thread 4 --gt-val"
```
