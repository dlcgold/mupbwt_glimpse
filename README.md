# mu-PBWT GLIMPSE

## Build

```shell
make              # build all tools; htslib 1.16 and boost 1.73 are downloaded and built locally if missing
make DNANEXUS=1   # build all tools against the system htslib/boost (static .a), no downloads
make phase        # build one tool only: phase, split_reference, build_ser_from_bin or build_ser_from_bcf
make upstream     # build only the upstream GLIMPSE tools: chunk, ligate, concordance, inspect
make clean        # remove objects and binaries
```

`make` builds `phase`, `split_reference`, `build_ser_from_bin`, `build_ser_from_bcf` and the upstream GLIMPSE2 `chunk`, `ligate`, `concordance`, `inspect`. The upstream tools are linked statically against the bundled boost and htslib; they need `libdeflate-dev`, and the SIMDe headers (`simde/`, e.g. in `/usr/include`) are needed to compile `phase`. The downloads need GNU make >= 4.3 (grouped targets).

## Folders

| Folder                                               | Content                                                                                                                             |
| ---------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------- |
| `phase`                                              | `GLIMPSE2_phase`, with the mu-PBWT selection (`--mupbwt`); without `--mupbwt` it is the stock GLIMPSE2 phasing                      |
| `split_reference`                                    | `GLIMPSE2_split_reference`; without `--mupbwt` it is the stock reference splitting, with `--mupbwt` it also writes the `.ser` index |
| `build_ser_from_bin`                                 | builds the mu-PBWT `.ser` from an existing `.bin`                                                                                   |
| `build_ser_from_bcf`                                 | builds the mu-PBWT `.ser` from a reference BCF                                                                                      |
| `chunk`, `ligate`, `concordance`, `inspect`          | upstream GLIMPSE2 master tools, unmodified; `GLIMPSE2_inspect` reads a `split_reference` `.bin`                                     |
| `common`, `versions`                                 | shared sources and version header used by the tools                                                                                 |
| `phase_legacy_index`, `split_reference_legacy_index` | legacy variants for old mupbwt indexes, not built by the root `makefile`                                                            |

## Quick start

```shell
# 1. simulate low-coverage BAMs
mkdir -p unmasked_10k/reads100_150_default
./GLIMPSE2_simulate_bams_static --input-vcf unmasked_10k/target_unmasked_100_msprime_sim.bcf \
  --read-length 150 --contig 1 -O unmasked_10k/reads100_150_default --thread 7

# 2. index the BAMs and build the list it expects
samtools index -M unmasked_10k/reads100_150_default/*.bam
realpath unmasked_10k/reads100_150_default/*.bam > unmasked_10k/reads100_150_default/all.txt

# 3. split + index the reference panel
./split_reference/bin/GLIMPSE2_split_reference \
  --reference 10k/ref_panel_10k_msprime_sim.bcf \
  --input-region 1:77-9999973 --output-region 1:77-9999973 \
  --output 10k/ref_panel_10k --mupbwt

# 3b. alternatively, build the .ser from the .bin made without --mupbwt
./build_ser_from_bin/bin/GLIMPSE2_build_ser_from_bin \
  --bin 10k/ref_panel_10k_1_77_9999973.bin --output 10k/ref_panel_10k_1_77_9999973.ser

# 3c. or from the reference BCF cut to the same region
bcftools view -r 1:77-9999973 10k/ref_panel_10k_msprime_sim.bcf -Ob -o 10k/ref_panel_10k_region.bcf
./build_ser_from_bcf/bin/GLIMPSE2_build_ser_from_bcf index \
  --input 10k/ref_panel_10k_region.bcf --output 10k/ref_panel_10k_1_77_9999973.ser --threads 4

# 4. phase / impute
./phase/bin/GLIMPSE2_phase \
  --bam-list unmasked_10k/reads100_150_default/all.txt \
  --reference 10k/ref_panel_10k_1_77_9999973.bin \
  --output unmasked_10k/imputed.bcf \
  --main 15 --burnin 5 --thread 5 --mupbwt

# 5. allele frequencies of the reference panel (sites only), used by concordance
bcftools view -G -c0 10k/ref_panel_10k_msprime_sim.bcf -Ob > unmasked_10k/frq.bcf
bcftools index unmasked_10k/frq.bcf

# 6. concordance: one line per chromosome, <chr> <allele-freq bcf> <truth bcf> <imputed bcf>
echo "1 unmasked_10k/frq.bcf unmasked_10k/target_unmasked_100_msprime_sim.bcf unmasked_10k/imputed.bcf" > unmasked_10k/con.txt
./concordance/bin/GLIMPSE2_concordance --input unmasked_10k/con.txt \
  --output unmasked_10k/concordance \
  --bins 0.00000 0.00100 0.00200 0.00500 0.01000 0.05000 0.10000 0.20000 0.50000 \
  --thread 4 --gt-val

python conc.py unmasked_10k/concordance.rsquare.grp.txt.gz unmasked_10k/plot.pdf
```

## mu-PBWT flags

### split_reference

use `--mupbwt` to build mu-PBWT alongside stock:

|        | without `--mupbwt` (default) | with `--mupbwt`                                         |
| ------ | ---------------------------- | ------------------------------------------------------- |
| Builds | stock PBWT only              | mu-PBWT index and stock required information            |
| Output | `<output>_<region>.bin` only | `<output>_<region>.bin` **and** `<output>_<region>.ser` |

### phase options for mu-PBWT

All of them only matter with `--mupbwt`. The reference must have been indexed (`.ser` next to the `.bin`).

| Option                       | Default | Meaning                                                                                                                          |
| ---------------------------- | ------- | -------------------------------------------------------------------------------------------------------------------------------- |
| `--mupbwt`                   | off     | select conditioning haplotypes with mu-PBWT instead of GLIMPSE's compressed PBWT                                                 |
| `--mupbwt-depth`             | 10      | selection depth; must be `<= --pbwt-depth` and `<= --max-depth`; bigger uses more memory                                         |
| `--mupbwt-chunk`             | 50      | site-scan chunk size; bigger uses more memory                                                                                    |
| `--mupbwt-max`               | 16      | multiplier of the per-target exact-step budget (`K * mupbwt-max`)                                                                |
| `--mupbwt-min-cm`            | 0.01    | matches shorter than this genetic length (cM) are not scored                                                                     |
| `--mupbwt-short-cm`          | 0.02    | boundary between short and medium matches (cM); must be `< --mupbwt-medium-cm`                                                   |
| `--mupbwt-medium-cm`         | 0.05    | boundary between medium and long matches (cM)                                                                                    |
| `--mupbwt-persistence`       | true    | carry matches over from the previous iteration                                                                                   |
| `--mupbwt-persistence-decay` | 0.9     | per-iteration decay applied to carried-over match scores                                                                         |
| `--mupbwt-persistence-floor` | 0.02    | carried-over matches whose decayed score falls below this are dropped                                                            |
| `--mupbwt-batch`             | 2       | maximum number of queries scanned together in one work item; must be `>= 1`. Bigger makes selection faster and uses more memory  |
| `--mupbwt-numa`              | true    | when several threads are used, interleave the mu-PBWT index over the machine's NUMA memory nodes; disable with `--mupbwt-numa 0` |

## Docker

`Dockerfile` builds all of the above into one image: the four mu-PBWT executables
(`phase`, `split_reference`, `build_ser_from_bin`, `build_ser_from_bcf`) and the four upstream ones, plus `samtools`.
`Dockerfile_standalone` does the same but clones the repository itself, so it needs this layout to be in
the repository branch it clones (`--build-arg REPO_URL=... --build-arg GIT_REF=...`).

```shell
docker build --load -t glimpse2-mupbwt .
```

On `PATH` inside the image: the four mu-PBWT executables and the upstream `GLIMPSE2_chunk`, `GLIMPSE2_ligate`,
`GLIMPSE2_concordance`, `GLIMPSE2_inspect`.

```shell
docker run --rm -v "$PWD:/data" -w /data glimpse2-mupbwt -c "
GLIMPSE2_split_reference --reference 10k/ref_panel_10k_msprime_sim.bcf \
  --input-region 1:77-9999973 --output-region 1:77-9999973 --output 10k/ref_panel_10k --mupbwt"

docker run --rm -v "$PWD:/data" -w /data glimpse2-mupbwt -c "
GLIMPSE2_build_ser_from_bin --bin 10k/ref_panel_10k_1_77_9999973.bin --output 10k/from_bin.ser"

docker run --rm -v "$PWD:/data" -w /data glimpse2-mupbwt -c "
GLIMPSE2_build_ser_from_bcf index --input 10k/ref_panel_10k_region.bcf --output 10k/from_bcf.ser --threads 4"

docker run --rm -v "$PWD:/data" -w /data glimpse2-mupbwt -c "
GLIMPSE2_phase --bam-list all.txt --reference 10k/ref_panel_10k_1_77_9999973.bin \
  --output imputed.bcf --main 15 --burnin 5 --thread 5 --mupbwt"
```
