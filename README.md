# mu-PBWT GLIMPSE

## Build

```shell
make                 # downloads htslib 1.16 and boost 1.73 if missing, then builds everything
make DNANEXUS=1      # servers with system boost/htslib (static .a): skips the downloads, upstream tools use their `system` target
make original        # only the upstream tools under original/
make build_ser_from_bcf
```

`make` builds `phase`, `split_reference`, `build_ser_from_bin`, `build_ser_from_bcf` and the five upstream tools. The upstream tools are linked statically against the bundled boost and htslib; they need `libdeflate-dev`.

## phase options for mu-PBWT

All of them only matter with `--mupbwt`. The reference must have been indexed (`.ser` next to the `.bin`).

| Option                       | Default | Meaning                                                                                                                         |
| ---------------------------- | ------- | ------------------------------------------------------------------------------------------------------------------------------- |
| `--mupbwt`                   | off     | select conditioning haplotypes with mu-PBWT instead of GLIMPSE's compressed PBWT                                                |
| `--mupbwt-depth`             | 10      | selection depth; must be `<= --pbwt-depth` and `<= --max-depth`; bigger uses more memory                                        |
| `--mupbwt-chunk`             | 50      | site-scan chunk size; bigger uses more memory                                                                                   |
| `--mupbwt-max`               | 16      | multiplier of the per-target exact-step budget (`K * mupbwt-max`)                                                               |
| `--mupbwt-min-cm`            | 0.01    | matches shorter than this genetic length (cM) are not scored                                                                    |
| `--mupbwt-short-cm`          | 0.02    | boundary between short and medium matches (cM); must be `< --mupbwt-medium-cm`                                                  |
| `--mupbwt-medium-cm`         | 0.05    | boundary between medium and long matches (cM)                                                                                   |
| `--mupbwt-persistence`       | true    | carry matches over from the previous iteration                                                                                  |
| `--mupbwt-persistence-decay` | 0.9     | per-iteration decay applied to carried-over match scores                                                                        |
| `--mupbwt-persistence-floor` | 0.02    | carried-over matches whose decayed score falls below this are dropped                                                           |
| `--mupbwt-batch`             | 2       | maximum number of queries scanned together in one work item; must be `>= 1`. Bigger makes selection faster and uses more memory |
| `--mupbwt-numa`              | true    | when several threads are used, interleave the mu-PBWT index over the machine's NUMA memory nodes                                |

The cM thresholds use the genetic map the reference was built with. If the `.bin` was built without `--map`, cM is 1 cM per Mb and they are physical distances.

Measurements so far, all on the test machine used for the benchmarks:

- 1000 Genomes chr20, 24 CEU samples: batch 2 with NUMA on gave no gain over batch 1 with NUMA off on this small panel.
- Earlier 16-thread sweep: batch 8 made selection 16% faster than batch 2, but total wall time about 4% slower and 1.6 GiB more RAM.

`--mupbwt-smems` and `--mupbwt-mpsc` are still accepted but do nothing.

## Index tools

```shell
split_reference/bin/GLIMPSE2_split_reference --reference ref.bcf --input-region 1:11-500000 \
  --output-region 1:11-500000 --output out/ref --mupbwt            # writes .bin and .ser
build_ser_from_bin/bin/GLIMPSE2_build_ser_from_bin --bin out/ref_1_11_500000.bin --output out/ref_1_11_500000.ser
build_ser_from_bcf/bin/GLIMPSE2_build_ser_from_bcf index --input ref_region.bcf --output out/ref_1_11_500000.ser --threads 4
```

All three give the same `.ser` bytes. `phase` reads the `.ser` next to the `.bin` (same name, `.ser` extension).

## Docker

`Dockerfile` builds all of the above into one image: the four mu-PBWT executables
(`phase`, `split_reference`, `build_ser_from_bin`, `build_ser_from_bcf`) and the five upstream ones, plus `samtools`.
`Dockerfile_standalone` does the same but clones the repository itself, so it needs this layout to be in
the repository branch it clones (`--build-arg REPO_URL=... --build-arg GIT_REF=...`).

```shell
docker build --load -t glimpse2-mupbwt .
```

On `PATH` inside the image: the four mu-PBWT executables and the upstream `GLIMPSE2_chunk`, `GLIMPSE2_ligate`,
`GLIMPSE2_concordance`. The upstream `phase` and `split_reference` are not on `PATH` (same names as the
mu-PBWT ones); call them as `/app/original/phase/bin/GLIMPSE2_phase` and `/app/original/split_reference/bin/GLIMPSE2_split_reference`.

```shell
docker run --rm -v "$PWD:/data" -w /data glimpse2-mupbwt -c "
GLIMPSE2_split_reference --reference 10k/ref_panel_10k_msprime_sim.bcf \
  --input-region 1:77-9999973 --output-region 1:77-9999973 --output 10k/ref_panel_10k --mupbwt"

docker run --rm -v "$PWD:/data" -w /data glimpse2-mupbwt -c "
GLIMPSE2_build_ser_from_bin --bin 10k/ref_panel_10k_1_77_9999973.bin --output 10k/from_bin.ser"

docker run --rm -v "$PWD:/data" -w /data glimpse2-mupbwt -c "
GLIMPSE2_build_ser_from_bcf index --input 10k/ref_panel_10k_msprime_sim.bcf --output 10k/from_bcf.ser --threads 4"

docker run --rm -v "$PWD:/data" -w /data glimpse2-mupbwt -c "
GLIMPSE2_phase --bam-list all.txt --reference 10k/ref_panel_10k_1_77_9999973.bin \
  --output imputed.bcf --main 15 --burnin 5 --thread 5 --mupbwt"
```

The build was tested with `make` in a clean copy of this folder; the Dockerfiles themselves have not been built (no Docker on the build machine).
