# build image
FROM debian:bookworm AS builder

# deps
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
  build-essential \
  g++ \
  make \
  git \
  ca-certificates \
  wget \
  autoconf \
  zlib1g-dev \
  libbz2-dev \
  liblzma-dev \
  libcurl4-openssl-dev \
  libssl-dev \
  && rm -rf /var/lib/apt/lists/*

# simde
ARG SIMDE_VERSION=0.8.2
RUN wget -q https://github.com/simd-everywhere/simde/archive/refs/tags/v${SIMDE_VERSION}.tar.gz -O /tmp/simde.tar.gz \
  && tar -xzf /tmp/simde.tar.gz -C /usr/include --strip-components=1 simde-${SIMDE_VERSION}/simde \
  && rm -f /tmp/simde.tar.gz

WORKDIR /app
COPY . .

# compile downloading htslib and boost
RUN make

# run image
FROM debian:bookworm-slim AS runtime

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
  libcurl4 \
  libssl3 \
  zlib1g \
  libbz2-1.0 \
  liblzma5 \
  libgomp1 \
  libstdc++6 \
  ca-certificates \
  && rm -rf /var/lib/apt/lists/*

WORKDIR /app

# copy exes
COPY --from=builder /app/boost ./boost
COPY --from=builder /app/htslib ./htslib
COPY --from=builder /app/chunk/bin ./chunk/bin
COPY --from=builder /app/concordance/bin ./concordance/bin
COPY --from=builder /app/ligate/bin ./ligate/bin
COPY --from=builder /app/phase/bin ./phase/bin
COPY --from=builder /app/split_reference/bin ./split_reference/bin
COPY --from=builder /app/GLIMPSE2_simulate_bams_static ./GLIMPSE2_simulate_bams_static

ENV PATH="/app/chunk/bin:/app/concordance/bin:/app/ligate/bin:/app/phase/bin:/app/split_reference/bin:${PATH}"

ENTRYPOINT ["/bin/bash"]
