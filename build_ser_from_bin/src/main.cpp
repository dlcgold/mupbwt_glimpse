#define _DECLARE_TOOLBOX_HERE
#include <utils/otools.h>

#include <containers/ref_haplotype_set.h>
#include <mupbwt/pbwt.h>

#include <boost/archive/binary_iarchive.hpp>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

static void print_usage() {
  std::cerr
      << "GLIMPSE2_build_mupbwt_index\n\n"
      << "Builds the mu-PBWT \".ser\" index from an existing stock \".bin\"\n"
      << "reference panel, without re-reading the reference BCF.\n\n"
      << "Usage:\n"
      << "  GLIMPSE2_build_mupbwt_index --bin <ref.bin> --output <out.ser>\n\n"
      << "Options:\n"
      << "  --bin <path>     Existing .bin reference panel (built by\n"
      << "                   GLIMPSE2_split_reference, --mupbwt or not).\n"
      << "  --output <path>  Output .ser path. Name it "
         "<prefix>_<region>.ser\n"
      << "                   to match what GLIMPSE2_phase --mupbwt expects\n"
      << "                   next to --reference <prefix>_<region>.bin.\n"
      << "  --help           Show this message.\n";
}

int main(int argc, char **argv) {
  std::string bin_path, output_path;

  for (int a = 1; a < argc; ++a) {
    std::string arg = argv[a];
    if (arg == "--bin" && a + 1 < argc) {
      bin_path = argv[++a];
    } else if (arg == "--output" && a + 1 < argc) {
      output_path = argv[++a];
    } else if (arg == "--help" || arg == "-h") {
      print_usage();
      return EXIT_SUCCESS;
    } else {
      std::cerr << "Unknown or incomplete argument: " << arg << "\n";
      print_usage();
      return EXIT_FAILURE;
    }
  }

  if (bin_path.empty() || output_path.empty()) {
    print_usage();
    return EXIT_FAILURE;
  }

  vrb.title("Building mu-PBWT index from an existing .bin reference panel:");
  vrb.bullet("Input .bin  : " + bin_path);
  vrb.bullet("Output .ser : " + output_path);

  ref_haplotype_set H;
  tac.clock();
  {
    std::ifstream ifs(bin_path, std::ios::binary | std::ios::in);
    if (!ifs) {
      vrb.error("Could not open .bin file: " + bin_path);
    }
    boost::archive::binary_iarchive ia(ifs);
    ia >> H;
  }
  vrb.bullet("Deserialized .bin [Nhaps=" + stb.str(H.n_ref_haps) +
             " Nsites=" + stb.str(H.n_tot_sites) +
             " Ncommon=" + stb.str(H.n_com_sites) +
             " Nrare=" + stb.str(H.n_rar_sites) + "] (" +
             stb.str(tac.rel_time() * 1.0 / 1000, 2) + "s)");

  if (H.n_ref_haps == 0 || H.n_tot_sites == 0) {
    vrb.error(".bin file carries no haplotypes/sites -- was it built "
              "correctly by GLIMPSE2_split_reference?");
  }
  if (H.n_ref_haps % 2 != 0) {
    vrb.error(
        "This tool assumes a fully diploid reference panel (n_ref_haps "
        "must be even); got n_ref_haps=" +
        stb.str(H.n_ref_haps) +
        ". The original split_reference --mupbwt path derives n_haps from "
        "n_ref_samples*2 directly and does not store n_ref_samples "
        "separately in the .bin, so a reference panel with haploid samples "
        "cannot be disambiguated here -- re-run split_reference --mupbwt "
        "from the BCF for that panel instead.");
  }
  const int n_ref_samples = (int)(H.n_ref_haps / 2);

  pbwt mupbwt;
  pbwt_build_state st;
  pbwt_build_af_init(&mupbwt, &st, n_ref_samples);

  tac.clock();
  std::vector<std::vector<uint32_t>> rare_deviators(H.n_tot_sites);
  for (unsigned int h = 0; h < H.n_ref_haps; ++h)
    for (int s : H.ShapRef[h])
      rare_deviators[(unsigned int)s].push_back(h);
  vrb.bullet("Inverted rare-site deviator index (" +
             stb.str(tac.rel_time() * 1.0 / 1000, 2) + "s)");

  std::vector<uint8_t> col(H.n_ref_haps);
  unsigned int i_common = 0;

  tac.clock();
  for (unsigned int s = 0; s < H.n_tot_sites; ++s) {
    if (H.flag_common[s]) {
      for (unsigned int h = 0; h < H.n_ref_haps; ++h)
        col[h] = H.HvarRef.get(i_common, h);
      ++i_common;
    } else {
      const uint8_t major = H.major_alleles[s] ? 1 : 0;
      const uint8_t minor = 1 - major;
      std::fill(col.begin(), col.end(), major);
      for (uint32_t h : rare_deviators[s])
        col[h] = minor;
    }
    pbwt_build_af_process_site_from_bits(&mupbwt, &st, col.data());
    vrb.progress("  * mu-PBWT reconstruction ",
                 std::min(1.0f, (s + 1.0f) / H.n_tot_sites));
  }
  pbwt_build_af_finalize(&mupbwt, &st);
  vrb.bullet("mu-PBWT reconstructed from .bin (" +
             stb.str(tac.rel_time() * 1.0 / 1000, 2) + "s)");

  FILE *out = fopen(output_path.c_str(), "wb");
  if (!out) {
    perror("fopen");
    vrb.error("Could not open output .ser file for writing: " + output_path);
  }
  if (pbwt_serialize(out, &mupbwt) != 0) {
    fclose(out);
    vrb.error("Serialization error writing: " + output_path);
  }
  fclose(out);
  pbwt_print_size(&mupbwt);
  free_pbwt(&mupbwt);

  vrb.bullet("Done. Wrote: " + output_path);
  return EXIT_SUCCESS;
}
