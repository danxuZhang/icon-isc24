#include <cstdlib>
#include <filesystem>
#include <gtest/gtest.h>
#include <iostream>
#include <regex>

#include <core/common/constants.hpp>
#include <core/properties/thermo.hpp>
#include <io/io.hpp>

double touched_cells(size_t &ke, size_t &ivend, array_1d_t<real_t> &t,
                     array_1d_t<real_t> &rho, array_1d_t<real_t> &qv,
                     array_1d_t<real_t> &qc, array_1d_t<real_t> &qi,
                     array_1d_t<real_t> &qr, array_1d_t<real_t> &qs,
                     array_1d_t<real_t> &qg) {
  size_t oned_vec_index;
  size_t jmx = 0;
  for (size_t i = ke - 1; i < ke; --i) {
    for (size_t j = 0; j < ivend; j++) {
      oned_vec_index = i * ivend + j;
      if ((std::max({qc[oned_vec_index], qr[oned_vec_index], qs[oned_vec_index],
                     qi[oned_vec_index], qg[oned_vec_index]}) >
           graupel_ct::qmin) or
          ((t[oned_vec_index] < graupel_ct::tfrz_het2) and
           (qv[oned_vec_index] >
            thermo::qsat_ice_rho(t[oned_vec_index], rho[oned_vec_index])))) {
        jmx = jmx + 1;
      }
    }
  }

  return static_cast<double>(jmx) / (ke * ivend);
}

TEST(IOTestSuite, CheckFile) {

  size_t ncells, nlev, itime = 0;
  array_1d_t<real_t> z, t, p, rho, qv, qc, qi, qr, qs, qg;
  std::vector<std::string> input_files;

  if (const char *env_file = std::getenv("TEST_FILE")) {
    input_files.push_back(env_file);
  } else {

    std::filesystem::path folder("input");
    const std::regex regex("[_[:alnum:]]+\\.nc",
                           std::regex::ECMAScript | std::regex::icase);

    if (!std::filesystem::is_directory(folder)) {
      throw std::runtime_error(folder.string() + " is not a folder");
    }

    for (const auto &entry : std::filesystem::directory_iterator(folder)) {
      const auto full_name = entry.path().string();

      if (entry.is_regular_file()) {
        const auto base_name = entry.path().filename().string();
        if (std::regex_match(base_name, regex)) {
          input_files.push_back(full_name);
        }
      }
    }
  }

  for (auto const &ifile : input_files) {
    io_muphys::read_fields(ifile, itime, ncells, nlev, z, t, p, rho, qv, qc, qi,
                           qr, qs, qg);

    double result = touched_cells(nlev, ncells, t, rho, qv, qc, qi, qr, qs, qg);
    std::cout << "CHECKING file " << ifile << ":";
    std::cout << " active cells fraction " << result << std::endl;
    // make sure at least 1 cell changed by the second loop of the algorithm
    EXPECT_GT(result, 0.);
  }
}
