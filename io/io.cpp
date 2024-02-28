// ICON
//
// ---------------------------------------------------------------
// Copyright (C) 2004-2024, DWD, MPI-M, DKRZ, KIT, ETH, MeteoSwiss
// Contact information: icon-model.org
//
// See AUTHORS.TXT for a list of authors
// See LICENSES/ for license information
// SPDX-License-Identifier: BSD-3-Clause
// ---------------------------------------------------------------
//
#include "io.hpp"
#include <algorithm>
#include <map>

static const int NC_ERR = 2;
static const std::string BASE_VAR = "zg";

void io_muphys::parse_args(string &file, size_t &itime, real_t &dt, real_t &qnc,
                           int argc, char **argv) {
  file = "aes-new-gr_moderate-dt30s_atm_3d_ml_20080801T000000Z.nc";
  itime = 0; /* default to the first timestep in the input file */
  dt = 30.0;
  qnc = 100.0;
  char *end = nullptr;

  if (argc > 1) {
    file = argv[1];
  }
  cout << "input file: " << file << "\n";

  if (argc > 2) {
    itime = stoi(argv[2]);
  }
  cout << "itime: " << itime << "\n";

  if (argc > 3) {
    dt = type_converter<real_t>(argv[3], &end);
  }
  cout << "dt: " << dt << "\n";

  if (argc > 4) {
    qnc = type_converter<real_t>(argv[4], &end);
  }
  cout << "qnc: " << qnc << endl;
}

/* read-in time-constant data fields without a time dimension */
void io_muphys::input_vector(NcFile &datafile, array_1d_t<real_t> &v,
                             const string input, size_t &ncells, size_t &nlev) {
  NcVar var;
  v.resize(ncells * nlev);
  /*  access the input variable */
  try {
    var = datafile.getVar(input);
  } catch (NcNotVar &e) {
    cout << "FAILURE in accessing " << input << " (no time dimension) *******"
         << endl;
    e.what();
    e.errorCode();
  }
  /*  read-in input field values */
  try {
    array_1d_t<size_t> startp = {0, 0};
    array_1d_t<size_t> count = {nlev, ncells};
    var.getVar(startp, count, v.data());
  } catch (NcNotVar &e) {
    cout << "FAILURE in reading values from " << input
         << " (no time dimensions) *******" << endl;
    e.what();
    e.errorCode();
  }
}

void io_muphys::input_vector(NcFile &datafile, array_1d_t<real_t> &v,
                             const string input, size_t &ncells, size_t &nlev,
                             size_t itime) {
  NcVar att = datafile.getVar(input);
  try {
    v.resize(ncells * nlev);
    if (att.isNull()) {
      throw NC_ERR;
    }
    array_1d_t<size_t> startp = {itime, 0, 0};
    array_1d_t<size_t> count = {1, nlev, ncells};
    att.getVar(startp, count, v.data());
  } catch (NcException &e) {
    e.what();
    cout << "FAILURE in reading " << input << " **************************"
         << endl;
    throw NC_ERR;
  }
}

void io_muphys::output_vector(NcFile &datafile, array_1d_t<NcDim> &dims,
                              const string output, array_1d_t<real_t> &v,
                              size_t &ncells, size_t &nlev,
                              int &deflate_level) {
  // fortran:column major while c++ is row major
  NCreal_t ncreal_t;
  netCDF::NcVar var = datafile.addVar(output, ncreal_t, dims);

  for (size_t i = 0; i < nlev; ++i) {
    var.putVar({i, 0}, {1, ncells}, &v[i * ncells]);
  }
  if (deflate_level > 0) {
    var.setCompression(true, false, deflate_level);
  }
}
void io_muphys::output_vector(NcFile &datafile, array_1d_t<NcDim> &dims,
                              const string output,
                              std::map<std::string, NcVarAtt> varAttributes,
                              array_1d_t<real_t> &v, size_t &ncells,
                              size_t &nlev, int &deflate_level) {
  // fortran:column major while c++ is row major
  NCreal_t ncreal_t;
  netCDF::NcVar var = datafile.addVar(output, ncreal_t, dims);

  for (size_t i = 0; i < nlev; ++i) {
    var.putVar({i, 0}, {1, ncells}, &v[i * ncells]);
  }
  if (deflate_level > 0) {
    var.setCompression(true, false, deflate_level);
  }
  /* Add given attribues to the output variables (string, only) */
  for (auto &attribute_name : {"standard_name", "long_name", "units",
                               "coordinates", "CDI_grid_type"}) {
    auto attribute = varAttributes[attribute_name];

    std::string dataValues = "default";
    attribute.getValues(dataValues);

    var.putAtt(attribute.getName(), attribute.getType(),
               attribute.getAttLength(), dataValues.c_str());
  }
}

void io_muphys::read_fields(const string input_file, size_t &itime,
                            size_t &ncells, size_t &nlev, array_1d_t<real_t> &z,
                            array_1d_t<real_t> &t, array_1d_t<real_t> &p,
                            array_1d_t<real_t> &rho, array_1d_t<real_t> &qv,
                            array_1d_t<real_t> &qc, array_1d_t<real_t> &qi,
                            array_1d_t<real_t> &qr, array_1d_t<real_t> &qs,
                            array_1d_t<real_t> &qg) {

  NcFile datafile(input_file, NcFile::read);

  /*  read in the dimensions from the base variable: zg
   *  1st) vertical
   *  2nd) horizontal
   */
  auto baseDims = datafile.getVar(BASE_VAR).getDims();
  nlev = baseDims[0].getSize();
  ncells = baseDims[1].getSize();

  io_muphys::input_vector(datafile, z, "zg", ncells, nlev);
  io_muphys::input_vector(datafile, t, "ta", ncells, nlev, itime);

  io_muphys::input_vector(datafile, p, "pfull", ncells, nlev, itime);
  io_muphys::input_vector(datafile, rho, "rho", ncells, nlev, itime);
  io_muphys::input_vector(datafile, qv, "hus", ncells, nlev, itime);
  io_muphys::input_vector(datafile, qc, "clw", ncells, nlev, itime);
  io_muphys::input_vector(datafile, qi, "cli", ncells, nlev, itime);
  io_muphys::input_vector(datafile, qr, "qr", ncells, nlev, itime);
  io_muphys::input_vector(datafile, qs, "qs", ncells, nlev, itime);
  io_muphys::input_vector(datafile, qg, "qg", ncells, nlev, itime);

  datafile.close();
}

void io_muphys::write_fields(string output_file, size_t &ncells, size_t &nlev,
                             array_1d_t<real_t> &t, array_1d_t<real_t> &qv,
                             array_1d_t<real_t> &qc, array_1d_t<real_t> &qi,
                             array_1d_t<real_t> &qr, array_1d_t<real_t> &qs,
                             array_1d_t<real_t> &qg) {
  NcFile datafile(output_file, NcFile::replace);
  NcDim ncells_dim = datafile.addDim("ncells", ncells);
  NcDim nlev_dim = datafile.addDim("height", nlev);
  array_1d_t<NcDim> dims = {nlev_dim, ncells_dim};
  int deflate_level = 0;

  io_muphys::output_vector(datafile, dims, "ta", t, ncells, nlev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims, "hus", qv, ncells, nlev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims, "clw", qc, ncells, nlev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims, "cli", qi, ncells, nlev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims, "qr", qr, ncells, nlev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims, "qs", qs, ncells, nlev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims, "qg", qg, ncells, nlev,
                           deflate_level);

  datafile.close();
}

static void copy_coordinate_variables(NcFile &datafile, NcFile &inputfile,
                                      array_1d_t<std::string> coordinates) {
  for (auto &coordinate_name : coordinates) {
    auto coordinate = inputfile.getVar(coordinate_name);
    /* copy possible new dimensions from input coordinates to the output
     * befor adding the related data variables */
    for (NcDim &dim : coordinate.getDims()) {
      auto currentDims = datafile.getDims();
      /* map.contains() would be better, but requires c++20 */
      if (auto search = currentDims.find(dim.getName());
          search == currentDims.end())
        datafile.addDim(dim.getName(), dim.getSize());
    }
    auto var = datafile.addVar(coordinate.getName(), coordinate.getType(),
                               coordinate.getDims());

    for (auto &[attribute_name, attribute] : coordinate.getAtts()) {
      std::string dataValues = "default";
      attribute.getValues(dataValues);
      var.putAtt(attribute.getName(), attribute.getType(),
                 attribute.getAttLength(), dataValues.c_str());
    }

    /* Calculate the total size of the varibale */
    auto dimensions = coordinate.getDims();
    struct Prod {
      void operator()(NcDim n) { prod *= n.getSize(); }
      size_t prod{1};
    };
    size_t totalSize =
        std::for_each(dimensions.cbegin(), dimensions.cend(), Prod()).prod;

    /* Create a one-dimensional vector to store its values */
    array_1d_t<real_t> oneDimensionalVariable(totalSize);

    /* Copy the original values to the output */
    coordinate.getVar(oneDimensionalVariable.data());
    var.putVar(oneDimensionalVariable.data());
  }
}

void io_muphys::write_fields(string output_file, string input_file,
                             size_t &ncells, size_t &nlev,
                             array_1d_t<real_t> &t, array_1d_t<real_t> &qv,
                             array_1d_t<real_t> &qc, array_1d_t<real_t> &qi,
                             array_1d_t<real_t> &qr, array_1d_t<real_t> &qs,
                             array_1d_t<real_t> &qg) {
  NcFile datafile(output_file, NcFile::replace);
  NcFile inputfile(input_file, NcFile::read);
  auto baseDims = inputfile.getVar(BASE_VAR).getDims();
  NcDim nlev_dim =
      datafile.addDim(baseDims[0].getName(), baseDims[0].getSize());
  NcDim ncells_dim =
      datafile.addDim(baseDims[1].getName(), baseDims[1].getSize());

  copy_coordinate_variables(datafile, inputfile,
                            {"clon", "clon_bnds", "clat", "clat_bnds",
                             baseDims[0].getName(), "height_bnds"});
  /*TODO  height_bnds might have a different name */

  array_1d_t<NcDim> dims = {nlev_dim, ncells_dim};
  int deflate_level = 0;

  io_muphys::output_vector(datafile, dims, "ta",
                           inputfile.getVar("ta").getAtts(), t, ncells, nlev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims, "hus",
                           inputfile.getVar("hus").getAtts(), qv, ncells, nlev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims, "clw",
                           inputfile.getVar("clw").getAtts(), qc, ncells, nlev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims, "cli",
                           inputfile.getVar("cli").getAtts(), qi, ncells, nlev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims, "qr",
                           inputfile.getVar("qr").getAtts(), qr, ncells, nlev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims, "qs",
                           inputfile.getVar("qs").getAtts(), qs, ncells, nlev,
                           deflate_level);
  io_muphys::output_vector(datafile, dims, "qg",
                           inputfile.getVar("qg").getAtts(), qg, ncells, nlev,
                           deflate_level);

  inputfile.close();
  datafile.close();
}

void io_muphys::log_time(int64_t value) {
  std::ofstream time_file;
  time_file.open("timing.csv", std::ios_base::app);
  time_file << value << "\n";
  time_file.close();
}
