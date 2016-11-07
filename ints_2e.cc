
// standard C++ headers
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include <chrono>
#include <thread>
#include <atomic>
#include <iterator>
#include <unordered_map>
#include <mutex>
#include <memory>

// Eigen matrix algebra library
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <Eigen/Cholesky>

// have BTAS library?
#ifdef LIBINT2_HAVE_BTAS
# include <btas/btas.h>
#endif // LIBINT2_HAVE_BTAS

// Libint Gaussian integrals library
#include <libint2.hpp>
#include <libint2/diis.h>


typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>
        Matrix;  // import dense, dynamically sized Matrix type from Eigen;
                 // this is a matrix with row-major storage (http://en.wikipedia.org/wiki/Row-major_order)
                 // to meet the layout of the integrals returned by the Libint integral library
typedef Eigen::DiagonalMatrix<double, Eigen::Dynamic, Eigen::Dynamic>
        DiagonalMatrix;

using libint2::Shell;
using libint2::BasisSet;
using libint2::Operator;
using libint2::BraKet;
using libint2::Engine;

struct Atom {
    int atomic_number;
    double x, y, z;
};

std::vector<Atom> read_geometry(const std::string& filename);

std::vector<libint2::Shell> make_sto3g_basis(const std::vector<Atom>& atoms);
std::vector<libint2::Shell> make_ucsto3g_basis(const std::vector<Atom>& atoms);
std::vector<libint2::Shell> make_pvdz_basis(const std::vector<Atom>& atoms);
std::vector<libint2::Shell> make_augpvdz_basis(const std::vector<Atom>& atoms);
std::vector<libint2::Shell> make_ucpvdz_basis(const std::vector<Atom>& atoms);
std::vector<libint2::Shell> make_ucaugpvdz_basis(const std::vector<Atom>& atoms);

size_t nbasis(const std::vector<libint2::Shell>& shells);

size_t max_nprim(const std::vector<libint2::Shell>& shells);

int max_l(const std::vector<libint2::Shell>& shells);

std::vector<size_t> map_shell_to_basis_function(const std::vector<libint2::Shell>& shells);

Matrix compute_1body_ints(const std::vector<libint2::Shell>& shells,
                          libint2::Operator t,
                          const std::vector<Atom>& atoms = std::vector<Atom>());

int main(int argc, char *argv[]) {

  using std::cout;
  using std::cerr;
  using std::endl;

  try {

    /*** =========================== ***/
    /*** initialize molecule         ***/
    /*** =========================== ***/

    // read geometry from a file; by default read from h2o.xyz, else take filename (.xyz) from the command line
    const auto filename = (argc > 1) ? argv[1] : "h2o.xyz";

    std::vector<Atom> atoms = read_geometry(filename);

    // count the number of electrons
    auto nelectron = 0;
    for (auto i = 0; i < atoms.size(); ++i)
      nelectron += atoms[i].atomic_number;
    const auto ndocc = nelectron / 2;
    cout << "# of electrons = " << nelectron << endl;

    // compute the nuclear repulsion energy
    auto enuc = 0.0;
    for (auto i = 0; i < atoms.size(); i++)
      for (auto j = i + 1; j < atoms.size(); j++) {
        auto xij = atoms[i].x - atoms[j].x;
        auto yij = atoms[i].y - atoms[j].y;
        auto zij = atoms[i].z - atoms[j].z;
        auto r2 = xij*xij + yij*yij + zij*zij;
        auto r = sqrt(r2);
        enuc += atoms[i].atomic_number * atoms[j].atomic_number / r;
      }
    cout << "Nuclear repulsion energy = " << std::setprecision(15) << enuc << endl;

    libint2::Shell::do_enforce_unit_normalization(false);

    cout << "Atomic Cartesian coordinates (a.u.):" << endl;
    for(const auto& a: atoms)
      std::cout << a.atomic_number << " " << a.x << " " << a.y << " " << a.z << std::endl;

    // read basis sets
    auto shells = make_pvdz_basis(atoms);
//    auto shells = make_ucpvdz_basis(atoms);
//    auto shells = make_augpvdz_basis(atoms);
//    auto shells = make_ucaugpvdz_basis(atoms);
//    auto shells = make_sto3g_basis(atoms);
//    auto shells = make_ucsto3g_basis(atoms);
    size_t nao = 0;
    for (auto s=0; s<shells.size(); ++s)
      nao += shells[s].size();
    cout << "number of ao: " << nao << endl;

    // initializes the Libint integrals library ... now ready to compute
    libint2::initialize();

    // compute overlap integrals
    auto S = compute_1body_ints(shells, Operator::overlap);
    cout << "\nOverlap Integrals:\n";
    cout << S << endl<< endl;

    //
    // compute electron repulsion integrals
    //
    const auto n = nbasis(shells);
    double* g = new double[n * n * n * n];

    // print 2e-integral into a file
    std::ofstream outfile;
    outfile.open("G_2eints.out");

    // construct the electron repulsion integrals engine
    Engine engine(Operator::coulomb, max_nprim(shells), max_l(shells), 0);

    auto shell2bf = map_shell_to_basis_function(shells);

    // matrix is symmetric, but skipping it here for simplicity
    for(auto s1=0; s1!=shells.size(); ++s1) {

      auto bf1_first = shell2bf[s1]; // first basis function in this shell
      auto n1 = shells[s1].size();

      for(auto s2=0; s2!=shells.size(); ++s2) {

        auto bf2_first = shell2bf[s2];
        auto n2 = shells[s2].size();

        // loop over shell pairs of the density matrix, {s3,s4}
        // again symmetry is not used for simplicity
        for(auto s3=0; s3!=shells.size(); ++s3) {

          auto bf3_first = shell2bf[s3];
          auto n3 = shells[s3].size();

          for(auto s4=0; s4!=shells.size(); ++s4) {

            auto bf4_first = shell2bf[s4];
            auto n4 = shells[s4].size();

            // Coulomb contribution to the Fock matrix is from {s1,s2,s3,s4} integrals
            const auto* buf_1234 = engine.compute(shells[s1], shells[s2], shells[s3], shells[s4]);

            // loop over every integral in the shell set (= nested loops over basis functions in each shell)
            for(auto f1=0, f1234=0; f1!=n1; ++f1) {
              const auto bf1 = f1 + bf1_first;
              for(auto f2=0; f2!=n2; ++f2) {
                const auto bf2 = f2 + bf2_first;
                for(auto f3=0; f3!=n3; ++f3) {
                  const auto bf3 = f3 + bf3_first;
                  for(auto f4=0; f4!=n4; ++f4, ++f1234) {
                    const auto bf4 = f4 + bf4_first;
                    g[((bf1*n + bf2)*n + bf3)*n + bf4] = buf_1234[f1234];
                  }
                }
              }
            }

          }
        }
      }
    }

    for(auto i=0, idx=0; i!=n; ++i) {
      for(auto j=0; j!=n; ++j) {
        for(auto k=0; k!=n; ++k) {
          for(auto l=0; l!=n; ++l, ++idx) {
//            printf("%10.8f   ", g[idx]);
            outfile << std::setprecision(12) << g[idx] << endl;
          }
          //cout << endl;
        }
        //cout << endl;
      }
      //cout << endl;
    }

    outfile.close();

    double ehf = 0;

    printf("** Hartree-Fock energy = %20.12f\n", ehf + enuc);

    libint2::finalize(); // done with libint

  } // end of try block; if any exceptions occurred, report them and exit cleanly

  catch (const char* ex) {
    cerr << "caught exception: " << ex << endl;
    return 1;
  }
  catch (std::string& ex) {
    cerr << "caught exception: " << ex << endl;
    return 1;
  }
  catch (std::exception& ex) {
    cerr << ex.what() << endl;
    return 1;
  }
  catch (...) {
    cerr << "caught unknown exception\n";
    return 1;
  }

  return 0;
}

// this reads the geometry in the standard xyz format supported by most chemistry software
std::vector<Atom> read_dotxyz(std::istream& is) {
  // line 1 = # of atoms
  size_t natom;
  is >> natom;
  // read off the rest of line 1 and discard
  std::string rest_of_line;
  std::getline(is, rest_of_line);

  // line 2 = comment (possibly empty)
  std::string comment;
  std::getline(is, comment);

  std::vector<Atom> atoms(natom);
  for (auto i = 0; i < natom; i++) {
    std::string element_label;
    double x, y, z;
    is >> element_label >> x >> y >> z;

    // .xyz files report element labels, hence convert to atomic numbers
    int Z;
    if (element_label == "H")
      Z = 1;
    else if (element_label == "C")
      Z = 6;
    else if (element_label == "N")
      Z = 7;
    else if (element_label == "O")
      Z = 8;
    else if (element_label == "F")
      Z = 9;
    else if (element_label == "S")
      Z = 16;
    else if (element_label == "Cl")
      Z = 17;
    else {
      std::cerr << "read_dotxyz: element label \"" << element_label << "\" is not recognized" << std::endl;
      throw "Did not recognize element label in .xyz file";
    }

    atoms[i].atomic_number = Z;

    // .xyz files report Cartesian coordinates in angstroms; convert to bohr
    const auto angstrom_to_bohr = 0.52917721092; // 2010 CODATA value
    atoms[i].x = x * angstrom_to_bohr;
    atoms[i].y = y * angstrom_to_bohr;
    atoms[i].z = z * angstrom_to_bohr;
//    atoms[i].x = x;
//    atoms[i].y = y;
//    atoms[i].z = z;
  }

  return atoms;
}

std::vector<Atom> read_geometry(const std::string& filename) {

  std::cout << "Will read geometry from " << filename << std::endl;
  std::ifstream is(filename);
  assert(is.good());

  // to prepare for MPI parallelization, we will read the entire file into a string that can be
  // broadcast to everyone, then converted to an std::istringstream object that can be used just like std::ifstream
  std::ostringstream oss;
  oss << is.rdbuf();
  // use ss.str() to get the entire contents of the file as an std::string
  // broadcast
  // then make an std::istringstream in each process
  std::istringstream iss(oss.str());

  // check the extension: if .xyz, assume the standard XYZ format, otherwise throw an exception
  if ( filename.rfind(".xyz") != std::string::npos)
    return read_dotxyz(iss);
  else
    throw "only .xyz files are accepted";
}

size_t nbasis(const std::vector<libint2::Shell>& shells) {
  size_t n = 0;
  for (const auto& shell: shells)
    n += shell.size();
  return n;
}

size_t max_nprim(const std::vector<libint2::Shell>& shells) {
  size_t n = 0;
  for (auto shell: shells)
    n = std::max(shell.nprim(), n);
  return n;
}

int max_l(const std::vector<libint2::Shell>& shells) {
  int l = 0;
  for (auto shell: shells)
    for (auto c: shell.contr)
      l = std::max(c.l, l);
  return l;
}

std::vector<libint2::Shell> make_sto3g_basis(const std::vector<Atom>& atoms) {

  using libint2::Shell;

  std::vector<Shell> shells;

  for(auto a=0; a<atoms.size(); ++a) {

    // STO-3G basis set
    // cite: W. J. Hehre, R. F. Stewart, and J. A. Pople, The Journal of Chemical Physics 51, 2657 (1969)
    //       doi: 10.1063/1.1672392
    // obtained from https://bse.pnl.gov/bse/portal
    switch (atoms[a].atomic_number) {
      case 1: // Z=1: hydrogen
        shells.push_back(
            {
              {3.425250910, 0.623913730, 0.168855400}, // exponents of primitive Gaussians
              {  // contraction 0: s shell (l=0), spherical=false, contraction coefficients
                {0, false, {0.15432897, 0.53532814, 0.44463454}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}   // origin coordinates
            }
        );
        break;

      case 6: // Z=6: carbon
        shells.push_back(
            {
              {71.616837000, 13.045096000, 3.530512200},
              {
                {0, false, {0.15432897, 0.53532814, 0.44463454}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {2.941249400, 0.683483100, 0.222289900},
              {
                {0, false, {-0.09996723, 0.39951283, 0.70011547}}
              }, 
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {2.941249400, 0.683483100, 0.222289900},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {0.15591627, 0.60768372, 0.39195739}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        break;

      case 7: // Z=7: nitrogen
        shells.push_back(
            {
              {99.106169000, 18.052312000, 4.885660200},
              {
                {0, false, {0.15432897, 0.53532814, 0.44463454}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {3.780455900, 0.878496600, 0.285714400},
              {
                {0, false, {-0.09996723, 0.39951283, 0.70011547}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
          {3.780455900, 0.878496600, 0.285714400},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {0.15591627, 0.60768372, 0.39195739}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        break;

      case 8: // Z=8: oxygen
        shells.push_back(
            {
              {130.709320000, 23.808861000, 6.443608300},
              {
                {0, false, {0.15432897, 0.53532814, 0.44463454}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {5.033151300, 1.169596100, 0.380389000},
              {
                {0, false, {-0.09996723, 0.39951283, 0.70011547}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {5.033151300, 1.169596100, 0.380389000},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {0.15591627, 0.60768372, 0.39195739}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        break;

      default:
        throw "do not know STO-3G basis for this Z";
    }

  }

  return shells;
}

std::vector<libint2::Shell> make_ucsto3g_basis(const std::vector<Atom>& atoms) {

  using libint2::Shell;

  std::vector<Shell> shells;

  for(auto a=0; a<atoms.size(); ++a) {

    // uncontracted STO-3G basis set
    // cite: W. J. Hehre, R. F. Stewart, and J. A. Pople, The Journal of Chemical Physics 51, 2657 (1969)
    //       doi: 10.1063/1.1672392
    // obtained from https://bse.pnl.gov/bse/portal
    switch (atoms[a].atomic_number) {
      case 1: // Z=1: hydrogen
        shells.push_back(
            {
              {3.425250910}, // exponents of primitive Gaussians
              {  // contraction 0: s shell (l=0), spherical=false, contraction coefficients
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}   // origin coordinates
            }
        );
        shells.push_back(
            {
              {0.623913730}, // exponents of primitive Gaussians
              {  // contraction 0: s shell (l=0), spherical=false, contraction coefficients
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}   // origin coordinates
            }
        );
        shells.push_back(
            {
              {0.168855400}, // exponents of primitive Gaussians
              {  // contraction 0: s shell (l=0), spherical=false, contraction coefficients
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}   // origin coordinates
            }
        );
        break;

      case 6: // Z=6: carbon
        shells.push_back(
            {
              {71.616837000, 13.045096000, 3.530512200},
              {
                {0, false, {0.15432897, 0.53532814, 0.44463454}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {2.941249400, 0.683483100, 0.222289900},
              {
                {0, false, {-0.09996723, 0.39951283, 0.70011547}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {2.941249400, 0.683483100, 0.222289900},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {0.15591627, 0.60768372, 0.39195739}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        break;

      case 7: // Z=7: nitrogen
        shells.push_back(
            {
              {99.106169000, 18.052312000, 4.885660200},
              {
                {0, false, {0.15432897, 0.53532814, 0.44463454}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {3.780455900, 0.878496600, 0.285714400},
              {
                {0, false, {-0.09996723, 0.39951283, 0.70011547}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
          {3.780455900, 0.878496600, 0.285714400},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {0.15591627, 0.60768372, 0.39195739}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        break;

      case 8: // Z=8: oxygen
        shells.push_back(
            {
              {130.709320000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {23.808861000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {6.443608300},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {5.033151300},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {1.169596100},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {0.380389000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {5.033151300},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {1.169596100},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {0.380389000},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        break;

      default:
        throw "do not know STO-3G basis for this Z";
    }

  }

  return shells;
}

std::vector<libint2::Shell> make_pvdz_basis(const std::vector<Atom>& atoms) {

  using libint2::Shell;

  std::vector<Shell> shells;

  for(auto a=0; a<atoms.size(); ++a) {

    // cc-pVDZ basis set for H atom
    // obtained from https://bse.pnl.gov/bse/portal
    switch (atoms[a].atomic_number) {
      case 1: // Z=1: hydrogen
        shells.push_back(
            {
              {13.0100000, 1.9620000, 0.4446000}, // exponents of primitive Gaussians
              {  // contraction 0: s shell (l=0), spherical=false, contraction coefficients
                {0, false, {0.0196850, 0.1379770, 0.4781480}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}   // origin coordinates
            }
        );
        shells.push_back(
            {
              {0.1220000},
              { // contraction 0: p shell (l=1), spherical=false
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {0.7270000},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        break;
      case 8: // Z=8: oxygen
        shells.push_back(
            {
              {11720.0000000,1759.0000000, 400.8000000, 113.7000000, 37.0300000, 13.2700000, 5.0250000, 1.0130000},
              {
                {0, false, {0.0007100, 0.0054700, 0.0278370, 0.1048000, 0.2830620, 0.4487190, 0.2709520, 0.0154580}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {11720.0000000, 1759.0000000, 400.8000000, 113.7000000, 37.0300000, 13.2700000, 5.0250000, 1.0130000},
              {
                {0, false, {-0.0001600, -0.0012630, -0.0062670, -0.0257160, -0.0709240, -0.1654110, -0.1169550, 0.5573680}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {0.3023000},
              { // contraction 0: p shell (l=1), spherical=false
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {17.7000000, 3.8540000, 1.0460000},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {0.0430180, 0.2289130, 0.5087280}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {0.2753000},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {1.1850000},
              { // contraction 0: d shell (l=2), spherical=false
                {2, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        break;
      default:
        throw "do not know cc-pVDZ basis for this Z";
    }

  }

  return shells;
}

std::vector<libint2::Shell> make_ucpvdz_basis(const std::vector<Atom>& atoms) {

  using libint2::Shell;

  std::vector<Shell> shells;

  for(auto a=0; a<atoms.size(); ++a) {

    // cc-pVDZ basis set for H atom
    // obtained from https://bse.pnl.gov/bse/portal
    switch (atoms[a].atomic_number) {
      case 1: // Z=1: hydrogen
        shells.push_back(
            {
              {13.0100000}, // exponents of primitive Gaussians
              {  // contraction 0: s shell (l=0), spherical=false, contraction coefficients
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}   // origin coordinates
            }
        );
        shells.push_back(
            {
              {1.9620000}, // exponents of primitive Gaussians
              {  // contraction 0: s shell (l=0), spherical=false, contraction coefficients
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}   // origin coordinates
            }
        );
        shells.push_back(
            {
              {0.4446000}, // exponents of primitive Gaussians
              {  // contraction 0: s shell (l=0), spherical=false, contraction coefficients
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}   // origin coordinates
            }
        );
        shells.push_back(
            {
              {0.1220000},
              { // contraction 0: p shell (l=1), spherical=false
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {0.7270000},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        break;
      case 8: // Z=8: oxygen
        shells.push_back(
            {
              {11720.0000000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {1759.0000000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {400.8000000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {113.7000000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {37.0300000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {13.2700000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {5.0250000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {1.0130000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {11720.0000000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {1759.0000000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {400.8000000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {113.7000000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {37.0300000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {13.2700000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {5.0250000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {1.0130000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {0.3023000},
              { // contraction 0: p shell (l=1), spherical=false
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {17.7000000},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {3.8540000},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {1.0460000},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {0.2753000},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {1.1850000},
              { // contraction 0: d shell (l=2), spherical=false
                {2, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        break;
      default:
        throw "do not know cc-pVDZ basis for this Z";
    }

  }

  return shells;
}

std::vector<libint2::Shell> make_augpvdz_basis(const std::vector<Atom>& atoms) {

  using libint2::Shell;

  std::vector<Shell> shells;

  for(auto a=0; a<atoms.size(); ++a) {

    // aug-cc-pVDZ basis set for O atom
    // obtained from https://bse.pnl.gov/bse/portal
    switch (atoms[a].atomic_number) {

      case 8: // Z=8: oxygen
        shells.push_back(
            {
              {11720.0000000,1759.0000000, 400.8000000, 113.7000000, 37.0300000, 13.2700000, 5.0250000, 1.0130000},
              {
                {0, false, {0.0007100, 0.0054700, 0.0278370, 0.1048000, 0.2830620, 0.4487190, 0.2709520, 0.0154580}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {11720.0000000, 1759.0000000, 400.8000000, 113.7000000, 37.0300000, 13.2700000, 5.0250000, 1.0130000},
              {
                {0, false, {-0.0001600, -0.0012630, -0.0062670, -0.0257160, -0.0709240, -0.1654110, -0.1169550, 0.5573680}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {0.3023000},
              { // contraction 0: p shell (l=1), spherical=false
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {17.7000000, 3.8540000, 1.0460000},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {0.0430180, 0.2289130, 0.5087280}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {0.2753000},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {1.1850000},
              { // contraction 0: d shell (l=2), spherical=false
                {2, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {0.3320000},
              { // contraction 0: d shell (l=2), spherical=false
                {2, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        break;

      default:
        throw "do not know aug-cc-pVDZ basis for this Z";
    }

  }

  return shells;
}

std::vector<libint2::Shell> make_ucaugpvdz_basis(const std::vector<Atom>& atoms) {

  using libint2::Shell;

  std::vector<Shell> shells;

  for(auto a=0; a<atoms.size(); ++a) {

    // uncontracted aug-cc-pVDZ basis set for O atom
    // obtained from https://bse.pnl.gov/bse/portal
    switch (atoms[a].atomic_number) {
       case 8: // Z=8: oxygen
        shells.push_back(
            {
              {11720.0000000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {1759.0000000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {400.8000000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {113.7000000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {37.0300000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {13.2700000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {5.0250000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {1.0130000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {11720.0000000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {1759.0000000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {400.8000000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {113.7000000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {37.0300000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {13.2700000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {5.0250000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {1.0130000},
              {
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {0.3023000},
              { // contraction 0: p shell (l=1), spherical=false
                {0, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {17.7000000},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {3.8540000},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {1.0460000},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {0.2753000},
              { // contraction 0: p shell (l=1), spherical=false
                {1, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {1.1850000},
              { // contraction 0: d shell (l=2), spherical=false
                {2, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        shells.push_back(
            {
              {0.3320000},
              { // contraction 0: d shell (l=2), spherical=false
                {2, false, {1.0}}
              },
              {{atoms[a].x, atoms[a].y, atoms[a].z}}
            }
        );
        break;
      default:
        throw "do not know cc-pVDZ basis for this Z";
    }

  }

  return shells;
}

std::vector<size_t> map_shell_to_basis_function(const std::vector<libint2::Shell>& shells) {
  std::vector<size_t> result;
  result.reserve(shells.size());

  size_t n = 0;
  for (auto shell: shells) {
    result.push_back(n);
    n += shell.size();
  }

  return result;
}

Matrix compute_1body_ints(const std::vector<libint2::Shell>& shells,
                          libint2::Operator obtype,
                          const std::vector<Atom>& atoms)
{
  using libint2::Shell;
  using libint2::Engine;
  using libint2::Operator;

  const auto n = nbasis(shells);
  Matrix result(n,n);

  // construct the overlap integrals engine
  Engine engine(obtype, max_nprim(shells), max_l(shells), 0);
  // nuclear attraction ints engine needs to know where the charges sit ...
  // the nuclei are charges in this case; in QM/MM there will also be classical charges
  if (obtype == Operator::nuclear) {
    std::vector<std::pair<double,std::array<double,3>>> q;
    for(const auto& atom : atoms) {
      q.push_back( {static_cast<double>(atom.atomic_number), {{atom.x, atom.y, atom.z}}} );
    }
    engine.set_params(q);
  }

  auto shell2bf = map_shell_to_basis_function(shells);

  // loop over unique shell pairs, {s1,s2} such that s1 >= s2
  // this is due to the permutational symmetry of the real integrals over Hermitian operators: (1|2) = (2|1)
  for(auto s1=0; s1!=shells.size(); ++s1) {

    auto bf1 = shell2bf[s1]; // first basis function in this shell
    auto n1 = shells[s1].size();

    for(auto s2=0; s2<=s1; ++s2) {

      auto bf2 = shell2bf[s2];
      auto n2 = shells[s2].size();

      // compute shell pair; return is the pointer to the buffer
      const auto* buf = engine.compute(shells[s1], shells[s2]);

      // "map" buffer to a const Eigen Matrix, and copy it to the corresponding blocks of the result
      Eigen::Map<const Matrix> buf_mat(buf, n1, n2);
      result.block(bf1, bf2, n1, n2) = buf_mat;
      if (s1 != s2) // if s1 >= s2, copy {s1,s2} to the corresponding {s2,s1} block, note the transpose!
      result.block(bf2, bf1, n2, n1) = buf_mat.transpose();

    }
  }

  return result;
}
