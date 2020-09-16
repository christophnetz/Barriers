#include <iostream>
#include <vector>
#include <random>
#include <cmath>
#include <chrono>
#include <numeric>
#include <fstream>

#include "rndutils.hpp"
#include "parameters.h"
#include "cmd_line.h"
using namespace std;

/////////////////////////////////////////////////////////////
//Global Parameters

std::mt19937_64 rng;
/////////////////////////////////////////////////////////////


class Individual
{
public:
  Individual(int genome_length) : genome(genome_length, 1.0) {};

  const double get_phenotype() {
    return genome[0];
  }

  void mutate(double m_rate, double m_shape, double tolerance, bool chain) {
    bernoulli_distribution mdist(m_rate);
    normal_distribution<double> sdist(0.0, m_shape);
    bool do_mutation = true;
    double change;
    for (int i = 0; i < genome.size(); ++i) {
      if (mdist(rng)) {

        change = sdist(rng);

        if (!chain) {
          if (i == 0) {
            for (int j = 1; j < genome.size(); ++j) {
              if (genome[j] + tolerance < change + genome[0]) {
                do_mutation = false;
              }
            }
          }
          else {
            if (genome[i] + change + tolerance < genome[0]) {
              do_mutation = false;
            }

          }
        }

        else {    //when the hierarchy is not wide
          if (i < genome.size() - 1) {
            if (genome[i + 1] + tolerance < genome[i] + change) {
              do_mutation = false;
            }
          }
          if (i > 0) {

            if (genome[i] + change + tolerance < genome[i - 1]) {
              do_mutation = false;
            }
          }
        }

        if (do_mutation) {
          genome[i] += change;

        }
      }
    }
  }


private:

  vector<double> genome;

};



void reproduction(vector<Individual>& pop, const vector<double>& phenotypes, const param_t params) {

  auto tmp_pop = pop;
  rndutils::mutable_discrete_distribution<int, rndutils::all_zero_policy_uni> rdist;


  rdist.mutate(phenotypes.cbegin(), phenotypes.cend());

  for (int i = 0; i < params.pop_size; ++i) {
    const int ancestor = rdist(rng);
    tmp_pop[i] = pop[ancestor];
    tmp_pop[i].mutate(params.mutation_rate, params.mutation_shape, params.tolerance, params.chain);
  }

  using std::swap;
  swap(pop, tmp_pop);
}


void run_sim(param_t params) {

  std::ofstream ofs;

  if (!params.outdir.empty()) {
    ofs.open(std::string(params.outdir + ".txt"));
    ofs << "g\tmean" << "\n";
  }

  unsigned seed =
    static_cast<unsigned>(std::chrono::high_resolution_clock::now().time_since_epoch().count());
  std::clog << "random_seed : " << seed << '\n';
  rng.seed(seed);

  Individual ind1(params.genome_length);
  vector<Individual> pop(params.pop_size, ind1);
  vector<double> phenotypes(params.pop_size);
  vector<double> fitness(params.pop_size);

  for (int g = 0; g < params.Gmax; g++) {

    for (int i = 0; i < params.pop_size; ++i) {
      phenotypes[i] = pop[i].get_phenotype();
    }

    double mean = std::accumulate(phenotypes.begin(), phenotypes.end(), 0.0) / phenotypes.size();
    if (ofs.is_open()) {
      ofs << g << "\t" << mean << "\n";
    }
    else if (!ofs.is_open()){
      cout << "Gen: \t" << g << "\t\tPop mean: \t" << mean << endl;
    }
    for (int i = 0; i < params.pop_size; ++i) {
      fitness[i] = 0.2 + (phenotypes[i] - mean);
      if (fitness[i] < 0.0)
        fitness[i] = 0.0;
    }

    reproduction(pop, fitness, params);
  }

  ofs.close();
}


int main(int argc, const char** argv) {

  try {
    auto param = param_t{};
    cmd::cmd_line_parser clp(argc, argv);
    auto config = clp.optional_val("config", std::string{});
    if (!config.empty()) clp.append(config_file_parser(config));
    clp.optional("Gmax", param.Gmax);
    clp.optional("pop_size", param.pop_size);
    clp.optional("mutation_rate", param.mutation_rate);
    clp.optional("mutation_shape", param.mutation_shape);
    clp.optional("tolerance", param.tolerance);
    clp.optional("chain", param.chain);
    clp.optional("genome_length", param.genome_length);
    clp.optional("outdir", param.outdir);

    auto unknown = clp.unrecognized();
    if (!unknown.empty()) {
      for (const auto& arg : unknown) {
        std::cerr << "unknown argument '" << arg << "'\n";
      }
      return 1;
    }

    run_sim(param);

    return 0;
  }

  catch (cmd::parse_error & err) {
    std::cerr << "\nParameter trouble: " << err.what() << '\n';
  }
  catch (std::exception & err) {
    std::cerr << "\nExeption caught: " << err.what() << '\n';
  }
  catch (...) {
    std::cerr << "\nUnknown exeption caught\n";
  }

  return 0;
}