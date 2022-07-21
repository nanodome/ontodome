#ifndef NANONETWORK_H
#define NANONETWORK_H

#include "../../ontodome.h"

class nanoCell : public SoftwareModel {

 private:

    /// cell simulation time
    double t = 0.;

 public:
    /// methods
    PBMFractalParticlePhase<PBMAggregate<Particle>>* pbm;
    NucleationTheory* nt;
    GasModels* gp;

    /// Species container
    std::vector<SingleComponentComposition> specs;

    /// Constructor
    nanoCell(std::vector<SingleComponentComposition> species, double p, double T, std::valarray<double> c) {

      // Get all the required inputs
      nt = findNearest<NucleationTheory>();
      gp = findNearest<GasModels>();
      pbm = findNearest<PBMFractalParticlePhase<PBMAggregate<Particle>>>();

      // initialize the species container
      specs = species;

      // Set the methods for each cell
      pbm->set_max_aggregates(500);
      pbm->set_min_aggregates(50);

      // Update the gas phase
      this->gp->update(p,T,c);
    }

    int print_aggregates_number();

    std::valarray<double> advance(double dt, double p, double T, std::valarray<double> cs);

};

class nanoNetwork {
  private:
    /// number of network cells
    int n_cells;

    /// methods containers for cells
    std::vector<std::shared_ptr<nanoCell>> cells;

    /// Network simulation time [s]
    double t = 0.;

    /// Network last timestep [s]
    double dt = 0.;

  public:

    /// Constructor
    nanoNetwork(std::vector<std::shared_ptr<nanoCell>> _cells){
      // set the number of cells
      n_cells = _cells.size();
      cells = _cells;
    };

    /// Prints the current number of aggregates for each cell
    void print_cells();

    /// Get network's time [s]
    double get_t() { return this->t; }

    /// Get network's current timestep [s]
    double get_dt() { return this->dt; }

    /// Get particles diameters for cell with index idx [m]
    std::valarray<double> get_cell_particles_diameters(int idx);

    /// Get aggregates diameters for cell with index idx [m]
    std::valarray<double> get_cell_aggregates_diameters(int idx);

    /// Get pbm volume for cell with index idx [m]
    double get_cell_pbm_volume(int idx) { return this->cells[idx]->pbm->get_volume(); }

    /// Get pbm mean fractal dimension for cell with index idx [m]
    double get_cell_mean_fractal_dimension(int idx) { return this->cells[idx]->pbm->get_mean_fractal_dimension(); }

    /// Advances the network using a synchronous timestep
    std::vector<std::valarray<double>> timestep(std::vector<double> p, std::vector<double> T, std::vector<double> vels, std::vector<std::valarray<double>> cs);

  private:

    double calc_dt ();
};

int nanoCell::print_aggregates_number() {
  return pbm->get_aggregates_number();
}

std::valarray<double> nanoCell::advance(double dt, double p, double T, std::valarray<double> cs) {

  // Update the gasphase
  gp->update(p,T,cs);

  // advance simulation time
  t += dt;

  // Advance simulation time using the Strang algorithm

  // Species source term for the gas phase
  double g_prec = 0.;

  // species consumption container
  std::valarray<double> cons(0.,int(specs.size()));

  // Strang first step
  gp->timestep(dt/2.0,cons);
  pbm->volume_expansion(dt/2.);

  // Strang second step
  // species source term for the condensing species in gas phase
  g_prec -= pbm->timestep(dt);

  cons[0] = g_prec;
  gp->timestep(dt,cons);

  // Strang third step
  cons[0] = 0.;
  gp->timestep(dt/2.0,cons);
  pbm->volume_expansion(dt/2.0);

  return(gp->get_c());
}


void nanoNetwork::print_cells() {
  for (int k=0; k < n_cells; k++) {
    std::cout << "Cell number: " << k << " with: " << cells[k]->print_aggregates_number() << " particles" << std::endl;
  }
}

/// Get particles diameters for cell with index idx [m]
std::valarray<double> nanoNetwork::get_cell_particles_diameters(int idx) {

  return this->cells[idx]->pbm->get_particles_sizes();
}

/// Get aggregates diameters for cell with index idx [m]
std::valarray<double> nanoNetwork::get_cell_aggregates_diameters(int idx) {

  return this->cells[idx]->pbm->get_aggregates_sizes();
}

/// Advances the network using a synchronous timestep
std::vector<std::valarray<double>> nanoNetwork::timestep(std::vector<double> p, std::vector<double> T, std::vector<double> vels, std::vector<std::valarray<double>> cs) {

  // the timestep is choosen as the minimum above all cells
  dt = calc_dt();

  // advance network simulation time
  t += dt;

  // Move particles from cell to cell based on fluxes;
  // Collect the aggregates to be removed until the flux value is reached
  // The aggregates are select randomly
  PBMFractalParticlePhase<PBMAggregate<Particle>> *dest_pbm = nullptr;

  for (int k=0; k < n_cells; k++) {

    if (k+1 < n_cells) {
      dest_pbm = cells[k+1]->pbm;
    }

    cells[k]->pbm->move_aggregates(dt,vels[k],dest_pbm);

  }

  // initialize the molar fractions container
  std::vector<std::valarray<double>> cs_o;

  // advance the simulations and store the results in the molar fractions container
  for (int k=0; k < n_cells; k++) {
    cs_o.push_back(cells[k]->advance(dt,p[k],T[k],cs[k]));
  }

  return cs_o;
}

double nanoNetwork::calc_dt () {

  double dt = 1e-8; // 1e-8 is the maximum time-step allowed by PBM

  for (int k=0; k < n_cells; k++) {

    // calculate the timestep using an exponential waiting time
    double R_tot = cells[k]->pbm->get_total_processes_rate();
    double rho = ndm::uniform_double_distr(ndm::rand_gen);

    // exponential waiting time
    double dt_t = -log(rho)/R_tot;
    dt = std::min(dt,dt_t);

  }

  return dt;
}

#endif // NANONETWORK_H
