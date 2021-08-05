/*
    NanoDome - H2020 European Project NanoDome, GA n.646121
    (www.nanodome.eu, https://github.com/nanodome/nanodome)
    e-mail: Emanuele Ghedini, emanuele.ghedini@unibo.it

    Copyright (C) 2018  Alma Mater Studiorum - Università di Bologna

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU Lesser General Public License as
    published by the Free Software Foundation, either version 3 of the
    License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#ifndef CONSTRAINEDLANGEVINPARTICLEPHASE_H
#define CONSTRAINEDLANGEVINPARTICLEPHASE_H

#include "../particlephase/dynamicparticlephase.h"
#include "../PBM/pbmfractalparticlephase.h"

#include <time.h>

//#define OMP
#include <omp.h>

template<typename A>
class ConstrainedLangevinParticlePhase : public DynamicParticlePhase<A> {

	/// Function for calculating the value of a 3d index a grid in a periodic domain
	///	\param	int c_idx: cell actual component position
	///	\param int SIDE: Number of cells in the dimension
	///	\param double& shift: shift to impose in the point coordinates at the correspondig component
	///	\param double domain_side: lenght of the control volume side
	int periodic_position(int c_idx, int SIDE, double& shift, double domain_side);

public:

        /// Constructor
        ConstrainedLangevinParticlePhase(double _volume):
                DynamicParticlePhase<A>(_volume){}

/*
	/// Constructor from PBM Particle phase
	///	\param PBMFractalParticlePhase<PBMAggregate<Particle>>& _fractal_pp PBM particle phase
	ConstrainedLangevinParticlePhase(PBMFractalParticlePhase<PBMAggregate<Particle>>& _fractal_pp, double _volume, double T) :
		DynamicParticlePhase<A>(_volume)
	{

		double PBM_Volume = _fractal_pp.get_volume();


		int N = _fractal_pp.get_aggregates_number();

		//if volumes are the same

		if (PBM_Volume == _volume) {
			// Spawn aggregates from PBM aggregates
			for (int i = 0; i < N; i++) {
				std::shared_ptr<PBMAggregate<Particle>> pbm_agg = _fractal_pp.get_aggregate(i);
				this->add_particle_lc(pbm_agg->get_n_monomers(), pbm_agg->get_species());

			}
		}
		else if (_volume < PBM_Volume) {

			// Choose ramdom aggregates from the PBM Particle Phase
//			double pbm_p_density = N / PBM_Volume;
			int n_particles = (N*_volume) / PBM_Volume;
			std::cout << "PBM Particles:" << N << std::endl;
			std::cout << "Langevin Particles" << n_particles << std::endl;

			srand(time(NULL));
			int count = 0;
			while (count < n_particles) {
				int rand_index = rand() % N;
				std::shared_ptr<PBMAggregate<Particle>> pbm_agg = _fractal_pp.get_aggregate(rand_index);
				this->add_particle_lc(pbm_agg->get_n_monomers(), pbm_agg->get_species());
				count++;
			}

		}
		else { // if the volume is bigger than the PBM one
			// TO BE IMPLEMENTED
		}

		// Initialize velocities
		this->initialize_velocities(T);
	}
*/
	/// Constructor for Particle Phase starting from a list of aggregates (hotstart/recovery)
	ConstrainedLangevinParticlePhase(std::list<std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& _aggregates, double _volume, double _T):
		DynamicParticlePhase<A>(_volume) {

		// Assign the list
		this->aggregates = _aggregates;

		// Initialize velocities
		this->initialize_velocities(_T);

		// Displace aggregates in the grid
		this->grid.displace_aggregates(this->aggregates);
	}

    /// Saves VTK data
    ///	\param int iteration: iteration in which the snapshot is saved
    ///	\param s_idx: int streamline index
    void save_vtk(int iteration, std::string root_path);

    /// Time step
    ///	\param double dt: timestep [sec]
    ///	\param const GasModels& gp: GasPhase object [by reference]
    ///	\param const NucleationTheory& nt: Nucleation Theory object [by reference]
    ///	\param double& T: Temperature[K] [by reference]
    double timestep(double dt, GasModels* gp, NucleationTheory* nt, double& T, SingleComponentComposition* sp);

    /// Linked-Cell time step
    ///	\param double dt: timestep [sec]
    ///	\param const GasModels& gp: GasPhase object [by reference]
    ///	\param const NucleationTheory& nt: Nucleation Theory object [by reference]
    ///	\param double& T: Temperature[K] [by reference]
    double timestep_lc(double dt, GasModels* gp, NucleationTheory* nt, SingleComponentComposition* sp, double& T);

        /// Returns cardinality of the active particles in the simulation
    int get_aggregates_cardinality();

    /// Updates graph constraints values
    void update_constraints();

    /// biggest spherical enclosure diameter
    double get_biggest_spherical_enclosure(int& _id) const;


protected:

    /// Motion
    ///	\param dt timestep [s]
    ///	\param gp Gas Phase object
    ///	\param T Temperature [K]
    void motion(double dt, GasModels* gp, double T);

    /// Motion (Linked Cell Version)
    ///	\param dt timestep [s]
    ///	\param gp Gas Phase object
    ///	\param T Temperature [K]
    void motion_lc(double dt, GasModels* gp, double T);

    ///Coagulation
    void coagulation();

    /// Coagulation (Linked Cell Version)
    void coagulation_lc();

    /// Creates a Particle in the simulation
    void add_particle(double n, SingleComponentComposition* s, NucleationTheory* nt);

    /// Creates a Particle in the simulation (Linked Cell Version)
    void add_particle_lc(double n, SingleComponentComposition* s, NucleationTheory* nt);

    /// Checks if some particles has coalesced
    void rearrange_aggregate();

    /// Keeps constant the number of aggregates in the simulation
    void aggregates_number_balance();

    /// Initialize Velocity
    ///	\param double T Temperature [K]
    void initialize_velocities(double T);

    /// Evaluate nucleation events for species s and return the molecules consumption
    /// [#/m3 s]
    /// \param dt timestep [s]
    /// \param J particle nucleation rate [#/m3 s]
    /// \param s condensing species
    double nucleation(double dt, double J, double j, SingleComponentComposition* s, NucleationTheory* nt);

    /// Evaluate nucleation events for species s and return the molecules consumption (Linked Cell Version)
    /// [#/m3 s]
    /// \param dt timestep [s]
    /// \param J particle nucleation rate [#/m3 s]
    /// \param s condensing species
    double nucleation_lc(double dt, double J, double j, SingleComponentComposition* s);

    /// Expands the volume and shifts the aggregates
    /// Correct volume due to gas expansion
    /// \param dt timestep [s]
    /// \param gp gas phase surrounding particles
    void volume_resize(double dt, GasModels* gp);


};

#include "constrainedlangevinparticlephase.h"
#include "../utilities/ndm_random.h"

// External Libs
#include "../../../../tools//tinyxml2.h"
#include "../../../../tools/utilities.h"

/*
 s t*d::vector<double> getParticleCoords(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble);
 std::vector<double> getParticleMass(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble);
 std::vector<int>	getParticleAggregate(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble);
 std::vector<double> getParticleDiameter(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble);
 std::vector<double>	getAggregatesDiameter(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble);
 std::vector<int>	getParticleID(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble);

 int getTotalEdges(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble);
 */

template<typename A>
int ConstrainedLangevinParticlePhase<A>::periodic_position(int c_idx, int SIDE, double& shift, double domain_side) {
	if (c_idx < 0) { shift = +domain_side; return SIDE - 1; }
	else if (c_idx > (SIDE - 1)) { shift = -domain_side; return 0; }
	else {
		return c_idx;
	}

}


// s t*d::vector<std::tuple<int, int, int>> getConnections(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble);

 // vect[0] -> couples
 // vect[1] -> types
 // vect[2] -> offsets
std::vector<std::vector<int>> GetAbsoluteIndices(std::vector<std::tuple<int, int, int>> graph_data, std::vector<int> p_list);

template<class dataType>std::string V_To_S(std::vector<dataType> vect);



template<typename A>
double ConstrainedLangevinParticlePhase<A>::timestep(double dt, GasModels *gp, NucleationTheory* nt, double& T, SingleComponentComposition* sp) {

	double g_si = 0.0; // Consumption

	// get Temperature from the Gas Phase
	T = gp->get_T();

	// nucleating species concentration and temperature for this timestep
//	double ns = gp->get_n();
	//double S = ns / species[0].n_sat(T);
//	double S = ns / gp->get_n_sat(sp);

	// nucleation parameters
	double J_si = nt->nucleation_rate();
	double j_si = nt->stable_cluster_diameter();

	// Condensation parameters
	double F_si = nt->condensation_rate();

	// Consumption of species from the gas phase and nucleation and condensation in the particle phase
	g_si += nucleation(dt, J_si, j_si, sp, nt);

	// Condensation
	g_si += this->condensation(dt, F_si);

	// Move particles
	motion(dt, gp, T);

	// Coagulation
	coagulation();

	// Sintering
	this->sintering(dt, T);

	// Check Coalescence and applies shake after sintering
	rearrange_aggregate();

	// Update volume from gas phase data
	this->volume_expansion(dt,gp);

	// adjust the aggregates number
	this->aggregates_number_balance();

	return g_si;
}

template<typename A>
double ConstrainedLangevinParticlePhase<A>::timestep_lc(double dt, GasModels* gp, NucleationTheory* nt, SingleComponentComposition* sp, double& T) {

	double g_si = 0.0; // Consumption

	// get Temperature from the Gas Phase
	T = gp->get_T();

	// nucleating species concentration and temperature for this timestep
//	double ns = gp->get_n();
	//double S = ns / species[0].n_sat(T);
//	double S = ns / gp->get_n_sat(sp);

	// nucleation parameters
	double J_si = nt->nucleation_rate();
	double j_si = nt->stable_cluster_diameter();

	// Condensation parameters
	double F_si = nt->condensation_rate();

	// Consumption of species from the gas phase and nucleation and condensation in the particle phase
	g_si += nucleation_lc(dt, J_si, j_si, sp);

	// Surface Condensation
	g_si += this->condensation(dt, F_si);

	// Particle Motion
	this->motion_lc(dt, gp, T);

	// Coagulation
	this->coagulation_lc();

	// Sintering
	this->sintering(dt, T);

	// Update constraints
	this->update_constraints();

	// Check Coalescence and applies shake after sintering
	this->rearrange_aggregate();

	// Update the grid
	int biggest_agg = -1;
	double biggest_s_enclosure = this->get_biggest_spherical_enclosure(biggest_agg);
	this->grid.grid_update(this->aggregates, this->volume, biggest_s_enclosure, biggest_agg);

	// Keeps constant the number of aggregates
	if (T < 2600)
		this->aggregates_number_balance();


	// Update volume from gas phase data
	//this->volume_expansion(dt, gp);
	//this->volume_resize(dt, gp);

	// Update the grid
	biggest_s_enclosure = this->get_biggest_spherical_enclosure(biggest_agg);
	this->grid.grid_update(this->aggregates, this->volume, biggest_s_enclosure, biggest_agg);


	// Update constraints
	//this->update_constraints();

	return g_si;

}

template<typename A>
void ConstrainedLangevinParticlePhase<A>::volume_resize(double dt, GasModels* gp) {

//	double volume_prev = this->volume;
	double scale_factor = exp(gp->get_gamma()*dt);
	this->volume = this->volume* scale_factor;

	double shift = pow(scale_factor, 1.0 / 3.0);

	std::valarray<double> x_shift = { shift,shift,shift };

	for (auto agg = this->aggregates.begin(); agg != this->aggregates.end(); agg++) {
		(*agg)->scale_coordinates(x_shift);
	}

}


template<typename A>
void ConstrainedLangevinParticlePhase<A>::rearrange_aggregate() {

	for (auto agg = this->aggregates.begin(); agg != this->aggregates.end(); agg++) {
		// Check if the aggregates has coalesced
		if ((*agg)->check_coalescence()) {
			// if yes, update graph structure
			(*agg)->update_graph();
			// Update actual cardinality
			(*agg)->update_particles_cardinality();
		}
		else {
			// Update the constraints with the new distances
			(*agg)->update_constraints();
		}
	}

	std::list<int> erased_agg;

	#ifdef OMP

	bool SHAKE_ERROR;

	#pragma omp parallel
	#pragma omp single nowait
	{
		for (auto agg = this->aggregates.begin(); agg != this->aggregates.end(); agg++) {
			#pragma omp task firstprivate(agg)
			SHAKE_ERROR = (*agg)->SHAKE();
			// Apply shake to the sintered aggregate
			#pragma omp critical
			if (!(SHAKE_ERROR)) {
				std::cout << "SHAKE DIVERGENCE IN REARRANGING AGGREGATES AFTER SINTERING" << std::endl;
				erased_agg.push_back((*agg)->get_id());
				this->volume *= (this->aggregates.size() - 1.0) / this->aggregates.size();
			}
		}
		#pragma omp taskwait
	}
	#else
	for (auto agg = this->aggregates.begin(); agg != this->aggregates.end(); agg++) {
		// Apply shake to the sintered aggregate
		if (!((*agg)->SHAKE())) {
			std::cout << "SHAKE DIVERGENCE IN REARRANGING AGGREGATES AFTER SINTERING" << std::endl;
			erased_agg.push_back((*agg)->get_id());
			this->volume *= (this->aggregates.size() - 1.0) / this->aggregates.size();
		}
	}
	#endif


	/*
		f *or (auto agg = this->aggregates.begin(); agg != this->aggregates.end(); agg++) {
		// Check if the aggregates has coalesced
		if ((*agg)->check_coalescence()) {
			// if yes, update graph structure
			(*agg)->update_graph();
			// Update actual cardinality
			(*agg)->update_particles_cardinality();
}
else {
	// Update the constraints with the new distances
	(*agg)->update_constraints();
}

// Apply shake to the sintered aggregate
if (!((*agg)->SHAKE())) {
	std::cout << "SHAKE DIVERGENCE IN REARRANGING AGGREGATES AFTER SINTERING" << std::endl;
	erased_agg.push_back((*agg)->get_id());
	this->volume *= (this->aggregates.size() - 1.0) / this->aggregates.size();
}
}

*/

	if (erased_agg.size() > 0) {
		// Erase from aggregate
		for(auto e = erased_agg.begin(); e != erased_agg.end(); e++){
			for(auto a = this->aggregates.begin(); a !=this->aggregates.end(); a++){
				if((*a)->get_id() == (*e)){
					this->aggregates.erase(a);
					break;
				}
			}
			// Erase from grid
			for(int c = 0; c < this->grid.get_n_cells(); c++){
				for(int a = 0; a < this->grid[c].get_n(); a++){
					if(this->grid[c][a]->get_id() == (*e)){
						this->grid[c].erase_n(a);
						break;
					}
				}
			}
		}
	}

}

template<typename A>
void ConstrainedLangevinParticlePhase<A>::add_particle(double n, SingleComponentComposition* s, NucleationTheory* nt) {

	// Control volume side
	double side = pow(this->volume, 1.0 / 3.0);

	// generate position of the particle
	std::valarray<double> x = { ndm::uniform_double_distr(ndm::rand_gen),
		ndm::uniform_double_distr(ndm::rand_gen),
		ndm::uniform_double_distr(ndm::rand_gen) };

		// scale coordinates on the control volume
		x = side*(x - 0.5);

		// create particle
		std::shared_ptr<DynamicParticle> data(new DynamicParticle(n, s, nt, x));
		std::shared_ptr<RATTLEAggregate<DynamicParticle>> agg(new RATTLEAggregate<DynamicParticle>(data));

		// Add particle to the ensamble.
		this->aggregates.push_back(agg);

}

template<typename A>
void ConstrainedLangevinParticlePhase<A>::add_particle_lc(double n, SingleComponentComposition* s, NucleationTheory* nt) {

	// Control volume side
	double side = pow(this->volume, 1.0 / 3.0);

	// generate position of the particle
	std::valarray<double> x = { ndm::uniform_double_distr(ndm::rand_gen),
				    ndm::uniform_double_distr(ndm::rand_gen),
				    ndm::uniform_double_distr(ndm::rand_gen) };

	// scale coordinates on the control volume
	x = side*(x - 0.5);

	// create particle
	std::shared_ptr<DynamicParticle> data(new DynamicParticle(n, s, nt, x));
	std::shared_ptr<RATTLEAggregate<DynamicParticle>> agg(new RATTLEAggregate<DynamicParticle>(data));

	// Add particle to the ensamble.
	this->aggregates.push_back(agg);

	// Update the grid
	int biggest_agg = -1;
	double biggest_s_enclosure = this->get_biggest_spherical_enclosure(biggest_agg);
	this->grid.grid_update(this->aggregates, this->volume, biggest_s_enclosure, biggest_agg);

	// Add pointer to particle in the grid
	this->grid.insert_aggregate(agg);

}

template<typename A>
void ConstrainedLangevinParticlePhase<A>::initialize_velocities(double T) {

	for (auto agg = this->aggregates.begin(); agg != this->aggregates.end(); agg++) {
		(*agg)->v_initialization(T);
	}

}

template<typename A>
void ConstrainedLangevinParticlePhase<A>::update_constraints() {

	for (auto agg = this->aggregates.begin(); agg != this->aggregates.end(); agg++) {
		(*agg)->update_constraints();
	}

}

/// biggest spherical enclosure diameter
template<typename A>
double ConstrainedLangevinParticlePhase<A>::get_biggest_spherical_enclosure(int& _id) const {

	double d;
	double d_max = DBL_MIN;
	for (auto& a : this->aggregates) {
		d = a->get_enclosing_sphere_diameter();
		if (d > d_max) {
			d_max = d;
			_id = a->get_id();
		}
	}

	return d_max;
}

template<typename A>
void ConstrainedLangevinParticlePhase<A>::save_vtk(int iter, std::string root_path) {

	using namespace tinyxml2;

	// All points (nanoparticles) coords
	std::vector<double> coords = getParticleCoords(this->aggregates);

	// All points (nanoparticles) mass
	std::vector<double> masses = getParticleMass(this->aggregates);

	// All points (nanoparticles) diameter
	std::vector<double> diameter = getParticleDiameter(this->aggregates);

	// All aggregates diameter
	std::vector<double> agg_diameter = getAggregatesDiameter(this->aggregates);

	// All points (nanoparticles) affiliation
	std::vector<int> affil = getParticleAggregate(this->aggregates);

	// All points IDs
	std::vector<int> p_IDs = getParticleID(this->aggregates);

	// Get the connectivity data
	std::vector <std::tuple<int,int,int>> graph_data = getConnections(this->aggregates);

	std::vector<std::vector<int>> abs_idx = GetAbsoluteIndices(graph_data, p_IDs);

	// Number of elementary nanoparticle in the system
	int nPart = masses.size();

	// Total number of edges in the simulation
	int nEdges = getTotalEdges(this->aggregates);

	XMLDocument xmlDoc;

	XMLElement *pRoot = xmlDoc.NewElement("VTKFile");
	pRoot->SetAttribute("type", "PolyData");
	pRoot->SetAttribute("version", "0.1");
	xmlDoc.InsertFirstChild(pRoot);

	{
		// "Polydata" Sub Tree
		XMLElement *pPolyData = xmlDoc.NewElement("PolyData");

		{
			// "Piece" Sub Tree
			XMLElement *pPiece = xmlDoc.NewElement("Piece");
			// "NumberOfPoints" attribute, indicating the number of points (particles) to draw
			pPiece->SetAttribute("NumberOfPoints", nPart); // Particles in the simulation
			pPiece->SetAttribute("NumberOfLines", nEdges); // Edges in the simulation

			{
				// Points(Nanoparticles) coordinates
				XMLElement *pPoints = xmlDoc.NewElement("Points");

				{
					XMLElement *pDataArray = xmlDoc.NewElement("DataArray");
					pDataArray->SetAttribute("NumberOfComponents", 3); // 3 Dimensions
					pDataArray->SetAttribute("type", "Float32");
					pDataArray->SetAttribute("format", "ascii");

					std::string s_coords = V_To_S(coords);
					pDataArray->SetText(s_coords.c_str());


					pPoints->InsertEndChild(pDataArray);
				}

				pPiece->InsertEndChild(pPoints);
			}

			{
				XMLElement *pPointData = xmlDoc.NewElement("PointData");
				// Mass scalar for each point (nanoparticle)
				pPointData->SetAttribute("Scalars", "particles data");

				{
					// Masses
					XMLElement *pDataArray = xmlDoc.NewElement("DataArray");
					pDataArray->SetAttribute("type", "Float32");
					pDataArray->SetAttribute("Name", "mass");
					pDataArray->SetAttribute("format", "ascii");

					std::string s_masses = V_To_S(masses);
					pDataArray->SetText(s_masses.c_str());

					pPointData->InsertEndChild(pDataArray);
				}

				{
					// particles diameter
					XMLElement *pDataArray = xmlDoc.NewElement("DataArray");
					pDataArray->SetAttribute("type", "Float32");
					pDataArray->SetAttribute("Name", "particles diameter");
					pDataArray->SetAttribute("format", "ascii");

					std::string s_diameter = V_To_S(diameter);
					pDataArray->SetText(s_diameter.c_str());

					pPointData->InsertEndChild(pDataArray);
				}

				{
					// aggregates diameter
					XMLElement *pDataArray = xmlDoc.NewElement("DataArray");
					pDataArray->SetAttribute("type", "Float32");
					pDataArray->SetAttribute("Name", "aggregates diameter");
					pDataArray->SetAttribute("format", "ascii");

					std::string s_agg_diameter = V_To_S(agg_diameter);
					pDataArray->SetText(s_agg_diameter.c_str());

					pPointData->InsertEndChild(pDataArray);
				}

				{
					// aggregate ID in which are included
					XMLElement *pDataArray = xmlDoc.NewElement("DataArray");
					pDataArray->SetAttribute("type", "Int32");
					pDataArray->SetAttribute("Name", "aggregateId");
					pDataArray->SetAttribute("format", "ascii");

					std::string s_affil = V_To_S(affil);
					pDataArray->SetText(s_affil.c_str());

					pPointData->InsertEndChild(pDataArray);
				}

				{
					// particle ID
					XMLElement *pDataArray = xmlDoc.NewElement("DataArray");
					pDataArray->SetAttribute("type", "Int32");
					pDataArray->SetAttribute("Name", "particleId");
					pDataArray->SetAttribute("format", "ascii");

					std::string s_p_IDs = V_To_S(p_IDs);
					pDataArray->SetText(s_p_IDs.c_str());

					pPointData->InsertEndChild(pDataArray);
				}

				pPiece->InsertEndChild(pPointData);
			}

			// Edges connecting the particles
			////////////////////////////////////////////////////////////////////////
			XMLElement *pCellData = xmlDoc.NewElement("CellData");
			{
				// DataArray "edge_type" indicating the type of the edge, 0 physical, 1 contraints
				{
					XMLElement *pDataArray = xmlDoc.NewElement("DataArray");
					pDataArray->SetAttribute("type", "Int32");
					pDataArray->SetAttribute("Name", "edge_type");
					pDataArray->SetAttribute("format", "ascii");

					std::string s_types = V_To_S(abs_idx[1]);
					pDataArray->SetText(s_types.c_str());

					pCellData->InsertEndChild(pDataArray);
				}

				pPiece->InsertEndChild(pCellData);
			}

			XMLElement *pLines = xmlDoc.NewElement("Lines");
			{
				// DataArray "Connectivity", points indices couples for indicating the edge extremes
				{
					XMLElement *pDataArray = xmlDoc.NewElement("DataArray");
					pDataArray->SetAttribute("type", "Int32");
					pDataArray->SetAttribute("Name", "connectivity");
					pDataArray->SetAttribute("format", "ascii");

					std::string s_edges = V_To_S(abs_idx[0]);
					pDataArray->SetText(s_edges.c_str());

					pLines->InsertEndChild(pDataArray);
				}
				// DataArray "offset"
				{
					XMLElement *pDataArray = xmlDoc.NewElement("DataArray");
					pDataArray->SetAttribute("type", "Float32");
					pDataArray->SetAttribute("Name", "offsets");
					pDataArray->SetAttribute("format", "ascii");

					std::string s_types = V_To_S(abs_idx[2]);
					pDataArray->SetText(s_types.c_str());

					pLines->InsertEndChild(pDataArray);
				}
				pPiece->InsertEndChild(pLines);
			}


			//////////////////////////////////////////////////////////////////////////////////////////
			pPolyData->InsertEndChild(pPiece);
		}

		pRoot->InsertEndChild(pPolyData);
	}


	std::string file_name = root_path + "nanodome" + std::to_string(iter) + ".vtp";
	xmlDoc.SaveFile(file_name.c_str());



}

template<typename A>
void ConstrainedLangevinParticlePhase<A>::motion(double dt, GasModels* gp, double T) {

	double cube_side = pow(this->volume, 1. / 3.);

	// For each aggregate in the simulation
	for (auto& agg : this->aggregates) {

		// update the particle positions
		//agg.vsv_update_x_lgv(p.dt, gamma, p.T);
		agg->langevin_x_update(dt,gp,T);

		// update the velocity
		//agg.vsv_update_v_lgv(p.dt);
		agg->langevin_v_update(dt);

		// check reflections on the wall
		agg->wall_reflection(cube_side);
	}

	// apply aggregate constraints
	for (auto& a : this->aggregates)
		a->SHAKE();

	return;
}

template<typename A>
void ConstrainedLangevinParticlePhase<A>::motion_lc(double dt, GasModels* gp, double T) {

//	double cube_side = pow(this->volume, 1. / 3.);

	std::list<int> erased_agg;

	#ifdef OMP

	omp_set_num_threads(8);

	bool SHAKE_ERROR;

	#pragma omp parallel
	#pragma omp single
	{
		for (auto agg = this->aggregates.begin(); agg != this->aggregates.end(); ++agg) {
			#pragma omp task firstprivate(agg)
			// Update position
			(*agg)->langevin_x_update(dt, gp, T);
			// Update velocities
			(*agg)->langevin_v_update(dt);

			SHAKE_ERROR = (*agg)->SHAKE();
			// SHAKE Algorithm
			#pragma omp critical
			if (!(SHAKE_ERROR)) { // if the algorithm diverges, the aggregate is erased from the list
				std::cout << "SHAKE DIVERGED IN PARTICLES MOTION" << std::endl;
				erased_agg.push_back((*agg)->get_id());
				this->volume *= (this->aggregates.size() - 1.0) / this->aggregates.size();
			}
		}
		#pragma omp taskwait
	}
	#else
	for (auto agg = this->aggregates.begin(); agg != this->aggregates.end(); ++agg) {
		// Update position
		(*agg)->langevin_x_update(dt, gp, T);
		// Update velocities
		(*agg)->langevin_v_update(dt);
		// SHAKE Algorithm
		if (!((*agg)->SHAKE())) { // if the algorithm diverges, the aggregate is erased from the list
			std::cout << "SHAKE DIVERGED IN PARTICLES MOTION" << std::endl;
			erased_agg.push_back((*agg)->get_id());
			this->volume *= (this->aggregates.size() - 1.0) / this->aggregates.size();
		}
	}
	#endif


	// Erase the diverged aggregates from the grid and aggregates, if there are any
	if (erased_agg.size() > 0) {
		// Erase from aggregate
		for(auto e = erased_agg.begin(); e != erased_agg.end(); e++){
			for(auto a = this->aggregates.begin(); a !=this->aggregates.end(); a++){
				if((*a)->get_id() == (*e)){
					this->aggregates.erase(a);
					break;
				}
			}
			// Erase from grid
			for(int c = 0; c < this->grid.get_n_cells(); c++){
				for(int a = 0; a < this->grid[c].get_n(); a++){
					if(this->grid[c][a]->get_id() == (*e)){
						this->grid[c].erase_n(a);
						break;
					}
				}
			}
		}

		/*
			* int N_cells = this->grid.get_n_cells();
			* for (int c = 0; c < N_cells; c++) {
			*	int N_agg = this->grid[c].get_n();
			*	if (N_agg == 0) continue;
			*	for (int a = 0; a < N_agg; a++) {
			*		int agg_ID = this->grid[c][a]->get_id();
			*		for (auto idx = erased_agg.begin(); idx != erased_agg.end(); idx++) {
			*			if (agg_ID == (*idx)) {
			*	//test
			*	this->grid[c].erase_n(a);
			*	std::cout<<"aggregate "<<agg_ID <<" "<<(*idx)<<"erased from grid linear location: "<<a<<std::endl;
			*	if(std::next(idx) == erased_agg.end()){
			*	    erased_agg.erase(idx);
			*	    break;
	}
	else {
		std::advance(idx, 1);
		erased_agg.erase(std::prev(idx));

	}
	// update cell index and cardinality
	a--;
	N_agg--;
	}
	}
	}
	}
	// Test
	std::cout<<"Aggregate erased from grid"<<std::endl;
	*/
	}

	// update the grid
	int biggest_agg = -1;
	double biggest_s_enclosure = this->get_biggest_spherical_enclosure(biggest_agg);
	this->grid.grid_update(this->aggregates, this->volume, biggest_s_enclosure, biggest_agg);

	// update constraits
	this->update_constraints();

	return;

}

template<typename A>
double ConstrainedLangevinParticlePhase<A>::nucleation(double dt, double J, double j, SingleComponentComposition* s, NucleationTheory* nt) {

	// average number of nucleation events in this timestep
	double n_events = J * this->volume * dt;

	// single nucleation events
	int n = int(n_events);

	// fractional part
	double n_frac = n_events - n;

	// add all the new particles according to J
	for (int i = 0; i<n; ++i)
		this->add_particle(j, s, nt);

	// the fraction < 1 is treated in a probabilistic way
	if (ndm::uniform_double_distr(ndm::rand_gen) <= n_frac) {
		this->add_particle(j, s, nt);
		n++;
	}

	// molecules consuption per unit volume and time
	return n*j / (this->volume*dt);
}

template<typename A>
double ConstrainedLangevinParticlePhase<A>::nucleation_lc(double dt, double J, double j, SingleComponentComposition* s) {

	// average number of nucleation events in this timestep
	double n_events = J * this->volume * dt;

	// single nucleation events
	int n = int(n_events);

	// fractional part
	double n_frac = n_events - n;

	// add all the new particles according to J
	for (int i = 0; i<n; ++i)
		this->add_particle_lc(j, s);

	// the fraction < 1 is treated in a probabilistic way
	if (ndm::uniform_double_distr(ndm::rand_gen) <= n_frac) {
		this->add_particle_lc(j, s);
		n++;
	}

	// molecules consuption per unit volume and time
	return n*j / (this->volume*dt);


}

template<typename A>
void ConstrainedLangevinParticlePhase<A>::aggregates_number_balance() {

	int N = this->aggregates.size();

	if (N>this->max_aggregates_number) {

//	    while (N > this->max_aggregates_number) {
		auto it = this->aggregates.begin();
		double rho = ndm::uniform_double_distr(ndm::rand_gen);

		std::advance(it, int(rho*(N - 1)));

		// ID erased aggregate
		int eraased_ID = (*it)->get_id();
		// center of mass
		std::valarray<double> com = (*it)->get_center_of_mass();

		this->aggregates.erase(it);

		// Erase from the grid
		// Get the cell coordinates
		std::valarray<int> cell_coords = this->grid.get_cell_index(com);
		// erase from the grid
		auto aggregates_cell = this->grid(cell_coords[0], cell_coords[1], cell_coords[2]).get_list();
		int agg_idx = 0;
		for (auto agg = aggregates_cell.begin(); agg != aggregates_cell.end(); agg++) {
			if ((*agg)->get_id() == eraased_ID) {
				this->grid(cell_coords[0], cell_coords[1], cell_coords[2]).erase_n(agg_idx);
			}
			else
				agg_idx++;
		}

		// Update volume and aggregates number
		this->volume *= (N - 1.0) / N;
//		N--

//	    }
	}

	if (N < this->min_aggregates_number) {

		while (N != this->max_aggregates_number) {

			double side = pow(this->volume, 1.0 / 3.0);

			auto it = this->aggregates.begin();
			double rho = ndm::uniform_double_distr(ndm::rand_gen);

			std::advance(it, int(rho*(N - 1)));

			// Using the copy contructor, create a shared pointer for the aggregate copy
			auto agg = std::make_shared<RATTLEAggregate<DynamicParticle>>(*(*it));

			// Debug
			std::cout << "ID Aggregate to copy: "
			<< (*it)->get_id()
			<< "ID Copied Aggregate: "
			<< agg->get_id() << std::endl;

			// get the enclosing sphere radius
			double r0 = agg->get_enclosing_sphere_diameter()*0.5;

			// get volume
			double v0 = agg->get_volume();

			// while guard
			bool placement = false;
			// new center of mass coordinates
			std::valarray<double> n_c_mass;

			while (!placement) {
				// get a random new position to place the aggregate
				n_c_mass = {
					ndm::uniform_double_distr(ndm::rand_gen),
					ndm::uniform_double_distr(ndm::rand_gen),
					ndm::uniform_double_distr(ndm::rand_gen) };

					//center of mass scaled in the domain
					n_c_mass = side*(n_c_mass - 0.5);

					placement = true;

					// Check if in the new position, can collide or is superposed with some aggregate

					// get the cell where the aggregate can be placed
					std::valarray<int> cell_idx = this->grid.get_cell_index(n_c_mass);

					std::vector<std::shared_ptr<RATTLEAggregate<DynamicParticle>>> agg_cell_list;

					agg_cell_list = this->grid(cell_idx[0], cell_idx[1], cell_idx[2]).get_list();

					//try {
					//	// Get the aggregates in the cell
					//	agg_cell_list = this->grid(cell_idx[0], cell_idx[1], cell_idx[2]).get_list();
					//}
					//catch(std::exception& e) {
					//	std::cout << n_c_mass[0] << " " << n_c_mass[1] << " " << n_c_mass[2] << std::endl;
					//	std::cout << cell_idx[0] << " " << cell_idx[1] << " " << cell_idx[2] << std::endl;
					//}

					for (auto a1 = agg_cell_list.begin(); a1 != agg_cell_list.end(); ++a1) {
						double r1 = (*a1)->get_enclosing_sphere_diameter()*0.5;
						std::valarray<double> diff = n_c_mass - (*a1)->get_center_of_mass();
						double dist = sqrt((diff*diff).sum());
						// if the collision is possible, new center
						if (dist - r0 - r1 <= 1.0e-10) {
							std::cout << "POSSIBLE COLLISION IN DUPLICATING AGGREGATE, NEW CENTER OF MASS" << std::endl;
							std::cout << n_c_mass[0] << " " << n_c_mass[1] << " " << n_c_mass[2]<< std::endl;
							std::cout << "Volume: " << v0 << std::endl;
							std::cout << cell_idx[0] << " " << cell_idx[1] << " " << cell_idx[2]
							<< " cell side: " << this->grid.get_cell_side() << " volume side: " << side << " grid side: " << this->grid.get_grid_side()
							<< std::endl << "Volume_side - grid side = " << abs(side - this->grid.get_grid_side())
							<< "Cells: " << this->grid.get_C() << std::endl;

							placement = false;
							break;
						}
					}

			}

			// Change the aggregate's center of mass
			agg->change_center_of_mass(n_c_mass);

			// Add the the aggregates
			this->aggregates.push_back(agg);

			// Insert in the grid
			this->grid.insert_aggregate(agg);

			// Update volume
			this->volume *= (N + 1.0) / N;

			// Update the grid
			int biggest_agg = -1;
			double biggest_s_enclosure = this->get_biggest_spherical_enclosure(biggest_agg);
			this->grid.grid_update(this->aggregates, this->volume, biggest_s_enclosure, biggest_agg);

			N++;
		}

	}
}


template<typename A>
void ConstrainedLangevinParticlePhase<A>::coagulation() {


	for (auto a0 = this->aggregates.begin(); a0 != this->aggregates.end(); ++a0) {
		for (auto a1 = std::next(a0); a1 != this->aggregates.end(); ++a1) {
			std::shared_ptr<DynamicParticle> p0;
			std::shared_ptr<DynamicParticle> p1;
			bool collision = (*a0)->check_collision((*a1), 1.0e-10, p0, p1);
			if (collision) {
				// Merge aggregates
				(*a0)->merge((*a1), p0, p1);
				// Update actual particles cardinality inside the agregate (for sintering)
				(*a0)->update_particles_cardinality();
				// Update graph of the new aggregate
				(*a0)->update_graph();
				// Erase graph of the included aggregate
				(*a1)->erase_graph();
				#ifdef VERBOSE
				std::cout << "CHECK INDICES:" << std::endl;
				std::cout << "ID 0: " << (*a0)->get_id() << " ID 1: " << (*a1)->get_id() << std::endl;
				#endif

				// trick for deleting the current bond without crash
				if (a1 == this->aggregates.begin()) {
					this->aggregates.remove((*a1));
					break;
				}
				else {
					a1 = std::prev(a1);
					this->aggregates.remove((*(std::next(a1))));
				}

			}


		}
	}

	return;
}

template<typename A>
void ConstrainedLangevinParticlePhase<A>::coagulation_lc() {

	// typedef for shortening the code
	typedef std::shared_ptr<RATTLEAggregate<DynamicParticle>> Aggregate_ptr; // Aggregate's shared pointer
	typedef std::shared_ptr<DynamicParticle> Particle_ptr; // particle's shared pointer

	// event: First aggregate, Second Aggregate, contact particle 1, contact particle 2, shift for periodic boundaries,
	// 3D cell index(x, y, z) and index in the cell
	typedef std::tuple<Aggregate_ptr, Aggregate_ptr,
	Particle_ptr, Particle_ptr,
	std::valarray<double>,
	int, int, int, int> collision_event;


	// Use the Grid
	int n_cells = this->grid.get_C(); ///< number of cells for each dimension

	// domain side
	double d_side = pow(this->volume, 1 / 3);

	int x1, x2, x3;
	int c_x1, c_x2, c_x3; // Neighbourg cell coordinates
	int c_idx, nc_index;

	// Inspect the near Cells
	int x1_start, x2_start, x3_start;
	int x1_end, x2_end, x3_end;

	//	// Coagulation -> Wall Reflection
	//	for (x3 = 0; x3 < n_cells; x3++) { // z // Coagulation
	//		for (x2 = 0; x2 < n_cells; x2++) { // y
	//			for (x1 = 0; x1 < n_cells; x1++) { // x
	//
	//				for (c_idx = 0; c_idx < grid(x1, x2, x3).get_n(); c_idx++) {
	//					// Get the aggregate in the cell
	//					std::shared_ptr<RATTLEAggregate<DynamicParticle>> a0 = grid(x1, x2, x3)[c_idx];
	//					// if the aggregate is not empty
	//					if (a0->get_particles_number() == 0)
	//						continue;
	//
	//
	//						// Set the boundaries (make it better)
	//						(x1 - 1 < 0) ? x1_start = 0 : x1_start = x1 - 1;
	//						(x2 - 1 < 0) ? x2_start = 0 : x2_start = x2 - 1;
	//						(x3 - 1 < 0) ? x3_start = 0 : x3_start = x3 - 1;
	//						(x1 + 2 < n_cells) ? x1_end = x1 + 2 : x1_end = n_cells;
	//						(x2 + 2 < n_cells) ? x2_end = x2 + 2 : x2_end = n_cells;
	//						(x3 + 2 < n_cells) ? x3_end = x3 + 2 : x3_end = n_cells;
	//
	//						for (c_x3 = x3_start; c_x3 < x3_end; c_x3++) {
	//							for (c_x2 = x2_start; c_x2 < x2_end; c_x2++) {
	//								for (c_x1 = x1_start; c_x1 < x1_end; c_x1++) {
	//
	//									// all the aggregates in the near cell
	//									int cell_N = grid(c_x1, c_x2, c_x3).get_n();
	//									for (nc_index = 0; nc_index < cell_N; nc_index++) {
	//										std::shared_ptr<RATTLEAggregate<DynamicParticle>> a1 = grid(c_x1, c_x2, c_x3)[nc_index];
	//
	//										// if the aggregate is not empty
	//										if (a1->get_particles_number() == 0)
	//											continue;
	//
	//										std::shared_ptr<DynamicParticle> p0;
	//										std::shared_ptr<DynamicParticle> p1;
	//										bool collision = a0->check_collision(a1, 1.0e-10, p0, p1);
	//										if (collision) {
	//											// Merge aggregates
	//											a0->merge(a1, p0, p1);
	//											// Update actual particles cardinality inside the agregate (for sintering)
	//											a0->update_particles_cardinality();
	//											// Update graph of the new aggregate
	//											a0->update_graph();
	//											// Erase graph of the included aggregate
	//											a1->erase_graph();
	//											// Erase the aggregate from the cell
	//											grid(c_x1, c_x2, c_x3).erase_n(nc_index);
	//											// Update indices
	//											nc_index--;
	//											cell_N--;
	//
	//#ifdef VERBOSE
	//											std::cout << "CHECK INDICES:" << std::endl;
	//											std::cout << "ID 0: " << a0->get_id()
	//													  << " ID 1: " << a1->get_id() << std::endl;
	//#endif
	//
	//										}
	//
	//									}
	//								}
	//
	//							}
	//						}
	//
	//				}
	//			}
	//		}
	//
	//	}

	// Coagulation periodic boundaries

	//	for (x3 = 0; x3 < n_cells; x3++) { // z // Coagulation
	//		for (x2 = 0; x2 < n_cells; x2++) { // y
	//			for (x1 = 0; x1 < n_cells; x1++) { // x
	//
	//                for (c_idx = 0; c_idx < this->grid(x1, x2, x3).get_n(); c_idx++) {
	//					// Get the aggregate in the cell
	//                    std::shared_ptr<RATTLEAggregate<DynamicParticle>> a0 = this->grid(x1, x2, x3)[c_idx];
	//					// if the aggregate is not empty due to collision with other aggregate
	//					if (a0->get_particles_number() == 0)
	//						continue;
	//
	//
	//					/* FULL NEIGHBOURGS
	//					x3_start = x3 - 1; x3_end = x3 + 1;
	//					x2_start = x2 - 1; x2_end = x2 + 1;
	//					x1_start = x1 - 1; x1_end = x1 + 1;
	//					*/
	//					// HALF NEIGHBOURGS
	//					x3_start = x3 - 1; x3_end = x3 + 1;
	//					x2_start = x2	 ; x2_end = x2 + 1;
	//					x1_start = x1 - 1; x1_end = x1 + 1;
	//
	//					for (c_x3 = x3_start; c_x3 < x3_end; c_x3++) {
	//						for (c_x2 = x2_start; c_x2 < x2_end; c_x2++) {
	//							for (c_x1 = x1_start; c_x1 < x1_end; c_x1++) {
	//
	//								std::valarray<double> shift = { 0.0,0.0,0.0 };
	//
	//								// Get the periodic position and the shift for the coordinats
	//								int p_c_x1 = periodic_position(c_x1, n_cells, shift[0], d_side);
	//								int p_c_x2 = periodic_position(c_x2, n_cells, shift[1], d_side);
	//								int p_c_x3 = periodic_position(c_x3, n_cells, shift[2], d_side);
	//
	//								// all the aggregates in the near cell
	//                                int cell_N = this->grid(p_c_x1, p_c_x2, p_c_x3).get_n();
	//								for (nc_index = 0; nc_index < cell_N; nc_index++) {
	//                                    auto a1 = this->grid(p_c_x1, p_c_x2, p_c_x3)[nc_index];
	//
	//									// if the aggregate is not empty due to collision with other aggregate
	//									if (a1->get_particles_number() == 0)
	//										continue;
	//
	//									std::shared_ptr<DynamicParticle> p0;
	//									std::shared_ptr<DynamicParticle> p1;
	//									bool collision = a0->check_collision_periodic(a1, 1.0e-10, shift, p0, p1);
	//									if (collision) {
	//										std::cout << "COLLISION!" << std::endl;
	//										// Merge aggregates
	//										a0->merge_periodic(a1, shift, p0, p1);
	//										// Update actual particles cardinality inside the aggregate (for sintering)
	//										a0->update_particles_cardinality();
	//										// Update graph of the new aggregate
	//										a0->update_graph();
	//										// Erase graph of the included aggregate
	//										a1->erase_graph();
	//										// Erase the aggregate from the cell
	//                                        this->grid(p_c_x1, p_c_x2, p_c_x3).erase_n(nc_index);
	//										// Update indices
	//										nc_index--;
	//										cell_N--;
	//#ifdef VERBOSE
	//										std::cout << "CHECK INDICES:" << std::endl;
	//										std::cout << "ID 0: " << a0->get_id()
	//											<< " ID 1: " << a1->get_id() << std::endl;
	//#endif
	//
	//									}
	//
	//								}
	//							}
	//
	//						}
	//					}
	//
	//				}
	//			}
	//		}
	//
	//	}
	//
	//
	//	/*!!!!!!!!!!!!!!!! DO IT BETTER!!!!!!!!!!!!!!!!!!!!!!*/
	//	// Scrub the coagulated aggregates and update volume
	//	for (auto agg = this->aggregates.begin(); agg != this->aggregates.end(); agg++) {
	//		if ((*agg)->get_particles_number() == 0) {
	//			if (agg == this->aggregates.begin()) {
	//				this->aggregates.remove((*agg));
	//
	//				int N = this->aggregates.size();
	//				if (N > 0)
	//					this->volume *= (N - 1.0) / N;
	//			}
	//			else {
	//				agg = std::prev(agg);
	//				this->aggregates.remove((*(std::next(agg))));
	//
	//				int N = this->aggregates.size();
	//				if (N > 0)
	//					this->volume *= (N - 1.0) / N;
	//			}
	//		}
	//
	//	}

	// Coagulation with coagulation list and periodic boundaries

	// collisions list
	std::list<collision_event> collision_List;

	// unordered set for avoiding multiple collisions
	std::unordered_set<int> coagulated_aggrgates_IDs;

	for (x3 = 0; x3 < n_cells; x3++) { // z // Coagulation
		for (x2 = 0; x2 < n_cells; x2++) { // y
			for (x1 = 0; x1 < n_cells; x1++) { // x

				for (c_idx = 0; c_idx < this->grid(x1, x2, x3).get_n(); c_idx++) {
					// Get the aggregate in the cell
					Aggregate_ptr a0 = this->grid(x1, x2, x3)[c_idx];

					// HALF NEIGHBOURGS
					x3_start = x3 - 1; x3_end = x3 + 1;
					x2_start = x2; x2_end = x2 + 1;
					x1_start = x1 - 1; x1_end = x1 + 1;

					for (c_x3 = x3_start; c_x3 < x3_end; c_x3++) {
						for (c_x2 = x2_start; c_x2 < x2_end; c_x2++) {
							for (c_x1 = x1_start; c_x1 < x1_end; c_x1++) {

								std::valarray<double> shift = { 0.0,0.0,0.0 };

								// Get the periodic position and the shift for the coordinats
								int p_c_x3 = periodic_position(c_x3, n_cells, shift[2], d_side);
								int p_c_x2 = periodic_position(c_x2, n_cells, shift[1], d_side);
								int p_c_x1 = periodic_position(c_x1, n_cells, shift[0], d_side);

								// all the aggregates in the near cell
								int cell_N = this->grid(p_c_x1, p_c_x2, p_c_x3).get_n();
								for (nc_index = 0; nc_index < cell_N; nc_index++) {
									Aggregate_ptr a1 = this->grid(p_c_x1, p_c_x2, p_c_x3)[nc_index];

									Particle_ptr p0;
									Particle_ptr p1;
									bool collision = a0->check_collision_periodic(a1, 1.0e-10, shift, p0, p1);

									if (collision) {

										int ID0 = a0->get_id();
										int ID1 = a1->get_id();

										if ((coagulated_aggrgates_IDs.find(ID0) == coagulated_aggrgates_IDs.end()) &&
											(coagulated_aggrgates_IDs.find(ID1) == coagulated_aggrgates_IDs.end())) {

											std::cout << "COLLISION EVENT!" << std::endl;
										collision_event evt = std::make_tuple(a0, a1, p0, p1, shift, p_c_x1, p_c_x2, p_c_x3, nc_index);
										collision_List.push_back(evt);

										// update coagulated aggregates ID
										coagulated_aggrgates_IDs.insert(ID0);
										coagulated_aggrgates_IDs.insert(ID1);

										std::cout << "CHECK INDICES:" << std::endl;
										std::cout << "ID 0: " << a0->get_id()
										<< " ID 1: " << a1->get_id() << std::endl;

											}

											//// Merge aggregates
											//a0->merge_periodic(a1, shift, p0, p1);
											//// Update actual particles cardinality inside the aggregate (for sintering)
											//a0->update_particles_cardinality();
											//// Update graph of the new aggregate
											//a0->update_graph();
											//// Erase graph of the included aggregate
											//a1->erase_graph();
											//// Erase the aggregate from the cell
											//this->grid(p_c_x1, p_c_x2, p_c_x3).erase_n(nc_index);
											//// Update indices
											//nc_index--;
											//cell_N--;
											#ifdef VERBOSE
											std::cout << "CHECK INDICES:" << std::endl;
											std::cout << "ID 0: " << a0->get_id()
											<< " ID 1: " << a1->get_id() << std::endl;
											#endif

									}

								}
							}

						}
					}

				}
			}
		}

	}



	// Update the colliding aggregates
	// event: First aggregate, Second Aggregate, contact particle 1, contact particle 2, shift for periodic boundaries,
	// 3D cell index(x, y, z) and index in the cell
	for (auto c_evt = collision_List.begin(); c_evt != collision_List.end(); c_evt++) {
		// Unpack the tuple
		auto a0 = std::get<0>((*c_evt)); auto a1 = std::get<1>((*c_evt));
		auto p0 = std::get<2>((*c_evt)); auto p1 = std::get<3>((*c_evt));
		auto shift = std::get<4>((*c_evt));
		auto c_x1 = std::get<5>((*c_evt)); auto c_x2 = std::get<6>((*c_evt)); auto c_x3 = std::get<7>((*c_evt));
		auto nc_index = std::get<8>((*c_evt));
		std::cout << "Contact Particles: " << p0->get_id() << " - " << p1->get_id()<<std::endl;

		// Merge aggregates
		a0->merge_periodic(a1, shift, p0, p1);
		// Update actual particles cardinality inside the aggregate (for sintering)
		a0->update_particles_cardinality();
		// Update graph of the new aggregate
		a0->update_graph();
		// Erase graph of the included aggregate
		a1->erase_graph();
		// Erase the aggregate a1 from the cell
		std::cout << "aggregate to erase from the grid:" << this->grid(c_x1, c_x2, c_x3)[nc_index]->get_id() << std::endl;
		this->grid(c_x1, c_x2, c_x3).erase_n(nc_index);


		// Erase a1 from the aggregates list and update volume
		for (auto it = this->aggregates.begin(); it != this->aggregates.end(); it++) {
			if ((*it)->get_id() == a1->get_id()) {
				// erase aggregate
				this->aggregates.erase(it);
				// update volume
				int N = this->aggregates.size();
				if (N > 0)
					this->volume *= (N - 1.0) / N;
				// exit from the cycle
				break;
			}
		}

	}

	// Update the grid
	int biggest_agg = -1;
	double biggest_s_enclosure = this->get_biggest_spherical_enclosure(biggest_agg);
	this->grid.grid_update(this->aggregates, this->volume, biggest_s_enclosure, biggest_agg);

	// update constraits
	this->update_constraints();

}

template<typename A>
int ConstrainedLangevinParticlePhase<A>::get_aggregates_cardinality() {
	int card = 0;
	for (auto a = this->aggregates.begin(); a != this->aggregates.end(); a++) {
		if ((*a)->get_particles_number() != 0)
			card++;
	}

	return card;
}

/*
	/ /* Auxiliary Function to get all the coordinates of each particle in the system
	std::vector<double> getParticleCoords(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble) {

	std::vector<double> points;

	for (auto agg = ensamble.begin(); agg != ensamble.end(); ++agg) {

		std::list<std::shared_ptr<DynamicParticle>> lp_tmp = (*agg)->get_particles();
		for (auto p = lp_tmp.begin(); p != lp_tmp.end(); ++p) {
			std::valarray<double> cTmp = (*p)->get_x();
			for (int d = 0; d < 3; d++)
				points.push_back(cTmp[d]);
			}

			}

			return points;

			}

			// Auxiliary Function to get the mass of each particle in the system
			std::vector<double> getParticleMass(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble) {

			std::vector<double> masses;

			for (auto agg = ensamble.begin(); agg != ensamble.end(); ++agg) {

				std::list<std::shared_ptr<DynamicParticle>> lp_tmp = (*agg)->get_particles();
				for (auto p = lp_tmp.begin(); p != lp_tmp.end(); ++p) {
					masses.push_back((*p)->get_mass());
					}
					}

					return masses;
					}

					// Auxiliary Function to get the ID of each nanoparticle in the system
					std::vector<int> getParticleAggregate(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble) {

					//Aggregate affiliation of particles
					std::vector<int> aggAffil;

					//Aggregate's affiliation id
					int ident;

					for (auto agg = ensamble.begin(); agg != ensamble.end(); ++agg) {
						ident = (*agg)->get_id();
						std::list<std::shared_ptr<DynamicParticle>> lp_tmp = (*agg)->get_particles();
						for (auto p = lp_tmp.begin(); p != lp_tmp.end(); ++p) {
							aggAffil.push_back(ident);
							}
							}
							return aggAffil;


							}

							// Auxiliary Function to retrieve the ID of each particle in the aggregate
							std::vector<int> getParticleID(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble) {

							std::vector<int> IDs;

							for (auto agg = ensamble.begin(); agg != ensamble.end(); ++agg) {
								std::vector<int> p_IDs;
								std::list<std::shared_ptr<DynamicParticle>> lp_tmp = (*agg)->get_particles();
								for (auto p = lp_tmp.begin(); p != lp_tmp.end(); ++p) {
									p_IDs.push_back((*p)->get_id());
									}
									IDs.insert(IDs.end(), p_IDs.begin(), p_IDs.end());
									}


									return IDs;

									}

									// Auxiliary Function to get the diameter of each nanoparticle in the system
									std::vector<double> getParticleDiameter(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble) {

									std::vector<double> diameter;

									for (auto agg = ensamble.begin(); agg != ensamble.end(); ++agg) {
										std::list<std::shared_ptr<DynamicParticle>> lp_tmp = (*agg)->get_particles();
										for (auto p = lp_tmp.begin(); p != lp_tmp.end(); ++p) {
											diameter.push_back((*p)->get_diameter());
											}
											}

											return diameter;
											}

											std::vector<double>	getAggregatesDiameter(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble) {

											std::vector<double> diameter;

											for (auto agg = ensamble.begin(); agg != ensamble.end(); ++agg) {
												std::list<std::shared_ptr<DynamicParticle>> lp_tmp = (*agg)->get_particles();
												for (auto p = lp_tmp.begin(); p != lp_tmp.end(); ++p) {
													diameter.push_back((*agg)->get_enclosing_sphere_diameter());
													}
													}

													return diameter;


													}

													int getTotalEdges(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble) {

													int tot_edges = 0;

													for (auto agg = ensamble.begin(); agg != ensamble.end(); agg++) {
														tot_edges += (*agg)->get_graph().get_n_edges();
														}

														return tot_edges;

														}

														// Auxiliary function to get the connectiob data of each aggregate
														std::vector<std::tuple<int, int, int>> getConnections(std::list < std::shared_ptr<RATTLEAggregate<DynamicParticle>>>& ensamble) {

														std::vector<std::tuple<int, int, int>> graph_data;

														std::vector<int> couples;
														std::vector<int> types;
														std::vector<int> offsets;

														int off = 0;

														for (auto agg = ensamble.begin(); agg != ensamble.end(); agg++) {

															// Get graph edges (extreme indices and types)
															std::vector<std::tuple<int,int,int>> tmp = (*agg)->get_graph().get_edges_indices();
															graph_data.insert(graph_data.end(), tmp.begin(), tmp.end());
															}

															return graph_data;


															}

*/

// Auxiliary function to map edges extremes indices with position in the array for vtk
std::vector<std::vector<int>> GetAbsoluteIndices(std::vector<std::tuple<int, int, int>> graph_data, std::vector<int> p_list) {

  std::vector<std::vector<int>> res;

  std::vector<int> couples;
  std::vector<int> types;
  std::vector<int> offsets;

  int off = 0;

  for (auto t = graph_data.begin(); t != graph_data.end(); t++) {

      //get first extreme
      int idx1 = std::get<0>((*t));
      int idx_abs1 = -1;
      //get second extreme
      int idx2 = std::get<1>((*t));
      int idx_abs2 = -1;
      //get type
      int type = std::get<2>((*t));
      // find the first index
      for (int i = 0; i < p_list.size(); i++) {
          if (p_list[i] == idx1) {
              idx_abs1 = i;
              off++;
              break;
            }
        }
      // find the second index
      for (int i = 0; i < p_list.size(); i++) {
          if (p_list[i] == idx2) {
              idx_abs2 = i;
              off++;
              break;
            }
        }
      // insert the absolute indices
      couples.push_back(idx_abs1); couples.push_back(idx_abs2);
      // insert the edge type
      types.push_back(type);
      // give offset
      offsets.push_back(off);
    }

  res.push_back(couples); res.push_back(types); res.push_back(offsets);

  return res;

}

template<class dataType>
std::string V_To_S(std::vector<dataType> vect) {

  std::string ret_str;

  for (auto it : vect) {

      std::stringstream ss;

      ss << std::setprecision(4) << it;

      std::string part;

      ss >> part;
      ret_str += part + " ";
    }

  return ret_str;
}

#endif // CONSTRAINEDLANGEVINPARTICLEPHASE_H
