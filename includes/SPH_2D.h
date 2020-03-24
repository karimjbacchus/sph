#pragma once

#include <vector>
#include <cmath>
#include <iostream>

using namespace std;

class SPH_main;

class SPH_particle
{
public:
	// constructor
	SPH_particle()
	{
		// by default a new particle is not a boundary particle
		boundary = false;
		is_boundary_close = false;

		// set initial velocity to be zero
		for (int i = 0; i < 2; i++){
			v[i] = 0;
			v_h[i] = 0;
		}
	}

	//position and velocity
	double x[2], v[2], x_h[2], v_h[2];

	//density and pressure
	double rho, P, rho_h, P_h;

	// derivative of density and acceleration, and mass
	double drho, a[2], m, drho_h, a_h[2], m_h;

	// variables for density smoothing
	double sum_W, sum_W_rho;

	//link to SPH_main class so that it can be used in calc_index
	static SPH_main *main_data;

	//index in neighbour finding array
	int list_num[2];

	// for XSPH algorithm 
	double dv[2];	

	bool boundary;

	bool is_boundary_close;

	void calc_index(void);
};


class SPH_main
{
public:
	SPH_main();

	void set_values(void);
	void initialise_grid(bool shoaling = false);

	void place_points(double *min, double *max, bool shoaling = false);
	// Start with a wave beginning in the centre of the domain
	void place_points_conf1(double* min, double* max, bool shoaling = false);
	// Add a rectangular obstacle to the right boundary
	void place_points_conf2(double* min, double* max, bool shoaling = false);
	// Start with a "perfect" wave and add a step to the right boundary
	void place_points_conf3(double* min, double* max, bool shoaling = false);
	// Increase the steepness of the right step
	void place_points_conf4(double* min, double* max, bool shoaling = false);
	// Same level everywhere apart from a hole in the middle
	void place_points_conf5(double* min, double* max, bool shoaling = false);

	//allocates all the points to the search grid (assumes that index has been appropriately updated)
	void allocate_to_grid();

	// calculate W(r_ij,h)
	double compute_W(double dist, double h);

	void neighbour_iterate(int i_c, int j_c, bool den_smoothing, bool half=false); 

	// update masses of each particled
	void initialise_values(SPH_particle *part);

	void update_mass_pressure(SPH_particle *part, bool shoaling = false, bool half=false);

	// add elements of ns equation for each particle
	void add_to_navier_stokes(SPH_particle *part, SPH_particle *other_part, double dn[2], double dist, bool half=false);

	// time-stepping using Euler method
	void Euler_time_stepping(SPH_particle *part);
	
	//time-stepping using PC method
	void PC_half_stepping(SPH_particle *part);
	void PC_full_stepping(SPH_particle *part);

	// Enforce boundary conditions
	void remove_close_particles(bool half = false);

	//smoothing length
	double h;
	double h_fac;
	//particle initial spacing
	double dx;
	//time stepping
	double dt;
	// particle initial density
	double rho0;
	// speed of sound
	double c0;
	// gamma
	double gamma;
	// B
	double B;

	//dimensions of simulation region
	double min_x[2], max_x[2];

	int max_list[2];

	//list of all the particles
	vector<SPH_particle> particle_list;

	//Outer 2 are the grid, inner vector is the list of pointers in each cell
	vector<vector<vector<SPH_particle*> > > search_grid;
};
