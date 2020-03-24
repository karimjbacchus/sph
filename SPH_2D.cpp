#include "../includes/SPH_2D.h"

#ifndef M_PI
#define M_PI 3.14159265
#endif
#define G_CONST 9.80665
#define mu 0.001

SPH_main *SPH_particle::main_data;

void SPH_particle::calc_index(void)
{
	for (int i = 0; i < 2; i++)
		list_num[i] = int((x[i] - main_data->min_x[i]) / (2.0 * main_data->h));
}

SPH_main::SPH_main()
{
	SPH_particle::main_data = this;
}

void SPH_main::set_values(void)
{
	min_x[0] = 0.0;
	min_x[1] = 0.0;

	// 20 x 10 domain of fluid
	max_x[0] = 20.0;
	max_x[1] = 10.0;

	dx = 0.2;
	c0 = 20.;
	h_fac = 1.3;
	h = dx * h_fac;
	dt = 0.1 * h / c0;
	rho0 = 1000.0;
	gamma = 7.0;
	B = rho0 * c0 * c0 / gamma;
}

void SPH_main::initialise_grid(bool shoaling)
{
	for (int i = 0; i < 2; i++)
	{
		//add buffer for virtual wall particles
		min_x[i] -= 2.0*h;
		max_x[i] += 2.0*h;												

		max_list[i] = int((max_x[i] - min_x[i]) / (2.0*h) + 1.0);
	}

	search_grid.resize(max_list[0]);
	for (int i = 0; i < max_list[0]; i++)
		search_grid[i].resize(max_list[1]);
}


// do this at the beginning of each timestep: update masses and zero out derivatives
void SPH_main::initialise_values(SPH_particle *part)
{
	// initial velocity and density
	part->rho = this->rho0;
	part->v[0] = 0.0;
	part->v[1] = 0.0;
	part->rho_h = this->rho0;
	part->v_h[0] = 0.0;
	part->v_h[1] = 0.0;

	// zero out derivatives in preparation for recalculation in neighbour_iterate
	part->a[0] = 0.0;
	part->a[1] = -G_CONST;
	part->drho = 0.0;
	part->a_h[0] = 0.0;
	part->a_h[1] = -G_CONST;
	part->drho_h = 0.0;
	
	// mass and pressure equation
	part->m = this->dx * this->dx * part->rho;
	part->P = this->B * (pow(part->rho / this->rho0, this->gamma) - 1.0);
	part->m_h = this->dx * this->dx * part->rho_h;
	part->P_h = this->B * (pow(part->rho_h / this->rho0, this->gamma) - 1.0);

	// density smoothing
	part->sum_W = 0.0;
	part->sum_W_rho = 0.0;
}

void SPH_main::update_mass_pressure(SPH_particle *part, bool shoaling, bool half)
{
	// zero out derivatives in preparation for recalculation in neighbour_iterate
	if (!half)
	{
		part->a[0] = 0.0;
		part->a[1] = -G_CONST;
		part->drho = 0.0;
	}
	else
	{
		part->a_h[0] = 0.0;
		part->a_h[1] = -G_CONST;
		part->drho_h = 0.0;
	}

	// mass and pressure equation
	if (!half)
	{
		part->m = this->dx * this->dx * part->rho;
		part->P = this->B * (pow(part->rho / this->rho0, this->gamma) - 1.0);
	}
	else
	{
		part->m_h = this->dx * this->dx * part->rho_h;
		part->P_h = this->B * (pow(part->rho_h / this->rho0, this->gamma) - 1.0);
	}
	
	// if (shoaling){
	// 	if (part->x[1] < 0.1*part->x[0] && part->boundary == false){
	// 		double below_dist = part->x[1] - (0.1*part->x[0] - 2.0 * h);
	// 		part->a[1] += 0.1*( pow(2.0 * h / below_dist, 12) - 2*pow(2.0 * h / below_dist, 6) );
	// 	}
	// }

	//force to contain from bottom
	if (!half){
		if (part->x_h[1] <= 0 && part->boundary == false){
			double below_dist = 2.0 * h + part->x_h[1];
			part->a[1] += 2*( pow(2.0 * h / below_dist, 12) - 2*pow(2.0 * h / below_dist, 6) );
		}

		// force to contain from left
		if (part->x_h[0] <= 0 && part->boundary == false){
			double below_dist = 2.0 * h + part->x_h[0];
			if (below_dist > 0)
				part->a[0] += 0.5*( pow(2.0 * h / below_dist, 12) - 2*pow(2.0 * h / below_dist, 6) );
		}

		// force to contain from right
		if (part->x_h[0] >= (this->max_x[0]-2.0 * h) && part->boundary == false){
			double below_dist = this->max_x[0] - part->x_h[0];
			if (below_dist > 0)
				part->a[0] -= 0.5*( pow(2.0 * h / below_dist, 12) - 2*pow(2.0 * h / below_dist, 6) );
		}
	}

	else{
		if (part->x[1] <= 0 && part->boundary == false){
			double below_dist = 2.0 * h + part->x[1];
			part->a_h[1] += 2*( pow(2.0 * h / below_dist, 12) - 2*pow(2.0 * h / below_dist, 6) );
		}

		// force to contain from left
		if (part->x[0] <= 0 && part->boundary == false){
			double below_dist = 2.0 * h + part->x[0];
			if (below_dist > 0)
				part->a_h[0] += 0.5*( pow(2.0 * h / below_dist, 12) - 2*pow(2.0 * h / below_dist, 6) );
		}

		// force to contain from right
		if (part->x[0] >= (this->max_x[0]-2.0 * h) && part->boundary == false){
			double below_dist = this->max_x[0] - part->x[0];
			if (below_dist > 0)
				part->a_h[0] -= 0.5*( pow(2.0 * h / below_dist, 12) - 2*pow(2.0 * h / below_dist, 6) );
		}
	}
	
	// density smoothing
	part->sum_W = 0.0;
	part->sum_W_rho = 0.0;

	part->calc_index();
}

// Places initial points - will need to be modified to include boundary points 
// and the specifics of where the fluid is in the domain
void SPH_main::place_points(double* min, double* max, bool shoaling)
{
	double x[2] = { min[0], min[1] };
	SPH_particle particle;

	for (x[0] = min[0]; x[0] <= max[0]; x[0] += dx) // across x axis
	{
		for (x[1] = min[1]; x[1] <= max[1]; x[1] += dx) // accross y axis
		{
			// skip empty space of domain
			if (x[0] >= 0 && x[0] < max[0] - 2.0 * h) // left/right boundaries
			{
				if (x[1] < max[1] - 2.0 * h) // top boundary 
				{
					// empty top rectangle
					if (x[1] > 5)
						continue;
					// empty lower right rectangle
					else if (x[1] > 2 && x[0] > 3)
						continue;
				}
			}

			for (int i = 0; i < 2; i++)
			{
				particle.x[i] = x[i];
				particle.x_h[i] = x[i];
			}

			// check if particle is a boundary point
			if (x[0] < 0 || x[1] < 0 || x[0] > max[0] - 2.0 * h || x[1] > max[1] - 2.0 * h)
				particle.boundary = true;
			else
				particle.boundary = false;

			if (shoaling) {
				if (x[1] < 0.1 * x[0])
					particle.boundary = true;
			}

			particle.calc_index();

			this->initialise_values(&particle);

			particle_list.push_back(particle);
		}
	}


}

// Start with a wave beginning in the centre of the domain
void SPH_main::place_points_conf1(double *min, double *max, bool shoaling)
{
	double x[2] = { min[0], min[1] };
	SPH_particle particle;

	for (x[0] = min[0]; x[0] <= max[0]; x[0] += dx) // across x axis
	{
		for (x[1] = min[1]; x[1] <= max[1]; x[1] += dx) // accross y axis
		{
			// skip empty space of domain
			if (x[0] >= 0 && x[0] < max[0] - 2.0*h) // left/right boundaries
			{ 
				if (x[1] < max[1] - 2.0*h) // top boundary 
				{ 
					if (x[1] > 9)
						continue; 
					else if (x[1] > 2 && (x[0] < 9 || x[0] > 11))
						continue;
				}
			}

			for (int i = 0; i < 2; i++)
			{
				particle.x[i] = x[i];
				particle.x_h[i] = x[i];
			}

			// check if particle is a boundary point
			if (x[0] < 0 || x[1] < 0 || x[0] > max[0] - 2.0*h || x[1] > max[1] - 2.0*h)
				particle.boundary = true;
			else
				particle.boundary = false;

			if (shoaling){
				if (x[1] < 0.1*x[0])
					particle.boundary = true;
			}

			particle.calc_index();
			
			this->initialise_values(&particle);

			particle_list.push_back(particle);
		}
	}
	
	
}

// Add a rectangular obstacle to the right boundary
void SPH_main::place_points_conf2(double* min, double* max, bool shoaling)
{
	double x[2] = { min[0], min[1] };
	SPH_particle particle;

	for (x[0] = min[0]; x[0] <= max[0]; x[0] += dx) // across x axis
	{
		for (x[1] = min[1]; x[1] <= max[1]; x[1] += dx) // accross y axis
		{
			// skip empty space of domain
			if (x[0] >= 0 && x[0] < max[0] - 2.0 * h) // left/right boundaries
			{
				if (x[1] < max[1] - 2.0 * h) // top boundary 
				{
					if (x[1] > 5)
						continue;
					else
					{
						if (x[1] > 2 && x[0] > 3 && x[0] < 12)
							continue;
						if (x[1] > 3 && x[0] >= 12)
							continue;
					}
				}
			}

			for (int i = 0; i < 2; i++)
			{
				particle.x[i] = x[i];
				particle.x_h[i] = x[i];
			}

			// check if particle is a boundary point
			if (x[0] < 0 || x[1] < 0 || x[0] > max[0] - 2.0 * h || x[1] > max[1] - 2.0 * h || (x[0] >= 12 && x[1] > 1 && x[1] < 3))
				particle.boundary = true;
			else
				particle.boundary = false;

			if (shoaling) {
				if (x[1] < 0.1 * x[0])
					particle.boundary = true;
			}

			particle.calc_index();

			this->initialise_values(&particle);

			particle_list.push_back(particle);
		}
	}


}

// Start with a "perfect" wave and add a step to the right boundary
void SPH_main::place_points_conf3(double* min, double* max, bool shoaling)
{
	double x[2] = { min[0], min[1] };
	SPH_particle particle;

	for (x[0] = min[0]; x[0] <= max[0]; x[0] += dx) // across x axis
	{
		for (x[1] = min[1]; x[1] <= max[1]; x[1] += dx) // accross y axis
		{
			// skip empty space of domain
			if (x[0] >= 0 && x[0] < max[0] - 2.0 * h) // left/right boundaries
			{
				if (x[1] < max[1] - 2.0 * h) // top boundary 
				{
					if (x[1] > 7)
						continue;
					if (x[0] > 13 && x[1] > (x[0] - 10))
						continue;
					if (x[1] > 3 && x[0] > 6 && x[0] <= 13)
						continue;
					if (x[1] > 5 && x[0] <= 6 && x[0] > 3)
						continue;
					if (x[1] > 1 && x[0] <= 3)
						continue;

				}
			}

			for (int i = 0; i < 2; i++)
			{
				particle.x[i] = x[i];
				particle.x_h[i] = x[i];
			}

			// check if particle is a boundary point
			if (x[0] < 0 || x[1] < 0 || x[0] > max[0] - 2.0 * h || x[1] > max[1] - 2.0 * h || (x[0] > 10 && x[1] < (x[0] - 10)))
				particle.boundary = true;
			else
				particle.boundary = false;

			if (shoaling) {
				if (x[1] < 0.1 * x[0])
					particle.boundary = true;
			}

			particle.calc_index();

			this->initialise_values(&particle);

			particle_list.push_back(particle);
		}
	}


}

// Increase the steepness of the right step
void SPH_main::place_points_conf4(double* min, double* max, bool shoaling)
{
	double x[2] = { min[0], min[1] };
	SPH_particle particle;

	for (x[0] = min[0]; x[0] <= max[0]; x[0] += dx) // across x axis
	{
		for (x[1] = min[1]; x[1] <= max[1]; x[1] += dx) // accross y axis
		{
			// skip empty space of domain
			if (x[0] >= 0 && x[0] < max[0] - 2.0 * h) // left/right boundaries
			{
				if (x[1] < max[1] - 2.0 * h) // top boundary 
				{
					if (x[1] > 6)
						continue;
					if (x[0] > 16.5 && x[1] > 2*(x[0] - 15))
						continue;
					if (x[1] > 3 && x[0] > 6 && x[0] <= 16.5)
						continue;
					if (x[1] > 5 && x[0] <= 6 && x[0] > 3)
						continue;
					if (x[1] > 1 && x[0] <= 3)
						continue;

				}
			}

			for (int i = 0; i < 2; i++)
			{
				particle.x[i] = x[i];
				particle.x_h[i] = x[i];
			}

			// check if particle is a boundary point
			if (x[0] < 0 || x[1] < 0 || x[0] > max[0] - 2.0 * h || x[1] > max[1] - 2.0 * h || (x[0] > 15 && x[1] < 2*(x[0] - 15)))
				particle.boundary = true;
			else
				particle.boundary = false;

			if (shoaling) {
				if (x[1] < 0.1 * x[0])
					particle.boundary = true;
			}

			particle.calc_index();

			this->initialise_values(&particle);

			particle_list.push_back(particle);
		}
	}


}

// Same level everywhere apart from a hole in the middle
void SPH_main::place_points_conf5(double* min, double* max, bool shoaling)
{
	double x[2] = { min[0], min[1] };
	SPH_particle particle;

	for (x[0] = min[0]; x[0] <= max[0]; x[0] += dx) // across x axis
	{
		for (x[1] = min[1]; x[1] <= max[1]; x[1] += dx) // accross y axis
		{
			// skip empty space of domain
			if (x[0] >= 0 && x[0] < max[0] - 2.0 * h) // left/right boundaries
			{
				if (x[1] < max[1] - 2.0 * h) // top boundary 
				{
					if (x[1] > 5)
						continue;
					if (x[1] > 2 && (x[0] < 10.5 && x[0] > 9.5))
						continue;

				}
			}

			for (int i = 0; i < 2; i++)
			{
				particle.x[i] = x[i];
				particle.x_h[i] = x[i];
			}

			// check if particle is a boundary point
			if (x[0] < 0 || x[1] < 0 || x[0] > max[0] - 2.0 * h || x[1] > max[1] - 2.0 * h || (x[0] < 10.5 && x[0] > 9.5))
				particle.boundary = true;
			else
				particle.boundary = false;

			if (shoaling) {
				if (x[1] < 0.1 * x[0])
					particle.boundary = true;
			}

			particle.calc_index();

			this->initialise_values(&particle);

			particle_list.push_back(particle);
		}
	}


}

//needs to be called each time that all the particles have their positions updated
void SPH_main::allocate_to_grid(void)				
{
	for (int i = 0; i < max_list[0]; i++)
		for (int j = 0; j < max_list[1]; j++)
			search_grid[i][j].clear();

	for (unsigned int cnt = 0; cnt < particle_list.size(); cnt++)
	{
		if ((particle_list[cnt].list_num[0] < -2.0*h || particle_list[cnt].list_num[0] >= max_list[0]) || (particle_list[cnt].list_num[1] < -2.0 * h || particle_list[cnt].list_num[1] >= max_list[1]))
			std::cout << "Particle out of bounds!"<<std::endl;
		search_grid[particle_list[cnt].list_num[0]][particle_list[cnt].list_num[1]].push_back(&particle_list[cnt]);
	}
	return;
}


// Calculate W(r_ij,h)
double SPH_main::compute_W(double dist, double h)
{
	double q = dist / h;
	double W_ij = 10.0 / (7.0 * M_PI * h * h);

	if (h < dist)
		W_ij *= 1.0 - 3.0 * q * q / 2.0 + 3.0 * q * q * q / 4.0;
	else
		W_ij *= (2.0 - q) * (2.0 - q) * (2.0 - q) / 4.0;

	return W_ij;
}

// add elements of ns equation for each particle, during iteraction
void SPH_main::add_to_navier_stokes(SPH_particle *part, SPH_particle *other_part, double dn[2], double dist, bool half)
{
	double dWdr;
	double q = dist / h;
	if (q < 1)
		dWdr = -3.0 * q + 9.0 * (q * q) / 4.0;
	else
		dWdr = -3.0 * (2.0 - q) * (2.0 - q) / 4.0;

	dWdr *= 10.0 / (7.0 * M_PI * h * h * h); // multiply by constant scaling factor

	double e[2];
	e[0] = dn[0] / dist;
	e[1] = dn[1] / dist;

	double v_ij[2];
	double const_1;
	double const_2;
	if (!half)
	{
		v_ij[0] = part->v[0] - other_part->v[0];
		v_ij[1] = part->v[1] - other_part->v[1];

		const_1 = other_part->m * (part->P / (part->rho * part->rho) + other_part->P / (other_part->rho * other_part->rho)) * dWdr;
		const_2 = other_part->m * (1.0 / (part->rho * part->rho) + 1.0 / (other_part->rho * other_part->rho)) * dWdr;
	
	}
	else
	{
		v_ij[0] = part->v_h[0] - other_part->v_h[0];
		v_ij[1] = part->v_h[1] - other_part->v_h[1];

		const_1 = other_part->m_h * (part->P_h / (part->rho_h * part->rho_h) + other_part->P_h / (other_part->rho_h * other_part->rho_h)) * dWdr;
		const_2 = other_part->m_h * (1.0 / (part->rho_h * part->rho_h) + 1.0 / (other_part->rho_h * other_part->rho_h)) * dWdr;
	}	
	
	double mu_term_0 = mu * const_2 * v_ij[0] / dist;
	double mu_term_1 = mu * const_2 * v_ij[1] / dist;

	if (!half)
	{
		part->a[0] += -const_1 * e[0] + mu_term_0;
		part->a[1] += -const_1 * e[1] + mu_term_1;

		part->drho += other_part->m * dWdr * (v_ij[0] * e[0] + v_ij[1] * e[1]);
	}
	else
	{
		part->a_h[0] += -const_1 * e[0] + mu_term_0;
		part->a_h[1] += -const_1 * e[1] + mu_term_1;

		part->drho_h += other_part->m_h * dWdr * (v_ij[0] * e[0] + v_ij[1] * e[1]);
	}

    // Symmetric part
    // const_1 and const_2 for symmetric updates are NOT the same: part should be switched to other_part and viceversa

	double const_1_symmetric, const_2_symmetric;
	if (!half)
	{
    	const_1_symmetric = part->m * (other_part->P / (other_part->rho * other_part->rho) + part->P / (part->rho * part->rho)) * dWdr;
    	const_2_symmetric = part->m * (1.0 / (other_part->rho * other_part->rho) + 1.0 / (part->rho * part->rho)) * dWdr;
	}
	else 
	{
    	const_1_symmetric = part->m_h * (other_part->P_h / (other_part->rho_h * other_part->rho_h) + part->P_h / (part->rho_h * part->rho_h)) * dWdr;
    	const_2_symmetric = part->m_h * (1.0 / (other_part->rho_h * other_part->rho_h) + 1.0 / (part->rho_h * part->rho_h)) * dWdr;
	}
    // mu_term_0 and mu_term_1 for symmetric updates are NOT the same: part should be switched to other_part and viceversa
    // (together with change in v sign)
    double mu_term_0_symmetric = mu * const_2_symmetric * v_ij[0] / dist;
    double mu_term_1_symmetric = mu * const_2_symmetric * v_ij[1] / dist;

	// Symmetric part: update corresponding acceleration of other part, v_ji = -v_ij
	// Symmetric part: update corresponding drho to other_part, v_ji = -v_ij
	if (!half)
	{
		other_part->a[0] += const_1_symmetric * e[0] - mu_term_0_symmetric;
		other_part->a[1] += const_1_symmetric * e[1] - mu_term_1_symmetric;
    	other_part->drho += part->m * dWdr * (v_ij[0] * e[0] + v_ij[1] * e[1]);
	}
	else 
	{
		other_part->a_h[0] += const_1_symmetric * e[0] - mu_term_0_symmetric;
		other_part->a_h[1] += const_1_symmetric * e[1] - mu_term_1_symmetric;
		other_part->drho_h += part->m_h * dWdr * (v_ij[0] * e[0] + v_ij[1] * e[1]);
	}

}


// Euler time-stepping
void SPH_main::Euler_time_stepping(SPH_particle *part)
{
	// ignore position and velocity changes to boundary particles
	if (!part->boundary){
		for (int i = 0; i < 2; i++)
		{
			part->x[i] += this->dt * part->v[i];
			part->v[i] += this->dt * part->a[i];
		}
	}
	part->rho += this->dt * part->drho;
}

// Predictor-Corrector time-stepping
void SPH_main::PC_half_stepping(SPH_particle *part)
{
	// half step
	if (!part->boundary){
		for (int i = 0; i < 2; i++)
		{
			part->x_h[i] = part->x[i] + 0.5 * this->dt * part->v[i];
			part->v_h[i] = part->v[i] + 0.5 * this->dt * part->a[i];
		}
	}
	part->rho_h = part->rho + 0.5 * this->dt * part->drho;
}

// full step
void SPH_main::PC_full_stepping(SPH_particle *part)
{
	if (!part->boundary){
		for (int i = 0; i < 2; i++)
		{
			part->x_h[i] = part->x[i] + 0.5 * this->dt * part->v_h[i];
			part->v_h[i] = part->v[i] + 0.5 * this->dt * part->a_h[i];
			part->x[i] = 2*part->x_h[i] - part->x[i];
			part->v[i] = 2*part->v_h[i] - part->v[i];
		}
	}
	part->rho_h = part->rho + 0.5 * this->dt * part->drho_h;
	part->rho = 2*part->rho_h - part->rho;

}

// deletes particles too close
void SPH_main::remove_close_particles(bool half)
{
	SPH_particle *part, *other_part;
	double dn[2], dist2, tol2 = dx/1000;
	int index[2];

	for (int i = 0; i < max_list[0]; i++)
		for (int j = 0; j < max_list[1]; j++){

			for (unsigned int cnt_1 = 0; cnt_1 < search_grid[i][j].size(); cnt_1++)
			{	
				endloop: ;
				part = search_grid[i][j][cnt_1];
				//part->calc_index();

				// interact with other particles in grid

				for (unsigned int cnt_2 = cnt_1 + 1; cnt_2 < search_grid[i][j].size(); cnt_2++)
				{
					endloop_2: ;
					other_part = search_grid[i][j][cnt_2];

					for (int n = 0; n < 2; n++)
						if (!half)
							dn[n] = part->x[n] - other_part->x[n];
						else 
							dn[n] = part->x_h[n] - other_part->x_h[n];

					dist2 = (dn[0] * dn[0] + dn[1] * dn[1]);

					if ( dist2 < tol2 ){
						// delete particle
						if (part->boundary == false){
							unsigned int z;
							// get other_particle's list number =: z
							for (z = 0; z < particle_list.size(); z++)
								if (part == &particle_list[z])
									break;
							// delete particle
							particle_list.erase(particle_list.begin() + z);
							search_grid[i][j].erase(search_grid[i][j].begin() + cnt_1);
							goto endloop; // restart part iteration from where we left off 
						}

						else if (other_part->boundary == false){
							unsigned int z;
							// get other_particle's list number =: z
							for (z = 0; z < particle_list.size(); z++)
								if (other_part == &particle_list[z])
									break;
							// delete other particle
							particle_list.erase(particle_list.begin() + z);
							search_grid[i][j].erase(search_grid[i][j].begin() + cnt_2);
							goto endloop_2; // restart other_part iteration from where we left off 
						}
					}
				}

				for (int k = 0; k <= 3; k++)
				{
					//using the stencil - this halves the number of operations
					if (k == 0) { index[0] = -1; index[1] = 1; }
					if (k == 1) { index[0] = 0; index[1] = 1; }
					if (k == 2) { index[0] = 1; index[1] = 1; }
					if (k == 3) { index[0] = 1; index[1] = 0; }
				
					int ii = i + index[0];
					int jj = j + index[1];
				
					if (ii >= 0 && ii < max_list[0] && jj >= 0 && jj < max_list[1])
					{
						for (unsigned int cnt_2 = 0; cnt_2 < search_grid[ii][jj].size(); cnt_2++)
						{
							endloop_3: ;
							other_part = search_grid[ii][jj][cnt_2];

							//Calculates the distance between potential neighbours
							for (int n = 0; n < 2; n++)
							{
								if (!half)
									dn[n] = part->x[n] - other_part->x[n];
								else
									dn[n] = part->x_h[n] - other_part->x_h[n];	
							}
							dist2 = dn[0] * dn[0] + dn[1] * dn[1];

							if ( dist2 < tol2 ){
								// delete particle
								if (part->boundary == false){
									unsigned int z;
									// get other_particle's list number =: z
									for (z = 0; z < particle_list.size(); z++)
										if (part == &particle_list[z])
											break;
									// delete particle
									particle_list.erase(particle_list.begin() + z);
									search_grid[i][j].erase(search_grid[i][j].begin() + cnt_1);
									goto endloop; // restart part iteration from where we left off 
 								}

								else if (other_part->boundary == false){
									unsigned int z;
									// get other_particle's list number =: z
									for (z = 0; z < particle_list.size(); z++)
										if (other_part == &particle_list[z])
											break;
									// delete other particle
									particle_list.erase(particle_list.begin() + z);
									search_grid[ii][jj].erase(search_grid[ii][jj].begin() + cnt_2);
									goto endloop_3; // restart other_part iteration from where we left off 
								}
							}

						}
					}

				}


			}

		}


}

//iterates over all particles within 2h of part
void SPH_main::neighbour_iterate(int i_c, int j_c, bool den_smoothing, bool half)
{
    SPH_particle *part, *other_part;
    double dist;                            //distance between particles
    double dn[2];                           //vector from 1st to 2nd particle
    int index[2];                           //to take advantage of symmetry
    //this->allocate_to_grid();
	
    for (unsigned int cnt_1 = 0; cnt_1 < search_grid[i_c][j_c].size(); cnt_1++)
    {
        part = search_grid[i_c][j_c][cnt_1];
        //part->calc_index();

		// interact with other particles in grid
        for (unsigned int cnt_2 = cnt_1 + 1; cnt_2 < search_grid[i_c][j_c].size(); cnt_2++)
        {
            other_part = search_grid[i_c][j_c][cnt_2];

            //other_part->calc_index();
			//Calculates the distance between potential neighbours
			for (int n = 0; n < 2; n++)
			{
				if (!half)
					dn[n] = part->x[n] - other_part->x[n];
				else
					dn[n] = part->x_h[n] - other_part->x_h[n];	
			}
            dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);

            if (dist < 2.0 * h) //only particle within 2h
            {
                if (den_smoothing)
                {
                    double W_ij = compute_W(dist, h);
                    part->sum_W += W_ij;
                    part->sum_W_rho += W_ij / other_part->rho;
                    //Symmetric part
                    other_part->sum_W += W_ij;
                    other_part->sum_W_rho += W_ij / part->rho;
                }
				// add terms to sum into the navier stokes equations for both part and other_part
				if (!half){
					add_to_navier_stokes(part, other_part, dn, dist);
				}
				else{
					add_to_navier_stokes(part, other_part, dn, dist, true);
				}
            }
        }	

		//Searching for neighbours only in the top-right neighbouring cells
		
		for (int k = 0; k <= 3; k++)
		{
			//using the stencil - this halves the number of operations
			if (k == 0) { index[0] = -1; index[1] = 1; }
			if (k == 1) { index[0] = 0; index[1] = 1; }
			if (k == 2) { index[0] = 1; index[1] = 1; }
			if (k == 3) { index[0] = 1; index[1] = 0; }
		
			int i = i_c + index[0];
			int j = j_c + index[1];
		
			if (i >= 0 && i < max_list[0] && j >= 0 && j < max_list[1])
			{
				for (unsigned int cnt_2 = 0; cnt_2 < search_grid[i][j].size(); cnt_2++)
				{
					other_part = search_grid[i][j][cnt_2];
					//other_part->calc_index();

					//Calculates the distance between potential neighbours
					for (int n = 0; n < 2; n++)
					{
						if (!half)
							dn[n] = part->x[n] - other_part->x[n];
						else
							dn[n] = part->x_h[n] - other_part->x_h[n];	
					}

					dist = sqrt(dn[0] * dn[0] + dn[1] * dn[1]);

					if (dist < 2.0 * h)                 //only particle within 2h
					{
						// add terms to sum into the navier stokes equations for both part and other_part
						if (!half){
							add_to_navier_stokes(part, other_part, dn, dist);
						}
						else{
							add_to_navier_stokes(part, other_part, dn, dist, true);
						}

						if (den_smoothing)
						{
							double W_ij = compute_W(dist, h);
							part->sum_W += W_ij;
							part->sum_W_rho += W_ij / other_part->rho;
							//Symmetric part
							other_part->sum_W += W_ij;
							other_part->sum_W_rho += W_ij / part->rho;
						}
					}
					
				}

			}
		}

		if (den_smoothing)
		{
			double W_ii = compute_W(0, h);
			part->sum_W += W_ii;
			part->sum_W_rho += W_ii / part->rho;
			part->rho = part->sum_W / part->sum_W_rho;
		}	

	}

}
