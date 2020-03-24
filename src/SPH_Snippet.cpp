#include "../includes/SPH_2D.h"
#include "../includes/file_writer.h"
// #include "SPH_2D.cpp"
// #include "file_writer.cpp"
#include <omp.h>

SPH_main domain;

int main(void)
{

	bool shoaling = false;

	domain.set_values();										//Set simulation parameters
	domain.initialise_grid();									//initialise simulation grid

	// Places initial points - will need to be modified to include boundary points 
	// and the specifics of where the fluid is in the domain
	domain.place_points(domain.min_x, domain.max_x, shoaling);
	//domain.place_points_conf1(domain.min_x, domain.max_x, shoaling);
	//domain.place_points_conf2(domain.min_x, domain.max_x, shoaling);
	//domain.place_points_conf3(domain.min_x, domain.max_x, shoaling);
	//domain.place_points_conf4(domain.min_x, domain.max_x, shoaling);
	//domain.place_points_conf5(domain.min_x, domain.max_x, shoaling);

	//needs to be called for each time step
	domain.allocate_to_grid();
	
	write_file(0, &domain.particle_list);
	
	double maxt = 30;
	int cnt = 1;

	for (double t = 0; t <= maxt; t += domain.dt)
	{
		// half step
		domain.allocate_to_grid();

		for (unsigned int i = 0; i < domain.particle_list.size(); i++){
			domain.update_mass_pressure(&domain.particle_list[i], shoaling, true);
		}

        for (int i = 0; i < domain.max_list[0]; i++)
            for (int j = 0; j < domain.max_list[1]; j++)
                domain.neighbour_iterate(i, j, false, true);


		for (unsigned int i = 0; i < domain.particle_list.size(); i++)
		{
			domain.PC_half_stepping(&domain.particle_list[i]);
		}

		// remove close particles
		domain.remove_close_particles(true);

		// full step
		domain.allocate_to_grid();

		for (unsigned int i = 0; i < domain.particle_list.size(); i++){
			domain.update_mass_pressure(&domain.particle_list[i], shoaling, false);
		}

        for (int i = 0; i < domain.max_list[0]; i++)
            for (int j = 0; j < domain.max_list[1]; j++)
            {
                //std::cout << i << " " << j << std::endl;
                if (cnt % 10 == 0) // pass in smoothing = true every 20 steps
                    domain.neighbour_iterate(i, j, true);
                else
                    domain.neighbour_iterate(i, j, false);
            }


		for (unsigned int i = 0; i < domain.particle_list.size(); i++)
		{
			domain.PC_full_stepping(&domain.particle_list[i]);
		}

		if (cnt % 32 == 0){ //output to file
			cerr << "time: " << t << endl;
			write_file(cnt, &domain.particle_list);
		}

		// remove close particles
		domain.remove_close_particles();

		cnt++;
	}

	write_file(cnt, &domain.particle_list);
	
	return 0;
}
