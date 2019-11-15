/*
###############################################################################
# If you use PhysiCell in your project, please cite PhysiCell and the version #
# number, such as below:                                                      #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1].    #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# See VERSION.txt or call get_PhysiCell_version() to get the current version  #
#     x.y.z. Call display_citations() to get detailed information on all cite-#
#     able software used in your PhysiCell application.                       #
#                                                                             #
# Because PhysiCell extensively uses BioFVM, we suggest you also cite BioFVM  #
#     as below:                                                               #
#                                                                             #
# We implemented and solved the model using PhysiCell (Version x.y.z) [1],    #
# with BioFVM [2] to solve the transport equations.                           #
#                                                                             #
# [1] A Ghaffarizadeh, R Heiland, SH Friedman, SM Mumenthaler, and P Macklin, #
#     PhysiCell: an Open Source Physics-Based Cell Simulator for Multicellu-  #
#     lar Systems, PLoS Comput. Biol. 14(2): e1005991, 2018                   #
#     DOI: 10.1371/journal.pcbi.1005991                                       #
#                                                                             #
# [2] A Ghaffarizadeh, SH Friedman, and P Macklin, BioFVM: an efficient para- #
#     llelized diffusive transport solver for 3-D biological simulations,     #
#     Bioinformatics 32(8): 1256-8, 2016. DOI: 10.1093/bioinformatics/btv730  #
#                                                                             #
###############################################################################
#                                                                             #
# BSD 3-Clause License (see https://opensource.org/licenses/BSD-3-Clause)     #
#                                                                             #
# Copyright (c) 2015-2018, Paul Macklin and the PhysiCell Project             #
# All rights reserved.                                                        #
#                                                                             #
# Redistribution and use in source and binary forms, with or without          #
# modification, are permitted provided that the following conditions are met: #
#                                                                             #
# 1. Redistributions of source code must retain the above copyright notice,   #
# this list of conditions and the following disclaimer.                       #
#                                                                             #
# 2. Redistributions in binary form must reproduce the above copyright        #
# notice, this list of conditions and the following disclaimer in the         #
# documentation and/or other materials provided with the distribution.        #
#                                                                             #
# 3. Neither the name of the copyright holder nor the names of its            #
# contributors may be used to endorse or promote products derived from this   #
# software without specific prior written permission.                         #
#                                                                             #
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" #
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE   #
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE  #
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE   #
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR         #
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF        #
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS    #
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN     #
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)     #
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE  #
# POSSIBILITY OF SUCH DAMAGE.                                                 #
#                                                                             #
###############################################################################
*/

#include "./Mitra_biotech_sim.h"
#include "../modules/PhysiCell_settings.h"

#include <cmath>
#include <iostream>
#include <random>

Cell_Definition immune_cell; 
Cell_Definition cancer_cell; 

void create_immune_cells( void )
{
	immune_cell = cell_defaults;
	
	immune_cell.name = "immune cell";
	immune_cell.type = 1;
	
	// turn off proliferation 
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 

	immune_cell.phenotype.secretion.uptake_rates[1] = 0.1;
	
	immune_cell.phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 0.0; 	
	
	return;
}

void create_cancer_cells( void )
{
	cancer_cell = cell_defaults;
	
	cancer_cell.phenotype.secretion.uptake_rates[1] = 0.1;
	
	cancer_cell.name = "cancer cell";
	cancer_cell.type = 2; 
	
	return;
}

void create_cell_types( void )
{
	// housekeeping 
	SeedRandom( parameters.ints( "random_seed" ) ); 
	initialize_default_cell_definition();	
	
	//setting cycle model to live
	cell_defaults.phenotype.cycle.sync_to_cycle_model( live ); 
	
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	
	double C0 = parameters.doubles("C0");
	cell_defaults.phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0.0000064812;
	
	// Make sure we're ready for 2D
	cell_defaults.functions.set_orientation = up_orientation;  
	cell_defaults.phenotype.geometry.polarity = 1.0; 
	cell_defaults.phenotype.motility.restrict_to_2D = true; 
	
	cell_defaults.phenotype.secretion.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.molecular.sync_to_microenvironment( &microenvironment );
	cell_defaults.phenotype.sync_to_functions( cell_defaults.functions ); 


 	// turn off death
	int apoptosis_index = cell_defaults.phenotype.death.find_death_model_index( PhysiCell_constants::apoptosis_death_model ); 
	cell_defaults.phenotype.death.rates[apoptosis_index] = 0.0; 
	
	// reduce cell velocity
	cell_defaults.phenotype.motility.migration_speed = 0.05;
	
	// add variables to track virus infection start time, length of time and amount of virus
	cell_defaults.custom_data.add_variable( "intracellular_virus_amount", "dimensionless", 0.0 );
	cell_defaults.custom_data.add_variable( "bias_phenotype", "dimensionless", 0.0 );
	cell_defaults.custom_data.add_variable( "replication_start_time", "dimensionless", 0.0 );
	cell_defaults.custom_data.add_variable( "infection_time_start", "dimensionless", -1.0 );
	cell_defaults.custom_data.add_variable( "infection_time_length", "dimensionless", -1.0 );
	cell_defaults.custom_data.add_variable( "time_of_injection", "dimensionless", 0.0 );
	
	// turn off secretion from these cells (oxygen and virus)
	cell_defaults.phenotype.secretion.secretion_rates[0] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[0] = 0;
	cell_defaults.phenotype.secretion.saturation_densities[0] = 10; 
	
	static int virus_index = microenvironment.find_density_index( "virus");
	
	cell_defaults.phenotype.secretion.secretion_rates[virus_index] = 0;
	cell_defaults.phenotype.secretion.uptake_rates[virus_index] = 0.1;
	cell_defaults.phenotype.secretion.saturation_densities[virus_index] = 10; 
	
	// update cell and phenotype based on virus dynamics only
	//cell_defaults.functions.update_phenotype = virus_dynamics_TRAIL_MODEL3; 
	
	//cell_defaults.functions.update_phenotype = virus_dynamics_TRAIL_MODEL3;//ECM_movement;//virus_dynamics; 
	cell_defaults.functions.update_phenotype = infection_dynamics;
	cell_defaults.phenotype.motility.is_motile = true; 
	
	//create the immune and cancer cell types
	create_immune_cells();
	create_cancer_cells();
		
	return; 
}

void setup_microenvironment( void )
{
	// set domain parameters
	
	// make sure not override and go back to 2D 
	if( default_microenvironment_options.simulate_2D == false )
	{
		std::cout << "Warning: overriding XML config option and setting to 2D!" << std::endl; 
		default_microenvironment_options.simulate_2D = true; 
	}
	initialize_microenvironment(); 	
	
	
	static int virus_index = microenvironment.find_density_index( "virus");

	for( int n = 0 ; n < microenvironment.mesh.voxels.size(); n++ )
	{	
		std::vector<double> ECMdense = microenvironment.mesh.voxels[n].center; 
		//assign random ECM density to microenvironment voxel
		if( ECMdense[0]<-970) //white
		{microenvironment(n)[virus_index] = 10;}
		else if( ECMdense[0]>970)
		{microenvironment(n)[virus_index] = 10;}
		else if(ECMdense[1]<-970)
		{microenvironment(n)[virus_index] = 10;}
		else if(ECMdense[1]>970)
		{microenvironment(n)[virus_index] = 10;}		
		else
		{microenvironment(n)[virus_index] = 0;}
	}
		
	return; 
}	
void setup_tissue_square( void )
{
	double C0 = 2*1e4;//parameters.doubles("C0");//5127;
		
	static int virus_index = microenvironment.find_density_index( "virus");
	double MOI = parameters.doubles("MOI");
	
	double start_of_boundary = -979;//6913/2;
	double end_of_boundary = 979;//6913/2;
	double interval_lengths = (end_of_boundary-start_of_boundary)/(sqrt(C0)+1);
	double start_of_xboundary = start_of_boundary+interval_lengths;
		
	double length_x = microenvironment.mesh.bounding_box[3] - 
		microenvironment.mesh.bounding_box[0]; 
		
	double length_y = microenvironment.mesh.bounding_box[4] - 
		microenvironment.mesh.bounding_box[1]; 
		
		
	Cell* pCell = NULL; 
	
	double x = 0.0;
	double y = 0.0;
	
	double xgrid = 0.0;
	double ygrid = 0.0;
		
	
	//placing vien cells and normal cells in a hexagonal grid
	for( int i=0; i<sqrt(C0); i++ )
	{	
		if( i % 2==0 )
		{x = start_of_xboundary+i*interval_lengths;}//18;
		else
		{x = start_of_xboundary+1/2*interval_lengths+i*interval_lengths;}
		
		for( int j=0; j<sqrt(C0); j++ )
		{
			y = start_of_boundary+j*interval_lengths;//18;
			if(UniformRandom() < 0.3679)
			{pCell = create_cell( cancer_cell );}
		    else 
			{pCell = create_cell( immune_cell );}
			double x = microenvironment.mesh.bounding_box[0] + UniformRandom() * length_x; 
			double y = microenvironment.mesh.bounding_box[1] + UniformRandom() * length_y; 
		
			pCell->assign_position( x , y , 0.0 );
			pCell->phenotype.molecular.internalized_total_substrates[virus_index] = 0; 
			double bias_phenotype_index = pCell->custom_data.find_variable_index( "bias_phenotype" );
		    pCell->custom_data.variables[bias_phenotype_index].value =  UniformRandom();
		}
		
	}

	return;
	
}

void infection_dynamics( Cell* pCell, Phenotype& phenotype, double dt )
{
	
	static int v_index = microenvironment.find_density_index( "virus");
	pCell->phenotype.molecular.internalized_total_substrates[v_index];
	static int apoptosis_model_index =	pCell->phenotype.death.find_death_model_index( "apoptosis" );
	
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	 
	double internal_virus = pCell->phenotype.molecular.internalized_total_substrates[v_index];
	double replication_rate = 10;
	
	// find intracellular virus index
	static double bias_phenotype_index = pCell->custom_data.find_variable_index( "bias_phenotype" );
	double bias_phenotype = pCell->custom_data.variables[bias_phenotype_index].value;
	
   if( bias_phenotype<0.8 && internal_virus>1)//infected so not replicating
   {
	    pCell->phenotype.cycle.data.transition_rate(cycle_start_index,cycle_end_index) = 0.0;
	   //std::cout<<bias_phenotype<<std::endl;
		pCell->phenotype.molecular.internalized_total_substrates[v_index] = internal_virus+replication_rate*6;//virus replicates at linear rate
		if(pCell->phenotype.molecular.internalized_total_substrates[v_index]>1e2) 
		{pCell->phenotype.secretion.uptake_rates[v_index] = 0;}
   
		if(pCell->phenotype.molecular.internalized_total_substrates[v_index]>1e4)
		{pCell->start_death( apoptosis_model_index );}	   
	
		pCell->custom_data.variables[bias_phenotype_index].value = 0;
   }
   
	return;
}


void virus_dynamics_TRAIL_MODEL3( Cell* pCell, Phenotype& phenotype, double dt )
{
	// PROLIFERATION ^^
	
	int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
    //static int pressure_store_index =pCell->phenotype.death.find_death_model_index( "pressure_store" );	
	//pCell->custom_data.variables[pressure_store_index].value = pCell->state.simple_pressure;
		double U = all_cells->size();
	pCell->phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0.00064812*(1-U/5127);	
	
	
	
	// PROLIFERATION ^^
	
	
	// find microenvironment virus index
	//static int virus_signal_index	= microenvironment.find_density_index( "virus" ); 
	
	// find microenvironment TRAIL index
	static int TRAIL_signal_index	= microenvironment.find_density_index( "TRAIL" ); 
	
	// sample amount of virus in microenvironment
	//double virus_amount = pCell->nearest_density_vector()[virus_signal_index];
	// sample amount of TRAIL in microenvironment
	double TRAIL_amount = pCell->nearest_density_vector()[TRAIL_signal_index];
	
	// find intracellular virus index
	static double intracellular_virus_index = pCell->custom_data.find_variable_index( "intracellular_virus_amount" );
	
	// find intracellular virus index
	static double intracellular_TRAIL_index = pCell->custom_data.find_variable_index( "intracellular_TRAIL_amount" );
	
	// find replication start time index
	static double replication_start_index = pCell->custom_data.find_variable_index( "replication_start_time" );
	
	// find death model index 
	static int apoptosis_model_index =	pCell->phenotype.death.find_death_model_index( "apoptosis" );
	
	// parameters 
	static double infection_density_capacity = parameters.doubles("infection_density_capacity");
	static double intrinsic_infection_rate = parameters.doubles("intrinsic_infection_rate"); //2; 
	static double viral_replication_rate = parameters.doubles("viral_replication_rate");	
	static double virus_replication_minimum = parameters.doubles("virus_replication_minimum");	
	static double TRAIL_secretion_rate = parameters.doubles("TRAIL_secretion_rate");
	static double TRAIL_generation_rate = parameters.doubles("TRAIL_generation_rate");
	static double TRAIL_killing_level = parameters.doubles("TRAIL_killing_level");
	static double M = parameters.doubles("M");
	
	double c_I = intrinsic_infection_rate;//rate of infection  0.1
	double K = infection_density_capacity;//capacity of cell
	double c_R = 0;//viral_replication_rate; // 0.1
	double c_T = TRAIL_generation_rate; // 0.1
	double n_Istar = virus_replication_minimum;//10
	//double rho_E = virus_amount;
	double s_T =  TRAIL_secretion_rate;
	int alpha = 1;
	double L = 3000;//K/2+1.3*K/2;
	double s_tau = parameters.doubles("TRAIL_secretion_time_commence");
	
	double V_voxel = microenvironment.mesh.voxels[1].volume;//volume of voxel
	
	double n_I = pCell->custom_data.variables[intracellular_virus_index].value;
	//double n_E = virus_amount*V_voxel;
	double T_I = pCell->custom_data.variables[intracellular_TRAIL_index].value;
	double replication_start_time = pCell->custom_data.variables[replication_start_index].value;
	double T_E = pCell->nearest_density_vector()[TRAIL_signal_index];
	double pstarT = phenotype.secretion.saturation_densities[TRAIL_signal_index]; 

	double n_T = infection_density_capacity;

	double n_I_zero = 0;
	if( n_I < n_T )
	{n_I_zero = n_I;}
	else
	{n_I_zero = n_T;}
	
	if( pCell->phenotype.death.dead == false )// cell not dead	
	{	
		if( replication_start_time = 0)
		{
			pCell->custom_data.variables[replication_start_index].value = PhysiCell_globals.current_time;
		}
		
		if( PhysiCell_globals.current_time>s_tau+replication_start_time )
		{
			pCell->custom_data.variables[intracellular_TRAIL_index].value = T_I+dt*(c_T*n_I_zero-s_T*T_I*(pstarT-T_E));
			pCell->phenotype.secretion.secretion_rates[TRAIL_signal_index] = s_T*T_I/V_voxel;
			//std::cout<<"here 2 "<<pCell->custom_data.variables[intracellular_TRAIL_index].value<<" "<<s_T*T_I/V_voxel<<" "<<c_T*n_I_zero-s_T*T_I*(pstarT-T_E)<<std::endl;
		}
		else
		{
			pCell->custom_data.variables[intracellular_TRAIL_index].value = T_I+dt*(c_T*n_I_zero);
			pCell->phenotype.secretion.secretion_rates[TRAIL_signal_index] = 0;
		}
		
		if(pCell->custom_data.variables[intracellular_TRAIL_index].value<0)
		{	
			std::cout<<"NEGATIVE HELP TRAIL7 "<<pCell->custom_data.variables[intracellular_TRAIL_index].value<<" "<<T_I<<" "<<c_T<<std::endl;	
			system("pause");
		}
	}
	else if( pCell->phenotype.death.dead==true && pCell->custom_data.variables[intracellular_TRAIL_index].value>0)
	{
		//std::cout<<"Cell died"<<std::endl;	
		pCell->phenotype.secretion.secretion_rates[TRAIL_signal_index] = 0;
		double TRAIL_amount_in_cell = pCell->custom_data.variables[intracellular_TRAIL_index].value;
		if(TRAIL_amount_in_cell/V_voxel+pCell->nearest_density_vector()[TRAIL_signal_index]>pstarT)
		{
			//std::cout<<TRAIL_amount_in_cell/V_voxel<<" "<<pCell->nearest_density_vector()[TRAIL_signal_index]<<std::endl;
			double amount_addble = pstarT-pCell->nearest_density_vector()[TRAIL_signal_index];
			pCell->nearest_density_vector()[TRAIL_signal_index] += amount_addble;
			pCell->custom_data.variables[intracellular_TRAIL_index].value = TRAIL_amount_in_cell-amount_addble*V_voxel;
		}
		else
		{pCell->nearest_density_vector()[TRAIL_signal_index] += TRAIL_amount_in_cell/V_voxel;}
	}
	
	/*	if( pCell->phenotype.death.dead==false && virus_amount>1e-2)
		{// cell becomes infected
	
			//std::cout<<n_I<<" "<<n_Istar<<std::endl;
			if( n_I > n_Istar && n_I< K )
			{
				
				//std::cout<<"NEGATIVE HELP0"<<" "<<n_I<<" "<<K<<" "<<n_Istar<<std::endl;
				
				
				if( replication_start_time = 0)
				{
					pCell->custom_data.variables[replication_start_index].value = PhysiCell_globals.current_time;
				}
				if( PhysiCell_globals.current_time>s_tau+replication_start_time )
				{
					pCell->custom_data.variables[intracellular_TRAIL_index].value = T_I*n_I_zero+dt*(c_T-s_T*T_I*(pstarT-T_E));
					pCell->phenotype.secretion.secretion_rates[TRAIL_signal_index] = s_T*T_I/V_voxel;
					
				if(pCell->custom_data.variables[intracellular_TRAIL_index].value<0)
				{
					std::cout<<"NEGATIVE HELP TRAIL"<<std::endl;
					system("pause");
				}
				}
				else
				{
					pCell->custom_data.variables[intracellular_TRAIL_index].value = T_I*n_I_zero+dt*(c_T);
					pCell->phenotype.secretion.secretion_rates[TRAIL_signal_index] = 0;
					
					if(pCell->custom_data.variables[intracellular_TRAIL_index].value<0)
					{
						std::cout<<"NEGATIVE HELP TRAIL7"<<T_I<<c_T<<std::endl;
						system("pause");
					}
				}
				n_I = pCell->custom_data.variables[intracellular_virus_index].value;
				pCell->phenotype.secretion.uptake_rates[virus_signal_index] = c_I*(1-n_I/K);
				
				if( pCell->custom_data.variables[intracellular_virus_index].value<0 )
				{
					
					std::cout<<"NEGATIVE HELP1"<<n_I+dt*(c_I*n_E*(1-n_I/K)+c_R)<<" "<<n_I<<" "<<K<<" "<<n_Istar<<" "<<n_E<<" "<<dt<<" "<<c_I<<" "<<c_R<<std::endl;
					system("pause");
				}
			}
			else if( n_I > K)
			{
				if( replication_start_time = 0)
				{
					pCell->custom_data.variables[replication_start_index].value = PhysiCell_globals.current_time;
				}
				if( PhysiCell_globals.current_time>s_tau+replication_start_time )
				{
					pCell->custom_data.variables[intracellular_TRAIL_index].value = T_I*n_I_zero+dt*(c_T-s_T*T_I*(pstarT-T_E));
					pCell->phenotype.secretion.secretion_rates[TRAIL_signal_index] = s_T*T_I/V_voxel;
					
					if(pCell->custom_data.variables[intracellular_TRAIL_index].value<0)
					{
						std::cout<<"NEGATIVE HELP TRAIL6"<<std::endl;
						system("pause");
					}
								}
								else
								{
									pCell->custom_data.variables[intracellular_TRAIL_index].value = T_I*n_I_zero+dt*(c_T);//-s_T*T_I);
									pCell->phenotype.secretion.secretion_rates[TRAIL_signal_index] = 0;//s_T*T_I/V_voxel;
									
					if(pCell->custom_data.variables[intracellular_TRAIL_index].value<0)
					{
						std::cout<<"NEGATIVE HELP TRAIL5"<<std::endl;
						system("pause");
					}
				}
				pCell->phenotype.secretion.uptake_rates[virus_signal_index] = 0;
				if( pCell->custom_data.variables[intracellular_virus_index].value<0 )
				{
					
					std::cout<<"NEGATIVE HELP2"<<" "<<pCell->custom_data.variables[intracellular_virus_index].value<<" "<<n_I+dt*(c_R)<<" "<<T_I+dt*(c_T-s_T*T_I)<<" "<<n_I<<" "<<K<<" "<<n_Istar<<std::endl;
					system("pause");
				}
			}
			else
			{
				if( pCell->custom_data.variables[intracellular_virus_index].value<0 )
				{
					
					std::cout<<"NEGATIVE HELP3"<<" "<<n_I<<" "<<K<<" "<<n_Istar<<std::endl;
					system("pause");
				}
				
			}
			}
		else if( n_I>0 && pCell->phenotype.death.dead==false )
		{// cell becomes infected
	
			
	
			if( n_I > n_Istar && n_I< K )
			{
				
				//std::cout<<"NEGATIVE HELP0"<<" "<<n_I<<" "<<K<<" "<<n_Istar<<std::endl;
				
				
				if( replication_start_time = 0)
				{
					pCell->custom_data.variables[replication_start_index].value = PhysiCell_globals.current_time;
				}
				if( PhysiCell_globals.current_time>s_tau+replication_start_time )
				{
					pCell->custom_data.variables[intracellular_TRAIL_index].value = T_I*n_I_zero+dt*(c_T-s_T*T_I*(pstarT-T_E));
					pCell->phenotype.secretion.secretion_rates[TRAIL_signal_index] = s_T*T_I/V_voxel;
					
					if(pCell->custom_data.variables[intracellular_TRAIL_index].value<0)
					{
						std::cout<<"NEGATIVE HELP TRAIL4"<<std::endl;
						system("pause");
					}
				}
				else
				{
					pCell->custom_data.variables[intracellular_TRAIL_index].value = T_I*n_I_zero+dt*(c_T);
					pCell->phenotype.secretion.secretion_rates[TRAIL_signal_index] = 0;
					
					if(pCell->custom_data.variables[intracellular_TRAIL_index].value<0)
					{
						std::cout<<"NEGATIVE HELP TRAIL3"<<T_I<<" "<<c_T<<" "<<pCell->custom_data.variables[intracellular_TRAIL_index].value<<std::endl;
						system("pause");
					}
				}
				n_I = pCell->custom_data.variables[intracellular_virus_index].value;
				pCell->phenotype.secretion.uptake_rates[virus_signal_index] = c_I*(1-n_I/K);
				
				if( pCell->custom_data.variables[intracellular_virus_index].value<0 )
				{
					
					std::cout<<"NEGATIVE HELP1"<<n_I+dt*(c_I*n_E*(1-n_I/K)+c_R)<<" "<<n_I<<" "<<K<<" "<<n_Istar<<" "<<n_E<<" "<<dt<<" "<<c_I<<" "<<c_R<<std::endl;
					system("pause");
				}
			}
			else if( n_I > K)
			{
				if( replication_start_time = 0)
				{
					pCell->custom_data.variables[replication_start_index].value = PhysiCell_globals.current_time;
				}
				if( PhysiCell_globals.current_time>s_tau+replication_start_time )
				{
					pCell->custom_data.variables[intracellular_TRAIL_index].value = T_I*n_I_zero+dt*(c_T-s_T*T_I*(pstarT-T_E));
					pCell->phenotype.secretion.secretion_rates[TRAIL_signal_index] = s_T*T_I/V_voxel;
					
					if(pCell->custom_data.variables[intracellular_TRAIL_index].value<0)
					{
						std::cout<<"NEGATIVE HELP TRAIL2"<<std::endl;
						system("pause");
					}
				}
				else
				{
					pCell->custom_data.variables[intracellular_TRAIL_index].value = T_I*n_I_zero+dt*(c_T);//-s_T*T_I);
					pCell->phenotype.secretion.secretion_rates[TRAIL_signal_index] = 0;//s_T*T_I/V_voxel;
					
					if(pCell->custom_data.variables[intracellular_TRAIL_index].value<0)
					{
						std::cout<<"NEGATIVE HELP TRAIL1"<<std::endl;
						system("pause");
					}
				}
				pCell->phenotype.secretion.uptake_rates[virus_signal_index] = 0;
				if( pCell->custom_data.variables[intracellular_virus_index].value<0 )
				{
					
					std::cout<<"NEGATIVE HELP2"<<" "<<pCell->custom_data.variables[intracellular_virus_index].value<<" "<<n_I+dt*(c_R)<<" "<<T_I+dt*(c_T-s_T*T_I)<<" "<<n_I<<" "<<K<<" "<<n_Istar<<std::endl;
					system("pause");
				}
			}
			else
			{
				n_I = pCell->custom_data.variables[intracellular_virus_index].value;
				pCell->phenotype.secretion.uptake_rates[virus_signal_index] = c_I*(1-n_I/K);
				if( pCell->custom_data.variables[intracellular_virus_index].value<0 )
				{
					
					std::cout<<"NEGATIVE HELP3"<<" "<<n_I<<" "<<K<<" "<<n_Istar<<std::endl;
					system("pause");
				}
				
			}
			
			/*if( n_I>L/5 && n_I<L )
			{
				pCell->phenotype.death.rates[apoptosis_model_index] = n_I*n_I*n_I/(L/2*L/2*L/2+n_I*n_I*n_I);
				
			}
			else if(n_I > L)
			{
				//pCell->phenotype.death.rates[apoptosis_model_index] = 1000;
				pCell->start_death( apoptosis_model_index );				
				
			}
			else
			{
				pCell->phenotype.death.rates[apoptosis_model_index]=0;
			}*/
		//}
		/*else if( T_E>1e-5 && pCell->phenotype.death.dead==false)
		{
			double M = 0.0005;
			pCell->phenotype.death.rates[apoptosis_model_index] = T_E*T_E/(M*M+T_E*T_E);
			//if (pCell->phenotype.death.rates[apoptosis_model_index]>0.01)
			//std::cout<<"apop rate"<<pCell->phenotype.death.rates[apoptosis_model_index] <<std::endl;
			
		}*/
		/*else if( pCell->phenotype.death.dead==true)
		{
			//std::cout<<"Cell died"<<std::endl;
			pCell->phenotype.secretion.uptake_rates[virus_signal_index] = 0;
			pCell->phenotype.secretion.secretion_rates[virus_signal_index] = 0;		
			pCell->phenotype.secretion.secretion_rates[TRAIL_signal_index] = 0;
			if( pCell->custom_data.variables[intracellular_virus_index].value>0)
				{
					double virus_amount_in_cell = pCell->custom_data.variables[intracellular_virus_index].value;
					//std::cout<<"Cell dies: " << virus_amount_in_cell<<std::endl;
					pCell->phenotype.secretion.saturation_densities[virus_signal_index] = 5;
					pstarV = 5;
					double TRAIL_amount_in_cell = pCell->custom_data.variables[intracellular_TRAIL_index].value;
					if(TRAIL_amount_in_cell/V_voxel+pCell->nearest_density_vector()[TRAIL_signal_index]>pstarT)
					{
						std::cout<<TRAIL_amount_in_cell/V_voxel<<" "<<pCell->nearest_density_vector()[TRAIL_signal_index]<<std::endl;
						double amount_addble = pstarT-pCell->nearest_density_vector()[TRAIL_signal_index];
						pCell->nearest_density_vector()[TRAIL_signal_index] += amount_addble;
						pCell->custom_data.variables[intracellular_TRAIL_index].value = TRAIL_amount_in_cell-amount_addble*V_voxel;
					}
					else
					{pCell->nearest_density_vector()[TRAIL_signal_index] += TRAIL_amount_in_cell/V_voxel;}
					
					if( virus_amount_in_cell/V_voxel+pCell->nearest_density_vector()[virus_signal_index]>pstarV )
					{
						//std::cout<<virus_amount_in_cell/V_voxel<<" "<<pCell->nearest_density_vector()[virus_signal_index]<<std::endl;
						double amount_addble = pstarV-pCell->nearest_density_vector()[virus_signal_index];
						pCell->nearest_density_vector()[virus_signal_index] += 0;//amount_addble;
						pCell->custom_data.variables[intracellular_virus_index].value = virus_amount_in_cell-amount_addble*V_voxel;
					
					}
					else
					{pCell->nearest_density_vector()[virus_signal_index] += 0;//virus_amount_in_cell/V_voxel;
					pCell->custom_data.variables[intracellular_virus_index].value=0;
					//pCell->phenotype.secretion.saturation_densities[virus_signal_index] = 5;
					pCell->custom_data.variables[intracellular_virus_index].value=0;
					}
				}
		}*/
	
	/*int cycle_start_index = live.find_phase_index( PhysiCell_constants::live ); 
	int cycle_end_index = live.find_phase_index( PhysiCell_constants::live ); 
	
	if (pCell->state.simple_pressure > 10)//30
	{
		pCell->phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0.0;
	}
	else 
	{
		pCell->phenotype.cycle.data.transition_rate( cycle_start_index , cycle_end_index ) = 0.0001428;
	}
	//std::cout<<pCell->state.simple_pressure <<std::endl;
	
	*/
	if(pCell->custom_data.variables[intracellular_TRAIL_index].value<0)
	{
		std::cout<<"NEGATIVE HELP TRAIL"<<std::endl;
	}
	
		
	 if( T_E>0.01 && pCell->phenotype.death.dead==false)
	{
		double M = 0.8;
		pCell->phenotype.death.rates[apoptosis_model_index] = (T_E*T_E*T_E/(M*M*M+T_E*T_E*T_E))*0.05;
	}
		
		
	
	return;
}

std::vector<std::string> colouring_by_intracellular_virus_amount( Cell* pCell )
{

	std::vector< std::string > output( 4, "darkgrey" ); 
	
	static int v_index = microenvironment.find_density_index( "virus");
	static int infection_density_capacity = parameters.doubles("infection_density_capacity");
	
	double p_min = 1;
	double p_max = 1e4;

	// find intracellular virus index
	static double bias_phenotype_index = pCell->custom_data.find_variable_index( "bias_phenotype" );
	double bias_phenotype = pCell->custom_data.variables[bias_phenotype_index].value;
	
	static double  V_0 = parameters.doubles("no_of_viruses_in_initial_injection")/221; //CIRCLE - no of vein cells
	
	
	double n_I = pCell->phenotype.molecular.internalized_total_substrates[v_index];//pCell->custom_data.find_variable_index( "intracellular_virus_amount" );

		if(pCell->type==1 && pCell->phenotype.death.dead==false)
		{
			int oncoprotein = (int) round( (1.0/(p_max-p_min)) * (pCell->phenotype.molecular.internalized_total_substrates[v_index]-p_min) * 255.0 ); 
			char szTempString [128]; // ceates a character array that can store 128
			sprintf( szTempString , "rgb(%u,%u,%u)", oncoprotein, oncoprotein, 255-oncoprotein ); // puts oncoprotein, oncoprotein and 255-oncoprotein in place of u u u
			
			
			if( bias_phenotype<0.9 && n_I>1)
			{
				
				int oncoprotein = (int) round( (1.0/(p_min-p_max)) * (pCell->phenotype.molecular.internalized_total_substrates[v_index]-p_max) * 230.0 ); 
				char szTempString [128]; // ceates a character array that can store 128
				sprintf( szTempString , "rgb(%u,%u,%u)", oncoprotein, oncoprotein,255); // puts oncoprotein, oncoprotein and 255-oncoprotein in place of u u u
				output[0].assign( szTempString );
				output[1]="blue";
				output[2].assign( szTempString );
				output[3]="blue";
							
				return output;
			}
			else if( pCell->phenotype.molecular.internalized_total_substrates[v_index] <0)
			{
				std::cout<<"NEGATIVE"<<std::endl;
				std::vector< std::string > output( 4, "darkgrey" ); 
				return output;
			}
			else
			{
					output[0] = "darkorange";
					output[1] = "darkorange";
					output[2] = "chocolate";
					output[3] = "chocolate)";
					return output; 
			}
		}
		
		if( pCell->type == 2 && pCell->phenotype.death.dead==false)
		{
			if(  bias_phenotype<0.8 && n_I>1)//n_I > 1 )
			{
				double p_min = 1;
				double p_max = 1e4;//V_0;
				
				int oncoprotein = (int) round( (1.0/(p_min-p_max)) * (pCell->phenotype.molecular.internalized_total_substrates[v_index]-p_max) * 230.0 ); 
				char szTempString [128]; // ceates a character array that can store 128
				sprintf( szTempString , "rgb(%u,%u,%u)", 255, oncoprotein, oncoprotein); // puts oncoprotein, oncoprotein and 255-oncoprotein in place of u u u
				output[0].assign( szTempString );
				output[1]="red";
				output[2].assign( szTempString );
				output[3]="red";
				
				return output;
			}
			else
			{
					output[0] = "orchid";//"rgb(255,230,230)";
					output[1] = "orchid";
					output[2] = "plum";//"rgb(255,230,230)";
					output[3] = "plum";

					return output; 
					
					return output; 
			}
		}
	if( pCell->phenotype.death.dead == true )
	{ 
			output[0] = "rgb(255, 255, 224)";
			output[1] = "rgb(255, 255, 224)";
			output[2] = "rgb(255, 228, 181)";
			output[3] = "rgb(255, 228, 181)";
			
		
		return output; 
	}
	 
	return output;
}

