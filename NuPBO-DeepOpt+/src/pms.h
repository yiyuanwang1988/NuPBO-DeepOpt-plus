#ifndef _PMS_H_
#define _PMS_H_

#include <cmath>
#include "basis_pms.h"
#include "deci.h"
#include <sstream>
#include <algorithm>
#include <iomanip>

Satlike::Satlike()
{
	problem_weighted = 1;
	partial = 1; // 1 if the instance has hard clauses, and 0 otherwise.

	max_clause_length = 0;
	min_clause_length = 100000000;

	// size of the instance
	num_vars = 0;	 // var index from 1 to num_vars
	num_clauses = 0; // clause index from 0 to num_clauses-1
	num_hclauses = 0;
	num_sclauses = 0;

	print_time = 240;
	cutoff_time = 300;
	ans = 0.2; // hhscore
	mark = 0;
	mark1 = 0;
	obj_300 = 0;
	obj_1800 = 0;
	obj_3600 = 0;
	num_vars_opb = 0;
}

void Satlike::settings()
{
	// steps
	max_tries = 100000000;
	// max_tries = 4;
	max_flips = 200000000;
	// max_flips = 10000;
	// max_flips = 30;
	max_non_improve_flip = 10000000;

	large_clause_count_threshold = 0;
	soft_large_clause_count_threshold = 0;

	rdprob = 0.01;
	hd_count_threshold = 15;
	rwprob = 0.1;
	smooth_probability = 0.01;

	h_inc = 1;
	softclause_weight_threshold = 1;
	max_clause_score = 1000;

	/*para*/
	MinStep = 100000;
	delta = 100;
	MinHard = 10;
	gamma = 0.05;
	MaxHard = 50;
	Lopt = 50;

	lambda = 0.9;
	backward_step = 21; // para d
	NumSample = 20;

	beta = 100;

	if (num_vars > 2000)
	{
		rdprob = 0.01;
		hd_count_threshold = 50; // 50
		rwprob = 0.1;
		smooth_probability = 0.0000001;
	}
}

void Satlike::allocate_memory()
{
	int malloc_var_length = num_vars + 10;
	int malloc_clause_length = num_clauses + 10;

	unit_clause = new lit[malloc_clause_length];

	var_lit = new lit *[malloc_var_length];
	var_lit_count = new int[malloc_var_length];
	clause_lit = new lit *[malloc_clause_length];
	clause_lit_count = new int[malloc_clause_length];
	clause_true_lit_thres = new int[malloc_clause_length];
	clause_visied_times = new int[malloc_clause_length];

	hhscore = new double[malloc_var_length]; // hhscore
	clause_max_weight = new int[malloc_clause_length];
	clause_total_sum = new long long[malloc_clause_length];
	gap1 = new int[malloc_clause_length]; // gap1
	best_soln_300 = new int[malloc_var_length];
	best_soln_1800 = new int[malloc_var_length];
	opb = new int[malloc_var_length];

	turb_hardunsat_stack = new int[malloc_clause_length]; // turb
	var_mark = new int[malloc_var_length];
	is_selected_clauses = new int[malloc_clause_length];
	selected_var = new int[malloc_var_length];
	turb_best_soln = new int[malloc_var_length];

	score_baocun = new double[10]; // fps
	sscore_baocun = new double[10];
	hhscore_baocun = new double[10];
	score2 = new double[malloc_var_length];
	sscore2 = new double[malloc_var_length];
	hhscore2 = new double[malloc_var_length];
	sat_count2 = new int[malloc_clause_length];
	goodvar_stack2 = new int[malloc_var_length];

	clause_score = new double[malloc_clause_length]; // MAB-s
	selected_clauses = new int[malloc_clause_length];
	selected_times = new int[malloc_clause_length];
	sampled_clauses = new int[malloc_clause_length];

	selected_clauses_hard = new int[malloc_clause_length]; // MAB-h
	clause_hard_score = new double[malloc_clause_length];
	selected_times_hard = new int[malloc_clause_length];
	sampled_clauses_hard = new int[malloc_clause_length];

	tune_soft_clause_weight = new double[malloc_clause_length]; // NuPBO
	tuned_degree_unit_weight = new double[malloc_clause_length];
	soft_clause_num_index = new int[malloc_clause_length];
	hard_clause_num_index = new int[malloc_clause_length];
	avg_clause_coe = new double[malloc_clause_length]();

	score = new double[malloc_var_length];
	sscore = new double[malloc_var_length];
	oscore = new long long[malloc_var_length];
	var_neighbor = new int *[malloc_var_length];
	var_neighbor_count = new int[malloc_var_length];
	time_stamp = new long long[malloc_var_length];
	neighbor_flag = new int[malloc_var_length];
	temp_neighbor = new int[malloc_var_length];

	org_clause_weight = new long long[malloc_clause_length];
	clause_weight = new long long[malloc_clause_length];
	unit_weight = new double[malloc_clause_length];
	org_unit_weight = new long long[malloc_clause_length];
	sat_count = new int[malloc_clause_length];
	sat_var = new int[malloc_clause_length];
	clause_selected_count = new long long[malloc_clause_length];
	best_soft_clause = new int[malloc_clause_length];

	hardunsat_stack = new int[malloc_clause_length];
	index_in_hardunsat_stack = new int[malloc_clause_length];
	softunsat_stack = new int[malloc_clause_length];
	index_in_softunsat_stack = new int[malloc_clause_length];

	unsatvar_stack = new int[malloc_var_length];
	index_in_unsatvar_stack = new int[malloc_var_length];
	unsat_app_count = new int[malloc_var_length];

	goodvar_stack = new int[malloc_var_length];
	already_in_goodvar_stack = new int[malloc_var_length];

	cur_soln = new int[malloc_var_length];
	best_soln = new int[malloc_var_length];
	local_opt_soln = new int[malloc_var_length];

	large_weight_clauses = new int[malloc_clause_length];
	soft_large_weight_clauses = new int[malloc_clause_length];
	already_in_soft_large_weight_stack = new int[malloc_clause_length];

	best_array = new int[malloc_var_length];
	temp_lit = new int[malloc_var_length];
}

void Satlike::free_memory()
{
	int i;
	for (i = 0; i < num_clauses; i++)
		delete[] clause_lit[i];

	for (i = 1; i <= num_vars; ++i)
	{
		delete[] var_lit[i];
		delete[] var_neighbor[i];
	}
	/*hhscore*/
	delete[] hhscore;
	delete[] clause_max_weight;
	delete[] clause_total_sum;
	delete[] gap1;
	delete[] best_soln_300;
	delete[] best_soln_1800;
	delete[] opb;
	// delete[] same_big_score;
	/*hhscore*/

	delete[] turb_hardunsat_stack; // turb
	delete[] var_mark;
	delete[] is_selected_clauses;
	delete[] selected_var;
	delete[] turb_best_soln;

	delete[] score_baocun; // fps
	delete[] sscore_baocun;
	delete[] hhscore_baocun;
	delete[] score2;
	delete[] sscore2;
	delete[] hhscore2;
	delete[] sat_count2;
	delete[] goodvar_stack2;

	delete[] tune_soft_clause_weight; // NuPBO
	delete[] tuned_degree_unit_weight;
	delete[] soft_clause_num_index;
	delete[] hard_clause_num_index;
	delete[] avg_clause_coe;

	delete[] clause_score; // MAB-s
	delete[] selected_clauses;
	delete[] selected_times;
	delete[] sampled_clauses;

	delete[] selected_clauses_hard; // MAB-h
	delete[] clause_hard_score;
	delete[] selected_times_hard;
	delete[] sampled_clauses_hard;

	delete[] var_lit;
	delete[] var_lit_count;
	delete[] clause_lit;
	delete[] clause_lit_count;
	delete[] clause_true_lit_thres;
	delete[] clause_visied_times;

	delete[] score;
	delete[] oscore;
	delete[] sscore;
	delete[] var_neighbor;
	delete[] var_neighbor_count;
	delete[] time_stamp;
	delete[] neighbor_flag;
	delete[] temp_neighbor;

	delete[] org_clause_weight;
	delete[] clause_weight;
	delete[] unit_weight;
	delete[] org_unit_weight;
	delete[] sat_count;
	delete[] sat_var;
	delete[] clause_selected_count;
	delete[] best_soft_clause;

	delete[] hardunsat_stack;
	delete[] index_in_hardunsat_stack;
	delete[] softunsat_stack;
	delete[] index_in_softunsat_stack;

	delete[] unsatvar_stack;
	delete[] index_in_unsatvar_stack;
	delete[] unsat_app_count;

	delete[] goodvar_stack;
	delete[] already_in_goodvar_stack;

	// delete [] fix;
	delete[] cur_soln;
	delete[] best_soln;
	delete[] local_opt_soln;

	delete[] large_weight_clauses;
	delete[] soft_large_weight_clauses;
	delete[] already_in_soft_large_weight_stack;

	delete[] best_array;
	delete[] temp_lit;
}

void Satlike::build_neighbor_relation()
{
	int i, j, count;
	int v, c, n;
	int temp_neighbor_count;

	for (v = 1; v <= num_vars; ++v)
	{
		neighbor_flag[v] = 1;
		temp_neighbor_count = 0;

		for (i = 0; i < var_lit_count[v]; ++i)
		{
			c = var_lit[v][i].clause_num;
			for (j = 0; j < clause_lit_count[c]; ++j)
			{
				n = clause_lit[c][j].var_num;
				if (neighbor_flag[n] != 1)
				{
					neighbor_flag[n] = 1;
					temp_neighbor[temp_neighbor_count++] = n;
				}
			}
		}

		neighbor_flag[v] = 0;

		var_neighbor[v] = new int[temp_neighbor_count];
		var_neighbor_count[v] = temp_neighbor_count;

		count = 0;
		for (i = 0; i < temp_neighbor_count; i++)
		{
			var_neighbor[v][count++] = temp_neighbor[i];
			neighbor_flag[temp_neighbor[i]] = 0;
		}
	}
}

void Satlike::build_instance(char *filename)
{
	istringstream iss;
	char line[1024];
	string line2;
	char tempstr1[10];
	char tempstr2[10];
	int cur_lit;
	int i, v, c;
	// int     temp_lit[MAX_VARS];

	ifstream infile(filename);
	if (!infile)
	{
		cout << "c the input filename " << filename << " is invalid, please input the correct filename." << endl;
		exit(-1);
	}

	/*** build problem data structures of the instance ***/

	getline(infile, line2);

	while (line2[0] != 'p')
	{
		getline(infile, line2);
	}
	for (i = 0; i < 1024; i++)
	{
		line[i] = line2[i];
	}
	int read_items;
	read_items = sscanf(line, "%s %s %d %d %lld", tempstr1, tempstr2, &num_vars, &num_clauses, &top_clause_weight);

	allocate_memory(); // 加入了hhscore

	for (c = 0; c < num_clauses; c++)
	{
		clause_lit_count[c] = 0;
		clause_true_lit_thres[c] = 1;
		clause_lit[c] = NULL;
	}

	for (v = 1; v <= num_vars; ++v)
	{
		var_lit_count[v] = 0;
		var_lit[v] = NULL;
		var_neighbor[v] = NULL;
	}

	num_hclauses = num_sclauses = 0;
	unit_clause_count = 0;
	// Now, read the clauses, one at a time.
	int lit_redundent, clause_redundent;
	int *temp_weight = new int[num_vars];
	int cur_weight;
	total_soft_weight = 0;

	c = 0;
	while (getline(infile, line2))
	{
		if (line2[0] == 'c')
			continue;
		else
		{
			iss.clear();
			iss.str(line2);
			iss.seekg(0, ios::beg);
		}
		clause_redundent = 0;
		clause_lit_count[c] = 0;

		int max_weight = 0; // find max_weight

		iss >> org_clause_weight[c];
		iss >> clause_true_lit_thres[c];

		if (clause_true_lit_thres[c] <= 0) // NuPBO
		{
			// cout << "c ----" << c << ": " << clause_true_lit_thres[c] << endl;
			num_clauses--;
			continue;
		}

		if (org_clause_weight[c] != top_clause_weight)
		{
			total_soft_weight += org_clause_weight[c];
			// num_sclauses++; //privious-DeepOpt-v1
			soft_clause_num_index[num_sclauses++] = c; // NuPBO
		}
		else
		{
			hard_clause_num_index[num_hclauses++] = c; // NuPBO
		}
		// num_hclauses++; //privious-DeepOpt-v1

		iss >> cur_weight;
		iss >> cur_lit;
		while (cur_weight != 0)
		{
			temp_weight[clause_lit_count[c]] = cur_weight;

			/*find max_weight*/
			if (cur_weight > max_weight)
			{
				max_weight = cur_weight;
				// clause_max_weight[c] = max_weight;
			}
			/*find max_weight*/

			temp_lit[clause_lit_count[c]] = cur_lit;
			clause_lit_count[c]++;
			//}
			iss >> cur_weight;
			iss >> cur_lit;
		}
		// sort(temp_weight, temp_weight + clause_lit_count[c], compare);//NuPBO

		clause_max_weight[c] = max_weight;
		// if(clause_redundent==0) //the clause is not tautology
		//{
		clause_lit[c] = new lit[clause_lit_count[c] + 1];

		for (i = 0; i < clause_lit_count[c]; ++i)
		{
			clause_lit[c][i].clause_num = c;
			clause_lit[c][i].var_num = abs(temp_lit[i]);
			clause_lit[c][i].weight = temp_weight[i];
			avg_clause_coe[c] += double(clause_lit[c][i].weight);

			if (temp_lit[i] > 0)
				clause_lit[c][i].sense = 1;
			else
				clause_lit[c][i].sense = 0;

			var_lit_count[clause_lit[c][i].var_num]++;
		}

		avg_clause_coe[c] = round(double(avg_clause_coe[c] / (double)clause_lit_count[c]));
		if (avg_clause_coe[c] < 1)
			avg_clause_coe[c] = 1;

		clause_lit[c][i].var_num = 0;
		clause_lit[c][i].clause_num = -1;
		clause_lit[c][i].weight = 0;

		if (clause_lit_count[c] == 1)
		{
			unit_clause[unit_clause_count++] = clause_lit[c][0];
		}

		if (clause_lit_count[c] > max_clause_length)
			max_clause_length = clause_lit_count[c];
		if (clause_lit_count[c] < min_clause_length)
			min_clause_length = clause_lit_count[c];

		c++;
	}
	infile.close();

	// creat var literal arrays
	for (v = 1; v <= num_vars; ++v)
	{
		var_lit[v] = new lit[var_lit_count[v] + 1];
		var_lit_count[v] = 0; // reset to 0, for build up the array
	}
	// scan all clauses to build up var literal arrays
	for (c = 0; c < num_clauses; ++c)
	{
		for (i = 0; i < clause_lit_count[c]; ++i)
		{
			v = clause_lit[c][i].var_num;
			var_lit[v][var_lit_count[v]] = clause_lit[c][i];
			++var_lit_count[v];
		}
		clause_visied_times[c] = 0;
	}
	for (v = 1; v <= num_vars; ++v)
		var_lit[v][var_lit_count[v]].clause_num = -1;

	build_neighbor_relation();

	best_soln_feasible = 0;
	opt_unsat_weight = total_soft_weight + 1;
}

void Satlike::init(vector<int> &init_solution)
{
	int v, c;
	int j;
	float initsoftw = 0;

	local_times = 0; // MAB
	if_exceed = 0;
	if_exceed_hard = 0;
	hard_unsat_weight = 0;

	local_soln_feasible = 0;
	local_opt_unsat_weight = top_clause_weight + 1;
	large_weight_clauses_count = 0;
	soft_large_weight_clauses_count = 0;

	ave_soft_weight = total_soft_weight / num_sclauses;
	ave_hard_weight = 0;
	inc_hard_weight = 0;
	// cout << "ave soft weight is " << ave_soft_weight << endl;

	double tmp_avg_soft_clause_weight = 0.0;
	tmp_avg_soft_clause_weight = round(double(top_clause_weight - 1) / num_sclauses);
	if (tmp_avg_soft_clause_weight < 1)
		tmp_avg_soft_clause_weight = 1;

	// Initialize clause information
	for (c = 0; c < num_clauses; c++)
	{
		selected_times[c] = 0;		// MAB-s
		clause_score[c] = 1;		// MAB-s
		selected_times_hard[c] = 0; // MAB-hselectd_times_hard
		clause_hard_score[c] = 1;	// MAB-h

		clause_visied_times[c] = 0;
		clause_selected_count[c] = 0;

		if (org_clause_weight[c] == top_clause_weight)
		{
			// clause_weight[c] = clause_true_lit_thres[c];
			org_unit_weight[c] = 1;
			unit_weight[c] = org_unit_weight[c];
			tuned_degree_unit_weight[c] = double(unit_weight[c]) / avg_clause_coe[c];
			long long total_sum = 0;
			for (int i = 0; i < clause_lit_count[c]; ++i)
			{
				total_sum += clause_lit[c][i].weight;
			}
			clause_weight[c] = total_sum / clause_lit_count[c];
			ave_hard_weight += clause_weight[c];
			clause_total_sum[c] = total_sum;
		}
		else
		{
			tune_soft_clause_weight[c] = double(org_clause_weight[c] / tmp_avg_soft_clause_weight);
			unit_weight[c] = initsoftw;

			clause_weight[c] = org_clause_weight[c];
			clause_total_sum[c] = org_clause_weight[c];
			org_unit_weight[c] = ceil((double)clause_weight[c] / (double)clause_true_lit_thres[c]);
			// unit_weight[c] = org_unit_weight[c];
		}
		/********min{k + amax, asum}**********/
		if (clause_true_lit_thres[c] + clause_max_weight[c] <= clause_total_sum[c])
		{
			gap1[c] = clause_true_lit_thres[c] + clause_max_weight[c];
		}
		else
			gap1[c] = clause_total_sum[c];
	}
	inc_hard_weight = ave_hard_weight % num_hclauses;
	ave_hard_weight /= num_hclauses;
	inc_hard_weight += ave_soft_weight;

	if (init_solution.size() == 0)
	{
		for (v = 1; v <= num_vars; v++)
		{
			cur_soln[v] = 0;
			time_stamp[v] = 0;
			unsat_app_count[v] = 0;
		}
	}
	else
	{
		for (v = 1; v <= num_vars; v++)
		{
			cur_soln[v] = init_solution[v];
			if (cur_soln[v] == 2)
				cur_soln[v] = rand() % 2;
			time_stamp[v] = 0;
			unsat_app_count[v] = 0;
		}
	}

	// init stacks
	hard_unsat_nb = 0;
	hardunsat_stack_fill_pointer = 0;
	softunsat_stack_fill_pointer = 0;
	unsatvar_stack_fill_pointer = 0;
	/* figure out sat_count, sat_var, soft_unsat_weight and init unsat_stack */
	soft_unsat_weight = 0;

	for (c = 0; c < num_clauses; ++c)
	{
		sat_count[c] = 0;
		for (j = 0; j < clause_lit_count[c]; ++j)
		{
			if (cur_soln[clause_lit[c][j].var_num] == clause_lit[c][j].sense)
			{
				sat_count[c] += clause_lit[c][j].weight;
				sat_var[c] = clause_lit[c][j].var_num;
			}
		}
		if (sat_count[c] < clause_true_lit_thres[c])
		{
			if (org_clause_weight[c] != top_clause_weight)
				soft_unsat_weight += (clause_true_lit_thres[c] - sat_count[c]) * org_unit_weight[c];
			else
				hard_unsat_weight += clause_true_lit_thres[c] - sat_count[c]; // zyj
			unsat(c);
		}
		// cout<<"soft_unsat_weight "<<soft_unsat_weight<<endl;
	}

	/*figure out score*/
	int sense, weight;

	for (v = 1; v <= num_vars; v++)
	{
		score[v] = 0;
		sscore[v] = 0;
		hhscore[v] = 0;
		for (int i = 0; i < var_lit_count[v]; ++i)
		{
			c = var_lit[v][i].clause_num;
			sense = var_lit[v][i].sense;
			weight = var_lit[v][i].weight;
			if (org_clause_weight[c] == top_clause_weight)
			{
				if (sat_count[c] < clause_true_lit_thres[c])
				{
					if (sense != cur_soln[v])
					{
						score[v] += double(tuned_degree_unit_weight[c] * min(clause_true_lit_thres[c] - sat_count[c], weight));
						// hhscore[v] += unit_weight[c] * max(weight - (clause_true_lit_thres[c] - sat_count[c]), 0);
						hhscore[v] += 1 * max(weight - (clause_true_lit_thres[c] - sat_count[c]), 0);
					}
					else
					{
						score[v] -= double(tuned_degree_unit_weight[c] * weight);
						hhscore[v] += 0;
					}
				}
				else if (sat_count[c] >= clause_true_lit_thres[c])
				{
					if (sat_count[c] <= gap1[c])
					{
						if (sense == cur_soln[v])
						{
							score[v] -= double(tuned_degree_unit_weight[c] * max(0, clause_true_lit_thres[c] - sat_count[c] + weight));
							// hhscore[v] -= unit_weight[c] * min(weight, sat_count[c] - clause_true_lit_thres[c]);
							hhscore[v] -= 1 * min(weight, sat_count[c] - clause_true_lit_thres[c]);
						}
						else
						{
							// hhscore[v] += unit_weight[c] * min(weight, gap1[c] - sat_count[c]);
							hhscore[v] += 1 * min(weight, gap1[c] - sat_count[c]);
						}
					}
					else if (sat_count[c] > gap1[c])
					{
						if (sense == cur_soln[v])
						{
							score[v] -= double(tuned_degree_unit_weight[c] * max(0, clause_true_lit_thres[c] - sat_count[c] + weight));
							// hhscore[v] -= unit_weight[c] * max(weight - (sat_count[c] - gap1[c]), 0);
							hhscore[v] -= 1 * max(weight - (sat_count[c] - gap1[c]), 0);
						}
					}
				}
			}
			else
			{
				if (sat_count[c] < clause_true_lit_thres[c])
				{
					if (sense != cur_soln[v])
					{
						sscore[v] += unit_weight[c] * tune_soft_clause_weight[c];
					}
					else
						sscore[v] -= unit_weight[c] * tune_soft_clause_weight[c];
				}
				else if (sat_count[c] >= clause_true_lit_thres[c])
				{
					if (sense == cur_soln[v])
					{
						sscore[v] -= unit_weight[c] * tune_soft_clause_weight[c];
					}
				}
			}
		}
	}

	// init goodvars stack
	goodvar_stack_fill_pointer = 0;

	for (v = 1; v <= num_vars; v++)
	{
		// if (score[v] + sscore[v] + ans * hhscore[v] > 0)
		if (score[v] + sscore[v] > 0)
		{
			already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
			mypush(v, goodvar_stack);
		}
		else
			already_in_goodvar_stack[v] = -1;
	}
}

int Satlike::pick_softc()
{
	int sel_c;

	sampled_clauses[0] = softunsat_stack[rand() % softunsat_stack_fill_pointer];
	double min = clause_score[sampled_clauses[0]], max = clause_score[sampled_clauses[0]];

	for (int i = 1; i < NumSample; i++)
	{
		sampled_clauses[i] = softunsat_stack[rand() % softunsat_stack_fill_pointer];

		if (clause_score[sampled_clauses[i]] < min)
			min = clause_score[sampled_clauses[i]];
		if (clause_score[sampled_clauses[i]] > max)
			max = clause_score[sampled_clauses[i]];
	}
	if (max == min)
	{
		sel_c = sampled_clauses[0];
		for (int i = 1; i < NumSample; i++)
		{
			if (selected_times[sampled_clauses[i]] < selected_times[sel_c])
				sel_c = sampled_clauses[i];
			else if (selected_times[sampled_clauses[i]] == selected_times[sel_c])
			{
				if (clause_weight[sampled_clauses[i]] > clause_weight[sel_c])
					sel_c = sampled_clauses[i];
			}
		}
	}
	else
	{
		double max_value = clause_score[sampled_clauses[0]] + 1.0 * sqrt((log(local_times + 1)) / ((double)(selected_times[sampled_clauses[0]] + 1)));
		sel_c = sampled_clauses[0];
		for (int i = 1; i < NumSample; i++)
		{
			double dtemp = clause_score[sampled_clauses[i]] + 1.0 * sqrt((log(local_times + 1)) / ((double)(selected_times[sampled_clauses[i]] + 1)));
			if (dtemp > max_value)
			{
				max_value = dtemp;
				sel_c = sampled_clauses[i];
			}
		}
	}
	selected_times[sel_c]++;
	selected_clauses[local_times % backward_step] = sel_c;
	if (local_times > 0)
	{
		long long s = pre_unsat_weight - soft_unsat_weight;
		update_clause_scores(s);
	}
	pre_unsat_weight = soft_unsat_weight;
	local_times++;

	return sel_c;
}

void Satlike::update_clause_scores(long long s)
{
	int i, j;
	double stemp;

	long long opt = opt_unsat_weight;
	if (soft_unsat_weight < opt)
		opt = soft_unsat_weight;

	stemp = ((double)s) / (pre_unsat_weight - opt + 1);

	double discount;
	if (local_times < backward_step)
	{
		for (i = 0; i < local_times; i++)
		{
			discount = pow(lambda, local_times - 1 - i);
			clause_score[selected_clauses[i]] += (discount * ((double)stemp));
			if (abs(clause_score[selected_clauses[i]]) > max_clause_score)
				if_exceed = 1;
		}
	}
	else
	{
		for (i = 0; i < backward_step; i++)
		{
			if (i == local_times % backward_step)
				continue;
			if (i < local_times % backward_step)
				discount = pow(lambda, local_times % backward_step - 1 - i);
			else
				discount = pow(lambda, local_times % backward_step + backward_step - 1 - i);
			clause_score[selected_clauses[i]] += (discount * ((double)stemp));
			if (abs(clause_score[selected_clauses[i]]) > max_clause_score)
				if_exceed = 1;
		}
	}
	if (if_exceed)
	{
		for (i = 0; i < num_clauses; i++)
			clause_score[i] = clause_score[i] / 2.0;
		if_exceed = 0;
	}
}

int Satlike::pick_hardc()
{
	int sel_c;

	sampled_clauses_hard[0] = hardunsat_stack[rand() % hardunsat_stack_fill_pointer];
	double min = clause_hard_score[sampled_clauses_hard[0]], max = clause_hard_score[sampled_clauses_hard[0]];

	for (int i = 1; i < NumSample; i++)
	{
		sampled_clauses_hard[i] = hardunsat_stack[rand() % hardunsat_stack_fill_pointer];

		if (clause_hard_score[sampled_clauses_hard[i]] < min)
			min = clause_hard_score[sampled_clauses_hard[i]];
		if (clause_hard_score[sampled_clauses_hard[i]] > max)
			max = clause_hard_score[sampled_clauses_hard[i]];
	}
	if (max == min)
	{
		sel_c = sampled_clauses_hard[0];
		for (int i = 1; i < NumSample; i++)
		{
			if (selected_times_hard[sampled_clauses_hard[i]] < selected_times_hard[sel_c])
				sel_c = sampled_clauses_hard[i];
			else if (selected_times_hard[sampled_clauses_hard[i]] == selected_times_hard[sel_c])
			{
				// if (clause_weight[sampled_clauses[i]] > clause_weight[sel_c])
				if (rand() % 100 < 50)
					sel_c = sampled_clauses_hard[i];
			}
		}
	}
	else
	{
		double max_value = clause_hard_score[sampled_clauses_hard[0]] + 1.0 * sqrt((log(local_times_hard + 1)) / ((double)(selected_times_hard[sampled_clauses_hard[0]] + 1)));
		sel_c = sampled_clauses_hard[0];
		for (int i = 1; i < NumSample; i++)
		{
			double dtemp = clause_hard_score[sampled_clauses_hard[i]] + 1.0 * sqrt((log(local_times_hard + 1)) / ((double)(selected_times_hard[sampled_clauses_hard[i]] + 1)));
			if (dtemp > max_value)
			{
				max_value = dtemp;
				sel_c = sampled_clauses_hard[i];
			}
		}
	}
	selected_times_hard[sel_c]++;
	selected_clauses_hard[local_times_hard % backward_step] = sel_c;
	if (local_times_hard > 0)
	{
		long long s = pre_hard_unsat_weight - hard_unsat_weight;
		update_clause_scores_hard(s);
	}
	pre_hard_unsat_weight = hard_unsat_weight;
	local_times_hard++;

	return sel_c;
}

void Satlike::update_clause_scores_hard(long long s)
{
	int i, j;
	double stemp;

	stemp = ((double)s) / (pre_hard_unsat_weight + 1);

	double discount;
	if (local_times_hard < backward_step)
	{
		for (i = 0; i < local_times_hard; i++)
		{
			discount = pow(lambda, local_times_hard - 1 - i);
			clause_hard_score[selected_clauses_hard[i]] += (discount * ((double)stemp));
			if (abs(clause_hard_score[selected_clauses_hard[i]]) > 1000)
				if_exceed_hard = 1;
		}
	}
	else
	{
		for (i = 0; i < backward_step; i++)
		{
			if (i == local_times_hard % backward_step)
				continue;
			if (i < local_times_hard % backward_step)
				discount = pow(lambda, local_times_hard % backward_step - 1 - i);
			else
				discount = pow(lambda, local_times_hard % backward_step + backward_step - 1 - i);
			clause_hard_score[selected_clauses_hard[i]] += (discount * ((double)stemp));
			if (abs(clause_hard_score[selected_clauses_hard[i]]) > 1000)
				if_exceed_hard = 1;
		}
	}
	if (if_exceed_hard)
	{
		for (i = 0; i < num_clauses; i++)
			clause_hard_score[i] = clause_hard_score[i] / 2.0;
		if_exceed_hard = 0;
	}
}

void Satlike::pick_vars()
{
	int i, v, j;
	int best_var;
	int rand_select;
	int mark_v1 = 0;
	int mark_v2 = 0;
	// int better_var1;

	update_clause_weights();
	int sel_c;
	int sel_c_set[110];
	lit *p;
	int mark_which = 0;
	if (hardunsat_stack_fill_pointer > 0)
	{
		// if (step < 10000 || best_soln_feasible > 0) // MAB-h
		if (best_soln_feasible > 0) // MAB-h
		{
			if (hardunsat_stack_fill_pointer < beta)
			{
				sel_c = hardunsat_stack[rand() % hardunsat_stack_fill_pointer];
				mark_which = 1;
			}
			else
			{
				sel_c = hardunsat_stack[rand() % hardunsat_stack_fill_pointer];
				for (i = 0; i < beta; i++)
				{
					sel_c_set[i] = hardunsat_stack[rand() % hardunsat_stack_fill_pointer];
				}
				mark_which = 2;
			}
		}
		else
		{
			if ((rand() % MY_RAND_MAX_INT) * BASIC_SCALE < rwprob)
			{
				sel_c = hardunsat_stack[rand() % hardunsat_stack_fill_pointer];
			}
			else
				sel_c = pick_hardc();
			mark_which = 1;
		}
	}
	else
	{
		// sel_c = softunsat_stack[rand() % softunsat_stack_fill_pointer];
		sel_c = pick_softc();
		mark_which = 3;
	}

	if ((rand() % MY_RAND_MAX_INT) * BASIC_SCALE < rwprob)
	{
		best_var = clause_lit[sel_c][rand() % clause_lit_count[sel_c]].var_num;
		flip(best_var);
		return;
	}

	if (clause_lit_count[sel_c] == 1)
	{
		best_var = clause_lit[sel_c][0].var_num;
		flip(best_var);
		return;
	}
	else if (clause_lit_count[sel_c] == 2)
	{
		int v1 = clause_lit[sel_c][0].var_num;
		int v2 = clause_lit[sel_c][1].var_num;
		if (rand() % 100 < 0)
		{
			if (rand() % 100 < 50)
				best_var = v1;
			else
				best_var = v2;
		}
		else
		{
			score_baocun[0] = score[v1];
			sscore_baocun[0] = sscore[v1];
			hhscore_baocun[0] = hhscore[v1];
			score_baocun[1] = score[v2];
			sscore_baocun[1] = sscore[v2];
			hhscore_baocun[1] = hhscore[v2];
			for (i = 0; i < num_vars; i++)
			{
				goodvar_stack2[i] = 0;
			}

			int best_v1_neighbor;
			int best_v2_neighbor;
			flip_fps(v1);
			if (goodvar_stack2_num > 0)
			{
				mark_v1 = 1;
				if (goodvar_stack2_num < 15)
				{
					best_v1_neighbor = goodvar_stack2[0];
					for (j = 1; j < goodvar_stack2_num; ++j)
					{
						v = goodvar_stack2[j];
						if (score2[v] + sscore2[v] > score2[best_v1_neighbor] + sscore2[best_v1_neighbor])
							best_v1_neighbor = v;
						else if (score2[v] + sscore2[v] == score2[best_v1_neighbor] + sscore2[best_v1_neighbor])
						{
							if (hhscore2[v] > hhscore2[best_v1_neighbor])
							{
								best_v1_neighbor = v;
							}
							else if (hhscore2[v] == hhscore2[best_v1_neighbor])
							{
								rand_select = rand() % 2;
								if (rand_select == 1)
								{
									best_v1_neighbor = v;
								}
							}
						}
					}
				}
				else
				{
					best_v1_neighbor = goodvar_stack2[rand() % goodvar_stack2_num];
					for (j = 1; j < 15; ++j)
					{
						v = goodvar_stack2[rand() % goodvar_stack2_num];
						if (score2[v] + sscore2[v] > score2[best_v1_neighbor] + sscore2[best_v1_neighbor])
							best_v1_neighbor = v;
						else if (score2[v] + sscore2[v] == score2[best_v1_neighbor] + sscore2[best_v1_neighbor])
						{
							if (hhscore2[v] > hhscore2[best_v1_neighbor])
							{
								best_v1_neighbor = v;
							}
							else if (hhscore2[v] == hhscore2[best_v1_neighbor])
							{
								rand_select = rand() % 2;
								if (rand_select == 1)
								{
									best_v1_neighbor = v;
								}
							}
						}
					}
				}
				score_baocun[0] += score2[best_v1_neighbor];
				sscore_baocun[0] += sscore2[best_v1_neighbor];
				hhscore_baocun[0] += hhscore2[best_v1_neighbor];
			}
			else
			{
				mark_v1 = 0;
			}

			for (i = 0; i < num_vars; i++)
			{
				goodvar_stack2[i] = 0;
			}

			flip_fps(v2);
			if (goodvar_stack2_num > 0)
			{
				mark_v2 = 1;
				if (goodvar_stack2_num < 15)
				{
					best_v2_neighbor = goodvar_stack2[0];
					for (j = 1; j < goodvar_stack2_num; ++j)
					{
						v = goodvar_stack2[j];
						if (score2[v] + sscore2[v] > score2[best_v2_neighbor] + sscore2[best_v2_neighbor])
							best_v2_neighbor = v;
						else if (score2[v] + sscore2[v] == score2[best_v2_neighbor] + sscore2[best_v2_neighbor])
						{
							if (hhscore2[v] > hhscore2[best_v2_neighbor])
							{
								best_v2_neighbor = v;
							}
							else if (hhscore2[v] == hhscore2[best_v2_neighbor])
							{
								rand_select = rand() % 2;
								if (rand_select == 1)
								{
									best_v2_neighbor = v;
								}
							}
						}
					}
				}
				else
				{
					best_v2_neighbor = goodvar_stack2[rand() % goodvar_stack2_num];
					for (j = 1; j < 15; ++j)
					{
						v = goodvar_stack2[rand() % goodvar_stack2_num];
						if (score2[v] + sscore2[v] > score2[best_v2_neighbor] + sscore2[best_v2_neighbor])
							best_v2_neighbor = v;
						else if (score2[v] + sscore2[v] == score2[best_v2_neighbor] + sscore2[best_v2_neighbor])
						{
							if (hhscore2[v] > hhscore2[best_v2_neighbor])
							{
								best_v2_neighbor = v;
							}
							else if (hhscore2[v] == hhscore2[best_v2_neighbor])
							{
								rand_select = rand() % 2;
								if (rand_select == 1)
								{
									best_v2_neighbor = v;
								}
							}
						}
					}
				}
				score_baocun[1] += score2[best_v2_neighbor];
				sscore_baocun[1] += sscore2[best_v2_neighbor];
				hhscore_baocun[1] += hhscore2[best_v2_neighbor];
			}
			else
			{
				mark_v2 = 0;
			}
			// int best_flip_neighbor;
			// int best_flip;
			if (mark_v1 == 1 && mark_v2 == 1)
			{
				if (score_baocun[0] + sscore_baocun[0] > 0 && score_baocun[1] + sscore_baocun[1] > 0)
				{
					if (score_baocun[0] + sscore_baocun[0] > score_baocun[1] + sscore_baocun[1])
					{
						best_flip_neighbor = best_v1_neighbor;
						best_var = v1;
					}
					else if (score_baocun[0] + sscore_baocun[0] < score_baocun[1] + sscore_baocun[1])
					{
						best_flip_neighbor = best_v2_neighbor;
						best_var = v2;
					}
					else
					{
						if (hhscore_baocun[0] > hhscore_baocun[1])
						{
							best_flip_neighbor = best_v1_neighbor;
							best_var = v1;
						}
						else if (hhscore_baocun[0] < hhscore_baocun[1])
						{
							best_flip_neighbor = best_v2_neighbor;
							best_var = v2;
						}
						else
						{
							rand_select = rand() % 2;
							if (rand_select == 1)
							{
								best_var = v1;
								best_flip_neighbor = best_v1_neighbor;
							}
							else
							{
								best_var = v2;
								best_flip_neighbor = best_v2_neighbor;
							}
						}
					}
					flip(best_var);
					flip(best_flip_neighbor);
					return;
				}
				else if (score_baocun[0] + sscore_baocun[0] > 0 && score_baocun[1] + sscore_baocun[1] < 0)
				{
					best_flip_neighbor = best_v1_neighbor;
					best_var = v1;
					flip(best_var);
					flip(best_flip_neighbor);
					return;
				}
				else if (score_baocun[0] + sscore_baocun[0] < 0 && score_baocun[1] + sscore_baocun[1] > 0)
				{
					best_flip_neighbor = best_v2_neighbor;
					best_var = v2;
					flip(best_var);
					flip(best_flip_neighbor);
					return;
				}
				else // new addflip
				{
					if (score[v1] + sscore[v1] < score_baocun[0] + sscore_baocun[0] && score[v2] + sscore[v2] < score_baocun[1] + sscore_baocun[1])
					{
						if (score_baocun[0] + sscore_baocun[0] > score_baocun[1] + sscore_baocun[1])
						{
							best_flip_neighbor = best_v1_neighbor;
							best_var = v1;
							flip(best_var);
							flip(best_flip_neighbor);
							return;
						}
						else if (score_baocun[0] + sscore_baocun[0] < score_baocun[1] + sscore_baocun[1])
						{
							best_flip_neighbor = best_v2_neighbor;
							best_var = v2;
							flip(best_var);
							flip(best_flip_neighbor);
							return;
						}
						else // 可以用hhscore or rand or teacher
						{
							if (hhscore_baocun[0] > hhscore_baocun[1])
							{
								best_flip_neighbor = best_v1_neighbor;
								best_var = v1;
							}
							else if (hhscore_baocun[0] < hhscore_baocun[1])
							{
								best_flip_neighbor = best_v2_neighbor;
								best_var = v2;
							}
							else
							{
								rand_select = rand() % 2;
								if (rand_select == 1)
								{
									best_var = v1;
									best_flip_neighbor = best_v1_neighbor;
								}
								else
								{
									best_var = v2;
									best_flip_neighbor = best_v2_neighbor;
								}
							}
							flip(best_var);
							flip(best_flip_neighbor);
							return;
						}
					}
					else
					{
						if (score[v1] + sscore[v1] > score[v2] + sscore[v2])
						{
							best_var = v1;
						}
						else if (score[v1] + sscore[v1] < score[v2] + sscore[v2])
						{
							best_var = v2;
						}
						else
						{
							if (hhscore[v1] > hhscore[v2])
							{
								best_var = v1;
							}
							else if (hhscore[v1] < hhscore[v2])
							{
								best_var = v2;
							}
							else
							{
								if (rand() % 100 < 50)
									best_var = v1;
								else
									best_var = v2;
							}
						}
						flip(best_var);
						return;
					}
				}
			}
			else if (mark_v1 == 1 && mark_v2 == 0)
			{
				if (score_baocun[0] + sscore_baocun[0] > 0)
				{
					best_var = v1;
					best_flip_neighbor = best_v1_neighbor;
					flip(best_var);
					flip(best_flip_neighbor);
					return;
				}
				else
				{
					if (score[v1] + sscore[v1] < score_baocun[0] + sscore_baocun[0] && score[v2] + sscore[v2] < score_baocun[1] + sscore_baocun[1])
					{
						if (score_baocun[0] + sscore_baocun[0] > score_baocun[1] + sscore_baocun[1])
						{
							best_flip_neighbor = best_v1_neighbor;
							best_var = v1;
							flip(best_var);
							flip(best_flip_neighbor);
							return;
						}
						else if (score_baocun[0] + sscore_baocun[0] < score_baocun[1] + sscore_baocun[1])
						{
							best_flip_neighbor = best_v2_neighbor;
							best_var = v2;
							flip(best_var);
							flip(best_flip_neighbor);
							return;
						}
						else
						{
							if (hhscore_baocun[0] > hhscore_baocun[1])
							{
								best_flip_neighbor = best_v1_neighbor;
								best_var = v1;
							}
							else if (hhscore_baocun[0] < hhscore_baocun[1])
							{
								best_flip_neighbor = best_v2_neighbor;
								best_var = v2;
							}
							else
							{
								rand_select = rand() % 2;
								if (rand_select == 1)
								{
									best_var = v1;
									best_flip_neighbor = best_v1_neighbor;
								}
								else
								{
									best_var = v2;
									best_flip_neighbor = best_v2_neighbor;
								}
							}
							flip(best_var);
							flip(best_flip_neighbor);
							return;
						}
					}
					else
					{
						if (score[v1] + sscore[v1] > score[v2] + sscore[v2])
						{
							best_var = v1;
						}
						else if (score[v1] + sscore[v1] < score[v2] + sscore[v2])
						{
							best_var = v2;
						}
						else
						{
							if (hhscore[v1] > hhscore[v2])
							{
								best_var = v1;
							}
							else if (hhscore[v1] < hhscore[v2])
							{
								best_var = v2;
							}
							else
							{
								if (rand() % 100 < 50)
									best_var = v1;
								else
									best_var = v2;
							}
						}
						flip(best_var);
						return;
					}
				}
			}
			else if (mark_v1 == 0 && mark_v2 == 1)
			{
				if (score_baocun[1] + sscore_baocun[1] > 0)
				{
					best_var = v2;
					best_flip_neighbor = best_v2_neighbor;
					flip(best_var);
					flip(best_flip_neighbor);
					return;
				}
				else
				{
					if (score[v1] + sscore[v1] < score_baocun[0] + sscore_baocun[0] && score[v2] + sscore[v2] < score_baocun[1] + sscore_baocun[1])
					{
						if (score_baocun[0] + sscore_baocun[0] > score_baocun[1] + sscore_baocun[1])
						{
							best_flip_neighbor = best_v1_neighbor;
							best_var = v1;
							flip(best_var);
							flip(best_flip_neighbor);
							return;
						}
						else if (score_baocun[0] + sscore_baocun[0] < score_baocun[1] + sscore_baocun[1])
						{
							best_flip_neighbor = best_v2_neighbor;
							best_var = v2;
							flip(best_var);
							flip(best_flip_neighbor);
							return;
						}
						else
						{
							if (hhscore_baocun[0] > hhscore_baocun[1])
							{
								best_flip_neighbor = best_v1_neighbor;
								best_var = v1;
							}
							else if (hhscore_baocun[0] < hhscore_baocun[1])
							{
								best_flip_neighbor = best_v2_neighbor;
								best_var = v2;
							}
							else
							{
								rand_select = rand() % 2;
								if (rand_select == 1)
								{
									best_var = v1;
									best_flip_neighbor = best_v1_neighbor;
								}
								else
								{
									best_var = v2;
									best_flip_neighbor = best_v2_neighbor;
								}
							}
							flip(best_var);
							flip(best_flip_neighbor);
							return;
						}
					}
					else
					{
						if (score[v1] + sscore[v1] > score[v2] + sscore[v2])
						{
							best_var = v1;
						}
						else if (score[v1] + sscore[v1] < score[v2] + sscore[v2])
						{
							best_var = v2;
						}
						else
						{
							if (hhscore[v1] > hhscore[v2])
							{
								best_var = v1;
							}
							else if (hhscore[v1] < hhscore[v2])
							{
								best_var = v2;
							}
							else
							{
								if (rand() % 100 < 50)
									best_var = v1;
								else
									best_var = v2;
							}
						}
						flip(best_var);
						return;
					}
				}
			}
			else
			{
				if (score[v1] + sscore[v1] < score_baocun[0] + sscore_baocun[0] && score[v2] + sscore[v2] < score_baocun[1] + sscore_baocun[1])
				{
					if (score_baocun[0] + sscore_baocun[0] > score_baocun[1] + sscore_baocun[1])
					{
						best_flip_neighbor = best_v1_neighbor;
						best_var = v1;
						flip(best_var);
						flip(best_flip_neighbor);
						return;
					}
					else if (score_baocun[0] + sscore_baocun[0] < score_baocun[1] + sscore_baocun[1])
					{
						best_flip_neighbor = best_v2_neighbor;
						best_var = v2;
						flip(best_var);
						flip(best_flip_neighbor);
						return;
					}
					else
					{
						if (hhscore_baocun[0] > hhscore_baocun[1])
						{
							best_flip_neighbor = best_v1_neighbor;
							best_var = v1;
						}
						else if (hhscore_baocun[0] < hhscore_baocun[1])
						{
							best_flip_neighbor = best_v2_neighbor;
							best_var = v2;
						}
						else
						{
							rand_select = rand() % 2;
							if (rand_select == 1)
							{
								best_var = v1;
								best_flip_neighbor = best_v1_neighbor;
							}
							else
							{
								best_var = v2;
								best_flip_neighbor = best_v2_neighbor;
							}
						}
						flip(best_var);
						flip(best_flip_neighbor);
						return;
					}
				}
				else
				{
					if (score[v1] + sscore[v1] > score[v2] + sscore[v2])
					{
						best_var = v1;
					}
					else if (score[v1] + sscore[v1] < score[v2] + sscore[v2])
					{
						best_var = v2;
					}
					else
					{
						if (hhscore[v1] > hhscore[v2])
						{
							best_var = v1;
						}
						else if (hhscore[v1] < hhscore[v2])
						{
							best_var = v2;
						}
						else
						{
							if (rand() % 100 < 50)
								best_var = v1;
							else
								best_var = v2;
						}
					}
					flip(best_var);
					return;
				}
			}
		}
	}
	else
	{
		/********BMS********/
		// if() {printf("wrong else\n");
		if (mark_which == 1 || mark_which == 3)
		{
			int bms;
			if (clause_lit_count[sel_c] <= 200)
			{
				bms = (clause_lit_count[sel_c] - 1) / 2 + 1;
			}
			else
			{
				bms = 100;
			}

			best_var = clause_lit[sel_c][rand() % clause_lit_count[sel_c]].var_num;
			for (i = 1; i < bms; i++)
			{
				v = clause_lit[sel_c][rand() % clause_lit_count[sel_c]].var_num;
				if (score[v] + sscore[v] > score[best_var] + sscore[best_var])
					best_var = v;
				else if (score[v] + sscore[v] == score[best_var] + sscore[best_var])
				{
					if (hhscore[v] > hhscore[best_var])
					{
						best_var = v;
					}
					else if (hhscore[v] == hhscore[best_var])
					{
						rand_select = rand() % 2;
						if (rand_select == 1)
						{
							best_var = v;
						}
					}
				}
			}
		}
		else
		{
			best_var = clause_lit[sel_c_set[0]][rand() % clause_lit_count[sel_c_set[0]]].var_num;
			for (j = 0; j < 100; j++)
			{
				int bms;
				if (clause_lit_count[sel_c] <= 200)
				{
					bms = (clause_lit_count[sel_c] - 1) / 2 + 1;
				}
				else
				{
					bms = 100;
				}
				for (i = 0; i < bms; i++)
				{
					v = clause_lit[sel_c_set[j]][rand() % clause_lit_count[sel_c_set[j]]].var_num;
					if (score[v] + sscore[v] > score[best_var] + sscore[best_var])
						best_var = v;
					else if (score[v] + sscore[v] == score[best_var] + sscore[best_var])
					{
						if (hhscore[v] > hhscore[best_var])
						{
							best_var = v;
						}
						else if (hhscore[v] == hhscore[best_var])
						{
							rand_select = rand() % 2;
							if (rand_select == 1)
							{
								best_var = v;
							}
						}
					}
				}
			}
		}
		flip(best_var);
		return;
	}
}

int Satlike::pick_var()
{
	int i, v, j;
	int best_var;
	int rand_select;

	if (goodvar_stack_fill_pointer > 0)
	{
		if ((rand() % MY_RAND_MAX_INT) * BASIC_SCALE < rdprob)
		{
			return goodvar_stack[rand() % goodvar_stack_fill_pointer];
		}

		if (goodvar_stack_fill_pointer < hd_count_threshold) // 15
		{
			best_var = goodvar_stack[0];
			for (i = 1; i < goodvar_stack_fill_pointer; ++i)
			{
				v = goodvar_stack[i];
				if (score[v] + sscore[v] > score[best_var] + sscore[best_var])
					best_var = v;
				else if (score[v] + sscore[v] == score[best_var] + sscore[best_var])
				{
					if (hhscore[v] > hhscore[best_var])
					{
						best_var = v;
					}
					else if (hhscore[v] == hhscore[best_var])
					{
						rand_select = rand() % 2;
						if (rand_select == 1)
						{
							best_var = v;
						}
					}
				}
			}
			return best_var;
		}
		else
		{
			int r = rand() % goodvar_stack_fill_pointer;
			best_var = goodvar_stack[r];

			for (i = 1; i < hd_count_threshold; ++i)
			{
				r = rand() % goodvar_stack_fill_pointer;
				v = goodvar_stack[r];
				if (score[v] + sscore[v] > score[best_var] + sscore[best_var])
					best_var = v;
				else if (score[v] + sscore[v] == score[best_var] + sscore[best_var])
				{
					if (hhscore[v] > hhscore[best_var])
					{
						best_var = v;
					}
					else if (hhscore[v] == hhscore[best_var])
					{
						rand_select = rand() % 2;
						if (rand_select == 1)
						{
							best_var = v;
						}
					}
				}
			}
			return best_var;
		}
	}

	update_clause_weights();

	int sel_c;
	int sel_c_set[100];
	lit *p;
	int mark_which = 0;

	if (hardunsat_stack_fill_pointer > 0)
	{
		if (hardunsat_stack_fill_pointer < 100)
		{
			sel_c = hardunsat_stack[rand() % hardunsat_stack_fill_pointer];
			mark_which = 1;
		}
		else
		{
			sel_c = hardunsat_stack[rand() % hardunsat_stack_fill_pointer];
			for (i = 0; i < 100; i++)
			{
				sel_c_set[i] = hardunsat_stack[rand() % hardunsat_stack_fill_pointer];
			}
			mark_which = 2;
		}
	}
	else
	{
		// sel_c = softunsat_stack[rand() % softunsat_stack_fill_pointer];
		sel_c = pick_softc(); // MAB
		mark_which = 3;
	}
	if ((rand() % MY_RAND_MAX_INT) * BASIC_SCALE < rwprob)
		return clause_lit[sel_c][rand() % clause_lit_count[sel_c]].var_num;

	if (clause_lit_count[sel_c] == 1)
	{
		best_var = clause_lit[sel_c][0].var_num;
	}
	else if (clause_lit_count[sel_c] == 2)
	{
		int v1 = clause_lit[sel_c][0].var_num;
		int v2 = clause_lit[sel_c][1].var_num;
		if (rand() % 100 < 0)
		{
			if (rand() % 100 < 50)
				best_var = v1;
			else
				best_var = v2;
		}
		else
		{

			if (score[v1] + sscore[v1] > score[v2] + sscore[v2])
			{
				best_var = v1;
			}
			else if (score[v1] + sscore[v1] < score[v2] + sscore[v2])
			{
				best_var = v2;
			}
			else
			{
				if (hhscore[v1] > hhscore[v2])
				{
					best_var = v1;
				}
				else if (hhscore[v1] < hhscore[v2])
				{
					best_var = v2;
				}
				else
				{
					int rand_select = rand() % 2;
					if (rand_select == 1)
					{
						best_var = v1;
					}
					else
					{
						best_var = v2;
					}
				}
			}
		}
	}
	else
	{
		/********BMS********/
		if (mark_which == 1 || mark_which == 3)
		{
			int bms;
			if (clause_lit_count[sel_c] <= 200)
			{
				bms = (clause_lit_count[sel_c] - 1) / 2 + 1;
			}
			else
			{
				bms = 100;
			}

			best_var = clause_lit[sel_c][rand() % clause_lit_count[sel_c]].var_num;
			for (i = 1; i < bms; i++)
			{
				v = clause_lit[sel_c][rand() % clause_lit_count[sel_c]].var_num;
				if (score[v] + sscore[v] > score[best_var] + sscore[best_var])
					best_var = v;
				else if (score[v] + sscore[v] == score[best_var] + sscore[best_var])
				{
					if (hhscore[v] > hhscore[best_var])
					{
						best_var = v;
					}
					else if (hhscore[v] == hhscore[best_var])
					{
						rand_select = rand() % 2;
						if (rand_select == 1)
						{
							best_var = v;
						}
					}
				}
			}
		}
		else
		{
			best_var = clause_lit[sel_c_set[0]][rand() % clause_lit_count[sel_c_set[0]]].var_num;
			for (j = 0; j < 100; j++)
			{
				int bms;
				if (clause_lit_count[sel_c] <= 200)
				{
					bms = (clause_lit_count[sel_c] - 1) / 2 + 1;
				}
				else
				{
					bms = 100;
				}
				for (i = 0; i < bms; i++)
				{
					v = clause_lit[sel_c_set[j]][rand() % clause_lit_count[sel_c_set[j]]].var_num;
					if (score[v] + sscore[v] > score[best_var] + sscore[best_var])
						best_var = v;
					else if (score[v] + sscore[v] == score[best_var] + sscore[best_var])
					{
						if (hhscore[v] > hhscore[best_var])
						{
							best_var = v;
						}
						else if (hhscore[v] == hhscore[best_var])
						{
							rand_select = rand() % 2;
							if (rand_select == 1)
							{
								best_var = v;
							}
						}
					}
				}
			}
		}
	}
	return best_var;
}

void Satlike::flip_fps(int flipvar)
{
	int i, v, c;
	int index;
	lit *clause_c;
	int weight;
	int gap;

	double org_flipvar_score = score[flipvar];
	double org_sscore = sscore[flipvar];
	double org_hhscore = hhscore[flipvar];
	cur_soln[flipvar] = 1 - cur_soln[flipvar];

	for (i = 0; i < var_neighbor_count[flipvar]; i++)
	{
		score2[var_neighbor[flipvar][i]] = score[var_neighbor[flipvar][i]];
		sscore2[var_neighbor[flipvar][i]] = sscore[var_neighbor[flipvar][i]];
		hhscore2[var_neighbor[flipvar][i]] = hhscore[var_neighbor[flipvar][i]];
	}

	for (i = 0; i < var_lit_count[flipvar]; i++)
	{
		c = var_lit[flipvar][i].clause_num;
		sat_count2[c] = sat_count[c];
	}

	for (i = 0; i < var_lit_count[flipvar]; ++i)
	{
		c = var_lit[flipvar][i].clause_num;
		clause_c = clause_lit[c];
		weight = var_lit[flipvar][i].weight;

		if (org_clause_weight[c] == top_clause_weight)
		{
			if (cur_soln[flipvar] == var_lit[flipvar][i].sense)
			{
				if (sat_count2[c] + weight <= clause_true_lit_thres[c])
				{
					gap = clause_true_lit_thres[c] - sat_count2[c];
					for (lit *p = clause_c; (v = p->var_num) != 0; p++)
					{
						if (p->sense != cur_soln[v])
						{
							score2[v] -= double((tuned_degree_unit_weight[c] * (min(gap, p->weight) - min(gap - weight, p->weight))));
							hhscore2[v] += (1 * (max(p->weight - gap + weight, 0) - max(p->weight - gap, 0)));
						}
					}
				}
				else if (sat_count2[c] <= clause_true_lit_thres[c])
				{
					gap = clause_true_lit_thres[c] - sat_count2[c];
					for (lit *p = clause_c; (v = p->var_num) != 0; p++)
					{
						if (p->sense != cur_soln[v])
						{
							score2[v] -= double((tuned_degree_unit_weight[c] * min(gap, p->weight)));
							// hhscore[v] += (unit_weight[c] * (min(p->weight, gap1[c] - sat_count[c] - weight) - max(p->weight - gap, 0)));//2 unsat
							hhscore2[v] += (1 * (min(p->weight, gap1[c] - sat_count2[c] - weight) - max(p->weight - gap, 0))); // 2 unsat
						}
						else
						{
							score2[v] += double(tuned_degree_unit_weight[c] * (p->weight - max(0, gap - weight + p->weight)));
							// hhscore[v] -= unit_weight[c] * min(p->weight, weight - gap);//2 sat
							hhscore2[v] -= 1 * min(p->weight, weight - gap); // 2 sat
						}
					}
				}
				else
				{
					/**************hhscore****************/
					if (sat_count2[c] <= gap1[c])
					{
						// gap = clause_true_lit_thres[c] - sat_count[c];
						if (sat_count2[c] + weight <= gap1[c])
						{
							gap = clause_true_lit_thres[c] - sat_count2[c];
							for (lit *p = clause_c; (v = p->var_num) != 0; p++)
							{
								if (p->sense == cur_soln[v])
								{
									score2[v] += double(tuned_degree_unit_weight[c] * (max(0, gap + p->weight) - max(0, gap - weight + p->weight)));
									// hhscore[v] += unit_weight[c] * (min(p->weight, -gap) - min(p->weight, weight - gap));//4
									hhscore2[v] += 1 * (min(p->weight, -gap) - min(p->weight, weight - gap)); // 4
								}
								else
								{
									// hhscore[v] += unit_weight[c] * (min(p->weight, gap1[c] - sat_count[c] - weight) - min(p->weight, gap1[c] - sat_count[c]));//3
									hhscore2[v] += 1 * (min(p->weight, gap1[c] - sat_count2[c] - weight) - min(p->weight, gap1[c] - sat_count2[c])); // 3
								}
							}
						}
						else if (sat_count2[c] + weight > gap1[c]) //(2)-(3)
						{
							gap = clause_true_lit_thres[c] - sat_count2[c];
							for (lit *p = clause_c; (v = p->var_num) != 0; p++)
							{
								if (p->sense == cur_soln[v])
								{
									score2[v] += double(tuned_degree_unit_weight[c] * (max(0, gap + p->weight) - max(0, gap - weight + p->weight)));
									// hhscore[v] += unit_weight[c] * (min(p->weight, -gap) - max(p->weight - weight - sat_count[c] + gap , 0));//5
									hhscore2[v] += 1 * (min(p->weight, -gap) - max(p->weight - weight - sat_count2[c] + gap, 0)); // 5
								}
								else
								{
									// hhscore[v] -= unit_weight[c] * (min(p->weight, gap1[c] - sat_count[c]));//6
									hhscore2[v] -= 1 * (min(p->weight, gap1[c] - sat_count2[c])); // 6
								}
							}
						}
					}
					else if (sat_count2[c] > gap1[c]) //(3) - (3)
					{
						gap = clause_true_lit_thres[c] - sat_count2[c];
						for (lit *p = clause_c; (v = p->var_num) != 0; p++)
						{
							if (p->sense == cur_soln[v])
							{
								score2[v] += double(tuned_degree_unit_weight[c] * (max(0, gap + p->weight) - max(0, gap - weight + p->weight)));
								// hhscore[v] += unit_weight[c] * (max(p->weight - sat_count[c] + gap1[c] , 0) - max(p->weight - sat_count[c] - weight + gap1[c], 0));//7
								hhscore2[v] += 1 * (max(p->weight - sat_count2[c] + gap1[c], 0) - max(p->weight - sat_count2[c] - weight + gap1[c], 0)); // 7
							}
						}
					}
				}
				sat_count2[c] += weight;
			}
			else
			{
				if (sat_count2[c] - weight >= clause_true_lit_thres[c])
				{
					if (sat_count2[c] > gap1[c]) //(3)
					{
						if (sat_count2[c] - weight > gap1[c]) //(3) - (3)
						{
							gap = clause_true_lit_thres[c] - sat_count2[c];
							for (lit *p = clause_c; (v = p->var_num) != 0; p++)
							{
								if (p->sense == cur_soln[v])
								{
									score2[v] -= double(tuned_degree_unit_weight[c] * (max(0, gap + weight + p->weight) - max(0, gap + p->weight)));
									// hhscore[v] += unit_weight[c] * (max(p->weight - sat_count[c] + gap1[c], 0) - max(p->weight - sat_count[c] + weight + gap1[c], 0));//8
									hhscore2[v] += 1 * (max(p->weight - sat_count2[c] + gap1[c], 0) - max(p->weight - sat_count2[c] + weight + gap1[c], 0)); // 8
								}
							}
						}
						else if (sat_count2[c] - weight <= gap1[c]) //(3) - (2)
						{
							gap = clause_true_lit_thres[c] - sat_count2[c];
							for (lit *p = clause_c; (v = p->var_num) != 0; p++)
							{
								if (p->sense == cur_soln[v])
								{
									score2[v] -= double(tuned_degree_unit_weight[c] * (max(0, gap + weight + p->weight) - max(0, gap + p->weight)));
									// hhscore[v] += unit_weight[c] * (max(p->weight - sat_count[c] + gap1[c], 0) - min(p->weight, sat_count[c] - weight - sat_count[c]));//9
									hhscore2[v] += 1 * (max(p->weight - sat_count2[c] + gap1[c], 0) - min(p->weight, sat_count2[c] - weight - sat_count2[c])); // 9
								}
								else
								{
									// hhscore[v] += unit_weight[c] * (min(p->weight, gap1[c] - sat_count[c] + weight));//10
									hhscore2[v] += 1 * (min(p->weight, gap1[c] - sat_count2[c] + weight)); // 10
								}
							}
						}
					}
					else if (sat_count2[c] <= gap1[c]) //(2) - (2)
					{
						gap = clause_true_lit_thres[c] - sat_count2[c];
						for (lit *p = clause_c; (v = p->var_num) != 0; p++)
						{
							if (p->sense == cur_soln[v])
							{
								score2[v] -= double(tuned_degree_unit_weight[c] * (max(0, gap + weight + p->weight) - max(0, gap + p->weight)));
								// hhscore[v] += unit_weight[c] * (min(p->weight, -gap) - min(p->weight, sat_count[c] - weight - clause_true_lit_thres[c]));//11
								hhscore2[v] += 1 * (min(p->weight, -gap) - min(p->weight, sat_count2[c] - weight - clause_true_lit_thres[c])); // 11
							}
							else
							{
								// hhscore[v] += unit_weight[c] * (min(p->weight, gap1[c] - sat_count[c] + weight) - min(p->weight, gap1[c] - sat_count[c]));//12
								hhscore2[v] += 1 * (min(p->weight, gap1[c] - sat_count2[c] + weight) - min(p->weight, gap1[c] - sat_count2[c])); // 12
							}
						}
					}
				}
				else if (sat_count2[c] >= clause_true_lit_thres[c]) //(2) -(1)
				{
					gap = clause_true_lit_thres[c] - sat_count2[c];
					for (lit *p = clause_c; (v = p->var_num) != 0; p++)
					{
						if (p->sense == cur_soln[v])
						{
							score2[v] -= double(tuned_degree_unit_weight[c] * (p->weight - max(0, gap + p->weight)));
							// hhscore[v] += unit_weight[c] * (p->weight, -gap);//13
							hhscore2[v] += 1 * (p->weight, -gap); // 13
						}
						else
						{
							score2[v] += double(tuned_degree_unit_weight[c] * min(p->weight, gap + weight));
							// hhscore[v] += unit_weight[c] * (min(p->weight - gap - weight, 0) - min(p->weight, gap1[c] - sat_count[c]));//14
							hhscore2[v] += 1 * (min(p->weight - gap - weight, 0) - min(p->weight, gap1[c] - sat_count2[c])); // 14
						}
					}
				}
				else //(1) -(1)
				{
					gap = clause_true_lit_thres[c] - sat_count2[c];
					for (lit *p = clause_c; (v = p->var_num) != 0; p++)
					{
						if (p->sense != cur_soln[v])
						{
							score2[v] += double(tuned_degree_unit_weight[c] * (min(p->weight, gap + weight) - min(p->weight, gap)));
							// hhscore[v] += unit_weight[c] * (max(p->weight - gap - weight, 0) - max(p->weight - gap, 0));//15
							hhscore2[v] += 1 * (max(p->weight - gap - weight, 0) - max(p->weight - gap, 0)); // 15
						}
					}
				}
				sat_count2[c] -= weight;
			}
		}
	}

	// update information of flipvar
	cur_soln[flipvar] = 1 - cur_soln[flipvar];
	score2[flipvar] = -org_flipvar_score;
	sscore2[flipvar] = -org_sscore;
	hhscore2[flipvar] = -org_hhscore;
	update_goodvarstack12(flipvar);
}

void Satlike::update_goodvarstack12(int flipvar)
{
	int v;
	goodvar_stack2_num = 0;
	// remove the vars no longer goodvar in goodvar stack
	// add goodvar
	for (int i = 0; i < var_neighbor_count[flipvar]; ++i)
	{
		v = var_neighbor[flipvar][i];
		if (score2[v] + sscore2[v] > 0)
		{
			goodvar_stack2[goodvar_stack2_num] = v;
			goodvar_stack2_num++;
		}
	}
}

void Satlike::update_goodvarstack1(int flipvar)
{
	int v;

	// remove the vars no longer goodvar in goodvar stack
	for (int index = goodvar_stack_fill_pointer - 1; index >= 0; index--)
	{
		v = goodvar_stack[index];
		if (score[v] + sscore[v] <= 0)
		// if (score[v] + sscore[v] + ans * hhscore[v] <= 0)
		{
			int top_v = mypop(goodvar_stack);
			goodvar_stack[index] = top_v;
			already_in_goodvar_stack[top_v] = index;
			already_in_goodvar_stack[v] = -1;
		}
	}

	// add goodvar
	for (int i = 0; i < var_neighbor_count[flipvar]; ++i)
	{
		v = var_neighbor[flipvar][i];
		if (score[v] + sscore[v] > 0)
		// if (score[v] + sscore[v] + ans * hhscore[v] > 0)
		{
			if (already_in_goodvar_stack[v] == -1)
			{
				already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
				mypush(v, goodvar_stack);
			}
		}
	}
}
void Satlike::update_goodvarstack2(int flipvar)
{
	if (score[flipvar] > 0 && already_in_goodvar_stack[flipvar] == -1)
	{
		already_in_goodvar_stack[flipvar] = goodvar_stack_fill_pointer;
		mypush(flipvar, goodvar_stack);
	}
	else if (score[flipvar] <= 0 && already_in_goodvar_stack[flipvar] != -1)
	{
		int index = already_in_goodvar_stack[flipvar];
		int last_v = mypop(goodvar_stack);
		goodvar_stack[index] = last_v;
		already_in_goodvar_stack[last_v] = index;
		already_in_goodvar_stack[flipvar] = -1;
	}
	int i, v;
	for (i = 0; i < var_neighbor_count[flipvar]; ++i)
	{
		v = var_neighbor[flipvar][i];
		if (score[v] > 0)
		{
			if (already_in_goodvar_stack[v] == -1)
			{
				already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
				mypush(v, goodvar_stack);
			}
		}
		else if (already_in_goodvar_stack[v] != -1)
		{
			int index = already_in_goodvar_stack[v];
			int last_v = mypop(goodvar_stack);
			goodvar_stack[index] = last_v;
			already_in_goodvar_stack[last_v] = index;
			already_in_goodvar_stack[v] = -1;
		}
	}
}

void Satlike::flip(int flipvar)
{
	int i, v, c;
	int index;
	lit *clause_c;
	int weight;
	int gap;

	double org_flipvar_score = score[flipvar];
	double org_sscore = sscore[flipvar];
	double org_hhscore = hhscore[flipvar];
	cur_soln[flipvar] = 1 - cur_soln[flipvar];

	// cout<<"filpvar = "<<flipvar<<endl;

	for (i = 0; i < var_lit_count[flipvar]; ++i)
	{
		c = var_lit[flipvar][i].clause_num;
		clause_c = clause_lit[c];
		weight = var_lit[flipvar][i].weight;
		if (org_clause_weight[c] == top_clause_weight)
		{
			if (cur_soln[flipvar] == var_lit[flipvar][i].sense) // SatL1 = SatL + weight
			{
				if (sat_count[c] < clause_true_lit_thres[c] && sat_count[c] + weight < clause_true_lit_thres[c])
				{
					// hard_unsat_weight -= double(avghard[c] * weight);
					hard_unsat_weight -= weight;
				}
				else if (sat_count[c] < clause_true_lit_thres[c] && sat_count[c] + weight >= clause_true_lit_thres[c])
				{
					// hard_unsat_weight -= double(avghard[c] * (clause_true_lit_thres[c] - sat_count[c]));
					hard_unsat_weight -= (clause_true_lit_thres[c] - sat_count[c]);
				}

				if (sat_count[c] + weight <= clause_true_lit_thres[c])
				{
					gap = clause_true_lit_thres[c] - sat_count[c];
					for (lit *p = clause_c; (v = p->var_num) != 0; p++)
					{
						if (p->sense != cur_soln[v])
						{
							score[v] -= double((tuned_degree_unit_weight[c] * (min(gap, p->weight) - min(gap - weight, p->weight))));
							// hhscore[v] += (unit_weight[c] * (max(p->weight - gap + weight, 0) - max(p->weight - gap , 0)));//1
							hhscore[v] += (1 * (max(p->weight - gap + weight, 0) - max(p->weight - gap, 0))); // 1
						}
					}
				}
				else if (sat_count[c] <= clause_true_lit_thres[c]) // sat_count[c]+weight > clause_true_lit_thres[c]
				{
					gap = clause_true_lit_thres[c] - sat_count[c];
					for (lit *p = clause_c; (v = p->var_num) != 0; p++)
					{
						if (p->sense != cur_soln[v])
						{
							score[v] -= double((tuned_degree_unit_weight[c] * min(gap, p->weight)));
							// hhscore[v] += (unit_weight[c] * (min(p->weight, gap1[c] - sat_count[c] - weight) - max(p->weight - gap, 0)));//2 unsat
							hhscore[v] += (1 * (min(p->weight, gap1[c] - sat_count[c] - weight) - max(p->weight - gap, 0))); // 2 unsat
						}
						else
						{
							score[v] += double(tuned_degree_unit_weight[c] * (p->weight - max(0, gap - weight + p->weight)));
							// hhscore[v] -= unit_weight[c] * min(p->weight, weight - gap);//2 sat
							hhscore[v] -= 1 * min(p->weight, weight - gap); // 2 sat
						}
					}
				}
				else // sat_count[c]+weight > clause_true_lit_thres[c], sat_count[c] > clause_true_lit_thres[c]
				{
					if (sat_count[c] <= gap1[c])
					{
						// gap = clause_true_lit_thres[c] - sat_count[c];
						if (sat_count[c] + weight <= gap1[c])
						{
							gap = clause_true_lit_thres[c] - sat_count[c];
							for (lit *p = clause_c; (v = p->var_num) != 0; p++)
							{
								if (p->sense == cur_soln[v])
								{
									score[v] += double(tuned_degree_unit_weight[c] * (max(0, gap + p->weight) - max(0, gap - weight + p->weight)));
									// hhscore[v] += unit_weight[c] * (min(p->weight, -gap) - min(p->weight, weight - gap));//4
									hhscore[v] += 1 * (min(p->weight, -gap) - min(p->weight, weight - gap)); // 4
								}
								else
								{
									// hhscore[v] += unit_weight[c] * (min(p->weight, gap1[c] - sat_count[c] - weight) - min(p->weight, gap1[c] - sat_count[c]));//3
									hhscore[v] += 1 * (min(p->weight, gap1[c] - sat_count[c] - weight) - min(p->weight, gap1[c] - sat_count[c])); // 3
								}
							}
						}
						else if (sat_count[c] + weight > gap1[c]) //(2)-(3)
						{
							gap = clause_true_lit_thres[c] - sat_count[c];
							for (lit *p = clause_c; (v = p->var_num) != 0; p++)
							{
								if (p->sense == cur_soln[v])
								{
									score[v] += double(tuned_degree_unit_weight[c] * (max(0, gap + p->weight) - max(0, gap - weight + p->weight)));
									// hhscore[v] += unit_weight[c] * (min(p->weight, -gap) - max(p->weight - weight - sat_count[c] + gap , 0));//5
									hhscore[v] += 1 * (min(p->weight, -gap) - max(p->weight - weight - sat_count[c] + gap, 0)); // 5
								}
								else
								{
									// hhscore[v] -= unit_weight[c] * (min(p->weight, gap1[c] - sat_count[c]));//6
									hhscore[v] -= 1 * (min(p->weight, gap1[c] - sat_count[c])); // 6
								}
							}
						}
					}
					else if (sat_count[c] > gap1[c]) //(3) - (3)
					{
						gap = clause_true_lit_thres[c] - sat_count[c];
						for (lit *p = clause_c; (v = p->var_num) != 0; p++)
						{
							if (p->sense == cur_soln[v])
							{
								score[v] += double(tuned_degree_unit_weight[c] * (max(0, gap + p->weight) - max(0, gap - weight + p->weight)));
								// hhscore[v] += unit_weight[c] * (max(p->weight - sat_count[c] + gap1[c] , 0) - max(p->weight - sat_count[c] - weight + gap1[c], 0));//7
								hhscore[v] += 1 * (max(p->weight - sat_count[c] + gap1[c], 0) - max(p->weight - sat_count[c] - weight + gap1[c], 0)); // 7
							}
						}
					}
				}
				if (sat_count[c] < clause_true_lit_thres[c] && sat_count[c] + weight >= clause_true_lit_thres[c])
				{
					sat(c);
				}
				sat_count[c] += weight;
			}
			else // cur_soln[flipvar] != cur_lit.sense
			{
				if (sat_count[c] >= clause_true_lit_thres[c] && sat_count[c] - weight < clause_true_lit_thres[c])
				{
					// hard_unsat_weight += double(avghard[c] * (clause_true_lit_thres[c] - sat_count[c] + weight));
					hard_unsat_weight += (clause_true_lit_thres[c] - sat_count[c] + weight);
				}
				else if (sat_count[c] < clause_true_lit_thres[c] && sat_count[c] - weight < clause_true_lit_thres[c])
				{
					// hard_unsat_weight += double(avghard[c] * weight);
					hard_unsat_weight += weight;
				}

				if (sat_count[c] - weight >= clause_true_lit_thres[c])
				{
					if (sat_count[c] > gap1[c]) //(3)
					{
						if (sat_count[c] - weight > gap1[c])
						{
							gap = clause_true_lit_thres[c] - sat_count[c];
							for (lit *p = clause_c; (v = p->var_num) != 0; p++)
							{
								if (p->sense == cur_soln[v])
								{
									score[v] -= double(tuned_degree_unit_weight[c] * (max(0, gap + weight + p->weight) - max(0, gap + p->weight)));
									// hhscore[v] += unit_weight[c] * (max(p->weight - sat_count[c] + gap1[c], 0) - max(p->weight - sat_count[c] + weight + gap1[c], 0));//8
									hhscore[v] += 1 * (max(p->weight - sat_count[c] + gap1[c], 0) - max(p->weight - sat_count[c] + weight + gap1[c], 0)); // 8
								}
							}
						}
						else if (sat_count[c] - weight <= gap1[c]) //(3) - (2)
						{
							gap = clause_true_lit_thres[c] - sat_count[c];
							for (lit *p = clause_c; (v = p->var_num) != 0; p++)
							{
								if (p->sense == cur_soln[v])
								{
									score[v] -= double(tuned_degree_unit_weight[c] * (max(0, gap + weight + p->weight) - max(0, gap + p->weight)));
									// hhscore[v] += unit_weight[c] * (max(p->weight - sat_count[c] + gap1[c], 0) - min(p->weight, sat_count[c] - weight - sat_count[c]));//9
									hhscore[v] += 1 * (max(p->weight - sat_count[c] + gap1[c], 0) - min(p->weight, sat_count[c] - weight - sat_count[c])); // 9
								}
								else
								{
									// hhscore[v] += unit_weight[c] * (min(p->weight, gap1[c] - sat_count[c] + weight));//10
									hhscore[v] += 1 * (min(p->weight, gap1[c] - sat_count[c] + weight)); // 10
								}
							}
						}
					}
					else if (sat_count[c] <= gap1[c]) //(2) - (2)
					{
						gap = clause_true_lit_thres[c] - sat_count[c];
						for (lit *p = clause_c; (v = p->var_num) != 0; p++)
						{
							if (p->sense == cur_soln[v])
							{
								score[v] -= double(tuned_degree_unit_weight[c] * (max(0, gap + weight + p->weight) - max(0, gap + p->weight)));
								// hhscore[v] += unit_weight[c] * (min(p->weight, -gap) - min(p->weight, sat_count[c] - weight - clause_true_lit_thres[c]));//11
								hhscore[v] += 1 * (min(p->weight, -gap) - min(p->weight, sat_count[c] - weight - clause_true_lit_thres[c])); // 11
							}
							else
							{
								// hhscore[v] += unit_weight[c] * (min(p->weight, gap1[c] - sat_count[c] + weight) - min(p->weight, gap1[c] - sat_count[c]));//12
								hhscore[v] += 1 * (min(p->weight, gap1[c] - sat_count[c] + weight) - min(p->weight, gap1[c] - sat_count[c])); // 12
							}
						}
					}
				}
				else if (sat_count[c] >= clause_true_lit_thres[c])
				{
					gap = clause_true_lit_thres[c] - sat_count[c];
					for (lit *p = clause_c; (v = p->var_num) != 0; p++)
					{
						if (p->sense == cur_soln[v])
						{
							score[v] -= double(tuned_degree_unit_weight[c] * (p->weight - max(0, gap + p->weight)));
							// hhscore[v] += unit_weight[c] * (p->weight, -gap);//13
							hhscore[v] += 1 * (p->weight, -gap); // 13
						}
						else
						{
							score[v] += double(tuned_degree_unit_weight[c] * min(p->weight, gap + weight));
							// hhscore[v] += unit_weight[c] * (min(p->weight - gap - weight, 0) - min(p->weight, gap1[c] - sat_count[c]));//14
							hhscore[v] += 1 * (min(p->weight - gap - weight, 0) - min(p->weight, gap1[c] - sat_count[c])); // 14
						}
					}
				}
				else //(1) -(1)
				{
					gap = clause_true_lit_thres[c] - sat_count[c];
					for (lit *p = clause_c; (v = p->var_num) != 0; p++)
					{
						if (p->sense != cur_soln[v])
						{
							score[v] += double(tuned_degree_unit_weight[c] * (min(p->weight, gap + weight) - min(p->weight, gap)));
							// hhscore[v] += unit_weight[c] * (max(p->weight - gap - weight, 0) - max(p->weight - gap, 0));//15
							hhscore[v] += 1 * (max(p->weight - gap - weight, 0) - max(p->weight - gap, 0)); // 15
						}
					}
				}
				if (sat_count[c] >= clause_true_lit_thres[c] && sat_count[c] - weight < clause_true_lit_thres[c])
				{
					unsat(c);
				}
				sat_count[c] -= weight;

			} // end else
		}
		else
		{
			if (cur_soln[flipvar] == var_lit[flipvar][i].sense) // flip better
			{
				soft_unsat_weight -= org_clause_weight[c];
				sat(c);
				sat_count[c] += weight;
			}
			else // flip worse
			{
				soft_unsat_weight += org_clause_weight[c];
				unsat(c);
				sat_count[c] -= weight;
			} // end else
		}
	}

	// update information of flipvar
	score[flipvar] = -org_flipvar_score;
	sscore[flipvar] = -org_sscore;
	hhscore[flipvar] = -org_hhscore; // hhscore
	update_goodvarstack1(flipvar);
}

void Satlike::local_search(vector<int> &init_solution)
{
	settings();
	max_flips = 200000000;
	init(init_solution);
	cout << "time is " << get_runtime() << endl;
	cout << "hard unsat nb is " << hard_unsat_nb << endl;
	cout << "soft unsat nb is " << soft_unsat_weight << endl;
	cout << "goodvar nb is " << goodvar_stack_fill_pointer << endl;
}

void Satlike::print_best_solution(char *filename, char *seed) // 传入文件名和seed
{
	if (best_soln_feasible == 1)
	{

		if (verify_sol())
			// cout << filename << '\t' << seed << '\t' << opt_unsat_weight << '\t' << opt_time << '\t' << tries << '\t' << step << endl;
			// cout << filename << " " << seed << " " << opt_unsat_weight_300 << " " << obj_300 << " " << opt_time_300 << " " << opt_unsat_weight_1800 << " " << obj_1800 << " " << opt_time_1800 << " " << opt_unsat_weight << " " << obj_3600 << " " << opt_time << endl;
			cout << filename << " " << seed << " " << obj_3600 << " " << opt_time << endl;
		else
			cout << "verify solion wrong " << endl;
	}
	else
		cout << filename << '\t' << seed << '\t' << "no feasible solution" << endl;

	ofstream ofile("solution.res");
	ofile << num_vars << " ";
	for (int i = 1; i <= num_vars; i++)
	{
		ofile << best_soln[i] << " ";
	}
	ofile << endl;
}

void Satlike::local_search_with_decimation(vector<int> &init_solution, char *inputfile)
{
	int step_count = 0;
	int cur_hard_unsat_nb = 0;
	int xishu;
	// int MABh_step = 0;

	if (num_vars < 1800000)
		xishu = 1;
	else
		xishu = 16;

	settings();
	for (tries = 1; tries < max_tries; ++tries)
	{
		init(init_solution);
		cur_hard_unsat_nb = hard_unsat_nb;
		for (step = 1; step < max_flips; ++step)
		{
			step_count++;
			// MABh_step++;
			if (hard_unsat_nb < cur_hard_unsat_nb)
			{
				cur_hard_unsat_nb = hard_unsat_nb;
				step_count = step_count / 2 + 1;
			}
			if (hard_unsat_nb == 0 && (soft_unsat_weight < opt_unsat_weight || best_soln_feasible == 0))
			{
				// cout << "unsat soft stack fill pointer" << softunsat_stack_fill_pointer << endl;
				if (soft_unsat_weight < top_clause_weight)
				{
					softclause_weight_threshold += 0;
					best_soln_feasible = 1;
					opt_unsat_weight = soft_unsat_weight;
					step_count = 0;
					opt_time = get_runtime();
					// printf("opt_unsat_weight = %lld, opt_time = %.2f, step = %lld\n", opt_unsat_weight, opt_time, step);
					for (int v = 1; v <= num_vars; ++v)
						best_soln[v] = cur_soln[v];

					if (opt_unsat_weight == 0)
					{
						return;
					}
				}
				cur_hard_unsat_nb = num_hclauses;
				if (num_vars < 1800000 && xishu != 1)
					xishu = xishu / 2;
				else if (num_vars >= 1800000 && xishu != 16)
					xishu = xishu / 2;
			}
			if (step % 1000 == 0)
			{
				double elapse_time = get_runtime();
				if (mark == 0 && elapse_time > 300)
				{
					if (best_soln_feasible == 1)
					{
						opt_unsat_weight_300 = opt_unsat_weight;
						opt_time_300 = opt_time;
						mark = 1;
						for (int v = 1; v <= num_vars; ++v)
						{
							best_soln_300[v] = best_soln[v];
						}
					}
					else
					{
						mark = 1;
						opt_unsat_weight_300 = 9223372036;
						opt_time_300 = -1;
					}
				}
				else if (mark1 == 0 && elapse_time > 1800)
				{
					if (best_soln_feasible == 1)
					{
						opt_unsat_weight_1800 = opt_unsat_weight;
						opt_time_1800 = opt_time;
						mark1 = 1;
						for (int v = 1; v <= num_vars; ++v)
						{
							best_soln_1800[v] = best_soln[v];
						}
					}
					else
					{
						mark1 = 1;
						opt_unsat_weight_1800 = 9223372036;
						opt_time_1800 = -1;
					}
				}
				else if (elapse_time >= cutoff_time)
				{
					return;
				}
			}

			if (goodvar_stack_fill_pointer > 0)
			{
				int flipvar = pick_var();
				flip(flipvar);
				time_stamp[flipvar] = step;
			}
			else
			{
				pick_vars();
			}
			if ((step_count % (MinStep * xishu)) == (MinStep * xishu - 1))
			{
				if (hard_unsat_nb <= MinHard && rand() % 100 < 50)
				{
					turb();
					if (xishu <= delta)
						xishu = 2 * xishu;
				}
				else
				{
					// printf("missing turb\n");
					if (best_soln_feasible == 1)
					{
						for (int v = 1; v <= num_vars; ++v)
							cur_soln[v] = best_soln[v];
						init_turb();
						turb();
					}
				}
				step_count = 0;
				cur_hard_unsat_nb = hard_unsat_nb;
			}
		}
	}
}

void Satlike::turb()
{
	/*select clauses*/
	int c, sel_c, v;
	int cur_hard_unsat_nb;
	int count_var_mark, limit_var_mark, flag = 0;
	int turb_flipvar;
	int step_sum = 50;
	selected_var_num = 0;
	// int fix_value = rand() % 2;

	cur_hard_unsat_nb = hard_unsat_nb; // init
	memset(turb_hardunsat_stack, 0, num_clauses * sizeof(int));
	memset(is_selected_clauses, 0, num_clauses * sizeof(int));
	memset(selected_var, 0, num_vars * sizeof(int));
	memset(var_mark, 0, num_vars * sizeof(int));

	for (c = 0; c < cur_hard_unsat_nb; c++)
	{
		turb_hardunsat_stack[c] = hardunsat_stack[c];
	}
	count_var_mark = 0;
	int theshold = num_vars * gamma + 1;
	if (theshold == 1)
		return;
	if (cur_hard_unsat_nb >= 10)
	{
		for (c = 0; c < 10; c++)
		{

			int index = rand() % cur_hard_unsat_nb;
			sel_c = turb_hardunsat_stack[index];
			while (is_selected_clauses[sel_c] != 0)
			{
				index = (index + 1) % cur_hard_unsat_nb;
				sel_c = turb_hardunsat_stack[index];
			}

			is_selected_clauses[sel_c] = 1;
			for (lit *p = clause_lit[sel_c]; (v = p->var_num) != 0; p++)
			{
				if (var_mark[v] == 0)
				{
					count_var_mark++;

					if (count_var_mark >= theshold)
					{
						flag = 1;
						break;
					}
					var_mark[v] = 1;
					selected_var[selected_var_num++] = v;

					if (hard_unsat_nb <= MaxHard && rand() % 100 < 50 && cur_soln[v] == 1)
					{
						turb_flipvar = v;
						flip(turb_flipvar);
						step_sum++;
					}
				}
			}
			if (flag == 1)
			{
				break;
			}
		}
	}
	else
	{
		int unsat_hard_num = cur_hard_unsat_nb;
		int ss = 0;

		for (c = 0; c < 10; c++)
		{
			if (unsat_hard_num > 0)
			{
				sel_c = turb_hardunsat_stack[c];
				is_selected_clauses[sel_c] = 1;
				unsat_hard_num--;
				ss = 1;
			}

			else
			{
				sel_c = rand() % num_hclauses;
				while (is_selected_clauses[sel_c] != 0)
				{
					sel_c = (sel_c + 1) % num_hclauses;
				}
				is_selected_clauses[sel_c] = 2;
				ss = 2;
			}

			for (lit *p = clause_lit[sel_c]; (v = p->var_num) != 0; p++)
			{
				if (var_mark[v] == 0)
				{
					count_var_mark++;

					if (count_var_mark >= theshold)
					{
						flag = 1;
						break;
					}
					var_mark[v] = 1;
					selected_var[selected_var_num++] = v;

					if (hard_unsat_nb <= MaxHard && ss == 1 && rand() % 100 < 50 && cur_soln[v] == 1)
					{
						turb_flipvar = v;
						flip(turb_flipvar);
						step_sum++;
					}
				}
			}
			if (flag == 1)
			{
				break;
			}
		}
	}

	int turb_step;
	int turb_hardunsat_nb = hard_unsat_nb;
	long long turb_opt_unsat_weight;
	int mark_soln = 0;

	int ttt = 0;
	int turb_best_var;
	for (turb_step = 1; turb_step <= Lopt; turb_step++)
	{
		ttt++;

		if (ttt == 1)
			turb_best_var = turb_pick_var(-1);
		else
			turb_best_var = turb_pick_var(turb_best_var);
		flip(turb_best_var);
		if (hard_unsat_nb < turb_hardunsat_nb)
		{
			turb_hardunsat_nb = hard_unsat_nb;
			turb_opt_unsat_weight = soft_unsat_weight;
			mark_soln = 1;
			turb_step = 1;
		}
		if (hard_unsat_nb == 0 && (soft_unsat_weight < opt_unsat_weight || best_soln_feasible == 0))

		{
			break;
		}
	}

	if (mark_soln == 1)
	{

		;
	}
}

void Satlike::init_turb()
{
	int v, c;
	int j, i;

	// init stacks
	hard_unsat_nb = 0;
	// soft_unsat_nb = 0;//turb
	hardunsat_stack_fill_pointer = 0;
	softunsat_stack_fill_pointer = 0;
	unsatvar_stack_fill_pointer = 0;

	/* figure out sat_count, sat_var, soft_unsat_weight and init unsat_stack */
	soft_unsat_weight = 0;
	hard_unsat_weight = 0;

	for (c = 0; c < num_clauses; ++c)
	{
		sat_count[c] = 0;
		for (j = 0; j < clause_lit_count[c]; ++j)
		{

			if (cur_soln[clause_lit[c][j].var_num] == clause_lit[c][j].sense)
			{
				sat_count[c] += clause_lit[c][j].weight;
				sat_var[c] = clause_lit[c][j].var_num;
			}
		}
		if (sat_count[c] < clause_true_lit_thres[c])
		{
			if (org_clause_weight[c] != top_clause_weight)
				soft_unsat_weight += (clause_true_lit_thres[c] - sat_count[c]) * org_unit_weight[c];
			else
				hard_unsat_weight += clause_true_lit_thres[c] - sat_count[c]; // zyj
			unsat(c);
		}
		// cout<<"soft_unsat_weight "<<soft_unsat_weight<<endl;
	}

	/*figure out score*/
	int sense, weight;

	for (v = 1; v <= num_vars; v++)
	{
		score[v] = 0;
		sscore[v] = 0;
		hhscore[v] = 0;
		for (int i = 0; i < var_lit_count[v]; ++i)
		{
			c = var_lit[v][i].clause_num;
			sense = var_lit[v][i].sense;
			weight = var_lit[v][i].weight;
			if (org_clause_weight[c] == top_clause_weight)
			{
				if (sat_count[c] < clause_true_lit_thres[c])
				{
					if (sense != cur_soln[v])
					{
						score[v] += double(tuned_degree_unit_weight[c] * min(clause_true_lit_thres[c] - sat_count[c], weight));
						// hhscore[v] += unit_weight[c] * max(weight - (clause_true_lit_thres[c] - sat_count[c]), 0);
						hhscore[v] += 1 * max(weight - (clause_true_lit_thres[c] - sat_count[c]), 0);
					}
					else
					{
						score[v] -= double(tuned_degree_unit_weight[c] * weight);
						hhscore[v] += 0;
					}
				}
				else if (sat_count[c] >= clause_true_lit_thres[c])
				{
					if (sat_count[c] <= gap1[c])
					{
						if (sense == cur_soln[v])
						{
							score[v] -= double(tuned_degree_unit_weight[c] * max(0, clause_true_lit_thres[c] - sat_count[c] + weight));
							// hhscore[v] -= unit_weight[c] * min(weight, sat_count[c] - clause_true_lit_thres[c]);
							hhscore[v] -= 1 * min(weight, sat_count[c] - clause_true_lit_thres[c]);
						}
						else
						{
							// hhscore[v] += unit_weight[c] * min(weight, gap1[c] - sat_count[c]);
							hhscore[v] += 1 * min(weight, gap1[c] - sat_count[c]);
						}
					}
					else if (sat_count[c] > gap1[c])
					{
						if (sense == cur_soln[v])
						{
							score[v] -= double(tuned_degree_unit_weight[c] * max(0, clause_true_lit_thres[c] - sat_count[c] + weight));
							// hhscore[v] -= unit_weight[c] * max(weight - (sat_count[c] - gap1[c]), 0);
							hhscore[v] -= 1 * max(weight - (sat_count[c] - gap1[c]), 0);
						}
					}
				}
			}
			else
			{
				if (sat_count[c] < clause_true_lit_thres[c])
				{
					if (sense != cur_soln[v])
					{
						sscore[v] += unit_weight[c] * tune_soft_clause_weight[c];
					}
					else
						sscore[v] -= unit_weight[c] * tune_soft_clause_weight[c];
				}
				else if (sat_count[c] >= clause_true_lit_thres[c])
				{
					if (sense == cur_soln[v])
					{
						sscore[v] -= unit_weight[c] * tune_soft_clause_weight[c];
					}
				}
			}
		}
	}

	// init goodvars stack
	goodvar_stack_fill_pointer = 0;

	for (v = 1; v <= num_vars; v++)
	{
		// if (score[v] + sscore[v] + ans * hhscore[v] > 0)//加入hhscore
		if (score[v] + sscore[v] > 0)
		{
			already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
			mypush(v, goodvar_stack);
		}
		else
			already_in_goodvar_stack[v] = -1;
	}
}

int Satlike::turb_pick_var(int last_flip_var)
{
	int turb_best_var;
	int i, v;
	long long pick_score[300] = {0};
	int pos_score_var[300] = {0};
	// long long *pick_score = new long long[30]();
	// int *pos_score_var = new int[30]();
	int bms_turb_pick;
	long long abs_min_pick_score, abs_min;
	int turb_size = 1;
	int pos_size = 1;
	long long min_pick_score = 200000000000000;
	int turb_seed;
	long long sum_pick_score;
	long long ssscore;

	bms_turb_pick = selected_var_num / 2 + 1;

	turb_best_var = selected_var[rand() % selected_var_num];

	for (i = 1; i < bms_turb_pick; i++)
	{
		v = selected_var[rand() % selected_var_num];

		// v = selected_var[i];
		if (last_flip_var == v)
			continue;
		// v = i;
		if (score[v] + sscore[v] > score[turb_best_var] + sscore[turb_best_var])
			turb_best_var = v;
		else if (score[v] + sscore[v] == score[turb_best_var] + sscore[turb_best_var])
		{
			if (hhscore[v] > hhscore[turb_best_var])
			{
				turb_best_var = v;
			}
			else if (hhscore[v] == hhscore[turb_best_var])
			{
				int rand_select = rand() % 2;
				if (rand_select == 1)
				{
					turb_best_var = v;
				}
			}
		}
	}
	return turb_best_var;

	for (i = 0; i < bms_turb_pick; i++)
	{
		v = selected_var[rand() % selected_var_num];
		ssscore = score[v] + sscore[v];
		pick_score[turb_size++] = ssscore;
		pos_score_var[pos_size++] = v;
		if (ssscore < min_pick_score)
		{
			min_pick_score = ssscore;
		}
	}

	if (min_pick_score <= 0)
	{
		abs_min = abs(min_pick_score);
		abs_min_pick_score = abs_min + 1;
		// pick_score[0] +=  abs_min_pick_score;
		for (i = 1; i <= bms_turb_pick; i++)
		{
			pick_score[i] = pick_score[i] + abs_min_pick_score + pick_score[i - 1];
		}
	}
	else
	{
		for (i = 1; i <= bms_turb_pick; i++)
		{
			pick_score[i] = pick_score[i] + pick_score[i - 1];
		}
	}
	sum_pick_score = pick_score[bms_turb_pick];

	turb_seed = rand() % sum_pick_score;
	for (i = 1; i <= bms_turb_pick; i++)
	{
		if (turb_seed < pick_score[i])
		{
			turb_best_var = pos_score_var[i];
			break;
		}
	}
	return turb_best_var;
}

void Satlike::get_obj(string opb_file_name)
{
	ifstream opb_ifile(opb_file_name.c_str());
	int num_hardclauses, i;
	istringstream ss;
	istringstream ss2;

	string line, coeff, var, other;
	getline(opb_ifile, line);

	ss.clear();
	ss.str(line);
	ss.seekg(0, ios::beg);
	ss >> other;
	ss >> other;
	ss >> num_vars_opb;
	ss >> other;
	ss >> num_hardclauses;

	int v;

	for (v = 1; v <= num_vars_opb; v++)
	{
		opb[v] = 0;
	}

	getline(opb_ifile, line);
	while (line[0] == '*')
	{
		getline(opb_ifile, line);
	}

	ss.clear();
	ss.str(line);
	ss.seekg(0, ios::beg);
	ss >> coeff;
	ss >> coeff >> var;
	int icoeff, ivar;
	while (coeff != ";")
	{
		ss2.clear();
		ss2.str(coeff);
		ss2.seekg(0, ios::beg);
		ss2 >> icoeff;

		ss2.clear();
		ss2.str(var.substr(1));
		ss2.seekg(0, ios::beg);
		ss2 >> ivar; // 1

		opb[ivar] = icoeff;

		ss >> coeff >> var;
	}
	if (num_vars == num_vars_opb)
	{
		for (i = 1; i <= num_vars; i++)
		{
			obj_300 += best_soln_300[i] * opb[i];
			obj_1800 += best_soln_1800[i] * opb[i];
			obj_3600 += best_soln[i] * opb[i];
		}
	}
	else
	{
		cout << " error : num_vars != num_vars_opb "
			 << " " << opb_file_name << endl;
		cout << "num_vars = " << num_vars << "num_vars_opb = " << num_vars_opb << endl;
	}
}

void Satlike::check_softunsat_weight()
{
	int c, j, flag;
	long long verify_unsat_weight = 0;

	for (c = 0; c < num_clauses; ++c)
	{
		flag = 0;
		int tem_clause_true_lit_count = 0;
		for (j = 0; j < clause_lit_count[c]; ++j)
		{
			if (cur_soln[clause_lit[c][j].var_num] == clause_lit[c][j].sense)
			{
				tem_clause_true_lit_count++;
			}
		}
		if (tem_clause_true_lit_count >= clause_true_lit_thres[c])
			flag = 1;

		if (flag == 0)
		{
			if (org_clause_weight[c] == top_clause_weight) // verify hard clauses
			{
				continue;
			}
			else
			{
				verify_unsat_weight += org_unit_weight[c] * (clause_true_lit_thres[c] - tem_clause_true_lit_count);
			}
		}
	}

	if (verify_unsat_weight != soft_unsat_weight)
	{
		cout << step << endl;
		cout << "verify unsat weight is" << verify_unsat_weight << " and soft unsat weight is " << soft_unsat_weight << endl;
	}
	// return 0;
}

bool Satlike::verify_sol()
{
	int c, j, flag;
	long long verify_unsat_weight = 0;

	for (c = 0; c < num_clauses; ++c)
	{
		flag = 0;
		int tem_clause_true_lit_count = 0;
		for (j = 0; j < clause_lit_count[c]; ++j)
		{
			if (best_soln[clause_lit[c][j].var_num] == clause_lit[c][j].sense)
			{
				tem_clause_true_lit_count += clause_lit[c][j].weight;
			}
		}
		if (tem_clause_true_lit_count >= clause_true_lit_thres[c])
			flag = 1;

		if (flag == 0)
		{
			if (org_clause_weight[c] == top_clause_weight) // verify hard clauses
			{
				// output the clause unsatisfied by the solution
				cout << "c Error: hard clause " << c << " is not satisfied" << endl;

				cout << "c ";
				for (j = 0; j < clause_lit_count[c]; ++j)
				{
					if (clause_lit[c][j].sense == 0)
						cout << "-";
					cout << clause_lit[c][j].var_num << " ";
				}
				cout << endl;
				cout << "c ";
				for (j = 0; j < clause_lit_count[c]; ++j)
					cout << best_soln[clause_lit[c][j].var_num] << " ";
				cout << endl;
				return 0;
			}
			else
			{
				verify_unsat_weight += org_unit_weight[c] * (clause_true_lit_thres[c] - tem_clause_true_lit_count);
				/*
				cout << "c wanning: soft clause " << c << " is not satisfied" << endl;

				cout << "c org clause weight " << org_clause_weight[c] << " " << endl;

				for (j = 0; j < clause_lit_count[c]; ++j)
				{
					if (clause_lit[c][j].sense == 0)
						cout << "-";
					cout << clause_lit[c][j].var_num << " ";
				}
				cout << endl;
				//cout << "c ";
				for (j = 0; j < clause_lit_count[c]; ++j)
					cout << best_soln[clause_lit[c][j].var_num] << " ";
				cout << endl;*/
				// return 0;
			}
		}
	}

	if (verify_unsat_weight == opt_unsat_weight)
		return 1;
	else
	{
		cout << "c Error: find opt=" << opt_unsat_weight << ", but verified opt=" << verify_unsat_weight << endl;
	}
	return 0;
}

void Satlike::simple_print()
{
	if (best_soln_feasible == 1)
	{
		if (verify_sol() == 1)
			cout << opt_unsat_weight << '\t' << opt_time << endl;
		else
			cout << "solution is wrong " << endl;
	}
	else
		cout << -1 << '\t' << -1 << endl;
}

void Satlike::increase_weights()
{
	int i, c, v;
	int weight;
	int flag = 0;

	for (i = 0; i < hardunsat_stack_fill_pointer; ++i)
	{
		c = hardunsat_stack[i];
		if (clause_visied_times[c] < clause_true_lit_thres[c] / clause_weight[c])
		{
			clause_visied_times[c]++;
			continue;
		}
		else
		{
			clause_visied_times[c] = 0;
		}

		inc_hard_weight += clause_weight[c];
		// clause_weight[c] += (h_inc * clause_true_lit_thres[c]);
		// cout << "c: " << c << endl;
		unit_weight[c] += h_inc;
		tuned_degree_unit_weight[c] = double(unit_weight[c]) / avg_clause_coe[c]; // Nupbo

		for (lit *p = clause_lit[c]; (v = p->var_num) != 0; p++)
		{
			weight = p->weight;
			if (p->sense != cur_soln[v])
			{
				score[v] += double(h_inc * min(clause_true_lit_thres[c] - sat_count[c], weight)) / avg_clause_coe[c];
				if (score[v] + sscore[v] > 0 && already_in_goodvar_stack[v] == -1)
				// if (score[v] + sscore[v] + ans * hhscore[v] > 0 && already_in_goodvar_stack[v] == -1)
				{
					already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
					mypush(v, goodvar_stack);
				}
			}
			else
			{
				score[v] -= double(h_inc * weight) / avg_clause_coe[c];
				if (already_in_goodvar_stack[v] != -1 && score[v] + sscore[v] <= 0)
				// if (already_in_goodvar_stack[v] != -1 && score[v] + sscore[v] + ans * hhscore[v] <= 0)
				{
					int top_v = mypop(goodvar_stack);
					goodvar_stack[already_in_goodvar_stack[v]] = top_v;
					already_in_goodvar_stack[top_v] = already_in_goodvar_stack[v];
					already_in_goodvar_stack[v] = -1;
				}
			}
		}
	}

	if (0 == hard_unsat_nb)
	{
		// flag = 1;
		ave_soft_weight += total_soft_weight / num_sclauses;
		inc_hard_weight += total_soft_weight / num_sclauses;
		for (c = 0; c < num_clauses; ++c)
		{
			if (org_clause_weight[c] == top_clause_weight)
				continue;

			// unit_weight[c] += org_unit_weight[c];
			unit_weight[c] += 1;

			if (sat_count[c] < clause_true_lit_thres[c])
			{
				for (lit *p = clause_lit[c]; (v = p->var_num) != 0; p++)
				{
					sscore[v] += tune_soft_clause_weight[c];
					// min(clause_true_lit_thres[c] - sat_count[c], weight);
					if (score[v] + sscore[v] > 0 && already_in_goodvar_stack[v] == -1)
					// if (score[v] + sscore[v] + ans * hhscore[v]> 0 && already_in_goodvar_stack[v] == -1)
					{
						already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
						mypush(v, goodvar_stack);
					}
				}
			}
			else if (sat_count[c] == 1)
			{
				for (lit *p = clause_lit[c]; (v = p->var_num) != 0; p++)
				{
					if (p->sense == cur_soln[v])
					{
						sscore[v] -= tune_soft_clause_weight[c];
						if (already_in_goodvar_stack[v] != -1 && score[v] + sscore[v] <= 0)
						{
							int top_v = mypop(goodvar_stack);
							goodvar_stack[already_in_goodvar_stack[v]] = top_v;
							already_in_goodvar_stack[top_v] = already_in_goodvar_stack[v];
							already_in_goodvar_stack[v] = -1;
						}
					}
				}
			}
		}
	}

	ave_hard_weight += (inc_hard_weight / num_hclauses);
	inc_hard_weight %= num_hclauses;
}

void Satlike::smooth_weights()
{
	int i, clause, v;
	int weight;
	if (soft_unsat_weight < opt_unsat_weight && ave_soft_weight > (total_soft_weight / num_sclauses))
	{
		ave_soft_weight -= (total_soft_weight / num_sclauses);
		inc_hard_weight -= (total_soft_weight / num_sclauses);
		if (inc_hard_weight < 0)
			inc_hard_weight = 0;
		for (int c = 0; c < num_clauses; ++c)
		{
			if (org_clause_weight[c] == top_clause_weight)
			{
				if (unit_weight[c] == 1 && sat_count[c] < clause_true_lit_thres[c])
					continue;

				unit_weight[c]--;
				for (lit *p = clause_lit[c]; (v = p->var_num) != 0; p++)
				{
					weight = p->weight;
					if (p->sense == cur_soln[v])
					{
						if (sat_count[c] - weight < clause_true_lit_thres[c])
						{
							score[v] += clause_true_lit_thres[c] - sat_count[c] + weight;
							if (score[v] + sscore[v] > 0 && already_in_goodvar_stack[v] == -1)
							// if (score[v] + sscore[v] + ans * hhscore[v] > 0 && already_in_goodvar_stack[v] == -1)
							{
								already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
								mypush(v, goodvar_stack);
							}
						}
					}
				}
			}
			else
			{
				unit_weight[c]--;
				for (lit *p = clause_lit[c]; (v = p->var_num) != 0; p++)
				{
					weight = p->weight;
					if (p->sense == cur_soln[v])
					{
						if (sat_count[c] - weight < clause_true_lit_thres[c])
						{
							sscore[v] += clause_true_lit_thres[c] - sat_count[c] + weight;
							if (score[v] + sscore[v] > 0 && already_in_goodvar_stack[v] == -1)
							// if (score[v] + sscore[v] + ans * hhscore[v] > 0 && already_in_goodvar_stack[v] == -1)
							{
								already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
								mypush(v, goodvar_stack);
							}
						}
					}
				}
			}
		}
	}
	else
	{
		for (int c = 0; c < num_clauses; ++c)
		{
			if (org_clause_weight[c] == top_clause_weight)
			{
				if (unit_weight[c] == 1 && sat_count[c] < clause_true_lit_thres[c])
					continue;

				unit_weight[c]--;
				for (lit *p = clause_lit[c]; (v = p->var_num) != 0; p++)
				{
					weight = p->weight;
					if (p->sense == cur_soln[v])
					{
						if (sat_count[c] - weight < clause_true_lit_thres[c])
						{
							score[v] += clause_true_lit_thres[c] - sat_count[c] + weight;
							if (score[v] + sscore[v] > 0 && already_in_goodvar_stack[v] == -1)
							// if (score[v] + sscore[v] + ans * hhscore[v] > 0 && already_in_goodvar_stack[v] == -1)
							{
								already_in_goodvar_stack[v] = goodvar_stack_fill_pointer;
								mypush(v, goodvar_stack);
							}
						}
					}
				}
			}
		}
	}
}

void Satlike::update_clause_weights()
{
	/*
	if (((rand() % MY_RAND_MAX_INT) * BASIC_SCALE) < smooth_probability)
	{
		smooth_weights();
	}
	else
	{*/
	increase_weights();
	//}
}

inline void Satlike::unsat(int clause)
{
	if (org_clause_weight[clause] == top_clause_weight)
	{
		index_in_hardunsat_stack[clause] = hardunsat_stack_fill_pointer;
		mypush(clause, hardunsat_stack);
		hard_unsat_nb++;
	}
	else
	{
		index_in_softunsat_stack[clause] = softunsat_stack_fill_pointer;
		mypush(clause, softunsat_stack);
		// soft_unsat_weight += org_clause_weight[clause];
	}
}

inline void Satlike::sat(int clause)
{
	int index, last_unsat_clause;

	if (org_clause_weight[clause] == top_clause_weight)
	{

		last_unsat_clause = mypop(hardunsat_stack);
		index = index_in_hardunsat_stack[clause];
		hardunsat_stack[index] = last_unsat_clause;
		index_in_hardunsat_stack[last_unsat_clause] = index;

		hard_unsat_nb--;
	}
	else
	{
		last_unsat_clause = mypop(softunsat_stack);
		index = index_in_softunsat_stack[clause];
		softunsat_stack[index] = last_unsat_clause;
		index_in_softunsat_stack[last_unsat_clause] = index;

		// soft_unsat_weight -= org_clause_weight[clause];
	}
}

void Satlike::check_new_score()
{
	long long tem_score = 0;
	long long tem_sscore = 0;
	int sense, c, v, i;
	int weight;
	for (v = 1; v <= num_vars; v++)
	{
		tem_score = 0;
		tem_sscore = 0;
		for (i = 0; i < var_lit_count[v]; ++i)
		{
			c = var_lit[v][i].clause_num;
			sense = var_lit[v][i].sense;
			weight = var_lit[v][i].weight;
			if (org_clause_weight[c] == top_clause_weight)
			{
				if (sat_count[c] < clause_true_lit_thres[c])
				{
					if (sense != cur_soln[v])
					{
						tem_score += unit_weight[c] * min(clause_true_lit_thres[c] - sat_count[c], weight);
					}
					else
						tem_score -= unit_weight[c] * weight;
				}
				else if (sat_count[c] >= clause_true_lit_thres[c])
				{
					if (sense == cur_soln[v])
					{
						tem_score -= unit_weight[c] * max(0, clause_true_lit_thres[c] - sat_count[c] + weight);
					}
				}
			}
			else
			{
				if (sat_count[c] < clause_true_lit_thres[c])
				{
					if (sense != cur_soln[v])
					{
						tem_sscore += unit_weight[c] * min(clause_true_lit_thres[c] - sat_count[c], weight);
					}
					else
						tem_sscore -= unit_weight[c] * weight;
				}
				else if (sat_count[c] >= clause_true_lit_thres[c])
				{
					if (sense == cur_soln[v])
					{
						tem_sscore -= unit_weight[c] * max(0, clause_true_lit_thres[c] - sat_count[c] + weight);
					}
				}
			}
		}
		if (tem_score != score[v] || tem_sscore != sscore[v])
		{

			cout << "score is worng in variable " << v << endl;
			cout << "tem_score is " << tem_score << endl;
			cout << "score function is " << score[v] << endl;
			cout << "flip num is " << step << endl;

			for (i = 0; i < var_lit_count[v]; ++i)
			{
				c = var_lit[v][i].clause_num;
				sense = var_lit[v][i].sense;
				weight = var_lit[v][i].weight;
				cout << c << " ";
			}
			cout << endl;
			exit(0);
			break;
		}
	}

	int tem_unsat_softweight = 0;
	for (int i = 0; i < num_clauses; ++i)
	{
		if (org_clause_weight[i] == top_clause_weight)
			continue;
		if (sat_count[i] < clause_true_lit_thres[i])
		{
			tem_unsat_softweight += (clause_true_lit_thres[i] - sat_count[i]);
		}
	}
	if (tem_unsat_softweight != soft_unsat_weight)
	{
		cout << "verify softunsat weight wrong " << endl;
		exit(0);
	}
}

#endif
