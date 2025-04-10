-- list of unique hits (standard gene name)
-- columns: effect-stage
-- 1: identified as a hit for particular effect-stage
-- 0: not considered a hit  for particular effect-stage
drop procedure if exists hc_microscopy_data_v2.p_hits_stages;
delimiter //
create procedure p_hits_stages()
begin
	select
		uh.hit_standard_name,
		-- esl.effect_stage_label,
		max(case when esl.effect_stage_label= 'decreased formation\r' then 1 else 0 end) as decreased_formation,
		max(case when esl.effect_stage_label= 'disrupted relocation & fusion\r' then 1 else 0 end) as disrupted_relocation_and_fusion,
		max(case when esl.effect_stage_label= 'slower clearance\r' then 1 else 0 end) as slower_clearance
	from
		hits_clusters as hc
	inner join 
		unique_hits as uh
	on
		hc.hit_systematic_name=uh.hit_systematic_name
	inner join
		effect_stage_labels as esl
	on
		hc.effect_stage_label_id=esl.effect_stage_label_id
	group by
		uh.hit_standard_name
	order by
		uh.hit_standard_name asc;
end //
delimiter ;
-- call p_hits_stages;



-- all alleles/mutation for each unique hit + average difference from their corresponding WT control 
-- wt controls averaged per plate 
-- mutations occuring multiple times averaged
drop procedure if exists hc_microscopy_data_v2.p_hit_alleles_percentage;
delimiter //
create procedure p_hit_alleles_percentage() 
begin
	with
	cte_unique_hits_systematic_name as -- unique hits
	(
	select
		hit_systematic_name
	from
		unique_hits
	),
	cte_all_hit_mutations as  -- all mutations corresponding to the uniqze hits
	(
	select distinct
		mutation
	from
		strains_and_conditions_main
	where
		mutated_gene_systematic_name in (select * from cte_unique_hits_systematic_name)
	order by
		mutation asc
	),
	cte_selected_experiments as -- list of relevant experiments (containing at least one of the selected mutations)
	(
	select distinct
		e.date_label
	from
		strains_and_conditions_main as sacm
	inner join 
		experiments as e
	on
		sacm.date_label=e.date_label
	inner join
		experiment_types as et
	on
		e.experiment_type_id=et.experiment_type_id
	where
		mutated_gene_systematic_name in (select * from cte_unique_hits_systematic_name) and
		e.data_quality= 'Good' and
		et.experiment_type= 'TS collection screening' and
		et.experiment_subtype= 'first round'
	),
	cte_control_data as -- averaged control data from relevant experiments (averaged per plate and timepoint), single control dataset for each experiment
	(
	select
		cac.date_label,
		cac.timepoint,
		avg(round(cac.number_of_cells_with_foci/cac.number_of_cells*100,2)) as percentage_control_data
	from
		experimental_data_sbw_cell_area_and_counts as cac
	inner join 
		strains_and_conditions_main as sacm
	on
		cac.date_label=sacm.date_label and
		cac.experimental_well_label= sacm.experimental_well_label
	where
		cac.date_label in (select * from cte_selected_experiments) and
		cac.number_of_cells >= 50 and
		sacm.mutation= 'wt control'
	group by
		cac.date_label,
		cac.timepoint
	),
	cte_microscopy_interval as -- microscopy interval for TS screen (1st round)
	(
	select distinct
		microscopy_interval_min
	from
		experiments
	where
		date_label in (select * from cte_selected_experiments)
	),
	cte_microscopy_initial_delay as -- microscopy initial delay for TS screen (1st round)
	(
	select distinct
		microscopy_initial_delay_min
	from
		experiments
	where
		date_label in (select * from cte_selected_experiments)
	)
	select
		*
	from
	(
	select -- pivot: mutation(rows)-stage(columns)-average difference WT-mutant(values)
		a.mutated_gene_standard_name,
		a.mutation,
		avg(case when a.stage= 'formation' then a.percentage_control_minus_mutant else null end) as formation,
		avg(case when a.stage= 'relocation & fusion' then a.percentage_control_minus_mutant else null end) as relocation_and_fusion,
		avg(case when a.stage= 'clearance' then a.percentage_control_minus_mutant else null end) as clearance
	from
	(
	select -- substract (WT-mutant), define stages
		sacm.mutated_gene_standard_name,
		sacm.mutation,
		round(ctrl.percentage_control_data, 2) - round(cac.number_of_cells_with_foci/cac.number_of_cells*100, 2) as percentage_control_minus_mutant,
		case
			when 
				cac.timepoint * (select * from cte_microscopy_interval) - ((select * from cte_microscopy_interval) - (select * from cte_microscopy_initial_delay)) <= 62 then 'formation'
			when 
				cac.timepoint * (select * from cte_microscopy_interval) - ((select * from cte_microscopy_interval) - (select * from cte_microscopy_initial_delay)) > 62 and 
				cac.timepoint * (select * from cte_microscopy_interval) - ((select * from cte_microscopy_interval) - (select * from cte_microscopy_initial_delay)) <= 300 then 'relocation & fusion'
			else 
				'clearance'
		end as stage
	from
		experimental_data_sbw_cell_area_and_counts as cac
	inner join 
		strains_and_conditions_main as sacm
	on
		cac.date_label=sacm.date_label and
		cac.experimental_well_label= sacm.experimental_well_label
	inner join
		cte_control_data as ctrl
	on
		cac.date_label= ctrl.date_label and
		cac.timepoint= ctrl.timepoint
	where
		sacm.mutation in (select * from cte_all_hit_mutations) and
		cac.number_of_cells >= 50
	order by 
		sacm.mutated_gene_standard_name asc,
		sacm.mutation asc
	) as a
	group by
		a.mutated_gene_standard_name,
		a.mutation
	) as b
	where -- remove null cells, mutants cell count of which goes below the threshold at one or more stages
		b.formation is not null and
		b.relocation_and_fusion is not null and
		b.clearance is not null;
end //
delimiter ;
-- call p_hit_alleles_percentage;



drop procedure if exists hc_microscopy_data_v2.p_most_affected_alleles_slower_clearance_hits;
delimiter //
-- for each hit from the 'slower clearance' group returns 3 of the most affected alleles (mutations)
-- calculated as an average difference (between a particular mutant and a corresponding control) of a percentage of cells containing aggregates in the clearance stage (time > 300 minutes)
-- control data averaged per each plate and assigned to each mutant based on date label
create procedure hc_microscopy_data_v2.p_most_affected_alleles_slower_clearance_hits()
#most affected alleles (top3) of the decreased-clearance hits, for follow-up
begin
	with
	cte_slower_clearance_hits as -- systematic names for all slower-clearance hits
	(
	select
		hc.hit_systematic_name
	from
		hits_clusters as hc
	inner join
		effect_stage_labels as esl
	on
		hc.effect_stage_label_id=esl.effect_stage_label_id
	where
		esl.effect_stage_label= 'slower clearance\r'
	),
	cte_selected_experiments as -- selected experiments (date_labels) that contain at least one of the slower-clearance hits (only TS screening 1st round, 'Good' data quality)
	(
	select distinct
		e.date_label
	from
		experiments as e
	inner join
		experiment_types as et
	on
		e.experiment_type_id=et.experiment_type_id
	inner join
		strains_and_conditions_main as sacm
	on
		sacm.date_label=e.date_label
	where
		et.experiment_type = 'TS collection screening' and
		et.experiment_subtype= 'first round' and
		e.data_quality= 'Good' and
		sacm.mutated_gene_systematic_name in (select * from cte_slower_clearance_hits)
	),
	cte_microscopy_interval as -- microscopy interval for TS screen (1st round)
	(
	select distinct
		microscopy_interval_min
	from
		experiments
	where
		date_label in (select * from cte_selected_experiments)
	),
	cte_microscopy_initial_delay as -- microscopy initial delay for TS screen (1st round)
	(
	select distinct
		microscopy_initial_delay_min
	from
		experiments
	where
		date_label in (select * from cte_selected_experiments)
	),
	cte_control_data as -- date label plus data on percentage of cells containing aggregates, only WT controls, averaged per plate (only clearance timepoints included)
	(
	select
		cac.date_label,
		cac.timepoint,
		cac.timepoint * (select * from cte_microscopy_interval) - ((select * from cte_microscopy_interval) - (select * from cte_microscopy_initial_delay)) as timepoint_minutes,
		round(avg(cac.number_of_cells_with_foci/cac.number_of_cells*100), 4) as percentage_control
	from
		experimental_data_sbw_cell_area_and_counts as cac
	inner join
		strains_and_conditions_main as sacm
	on
		sacm.date_label= cac.date_label and
		sacm.experimental_well_label= cac.experimental_well_label
	where
		sacm.mutation= 'wt control' and
		cac.date_label in (select * from cte_selected_experiments) and
		cac.timepoint * (select * from cte_microscopy_interval) - ((select * from cte_microscopy_interval) - (select * from cte_microscopy_initial_delay)) > 300
	group by
		cac.date_label,
		cac.timepoint
	order by 
		cac.date_label asc
	)
	select -- pivot table
		a.mutated_gene_standard_name,
		ifnull(max(case when a.allele_index = 1 then a.mutation else null end), "-") as most_affected_allele,
		ifnull(max(case when a.allele_index = 2 then a.mutation else null end), "-") as second_most_affected_allele,
		ifnull(max(case when a.allele_index = 3 then a.mutation else null end), "-") as third_most_affected_allele
	from
	(
	select -- control minus mutant data
		sacm.mutated_gene_standard_name,
		sacm.mutation,
		row_number() over (partition by sacm.mutated_gene_standard_name order by round(avg(ctrl.percentage_control - cac.number_of_cells_with_foci/cac.number_of_cells*100),2) asc) as allele_index,
		round(avg(ctrl.percentage_control - cac.number_of_cells_with_foci/cac.number_of_cells*100),2) as control_minus_hit_percentage
	from
		experimental_data_sbw_cell_area_and_counts as cac
	inner join
		strains_and_conditions_main as sacm
	on
		sacm.date_label= cac.date_label and
		sacm.experimental_well_label= cac.experimental_well_label
	inner join
		cte_control_data as ctrl
	on
		cac.date_label = ctrl.date_label and 
		cac.timepoint=ctrl.timepoint
	where
		sacm.mutated_gene_systematic_name in (select * from cte_slower_clearance_hits) and
		cac.date_label in (select * from cte_selected_experiments) and
		cac.timepoint * (select * from cte_microscopy_interval) - ((select * from cte_microscopy_interval) - (select * from cte_microscopy_initial_delay)) > 300 -- only clearance-stage timepoints
	group by
		sacm.mutated_gene_standard_name,
		sacm.mutation -- averaged per allele (mutation), in case multiple similar alleles present (duplicate entries present in the TS collection in high frequency)
	order by
		sacm.mutated_gene_standard_name asc
	) as a
	where
		a.control_minus_hit_percentage < 0 -- only negative values (mutants where the clearance is negatively affected, percentage of cells with aggregates higher in the clearance stage, compared to control)
	group by
		a.mutated_gene_standard_name;
end //
delimiter ;
-- call p_most_affected_alleles_slower_clearance_hits();
