--
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



-- all alleles for each unique hit + average difference from their corresponding WT control (averaged per plate), if mutation occurs multiple times- averaged
drop procedure if exists hc_microscopy_data_v2.p_hit_alleles_percentage;
delimiter //
create procedure p_hit_alleles_percentage() 
begin
	with
	cte_unique_hits_systematic_name as
	(
	select
		hit_systematic_name
	from
		unique_hits
	),
	cte_all_hit_mutations as
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
	cte_selected_experiments as
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
	cte_control_data as
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
	select
		a.mutated_gene_standard_name,
		a.mutation,
		avg(case when a.stage= 'formation' then a.percentage_control_minus_mutant else null end) as formation,
		avg(case when a.stage= 'relocation & fusion' then a.percentage_control_minus_mutant else null end) as relocation_and_fusion,
		avg(case when a.stage= 'clearance' then a.percentage_control_minus_mutant else null end) as clearance
	from
	(
	select
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
