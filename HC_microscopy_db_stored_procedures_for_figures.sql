use hc_microscopy_data_v2;

-- FIGURE 1, S1 and S2: WT characterisation --
# stored procedure p_wt_characterisation_data called in python scripts 'Figure1_S1_wt_characterisation' and 'FigureS2_As_pre-treatment'
# returns data from WT characterisation experiment: well labels, timepoint fields, number of cells, number of cells with foci, percentage of cells with foci, average number of foci per cell and average size of a single focus
# data on both As-exposed and control cells
# arguments: 'p_initial_timepoints_skipped': int, number of initial timepoints skipped (generally low quality data from initital timepoint, use at least 1)
#			 'p_experiment_subtype': varchar, 'basic' for Figure 1/S1 and 'pretreatment' for Figure 2
drop procedure if exists p_wt_characterisation_data;
delimiter //
create procedure p_wt_characterisation_data(in p_initial_timepoints_skipped int, in p_experiment_subtype varchar(15))
begin
	with 
	cte_wt_characterisation_date_label as -- relevant experiments
	(
	select distinct
		date_label
	from
		experiments as e
	inner join
		experiment_types as et
	on
		e.experiment_type_id=et.experiment_type_id
	where
		et.experiment_type= 'WT characterisation' and
		et.experiment_subtype= p_experiment_subtype
	),
	cte_microscopy_interval as -- microscopy interval for WT characterisation experiments
	(
	select distinct
		microscopy_interval_min
	from
		experiments
	where
		date_label = (select * from cte_wt_characterisation_date_label)
	),
	cte_microscopy_initial_delay as -- microscopy initial delay for for WT characterisation experiments
	(
	select distinct
		microscopy_initial_delay_min
	from
		experiments
	where
		date_label = (select * from cte_wt_characterisation_date_label)
	)
	select -- fields
		cac.experimental_well_label,
		cac.timepoint,
        (cac.timepoint * (select * from cte_microscopy_interval) - ((select * from cte_microscopy_interval) - (select * from cte_microscopy_initial_delay)))/60 as timepoint_hours,
		cac.timepoint * (select * from cte_microscopy_interval) - ((select * from cte_microscopy_interval) - (select * from cte_microscopy_initial_delay)) as timepoint_minutes,
		cac.number_of_cells,
		cac.number_of_cells_with_foci,
		cac.number_of_cells_with_foci/cac.number_of_cells*100 as percentage_of_cells_with_foci,
		nas.avg_number_of_foci_per_cell,
		nas.avg_size_single_focus
	from
		experimental_data_sbw_cell_area_and_counts as cac
	inner join -- join
		experimental_data_sbw_foci_number_and_size as nas
	on
		cac.date_label= nas.date_label and
		cac.experimental_well_label= nas.experimental_well_label and
		cac.timepoint= nas.timepoint
	where -- filtering
		cac.date_label= (select * from cte_wt_characterisation_date_label) and
		cac.timepoint > p_initial_timepoints_skipped;
end//
delimiter ;
-- call p_wt_characterisation_data(1, "basic");


-- FIGURE 2: number of foci per cell (single cell data)
# stored procedure p_number_of_foci_single_cell_data_two_timepoints called in python script 'Figure2_single_cell_agg_counts' 
# returns data from WT characterisation experiment; data on number of foci per single cell (not averaged per well), 2 timepoints: well labels, timepoint fields, cell id (based on FOV-object_id), number of foci, conditions fields
# data on both As-exposed and control cells
# arguments: 'p_reference_timepoint': decimal, reference timepoint (late formation stage)
#			 'p_selected_timepoint': decimal, selected timepoint (late relocation & fusion stage)
drop procedure if exists hc_microscopy_data_v2.p_number_of_foci_single_cell_data_two_timepoints;
delimiter //
create procedure hc_microscopy_data_v2.p_number_of_foci_single_cell_data_two_timepoints(in p_reference_timepoint decimal(4,1), in p_selected_timepoint decimal(4,1))
begin
	select
		*
    	from
    	(
	with 
	cte_experiment_id as -- relevant experiments
	(
	select distinct
		e.date_label
	from 
		experiments as e
	inner join
		experiment_types as et
	on
		e.experiment_type_id=et.experiment_type_id
	where
		et.experiment_type= 'WT characterisation' and
		et.experiment_subtype= 'basic'
	),
	cte_microscopy_interval as -- microscopy interval for WT characterisation experiments
	(
	select distinct
		microscopy_interval_min
	from
		experiments
	where
		date_label in (select * from cte_experiment_id)
	),
	cte_microscopy_initial_delay as -- microscopy initial delay for for WT characterisation experiments
	(
	select distinct
		microscopy_initial_delay_min
	from
		experiments
	where
		date_label in (select * from cte_experiment_id)
	)
	select -- selected fields
		scd_fna.experimental_well_label,
		scd_fna.timepoint,
		scd_fna.timepoint * (select * from cte_microscopy_interval) - ((select * from cte_microscopy_interval) - (select * from cte_microscopy_initial_delay)) as timepoint_minutes,
		scd_fna.fov_cell_id,
		scd_fna.number_of_foci,
		e.tested_metal,
		sacm.metal_concentration,
		sacm.metal_concentration_unit
	from
		experimental_data_scd_foci_number_and_area as scd_fna
	inner join
		experiments as e
	on 
		e.date_label= scd_fna.date_label
	inner join 
		strains_and_conditions_main as sacm
	on
		scd_fna.date_label= sacm.date_label and
		scd_fna.experimental_well_label= sacm.experimental_well_label
	where
		scd_fna.date_label in (select * from cte_experiment_id)
	) as a
    	where -- filter down to 2 timepoints (reference and selected)
		a.timepoint_minutes = p_reference_timepoint or
        a.timepoint_minutes = p_selected_timepoint;
end //
delimiter ;
-- call p_number_of_foci_single_cell_data_two_timepoints(98, 315);


-- Figure 3 (A): TS mutant library screening 1st round
# stored procedure p_ts_screen_first_round_all_data called in python script 'Figure3_TS_screening_first_round'
# returns all data from the 1st round of TS mutant library screeniing (only 'Good' quality data): date label, plate label, well label, mutant label fields, timepoint fields, number of cells, number of cells with foci, percentage of cells with foci, size and number of foci 
# data on As-exposed cells
# arguments: 'p_initial_timepoints_skipped': int, number of initial timepoints skipped (generally low quality data from initital timepoint, use at least 1)
drop procedure if exists hc_microscopy_data_v2.p_ts_screen_first_round_all_data;
delimiter //
create procedure hc_microscopy_data_v2.p_ts_screen_first_round_all_data(in p_initital_skipped_timepoints int)
begin
	with
	cte_TS_screen_first_round_dates as -- date labels for selected experiments (based on experiment type/subtype: 'TS collection screening'/'first round' and data quality: 'Good')
	(
	select distinct
		e.date_label
	from 
		experiments as e
	inner join 
		experiment_types as et
	on
		e.experiment_type_id= et.experiment_type_id
	where
		e.data_quality = 'Good' and
		et.experiment_type= 'TS collection screening' and
		et.experiment_subtype= 'first round'
	),
	cte_microscopy_interval_min as -- microscopy interval (in minutes) based on data quality and experiment type/subtype
	(
	select distinct
		microscopy_interval_min
	from 
		experiments
	where
		date_label in (select * from cte_TS_screen_first_round_dates)
	),
	cte_initital_delay_min as -- microscopy initial delay (in minutes) based on data quality and experiment type/subtype
	(
	select distinct
		microscopy_initial_delay_min
	from 
		experiments
	where
		date_label in (select * from cte_TS_screen_first_round_dates)
	)
	select -- selected columns and calculations
		scm.date_label,
		scm.collection_plate_label,
		scm.corresponding_collection_well_label,
		scm.mutated_gene_systematic_name,
		scm.mutated_gene_standard_name,
		scm.mutation,
        cac.timepoint,
		cac.timepoint * (select * from cte_microscopy_interval_min) - ((select * from cte_microscopy_interval_min) - (select * from cte_initital_delay_min)) as timepoint_minutes,
        (cac.timepoint * (select * from cte_microscopy_interval_min) - ((select * from cte_microscopy_interval_min) - (select * from cte_initital_delay_min)))/60 as timepoint_hours,
		cac.number_of_cells,
        cac.number_of_cells_with_foci,
		cac.number_of_cells_with_foci/cac.number_of_cells*100 as percentage_of_cells_with_foci,
		nas.avg_number_of_foci_per_cell,
		nas.avg_size_single_focus
	from
		strains_and_conditions_main as scm
	inner join
		experimental_data_sbw_cell_area_and_counts as cac
	on
		scm.date_label= cac.date_label and 
		scm.experimental_well_label= cac.experimental_well_label
	inner join 
		experimental_data_sbw_foci_number_and_size as nas
	on
		cac.date_label= nas.date_label and 
		cac.experimental_well_label= nas.experimental_well_label and
		cac.timepoint=nas.timepoint
	where -- filtering to data from TS screening 1st round
		scm.date_label in (select * from cte_TS_screen_first_round_dates) and
        cac.timepoint > p_initital_skipped_timepoints;
end //
delimiter ;
-- call p_ts_screen_first_round_all_data(3);


-- Figure 3 (B): data on clusters and enrichments
# stored procedure p_ts_screen_first_round_all_data called in python script 'Figure3_enrichments'
# returns data on all clusters and corresponding enrichments based on hits from TS screening 1st round
# arguments: 'p_strength_threshold': minimal enrichment strength, 'p_fdr_threshold': maximal false discovery rate (fdr)
drop procedure if exists hc_microscopy_data_v2.p_clusters_enrichments;
delimiter //
create procedure hc_microscopy_data_v2.p_clusters_enrichments(in p_strength_threshold decimal(3,2), in p_fdr_threshold decimal(7,6))
begin

	with
	cte_unique_effect_stage_cluster as -- unique list of effect-stage-cluster combinations 
	(
	select distinct
		sl.effect_stage_label_id,
		sl.effect_stage_label,
		c.cluster_id,
		c.cluster_label
	from
		hits_clusters as hc
	inner join
		clusters as c
	on
		c.cluster_id= hc.cluster_id
	inner join
		effect_stage_labels as sl
	on
		hc.effect_stage_label_id= sl.effect_stage_label_id
	where
		cluster_label != 'no cluster\r'
	),
	cte_enrichments_per_cluster as -- list of clusters with corresponding number of enrichments
	(
	select
		cluster_id,
		count(distinct enrichment_id) as enrichments_per_cluster
	from
		cluster_enrichment
	group by
		cluster_id
	),
	cte_nodes_per_cluster as -- list of clusters with corresponding number of nodes/genes
	(
	select
		cluster_id,
		count(distinct hit_systematic_name) as nodes_per_cluster
	from
		hits_clusters
	group by
		cluster_id
	)
	select -- selected fields
		cte1.effect_stage_label,
		replace(replace(replace(replace(replace(cte1.cluster_label, '_', ' '), 'dna', 'DNA'), 'mrna', 'mRNA'), 'rna', 'RNA'), 'sec62 63', 'Sec62/Sec63') as cluster_label,
		case
			when e.go_enrichment_category= 'GO Process' then 'GO Biological Process'
			else 'GO Cellular Compartment'
		end as go_enrichment_category,
		e.go_id as go_enrichment_id,
		e.enrichment_description,
		ce.strength,
		round(ce.fdr, 4) as fdr,
		cte3.nodes_per_cluster,
		cte2.enrichments_per_cluster
	from
		cluster_enrichment as ce
	inner join
		enrichments as e
	on 
		ce.enrichment_id=e.enrichment_id
	inner join 
		cte_unique_effect_stage_cluster as cte1
	on
		cte1.cluster_id=ce.cluster_id
	inner join
		cte_enrichments_per_cluster as cte2
	on
		cte2.cluster_id= ce.cluster_id
	inner join
		cte_nodes_per_cluster as cte3
	on
		cte3.cluster_id= ce.cluster_id
	where
		ce.strength >=  p_strength_threshold and
        	ce.fdr < p_fdr_threshold
	order by
		cte1.effect_stage_label asc,
		cte3.nodes_per_cluster desc,
        	cluster_label asc,
        	go_enrichment_category asc,
        	ce.strength desc;
        
end //
delimiter ;

-- call hc_microscopy_data_v2.p_clusters_enrichments(1.25, 0.05);
