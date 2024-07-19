from pathlib import Path
import os

workdir: config["workdir"]

rule all:
	input:
		# expand("{patient}/opt_nclones/optimal_cell_assignments.csv", patient=config["patients"]),
		# expand("{patient}/clonal_refinement/refined_cell_assignments.csv", patient=config["patients"]),
		expand("condor_inputs/{patient}/character_vaf_matrix.csv", patient=config["patients"]),
		expand("condor_outputs/{patient}/out_tree.newick", patient=config["patients"]),
		expand("condor_outputs/{patient}/heatmaps/condor_solution_heatmap.png", patient=config["patients"]),
		expand("condor_downstream/{patient}/{patient}_ETE_tree.refined.subclonal_snvs.png", patient=config["patients"]),
		expand("condor_downstream/{patient}/{patient}_final_sc_clone_assignment.csv", patient=config["patients"]),
		expand("pre_condor_sc_heatmaps/{patient}_DNA_heatmap.pdf", patient=config["patients"]),
		expand("post_condor_sc_heatmaps/{patient}_DNA_heatmap.pdf", patient=config["patients"]),

# sanity check first
patient_subclonal_snv_yaml = {}
for patient_i in config["datasets"]:
	if len(list(Path(config["raw_data_directory"]).glob(f"{patient_i}*.h5"))) != 1:
		raise ValueError(f"Expected 1 h5 file for {patient_i}, found {len(list(Path(config['raw_data_directory']).glob(f'{patient_i}*.h5')))}")
	if len(list(Path(config["falcon_solutions"]).glob(f"{patient_i}*assignment*.csv"))) != 1:
		raise ValueError(f"Expected 1 clone-assignment file for {patient_i}, found {len(list(Path(config['falcon_solutions']).glob(f'{patient_i}*assignment*.csv')))}")
	if len(list(Path(config["falcon_solutions"]).glob(f"{patient_i}*clone_profile*.csv"))) != 1:
		raise ValueError(f"Expected 1 clone-profile file for {patient_i}, found {len(list(Path(config['falcon_solutions']).glob(f'{patient_i}*clone_profile*.csv')))}")
	snv_list = list(Path(config["annotated_mutations"]).glob(f"{patient_i}*voi*.txt"))
	if len(snv_list) != 1:
		raise ValueError(f"Expected 1 annotated mutation file for {patient_i}, found {len(snv_list)}")
	subclonal_snv_yaml = list(Path(config["subclonal_mutations"]).glob(f"{patient_i}*subclonal_mutations*.yaml"))
	if len(subclonal_snv_yaml) >0:
		print(f"Found {len(subclonal_snv_yaml)} subclonal mutation files for {patient_i}")
		if len(subclonal_snv_yaml) > 1:
			raise ValueError(f"Expected 1 subclonal mutation file for {patient_i}, found {len(subclonal_snv_yaml)}")
		patient_subclonal_snv_yaml[patient_i] = subclonal_snv_yaml[0]
	else:
		patient_subclonal_snv_yaml[patient_i] = None

def __get_input_h5(wildcards):
	return list(Path(config["raw_data_directory"]).glob(f"{wildcards.patient}*.h5"))[0]
def __get_falcon_solution_cell_assignment(wildcards):
	return list(Path(config["falcon_solutions"]).glob(f"{wildcards.patient}*assignment*.csv"))[0]
def __get_falcon_solution_cn_profiles(wildcards):
	return list(Path(config["falcon_solutions"]).glob(f"{wildcards.patient}*clone_profile*.csv"))[0]
def __get_annotated_mutations(wildcards):
    return list(Path(config["annotated_mutations"]).glob(f"{wildcards.patient}*voi*.txt"))[0]


rule condor_inputs:
	input:
		h5 = __get_input_h5,
		refined_clone_assignment = __get_falcon_solution_cell_assignment,
		# "opt_nclones/{patient}.sample_sc_clone_assignment.updated.csv",
		annotated_mutations = __get_annotated_mutations,
	output:
		character_vaf="condor_inputs/{patient}/character_vaf_matrix.csv",
		character_mat="condor_inputs/{patient}/character_bin_matrix.csv",
		alt_readcounts="condor_inputs/{patient}/alt_readcounts.csv",
		total_readcounts="condor_inputs/{patient}/total_readcounts.csv",
		germline_mutations="condor_inputs/{patient}/germline_mutations.txt",
		somatic_mutations="condor_inputs/{patient}/somatic_mutations.txt",
	params: 
		condor_input_script = Path(config["condor_pipeline_scripts_dir"]) / "1a_generate_condor_input.py",
	log:
		std="condor_inputs/{patient}/condor_inputs.log",
		err="condor_inputs/{patient}/condor_inputs.err.log",
	conda: 
		# "envs/mosaic-custom.yaml",
		"mosaic",
	threads: lambda wildcards, attempt: attempt * 4
	resources:
		mem_mb = 8000,
		time_min = lambda wildcards, attempt: attempt * 59,
	retries: 2
	shell:
		"""
		python {params.condor_input_script} \
			-d {wildcards.patient} \
			-h5 {input.h5} \
			-i {input.refined_clone_assignment} \
			-snvs {input.annotated_mutations} \
			-v {output.character_vaf} \
			-m {output.character_mat} \
			-a {output.alt_readcounts} \
			-t {output.total_readcounts} \
			-g {output.germline_mutations} \
			-s {output.somatic_mutations} \
			1> {log.std} 2> {log.err}
		"""

rule generate_pre_condor_sc_heatmaps:
	input:
		h5 = __get_input_h5,
		refined_clone_assignment = __get_falcon_solution_cell_assignment,
		# "opt_nclones/{patient}.sample_sc_clone_assignment.updated.csv",
		annotated_mutations = __get_annotated_mutations,
	output:
		pre_condor_sc_heatmap="pre_condor_sc_heatmaps/{patient}_DNA_heatmap.pdf",
	params:
		sc_heatmap_script = Path(config["condor_pipeline_scripts_dir"]) / "1b_generate_sc_heatmaps.py",
		output_dir = "pre_condor_sc_heatmaps",
	conda: "mosaic"
	log:
		std="pre_condor_sc_heatmaps/logs/{patient}_heatmap.log",
		err="pre_condor_sc_heatmaps/logs/{patient}_heatmap.err.log",
	threads: lambda wildcards, attempt: attempt * 4
	resources:
		mem_mb = 8000,
		time_min = lambda wildcards, attempt: attempt * 59,
	retries: 1
	shell:
		"""
		python {params.sc_heatmap_script} \
			--patient_name {wildcards.patient} \
			--patient_h5 {input.h5} \
			--snv_f {input.annotated_mutations} \
			--clone_assignment {input.refined_clone_assignment} \
			--output_dir {params.output_dir} \
			1> {log.std} 2> {log.err}
		"""

rule fast_condor:
	input:
		character_mat="condor_inputs/{patient}/character_bin_matrix.csv",
		alt_readcounts="condor_inputs/{patient}/alt_readcounts.csv",
		total_readcounts="condor_inputs/{patient}/total_readcounts.csv",
		germline_mutations="condor_inputs/{patient}/germline_mutations.txt",
		somatic_mutations="condor_inputs/{patient}/somatic_mutations.txt",
		# refined_clone_assignment = lambda wildcards: Path(config["falcon_solutions"]) / f"{wildcards.patient}.sample_sc_clone_assignment.updated.csv",
		refined_clone_profiles = __get_falcon_solution_cn_profiles,
	output:
		newick_tree_file = "condor_outputs/{patient}/out_tree.newick",
		ete_tree_pickle = "condor_outputs/pickle_files/{patient}_self.solT_cell",
	params:
		fast_condor_script = config["fast_condor_script"],
		amplicon_coordinates_file = config["amplicon_coordinates_file"],
		output_prefix = "condor_outputs/{dataset}/out",
		subclonal_mutations = lambda wildcards: patient_subclonal_snv_yaml[wildcards.dataset],
	log:
		std="condor_outputs/{patient}/condor_outputs.log",
		err="condor_outputs/{patient}/condor_outputs.err.log",
	conda: "condor"
	threads: lambda wildcards, attempt: (attempt**2) * 4
	resources:
		mem_gb = 16,
		time_min = lambda wildcards, attempt: attempt * 119,
	retries: 2
	shell:
		"""
		python {params.fast_condor_script} \
			-i {input.character_mat} \
			-r {input.total_readcounts} \
			-v {input.alt_readcounts} \
			-s {input.germline_mutations} \
			-s2 {input.somatic_mutations} \
			-o {params.output_prefix} \
			-m {params.amplicon_coordinates_file} \
			-d {wildcards.patient} \
			--scr \
			--cnp {input.refined_clone_profiles} \
			--subclonal_mutations {params.subclonal_mutations} \
			1> {log.std} 2> {log.err}
		"""

rule generate_solution_heatmaps:
	input:
		expand("condor_outputs/{patient}/condor_outputs.log", patient=config["patients"]),
		vaf_matrix="condor_inputs/{patient}/character_vaf_matrix.csv",
		alt_readcounts="condor_inputs/{patient}/alt_readcounts.csv",
		total_readcounts="condor_inputs/{patient}/total_readcounts.csv",
		character_mat="condor_inputs/{patient}/character_bin_matrix.csv",
		germline_mutations="condor_inputs/{patient}/germline_mutations.txt",
		somatic_mutations="condor_inputs/{patient}/somatic_mutations.txt",
		annotated_mutations = __get_annotated_mutations,
	output:
		vaf_heatmap="condor_outputs/{patient}/heatmaps/vaf_heatmap.png",
		solution_heatmap="condor_outputs/{patient}/heatmaps/condor_solution_heatmap.png",
	params:
		solution_heatmap_script = Path(config["condor_pipeline_scripts_dir"]) / "3a_generate_solution_heatmaps.py",
		condor_solution="condor_outputs/{patient}/out_B.csv",
		amplicon_coordinates_file=config["amplicon_coordinates_file"],
	conda: "condor"
	log:
		std="condor_outputs/{patient}/heatmaps/heatmaps.log",
		err="condor_outputs/{patient}/heatmaps/heatmaps.err.log",
	threads: lambda wildcards, attempt: attempt * 4
	resources:
		mem_mb = 8000,
		time_min = lambda wildcards, attempt: attempt * 59,
	retries: 1
	shell:
		"""
		python {params.solution_heatmap_script} \
			-d {wildcards.patient} \
			-c {input.character_mat} \
			-v {input.vaf_matrix} \
			-s {params.condor_solution} \
			-g {input.germline_mutations} \
			-m {input.somatic_mutations} \
			-t {input.total_readcounts} \
			-a {input.alt_readcounts} \
			-o {output.solution_heatmap} \
			-p {output.vaf_heatmap} \
			-i {params.amplicon_coordinates_file} \
			-snvs {input.annotated_mutations} \
			1> {log.std} 2> {log.err}
		"""

rule generate_ete_trees:
	input:
		ete_tree_pickle = "condor_outputs/pickle_files/{patient}_self.solT_cell",
		snv_ann_f = __get_annotated_mutations,
		amplicon_gene_mapping_f = config["amplicon_coordinates_file"],
		sample_name_map_f = config["sample_name_map"],
	output:
		ete_tree_refined = "condor_downstream/{patient}/{patient}_ETE_tree.refined.subclonal_snvs.png",
		final_clone_assignment = "condor_downstream/{patient}/{patient}_final_sc_clone_assignment.csv",
	params:
		condor_tree_to_ete_tree_script = os.path.join(config["condor_downstream_scripts_dir"], "4a_make_ete_tree_with_subclonal_snvs.py"),
		output_dir = "condor_downstream",
	conda: "condor"
	log:
		std="condor_downstream/logs/{patient}_ete_tree.log",
		err="condor_downstream/logs/{patient}_ete_tree.err.log",
	threads: lambda wildcards, attempt: attempt * 4
	resources:
		mem_mb = 8000,
		time_min = lambda wildcards, attempt: attempt * 59,
	retries: 1
	shell:
		"""
		python {params.condor_tree_to_ete_tree_script} \
			--amp_gene_map {input.amplicon_gene_mapping_f} \
			--sample_name_map {input.sample_name_map_f} \
			--patient_name {wildcards.patient} \
			--condor_tree_pickle {input.ete_tree_pickle} \
			--snv_ann_f {input.snv_ann_f} \
			--output_dir {params.output_dir} \
			1> {log.std} 2> {log.err}
		"""


rule generate_post_condor_sc_heatmaps:
	input:
		h5 = __get_input_h5,
		post_condor_clone_assignment = "condor_downstream/{patient}/{patient}_final_sc_clone_assignment.csv",
		annotated_mutations = __get_annotated_mutations,
	output:
		post_condor_sc_heatmap="post_condor_sc_heatmaps/{patient}_DNA_heatmap.pdf",
	params:
		sc_heatmap_script = Path(config["condor_pipeline_scripts_dir"]) / "1b_generate_sc_heatmaps.py",
		output_dir = "post_condor_sc_heatmaps",
	conda: "mosaic"
	log:
		std="post_condor_sc_heatmaps/logs/{patient}_heatmap.log",
		err="post_condor_sc_heatmaps/logs/{patient}_heatmap.err.log",
	threads: lambda wildcards, attempt: attempt * 4
	resources:
		mem_mb = 8000,
		time_min = lambda wildcards, attempt: attempt * 59,
	retries: 1
	shell:
		"""
		python {params.sc_heatmap_script} \
			--patient_name {wildcards.patient} \
			--patient_h5 {input.h5} \
			--snv_f {input.annotated_mutations} \
			--clone_assignment {input.post_condor_clone_assignment} \
			--output_dir {params.output_dir} \
			--post_condor \
			1> {log.std} 2> {log.err}
		"""
	


