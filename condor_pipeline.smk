from pathlib import Path

workdir: config["workdir"]

rule all:
	input:
		# expand("{dataset}/opt_nclones/optimal_cell_assignments.csv", dataset=config["datasets"]),
		# expand("{dataset}/clonal_refinement/refined_cell_assignments.csv", dataset=config["datasets"]),
		expand("condor_inputs/{dataset}/character_vaf_matrix.csv", dataset=config["datasets"]),
		expand("condor_outputs/{dataset}/out_tree.newick", dataset=config["datasets"]),
		expand("condor_outputs/{dataset}/heatmaps/condor_solution_heatmap.png", dataset=config["datasets"]),

rule condor_inputs:
	input:
		refined_clone_assignment = lambda wildcards: f"{config['falcon_solutions']}/{wildcards.dataset}.sample_sc_clone_assignment.updated.csv",
		# "opt_nclones/{dataset}.sample_sc_clone_assignment.updated.csv",
	output:
		character_vaf="condor_inputs/{dataset}/character_vaf_matrix.csv",
		character_mat="condor_inputs/{dataset}/character_bin_matrix.csv",
		alt_readcounts="condor_inputs/{dataset}/alt_readcounts.csv",
		total_readcounts="condor_inputs/{dataset}/total_readcounts.csv",
		germline_mutations="condor_inputs/{dataset}/germline_mutations.txt",
		somatic_mutations="condor_inputs/{dataset}/somatic_mutations.txt",
	params: 
		condor_input_script = Path(config["condor_pipeline_scripts_dir"]) / "generate_condor_input.py",
		hdf5_directory=config["raw_data_directory"],
		annotated_mutations=config["annotated_mutations"],
	log:
		std="condor_inputs/{dataset}/condor_inputs.log",
		err="condor_inputs/{dataset}/condor_inputs.err.log",
	conda: 
		"envs/mosaic-custom.yaml",
	threads: lambda wildcards, attempt: attempt * 4
	resources:
		mem_mb = 8000,
		time_min = lambda wildcards, attempt: attempt * 59,
	retries: 2
	shell:
		"""
		python {params.condor_input_script} \
			-d {wildcards.dataset} \
			-l {params.hdf5_directory} \
			-i {input.refined_clone_assignment} \
			-snvs {params.annotated_mutations}{wildcards.dataset}-patient-all_vars-voi.hz_curated.txt \
			-v {output.character_vaf} \
			-m {output.character_mat} \
			-a {output.alt_readcounts} \
			-t {output.total_readcounts} \
			-g {output.germline_mutations} \
			-s {output.somatic_mutations} \
			1> {log.std} 2> {log.err}
		"""

rule fast_condor:
	input:
		character_mat="condor_inputs/{dataset}/character_bin_matrix.csv",
		alt_readcounts="condor_inputs/{dataset}/alt_readcounts.csv",
		total_readcounts="condor_inputs/{dataset}/total_readcounts.csv",
		germline_mutations="condor_inputs/{dataset}/germline_mutations.txt",
		somatic_mutations="condor_inputs/{dataset}/somatic_mutations.txt",
	output:
		newick_tree_file = "condor_outputs/{dataset}/out_tree.newick",
	params:
		fast_condor_script = config["fast_condor_script"],
		amplicon_coordinates_file = config["amplicon_coordinates_file"],
		output_prefix = "condor_outputs/{dataset}/out",
		hdf5_directory = config["raw_data_directory"],
	log:
		std="condor_outputs/{dataset}/condor_outputs.log",
		err="condor_outputs/{dataset}/condor_outputs.err.log",
	conda: "condor"
	threads: lambda wildcards, attempt: attempt * 4
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
			-c {params.hdf5_directory} \
			-d {wildcards.dataset} \
			1> {log.std} 2> {log.err}
		"""

rule generate_heatmaps:
	input:
		expand("condor_outputs/{dataset}/condor_outputs.log", dataset=config["datasets"]),
		vaf_matrix="condor_inputs/{dataset}/character_vaf_matrix.csv",
		alt_readcounts="condor_inputs/{dataset}/alt_readcounts.csv",
		total_readcounts="condor_inputs/{dataset}/total_readcounts.csv",
		character_mat="condor_inputs/{dataset}/character_bin_matrix.csv",
		germline_mutations="condor_inputs/{dataset}/germline_mutations.txt",
		somatic_mutations="condor_inputs/{dataset}/somatic_mutations.txt",
	output:
		vaf_heatmap="condor_outputs/{dataset}/heatmaps/vaf_heatmap.png",
		solution_heatmap="condor_outputs/{dataset}/heatmaps/condor_solution_heatmap.png",
	params:
		heatmap_script = Path(config["condor_pipeline_scripts_dir"]) / "generate_heatmaps.py",
		hdf5_directory=config["raw_data_directory"],
		condor_solution="condor_outputs/{dataset}/out_B.csv",
		amplicon_coordinates_file=config["amplicon_coordinates_file"],
	conda: "condor"
	log:
		std="condor_outputs/{dataset}/heatmaps/heatmaps.log",
		err="condor_outputs/{dataset}/heatmaps/heatmaps.err.log",
	threads: lambda wildcards, attempt: attempt * 4
	resources:
		mem_mb = 8000,
		time_min = lambda wildcards, attempt: attempt * 59,
	retries: 2
	shell:
		"""
		python {params.heatmap_script} \
			-d {wildcards.dataset} \
			-c {input.character_mat} \
			-v {input.vaf_matrix} \
			-s {params.condor_solution} \
			-l {params.hdf5_directory} \
			-g {input.germline_mutations} \
			-m {input.somatic_mutations} \
			-t {input.total_readcounts} \
			-a {input.alt_readcounts} \
			-o {output.solution_heatmap} \
			-p {output.vaf_heatmap} \
			-i {params.amplicon_coordinates_file} \
			1> {log.std} 2> {log.err}
		"""

	


