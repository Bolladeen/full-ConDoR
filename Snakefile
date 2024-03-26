from pathlib import Path
configfile: "config.yaml"
workdir: config["workdir"]

rule all:
  input:
    expand("results/opt_nclones/{dataset}/optimal_cell_assignments.csv", dataset=config["datasets"]),
		# expand("results/{dataset}/clonal_refinement/refined_cell_assignments.csv", dataset=config["datasets"]),
    expand("results/condor_inputs/{dataset}/character_vaf_matrix.csv", dataset=config["datasets"]),
    expand("results/condor_outputs/{dataset}/condor_outputs.log", dataset=config["datasets"]),
    expand("results/condor_outputs/{dataset}/heatmaps/condor_solution_heatmap.png", dataset=config["datasets"]),

rule opt_nclones:
  params:
    falcon_solutions=config['falcon_solutions'],
    nclones=config['nclones'],
    output_dir='results/opt_nclones/{dataset}/'
  output:
    result_clone_assignment='results/opt_nclones/{dataset}/optimal_cell_assignments.csv',
    result_cn_profile='results/opt_nclones/{dataset}/optimal_clone_profiles.csv',
  log:
    std='results/opt_nclones/{dataset}/opt_nclones.log',
    err='results/opt_nclones/{dataset}/opt_nclones.err.log',
  shell:
      """
      python scripts/select_optimal_nclones.py \
        -d {wildcards.dataset} \
        -c {params.nclones} \
        -i {params.falcon_solutions}{wildcards.dataset}/solutions/  \
        -a {output.result_clone_assignment}  \
        -p {output.result_cn_profile} \
        --output_dir {params.output_dir} \
        1> {log.std} 2> {log.err}
      """
rule condor_inputs:
  input:
    refined_clone_assignment='results/opt_nclones/{dataset}/optimal_cell_assignments.csv',
	output:
		character_vaf="results/condor_inputs/{dataset}/character_vaf_matrix.csv",
		character_mat="results/condor_inputs/{dataset}/character_bin_matrix.csv",
		alt_readcounts="results/condor_inputs/{dataset}/alt_readcounts.csv",
		total_readcounts="results/condor_inputs/{dataset}/total_readcounts.csv",
		germline_mutations="results/condor_inputs/{dataset}/germline_mutations.txt",
		somatic_mutations="results/condor_inputs/{dataset}/somatic_mutations.txt",
	params: 
		condor_input_script = Path(config["condor_pipeline_scripts_dir"]) / "generate_condor_input.py",
		hdf5_directory=config["raw_data_directory"],
		annotated_mutations=config["annotated_mutations"],
	log:
		std="results/condor_inputs/{dataset}/condor_inputs.log",
		err="results/condor_inputs/{dataset}/condor_inputs.err.log",
	conda: 
		# "envs/mosaic-custom.yaml",
		"mosaic"
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
			-snvs {params.annotated_mutations}{wildcards.dataset}-patient-all_vars-voi.hz.txt \
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
    character_mat="results/condor_inputs/{dataset}/character_bin_matrix.csv",
    alt_readcounts="results/condor_inputs/{dataset}/alt_readcounts.csv",
    total_readcounts="results/condor_inputs/{dataset}/total_readcounts.csv",
    germline_mutations="results/condor_inputs/{dataset}/germline_mutations.txt",
    somatic_mutations="results/condor_inputs/{dataset}/somatic_mutations.txt",
    cn_profiles='results/opt_nclones/{dataset}/optimal_clone_profiles.csv',
	params:
		fast_condor_script = config["fast_condor_script"],
		amplicon_parameters = config["amplicon_parameters"],
		output_prefix = "results/condor_outputs/{dataset}/out",
		hdf5_directory = config["raw_data_directory"],
	log:
		std="results/condor_outputs/{dataset}/condor_outputs.log",
		err="results/condor_outputs/{dataset}/condor_outputs.err.log",
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
			-m {params.amplicon_parameters} \
			-c {params.hdf5_directory} \
			-d {wildcards.dataset} \
      --scr \
      --cnp {input.cn_profiles} \
			1> {log.std} 2> {log.err}
		"""

rule generate_heatmaps:
	input:
		expand("results/condor_outputs/{dataset}/condor_outputs.log", dataset=config["datasets"]),
		vaf_matrix="results/condor_inputs/{dataset}/character_vaf_matrix.csv",
		alt_readcounts="results/condor_inputs/{dataset}/alt_readcounts.csv",
		total_readcounts="results/condor_inputs/{dataset}/total_readcounts.csv",
		character_mat="results/condor_inputs/{dataset}/character_bin_matrix.csv",
		germline_mutations="results/condor_inputs/{dataset}/germline_mutations.txt",
		somatic_mutations="results/condor_inputs/{dataset}/somatic_mutations.txt",
	output:
		vaf_heatmap="results/condor_outputs/{dataset}/heatmaps/vaf_heatmap.png",
		solution_heatmap="results/condor_outputs/{dataset}/heatmaps/condor_solution_heatmap.png",
	params:
		heatmap_script = Path(config["condor_pipeline_scripts_dir"]) / "generate_heatmaps.py",
		hdf5_directory=config["raw_data_directory"],
		condor_solution="results/condor_outputs/{dataset}/out_B.csv",
		amplicon_parameters=config["amplicon_parameters"],
	conda: "condor"
	log:
		std="results/condor_outputs/{dataset}/heatmaps/heatmaps.log",
		err="results/condor_outputs/{dataset}/heatmaps/heatmaps.err.log",
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
			-i {params.amplicon_parameters} \
			1> {log.std} 2> {log.err}
		"""

	


