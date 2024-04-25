mamba create \
	-n condor \
	-c conda-forge \
	ete3=3.1.3=pyhd8ed1ab_0 python=3.11
mamba activate condor
mamba install \
	-c conda-forge -c bioconda \
	numpy pandas networkx "snakemake<8" biopython seaborn openpyxl ipykernel
mamba install -c plotly plotly
pip install -U kaleido # for plotly plotting
conda config --add channels https://conda.anaconda.org/gurobi
pip install gurobipy==10.0.3 # version 11 is not compatible with condor

# obtain license and put it to ~/gurobi.lic or set up $GRB_LICENSE_FILE env variable

# need gurobi version > 10.0.0
# https://support.gurobi.com/hc/en-us/articles/4416277022353-Error-10024-Web-license-service-only-available-for-container-environments

# pip install /home/zhangh5/work/Tapestri_project/TapVarCallSmk/resources/mosaic
# mamba install -c conda-forge "h5py>=2.10.0"
# pip install umap-learn>=0.4.6