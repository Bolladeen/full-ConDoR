# %% [markdown]
# # Draw pie charts

# %%
cn_clone_assignment_df_updated_for_pie = cn_clone_assignment_df_updated.set_index('sample')

# %%
import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from pathlib import Path


# load clone assignment df
ete_tree_dir = Path("/juno/work/iacobuzc/haochen/Tapestri_batch2/analysis/condor_downstream/HZ_ete_trees_refined_subclonal_snvs")

# load MASTER df
sample_sheet = Path("../references/Tapestri_batch2_samples_MASTER.xlsx")
sample_sheet_df = pd.read_excel(sample_sheet)
sample_organ_map = dict(zip(sample_sheet_df['sample'], sample_sheet_df['HZ_official_site_name']))

# %% Load clone assignment df
patient_names = ["RA16_29", "RA17_13", "RA17_22"]
for patient_name in patient_names:
    wd = ete_tree_dir / patient_name / "sample_pie_charts"
    wd.mkdir(exist_ok=True, parents=True)
    clone_compo = ete_tree_dir / patient_name / f"{patient_name}_clone_compo.csv"
    clone_compo_df = pd.read_csv(clone_compo, index_col=0)
    unique_cluster_ids_sorted = np.sort(clone_compo_df.columns.astype(int))

    cn_clone_palette = dict(zip(unique_cluster_ids_sorted, np.array(px.colors.qualitative.Pastel)[unique_cluster_ids_sorted]))

    # # %% Make a pie chart for each sample

    for sample_i in clone_compo_df.index.unique():
        sample_name_mapped = sample_organ_map[sample_i]
        # filter columns that have 0
        filtered_clone_compo_df = clone_compo_df.loc[sample_i].loc[clone_compo_df.loc[sample_i] > 0]
        fig = go.Figure(data=[go.Pie(
            labels=filtered_clone_compo_df.index, 
            values=filtered_clone_compo_df.values, 
            hole=.3
        )])
        colors = [cn_clone_palette[int(i)] for i in filtered_clone_compo_df.index]
        fig.update_traces(
            # hoverinfo='label+percent', 
            # textinfo='value', 
            # textinfo = 'label+percent',
            textfont_size=20,
            marker=dict(colors=colors, 
            line=dict(color='#000000', width=2))
        )
        fig.update_layout(
            paper_bgcolor='rgba(0,0,0,0)',
            plot_bgcolor='rgba(0,0,0,0)',
            showlegend=False,
        )
        fig.write_image(str(wd /f"{sample_name_mapped}_clone_compo_pie.png"), scale=2)




# %%
