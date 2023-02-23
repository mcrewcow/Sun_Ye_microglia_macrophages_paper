#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().system('pip install scFates')


# In[2]:


import sys
get_ipython().system('{sys.executable} -m pip -q install palantir fa2')


# In[3]:


import warnings
warnings.filterwarnings("ignore")
from anndata import AnnData
import numpy as np
import pandas as pd
import scanpy as sc
import scFates as scf
import palantir
import matplotlib.pyplot as plt
sc.settings.verbosity = 3
sc.settings.logfile = sys.stdout
## fix palantir breaking down some plots
import seaborn
seaborn.reset_orig()
get_ipython().run_line_magic('matplotlib', 'inline')

sc.set_figure_params()
scf.set_figure_pubready()


# In[68]:


sc.pp.filter_genes(adata1,min_cells=3)
sc.pp.normalize_total(adata1)
sc.pp.log1p(adata1,base=10)
sc.pp.highly_variable_genes(adata1)


# In[145]:


sc.pp.pca(adatamergedsubset)
pca_projections = pd.DataFrame(adatamergedsubset.obsm["X_pca"],index=adatamergedsubset.obs_names)


# In[261]:


pca_projections = pd.DataFrame(adatamerged.obsm["X_pca"],index=adatamerged.obs_names)


# In[284]:


pca_projections


# In[146]:


dm_res = palantir.utils.run_diffusion_maps(pca_projections)
ms_data = palantir.utils.determine_multiscale_space(dm_res,n_eigs=4)


# In[147]:


# generate neighbor draph in multiscale diffusion space
adatamergedsubset.obsm["X_palantir"]=ms_data.values
sc.pp.neighbors(adatamergedsubset,n_neighbors=30,use_rep="X_palantir")


# In[148]:


# draw ForceAtlas2 embedding using 3 first PCs as initial positions
adatamergedsubset.obsm["X_pca2d"]=adatamergedsubset.obsm["X_pca"][:,:3]
sc.tl.draw_graph(adatamergedsubset, layout = 'fa',init_pos='X_pca2d')


# In[334]:


adatamerged.obsm["X_fates"]


# In[318]:


adatamerged.obsm["X_fates"]=np.concatenate([adatamerged.obsm["X_draw_graph_fa"], adatamerged.obsm["X_pca"][:,2:3]], axis=1)


# In[351]:


adataWT.obsm["X_fates"]=np.concatenate([adataWT.obsm["X_draw_graph_fa"], adataWT.obs["velocity_pseudotime"].values.reshape(-1,1)], axis=1)


# In[350]:


adataWT.obs["velocity_pseudotime"]


# In[325]:


adatamerged.obs


# In[211]:


sc.tl.umap(adatamerged,n_components=3)


# In[102]:


sc.set_figure_params()
sc.pl.draw_graph(adatamerged,color=['Cx3cr1', 'labels1'], layout = 'fa')


# In[149]:


sc.set_figure_params()
sc.pl.draw_graph(adatamergedsubset,color=['Cx3cr1', 'labels1'], layout = 'fa')


# In[191]:


sc.set_figure_params()
sc.pl.draw_graph(adatamerged,color='Spp1', vmax = 2, projection = '2d', layout = 'fa', save = 'mergedSpp1.pdf')


# In[340]:


sc.set_figure_params()
sc.pl.draw_graph(adatamerged,color=['Cx3cr1', 'labels1'], projection = '2d', layout = 'fa')


# In[1]:


sc.set_figure_params()
sc.pl.draw_graph(adatamerged,color=['Cx3cr1', 'labels1'], projection = '2d', layout = 'fa')


# In[327]:


sc.set_figure_params()
sc.pl.embedding(adatamerged,color=['labels1'], projection = '3d', basis = 'fates')


# In[231]:


color = dict(zip(range(0,16), plt.cm.tab20(range(0,16))))


# In[269]:


adatamerged.obsm


# In[306]:


draw_graph = adatamerged.obsm['X_draw_graph_fr']
draw_graph


# In[283]:


umap = adatamerged.obsm['X_pca']
umap


# In[242]:


for i in range(0,360,2):
    fig = plt.figure(figsize = (10,10))
    ax = fig.add_subplot(projection = '3d')
    ax.scatter(umap[:,0],umap[:,1],umap[:,2], c = adatamerged.obs.seurat_clusters.astype('int').map(color))

    ax.axis('off')
    ax.view_init(20, i)
    
    plt.savefig(f'figs/{i:003}.png',dpi = 300, facecolor = 'white')
    plt.show()


# In[241]:


get_ipython().system('ls')


# In[233]:


sc.pl.umap(adatamerged,color=['Cx3cr1', 'labels1'])


# In[ ]:


get_ipython().system('convert -delay 5 figs/*.png umap.gif')


# In[235]:


sc.pl.umap(adatamerged, color = ['seurat_clusters'])


# In[223]:


sc.pl.umap(adatamerged,color=['Cx3cr1', 'Cell.ident'])


# In[96]:


adatamerged = AnnData.concatenate(adata1,adata2,adata3)


# In[97]:


adatamerged.obs


# In[118]:


scv.pl.proportions(adatamergedvelo)


# In[112]:


scv.pl.velocity_embedding_stream(adatanormalsubset, basis='draw_graph_fa', color = 'labels1', dpi = 300, add_margin = 0.2, figsize = (8,8), arrow_size = 2, density = 6)


# In[56]:


sc.set_figure_params(scanpy = True,format = 'svg')


scv.pl.velocity_embedding_stream(adataWTOIRsubset, basis='draw_graph_fa', color = 'labels1', dpi = 300, add_margin = 0.2, legend_loc = 'right_margin',save="figures/microglia_WTOIRsubset_legends_lessarrows.svg")


# In[57]:


sc.set_figure_params(scanpy = True,format = 'svg')


scv.pl.velocity_embedding_stream(adataKOOIRsubset, basis='draw_graph_fa', color = 'labels1', dpi = 300, add_margin = 0.2, legend_loc = 'right_margin',save="figures/microglia_KOOIRsubset_legends_lessarrows.svg")


# In[58]:


sc.set_figure_params(scanpy = True,format = 'svg')


scv.pl.velocity_embedding_stream(adatanormalsubset, basis='draw_graph_fa', color = 'labels1', dpi = 300, add_margin = 0.2, legend_loc = 'right_margin',save="figures/microglia_normalsubset_legends_lessarrows.svg")


# In[50]:


import scvelo as scv


# In[122]:


scv.pl.velocity_embedding_grid(adataWTOIR, basis='draw_graph_fa', arrow_length = 4, arrow_size = 3, color = 'labels1', dpi = 300, add_margin = 0.2,save="figures/microglia_WTOIR.svg")


# In[124]:


scv.pl.velocity_embedding_grid(adataKOOIR, basis='draw_graph_fa', arrow_length = 4, arrow_size = 3, color = 'labels1', dpi = 300, add_margin = 0.2,save="figures/microglia_KOOIR.svg")


# In[148]:


sc.pl.draw_graph(adataKOOIRsubset,color=["Cx3cr1",'Tmem119'],save="KOOIRcxtmem.svg")


# In[156]:


adataKOOIRsubset.uns['neighbors']['distances'] = adataKOOIRsubset.obsp['distances']
adataKOOIRsubset.uns['neighbors']['connectivities'] = adataKOOIRsubset.obsp['connectivities']


# In[157]:


scv.tl.paga(adataKOOIRsubset, groups='labels1')
df = scv.get_df(adataKOOIRsubset, 'paga/transitions_confidence', precision=2).T
df.style.background_gradient(cmap='Blues').format('{:.2g}')


# In[158]:


scv.pl.paga(adataKOOIRsubset, basis='draw_graph_fa', size=50, alpha=.1,
            min_edge_width=2, node_size_scale=1.5, save = 'figures/KOOIRsubsetpaga.pdf')


# In[169]:


scv.tl.recover_dynamics(adatanormalsubset)
scv.tl.latent_time(adatanormalsubset)


# In[170]:


scv.pl.scatter(adatanormalsubset, color='latent_time', color_map='gnuplot', size=80, basis = 'draw_graph_fa', save = 'figures/latenttimenormal.pdf')


# In[171]:


top_genes = adatanormalsubset.var['fit_likelihood'].sort_values(ascending=False).index[:300]
scv.pl.heatmap(adatanormalsubset, var_names=top_genes, sortby='latent_time', col_color='labels1', n_convolve=100, save = 'figures/latenttimeheatmapnormal.pdf')


# In[173]:


top_genes = adatanormalsubset.var['fit_likelihood'].sort_values(ascending=False).index
scv.pl.scatter(adatanormalsubset, basis=top_genes[:15], ncols=5, frameon=False)


# In[180]:


adataKOOIRsubset.var


# In[74]:


adata1 = sc.read_h5ad(r"/mnt/c/Users/Emil/10X/scretina/sunyenorm.h5ad")


# In[75]:


adata2 = sc.read_h5ad(r"/mnt/c/Users/Emil/10X/scretina/sunyeOIRKO.h5ad")


# In[76]:


adata3 = sc.read_h5ad(r"/mnt/c/Users/Emil/10X/scretina/sunyeOIRWT.h5ad")


# In[79]:


adata1


# In[94]:


adata1.obs


# In[81]:


sc.pl.umap(adata1,color=['seurat_clusters'])


# In[77]:


sunyenormal_obs = pd.read_csv(r"//mnt//c//Users//Emil//10X//scretina//sunyenormal_metadata.csv")


# In[78]:


sunyeOIRKO_obs = pd.read_csv(r"//mnt//c//Users//Emil//10X//scretina//sunyeOIRKO_metadata.csv")


# In[79]:


sunyeOIRWT_obs = pd.read_csv(r"//mnt//c//Users//Emil//10X//scretina//sunyeOIRWT_metadata.csv")


# In[80]:


sunyenormal_obs


# In[92]:


sunyeOIRWT_obs = sunyeOIRWT_obs.rename(columns = {'Unnamed: 0':'CellID'})


# In[93]:


sunyeOIRWT_obs = sunyeOIRWT_obs.set_index('CellID')


# In[88]:


sunyeOIRKO_obs


# In[94]:


adata3.obs['labels1'] = sunyeOIRWT_obs['labels']


# In[84]:


sunyenormal_obs['labels']


# In[95]:


adata3.obs


# In[28]:


import scvelo as scv


# In[109]:


WTOIR = scv.read("/mnt/c/Users/Emil/10X/sunye/WTOIR/dataset1.loom", cache = True)


# In[110]:


KOOIR = scv.read("/mnt/c/Users/Emil/10X/sunye/KOOIR/dataset1.loom", cache = True)


# In[111]:


normal = scv.read("/mnt/c/Users/Emil/10X/sunye/normal/dataset1.loom", cache = True)


# In[37]:


adatamerged.obs


# In[195]:


adataWTOIRsubsettest = adatamergedsubset[adatamergedsubset.obs['group'].isin(['OIR_WT'])]


# In[153]:


adataKOOIRsubset = adatamergedsubset[adatamergedsubset.obs['group'].isin(['OIR_KO'])]


# In[150]:


adatanormalsubset = adatamergedsubset[adatamergedsubset.obs['group'].isin(['Normal'])]


# In[156]:


adataWTOIRsubset


# In[40]:


adataWTOIR = scv.utils.merge(adataWTOIR, WTOIR)


# In[41]:


scv.pl.proportions(adataWTOIR)


# In[155]:


adataWTOIRsubset = scv.utils.merge(adataWTOIRsubset, WTOIR)


# In[157]:


adataKOOIRsubset = scv.utils.merge(adataKOOIRsubset, KOOIR)


# In[108]:


adatanormal


# In[158]:


adatanormalsubset = scv.utils.merge(adatanormalsubset, normal)


# In[59]:


WTOIR = scv.utils.merge(WTOIR, KOOIR, normal)


# In[136]:


adatamergedvelo = adataWTOIRsubset.concatenate(adataKOOIRsubset,adatanormalsubset)


# In[181]:


adatamerged = adataWTOIR.concatenate(adataKOOIR,adatanormal)


# In[135]:


adatamergedvelo


# In[287]:


sc.set_figure_params(figsize = [6,6])
sc.pl.draw_graph(adatamerged,color=['labels1','Spp1'],vmax = 4,
                 legend_loc = 'on data', edges = True)


# In[209]:


sc.set_figure_params(figsize = [8,8])
sc.pl.umap(adatamerged,color='labels1',
                 legend_loc = 'on data')


# In[210]:


sc.set_figure_params(figsize = [8,8])
sc.pl.tsne(adatamerged,color='labels1',
                 legend_loc = 'on data')


# In[288]:


from adjustText import adjust_text


# In[291]:


adatamerged.obsm


# In[301]:


def gen_mpl_labels(
    adata, groupby, exclude=(), ax=None, adjust_kwargs=None, text_kwargs=None
):
    if adjust_kwargs is None:
        adjust_kwargs = {"text_from_points": False}
    if text_kwargs is None:
        text_kwargs = {}

    medians = {}

    for g, g_idx in adata.obs.groupby(groupby).groups.items():
        if g in exclude:
            continue
        medians[g] = np.median(adata[g_idx].obsm["X_draw_graph_fa"], axis=0)

    if ax is None:
        texts = [
            plt.text(x=x, y=y, s=k, **text_kwargs) for k, (x, y) in medians.items()
        ]
    else:
        texts = [ax.text(x=x, y=y, s=k, **text_kwargs) for k, (x, y) in medians.items()]

    adjust_text(texts, **adjust_kwargs)

with plt.rc_context({"figure.figsize": (6, 6), "figure.dpi": 300, "figure.frameon": False}):
    ax = sc.pl.draw_graph(adatamerged,color='labels1', show=False,edges =  True, legend_loc='None', frameon=True)
    gen_mpl_labels(
        adatamerged,
        "labels1",
        exclude=("None",),  # This was before we had the `nan` behaviour
        ax=ax,
        adjust_kwargs=dict(arrowprops=dict(arrowstyle='-', color='black')),
        text_kwargs=dict(fontsize=12),
    )
    fig = ax.get_figure()
    fig.tight_layout()
    plt.show()


# In[214]:


sc.set_figure_params(figsize = [6,6])
sc.pl.draw_graph(adataWTOIR,color=['labels1','Spp1'],vmax = 4,
                 legend_loc = 'on data')


# In[215]:


sc.set_figure_params(figsize = [6,6])
sc.pl.draw_graph(adataKOOIR,color=['labels1','Spp1'],vmax = 4,
                 legend_loc = 'on data')


# In[216]:


sc.set_figure_params(figsize = [6,6])
sc.pl.draw_graph(adatanormal,color=['labels1','Spp1'],vmax = 4,
                 legend_loc = 'on data')


# In[143]:


adatamergedsubset = adatamerged[adatamerged.obs['labels1'].isin(['macrophage_2','microglia_1','microglia_4','macrophage_1','microglia_3','microglia_2'])]


# In[213]:


sc.set_figure_params(figsize = [6,6])
sc.pl.draw_graph(adatamergedsubset,color=['labels1','Spp1'], vmax = 4, legend_loc = 'on data')


# In[48]:


sc.set_figure_params(figsize = [6,6])
sc.pl.draw_graph(adatamergedsubset,color=['labels1','Spp1'], vmax = 4, legend_loc = 'on data')


# In[217]:


sc.set_figure_params(figsize = [6,6])
sc.pl.draw_graph(adataWTOIRsubset,color=['labels1','Spp1'],vmax = 4,
                 legend_loc = 'on data')


# In[218]:


sc.set_figure_params(figsize = [6,6])
sc.pl.draw_graph(adataKOOIRsubset,color=['labels1','Spp1'],vmax = 4,
                 legend_loc = 'on data')


# In[219]:


sc.set_figure_params(figsize = [6,6])
sc.pl.draw_graph(adatanormalsubset,color=['labels1','Spp1'],vmax = 4,
                 legend_loc = 'on data')


# In[177]:


AnnData.write_h5ad(adatanormal, "/mnt/c/Users/Emil/10X/sunye/adatanormal.h5ad")


# In[178]:


AnnData.write_h5ad(adatanormalsubset, "/mnt/c/Users/Emil/10X/sunye/adatanormalsubset.h5ad")


# In[179]:


AnnData.write_h5ad(adataWTOIR, "/mnt/c/Users/Emil/10X/sunye/adataWTOIR.h5ad")


# In[196]:


AnnData.write_h5ad(adataWTOIRsubsettest, "/mnt/c/Users/Emil/10X/sunye/adataWTOIRsubsettest.h5ad")


# In[181]:


AnnData.write_h5ad(adataKOOIR, "/mnt/c/Users/Emil/10X/sunye/adataKOOIR.h5ad")


# In[53]:


adatanormal = sc.read_h5ad("/mnt/c/Users/Emil/10X/sunye/adatanormal1.h5ad")
adatanormalsubset = sc.read_h5ad("/mnt/c/Users/Emil/10X/sunye/adatanormalsubset.h5ad")
adataWTOIR = sc.read_h5ad("/mnt/c/Users/Emil/10X/sunye/adataWTOIR.h5ad")
adataWTOIRsubset = sc.read_h5ad("/mnt/c/Users/Emil/10X/sunye/adataWTOIRsubset.h5ad")
adataKOOIR = sc.read_h5ad("/mnt/c/Users/Emil/10X/sunye/adataKOOIR.h5ad")
adataKOOIRsubset = sc.read_h5ad("/mnt/c/Users/Emil/10X/sunye/adataKOOIRsubset.h5ad")


# In[107]:


adatanormal


# In[194]:


AnnData.write_csvs(adataKOOIRsubset, "/mnt/c/Users/Emil/10X/sunye/adatanormal.h5ad")


# In[192]:


adataKOOIR


# In[206]:


scv.tl.velocity_pseudotime(adataWTOIRsubset)
scv.pl.scatter(adataWTOIRsubset, color='velocity_pseudotime',basis='draw_graph_fa', cmap='gnuplot')


# In[228]:


pip list


# In[ ]:




