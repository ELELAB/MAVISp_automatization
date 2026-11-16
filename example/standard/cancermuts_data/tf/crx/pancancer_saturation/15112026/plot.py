import pandas as pd
from cancermuts.table import Table

df = pd.read_csv('metatable_pancancer_MLH1.csv')
tbl=Table()
ax, fig = tbl.plot_metatable(df, fname='plotall.pdf', section_size=20, figsize=(15,100))
#ax, fig = tbl.plot_metatable(df, fname='mutations.pdf', section_size=20, figsize=(15,90), mutations=True, elm=False, ptms=False, mutations_revel=True, structure=True, structure_mobidb=False)
#ax, fig = tbl.plot_metatable(df, fname='ptms.pdf', section_size=50, figsize=(15,70), mutations=False, elm=False, ptms=True, mutations_revel=False, structure=True, structure_mobidb=False)
