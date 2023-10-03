import knime.scripting.io as knio

# This example script simply outputs the node's input table.
import geopandas as gp
import pandas as pd
import numpy as np
# Non-liner
from scipy.optimize import curve_fit
import statsmodels.api as sm
from statsmodels.formula.api import ols, wls

input_table = knio.input_tables[0].to_pandas()

lnDr=input_table['lnpopden']
r=input_table['Distance']
w=input_table['AREA']
# Weighted Regression
df = pd.concat([lnDr, r], axis=1)  

# WLS
weighted_model = wls('lnDr ~ r', weights = w, data = df).fit()
weighted_model_summary = weighted_model.summary2().tables[1]
weighted_model_summary.reset_index(drop=True, inplace=True)
# Copy input to output
 
print(weighted_model.summary())
weighted_model_summary2 = weighted_model.summary2().tables[0]
weighted_model_summary2.reset_index(drop=True, inplace=True)
knio.output_tables[0] = knio.Table.from_pandas(weighted_model_summary)
knio.output_tables[1] = knio.Table.from_pandas(weighted_model_summary2)
