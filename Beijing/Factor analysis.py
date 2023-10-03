import knime.scripting.io as knio

import numpy as np
import pandas as pd

df = knio.input_tables[0].to_pandas()
ref_id_col=pd.Series(df['ref_id'].tolist(),name='ref_id')
df_final = df.drop(['ref_id'],axis=1)

 
from factor_analyzer import FactorAnalyzer
fa = FactorAnalyzer(4, rotation='varimax',method='principal')
fa.fit(df_final)
ev, v = fa.get_eigenvalues() 
var = fa.get_factor_variance()
 

# factor loading
df_loading = pd.DataFrame(fa.loadings_,index=df_final.columns.tolist(),columns=['factor1','factor2','factor3','factor4'])
#output_table_1
knio.output_tables[0]  = knio.Table.from_pandas(df_loading)


# factor score
factscore = pd.DataFrame(fa.transform(df_final), columns=['factor1','factor2','factor3','factor4'])
df=pd.concat([ref_id_col,factscore],axis=1)

#Factor scores are first weighted by their relative importance measured by variance portions accounted for (based on FA) 
df['factor1']=0.3516*df['factor1']
df['factor2']=0.1542*df['factor1']
df['factor3']=0.1057*df['factor1']
df['factor4']=0.0922*df['factor1']

#output_table_2
knio.output_tables[1] = knio.Table.from_pandas(df)



 