import knime.scripting.io as knio
import numpy as np
import pandas as pd
from scipy.linalg import lu_factor, lu_solve

A = knio.input_tables[0].to_pandas()
B =knio.input_tables[1].to_pandas()

# convert dataframes to matrices
IGT= A.values
PB= B.values 
N = len(IGT)
# 4.1 Factor AA matrix, print out INFO (indicating if the solution exists).
lu, piv = lu_factor(IGT.T)
x = lu_solve((lu, piv), PB )
print(np.allclose(IGT.T @ x - PB, np.zeros((N,))))

# convert matrix to dataframe
df = pd.DataFrame(x)

knio.output_tables[0] = knio.Table.from_pandas(df)
