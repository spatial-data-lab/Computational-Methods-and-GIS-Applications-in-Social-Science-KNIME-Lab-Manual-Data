import knime.scripting.io as knio
import numpy as np
import pandas as pd

A = knio.input_tables[0].to_pandas()
B =knio.input_tables[1].to_pandas()

# convert dataframes to matrices
T= A.values
G= B.values

# define N
N = len(A)

# calculate IAB
IGT = np.identity(N) - np.dot(G,T )

# convert IAB to dataframe
df = pd.DataFrame(IGT)

knio.output_tables[0] = knio.Table.from_pandas(df)
