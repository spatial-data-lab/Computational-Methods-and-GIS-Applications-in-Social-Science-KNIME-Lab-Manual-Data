import knime.scripting.io as knio
import numpy as np
import pandas as pd

A = knio.input_tables[0].to_pandas()
B =knio.input_tables[1].to_pandas()

# convert dataframes to matrices
G= A.values
B= B.values
# calculate IAB
GB = np.dot(G, B)

# convert IAB to dataframe
df = pd.DataFrame(GB)

knio.output_tables[0] = knio.Table.from_pandas(df)
