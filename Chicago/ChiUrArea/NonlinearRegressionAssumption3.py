import knime.scripting.io as knio
import geopandas as gp
import pandas as pd
import numpy as np
from scipy.optimize import curve_fit

df = knio.input_tables[0].to_pandas()


def func_Assum3(r, a1, a2, b1,b2):
    return a1*pow(np.exp(1), b1*r[:,0]) + a2*pow(np.exp(1), b2*r[:,1])

X = df[['15','5']].values
popt, pcov = curve_fit(func_Assum3, X, df.popden.values,bounds = ([-5000,-10000,-0.5,-0.5],[15000,10000,0,0]))




def get_indexes(y_predict, y_data):
    n = y_data.size
    SSE = ((y_data - y_predict)**2).sum()
    MSE = SSE / n
    RMSE = np.sqrt(MSE)
    u = y_data.mean()
    SST = ((y_data - u)**2).sum()
    SSR = SST - SSE
    R_square = SSR / SST
    return SSE, MSE, RMSE, R_square
    
y_predict = func_Assum3(X, *popt)
indexes = get_indexes(y_predict, df['popden'])
R2 = 1-indexes[3]

# Copy input to output
out_dict = {'param_name':['a15','a5','b15','b5','R2'],'Value':[popt[0],popt[1],popt[2],popt[3],R2]}
df1 = pd.DataFrame(out_dict)


df1.reset_index(drop=True, inplace=True)
knio.output_tables[0] = knio.Table.from_pandas(df1)

