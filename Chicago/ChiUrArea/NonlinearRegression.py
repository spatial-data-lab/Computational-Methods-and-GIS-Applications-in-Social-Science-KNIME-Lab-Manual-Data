import knime.scripting.io as knio

# This example script simply outputs the node's input table.
import geopandas as gp
import pandas as pd
import numpy as np
# Non-liner
from scipy.optimize import curve_fit

input_table = knio.input_tables[0].to_pandas()

y=input_table['popden']
x=input_table['Distance']

# Dr = a*r^b
def func1(r, a, b):  
    return a * pow(r, b)
# Dr = a*e^(b*r)
def func2(r, a, b):
    return (r * b).apply(np.exp) * a        
# Dr = a*e^(b*r*r)
def func3(r, a, b):
    return a * pow(np.exp(1), (b*r*r))   
 # Dr = a*e^(b*r + c*r^2)
def func4(r, a, b, c):
    return (b*r+c*r*r).apply(np.exp) * a


def nls(func, x, y,p0):
    popt, pcov = curve_fit(func, x, y,p0)
    y_predict = func(x, *popt)
    n = y.size
    SSE = ((y - y_predict)**2).sum()
    MSE = SSE / n
    RMSE = np.sqrt(MSE)
    u = y.mean()
    SST = ((y - u)**2).sum()
    SSR = SST - SSE
    R_square = SSR / SST
    b = ["b" + str(i) for i in range(1,len(popt))]
    namelist=['Intercept']+b+['R2']
    valuelist=list(popt)+[R_square]    
    out_dict = {'param_name':namelist,'Value':valuelist}
    df = pd.DataFrame(out_dict)
    return df

df1=nls(func1, x, y,p0 = [10000,-0.1])
df2=nls(func2, x, y,p0 = [10000,-0.01])
df3=nls(func3, x, y,p0 = [10000,-0.01])
df4=nls(func4, x, y,p0 = [1000, 0.005, -0.005])

df = pd.concat([df1, df2,df3,df4])
df.reset_index(drop=True, inplace=True)
knio.output_tables[0] = knio.Table.from_pandas(df)

