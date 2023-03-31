import numpy as np
from scipy.stats import rankdata
from scipy.stats import friedmanchisquare
from scipy.stats import chi2
import matplotlib.pyplot as plt
import pandas as pd
filename = 'E:/Github/AEP-anaylsis/code-Python/audiogramdataonesheet_2.xlsx'
df = pd.read_excel(filename)

df_MLR_org =  df.loc[(df['Type'] == 'LLR') & (df['stim'] == 'BMLD')]
# df_LLR =  df.loc[(df['Type'] == 'LLR') & (df['stim'] == 'BMLD')]
# df_ABR =  df.loc[(df['Type'] == 'ABR') & (df['stim'] == 'BMLD')]
df_MLR = df_MLR_org.drop(columns=['Subjects','Type', 'stim'])
print(df_MLR.head)




# Define the data as a pandas DataFrame
# data = pd.DataFrame({
#     'Group 1': [8, 9, 7, 6, 5],
#     'Group 2': [9, 9, 6, 7, 5],
#     'Group 3': [8, 8, 8, 7, 4],
#     'Group 4': [7, 9, 6, 7, 5],
#     'Group 5': [8, 8, 6, 7, 6]
# })
data=df_MLR
# Rank the data
ranks = rankdata(-data.values, method="average", axis=1)
mean_ranks = np.mean(ranks, axis=0)
print(data)
n, k = data.shape
F, p = friedmanchisquare(*data.values.T)
alpha = 0.05

# Compute the critical difference
cd = np.sqrt(k*(k+1)/(6*n)*chi2.ppf(1-alpha, k))
print(cd)

# Plot the Nemenyi plot
# Plot the Nemenyi plot
plt.errorbar(x=np.arange(k), y=mean_ranks, yerr=cd, fmt="o")
plt.xticks(np.arange(k), data.columns)
plt.ylabel("Mean rank")
plt.xlabel("Frequency")
plt.title("Nemenyi plot signal duration 48ms and CD = 1.25")
plt.show()
print(mean_ranks)



