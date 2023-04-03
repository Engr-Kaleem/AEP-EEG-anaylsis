import numpy as np
from scipy.stats import rankdata
from scipy.stats import friedmanchisquare
from scipy.stats import chi2
import matplotlib.pyplot as plt
import pandas as pd
filename = 'E:/Github/AEP-anaylsis/code-Python/audiogramdataonesheet_2.xlsx'
df = pd.read_excel(filename)

df_MLR_org =  df.loc[(df['Type'] == 'MLR') & (df['stim'] == 'BMLD')]
# df_LLR =  df.loc[(df['Type'] == 'LLR') & (df['stim'] == 'BMLD')]
# df_ABR =  df.loc[(df['Type'] == 'ABR') & (df['stim'] == 'BMLD')]
df_MLR = df_MLR_org.drop(columns=['Subjects','Type', 'stim'])
print(df_MLR.head)
# freq_125 = df['125'].values
# freq_250 = df['250'].values
# freq_500 = df['500'].values
# freq_750 = df['750'].values
# freq_1000 = df['1000'].values

# Run the Friedman rank test
# data = [freq_125, freq_250, freq_500, freq_750, freq_1000]
# friedman_stat, p_value = friedmanchisquare(*data)

# # Print the results
# print('Friedman test statistic:', friedman_stat)
# print('p-value:', p_value)



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
print(n,k)
F, p = friedmanchisquare(*data.values.T)
alpha = 0.05

# Compute the critical difference
cd = np.sqrt(k*(k+1)/(6*n)*chi2.ppf(1-alpha, k))
print(cd)

# Plot the Nemenyi plot
# Plot the Nemenyi plot
plt.errorbar(x=np.arange(k), y=mean_ranks, yerr=cd, fmt="o")
plt.plot(np.arange(k), mean_ranks)
plt.xticks(np.arange(k), data.columns)
plt.ylabel("Mean rank")
plt.xlabel("Frequency")
plt.title(f"Nemenyi plot signal duration 3ms and CD ={cd}")

plt.show()
print(mean_ranks)



