import numpy as  np
import seaborn as sns
import pandas as pd 
import matplotlib.pyplot as plt
import scipy.stats as stats
from scipy.stats import binom_test
from statsmodels.stats.contingency_tables import mcnemar
import statsmodels.api as sm
from statsmodels.formula.api import ols
from scipy.stats import friedmanchisquare,rankdata
from matplotlib.backends.backend_pdf import PdfPages
filename = 'E:/Github/AEP-anaylsis/code-Python/audiogramdataonesheet_2.xlsx'



df = pd.read_excel(filename)

df_500_org =  df.loc[(df['stim'] == 'BMLD')]
# df_LLR =  df.loc[(df['Type'] == 'LLR') & (df['stim'] == 'BMLD')]
# df_ABR =  df.loc[(df['Type'] == 'ABR') & (df['stim'] == 'BMLD')]
df_MLR_long = df_500_org.iloc[:, [0,5,7]]




from scipy.stats import f 
# Set the significance level
alpha = 0.05

# Set the number of treatments and degrees of freedom
k = 5
df = k - 1

# Set the ranks of the treatments
ranks = np.arange(1, k+1)

# Set the number of observations
n = 35

# Calculate the critical value of the Studentized range distribution
q = f.ppf(1 - alpha / (k * (k - 1)), k, (k - 1) * (n - 1))

# Calculate the critical difference value
CD = 2.343 * np.sqrt(k * (k+1) / (6 * n))

print("Critical difference value:", CD,q)


df_MLR_long=df_MLR_long.rename(columns={1000: 'BMLD'})
print(df_MLR_long)
df_MLR_long['Type'] = df['Type'].replace({'MLR': '18ms', 'LLR': '48ms','ABR':'3ms'})
# Specify the dependent variable as "BMLD"
dependent_var = 'BMLD'
factor = 'Type'
levels = ['3ms', '18ms', '48ms']



f_value, p_value = friedmanchisquare(*[df_MLR_long[df_MLR_long[factor] == level][dependent_var].values for level in df_MLR_long[factor].unique()])





# Compute the ranks of the BMLD values for each Type level
df_MLR_long['Rank'] = rankdata(-df_MLR_long['BMLD'].values).astype(int)
ranks = np.zeros((len(df_MLR_long['Type'].unique()), len(df_MLR_long)))
for i, freq in enumerate(df_MLR_long['Type'].unique()):
    ranks[i, df_MLR_long['Type'] == freq] = df_MLR_long[df_MLR_long['Type'] == freq]['Rank'].values

# Compute the mean ranks and the rank differences between each Type pair
mean_ranks = np.mean(ranks, axis=1)
pairwise_differences = np.abs(np.subtract.outer(mean_ranks, mean_ranks))

# Compute the critical value for the Nemenyi test
n_groups = len(df_MLR_long['Type'].unique())
n_samples = len(df_MLR_long)
q_alpha = 2.343 * np.sqrt(n_groups * (n_groups + 1) / (6 * n_samples))
print(f' CD={q_alpha}')
# Perform the post-hoc Nemenyi test for all Type pairs and save the results to a DataFrame
results = pd.DataFrame(columns=['Type 1', 'Type 2', 'Difference in mean ranks', 'Z-value', 'p-value', 'Significant'])
index = 0
for i in range(n_groups):
    for j in range(i+1, n_groups):
        diff = mean_ranks[i] - mean_ranks[j]
        q = diff / np.sqrt((n_groups * (n_groups + 1)) / (6 * n_samples))
        p = 1 - np.sum(pairwise_differences > np.abs(diff)) / (n_groups * (n_groups - 1) / 2)
        if np.abs(q) > q_alpha:
            significant = "Yes"
        else:
            significant = "No"
        results.loc[index] = [df_MLR_long['Type'].unique()[i], df_MLR_long['Type'].unique()[j], diff, q, p, significant]
        index += 1

# Save the results DataFrame to a CSV file
results.to_csv('bml_statistics_1000.csv', index=False)

# Print the test result
print(f"Friedman test: F = {f_value:.2f}, p = {p_value:.4f}")

# Print the test result
print(f"Friedman test: F = {f_value}, p = {p_value}")


# Create a boxplot to visualize the BMLD values for each frequency level
f1=plt.figure()
sns.boxplot(x='Type', y='BMLD', data=df_MLR_long)
plt.title('Boxplot for  BMLD values for each duration @1000hz')
plt.xlabel('duration')
plt.ylabel('BMLD')
# plt.show()



f2=plt.figure()
# Create a heatmap to visualize the pairwise differences in mean ranks between each frequency pair
sns.heatmap(pairwise_differences, xticklabels=levels, yticklabels=levels, cmap='coolwarm', annot=True, fmt='.2f')
plt.title('Pairwise differences in mean ranks between each duration pair for 1000hz')
plt.xlabel('duration')
plt.ylabel('duration')
# plt.show()





f3=plt.figure()
# Create a bar chart to visualize the results of the Nemenyi test
results['Type Pair'] = results['Type 1'].astype(str) + '-' + results['Type 2'].astype(str)
sns.barplot(x='Type Pair', y='Difference in mean ranks', data=results, hue='Significant')
plt.title('Mean rank differences between each duration pair for 1000hz')
plt.xlabel('duration pair')
plt.ylabel('Difference in mean ranks')
# plt.show()



def save_image(filename):
    
    # PdfPages is a wrapper around pdf 
    # file so there is no clash and create
    # files with no error.
    p = PdfPages(filename)
      
    # get_fignums Return list of existing 
    # figure numbers
    fig_nums = plt.get_fignums()  
    figs = [plt.figure(n) for n in fig_nums]
      
    # iterating over the numbers in list
    for fig in figs: 
        
        # and saving the files
        fig.savefig(p, format='pdf') 
      
    # close the object
    p.close()  
  
# name your Pdf file
filename = "plots_1000.pdf"  
  
# call the function
save_image(filename) 