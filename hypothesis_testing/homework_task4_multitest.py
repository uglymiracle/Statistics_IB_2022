import pandas as pd
import numpy as np
import scipy.stats as st
from statsmodels.stats.weightstats import ztest
import statsmodels.stats.multitest


def check_intervals_intersect(first_ci, second_ci):   
    
    first_intervals = st.t.interval(alpha=0.95, 
                                    df=len(first_ci) - 1, 
                                    loc=np.mean(first_ci), 
                                    scale=st.sem(first_ci))

    second_intervals = st.t.interval(alpha=0.95, 
                                    df=len(second_ci) - 1, 
                                    loc=np.mean(second_ci), 
                                    scale=st.sem(second_ci))
    
    overlap = max(0, min(first_intervals[1], second_intervals[1]) 
                  - max(first_intervals[0], second_intervals[0]))
    
    return overlap != 0


def check_dge_with_ci(first_table, second_table):
    
    ci_test_results = []
    genes = []
    
    for i in first_table:
    
        if i == 'Cell_type':
            continue
        
        overlap = check_intervals_intersect(first_table[i], second_table[i])
        ci_test_results.append(not overlap)
        genes.append(i)
        
    table = {'gene': genes,
             'ci_test_results': ci_test_results}
    ci_test_results = pd.DataFrame(table)

    return ci_test_results


def check_dge_with_ztest(first_table, second_table):
    
    z_test_results = []
    genes = []
    z_test_p_values = []
    
    for i in first_table:
        if i == 'Cell_type':
            continue
        p_val = ztest(first_table[i], second_table[i])[1]

        if p_val < 0.05:
            z_test_result = True
        else:
            z_test_result = False
    
        z_test_results.append(z_test_result)
        z_test_p_values.append(p_val)
        genes.append(i)
    
    table = {'gene': genes,
             'z_test_results': z_test_results,
            'z_test_p_values': z_test_p_values}
    z_test_results = pd.DataFrame(table)
    return z_test_results


def mean_diff(first_table, second_table):
    
    mean_diffs = []
    genes = []
    
    for i in first_table:
        if i == 'Cell_type':
            continue
            
        mean_d = np.mean(first_table[i]) - np.mean(second_table[i])
        mean_diffs.append(mean_d)
        genes.append(i)
    
    table = {'gene': genes,
             'mean_diff': mean_diffs}
    mean_diff_results = pd.DataFrame(table)
    
    return mean_diff_results


def multitest(results_table, method, alpha):

    
    pvals = results_table['z_test_p_values'].values.tolist()
    multitest_table = statsmodels.stats.multitest.multipletests(pvals, alpha=alpha, method=method)

    z_test_adj_results = multitest_table[0]
    pvals_adjusted = multitest_table[1]
    
    pvals_adj_table = {'z_test_p_values_adjusted': pvals_adjusted,
                       'z_test_adjusted_result': z_test_adj_results}
    table = pd.DataFrame(pvals_adj_table)
    
    return table


def main():
    
    first_cell_type_expressions_path = input('Path to first table: ') 
    second_cell_type_expressions_path = input('Path to second table: ') 
    save_results_table = input('Name of results table: ')
    method = input('Method used for testing and adjustment of pvalues: ')
    if method != '':
        alpha = float(input('Alpha: '))

    first_table = pd.read_csv(first_cell_type_expressions_path)
    second_table = pd.read_csv(second_cell_type_expressions_path)   
    
    ci_table = check_dge_with_ci(first_table, second_table)
    ztest_table = check_dge_with_ztest(first_table, second_table)
    mean_diff_table = mean_diff(first_table, second_table)

    results = pd.merge(ci_table, ztest_table, on=('gene'))
    results = pd.merge(results, mean_diff_table, on=('gene'))
    results = results.iloc[1: , : ]

    if method != '':

        adj_pvals_table = multitest(results, method, alpha)
        adj_pvals_table = adj_pvals_table.iloc[1: , : ]
        results = pd.concat([results, adj_pvals_table], axis=1)
        results = results.reindex(columns=['gene',
                                           'ci_test_results',
                                           'z_test_results',
                                           'z_test_p_values',
                                           'z_test_p_values_adjusted',
                                           'z_test_adjusted_result',
                                           'mean_diff'])
    
    results.to_csv(save_results_table, index=False)


if __name__ == '__main__':
    main()
