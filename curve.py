#! /usr/bin/env python

import matplotlib
# matplotlib.use('Agg')
import matplotlib.pyplot as plt

import sys
import numpy as np
import pandas as pd

from tqdm import tqdm
from itertools import islice
from sklearn.metrics import r2_score
from scipy.optimize import curve_fit

import uno_data as ud


HS_BOUNDS_ORIG = ([0, 10**-12, 0], [1, 1, 4])

def hs_response_curve_original(x, einf, ec50, hs):
    """ from PharmacoDB supp. https://doi.org/10.1093/nar/gkx911
        bounds:
          einf: [0, 1]       # fraction of cells not susceptible to drug
          ec50: [10^-12, 1]  # concentration to have half target receptors bound: [1pM, 1M]
          hs:   [0, 4]       # hill slope binding cooperativity
    """
    return einf + (1 - einf) / (1 + np.power(x/ec50, hs))


HS_BOUNDS = ([0, 0, 0], [1, 12, 4])

def response_curve(x, einf, ec50, hs):
    """ transformed the original function with ec50 in -log10(M) instead of M
    """
    return einf + (1 - einf) / (1 + 10 ** ((ec50 - x) * hs))


def response_integral(x, einf, ec50, hs):
    return (1 - einf) * np.log10(1 + 10 ** ((ec50 - x) * hs)) / hs + x


def compute_area(x1, x2, einf, ec50, hs, mode='trapz'):
    popt = (einf, ec50, hs)
    if mode == 'trapz':
        # trapezoidal numerical integrationcollapse
        xx = np.linspace(x1, x2, 100)
        yy = response_curve(xx, *popt)
        area = np.trapz(yy, xx, dx=0.01)
    else:
        # the integral function can be expressed analytically
        # but sometimes less accurate due to float precision issues
        area = response_integral(x2, *popt) - response_integral(x1, *popt)
    return area


def compute_fit_metrics(xdata, ydata, popt, pcov, d1=4, d2=10):
    if popt is None:
        cols = 'AUC IC50 EC50 EC50se R2fit Einf HS AAC1 AUC1 DSS1'.split(' ')
        return pd.Series([np.nan] * len(cols), index=cols)

    einf, ec50, hs = popt
    perr = np.sqrt(np.diag(pcov))
    ec50se = perr[1]

    xmin = xdata.min()
    xmax = xdata.max()

    ypred = response_curve(xdata, *popt)
    r2 = r2_score(ydata, ypred)

    auc1 = compute_area(xmin, xmax, *popt) / (xmax - xmin)
    aac1 = 1 - auc1

    ic50 = ec50 - np.log10(0.5/(0.5-einf)) / hs if einf < 0.5 else np.nan
    ic90 = ec50 - np.log10(0.9/(0.1-einf)) / hs if einf < 0.1 else np.nan
    ic10 = ec50 - np.log10(0.1/(0.9-einf)) / hs if einf < 0.9 else np.nan

    ic10x = min(ic10, xmax)
    int10x = compute_area(xmin, ic10x, *popt)
    dss1 = (0.9 * (ic10x - xmin) - int10x) / (0.9 * (xmax - xmin)) if xmin < ic10x else 0
    auc = (response_integral(d2, *popt) - response_integral(d1, *popt)) / (d2 - d1)

    metrics = pd.Series({'AUC':auc, 'IC50':ic50, 'EC50':ec50,
                         'EC50se':ec50se, 'R2fit':r2, 'Einf':einf, 'HS':hs,
                         'AAC1':aac1, 'AUC1':auc1, 'DSS1':dss1}).round(4)

    return metrics


def response_curve_fit(xdata, ydata, bounds=HS_BOUNDS):
    ydata = ydata.clip(lower=0, upper=1.0)
    popt, pcov = None, None
    nfev = 100 * 3
    while popt is None and nfev < 10000:
        # print(nfev)
        try:
            popt, pcov = curve_fit(response_curve, xdata, ydata, bounds=bounds, max_nfev=nfev)
            # popt, pcov = curve_fit(response_curve, xdata, ydata, bounds=bounds, max_nfev=nfev, method='dogbox')
        except RuntimeError:
            pass
        nfev *= 2
    return popt, pcov


def fit_exp(df_exp, title=None, dmin=None, dmax=None):
    print(df_exp)
    xdata = df_exp.DOSE.astype(np.float)
    ydata = df_exp.GROWTH.astype(np.float)
    # ydata = df_exp.GROWTH.clip(lower=0, upper=1.0).astype(np.float)

    # print(xdata)
    # print(ydata)

    popt, pcov = response_curve_fit(xdata, ydata)
    metrics = compute_fit_metrics(xdata, ydata, popt, pcov)

    if popt is None:
        return metrics

    dmin = dmin or xdata.min()
    dmax = dmax or xdata.max()
    xx = np.linspace(dmin, dmax, 100)
    yy = response_curve(xx, *popt)

    plt.xlim(dmax, dmin)
    plt.ylim(0, np.max([105, np.max(yy)]))
    plt.plot(xx, yy*100, 'r-', label='fit: Einf=%.3f, EC50=%.3f, HS=%.3f' % tuple(popt))
    plt.plot(xdata, ydata.clip(lower=0, upper=1.0)*100, 'b*', label='')
    plt.xlabel('Dose (-log10(M))')
    plt.ylabel('Growth%')
    plt.title(title)
    plt.legend()
    plt.show()

    return metrics.to_frame(name='metrics').T


def fit_response(df_all, cell, drug, source, study=None):
    cell_ids = ud.cell_name_to_ids(cell) or [cell]
    drug_ids = ud.drug_name_to_ids(drug) or [drug]

    df_exp = df_all[df_all.CELL.isin(cell_ids) & df_all.DRUG.isin(drug_ids)].copy()
    df_exp.GROWTH = (df_exp.GROWTH/2 + 0.5)
    df_exp = df_exp[df_exp.SOURCE == source]

    title = f'{cell} treated with {drug} in {source}'

    studies = df_exp.STUDY.unique()
    if len(studies) > 1:
        study = studies[study] if type(study) == int else study or studies[0]
        title += f' study {study}'
        df_exp = df_exp[df_exp.STUDY == study]

    return fit_exp(df_exp, title)


def show_dose_distribution(df_all):
    sources = df_all.SOURCE.unique()
    qs = [0, 0.02, 0.05, 0.1, 0.2, 0.5, 0.8, 0.9, 0.95, 0.98, 1]
    series = []
    for src in sources:
        s = df_all[df_all.SOURCE == src].DOSE.quantile(qs)
        s.name = src
        series.append(s)
    df_dose = pd.concat(series, axis=1)
    return df_dose


def process_df(df, fname, sep='\t', ngroups=None):
    # df = df1.copy()
    i = 0
    header = None
    cols = ['SOURCE', 'CELL', 'DRUG', 'STUDY']
    groups = df.groupby(cols)
    f = open(fname, 'w')
    for name, group in tqdm(groups):
        # print(name)
        xdata = group.DOSE.astype(np.float)
        ydata = group.GROWTH.clip(lower=0, upper=1.0).astype(np.float)
        popt, pcov = response_curve_fit(xdata, ydata)
        metrics = compute_fit_metrics(xdata, ydata, popt, pcov)
        if header is None:
            header = cols + metrics.index.tolist()
            print(sep.join(header), file=f)
        print(sep.join(name), end=sep, file=f)
        print(sep.join([f'{x:.4g}' for x in metrics]), file=f)
        i += 1
        if ngroups and i >= ngroups:
            break
    f.close()


def process_df_part(df, fname, sep='\t', start=0, count=None):
    header = None
    cols = ['SOURCE', 'CELL', 'DRUG', 'STUDY']
    groups = df.groupby(cols)
    # count = count or (len(groups) - start)
    count = count or (4484081 - start)
    groups = islice(groups, start, start+count)
    f = open(f'{fname}.{start}', 'w')
    for name, group in tqdm(groups):
        # print(name)
        xdata = group.DOSE.astype(np.float)
        ydata = group.GROWTH.clip(lower=0, upper=1.0).astype(np.float)
        popt, pcov = response_curve_fit(xdata, ydata)
        metrics = compute_fit_metrics(xdata, ydata, popt, pcov)
        if start == 0 and header is None:
            header = cols + metrics.index.tolist()
            print(sep.join(header), file=f)
        print(sep.join(name), end=sep, file=f)
        print(sep.join([f'{x:.4g}' for x in metrics]), file=f)
    f.close()


def test():
    df0 = ud.load_single_dose_response(fraction=True)

    cell_ids = ud.cell_name_to_ids('LOXIMVI')
    drug_ids = ud.drug_name_to_ids('paclitaxel')

    df1 = df0[df0.CELL.isin(cell_ids) & df0.DRUG.isin(drug_ids)].copy()
    df1.GROWTH = df1.GROWTH/2 + 0.5
    df2 = df1[df1.SOURCE == 'NCI60']

    fit_exp(df2)


def process_chem_partner_data():
    df_cp = pd.read_csv('curve/ChemPartner_dose_response', sep='\t')
    df_cp = df_cp[df_cp.DRUG2.isnull() & df_cp.DOSE2.isnull()].drop(['DRUG2', 'DOSE2'], axis=1)
    df_cp = df_cp.rename(columns={'DRUG1':'DRUG', 'DOSE1':'DOSE'})
    df_cp.DOSE = -df_cp.DOSE
    # df_cp.GROWTH = df_cp.GROWTH/100
    df_cp.GROWTH = df_cp.GROWTH/200 + 0.5

    # process_df(df_cp, 'curve/ChemPartner_single_response_agg', ngroups=10)
    process_df(df_cp, 'curve/ChemPartner_single_response_agg.new')


def fix_auc_gt_one():
    dfx = pd.read_table('curve/combined_single_response_agg.0', engine='c', low_memory=False)


def notebook():
    d = pd.read_csv('./curve/combined_single_response_agg', engine='c', sep='\t', low_memory=False)
    m = 1e-3
    print(d[(d.AUC < m) & (d.R2fit < m) & (d.EC50se < m)].shape)
    print(d[(d.AUC < m) & (d.R2fit < m) & (d.EC50se < m)].head())


def main():
    df_all = ud.load_single_dose_response(fraction=True)
    df_all.GROWTH = df_all.GROWTH/2 + 0.5
    process_df(df_all, 'curve/combined_single_response_agg')
    # 4484081
    # process_df_part(df_all, 'curve/combined_single_response_agg.debug', start=0, count=270)
    # process_df_part(df_all, 'curve/combined_single_response_agg', start=int(sys.argv[1]), count=int(sys.argv[2]))


if __name__ == '__main__':
    main()
