import scipy as sc
import numpy as np
import pandas as pd
import plotly.express as px

import plotly.graph_objects as go
from scipy.optimize import curve_fit
from .fit import gauss

def gauss_fit_distribution(df,variable,color,color_map,xlim=(-100,100),acc=100):
    plot_array = []; fit_labels = ["AMP", "MEAN", "STD"]
    for idx,filter in enumerate(df[color].unique()):
        filter_df = df[df[color]==filter].copy()
        hist, bin_edges = np.histogram(filter_df[variable], bins=acc, range=xlim)
        bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2
        popt, pcov = curve_fit(gauss, bin_centers, hist, p0=[np.max(hist), 0, 10], bounds=([0, -np.inf, 0], [np.inf, np.inf, np.inf]))
        perr = np.sqrt(np.diag(pcov))
        print("----------------"+filter+"----------------")
        for i in range(len(popt)):
            print(fit_labels[i]+'\t= %.2E\t± %.2E' % (popt[i], perr[i]))
        print("-------------------------------------")
        fit = gauss(bin_centers, *popt)

        plot_dict = {"filter": filter, "hist":hist, "bin_centers":bin_centers, "fit": fit}
        plot_array.append(plot_dict)

    plot_df = pd.DataFrame(plot_array).explode(["hist","bin_centers","fit"])
    fig = px.bar(plot_df,
                x="bin_centers",
                y="hist", 
                color="filter", 
                color_discrete_sequence=pd.Series(plot_df["filter"].unique()).map(color_map),
                barmode="overlay", 
                template="presentation",
                )
    for idx,filter in enumerate(df[color].unique()):
        fig.add_trace(go.Scatter(name="FIT "+filter,x=plot_df[plot_df["filter"] == filter]["bin_centers"], y=plot_df[plot_df["filter"] == filter]["fit"],line=dict(color=color_map[filter])))
        # fig.add_vline(x=np.mean(df[df[color] == filter][variable]))

    return fig