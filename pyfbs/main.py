import data
import plotting
import stability_curve as scc

import numpy as np
import matplotlib.pyplot as plt  # for testing purposes. remove later

if __name__ == "__main__":
    # -------------------------------------------
    # Set name of input files:
    filename1 = "../plots/NbNf_test1.txt"
    filename2 = "../plots/tlncurve_mu2-lambda-200.txt"
    filename3 = "../plots/tlncurve_test3dd2_no_re-computing-NSvalues.txt"
    filename4 = "../plots/twofluid_MRphi-diagram3.txt"
    df, indices = data.load_file(filename4)

    # print(df)

    # -------------------------------------------
    # Compute the stability curve:
    stab_curve = scc.calc_stability_curve(df, indices, debug=False)

    # -------------------------------------------
    # Set labels, colour and corresponding levels and contour lines for 2D Mass-plot:
    contourlines_colour = ["purple", "brown", "red", "green", "orange"]
    # contourlines_colour = ['orange', 'yellow', 'red', 'green', 'purple']
    contourlines_levels = [0.62, 1.2, 1.6, 2.0, 2.3]
    # contourlines_levels = [10, 200, 700, 1000, 3000]
    plot_colorbar_label = "Total Gravitational Mass [M$_\odot$]"
    # plot_colorbar_label = '$\Lambda_{tidal}$'
    # plot_colorbar_label = 'Dark Matter fraction'

    # plot the rho-phi-diagram
    plotting.plotRhoPhi(
        np.clip(df, 0.0, 5000),
        indices,
        "M_T",
        stabCurve=stab_curve,
        contourColors=contourlines_colour,
        contourLevels=contourlines_levels,
        cbarLabel=plot_colorbar_label,
    )

    # -------------------------------------------
    # Set labels, colour and corresponding levels for MR-plot:
    my_plot_filename = "plots/tldtest_mu2-lambda-200_MR_plottest.png"
    my_plot_filename = "plots/effective_EOS-Dm-fraction2.png"
    plot_title = "stab region: DD2 EOS mu=2 lamb=-200 tln solution"
    plot_title = "DD2 EOS + effective boson EOS $\mu=1, \lambda=300$"
    plot_colorbar_label = "log10 $\Lambda_{tidal}$"
    plot_colorbar_label = "Dark Matter fraction"
    # sets the number of bins in the colorbar. Can be int or array. if it is an array, all bins have to be set manually in ascending order
    colorbar_levels = 100

    # filter the stars inside the stability region
    filtered_data_arr = scc.filter_stab_curve_data(df, indices, stab_curve)

    # plot the stars in an MR-diagram with a specified quantity as values in the colorbar
    plotting.plot_interpolate_stability_region(
        np.clip(filtered_data_arr, 0.0, 5000.0),
        indices,
        "DM-fraction",
        20.0,
        2.5,
        0.55,
        False,
        plot_title,
        my_plot_filename,
        plot_colorbar_label,
        False,
        colorbar_levels,
    )

    # -------------------------------------------
    # Obtain the stability curve as a function of Mass and Radius:
    my_stabcurveMR = scc.stab_curve_to_MR(df, indices, stab_curve)
    print(my_stabcurveMR)

    ##----------
    # (for testing purposes. remove later):
    # plt.plot(my_stabcurveMR[:,1], my_stabcurveMR[:,0], c='black', linewidth=2)
    plt.show()
    ##----------

    exit()
