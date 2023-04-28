from HARK.core import _log
import matplotlib.pyplot as plt
import numpy as np
from copy import deepcopy
import ipywidgets as widgets
import logging
import warnings

from HARK.ConsumptionSaving.ConsIndShockModel \
    import (IndShockConsumerType, init_idiosyncratic_shocks)

warnings.filterwarnings("ignore")

base_params = deepcopy(init_idiosyncratic_shocks)

# Set the parameters for the baseline results in the paper
base_params['PermGroFac'] = [1.03]  # Permanent income growth factor
base_params['Rfree'] = Rfree = 1.04  # Interest factor on assets
base_params['DiscFac'] = DiscFac = 0.96  # Time Preference Factor
base_params['CRRA'] = CRRA = 2.00  # Coefficient of relative risk aversion
# Probability of unemployment (e.g. Probability of Zero Income in the paper)
base_params['UnempPrb'] = UnempPrb = 0.005
base_params['IncUnemp'] = IncUnemp = 0.0   # Induces natural borrowing constraint
base_params['PermShkStd'] = [0.1]   # Standard deviation of log permanent income shocks
base_params['TranShkStd'] = [0.1]   # Standard deviation of log transitory income shocks
# %%
# Uninteresting housekeeping and details
# Make global variables for the things that were lists above -- uninteresting housekeeping
PermGroFac, PermShkStd, TranShkStd = base_params['PermGroFac'][0], base_params['PermShkStd'][0], base_params['TranShkStd'][0]

# Some technical settings that are not interesting for our purposes
base_params['LivPrb'] = [1.0]   # 100 percent probability of living to next period
base_params['CubicBool'] = True    # Use cubic spline interpolation
base_params['BoroCnstArt'] = None    # No artificial borrowing constraint

# Settings to speed up the calcs for the widgets (at the cost of accuracy)
base_params['TranShkCount'] = 2    # 2 shocks instead of 7 speeds things up a lot!
base_params['PermShkCount'] = 2    #
base_params['aXtraCount'] = 20     # not very many gridpoints; solution not too accurate

fssml, fsmid, fsbig = 18, 22, 26

# Define a slider for the discount factor
DiscFac_widget = [
    widgets.FloatSlider(
        min=0.90,
        max=0.99,
        step=0.015,
        value=DiscFac,  # Default value
        continuous_update=False,
        readout_format=".3f",
        description="\u03B2",
    )
    for i in range(5)
]  # beta unicode

DiscFac_growth_widget = [
    widgets.FloatSlider(
        min=0.94,
        max=0.99,
        step=0.01,
        value=0.98,  # Default value
        continuous_update=False,
        readout_format=".3f",
        description="\u03B2",
    )
    for i in range(5)
]  # beta unicode

# Define a slider for relative risk aversion
CRRA_widget = [
    widgets.FloatSlider(
        min=1.1,
        max=5.0,
        step=0.1,
        value=CRRA,  # Default value
        continuous_update=False,
        readout_format=".2f",
        description="\u03C1",
    )
    for i in range(5)
]  # rho unicode

# Define a slider for the interest factor
Rfree_widget = [
    widgets.FloatSlider(
        min=1.01,
        max=1.08,
        step=0.01,
        value=Rfree,  # Default value
        continuous_update=False,
        readout_format=".2f",
        description="R",
    )
    for i in range(5)
]


# Define a slider for permanent income growth
PermGroFac_widget = [
    widgets.FloatSlider(
        min=1.00,
        max=1.03,
        step=0.01,
        value=PermGroFac,  # Default value
        continuous_update=False,
        readout_format=".2f",
        description="\u03D5",
    )
    for i in range(5)
]

# Define a slider for permanent income growth
PermGroFac_growth_widget = [
    widgets.FloatSlider(
        min=1.00,
        max=1.04,
        step=0.005,
        value=1.02,  # Default value
        continuous_update=False,
        readout_format=".3f",
        description="\u03D5",
    )
    for i in range(5)
]

# change default value for GIC fail figure
PermGroFac_widget[1].value = 1.0

# change default value for Bounds fig
Rfree_widget[3].min = 1.0
Rfree_widget[3].max = 1.0
Rfree = 1.0
PermGroFac_widget[3].min = Rfree_widget[3].min-0.03
PermGroFac_widget[3].max = Rfree_widget[3].max
DiscFacMin = 0.92
DiscFacMax = ((Rfree_widget[3].max)**(CRRA-1))/Rfree_widget[3].max - 0.01

# Define a slider for the discount factor
DiscFac_widget = [
    widgets.FloatSlider(
        min=0.92,
        max=DiscFacMax,
        step=0.01,
        value=DiscFac,  # Default value
        continuous_update=False,
        readout_format=".2f",
        description="\u03B2",
    )
    for i in range(5)
]  # beta unicode

# Define a slider for unemployment probability
UnempPrb_widget = [
    widgets.FloatSlider(
        min=0.001,
        max=0.05,  # Go up to twice the default value
        step=0.001,
        value=UnempPrb,
        continuous_update=False,
        readout_format=".3f",
        description="℘",
    )
    for i in range(5)
]

# Define a slider for unemployment income
IncUnemp_widget = [
    widgets.FloatSlider(
        min=0.001,
        max=0.01,  # Go up to twice the default value
        step=0.001,
        value=IncUnemp,
        continuous_update=False,
        readout_format=".3f",
        description="$\\mho$",
    )
    for i in range(5)
]

# Define a slider for PermShkStd
PermShkStd_widget = [
    widgets.FloatSlider(
        min=0.01,
        max=0.30,  # Go up to twice the default value
        step=0.01,
        value=PermShkStd,
        continuous_update=False,
        readout_format=".2f",
        description="$\sigma_\psi$",
    )
    for i in range(5)
]

# Define an alternative slider for PermShkStd
PermShkStd_alt_start_widget = [
    widgets.FloatSlider(
        min=0.10,
        max=0.20,
        step=0.02,
        value=0.2,
        continuous_update=False,
        readout_format=".2f",
        description="$\sigma_\psi$",
    )
    for i in range(5)
]

# Define a slider for the std of the transitory shock
TranShkStd_widget = [
    widgets.FloatSlider(
        min=0.00,
        max=0.30,  # Go up to twice the default value
        step=0.025,
        value=TranShkStd,
        continuous_update=False,
        readout_format=".2f",
        description="$\sigma_θ$",
    )
    for i in range(5)
]


def makeConvergencePlot(DiscFac, CRRA, Rfree, PermShkStd):
    # Construct finite horizon agent with baseline parameters
    baseAgent_Fin = IndShockConsumerType(
        quietly=True, messaging_level=logging.CRITICAL, **base_params)
    baseAgent_Fin.DiscFac = DiscFac
    baseAgent_Fin.CRRA = CRRA
    baseAgent_Fin.Rfree = Rfree
    baseAgent_Fin.PermShkStd = [PermShkStd]
    baseAgent_Fin.cycles = 100
    baseAgent_Fin.solve(quietly=True, messaging_level=logging.CRITICAL)
    baseAgent_Fin.unpack('cFunc')

    # figure plot limits
    mMin, mMax, cPlotMin, cPlotMax = 0, 7, 0, 7

    mPlotMin = 0
    mLocCLabels = 5.6  # Defines horizontal limit of figure
    mPlotTop = 6.5  # 3.5 # 6.5    # Defines maximum m value where functions are plotted
    mPts = 100      # Number of points at which functions are evaluated

    plt.figure(figsize=(12, 8))
    plt.ylim([cPlotMin, cPlotMax])
    plt.xlim([mMin, mMax])

    mBelwLabels = np.linspace(mPlotMin, mLocCLabels-0.1, mPts)  # Range of m below loc of labels
    m_FullRange = np.linspace(mPlotMin, mPlotTop, mPts)        # Full plot range
    # c_Tm0  defines the last period consumption rule (c=m)
    c_Tm0 = m_FullRange
    # c_Tm1 defines the second-to-last period consumption rule
    c_Tm1 = baseAgent_Fin.cFunc[-2](mBelwLabels)
    c_Tm5 = baseAgent_Fin.cFunc[-6](mBelwLabels)  # c_Tm5 defines the T-5 period consumption rule
    # c_Tm10 defines the T-10 period consumption rule
    c_Tm10 = baseAgent_Fin.cFunc[-11](mBelwLabels)
    # c_Limt defines limiting inﬁnite-horizon consumption rule
    c_Limt = baseAgent_Fin.cFunc[0](mBelwLabels)

    plt.plot(mBelwLabels, c_Limt, label="$c(m)$")
    plt.plot(mBelwLabels, c_Tm1, label="$c_{T-1}(m)$")
    plt.plot(mBelwLabels, c_Tm5, label="$c_{T-5}(m)$")
    plt.plot(mBelwLabels, c_Tm10, label="$c_{T-10}(m)$")
    plt.plot(m_FullRange, c_Tm0, label="$c_{T}(m) = 45$ degree line")
    plt.legend(fontsize='x-large')
    plt.tick_params(
        labelbottom=False, labelleft=False, left="off", right="off",
        bottom="off", top="off")

    plt.show()
    return None


def makeGICFailExample(DiscFac, PermShkStd, UnempPrb):
    # replace default params with passed
    GIC_fails_dict = deepcopy(base_params)
    GIC_fails_dict['DiscFac'] = DiscFac
    GIC_fails_dict['PermShkStd'] = [PermShkStd]
    GIC_fails_dict['UnempPrb'] = UnempPrb

    GIC_fails_dict['PermShkCount'] = 7  # Need more accuracy

    GICFailsExample = IndShockConsumerType(
        **GIC_fails_dict, quietly=True, messaging_level=logging.WARNING)

    # Prior command set up the problem but did not solve it
    GICFailsExample.solve(quietly=True, messaging_level=logging.WARNING)

    # Shortcuts/aliases
    soln = GICFailsExample.solution[0]  # solution
    cFunc, Bilt, E_Next_ = soln.cFunc, soln.Bilt, soln.E_Next_
    mBalLvl, mTrgNrm = Bilt.mBalLvl, Bilt.mTrgNrm

    E_d_mtp1_0 = E_Next_.c_where_E_Next_m_tp1_minus_m_t_eq_0
    E_MGro_Bal = E_Next_.c_where_E_Next_PermShk_tp1_times_m_tp1_minus_m_t_eq_0

    mPlotMin, mPlotMax = 0, 20
    cPlotMin, cPlotMax = 0, 1.1 * E_d_mtp1_0(mPlotMax)

    mPtsCount = 50
    mVec = np.linspace(mPlotMin, mPlotMax, mPtsCount)
    c_Limt = cFunc(mVec)
    c_E_d_mtp1_0 = E_d_mtp1_0(mVec)  # "sustainable" consumption ratio
    c_E_MGro_Bal = E_MGro_Bal(mVec)  # "M-Balanced-Growth" consumption ratio

    plt.figure(figsize=(12, 8))
    plt.plot(mVec, c_Limt, label="$c(m_{t})$", color="black")
    plt.plot(mVec, c_E_d_mtp1_0, label="$\mathsf{E}_{t}[\Delta m_{t+1}] = 0$", color="orange")
    plt.plot(mVec, c_E_MGro_Bal, label="$\mathsf{E}_{t}[M_{t+1}/M_{t}] = {\Phi}$", color="green")

    plt.xlim(mPlotMin, mPlotMax)
    plt.ylim(cPlotMin, cPlotMax)

    if Bilt.GICLiv:  # Growth impatience condition (including mortality)
        plt.axvline(Bilt.mBalLvl, label="mBalLvl exists", color="orange", linestyle="dashed")

    if Bilt.GICMod:  # Normalized GIC
        plt.axvline(Bilt.mTrgNrm, label="mTrgNrm exists", color="green", linestyle="dashed")

    plt.tick_params(
        labelbottom=False, labelleft=False, left="on", right="off",  # bottom="off",
        top="off"
    )

    _log.critical(f'mTrgNrm: {mTrgNrm:.3f}')
    _log.critical(f'mBalLvl: {mBalLvl:.3f}')

    plt.legend(fontsize='medium')
    plt.show()
    return None


def cNrmTargetFig_make(PermGroFac, DiscFac):

    gro_params = deepcopy(base_params)
    gro_params['PermGroFac'] = [PermGroFac]
    gro_params['DiscFac'] = DiscFac
    gro_params['aXtraGrid'] = 40
    gro_params['aXtraMax'] = 100

    baseAgent_Inf = IndShockConsumerType(
        **gro_params, quietly=True, messaging_level=logging.WARNING)  # construct it silently

    baseAgent_Inf.tolerance = baseAgent_Inf.tolerance/10

    baseAgent_Inf.solve(
        quietly=True, messaging_level=logging.WARNING)

    soln = baseAgent_Inf.solution[0]        # shorthand
    Bilt, Pars, E_Next_, cFunc = soln.Bilt, soln.Pars, soln.E_Next_, soln.cFunc

    Rfree, DiscFac, CRRA, G = Pars.Rfree, Pars.DiscFac, Pars.CRRA, Pars.PermGroFac

    color_cons, color_mrktLev, color_mrktRat, color_perm = "blue", "red", "green", "orange"

    mPlotMin, mCalcMax, mPlotMax = 0.3, 50, 8

    # Get StE and target values
    mBalLvl, mTrgNrm = Bilt.mBalLvl, Bilt.mTrgNrm

    pts_num = 200  # Plot this many points

    m_pts = np.linspace(1, mPlotMax, pts_num)   # values of m for plot
    c_pts = soln.cFunc(m_pts)                   # values of c for plot
    a_pts = m_pts - c_pts                       # values of a

    Ex_cLvl_tp1_Over_pLvl_t = [
        E_Next_.cLvl_tp1_Over_pLvl_t_from_a_t(a) for a in a_pts]
    Ex_mLvl_tp1_Over_pLvl_t = [
        E_Next_.mLvl_tp1_Over_pLvl_t_from_a_t(a) for a in a_pts]
    Ex_m_tp1_from_a_t = [
        E_Next_.m_tp1_from_a_t(a) for a in a_pts]

    Ex_cLvlGro = np.array(Ex_cLvl_tp1_Over_pLvl_t)/c_pts
    Ex_mLvlGro = np.array(Ex_mLvl_tp1_Over_pLvl_t)/m_pts
    Ex_mNrmGro = np.array(Ex_m_tp1_from_a_t)/m_pts

    # Absolute Patience Factor = lower bound of consumption growth factor
    APFac = (Rfree*DiscFac)**(1.0/CRRA)

    fig, ax = plt.subplots(figsize=(12, 4))

    # Plot the Absolute Patience Factor line
    ax.plot([0, mPlotMax], [APFac, APFac], color=color_cons, linestyle='dashed')

    # Plot the Permanent Income Growth Factor line
    ax.plot([0, mPlotMax], [G, G], color=color_perm)

    # Plot the expected consumption growth factor
    ax.plot(m_pts, Ex_cLvlGro, color=color_cons,
            label=r'$\mathbf{c}$-Level-Growth: $\mathbb{E}_{t}[{\mathbf{c}}_{t+1}/{\mathbf{c}}_{t}]$'
            )

    # Plot the expect growth for the level of market resources
    mLvlGro_lbl, = ax.plot(m_pts, Ex_mLvlGro, color=color_mrktLev,
                           label=r'$\mathbf{m}$-Level-Growth: $\mathbb{E}_{t}[{\mathbf{m}}_{t+1}/{\mathbf{m}}_{t}]$')

    # Plot the expect growth for the market resources ratio
    mNrmGro_lbl, = ax.plot(m_pts, Ex_mNrmGro, color=color_mrktRat,
                           label=r'$m$-ratio Growth: $\mathbb{E}_{t}[m_{t+1}/m_{t}]$')

    # Axes limits
    GroFacMin, GroFacMax, xMin = 0.96, 1.08, 1.1

    if mBalLvl and mBalLvl < mPlotMax:
        ax.plot(mBalLvl, PermGroFac, marker=".", markersize=15, color="black")

#    mLvlGro, = ax.plot([mBalLvl, mBalLvl], [0, GroFacMax], color=color_mrktLev,
#                       label=r'$\mathbf{m}$-Level-Growth: $\mathbb{E}_{t}[{\mathbf{m}}_{t+1}/{\mathbf{m}}_{t}]$')
#    ax.legend(handles=[mLvlGro])

    ax.set_xlim(xMin, mPlotMax * 1.2)
    ax.set_ylim(GroFacMin, GroFacMax)

#    mNrmGro, = ax.plot([mTrgNrm, mTrgNrm], [0, GroFacMax], color=color_mrktRat,
    ax.legend(handles=[mLvlGro_lbl, mNrmGro_lbl])
    ax.legend(prop=dict(size=fssml))

    ax.text(mPlotMax+0.01, PermGroFac,
            r"${\Phi}$", fontsize=fssml, fontweight='bold')

#    ax.text(mPlotMax+0.01, Ex_cLvlGro[-1],
#            r"$\mathsf{E}_{t}[\mathbf{c}_{t+1}/\mathbf{c}_{t}]$", fontsize=fssml, fontweight='bold')
    ax.text(mPlotMax+0.01, APFac-0.003,
            r'$(R\beta)^{1/\rho}$', fontsize=fssml, fontweight='bold')

    # Ticks
#    ax.tick_params(labelbottom=False, labelleft=True, left='off',
    ax.tick_params(labelbottom=True, labelleft=True, left='off',
                   right='on', bottom='on', top='off')
    plt.setp(ax.get_yticklabels(), fontsize=fssml)

    _log.critical(f'mTrgNrm: {mTrgNrm:.3f}')
    _log.critical(f'mBalLvl: {mBalLvl:.3f}')

    plt.legend(fontsize='medium')
    ax.set_ylabel('Growth Factors', fontsize=fsmid)

    plt.show()
    return None


def makeBoundsFigure(UnempPrb, PermShkStd, TranShkStd, DiscFac, CRRA):
    inf_hor = IndShockConsumerType(quietly=True, messaging_level=logging.CRITICAL,
                                   **base_params)
    inf_hor.UnempPrb = UnempPrb
    inf_hor.PermShkStd = [PermShkStd]
    inf_hor.TranShkStd = [TranShkStd]
    inf_hor.DiscFac = DiscFac
    inf_hor.CRRA = CRRA
    inf_hor.update_income_process()

    inf_hor.solve(quietly=True, messaging_level=logging.CRITICAL)
    soln = inf_hor.solution[0]
    Bilt, Pars = soln.Bilt, soln.Pars

    cFunc = soln.cFunc
    mPlotMin, mPlotMax = 0, 25
    inf_hor.aXtraMax = mPlotMax

    # Retrieve parameters (makes code more readable)
    Rfree, EPermGroFac = Pars.Rfree, Pars.PermGroFac

    κ_Min = 1.0-(Rfree**(-1.0))*(Rfree*DiscFac)**(1.0/CRRA)
    h_inf = (1.0/(1.0-EPermGroFac/Rfree))
    def cFunc_Uncnst(mVec): return (h_inf - 1) * κ_Min + κ_Min*mVec
    def cFunc_TopBnd(mVec): return (1 - UnempPrb ** (1.0/CRRA)
                                    * (Rfree*DiscFac)**(1.0/CRRA)/Rfree)*mVec

    def cFunc_BotBnd(mVec): return (1 - (Rfree*DiscFac)**(1.0/CRRA)/Rfree) * mVec

    # Plot the consumption function and its bounds
    cPlotMaxLabel = r"c̅$(m) = (m-1+h)κ̲$"  # Use unicode kludge
    cPlotMinLabel = r"c̲$(m)= (1-\Phi_{R})m = κ̲ m$"

    # mKnk is point where the two upper bounds meet
    mKnk = ((h_inf-1) * κ_Min)/((1 - UnempPrb**(1.0/CRRA)*(Rfree*DiscFac)**(1.0/CRRA)/Rfree)-κ_Min)
    mBelwKnkPts, mAbveKnkPts = 50, 100
    mBelwKnk = np.linspace(mPlotMin, mKnk, mBelwKnkPts)
    mAbveKnk = np.linspace(mKnk, mPlotMax, mAbveKnkPts)
    mFullPts = np.linspace(mPlotMin, mPlotMax, mBelwKnkPts+mAbveKnkPts)

    plt.figure(figsize=(12, 8))
    plt.plot(mFullPts, cFunc(mFullPts), label=r'$c(m)$')
    plt.plot(mBelwKnk, cFunc_Uncnst(mBelwKnk), label=cPlotMaxLabel, linestyle="--")
    plt.plot(mAbveKnk, cFunc_Uncnst(mAbveKnk),
             label=r'Upper Bound $ = $ Min $[\overline{\overline{c}}(m),\overline{c}(m)]$', linewidth=2.5, color='black')
    plt.plot(mBelwKnk, cFunc_TopBnd(mBelwKnk), linewidth=2.5, color='black')
    plt.plot(mAbveKnk, cFunc_TopBnd(mAbveKnk), linestyle="--",
             label=r"$\overline{\overline{c}}(m) = κ̅m = (1 - ℘^{1/ρ}Φᵣ)m$")
    plt.plot(mBelwKnk, cFunc_BotBnd(mBelwKnk), color='red', linewidth=2.5)
    plt.plot(mAbveKnk, cFunc_BotBnd(mAbveKnk), color='red', label=cPlotMinLabel, linewidth=2.5)
    plt.tick_params(labelbottom=False, labelleft=False, left='off',
                    right='off', bottom='off', top='off')
    plt.xlim(mPlotMin, mPlotMax)
    plt.ylim(mPlotMin, 1.12*cFunc_Uncnst(mPlotMax))
    plt.text(mPlotMin, 1.12*cFunc_Uncnst(mPlotMax)+0.05, "$c$", fontsize=22)
    plt.text(mPlotMax+0.1, mPlotMin, "$m$", fontsize=22)
    plt.legend(fontsize='x-large')
    plt.show()
    return None


def makeTargetMfig(Rfree, DiscFac, CRRA, PermShkStd, TranShkStd):
    inf_hor = IndShockConsumerType(quietly=True, **base_params)
    inf_hor.Rfree = Rfree
    inf_hor.DiscFac = DiscFac
    inf_hor.CRRA = CRRA
    inf_hor.PermShkStd = [PermShkStd]
    inf_hor.TranShkStd = [TranShkStd]
    inf_hor.update_income_process()
    mPlotMin = 0
    mPlotMax = 250
    inf_hor.aXtraMax = mPlotMax
    inf_hor.solve(quietly=True, messaging_level=logging.CRITICAL)
    soln = inf_hor.solution[0]
    Bilt, cFunc = soln.Bilt, soln.cFunc
    cPlotMin = 0, cFunc(mPlotMax)

    if Bilt.GICMod:  # tattle
        soln.check_GICMod(soln, quietly=False, messaging_level=logging.WARNING)

    mBelwStE = np.linspace(mPlotMin, mPlotMax, 1000)
    EPermGroFac = inf_hor.PermGroFac[0]
    def EmDelEq0(mVec): return (EPermGroFac/Rfree)+(1.0-EPermGroFac/Rfree)*mVec
    cBelwStE_Best = cFunc(mBelwStE)  # "best" = optimal c
    cBelwStE_Sstn = EmDelEq0(mBelwStE)               # "sustainable" c
    mBalLvl = Bilt.mBalLvl

    plt.figure(figsize=(12, 8))
    plt.plot(mBelwStE, cBelwStE_Best, label="$c(m_{t})$")
    plt.plot(mBelwStE, cBelwStE_Sstn, label="$\mathsf{E}_{t}[\Delta m_{t+1}] = 0$")
    plt.xlim(mPlotMin, mPlotMax)
    plt.ylim(cPlotMin, cFunc(mPlotMax))
    plt.plot(
        [mBalLvl, mBalLvl],
        [0, 2.5],
        color="black",
        linestyle="--",
    )
    plt.tick_params(
        #        labelbottom=False,
        #        labelleft=False,
        #        left="off",
        right="off",
        #        bottom="off",
        top="off",
    )
    plt.text(0, 1.47, r"$c$", fontsize=26)
    plt.text(3.02, 0, r"$m$", fontsize=26)
    plt.text(mBalLvl - 0.05, -0.1, "m̌", fontsize=26)
    plt.legend(fontsize='x-large')
    plt.show()
    return None
