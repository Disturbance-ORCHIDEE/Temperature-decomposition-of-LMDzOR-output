# Temperature-decomposition-of-LMDzOR-output
Temperature decomposition uses a meta model of the processes resulting in a temperature change to quantify the contribution of each process of from the difference between only a control and treatment simulation rather than a full factorial experiment.

# Description of the meta model
In ecology, the temperature of air near the surface of the Earth, i.e., at 2\,m above the ground, is among the most relevant climate variables, as it represents the temperature experienced by plants, animals, and human society. The surface energy budget conceptualizes the land-atmosphere interface as a coupled system: a surface layer overlain by an atmospheric layer. To quantitatively evaluate this system, we adopt the theoretical foundations laid by previous research \citep{Juang2007, naudts2016}. After ignoring heat storage, heat advection, and bond energy in the organic molecules that make up the ecosystem, the simplified energy budget of the surface layer can be written as:

\begin{equation}
    \label{eq:energy_budget}
    R_{net} = LE + H + G 
\end{equation}

where $R_{net}$ is the net radiation, $H$ sensible heat, $LE$ latent heat and $G$ ground heat flux all expressed in W m$^{-2}$. The $R_{net}$ flux is calculated from the longwave incoming radiation ($R_{li}$), longwave outgoing radiation ($R_{lo}$), shortwave incoming radiation($R_{si}$), and shortwave outgoing radiation ($R_{so}$):

\begin{equation}
    \label{eq:rn_flux}
    R_{net} = R_{si} - R_{so} + R_{li} - R_{lo} 
\end{equation}

where the shortwave outgoing radiation ($R_{so}$) can be expressed as a function of surface albedo ($\alpha$) and shortwave incoming radiation ($R_{si}$):

\begin{equation}
    \label{eq:rso_calculation}
    R_{so} = \alpha R_{si}
\end{equation}

According to the law of Stefan-Boltzmann, longwave radiation can be written as a function of emissivity, temperature, and the Stefan-Boltzmann constant $\sigma$ ($5.670367 \times 10^{-8}$ W m$^{-2}$ K$^{-4}$). The surface layer is characterized by its emissivity ($\varepsilon_s$), albedo ($\alpha$) and temperature ($T_s$) and is overlain by an atmospheric layer with an atmospheric emissivity ($\varepsilon_a$) and air temperature ($T_a$). The longwave radiation leaving the surface layer ($R_{lo}$) can therefore be written as a function of the surface layer emissivity ($\varepsilon_s$), its temperature, the Stefan-Boltzmann constant $\sigma$ ($5.67 \times 10^{-8}$ W m$^{-2}$ K$^{-4}$) and a term accounting for the incoming longwave radiation that will be reflected by the surface layer:

\begin{equation}
    \label{eq:rlo_stefan_boltzmann}
    R_{lo} = \varepsilon_s \sigma T_s^4 + (1 - \varepsilon_s) \varepsilon_a \sigma T_a^4
\end{equation}

Following a similar reasoning, the incoming longwave radiation reaching the surface layer ($R_{li}$) is a function of the atmospheric emissivity and air temperature:

\begin{equation}
    \label{eq:rli_calculation}
    R_{li} = \varepsilon_a  \sigma  T_a^4
\end{equation}

Hence, the energy budget of the surface layer can thus be written as:

\begin{equation}
    \label{eq:simplified_energy_budget}
    (1 - \alpha)R_{si} - (\varepsilon_s \sigma T_s^4) + (\varepsilon_s  \varepsilon_a  \sigma T_a^4) = LE + H + G
\end{equation}

To quantify the contribution of each component of the energy budget to a change in air temperature, the equation is reorganized to obtain a function of air temperature ($T_a$):

\begin{equation}
    \label{eq:reorganized_energy_budget}
    \varepsilon_a \sigma T_a^4 = \frac{1}{\varepsilon_s} \left( LE + H + G + (\varepsilon_s  \sigma  T_s^4)  - (1 - \alpha) R_{si}\right)
\end{equation}

To attribute the change in $T_a$ to its components, i.e., $R_{si}$, $\alpha$, $T_s$, $\varepsilon_a$, $H$, $LE$ and $G$, each component is treated as a mathematically independent contributor to the air temperature response \citep{naudts2016, zeng2017} and surface emissivity is assumed constant. In the LMDzOR model surface emissivity is not only constant, it is also set to one, assuming the land surface is a perfect black body for long wave radiation. We take the partial derivatives of Eq. \ref{eq:reorganized_energy_budget} with respect to each component. The left side can be written as :

\begin{equation}
    \label{eq:left_side}
    f_l = \frac{\partial (\varepsilon_a \sigma T_a^4)}{\partial T_a} \Delta T_a + \frac{\partial (\varepsilon_a \sigma T_a^4)}{\partial \varepsilon_a} \Delta \varepsilon_a = 4 T_a^3 \varepsilon_a \sigma \Delta T_a + \sigma T_a^4 \Delta \varepsilon_a
\end{equation}

The right side can be written as:

\begin{equation}
    \label{eq:right_side}
    f_r = \Bigl( \Delta LE + \Delta H + \Delta G + R_{si}\Delta\alpha 
    - (1 - \alpha)\Delta R_{si} + 4 \sigma T_s^3 \Delta T_s \Bigr)
\end{equation}

Combining Eqs. \ref{eq:left_side} and \ref{eq:right_side} then gives:

\begin{equation}
    \label{eq:final_simplified}
    \Delta T_a = \frac{1}{4 \varepsilon_a \sigma T_a^3} \left(
    \Delta LE + \Delta H + \Delta G + R_{si} \Delta \alpha 
    - (1 - \alpha)\Delta R_{si} + 4\sigma T_s^3 \Delta T_s 
    - \sigma T_a^4 \Delta \varepsilon_a  
    \right)
\end{equation}

This temperature decomposition equation can be applied to any two experiment scenarios with identical temporal and spatial domains. The $\Delta T_a$ on the left-hand side of Eq. \ref{eq:final_simplified}, representing the air temperature difference between a control scenario and a reference scenario, can be estimated from a linearized surface energy budget. The $\Delta LE$, $\Delta H$, $\Delta G$, $\Delta R_{si}$, $\Delta T_s$, $\Delta\alpha$, and $\Delta\varepsilon_a$ on the right-hand side of Eq. \ref{eq:final_simplified} are calculated by subtracting the control experiment from the treatment experiment (i.e., $\Delta X = X_{\text{treatment}} - X_{\text{control}}$). The variables $R_{si}$, $T_a$, $T_s$, $\alpha$, and $\varepsilon_a$ on the right-hand side of Eq. \ref{eq:final_simplified} are then taken from the control experiment. For both the treatment and control experiments, $\alpha$ is calculated as the ratio of $R_{so}$ to $R_{si}$ (Eq. \ref{eq:rso_calculation}), and $\varepsilon_a$ is derived from $R_{li}$ by inverting the Stefan--Boltzmann law (Eq. \ref{eq:rli_calculation}).

\(\Delta\)$T_a$ can now be obtained from the actual simulated change between the control and treatment experiments, or calculated from the decomposition approach given in Eq.~\ref{eq:final_simplified}. The actual simulated temperature change includes all the processes accounted for in the climate model, whereas the decomposed temperature change only includes the major processes accounted for in the meta-model used in the decomposition. Therefore, a residual term needs to be added to ensure equality between both approaches: 

\begin{equation}
    \label{eq:decomposition_vs_actual}
    \Delta T_a{^\text{ actual}} = \Delta T_a^{\text{decomposed}} + \epsilon
\end{equation}

The residual term ($\epsilon$) includes processes accounted for in the climate model but omitted from the meta-model used in the decomposition, i.e., heat transport through vertical and horizontal advection, as well as the higher-order terms that were neglected in the decomposition \citep{zeng2017, naudts2016}. When assuming that the residual term is dominated by the vertical and horizontal advection $\Delta T_a^{\text{advection}}$ can be estimated as:

\begin{equation}
    \label{eq:circulation}
    \Delta T_a^{\text{advection}} = \Delta T_a^{\text{actual}} - \Delta T_a^{\text{decomposed}}
\end{equation}

# Dimensional analysis
Dimensional analysis for Eq. \ref{eq:final_simplified} was used to show that the units on both sides of the temperature decomposition equation are consistent: 

\begin{equation}
\label{eq:unit_check_temperature_1}
\underbrace{\mathrm{K}}_{\Delta T_a}=
\underbrace{\displaystyle \frac{1}{1\cdot \mathrm{W\,m^{-2}K^{-4}} \cdot \mathrm{K^3}}}_{\frac{1}{4\varepsilon_a \sigma T_a^3}}
\left[
\begin{aligned}
& \underbrace{\mathrm{W\,m^{-2}}}_{\Delta LE} \\
& + \underbrace{\mathrm{W\,m^{-2}}}_{\Delta H} \\
& + \underbrace{\mathrm{W\,m^{-2}}}_{\Delta G} \\

# Known shortcomings

# License
& + \underbrace{\mathrm{W\,m^{-2}} \cdot 1}_{R_{si}\Delta\alpha} \\
& - \underbrace{1 \cdot \mathrm{W\,m^{-2}}}_{(1 - \alpha)\Delta R_{si}} \\
& + \underbrace{\mathrm{W\,m^{-2}K^{-4}} \cdot \mathrm{K^3} \cdot \mathrm{K}}_{4\sigma T_s^3 \Delta T_s} \\
