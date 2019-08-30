# PC-SAFT Code

### Brief overview of PC-SAFT

PC-SAFT was one of the first variant of the SAFT equation to pop-up and, due to it's slightly simpler configuration, is more commonly used (like in ASPEN). The main distinction between PC-SAFT and it's predecessor is that the reference fluid is a hard-chain rather than a hard sphere as it is in the original SAFT (and in other variant such as SAFT-$\gamma$ Mie). The overall form is:

$$ \tilde{a}_{\mathrm{res}}=\tilde{a}_{\mathrm{HC}}+\tilde{a}_{\mathrm{Disp}}$$

Where:

$$ \tilde{a}=\frac{A}{Nk_{b}T}$$

In terms of required inputs:

* Number of Particles ($N$): Usually just set to $6.022\times 10^{23}$

* Volume ($V$): In $m^{3}$

* Temperature ($T$): In $K$

* Composition ($x$): Must be the molar amount

* Number of segments ($m$): This is the number of beads required to represent a particular component (need not be an integer)

* Molar mass ($Mr$) : Simply to convert between molar and mass units

* Segment diameter ($\sigma$ ): Measured in $m$ , this acts as the 'diameter' of the bead used to represent the component

* Potential well depth ($\epsilon$) : Measured in $K$ although, it is later converted to be in Joules within the code

* Cross interaction parameter ($k$): A key parameter that determines whether the interactions between two components are favourable. It is dependent on the two groups in question

  Within my code, I have tried to maintain consistent use of symbols to avoid it being cluttered. The required inputs above are the same that I use within the code

### Obtaining the Gibbs Free Energy of Mixing for polymer blends

As shown above, the PC-SAFT code only provides the residual helmholtz free energy but to obtain the Gibbs free energy of mixing, we also need the ideal contribution. Unfortunately, due to the vibrational and rotatioanal contributions being hard to quantify in a polymer, such a term does not exist. However, within the Gibbs free energy of mixing, these terms do cancel out, thus:

$$\Delta G_{\mathrm{mix}}=\Delta G_{\mathrm{Ideal}}+\Delta A_{\mathrm{res}}-Nk_{B}T\Delta Z_{\mathrm{res}}$$

### Changing from the $N,V,T$ ensemble to $N,P,T$ 

One issue with SAFT in general is that you can't specify the pressure due to the ensemble it is based on, however, this can be circumvented using a non-linear equation solver. Based on Kezheng's advice, I used a least-squares solver in scipy. As we will always be working in the liquid phase and since polymers are so incompressible, the corresponding volume at a particular pressure will be relatively constant near a packing fraction of 1. Thus, I have provided an initial guess for a packing fraction of 0.9 and based it on the Baker-Henderson diameter (different to the segment diameter) to make it a temperature dependent initial guess.

## Using the PC-SAFT Code for Gibbs Free Energy of Mixing

Once the two classes are within a particular code, all that is needed is to call GibbsMixing as such:

r=GibbsMixing(NSegment,Mr,Epsilon,Sigma,NParticles,Temperature,Pressure,xComp,k)

The Gibbs Free energy of mixing is then obtained by:

r.GibbsFreeEnergyMixing()