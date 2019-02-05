
#Import Plotting Packages

using Pkg; Pkg.add("Plots")
using Plots
Pkg.add("PyPlot")

gr()

#PART A and B

#Biological Constants 

#Closed to open complex = k2 = .024 s-1 (McClure Paper) 
kI = 0.024; #s-1

# Elongation rate = e_X =42 nts/sec

# Note for all links, they are too long to be pasted, so they are 
# put onto multiple lines. All are official cited in a separate
# doc.

#All citations also included in constants document

# https://bionumbers.hms.harvard.edu/bionumber.aspx?id=108488&ve
# r=3&trm=elongation+rate+in+E.+coli&org=
e_X = 42; # nt/sec
Lj = 3075; #nt/gene
kej = e_X/Lj #s-1

#Find K_x,j = McClureSlope*kI (see paper), converted to units of [mM]
McClureSlope = 1.04*10^-3; #mM*s
K_xj = McClureSlope*kI #mM

#RNAP concentration = 30 nM
# https://bionumbers.hms.harvard.edu/bionumber.aspx?id
# =100194&ver=8&trm=rnap+e+coli+M&org=
R_xt=30*10^-6; #mM

#Find the gene concentration 
Gj_initial = 2500; # copies/cell 
Gj_2 = Gj_initial / (6.02*10^23); #mol/cell
#Bionumbers: volume per cell = 6.7E-10 μL/cell
# https://bionumbers.hms.harvard.edu/bionumber.aspx?id
# =108815&ver=1&trm=volume+of+e+coli+cell&org=
Gj_3= Gj_2 / (6.7*10^-16) #mol/L = M
Gj = Gj_3*1000 # mM

tau = kej/kI
#also PART B

rxj = kej*R_xt*(Gj/(tau*K_xj+(tau+1)*Gj)) #mM/s
#ANSWER TO rxj in PART A

#PART C

pyplot()

#Constants
W1 = 0.26;
W2 = 300;
n = 1.5; 
K_c = 0.30; #mM
#find f(I)
f(I) = (I^n)/(K_c^n+(I^n));
#find u(I)
u(I) = (W1 + W2*f(I))/(1 + W1 + W2*f(I))

#growth rate, doubling time = 30 min
grow = 30; #min
#dilution is 1/grow, also convert to seconds
B_term = 1/(30*60) #s-1

#global half life of mRNA = 5 min
# https://bionumbers.hms.harvard.edu/bionumber.aspx?id=111927&ver
# =2&trm=mrna+half+life+e+coli&org=
halflife = 5;  #min
#convert to seconds
halflife2 = 5*60; #sec
#convert to degradation rate, assuming first order kinetics
kdeg = .693/halflife2 #s-1


m_j(I) = (rxj*u(I)/(kdeg+B_term)) #units of mM
#convert to uM for the plot
m_j_2(I) = m_j(I)*1000 #uM

#lower bound
a = m_j_2(0.0001)

#upper bound
b = m_j_2(10)

#Plotted the entire range incrementally 
#(will not plot all at once, possible due to Julia Box limitations)
plot(m_j_2,0.0001,0.001,xaxis=:log,xlims=(0.0001,10),ylims=(1*10^-2,10*10^-2), 
    xlabel = "I (mM)",ylabel="mRNA (uM)")
plot!(m_j_2,0.0009,0.01,xaxis=:log,xlims=(0.0001,10),ylims=(1*10^-2,10*10^-2))
plot!(m_j_2,0.001,1,xaxis=:log,xlims=(0.0001,10),ylims=(1*10^-2,10*10^-2))
plot!(m_j_2,0.2,10,xaxis=:log,xlims=(0.0001,10),ylims=(1*10^-2,10*10^-2))


