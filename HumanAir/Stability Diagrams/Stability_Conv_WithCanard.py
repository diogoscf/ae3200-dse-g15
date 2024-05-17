import HumanAir.LoadingDiagram.Plotting as plt
import matplotlib.pyplot as plot
import numpy as np

SSh=0.2257

xlst=np.arange(-1.4,1.3, 0.1)
Stability=[0.63824447271722,0.598505635363738,0.558766798010256,0.519027960656774,0.479289123303292,0.43955028594981,0.399811448596329,0.360072611242846,0.320333773889365,0.280594936535883,0.240856099182401,0.201117261828919,0.161378424475437,0.121639587121955,0.0819007497684728,0.0421619124149908,0.00242307506150882,-0.0373157622919731,-0.0770545996454551,-0.116793436998937,-0.156532274352419,-0.196271111705901,-0.236009949059383,-0.275748786412865,-0.315487623766347,-0.355226461119829,-0.394965298473311]
Stability_NoMargin=[0.658113891393961,0.618375054040479,0.578636216686997,0.538897379333515,0.499158541980033,0.459419704626551,0.419680867273069,0.379942029919587,0.340203192566105,0.300464355212624,0.260725517859142,0.22098668050566,0.181247843152178,0.141509005798696,0.101770168445214,0.0620313310917318,0.0222924937382498,-0.0174463436152321,-0.0571851809687141,-0.0969240183221961,-0.136662855675678,-0.17640169302916,-0.216140530382642,-0.255879367736124,-0.295618205089606,-0.335357042443088,-0.37509587979657]
Controllability=[0.347979972566948,0.328574809701961,0.309169646836974,0.289764483971986,0.270359321106999,0.250954158242012,0.231548995377025,0.212143832512038,0.192738669647051,0.173333506782063,0.153928343917076,0.134523181052089,0.115118018187102,0.0957128553221148,0.0763076924571277,0.0569025295921405,0.0374973667271534,0.0180922038621662,-0.00131295900282094,-0.0207181218678081,-0.0401232847327952,-0.0595284475977824,-0.0789336104627696,-0.0983387733277567,-0.117743936192744,-0.137149099057731,-0.156554261922718]

Min=-0.77
Max=-0.39
ylst=np.arange(SSh-0.02,SSh+0.02, 0.005)

xline=np.linspace(Min, Max, 2)
yline=np.ones(len(xline))*SSh

plot.figure()

plt.Ploty(xlst, Stability, "Stability", colour="green", linestyle="solid")
plt.Ploty(xlst, Stability_NoMargin, "Stability_NoMargin", colour="green", linestyle="dashed")
plt.Ploty(xlst, Controllability, "Controllability", colour='blue', linestyle='solid')
plt.Plotx(Min, ylst, "CG Excursion", colour='red', linestyle='solid')
plt.Plotx(Max, ylst, "", colour='red', linestyle='solid')
plot.plot(xline,yline, color="red", linestyle='solid')
plot.fill_between(xlst,  Stability,1, color='red', alpha=.1)
plot.fill_between(xlst,0, Controllability, color='red', alpha=.1)

plot.xlabel(r'$X_{cg}/MAC$', fontsize=12, loc='right')
plot.ylabel(r'$S_h/S$', fontsize=12, loc='top')

plot.axhline(color='black', lw=0.5)
plot.axvline(color='black', lw=0.5)
plot.legend()
plot.xlim((-1,0.5))
plot.ylim((0,0.5))

plot.title("Stability Diagram - Conventional No Canard")
plot.show()