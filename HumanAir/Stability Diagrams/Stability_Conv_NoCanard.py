import HumanAir.LoadingDiagram.Plotting as plt
import matplotlib.pyplot as plot
import numpy as np

SSh=0.2146

xlst=np.arange(-0.3,1.1, 0.1)
Stability=[-0.210404306354706,-0.172171881752751,-0.133939457150796,-0.0957070325488414,-0.0574746079468865,-0.0192421833449316,0.0189902412570234,0.0572226658589783,0.0954550904609332,0.133687515062888,0.171919939664843,0.210152364266798,0.248384788868753,0.286617213470708]
Stability_NoMargin=[-0.229520518655684,-0.191288094053729,-0.153055669451774,-0.114823244849819,-0.0765908202478639,-0.038358395645909,-0.000125971043954093,0.0381064535580008,0.0763388781599558,0.114571302761911,0.152803727363866,0.191036151965821,0.229268576567775,0.26750100116973]
Controllability=[0.903234936729755,0.778030530617833,0.652826124505911,0.527621718393989,0.402417312282067,0.277212906170145,0.152008500058223,0.0268040939463015,-0.0984003121656204,-0.223604718277542,-0.348809124389464,-0.474013530501386,-0.599217936613308,-0.72442234272523]

Min=0.25
Max=0.75
ylst=np.arange(SSh-0.02,SSh+0.03, 0.01)

xline=np.linspace(Min, Max, 2)
yline=np.ones(len(xline))*SSh

plot.figure()

plt.Ploty(xlst, Stability, "Stability",colour='green', linestyle='solid')
plt.Ploty(xlst, Stability_NoMargin, "Stability_NoMargin",colour='green', linestyle='dashed')
plt.Ploty(xlst, Controllability, "Controllability",colour='blue', linestyle='solid')
plt.Plotx(Min, ylst, "CG Excursion", colour='red', linestyle='solid')
plt.Plotx(Max, ylst, "",colour='red', linestyle='solid')
plot.plot(xline,yline, color="red", linestyle='solid')
plot.fill_between(xlst, 0, Stability, color='red', alpha=.1)
plot.fill_between(xlst,0, Controllability, color='red', alpha=.1)

plot.xlabel(r'$X_{cg}/MAC$', fontsize=12, loc='right')
plot.ylabel(r'$S_h/S$', fontsize=12, loc='top')

plot.axhline(color='black', lw=0.5)
plot.axvline(color='black', lw=0.5)
plot.legend()
plot.xlim((0,1))
plot.ylim((0,0.5))
plot.title("Stability Diagram - Conventional No Canard")
plot.show()