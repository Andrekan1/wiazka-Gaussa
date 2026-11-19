# -*- coding: utf-8 -*-
# Schemat modelu COMSOL (przekrój r–z, nie w skali). Uruchomienie: python schemat_modelu.py
import os, numpy as np, matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, FancyArrowPatch
HERE=os.path.dirname(os.path.abspath(__file__)); OUT=os.path.join(HERE,'..','img','schemat_modelu.png')
SURF='#fcfcfb'; T1='#0b0b0b'; T2='#52514e'; BLUE='#2a78d6'; ORANGE='#eb6834'; GRID='#c9c8c3'
plt.rcParams.update({'font.family':'DejaVu Sans','font.size':14})
fig,ax=plt.subplots(figsize=(16,9),dpi=100); fig.patch.set_facecolor(SURF); ax.set_facecolor(SURF)
L,H=10.0,5.6   # rysunek: długość (z) i wysokość (przekrój średnicowy), nie w skali
x0,y0=2.0,2.4
ax.add_patch(Rectangle((x0,y0),L,H,fill=False,ec=T2,lw=2))
for i in range(1,5): ax.plot([x0+i*L/5]*2,[y0,y0+H],color=GRID,lw=1,ls='--')   # warstwy siatki
yc=y0+H/2; ax.plot([x0-0.6,x0+L+0.8],[yc,yc],color=T2,lw=1,ls='-.')
yy=np.linspace(-1.5,1.5,200)
ax.fill_betweenx(yc+yy,x0,x0+1.4*np.exp(-(yy/0.75)**2),color=BLUE,alpha=0.35,lw=0)
ax.plot(x0+1.4*np.exp(-(yy/0.75)**2),yc+yy,color=BLUE,lw=2.5)
ax.fill_between([x0,x0+L],[yc-0.75]*2,[yc+0.75]*2,color=BLUE,alpha=0.10,lw=0)
ax.add_patch(FancyArrowPatch((x0+2.2,yc),(x0+L-0.5,yc),arrowstyle='-|>',mutation_scale=25,color=BLUE,lw=2.5))
ax.plot([x0,x0],[y0,y0+H],color=ORANGE,lw=5); ax.plot([x0+L,x0+L],[y0,y0+H],color=ORANGE,lw=5)
ax.plot([x0,x0+L],[y0+H,y0+H],color=T1,lw=5); ax.plot([x0,x0+L],[y0,y0],color=T1,lw=5)
def lab(x,y,t,ha='left',c=T1,fs=14,w='normal',va='center'): ax.text(x,y,t,ha=ha,va=va,color=c,fontsize=fs,fontweight=w)
lab(x0+L/2,y0+H+0.4,'powierzchnia boczna r = 5w₀ = 50 µm: PEC (pole ≈ 0)',ha='center',fs=14)
lab(x0-0.1,y0-0.35,'z = 0 (brzeg 3): Matched Boundary Condition\n+ pole padające E0i = profil modu w przewężeniu',ha='left',fs=13,va='top',c=ORANGE)
lab(x0+L+0.1,y0-0.35,'z = h (brzeg 4): Matched Boundary Condition\nbez pola padającego – wiązka wychodzi',ha='right',fs=13,va='top',c=ORANGE)
lab(x0+L/2+0.6,yc+2.0,'E = E₁(r)·exp(−i k₁·r),  k₁ = (0, 0, k₀)',ha='center',c=BLUE,fs=15,w='bold')
lab(x0+L/2+0.6,yc+1.45,'COMSOL liczy tylko wolnozmienną obwiednię E₁',ha='center',c=BLUE,fs=14)
lab(x0+L/2+0.6,yc-1.6,'siatka „Swept”: trójkąty z wejścia przeciągnięte wzdłuż z w N_L warstwach',ha='center',c=T2,fs=13)
lab(x0+L/2+0.6,yc-2.05,'(Nl = 5 w modelu HG, 10 w modelu LG)',ha='center',c=T2,fs=13)
lab(x0+L/2,0.55,'h = 3 µm = 3λ   ≪   z₀ = π·w₀²/λ = 314 µm        (rysunek nie w skali)',ha='center',c=T2,fs=14)
ax.set_xlim(0.5,13.5); ax.set_ylim(0,9.4); ax.axis('off')
fig.suptitle('Model COMSOL (HG i LG): geometria, warunki brzegowe, metoda obwiedni',x=0.06,y=0.97,ha='left',fontsize=21,fontweight='bold',color=T1)
fig.savefig(OUT,facecolor=SURF); print('ok',OUT)
