import sqlite3 as lite
import numpy as np
from scipy import stats
from scipy.optimize import curve_fit


##fitting function
def func(x,a,b,c,d):
	return b*np.exp(-c*x)/((d*x+1)**d)+a
	
path="/media/Data/PLOS_Biol/sqlite_S/PS06_S_E001_R50.sqlite"
con = lite.connect(path)
with con:
		cur = con.cursor()
		cur.execute("SELECT duration, infection_id FROM sampled_duration")
		rows = cur.fetchall()
		

xd = np.array([x[1] for x in rows])-1  #number of infections
yd = np.array([x[0] for x in rows])  #infection duration

#get range of infection times
newxd=range(xd.min(),xd.max()+2)
# mean of infection duration as a function of infection times
xys = stats.binned_statistic(xd,yd,'mean',bins=newxd)
bc = np.bincount(xd)
flt = np.where(bc>1)  #take only bins that have more than 1 incidents
x = flt[0]
y = xys[0][flt[0]]-14 
infectionTimesToImmune = max(x)
if len(yd[xd>max(x)])>0:
	clearanceRateConstantImmune = 1/(np.mean(yd[xd>max(x)]) - 14)
else:
	clearanceRateConstantImmune = 1/y[-1]

while True:
	try:
		popt,pcov = curve_fit(func,x,y,p0 = (0.01,50,0.0017,0.8))
		break
	except:
		pass
		
		
generalImmunityParams = list(popt)
print list(popt)
