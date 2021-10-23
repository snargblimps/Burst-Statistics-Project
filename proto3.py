from scipy.fft import fft
import numpy as np
import pandas as pd
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
%matplotlib notebook
import statistics as st
import pandas as pd
import time
from tabulate import tabulate

## DATA ##

tInit = time.time()
kernel = 100
infile = r"C:\Users\mails\Desktop\Et sidus oritur\ARate.dat.gz"
rawdata2 = (np.loadtxt(infile)[:,10])
rawdata = (np.loadtxt(infile)[:,9])
combin = rawdata + rawdata2
times = np.loadtxt(infile)[:,0]

## LOOP VARIABLES ##

prev = 0
nexte = 1
current = 0
j = 0
serial = 1

## MOVING AVERAGE ##

def movingAverage(x,y,N=1000):                                # function calculating and storing moving average 
    arrSize = len(x)
    m=np.array(np.zeros(arrSize))
    for i in range(N):
        np.put(m,i, (np.sum((y[0:i+N]))/(i+N)))
    for i in range(N,arrSize-N):
        np.put(m,i,(np.sum(y[i-N:i+N+1])/(2*N+1)))
    for i in range(arrSize-N,arrSize):
        np.put(m,i, (np.sum((y[i-N:arrSize]))/(arrSize-i+N)))
    return(m)
m = movingAverage(times,combin,kernel) # function call with the parameters being the spliced datasets and the kernel that was set
filt = combin              # Currently not using the moving average, just using the unfiltered combination of row 10 and 11


## REMOVING INITIAL BURSTS ##

prelimpeaks, _ = find_peaks(filt, height= 50, distance=100) #finding peaks that are within 100 indices, from one peak to the next
prelimpeaks = np.array(prelimpeaks) #creating an array of all the peaks
while ((prelimpeaks[nexte]-prelimpeaks[current])< 5000):#loop to find the threshold at which the peaks start becoming 5000 indices apart
    threshold = prelimpeaks[nexte]
    nexte = nexte+1
    current = current+1
noiseless_filt = filt[threshold:] #splicing the data with respect to the threshold

## IDENTIFY PEAKS AND PROPERTIES AGAIN ##

peaks,properties = find_peaks(noiseless_filt, height= 100, distance=100, width=(None,None),plateau_size = (None,None)
                              ,prominence=(None,None)) #finding peaks again with the spliced data, doesn't really matter that the distance is 100 again
properties["prominences"], properties["widths"] #calling prominences and widths to use in finding, width, and plotting the structure etc.
filtered, = plt.plot(times[threshold:],noiseless_filt,label='filtered data') #plotting the spliced version of the original data
plt.plot((peaks*5e-6 + (threshold*5e-6)),noiseless_filt[peaks], "x")# plotting the peaks on a time scale instead of against indices
# scaling the peaks (index numbers) to time in Million years offsetting it by the threshold also converted into Myrs
plt.xlabel("M yrs")

plt.vlines(x=(peaks*5e-6 + (threshold*5e-6)), ymin=noiseless_filt[peaks] - properties["prominences"],ymax = noiseless_filt[peaks],color = "C1")
# plotting vertical lines along each peak, from the base line to the y value of each peak
plt.hlines(y=properties["width_heights"], xmin=((properties["left_ips"])*5e-6 + (threshold*5e-6)),
           xmax=((properties["right_ips"])*5e-6 + (threshold*5e-6)), color = "C3")
#plotting horizontal lines from the left most point of the peak uunder consideration to the right most point of that peak

#just converting the property data into numpy arrays for ease of use later on, this section may not be needed
print(peaks)
peaks = np.array(peaks)
properties["peak_heights"] = np.array(properties["peak_heights"])
properties["widths"] = np.array(properties["widths"])
properties["prominences"] = np.array(properties["prominences"])
properties["left_ips"] = np.array(properties["left_ips"])
properties["right_ips"] = np.array(properties["right_ips"])

## FILL PROPERTY TABLES ##
newwidths = []
newheights = []
newpeaks = []
newprom = []
newleft = []
newright = []
rise = []
#adding the first element as the name of the column
newwidths.append("width (yrs)")
newheights.append("height")
newprom.append("prominences")
newleft.append("start time")
newright.append("end time")
rise.append("rise time (yrs)")
Sno = []
Sno.append("peaks")

# looping and filling the property tables, after formatting to two decimal places and converting to units of yrs
while (j<(peaks.size)):
    newpeaks.append(peaks[j])
    newheights.append("{:.2f}".format(((properties["peak_heights"])[j])))
    newwidths.append("{:.2f}".format(((properties["widths"])[j])*5))
    newprom.append("{:.2f}".format(((properties["prominences"])[j])*5))
    newleft.append("{:.2f}".format(((properties["left_ips"])[j])*5))
    newright.append("{:.2f}".format(((properties["right_ips"])[j])*5))
    rise.append((peaks[j]-((properties["left_ips"])[j]))*5)
    Sno.append(serial)
    serial = serial +1;
        
    j= j+1

newpeaks.append(np.amax(peaks))
newpeaks = np.array(newpeaks)    

## OUTPUT PLOTS AND TABLE ##

plt.savefig('Noisy Peaks removed.pdf')
print(noiseless_filt[peaks])
table = [Sno,newheights,newwidths,newprom,newleft,newright,rise]
print(tabulate(table, tablefmt='fancy_grid'))
print(threshold)
print(peaks.shape)
