import seaborn as sns
import matplotlib.pyplot as plt
import os, sys
import csv
import pandas
import scipy
import plotting
reload (plotting)

root = '/hpc/mosa004/PCA/lobe/results'

def calculate_EIM1(age,fvc,dlco):
    eim1 = 1.3773 + (-0.0364*age) + (0.3754*fvc) + (-0.0409*dlco)
            
    return eim1
    
def calculate_EIM2(bmi,rvtlc):
    eim2 = 3.4902 + (-0.1642*bmi) + (0.0167*rvtlc)
            
    return eim2
    
def calculate_EIM3(age,tlc,dlco):
    eim3 = 4.8953 + (-0.018*age) + (-0.4557*tlc) + (-0.0506*dlco)
            
    return eim3
    
def calculate_EEM1(age,frc,rv):
    eem1 = 0.7802 + (-0.0257*age) + (0.4243*frc) + (-0.359*rv)
            
    return eem1
    
def calculate_EEM2(rvtlc,dlco,frc,fev1,pefr):
    eem2 = 4.1438 + (-0.0383*rvtlc) + (-0.093*dlco) + (-0.5009*frc) + (0.572*fev1) + (-0.1269*pefr)
            
    return eem2
    
def calculate_EEM3(age,vc,rv):
    eem3 = 3.3263 + (-0.0317*age) + (-0.5316*vc) + (0.2464*rv)
            
    return eem3
    
#~ files = ['AGING001','AGING004','AGING006',
#~ 'AGING010','AGING012','AGING013','AGING015','AGING017',
#~ 'AGING020','AGING021','AGING022','AGING023','AGING024',
#~ 'AGING025','AGING026','AGING027','AGING028','AGING030',
#~ 'AGING031','AGING040','AGING041','AGING043','AGING045',
#~ 'AGING046','AGING048','AGING050','AGING051','AGING052',
#~ 'AGING054','AGING055']
#~ files = ['H18','H653','H673','H682','H684',
#~ 'H1335','H5977','H6229','H6310','H6785','H11303','H11750',
#~ 'H11779','H12076','H12816','H670','H245','H817','H6508','H7311','H7351',
#~ 'H8753','H10769','H11194','H12051']

fname = root+'/data_for_analysis.csv'
#~ for i in range(len(files)):
with open(fname) as csvfile:
    next(csvfile)
    sub = 'AGING031'
    readCSV = csv.reader(csvfile, delimiter=',')
    for row in readCSV:
        s = row[1]
        if s == sub:
            age = int(row[2])
            bmi = float(row[4])
            fvc = float(row[5])
            fev1 = float(row[6])
            frc = float(row[8])
            tlc = float(row[9])
            vc = float(row[10])
            rv = float(row[11])
            pefr = float(row[13])
            rvtlc = float(row[14])
            dlco = float(row[15])

#~ print "Subject = %s" %sub
eim1 = calculate_EIM1(age,fvc,dlco)
#~ print "EIM1 = %.2f" %eim1
eim2 = calculate_EIM2(bmi,rvtlc)
#~ print "EIM2 = %.2f" %eim2
eim3 = calculate_EIM3(age,tlc,dlco)
#~ print "EIM3 = %.2f" %eim3
eem1 = calculate_EEM1(age,frc,rv)
#~ print "EEM1 = %.2f" %eem1
eem2 = calculate_EEM2(rvtlc,dlco,frc,fev1,pefr)
#~ print "EEM2 = %.2f" %eem2
eem3 = calculate_EEM3(age,vc,rv)
#~ print "EEM3 = %.2f" %eem3

print "%.2f,%.2f,%.2f,%.2f,%.2f,%.2f" %(eim1,eim2,eim3,eem1,eem2,eem3)
