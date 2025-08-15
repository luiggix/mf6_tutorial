time = 0.0
year = 2018
station = 260
file = open("METE_SVAT.INP", "w")
with open('KNMI_20181231.txt', 'r') as f:
    for n, line in enumerate(f):
        if line.startswith('#'):
            continue
        dummy_station, dummy_date, rh, ev24 = line.split(',')
        daily_rate = max(0.0, float(rh.strip())*0.1)
        evapo_ref = float(ev24.strip())*0.1
        file.write('{:15.3f}{:5d}{:10.2f}{:10.2f}{:10d}\n'.format(time,year,daily_rate,evapo_ref,station))
        
        time = time+1.0
        
file.close()