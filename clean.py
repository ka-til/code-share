import numpy as np
import matplotlib.pylab as plt
import re
from tables import open_file
from statistics import median, mode
from scipy.optimize import curve_fit, minimize
from math import acos, degrees, log
from matplotlib.colors import LogNorm

class clean(object):
    def __init__(self, file_name):
        self.file_name = file_name
        self.file_path = '/Users/cosmos2299/documents/data/'
        l = open_file(self.file_path + self.file_name)
        list(l)
    
        self.atstamp = l.root.absolute_timestamp[150:].flatten()
        self.r_0 = l.root.rising_0[150:].flatten()
        self.r_1 = l.root.rising_1[150:].flatten()
        self.r_2 = l.root.rising_2[150:].flatten()
        self.r_3 = l.root.rising_3[150:].flatten()
        self.f_0 = l.root.falling_0[150:].flatten()
        self.f_1 = l.root.falling_1[150:].flatten()
        self.f_2 = l.root.falling_3[150:].flatten()
        self.f_3 = l.root.falling_3[150:].flatten()
    
        self.p_jumps = ((self.atstamp.size - self.atstamp[(self.atstamp < 1e11) & (self.atstamp > -1e7)].size)/
                     self.atstamp.size) * 100
        print('percentage of high jumps in the file -', self.p_jumps)
    
    #POCAM, sDOM, PMT, LED and voltage used
    def P_S_used(self):                                         
        POCAM_re = re.compile('P[1-2]')
        SDOM_re = re.compile('SDOM[1-5]')
        frequency_re = re.compile('[0-9][0-9][0-9][0-9]Hz')
        voltage_re = re.compile('[0-9][0-9]V')
        flash_time_re = re.compile('[0-9][0-9]s')
        LED_re = re.compile('P[0-9]_[a-z]')
        PMT_re = re.compile('hld_[a-z]')

        self.SDOM_num = SDOM_re.findall(self.file_name)
        self.POCAM_num = POCAM_re.findall(self.file_name)
        self.frequency = frequency_re.findall(self.file_name)
        voltage = voltage_re.findall(self.file_name)
        flash_time = flash_time_re.findall(self.file_name)
        LED = LED_re.findall(self.file_name)
        PMT = PMT_re.findall(self.file_name)
    
        if PMT == ['hld_d']:
            self.PMT = 'down'
        else:
            self.PMT = 'up'
    
        if LED == ['P1_b']:
            LED = 'blue'
        if LED == ['P1_v']:
            LED = 'violet'
    
        if LED == ['P2_b']:
            LED = 'blue'
        if LED == ['P2_v']:
            LED = 'violet'
        if LED == ['P1_o']:
        	LED = 'orange'
        if LED == ['P1_u']:
        	LED == 'uv'
        
        graph_title = [self.POCAM_num, self.SDOM_num, self.PMT, LED, voltage, self.frequency]
        values = ','.join(str(v) for v in graph_title)
         
    #timestamp graph
        plt.figure(figsize=(10,9))
        plt.title(values, fontsize = 22)
        plt.ylabel('absolute_timestamps(ns)', fontsize = 19)
        plt.xlabel('index', fontsize = 19)
        plt.plot(self.atstamp, '.')
        plt.savefig(self.file_path + '/graphs/' + values + 'high_jumps.jpeg', dpi = 200)
    
    #cleaning large jumps
        elim_h_jumps = (self.atstamp < 1e11) & (self.atstamp > -1e7)
        abs_elim = self.atstamp[elim_h_jumps]
        rising_0_elim = self.r_0[elim_h_jumps]
        rising_1_elim = self.r_1[elim_h_jumps]
        rising_2_elim = self.r_2[elim_h_jumps]
        rising_3_elim = self.r_3[elim_h_jumps]
        falling_0_elim = self.f_0[elim_h_jumps]
        falling_1_elim = self.f_1[elim_h_jumps]
        falling_2_elim = self.f_2[elim_h_jumps]
        falling_3_elim = self.f_3[elim_h_jumps]
    
        plt.figure(figsize=(10,9))
        plt.title(values, fontsize = 22)
        plt.ylabel('absolute_timestamp(ns)', fontsize = 19)
        plt.xlabel('index', fontsize = 19)
        plt.plot(abs_elim, '.')
        plt.savefig(self.file_path + '/graphs/' + values + 'high_jumps_cleaned.jpeg', dpi = 200)
        plt.show()
   
        plt.figure(figsize=(10,9))
        plt.title(values + ' Negative Timestamps', fontsize = 18)
        plt.ylabel('absolute_timestamp(ns)', fontsize = 19)
        plt.plot(abs_elim, '.')
        plt.ylim(-1e7, 0)
        plt.savefig(self.file_path + '/graphs/' + values + 'negative_values.jpeg', dpi = 200)
        plt.show()

    #time difference graph
        abs_elim_diff = abs_elim[1:] - abs_elim[:-1]
    
        abs_elim_bool = abs_elim[:-1][abs_elim_diff < 0]
    
        plt.figure(figsize=(10,9))
        plt.title(values, fontsize = 22)
        plt.ylabel('absolute_timestamp difference(ns)', fontsize = 19)
        plt.xlabel('index', fontsize = 19)
        plt.plot(abs_elim_diff, '.')
        plt.savefig(self.file_path + '/graphs/' + values + 'timestamp_differences.jpeg', dpi = 200)
        plt.show()
    
        if abs_elim_bool.size != 0:
            plt.figure(figsize=(10,9))
            plt.ylabel('absolute_timestamp(ns)', fontsize = 19)
            plt.xlabel('index', fontsize = 19)
            plt.title(values + ' Jumps in Timestamps', fontsize = 16)
            plt.plot(abs_elim[abs_elim_diff.argmin()-10:abs_elim_diff.argmin()+10], '.')
            plt.savefig(self.file_path + '/graphs/' + values + 'small_jumps.jpeg', dpi = 200)
            plt.show()
    
    #cleaning small jumps
        abs_elim_diff_2 = abs_elim[1:] - abs_elim[:-1]
    
        s_jump_index = []
        s_jump_1 = []
        s_jump_2 = []
        list_1 = np.array([])

        for r in range(0, abs_elim_diff_2.size):
            if abs_elim_diff_2[r] < 0:
                s_jump_index.append(r)
                s_jump_1.append(abs_elim_diff_2[r])

        for t in range(0, len(s_jump_index)):
            select = abs_elim_diff_2[s_jump_index[t] - 10:s_jump_index[t]]
            x = s_jump_index[t] - (10 - (np.abs(select+abs_elim_diff_2[s_jump_index[t]])).argmin())
            jump_length = s_jump_index[t] - x
    
            if jump_length == 1:
                list_1 = np.append(list_1, [x+1])
            if jump_length == 2:
                list_1 = np.append(list_1,[x+1, x+2])
            if jump_length == 3:
                list_1 = np.append(list_1,[x+1, x+2, x+3])
            if jump_length == 4:
                list_1 = np.append(list_1,[x+1, x+2, x+3, x+4])
            if jump_length == 5:
                list_1 = np.append(list_1,[x+1, x+2, x+3, x+4, x+5])
            if jump_length == 6:
                list_1 = np.append(list_1,[x+1, x+2, x+3, x+4, x+5, x+6])
            if jump_length == 7:
                list_1 = np.append(list_1,[x+1, x+2, x+3, x+4, x+5, x+6, x+7])
            if jump_length == 8:
                list_1 = np.append(list_1,[x+1, x+2, x+3, x+4, x+5, x+6, x+7, x+8])
            if jump_length == 9:
                list_1 = np.append(list_1,[x+1, x+2, x+3, x+4, x+5, x+6, x+7, x+8, x+9])
            if jump_length == 10:
                list_1 = np.append(list_1,[x+1, x+2, x+3, x+4, x+5, x+6, x+7, x+8, x+9, x+10])
   
        print(list_1)
    
        abs_elim_3e = np.delete(abs_elim, list_1)
        self.abs_elim_diff_3 = abs_elim_3e[1:] - abs_elim_3e[:-1]
        
        rising_0_elim_3e = np.delete(rising_0_elim, list_1)
        rising_1_elim_3e = np.delete(rising_1_elim, list_1)
        rising_2_elim_3e = np.delete(rising_2_elim, list_1)
        rising_3_elim_3e = np.delete(rising_3_elim, list_1)
        falling_0_elim_3e = np.delete(falling_0_elim, list_1)
        falling_1_elim_3e = np.delete(falling_1_elim, list_1)
        falling_2_elim_3e = np.delete(falling_2_elim, list_1)
        falling_3_elim_3e = np.delete(falling_3_elim, list_1)
        
        plt.figure(figsize=(10,9))
        plt.title(graph_title, fontsize = 22)
        plt.ylabel('absolute_timestamp difference(ns)', fontsize = 19)
        plt.xlabel('index', fontsize = 19)
        #plt.ylim(-10000, 0)
        plt.plot(self.abs_elim_diff_3, '.')
        plt.savefig(self.file_path + '/graphs/' + values + 'small_jumps_cleaned.jpeg', dpi = 200)
        plt.show()
        
        dt_mean = (abs_elim_3e[1:] - abs_elim_3e[:-1]).mean()
        
        #cleaning faulty falling time events
        boolean = ((falling_0_elim_3e != 0) & ((falling_1_elim_3e - rising_1_elim_3e) >= 0) & 
                   ((falling_2_elim_3e - rising_2_elim_3e) >= 0) & ((falling_3_elim_3e - rising_3_elim_3e) >= 0))
        abs_elim_3 = abs_elim_3e[boolean]
        rising_0_elim_3 = rising_0_elim_3e[boolean]
        rising_1_elim_3 = rising_1_elim_3e[boolean]
        rising_2_elim_3 = rising_2_elim_3e[boolean]
        rising_3_elim_3 = rising_3_elim_3e[boolean]
        falling_0_elim_3 = falling_0_elim_3e[boolean] 
        falling_1_elim_3 = falling_1_elim_3e[boolean]
        falling_2_elim_3 = falling_2_elim_3e[boolean]
        falling_3_elim_3 = falling_3_elim_3e[boolean]
        
        #cleaning faulty falling time events
        delete = []
        for i in range(0, abs_elim_3.size):
            if ((falling_1_elim_3[i] != 0 and rising_1_elim_3[i] == 0)
                or (falling_2_elim_3[i] != 0 and rising_2_elim_3[i] == 0)
                or (falling_3_elim_3[i] != 0 and rising_3_elim_3[i] == 0)
               ):
                delete.append(i)
        
        self.abs_elim_3 = np.delete(abs_elim_3, delete)
        self.rising_0_elim_3 = np.delete(rising_0_elim_3, delete)
        self.rising_1_elim_3 = np.delete(rising_1_elim_3, delete)
        self.rising_2_elim_3 = np.delete(rising_2_elim_3, delete)
        self.rising_3_elim_3 = np.delete(rising_3_elim_3, delete)
        self.falling_0_elim_3 = np.delete(falling_0_elim_3, delete)
        self.falling_1_elim_3 = np.delete(falling_1_elim_3, delete)
        self.falling_2_elim_3 = np.delete(falling_2_elim_3, delete)
        self.falling_3_elim_3 = np.delete(falling_3_elim_3, delete)  
        
        f_r_percent_error = (abs_elim_3e.size - self.abs_elim_3.size)/(self.abs_elim_3.size) * 100
        print(f_r_percent_error)
        return (self.abs_elim_3,self.rising_0_elim_3, self.rising_1_elim_3, self.rising_2_elim_3, self.rising_3_elim_3,
                self.falling_0_elim_3, self.falling_1_elim_3, self.falling_2_elim_3, self.falling_3_elim_3, 
                self.POCAM_num, values, self.atstamp, self.p_jumps, dt_mean, f_r_percent_error, values,
                self.file_path, self.SDOM_num, self.PMT)