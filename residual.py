import numpy as np
import matplotlib.pylab as plt
import re
from tables import open_file
from statistics import median, mode
from scipy.optimize import curve_fit, minimize
from math import acos, degrees, log
from matplotlib.colors import LogNorm

class residual(object):
    def __init__(self, abs_elim_3, rising_a_elim_3, lspace_lower, lspace_upper, lspace, lower_bound, upper_bound,
                xlim_min, xlim_max, ylim_min, ylim_max, file_path, values, POCAM_num):
    
        self.t_threshold = abs_elim_3[rising_a_elim_3 > 0] + rising_a_elim_3[rising_a_elim_3 > 0]
        t_diff = self.t_threshold[1:] - self.t_threshold[:-1]
        self.events_in_peak = self.t_threshold[:-1][(t_diff[:] > lower_bound) & (t_diff[:] < upper_bound)] #selecting events in peak
    
        plt.figure(figsize=(10,9))
        _ = plt.hist(t_diff, np.linspace(lspace_lower,lspace_upper, lspace), log = True)
        plt.title(values + '- Time Difference', fontsize = 19)
        plt.xlabel('time difference(ns)', fontsize = 16)

        if POCAM_num ==['P1']:
            self.estimate_peak = 200100.71     #for POCAM_1
        if POCAM_num == ['P2']:
            self.estimate_peak = 200100.33     #for POCAM_2
    
        self.estimate_residual = self.events_in_peak%self.estimate_peak
    
        plt.figure(figsize=(10,9))
        plt.plot(self.events_in_peak, self.estimate_residual, '.')
        plt.title(values + '-time residuals of threshold 1', fontsize = 19)
        plt.xlabel('timestamps', fontsize = 16)
        plt.ylabel('time residuals', fontsize = 16)
        plt.xlim(xlim_min, xlim_max)
        plt.ylim(ylim_min,ylim_max)
        
        self.abs_elim_3 = abs_elim_3
        self.rising_1_elim_3 = rising_a_elim_3
        self.values = values
        self.file_path = file_path
        
    def minimizer(self, select_min, select_max):
        
        #Minimizer
        tmin = select_min
        tmax = select_max
        def dfunc2(delta_t):
            x = self.events_in_peak
            tres = x%delta_t
            selection = (self.events_in_peak>tmin)*(self.events_in_peak<tmax)
            return np.sum((tres[selection]-tres[selection].mean())**2)
        from scipy.optimize import minimize
        m = minimize(dfunc2, [200100.], method='Powell')
        #print(m)
        print(m.x)
        self.gaus_peak = m.x
        return self.gaus_peak
    
    def res(self, greater_than, med_bound):
        print(self.events_in_peak.size)
        est_res_diff = abs(self.estimate_residual[1:] - self.estimate_residual[:-1])
        jump_index = ([])
        est_res_diff_ji = ([])
        jump_index = np.append(jump_index, 0)
        for r in range(0, est_res_diff.size):
            if est_res_diff[r] > greater_than:
                jump_index = np.append(jump_index, r)
            
        JumpIndex = jump_index.astype(int)
            
        print('jump_index', jump_index)
        print('jump_index size', jump_index.size)
    
        self.peak_1 = [] 
        t_res_all = np.array([])
        t_res_all_all = np.array([])
        abs_elim = self.abs_elim_3[self.rising_1_elim_3 > 0] + self.rising_1_elim_3[self.rising_1_elim_3 > 0]
        
        for p in range(0, len(JumpIndex)):      

            if p == 0:
                print("index_p", p)

                v = self.events_in_peak[:][(self.events_in_peak[:] >= self.events_in_peak[0]) & 
                                      (self.events_in_peak[:] < self.events_in_peak[JumpIndex[p+1]+1])]
                print('Jump Indices',0,JumpIndex[p+1]+1)
                if v.size == 0:
                    continue
                print(v.size)
                b = self.events_in_peak[JumpIndex[1]+1]
                a = min(v)
                run_time = abs_elim[(abs_elim >= 0) & (abs_elim < b)]
            
            elif p == JumpIndex.size - 1:
                print("index_p", p)
                print('Jump Index size =', JumpIndex.size, self.events_in_peak.size - 1)
                print([JumpIndex[p]+1])
                v = self.events_in_peak[:][(self.events_in_peak[:] >= self.events_in_peak[JumpIndex[p]+1]) & 
                                      (self.events_in_peak[:] <= self.events_in_peak[self.events_in_peak.size - 1])]
            
                print('Jump Indices', JumpIndex[p]+1, self.events_in_peak.size - 1)
                a = min(v)
                run_time = abs_elim[(abs_elim >= a) & (abs_elim <= abs_elim[abs_elim.size-1])]
        
            else: 
                print("index_p", p)
                v = self.events_in_peak[:][(self.events_in_peak[:] >= self.events_in_peak[JumpIndex[p]+1]) & 
                                      (self.events_in_peak[:] < self.events_in_peak[JumpIndex[p+1]+1])]
                print('Jump Indices', JumpIndex[p]+1, JumpIndex[p+1]+1)
                print(v.size)
                b = self.events_in_peak[JumpIndex[p+1]+1]
                a = min(v)
                run_time = abs_elim[(abs_elim >= a) & (abs_elim < b)]
        
            time_res = v%self.gaus_peak
            time_res_all = run_time%self.gaus_peak
            #print('time_res length', time_res.size)
        
            plt.figure(figsize=(10,9))            #plotting time residual graph for individual runs
            plt.plot(time_res, '.')
            plt.ylabel('time_residual')
            plt.show()
    
            #Gaussian fit
            med = median(time_res)                    
            med_all = median(time_res_all)
            peak = time_res[(time_res >= med - med_bound) & (time_res <= med + med_bound)]
            peak_all = time_res_all[(time_res_all >= med_all - med_bound) & (time_res_all <= med_all + med_bound)]
    
            def gaussian(x, mean, amplitude, standard_deviation):
                return amplitude * np.exp( - ((x - mean) / standard_deviation) ** 2)

            bins = np.linspace(med-med_bound, med + med_bound, 11)
            bins_all = np.linspace(med-med_bound, med + med_bound, 11)
            data_entries_1, bins_1, _ = plt.hist(peak, bins, alpha = 0.5)
            #plt.show()
            data_entries_1_all, bins_1_all, _ = plt.hist(peak_all, bins_all, alpha = 0.5)
            #plt.show()
    
            data = peak
            data_all = peak_all
            bincenters = ((bins[:-1]+bins[1:])/2)
            bincenters_all = ((bins_all[:-1]+bins_all[1:])/2)
    
            from scipy.optimize import curve_fit
            data_entries = data_entries_1
            popt, pcov = curve_fit(gaussian, xdata = bincenters, 
                                    ydata = data_entries,  
                                    absolute_sigma = True, 
                                    p0 = (med, 10, 5),
                                    sigma = np.sqrt(data_entries))
            data_entries_all = data_entries_1_all
    
            popt_all, pcov_all = curve_fit(gaussian, xdata = bincenters_all, 
                                            ydata = data_entries_all,  
                                            absolute_sigma = True, 
                                            p0 = (med, 10, 5),
                                            sigma = np.sqrt(data_entries_all))
    
            time_res_sub = time_res - popt[0]
            time_res_sub_all = time_res_all - popt_all[0]
    
            self.peak_1.append(time_res_all) #peak_1 is for function CheckPeak

            #print(popt)
            #print(popt_all)
            t_res_all =np.append(t_res_all, time_res_sub)
            t_res_all_all =np.append(t_res_all_all, time_res_sub_all)
            
        plt.figure(figsize=(10,9))
        n, bins, patches = plt.hist(t_res_all, 480,
                                    #np.linspace(-10, 40, 100), 
                                    log = True)
        plt.title(self.values + '-time residuals of threshold 1', fontsize = 19)
        plt.xlabel('time_ns', fontsize = 16)
        plt.ylabel('bincount', fontsize = 16)
        plt.axvline(color = 'r')
    
        plt.figure(figsize=(10,9))
        n, bins, patches = plt.hist(t_res_all_all, 
                    #500,
                    np.linspace(-100, 200, 480), 
                    log = True)
        
        plt.title(self.values + '-time residuals of threshold 1', fontsize = 19)
        plt.xlabel('time_ns', fontsize = 16)
        plt.ylabel('bincount', fontsize = 16)
        plt.axvline(x = 0, color = 'r')
        plt.axvline(x = -10, color = 'k')
        plt.axvline(x = 10, color = 'k')
        plt.savefig(self.file_path + '/graphs/' + self.values + 'time resi graph', dpi = 200)
        
        plt.figure(figsize=(10,9))
        n, bins, patches = plt.hist(t_res_all_all, 
                                    #500,
                                    np.linspace(-1000, 20000, 480), 
                                    log = True)
        plt.figure(figsize=(10,9))
        n, bins, patches = plt.hist(t_res_all_all, 
                                    500,
                                    #np.linspace(-1000, 20000, 480), 
                                    log = True)
    
        plt.figure(figsize=(10,9))
        plt.axhline(t_res_all[(t_res_all<25000)].mean(), 0, 1, color='k')
        plt.plot(t_res_all, '.')
    
        what_peak = abs_elim[(t_res_all_all > 24) & (t_res_all_all < 40)]
        num_events = t_res_all_all[(t_res_all_all > -10) & (t_res_all_all < 10)]
        noise_events = t_res_all_all[(t_res_all_all < -10000) ^ (t_res_all_all > 10000)]
        return t_res_all_all, num_events.size, noise_events.size
    
    def HIST2D(self, BinsHist, SDOM_num):
    	plt.figure(figsize=(12,12))
    	self.abs_elim = self.abs_elim_3[self.rising_1_elim_3 > 0] + self.rising_1_elim_3[self.rising_1_elim_3 > 0]
    	x = self.abs_elim
    	y = self.abs_elim % self.gaus_peak
    	h, self.xedges, self.yedges, img = plt.hist2d(x, y,
                            				BinsHist,
                            				#[np.linspace(0.0e10,0.2e10, 150), np.linspace(97620, 102000, 150)],
                            				#[np.linspace(2.68e10, 3e10, 150), np.linspace(13410, 13500, 150)],
                            				#cmin = 4 , 
                            				#norm = LogNorm() 
                            				)
    	cb = plt.colorbar()
    
    	self.POCAM_bins = ([])
    	for j in range (0, BinsHist):
        	#print(j, j+1)
        	bins = h[j:j+1, 0:].flatten()
        	max_ind = np.argmax(bins)
        	self.POCAM_bins = np.append(self.POCAM_bins, max_ind)
        
    	POCAM_diff = abs(self.POCAM_bins[1:] - self.POCAM_bins[:-1])
    
    	if SDOM_num == ['SDOM5']:
        	Mode = mode(POCAM_diff[POCAM_diff > 1])
        	#Mode = 10
        	jump_index = np.where((POCAM_diff > Mode - 25) * (POCAM_diff < Mode + 25))
    	if SDOM_num == ['SDOM1']:
        	Mode = 1
        	jump_index = np.where((POCAM_diff > Mode - 5) * (POCAM_diff < Mode + 5))
    	if SDOM_num == ['SDOM2']:
        	Mode = 1
        	jump_index = np.where((POCAM_diff > Mode - 25) * (POCAM_diff < Mode + 25))
    	if SDOM_num == ['SDOM3']:
        	Mode = 1
        	jump_index = np.where((POCAM_diff > Mode - 25) * (POCAM_diff < Mode + 25))
        
    	self.JumpIndex = (np.array(jump_index).flatten()) + 1
    	print(self.JumpIndex)    
    
    	plt.plot(self.POCAM_bins, '.')
    	return self.abs_elim, BinsHist, self.JumpIndex, self.xedges, self.yedges, self.POCAM_bins, POCAM_diff
    		
    def calc_res(self, BinsHist, med_bound, yaxis_lbound, yaxis_ubound):
    	#gaus_peak = 200100.33417353552
    	a = 0
    	Min = 0
    	Max = self.JumpIndex[0]
    	t_res_all = np.array([])
    	self.peak_1 = []
    	for d in range(0, self.JumpIndex.size + 1):
        	print('run#', d)
        	if d == self.JumpIndex.size:
        		b = BinsHist
        		Max = BinsHist
        		print('b', b, 'Max', Max)
        		select = self.abs_elim[(self.abs_elim >= a) & (self.abs_elim <= self.xedges[b])] % self.gaus_peak
        		y_axis = int(self.POCAM_bins[Min:Max].mean())
        		print('yaxis - ', y_axis)
        		lower_bound = y_axis -  yaxis_lbound
        		upper_bound = y_axis + yaxis_ubound
        		time_res = select[(select >= self.yedges[lower_bound]) & (select <= self.yedges[upper_bound])]
        	else:
        		b = self.JumpIndex[d] 
        		Max = self.JumpIndex[d]
        		print('b', b, 'Max', Max)
        		select = self.abs_elim[(self.abs_elim >= a) & (self.abs_elim < self.xedges[b])] % self.gaus_peak
        		print('select size', select.size)
        		y_axis = int(self.POCAM_bins[Min:Max].mean())
        		lower_bound = y_axis -  yaxis_lbound
        		upper_bound = y_axis + yaxis_ubound
        		time_res = select[(select >= self.yedges[lower_bound]) & (select <= self.yedges[upper_bound])]
        		
        	if time_res.size == 0:
        		continue
        	else:
        		med = median(time_res)
        		print('median - ', med)
        		peak = time_res[(time_res >= med - med_bound) & (time_res <= med + med_bound)]
        		
        		def gaussian(x, mean, amplitude, standard_deviation):
        			return amplitude * np.exp( - ((x - mean) / standard_deviation) ** 2)
        			
        		bins = np.linspace(med - med_bound, med + med_bound, 11)
        		data_entries_1, bins_1, _ = plt.hist(peak, bins, alpha = 0.5)
        		
        		data = peak
        		bincenters = ((bins[:-1]+bins[1:])/2)
        		
        		from scipy.optimize import curve_fit
        		data_entries = data_entries_1
        		popt, pcov = curve_fit(gaussian, xdata = bincenters, 
        								ydata = data_entries,
        								absolute_sigma = True,
        								p0 = (med, 10, 5),
        								sigma = np.sqrt(data_entries))
        								
        		time_res_sub = select - popt[0]
        		
        		self.peak_1.append(time_res_sub)
        		
        		t_res_all =np.append(t_res_all, time_res_sub)
        		
        	a = self.xedges[b]
        	Min = Max
        	print('a', b, Min, 'min')
        
    	plt.figure(figsize=(10,9))
    	_ = plt.hist(t_res_all, 
                	#500,
                	np.linspace(-100, 200, 480), 
                	log = True)
    	plt.axvline(x = 0, color = 'r')
    	plt.axvline(x = -10, color = 'k')
    	plt.axvline(x = 10, color = 'k')
    
        
    	plt.figure(figsize=(10,9))
    	_ = plt.hist(t_res_all,
    				#500,
    				np.linspace(-1000, 20000, 480),
    				log = True)
    	plt.show()
    
    	plt.figure(figsize=(10,9))
    	_ = plt.hist(t_res_all, 
                	500,
                	#np.linspace(-100500, -99500, 480), 
                	log = True)
    	plt.show()
    
    	num_events = t_res_all[(t_res_all > -10) & (t_res_all < 10)]
    	noise_events = t_res_all[(t_res_all < -10000) ^ (t_res_all > 10000)]
    
    	if t_res_all.size != self.abs_elim.size:
        	print('TRUE')

    	return t_res_all,num_events.size, noise_events.size
    	
        
    def CheckPeak(self):
        plt.figure(figsize=(10,9))
        for b in range(0, len(self.peak_1)):
            plt.figure(figsize=(10,9))
            n, bins, patches = plt.hist(self.peak_1[b], 500, alpha = 0.4, log = True)
        plt.title(self.values + '-time residuals of threshold 1', fontsize = 19)
        plt.xlabel('time_ns', fontsize = 16)
        plt.ylabel('bincount', fontsize = 16)
    		