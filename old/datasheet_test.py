'''
def plot_cycle_times(list_of_files):
    if list_of_files is None:
        list_of_files = glob.glob('./*.raw')
    
    num_files = len(list_of_files)
    num_pages = num_files//15
    
    for page in num_pages:
        fig, fig_axes = plt.subplots(min(15, num_files))
        for file_name, file_index in list_of_file[page*15:(page+1)*15]:
            rawfile = MSFileReader.ThermoRawfile(file_name)
            ms1_scans = get_all_scan_numbers(file_name, ms_level=1)
            fig_axes[file_index] = calc_cycle_time(rawfile, ms1_scans, fig_axis=fig_axes[file_index])
            fig_axes[file_index].set_title(file_name)

def calc_cycle_time(rawfile, ms1_scans, fig_axis=None):
    if fig_axis is None:
        fig, fig_axis = plt.subplots(1)
    cycle_times = []
    rt_times = []
    prior_time = rawfile.RTFromScanNum(ms1_scans[0])
    for scan_num in ms1_scans[1:]:
        current_time = rawfile.RTFromScanNum(scan_num)
        rt_times.append(current_time)
        cycle_times.append((current_time - prior_time)*60.0)
        prior_time = current_time
    cycle_times = np.array(cycle_times)
    rt_times = np.array(rt_times)
    fig_axis.plot(rt_times, cycle_times)
    fig_axis.set_xlabel('retention time (minutes)')
    fig_axis.set_ylabel('cycles time (seconds)')
    return fig_axis

def plot_tic(rawfile, ms1_scans, fig_axis=None):
    if fig_axis is None:
        fig, fig_axes = plt.subplots(2)

    total_ion_current = []
    base_peak_current = []
    rt_times = []
    total_tic = 0.0
    total_bpc = 0.0
    for scan_num in ms1_scans:
        current_scan = rawfile.GetScanHeaderInfoForScanNum(scan_num)
        rt_times.append(current_scan['StartTime'])
        total_ion_current.append(current_scan['TIC'])
        total_tic+=current_scan['TIC']
        
        base_peak_current.append(current_scan['BasePeakIntensity'])
        total_bpc+=current_scan['BasePeakIntensity']
    
    total_ion_current = np.array(total_ion_current)
    base_peak_current = np.array(base_peak_current)    
    rt_times = np.array(rt_times)
    
    fig_axes[0].plot(rt_times, total_ion_current)
    fig_axes[1].plot(rt_times, base_peak_current, 'g-')
    #fig_axes[0].set_xlabel('retention time (minutes)')
    fig_axes[1].set_xlabel('retention time (minutes)')
    fig_axes[0].set_ylabel('total ion current')
    fig_axes[1].set_ylabel('base peak intensity')
    fig_axes[0].set_title('total ion current = ' + '{:.2e}'.format(total_tic))
    fig_axes[1].set_title('total base peak = ' + '{:.2e}'.format(total_bpc))
    #plt.tight_layout()
    return fig_axes, total_tic, total_bpc
    
def get_tics(rawfile, scan_list):
    total_tic = 0.0
    total_bpc = 0.0
    for scan_num in scan_list:
        current_scan = rawfile.GetScanHeaderInfoForScanNum(scan_num)
        total_tic+=current_scan['TIC']
        total_bpc+=current_scan['BasePeakIntensity']
    return total_tic, total_bpc

def get_total_ion_currents(list_of_files, ms_level=1):
    entries = {}
    for file_name in list_of_files:
        
        rawfile = MSFileReader.ThermoRawfile(file_name)
        ms_scans = get_all_scan_numbers(rawfile, ms_level=ms_level)
        total_tic, total_bpc = get_tics(rawfile, ms_scans)
        entries[file_name] = {'total_tic': '{:.3e}'.format(total_tic), 'total_bpc':'{:.3e}'.format(total_bpc)}
    return entries

def plot_ms2_dia_tics(rawfile):
    all_scan_types = rawfile.GetFilters()
    all_ms2_scan_types = [i.split('@')[0].split(' ')[-1] for i in all_scan_types if "ms2" in i]
''' 