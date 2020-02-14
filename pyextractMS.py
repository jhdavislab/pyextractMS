###################################
#Library of functions to extract and 
#plot data from Thermo RAW files
#relies on MSFileReader python 
#bindings from Francois Allen
#Joey Davis
#jhdavislab.com
#jhdavis@mit.edu
###################################

import MSFileReader
import numpy as np
from matplotlib import pylab as plt
from matplotlib import cm
from scipy import interpolate
import pysoquant_report_defs_v1 as defs

__VERSION__='0.1.2'

###########################################
#######Thermo Raw file utilities###########
###########################################
def get_all_scan_numbers(rawfile, ms_level=1, rt_range=None):
    '''Get scan numbers from rawfile at a given ms_level

    Args:
        rawfile (string to thermo rawfile): string pointing to a Thermo rawfile or on opened rawfile
        ms_level (int): MS level of interest
        rt_range (list of floats): a range of retention times to consider
    
    Returns:
        numpy array with the relevant scan numbers
        
    Usage:
        all_ms1_scans = get_all_scan_numbers(rawfile, ms_level=1)
    '''
    
    if type(rawfile) is str:
        rawfile = MSFileReader.ThermoRawfile(rawfile)
    
    if rt_range is None:
        start_scan = rawfile.FirstSpectrumNumber
        end_scan = rawfile.LastSpectrumNumber
    else:
        try:
            start_scan = rawfile.ScanNumFromRT(rt_range[0])
            end_scan = rawfile.ScanNumFromRT(rt_range[1])
        except OSError:
            start_scan = rawfile.FirstSpectrumNumber
            end_scan = rawfile.LastSpectrumNumber
            start_RT = rawfile.RTFromScanNum(start_scan)
            end_RT = rawfile.RTFromScanNum(end_scan)
            print("provided retention time range (" + str(rt_range) + " is not valid for this .raw file. Reverting to use the full scan range (" + str([start_RT, end_RT]) + ").")
        
    scan_list = []
    for i in range(start_scan, end_scan+1):
        scan_type = rawfile.GetMSOrderForScanNum(i)
        if scan_type == ms_level:
            scan_list.append(i)
    
    return np.array(scan_list)

def get_ms2_scans(rawfile, precursor_mz, rt_range=None):
    '''Get scan numbers from rawfile that contain fragmented the precursor of interest

    Args:
        rawfile (string to thermo rawfile): string pointing to a Thermo rawfile or on opened rawfile
        precursor_mz (float): a float with the m/z of the precursor of interest
        rt_range (list of floats: Optional): a range of retention times to consider
    
    Returns:
        numpy array with the relevant scan numbers
        
    Usage:
        mz_ms2_scans = get_ms2_scans(rawfile, 695.8324, rt_range=[40,52], )
    '''

    if type(rawfile) is str:
        rawfile = MSFileReader.ThermoRawfile(rawfile)
    
    allms2 = get_all_scan_numbers(rawfile, ms_level=2, rt_range=rt_range)
    filled_scans = []
    for scan_num in allms2:
        current_trail = rawfile.GetTrailerExtraForScanNum(scan_num)
        mono_mz = rawfile.GetPrecursorMassForScanNum(scan_num, 2)
        current_mz_range_min = mono_mz - current_trail['MS2 Isolation Width']/2.0
        current_mz_range_max = mono_mz + current_trail['MS2 Isolation Width']/2.0
        if (precursor_mz >= current_mz_range_min) and (precursor_mz <= current_mz_range_max):
            filled_scans.append(scan_num)
    return np.array(filled_scans)

def get_xic(rawfile, rt_range, mz_range, scan_filter = "Full ms "):
    '''Extract ion chromatograms from a Thermo RawFile given retention time range and m/z range

        Args:
            rawfile (string to thermo rawfile): string pointing to a Thermo rawfile or on opened rawfile
            rt_range (list of floats): List of the start[0] and end[1] retention times - typically in minutes
            mz_range (list of floats): List of the m/z_start_range[0] and m/z_end_range[1]
            scan_filter (string: optional): String with the types of scans to inspect for the xic. Default to "Full ms ", which is MS1.

        Returns:
            numpy array with tuples of [retention_time, ion_intensity]
        
        Usage:
            example_xic = get_xic('0.raw', [50,60], [666.6, 666.7])
            example_xic_direct = get_xic(open_rawfile, [50,60], [666.6, 666.7])
    '''
    if type(rawfile) is str:
        rawfile = MSFileReader.ThermoRawfile(rawfile)
    
    try:
        time_array, intensity_array = rawfile.GetChroData(startTime=rt_range[0],
                                       endTime=rt_range[1],
                                       massRange1="{}-{}".format(mz_range[0], mz_range[1]),
                                       scanFilter="Full ms ")[0]
    except OSError:
        start_scan = rawfile.FirstSpectrumNumber
        end_scan = rawfile.LastSpectrumNumber
        start_RT = rawfile.RTFromScanNum(start_scan)
        end_RT = rawfile.RTFromScanNum(end_scan)
        
        print("provided retention time range (" + str(rt_range) + ") is not valid for this .raw file for the XICs. Reverting to use the full scan range (" + str([start_RT,end_RT])+").")
        time_array, intensity_array = rawfile.GetChroData(startTime=start_RT, endTime=end_RT,
                                                          massRange1="{}-{}".format(mz_range[0], mz_range[1]),
                                                          scanFilter="Full ms ")[0]
    return np.array(list(zip(time_array, intensity_array)))

def extract_ms_spectra(rawfile, scan_num, mz_range=[-1,10000]):
    '''Extract a mass spectra (m/z vs intensity) from a rawfile at a given scan number

    Args:
        rawfile (string to thermo rawfile): string pointing to a Thermo rawfile or on opened rawfile    
        scan_num (float): the scan number of the desired scan
        mz_range (tuple of floats): mz_range[0] is mz_start, mz_range[1] is mz_end (inclusive)
    
    Returns:
        numpy array with tuples of (m/z, intensity)
        
    Usage:
        ms_spectra = extract_ms_spectra(rawfile, 5555, [443.7,443.85])
    '''
    
    if type(rawfile) is str:
        rawfile = MSFileReader.ThermoRawfile(rawfile)

    mz_array, int_array = rawfile.GetMassListFromScanNum(scan_num)[0]
    masslist = np.column_stack((mz_array, int_array))
    return masslist[(masslist[:, 0] > mz_range[0]) & (masslist[:, 0] < mz_range[1])] #select only the values within the desired mz range

def interp_ms_spectra(ms_spectra, mz_range, sig_dec=3):
    '''Interpolate a mass spectra between a mz_range using a given number of significant digits.
    
    Args:
        ms_spectra (2D numpy array): should be the output of extract_ms_spectra with two columns (m/z, intensity)
        mz_range (tuple of floats): mz_range[0] is mz_start, mz_range[1] is mz_end (inclusive). If larger than the input spectra, will be padded with zeros.
        sig_dec (int: optoinal - default=3): specifies the number of significant decimal places to consider when summing spectra
    
    Returns:
        numpy array with tuples of (m/z, intensity)
        
    Usage:
        ms_spectra = interp_ms_spectra(ms_spectra, [400,800], 3)
    '''
    mz_interp_axis = np.arange(mz_range[0], mz_range[1], 1/10**sig_dec)
    if mz_range[0] < ms_spectra[0,0]:
        ms_spectra_full = np.concatenate((np.array([[mz_range[0], 0.0]]), ms_spectra), axis=0)
    if mz_range[1] > ms_spectra_full[-1,0]:
        ms_spectra_full = np.concatenate((ms_spectra_full, np.array([[mz_range[1], 0.0]])), axis=0)
    interp_func = interpolate.interp1d(ms_spectra_full[:,0], ms_spectra_full[:,1])
    return np.column_stack((mz_interp_axis, interp_func(mz_interp_axis)))

def sum_ms_spectra(rawfile, scan_list, mz_range, sig_dec=3):
    '''Sum ms spectra (m/z vs intensity). Sum a series of scans over a given mz_range. Spectra will first be interpolated 
    using the given sig_dec values so they can be more readily summed

    Args:
        rawfile (string to thermo rawfile): string pointing to a Thermo rawfile or on opened rawfile.
        scan_list (list of ints): scan numbers to sum.
        mz_range (tuple of floats): mz_range[0] is mz_start, mz_range[1] is mz_end (inclusive)
        sig_dec (int: optional - default=3): specifies the number of significant digits to consider when summing spectra
    
    Returns:
        numpy array with tuples of (m/z, intensity)
        
    Usage:
        summed_ms_spectra = sum_ms_spectra(rawfile, [500, 505, 510, 515], [400,1200], sig_dec=3)
    '''
    if type(rawfile) is str:
        rawfile = MSFileReader.ThermoRawfile(rawfile)

    assert mz_range[1]-mz_range[0] > 1/10**sig_dec, 'Please pass a mz_range as follows [low_value, high_value] where high_value-low_value > sig_dec'
    
    summed_mz = np.arange(mz_range[0], mz_range[1], 1/10**sig_dec)
    summed_int = np.zeros(summed_mz.shape)
    for scan in scan_list:
        spectra = extract_ms_spectra(rawfile, scan, mz_range)
        try:
            int_spectra = interp_ms_spectra(spectra, mz_range, sig_dec=sig_dec)
            summed_int+=int_spectra[:,1]
        except IndexError:
            pass
    return np.column_stack((summed_mz, summed_int))

###########################################
#######_ms_ting utilities###
###########################################
def plot_tic(rawfile, fig_axis=None, rt_range=None, color='black'):
    '''Plot total ion chromatograms from a Thermo RawFile 

        Args:
            rawfile (string to thermo rawfile): string pointing to a Thermo rawfile or on opened rawfile
            fig_axis (matplotlib figure axis: optional): if provided, data will be plotted on this axis,
                                                       otherwise a new figure axis is created.
            rt_range (list of floats: optional): List of the start[0] and end[1] retention times - typically in minutes
            color (string: optional): color of the plotted data. default='black'

        Returns:
            matplotlib figure axis
        
        Usage:
            tic_figure = plot_tic('0.raw', rt_range=[20,80])
    '''
    #Calc TIC
    if type(rawfile) is str:
        rawfile = MSFileReader.ThermoRawfile(rawfile)
    
    ms1_scans = get_all_scan_numbers(rawfile, ms_level=1, rt_range=rt_range)
    total_ion_current = []
    rt_times = []
    total_tic = 0.0
    for scan_num in ms1_scans:
        current_scan = rawfile.GetScanHeaderInfoForScanNum(scan_num)
        rt_times.append(current_scan['StartTime'])
        total_ion_current.append(current_scan['TIC'])
        total_tic+=current_scan['TIC']
            
    total_ion_current = np.array(total_ion_current)
    rt_times = np.array(rt_times)

    #Plot TIC
    if fig_axis is None:
        fig, fig_axis = plt.subplots(1)

    fig_axis.plot(rt_times, total_ion_current, color=color, linewidth=.2)
    fig_axis.set_xlabel('RT (mins)')
    fig_axis.set_ylabel('TIC')
    fig_axis.text(0.1, 0.9, 'total ion current = ' + '{:.2e}'.format(total_tic), transform=fig_axis.transAxes)
    fig_axis.grid()
    return fig_axis

def plot_bpc(rawfile, fig_axis=None, rt_range=None, color='black'):
    '''Plot base peakchromatograms from a Thermo RawFile 

        Args:
            rawfile (string to thermo rawfile): string pointing to a Thermo rawfile or on opened rawfile
            fig_axis (matplotlib figure axis: optional): if provided, data will be plotted on this axis,
                                                       otherwise a new figure axis is created.
            rt_range (list of floats: optional): List of the start[0] and end[1] retention times - typically in minutes
            color (string: optional): color of the plotted data. default='black'

        Returns:
            matplotlib figure axis
        
        Usage:
            bpc_figure = plot_bpc('0.raw', rt_range=[20,80])
    '''
    #Calc BPC
    if type(rawfile) is str:
        rawfile = MSFileReader.ThermoRawfile(rawfile)
        
    ms1_scans = get_all_scan_numbers(rawfile, ms_level=1, rt_range=rt_range)
    base_peak_current = []
    rt_times = []
    total_bpc = 0.0
    for scan_num in ms1_scans:
        current_scan = rawfile.GetScanHeaderInfoForScanNum(scan_num)
        rt_times.append(current_scan['StartTime'])
        base_peak_current.append(current_scan['BasePeakIntensity'])
        total_bpc+=current_scan['BasePeakIntensity']
    
    base_peak_current = np.array(base_peak_current)    
    rt_times = np.array(rt_times)

    #Plot BPC
    if fig_axis is None:
        fig, fig_axis = plt.subplots(1)
    fig_axis.plot(rt_times, base_peak_current, color=color, linewidth=.2)
    fig_axis.set_xlabel('RT (mins)')
    fig_axis.set_ylabel('BPC')
    fig_axis.text(0.1, 0.9, 'total base current = ' + '{:.2e}'.format(total_bpc), transform=fig_axis.transAxes)
    fig_axis.grid()
    return fig_axis

def plot_xic(xic, title='XIC', fig_axis=None, figsize=(5,5), color='black', x_label='retention time', y_label='intensity', circle_size=2):
    '''Plot extracted ion chromatograms - input should be output of get_xic()

        Args:
            xic (numpy array with tuples of RT/intensity): xic output is expected input
            title (string): plot title
            fig_axis (matplotlib figure axis: optional): if provided, data will be plotted on this axis,
                                                       otherwise a new figure axis is created.
            figsize (tuple of floats: optional): figure size (x_dim, y_dim). Default = (5,5)
            color (string: optional): color of the plotted data. default='black'
            x_label (string: optional): x axis label. default='retention time'
            y_label (string: optional): y axis label. default='intensity'
            circle_size (int: optional): size of the circles in the plotted data. default=2
            
        Returns:
            matplotlib figure axis
        
        Usage:
            new_fig_axis = plot_xic(example_xic, title='test_title', figsize=(8,8))
    '''
    plt.rcParams['grid.alpha'] = 0.5
    plt.rcParams['grid.color'] = "grey"
    if fig_axis is None:
        fig = plt.figure(figsize=figsize)
        fig_axis = fig.add_subplot(111)
    
    fig_axis.plot(xic[:,0], xic[:,1], 'o-', color=color, markersize=circle_size)
    fig_axis.set_xlabel(x_label)
    fig_axis.set_ylabel(y_label)
    fig_axis.set_title(title)
    fig_axis.grid()
    return fig_axis

def plot_cycle_time(rawfile, fig_axis=None, rt_range=None, color='black'):
    '''Plot ms1-to-ms1 cycle time

        Args:
            rawfile (string to thermo rawfile): string pointing to a Thermo rawfile or on opened rawfile
            fig_axis (matplotlib figure axis: optional): if provided, data will be plotted on this axis,
                                                       otherwise a new figure axis is created.
            rt_range (tuple of floats): the start ([0]) and end([1]) range to extract (nearest scans will be taken)
            color (string: optional): color of the plotted data. default='black'
            
        Returns:
            matplotlib figure axis
        
        Usage:
            new_fig_axis = plot_cycle_time(file_name)
    '''
    ###Calc result##
    if type(rawfile) is str:
        rawfile = MSFileReader.ThermoRawfile(rawfile)
        
    ms1_scans = get_all_scan_numbers(rawfile, rt_range=rt_range, ms_level=1)
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
    
    ###Plot result##
    if fig_axis is None:
        fig, fig_axis = plt.subplots(1)

    fig_axis.scatter(rt_times, cycle_times, alpha=0.2, s=1, color=color)
    fig_axis.set_xlabel('RT (mins)')
    fig_axis.set_ylabel('cycles time (secs)')
    fig_axis.grid()
    return fig_axis

def plot_inj_times(rawfile, ms_level, fig_axis=None, rt_range=None, num_bins=50, color='grey'):
    '''Plot scan injection time histogram

        Args:
            rawfile (string to thermo rawfile): string pointing to a Thermo rawfile or on opened rawfile
            ms_level (int): chose which MS level to analyze (e.g. 1, 2)
            fig_axis (matplotlib figure axis: optional): if provided, data will be plotted on this axis,
                                                       otherwise a new figure axis is created.
            rt_range (tuple of floats): the start ([0]) and end([1]) range to extract (nearest scans will be taken)
            num_bins (int): number of bins in the histogram
            color (string: optional): color of the plotted data. default='grey'
            
        Returns:
            matplotlib figure axis
        
        Usage:
            new_fig_axis = plot_inj_times(file_name, 1)
    '''
    if type(rawfile) is str:
        rawfile = MSFileReader.ThermoRawfile(rawfile)

    if fig_axis is None:
        fig, fig_axis = plt.subplots(1)    
        
    scans = get_all_scan_numbers(rawfile, ms_level=ms_level, rt_range=rt_range)
    inj_times = [rawfile.GetTrailerExtraForScanNum(scan_num)['Ion Injection Time (ms)'] for scan_num in scans]
    fig_axis.hist(inj_times, bins=num_bins, color=color)
    fig_axis.set_xlabel('injection time (ms)')
    fig_axis.set_ylabel('# scans')
    if rt_range is None:
        rt_range_string = 'all'
    else:
        rt_range_string = str(rt_range[0])+'-'+str(rt_range[1])
    fig_axis.text(0.1, 0.9, 'MS'+str(ms_level)+' Inj. Times | RT='+rt_range_string, transform=fig_axis.transAxes)

def plot_ms_spectra(ms_spectra, title='mass spectra', fig_axis=None, figsize=(5,5), 
                    color='black', x_label='mass/charge', y_label='intensity', circle_size=2,
                    linewidth=0.2):
    '''Plot mass spectra - input should be ouptut of extract_ms1_spectra or extract_ms2_spectra

        Args:
            ms_spectra (numpy array with tuples of m/z, intensity): extract_ms1 or extract_ms2 output is expected input
            title (string: optional): plot title. default='mass spectra'
            fig_axis (matplotlib figure axis: optional): if provided, data will be plotted on this axis, otherwise a new figure axis is created.
            figsize (tuple of floats: optional): figure size (x_dim, y_dim). Default = (5,5)
            color (string: optional): color of the plotted data. default='black'
            x_label (string: optional): x axis label. default='mass/charge'
            y_label (string: optional): y axis label. default='intensity'
            circle_size (int: optional): size of the circles in the plotted data. default=2
            linewidth (float: optional): size of the line in the plotted data. default=0.2
            
            
        Returns:
            matplotlib figure axis
        
        Usage:
            ax = plot_ms_spectra(ms_spectra)
    '''
    
    plt.rcParams['grid.alpha'] = 0.5
    plt.rcParams['grid.color'] = "grey"
    if fig_axis is None:
        fig = plt.figure(figsize=figsize)
        fig_axis = fig.add_subplot(111)
    
    fig_axis.plot(*zip(*ms_spectra),'o-', color=color, markersize=circle_size, linewidth=linewidth)
    fig_axis.set_xlabel(x_label)
    fig_axis.set_ylabel(y_label)
    fig_axis.set_title(title)
    fig_axis.grid()
    return fig_axis

def plot_xics(rawfile, fig_axis=None, colors=cm.get_cmap(name='plasma'), mass_accuracy=0.01, rt_window=7.5, 
              peptides={'DIPVPKPK':[451.2835,43.5], \
                    'IGDYAGIK':[422.7364,47.5], \
                    'TASEFDSAIAQDK':[695.8324,53.5], \
                    'SAAGAFGPELSR':[586.8003,56.5], \
                    'SFANQPLEVVYSK':[745.3925,67.5], \
                    'ELASGLSFPVGFK':[680.3736,81.5], \
                    'LSSEAPALFQFDLK':[787.4212,85.5]}):

    '''Plot a series of extracted ion chromatograms from a single injection - can be used to inspect standard peptides

        Args:
            rawfile (string to thermo rawfile): string pointing to a Thermo rawfile or on opened rawfile
            fig_axis (matplotlib figure axis: optional): if provided, data will be plotted on this axis,
                                                       otherwise a new figure axis is created.
            colors (matplotlib color map: optional): color map to use (defaults to plasma)
            mass_accuracy (float: optional): m/z width to extract around given exact m/z
            rt_window (float: optional): retention time range to highlight/zoom around given RT. Default=5
            peptide (hash of pepties: optional): has with keys of peptide sequence, and values of tuples of [exact_mass, retention_time]. Defaults to pierce iRTs
            color (string: optional): color of the plotted data. default='black'
            
        Returns:
            matplotlib figure axis
        
        Usage:
            new_fig_axis = plot_xics(file_name)
    '''
    if type(rawfile) is str:
        rawfile = MSFileReader.ThermoRawfile(rawfile)
    
    if fig_axis is None:
        fig, fig_axis = plt.subplots(1)
        
    total_peptides = len(peptides.keys())
    all_rts = []
    for peptide_num, peptide in enumerate(peptides.keys()):
        m_z = peptides[peptide][0]
        rt = peptides[peptide][1]
        color = colors(peptide_num/total_peptides)
        full_xic = get_xic(mz_range=[m_z-mass_accuracy,m_z+mass_accuracy], rawfile=rawfile,
                                   rt_range=[rt-rt_window, rt+rt_window])
        peak_max_rt, peak_max_intensity = get_max_rt(full_xic)
        all_rts.append(peak_max_rt)
    
        target_xic = get_xic(mz_range=[m_z-mass_accuracy,m_z+mass_accuracy], rawfile=rawfile, 
                                    rt_range=[peak_max_rt-rt_window/4, peak_max_rt+rt_window/4])
        
        fig_axis.plot(target_xic[:,0], target_xic[:,1], color=color, linewidth=0.75, label=peptide)
        fig_axis.axvline(peak_max_rt, linestyle='dotted', color=color, alpha=0.75, linewidth=1)
    fig_axis.set_xlabel('retention time (mins)')
    fig_axis.set_ylabel('intensity')
    all_rts_np = np.array(all_rts)
    xmin = all_rts_np.min()-rt_window
    xmax= all_rts_np.max()+((all_rts_np.max()-all_rts_np.min())/2)
    fig_axis.set_xlim(xmin,xmax)
    fig_axis.text(0.1, 0.9, 'XICs', transform=fig_axis.transAxes)
    fig_axis.legend(loc='center left', bbox_to_anchor=(0.75, 0.5))
    return fig_axis

def plot_pressure_traces(rawfile, fig_axis=None, rt_range=None, colors=['red', 'blue']):
    '''Plot a series of extracted ion chromatograms from a single injection - can be used to inspect standard peptides
    Args:
        rawfile (string to thermo rawfile): string pointing to a Thermo rawfile or on opened rawfile
        fig_axis (matplotlib figure axis: optional): if provided, data will be plotted on this axis,
                                                       otherwise a new figure axis is created.
        colors (list of strings): colors to plot for [0]=gradient_pump; [1]:loading_pump (defaults to red, blue)
        rt_range (list of floats: optional): retention time range to highlight/zoom around given RT. Default=full range
            
        Returns:
            matplotlib figure axis
        
        Usage:
            new_fig_axis = plot_pressure_traces(file_name)
    '''
    if type(rawfile) is str:
        rawfile = MSFileReader.ThermoRawfile(rawfile)
    
    if fig_axis is None:
        fig, fig_axis = plt.subplots(1)

    rawfile.SetCurrentController('A/D card', 1)
    loading_pump_chro_data = np.array(rawfile.GetChroData())[0]
    rawfile.SetCurrentController('A/D card', 2)
    gradient_pump_chro_data = np.array(rawfile.GetChroData())[0]
    rawfile.SetCurrentController('MS', 1)
    fig_axis.plot(gradient_pump_chro_data[0], gradient_pump_chro_data[1], color=colors[0], linewidth=.5)
    fig_axis.set_xlabel('RT (mins)')
    fig_axis.set_ylabel('grad pressure', color=colors[0])
    fig_axis2 = fig_axis.twinx()
    fig_axis2.plot(loading_pump_chro_data[0], loading_pump_chro_data[1], color=colors[1], linewidth=.5)
    fig_axis2.set_ylabel('load pressure', color=colors[1])
    fig_axis.grid()
    return fig_axis
    
###########################################
#######  Helper functions  ################
###########################################

def calc_mz_range(peptide_seq, mz, labeling='N15', charge=2, offset=2):
    num_nit=num_nit=np.array([defs.AANITROGENS[i] for i in peptide_seq]).sum()
    num_car=num_nit=np.array([defs.AANITROGENS[i] for i in peptide_seq]).sum() #NEED TO CHANGE FOR NUMBER OF CARBONS
    num_lys=peptide_seq.upper().count('K')
    num_arg=peptide_seq.upper().count('R')
    
    if labeling=='K8R10':
        return [mz-1, mz+num_lys*8/charge+num_arg*10/charge+offset]
    if labeling=='K6R4':
        return [mz-1, mz+num_lys*6/charge+num_arg*4/charge+offset]
    if labeling=='N15':
        return [mz-1, mz+num_nit/2+offset]
    if labeling=='C13':
        return [mz-1, mz+num_car/2+offset]
    if labeling=='None':
        return [mz-1, mz+offset]
    else:
        raise Exception('labeling type: "'+labeling + '" not defined, please add this definition to pyextractMS if you want to use it') 

def find_nearest(array, value):
    '''Helper function to find nearest value in an array

    Args:
        array (numpy array): array with values to search
        value (float or int): value to search for
    
    Returns:
        the value in the array closest to the value passed
        
    Usage:
        scan_num = find_nearest(ms1_scans, rawfile.ScanNumFromRT(retention_time))
    '''
    idx = (np.abs(array - value)).argmin()
    return array[idx]

def get_max_rt(xic, rt_range=None):
    '''get_max_rt finds the maximum peak location in a provided xic, which should be the output of get_xic

    Args:
        xic (numpy array with tuples of [retention_time, ion_intensity]): Should be the output of get_xic
        rt_range (list of floats: optional): List of start[0] and end[1] retention times, which can be used to limit the search space.

    Returns:
        tuple of floats [0]=retention time of max intensity; [1]=max intensity
        
    Usage:
        example_xic = get_xic('0.raw', [50,60], [666.6, 666.7])
        peak_max_rt, peak_max_intensity = get_max_rt(example_xic)
    '''
    if rt_range is None:
        rt_range_idx = [0, -1]
    else:
        rt_range_idx = [(np.abs(xic[:,0] - rt_range[0])).argmin(), (np.abs(xic[:,0] - rt_range[1])).argmin()]
    
    peak_max_idx = np.argmax(xic[rt_range_idx[0]:rt_range_idx[1],1])
    return xic[peak_max_idx]

##########################################
############ Example Usage ###############
##########################################

if __name__ == "__main__":
    rawfile = MSFileReader.ThermoRawfile('0.raw')

    #plot xic
    rt_range = [46,50]
    mz_range = [695.7, 695.9]
    example_xic = get_xic(rawfile, rt_range, mz_range)
    xic_axis = plot_xic(example_xic, title='xic_'+str(mz_range), figsize=(8,8))

    #plot single spectra
    ms1_peptide = 'TASEFDSAIAQDK'
    ms1_charge = 2
    ms1_mz = 695.8324
    ms1_rt = 46.75
    ms1_mass_accuracy=0.01
    ms1_rt_window=5
    full_xic = get_xic(rawfile, [ms1_rt-ms1_rt_window, ms1_rt+ms1_rt_window], [ms1_mz-ms1_mass_accuracy,ms1_mz+ms1_mass_accuracy])
    ms1_peak_max_rt, ms1_peak_max_intensity = get_max_rt(full_xic)
    ms1_peak_max_scan = rawfile.ScanNumFromRT(ms1_peak_max_rt)
    ms1_mz_range = calc_mz_range(ms1_peptide, ms1_mz, labeling='None', charge=2, offset=2)
    ms1 = extract_ms_spectra(rawfile, ms1_peak_max_scan, ms1_mz_range)
    single_ax = plot_ms_spectra(ms1, title='Single spectra: '+str(ms1_peak_max_rt))
    single_ax.axvline(ms1_mz, linestyle='--', color='red', alpha=0.75, linewidth=1)

    #plot interpolated spectra
    ms1_interp = interp_ms_spectra(ms1, ms1_mz_range, sig_dec=3)
    interp_ax = plot_ms_spectra(ms1_interp, title='Interp spectra: '+str(ms1_peak_max_rt))
    interp_ax.axvline(ms1_mz, linestyle='--', color='red', alpha=0.75, linewidth=1)
    
    #plot integrated spectra
    rt_range = [47.45, 47.75]
    ms1_scans = get_all_scan_numbers(rawfile, ms_level=1, rt_range=rt_range)
    ms_summed_spectra = sum_ms_spectra(rawfile, ms1_scans, ms1_mz_range, sig_dec=3)
    summed_ax = plot_ms_spectra(ms_summed_spectra, title='summed spectra: '+str(rt_range[0])+'-'+str(rt_range[1]))
    summed_ax.axvline(ms1_mz, linestyle='--', color='red', alpha=0.75, linewidth=1)

    xic_axis.axvline(ms1_peak_max_rt, linestyle='--', color='red', alpha=0.75, linewidth=1)
    xic_axis.axvline(rt_range[0], linestyle='--', color='blue', alpha=0.75, linewidth=1)
    xic_axis.axvline(rt_range[1], linestyle='--', color='blue', alpha=0.75, linewidth=1)
    plt.show()