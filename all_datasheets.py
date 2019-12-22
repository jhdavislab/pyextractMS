import argparse
import datasheet
import glob
from matplotlib import cm

if __name__ =='__main__':
    parser = argparse.ArgumentParser(description='Plot a single datasheet for a each raw file in a given directory.')
    parser.add_argument('directory', type=str,
                       help='path to the directory to analyze')
    parser.add_argument('--display', default=True, action='store_false',
                        help='just display the plots, but do not save them')
    parser.add_argument('--extension', default='.pdf', type=str,
                        help='string for figure filetype (e.g. .pdf, .png)')
    args = parser.parse_args()

    directory= vars(args)['directory']
    savefig= vars(args)['display']
    fig_extension = vars(args)['extension']
    
    datasheet.set_style_single_inj()
    all_files = glob.glob(directory+'*.raw')
    print('analyzing the following files:')
    print(all_files)
    for file_name in all_files:
        datasheet.plot_datapage(file_name, savefig=savefig, fig_extension=fig_extension, colors=cm.get_cmap(name='plasma'))
    
