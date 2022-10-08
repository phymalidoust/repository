#!/usr/bin/env python

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import colors
import matplotlib.pyplot as plt
import numpy as np
from math import sqrt
import re
import h5py

# Create a single 1D line plot of a property
def plot_real_1D(   data,
                    scaling=1.0,
                    xlabel='Y',
                    ylabel='Z',
                    series=None,
                    title='Title'):
    m = len(data)
    n = len(data[0])
    x = scaling*(np.arange(n))

    fig = plt.figure()
    if(m > 1):
        if (series == None):
            series = ['Series ' + str(i) for i in range(m)]
        if(len(series) < m):
            for i in range(m - len(series)):
                series.append('Series ' + str(i))
            
        for i in range(m):
            plt.plot(x,data[i], label=series[i])

        plt.legend()
    else:
        plt.plot(x, data[0])

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    fig.suptitle(title)
    fig.set_size_inches(11, 8.5)
    return fig

# Create a single rectangular plot with real data
def plot_real_rect( data,
                    transpose=True,
                    scaling=(1.0,1.0),
                    xlabel='Y',
                    ylabel='Z',
                    title='Title'):
    values = data

    try:
        values = data[()]
    except AttributeError:
        values = data

    if(transpose):
        values = values.transpose()

    rows = len(values)
    cols = len(values[0])

# Plot the Density of States data
def plot_DOS( up, down, total, scaling, eps):
    up = up[()].transpose()
    down = down[()].transpose()
    total = total[()].transpose()

    rows = len(up)
    cols = len(up[0])

    fig, axs = plt.subplots(3, 1)
    fig.suptitle(f'Density of States, Epsilon = {eps}')

    images = []
    images.append(axs[0].imshow(up, cmap='plasma', origin='lower',
                                extent=[0,cols*scaling[1], 0, rows*scaling[0]]))
    axs[0].label_outer()
    axs[0].set_title('Up')
    axs[0].set_ylabel('Z')
    axs[0].set_xlabel('Y')

    images.append(axs[1].imshow(down, cmap='plasma', origin='lower',
                                extent=[0,cols*scaling[1], 0, rows*scaling[0]]))
    axs[1].label_outer()
    axs[1].set_title('Down')
    axs[1].set_ylabel('Z')
    axs[1].set_xlabel('Y')

    images.append(axs[2].imshow(total, cmap='plasma', origin='lower',
                                extent=[0,cols*scaling[1], 0, rows*scaling[0]]))
    axs[2].label_outer()
    axs[2].set_title('Total')
    axs[2].set_ylabel('Z')
    axs[2].set_xlabel('Y')

    mn = min(up.min(), down.min(), total.min())
    mx = max(up.max(), down.max(), total.max())
    norm = colors.Normalize(mn, mx)

    images[0].set_norm(norm)
    images[1].set_norm(norm)

    fig.colorbar(images[0], ax=axs, orientation='vertical', fraction=.1)
    fig.set_size_inches(11, 8.5)

    return fig

# Create a "single" rectangular plot of the domain with complex data
def plot_complex_rect(  data,
                        transpose=True,
                        scaling=(1.0,1.0),
                        xlabel='Y',
                        ylabel='Z',
                        title='Title'):
    re_values = data[0,:,:]
    im_values = data[1,:,:]
    if(transpose):
        re_values = re_values.transpose()
        im_values = im_values.transpose()

    
    (rows,cols) = re_values.shape
    print('complex, rows =', rows, 'cols =', cols)
    fig,axs = plt.subplots(2, 1)
    fig.suptitle(title)

    images = []
    images.append(axs[0].imshow(re_values, cmap='plasma',
                                origin='lower', extent=[0,cols*scaling[1], 0, rows*scaling[0]]))
    axs[0].label_outer()
    axs[0].set_title('Real')
    axs[0].set_ylabel(ylabel)
    axs[0].set_xlabel(xlabel)

    images.append(axs[1].imshow(im_values, cmap='plasma',
                                origin='lower', extent=[0,cols*scaling[1], 0, rows*scaling[0]]))
    axs[1].label_outer()
    axs[1].set_title('Imaginary')
    axs[1].set_ylabel(ylabel)
    axs[1].set_xlabel(xlabel)

    mn = min(re_values.min(), im_values.min())
    mx = max(re_values.max(), im_values.max())
    norm = colors.Normalize(mn, mx)

    images[0].set_norm(norm)
    images[1].set_norm(norm)

    fig.colorbar(images[0], ax=axs, orientation='vertical', fraction=.1)
    fig.set_size_inches(11, 8.5)

    return fig

# Create a series of rectangular domain plots for each component of a property
def plot_real_xyz_rects(data,
                        transpose=True,
                        scaling=(1.0,1.0),
                        xlabel='Y',
                        ylabel='Z',
                        title='Title',
                        subtitles=['X Component', 'Y Component', 'Z Component']):

    x = data[0]
    y = data[1]
    z = data[2]

    if(transpose):
        try:
            x = x.transpose()
            y = y.transpose()
            z = z.transpose()
        except AttributeError:
            x = x[()].transpose()
            y = y[()].transpose()
            z = z[()].transpose()

    (rows,cols) = x.shape
    print('xyz, rows =', rows, 'cols =', cols)

    fig,axs = plt.subplots(2, 2)
    fig.suptitle(title)

    images = []
    ax = axs[0,0]
    images.append(ax.imshow(x, cmap='plasma',
                                origin='lower', extent=[0,cols*scaling[1], 0, rows*scaling[0]]))
    ax.set_title(subtitles[0])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    ax = axs[0,1]
    images.append(ax.imshow(y, cmap='plasma',
                                origin='lower', extent=[0,cols*scaling[1], 0, rows*scaling[0]]))
    ax.set_title(subtitles[1])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    ax = axs[1,0]
    images.append(ax.imshow(z, cmap='plasma',
                                origin='lower', extent=[0,cols*scaling[1], 0, rows*scaling[0]]))
    ax.set_title(subtitles[2])
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    mn = min(x.min(), y.min(), z.min())
    mx = max(x.max(), y.max(), z.max())
    norm = colors.Normalize(mn, mx)

    images[0].set_norm(norm)
    images[1].set_norm(norm)
    images[2].set_norm(norm)

    fig.delaxes(axs[1,1])

    fig.colorbar(images[0], ax=axs, orientation='horizontal', fraction=.1)
    fig.set_size_inches(11, 8.5)

    return fig

# Create a single rectangular plot with the locations of the impurities
def plot_imp_centers( datay,
                      dataz,
                      transpose=True,
                      scaling=(1.0,1.0),
                      xlabel='Y',
                      ylabel='Z',
                      title='Title'):

    rows = len(dataz[0])
    cols = len(datay[0])
    print('real, rows =', rows, 'cols =', cols)
    
    fig = plt.figure()
    plt.plot(datay[0],dataz[0],marker='.', linestyle='')
    plt.ylabel(ylabel)
    plt.xlabel(xlabel)
    fig.suptitle(title)
    fig.set_size_inches(11, 8.5)
    return fig

# Create a single 1D line plot of a property
def plot_DOS_1D(   data,
                   energy_range,
                   xlabel='Y',
                   ylabel='Z',
                   series=None,
                   title='Title'):
    m = len(data)
    n = len(data[0])
    x = energy_range[0:n]

    fig = plt.figure()
    if(m > 1):
        if (series == None):
            series = ['Series ' + str(i) for i in range(m)]
        if(len(series) < m):
            for i in range(m - len(series)):
                series.append('Series ' + str(i))
            
        for i in range(m):
            plt.plot(x,data[i], label=series[i])

        plt.legend()
    else:
        plt.plot(x, data[0])

    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    fig.suptitle(title)
    fig.set_size_inches(11, 8.5)
    return fig

# Main program unit
if __name__ == '__main__':
    # Initialize the report PDF
    pdf = PdfPages('Spin_flip_report.pdf')
    # Open the data HDF5 from spin_flip_2d
    f = h5py.File('spin_flip.h5', 'r')
    
    # Establish the domain from the HDF5
    # y_coord hold the y coordinate of all "grid points" in the domain
    # z_coord is the same but for the other z-direction
    y_coord = f['y_coord']
    z_coord = f['z_coord']
    scale = (z_coord[1], y_coord[1])

    # Delta
    # For itmax = 1, delta isn't written to HDF5. If not present, pass
    try:
        data = f['delta']['values']
        
        fig = plot_complex_rect(data, scaling=scale, title='Delta')
        pdf.savefig(fig)
    
        plt.close(fig)
    except:
        pass

    # djspin
    data = []
    data.append(f['divx'])
    data.append(f['divy'])
    data.append(f['divz'])
    fig = plot_real_xyz_rects(data, scaling=scale, title='Derivative of J spin')
    pdf.savefig(fig)

    # J[y] Charge
    data = []
    data.append(f['jytot']['values'][0,:])
    data.append(f['jyup']['values'][0,:])
    data.append(f['jydown']['values'][0,:])
    fig = plot_real_xyz_rects(data, scaling=scale, title='J Charge [Y Component]',
                                subtitles=['Total charge', 'J[y,up]', 'J[y,down]'])
    pdf.savefig(fig)
    plt.close(fig)

    # J[z] Charge
    data = []
    data.append(f['jztot']['values'][0,:])
    data.append(f['jzup']['values'][0,:])
    data.append(f['jzdown']['values'][0,:])

    fig = plot_real_xyz_rects(data, scaling=scale, title='J Charge [Z Component]',
                                subtitles=['Total charge', 'J[z,up]', 'J[z,down]'])
    pdf.savefig(fig)
    plt.close(fig)

    # J[y] Spin
    data = []
    data.append(f['jyxspin']['values'][0,:])
    data.append(f['jyyspin']['values'][0,:])
    data.append(f['jyzspin']['values'][0,:])
    fig = plot_real_xyz_rects(data, scaling=scale, title='J[Y] Spin')
    pdf.savefig(fig)
    plt.close(fig)

    # J[z] Spin
    data = []
    data.append(f['jzxspin']['values'][0,:])
    data.append(f['jzyspin']['values'][0,:])
    data.append(f['jzzspin']['values'][0,:])
    fig = plot_real_xyz_rects(data, scaling=scale, title='J[Z] Spin')
    pdf.savefig(fig)
    plt.close(fig)

    # Cooper 0, 1, 2, 3
    fig = plot_complex_rect(f['cooper0']['values'], scaling=scale, title='Cooper Pair 0 Interaction')
    pdf.savefig(fig)
    plt.close(fig)
    fig = plot_complex_rect(f['cooper1']['values'], scaling=scale, title='Cooper Pair 1 Interaction')
    pdf.savefig(fig)
    plt.close(fig)
    fig = plot_complex_rect(f['cooper2']['values'], scaling=scale, title='Cooper Pair 2 Interaction')
    pdf.savefig(fig)
    plt.close(fig)
    fig = plot_complex_rect(f['cooper3']['values'], scaling=scale, title='Cooper Pair 3 Interaction')
    pdf.savefig(fig)
    plt.close(fig)

    # Rotated Cooper 0, 1, 2
    fig = plot_complex_rect(f['cooper0_r']['values'], scaling=scale, title='Spin-Rotated Cooper Pair 0 Interaction')
    pdf.savefig(fig)
    plt.close(fig)
    fig = plot_complex_rect(f['cooper1_r']['values'], scaling=scale, title='Spin-Rotated Cooper Pair 1 Interaction')
    pdf.savefig(fig)
    plt.close(fig)
    fig = plot_complex_rect(f['cooper2_r']['values'], scaling=scale, title='Spin-Rotated Cooper Pair 2 Interaction')
    pdf.savefig(fig)
    plt.close(fig)

    # Magnitization
    mm = None
    data = []
    data.append(f['mx']['values'][0,:])
    data.append(f['my']['values'][0,:])
    data.append(f['mz']['values'][0,:])
    fig = plot_real_xyz_rects(data, scaling=scale, title='Magnetization')
    pdf.savefig(fig)
    plt.close(fig)

    # Torque
    data = []
    data.append(f['tx'])
    data.append(f['ty'])
    data.append(f['tz'])
    fig = plot_real_xyz_rects(data, scaling=scale, title='Torque')
    pdf.savefig(fig)
    plt.close(fig)

    # J Charge on Z centerline
    data = []
    data.append(f['jytot']['y_scan'][0,:])
    data.append(f['jztot']['y_scan'][0,:])
    fig = plot_real_1D( data, scale[0], 'Y', 'J Charge',
                        ['J[y]', 'J[z]'], 'J Charge along Z Centerline')
    pdf.savefig(fig)
    plt.close(fig)

    data = []
    data.append(f['jytot']['z_scan'][0,:])
    data.append(f['jztot']['z_scan'][0,:])
    fig = plot_real_1D( data, scale[1], 'Z', 'J Charge',
                        ['J[y]', 'J[z]'], 'J Charge along Y Centerline')
    pdf.savefig(fig)
    plt.close(fig)

    # J spin y scan
    data = []
    data.append(f['jyxspin']['y_scan'][0,:])
    data.append(f['jyyspin']['y_scan'][0,:])
    data.append(f['jyzspin']['y_scan'][0,:])
    data.append(f['jzxspin']['y_scan'][0,:])
    data.append(f['jzyspin']['y_scan'][0,:])
    data.append(f['jzzspin']['y_scan'][0,:])
    fig = plot_real_1D( data, scale[0], 'Y', 'J Spin',
                        ['J spin[y,x]', 'J Spin[y,y]', 'J Spin[y,z]',
                         'J spin[z,x]', 'J Spin[z,y]', 'J Spin[z,z]'],
                        'J Spin along Z Centerline')
    pdf.savefig(fig)
    plt.close(fig)

    # J spin z scan
    data = []
    data.append(f['jyxspin']['z_scan'][0,:])
    data.append(f['jyyspin']['z_scan'][0,:])
    data.append(f['jyzspin']['z_scan'][0,:])
    data.append(f['jzxspin']['z_scan'][0,:])
    data.append(f['jzyspin']['z_scan'][0,:])
    data.append(f['jzzspin']['z_scan'][0,:])
    fig = plot_real_1D( data, scale[1], 'Z', 'J Spin',
                        ['J spin[y,x]', 'J Spin[y,y]', 'J Spin[y,z]',
                         'J spin[z,x]', 'J Spin[z,y]', 'J Spin[z,z]'],
                        'J Spin along Y Centerline')
    pdf.savefig(fig)
    plt.close(fig)

    # Magnetization y scan
    data = []
    data.append(f['mx']['y_scan'][0,:])
    data.append(f['my']['y_scan'][0,:])
    data.append(f['mz']['y_scan'][0,:])
    fig = plot_real_1D( data, scale[0], 'Y', 'Magnetization',
                        ['M[x]', 'M[y]', 'M[z]'], 'Magnetization along Z Centerline')
    pdf.savefig(fig)
    plt.close(fig)

    # Magnetization z scan
    data = []
    data.append(f['mx']['z_scan'][0,:])
    data.append(f['my']['z_scan'][0,:])
    data.append(f['mz']['z_scan'][0,:])
    fig = plot_real_1D( data, scale[1], 'Y', 'Magnetization',
                        ['M[x]', 'M[y]', 'M[z]'], 'Magnetization along Y Centerline')
    pdf.savefig(fig)
    plt.close(fig)

    # Cooper Interaction y scan
    data = []
    data.append(f['cooper0']['y_scan'][0,:])
    data.append(f['cooper0']['y_scan'][1,:])
    data.append(f['cooper1']['y_scan'][0,:])
    data.append(f['cooper1']['y_scan'][1,:])
    data.append(f['cooper2']['y_scan'][0,:])
    data.append(f['cooper2']['y_scan'][1,:])
    data.append(f['cooper3']['y_scan'][0,:])
    data.append(f['cooper3']['y_scan'][1,:])
    fig = plot_real_1D( data, scale[0], 'Y', 'Interaction',
                        [   'Cooper 0 (real)','Cooper 0 (imag)',
                            'Cooper 1 (real)','Cooper 1 (imag)',
                            'Cooper 2 (real)','Cooper 2 (imag)',
                            'Cooper 3 (real)','Cooper 3 (imag)'],
                        'Cooper Interaction along Z Centerline')
    pdf.savefig(fig)
    plt.close(fig)

    # Cooper Interaction z scan
    data = []
    data.append(f['cooper0']['z_scan'][0,:])
    data.append(f['cooper0']['z_scan'][1,:])
    data.append(f['cooper1']['z_scan'][0,:])
    data.append(f['cooper1']['z_scan'][1,:])
    data.append(f['cooper2']['z_scan'][0,:])
    data.append(f['cooper2']['z_scan'][1,:])
    data.append(f['cooper3']['z_scan'][0,:])
    data.append(f['cooper3']['z_scan'][1,:])
    fig = plot_real_1D( data, scale[1], 'Z', 'Interaction',
                        [   'Cooper 0 (real)','Cooper 0 (imag)',
                            'Cooper 1 (real)','Cooper 1 (imag)',
                            'Cooper 2 (real)','Cooper 2 (imag)',
                            'Cooper 3 (real)','Cooper 3 (imag)'],
                        'Cooper Interaction along Z Centerline')
    pdf.savefig(fig)
    plt.close(fig)

    # Spin-rotated Cooper Interaction y scan
    data = []
    data.append(f['cooper0_r']['y_scan'][0,:])
    data.append(f['cooper0_r']['y_scan'][1,:])
    data.append(f['cooper1_r']['y_scan'][0,:])
    data.append(f['cooper1_r']['y_scan'][1,:])
    data.append(f['cooper2_r']['y_scan'][0,:])
    data.append(f['cooper2_r']['y_scan'][1,:])
    fig = plot_real_1D( data, scale[0], 'Y', 'Interaction',
                        [   'Cooper 0 (real)','Cooper 0 (imag)',
                            'Cooper 1 (real)','Cooper 1 (imag)',
                            'Cooper 2 (real)','Cooper 2 (imag)'],
                        'Spin-Rotated Cooper Interaction along Z Centerline')
    pdf.savefig(fig)
    plt.close(fig)

    # Cooper Interaction z scan
    data = []
    data.append(f['cooper0_r']['z_scan'][0,:])
    data.append(f['cooper0_r']['z_scan'][1,:])
    data.append(f['cooper1_r']['z_scan'][0,:])
    data.append(f['cooper1_r']['z_scan'][1,:])
    data.append(f['cooper2_r']['z_scan'][0,:])
    data.append(f['cooper2_r']['z_scan'][1,:])
    fig = plot_real_1D( data, scale[1], 'Z', 'Interaction',
                        [   'Cooper 0 (real)','Cooper 0 (imag)',
                            'Cooper 1 (real)','Cooper 1 (imag)',
                            'Cooper 2 (real)','Cooper 2 (imag)'],
                        'Spin-Rotated Cooper Interaction along Y Centerline')
    pdf.savefig(fig)
    plt.close(fig)

    # Torque y scan
    z_center = int(len(z_coord) / 2)
    data = []
    data.append(f['tx'][0:-1,z_center])
    data.append(f['ty'][0:-1,z_center])
    data.append(f['tz'][0:-1,z_center])
    data.append(f['divx'][:,z_center])
    data.append(f['divy'][:,z_center])
    data.append(f['divz'][:,z_center])

    fig = plot_real_1D( data, scale[0], 'Y', 'Torque',
                        [   'Torque [x]', 'Torque [Y]', 'Torque [Z]',
                            'Div [x]', 'Div [Y]', 'Div [Z]'],
                        'Torque along Z Centerline')
    pdf.savefig(fig)
    plt.close(fig)

    # Density of States field plots
    fields = f.keys()
    for k in fields:
        if(('DOS' in k) and ('up' in k)):
            m = re.match(r'DOS(.*)-up', k)
            eps = m.group(1)
            print(f"eps = {eps}")
            
            up = f[k]
            down = f[f'DOS{eps}-down']
            total = f[f'DOS{eps}-total']
            fig = plot_DOS(up, down, total, scale, eps)
            pdf.savefig(fig)
            plt.close(fig)            
            
    # DOS location-specific energy range 
    for k in fields:
        if ('DOS_y' in k):
            m = re.match(r'DOS_y(.*)_z(.*)',k)
            y = m.group(1)
            z = m.group(2)
            print(f'y = {y}')
            print(f'z = {z}')

            data = []
            data.append(f[k][1,:])
            data.append(f[k][2,:])
            data.append(f[k][3,:])
            energy_range = f[k][0,:]
            fig = plot_DOS_1D( data, energy_range, 'Energy', 'DOS',
                                [   'DOS Up', 'DOS Down', 'Total DOS'],
                                'DOS at y = {y}, z = {z}')
            pdf.savefig(fig)
            plt.close(fig)

    # Impurity locations
    for k in fields:
        if ('imp_loc_y' in k):
            data_y = []
            data_y.append(f['imp_loc_y'][:])
            data_z = []
            data_z.append(f['imp_loc_z'][:])
            fig = plot_imp_centers( data_y,data_z,
                                    scaling=scale,
                                    xlabel='Y',
                                    ylabel='Z',
                                    title='Impurity Locations')
            pdf.savefig(fig)
            plt.close(fig)            

    # Wmin_N
    fig = plot_real_1D([f['initial_eigen_values-min'],], 1.0, 'nKx', 'Initial Eigenvalues',
                        ['Eigen Values', ], 'Minimum Eigenvalue After Normalstate')
    pdf.savefig(fig)
    plt.close(fig)

    # Wmin_S

    fig = plot_real_1D([f['final_eigen_values-min'],], 1.0, 'nKx', 'Final Eigenvalue',
                        ['Eigen Values', ], 'Minimum Eigenvalue After Convergence')
    pdf.savefig(fig)

    pdf.close()
    plt.close(fig)
    
