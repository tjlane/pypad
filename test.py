
import tables

from autogeom import optimize as opt
from autogeom import cspad
from autogeom import score
from autogeom import utils

import numpy as np
from scipy.ndimage import filters
import matplotlib.pyplot as plt


def get_avg_from_hdf5(hd5_image_file, calibration_path, dsname, start_event,
                      n_events_to_avg=5):
    """
    Extract an average image from a psana hdf5 file.    
    """

    print 'Getting raw CSPad event %d from file %s \ndataset %s' % (event, fname, dsname)
    ds1ev = hm.getAverageCSPadEvent( fname, dsname, start_event, nevents=n_events_to_avg )
    if not ds1ev.shape == (32, 185, 388):
        print 'WARNING: ds1ev.shape =', ds1ev.shape, "should be (32, 185, 388)"

    return ds1ev
    
    
def test_metrology():
    raw_image = get_event_from_npz('../data/test_images/cxi64813_r58_evt1.npz')
    d = CSPad.from_dir('example_calibration_dir')
    x,y,z = d.coordinate_map(metrology_file="../CSPad/cspad_2011-08-10-Metrology.txt")
    
    plt.imshow(x[0,0,:,:])
    plt.show()
    
    plt.imshow(x[0,1,:,:])
    plt.show()
    
    return


def cheetah_to_psana(cheetah_image):
    """
    Takes a raw cheetah image (2D) and returns it in psana format (3D)
    """
    
    psana_image = np.zeros((32, 185, 388))
    assert cheetah_image.shape == (1480, 1552)
    
    for i in range(8):
        for j in range(4):
            x_start = 185 * i
            x_stop  = 185 * (i+1)
            y_start = 388 * j
            y_stop  = 388 * (j+1)
            psind = i + j * 8 # confirmed visually
            psana_image[psind,:,:] = cheetah_image[x_start:x_stop,y_start:y_stop]
    
    return psana_image

    
def load_AgBe():
    
    f = tables.File('data/test_images/AgBe/r0003-RawSum.h5')
    cheetah_agbe = f.root.data.data.read()
    psana_agbe = cheetah_to_psana(cheetah_agbe)
    f.close()
    
    return psana_agbe

    
def test_cheetah_conv():
    raw_image = load_AgBe()
    ai = assemble.assemble_image_from_dir(raw_image, 'example_calibration_dir')
    assemble.plot_assembled_image(ai)
    return

    
def agbe_score(real_peak_locations):
    
    # from jonas // cheetah
    path_length = 129.0148
    energy      = 9394.363725

    sref = score.PowderReference.agbe(real_peak_locations, energy, path_length)

    sref.score()

    return
    
    
def test_agbe_assembly():
    
    params_to_opt = ['offset_corr']
    
    cal_image = load_AgBe()
    init_cspad = cspad.CSPad.from_dir('data/ex_params')
    opter = opt.Optimizer(initial_cspad=init_cspad, params_to_optimize=['offset_corr'])
        
    opt_cspad, maxima = opter(cal_image, return_maxima_locations=True)
    
    agbe_score(maxima)
    
    plt.imshow( opt_cspad(cal_image).T )
    plt.show()
    
    return


def test_filter(threshold=0.025):
    
    image = load_AgBe()
    print image.shape
    
    image = np.abs(filters.sobel(image, 0)) + np.abs(filters.sobel(image, 1))
    
    for i in range(32):
        image[i,:,:] -= image[i,:,:].min()
    
    image = (image > (image.max() * threshold)).astype(np.bool)
    
    plt.imshow( utils.flatten_2x1s(image).T )
    plt.show()
    
    return
    
    
if __name__ == '__main__':
    test_agbe_assembly()
    #test_filter()
    