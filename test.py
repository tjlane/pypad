
import tables

import optimize as opt
import cspad

import numpy as np
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


def get_event_from_npz(npz_image_file):
    return np.load(npz_image_file)['arr_0']
    
    
def show_assembled_image(image):
    plt.imshow(image.T)
    plt.show()
    return

    
def pyana_assembly(raw_image, calibration_path, run_number=0):
    
    print 'Loading calibration parameters from: %s' % calibration_path
    calp.calibpars.setCalibParsForPath( run=run_number, path=calibration_path )
    calp.calibpars.printCalibPars()
    cpe.cpeval.printCalibParsEvaluatedAll()
    
    print 'Constructing the CSPad image from raw array'
    cspadimg = cip.CSPadImageProducer(rotation=0, tiltIsOn=True )
    image = cspadimg.getCSPadImage( raw_image )
    
    return image


def test_assembly_from_dir():
    
    # this one
    raw_image = get_event_from_npz('../test_data/cxi64813_r58_evt1.npz')
    d = CSPad.from_dir('../ex_params')
    show_assembled_image( d(raw_image) )
    
    # should be the same as this one
    # ai = pyana_assembly(raw_image, 'example_calibration_dir')
    # show_assembled_image(ai)
    
    return
    
    
def test_metrology():
    raw_image = get_event_from_npz('../test_data/cxi64813_r58_evt1.npz')
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
    
    f = tables.File('../test_data/AgBe/r0003-RawSum.h5')
    cheetah_agbe = f.root.data.data.read()
    psana_agbe = cheetah_to_psana(cheetah_agbe)
    
    return psana_agbe

    
def test_cheetah_conv():
    raw_image = load_AgBe()
    ai = assemble.assemble_image_from_dir(raw_image, 'example_calibration_dir')
    assemble.plot_assembled_image(ai)
    return

    
def test_agbe_assembly():
    
    params_to_opt = ['offset_corr']
    
    cal_image = load_AgBe()
    init_cspad = cspad.CSPad.from_dir('../ex_params')
    opter = opt.Optimizer(initial_cspad=init_cspad, params_to_optimize=['offset_corr'])
    
    opt_cspad = opter(cal_image)
    plt.imshow( opt_cspad(cal_image).T )
    
    return

    
if __name__ == '__main__':
    test_agbe_assembly()
