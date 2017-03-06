#ifndef __CONF_H__ // __CONF_H__
#define __CONF_H__

#define CFG_RECONSTR           0     // reconstruction mode
#define CFG_CALIBRATE          1     // calibration mode
#define CFG_FF_ONLY            2     // only flatfield and clip, then exit

#define CFG_GRADIENT_DIFF      1     // finite difference
#define CFG_GRADIENT_VOGEL     2     // Vogel

#define CFG_GETSTEP_SDSC       1     // steepest descent
#define CFG_GETSTEP_CNJG       2     // conjugate gradient
#define CFG_GETSTEP_BFGS_inv   3     // inverse BFGS?
#define CFG_GETSTEP_BFGS       4     // BFGS?
#define CFG_GETSTEP_NAMES {"unknown","steepest descent","conjugate_gradient","inverse BFGS","unknown"}

#define CFG_FPM_MEDIAN         1
#define CFG_FPM_INVDISTWEIGHT  2
#define CFG_FPMETHOD_NAMES {"unknown","median","inverse distance weighted","unknown"}

#define CFG_DEFAULT           0
#define CFG_ZERNIKE           1
#define CFG_KARHUNEN_LOEVE    2

#define CFG_FT_NONE           0
#define CFG_FT_ANA            1
#define CFG_FT_FITS           2
#define CFG_FT_OMASK         (CFG_FT_ANA|CFG_FT_FITS)
#define CFG_FT_MOMFBD         4
#define CFG_FT_NMASK         (CFG_FT_MOMFBD)
#define CFG_FT_NAMES         {"NONE","ANA","FITS","MOMFBD"}
#define CFG_FT_EXT           {"NONE","f0","fits","momfbd"}

#define CFG_PSF_ALL           1
#define CFG_PSF_AVG           2

#endif             // __CONF_H__
