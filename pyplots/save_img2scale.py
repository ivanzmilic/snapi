from matplotlib import pyplot as plt
import matplotlib
import scipy.ndimage as ndimage


def save_img2scale(img, path, rng=[3, 3], resize=-1, intord=0,
                   title='No title defined', xlabel='Pixel', ylabel='Pixel',
                   cblabel='Signal [DN]', fontsz=25):
    """
    Francisco  A. Iglesias (franciscoaiglesias@hotmail.com), 2017.04.25

    Saves as an image the plot of the input numpy array keeping the original size
    and adding a color bar and pixel axes. Note that:
        -No interpolation is done by default
        -No resizeing is done by default
        -The output image is larger to fit the color bar and the axes' labels
         However, all the pixels in the input all plotted, no more no less (if resize=-1)

    INPUTS
    img(2D Numpy arr): Array to plot
    path(string): Full path to the output file. The file name defines the output image
                  format. Options are: '.png' (tested), '.eps', '.pdf', etc

    KEYWORDS
    rng([int,int]): Range of the color bar in fractions of the spatial
                    std deviaiton [min, max].
    resize(int): Set to resize the img array before plotting using nearest interpolation.
               Use -1 (default) to avoid resizeing.
    intord(int): Order of the spline used for interpolating when resizeing.
    title (string): Main title. Note that a new line showing the spatial mean
                    and stddev is alwways added below title
    xlable, ylable and cblabel (string):
        To include as x and y axes labels, and color bar label
    fontsz(int): Font size. When resize!=-1, then efective font used is fontsz*resize

    OUTPUTS
    The plot of img is saved in path
    """

    # CONSTANTS
    DPI = 96.
    MARGIN = [0.2, 0.17]  # Margins to allocate cbar and axes (not symetric)
    SCIFMT = '{:6.2e}'  # format string for mean and stddev print

    # MAIN
    if resize != -1:
        im_data = ndimage.zoom(img, resize, order=intord)
        fontsz *= resize
    else:
        im_data = img

    height, width = im_data.shape

    # figure size in inches
    figsize = width / float(DPI) / (1. - MARGIN[0]), height / float(DPI) / (1. - MARGIN[1])

    # img mean and stdddev
    mimg = im_data.mean()
    stdimg = im_data.std()

    # Create a figure of the right size with one axis that takes up the full figure
    fig = plt.figure(figsize=figsize)
    matplotlib.rcParams.update({'font.size': fontsz})
    ax = fig.add_axes([MARGIN[0] / 3., MARGIN[1] / 3., (1. - MARGIN[0]), (1. - MARGIN[1])])

    # Display the image.
    im = ax.imshow(im_data, cmap='Greys_r', vmin=mimg - rng[0] * stdimg,
                   vmax=mimg + rng[1] * stdimg, interpolation='none')

    # Add axis for the color bar and plot it
    box = ax.get_position()
    cax = plt.axes([box.x0 + box.width * 1.01, box.y0, MARGIN[0] / 8., box.height])
    cbar = plt.colorbar(im, cax=cax, orientation="vertical")
    cbar.set_label(cblabel)
    ax.set_title(title + '\n mean=' + SCIFMT.format(mimg) + '; std=' + SCIFMT.format(stdimg))
    ax.set_xlabel(ylabel)
    ax.set_ylabel(xlabel)

    # Ensure we're displaying with square pixels and the right extent.
    # ax.set(xlim=[0, width], ylim=[height, 0], aspect=1)

    plt.savefig(path, dpi=DPI)
    plt.close()
